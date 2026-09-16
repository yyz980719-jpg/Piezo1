from pathlib import Path
import json,hashlib,datetime,gzip,itertools,sys
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];A=R/'audit/p21_external_effects';A.mkdir(exist_ok=True)
plan=R/'protocol/P21_EXTERNAL_EFFECTS.md'
receipt=A/'plan_receipt.json'
if not receipt.exists():receipt.write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(plan.read_bytes()).hexdigest(),'external_effects_seen_in_workflow':False},indent=2))
else:assert json.loads(receipt.read_text())['sha256']==hashlib.sha256(plan.read_bytes()).hexdigest()
M='COL2A1 SOX9 ACAN CHAD HAPLN1'.split();H='COL10A1 RUNX2 IBSP ALPL MMP13 SPP1 POSTN'.split();genes=['PIEZO1','TRPV4']+M+H
cross=pd.read_csv(R/'audit/p19_admission/donor_library_run_crosswalk.tsv',sep='\t');cohort=cross[cross.GSE.eq('GSE220243')];assert len(cohort)==6 and cohort.donor_key.is_unique and cohort.mapping_status.eq('SOURCE_REPORTED_DONOR').all()
manifest={x['GSM']:x for x in json.loads((R/'audit/p18_external/matrix_download_manifest.json').read_text())}
rows=[];sources=[]
for donor in cohort.itertuples():
 gsm=donor.GSM;print('READ',gsm,flush=True)
 f=next((R/'audit/p17_external').glob(gsm+'*features.tsv.gz'));g=pd.read_csv(f,sep='\t',header=None);symbols=g[1].tolist()
 if gsm=='GSM6797154':
  record=json.loads((R/'audit/p20_qc/download_receipt.json').read_text());path=Path(record['path']);d=pd.read_csv(R/'audit/p20_qc/CartOA1_screened_cells.tsv.gz',sep='\t');ix=d.source_index.to_numpy()-1
 else:
  record=manifest[gsm];path=Path(record['path']);all_d=pd.read_csv(R/f'audit/p19_admission/{gsm}_cell_admission.tsv.gz',sep='\t');ix=np.flatnonzero(all_d.QC_pass);d=all_d.loc[all_d.QC_pass].reset_index(drop=True)
 assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
 bc=gzip.open(next((R/'audit/p17_external').glob(gsm+'*barcodes.tsv.gz')),'rt').read().splitlines();assert np.array_equal(np.array(bc)[ix],d.barcode.to_numpy())
 x=mmread(path).tocsc();assert x.shape==(len(g),len(bc));x=x[:,ix];assert np.all(x.data>=0)
 full=np.asarray(x[g[2].eq('Gene Expression').to_numpy()].sum(axis=0)).ravel();assert np.array_equal(full,d.total_UMI.to_numpy()) and np.all(full>0)
 assert all(symbols.count(v)==1 and g.iloc[symbols.index(v),2]=='Gene Expression' for v in genes)
 counts=x[[symbols.index(v) for v in genes],:].toarray().T;expr=np.log1p(1e4*counts/full[:,None]);del x
 calls=pd.read_csv(R/f'audit/p20_qc/{gsm}_crosscheck_cell_decisions.tsv.gz',sep='\t');assert calls.barcode.is_unique
 calls=calls.set_index('barcode').loc[d.barcode].reset_index();assert np.array_equal(calls.barcode,d.barcode)
 neg=calls.negative_lineage.to_numpy();scr=calls.Scrublet_doublet.to_numpy();sdb=calls.scDblFinder_class.eq('doublet').to_numpy()
 masks={'primary_P19':~scr&~neg,'QC_only':np.ones(len(d),bool),'scDblFinder':~sdb&~neg,'union':~scr&~sdb&~neg}
 cell=pd.DataFrame(expr,columns=genes);cell.insert(0,'barcode',d.barcode.to_numpy());cell['full_UMI']=full
 for name,mask in masks.items():cell[name]=mask
 cell.to_csv(A/(gsm+'_normalized_targets_and_masks.tsv.gz'),sep='\t',index=False)
 for name,mask in masks.items():
  v=expr[mask];ct=counts[mask];n=len(v);base={'GSM':gsm,'donor_key':donor.donor_key,'selection':name,'n_cells':n,'P1_detected':int((ct[:,0]>0).sum()),'T4_detected':int((ct[:,1]>0).sum())}
  vals={}
  for program,js in [('M',range(2,7)),('H',range(7,14))]:
   score=v[:,list(js)].mean(axis=1);eligible=n>=200 and base['P1_detected']>=20 and base['T4_detected']>=20 and np.var(score)>0 and np.var(v[:,0])>0 and np.var(v[:,1])>0
   rhos=[float(stats.spearmanr(v[:,k],score).statistic) for k in range(2)] if eligible else [np.nan,np.nan]
   eligible=eligible and all(np.isfinite(z) and abs(z)<1 for z in rhos)
   vals[program]=np.arctanh(rhos[0])-np.arctanh(rhos[1]) if eligible else np.nan
   base.update({f'{program}_eligible':eligible,f'{program}_score_variance':float(np.var(score)),f'P1_{program}_rho':rhos[0],f'T4_{program}_rho':rhos[1]})
  base.update(DeltaM=vals['M'],DeltaH=vals['H'],Gamma=vals['M']-vals['H']);rows.append(base)
 sources.append({'GSM':gsm,'matrix_sha256':record['sha256'],'matrix_path':str(path),'feature_sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
out=pd.DataFrame(rows);out.to_csv(A/'donor_effects.tsv',sep='\t',index=False);(A/'input_manifest.json').write_text(json.dumps(sources,indent=2))
# Diagnostic points and normality checks precede inferential summaries; no data deletion or test switching.
sys.path.insert(0,'C:/Users/Administrator/.codex/skills/statistical-analysis/scripts')
from assumption_checks import check_normality
primary=out[out.selection.eq('primary_P19')];diagnostics=[]
fig,axs=plt.subplots(1,2,figsize=(9,4))
for ax,endpoint in zip(axs,['DeltaM','Gamma']):
 v=primary[endpoint].dropna().to_numpy();ax.scatter(np.arange(1,len(v)+1),v,color='black');ax.axhline(0,color='gray',linestyle='--');ax.set(xlabel='Source donor index',ylabel=endpoint,title=f'{endpoint}: n={len(v)} donors')
 if len(v)>=3:
  z=check_normality(v,name=endpoint,plot=False);diagnostics.append({'endpoint':endpoint,'result':z,'limitation':'n<=6; normality/non-rejection does not establish symmetry or independence; no test switch'})
fig.tight_layout();fig.savefig(A/'donor_points_diagnostic.png',dpi=160);plt.close(fig)
(A/'assumption_diagnostics.json').write_text(json.dumps(diagnostics,indent=2,default=lambda z:z.item() if isinstance(z,np.generic) else str(z)))
def flip(v):
 v=np.asarray(v);null=np.array([np.mean(v*np.array(s)) for s in itertools.product([-1,1],repeat=len(v))]);return float(np.mean(np.abs(null)>=abs(v.mean())-1e-14))
assert flip([1,1,1])==.25 and flip([1,-1])==1 and flip([1]*6)==.03125
summaries=[]
for selection,grp in out.groupby('selection',sort=False):
 for endpoint in ['DeltaM','Gamma']:
  v=grp[endpoint].dropna().to_numpy();n=len(v);mean=float(v.mean()) if n else np.nan;sd=float(v.std(ddof=1)) if n>1 else np.nan;half=stats.t.ppf(.975,n-1)*sd/np.sqrt(n) if n>=3 else np.nan
  summaries.append({'cohort':'GSE220243','selection':selection,'endpoint':endpoint,'n_donors':n,'mean':mean,'sd':sd,'median':float(np.median(v)) if n else np.nan,'negative_donors':int((v<0).sum()),'ci95_low':mean-half,'ci95_high':mean+half,'df':n-1,'p_signflip':flip(v) if selection=='primary_P19' and n>=3 else np.nan,'intent':'primary' if selection=='primary_P19' else 'descriptive sensitivity'})
summary=pd.DataFrame(summaries);summary.to_csv(A/'cohort_summary.tsv',sep='\t',index=False)
fam=[]
for gse in ['GSE220243','GSE169454']:
 for ep in ['DeltaM','Gamma']:
  val=summary.loc[summary.selection.eq('primary_P19')&summary.endpoint.eq(ep),'p_signflip'].iloc[0] if gse=='GSE220243' else np.nan
  fam.append({'cohort':gse,'endpoint':ep,'observed_p':val,'p_for_BH4':val if np.isfinite(val) else 1.,'status':'TESTED' if np.isfinite(val) else 'NOT_AVAILABLE_PATIENT_MAPPING' if gse=='GSE169454' else 'INSUFFICIENT_ELIGIBLE_DONORS'})
fam=pd.DataFrame(fam);fam['BH4_accounting_q']=stats.false_discovery_control(fam.p_for_BH4,method='bh');fam['reported_q']=np.where(fam.observed_p.notna(),fam.BH4_accounting_q,np.nan);fam.to_csv(A/'primary_BH4_family.tsv',sep='\t',index=False)
print(summary.to_string(index=False));print(fam.to_string(index=False))
