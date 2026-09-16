from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,hashlib,datetime
import numpy as np,pandas as pd
from scipy import stats
import scipy
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parents[1]
B=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P/'replay_output'
B.mkdir(parents=True,exist_ok=True)
import shutil
shutil.copy2(P/'V34_PLAN.md',B/'PLAN.md')
R=P/'legacy'
M='COL2A1 SOX9 ACAN CHAD HAPLN1'.split(); H='COL10A1 RUNX2 IBSP ALPL MMP13 SPP1 POSTN'.split(); G=['PIEZO1','TRPV4']+M+H
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
guards=[R/'audit/p15_components_v2/core_tests_BH2.tsv',R/'audit/p21_external_effects/primary_BH4_family.tsv',R/'audit/p21_external_effects/donor_effects.tsv']
receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),intent='POST_HOC_DESCRIPTIVE',plan_sha256=sha(B/'PLAN.md'),code_sha256=sha(__file__),guards={str(p):sha(p) for p in guards})
assert not (B/'execution_receipt.json').exists(),'Do not overwrite an existing execution'
(B/'execution_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
inputs=[]; datasets=[]
def read(p):
 inputs.append(dict(path=str(p),sha256=sha(p)));return pd.read_csv(p,sep='\t')
meta=read(R/'audit/p07_prefreeze/candidate_design_16x9.tsv')
for m in meta.itertuples():
 x=read(R/f'audit/p07_execution/E3_cell_input_{m.legacy_id}.tsv.gz');c=read(R/f'audit/p09_exploration/channel_counts_{m.legacy_id}.tsv.gz')
 assert x.cell_id.equals(c.cell_id) and x.cell_id.is_unique
 counts=np.column_stack([c[['PIEZO1','TRPV4']].to_numpy(),x[M+H].to_numpy()]);full=x.total_counts.to_numpy();expr=np.log1p(1e4*counts/full[:,None])
 datasets.append(dict(cohort='discovery',donor=str(m.donor_id),sample=m.GSM,region=m.region,ids=x.cell_id.to_numpy(),counts=counts,full=full,expr=expr))
cross=read(R/'audit/p19_admission/donor_library_run_crosswalk.tsv');cross=cross[cross.GSE.eq('GSE220243')]
for m in cross.itertuples():
 x=read(R/f'audit/p21_external_effects/{m.GSM}_normalized_targets_and_masks.tsv.gz');x=x.loc[x.primary_P19].reset_index(drop=True)
 assert x.barcode.is_unique
 expr=x[G].to_numpy();full=x.full_UMI.to_numpy();recon=np.expm1(expr)*full[:,None]/1e4;counts=np.rint(recon)
 assert np.max(abs(recon-counts))<1e-6
 datasets.append(dict(cohort='external',donor=str(m.donor_key),sample=m.GSM,region='sample',ids=x.barcode.to_numpy(),counts=counts,full=full,expr=expr))
assert sum(len(d['full']) for d in datasets if d['cohort']=='discovery')==101702
assert sum(len(d['full']) for d in datasets if d['cohort']=='external')==35389
(B/'input_hashes.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')
rows=[];parts=[];des=[];errors=[]
def corr(x,y):
 if len(x)<2 or np.var(x)==0 or np.var(y)==0:return np.nan
 return float(np.corrcoef(x,y)[0,1])
def fit(d,mode,expr,full,draw=-1):
 n=len(full);score=expr[:,2:7].mean(1);ch=expr[:,:2];det=(ch>0).sum(0)
 r=[np.nan,np.nan];status='PASS';err=np.nan;rank=0
 if n<200 or min(det)<20 or not np.isfinite(expr).all() or np.any(full<=0):status='CELL_DETECTION_OR_FINITE_GATE'
 else:
  ranks=np.column_stack([stats.rankdata(ch[:,0]),stats.rankdata(ch[:,1]),stats.rankdata(score)])
  if 'depth_cubic' in mode:
   z=np.log(full);z=(z-z.mean())/z.std();X=np.column_stack([np.ones(n),z,z*z,z*z*z]);rank=np.linalg.matrix_rank(X)
   if rank<4:status='RANK_DEFICIENT'
   else:
    resid=ranks-X@np.linalg.lstsq(X,ranks,rcond=None)[0];Q=np.linalg.qr(X,mode='reduced')[0];other=ranks-Q@(Q.T@ranks)
    r=[corr(resid[:,j],resid[:,2]) for j in [0,1]];err=max(abs(r[j]-corr(other[:,j],other[:,2])) for j in [0,1]);errors.append(err);assert err<1e-10
  else:r=[corr(ranks[:,j],ranks[:,2]) for j in [0,1]]
  if not all(np.isfinite(a) and abs(a)<1 for a in r):status='VARIANCE_OR_BOUNDARY_GATE'
 row={k:d[k] for k in ['cohort','donor','sample','region']};row.update(mode=mode,draw=draw,n_cells=n,PIEZO1_detected=int(det[0]),TRPV4_detected=int(det[1]),status=status,PIEZO1_rho=r[0],TRPV4_rho=r[1],contrast_z=float(np.arctanh(r[0])-np.arctanh(r[1])) if status=='PASS' else np.nan,rank=rank,qr_error=err)
 rows.append(row)
for di,d in enumerate(datasets):
 key={k:d[k] for k in ['cohort','donor','sample','region']};expr=d['expr'];ct=d['counts'];full=d['full'];n=len(full)
 assert np.isfinite(ct).all() and (ct>=0).all() and (ct==np.floor(ct)).all() and (full==np.floor(full)).all() and (full>=ct.sum(1)).all()
 score=expr[:,2:7].mean(1); V=stats.rankdata(score)
 des.append(dict(**key,n_cells=n,depth_min=float(full.min()),depth_median=float(np.median(full)),depth_mean=float(full.mean()),depth_sd=float(full.std(ddof=1)),depth_q10=float(np.quantile(full,.1)),joint_positive=int(((ct[:,:2]>0).all(1)).sum()),depth1000_cells=int((full>=1000).sum())))
 for j,ch in enumerate(G[:2]):
  y=expr[:,j];U=stats.rankdata(y);pos=ct[:,j]>0;p=pos.mean();den=U.std()*V.std();cov=lambda a,b:float(np.mean((a-a.mean())*(b-b.mean())))
  total=cov(U,V)/den;positive=p*cov(U[pos],V[pos])/den if pos.any() else 0.;detection=p*(1-p)*(U[pos].mean()-U[~pos].mean())*(V[pos].mean()-V[~pos].mean())/den if pos.any() and (~pos).any() else 0.
  assert abs(total-positive-detection)<1e-10
  parts.append(dict(**key,channel=ch,cells=n,detected=int(pos.sum()),detection_fraction=p,total_rho=total,positive_contribution=positive,detection_contribution=detection,identity_error=abs(total-positive-detection),detection_M_rho=corr(pos.astype(float),V),positive_only_reranked_rho=float(stats.spearmanr(y[pos],score[pos]).statistic) if pos.sum()>=20 else np.nan))
 fit(d,'all_cells',expr,full);fit(d,'all_cells_depth_cubic',expr,full)
 joint=(ct[:,:2]>0).all(1);fit(d,'joint_positive',expr[joint],full[joint]);fit(d,'joint_positive_depth_cubic',expr[joint],full[joint])
 eligible=full>=1000;c=ct[eligible].astype('int64');f=full[eligible].astype('int64');fit(d,'depth1000_unthinned',expr[eligible],f)
 remain=f-c.sum(1);prob=1000/f;rng=np.random.default_rng(np.random.SeedSequence([20260912,di]))
 for rep in range(20):
  thin=rng.binomial(c,prob[:,None]);total=thin.sum(1)+rng.binomial(remain,prob)
  val=np.log1p(1e4*thin/total[:,None]);fit(d,'depth1000_thinned',val,total,rep)
 print(d['cohort'],d['donor'],d['region'],n,'done',flush=True)
sample=pd.DataFrame(rows);sample.to_csv(B/'sample_estimates.csv',index=False);pd.DataFrame(parts).to_csv(B/'decomposition_sample.csv',index=False);pd.DataFrame(des).to_csv(B/'sample_descriptives.csv',index=False)
donors=[]
for (co,dn,mo,rep),v in sample.groupby(['cohort','donor','mode','draw'],sort=False):
 need=2 if co=='discovery' else 1;ok=len(v)==need and v.status.eq('PASS').all()
 donors.append(dict(cohort=co,donor=dn,mode=mo,draw=rep,status='PASS' if ok else 'INCOMPLETE_SAMPLE_SET',contrast_z=v.contrast_z.mean() if ok else np.nan))
don=pd.DataFrame(donors);don.to_csv(B/'donor_estimates_all_draws.csv',index=False)
thin=[]
for (co,dn),v in don[don['mode'].eq('depth1000_thinned')].groupby(['cohort','donor']):
 ok=len(v)==20 and v.status.eq('PASS').all();thin.append(dict(cohort=co,donor=dn,mode='depth1000_thinned_mean20',draw=-1,status='PASS' if ok else 'INCOMPLETE_DRAWS',contrast_z=v.contrast_z.mean() if ok else np.nan,valid_draws=int(v.status.eq('PASS').sum()),MC_min=v.contrast_z.min(),MC_max=v.contrast_z.max(),MC_sd=v.contrast_z.std()))
summdon=pd.concat([don[don.draw==-1],pd.DataFrame(thin)],ignore_index=True);summdon.to_csv(B/'donor_estimates.csv',index=False)
# Reproduce frozen main estimates before interpreting any new output.
base=summdon[summdon['mode'].eq('all_cells')]
assert abs(base.loc[base.cohort=='discovery','contrast_z'].mean()-(-.2939357581))<1e-8
old=pd.read_csv(R/'audit/p21_external_effects/donor_effects.tsv',sep='\t');old=old[old.selection.eq('primary_P19')]
got=base[base.cohort=='external'].set_index('donor').contrast_z;expected=old.set_index('donor_key').DeltaM
assert np.max(abs(got.sort_index()-expected.sort_index()))<1e-10
# Optional historical Shapiro diagnostics are archived; no test selection occurs here.
def summary(values,seed):
 a=np.asarray(values,dtype=float);a=a[np.isfinite(a)];n=len(a)
 if not n:return dict(n=0,mean=None)
 boot=a[np.random.default_rng(seed).integers(0,n,(10000,n))].mean(1)
 return dict(n=n,mean=float(a.mean()),median=float(np.median(a)),sd=float(a.std(ddof=1)) if n>1 else np.nan,min=float(a.min()),max=float(a.max()),negative=int((a<0).sum()),positive=int((a>0).sum()),bootstrap95_low=float(np.quantile(boot,.025)) if n>1 else np.nan,bootstrap95_high=float(np.quantile(boot,.975)) if n>1 else np.nan)
ss=[dict(cohort=co,mode=mo,**summary(v.contrast_z,20260912)) for (co,mo),v in summdon.groupby(['cohort','mode'])]
pd.DataFrame(ss).to_csv(B/'summary.csv',index=False)
decp=pd.DataFrame(parts);dc=[]
for (co,dn),v in decp.groupby(['cohort','donor']):
 a=v[v.channel=='PIEZO1'];b=v[v.channel=='TRPV4']
 for term in ['total_rho','positive_contribution','detection_contribution']:
  dc.append(dict(cohort=co,donor=dn,term=term,value=a[term].mean()-b[term].mean()))
dc=pd.DataFrame(dc);dc.to_csv(B/'decomposition_donor.csv',index=False)
pd.DataFrame([dict(cohort=co,term=t,**summary(v.value,20260912)) for (co,t),v in dc.groupby(['cohort','term'])]).to_csv(B/'decomposition_summary.csv',index=False)
assert all(sha(p)==h for p,h in receipt['guards'].items())
(B/'validation.json').write_text(json.dumps(dict(status='COMPLETED',new_hypothesis_tests=0,old_results_unchanged=True,max_qr_error=max(errors),baseline_external_max_error=float(np.max(abs(got.sort_index()-expected.sort_index()))),software=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)),indent=2),encoding='utf-8')
print(pd.DataFrame(ss).to_string(index=False))
