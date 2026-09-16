from pathlib import Path
import argparse,json,hashlib,datetime,sys
import numpy as np,pandas as pd,scipy
from scipy.stats import spearmanr,rankdata
P=Path(__file__).resolve().parent;I=P/'inputs'
a=argparse.ArgumentParser();a.add_argument('--output',type=Path,default=P/'results');O=a.parse_args().output;O.mkdir(parents=True,exist_ok=True)
assert not (O/'RUN_RECEIPT.json').exists()
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
assert sha(P/'PLAN.md')==json.loads((P/'FREEZE_RECEIPT.json').read_text())['plan_sha256']
for x in json.loads((I/'INPUT_MANIFEST.json').read_text())['files']:assert sha(I/x['name'])==x['sha256']
meta=pd.read_csv(I/'metadata.csv');oldM=pd.read_csv(I/'prior_M5_states.csv');oldG=pd.read_csv(I/'prior_members.csv')
states=sorted(set(oldM.state)-{'ALL_CELLS'});assert len(states)==11
targets=['PRG4','OGN','M5'];mgenes=['COL2A1','SOX9','ACAN','CHAD','HAPLN1'];cache={};context=[];gates=[];checks=[]
for m in meta.itertuples():
 d=np.load(I/(m.GSM+'.npz'));names=d['genes'].tolist();E=d['expression'];assert np.allclose(E,np.log1p(1e4*d['counts']/d['full_UMI'][:,None]),atol=1e-12,rtol=0)
 X=np.column_stack([E[:,names.index(g)] for g in ['PIEZO1','TRPV4','PRG4','OGN']]+[E[:,[names.index(g) for g in mgenes]].mean(1)])
 cache[m.GSM]=(X,d['states'])
 for state in ['ALL_CELLS']+states:
  mask=np.ones(len(X),bool) if state=='ALL_CELLS' else d['states']==state;Y=X[mask];n=len(Y)
  for j,g in enumerate(['PIEZO1','TRPV4']+targets):
   context.append(dict(GSM=m.GSM,donor=m.donor_id,region=m.region,state=state,gene=g,cells=n,detected=int((Y[:,j]>0).sum()),fraction=float((Y[:,j]>0).mean()) if n else np.nan,mean_including_zero=float(Y[:,j].mean()) if n else np.nan))
  for j,g in enumerate(targets,2):
   nd=(Y[:,[0,1,j]]>0).sum(0);bad=[]
   if n<200:bad.append('CELLS_LT_200')
   for label,v in zip(['PIEZO1','TRPV4',g],nd):
    if v<20:bad.append(label+'_DETECTED_LT_20')
   if not n or np.any(np.var(Y[:,[0,1,j]],axis=0)==0):bad.append('VARIANCE_GATE')
   gates.append(dict(GSM=m.GSM,donor=m.donor_id,region=m.region,state=state,target=g,cells=n,PIEZO1_detected=int(nd[0]),TRPV4_detected=int(nd[1]),target_detected=int(nd[2]),gate='PASS' if not bad else ';'.join(bad)))
gate=pd.DataFrame(gates);gate.to_csv(O/'state_region_admission.csv',index=False);pd.DataFrame(context).to_csv(O/'state_expression_context.csv',index=False)
# Admission file is written before new correlation effects.
rows=[]
for row in gate.to_dict('records'):
 row.update(PIEZO1_z=np.nan,TRPV4_z=np.nan,delta=np.nan,status=row['gate'])
 if row['gate']=='PASS':
  X,labels=cache[row['GSM']];Y=X if row['state']=='ALL_CELLS' else X[labels==row['state']];j=targets.index(row['target'])+2
  r=np.array([spearmanr(Y[:,k],Y[:,j]).statistic for k in [0,1]]);ind=np.corrcoef(rankdata(Y[:,[0,1,j]],axis=0),rowvar=False)[:2,2]
  err=float(max(abs(r-ind)));assert err<1e-10;checks.append(dict(GSM=row['GSM'],state=row['state'],target=row['target'],error=err))
  if np.isfinite(r).all() and (abs(r)<1).all():
   z=np.arctanh(r);row.update(PIEZO1_z=z[0],TRPV4_z=z[1],delta=z[0]-z[1],status='ESTIMABLE')
  else:row['status']='CORRELATION_GATE'
 rows.append(row)
rg=pd.DataFrame(rows);rg.to_csv(O/'state_region_estimates.csv',index=False)
replay=[]
for row in rg.itertuples():
 old=None
 if row.target=='M5':
  old=oldM[(oldM.GSM==row.GSM)&(oldM.state==row.state)].iloc[0];assert (row.status=='ESTIMABLE')==(old.status=='ESTIMABLE');oc='contrast_z'
 elif row.state=='ALL_CELLS':
  legacy=meta.set_index('GSM').loc[row.GSM,'legacy_id'];old=oldG[(oldG['sample']==legacy)&(oldG.member==row.target)&(oldG.scenario=='primary')].iloc[0];oc='delta'
 if old is not None and row.status=='ESTIMABLE':
  err=max(abs(row.PIEZO1_z-old.PIEZO1_z),abs(row.TRPV4_z-old.TRPV4_z),abs(row.delta-old[oc]));assert err<1e-10;replay.append(dict(GSM=row.GSM,state=row.state,target=row.target,error=err))
pd.DataFrame(replay).to_csv(O/'historical_replay.csv',index=False);pd.DataFrame(checks).to_csv(O/'rank_checks.csv',index=False)
dn=[]
for (state,donor,target),x in rg.groupby(['state','donor','target']):
 assert len(x)==2 and set(x.region)=={'WB','NWB'};ok=x.status.eq('ESTIMABLE').all()
 row=dict(state=state,donor=donor,target=target,eligible=ok,failures='|'.join(x.loc[x.status.ne('ESTIMABLE'),'region']+':'+x.loc[x.status.ne('ESTIMABLE'),'status']))
 for k in ['PIEZO1_z','TRPV4_z','delta']:row[k]=float(x[k].mean()) if ok else np.nan
 dn.append(row)
D=pd.DataFrame(dn);D.to_csv(O/'donor_state_estimates.csv',index=False)
support=[];matched=[];summary=[]
def summ(v):
 v=np.asarray(v,float);n=len(v);return dict(n=n,mean=float(v.mean()) if n else np.nan,sd=float(v.std(ddof=1)) if n>1 else np.nan,median=float(np.median(v)) if n else np.nan,min=float(v.min()) if n else np.nan,max=float(v.max()) if n else np.nan,positive=int((v>0).sum()),negative=int((v<0).sum()))
for state in states:
 sets={g:set(D[(D.state==state)&(D.target==g)&D.eligible].donor) for g in targets};common=sorted(set.intersection(*sets.values()));ok=len(common)>=6
 support.append(dict(state=state,**{g+'_eligible_donors':len(sets[g]) for g in targets},common_donors=';'.join(common),n_common=len(common),supported=ok))
 for g in targets:
  for donor in common:
   now=D[(D.state==state)&(D.target==g)&(D.donor==donor)].iloc[0];base=D[(D.state=='ALL_CELLS')&(D.target==g)&(D.donor==donor)].iloc[0];assert base.eligible
   for k in ['PIEZO1_z','TRPV4_z','delta']:matched.append(dict(state=state,target=g,donor=donor,quantity=k,supported=ok,within_state=now[k],matched_overall=base[k],paired_change=now[k]-base[k]))
pd.DataFrame(support).to_csv(O/'state_common_admission.csv',index=False)
MD=pd.DataFrame(matched);MD.to_csv(O/'matched_donor_comparisons.csv',index=False)
for (state,g,k),x in MD.groupby(['state','target','quantity']):
 for what in ['within_state','matched_overall','paired_change']:summary.append(dict(state=state,target=g,quantity=k,estimate=what,supported=bool(x.supported.iloc[0]),donors=';'.join(x.donor),**summ(x[what])))
S=pd.DataFrame(summary);S.to_csv(O/'matched_descriptive_summaries.csv',index=False)
receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),plan_sha256=sha(P/'PLAN.md'),code_sha256=sha(__file__),input_manifest_sha256=sha(I/'INPUT_MANIFEST.json'),rank_max_error=max(x['error'] for x in checks),historical_max_error=max(x['error'] for x in replay),state_regions=len(gate),supported_states=[x['state'] for x in support if x['supported']],new_tests=0,manuscript_modified=False,python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
(O/'RUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(pd.DataFrame(support).to_string(index=False));print(S.query("supported and quantity=='delta' and estimate=='within_state'").to_string(index=False))
