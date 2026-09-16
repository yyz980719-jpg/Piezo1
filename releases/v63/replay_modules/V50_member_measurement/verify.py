from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata,t
P=Path(__file__).resolve().parent; R=P/'results'; errs=[]
s=pd.read_csv(R/'member_sample_scenarios.csv'); d=pd.read_csv(R/'member_donor_scenarios.csv')
assert len(s)==800 and len(d)==544 and s.status.eq('ESTIMABLE').all()
for sample,rows in s.groupby('sample'):
 z=np.load(P/'inputs'/(sample+'.npz')); names=z['genes'].tolist(); E=z['expression'].copy()
 if sample.startswith('GSM'):
  for j,g in enumerate(z['archived_genes'].tolist()): E[:,names.index(g)]=z['archived_expression'][:,j]
 for sc,ss in rows.groupby('scenario'):
  mask=np.ones(len(E),bool) if sc!='union' else pd.read_csv(P/'inputs'/(sample+'_fixed_masks.csv'))['union'].to_numpy()
  X=rankdata(E[mask],axis=0)
  if sc=='primary_cubic_depth':
   a=np.log(z['full_UMI'][mask]);a=(a-a.mean())/a.std(); Q=np.linalg.qr(np.column_stack([a**k for k in range(4)]))[0];X=X-Q@(Q.T@X)
  C=np.corrcoef(X,rowvar=False)
  for row in ss.itertuples():
   for ch in ['PIEZO1','TRPV4']:errs.append(abs(C[names.index(ch),names.index(row.member)]-getattr(row,ch+'_rho')))
for row in d.itertuples():
 ss=s[(s.cohort==row.cohort)&(s.donor==row.donor)&(s.member==row.member)&(s.scenario==row.scenario)]
 assert len(ss)==(2 if row.cohort=='discovery' else 1)
 for k in ['PIEZO1_z','TRPV4_z','delta']:errs.append(abs(ss[k].sum()/len(ss)-getattr(row,k)))
for row in pd.read_csv(R/'source_scenario_summaries.csv').itertuples():
 v=d[(d.cohort==row.cohort)&(d.member==row.member)&(d.scenario==row.scenario)][row.quantity].to_numpy()
 assert len(v)==row.n and sum(v>0)==row.positive and sum(v<0)==row.negative
 half=t.ppf(.975,len(v)-1)*np.std(v,ddof=1)/np.sqrt(len(v));errs.extend([abs(v.mean()-row.mean),abs(v.mean()-half-row.ci_low),abs(v.mean()+half-row.ci_high)])
m=pd.read_csv(R/'matched_donor_comparisons.csv');errs.append(float(abs(m.alternative-m.matched_primary-m.paired_change).max()))
pt=pd.read_csv(R/'detection_partition_donors.csv');errs.append(float(abs(pt.total_rho-pt.between_detection-pt.within_total).max()))
assert max(errs)<1e-10
files=[]
for f in R.glob('*.csv'):
 assert f.read_bytes()==(P/'qa_replay'/f.name).read_bytes();files.append(f.name)
mask=pd.read_csv(R/'union_mask_audit.csv').sort_values('sample');assert mask.union_cells.tolist()==[4320,4485,6697,5739,4835,5318]
guards={'PIEZO1_OA_Main_V49_REVIEW.docx':'074ce6c4352b004af256ea3c0f84d107399c40160fac67f59026343b8eacd0fa','PIEZO1_OA_Supplement_V49_REVIEW.docx':'339dd65cec824121a94ce34eeddacbb42ab5b8f34e760bcbe49ab39d83d9b85b'}
checked={}
for name,h in guards.items():
 f=P.parent/'V49_State_Integrated_20260915/manuscript'/name
 if f.exists():
  actual=hashlib.sha256(f.read_bytes()).hexdigest();assert actual==h;checked[name]=actual
report=dict(independent_implementation_max_error=max(errs),replay_byte_identical_csv=files,union_cells=int(mask.union_cells.sum()),sample_rows=len(s),donor_rows=len(d),manuscript_hashes=checked,scope='Independent numerical implementation, not independent human or biological validation')
(P/'VERIFICATION.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
print(d.query("cohort=='external' and member=='OGN'").to_string(index=False))
print(pd.read_csv(R/'source_scenario_summaries.csv').query("member in ['PRG4','OGN']")[['cohort','member','scenario','quantity','mean']].to_string(index=False))
print(pd.read_csv(R/'detection_partition_source_summaries.csv').query("member in ['PRG4','OGN'] and channel=='P1_minus_T4' and term in ['total_rho','between_detection','within_total']")[['cohort','member','partition','term','mean']].to_string(index=False))
