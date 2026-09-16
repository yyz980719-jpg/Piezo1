from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent;R=P/'results';errors=[]
S=pd.read_csv(R/'state_region_estimates.csv');D=pd.read_csv(R/'donor_state_estimates.csv');A=pd.read_csv(R/'state_common_admission.csv');M=pd.read_csv(R/'matched_donor_comparisons.csv')
for gsm,rows in S.groupby('GSM'):
 z=np.load(P/'inputs'/(gsm+'.npz'));gs=z['genes'].tolist();E=z['expression'];X=pd.DataFrame({g:E[:,gs.index(g)] for g in ['PIEZO1','TRPV4','PRG4','OGN']});X['M5']=E[:,[gs.index(g) for g in ['COL2A1','SOX9','ACAN','CHAD','HAPLN1']]].mean(1)
 for row in rows.itertuples():
  x=X if row.state=='ALL_CELLS' else X.loc[z['states']==row.state];cols=['PIEZO1','TRPV4',row.target];eligible=len(x)>=200 and (x[cols].gt(0).sum()>=20).all() and (x[cols].var(ddof=0)>0).all()
  assert eligible==(row.gate=='PASS')
  if eligible:
   c=x[cols].rank(method='average').corr().iloc[:2,2].to_numpy(); zz=np.arctanh(c);errors.extend(abs(zz-[row.PIEZO1_z,row.TRPV4_z]));errors.append(abs(zz[0]-zz[1]-row.delta))
for row in D.itertuples():
 s=S[(S.state==row.state)&(S.donor==row.donor)&(S.target==row.target)];assert len(s)==2;assert row.eligible==s.status.eq('ESTIMABLE').all()
 if row.eligible:
  for q in ['PIEZO1_z','TRPV4_z','delta']:errors.append(abs(s[q].mean()-getattr(row,q)))
for row in A.itertuples():
 ids=set.intersection(*[set(D[(D.state==row.state)&(D.target==g)&D.eligible].donor) for g in ['PRG4','OGN','M5']]);assert len(ids)==row.n_common and row.supported==(len(ids)>=6)
for row in M.itertuples():
 b=D[(D.state=='ALL_CELLS')&(D.donor==row.donor)&(D.target==row.target)].iloc[0];n=D[(D.state==row.state)&(D.donor==row.donor)&(D.target==row.target)].iloc[0]
 errors.extend([abs(row.matched_overall-b[row.quantity]),abs(row.within_state-n[row.quantity]),abs(row.paired_change-(n[row.quantity]-b[row.quantity]))])
for row in pd.read_csv(R/'matched_descriptive_summaries.csv').itertuples():
 vals=M[(M.state==row.state)&(M.target==row.target)&(M.quantity==row.quantity)][row.estimate].to_numpy();assert row.n==len(vals) and row.positive==sum(vals>0) and row.negative==sum(vals<0);errors.append(abs(row.mean-vals.mean()))
assert max(errors)<1e-10
files=[]
for f in R.glob('*.csv'):assert f.read_bytes()==(P/'qa_replay'/f.name).read_bytes();files.append(f.name)
guards={'PIEZO1_OA_Main_V49_REVIEW.docx':'074ce6c4352b004af256ea3c0f84d107399c40160fac67f59026343b8eacd0fa','PIEZO1_OA_Supplement_V49_REVIEW.docx':'339dd65cec824121a94ce34eeddacbb42ab5b8f34e760bcbe49ab39d83d9b85b'};verified={}
for name,h in guards.items():
 f=P.parent/'V49_State_Integrated_20260915/manuscript'/name
 if f.exists():assert hashlib.sha256(f.read_bytes()).hexdigest()==h;verified[name]=h
qa=dict(independent_numeric_max_error=float(max(errors)),byte_identical_tables=files,manuscript_unchanged_hashes=verified,independent_human_review=False)
(P/'VERIFICATION.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa,indent=2))
