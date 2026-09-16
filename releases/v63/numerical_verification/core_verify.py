from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import t
from itertools import product
import json, hashlib, shutil, sys

P=Path(__file__).resolve().parent
O=P.parent
for d in ['audit','source_tables','manuscript','qa','provenance']:(P/d).mkdir(exist_ok=True)
checks=[]; sources={}; results=[]
def read(path):
 path=Path(path); key=path.parent.parent.name+'__'+path.name
 actual=P/'source_tables'/key if '--archived' in sys.argv or not path.exists() else path
 if actual==path:shutil.copy2(path,P/'source_tables'/key)
 sources[str(path)]=hashlib.sha256(actual.read_bytes()).hexdigest()
 return pd.read_csv(actual,sep='\t' if path.suffix=='.tsv' else ',')
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
def near(name,a,b,tol=1e-12):
 err=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
 check(name,np.isfinite(err) and err<=tol,err)
def exact(x,tol=1e-12):
 x=np.array(x); z=np.array(list(product([-1,1],repeat=len(x))))@x/len(x)
 return float(np.mean(abs(z)>=abs(x.mean())-tol))
def adjust(p,kind):
 p=np.asarray(p); order=np.argsort(p); n=len(p); out=np.empty(n)
 v=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1] if kind=='BH' else np.maximum.accumulate(p[order]*np.arange(n,0,-1))
 out[order]=np.minimum(v,1); return out
def summary(label,x,test=False):
 x=np.asarray(x,float); n=len(x); m=x.mean(); h=t.ppf(.975,n-1)*x.std(ddof=1)/np.sqrt(n)
 r=dict(endpoint=label,n=n,mean=m,ci95_low=m-h,ci95_high=m+h,negative=int(sum(x<0)),positive=int(sum(x>0)))
 if test:r['raw_p']=exact(x)
 results.append(r);return r
def identities(df,label):
 q=df.dropna(subset=['PIEZO1_z','TRPV4_z','delta'])
 near(label+' delta identity',q.delta,q.PIEZO1_z-q.TRPV4_z)
def pairing(samples,donors,keys,label):
 q=samples[samples.status=='ESTIMABLE']; counts=q.groupby(keys).size()
 check(label+' regional multiplicity',all(counts.loc[[i for i in counts.index if i[0]=='discovery']]==2))
 means=q.groupby(keys)[['PIEZO1_z','TRPV4_z','delta']].mean()
 merged=donors.merge(means.reset_index(),on=keys,suffixes=('_saved','_recomputed'),validate='one_to_one')
 check(label+' donor coverage',len(merged)==len(donors))
 for col in ['PIEZO1_z','TRPV4_z','delta']:near(label+' equal-region '+col,merged[col+'_saved'],merged[col+'_recomputed'])

B=Path('E:/codex/2026-09-13/new-chat-2/outputs/V44_Matrix_Definition_Bridge/results')
d=read(B/'donor_estimates.csv'); s=read(B/'sample_estimates.csv')
check('bridge unique donor keys',not d.duplicated(['cohort','selection','donor','definition']).any())
identities(d,'bridge donor'); identities(s,'bridge sample')
pairing(s,d,['cohort','selection','donor','definition'],'bridge')
ps=[]
for cohort in ['discovery','external']:
 q=d[(d.selection=='primary')&(d.cohort==cohort)]
 pivot=q.pivot(index='donor',columns='definition',values='delta')
 for target in ['M5','M14']:summary(cohort+' '+target,pivot[target],True)
 summary(cohort+' D',pivot.M14-pivot.M5,True)
 ps.append(exact(pivot.M14))
ps += [next(r['raw_p'] for r in results if r['endpoint']==c+' D') for c in ['discovery','external']]
fam=read(B/'fixed_four_test_family.csv')
near('Holm4 raw P reconstruction',ps,fam.raw_p)
near('Holm4 adjustment',adjust(ps,'Holm'),fam.Holm4_adjusted_p)
L=O/'PIEZO1_OA_Round1_Revision_20260908/audit'
bh2=read(L/'p15_components_v2/core_tests_BH2.tsv'); bh4=read(L/'p21_external_effects/primary_BH4_family.tsv')
near('BH2 original family adjustment',adjust(bh2.p,'BH'),bh2.q_BH2)
near('BH4 original family adjustment including unavailable slots',adjust(bh4.p_for_BH4,'BH'),bh4.BH4_accounting_q)
near('M5 discovery independent P',results[0]['raw_p'],bh2.loc[bh2.target=='M','p'].iloc[0])
near('M5 external independent P',next(r['raw_p'] for r in results if r['endpoint']=='external M5'),bh4.loc[bh4.endpoint=='DeltaM','observed_p'].dropna().iloc[0])
for c in ['discovery','external']:
 r=next(r for r in results if r['endpoint']==c+' M5')
 expected=(-.2939357580516357,8) if c=='discovery' else (-.32908,6)
 near(c+' manuscript M5 rounding',r['mean'],expected[0],5e-6)
 check(c+' donor count and all-negative M5',r['n']==r['negative']==expected[1])

B=O/'V46_Bounded_Extension_20260915/results'
m=read(B/'member_donor.csv'); ms=read(B/'member_sample.csv')
identities(m,'all16 donor'); identities(ms,'all16 sample')
check('member unique keys',not m.duplicated(['cohort','donor','member']).any())
pairing(ms,m,['cohort','donor','member'],'all16')
classification=[]
for (cohort,member),g in m.groupby(['cohort','member']):
 r=summary(cohort+' member '+member,g.delta);classification.append(r)
for c in ['discovery','external']:
 q=[r for r in classification if r['endpoint'].startswith(c+' ') and not r['endpoint'].endswith(' U11')]
 check(c+' complete 16-member map',len(q)==16)
 check(c+' exactly 12 uniformly negative members',sum(r['negative']==r['n'] for r in q)==12)
 check(c+' only PRG4 and OGN uniformly positive',{r['endpoint'].split()[-1] for r in q if r['positive']==r['n']}=={'PRG4','OGN'})
reg=read(B/'regional_donor.csv'); r=reg[reg.definition=='M5']
near('R5 within-donor WB minus NWB',r.delta_change,r.delta_WB-r.delta_NWB)
rr=summary('R5 WB minus NWB',r.delta_change,True);near('R5 manuscript P',rr['raw_p'],.78125)
near('R5 manuscript mean',rr['mean'],.01237,5e-6)

B=O/'V50_Member_Measurement_20260915/results'
me=read(B/'member_donor_scenarios.csv');identities(me,'member measurement')
check('measurement unique keys',not me.duplicated(['cohort','donor','member','scenario']).any())
base=me[me.scenario=='primary'].merge(m,on=['cohort','donor','member'],suffixes=('_new','_old'),validate='one_to_one')
near('all16 primary replay across versions',base.delta_new,base.delta_old)
for (c,g,sc),q in me.groupby(['cohort','member','scenario']):
 check(f'{c} {g} {sc} complete donor support',len(q)==(8 if c=='discovery' else 6) and q.status.eq('ESTIMABLE').all())
 if g in ['PRG4','OGN']:summary(c+' '+g+' '+sc,q.delta)
pr=me[me.member=='PRG4'];check('PRG4 every examined measurement contrast positive',pr.delta.gt(0).all())
og=me[(me.member=='OGN')&(me.cohort=='external')]
check('OGN external depth 5 of 6 positive',sum(og[og.scenario=='primary_cubic_depth'].delta>0)==5)
check('OGN external union 4 of 6 positive',sum(og[og.scenario.str.contains('union')].delta>0)==4)

B=O/'V51_Member_State_20260915/results'
st=read(B/'matched_donor_comparisons.csv'); ds=read(B/'donor_state_estimates.csv'); sr=read(B/'state_region_estimates.csv')
identities(ds,'state donor'); identities(sr,'state region')
check('state estimate unique keys',not ds.duplicated(['state','donor','target']).any())
valid=ds[ds.eligible==True]; means=sr[sr.status=='ESTIMABLE'].groupby(['state','donor','target']).agg(delta=('delta','mean'),n=('delta','size')).reset_index()
joined=valid.merge(means,on=['state','donor','target'],validate='one_to_one',suffixes=('_saved','_recomputed'))
check('state paired-region completeness',len(joined)==len(valid) and joined.n.eq(2).all())
near('state equal-region donor reconstruction',joined.delta_saved,joined.delta_recomputed)
q=st[(st.supported==True)&(st.quantity=='delta')]
check('five admitted states',set(q.state)=={'HomC','RegC','RepC','preFC','preHTC'})
for state,g in q.groupby('state'):
 sets=[set(x.donor) for _,x in g.groupby('target')]
 check(state+' common target donor identities',len(sets)==3 and sets[0]==sets[1]==sets[2] and len(sets[0])>=6)
 for target,z in g.groupby('target'):
  summary(state+' '+target+' within',z.within_state)
  summary(state+' '+target+' matched overall',z.matched_overall)
  original=ds[(ds.state=='ALL_CELLS')&(ds.target==target)].set_index('donor').loc[z.donor,'delta']
  near(state+' '+target+' matched overall identities',z.matched_overall,original)
rep=q[(q.state=='RepC')&(q.target=='PRG4')]
check('PRG4 RepC six opposite donor orderings',len(rep)==6 and rep.within_state.lt(0).all() and rep.matched_overall.gt(0).all())
near('PRG4 RepC manuscript mean',rep.within_state.mean(),-.04603,5e-6)
near('PRG4 matched overall manuscript mean',rep.matched_overall.mean(),.07905,5e-6)

B=O/'V48_Region_State_20260915/results'
v=read(B/'donor_state_regional_changes.csv')
near('regional state R orientation',v.R,v.contrast_z_WB-v.contrast_z_NWB)
mat=v.pivot(index='donor',columns='state',values='R').values
check('regional state 7 by 7 completeness',mat.shape==(7,7) and np.isfinite(mat).all())
center=mat-mat.mean(axis=1,keepdims=True); obs=np.mean(center.mean(axis=0)**2)
null=np.array([np.mean((center*np.array(sign)[:,None]).mean(axis=0)**2) for sign in product([-1,1],repeat=len(mat))])
p=float(np.mean(null>=obs-1e-12));near('regional state exact omnibus P',p,.375)
near('regional state RMS',np.sqrt(obs),.038003682029611977)
report=dict(scope='Independent recomputation from archived sample and donor estimates; not a raw-read, cell-calling, doublet-model or cell-level correlation replay.',new_analyses=0,new_hypotheses=0,checks=len(checks),passed=sum(x['passed'] for x in checks),failed=[x for x in checks if not x['passed']],omnibus_p=p,sources=sources)
pd.DataFrame(checks).to_csv(P/'audit/checks.csv',index=False)
pd.DataFrame(results).to_csv(P/'audit/recomputed_results.csv',index=False)
(P/'audit/verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='sources'},indent=2))
assert not report['failed'],report['failed']
