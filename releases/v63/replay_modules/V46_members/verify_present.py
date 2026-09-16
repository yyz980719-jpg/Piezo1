from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,itertools,zipfile,datetime
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
P=Path(__file__).resolve().parent;I=P/'inputs';O=P/'results';F=P/'figures'
F.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for r in json.loads((P/'input_hashes.json').read_text()):assert sha(I/r['copy'])==r['sha256']
assert sha(P/'PLAN.md')==json.loads((P/'FREEZE_RECEIPT.json').read_text())['plan_sha256']
rg=pd.read_csv(O/'regional_donor.csv');pr=pd.read_csv(I/'prior_sample.csv');pr=pr[pr.selection.eq('primary')];ss=pd.read_csv(O/'member_sample.csv');dd=pd.read_csv(O/'member_donor.csv');summary=pd.read_csv(O/'member_summary.csv')
for r in rg.itertuples():
 z=pr[pr.cohort.eq('discovery')&pr.donor.eq(r.donor)&pr.definition.eq(r.definition)].set_index('region')
 for k in ['delta','PIEZO1_z','TRPV4_z']:assert abs(getattr(r,k+'_change')-(z.loc['WB',k]-z.loc['NWB',k]))<1e-12
v=rg[rg.definition.eq('M5')].delta_change.to_numpy();observed=float(v.mean());extreme=0
for code in range(2**len(v)):
 value=sum((1 if code&(1<<j) else -1)*float(x) for j,x in enumerate(v))/len(v);extreme+=abs(value)>=abs(observed)-1e-12
assert extreme/(2**len(v))==pd.read_csv(O/'regional_primary_inference.csv').raw_p.iloc[0]
q1,q3=np.quantile(v,[.25,.75]);flags=((v<q1-1.5*(q3-q1))|(v>q3+1.5*(q3-q1)));iqr={'n_flagged':int(flags.sum()),'donors':rg[rg.definition.eq('M5')].loc[flags,'donor'].tolist(),'action':'retain all donors'}
# Independently replay all 374 gene/score sample correlations from packaged inputs.
U11='COL9A1 COL9A2 COL9A3 COL11A1 COL11A2 COMP MATN3 PRG4 CILP CILP2 OGN'.split();errs=[]
for sample,t in ss.groupby('sample'):
 z=np.load(I/(sample+'.npz'));g=z['genes'].tolist();a=z['expression'].copy();raw=a.copy()
 if 'archived_expression' in z:
  for j,gene in enumerate(z['archived_genes'].tolist()):a[:,g.index(gene)]=z['archived_expression'][:,j]
 for r in t.itertuples():
  x=raw[:,[g.index(k) for k in U11]].mean(1) if r.member=='U11' else a[:,g.index(r.member)]
  assert len(x)==r.cells and int((x>0).sum())==r.member_detected
  rho=[]
  for ch in ['PIEZO1','TRPV4']:
   q=pd.Series(a[:,g.index(ch)]).rank(method='average').to_numpy();w=pd.Series(x).rank(method='average').to_numpy();q=q-q.mean();w=w-w.mean();val=float(q@w/np.sqrt((q@q)*(w@w)));rho.append(val);errs.append(abs(val-getattr(r,ch+'_rho')))
  assert abs(np.arctanh(rho[0])-np.arctanh(rho[1])-r.delta)<1e-10
assert max(errs)<1e-10 and len(ss)==22*17 and ss.status.eq('ESTIMABLE').all()
for r in dd.itertuples():
 z=ss[ss.cohort.eq(r.cohort)&ss.donor.eq(r.donor)&ss.member.eq(r.member)]
 assert len(z)==(2 if r.cohort=='discovery' else 1)
 for k in ['PIEZO1_z','TRPV4_z','delta']:assert abs(z[k].mean()-getattr(r,k))<1e-12
st=pd.read_csv(I/'state_regions.csv');pa=pd.read_csv(O/'state_pair_admission.csv')
for r in pa.itertuples():
 z=st[st.state.eq(r.state)&st.donor.eq(r.donor)];ok=len(z)==2 and set(z.region)=={'WB','NWB'} and z.status.eq('ESTIMABLE').all();assert ok==r.eligible
 if ok:assert z.cells.ge(200).all() and z.P1_detected.ge(20).all() and z.T4_detected.ge(20).all()
# Replace only presentation artifacts; immutable analysis code/results remain unchanged.
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
def save(fig,n):
 fig.savefig(F/(n+'.png'),dpi=300,facecolor='white');fig.savefig(F/(n+'.svg'),facecolor='white');plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(11,4.3),layout='constrained')
for ax,definition in zip(axes[:2],['M5','M14']):
 z=rg[rg.definition.eq(definition)].sort_values('donor');y=np.arange(len(z));ax.hlines(y,z.delta_NWB,z.delta_WB,color='#B7C1CD',lw=.9)
 ax.scatter(z.delta_NWB,y,s=21,facecolors='white',edgecolors='#263447',marker='o',label='NWB');ax.scatter(z.delta_WB,y,s=21,color='#263447',marker='s',label='WB');ax.axvline(0,color='.6',ls=':');ax.set(yticks=y,yticklabels=z.donor,xlabel='Channel contrast (Fisher z)',title=definition,xlim=(-.58,.025));ax.invert_yaxis();ax.legend(frameon=False,loc='lower left')
z=rg[rg.definition.eq('M5')].sort_values('donor');ax=axes[2];ax.scatter(z.delta_change,range(8),color='#263447',s=22);ax.axvline(0,color='.6',ls=':');ax.set(yticks=range(8),yticklabels=z.donor,xlabel='R5 = ΔWB − ΔNWB',title='Paired regional difference');ax.invert_yaxis()
save(fig,'A_regional_all_donors')
G='COL2A1 ACAN HAPLN1 SOX9 CHAD COL9A1 COL9A2 COL9A3 COL11A1 COL11A2 COMP MATN3 PRG4 CILP CILP2 OGN U11'.split()
fig,axes=plt.subplots(1,2,figsize=(10,7.2),layout='constrained')
for ax,co in zip(axes,['discovery','external']):
 for j,ch,c,m in [(0,'PIEZO1_z','#397BB6','o'),(1,'TRPV4_z','#C34F83','s')]:
  z=summary[summary.cohort.eq(co)&summary.quantity.eq(ch)].set_index('member').reindex(G);ax.scatter(z['mean'],np.arange(17)+(j-.5)*.18,color=c,marker=m,s=15,label=ch.replace('_z',''));ax.hlines(np.arange(17)+(j-.5)*.18,z.ci_low,z.ci_high,color=c,lw=.7)
 ax.axvline(0,color='.5',ls=':');ax.set(yticks=range(17),yticklabels=G,xlabel='Mean donor Fisher z (unadjusted t95 CI)',title=co,xlim=(-.11,.45));ax.invert_yaxis();ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,1.13),ncol=2)
save(fig,'B_absolute_associations')
meta=[]
for f in F.glob('*.png'):
 im=Image.open(f);assert abs(im.info['dpi'][0]-300)<.1;meta.append({'file':f.name,'pixels':list(im.size),'dpi':im.info['dpi']})
receipt={'independent_sample_correlations':len(errs),'max_absolute_error':max(errs),'all_374_sample_estimates_estimable':True,'all_238_donor_estimates_verified':len(dd)==238,'regional_algebra':'PASS','sign_flip_extreme_configurations':extreme,'sign_flip_configurations':2**len(v),'R5_IQR':iqr,'state_pair_gates':'PASS','input_hashes':'PASS','figure_metadata':meta,'presentation_repairs':'Replaced colliding endpoint labels by donor-row paired points; moved legends outside data; shared absolute-association x limits. No data edits.','new_test_count':1,'state_analysis':'Admission only, no new state effects'}
(O/'INDEPENDENT_QA.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
