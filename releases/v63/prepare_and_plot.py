"""Current descriptive aggregation and figures. No new hypothesis tests.
Run with Python 3, numpy, pandas, scipy, matplotlib and Pillow.
After first staging, --portable uses only the accompanying source_data directory.
"""
from pathlib import Path
import sys, shutil, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import t
P=Path(__file__).resolve().parent
S=P/'source_data'; F=P/'figures'; Q=P/'qa'
for d in [S,F,Q]: d.mkdir(exist_ok=True)
B=P.parent
source_map={
'bridge.csv':B/'V55_Core_Verification_Rewrite_20260915/source_tables/V44_Matrix_Definition_Bridge__donor_estimates.csv',
'members.csv':B/'V55_Core_Verification_Rewrite_20260915/source_tables/V46_Bounded_Extension_20260915__member_donor.csv',
'states.csv':B/'V55_Core_Verification_Rewrite_20260915/source_tables/V51_Member_State_20260915__donor_state_estimates.csv',
'state_matched.csv':B/'V55_Core_Verification_Rewrite_20260915/source_tables/V51_Member_State_20260915__matched_donor_comparisons.csv',
'sensitivity.csv':B/'V38_BluePink_20260913/source_data/sensitivity.csv',
'thinning_samples.csv':B/'V35_review_revision_20260913/Review_bundle/sample_estimates.csv',
'thinning_diagnostics.csv':B/'V35_review_revision_20260913/Review_bundle/review_outputs/thinning_before_after_measurement_diagnostics.csv',
'donor_crosswalk.tsv':B/'PIEZO1_OA_Round1_Revision_20260908/audit/p33_review/S34_source_donor_run_crosswalk.tsv',
'M14_dictionary.csv':Path('E:/codex/2026-09-13/new-chat-2/outputs/V44_Matrix_Definition_Bridge/audit/source_module_dictionary.csv'),
}
if False:  # V60 always uses bundled saved source tables
 for n,p in source_map.items(): shutil.copy2(p,S/n)
 (S/'source_provenance.json').write_text(json.dumps({n:{'original':str(p),'sha256':hashlib.sha256((S/n).read_bytes()).hexdigest()} for n,p in source_map.items()},indent=2),encoding='utf-8')
bridge=pd.read_csv(S/'bridge.csv'); bridge=bridge[bridge.selection.eq('primary')]
members=pd.read_csv(S/'members.csv'); matched=pd.read_csv(S/'state_matched.csv')
sens=pd.read_csv(S/'sensitivity.csv'); samples=pd.read_csv(S/'thinning_samples.csv')
diag=pd.read_csv(S/'thinning_diagnostics.csv')
EIDS=['2019-178','2019-211','2020-16','2021-092','2021-097','2021-101']
def display(x):
 x=str(x).replace('GSE220243:','')
 return 'CartOA'+str(EIDS.index(x)+1) if x in EIDS else x
for d in [bridge,members,matched,sens,samples,diag]:
 if 'donor' in d: d['display_donor']=d.donor.map(display)
pd.DataFrame([{'cohort':'discovery','source_donor':f'OA{i}','display_donor':f'OA{i}'} for i in range(1,9)]+[{'cohort':'external','source_donor':x,'display_donor':display(x)} for x in EIDS]).to_csv(S/'display_crosswalk.csv',index=False)
# Reconstruct absolute donor z from saved regional rho, then average 20 draws.
a=samples[samples['mode'].isin(['depth1000_unthinned','depth1000_thinned']) & samples.status.eq('PASS')].copy()
for ch in ['PIEZO1','TRPV4']: a[ch+'_z']=np.arctanh(a[ch+'_rho'])
az=a.groupby(['cohort','display_donor','mode','draw'],as_index=False)[['PIEZO1_z','TRPV4_z','contrast_z']].mean()
assert np.allclose(az.PIEZO1_z-az.TRPV4_z,az.contrast_z,atol=1e-12)
az=az.groupby(['cohort','display_donor','mode'],as_index=False)[['PIEZO1_z','TRPV4_z','contrast_z']].mean()
az.to_csv(S/'thinning_absolute_donor.csv',index=False)
az.groupby(['cohort','mode'],as_index=False)[['PIEZO1_z','TRPV4_z','contrast_z']].mean().to_csv(S/'thinning_absolute_summary.csv',index=False)
diag['condition']=np.where(diag.draw.eq(-1),'original','thinned')
dc=['depth_mean','PIEZO1_detection','TRPV4_detection','joint_positive_fraction','M_all_zero_fraction']
dg=diag.groupby(['cohort','display_donor','condition','draw'],as_index=False)[dc].mean().groupby(['cohort','display_donor','condition'],as_index=False)[dc].mean()
dg.to_csv(S/'thinning_diagnostic_donor.csv',index=False)
dg.groupby(['cohort','condition'],as_index=False)[dc].mean().to_csv(S/'thinning_diagnostic_summary.csv',index=False)
cw=pd.read_csv(S/'donor_crosswalk.tsv',sep='\t'); cw=cw[cw.GSE.eq('GSE220243')]
clinical=[{'cohort':'discovery','donor':f'OA{i+1}','age':a,'sex':s,'GEO':'GSE255460'} for i,(a,s) in enumerate(zip([58,51,61,73,61,70,75,66],['F','F','M','M','F','M','F','M']))]
for r in cw.itertuples():
 j=json.loads(r.source_attributes);clinical.append({'cohort':'external','donor':display(r.source_patient_id),'age':int(j['age']),'sex':j['sex'][0].upper(),'GEO':r.GSM})
pd.DataFrame(clinical).to_csv(S/'clinical_donors.csv',index=False)
BLUE='#397BB6'; PINK='#C34F83'; INK='#263447'; GREY='#8F99A4'; TEAL='#197C80'; ORANGE='#C77B32'
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':7,'xtick.labelsize':6,'ytick.labelsize':6,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,'xtick.major.width':.5,'ytick.major.width':.5,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'savefig.facecolor':'white'})
WIDTH=183/25.4
def panel(ax,letter,title):
 ax.text(-.12,1.09,letter,transform=ax.transAxes,fontweight='bold',fontsize=8,va='bottom'); ax.set_title(title,loc='left',pad=11)
def zero(ax,vertical=False):
 (ax.axvline if vertical else ax.axhline)(0,color='#CAD0D6',lw=.65,zorder=0)
def save(fig,n):
 for ext in ['pdf','svg','png','tiff']:
  kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
  fig.savefig(F/f'Figure_{n}.{ext}',dpi=600,**kw)
 fig.savefig(Q/f'Figure_{n}_preview.png',dpi=160);plt.close(fig)
def mean_ci(v):
 v=np.asarray(v);return v.mean(),t.ppf(.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))
def paired_channels(ax,d):
 for r in d.itertuples(): ax.plot([0,1],[r.PIEZO1_z,r.TRPV4_z],color=GREY,lw=.65,alpha=.6,zorder=1)
 for j,ch in enumerate(['PIEZO1','TRPV4']):
  ax.scatter(np.repeat(j,len(d)),d[ch+'_z'],s=13,color=[BLUE,PINK][j],alpha=.75,zorder=2)
  ax.scatter(j,d[ch+'_z'].mean(),marker='D',s=34,color=[BLUE,PINK][j],edgecolors='white',linewidth=.6,zorder=3)
 ax.set_xticks([0,1],['PIEZO1','TRPV4']);ax.set_xlim(-.4,1.4);zero(ax);ax.set_ylabel('Association (Fisher z)')
# F1: source separation and common absolute scale.
fig,axes=plt.subplots(2,2,figsize=(WIDTH,5.1),gridspec_kw={'width_ratios':[1,1.5]});fig.subplots_adjust(left=.10,right=.98,bottom=.09,top=.92,hspace=.65,wspace=.42)
for row,(cohort,label) in enumerate([('discovery','Discovery · 8 donors'),('external','Second source · 6 donors')]):
 d=bridge[bridge.cohort.eq(cohort)&bridge.definition.eq('M5')].sort_values('display_donor');paired_channels(axes[row,0],d);axes[row,0].set_ylim(-.09,.62);panel(axes[row,0],chr(97+2*row),label)
 ax=axes[row,1];y=np.arange(len(d));ax.scatter(d.delta,y,s=20,color=TEAL);ax.set_yticks(y,d.display_donor);m,e=mean_ci(d.delta);ax.errorbar(m,-1.3,xerr=e,fmt='D',color=INK,capsize=3,ms=4,lw=1);ax.text(-.01,-1.3,'Mean',fontsize=6,ha='right',va='center');zero(ax,True);ax.set_xlim(-.6,.01);ax.set_ylim(-2,len(d));ax.set_xlabel('ΔM5  (PIEZO1 − TRPV4, Fisher z)');panel(ax,chr(98+2*row),'Donor contrasts');ax.text(.99,.97,'Unadjusted\n95% t interval',transform=ax.transAxes,fontsize=6,ha='right',va='top');ax.text(.99,1.04,'BH2 q = 0.0078' if row==0 else 'BH4 q = 0.0625',transform=ax.transAxes,ha='right',fontsize=6)
save(fig,1)
# F2: all fixed genes, not a selected positive subset; score contrasts use same donors.
order=['COL2A1','ACAN','HAPLN1','SOX9','CHAD','COL9A1','COL9A2','COL9A3','COL11A1','COL11A2','COMP','MATN3','PRG4','CILP','CILP2','OGN']
cmap=LinearSegmentedColormap.from_list('contrast',[TEAL,'#FAFAFA',ORANGE]);norm=TwoSlopeNorm(vmin=-.65,vcenter=0,vmax=.65)
fig=plt.figure(figsize=(WIDTH,5.4));gs=fig.add_gridspec(2,2,height_ratios=[3.7,1.1],left=.13,right=.91,top=.91,bottom=.09,wspace=.23,hspace=.53)
for j,cohort in enumerate(['discovery','external']):
 ax=fig.add_subplot(gs[0,j]);ds=[f'OA{i}' for i in range(1,9)] if j==0 else [f'CartOA{i}' for i in range(1,7)]
 a=members[members.cohort.eq(cohort)].pivot(index='member',columns='display_donor',values='delta').loc[order,ds]
 im=ax.imshow(a,cmap=cmap,norm=norm,aspect='auto');ax.set_yticks(range(16),order if j==0 else []);ax.set_xticks(range(len(ds)),ds,rotation=45,ha='right');panel(ax,chr(97+j),'Discovery' if j==0 else 'Second source')
 for k in [2.5,4.5]:ax.axhline(k,color='white',lw=1.5)
 for spine in ax.spines.values():spine.set_visible(False)
 ax=fig.add_subplot(gs[1,j]);
 scorewide=pd.DataFrame({score:(bridge[bridge.cohort.eq(cohort)&bridge.definition.eq(score)] if score!='U11' else members[members.cohort.eq(cohort)&members.member.eq('U11')]).set_index('display_donor').delta for score in ['M5','M14','U11']}).loc[ds]
 assert scorewide.notna().all().all()
 for ii,(_,vals) in enumerate(scorewide.iterrows()):ax.plot(np.arange(3)+ii*.025-.08,vals.to_numpy(),color='#C8CDD2',lw=.5,zorder=0)
 for k,score in enumerate(['M5','M14','U11']):
  d=bridge[bridge.cohort.eq(cohort)&bridge.definition.eq(score)] if score!='U11' else members[members.cohort.eq(cohort)&members.member.eq('U11')]
  d=d.set_index('display_donor').loc[ds].reset_index()
  ax.scatter(np.arange(len(d))*.025+k-.08,d.delta,s=11,color=GREY);ax.plot([k-.15,k+.15],[d.delta.mean()]*2,color=INK,lw=2)
 ax.set_xticks(range(3),['M5','M14','U11']);zero(ax);ax.set_ylim(-.6,.025);ax.set_ylabel('Score contrast Δ');panel(ax,chr(99+j),'Fixed score definitions')
cax=fig.add_axes([.94,.39,.018,.46]);fig.colorbar(im,cax=cax,ticks=[-.6,-.3,0,.3,.6]);cax.set_title('Δz',fontsize=7,pad=8);save(fig,2)
# F3: sensitivity plus both absolute channel associations before/after thinning.
fig=plt.figure(figsize=(WIDTH,6.25));gs=fig.add_gridspec(3,2,left=.10,right=.98,bottom=.075,top=.94,wspace=.34,hspace=.62)
axes=np.array([[fig.add_subplot(gs[i,j]) for j in range(2)] for i in range(3)])
for j in range(2):
 pos=axes[0,j].get_position();axes[0,j].set_position([pos.x0+.12,pos.y0,pos.width-.12,pos.height])
modes=['all_cells','matched','all_cells_depth_cubic','joint_positive','joint_positive_depth_cubic','depth1000_unthinned','depth1000_thinned_mean20']
labels=['All cells','Matched donors','Depth adjusted','Jointly positive','Jointly positive + depth','Depth eligible','Thinned · 20 draws']
for j,cohort in enumerate(['discovery','external']):
 ax=axes[0,j];d=sens[sens.cohort.eq(cohort)&sens.status.eq('PASS')]
 for k,mode in enumerate(modes):
  ds=d[d['mode'].eq('all_cells') & d.donor.isin(d[d['mode'].eq('joint_positive')].donor)] if mode=='matched' else d[d['mode'].eq(mode)]
  ax.scatter(ds.contrast_z,np.repeat(k,len(ds)),s=11,color=GREY);ax.plot(ds.contrast_z.mean(),k,marker='|',ms=9,color=INK);ax.text(.015,k,str(len(ds)),fontsize=6,va='center')
 ax.set_yticks(range(7),labels if j==0 else []);ax.invert_yaxis();zero(ax,True);ax.set_xlim(-.6,.05);ax.set_xlabel('ΔM5');panel(ax,chr(97+j),'Discovery' if j==0 else 'Second source')
 for k,ch in enumerate(['PIEZO1','TRPV4']):
  ax=axes[k+1,j];p=az[az.cohort.eq(cohort)].pivot(index='display_donor',columns='mode',values=ch+'_z');cols=['depth1000_unthinned','depth1000_thinned'];co=[BLUE,PINK][k]
  for _,r in p.iterrows():ax.plot([0,1],r[cols],color=co,alpha=.45,lw=.7,marker='o',ms=2.5)
  ax.plot([0,1],p[cols].mean(),color=co,lw=2,marker='D',ms=4);ax.set_xticks([0,1],['Original','Thinned']);ax.set_xlim(-.3,1.3);ax.set_ylim(-.12,.62);zero(ax);ax.set_ylabel('Association (Fisher z)');panel(ax,chr(99+2*k+j),ch+' · '+('Discovery' if j==0 else 'Second source'))
save(fig,3)
# F4: matched donor contrasts, all five supported states and all three targets.
states=['HomC','RegC','RepC','preFC','preHTC'];targets=['PRG4','OGN','M5']
fig,axes=plt.subplots(1,3,figsize=(WIDTH,4.4));fig.subplots_adjust(left=.12,right=.985,bottom=.14,top=.86,wspace=.27)
for j,target in enumerate(targets):
 ax=axes[j]
 for k,state in enumerate(states):
  d=matched[matched.supported.eq(True)&matched.target.eq(target)&matched.state.eq(state)&matched.quantity.eq('delta')].sort_values('display_donor')
  for o,r in enumerate(d.itertuples()):
   y=k+(o-(len(d)-1)/2)*.065;ax.plot([r.matched_overall,r.within_state],[y,y],color='#CBD0D5',lw=.6);ax.scatter(r.matched_overall,y,s=10,facecolors='white',edgecolors=GREY,marker='s',linewidth=.7);ax.scatter(r.within_state,y,s=11,color=INK)
  ax.scatter(d.matched_overall.mean(),k+.34,marker='D',facecolors='white',edgecolors=GREY,s=20);ax.scatter(d.within_state.mean(),k+.34,marker='D',color=INK,s=20)
 ax.set_yticks(range(5),[f'{s}  (n={7 if s in ["HomC","RegC"] else 6 if s=="RepC" else 8})' for s in states] if j==0 else []);ax.invert_yaxis();ax.set_ylim(4.6,-.5);ax.set_xlim((-.65,.16) if target=='M5' else (-.2,.21));zero(ax,True);ax.set_xlabel('Channel contrast Δ (Fisher z)');panel(ax,chr(97+j),target)
axes[2].text(1,1.045,'Wider x range',transform=axes[2].transAxes,ha='right',fontsize=6.5,fontweight='bold')
fig.text(.5,.96,'Open squares: overall     Filled circles: within state     Diamonds: means',ha='center',fontsize=7);save(fig,4)
# F5: same six donors, absolute correlations prevent an exaggerated reversal claim.
fixed=pd.read_csv(P/'analysis/reference_v61/donor_estimates.csv')
fig,axes=plt.subplots(2,3,figsize=(WIDTH,4.7));fig.subplots_adjust(left=.10,right=.98,bottom=.12,top=.84,wspace=.40,hspace=.72)
for i,mode in enumerate(['unadjusted','cubic_depth']):
 for j,(quantity,title) in enumerate([('PIEZO1_z','PIEZO1–PRG4'),('TRPV4_z','TRPV4–PRG4'),('delta','Relative ordering')]):
  d=fixed[fixed['mode'].eq(mode)].pivot(index='donor',columns='population',values=quantity).sort_index();assert len(d)==6 and d.notna().all().all()
  ax=axes[i,j];color=[BLUE,PINK,INK][j]
  for _,r in d.iterrows():ax.plot([0,1],[r.ALL_CELLS,r.RepC],color=color,lw=.7,alpha=.5,marker='o',ms=3)
  ax.plot([0,1],[d.ALL_CELLS.mean(),d.RepC.mean()],color=color,lw=2,marker='D',ms=5);zero(ax);ax.set_xticks([0,1],['Overall','RepC']);ax.set_xlim(-.25,1.25);ax.set_ylim(-.15,.2);ax.set_ylabel('Fisher z' if j<2 else 'ΔPRG4 (Fisher z)');panel(ax,chr(97+3*i+j),title)
fig.text(.5,.975,'PRG4 · the same six donors in every panel',ha='center',fontsize=8)
fig.text(.10,.925,'Unadjusted',fontsize=8,fontweight='bold')
fig.text(.10,.485,'Cubic depth adjusted',fontsize=8,fontweight='bold')
save(fig,5)
checks={'primary_D':float(bridge.query('cohort=="discovery" and definition=="M5"').delta.mean()),'primary_E':float(bridge.query('cohort=="external" and definition=="M5"').delta.mean()),'thinning_summary':az.groupby(['cohort','mode']).contrast_z.mean().to_dict().__repr__(),'new_hypothesis_tests':0,'n_figures':5,'all_fixed_members':len(order)}
assert abs(checks['primary_D']+.29394)<.00001 and abs(checks['primary_E']+.32908)<.00001
(Q/'plot_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2));print(pd.read_csv(S/'thinning_absolute_summary.csv').to_string(index=False))
