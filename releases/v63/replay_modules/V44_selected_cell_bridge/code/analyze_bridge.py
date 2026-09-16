from audit_inputs import *
from scipy import stats
import itertools,sys,datetime,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# Portable release: unused local skill diagnostics removed. The reported
# Shapiro and IQR diagnostics below remain calculated directly with SciPy/NumPy.

def main():
 receipt=json.loads((P/'FREEZE_RECEIPT.json').read_text());assert sha(P/'FROZEN_PLAN.md')==receipt['plan_sha256']
 assert not (P/'results/RUN_RECEIPT.json').exists(),'Completed run must not be silently replaced'
 dictionary=OTHER/'outputs/phase81_external_review/package/tables/Table_S12_module_dictionary.csv'
 d=read(dictionary);row=d[d.genes.str.split(';').map(lambda x:set(x)==set(M14))];assert len(row)==1 and row.definition_n.iloc[0]==14
 shutil.copy2(dictionary,P/'audit/source_module_dictionary.csv')
 records=read(P/'audit/discovery_sample_audit.csv').to_dict('records')
 ext=read(B/'reporting/Table_S55_External_Donor_Metadata.csv')
 records += [dict(cohort='external',sample=r.GSM,donor=r.source_patient_id,region='source_OA') for r in ext.itertuples()]
 rows=[];qa=[];coverage=[]
 for rec in records:
  f=P/f"inputs/{rec['sample']}.npz";z=np.load(f);genes=z['genes'].tolist();a=z['expression'].copy();counts=z['counts'];old=a.copy()
  if rec['cohort']=='external':
   for j,g in enumerate(z['archived_genes'].tolist()):a[:,genes.index(g)]=z['archived_expression'][:,j]
  for g in GENES:coverage.append(dict(sample=rec['sample'],gene=g,detected=int((counts[:,genes.index(g)]>0).sum()),cells=len(a),identity_present=True))
  for selection in (['primary','intersection'] if rec['cohort']=='discovery' else ['primary']):
   mask=np.ones(len(a),dtype=bool) if selection=='primary' else z['intersection'];aa=a[mask];raw_normalized=old[mask]
   for name,members in [('M5',M5),('M14',M14)]:
    # Keep M14 on one coherent raw-count transformation, M5 on its legacy representation.
    score=(aa if name=='M5' else raw_normalized)[:,[genes.index(g) for g in members]].mean(1)
    channels=aa[:,:2];det=(channels>0).sum(0);status='ESTIMABLE'
    out={k:rec[k] for k in ['cohort','sample','donor','region']};out.update(selection=selection,definition=name,cells=len(aa),P1_detected=int(det[0]),T4_detected=int(det[1]),PIEZO1_rho=np.nan,TRPV4_rho=np.nan,PIEZO1_z=np.nan,TRPV4_z=np.nan,delta=np.nan)
    if len(aa)<200 or min(det)<20:status='CELL_OR_DETECTION_GATE'
    elif np.var(score)==0 or np.any(np.var(channels,axis=0)==0):status='VARIANCE_GATE'
    else:
     rho=np.array([stats.spearmanr(channels[:,j],score).statistic for j in range(2)])
     alternate=pd.DataFrame(dict(P1=channels[:,0],T4=channels[:,1],score=score)).rank().corr().loc[['P1','T4'],'score'].to_numpy()
     err=float(np.max(abs(rho-alternate)));assert err<1e-10;qa.append(dict(sample=rec['sample'],selection=selection,definition=name,rank_check_max_error=err))
     if not (np.isfinite(rho).all() and (abs(rho)<1).all()):status='CORRELATION_GATE'
     else:out.update(PIEZO1_rho=rho[0],TRPV4_rho=rho[1],PIEZO1_z=np.arctanh(rho[0]),TRPV4_z=np.arctanh(rho[1]),delta=np.arctanh(rho[0])-np.arctanh(rho[1]))
    out['status']=status;rows.append(out)
 sample=pd.DataFrame(rows);sample.to_csv(P/'results/sample_estimates.csv',index=False);pd.DataFrame(coverage).to_csv(P/'audit/gene_detection_coverage.csv',index=False)
 donors=[]
 for (cohort,selection,donor,definition),group in sample.groupby(['cohort','selection','donor','definition']):
  valid=group.status.eq('ESTIMABLE').all() and len(group)==(2 if cohort=='discovery' else 1)
  out=dict(cohort=cohort,selection=selection,donor=donor,definition=definition,status='ESTIMABLE' if valid else 'UNAVAILABLE_SAMPLE',cells=int(group.cells.sum()))
  for k in ['PIEZO1_z','TRPV4_z','delta']:out[k]=group[k].mean() if valid else np.nan
  donors.append(out)
 donor=pd.DataFrame(donors);donor.to_csv(P/'results/donor_estimates.csv',index=False)
 paired=donor.pivot(index=['cohort','selection','donor'],columns='definition',values='delta').reset_index();paired['D']=paired.M14-paired.M5
 paired.to_csv(P/'results/paired_definition_change.csv',index=False)
 olddisc=read(I/'audit/p15_components_v2/channel_absolute_and_contrast.tsv');olddisc=olddisc[olddisc.stratum.eq('all_cells')&olddisc.target.eq('M')].set_index('donor_id')
 oldext=read(I/'audit/p21_external_effects/donor_effects.tsv');oldext=oldext[oldext.selection.eq('primary_P19')].set_index('GSM')
 replay=[]
 for r in donor[(donor.selection=='primary')&(donor.definition=='M5')].itertuples():
  if r.cohort=='discovery':expected=olddisc.loc[r.donor,'contrast_z']
  else:expected=oldext.loc[ext.set_index('source_patient_id').loc[r.donor,'GSM'],'DeltaM']
  err=abs(r.delta-expected);assert err<1e-10;replay.append(dict(cohort=r.cohort,donor=r.donor,expected=expected,replay=r.delta,error=err))
 pd.DataFrame(replay).to_csv(P/'qa/M5_historical_replay.csv',index=False);pd.DataFrame(qa).to_csv(P/'qa/independent_rank_checks.csv',index=False)
 summaries=[]
 for (cohort,selection),group in paired.groupby(['cohort','selection']):
  for endpoint in ['M5','M14','D']:
   vals=group[endpoint].dropna().to_numpy();n=len(vals);valid=n>=(6 if cohort=='discovery' else 3)
   mean=vals.mean() if n else np.nan;sd=vals.std(ddof=1) if n>1 else np.nan;half=stats.t.ppf(.975,n-1)*sd/np.sqrt(n) if n>1 else np.nan
   summaries.append(dict(cohort=cohort,selection=selection,endpoint=endpoint,n=n,mean=mean,sd=sd,ci95_low=mean-half,ci95_high=mean+half,negative=int((vals<0).sum()),positive=int((vals>0).sum()),status='ESTIMABLE' if valid else 'UNAVAILABLE_DONOR_GATE'))
 summary=pd.DataFrame(summaries);summary.to_csv(P/'results/summary_estimates.csv',index=False)
 tests=[];diagnostics=[]
 figq,axq=plt.subplots(2,2,figsize=(9,8))
 for slot,(cohort,endpoint) in enumerate([('discovery','M14'),('external','M14'),('discovery','D'),('external','D')],1):
  vals=paired[(paired.cohort==cohort)&(paired.selection=='primary')][endpoint].dropna().to_numpy();valid=len(vals)>=(6 if cohort=='discovery' else 3)
  p=np.nan
  if valid:
   null=np.array(list(itertools.product([-1,1],repeat=len(vals))))@vals/len(vals);p=float(np.mean(abs(null)>=abs(vals.mean())-1e-12))
   q1,q3=np.quantile(vals,[.25,.75]);out=int(((vals<q1-1.5*(q3-q1))|(vals>q3+1.5*(q3-q1))).sum())
   diagnostics.append(dict(cohort=cohort,endpoint=endpoint,n=len(vals),shapiro_W=float(stats.shapiro(vals).statistic),shapiro_p=float(stats.shapiro(vals).pvalue),iqr_flags=out,action='Retain all donors and frozen test; small n cannot establish symmetry or normality.'))
   stats.probplot(vals,dist='norm',plot=axq.flat[slot-1]);axq.flat[slot-1].set_title(cohort+' '+endpoint)
  tests.append(dict(slot=slot,cohort=cohort,endpoint=endpoint,n=len(vals),raw_p=p,status='ESTIMABLE' if valid else 'UNAVAILABLE_DONOR_GATE',p_for_family=p if valid else 1.0))
 test=pd.DataFrame(tests);order=np.argsort(test.p_for_family);adjust=np.maximum.accumulate(np.minimum(1,test.p_for_family.to_numpy()[order]*np.arange(4,0,-1)));test['Holm4_adjusted_p']=np.nan;test.loc[order,'Holm4_adjusted_p']=adjust
 test.loc[test.status!='ESTIMABLE','Holm4_adjusted_p']=np.nan;test.to_csv(P/'results/fixed_four_test_family.csv',index=False)
 pd.DataFrame(diagnostics).to_csv(P/'qa/distribution_diagnostics.csv',index=False);figq.tight_layout();figq.savefig(P/'qa/donor_QQ_diagnostics.png',dpi=180);plt.close(figq)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(1,3,figsize=(13,5),gridspec_kw={'width_ratios':[1,1,1.3]})
 for ax,cohort,title in zip(axes[:2],['discovery','external'],['A  Discovery: 8 donors','B  Second source: 6 donors']):
  g=paired[(paired.cohort==cohort)&(paired.selection=='primary')]
  for j,r in enumerate(g.itertuples()):ax.plot([0,1],[r.M5,r.M14],'-o',alpha=.75,lw=1.2,ms=4,label=r.donor)
  ax.axhline(0,color='.5',ls=':',lw=1);ax.set_xticks([0,1],['M5','M14']);ax.set_xlim(-.2,1.2);ax.set_ylabel('PIEZO1 − TRPV4 association (Fisher z)');ax.set_title(title,loc='left',fontsize=11);ax.legend(fontsize=7,ncol=2,loc='lower left')
 lim=[min(ax.get_ylim()[0] for ax in axes[:2]),max(ax.get_ylim()[1] for ax in axes[:2])]
 for ax in axes[:2]:ax.set_ylim(lim)
 ax=axes[2];ax.set_title('C  Paired change of definition',loc='left',fontsize=11)
 for y,cohort in enumerate(['discovery','external']):
  vals=paired[(paired.cohort==cohort)&(paired.selection=='primary')].D.to_numpy();s=summary[(summary.cohort==cohort)&(summary.selection=='primary')&(summary.endpoint=='D')].iloc[0]
  ax.scatter(vals,np.linspace(y-.12,y+.12,len(vals)),s=22,color=['#21618c','#a04000'][y],alpha=.8)
  ax.errorbar(s['mean'],y+.23,xerr=[[s['mean']-s.ci95_low],[s.ci95_high-s['mean']]],fmt='D',color='black',capsize=3)
 ax.axvline(0,color='.5',ls=':',lw=1);ax.set_yticks([0,1],['Discovery','Second source']);ax.set_xlabel('D = Δ14 − Δ5 (Fisher z)');ax.set_ylim(-.4,1.55)
 fig.text(.5,.01,'Each line/dot represents one donor. Black diamonds: mean with unadjusted t-based 95% CI. Shared genes: 3.',ha='center',fontsize=9)
 fig.tight_layout(rect=[0,.055,1,1]);fig.savefig(P/'figures/Figure_7_Matrix_Definition_Bridge.png',dpi=240);fig.savefig(P/'figures/Figure_7_Matrix_Definition_Bridge.pdf');plt.close(fig)
 run=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),plan_sha256=receipt['plan_sha256'],code_sha256=sha(__file__),all_samples_estimable=bool(sample.status.eq('ESTIMABLE').all()),M5_replay_max_error=max(r['error'] for r in replay),dictionary_sha256=sha(dictionary),no_new_background_analysis=True)
 (P/'results/RUN_RECEIPT.json').write_text(json.dumps(run,indent=2));print(summary.to_string(index=False));print(test.to_string(index=False));print(json.dumps(run,indent=2))
if __name__=='__main__':main()
