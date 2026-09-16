from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys
import numpy as np,pandas as pd
from scipy import stats
PK=Path(__file__).resolve().parents[1]
OLD=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else PK/'replay_output'
R=PK/'legacy'
DATA=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else PK/'review_replay_output'
DATA.mkdir(parents=True,exist_ok=True)
M='COL2A1 SOX9 ACAN CHAD HAPLN1'.split();H='COL10A1 RUNX2 IBSP ALPL MMP13 SPP1 POSTN'.split();G=['PIEZO1','TRPV4']+M+H
def save(df,name):df.to_csv(DATA/name,index=False)
# R2.1: matched donors, absolute z and back-transformed donor r.
s=pd.read_csv(OLD/'sample_estimates.csv');d=pd.read_csv(OLD/'donor_estimates.csv');desc=pd.read_csv(OLD/'sample_descriptives.csv')
sd=s[s.draw.eq(-1)].merge(desc[['cohort','donor','sample','joint_positive','n_cells']],on=['cohort','donor','sample'],validate='many_to_one',suffixes=('','_primary'))
sd['joint_positive_fraction']=sd.joint_positive/sd.n_cells_primary
save(sd,'sample_absolute_correlations_and_admission.csv')
ad=[]
for (co,dn,mo),v in s[s.draw.eq(-1)].groupby(['cohort','donor','mode'],sort=False):
 ok=len(v)==(2 if co=='discovery' else 1) and v.status.eq('PASS').all()
 row=dict(cohort=co,donor=dn,mode=mo,status='PASS' if ok else 'INCOMPLETE_SAMPLE_SET',n_cells=int(v.n_cells.sum()))
 for ch in G[:2]:
  z=float(np.arctanh(v[ch+'_rho']).mean()) if ok else np.nan
  row[ch+'_mean_z']=z;row[ch+'_donor_r']=np.tanh(z)
 row['contrast_z']=row['PIEZO1_mean_z']-row['TRPV4_mean_z'];ad.append(row)
ad=pd.DataFrame(ad);save(ad,'donor_absolute_correlations.csv')
matched=[]
for co in ['discovery','external']:
 ids=ad.loc[ad.cohort.eq(co)&ad['mode'].eq('joint_positive')&ad.status.eq('PASS'),'donor']
 for mo in ['all_cells','all_cells_depth_cubic','joint_positive','joint_positive_depth_cubic']:
  v=ad[ad.cohort.eq(co)&ad.donor.isin(ids)&ad['mode'].eq(mo)]
  matched.append(dict(cohort=co,mode=mo,n=len(v),donors=';'.join(ids),mean_contrast_z=v.contrast_z.mean(),PIEZO1_mean_z=v.PIEZO1_mean_z.mean(),TRPV4_mean_z=v.TRPV4_mean_z.mean(),PIEZO1_backtransformed_mean_z=np.tanh(v.PIEZO1_mean_z.mean()),TRPV4_backtransformed_mean_z=np.tanh(v.TRPV4_mean_z.mean())))
save(pd.DataFrame(matched),'matched_donor_summary.csv')

# R2.2: retain all draws, distinguish MC error from donor sampling uncertainty.
draw=pd.read_csv(OLD/'donor_estimates_all_draws.csv');draw=draw[draw['mode'].eq('depth1000_thinned')].copy()
assert draw.status.eq('PASS').all()
baseline=d[d['mode'].eq('depth1000_unthinned')][['cohort','donor','contrast_z']].rename(columns={'contrast_z':'unthinned_z'})
draw=draw.merge(baseline,on=['cohort','donor'],validate='many_to_one');draw['paired_change']=draw.contrast_z-draw.unthinned_z
draw=draw.sort_values(['cohort','donor','draw']);draw['cumulative_mean']=draw.groupby(['cohort','donor']).contrast_z.transform(lambda a:a.expanding().mean())
save(draw,'all_thinning_draws_and_cumulative_means.csv')
mc=[]
for (co,dn),v in draw.groupby(['cohort','donor']):
 assert len(v)==20
 mc.append(dict(cohort=co,donor=dn,n_draws=len(v),unthinned_z=v.unthinned_z.iloc[0],thinned_mean_z=v.contrast_z.mean(),paired_change=v.paired_change.mean(),MC_min=v.contrast_z.min(),MC_max=v.contrast_z.max(),MC_sd=v.contrast_z.std(),MC_SE=v.contrast_z.std()/np.sqrt(20),negative_draws=int((v.contrast_z<0).sum()),cumulative_10=v.cumulative_mean.iloc[9],cumulative_20=v.cumulative_mean.iloc[19]))
mc=pd.DataFrame(mc);save(mc,'thinning_paired_changes_and_MC.csv')

# Exact replay of the existing RNG draws solely to add omitted diagnostics.
datasets=[];meta=pd.read_csv(R/'audit/p07_prefreeze/candidate_design_16x9.tsv',sep='\t')
for m in meta.itertuples():
 x=pd.read_csv(R/f'audit/p07_execution/E3_cell_input_{m.legacy_id}.tsv.gz',sep='\t');c=pd.read_csv(R/f'audit/p09_exploration/channel_counts_{m.legacy_id}.tsv.gz',sep='\t');assert x.cell_id.equals(c.cell_id)
 datasets.append(dict(cohort='discovery',donor=m.donor_id,sample=m.GSM,region=m.region,ct=np.column_stack([c[G[:2]].to_numpy(),x[M+H].to_numpy()]),full=x.total_counts.to_numpy(),ids=x.cell_id.to_numpy()))
cross=pd.read_csv(R/'audit/p19_admission/donor_library_run_crosswalk.tsv',sep='\t');cross=cross[cross.GSE.eq('GSE220243')]
for m in cross.itertuples():
 x=pd.read_csv(PK/'verified_external_counts'/f'{m.GSM}_counts.tsv.gz',sep='\t');x=x[x.primary_P19]
 datasets.append(dict(cohort='external',donor=m.donor_key,sample=m.GSM,region='sample',ct=x[G].to_numpy(),full=x.full_UMI.to_numpy(),ids=x.barcode.to_numpy()))
diag=[];maps=[];maxerr=0.
for di,ds in enumerate(datasets):
 key={k:ds[k] for k in ['cohort','donor','sample','region']};ct=ds['ct'];full=ds['full'];el=full>=1000;ct=ct[el].astype('int64');f=full[el].astype('int64')
 maps.append(dict(**key,sample_index=di,primary_cells=len(full),depth_eligible_cells=len(f),selected_cell_ids_in='legacy/audit/p07_execution/' if ds['cohort']=='discovery' else 'verified_external_counts/'))
 rng=np.random.default_rng(np.random.SeedSequence([20260912,di]));prob=1000/f;remain=f-ct.sum(1)
 for rep in range(-1,20):
  if rep==-1:thin=ct;total=f
  else:thin=rng.binomial(ct,prob[:,None]);total=thin.sum(1)+rng.binomial(remain,prob)
  val=np.log1p(1e4*thin/total[:,None]);score=val[:,2:7].mean(1)
  rr=[float(stats.spearmanr(val[:,j],score).statistic) for j in [0,1]];delta=np.arctanh(rr[0])-np.arctanh(rr[1])
  old=s[s['sample'].eq(ds['sample'])&s.draw.eq(rep)&s['mode'].eq('depth1000_unthinned' if rep==-1 else 'depth1000_thinned')]
  assert len(old)==1;err=abs(delta-old.contrast_z.iloc[0]);assert err<1e-12;maxerr=max(maxerr,float(err))
  diag.append(dict(**key,draw=rep,n_cells=len(f),depth_mean=total.mean(),depth_median=np.median(total),depth_sd=total.std(ddof=1),depth_min=total.min(),depth_max=total.max(),PIEZO1_detection=(thin[:,0]>0).mean(),TRPV4_detection=(thin[:,1]>0).mean(),joint_positive_fraction=(thin[:,:2]>0).all(1).mean(),M_all_zero_fraction=(thin[:,2:7]==0).all(1).mean(),contrast_z=delta))
 print('REPLAY DIAGNOSTICS',ds['sample'],flush=True)
save(pd.DataFrame(diag),'thinning_before_after_measurement_diagnostics.csv');save(pd.DataFrame(maps),'sample_order_and_mapping.csv')
pd.DataFrame(diag).groupby(['cohort','draw']).agg(n_samples=('sample','size'),mean_sample_depth=('depth_mean','mean'),mean_sample_PIEZO1_detection=('PIEZO1_detection','mean'),mean_sample_TRPV4_detection=('TRPV4_detection','mean'),mean_sample_joint_positive=('joint_positive_fraction','mean'),mean_sample_M_zero=('M_all_zero_fraction','mean')).reset_index().to_csv(DATA/'diagnostic_sample_mean_by_draw.csv',index=False)

