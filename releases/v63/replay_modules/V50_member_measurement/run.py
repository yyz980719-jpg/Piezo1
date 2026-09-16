from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,hashlib,datetime,shutil,sys,argparse
import numpy as np,pandas as pd,scipy
from scipy import stats
P=Path(__file__).resolve().parent;I=P/'inputs';I.mkdir(exist_ok=True)
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=P/'results');args=parser.parse_args();O=args.output;O.mkdir(exist_ok=True,parents=True)
assert not (O/'RUN_RECEIPT.json').exists(),'Completed run frozen; use a fresh --output directory'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(P/'PLAN.md')==json.loads((P/'FREEZE_RECEIPT.json').read_text())['plan_sha256']
V46=Path('E:/codex/2026-08-31/yu/outputs/V46_Bounded_Extension_20260915')
R=Path('E:/codex/2026-08-31/yu/outputs/PIEZO1_OA_Round1_Revision_20260908/audit/p21_external_effects')
G='COL2A1 ACAN HAPLN1 SOX9 CHAD COL9A1 COL9A2 COL9A3 COL11A1 COL11A2 COMP MATN3 PRG4 CILP CILP2 OGN'.split();M='COL2A1 SOX9 ACAN CHAD HAPLN1'.split();CH=['PIEZO1','TRPV4']
if not (I/'INPUT_MANIFEST.json').exists():
    source=[]
    def copy(f,n):
        shutil.copy2(f,I/n);source.append({'source':str(f),'copy':n,'sha256':sha(f)})
    for f in (V46/'inputs').glob('*.npz'):copy(f,f.name)
    for f,n in [(V46/'inputs/prior_sample.csv','prior_samples.csv'),(V46/'inputs/crosswalk.csv','crosswalk.csv'),(V46/'results/member_sample.csv','prior_member_samples.csv'),(V46/'inputs/measurement_samples.csv','prior_M5_measurement.csv'),(R/'donor_effects.tsv','prior_external_M5.tsv')]:copy(f,n)
    for f in sorted(I.glob('GSM*.npz')):
        src=R/(f.stem+'_normalized_targets_and_masks.tsv.gz');mask=pd.read_csv(src,sep='\t',usecols=['barcode','full_UMI','primary_P19','union']);assert mask.barcode.is_unique
        assert (~mask['union']|mask.primary_P19).all();z=np.load(f,allow_pickle=False);sub=mask.set_index('barcode').loc[z['barcodes']].reset_index();assert sub.primary_P19.all() and len(sub)==int(mask.primary_P19.sum()) and np.array_equal(sub.full_UMI,z['full_UMI'])
        sub.to_csv(I/(f.stem+'_fixed_masks.csv'),index=False);source.append({'source':str(src),'source_sha256':sha(src),'operation':'exact primary-cell barcode subset of frozen union mask; no flag change'})
    manifest={'sources':source,'files':[{'name':f.name,'sha256':sha(f)} for f in sorted(I.iterdir()) if f.is_file()]};(I/'INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
manifest=json.loads((I/'INPUT_MANIFEST.json').read_text())
for f in manifest['files']:assert sha(I/f['name'])==f['sha256']
prior=pd.read_csv(I/'prior_samples.csv');prior=prior[prior.selection.eq('primary')];records=prior.drop_duplicates('sample')[['cohort','sample','donor','region']].to_dict('records')
pm=pd.read_csv(I/'prior_member_samples.csv');ms=pd.read_csv(I/'prior_M5_measurement.csv');um=pd.read_csv(I/'prior_external_M5.tsv',sep='\t');cross=pd.read_csv(I/'crosswalk.csv').set_index('sample')
samples=[];parts=[];context=[];replay=[];checks=[];controls=[];mask_audit=[]
def corr(u,v):
    u=u-u.mean();v=v-v.mean();den=np.sqrt(np.dot(u,u)*np.dot(v,v));return float(np.dot(u,v)/den) if den>0 else np.nan
def fit(rec,target,A,y,full,depth=False,member_gate=True):
    n=len(y);det=(A>0).sum(0);td=int((y>0).sum());bad=[]
    if n<200:bad.append('CELLS_LT_200')
    if min(det)<20:bad.append('CHANNEL_DETECTION_LT_20')
    if member_gate and td<20:bad.append('MEMBER_DETECTION_LT_20')
    if n==0 or np.var(y)==0 or np.any(np.var(A,axis=0)==0):bad.append('VARIANCE_GATE')
    out=dict(**rec,member=target,cells=n,P1_detected=int(det[0]),T4_detected=int(det[1]),member_detected=td,PIEZO1_rho=np.nan,TRPV4_rho=np.nan,PIEZO1_z=np.nan,TRPV4_z=np.nan,delta=np.nan,design_rank=0)
    if bad:out['status']=';'.join(bad);return out
    ranks=np.column_stack([stats.rankdata(A[:,0]),stats.rankdata(A[:,1]),stats.rankdata(y)])
    if depth:
        x=np.log(full);sd=x.std()
        if not sd>0:out['status']='DEPTH_VARIANCE_GATE';return out
        x=(x-x.mean())/sd;design=np.column_stack([np.ones(n),x,x*x,x*x*x]);rank=int(np.linalg.matrix_rank(design));out['design_rank']=rank
        if rank!=4:out['status']='DESIGN_RANK_GATE';return out
        rr=ranks-design@np.linalg.lstsq(design,ranks,rcond=None)[0];Q=np.linalg.qr(design,mode='reduced')[0];alt=ranks-Q@(Q.T@ranks)
        rho=np.array([corr(rr[:,j],rr[:,2]) for j in range(2)]);rhob=np.array([corr(alt[:,j],alt[:,2]) for j in range(2)])
    else:
        rho=np.array([stats.spearmanr(A[:,j],y).statistic for j in range(2)]);rhob=np.array([corr(ranks[:,j],ranks[:,2]) for j in range(2)])
    err=float(np.nanmax(abs(rho-rhob)));assert err<1e-10
    checks.append(dict(**rec,member=target,check='QR_vs_lstsq' if depth else 'Pearson_on_ranks',max_error=err))
    if not (np.isfinite(rho).all() and (abs(rho)<1).all()):out['status']='RESIDUAL_VARIANCE_OR_RHO_GATE';return out
    zz=np.arctanh(rho);out.update(status='ESTIMABLE',PIEZO1_rho=rho[0],TRPV4_rho=rho[1],PIEZO1_z=zz[0],TRPV4_z=zz[1],delta=zz[0]-zz[1]);return out
def partition(rec,g,ch,u,v,group,label):
    uc=u-u.mean();vc=v-v.mean();den=np.sqrt(np.mean(uc**2)*np.mean(vc**2));total=float(np.mean(uc*vc)/den);within=[];EU=np.zeros(len(u));EV=np.zeros(len(u))
    for val in [False,True]:
        mask=group==val;p=mask.mean()
        if mask.any():
            a=u[mask];b=v[mask];EU[mask]=a.mean();EV[mask]=b.mean();within.append(float(p*np.mean((a-a.mean())*(b-b.mean()))/den))
        else:within.append(0.)
    between=float(np.mean((EU-u.mean())*(EV-v.mean()))/den);altwithin=float(np.mean((u-EU)*(v-EV))/den);err=max(abs(total-between-sum(within)),abs(altwithin-sum(within)));assert err<1e-10
    parts.append(dict(**rec,member=g,channel=ch,partition=label,total_rho=total,between_detection=between,within_nondetected=within[0],within_detected=within[1],within_total=sum(within),detected_fraction=float(group.mean()),identity_error=err))
for rec in records:
    d=np.load(I/(rec['sample']+'.npz'),allow_pickle=False);genes=d['genes'].tolist();ct=d['counts'];raw=d['expression'];full=d['full_UMI'];bars=d['barcodes'];assert len(bars)==len(set(bars)) and len(genes)==len(set(genes))
    assert (full>0).all() and np.isfinite(full).all() and (ct>=0).all() and np.equal(ct,np.floor(ct)).all() and np.allclose(raw,np.log1p(1e4*ct/full[:,None]),atol=1e-12,rtol=0)
    E=raw.copy();gsm=rec['sample']
    if rec['cohort']=='external':
        for j,g in enumerate(d['archived_genes'].tolist()):assert np.allclose(E[:,genes.index(g)],d['archived_expression'][:,j],atol=1e-12,rtol=0);E[:,genes.index(g)]=d['archived_expression'][:,j]
        ma=pd.read_csv(I/(gsm+'_fixed_masks.csv'));assert np.array_equal(ma.barcode,bars) and ma.primary_P19.all() and np.array_equal(ma.full_UMI,full);union=ma['union'].to_numpy();mask_audit.append(dict(**rec,primary_cells=len(full),union_cells=int(union.sum()),subset=True,unique_barcodes=True,full_UMI_match=True))
    else:
        info=cross.loc[gsm];assert info.donor==rec['donor'] and info.region==rec['region'];gsm=info.GSM
    A=E[:,[genes.index(x) for x in CH]]
    for g in CH+G:
        y=E[:,genes.index(g)];pos=y[y>0];ctx=dict(**rec,gene=g,cells=len(y),detected=len(pos),fraction=len(pos)/len(y),mean_including_zero=float(y.mean()),positive_mean=float(pos.mean()) if len(pos) else np.nan,positive_distinct=len(np.unique(pos)))
        for k,q in [('min',0),('q10',.1),('median',.5),('q90',.9),('max',1)]:ctx['positive_'+k]=float(np.quantile(pos,q)) if len(pos) else np.nan
        context.append(ctx)
    for g in G:
        y=E[:,genes.index(g)];base=fit(dict(**rec,scenario='primary'),g,A,y,full);samples.append(base)
        previous=pm[pm['sample'].eq(rec['sample'])&pm.member.eq(g)].iloc[0];assert base['status']==previous.status
        for c in ['PIEZO1_rho','TRPV4_rho','delta']:
            err=abs(base[c]-previous[c]);assert err<1e-10;replay.append(dict(**rec,member=g,quantity=c,error=err))
        samples.append(fit(dict(**rec,scenario='primary_cubic_depth'),g,A,y,full,depth=True))
        if rec['cohort']=='external':samples.append(fit(dict(**rec,scenario='union'),g,A[union],y[union],full[union]))
        if base['status']=='ESTIMABLE':
            v=stats.rankdata(y)
            for j,ch in enumerate(CH):
                u=stats.rankdata(A[:,j]);partition(rec,g,ch,u,v,A[:,j]>0,'channel_detected');partition(rec,g,ch,u,v,y>0,'member_detected')
    score=E[:,[genes.index(g) for g in M]].mean(1)
    for mo,depth in [('all_cells',False),('all_cells_depth_cubic',True)]:
        val=fit(dict(**rec,scenario='CONTROL_'+mo),'M5',A,score,full,depth=depth,member_gate=False);old=ms[ms['sample'].eq(gsm)&ms['mode'].eq(mo)].iloc[0]
        assert val['status']=='ESTIMABLE' and old.status=='PASS';err=abs(val['delta']-old.contrast_z);assert err<1e-10;controls.append(dict(**rec,control=mo,error=err))
    if rec['cohort']=='external':
        val=fit(dict(**rec,scenario='CONTROL_union'),'M5',A[union],score[union],full[union],member_gate=False);old=um[um.GSM.eq(gsm)&um.selection.eq('union')].iloc[0];err=abs(val['delta']-old.DeltaM);assert err<1e-10;controls.append(dict(**rec,control='M5_union',error=err))
    print('Verified',rec['cohort'],rec['sample'],flush=True)
sp=pd.DataFrame(samples);sp.to_csv(O/'member_sample_scenarios.csv',index=False);pd.DataFrame(context).to_csv(O/'primary_detection_expression_ranges.csv',index=False);pd.DataFrame(parts).to_csv(O/'detection_partition_samples.csv',index=False);pd.DataFrame(replay).to_csv(O/'baseline_replay.csv',index=False);pd.DataFrame(checks).to_csv(O/'numerical_crosschecks.csv',index=False);pd.DataFrame(controls).to_csv(O/'historical_M5_controls.csv',index=False);pd.DataFrame(mask_audit).to_csv(O/'union_mask_audit.csv',index=False)
don=[]
for (co,d,g,sc),t in sp.groupby(['cohort','donor','member','scenario']):
    ok=len(t)==(2 if co=='discovery' else 1) and t.status.eq('ESTIMABLE').all();row=dict(cohort=co,donor=d,member=g,scenario=sc,status='ESTIMABLE' if ok else 'INCOMPLETE_SAMPLE',failures='|'.join(t.loc[t.status.ne('ESTIMABLE'),'status']))
    for k in ['PIEZO1_z','TRPV4_z','delta']:row[k]=float(t[k].mean()) if ok else np.nan
    don.append(row)
dn=pd.DataFrame(don);dn.to_csv(O/'member_donor_scenarios.csv',index=False)
def summary(v,co):
    v=np.asarray(v,float);assert np.isfinite(v).all();n=len(v);mean=float(v.mean()) if n else np.nan;sd=float(v.std(ddof=1)) if n>1 else np.nan;gate=6 if co=='discovery' else 3;half=stats.t.ppf(.975,n-1)*sd/np.sqrt(n) if n>=gate else np.nan
    return dict(n=n,mean=mean,sd=sd,ci_low=mean-half,ci_high=mean+half,negative=int((v<0).sum()),positive=int((v>0).sum()),support='ADEQUATE_FOR_DESCRIPTIVE_INTERVAL' if n>=gate else 'LIMITED_SUPPORT_NO_CI')
su=[]
for (co,g,sc),t in dn.groupby(['cohort','member','scenario']):
    valid=t[t.status.eq('ESTIMABLE')]
    for k in ['PIEZO1_z','TRPV4_z','delta']:su.append(dict(cohort=co,member=g,scenario=sc,quantity=k,eligible_donors=';'.join(valid.donor),**summary(valid[k],co)))
pd.DataFrame(su).to_csv(O/'source_scenario_summaries.csv',index=False)
matched=[];mts=[];common=[]
for co,t in dn.groupby('cohort'):
    for sc in sorted(set(t.scenario)-{'primary'}):
        sets=[]
        for g in G:
            b=t[t.member.eq(g)&t.scenario.eq('primary')&t.status.eq('ESTIMABLE')].set_index('donor');s=t[t.member.eq(g)&t.scenario.eq(sc)&t.status.eq('ESTIMABLE')].set_index('donor');ids=sorted(set(b.index)&set(s.index));sets.append(set(ids))
            for k in ['PIEZO1_z','TRPV4_z','delta']:
                bv=b.loc[ids,k].to_numpy();sv=s.loc[ids,k].to_numpy()
                for d,x,y in zip(ids,bv,sv):matched.append(dict(cohort=co,member=g,scenario=sc,quantity=k,donor=d,matched_primary=x,alternative=y,paired_change=y-x))
                for what,v in [('matched_primary',bv),('alternative',sv),('paired_change',sv-bv)]:mts.append(dict(cohort=co,member=g,scenario=sc,quantity=k,estimate=what,eligible_donors=';'.join(ids),**summary(v,co)))
        ids=sorted(set.intersection(*sets));common.append(dict(cohort=co,scenario=sc,common_all16_donors=';'.join(ids),n_common=len(ids)))
pd.DataFrame(matched).to_csv(O/'matched_donor_comparisons.csv',index=False);pd.DataFrame(mts).to_csv(O/'matched_source_comparisons.csv',index=False);pd.DataFrame(common).to_csv(O/'all16_common_support.csv',index=False)
pt=pd.DataFrame(parts);pdons=[]
for keys,t in pt.groupby(['cohort','donor','member','channel','partition']):
    co,d,g,ch,part=keys;assert len(t)==(2 if co=='discovery' else 1)
    pdons.append(dict(cohort=co,donor=d,member=g,channel=ch,partition=part,**{k:float(t[k].mean()) for k in ['total_rho','between_detection','within_nondetected','within_detected','within_total']}))
dd=pd.DataFrame(pdons);dd.to_csv(O/'detection_partition_donors.csv',index=False);ps=[]
for (co,g,part),t in dd.groupby(['cohort','member','partition']):
    p=t[t.channel.eq('PIEZO1')].set_index('donor');q=t[t.channel.eq('TRPV4')].set_index('donor');ids=sorted(set(p.index)&set(q.index))
    for term in ['total_rho','between_detection','within_nondetected','within_detected','within_total']:
        for channel,values in [('PIEZO1',p.loc[ids,term]),('TRPV4',q.loc[ids,term]),('P1_minus_T4',p.loc[ids,term]-q.loc[ids,term])]:ps.append(dict(cohort=co,member=g,partition=part,channel=channel,term=term,eligible_donors=';'.join(ids),**summary(values,co)))
pd.DataFrame(ps).to_csv(O/'detection_partition_source_summaries.csv',index=False)
qa=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),plan_sha256=sha(P/'PLAN.md'),code_sha256=sha(__file__),input_manifest_sha256=sha(I/'INPUT_MANIFEST.json'),baseline_max_error=max(x['error'] for x in replay),M5_control_max_error=max(x['error'] for x in controls),crosscheck_max_error=max(x['max_error'] for x in checks),partition_identity_max_error=float(pt.identity_error.max()),sample_status_counts=sp.status.value_counts().to_dict(),new_tests=0,state_analysis=False,manuscript_modified=False,python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
(O/'RUN_RECEIPT.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa,indent=2));print(pd.DataFrame(su).query("member in ['PRG4','OGN'] and quantity == 'delta'").to_string(index=False))
