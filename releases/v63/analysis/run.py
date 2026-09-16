"""Fixed, nested cell bootstrap. No donor resampling or new hypothesis tests."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import argparse,hashlib,json,sys,datetime
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def estimate(Y,F):
    n=len(Y);det=(Y>0).sum(0);reasons=[]
    if n<200:reasons.append('CELLS_LT_200')
    for g,d in zip(['PIEZO1','TRPV4','PRG4'],det):
        if d<20:reasons.append(g+'_DETECTED_LT_20')
    r=rankdata(Y,axis=0)
    if (r.var(0)==0).any():reasons.append('CONSTANT_RANK')
    results=[]
    for mode in ['unadjusted','cubic_depth']:
        bad=list(reasons);a=r.copy()
        if mode=='cubic_depth':
            x=np.log(F);sd=x.std()
            if sd<=0:bad.append('DEPTH_VARIANCE')
            else:
                x=(x-x.mean())/sd;D=np.column_stack([np.ones(n),x,x*x,x*x*x])
                if np.linalg.matrix_rank(D)!=4:bad.append('DESIGN_RANK')
                else:a=r-D@np.linalg.lstsq(D,r,rcond=None)[0]
        if (a.var(0)<=0).any():bad.append('RESIDUAL_VARIANCE')
        rr=dict(mode=mode,status='|'.join(bad) or 'ESTIMABLE',cells=n,PIEZO1_detected=int(det[0]),TRPV4_detected=int(det[1]),PRG4_detected=int(det[2]),PIEZO1_z=np.nan,TRPV4_z=np.nan,delta=np.nan)
        if not bad:
            a-=a.mean(0);den=np.sqrt((a*a).sum(0));rho=(a[:,:2].T@a[:,2])/(den[:2]*den[2])
            if not np.isfinite(rho).all() or (abs(rho)>=1).any():rr['status']='CORRELATION'
            else:
                z=np.arctanh(rho);rr.update(PIEZO1_z=z[0],TRPV4_z=z[1],delta=z[0]-z[1])
        results.append(rr)
    return results

def sample_job(m,seed,draws,out):
    z=np.load(P/'inputs'/(m['GSM']+'.npz'),allow_pickle=False)
    genes=z['genes'].tolist();A=z['expression'][:,[genes.index(g) for g in ['PIEZO1','TRPV4','PRG4']]]
    F=z['full_UMI'];isrep=z['states']=='RepC';rep=np.flatnonzero(isrep);other=np.flatnonzero(~isrep)
    assert str(z['donor'])==m['donor_id'] and str(z['region'])==m['region']
    assert np.allclose(z['expression'],np.log1p(z['counts']*1e4/F[:,None]),rtol=0,atol=1e-12)
    records=[];saved={};rng=np.random.Generator(np.random.PCG64(seed))
    for b in range(-1,draws):
        ri=rep if b<0 else rng.choice(rep,len(rep),replace=True)
        oi=other if b<0 else rng.choice(other,len(other),replace=True)
        idx=np.concatenate([ri,oi])
        assert isrep[ri].all() and not isrep[oi].any() and np.array_equal(idx[isrep[idx]],ri)
        if b<5 and b>=0:saved[f'overall_{b}']=idx;saved[f'repc_{b}']=ri
        for population,ii in [('ALL_CELLS',idx),('RepC',ri)]:
            distinct=len(np.unique(ii))
            for r in estimate(A[ii],F[ii]):records.append(dict(GSM=m['GSM'],donor=m['donor_id'],region=m['region'],draw=b,population=population,distinct_original_cells=distinct,**r))
    df=pd.DataFrame(records);df.to_csv(out/'sample_draws'/(m['GSM']+'.csv'),index=False)
    np.savez_compressed(out/'verification_indices'/(m['GSM']+'.npz'),**saved)
    return m['GSM'],df.status.value_counts().to_dict()

def summary_row(values,observed,attempted):
    v=np.asarray(values,float);v=v[np.isfinite(v)];valid=len(v)
    same=int((v*np.sign(observed)>0).sum()) if observed!=0 else int((v==0).sum())
    lo,med,hi=np.quantile(v,[.025,.5,.975]) if valid else [np.nan]*3
    return dict(observed=observed,attempted=attempted,valid=valid,invalid=attempted-valid,p025=lo,median=med,p975=hi,same_sign=same,same_sign_valid=same/valid if valid else np.nan,same_sign_all=same/attempted,unresolved_upper=(same+attempted-valid)/attempted,max_mcse_valid=.5/np.sqrt(valid) if valid else np.nan)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=P/'results');args=ap.parse_args();O=args.output
    assert not O.exists(),'Use a fresh directory.'
    O.mkdir(parents=True);(O/'sample_draws').mkdir();(O/'verification_indices').mkdir()
    freeze=json.loads((P/'FREEZE_RECEIPT.json').read_text());assert sha(P/'PLAN.md')==freeze['plan_sha256'];assert sha(P/'INPUT_MANIFEST.json')==freeze['input_manifest_sha256']
    for f in json.loads((P/'INPUT_MANIFEST.json').read_text()):assert sha(P/f['path'])==f['sha256']
    meta=pd.read_csv(P/'inputs/metadata.csv');meta=meta[meta.donor_id.isin(freeze['donors'])].sort_values('GSM')
    assert len(meta)==12 and meta.GSM.is_unique and not meta.duplicated(['donor_id','region']).any()
    seeds=np.random.SeedSequence(freeze['seed']).spawn(12)
    with ProcessPoolExecutor(max_workers=4) as ex:
        jobs=[ex.submit(sample_job,m,seeds[i],freeze['draws'],O) for i,m in enumerate(meta.to_dict('records'))]
        for job in as_completed(jobs):print(job.result(),flush=True)
    S=pd.concat([pd.read_csv(f) for f in sorted((O/'sample_draws').glob('*.csv'))],ignore_index=True)
    baseline=S[S.draw==-1];old=pd.read_csv(P/'reference_v61/sample_estimates.csv')
    chk=baseline.merge(old,on=['GSM','donor','region','population','mode'],suffixes=('_new','_old'),validate='one_to_one')
    errs=[float(abs(chk[k+'_new']-chk[k+'_old']).max()) for k in ['PIEZO1_z','TRPV4_z','delta']]
    assert len(chk)==48 and max(errs)<1e-10
    ds=[]
    for keys,x in S.groupby(['donor','draw','population','mode'],sort=True):
        assert len(x)==2 and set(x.region)=={'WB','NWB'}
        valid=x.status.eq('ESTIMABLE').all();r=dict(zip(['donor','draw','population','mode'],keys));r.update(valid=bool(valid),failure=';'.join(x.loc[x.status!='ESTIMABLE','region']+':'+x.loc[x.status!='ESTIMABLE','status']))
        for k in ['PIEZO1_z','TRPV4_z','delta']:r[k]=float(x[k].mean()) if valid else np.nan
        ds.append(r)
    D=pd.DataFrame(ds);D.to_csv(O/'donor_draws.csv',index=False)
    summary=[]
    for (donor,pop,mode),x in D.groupby(['donor','population','mode']):
        for k in ['PIEZO1_z','TRPV4_z','delta']:
            summary.append(dict(donor=donor,population=pop,mode=mode,quantity=k,**summary_row(x.loc[x.draw>=0,k],x.loc[x.draw==-1,k].item(),freeze['draws'])))
    pd.DataFrame(summary).to_csv(O/'donor_summary.csv',index=False)
    group=[]
    for (b,pop,mode),x in D.groupby(['draw','population','mode']):
        assert set(x.donor)==set(freeze['donors']);valid=bool(x.valid.all());r=dict(draw=b,population=pop,mode=mode,valid=valid)
        for k in ['PIEZO1_z','TRPV4_z','delta']:r[k]=float(x[k].mean()) if valid else np.nan
        group.append(r)
    G=pd.DataFrame(group);G.to_csv(O/'group_draws.csv',index=False)
    gs=[]
    for (pop,mode),x in G.groupby(['population','mode']):
        for k in ['PIEZO1_z','TRPV4_z','delta']:gs.append(dict(population=pop,mode=mode,quantity=k,**summary_row(x.loc[x.draw>=0,k],x.loc[x.draw==-1,k].item(),freeze['draws'])))
    pd.DataFrame(gs).to_csv(O/'group_summary.csv',index=False)
    pairrows=[];pairs=[]
    both=pd.concat([D[['donor','draw','population','mode','delta']],G.assign(donor='FIXED_SIX_MEAN')[['donor','draw','population','mode','delta']]])
    for (donor,mode),x in both.groupby(['donor','mode']):
        w=x.pivot(index='draw',columns='population',values='delta');diff=w.RepC-w.ALL_CELLS;v=w.loc[w.index>=0].dropna();n=len(v);event=int(((v.ALL_CELLS>0)&(v.RepC<0)).sum())
        pairrows.append(dict(donor=donor,mode=mode,quantity='RepC_minus_overall',**summary_row(diff.loc[diff.index>=0],diff.loc[-1],freeze['draws']),joint_order_count=event,joint_order_fraction_valid=event/n if n else np.nan,joint_order_fraction_all=event/freeze['draws'],joint_order_unresolved_upper=(event+freeze['draws']-n)/freeze['draws']))
        for b,row in w.iterrows():pairs.append(dict(donor=donor,mode=mode,draw=b,overall=row.ALL_CELLS,RepC=row.RepC,RepC_minus_overall=row.RepC-row.ALL_CELLS))
    pd.DataFrame(pairrows).to_csv(O/'paired_summary.csv',index=False);pd.DataFrame(pairs).to_csv(O/'paired_draws.csv',index=False)
    S[S.draw>=0].groupby(['GSM','donor','region','population','mode','status']).size().rename('draws').reset_index().to_csv(O/'sample_support_failures.csv',index=False)
    receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),plan_sha256=sha(P/'PLAN.md'),code_sha256=sha(Path(__file__)),input_manifest_sha256=sha(P/'INPUT_MANIFEST.json'),baseline_estimates_checked=len(chk),baseline_max_error=max(errs),draws=freeze['draws'],seed=freeze['seed'],donors=freeze['donors'],sample_draw_estimates=int((S.draw>=0).sum()),new_tests=0,python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
    (O/'RUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(pd.DataFrame(gs).query("quantity=='delta'").to_string(index=False));print(pd.DataFrame(pairrows).to_string(index=False))
if __name__=='__main__':main()
