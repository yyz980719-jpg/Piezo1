"""Frozen six-donor PRG4/RepC diagnostic, no hypothesis tests or new selection."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse, json, hashlib, datetime, sys
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, spearmanr
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=P/'results');O=ap.parse_args().output
O.mkdir(parents=True,exist_ok=True)
assert not (O/'RUN_RECEIPT.json').exists(),'Use a fresh output folder.'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
freeze=json.loads((P/'FREEZE_RECEIPT.json').read_text())
assert sha(P/'PLAN.md')==freeze['plan_sha256']
assert sha(P/'INPUT_MANIFEST.json')==freeze['input_manifest_sha256']
for f in json.loads((P/'INPUT_MANIFEST.json').read_text()):assert sha(P/f['path'])==f['sha256']
I=P/'inputs';R=P/'reference';donors=freeze['donors']
assert donors==['OA1','OA2','OA4','OA5','OA7','OA8']
meta=pd.read_csv(I/'metadata.csv');meta=meta[meta.donor_id.isin(donors)].copy()
assert len(meta)==12 and meta.GSM.is_unique and not meta.duplicated(['donor_id','region']).any()
old=pd.read_csv(R/'V51_state_region_estimates.csv')
olddep=pd.read_csv(R/'V50_member_sample_scenarios.csv')
oldctx=pd.read_csv(R/'V51_state_expression_context.csv')
oldmatch=pd.read_csv(R/'V51_matched_donor_comparisons.csv')
admission=pd.read_csv(R/'V51_state_common_admission.csv')
assert admission.loc[admission.state.eq('RepC'),'common_donors'].item().split(';')==donors
manifest={x['name']:x['sha256'] for x in json.loads((I/'V51_INPUT_MANIFEST.json').read_text())['files']}
rows=[];contexts=[];joint=[];numerics=[];replays=[];identities=[];cellset=set()

def corr(u,v):
    u=u-u.mean();v=v-v.mean();d=np.sqrt(u@u*(v@v))
    return float(u@v/d) if d>0 else np.nan

def distribution(x,prefix):
    out={prefix+'mean':float(np.mean(x)),prefix+'sd':float(np.std(x,ddof=1)),prefix+'distinct':len(np.unique(x))}
    for name,q in [('min',0),('q10',.1),('q25',.25),('median',.5),('q75',.75),('q90',.9),('max',1)]:out[prefix+name]=float(np.quantile(x,q))
    return out

for m in meta.itertuples():
    f=I/(m.GSM+'.npz');assert sha(f)==manifest[f.name]
    z=np.load(f,allow_pickle=False);genes=z['genes'].tolist();C=z['counts'];E=z['expression'];F=z['full_UMI'];labels=z['states'];ids=z['cell_ids'].tolist()
    assert str(z['donor'])==m.donor_id and str(z['region'])==m.region
    assert len(genes)==len(set(genes)) and len(ids)==len(set(ids)) and cellset.isdisjoint(ids)
    cellset.update(ids)
    assert np.isfinite(C).all() and (C>=0).all() and np.equal(C,np.floor(C)).all()
    assert np.isfinite(F).all() and (F>0).all() and np.equal(F,np.floor(F)).all() and (C.sum(1)<=F).all()
    assert np.isfinite(E).all() and np.allclose(E,np.log1p(C*1e4/F[:,None]),rtol=0,atol=1e-12)
    A=E[:,[genes.index(g) for g in ['PIEZO1','TRPV4','PRG4']]]
    identities.append(dict(GSM=m.GSM,donor=m.donor_id,region=m.region,cells=len(ids),RepC_cells=int((labels=='RepC').sum()),input_sha256=sha(f),cell_ids_sha256=hashlib.sha256('\n'.join(ids).encode()).hexdigest()))
    for pop in ['ALL_CELLS','RepC']:
        mask=np.ones(len(F),bool) if pop=='ALL_CELLS' else labels=='RepC';Y=A[mask];full=F[mask];n=len(Y)
        rec=dict(GSM=m.GSM,donor=m.donor_id,region=m.region,population=pop,cells=n)
        nd=(Y>0).sum(0);bad=[]
        if n<200:bad.append('CELLS_LT_200')
        for g,count in zip(['PIEZO1','TRPV4','PRG4'],nd):
            if count<20:bad.append(g+'_DETECTED_LT_20')
        if not n or np.any(Y.var(0)==0):bad.append('VARIANCE_GATE')
        assert not bad,'Frozen baseline sample unexpectedly failed: '+str(rec)+str(bad)
        for j,g in enumerate(['PIEZO1','TRPV4','PRG4']):
            v=Y[:,j];positive=v[v>0]
            contexts.append(dict(**rec,gene=g,detected=int(nd[j]),fraction=float(nd[j]/n),**distribution(v,'expression_'),**distribution(positive,'positive_')))
            prior=oldctx[(oldctx.GSM==m.GSM)&(oldctx.state==pop)&(oldctx.gene==g)]
            assert len(prior)==1
            for k,vv in [('cells',n),('detected',nd[j]),('fraction',nd[j]/n),('mean_including_zero',v.mean())]:
                err=abs(float(prior.iloc[0][k])-vv);assert err<1e-10;replays.append(dict(GSM=m.GSM,population=pop,check='V51_'+g+'_'+k,error=err))
        joint.append(dict(**rec,channel_joint_n=int(((Y[:,:2]>0).all(1)).sum()),channel_joint_fraction=float((Y[:,:2]>0).all(1).mean()),triple_n=int((Y>0).all(1).sum()),triple_fraction=float((Y>0).all(1).mean()),**distribution(full,'full_UMI_')))
        ranks=rankdata(Y,axis=0)
        for mode in ['unadjusted','cubic_depth']:
            row=dict(**rec,mode=mode,PIEZO1_detected=int(nd[0]),TRPV4_detected=int(nd[1]),PRG4_detected=int(nd[2]),status='ESTIMABLE',design_rank=0,design_condition=np.nan,PIEZO1_rho=np.nan,TRPV4_rho=np.nan,PIEZO1_z=np.nan,TRPV4_z=np.nan,delta=np.nan)
            if mode=='unadjusted':
                rho=np.array([spearmanr(Y[:,j],Y[:,2]).statistic for j in [0,1]])
                alt=np.array([corr(ranks[:,j],ranks[:,2]) for j in [0,1]])
            else:
                x=np.log(full);sd=x.std()
                if not sd>0:row['status']='DEPTH_VARIANCE_GATE';rows.append(row);continue
                x=(x-x.mean())/sd;D=np.column_stack([np.ones(n),x,x*x,x*x*x]);rank=int(np.linalg.matrix_rank(D));row.update(design_rank=rank,design_condition=float(np.linalg.cond(D)))
                if rank!=4:row['status']='DESIGN_RANK_GATE';rows.append(row);continue
                resid=ranks-D@np.linalg.lstsq(D,ranks,rcond=None)[0]
                Q=np.linalg.qr(D,mode='reduced')[0];other=ranks-Q@(Q.T@ranks)
                if np.any(np.var(resid,axis=0)<=0):row['status']='RESIDUAL_VARIANCE_GATE';rows.append(row);continue
                rho=np.array([corr(resid[:,j],resid[:,2]) for j in [0,1]])
                alt=np.array([corr(other[:,j],other[:,2]) for j in [0,1]])
                for j,g in enumerate(['PIEZO1','TRPV4','PRG4']):row[g+'_residual_sd']=float(resid[:,j].std())
            error=float(max(abs(rho-alt)));assert error<1e-10
            numerics.append(dict(**rec,mode=mode,error=error))
            if not np.isfinite(rho).all() or not (abs(rho)<1).all():row['status']='CORRELATION_GATE';rows.append(row);continue
            zz=np.arctanh(rho);row.update(PIEZO1_rho=rho[0],TRPV4_rho=rho[1],PIEZO1_z=zz[0],TRPV4_z=zz[1],delta=zz[0]-zz[1]);rows.append(row)
            prior=None
            if mode=='unadjusted':prior=old[(old.GSM==m.GSM)&(old.state==pop)&(old.target=='PRG4')]
            elif pop=='ALL_CELLS':prior=olddep[(olddep['sample']==m.legacy_id)&(olddep.member=='PRG4')&(olddep.scenario=='primary_cubic_depth')]
            if prior is not None:
                assert len(prior)==1 and prior.iloc[0].status=='ESTIMABLE'
                for k in ['PIEZO1_z','TRPV4_z','delta']:
                    err=abs(row[k]-prior.iloc[0][k]);assert err<1e-10;replays.append(dict(GSM=m.GSM,population=pop,check=mode+'_'+k,error=err))

S=pd.DataFrame(rows);S.to_csv(O/'sample_estimates.csv',index=False)
pd.DataFrame(contexts).to_csv(O/'expression_detection_ranges.csv',index=False)
pd.DataFrame(joint).to_csv(O/'depth_joint_diagnostics.csv',index=False)
pd.DataFrame(identities).to_csv(O/'identity_checks.csv',index=False)
dn=[]
for (donor,pop,mode),x in S.groupby(['donor','population','mode'],sort=True):
    assert len(x)==2 and set(x.region)=={'WB','NWB'}
    ok=x.status.eq('ESTIMABLE').all();r=dict(donor=donor,population=pop,mode=mode,eligible=bool(ok),failure='|'.join(x.loc[x.status!='ESTIMABLE','status']))
    for k in ['PIEZO1_z','TRPV4_z','delta']:r[k]=float(x[k].mean()) if ok else np.nan
    dn.append(r)
D=pd.DataFrame(dn);D.to_csv(O/'donor_estimates.csv',index=False)
for r in D.query("mode=='unadjusted'").itertuples():
    for k in ['PIEZO1_z','TRPV4_z','delta']:
        p=oldmatch[(oldmatch.state=='RepC')&(oldmatch.target=='PRG4')&(oldmatch.donor==r.donor)&(oldmatch.quantity==k)]
        assert len(p)==1;val=p.iloc[0]['within_state' if r.population=='RepC' else 'matched_overall'];err=abs(getattr(r,k)-val);assert err<1e-10
        replays.append(dict(GSM='donor:'+r.donor,population=r.population,check=k,error=err))
summ=[]
for (pop,mode),x in D.groupby(['population','mode']):
    complete=bool(x.eligible.all()) and set(x.donor)==set(donors)
    for k in ['PIEZO1_z','TRPV4_z','delta']:
        vals=x[k].to_numpy();r=dict(population=pop,mode=mode,quantity=k,n_fixed=6,n_eligible=int(x.eligible.sum()),complete_fixed_set=complete)
        r.update({n:np.nan for n in ['mean','sd','median','min','max','positive','negative']})
        if complete:r.update(mean=float(vals.mean()),sd=float(vals.std(ddof=1)),median=float(np.median(vals)),min=float(vals.min()),max=float(vals.max()),positive=int((vals>0).sum()),negative=int((vals<0).sum()))
        summ.append(r)
pd.DataFrame(summ).to_csv(O/'summary.csv',index=False)
wide=D.pivot(index=['donor','population'],columns='mode',values=['PIEZO1_z','TRPV4_z','delta'])
wide.columns=['_'.join(x) for x in wide.columns];wide=wide.reset_index()
for k in ['PIEZO1_z','TRPV4_z','delta']:wide[k+'_adjusted_minus_unadjusted']=wide[k+'_cubic_depth']-wide[k+'_unadjusted']
wide.to_csv(O/'matched_adjustment_changes.csv',index=False)
context=pd.DataFrame(contexts);dc=context.groupby(['donor','population','gene'],as_index=False)[['fraction','expression_mean','expression_sd','positive_mean']].mean()
dc.to_csv(O/'donor_detection_expression.csv',index=False)
context_summary=dc.groupby(['population','gene']).agg(n=('donor','size'),mean_detection_fraction=('fraction','mean'),min_donor_detection_fraction=('fraction','min'),max_donor_detection_fraction=('fraction','max'),zero_inclusive_expression_mean=('expression_mean','mean'),positive_expression_mean=('positive_mean','mean')).reset_index()
context_summary.to_csv(O/'detection_summary.csv',index=False)
pd.DataFrame(numerics).to_csv(O/'numerical_crosschecks.csv',index=False)
pd.DataFrame(replays).to_csv(O/'historical_replay_checks.csv',index=False)
receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),plan_sha256=sha(P/'PLAN.md'),code_sha256=sha(Path(__file__)),input_manifest_sha256=sha(P/'INPUT_MANIFEST.json'),donors=donors,samples=len(meta),unique_selected_cells=len(cellset),sample_status_counts=S.status.value_counts().to_dict(),fixed_set_complete=bool(D.eligible.all()),historical_checks=len(replays),historical_max_error=max(r['error'] for r in replays),numerical_max_error=max(r['error'] for r in numerics),new_tests=0,python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
(O/'RUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(pd.DataFrame(summ).to_string(index=False));print(context_summary.to_string(index=False));print(json.dumps(receipt,indent=2))
