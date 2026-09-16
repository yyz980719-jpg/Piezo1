from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parent;A=P/'analysis';O=A/'results'
checked=0;maxerr=0;gates=0
for f in sorted((O/'verification_indices').glob('*.npz')):
    ix=np.load(f);z=np.load(A/'inputs'/f.name);genes=z['genes'].tolist();F=z['full_UMI'];Y=np.log1p(z['counts'][:,[genes.index(g) for g in ['PIEZO1','TRPV4','PRG4']]]*1e4/F[:,None]);label=z['states']=='RepC'
    saved=pd.read_csv(O/'sample_draws'/(f.stem+'.csv'))
    for b in range(5):
        allix=ix[f'overall_{b}'];repix=ix[f'repc_{b}'];assert np.array_equal(allix[label[allix]],repix)
        for pop,idx in [('ALL_CELLS',allix),('RepC',repix)]:
            X=Y[idx];R=pd.DataFrame(X).rank(method='average').to_numpy();depth=np.log(F[idx]);depth=(depth-depth.mean())/depth.std()
            for mode in ['unadjusted','cubic_depth']:
                row=saved[(saved.draw==b)&(saved.population==pop)&(saved['mode']==mode)].iloc[0]
                gate=((X>0).sum(0)<20).any() or len(X)<200
                if gate:assert row.status!='ESTIMABLE';gates+=1;continue
                assert row.status=='ESTIMABLE'
                if mode=='cubic_depth':
                    coef=np.polynomial.polynomial.polyfit(depth,R,3)
                    resid=R-np.polynomial.polynomial.polyval(depth,coef).T
                else:resid=R
                rr=np.corrcoef(resid.T);zz=np.arctanh([rr[0,2],rr[1,2]]);expected=np.array([zz[0],zz[1],zz[0]-zz[1]])
                err=float(np.max(abs(expected-row[['PIEZO1_z','TRPV4_z','delta']].to_numpy(float))))
                assert err<1e-10;maxerr=max(maxerr,err);checked+=1
replay=P/'qa/replay';same=[]
for f in sorted(O.rglob('*.csv')):
    other=replay/f.relative_to(O);assert other.exists()
    ok=hashlib.sha256(f.read_bytes()).digest()==hashlib.sha256(other.read_bytes()).digest();assert ok
    same.append(f.relative_to(O).as_posix())
receipt=dict(independent_reconstructed_estimates=checked,independent_gate_failures=gates,max_absolute_error=maxerr,nesting_checked_samples=12,nesting_checked_draws_per_sample=5,repeated_identical_csvs=same,scope='Internal computational crosscheck, not external independent replication.')
(P/'audit/bootstrap_verification.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt,indent=2))
