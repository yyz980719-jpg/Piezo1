from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import json
P=Path(__file__).resolve().parent
F=P/'figures';Q=P/'qa';Q.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'svg.fonttype':'none','pdf.fonttype':42})
blue='#397FB5';pink='#C46498';ink='#263648';grey='#89949E'
fixed=pd.read_csv(P/'analysis/reference_v61/donor_estimates.csv')
don=pd.read_csv(P/'analysis/results/donor_summary.csv')
grp=pd.read_csv(P/'analysis/results/group_summary.csv')
paired=pd.read_csv(P/'analysis/results/paired_summary.csv')
ids=['OA1','OA2','OA4','OA5','OA7','OA8']
fig=plt.figure(figsize=(183/25.4,173/25.4))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.8],left=.15,right=.98,bottom=.17,top=.94,hspace=.62,wspace=.38)
source=[]
for k,mode in enumerate(['unadjusted','cubic_depth']):
    ax=fig.add_subplot(gs[0,k]);d=fixed[fixed['mode'].eq(mode)]
    for j,ch in enumerate(['PIEZO1','TRPV4']):
        co=[blue,pink][j];x=np.array([0,1])+j*2.6
        w=d.pivot(index='donor',columns='population',values=ch+'_z').loc[ids,['ALL_CELLS','RepC']]
        for donor,r in w.iterrows():
            ax.plot(x,r.to_numpy(),color=co,lw=.6,alpha=.5,marker='o',ms=2.5)
            for pop,val in r.items():source.append(dict(panel=chr(97+k),donor=donor,mode=mode,population=pop,quantity=ch+'_z',observed=val))
        ax.plot(x,w.mean(),color=co,lw=1.6,marker='D',ms=4)
        ax.text(x.mean(),.205,ch,ha='center',color=co,fontsize=7)
    ax.axhline(0,color='#CDD3D9',lw=.7,zorder=0);ax.set_ylim(-.15,.235);ax.set_xlim(-.45,4.05)
    ax.set_xticks([0,1,2.6,3.6],['Overall','RepC','Overall','RepC']);ax.set_ylabel('Absolute association (Fisher z)')
    ax.set_title('Unadjusted' if k==0 else 'Cubic depth adjusted',loc='left',pad=14)
    ax.text(-.20,1.12,chr(97+k),transform=ax.transAxes,fontweight='bold',fontsize=9)
def getrow(donor,mode,ispair):
    if ispair:q=paired[paired.donor.eq(donor)&paired['mode'].eq(mode)]
    else:
        q=grp if donor=='FIXED_SIX_MEAN' else don[don.donor.eq(donor)]
        q=q[q.population.eq('RepC')&q['mode'].eq(mode)&q.quantity.eq('delta')]
    assert len(q)==1
    return q.iloc[0]
for j in range(2):
    ax=fig.add_subplot(gs[1,j]);units=ids+['FIXED_SIX_MEAN']
    for i,unit in enumerate(units):
        for k,mode in enumerate(['unadjusted','cubic_depth']):
            r=getrow(unit,mode,j==1);y=i+(-.13 if k==0 else .13)
            co=ink if k==0 else grey;mark='o' if k==0 else 's'
            ax.plot([r.p025,r.p975],[y,y],color=co,lw=1.2 if i==6 else .8)
            ax.scatter(r.observed,y,s=22 if i==6 else 15,marker=mark,facecolor=co if k==0 else 'white',edgecolor=co,linewidth=.8,zorder=3)
            source.append(dict(panel=chr(99+j),donor=unit,mode=mode,population='RepC' if j==0 else 'RepC_minus_overall',quantity='delta',observed=r.observed,p025=r.p025,p975=r.p975,valid=r.valid,attempted=r.attempted))
    ax.axvline(0,color='#AAB3BC',lw=.7,ls='--',zorder=0);ax.axhline(5.5,color='#D9DEE2',lw=.65)
    ax.set_yticks(range(7),[x if x!='FIXED_SIX_MEAN' else 'Fixed-six mean' for x in units] if j==0 else [f'{1254 if x in ["OA1","FIXED_SIX_MEAN"] else 2000}/2000' for x in units])
    ax.set_ylim(6.55,-.7);ax.set_xlim((-.22,.12) if j==0 else (-.36,.10))
    ax.set_xlabel('RepC contrast Δ (Fisher z)' if j==0 else 'RepC − overall contrast (Fisher z)')
    ax.set_title('RepC contrast stability' if j==0 else 'Paired change within donor',loc='left',pad=14)
    ax.text(-.20,1.07,chr(99+j),transform=ax.transAxes,fontweight='bold',fontsize=9)
    if j==1:ax.text(-.04,1.005,'Eligible draws',transform=ax.transAxes,ha='right',fontsize=6.5)
fig.legend(handles=[Line2D([0],[0],color=ink,marker='o',lw=1,label='Unadjusted'),Line2D([0],[0],color=grey,marker='s',markerfacecolor='white',lw=1,label='Depth adjusted')],loc='lower center',bbox_to_anchor=(.56,.093),ncol=2,frameon=False,fontsize=7)
fig.text(.15,.065,'Ranges: 2.5th–97.5th percentiles of eligible cell-resampling draws.',fontsize=7)
fig.text(.15,.042,'746/2000 attempts (37.3%) failed the OA1 PRG4 detection gate.',fontsize=7)
fig.text(.15,.019,'Conditional on fixed cells and labels; not donor-population confidence intervals.',fontsize=6.5)
for ext in ['pdf','svg','png','tiff']:
    kwargs={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(F/f'Figure_5.{ext}',dpi=600,**kwargs)
fig.savefig(Q/'Figure_5_preview.png',dpi=170)
pd.DataFrame(source).to_csv(P/'source_data/Figure_5_source.csv',index=False)
assert sum(getrow(x,'unadjusted',False).p025<0<getrow(x,'unadjusted',False).p975 for x in ids)==4
assert sum(getrow(x,'cubic_depth',False).p025<0<getrow(x,'cubic_depth',False).p975 for x in ids)==4
(P/'audit/figure5_checks.json').write_text(json.dumps(dict(donors=6,draws=2000,valid_complete=1254,failed=746,donor_ranges_crossing_zero_each_mode=4,source_rows=len(source),new_analysis=False),indent=2),encoding='utf-8')
print('Figure 5 and exact source rows exported')
