from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parents[1]
FIG=BASE/'figures'; FIG.mkdir(exist_ok=True)

# Figure: weak-feedback criterion and nonlinear validation
wf=pd.read_csv(BASE/'data'/'weak_feedback_sensitivity.csv')
sim=pd.read_csv(BASE/'data'/'sde_smallnoise_summary.csv')
sel=pd.read_csv(BASE/'data'/'selected_cases.csv')

fig,axs=plt.subplots(1,2,figsize=(10.6,4.2))
ax=axs[0]
x=wf['phi']/np.pi
ax.plot(x,wf['g_phase'],label=r'linewidth sensitivity $g_{\phi}$')
ax.plot(x,wf['g_RIN'],label=r'RIN sensitivity $g_{\rm RIN}$')
mask=(wf['g_phase']<0)&(wf['g_RIN']<0)
ax.fill_between(x,0,1,where=mask,transform=ax.get_xaxis_transform(),alpha=.16,label='first-order simultaneous suppression')
ax.axhline(0,lw=.9,ls='--',c='k')
ax.set_xlim(-.5,.5); ax.set_xlabel(r'feedback phase $\phi/\pi$'); ax.set_ylabel(r'first-order logarithmic sensitivity')
ax.set_title('(a) Weak-feedback phase criterion'); ax.legend(fontsize=8,loc='lower left')

ax=axs[1]
metrics=[('D_ratio_sim','D_ratio_lin','linewidth'),('RIN_ratio_sim','RIN_ratio_lin','RIN'),('N_ratio_sim','N_ratio_lin','carrier')]
cases=['tradeoff','win_win']; labels=['Trade-off','Win-win']
xx=np.arange(len(metrics)); width=.18
for j,(case,label) in enumerate(zip(cases,labels)):
    row=sim.loc[sim.case==case].iloc[0]
    base=sim.loc[sim.case=='solitary'].iloc[0]
    vals=[row[m[0]] for m in metrics]
    # approximate propagated 2SE for ratios; phase SE and stationary variances
    errs=[]
    errs.append(2*row.Dse/base.Dhat)
    errs.append(2*row.RINse/base.RINhat)
    errs.append(2*row.Nse/base.Nhat)
    theory=[row[m[1]] for m in metrics]
    off=(-.27 if j==0 else .09)
    ax.bar(xx+off,vals,width,yerr=errs,capsize=2,label=f'{label}: nonlinear DDE')
    ax.scatter(xx+off,theory,marker='D',s=28,facecolors='white',edgecolors='black',zorder=4,label=f'{label}: linear resolvent' if j==0 else None)
ax.axhline(1,lw=.9,ls='--',c='k')
ax.set_xticks(xx, [m[2] for m in metrics]); ax.set_ylabel('ratio to solitary laser')
ax.set_title('(b) Direct stochastic Lang-Kobayashi validation')
handles,labels0=ax.get_legend_handles_labels();
# remove duplicate linear label ambiguity
ax.legend(handles,labels0,fontsize=7,loc='upper right')
fig.tight_layout()
for ext in ['png','pdf']:
    fig.savefig(FIG/f'fig3_phase_criterion_validation.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)

# Figure: practical robustness and device-scale illustration
noise=pd.read_csv(BASE/'data'/'noise_model_ratios.csv')
bw=pd.read_csv(BASE/'data'/'bandwidth_ratios.csv')
rob=pd.read_csv(BASE/'data'/'extended_robustness.csv')
dev=pd.read_csv(BASE/'data'/'device_scale_illustration.csv')
fix=pd.read_csv(BASE/'data'/'fixed_output_power.csv')
fig,axs=plt.subplots(2,2,figsize=(10.8,8.0))
ax=axs[0,0]
mods=['field_spontaneous','balanced_mixed']; mlab=['field Langevin source','field + carrier source']
xx=np.arange(2); w=.18
for i,mod in enumerate(mods):
    for j,case in enumerate(['tradeoff','win_win']):
        val=noise[(noise.noise_model==mod)&(noise.case==case)].integrated_RIN_ratio.iloc[0]
        ax.bar(i+(j-.5)*w,val,w,label=('Trade-off' if j==0 else 'Win-win') if i==0 else None)
ax.axhline(1,c='k',ls='--',lw=.9); ax.set_xticks(xx,mlab); ax.set_ylabel('integrated RIN ratio'); ax.set_title('(a) Source-model robustness'); ax.legend(fontsize=8)

ax=axs[0,1]
for case,label in [('tradeoff','Trade-off'),('win_win','Win-win')]:
    d=bw[bw.case==case]
    ax.plot(d.omega_max,d.RIN_ratio,marker='o',label=label)
ax.axhline(1,c='k',ls='--',lw=.9); ax.set_xscale('log'); ax.set_xlabel(r'detection bandwidth $\omega_{\max}$'); ax.set_ylabel('band-limited RIN ratio'); ax.set_title('(b) Detection-bandwidth robustness'); ax.legend(fontsize=8)

ax=axs[1,0]
x=np.arange(len(rob));
ax.plot(x,rob.best_Dratio,marker='o',label='linewidth ratio')
ax.plot(x,rob.best_RINratio,marker='s',label='RIN ratio')
ax.axhline(1,c='k',ls='--',lw=.9); ax.set_xticks(x,rob.configuration,rotation=35,ha='right'); ax.set_ylim(0.35,1.02); ax.set_ylabel('best stable simultaneous-suppression ratio'); ax.set_title('(c) Pump, carrier time, and gain compression'); ax.legend(fontsize=8)

ax=axs[1,1]
xx=np.arange(3)
linewidth=dev.linewidth_Hz/1e6
rin=dev.RIN_peak_dBc_per_Hz
ax.bar(xx-.18,linewidth,.36,label='Lorentzian linewidth (MHz)')
ax.set_ylabel('linewidth (MHz)'); ax.set_xticks(xx,['Solitary','Trade-off','Win-win']); ax2=ax.twinx(); ax2.plot(xx+.18,rin,marker='o',ls='none',label='RIN peak (dBc/Hz)'); ax2.set_ylabel('RIN peak (dBc/Hz)')
ax.set_title('(d) Device-scale illustration')
lines,lab=ax.get_legend_handles_labels(); lines2,lab2=ax2.get_legend_handles_labels(); ax.legend(lines+lines2,lab+lab2,fontsize=7,loc='best')
fig.tight_layout()
for ext in ['png','pdf']:
    fig.savefig(FIG/f'fig5_practical_robustness.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)

# A compact fixed-pump/fixed-output comparison for the supplement
fig,ax=plt.subplots(figsize=(6.4,3.8))
metrics=['phase_diffusion_ratio','RIN_ratio','carrier_ratio']; labs=['linewidth','RIN','carrier']
xx=np.arange(3); w=.18
for i,(case,clab) in enumerate([('tradeoff','Trade-off'),('win_win','Win-win')]):
    for j,(comp,style) in enumerate([('fixed_pump','fixed pump'),('fixed_output','fixed output')]):
        row=fix[(fix.case==case)&(fix.comparison==comp)].iloc[0]
        vals=[row[m] for m in metrics]
        ax.bar(xx+(i*2+j-1.5)*w,vals,w,label=f'{clab}, {style}')
ax.axhline(1,c='k',ls='--',lw=.9); ax.set_xticks(xx,labs); ax.set_ylabel('ratio to solitary laser'); ax.legend(fontsize=7,ncol=2); ax.set_title('Fixed-pump and fixed-output-power comparisons')
fig.tight_layout()
for ext in ['png','pdf']:
    fig.savefig(FIG/f'figS_fixed_output.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)
