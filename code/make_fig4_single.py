import os, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
BASE=str(Path(__file__).resolve().parents[1]); FIG=os.path.join(BASE,'figures'); os.makedirs(FIG,exist_ok=True)
r=pd.read_csv(os.path.join(BASE,'data','robustness_selected.csv'))
plt.rcParams.update({'font.size':8.5,'axes.labelsize':9,'axes.titlesize':9,'legend.fontsize':7.3,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'pdf.fonttype':42,'ps.fonttype':42})
markers={10:'o',20:'s',40:'^'}; cols={2:'#56B4E9',4:'#E69F00',6:'#CC79A7'}
fig,ax=plt.subplots(figsize=(3.35,3.05),constrained_layout=True)
for _,row in r.iterrows():
    ax.scatter(row.Dphase_ratio,row.RIN_ratio,s=52,marker=markers[int(row.tau)],color=cols[int(row.alpha)],edgecolor='k',linewidth=.3)
ax.axvline(1,color='k',lw=.7,ls='--');ax.axhline(1,color='k',lw=.7,ls='--')
ax.set_xlim(.5,.9);ax.set_ylim(.45,.82);ax.set_xlabel('phase-diffusion ratio');ax.set_ylabel('integrated RIN ratio')
leg1=[Line2D([0],[0],marker=markers[t],color='none',markerfacecolor='white',markeredgecolor='k',label=f'$\\tau={t}$') for t in [10,20,40]]
leg2=[Line2D([0],[0],marker='o',color='none',markerfacecolor=cols[a],markeredgecolor='k',label=f'$\\alpha={a}$') for a in [2,4,6]]
L1=ax.legend(handles=leg1,frameon=False,loc='upper left');ax.add_artist(L1);ax.legend(handles=leg2,frameon=False,loc='lower right')
fig.savefig(os.path.join(FIG,'fig4_robustness_single.pdf'),bbox_inches='tight')
fig.savefig(os.path.join(FIG,'fig4_robustness_single.png'),bbox_inches='tight',dpi=220)
