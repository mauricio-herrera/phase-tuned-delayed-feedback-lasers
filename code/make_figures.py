import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.colors import LogNorm, TwoSlopeNorm, ListedColormap, BoundaryNorm
from lk_model import Params, solve_ecm, output_variance_grid, phase_diffusion
from lk_analysis import characteristic_roots, pole_residue
from pair_formula import modal_variance

from pathlib import Path
BASE=str(Path(__file__).resolve().parents[1])
FIG=os.path.join(BASE,'figures'); os.makedirs(FIG,exist_ok=True)
plt.rcParams.update({
    'font.size':8.5,'axes.labelsize':9,'axes.titlesize':9,'legend.fontsize':7.5,
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'figure.dpi':160,
    'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':0.8,
})
colors={'core':'#4C4C4C','trade':'#D55E00','win':'#0072B2','pair':'#009E73'}

# Selected cases
cases=[('Solitary',Params(kappa=0,tau=20,phi=0),colors['core']),
       ('Trade-off',Params(kappa=.008,tau=20,phi=0),colors['trade']),
       ('Win-win',Params(kappa=.009,tau=20,phi=np.pi/4),colors['win'])]

# FIG 1: schematic + spectra + ratios
fig=plt.figure(figsize=(7.05,5.4))
gs=fig.add_gridspec(2,2,height_ratios=[0.85,1.15],hspace=.35,wspace=.32)
ax=fig.add_subplot(gs[0,0]); ax.axis('off')
# --- schematic (redesigned: all labels sized to fit, no overflow) ---
from matplotlib.patches import Rectangle
lx0,lx1,ly0,ly1=.04,.40,.42,.74            # laser box
mx0,mx1,my0,my1=.86,.905,.38,.78           # mirror bar
laser=FancyBboxPatch((lx0,ly0),lx1-lx0,ly1-ly0,boxstyle='round,pad=.012',
                     fc='#E8EEF7',ec='black',lw=.9)
ax.add_patch(laser)
ax.text((lx0+lx1)/2,(ly0+ly1)/2+.055,'diode laser',ha='center',va='center',
        weight='bold',fontsize=8)
ax.text((lx0+lx1)/2,(ly0+ly1)/2-.065,r'$E(t),\ N(t)$',ha='center',va='center',fontsize=8)
ax.text((lx0+lx1)/2,ly0-.085,r'$\alpha,\ \epsilon,\ P$',ha='center',va='top',fontsize=7.5)
# external mirror (hatched bar)
ax.add_patch(Rectangle((mx0,my0),mx1-mx0,my1-my0,fc='#D9D9D9',ec='black',
                       lw=.9,hatch='////'))
ax.text((mx0+mx1)/2,my0-.085,'reflector',ha='center',va='top',fontsize=7.5)
# forward (output) beam
ax.add_patch(FancyArrowPatch((lx1+.017,.655),(mx0-.006,.655),arrowstyle='->',
                             mutation_scale=9,lw=1.2,color='#333333'))
ax.text((lx1+mx0)/2,.695,'output field',ha='center',va='bottom',fontsize=7.5)
# delayed return beam
ax.add_patch(FancyArrowPatch((mx0-.006,.50),(lx1+.017,.50),arrowstyle='->',
                             mutation_scale=9,lw=1.2,color='#B03A2E'))
ax.text((lx1+mx0)/2,.535,r'$\kappa\,e^{-i\phi}$',ha='center',va='bottom',
        fontsize=8,color='#B03A2E')
ax.text((lx1+mx0)/2,.455,r'delay $\tau$',ha='center',va='top',fontsize=7.5,
        color='#B03A2E')
# Langevin noise input
ax.add_patch(FancyArrowPatch((.115,.925),(.155,ly1+.02),arrowstyle='->',
                             mutation_scale=9,lw=1.0,color='#555555'))
ax.text(.045,.955,r'spontaneous-emission noise $\xi_E(t)$',ha='left',va='center',fontsize=7.2)
# observables strip
ax.text(.5,.10,r'observables:  phase diffusion $D_\varphi$   $\cdot$   RIN   $\cdot$   carrier noise',
        ha='center',va='center',fontsize=7.5)
ax.set_xlim(0,1);ax.set_ylim(0,1)
ax.set_title('(a) Delayed-feedback laser and observables',loc='left')

# calculate spectra
omega=np.concatenate([np.geomspace(1e-4,0.02,300),np.linspace(.02,.25,1400),np.geomspace(.25,3,500)])
axI=fig.add_subplot(gs[0,1])
axN=fig.add_subplot(gs[1,0])
ratios=[]
for label,p,c in cases:
    e=solve_ecm(p)
    vi,om,psdI=output_variance_grid(p,e,[2*e.A,0,0],omega=omega,return_psd=True)
    vn,_,psdN=output_variance_grid(p,e,[0,0,1],omega=omega,return_psd=True)
    rin=psdI/e.A**4
    axI.semilogy(om,rin,label=label,color=c,lw=1.5)
    axN.semilogy(om,psdN,label=label,color=c,lw=1.5)
    ratios.append((label,p,e,vi/e.A**4,vn))
axI.set_xlim(0.02,.22);axI.set_xlabel(r'angular frequency $\omega$');axI.set_ylabel('RIN spectrum (arb. units)')
axI.set_title('(b) Relative-intensity noise',loc='left');axI.legend(frameon=False)
axN.set_xlim(0.02,.22);axN.set_xlabel(r'angular frequency $\omega$');axN.set_ylabel(r'$S_N(\omega)$ (arb. units)')
axN.set_title('(c) Carrier-noise spectrum',loc='left')

axB=fig.add_subplot(gs[1,1])
core_p=cases[0][1];core_e=solve_ecm(core_p);D0,_=phase_diffusion(core_p,core_e)
core_rin=ratios[0][3];core_n=ratios[0][4]
vals=[]
for label,p,e,rin,nv in ratios[1:]:
    D,_=phase_diffusion(p,e)
    vals.append([D/D0,rin/core_rin,nv/core_n])
x=np.arange(3);w=.32
axB.bar(x-w/2,vals[0],width=w,label='Trade-off',color=colors['trade'])
axB.bar(x+w/2,vals[1],width=w,label='Win-win',color=colors['win'])
axB.axhline(1,color='k',lw=.8,ls='--')
axB.set_xticks(x,['phase\ndiffusion','integrated\nRIN','carrier\nvariance'])
axB.set_ylabel('ratio to solitary laser');axB.set_title('(d) Phase-controlled noise partition',loc='left')
axB.legend(frameon=False,ncol=2,loc='upper center')
axB.set_ylim(0,1.75)
fig.savefig(os.path.join(FIG,'fig1_model_spectra.pdf'),bbox_inches='tight')
fig.savefig(os.path.join(FIG,'fig1_model_spectra.png'),bbox_inches='tight',dpi=300)
plt.close(fig)

# FIG 2 map
D=np.load(os.path.join(BASE,'data','map_tau20_full.npz'))
k=D['kappas'];phi=D['phis'];g=k*20
X,Y=np.meshgrid(phi/np.pi,g)
st=D['margin_exact']>0
fig,axs=plt.subplots(2,2,figsize=(7.05,5.5),sharex=True,sharey=True,constrained_layout=True)
# stability
A=np.where(st,D['margin_exact'],np.nan)
im=axs[0,0].pcolormesh(X,Y,A,shading='gouraud',cmap='viridis',rasterized=True)
fig.colorbar(im,ax=axs[0,0],label='stability margin')
axs[0,0].contour(X,Y,D['margin_exact'],levels=[0],colors='k',linewidths=.9)
axs[0,0].set_title('(a) Stability of the selected ECM',loc='left')
# phase ratio log
A=np.where(st,D['Dphys'],np.nan)
im=axs[0,1].pcolormesh(X,Y,A,shading='gouraud',cmap='coolwarm',norm=LogNorm(vmin=.35,vmax=3),rasterized=True)
fig.colorbar(im,ax=axs[0,1],label=r'$D_\varphi/D_\varphi^{(0)}$')
axs[0,1].contour(X,Y,A,levels=[1],colors='k',linewidths=.9)
axs[0,1].set_title('(b) Physical phase diffusion',loc='left')
# RIN
A=np.where(st,D['RIN'],np.nan)
im=axs[1,0].pcolormesh(X,Y,A,shading='gouraud',cmap='coolwarm',norm=LogNorm(vmin=.25,vmax=4),rasterized=True)
fig.colorbar(im,ax=axs[1,0],label='integrated RIN ratio')
axs[1,0].contour(X,Y,A,levels=[1],colors='k',linewidths=.9)
axs[1,0].set_title('(c) Integrated relative-intensity noise',loc='left')
# classification; a small tolerance prevents the no-feedback baseline from
# being assigned to a gain/cost quadrant by roundoff.
K=np.broadcast_to(k[:,None],st.shape)
baseline=st&(K==0)
mfb=st&~baseline
cat=np.full(st.shape,4)  # unstable
cat[mfb&(D['Dphys']<1)&(D['RIN']<1)]=0
cat[mfb&(D['Dphys']<1)&(D['RIN']>=1)]=1
cat[mfb&(D['Dphys']>=1)&(D['RIN']<1)]=2
cat[mfb&(D['Dphys']>=1)&(D['RIN']>=1)]=3
cat[baseline]=5
cmap=ListedColormap(['#009E73','#E69F00','#56B4E9','#CC79A7','#B0B0B0','#FFFFFF'])
norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5,4.5,5.5],cmap.N)
im=axs[1,1].pcolormesh(X,Y,cat,shading='auto',cmap=cmap,norm=norm,rasterized=True)
cb=fig.colorbar(im,ax=axs[1,1],ticks=[0,1,2,3,4,5])
cb.ax.set_yticklabels(['win-win','linewidth gain / RIN cost','phase cost / RIN gain','bad-bad','unstable','no-feedback baseline'])
axs[1,1].set_title('(d) Noise-partition regimes',loc='left')
for ax in axs.flat:
    ax.set_xlabel(r'feedback phase $\phi/\pi$');ax.set_ylabel(r'$\kappa\tau$')
    ax.plot(0,.008*20,'o',ms=5,mec='k',mfc='white')
    ax.plot(.25,.009*20,'s',ms=5,mec='k',mfc='white')
fig.savefig(os.path.join(FIG,'fig2_partition_map.pdf'),bbox_inches='tight')
fig.savefig(os.path.join(FIG,'fig2_partition_map.png'),bbox_inches='tight',dpi=300)
plt.close(fig)

# FIG 3 residues
from lk_model import matrices,delta
from lk_analysis import pole_residue, characteristic_roots

def exact_H(om,p,e,B):
    A0,A1=matrices(p,e)
    DD=1j*om[:,None,None]*np.eye(3)[None]-A0[None]-np.exp(-1j*om*p.tau)[:,None,None]*A1[None]
    return np.linalg.solve(DD,np.broadcast_to(B,(len(om),)+B.shape))
def pair_H(om,p,e,roots,B):
    H=np.zeros((len(om),3,B.shape[1]),complex)
    for lam in roots:
        RB,_=pole_residue(lam,p,e,B)
        H+=RB[None]/(1j*om[:,None,None]-lam)
    return H
def out_psd(H,C):
    G=np.einsum('ab,wbc->wac',np.atleast_2d(C),H)
    return np.real(np.einsum('wac,wac->w',G,np.conj(G)))

fig,axs=plt.subplots(2,2,figsize=(7.05,5.2),constrained_layout=True)
B=np.array([[1.,0.],[0.,1.],[0.,0.]])
om=np.linspace(.02,.22,1800)
for col,(label,p,c) in enumerate(cases[1:]):
    e=solve_ecm(p)
    roots=characteristic_roots(p,e,M=180,nseeds=180,real_cut=-1,imag_cut=2)
    roots=np.array([z for z in roots if abs(z)>1e-7]);roots=roots[np.argsort(roots.real)[::-1]]
    pair=roots[:2]
    Hex=exact_H(om,p,e,B);H2=pair_H(om,p,e,pair,B)
    sI=out_psd(Hex,[2*e.A,0,0])/e.A**4;sI2=out_psd(H2,[2*e.A,0,0])/e.A**4
    axs[0,col].semilogy(om,sI,color=c,lw=1.6,label='full resolvent')
    axs[0,col].semilogy(om,sI2,color=colors['pair'],lw=1.2,ls='--',label='critical pair')
    axs[0,col].set_title(f'({chr(97+col)}) {label}: RIN spectrum',loc='left')
    axs[0,col].set_xlabel(r'$\omega$');axs[0,col].set_ylabel('RIN PSD')
    axs[0,col].legend(frameon=False)
    # roots
    axs[1,col].scatter(roots[:20].real,roots[:20].imag,s=16,facecolors='none',edgecolors=c)
    axs[1,col].scatter(pair.real,pair.imag,s=32,color=colors['pair'],marker='x')
    axs[1,col].axvline(0,color='k',lw=.7)
    axs[1,col].set_xlim(-.3,.03);axs[1,col].set_ylim(-1.3,1.3)
    axs[1,col].set_xlabel(r'$\mathrm{Re}\,\lambda$');axs[1,col].set_ylabel(r'$\mathrm{Im}\,\lambda$')
    axs[1,col].set_title(f'({chr(99+col)}) Delay spectrum and retained pair',loc='left')
fig.savefig(os.path.join(FIG,'fig3_residue_reduction.pdf'),bbox_inches='tight')
fig.savefig(os.path.join(FIG,'fig3_residue_reduction.png'),bbox_inches='tight',dpi=300)
plt.close(fig)

# FIG 4 robustness selected win-win points
r=pd.read_csv(os.path.join(BASE,'data','robustness_selected.csv'))
fig,axs=plt.subplots(1,2,figsize=(7.05,3.0),constrained_layout=True)
markers={10:'o',20:'s',40:'^'}
alpha_colors={2:'#56B4E9',4:'#E69F00',6:'#CC79A7'}
for tau,grp in r.groupby('tau'):
    for _,row in grp.iterrows():
        axs[0].scatter(row.Dphase_ratio,row.RIN_ratio,s=60,marker=markers[int(tau)],
                       color=alpha_colors[int(row.alpha)],edgecolor='k',linewidth=.3)
axs[0].axvline(1,color='k',lw=.7,ls='--');axs[0].axhline(1,color='k',lw=.7,ls='--')
axs[0].set_xlim(.5,.9);axs[0].set_ylim(.45,.82);axs[0].set_xlabel('phase-diffusion ratio');axs[0].set_ylabel('integrated RIN ratio')
axs[0].set_title('(a) Win-win points persist across parameters',loc='left')
# custom legends
from matplotlib.lines import Line2D
leg1=[Line2D([0],[0],marker=markers[t],color='none',markerfacecolor='white',markeredgecolor='k',label=f'$\\tau={t}$') for t in [10,20,40]]
leg2=[Line2D([0],[0],marker='o',color='none',markerfacecolor=alpha_colors[a],markeredgecolor='k',label=f'$\\alpha={a}$') for a in [2,4,6]]
L1=axs[0].legend(handles=leg1,frameon=False,loc='upper left');axs[0].add_artist(L1)
axs[0].legend(handles=leg2,frameon=False,loc='lower right')
for tau,grp in r.groupby('tau'):
    axs[1].plot(grp['alpha'],grp['phi']/np.pi,marker=markers[int(tau)],label=f'$\\tau={int(tau)}$')
axs[1].set_xlabel(r'Henry factor $\alpha$');axs[1].set_ylabel(r'selected feedback phase $\phi/\pi$')
axs[1].set_title('(b) Phase shifts with amplitude-phase coupling',loc='left')
axs[1].legend(frameon=False)
fig.savefig(os.path.join(FIG,'fig4_robustness.pdf'),bbox_inches='tight')
fig.savefig(os.path.join(FIG,'fig4_robustness.png'),bbox_inches='tight',dpi=300)
plt.close(fig)
print('Figures written to',FIG)
