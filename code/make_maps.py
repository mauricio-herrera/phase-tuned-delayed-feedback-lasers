import numpy as np, os, time
from scipy.linalg import eigvals
from lk_model import Params, solve_ecm, matrices, output_variance_grid, phase_diffusion, default_omega_grid
from lk_analysis import chain_matrix

from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
out=str(BASE/'data')
os.makedirs(out,exist_ok=True)
alpha=4.;eps=.005;P=1.;tau=20.
kappas=np.linspace(0,0.009,46)
phis=np.linspace(-np.pi,np.pi,145)
omega=default_omega_grid()
shape=(len(kappas),len(phis))
arrs={name:np.full(shape,np.nan) for name in ['theta','A','Ns','omega_ecm','margin','Dgeom','Dphys','varI','RIN','varN','peakRIN','peakfreq','ROarea']}
# core refs
p0=Params(alpha=alpha,eps=eps,P=P,kappa=0,tau=tau,phi=0)
e0=solve_ecm(p0)
D0,_=phase_diffusion(p0,e0)
# physical phase with vphys gives D0 same 8.5. Our phase_diffusion currently uses vphys.
varI0,om,psdI0=output_variance_grid(p0,e0,[2*e0.A,0,0],omega=omega,return_psd=True)
varN0=output_variance_grid(p0,e0,[0,0,1],omega=omega)
RIN0=varI0/e0.A**4
# define relaxation band around core peak excluding low freq
mask=(om>0.02)&(om<0.3)
ipeak=np.argmax(psdI0[mask]); wpk0=om[mask][ipeak]
band=(om>max(0.005,wpk0-0.04))&(om<wpk0+0.04)
area0=np.trapezoid(psdI0[band]/e0.A**4,om[band])/np.pi
print('core',D0,varI0,RIN0,varN0,wpk0,area0)
t0=time.time()
for ik,k in enumerate(kappas):
  for jp,phi in enumerate(phis):
    p=Params(alpha=alpha,eps=eps,P=P,kappa=float(k),tau=tau,phi=float(phi))
    try:
      e=solve_ecm(p)
      arrs['theta'][ik,jp]=e.theta;arrs['A'][ik,jp]=e.A;arrs['Ns'][ik,jp]=e.Ns;arrs['omega_ecm'][ik,jp]=e.omega
      # stability via chain M=45
      vals=eigvals(chain_matrix(p,e,M=45))
      vals=vals[np.abs(vals)>1e-6]
      lead=vals[np.argmax(vals.real)]
      arrs['margin'][ik,jp]=-lead.real
      # geometric chirp factor
      F=1+k*tau*(np.cos(e.theta)-alpha*np.sin(e.theta))
      arrs['Dgeom'][ik,jp]=1/(F*F)
      D,_=phase_diffusion(p,e)
      arrs['Dphys'][ik,jp]=D/D0
      vi,om,psd=output_variance_grid(p,e,[2*e.A,0,0],omega=omega,return_psd=True)
      vn=output_variance_grid(p,e,[0,0,1],omega=omega)
      arrs['varI'][ik,jp]=vi/varI0
      arrs['RIN'][ik,jp]=(vi/e.A**4)/RIN0
      arrs['varN'][ik,jp]=vn/varN0
      mm=(om>0.02)&(om<0.3)
      idx=np.argmax(psd[mm]/e.A**4)
      wpk=om[mm][idx]
      arrs['peakfreq'][ik,jp]=wpk
      arrs['peakRIN'][ik,jp]=(psd[mm][idx]/e.A**4)/(psdI0[mask][ipeak]/e0.A**4)
      b=(om>max(0.005,wpk-0.04))&(om<wpk+0.04)
      area=np.trapezoid(psd[b]/e.A**4,om[b])/np.pi
      arrs['ROarea'][ik,jp]=area/area0
    except Exception as ex:
      pass
  print(ik,k,'elapsed',time.time()-t0)
np.savez(os.path.join(out,'map_tau20.npz'),kappas=kappas,phis=phis,omega=omega,core=np.array([D0,varI0,RIN0,varN0,wpk0,area0]),**arrs)
# summary stable points
stable=arrs['margin']>0
print('stable fraction',np.nanmean(stable))
for name in ['Dgeom','Dphys','RIN','varN','peakRIN','ROarea']:
 a=arrs[name]
 print(name,'stable min/max',np.nanmin(a[stable]),np.nanmax(a[stable]))
# classify stable based Dphys and RIN
print('winwin',np.sum(stable&(arrs['Dphys']<1)&(arrs['RIN']<1)))
print('trade',np.sum(stable&(arrs['Dphys']<1)&(arrs['RIN']>1)))
print('badbad',np.sum(stable&(arrs['Dphys']>1)&(arrs['RIN']>1)))
print('phasepen',np.sum(stable&(arrs['Dphys']>1)&(arrs['RIN']<1)))
# top candidates
score=(1-arrs['Dphys'])+(1-arrs['RIN'])
score[~stable]=np.nan
for idx in np.dstack(np.unravel_index(np.argsort(np.nan_to_num(score,nan=-1e9).ravel())[::-1][:10],shape))[0]:
 ik,jp=idx
 print('best',kappas[ik],phis[jp],arrs['Dphys'][ik,jp],arrs['RIN'][ik,jp],arrs['varN'][ik,jp],arrs['margin'][ik,jp])
