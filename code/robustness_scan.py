import numpy as np, os
from scipy.linalg import eigvals
from lk_model import Params,solve_ecm,phase_diffusion,output_variance_grid,default_omega_grid
from lk_analysis import chain_matrix

from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
out=str(BASE/'data');os.makedirs(out,exist_ok=True)
omega=default_omega_grid()
configs=[]
for alpha in [2.,4.,6.]:
 for tau in [10.,20.,40.]:
  # fixed max g=kappa*tau=.18
  kappas=np.linspace(0,.18/tau,13)
  phis=np.linspace(-np.pi,np.pi,25)
  p0=Params(alpha=alpha,eps=.005,P=1,kappa=0,tau=tau,phi=0);e0=solve_ecm(p0)
  D0,_=phase_diffusion(p0,e0)
  vi0=output_variance_grid(p0,e0,[2,0,0],omega=omega);rin0=vi0
  vn0=output_variance_grid(p0,e0,[0,0,1],omega=omega)
  count={'stable':0,'winwin':0,'trade':0,'badbad':0,'phasepen':0}
  best=None
  for k in kappas:
   for phi in phis:
    p=Params(alpha=alpha,eps=.005,P=1,kappa=float(k),tau=tau,phi=float(phi))
    try:
     e=solve_ecm(p)
     vals=eigvals(chain_matrix(p,e,M=35));vals=vals[np.abs(vals)>1e-6]
     margin=-np.max(vals.real)
     if margin<=0:continue
     D,_=phase_diffusion(p,e);dr=D/D0
     vi=output_variance_grid(p,e,[2*e.A,0,0],omega=omega);rr=(vi/e.A**4)/rin0
     vn=output_variance_grid(p,e,[0,0,1],omega=omega)/vn0
     count['stable']+=1
     if dr<1 and rr<1:cat='winwin'
     elif dr<1 and rr>=1:cat='trade'
     elif dr>=1 and rr>=1:cat='badbad'
     else:cat='phasepen'
     count[cat]+=1
     score=(1-dr)+(1-rr)
     if best is None or score>best[0]:best=(score,k,phi,dr,rr,vn,margin)
    except Exception:pass
  configs.append((alpha,tau,count,best))
  print(alpha,tau,count,best)
np.save(str(BASE/'data'/'robustness.npy'),np.array(configs,dtype=object),allow_pickle=True)

# Save the selected representative point for each (alpha,tau) pair.
import pandas as pd
rows=[]
for alpha,tau,count,best in configs:
    if best is None:
        continue
    score,k,phi,dr,rr,vn,margin=best
    rows.append(dict(alpha=alpha,tau=tau,kappa=k,phi=phi,
                     Dphase_ratio=dr,RIN_ratio=rr,carrier_ratio=vn,margin=margin))
pd.DataFrame(rows).to_csv(str(BASE/'data'/'robustness_selected.csv'),index=False)
