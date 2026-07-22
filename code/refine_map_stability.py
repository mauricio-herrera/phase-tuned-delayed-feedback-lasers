import numpy as np, os, time
from scipy.linalg import eigvals
from lk_model import Params,solve_ecm
from lk_analysis import chain_matrix,refine_root
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
f=str(BASE/'data'/'map_tau20.npz')
d=np.load(f)
kappas=d['kappas'];phis=d['phis']
margin=np.full((len(kappas),len(phis)),np.nan,dtype=float)
leadre=np.full_like(margin,np.nan);leadim=np.full_like(margin,np.nan)
t0=time.time()
for ik,k in enumerate(kappas):
 for jp,phi in enumerate(phis):
  try:
   p=Params(kappa=float(k),tau=20,phi=float(phi));e=solve_ecm(p)
   vals=eigvals(chain_matrix(p,e,M=55));vals=vals[np.abs(vals)>1e-6]
   vals=vals[np.argsort(vals.real)[::-1]][:8]
   roots=[]
   for seed in vals:
    lam,ok,res=refine_root(seed,p,e)
    if ok and res<1e-7 and abs(lam)>1e-7 and all(abs(lam-r)>1e-6 for r in roots):roots.append(lam)
   if roots:
    lead=max(roots,key=lambda z:z.real)
    margin[ik,jp]=-lead.real;leadre[ik,jp]=lead.real;leadim[ik,jp]=lead.imag
  except Exception:pass
 print(ik,time.time()-t0)
# rewrite new npz preserving data
out={k:d[k] for k in d.files}
out['margin_exact']=margin;out['leadre']=leadre;out['leadim']=leadim
np.savez(str(BASE/'data'/'map_tau20_exact.npz'),**out)
print('diff',np.nanmax(np.abs(margin-d['margin'])))
