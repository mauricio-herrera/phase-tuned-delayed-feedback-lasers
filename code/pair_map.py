import numpy as np
from lk_model import Params,solve_ecm
from pair_formula import modal_variance
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
f=str(BASE/'data'/'map_tau20_exact.npz')
d=np.load(f)
kappas=d['kappas']; phis=d['phis']; marg=d['margin_exact']
pfi=np.full_like(marg,np.nan); pfn=np.full_like(marg,np.nan)
for i,k in enumerate(kappas):
 for j,phi in enumerate(phis):
  if not np.isfinite(marg[i,j]) or marg[i,j]<=0: continue
  im=abs(d['leadim'][i,j]); re=d['leadre'][i,j]
  if im<1e-8: continue
  p=Params(kappa=float(k),tau=20,phi=float(phi));e=solve_ecm(p)
  roots=np.array([re+1j*im,re-1j*im])
  vi_pair=modal_variance(roots,p,e,[2*e.A,0,0])
  vn_pair=modal_variance(roots,p,e,[0,0,1])
  # exact ratios times core values
  vi_exact=d['varI'][i,j]*136.32697994316268
  vn_exact=d['varN'][i,j]*0.33333330315847454
  pfi[i,j]=vi_pair/vi_exact
  pfn[i,j]=vn_pair/vn_exact
out={k:d[k] for k in d.files};out['pairfracI']=pfi;out['pairfracN']=pfn
np.savez(str(BASE/'data'/'map_tau20_full.npz'),**out)
st=marg>0
print('I frac percentiles',np.nanpercentile(pfi[st],[0,5,25,50,75,95,100]))
print('N frac percentiles',np.nanpercentile(pfn[st],[0,5,25,50,75,95,100]))
