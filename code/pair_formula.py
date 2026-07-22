import numpy as np
from lk_analysis import pole_residue

def modal_variance(roots,p,e,C,B=None):
    if B is None:
        B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    roots=np.asarray(roots,dtype=complex)
    C=np.atleast_2d(np.asarray(C,dtype=complex))
    Ms=[]
    for lam in roots:
        RB,_=pole_residue(lam,p,e,B)
        Ms.append(C@RB)
    V=0j
    for j,lj in enumerate(roots):
        for k,lk in enumerate(roots):
            V += -np.trace(Ms[j]@Ms[k].conj().T)/(lj+np.conj(lk))
    return float(V.real)
