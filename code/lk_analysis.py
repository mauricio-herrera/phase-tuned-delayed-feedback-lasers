import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import root
from lk_model import Params, solve_ecm, matrices, delta, null_vectors, delta_prime, phase_diffusion, output_variance_quad


def chain_matrix(p,ecm,M=80):
    A0,A1=matrices(p,ecm)
    n=3+2*M
    A=np.zeros((n,n),float)
    A[:3,:3]=A0
    A[:3,3+2*(M-1):3+2*M]=A1[:,:2]
    r=M/p.tau
    # first delay cell driven by current field
    A[3:5,0:2]=r*np.eye(2)
    A[3:5,3:5]=-r*np.eye(2)
    for j in range(1,M):
        rows=slice(3+2*j,3+2*(j+1))
        prev=slice(3+2*(j-1),3+2*j)
        A[rows,prev]=r*np.eye(2)
        A[rows,rows]=-r*np.eye(2)
    return A


def refine_root(seed,p,ecm):
    A0,A1=matrices(p,ecm)
    def fun(z):
        lam=z[0]+1j*z[1]
        d=np.linalg.det(delta(lam,A0,A1,p.tau))
        return [d.real,d.imag]
    sol=root(fun,[seed.real,seed.imag],tol=1e-11)
    lam=sol.x[0]+1j*sol.x[1]
    res=np.linalg.norm(fun(sol.x))
    return lam,sol.success,res


def characteristic_roots(p,ecm,M=80,nseeds=40,real_cut=-1.0,imag_cut=10.0):
    vals=eigvals(chain_matrix(p,ecm,M=M))
    vals=vals[(vals.real>real_cut)&(np.abs(vals.imag)<imag_cut)]
    vals=vals[np.argsort(vals.real)[::-1]]
    roots=[]
    for seed in vals[:nseeds]:
        lam,ok,res=refine_root(seed,p,ecm)
        if not ok or res>1e-7:
            continue
        if lam.real<real_cut-0.1 or abs(lam.imag)>imag_cut+1:
            continue
        if all(abs(lam-r)>1e-6 for r in roots):
            roots.append(lam)
    roots=sorted(roots,key=lambda z:z.real,reverse=True)
    return np.array(roots)


def leading_nonzero_root(p,ecm,M=80):
    rr=characteristic_roots(p,ecm,M=M,nseeds=60,real_cut=-0.5,imag_cut=2.0)
    rr=[z for z in rr if abs(z)>1e-7]
    if not rr:
        return np.nan+1j*np.nan, np.array([])
    return rr[0],np.array(rr)


def pole_residue(lam,p,ecm,B=None):
    A0,A1=matrices(p,ecm)
    if B is None:
        B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    u,v,res=null_vectors(lam,A0,A1,p.tau)
    d=np.vdot(u,delta_prime(lam,A1,p.tau)@v)
    R=np.outer(v,np.conj(u))/d
    return R@B, dict(u=u,v=v,d=d,res=res)


def multipole_transfer(omega,roots,p,ecm,B=None):
    if B is None:
        B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    H=np.zeros((3,B.shape[1]),complex)
    for lam in roots:
        RB,_=pole_residue(lam,p,ecm,B)
        H += RB/(1j*omega-lam)
    return H

if __name__=='__main__':
    p=Params(kappa=.008,tau=20,phi=0)
    e=solve_ecm(p)
    rr=characteristic_roots(p,e,M=120,nseeds=80,real_cut=-1,imag_cut=2)
    print(rr[:20])
