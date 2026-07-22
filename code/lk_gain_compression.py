import numpy as np
from dataclasses import dataclass
from scipy.optimize import root, root_scalar
from lk_model import Params, ECM

@dataclass
class GCParams(Params):
    gain_compression: float = 0.0

def solve_ecm_gc(p:GCParams,theta_guess=None,branch_m=0):
    base=p.phi+2*np.pi*branch_m
    def ft(th):return th-base+p.kappa*p.tau*(p.alpha*np.cos(th)+np.sin(th))
    if theta_guess is None:theta_guess=base
    sol=root(lambda z:[ft(z[0])],[theta_guess],tol=1e-12)
    if not sol.success or abs(ft(sol.x[0]))>1e-9:raise RuntimeError('theta root')
    th=float(sol.x[0]); c=np.cos(th); eta=p.gain_compression
    if eta==0:
        Ns=-p.kappa*c
        I=(p.P-Ns)/(1+2*Ns)
    else:
        def fi(I):
            N=-p.kappa*c*(1+eta*I)
            return p.P-N-(1+2*N)*I/(1+eta*I)
        xs=np.geomspace(1e-8,max(100,10*(p.P+1)),400)
        vals=[fi(x) for x in xs]
        idx=[i for i in range(len(xs)-1) if vals[i]*vals[i+1]<=0]
        if not idx: raise RuntimeError('intensity root')
        I=root_scalar(fi,bracket=[xs[idx[0]],xs[idx[0]+1]],xtol=1e-13).root
        Ns=-p.kappa*c*(1+eta*I)
    if I<=0:raise RuntimeError('nonphysical')
    A=np.sqrt(I);omega=-p.kappa*(p.alpha*c+np.sin(th))
    return ECM(theta=th,omega=omega,Ns=Ns,A=A)

def matrices_gc(p:GCParams,e:ECM):
    eta=p.gain_compression;I=e.A**2;den=1+eta*I
    G=e.Ns/den; GN=1/den; GI=-eta*e.Ns/den**2
    cr=p.kappa*np.cos(e.theta); s=p.kappa*np.sin(e.theta)
    A0=np.zeros((3,3),float)
    hx=2*e.A**2*GI
    A0[0,0]=-cr+hx; A0[0,1]=-s; A0[0,2]=e.A*GN
    A0[1,0]=s+p.alpha*hx; A0[1,1]=-cr; A0[1,2]=p.alpha*e.A*GN
    A0[2,0]=-2*p.eps*e.A*(1+2*e.Ns)/den**2
    A0[2,2]=p.eps*(-1-2*I/den)
    A1=np.zeros((3,3),float)
    A1[0,0]=cr;A1[0,1]=s;A1[1,0]=-s;A1[1,1]=cr
    return A0,A1
