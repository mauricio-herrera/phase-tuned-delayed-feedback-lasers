"""Independent, documented numerical routines for revision JOSA B 612251.

This module evaluates the analytic LK delay resolvent; no values are read from
manuscript tables. Convention: theta = phi + Omega*tau,
phi(theta) = theta + kappa*tau*(alpha*cos(theta)+sin(theta)).
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, root
from lk.rev_core import Model, ecm_theta, mats, variance, wrap

HERE=Path(__file__).resolve().parent
OUT=HERE/'revision_outputs'
OUT.mkdir(exist_ok=True)

def export_csv(name,rows):
    if not rows: raise ValueError('No rows')
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('Wrote', OUT/name, 'rows',len(rows))

def phase_edges(model,kappa):
    """Closest RIN-left and linewidth-right crossings around theta=0.

    Nonmonotonic ratios outside the joint window can have additional roots;
    using one global bracket can therefore miss the physical central window.
    """
    tau=model.tau
    def nearest_edge(which, left, right):
        grid=np.linspace(left,right,49)
        vals=[]
        for t in grid:
            st=model.state(float(t),kappa,stability=False)
            vals.append(st[which]-1)
        ix=np.flatnonzero(np.asarray(vals[:-1])*np.asarray(vals[1:])<0)
        if not len(ix):
            raise RuntimeError(f'No edge {which} tau={tau}, k={kappa}; endpoint signs {vals[0]},{vals[-1]}')
        # On negative side, select last transition; on positive, first.
        j=ix[-1] if right<=0 else ix[0]
        return brentq(lambda t:model.state(float(t),kappa,stability=False)[which]-1,grid[j],grid[j+1],xtol=1e-10)
    lo=nearest_edge('Rratio',-.30*np.pi,0.)
    hi=nearest_edge('Dratio',0.,.30*np.pi)
    below,above=model.state(lo,kappa,stability=True),model.state(hi,kappa,stability=True)
    mid=model.state((lo+hi)*.5,kappa,stability=True)
    if any(x['nunst']>0 or x['F']<=0 for x in (below,above,mid)):
        raise RuntimeError(f'Window edges not stable for tau={tau},k={kappa}')
    return lo,hi

def phase_map(theta,kt,alpha=4.):
    return theta+kt*(alpha*np.cos(theta)+np.sin(theta))

def all_theta(phi,kt,alpha=4.,n=5001):
    """Enumerate ECMs at phi modulo 2pi. Each root in [-pi,pi)."""
    grid=np.linspace(-np.pi,np.pi,n)
    q=grid+kt*(alpha*np.cos(grid)+np.sin(grid))
    low=int(np.floor((min(q)-phi)/(2*np.pi)))-1
    high=int(np.ceil((max(q)-phi)/(2*np.pi)))+1
    result=[]
    for m in range(low,high+1):
        a=q-(phi+2*np.pi*m)
        indices=np.nonzero(a[:-1]*a[1:]<=0)[0]
        for j in indices:
            z=brentq(lambda t:phase_map(t,kt,alpha)-(phi+2*np.pi*m),grid[j],grid[j+1],xtol=1e-12)
            if not any(abs(wrap(z-v))<1e-9 for v in result):result.append(z)
    return np.array(sorted(result))

def ecm_continuation(phi,kappa,tau,alpha=4.):
    """Return principal theta by integrating dtheta/d(kt) at fixed phi.

    Differential continuation avoids Newton basin jumps near multistable folds.
    Returns NaN if the branch terminates at F=0 before target kappa.
    """
    from scipy.integrate import solve_ivp
    def F(t,y):return 1+t*(np.cos(y[0])-alpha*np.sin(y[0]))
    def dy(t,y):return [-(alpha*np.cos(y[0])+np.sin(y[0]))/F(t,y)]
    def fold(t,y):return F(t,y)
    fold.terminal=True;fold.direction=-1
    sol=solve_ivp(dy,(0.,kappa*tau),[float(wrap(phi))],events=fold,
                  rtol=1e-10,atol=1e-12,max_step=.004)
    return float(wrap(sol.y[0,-1])) if sol.status==0 else float('nan')

def characteristic(lam,A0,A1,tau):
    return np.linalg.det(lam*np.eye(3)-A0-A1*np.exp(-lam*tau))

def relaxation_residue(tau,kappa,theta,alpha=4.,eps=.005,P=1.,guess=None):
    e=ecm_theta(theta,kappa,tau,alpha,P)
    A0,A1=mats(e,kappa,alpha,eps)
    if guess is None:
        ee=ecm_theta(0.,0.,tau,alpha,P)
        M0,_=mats(ee,0.,alpha,eps)
        ev=np.linalg.eigvals(M0)
        z=ev[np.argmax(ev.imag)]
    else:z=guess
    def f(t):
        zz=t[0]+1j*t[1]
        v=characteristic(zz,A0,A1,tau)
        return (v.real,v.imag)
    opt=root(f,[z.real,z.imag],tol=2e-12)
    lam=opt.x[0]+1j*opt.x[1]
    if abs(characteristic(lam,A0,A1,tau))>1e-8 or lam.imag<1e-5:
        raise RuntimeError(f'Relaxation pair refinement failed {tau} {kappa} {theta} {lam}')
    U,s,Vh=np.linalg.svd(lam*np.eye(3)-A0-A1*np.exp(-lam*tau))
    v=Vh.conj().T[:,-1];u=U[:,-1]
    DeltaP=np.eye(3)+tau*A1*np.exp(-lam*tau)
    d=np.vdot(u,DeltaP@v)
    B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    C=np.array([2/e['A'],0.,0.])
    q=(C@v)*(u.conj()@B)/d
    Q=np.conj(q)
    diag=float(np.sum(abs(q)**2)/abs(lam.real))
    pair=0j
    for l1,q1 in [(lam,q),(lam.conjugate(),Q)]:
        for l2,q2 in [(lam,q),(lam.conjugate(),Q)]:
            pair+=-np.vdot(q2,q1)/(l1+l2.conjugate())
    return dict(lam=lam,d=d,u=u,v=v,C=C,B=B,q=q,diagonal=diag,pair=float(pair.real))

def fo_rin_coefficients(tau,alpha=4.,eps=.005,P=1.):
    """Exact first derivative of the integrated power spectrum at kappa=0.

    Differentiates the resolvent, output row and ECM amplitude analytically.
    Uses the same finite numerical frequency quadrature as the full theory.
    """
    m=Model(tau,alpha,eps,P)
    om=m.om;A=np.sqrt(P);O=om
    e0=ecm_theta(0,0,tau,alpha,P)
    A0,A1=mats(e0,0,alpha,eps)
    D=1j*O[:,None,None]*np.eye(3)-A0
    H=np.linalg.solve(D,np.broadcast_to(np.array([[1.,0.],[0.,1.],[0.,0.]]),(len(O),3,2)))
    C=np.array([2/A,0.,0.]);G=H[:,:, :]
    G0=np.einsum('a,wab->wb',C,H)
    def deriv(phi):
        c=np.cos(phi);s=np.sin(phi)
        Ap=(1+2*P)*c/(2*A)
        dA0=np.array([[-c,-s,Ap],[s,-c,alpha*Ap],[-2*eps*(Ap-2*A*c),0.,-2*eps*(1+2*P)*c]])
        dA1=np.array([[c,s,0.],[-s,c,0.],[0.,0.,0.]])
        kk=dA0[None,:,:]+np.exp(-1j*O*tau)[:,None,None]*dA1
        dH=np.linalg.solve(D,np.einsum('wab,wbc->wac',kk,H))
        dC=np.array([-2*Ap/A**2,0,0])
        dg=np.einsum('a,wab->wb',dC,H)+np.einsum('a,wab->wb',C,dH)
        integrand=2*np.real(np.sum(dg*G0.conj(),axis=1))
        return np.trapezoid(integrand,O)/np.pi/m.R0
    return deriv(0.),deriv(np.pi/2)
