#!/usr/bin/env python3
"""Nine alpha-delay ECM survey, reproducing supplement S7 Table S8.

Default full: kappa=0.001,...,0.009 and 240 total-phase samples,
consistent with the supplement. --quick uses 60 phases for a smoke test.
No table numbers are used as inputs.
"""
import argparse
import numpy as np
from scipy.integrate import solve_ivp
from lk.rev_core import Model,wrap,mats,count_right
from rev_analysis import export_csv,all_theta,phase_map

p=argparse.ArgumentParser();p.add_argument('--quick',action='store_true');p.add_argument('--single',nargs=2,type=float,metavar=('ALPHA','TAU'));a=p.parse_args()
rows=[];details=[]
cases=[(float(al),float(tau)) for al in (2,4,6) for tau in (10,20,40)]
if a.single:cases=[tuple(a.single)]

def isprincipal(phi,theta,k,alpha,tau):
    """ODE continuation at constant laboratory phase, independent of Newton basin.
    dtheta/d(kt)=-(alpha*cos(theta)+sin(theta))/F.
    Terminates at a saddle-node when F=0.
    """
    def denom(t,y):return 1+t*(np.cos(y[0])-alpha*np.sin(y[0]))
    def ode(t,y):return [-(alpha*np.cos(y[0])+np.sin(y[0]))/denom(t,y)]
    def fold(t,y):return denom(t,y)
    fold.terminal=True;fold.direction=-1
    sol=solve_ivp(ode,(0,k*tau),[phi],events=fold,rtol=1e-10,atol=1e-12,max_step=.004)
    return bool(sol.status==0 and abs(wrap(sol.y[0,-1]-theta))<1e-6)

for alpha,tau in cases:
    model=Model(tau,alpha=alpha)
    nphase=60 if a.quick else 240
    ths=np.linspace(-np.pi,np.pi,nphase,endpoint=False)
    total=unique=principal=multi=multistable=0
    best=None
    for k in np.arange(1,10)*.001:
        kt=k*tau;C=kt*np.sqrt(1+alpha*alpha)
        for th in ths:
            e=model.state(th,k,stability=False)
            # Filter before expensive RIN quadrature: physical ECM and narrower linewidth.
            if e['F']<=0 or e['Dratio']>=1:continue
            if e['Rratio']>=1:continue
            A0,A1=mats(e,k,alpha)
            if count_right(A0,A1,tau,0.)!=0:continue
            total+=1
            if C<1:unique+=1
            phi=float(wrap(e['phi']))
            if C<1:
                isp=True; nsols=1; nstable=1
            else:
                sols=all_theta(phi,kt,alpha,n=1501)
                nsols=len(sols)
                isp=isprincipal(phi,th,k,alpha,tau)
                nstable=1
                for other in sols:
                    if abs(wrap(other-th))<1e-6:continue
                    eo=model.state(other,k,stability=False)
                    if eo['F']<=0:continue
                    AA,BB=mats(eo,k,alpha)
                    if count_right(AA,BB,tau,0.)==0:nstable+=1
            principal+=int(isp)
            multi+=int(nsols>1)
            multistable+=int(nstable>1)
            objective=2-e['Dratio']-e['Rratio']
            if C<1 and (best is None or objective>best['score']):
                best=dict(score=objective,kappa=k,phi_pi=phi/np.pi,linewidth=e['Dratio'],rin=e['Rratio'])
            details.append(dict(alpha=alpha,tau=tau,kappa=k,theta_pi=th/np.pi,phi_pi=phi/np.pi,C=C,
                        D=e['Dratio'],R=e['Rratio'],principal=int(isp),multi_ecm=int(nsols>1),multistable=int(nstable>1)))
    row=dict(alpha=alpha,tau=tau,Cmax=.009*tau*np.sqrt(1+alpha**2),total=total,C_lt_1=unique,principal=principal,
             multi_ecm=multi,multistable=multistable,**(best or {}))
    rows.append(row);print('ECM_AUDIT',row,flush=True)
suffix=(f'_a{int(a.single[0])}_t{int(a.single[1])}' if a.single else '') + ('_quick' if a.quick else '')
export_csv('ecm_audit_summary'+suffix+'.csv',rows)
export_csv('ecm_audit_states'+suffix+'.csv',details)
