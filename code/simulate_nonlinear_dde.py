import numpy as np
from pathlib import Path
from dataclasses import asdict
from lk_model import Params, solve_ecm, phase_diffusion, output_variance_grid, default_omega_grid

BASE=Path(__file__).resolve().parents[1]

def simulate_case(p, sigma=1e-3, carrier_sigma=0.0, n_ens=32, dt=0.02, T=2000.0, burn=400.0, record_dt=0.1, seed=0):
    e=solve_ecm(p)
    ndelay=int(round(p.tau/dt))
    if abs(ndelay*dt-p.tau)>1e-10: raise ValueError('tau/dt must be integer')
    nsteps=int(round(T/dt)); stride=int(round(record_dt/dt))
    rng=np.random.default_rng(seed)
    # ring history in co-rotating frame, steady ECM
    hist=np.full((ndelay+1,n_ens),e.A+0j,dtype=np.complex128)
    z=np.full(n_ens,e.A+0j,dtype=np.complex128)
    N=np.full(n_ens,e.Ns,dtype=float)
    theta=e.theta; fb=p.kappa*np.exp(-1j*theta)
    times=[]; zs=[]; Ns=[]
    ptr=ndelay
    for step in range(1,nsteps+1):
        zdel=hist[(ptr+1)%(ndelay+1)]
        drift=((1+1j*p.alpha)*N-1j*e.omega)*z + fb*zdel
        dW=np.sqrt(dt)*(rng.standard_normal(n_ens)+1j*rng.standard_normal(n_ens))
        znew=z+dt*drift+sigma*dW
        nd=p.eps*(p.P-N-(1+2*N)*np.abs(z)**2)
        Nnew=N+dt*nd
        if carrier_sigma:
            Nnew += carrier_sigma*np.sqrt(dt)*rng.standard_normal(n_ens)
        ptr=(ptr+1)%(ndelay+1)
        hist[ptr]=znew
        z,N=znew,Nnew
        if step%stride==0:
            times.append(step*dt); zs.append(z.copy()); Ns.append(N.copy())
    t=np.asarray(times); Z=np.asarray(zs); NN=np.asarray(Ns)
    q=np.unwrap(np.angle(Z),axis=0)
    varq=np.var(q,axis=1,ddof=1)
    # Diffusion estimate from time- and ensemble-averaged phase increments.
    # This estimator uses the stationary increment process and is substantially
    # less noisy than fitting the ensemble variance from a single time origin.
    stat=t>=burn
    qs=q[stat]
    dtrec=float(record_dt)
    max_lag_time=min(1500.0,0.30*(T-burn))
    min_lag_time=max(150.0,10*dtrec)
    lag_times=np.geomspace(min_lag_time,max_lag_time,18)
    lag_idx=np.unique(np.maximum(1,np.rint(lag_times/dtrec).astype(int)))
    lag_idx=lag_idx[lag_idx < qs.shape[0]//2]
    lag_times=lag_idx*dtrec
    msd=np.array([np.mean((qs[L:]-qs[:-L])**2) for L in lag_idx])
    slope,inter=np.polyfit(lag_times,msd,1)
    Dhat=slope/2
    # realization-level increment fits provide an ensemble standard error
    Dj=[]
    for j in range(n_ens):
        mj=np.array([np.mean((qs[L:,j]-qs[:-L,j])**2) for L in lag_idx])
        Dj.append(np.polyfit(lag_times,mj,1)[0]/2)
    Dj=np.asarray(Dj)
    Dse=float(np.std(Dj,ddof=1)/np.sqrt(n_ens))
    I=np.abs(Z[stat])**2
    rin=(I-I.mean(axis=0,keepdims=True))/(e.A**2)
    # realization-level time variances for uncertainty
    rin_var_j=np.var(rin,axis=0,ddof=1)
    ndev=NN[stat]-NN[stat].mean(axis=0,keepdims=True)
    n_var_j=np.var(ndev,axis=0,ddof=1)
    return dict(params=asdict(p),ecm=asdict(e),sigma=sigma,carrier_sigma=carrier_sigma,n_ens=n_ens,dt=dt,T=T,burn=burn,
                Dhat=Dhat,Dse=Dse,Dci=np.array([Dhat-1.96*Dse,Dhat+1.96*Dse]),
                lag_times=lag_times,phase_msd=msd,D_per_realization=Dj,
                RINhat=float(np.mean(rin_var_j)),RINse=float(np.std(rin_var_j,ddof=1)/np.sqrt(n_ens)),
                Nhat=float(np.mean(n_var_j)),Nse=float(np.std(n_var_j,ddof=1)/np.sqrt(n_ens)),
                t=t,varq=varq,rin_var_j=rin_var_j,n_var_j=n_var_j,
                q_sample=q[:,0],rin_sample=rin[:,0],N_sample=ndev[:,0])

def run(dt=0.02,n_ens=32,T=2000.0,tag='main',sigma=1e-3):
    cases=[('solitary',Params(kappa=0,tau=20,phi=0)),
           ('tradeoff',Params(kappa=.008,tau=20,phi=0)),
           ('win_win',Params(kappa=.009,tau=20,phi=np.pi/4))]
    rows=[]
    omega=default_omega_grid()
    for i,(name,p) in enumerate(cases):
        print('running',name,'dt',dt)
        r=simulate_case(p,sigma=sigma,n_ens=n_ens,dt=dt,T=T,seed=1729+100*i)
        e=solve_ecm(p)
        D,_=phase_diffusion(p,e)
        R=output_variance_grid(p,e,[2/e.A,0,0],omega=omega)
        Nv=output_variance_grid(p,e,[0,0,1],omega=omega)
        r['Dlin']=D*r['sigma']**2; r['RINlin']=R*r['sigma']**2; r['Nlin']=Nv*r['sigma']**2
        np.savez_compressed(BASE/'data'/f'sde_{tag}_{name}.npz',**{k:v for k,v in r.items() if k not in ['params','ecm']})
        rows.append((name,r))
    # summary relative ratios
    out=[]
    base=rows[0][1]
    for name,r in rows:
        out.append([name,r['Dhat'],r['Dse'],r['Dlin'],r['Dhat']/base['Dhat'],r['Dlin']/base['Dlin'],
                    r['RINhat'],r['RINse'],r['RINlin'],r['RINhat']/base['RINhat'],r['RINlin']/base['RINlin'],
                    r['Nhat'],r['Nse'],r['Nlin'],r['Nhat']/base['Nhat'],r['Nlin']/base['Nlin']])
    import csv
    path=BASE/'data'/f'sde_{tag}_summary.csv'
    with open(path,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['case','Dhat','Dse','Dlin','D_ratio_sim','D_ratio_lin','RINhat','RINse','RINlin','RIN_ratio_sim','RIN_ratio_lin','Nhat','Nse','Nlin','N_ratio_sim','N_ratio_lin']);w.writerows(out)
    print(open(path).read())

if __name__=='__main__':
    run()
