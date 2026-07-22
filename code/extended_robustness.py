import numpy as np,csv
from pathlib import Path
from scipy.linalg import eigvals
from lk_gain_compression import GCParams,solve_ecm_gc,matrices_gc
BASE=Path(__file__).resolve().parents[1]

def chain(A0,A1,tau,M=30):
    n=3+2*M;A=np.zeros((n,n));A[:3,:3]=A0;A[:3,3+2*(M-1):3+2*M]=A1[:,:2];r=M/tau
    A[3:5,:2]=r*np.eye(2);A[3:5,3:5]=-r*np.eye(2)
    for j in range(1,M):
        rr=slice(3+2*j,3+2*j+2);pp=slice(3+2*(j-1),3+2*j);A[rr,pp]=r*np.eye(2);A[rr,rr]=-r*np.eye(2)
    return A

def delta(lam,A0,A1,tau):return lam*np.eye(3)-A0-A1*np.exp(-lam*tau)
def phaseD(p,e,A0,A1,Q=np.diag([1.,1.,0.])):
    D=delta(0j,A0,A1,p.tau);_,_,vh=np.linalg.svd(D.conj().T);u=vh.conj().T[:,-1]
    v=np.array([0,e.A,0],complex);dp=np.eye(3)+p.tau*A1
    return .5*np.vdot(u,Q@u).real/abs(np.vdot(u,dp@v))**2

def rin_var(p,e,A0,A1,omega):
    B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    om=omega;DD=1j*om[:,None,None]*np.eye(3)[None]-A0[None]-np.exp(-1j*om*p.tau)[:,None,None]*A1[None]
    H=np.linalg.solve(DD,np.broadcast_to(B,(len(om),)+B.shape));C=np.array([[2/e.A,0,0.]])
    G=np.einsum('ab,wbc->wac',C,H);ps=np.real(np.einsum('wac,wac->w',G,np.conj(G)))
    return np.trapezoid(ps,om)/np.pi

def scan(label,P=1,eps=.005,eta=0):
    alpha=4.;tau=20.;omega=np.concatenate([np.linspace(1e-6,.5,1001),np.geomspace(.5,30,300)[1:]])
    p0=GCParams(alpha=alpha,eps=eps,P=P,kappa=0,tau=tau,phi=0,gain_compression=eta);e0=solve_ecm_gc(p0);a0,a1=matrices_gc(p0,e0);D0=phaseD(p0,e0,a0,a1);R0=rin_var(p0,e0,a0,a1,omega)
    kappas=np.linspace(0,.009,13);phis=np.linspace(-np.pi,np.pi,37)
    stable=win=0;best=None
    for k in kappas:
      for ph in phis:
       if k==0:continue
       p=GCParams(alpha=alpha,eps=eps,P=P,kappa=float(k),tau=tau,phi=float(ph),gain_compression=eta)
       try:
        e=solve_ecm_gc(p);A0,A1=matrices_gc(p,e);vals=eigvals(chain(A0,A1,tau));vals=vals[np.abs(vals)>1e-5];m=-np.max(vals.real)
        if m<=0:continue
        stable+=1;dr=phaseD(p,e,A0,A1)/D0;rr=rin_var(p,e,A0,A1,omega)/R0
        if dr<1 and rr<1:
          win+=1;score=(1-dr)+(1-rr)
          if best is None or score>best[0]:best=(score,k,ph,dr,rr,m)
       except Exception:pass
    return [label,P,eps,eta,stable,win,win/stable if stable else np.nan]+([] if best is None else list(best[1:]))

def main():
    cfg=[('P=0.7',.7,.005,0),('P=1.0',1,.005,0),('P=1.5',1.5,.005,0),('eps=0.003',1,.003,0),('eps=0.008',1,.008,0),('gain compression=0.05',1,.005,.05),('gain compression=0.10',1,.005,.10)]
    rows=[]
    for c in cfg:
      print('scan',c);r=scan(*c);rows.append(r);print(r)
    hdr=['configuration','P','eps','gain_compression','stable_nodes','winwin_nodes','winwin_fraction','best_kappa','best_phi','best_Dratio','best_RINratio','best_margin']
    with open(BASE/'data'/'extended_robustness.csv','w',newline='') as f:w=csv.writer(f);w.writerow(hdr);w.writerows(rows)
    print(open(BASE/'data'/'extended_robustness.csv').read())
if __name__=='__main__':main()
