import numpy as np
from pathlib import Path
from scipy.integrate import simpson
from lk_model import Params, solve_ecm, matrices, delta, default_omega_grid, phase_diffusion

BASE=Path(__file__).resolve().parents[1]
(BASE/'new_results').mkdir(exist_ok=True)
(BASE/'data').mkdir(exist_ok=True)

def first_order_matrices(alpha=4.0, eps=0.005, P=1.0, tau=20.0, phi=0.0):
    c=np.cos(phi); s=np.sin(phi)
    A=np.sqrt(P)
    Ns_p=-c
    A2_p=(1+2*P)*c
    A_p=A2_p/(2*A)
    A0=np.zeros((3,3),float)
    A0[0,2]=A
    A0[1,2]=alpha*A
    A0[2,0]=-2*eps*A
    A0[2,2]=-eps*(1+2*P)
    A1=np.zeros((3,3),float)
    A0p=np.zeros((3,3),float)
    # cr'=cos phi; ci'=-sin phi
    A0p[0,0]=-c; A0p[0,1]=-s; A0p[0,2]=A_p
    A0p[1,0]=s; A0p[1,1]=-c; A0p[1,2]=alpha*A_p
    # derivative of -2 eps (1+2 Ns) A
    A0p[2,0]=-2*eps*(2*Ns_p*A + A_p)
    A0p[2,2]=-4*eps*A*A_p
    A1p=np.zeros((3,3),float)
    A1p[0,0]=c; A1p[0,1]=s
    A1p[1,0]=-s; A1p[1,1]=c
    C=np.array([[2/A,0,0]],float)
    Cp=np.array([[-2*A_p/(A*A),0,0]],float)
    return A0,A1,A0p,A1p,C,Cp,A_p

def rin_sensitivity(phi,alpha=4.0,eps=0.005,P=1.0,tau=20.0,omega=None,B=None):
    if omega is None: omega=default_omega_grid()
    if B is None: B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    A0,A1,A0p,A1p,C,Cp,_=first_order_matrices(alpha,eps,P,tau,phi)
    om=np.asarray(omega)
    I=np.eye(3)
    D=1j*om[:,None,None]*I[None]-A0[None]
    BB=np.broadcast_to(B,(len(om),)+B.shape)
    R=np.linalg.inv(D)
    H=np.einsum('wij,wjk->wik',R,BB)
    K=A0p[None]+np.exp(-1j*om*tau)[:,None,None]*A1p[None]
    Hp=np.einsum('wij,wjk,wkl->wil',R,K,H)
    G=np.einsum('ab,wbc->wac',C,H)
    Gp=np.einsum('ab,wbc->wac',Cp,H)+np.einsum('ab,wbc->wac',C,Hp)
    psd=np.real(np.einsum('wac,wac->w',G,np.conj(G)))
    psdp=2*np.real(np.einsum('wac,wac->w',Gp,np.conj(G)))
    V=np.trapezoid(psd,om)/np.pi
    Vp=np.trapezoid(psdp,om)/np.pi
    return Vp/V,V,Vp

def phase_sensitivity(phi,alpha=4.0,P=1.0,tau=20.0):
    c=np.cos(phi); s=np.sin(phi)
    # derivative of log[D_phi(kappa)/D_phi(0)] at kappa=0
    return -((1+2*P)*c/P + 2*tau*(c-alpha*s))

def main():
    phis=np.linspace(-np.pi,np.pi,721)
    omega=default_omega_grid()
    g0=rin_sensitivity(0.0,omega=omega)[0]
    g90=rin_sensitivity(np.pi/2,omega=omega)[0]
    gr=g0*np.cos(phis)+g90*np.sin(phis)
    gp=np.array([phase_sensitivity(x) for x in phis])
    # harmonic fit, should be first harmonic to numerical precision
    X=np.column_stack([np.cos(phis),np.sin(phis),np.ones_like(phis)])
    br=np.linalg.lstsq(X,gr,rcond=None)[0]
    bp=np.linalg.lstsq(X,gp,rcond=None)[0]
    # intervals where both first derivatives negative
    mask=(gr<0)&(gp<0)
    # contiguous runs
    runs=[]
    start=None
    for i,m in enumerate(mask):
        if m and start is None:start=i
        if start is not None and (not m or i==len(mask)-1):
            end=i if (m and i==len(mask)-1) else i-1
            runs.append((phis[start],phis[end]))
            start=None
    out=np.column_stack([phis,gp,gr])
    np.savetxt(BASE/'data'/'weak_feedback_sensitivity.csv',out,delimiter=',',header='phi,g_phase,g_RIN',comments='')

    # finite-difference audit of the first-order formulas
    from lk_model import output_variance_grid
    ktest=1.0e-5
    pbase=Params(alpha=4.0,eps=0.005,P=1.0,kappa=0.0,tau=20.0,phi=0.0)
    ebase=solve_ecm(pbase)
    Db,_=phase_diffusion(pbase,ebase)
    Rb=output_variance_grid(pbase,ebase,[2/ebase.A,0,0],omega=omega)
    audit=[]
    for ph in [-0.25,0.0,0.20,0.50]:
        pp=Params(alpha=4.0,eps=0.005,P=1.0,kappa=ktest,tau=20.0,phi=ph)
        ee=solve_ecm(pp)
        Dp,_=phase_diffusion(pp,ee)
        Rp=output_variance_grid(pp,ee,[2/ee.A,0,0],omega=omega)
        gp_ex=(Dp/Db-1)/ktest
        gr_ex=(Rp/Rb-1)/ktest
        audit.append([ph,phase_sensitivity(ph),gp_ex,rin_sensitivity(ph,omega=omega)[0],gr_ex])
    np.savetxt(BASE/'data'/'weak_feedback_audit.csv',np.asarray(audit),delimiter=',',
               header='phi,g_phase_formula,g_phase_finite,g_RIN_formula,g_RIN_finite',comments='')
    with open(BASE/'new_results'/'weak_feedback_summary.txt','w') as f:
        f.write(f'g_phase fit = {bp[0]:.10f} cos(phi) + {bp[1]:.10f} sin(phi) + {bp[2]:.3e}\n')
        f.write(f'g_RIN fit = {br[0]:.10f} cos(phi) + {br[1]:.10f} sin(phi) + {br[2]:.3e}\n')
        f.write(f'max harmonic residual RIN = {np.max(np.abs(gr-X@br)):.3e}\n')
        f.write('win-win first-order intervals:\n')
        for a,b in runs:f.write(f'[{a:.8f},{b:.8f}] rad = [{a/np.pi:.6f},{b/np.pi:.6f}] pi\n')
    print(open(BASE/'new_results'/'weak_feedback_summary.txt').read())

if __name__=='__main__':main()
