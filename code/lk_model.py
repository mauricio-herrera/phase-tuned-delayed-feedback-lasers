import numpy as np
from dataclasses import dataclass
from scipy.optimize import root_scalar, root

@dataclass
class Params:
    alpha: float = 4.0
    eps: float = 0.005
    P: float = 1.0
    kappa: float = 0.0
    tau: float = 20.0
    phi: float = 0.0

@dataclass
class ECM:
    theta: float
    omega: float
    Ns: float
    A: float


def solve_ecm(p: Params, theta_guess=None, branch_m=0):
    # theta = phi + omega*tau + 2*pi*m convention. Choose branch near phi+2*pi*m.
    base = p.phi + 2*np.pi*branch_m
    def f(th):
        return th - base + p.kappa*p.tau*(p.alpha*np.cos(th) + np.sin(th))
    if theta_guess is None:
        theta_guess = base
    sol = root(lambda z: np.array([f(z[0])]), np.array([theta_guess]), method='hybr', tol=1e-12)
    if not sol.success or abs(f(sol.x[0]))>1e-9:
        # broad bracket around branch
        xs = np.linspace(base-4*np.pi, base+4*np.pi, 4001)
        vals = np.array([f(x) for x in xs])
        idx = np.where(vals[:-1]*vals[1:]<=0)[0]
        if len(idx)==0:
            raise RuntimeError(f'ECM root failed: {sol.message}')
        # choose closest root to guess
        candidates=[]
        for i in idx:
            try:
                r=root_scalar(f, bracket=[xs[i],xs[i+1]], xtol=1e-13).root
                candidates.append(r)
            except Exception:
                pass
        th=min(candidates,key=lambda x:abs(x-theta_guess))
    else:
        th=float(sol.x[0])
    Ns=-p.kappa*np.cos(th)
    denom=1+2*Ns
    num=p.P-Ns
    if denom<=0 or num<=0:
        raise RuntimeError('Nonphysical ECM amplitude')
    A=np.sqrt(num/denom)
    omega=-p.kappa*(p.alpha*np.cos(th)+np.sin(th))
    return ECM(theta=th, omega=omega, Ns=Ns, A=A)


def matrices(p: Params, ecm: ECM):
    th=ecm.theta; A=ecm.A; Ns=ecm.Ns
    cr=p.kappa*np.cos(th)
    ci=-p.kappa*np.sin(th)  # c = cr + i ci = k exp(-i theta)
    # instantaneous field coefficient is -c by equilibrium
    A0=np.zeros((3,3),float)
    A0[0,0]=-cr; A0[0,1]=ci;  A0[0,2]=A
    A0[1,0]=-ci; A0[1,1]=-cr; A0[1,2]=p.alpha*A
    A0[2,0]=-2*p.eps*(1+2*Ns)*A
    A0[2,1]=0.0
    A0[2,2]=-p.eps*(1+2*A*A)
    A1=np.zeros((3,3),float)
    A1[0,0]=cr; A1[0,1]=-ci
    A1[1,0]=ci; A1[1,1]=cr
    return A0,A1


def delta(lam, A0,A1,tau):
    return lam*np.eye(3)-A0-A1*np.exp(-lam*tau)

def delta_prime(lam,A1,tau):
    return np.eye(3)+tau*A1*np.exp(-lam*tau)


def null_vectors(lam,A0,A1,tau):
    D=delta(lam,A0,A1,tau)
    U,s,Vh=np.linalg.svd(D)
    v=Vh.conj().T[:,-1]
    # left null u satisfies u^* D=0 => D^* u=0
    U2,s2,Vh2=np.linalg.svd(D.conj().T)
    u=Vh2.conj().T[:,-1]
    return u,v,s[-1]


def phase_diffusion(p:Params,ecm:ECM,Q=None):
    if Q is None:
        Q=np.diag([1.0,1.0,0.0])
    A0,A1=matrices(p,ecm)
    u,v,res=null_vectors(0.0+0j,A0,A1,p.tau)
    # set physical phase mode v0=(0,A,0), then scale u arbitrarily handled by quotient
    vphys=np.array([0.0,ecm.A,0.0],complex)
    # project u from SVD; denominator uses vphys
    d=np.vdot(u, delta_prime(0.0+0j,A1,p.tau)@vphys)
    num=np.vdot(u,Q@u).real
    return 0.5*num/(abs(d)**2), dict(u=u,v=vphys,d=d,res=res)


def transfer(iw,A0,A1,tau,B=None):
    if B is None:
        B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    return np.linalg.solve(delta(iw,A0,A1,tau),B)


def output_psd(omega,p:Params,ecm:ECM,C,B=None):
    A0,A1=matrices(p,ecm)
    H=transfer(1j*omega,A0,A1,p.tau,B)
    G=np.atleast_2d(C)@H
    return float(np.real(np.trace(G@G.conj().T)))

if __name__=='__main__':
    p=Params(kappa=0.008,tau=20,phi=0)
    e=solve_ecm(p)
    print(e)
    A0,A1=matrices(p,e)
    print('A0',A0,'A1',A1,sep='\n')
    D,info=phase_diffusion(p,e)
    p0=Params(kappa=0,tau=20,phi=0)
    e0=solve_ecm(p0)
    D0,_=phase_diffusion(p0,e0)
    print(D0,D,D/D0, info['d'])

def output_variance_quad(p:Params,ecm:ECM,C,B=None,epsabs=1e-8,epsrel=1e-8):
    from scipy.integrate import quad
    A0,A1=matrices(p,ecm)
    C=np.atleast_2d(np.asarray(C,float))
    if B is None:
        B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    def ftheta(x):
        # omega=tan(pi*x/2), x in (0,1)
        om=np.tan(np.pi*x/2)
        jac=(np.pi/2)*(1+om*om)
        H=np.linalg.solve(delta(1j*om,A0,A1,p.tau),B)
        G=C@H
        psd=np.real(np.trace(G@G.conj().T))
        return psd*jac/np.pi # positive freq factor 1/pi
    val,err=quad(ftheta,0,1,epsabs=epsabs,epsrel=epsrel,limit=1000,points=[0.001,0.01,0.03,0.06,0.1,0.2,0.5,0.8])
    return val,err

def default_omega_grid():
    # dense around relaxation frequency and first delay resonances; long tail for white forcing
    w1=np.linspace(1e-7,0.5,2201)
    w2=np.linspace(0.5,3.0,1201)[1:]
    w3=np.geomspace(3.0,200.0,900)[1:]
    return np.concatenate([w1,w2,w3])


def output_variance_grid(p:Params,ecm:ECM,C,B=None,omega=None,return_psd=False):
    if omega is None:
        omega=default_omega_grid()
    if B is None:
        B=np.array([[1.,0.],[0.,1.],[0.,0.]])
    A0,A1=matrices(p,ecm)
    om=np.asarray(omega)
    D=1j*om[:,None,None]*np.eye(3)[None,:,:]-A0[None,:,:]-np.exp(-1j*om*p.tau)[:,None,None]*A1[None,:,:]
    BB=np.broadcast_to(B,(len(om),)+B.shape)
    H=np.linalg.solve(D,BB)
    C2=np.atleast_2d(np.asarray(C,float))
    G=np.einsum('ab,wbc->wac',C2,H)
    psd=np.real(np.einsum('wac,wac->w',G,np.conj(G)))
    val=np.trapezoid(psd,om)/np.pi
    if return_psd:
        return val,om,psd
    return val
