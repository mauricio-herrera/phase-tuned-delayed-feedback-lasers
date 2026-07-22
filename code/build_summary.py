import numpy as np, pandas as pd
from lk_model import Params,solve_ecm,phase_diffusion,output_variance_quad
from lk_analysis import leading_nonzero_root, characteristic_roots, pole_residue
from pair_formula import modal_variance

cases=[('solitary',0.0,0.0),('tradeoff',0.008,0.0),('win_win',0.009,np.pi/4)]
B=np.array([[1.,0.],[0.,1.],[0.,0.]])
omega_band=np.concatenate([np.geomspace(1e-5,0.02,500),np.linspace(.02,.25,3000),np.geomspace(.25,5,1200)])
def pair_band_fraction(p,e,pair):
    vi,om,psd=__import__('lk_model').output_variance_grid(p,e,[2*e.A,0,0],omega=omega_band,return_psd=True)
    H=np.zeros((len(om),3,2),complex)
    for lam in pair:
        RB,_=pole_residue(lam,p,e,B)
        H += RB[None]/(1j*om[:,None,None]-lam)
    C=np.atleast_2d([2*e.A,0,0])
    G=np.einsum('ab,wbc->wac',C,H)
    pair_psd=np.real(np.einsum('wac,wac->w',G,np.conj(G)))
    rin=psd/e.A**4; pair_rin=pair_psd/e.A**4
    m=(om>0.02)&(om<0.25); pk=om[m][np.argmax(rin[m])]
    band=(om>pk-0.04)&(om<pk+0.04)
    return np.trapezoid(pair_rin[band],om[band])/np.trapezoid(rin[band],om[band])

p0=Params(kappa=0,tau=20,phi=0);e0=solve_ecm(p0)
D0,_=phase_diffusion(p0,e0)
vi0=output_variance_quad(p0,e0,[2*e0.A,0,0])[0]
vn0=output_variance_quad(p0,e0,[0,0,1])[0]
rows=[]
for label,k,phi in cases:
 p=Params(kappa=k,tau=20,phi=phi);e=solve_ecm(p)
 D,_=phase_diffusion(p,e)
 F=1+k*p.tau*(np.cos(e.theta)-p.alpha*np.sin(e.theta))
 vi=output_variance_quad(p,e,[2*e.A,0,0])[0]
 vn=output_variance_quad(p,e,[0,0,1])[0]
 lead,roots=leading_nonzero_root(p,e,M=160)
 pair=np.array([lead.real+1j*abs(lead.imag),lead.real-1j*abs(lead.imag)])
 pvi=modal_variance(pair,p,e,[2*e.A,0,0])
 pvn=modal_variance(pair,p,e,[0,0,1])
 bandfrac=1.0 if k==0 else pair_band_fraction(p,e,pair)
 rows.append(dict(case=label,kappa=k,phi=phi,theta=e.theta,Omega=e.omega,Ns=e.Ns,A=e.A,
                  lead_real=lead.real,lead_imag=abs(lead.imag),margin=-lead.real,
                  chirp_factor=F,chirp_reduction=1/F**2,
                  phase_diffusion=D,phase_diffusion_ratio=D/D0,
                  intensity_variance=vi,intensity_variance_ratio=vi/vi0,
                  RIN_variance=vi/e.A**4,RIN_ratio=(vi/e.A**4)/(vi0/e0.A**4),
                  carrier_variance=vn,carrier_ratio=vn/vn0,
                  pair_fraction_intensity=pvi/vi,pair_fraction_RIN_band=bandfrac,pair_fraction_carrier=pvn/vn))
df=pd.DataFrame(rows)
df.to_csv(str((__import__('pathlib').Path(__file__).resolve().parents[1]/'data'/'selected_cases.csv')),index=False)
print(df.to_string(index=False))
