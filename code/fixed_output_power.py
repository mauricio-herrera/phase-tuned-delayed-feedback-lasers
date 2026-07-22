import csv
import numpy as np
from pathlib import Path
from lk_model import Params, solve_ecm, phase_diffusion, output_variance_grid, default_omega_grid
BASE=Path(__file__).resolve().parents[1]

def main():
    cases=[('tradeoff',0.008,0.0),('win_win',0.009,np.pi/4)]
    omega=default_omega_grid()
    p0=Params(P=1,kappa=0,tau=20,phi=0);e0=solve_ecm(p0)
    D0,_=phase_diffusion(p0,e0)
    R0=output_variance_grid(p0,e0,[2/e0.A,0,0],omega=omega)
    N0=output_variance_grid(p0,e0,[0,0,1],omega=omega)
    rows=[]
    for name,k,phi in cases:
        # first solve the fixed-pump state to get theta and Ns, which do not depend on P
        pf=Params(P=1,kappa=k,tau=20,phi=phi); ef=solve_ecm(pf)
        # choose P so A^2=(P-Ns)/(1+2Ns)=1
        Pfix=1+3*ef.Ns
        pp=Params(P=Pfix,kappa=k,tau=20,phi=phi); ep=solve_ecm(pp)
        Dfix,_=phase_diffusion(pp,ep)
        Dpump,_=phase_diffusion(pf,ef)
        R=output_variance_grid(pp,ep,[2/ep.A,0,0],omega=omega)
        Nv=output_variance_grid(pp,ep,[0,0,1],omega=omega)
        rows.append([name,'fixed_pump',1.0,ef.A,Dpump/D0,
                     output_variance_grid(pf,ef,[2/ef.A,0,0],omega=omega)/R0,
                     output_variance_grid(pf,ef,[0,0,1],omega=omega)/N0])
        rows.append([name,'fixed_output',Pfix,ep.A,Dfix/D0,R/R0,Nv/N0])
    path=BASE/'data'/'fixed_output_power.csv'
    with open(path,'w',newline='') as f:
        w=csv.writer(f); w.writerow(['case','comparison','pump_P','amplitude_A','phase_diffusion_ratio','RIN_ratio','carrier_ratio']); w.writerows(rows)
    print(path.read_text())
if __name__=='__main__': main()
