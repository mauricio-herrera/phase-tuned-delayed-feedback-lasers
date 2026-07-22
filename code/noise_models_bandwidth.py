import numpy as np, csv
from pathlib import Path
from lk_model import Params,solve_ecm,phase_diffusion,output_variance_grid,default_omega_grid
BASE=Path(__file__).resolve().parents[1]
(BASE/'new_results').mkdir(exist_ok=True)
(BASE/'data').mkdir(exist_ok=True)

def integrate_band(om,psd,wmax):
    m=om<=wmax
    return np.trapezoid(psd[m],om[m])/np.pi

def main():
    cases=[('solitary',Params(kappa=0,tau=20,phi=0)),('tradeoff',Params(kappa=.008,tau=20,phi=0)),('win_win',Params(kappa=.009,tau=20,phi=np.pi/4))]
    omega=default_omega_grid()
    Bopt=np.array([[1.,0.],[0.,1.],[0.,0.]])
    BN=np.array([[0.],[0.],[1.]])
    p0=cases[0][1];e0=solve_ecm(p0)
    vopt=output_variance_grid(p0,e0,[2/e0.A,0,0],B=Bopt,omega=omega)
    vN=output_variance_grid(p0,e0,[2/e0.A,0,0],B=BN,omega=omega)
    rho=vopt/vN
    Bmix=np.column_stack([Bopt,np.sqrt(rho)*BN])
    models=[('field_spontaneous',Bopt),('carrier_current',BN),('balanced_mixed',Bmix)]
    rows=[]
    for model,B in models:
      Q=B@B.T
      refs={}
      for name,p in cases:
        e=solve_ecm(p)
        D,_=phase_diffusion(p,e,Q=Q)
        V,om,psd=output_variance_grid(p,e,[2/e.A,0,0],B=B,omega=omega,return_psd=True)
        # RIN output is C=(2/A,0,0), so no extra A^4 factor is needed
        refs[name]=(D,V,om,psd)
      for name,p in cases:
        D,V,om,psd=refs[name]
        rows.append([model,name,D/refs['solitary'][0],V/refs['solitary'][1]])
    with open(BASE/'data'/'noise_model_ratios.csv','w',newline='') as f:
      w=csv.writer(f);w.writerow(['noise_model','case','phase_diffusion_ratio','integrated_RIN_ratio']);w.writerows(rows)
    # bandwidth ratios for field noise
    refs={}
    for name,p in cases:
      e=solve_ecm(p); V,om,psd=output_variance_grid(p,e,[2/e.A,0,0],B=Bopt,omega=omega,return_psd=True)
      refs[name]=(V,om,psd)
    bands=[0.15,0.25,0.5,3.0,200.0]
    brows=[]
    for wmax in bands:
      base=integrate_band(refs['solitary'][1],refs['solitary'][2],wmax)
      for name,p in cases:
        val=integrate_band(refs[name][1],refs[name][2],wmax)
        brows.append([wmax,name,val/base])
    with open(BASE/'data'/'bandwidth_ratios.csv','w',newline='') as f:
      w=csv.writer(f);w.writerow(['omega_max','case','RIN_ratio']);w.writerows(brows)
    # dimensional illustration calibrated to a 2 MHz solitary linewidth and photon lifetime 2.5 ps
    tau_p=2.5e-12; linewidth0=2e6
    D0,_=phase_diffusion(p0,e0)
    sigma2=np.pi*linewidth0*tau_p/D0
    # relaxation peak and RIN dBc/Hz for each case
    drows=[]
    for name,p in cases:
      e=solve_ecm(p)
      D,_=phase_diffusion(p,e)
      V,om,psd=output_variance_grid(p,e,[2/e.A,0,0],B=Bopt,omega=omega,return_psd=True)
      m=(om>0.02)&(om<0.3);idx=np.argmax(psd[m]); wpk=om[m][idx]; spk=psd[m][idx]
      fpk=wpk/(2*np.pi*tau_p)
      rin_hz=tau_p*sigma2*spk
      drows.append([name,D/D0*linewidth0,fpk,10*np.log10(rin_hz),V*sigma2])
    with open(BASE/'data'/'device_scale_illustration.csv','w',newline='') as f:
      w=csv.writer(f);w.writerow(['case','linewidth_Hz','RIN_peak_frequency_Hz','RIN_peak_dBc_per_Hz','integrated_RIN_variance']);w.writerows(drows)
    with open(BASE/'new_results'/'noise_bandwidth_summary.txt','w') as f:
      f.write(f'carrier scaling rho for equal solitary integrated RIN = {rho:.10e}\n')
      f.write(open(BASE/'data'/'noise_model_ratios.csv').read()+'\n')
      f.write(open(BASE/'data'/'bandwidth_ratios.csv').read()+'\n')
      f.write(f'dimensional calibration sigma^2={sigma2:.10e}, tau_p={tau_p}, linewidth0={linewidth0}\n')
      f.write(open(BASE/'data'/'device_scale_illustration.csv').read())
    print(open(BASE/'new_results'/'noise_bandwidth_summary.txt').read())
if __name__=='__main__':main()
