#!/usr/bin/env python3
"""Reproduce Table S2/Supp. S2 weak-feedback harmonic coefficients without hard-coded outputs."""
import numpy as np
from lk.rev_core import Model
from rev_analysis import export_csv,fo_rin_coefficients
rows=[]
for tau in (10,20,40,200,1000):
 m=Model(tau);k=1e-7
 aR,bR=fo_rin_coefficients(tau)
 aD=-3.-2*tau;bD=8*tau
 low=-np.arctan(aR/bR)/np.pi
 high=np.arctan(-aD/bD)/np.pi
 # second-order finite difference checks (not the analytic derivative)
 f0=m.state(0.,k,stability=False); f1=m.state(np.pi/2,k,stability=False)
 row=dict(tau=tau,a_phi=aD,b_phi=bD,a_R=aR,b_R=bR,low_pi=low,high_pi=high,
          fd_R_cos=(f0['Rratio']-1)/k,fd_R_sin=(f1['Rratio']-1)/k)
 rows.append(row)
 print(tau, 'a_R,b_R',round(aR,3),round(bR,3),'interval',round(low,5),round(high,5))
export_csv('firstorder.csv',rows)
