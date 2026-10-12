#!/usr/bin/env python3
"""Find full-resolvent joint-suppression boundaries in theta and laboratory phase."""
import numpy as np
from lk.rev_core import Model
from rev_analysis import export_csv,phase_edges,phase_map,fo_rin_coefficients
rows=[]
for tau in [20,200,1000]:
 m=Model(tau)
 a,b=fo_rin_coefficients(tau)
 low0=-np.arctan(a/b); high0=np.arctan((3+2*tau)/(8*tau))
 for kt in [.002,.02,.05,.10,.15,.18,.20,.24]:
  k=kt/tau
  low,high=phase_edges(m,k)
  row=dict(tau=tau,kt=kt,theta_low_pi=low/np.pi,theta_high_pi=high/np.pi,
      phi_low_pi=phase_map(low,kt)/np.pi,phi_high_pi=phase_map(high,kt)/np.pi,
      mapped_phi_low_pi=phase_map(low0,kt)/np.pi,mapped_phi_high_pi=phase_map(high0,kt)/np.pi,
      theta_lower_error_pi=(low-low0)/np.pi)
  rows.append(row)
  print(tau,kt,round(low/np.pi,5),round(high/np.pi,5),round(row['phi_low_pi'],5),round(row['phi_high_pi'],5))
export_csv('windows.csv',rows)
for tau in [20,200,1000]:
 export_csv(f'windows_tau{tau}.csv',[row for row in rows if row['tau']==tau])
