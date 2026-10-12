#!/usr/bin/env python3
"""Reproduce Table S11/S12: phase/path tolerance and principal ECM continuation."""
import numpy as np
from scipy.optimize import brentq
from lk.rev_core import Model
from rev_analysis import export_csv,phase_edges,phase_map
m=Model(20);kt=.18;phi=.25*np.pi
low,high=phase_edges(m,kt/20)
rows=[]
for edge,th in [('RIN lower edge',low),('linewidth upper edge',high)]:
 dphi=phase_map(th,kt)-phi
 for wl in (850,1310,1550):
  rows.append(dict(edge=edge,wavelength_nm=wl,delta_phi_pi=dphi/np.pi,
                   mirror_displacement_nm=abs(wl*dphi/(4*np.pi))))
export_csv('tolerances.csv',rows)
rows=[]
for kt in [.04,.08,.12,.16,.18,.20,.24,.28,.32,.36,.40]:
 k=kt/20
 theta=brentq(lambda t:phase_map(t,kt)-phi,-.7*np.pi,.7*np.pi)
 s=m.state(theta,k,stability=True,carrier=True)
 from lk.rev_core import mats,margin
 A0,A1=mats(s,k)
 marg,capped=margin(A0,A1,20)
 rows.append(dict(kt=kt,theta_pi=theta/np.pi,F=s['F'],linewidth_ratio=s['Dratio'],rin_ratio=s['Rratio'],carrier_ratio=s['Nratio'],stability_margin=marg,capped=capped,nunst=s['nunst']))
export_csv('kappa_continuation.csv',rows)
print('Central phase tolerance phi/pi',phase_map(low,.18)/np.pi,phase_map(high,.18)/np.pi)
