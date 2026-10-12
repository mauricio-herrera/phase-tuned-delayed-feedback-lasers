#!/usr/bin/env python3
"""Reproduce S13: count stable ECMs at feedback phase of theta=0; full window."""
import numpy as np
from lk.rev_core import Model
from rev_analysis import export_csv,all_theta,phase_edges,phase_map
rows=[]
for tau,ks in [(200,[.001,.002,.005,.009]),(1000,[.0005,.001,.002])]:
 m=Model(tau)
 for k in ks:
  kt=k*tau
  target_phi=phase_map(0.,kt)
  sols=all_theta(target_phi,kt)
  sts=[m.state(t,k,stability=True) for t in sols]
  stable=sum(s['nunst']==0 and s['F']>0 for s in sts)
  low,high=phase_edges(m,k)
  thgrid=np.linspace(-np.pi,np.pi,1440,endpoint=False)
  # Table S13 reports sampled-grid boundaries rather than refined continuous roots.
  sampled=thgrid[(thgrid>low)&(thgrid<high)]
  gridlow,gridhigh=sampled[0]/np.pi,sampled[-1]/np.pi
  v=m.state(0.,k,stability=True)
  row=dict(tau=tau,kappa=k,kt=kt,C=kt*np.sqrt(17),ecms=len(sols),stable_ecms=stable,
         window_low_pi=low/np.pi,window_high_pi=high/np.pi,window_low_grid_pi=gridlow,window_high_grid_pi=gridhigh,linewidth_ratio=v['Dratio'],rin_ratio=v['Rratio'],
         n_unstable_inphase=v['nunst'])
  print(row);rows.append(row)
export_csv('longcavity.csv',rows)
