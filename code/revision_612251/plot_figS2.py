#!/usr/bin/env python3
"""Regenerate Supplementary Fig. S2 from cached full-resolvent data.

python plot_figS2.py                    # redraw with included 1440-point CSVs (seconds)
python plot_figS2.py --recompute        # recompute all ECMs from rev_core.py (~3 min)

The two cases match the representative nonlinear SDE tests: (tau,kappa)=(200,.005)
and (1000,.001). Stable ECMs require F>0 and count_right=0. The plot is a
new reconstruction because the original submitted files lacked Fig. S2.
"""
from pathlib import Path
import argparse,csv,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'FigS2_reconstruction'
p=argparse.ArgumentParser();p.add_argument('--recompute',action='store_true');a=p.parse_args()
OUT.mkdir(exist_ok=True)
fig,axes=plt.subplots(1,2,figsize=(11.7,4.4),sharey=True)
for ax,(tau,k,lo,hi) in zip(axes,[(200,.005,-.0893,.0785),(1000,.001,-.0836,.0781)]):
 fn=OUT/f'figS2_points_tau{tau}.csv'
 if a.recompute:
  from lk.rev_core import Model
  m=Model(tau);th=np.linspace(-np.pi,np.pi,1440,endpoint=False)
  rows=[]
  for t in th:
   s=m.state(t,k,stability=True)
   rows.append((tau,k,t/np.pi,s['phi']/np.pi,s['F'],s['nunst'],s['Dratio'],s['Rratio']))
  with fn.open('w',newline='') as f:
   w=csv.writer(f);w.writerow(['tau','kappa','theta_over_pi','phi_over_pi','F','nunst','linewidth_ratio','rin_ratio']);w.writerows(rows)
  arr=np.array(rows)
 else:
  if not fn.exists():sys.exit(f'Missing {fn}; use --recompute or include the released data.')
  with fn.open() as f: arr=np.array([[float(x) for x in r.values()] for r in csv.DictReader(f)])
 valid=(arr[:,4]>0)&(arr[:,5]==0)
 ax.axvspan(lo,hi,alpha=.12,color='tab:green',label='First-order joint window')
 ax.scatter(arr[valid,2],arr[valid,6],s=9,marker='.',alpha=.8,label='Linewidth ratio')
 ax.scatter(arr[valid,2],arr[valid,7],s=10,marker='x',alpha=.72,label='RIN ratio')
 ax.axhline(1,color='.3',lw=.9,ls='--')
 ax.set(xlim=(-1,1),ylim=(0,2.2),xlabel=r'Total return phase $\theta/\pi$',title=rf'$\tau={tau}$, $\kappa={k:g}$ ($\kappa\tau={tau*k:g}$)')
 ax.grid(alpha=.13)
axes[0].set_ylabel('Ratio to solitary laser (linear delay resolvent)')
axes[1].legend(loc='upper right',fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'figS2_long_cavity.pdf',bbox_inches='tight')
print('Regenerated',OUT/'figS2_long_cavity.pdf')
