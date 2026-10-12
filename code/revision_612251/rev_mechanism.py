#!/usr/bin/env python3
"""Reproduce relaxation-residue factorization for Table S3."""
import numpy as np
from rev_analysis import export_csv,relaxation_residue
base=relaxation_residue(20,0,0)
rows=[]
for name,k,theta in [('tradeoff',.008,-.15584*np.pi),('win_win',.009,.01795*np.pi)]:
 r=relaxation_residue(20,k,theta)
 damping=base['lam'].real/r['lam'].real
 normal=abs(base['d']/r['d'])**2
 def proj(s): return abs(s['C']@s['v'])**2*np.linalg.norm(s['u'].conj()@s['B'])**2
 projection=proj(r)/proj(base)
 cross=(r['pair']/r['diagonal'])/(base['pair']/base['diagonal'])
 row=dict(case=name,lambda_real=r['lam'].real,lambda_imag=r['lam'].imag,damping=damping,normalization=normal,
          projection=projection,cross=cross,pair_ratio=r['pair']/base['pair'],product=damping*normal*projection*cross)
 rows.append(row);print(name,row)
export_csv('mechanism.csv',rows)
