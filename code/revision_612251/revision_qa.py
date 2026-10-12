#!/usr/bin/env python3
"""Quality gate vs journal Table S2,S3,S8,S10-S13 reference values.
Expected values are stored only here, never used in analytical calculations.
"""
import csv,math
from pathlib import Path
O=Path(__file__).resolve().parent/'revision_outputs'
def read(file):
    with (O/file).open() as f:return list(csv.DictReader(f))
def close(a,b,tol,label):
    if not math.isclose(float(a),float(b),rel_tol=0.,abs_tol=tol):raise AssertionError(f'{label}: {a} vs {b} tol {tol}')
# Supplement Table S8 exact counts (all nine parameter combinations).
expect={(2,10):(435,435,435,0,0),(2,20):(374,374,374,0,0),(2,40):(313,313,313,0,0),
 (4,10):(245,245,245,0,0),(4,20):(208,208,208,0,0),(4,40):(162,108,162,6,5),
 (6,10):(171,171,171,0,0),(6,20):(137,121,137,0,0),(6,40):(108,48,80,41,5)}
audit=read('ecm_audit_summary.csv')
assert len(audit)==9
for row in audit:
 k=(int(float(row['alpha'])),int(float(row['tau'])))
 actual=tuple(int(row[n]) for n in ('total','C_lt_1','principal','multi_ecm','multistable'))
 assert actual==expect[k],(k,actual,expect[k])
# Table S2 analytic harmonic coefficients
coefs={10:(-49.603,-116.042),20:(-106.862,-344.420),40:(-96.128,-408.912),200:(-67.746,-235.109),1000:(-71.606,-266.454)}
for r in read('firstorder.csv'):
 a,b=coefs[int(r['tau'])];close(r['a_R'],a,.003,'gR cos');close(r['b_R'],b,.003,'gR sin')
for r in read('mechanism.csv'):
 d,n,p,c,prod=((1.902,.773,1.057,.986,1.533) if r['case']=='tradeoff' else (.486,1.131,.830,.974,.444))
 for k,target in [('damping',d),('normalization',n),('projection',p),('cross',c),('product',prod)]:close(r[k],target,.001,'mech '+k)
# Supplement S10, tolerance and kappa continuation
for r in read('windows_tau20.csv'):
 if abs(float(r['kt'])-.18)<1e-9:
  close(r['phi_low_pi'],.0792,.0001,'phi lower');close(r['phi_high_pi'],.3198,.0001,'phi upper')
for r in read('tolerances.csv'):
 if int(r['wavelength_nm'])==1550:
  close(r['mirror_displacement_nm'],(66.2 if r['edge']=='RIN lower edge' else 27.0),.2,'mirror tolerance')
# Distinguish S13 grid edges from refined roots.
exp_long={(200,.001):(1,.692,.938),(200,.002):(2,.507,.884),(200,.005):(1,.246,.761),(200,.009):(1,.124,.650),
(1000,.0005):(2,.444,.966),(1000,.001):(2,.249,.935),(1000,.002):(3,.110,.880)}
for r in read('longcavity.csv'):
 t=int(r['tau']);k=float(r['kappa']);a,b,c=exp_long[(t,k)]
 assert int(r['stable_ecms'])==a,(t,k)
 close(r['linewidth_ratio'],b,.001,'long linewidth');close(r['rin_ratio'],c,.001,'long rin')
print('PASS: first-order (5), mechanism (2), ECM audit (9), tau20 boundaries/tolerances and long cavity (7).')
