#!/usr/bin/env python3
"""Audit independently supplied macOS SDE result files (JOSA B 612251).

This verifies every raw batch, seed uniqueness, ensemble sizes, finite results,
the published linear expectations and the stated acceptance conditions. It does
NOT generate new Monte Carlo trajectories; use run_sde.py for simulation and
aggregate.py for independent summary recomputation.
"""
import argparse,json,csv,math,hashlib
from pathlib import Path
import numpy as np
from lk.rev_core import Model
from rev_analysis import OUT
p=argparse.ArgumentParser();p.add_argument('--results',default='results');a=p.parse_args()
R=Path(a.results); report=[]; nraw=0
if not R.exists():raise SystemExit(f'Missing {R}; extract included Mac results')
for preset in ('tau20','long200','long1000'):
 raw=R/'raw'/preset
 config=json.loads((raw/'config_used.json').read_text())
 model=Model(config['tau']);sfn=R/f'summary_{preset}.csv'
 with sfn.open() as f:rows={r['case']:r for r in csv.DictReader(f)}
 assert set(rows)==set(config['cases']),preset
 for case,values in config['cases'].items():
  fs=sorted(raw.glob(f'{case}_b*.npz'))
  if len(fs)!=config['batches']:raise AssertionError((preset,case,'wrong batches',len(fs)))
  seeds=[];n=0
  for fn in fs:
   with np.load(fn) as z:
    if abs(float(z['tau'])-config['tau'])>1e-10 or abs(float(z['kappa'])-values[0])>1e-10:
     raise AssertionError((fn,'parameter mismatch'))
    if abs(float(z['theta'])/np.pi-values[1])>1e-7:raise AssertionError((fn,'theta mismatch'))
    for key in ('rin','varN','msd','mean_inc','meanI_rel'):
     if not np.isfinite(z[key]).all():raise AssertionError((fn,key,'nonfinite'))
    zz=[int(x) for x in z['seeds']];n+=len(zz);seeds+=zz
    if len(zz)!=config['batch_size']:raise AssertionError((fn,'batch size mismatch'))
   nraw+=1
  if len(set(seeds))!=len(seeds):raise AssertionError((preset,case,'duplicate seeds'))
  row=rows[case]
  if n!=int(float(row['n_real'])) or n!=config['batches']*config['batch_size']:
   raise AssertionError((preset,case,'inconsistent ensemble size'))
  e=model.state(values[1]*np.pi,values[0],carrier=True,stability=False)
  for obs,key in [('D','Dratio'),('RIN','Rratio'),('N','Nratio')]:
   if abs(float(row[f'{obs}_lin'])-e[key])>2e-4:raise AssertionError((preset,case,obs,'theory not reproduced'))
   z=float(row[f'{obs}_z'])
   if not math.isfinite(z) or abs(z)>3:raise AssertionError((preset,case,obs,'z out of range',z))
   if abs(z)>2:print(f'NOTE: {preset}/{case}/{obs} z={z:.2f} exceeds two SE')
  if float(row['max_meanI_drift'])>1e-3:raise AssertionError((preset,case,'drift'))
  report.append({'preset':preset,'case':case,'replicates':n,'batches':len(fs),
                 'D_sim':float(row['D_sim']),'RIN_sim':float(row['RIN_sim']),
                 'D_lin':float(row['D_lin']),'RIN_lin':float(row['RIN_lin'])})
 if len(rows)!=3:raise AssertionError('Expected three cases per preset')
assert nraw==72
(OUT/'sde_validation.json').write_text(json.dumps({'raw_batches':nraw,'cases':report},indent=2))
print(f'PASS: {nraw} production files, {len(report)} states, 128 realizations each, consistent seeds and linear predictions')
