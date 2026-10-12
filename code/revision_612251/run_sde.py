#!/usr/bin/env python3
"""Long nonlinear stochastic Lang-Kobayashi runs for the JOSA B revision (manuscript 612251).

Each batch of realizations is written to results/raw/<preset>/<case>_bNNN.npz and is never
recomputed, so a run can be interrupted and resumed, or extended by raising --batches.
Raw files store, per realization, the phase mean-square displacement on a dense lag grid
plus the intensity/carrier variances, so the lag window can be changed in aggregate.py
without rerunning.

Examples
  python run_sde.py --preset smoke
  python run_sde.py --preset long200
  python run_sde.py --preset long1000 --batches 16          # extend an existing run
  python run_sde.py --preset long200 --dt 0.01 --tag dt001  # time-step check
"""
import argparse, json, hashlib, time, sys, os
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'lk'))
from rev_core import ecm_theta
from sde_kernel import integrate, HAVE_NUMBA

ALPHA, EPS, P, SIGMA = 4.0, 0.005, 1.0, 5e-4


def seed_for(*parts):
    return int(hashlib.sha256('|'.join(map(str, parts)).encode()).hexdigest()[:8], 16)


def lag_grid(rec, T):
    return np.unique(np.rint(np.geomspace(rec, 0.25 * T, 60) / rec).astype(int)) * rec


def run_batch(cfg, preset, case, b, dt, backend, tag, sigma=SIGMA):
    kappa, theta_pi = cfg['cases'][case]
    tau = cfg['tau']; th = theta_pi * np.pi
    e = ecm_theta(th, kappa, tau, ALPHA, P)
    seeds = [seed_for(preset, tag, case, b, j) for j in range(cfg['batch_size'])]
    t0 = time.time()
    q, st, rdt = integrate(e, kappa, tau, ALPHA, EPS, P, sigma, dt, cfg['burn'], cfg['T'], cfg['rec'], seeds, backend)
    lags = lag_grid(rdt, cfg['T'])
    Ls = np.rint(lags / rdt).astype(int)
    # one lag at a time to keep memory at ~2x the phase record (long records would otherwise
    # need tens of GB). The mean increment per lag measures a residual mean-frequency offset,
    # which adds (dw*t)^2 to the MSD; aggregate.py removes it.
    msd = np.empty((len(Ls), q.shape[1])); mean_inc = np.empty_like(msd)
    for i, l in enumerate(Ls):
        x = q[l:] - q[:-l]
        msd[i] = np.mean(x * x, axis=0); mean_inc[i] = np.mean(x, axis=0)
    del x
    n = st[:, 4]
    varI = st[:, 1] / n - (st[:, 0] / n) ** 2
    varN = st[:, 3] / n - (st[:, 2] / n) ** 2
    # drift check: mean of |E|^2 relative to the ECM value (mode hopping shows up here)
    meanI_rel = (st[:, 0] / n) / e['A'] ** 2
    return dict(lags=lags, msd=msd, mean_inc=mean_inc, rin=varI / e['A'] ** 4, varN=varN, meanI_rel=meanI_rel,
                A=e['A'], tau=tau, kappa=kappa, theta=th, dt=dt, T=cfg['T'], burn=cfg['burn'],
                sigma=sigma, seeds=np.array(seeds), seconds=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--preset', required=True)
    ap.add_argument('--config', default=str(HERE / 'config' / 'cases.json'))
    ap.add_argument('--batches', type=int, help='override number of batches')
    ap.add_argument('--batch-size', type=int, help='override realizations per batch')
    ap.add_argument('--T', type=float, help='override record length')
    ap.add_argument('--dt', type=float, default=0.02)
    ap.add_argument('--sigma', type=float, default=SIGMA, help='noise amplitude (default 5e-4, as in the paper)')
    ap.add_argument('--cases', nargs='*', help='subset of cases')
    ap.add_argument('--backend', default='auto', choices=['auto', 'numba', 'numpy'])
    ap.add_argument('--tag', default='main', help='label for variant runs (e.g. dt001); results go to results/raw/<preset>_<tag>')
    ap.add_argument('--out', default=str(HERE / 'results' / 'raw'))
    a = ap.parse_args()
    cfg = json.load(open(a.config))[a.preset]
    for k, v in [('batches', a.batches), ('batch_size', a.batch_size), ('T', a.T)]:
        if v is not None: cfg[k] = v
    if a.backend == 'auto' and not HAVE_NUMBA:
        print('AVISO: numba no esta instalado; el backend numpy es ~30-100x mas lento.')
    out = Path(a.out) / (a.preset if a.tag == 'main' else f'{a.preset}_{a.tag}'); out.mkdir(parents=True, exist_ok=True)
    used = dict(cfg, dt=a.dt, tag=a.tag, sigma=a.sigma)
    prev = out / 'config_used.json'
    if prev.exists():
        old = json.load(open(prev))
        for k in ['tau', 'burn', 'T', 'rec', 'batch_size', 'dt', 'sigma']:
            if k == 'sigma' and k not in old: old[k] = SIGMA
            if old.get(k) != used.get(k):
                sys.exit(f'ERROR: {out} ya contiene lotes con {k}={old.get(k)}; ahora {k}={used.get(k)}. '
                         f'No se pueden mezclar. Use --tag para una variante (p.ej. --tag T2e5).')
    json.dump(used, open(prev, 'w'), indent=1)
    cases = a.cases or list(cfg['cases'])
    for b in range(cfg['batches']):
        for case in cases:
            f = out / f'{case}_b{b:03d}.npz'
            if f.exists():
                continue
            r = run_batch(cfg, a.preset, case, b, a.dt, a.backend, a.tag, a.sigma)
            np.savez_compressed(f, **r)
            print(f'[{time.strftime("%H:%M:%S")}] {a.preset}/{case} batch {b+1}/{cfg["batches"]}: '
                  f'{r["seconds"]:.0f} s, RIN/A^4 mean={r["rin"].mean():.3e}, max|<I>/A^2-1|={np.abs(r["meanI_rel"]).max():.2e}', flush=True)
    print('listo:', out)


if __name__ == '__main__':
    main()
