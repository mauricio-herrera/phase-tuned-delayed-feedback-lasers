#!/usr/bin/env python3
"""Aggregate long stochastic runs and compare with the linear delay resolvent.

Writes to results/:
  summary_<preset>.csv      ratios to the simulated solitary laser, 2SE, linear predictions, z-scores
  lag_sensitivity.csv       linewidth ratios for three phase-fit windows (robustness of the estimator)
  manuscript_numbers.json   every number that update_manuscript.py inserts in the paper
  tableS14_rows.tex         rows for Supplementary Table S14
  diag_<preset>.png         local slope of the phase MSD / (2 D_lin) versus lag (should scatter about 1 inside the fit window)

Usage: python aggregate.py [--presets tau20 long200 long1000] [--raw results/raw]
"""
import argparse, json, glob, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'lk'))
from rev_core import Model
from lk_model import Params, solve_ecm, phase_diffusion

SIGMA = 5e-4
DETREND = True


def load(raw, preset):
    d = Path(raw) / preset
    if not d.exists(): return None, None
    cfg = json.load(open(d / 'config_used.json'))
    if cfg.get('tag', 'main') != 'main':
        print(f'{preset}: variante "{cfg["tag"]}" (dt={cfg["dt"]}, T={cfg["T"]})')
    out = {}
    for case in cfg['cases']:
        fs = sorted(glob.glob(str(d / f'{case}_b*.npz')))
        if not fs: continue
        z = [np.load(f) for f in fs]
        has = all('mean_inc' in x for x in z)
        if has and DETREND:
            # Remove a residual mean-frequency offset dw common to all realizations (noise-induced
            # or O(dt) shift of the ECM frequency), which would add (dw t)^2 to the MSD. The offset is
            # estimated from the ensemble-averaged mean increment, so the per-realization phase
            # diffusion itself is not subtracted.
            mi = np.concatenate([x['mean_inc'] for x in z], axis=1).mean(axis=1, keepdims=True)
            msd = np.concatenate([x['msd'] for x in z], axis=1) - mi ** 2
        else:
            if DETREND: print(f'AVISO: {preset}/{case}: lotes sin mean_inc (version antigua); MSD sin corregir.')
            msd = np.concatenate([x['msd'] for x in z], axis=1)
        drift = np.concatenate([x['mean_inc'][-1] / x['lags'][-1] for x in z]) if has else np.array([np.nan])
        out[case] = dict(lags=z[0]['lags'], msd=msd, freq_offset=drift,
                         rin=np.concatenate([x['rin'] for x in z]), varN=np.concatenate([x['varN'] for x in z]),
                         meanI=np.concatenate([x['meanI_rel'] for x in z]), kappa=float(z[0]['kappa']),
                         theta=float(z[0]['theta']), tau=float(z[0]['tau']), nb=len(fs),
                         sigma=float(z[0]['sigma']))
    return cfg, out


def fitD(lags, msd, lo, hi):
    m = (lags >= lo) & (lags <= hi)
    if m.sum() < 3: raise ValueError(f'too few lags in window [{lo},{hi}]')
    X = np.vstack([lags[m], np.ones(m.sum())]).T
    return np.linalg.lstsq(X, msd[m], rcond=None)[0][0] / 2


def mse(x):
    return x.mean(), x.std(ddof=1) / np.sqrt(len(x))


def ratio(a, sa, b, sb):
    r = a / b; return r, r * np.hypot(sa / a, sb / b)


def analyze(cfg, data, window_tau):
    tau = cfg['tau']; M = Model(tau)
    sig = next(iter(data.values()))['sigma']
    D0 = phase_diffusion(Params(kappa=0.0, tau=tau), solve_ecm(Params(kappa=0.0, tau=tau)))[0] * sig ** 2
    base = {}
    rows = []
    lo, hi = window_tau[0] * tau, window_tau[1] * tau
    for case, c in data.items():
        Dj = fitD(c['lags'], c['msd'], lo, hi)
        s = M.state(c['theta'], c['kappa'], carrier=True, stability=False)
        lin = dict(D=s['Dratio'], RIN=s['Rratio'], N=s['Nratio'])
        stats = dict(D=mse(Dj), RIN=mse(c['rin']), N=mse(c['varN']))
        absl = dict(D=D0 * s['Dratio'], RIN=sig ** 2 * M.R0 * s['Rratio'], N=sig ** 2 * M.N0 * s['Nratio'])
        if case == 'solitary': base = stats
        rows.append(dict(freq_offset=float(np.nanmean(c['freq_offset'])), case=case, kappa=c['kappa'], theta_pi=c['theta'] / np.pi, n=len(Dj), stats=stats, lin=lin,
                         absl=absl, maxdrift=float(np.abs(c['meanI']).max())))
    if not base:  # variant without solitary reference: only absolute comparisons are meaningful
        base = {q: (np.nan, np.nan) for q in ['D', 'RIN', 'N']}
    out = []
    for r in rows:
        o = dict(case=r['case'], kappa=r['kappa'], theta_pi=round(r['theta_pi'], 6), n_real=r['n'], max_meanI_drift=r['maxdrift'],
                 mean_freq_offset=r['freq_offset'])
        for q in ['D', 'RIN', 'N']:
            v, se = ratio(r['stats'][q][0], r['stats'][q][1], base[q][0], base[q][1])
            o[f'{q}_sim'], o[f'{q}_2se'], o[f'{q}_lin'] = v, 2 * se, r['lin'][q]
            o[f'{q}_z'] = (v - r['lin'][q]) / se if se > 0 else np.nan
            o[f'{q}_abs_sim_over_lin'] = r['stats'][q][0] / r['absl'][q]
            o[f'{q}_abs_2se'] = 2 * r['stats'][q][1] / r['absl'][q]
            o[f'{q}_abs_z'] = (o[f'{q}_abs_sim_over_lin'] - 1) / (o[f'{q}_abs_2se'] / 2)
        out.append(o)
    return out


def diag_plot(cfg, data, preset, window_tau):
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    except ImportError:
        return
    tau = cfg['tau']; M = Model(tau)
    sig = next(iter(data.values()))['sigma']
    D0 = phase_diffusion(Params(kappa=0.0, tau=tau), solve_ecm(Params(kappa=0.0, tau=tau)))[0] * sig ** 2
    fig, ax = plt.subplots(figsize=(6, 4))
    for case, c in data.items():
        s = M.state(c['theta'], c['kappa'], stability=False)
        # local slope of the MSD; equals 2 D at lags beyond the relaxation of all non-neutral modes
        y = np.gradient(c['msd'].mean(1), c['lags']) / (2 * D0 * s['Dratio'])
        ax.semilogx(c['lags'] / tau, y, label=case)
    ax.axvspan(*window_tau, color='0.9'); ax.axhline(1, c='k', lw=.8)
    ax.set_xlabel(r'lag / $\tau$'); ax.set_ylabel(r'$d\,\mathrm{MSD}/dt\;/\;2D_{\rm lin}$'); ax.set_ylim(0.5, 1.5)
    ax.set_title(f'{preset}: local phase-diffusion rate / linear theory'); ax.legend()
    fig.tight_layout(); fig.savefig(HERE / 'results' / f'diag_{preset}.png', dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--presets', nargs='*', default=['tau20', 'long200', 'long1000'])
    ap.add_argument('--raw', default=str(HERE / 'results' / 'raw'))
    ap.add_argument('--no-detrend', action='store_true', help='use the raw MSD (includes residual frequency drift)')
    a = ap.parse_args()
    global DETREND
    DETREND = not a.no_detrend
    res = HERE / 'results'; res.mkdir(exist_ok=True)
    nfile = res / 'manuscript_numbers.json'
    numbers = json.load(open(nfile)) if nfile.exists() else {}   # merge with earlier presets
    sfile = res / 'lag_sensitivity.csv'
    sens = []
    if sfile.exists():
        import csv
        sens = [list(r.values()) for r in csv.DictReader(open(sfile)) if r['preset'] not in a.presets]
    for p in a.presets:
        cfg, data = load(a.raw, p)
        if not data:
            print(f'{p}: sin datos, se omite'); continue
        w = cfg['lag_window_tau']
        out = analyze(cfg, data, w)
        keys = list(out[0].keys())
        with open(res / f'summary_{p}.csv', 'w') as f:
            f.write(','.join(keys) + '\n')
            for o in out: f.write(','.join(str(o[k]) for k in keys) + '\n')
        numbers[p] = dict(tau=cfg['tau'], T=cfg['T'], burn=cfg['burn'], window_tau=w,
                          n_real=out[0]['n_real'], cases={o['case']: o for o in out})
        print(f'\n== {p} (tau={cfg["tau"]:g}, T={cfg["T"]:g}, {out[0]["n_real"]} realizaciones, ajuste de fase en {w} tau) ==')
        print('razones respecto al solitario simulado: valor ± 2SE (prediccion lineal)')
        print(f'{"caso":<10}{"linewidth":>22}{"RIN":>22}{"carrier":>22}   z(D,RIN,N)   deriva<I>')
        for o in out:
            print(f'{o["case"]:<10}' + ''.join(f'{o[q+"_sim"]:>9.3f}±{o[q+"_2se"]:.3f} ({o[q+"_lin"]:.3f})' for q in ['D', 'RIN', 'N'])
                  + f'   {o["D_z"]:+.1f},{o["RIN_z"]:+.1f},{o["N_z"]:+.1f}   {o["max_meanI_drift"]:.1e}')
        print('valores absolutos simulacion/teoria lineal: valor ± 2SE')
        for o in out:
            print(f'{o["case"]:<10}' + ''.join(f'{o[q+"_abs_sim_over_lin"]:>14.3f}±{o[q+"_abs_2se"]:.3f}' for q in ['D', 'RIN', 'N']))
        for ww in ([[2, 8], [3, 12], [5, 20]] if cfg['tau'] >= 100 else [[5, 25], [7.5, 40], [15, 75]]):
            try:
                oo = analyze(cfg, data, ww)
            except ValueError:
                continue
            for o in oo:
                sens.append([p, f'{ww[0]}-{ww[1]}', o['case'], o['D_sim'], o['D_2se'], o['D_lin'], o['D_abs_sim_over_lin'], o['D_abs_2se']])
        # phase-diffusion rate in successive lag windows (absolute, relative to linear theory):
        # a value drifting away from 1 at long lags signals slow phase wander beyond linear theory
        M_ = Model(cfg['tau']); tau_ = cfg['tau']
        sig_ = next(iter(data.values()))['sigma']
        D0_ = phase_diffusion(Params(kappa=0.0, tau=tau_), solve_ecm(Params(kappa=0.0, tau=tau_)))[0] * sig_ ** 2
        wins = [(3, 12), (12, 30), (30, 60), (60, 100), (100, 200)] if tau_ >= 100 else [(7.5, 40), (40, 100), (100, 300)]
        print('tasa de difusion de fase por ventana de lags (absoluta / lineal, ±2SE):')
        with open(res / f'lag_windows_{p}.csv', 'w') as fw:
            fw.write('case,window_tau,D_abs_over_lin,two_se\n')
            for case, c in data.items():
                s_ = M_.state(c['theta'], c['kappa'], stability=False); line = f'  {case:<10}'
                for lo_, hi_ in wins:
                    if hi_ * tau_ > c['lags'].max(): continue
                    try: Dj = fitD(c['lags'], c['msd'], lo_ * tau_, hi_ * tau_) / (D0_ * s_['Dratio'])
                    except ValueError: continue
                    v, e_ = Dj.mean(), 2 * Dj.std(ddof=1) / np.sqrt(len(Dj))
                    fw.write(f'{case},{lo_}-{hi_},{v},{e_}\n'); line += f'  {lo_:g}-{hi_:g}τ: {v:.2f}±{e_:.2f}'
                print(line)
        diag_plot(cfg, data, p, w)
    with open(res / 'lag_sensitivity.csv', 'w') as f:
        f.write('preset,window_tau,case,D_ratio_sim,D_2se,D_lin,D_abs_sim_over_lin,D_abs_2se\n')
        for r in sens: f.write(','.join(str(x) for x in r) + '\n')
    json.dump(numbers, open(nfile, 'w'), indent=1, default=float)
    tex = []
    for p in ['long200', 'long1000']:
        if p not in numbers: continue
        for case, o in numbers[p]['cases'].items():
            if case == 'solitary': continue
            lab = {'tradeoff': r'trade-off, $\theta=-0.156\pi$', 'win_win': r'in-phase, $\theta=0$'}.get(case, case)
            tex.append(f"{int(numbers[p]['tau'])} & {o['kappa']:g} & {lab} & ${o['D_sim']:.3f}\\pm{o['D_2se']:.3f}$ & {o['D_lin']:.3f} & "
                       f"${o['RIN_sim']:.3f}\\pm{o['RIN_2se']:.3f}$ & {o['RIN_lin']:.3f} & ${o['N_sim']:.3f}\\pm{o['N_2se']:.3f}$ & {o['N_lin']:.3f}\\\\")
    open(res / 'tableS14_rows.tex', 'w').write('\n'.join(tex) + '\n')
    print('\nArchivos escritos en', res)


if __name__ == '__main__':
    main()
