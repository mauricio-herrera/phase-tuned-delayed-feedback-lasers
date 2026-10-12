#!/usr/bin/env python3
"""Insert the long-run simulation results into the revised manuscript, supplement and letter.

Reads results/manuscript_numbers.json (from aggregate.py) and rewrites exactly three passages:
  main.tex             Sec. 7, the sentence block on nonlinear stochastic simulations
  supplement.tex       Sec. S10, paragraph "Nonlinear simulations" and the rows of Table S14
  response_letter.tex  Reviewer 3, item "Nonlinear validation"
Originals are never overwritten: files are written as *_updated.tex and a diff is printed.

Usage:
  python update_manuscript.py --main ../JOSAB_revised_source/main.tex \
        --supp ../Supplement_source/supplement.tex --letter ../Response_to_reviewers.tex
"""
import argparse, json, re, difflib, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent


def fmt(v, se):
    """value +/- 2SE; the uncertainty keeps two significant digits if its leading digit is 1 or 2,
    otherwise one (min. 2 decimals)."""
    import math
    if se <= 0:
        d = 2
    else:
        e = int(math.floor(math.log10(se))); lead = int(se / 10 ** e)
        d = max(2, (1 if lead < 3 else 0) - e)
    return f'${v:.{d}f}\\pm{se:.{d}f}$'


def num(x):
    """1e5 -> $10^{5}$, 4e5 -> $4\\times10^{5}$, 2000 -> 2000."""
    import math
    if x < 1e5: return f'{x:g}'
    e = int(math.floor(math.log10(x))); m = x / 10 ** e
    return f'${m:g}\\times10^{{{e}}}$' if abs(m - 1) > 1e-9 else f'$10^{{{e}}}$'


def triple(c):
    return ', '.join(fmt(c[q + '_sim'], c[q + '_2se']) for q in ['D', 'RIN', 'N'][:2]) + ', and ' + fmt(c['N_sim'], c['N_2se'])


def lin(c):
    return f"${c['D_lin']:.3f}$, ${c['RIN_lin']:.3f}$, ${c['N_lin']:.3f}$"


def checks(nums):
    ok = True
    for p, d in nums.items():
        for case, c in d['cases'].items():
            for q in ['D', 'RIN', 'N']:
                z = c.get(q + '_z')
                if z is not None and abs(z) > 3:
                    print(f'AVISO: {p}/{case}/{q}: |z|={abs(z):.1f} > 3 (simulacion vs teoria lineal). Revise antes de insertar.'); ok = False
            zabs = c.get('D_abs_z')
            if zabs is not None and abs(zabs) > 3:
                print(f'NOTA: {p}/{case}: ancho de linea absoluto simulacion/teoria = {c["D_abs_sim_over_lin"]:.3f} (|z|={abs(zabs):.1f}); '
                      'vea EXPECTED_RESULTS.md, seccion "Punto abierto".')
            if c['max_meanI_drift'] > 1e-3:
                print(f'AVISO: {p}/{case}: deriva de <|E|^2> = {c["max_meanI_drift"]:.1e}; posible salto de modo.'); ok = False
    return ok


def sub_between(text, start, ends, new, name):
    """Replace from `start` to the first of the candidate `ends` (old or already-updated text)."""
    i = text.find(start)
    if i < 0: sys.exit(f'No se encontro el inicio del pasaje "{name}". ¿Se edito ese texto a mano?')
    # candidates are tried in priority order: first the original ending, then the ending
    # written by a previous run of this script
    hits = [(text.find(e, i), e) for e in ends if text.find(e, i) >= 0]
    if not hits: sys.exit(f'No se encontro el final del pasaje "{name}".')
    j, e = hits[0]
    return text[:i] + new + text[j + len(e):]


def lag_statement(resdir):
    import csv
    rows = list(csv.DictReader(open(resdir / 'lag_sensitivity.csv')))
    worst = 0.0
    for p in ['long200', 'long1000']:
        ref = {r['case']: r for r in rows if r['preset'] == p and r['window_tau'] == '3-12'}
        for r in rows:
            if r['preset'] == p and r['case'] in ref:
                d = abs(float(r['D_ratio_sim']) - float(ref[r['case']]['D_ratio_sim']))
                worst = max(worst, d / float(ref[r['case']]['D_2se']))
    if worst <= 1:
        return r"changing this window to $2\tau$--$8\tau$ or $5\tau$--$20\tau$ changes the linewidth ratios by less than their two-standard-error uncertainty (file \texttt{lag\_sensitivity.csv})."
    print(f'AVISO: la sensibilidad a la ventana de ajuste es {worst:.1f} x 2SE; revise lag_sensitivity.csv y el texto del suplemento.')
    return rf"changing this window to $2\tau$--$8\tau$ or $5\tau$--$20\tau$ changes the linewidth ratios by up to {worst:.1f} times their two-standard-error uncertainty (file \texttt{{lag\_sensitivity.csv}})."


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--numbers', default=str(HERE / 'results' / 'manuscript_numbers.json'))
    ap.add_argument('--main'); ap.add_argument('--supp'); ap.add_argument('--letter')
    a = ap.parse_args()
    nums = json.load(open(a.numbers))
    resdir = Path(a.numbers).resolve().parent   # tables are read from the same folder
    for p in ['long200', 'long1000']:
        if p not in nums: sys.exit(f'Faltan resultados de {p}; ejecute run_sde.py y aggregate.py.')
    hop_ok = checks(nums)
    L2, L1 = nums['long200'], nums['long1000']
    t2, w2 = L2['cases']['tradeoff'], L2['cases']['win_win']
    t1, w1 = L1['cases']['tradeoff'], L1['cases']['win_win']
    nreal = min(L2['n_real'], L1['n_real'])
    hop = ('The stationary intensity shows no sign of mode hopping during the simulated records.' if hop_ok else
           'CHECK: possible mode hopping or disagreement with linear theory; see update_manuscript.py output.')

    main_new = (r"Nonlinear stochastic simulations for $\tau=200$, $\kappa=0.005$ (Supplementary Sec.~S10) confirm the classification: "
                rf"for the trade-off ECM at $\theta=-0.156\pi$ the simulated linewidth, RIN, and carrier ratios are {triple(t2)} (linear: {lin(t2)}), "
                rf"and for the in-phase ECM they are {triple(w2)} (linear: {lin(w2)}). "
                rf"For $\tau=1000$ and $\kappa=0.001$ the corresponding ratios are {triple(t1)} (linear: {lin(t1)}) and {triple(w1)} (linear: {lin(w1)}). "
                rf"Uncertainties are two standard errors over {nreal} realizations. {hop}")
    supp_par = (r"\paragraph{Nonlinear simulations.} The simulations of Sec.~S5 were repeated with the ECM specified by $\theta$, so that a selected mode is followed even when other ECMs coexist. "
                rf"The history was initialized on the selected ECM, $\Delta t=0.02$, $\sigma=5\times10^{{-4}}$, and {nreal} independent realizations were simulated per state. "
                rf"Burn-in and record length were {num(L2['burn'])} and {num(L2['T'])} for $\tau=200$ and {num(L1['burn'])} and {num(L1['T'])} for $\tau=1000$. "
                rf"Phase diffusion was obtained from the slope of the phase mean-square displacement over lags from {L2['window_tau'][0]:g}$\tau$ to {L2['window_tau'][1]:g}$\tau$, after removing a small residual offset of the mean optical frequency that is common to all realizations and would otherwise add a term quadratic in the lag; "
                + lag_statement(resdir) + " "
                r"Intensity and carrier variances were accumulated at every time step after burn-in. Table~\ref{tab:sdelong} compares the results with the linear resolvent; the solitary reference was simulated with the same settings.")
    rows = open(resdir / 'tableS14_rows.tex').read().strip()
    letter_new = (r"Nonlinear stochastic simulations that follow a selected ECM confirm the classification in the multi-ECM regime. "
                  rf"For $\tau=200$ and $\kappa=0.005$ ($C=4.1$), the trade-off ECM gives simulated linewidth, RIN, and carrier ratios of {triple(t2)} (linear: {lin(t2)}), "
                  rf"and the in-phase ECM gives {triple(w2)} (linear: {lin(w2)}). "
                  rf"For $\tau=1000$ and $\kappa=0.001$, the corresponding values are {triple(t1)} and {triple(w1)} (linear: {lin(t1)} and {lin(w1)}). "
                  rf"Uncertainties are two standard errors over {nreal} realizations with records of up to {num(L1['T'])} normalized time units.")

    jobs = []
    if a.main:
        jobs.append((a.main, [("Nonlinear stochastic simulations for $\\tau=200$, $\\kappa=0.005$ (Supplementary Sec.~S10)",
                               ["during the simulated records.", "see update_manuscript.py output."], main_new, 'main: Sec. 7')]))
    if a.supp:
        jobs.append((a.supp, [("\\paragraph{Nonlinear simulations.}", ["must exceed several delays.", "simulated with the same settings."], supp_par, 'supp: S10 paragraph'),
                              ("ECM & linewidth & lin. & RIN & lin. & carrier & lin.\\\\\n\\midrule\n", ["\\bottomrule"],
                               "ECM & linewidth & lin. & RIN & lin. & carrier & lin.\\\\\n\\midrule\n" + rows + "\n\\bottomrule", 'supp: Table S14')]))
    if a.letter:
        jobs.append((a.letter, [("Nonlinear stochastic simulations that follow a selected ECM", ["compared with $1.069$ and $0.935$.", "normalized time units."], letter_new, 'letter: R3')]))
    if not jobs:
        print(main_new, '\n\n', supp_par, '\n\n', rows, '\n\n', letter_new); return
    for path, subs in jobs:
        old = open(path).read(); new = old
        for s, e, txt, name in subs:
            new = sub_between(new, s, e, txt, name)
        out = Path(path).with_name(Path(path).stem + '_updated.tex')
        out.write_text(new)
        print(f'\n===== {path} -> {out} =====')
        sys.stdout.writelines(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=0))


if __name__ == '__main__':
    main()
