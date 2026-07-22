#!/usr/bin/env python3
"""Reproduce publication figures and selected numerical outputs.

The default mode regenerates the five main figures and the supplementary
fixed-output-power figure from the numerical files shipped in ``data/``.
The full mode reruns the deterministic scans before regenerating figures.
Nonlinear stochastic simulations are intentionally opt-in because they are
substantially more expensive.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CODE = ROOT / "code"
FIGURES = ROOT / "figures"
DATA = ROOT / "data"
RESULTS = ROOT / "new_results"


def run(script: str) -> None:
    print(f"\n>>> {script}", flush=True)
    subprocess.run([sys.executable, str(CODE / script)], cwd=ROOT, check=True)


def prepare() -> None:
    FIGURES.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)


def deterministic_pipeline() -> None:
    # These scripts rebuild the deterministic tables and delay-resolvent maps.
    for script in [
        "make_maps.py",
        "refine_map_stability.py",
        "pair_map.py",
        "weak_feedback.py",
        "robustness_scan.py",
        "build_summary.py",
        "root_audit.py",
        "noise_models_bandwidth.py",
        "fixed_output_power.py",
        "extended_robustness.py",
    ]:
        run(script)


def nonlinear_pipeline() -> None:
    # Match the simulation configuration used for the main validation figure.
    command = (
        "from simulate_nonlinear_dde import run; "
        "run(dt=0.02, n_ens=64, T=3000.0, tag='smallnoise', sigma=5e-4)"
    )
    print("\n>>> nonlinear stochastic simulations (small-noise validation)", flush=True)
    subprocess.run([sys.executable, "-c", command], cwd=CODE, check=True)


def figure_pipeline() -> None:
    run("make_figures.py")
    run("make_josab_revision_figures.py")

    # The original general-purpose script predates the final manuscript numbering.
    # Preserve its output and also expose the final JOSAB filename.
    for ext in ("pdf", "png"):
        source = FIGURES / f"fig3_residue_reduction.{ext}"
        target = FIGURES / f"fig4_residue_reduction.{ext}"
        if source.exists():
            shutil.copy2(source, target)

    expected = [
        "fig1_model_spectra.pdf",
        "fig2_partition_map.pdf",
        "fig3_phase_criterion_validation.pdf",
        "fig4_residue_reduction.pdf",
        "fig5_practical_robustness.pdf",
        "figS_fixed_output.pdf",
    ]
    missing = [name for name in expected if not (FIGURES / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing expected figure files: {missing}")
    print("\nReproduction complete. Final figures are in figures/.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="rerun deterministic scans and tables before making figures",
    )
    parser.add_argument(
        "--simulate",
        action="store_true",
        help="rerun the nonlinear small-noise stochastic simulations",
    )
    args = parser.parse_args()

    prepare()
    if args.full:
        deterministic_pipeline()
    if args.simulate:
        nonlinear_pipeline()
    figure_pipeline()


if __name__ == "__main__":
    main()
