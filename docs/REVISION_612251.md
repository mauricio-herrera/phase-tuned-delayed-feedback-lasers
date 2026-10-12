# JOSA B manuscript 612251 — reproducibility of the revised analysis

Revision prepared 11–12 October 2026 for *Phase-tuned delayed feedback for simultaneous linewidth and relative intensity noise suppression in semiconductor lasers*.

This update adds **seven numerical revision entrypoints** under [`code/revision_612251/`](../code/revision_612251/), the source model (`lk/`), the verified nonlinear stochastic integration scripts, and tabulated outputs under [`data/revision_612251/`](../data/revision_612251/). It supplements, rather than overwrites, the initial submission reproduction code.

The complete archived results, including the 72 production and 3 smoke stochastic `.npz` batches, are in [`JOSAB_612251_Revision_Reproducibility.zip`](../data/revision_612251/JOSAB_612251_Revision_Reproducibility.zip). The analytic audit and processed CSVs are provided separately for convenient inspection.

## Reproduction

```bash
cd code/revision_612251
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --only-binary=:all: -r requirements.txt
# Regenerate first-order, finite-feedback, mechanisms, and ECM checks
bash run_revision_analyses.sh
python revision_qa.py
# Reproduce a short stochastic smoke test
python run_sde.py --preset smoke --backend numba
```

Full SDE reruns are defined by `run_all.sh` and `config/cases.json`. Numerical comparisons and expected uncertainty are described in [`EXPECTED_RESULTS_612251.md`](EXPECTED_RESULTS_612251.md). Note that reproducibility of seeded computational trajectories is not an independent experimental measurement.

## Availability and provenance

- The existing project's original Zenodo record remains an earlier archive; a **new Zenodo version has not yet been verified**. Do not cite an unconfirmed new-version DOI.
- Supplemental Figure S2 is available as `data/revision_612251/figS2_long_cavity.pdf`.
- See [`REVIEWER_RESPONSE_CROSSCHECK_612251.md`](REVIEWER_RESPONSE_CROSSCHECK_612251.md) for the analysis-to-reviewer cross-check.

Scientific repository: https://github.com/mauricio-herrera/phase-tuned-delayed-feedback-lasers
