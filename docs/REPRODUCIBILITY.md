# Reproducibility map

## Final publication figures

| Output | Generating script | Principal inputs |
|---|---|---|
| `fig1_model_spectra.pdf` | `code/make_figures.py` | model equations and selected cases |
| `fig2_partition_map.pdf` | `code/make_figures.py` | `data/map_tau20_full.npz` |
| `fig3_phase_criterion_validation.pdf` | `code/make_josab_revision_figures.py` | weak-feedback and SDE summary CSV files |
| `fig4_residue_reduction.pdf` | `code/make_figures.py`, renamed by `reproduce.py` | full-resolvent and characteristic-root calculations |
| `fig5_practical_robustness.pdf` | `code/make_josab_revision_figures.py` | robustness, bandwidth, source-model, and device-scale CSV files |
| `figS_fixed_output.pdf` | `code/make_josab_revision_figures.py` | `data/fixed_output_power.csv` |

## Main numerical stages

1. `make_maps.py` creates the coarse delay-resolvent map.
2. `refine_map_stability.py` refines characteristic roots and stability margins.
3. `pair_map.py` adds relaxation-pair variance fractions.
4. `weak_feedback.py` computes the first-order phase sensitivities and finite-difference audit.
5. `robustness_scan.py` scans linewidth-enhancement factor and delay.
6. `build_summary.py` creates the selected-case table.
7. `root_audit.py` performs the transport-chain and argument-principle audit.
8. `noise_models_bandwidth.py`, `fixed_output_power.py`, and `extended_robustness.py` produce the practical robustness tables.
9. `simulate_nonlinear_dde.py` generates the nonlinear stochastic delay simulations.

The repository includes the numerical outputs used in the manuscript, so the default command regenerates figures without repeating the most expensive scans and stochastic simulations.

## Nonlinear simulation configurations represented in `data/`

| Tag | Noise amplitude | Ensembles | Time step | Total time |
|---|---:|---:|---:|---:|
| `main` | 0.001 | 32 | 0.02 | 2000 |
| `smallnoise` | 0.0005 | 64 | 0.02 | 3000 |
| `final` | 0.001 | 64 | 0.02 | 3000 |
| `long` | 0.001 | 24 | 0.02 | 8000 |

The main paper's direct-validation figure uses the `smallnoise` set.
