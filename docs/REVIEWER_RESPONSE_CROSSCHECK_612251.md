# JOSA B 612251 — reviewer-response crosscheck (11 October 2026)

## Scientific responses (7 reviewer points)

| Reviewer | Question | Revised location | Independent executable evidence | Finding |
|---|---|---|---|---|
| R1.1 | Experimental phase/strength tolerance, path length, stability | Main Sec. 4.2, Fig. 4(c), Supp. S9, Tables S11–S12 | `rev_tolerance.py`, `tolerances.csv`, `kappa_continuation.csv` | Phase interval [0.079226π, 0.319760π] at κτ=.18; margins and displacement independently reproduced |
| R1.2 | Physical origin of tradeoff-to-joint transition | Main Sec. 5, Fig. 4(d), Supp. S3 Table S3 | `rev_mechanism.py`, `mechanism.csv` | Damping factors 1.9019 (tradeoff), 0.4861 (joint); residue factors agree within rounding |
| R1.3 | First order vs finite feedback | Main Secs. 3.3, 4.1; Fig. 4(a–b), Supp. S2/S9 | `rev_firstorder.py`, `rev_windows.py` | θ-window mapped to φ and finite-κ drift reproduced |
| R2.1 | Novelty relative to separate theories | Main Introduction, Discussion, Conclusion | Argument-level review against shared operator formulas | Joint sign criterion and distinct neutral/relaxation boundaries are emphasized; theoretical novelty remains an editorial judgement |
| R2.2 | Quantitative range of weak-feedback criterion | Main Sec. 4.1, Supp. S9 Table S10 | `rev_windows.py` | All τ20 listed values reproduced to shown precision; τ200/1000 first-order agreement confirmed |
| R2.3 | ECM selection, principal/multiple branches | Main Secs. 2, 6, 8; Supp. S1, S4, S7 Table S8 | `rev_ecm_audit.py` and `consolidate_ecm.py` | **All nine rows of Table S8 reproduced exactly**, including 80 principal for α=6,τ=40 via differential continuation |
| R3 | Long-cavity validation | Main Sec. 7/Table 2; Supp. S10/Table S13/Fig S2/Table S14 | `rev_longcavity_strong.py`, `rev_sde.py`, Mac `results.zip`, `check_results.py` | Seven long-cavity configurations and stochastic checks reproduced |

## Checks and editorial caveats

- The S13 intervals are **phase-grid** boundaries from 1440 samples, whereas `rev_longcavity_strong.py` also records refined continuous crossings. These can differ by up to one grid spacing (~0.00139π); do not substitute one for the other when checking the paper.
- Main Sec. 7 says "shows no sign of mode hopping": the numerical records support no measurable persistent mean-intensity drift in initialized ECM trajectories; they do not exclude every transient switch. The subsequent sentence correctly treats device mode selection as a limitation.
- Mac original 128-realization results revalidated; τ20 joint linewidth z≈2.39 is documented. Do not describe all comparisons as within 2 SE.
- Comparison made numerically against supplied journal material. It does not constitute experimental validation, independent peer review or a guarantee of acceptance.
- Funding assertion in response changed from previously "checked" to "will be verified" until the user checks the Prism form. The claim that no external funds supported this *work* must be factually confirmed by the author.
- Main availability text distinguishes the archived original GitHub/Zenodo files and the **attached revision archive**. Before submitting, upload that code archive as a supplemental code/data file and then version the public deposits.
- Check original institutional affiliation formatting and address in Prism, and ensure author ORCID is linked if desired.
- The portal deadline should be confirmed inside Prism; 5 Oct +21 days suggests 26 Oct 2026.
