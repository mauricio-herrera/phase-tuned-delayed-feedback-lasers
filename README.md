# Phase-tuned delayed feedback in semiconductor lasers

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23320753.svg)](https://doi.org/10.5281/zenodo.23320753)

Reproducibility package for the manuscript:

> **Phase-tuned delayed feedback for simultaneous linewidth and relative intensity noise suppression in semiconductor lasers**

This repository contains the Python code, processed numerical outputs, supplementary material, and publication figures used to study the simultaneous suppression of optical linewidth and relative intensity noise (RIN) in semiconductor lasers with delayed optical feedback.

The archived release is permanently available on Zenodo:

**DOI:** [10.5281/zenodo.23320753](https://doi.org/10.5281/zenodo.23320753)

---

## Scientific scope

The repository implements a stochastic Lang–Kobayashi delay model and includes calculations and numerical experiments for:

- phase diffusion and linewidth reduction;
- integrated relative intensity noise;
- carrier-noise contributions;
- delay-spectrum residues;
- stability and characteristic-root audits;
- weak-feedback phase criteria;
- nonlinear stochastic delay simulations;
- robustness to parameter changes;
- bandwidth and noise-source sensitivity;
- fixed-output-power comparisons.

The package is intended to support transparent verification and reproduction of the results reported in the manuscript and its supplementary material.

---

## Repository structure

```text
phase-tuned-delayed-feedback-lasers/
├── code/                  Python analysis and simulation scripts
├── data/                  Processed CSV and compressed NumPy outputs
├── figures/               Final manuscript and supplementary figures
├── supplement/            Supplementary material PDF
├── docs/                  Reproducibility documentation
├── new_results/           Output directory for recomputed results
├── reproduce.py           Main reproduction entry point
├── requirements.txt       Minimal Python dependencies
├── requirements-tested.txt
├── CITATION.cff            Citation metadata for GitHub
├── LICENSE                MIT license for source code
├── LICENSES.md            Licensing information for mixed content
├── CHANGELOG.md
└── SHA256SUMS.txt         File-integrity checksums
```

A detailed figure-to-script and result-to-script map is available in:

```text
docs/REPRODUCIBILITY.md
```

---

## Quick start

Python 3.10 or newer is recommended. NumPy 2 or newer is required because the scripts use `numpy.trapezoid`.

Clone the repository and enter its directory:

```bash
git clone https://github.com/mauricio-herrera/phase-tuned-delayed-feedback-lasers.git
cd phase-tuned-delayed-feedback-lasers
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

On Windows PowerShell, use:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

Install the required packages:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Regenerate the publication figures from the archived numerical outputs:

```bash
python reproduce.py
```

The regenerated files are written to the `figures/` directory.

---

## Reproduction modes

### Figure-only reproduction

```bash
python reproduce.py
```

This is the recommended first test. It uses the archived numerical outputs and regenerates the five main figures together with the supplementary fixed-output-power figure.

### Full deterministic recomputation

```bash
python reproduce.py --full
```

This reruns the principal deterministic analyses before regenerating the figures, including:

- the delay-resolvent parameter map;
- refined stability calculations;
- the residue-pair map;
- weak-feedback sensitivity;
- characteristic-root audits;
- robustness calculations.

This option is slower than the default figure-only workflow.

### Nonlinear stochastic validation

```bash
python reproduce.py --simulate
```

This reruns the 64-realization nonlinear stochastic delay simulation used for direct validation and then regenerates the figures. Precomputed simulation outputs are already included in `data/`.

---

## Baseline parameters and representative states

The baseline dimensionless parameters are:

```text
alpha   = 4
epsilon = 0.005
P       = 1
tau     = 20
```

The representative feedback states are:

```text
Trade-off state:
kappa = 0.008
phi   = 0

Simultaneous-suppression state:
kappa = 0.009
phi   = pi/4
```

---

## Main outputs

The repository includes the numerical data and scripts associated with:

```text
fig1_model_spectra.pdf
fig2_partition_map.pdf
fig3_phase_criterion_validation.pdf
fig4_residue_reduction.pdf
fig4_robustness.pdf
fig5_practical_robustness.pdf
figS_fixed_output.pdf
```

Some retained filenames reflect earlier internal figure numbering. The reproducibility documentation identifies the final manuscript correspondence.

---

## Data and numerical outputs

The `data/` directory contains processed outputs in CSV and compressed NumPy (`.npz`) formats. These files allow the figures and summary results to be reproduced without rerunning every computationally intensive stochastic simulation.

Examples include:

- parameter maps;
- selected operating points;
- weak-feedback audits;
- stability and root audits;
- robustness summaries;
- bandwidth and noise-model comparisons;
- nonlinear stochastic simulation ensembles.

The archive does not require access to proprietary experimental data.

---

## Integrity verification

The archive includes SHA-256 checksums in:

```text
SHA256SUMS.txt
```

On macOS or Linux, the files can be checked with:

```bash
sha256sum -c SHA256SUMS.txt
```

On macOS systems where `sha256sum` is unavailable, install GNU core utilities or verify individual files using:

```bash
shasum -a 256 path/to/file
```

---

## Citation

When using the code, numerical outputs, or figures, please cite the archived Zenodo release:

```text
Herrera-Marín, Mauricio. Phase-tuned delayed feedback for simultaneous
linewidth and relative intensity noise suppression in semiconductor lasers:
reproducibility package. Version 2.0.0. Zenodo.
https://doi.org/10.5281/zenodo.23320753
```

A BibTeX entry is:

```bibtex
@software{herrera_marin_phase_tuned_feedback_2026,
  author    = {Herrera-Marín, Mauricio},
  title     = {Phase-tuned delayed feedback for simultaneous linewidth and
               relative intensity noise suppression in semiconductor lasers:
               reproducibility package},
  version   = {2.0.0},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23320753},
  url       = {https://doi.org/10.5281/zenodo.23320753}
}
```

The repository also contains a `CITATION.cff` file, allowing GitHub to display the **Cite this repository** control.

After the associated journal article is published, its DOI should also be cited and added to the Zenodo record as a related identifier.

---

## Suggested data and code availability statement

```text
The code, processed numerical outputs, and scripts used to reproduce all
figures and supplementary analyses are publicly available in the GitHub
repository https://github.com/mauricio-herrera/phase-tuned-delayed-feedback-lasers
and in its archived Zenodo release at
https://doi.org/10.5281/zenodo.23320753.
```

LaTeX version:

```latex
\section*{Data and code availability}
The code, processed numerical outputs, and scripts used to reproduce all
figures and supplementary analyses are publicly available in the GitHub
repository
\url{https://github.com/mauricio-herrera/phase-tuned-delayed-feedback-lasers}
and in its archived Zenodo release,
\url{https://doi.org/10.5281/zenodo.23320753}.
```

---

## Licensing

The source code is distributed under the MIT License.

Author-created numerical outputs and figures are available under the Creative Commons Attribution 4.0 International license (CC BY 4.0).

See `LICENSES.md` for the treatment of the supplementary PDF and other mixed-license content.

---

## Contact

**Mauricio Herrera-Marín**  
Faculty of Engineering  
Universidad del Desarrollo  
Santiago, Chile  

Email: `mherrera@udd.cl`

Repository:  
https://github.com/mauricio-herrera/phase-tuned-delayed-feedback-lasers

Zenodo archive:  
https://doi.org/10.5281/zenodo.23320753

---

## JOSA B manuscript 612251 — October 2026 revised results

The **revised** phase-resolved analysis, including the seven new numerical validation scripts,
ECM branch enumeration, finite-feedback windows, parameter tolerances, mechanism factorization,
and long-cavity noise results, is documented in
[docs/REVISION_612251.md](docs/REVISION_612251.md).

- [Numerical code and configurations](code/revision_612251/)
- [Analytic tables, stochastic summaries, and full raw-data archive](data/revision_612251/)
- [Revised expected-results comparison](docs/EXPECTED_RESULTS_612251.md)
- [Re-review response cross-check](docs/REVIEWER_RESPONSE_CROSSCHECK_612251.md)

The revised Zenodo software release **v2.0.0** was published on **12 October 2026**, with version-specific DOI [10.5281/zenodo.23320753](https://doi.org/10.5281/zenodo.23320753). The all-versions DOI is [10.5281/zenodo.21498282](https://doi.org/10.5281/zenodo.21498282); the initial v1.0.0 release was [10.5281/zenodo.21498283](https://doi.org/10.5281/zenodo.21498283). The current DOI identifies the ZIP containing revision analyses and raw stochastic simulations.
