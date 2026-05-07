# Python Prototype for Incremental SLAM

This directory contains a **Python prototype** of the incremental SLAM backend associated with the paper:

> Reza Arablouei, "Efficient Incremental SLAM via Information-Guided Gating and Selective Partial Optimization," *Robotics*, 2026, 15(5), 87. https://doi.org/10.3390/robotics15050087

The implementation focuses on readable, testable versions of the main IGG/SPO algorithmic ideas:

- **Information-Guided Gating (IGG)**: compute the log-determinant information surrogate
  \(\eta_t = \sum_i \log |\rho_{t,i}| = \frac{1}{2}\log\det(J_t^T J_t)\), then use the detrended increment
  \(\Delta\eta_t = \eta_t - (N_{t-1}/N_t)\eta_{t-1}\) to decide whether an increment should trigger a global active set.
- **Selective Partial Optimization (SPO)**: solve/update only the active variables, prune converged node blocks, and expand the active set through incident graph edges.
- **Paper-style accuracy metrics**: normalized chi-squared and ATE utilities are included for experimentation and regression checks.

The **MATLAB/Octave code in the repository root remains the reference implementation for reproducing the published paper results**. This Python directory is provided for readability, experimentation, and unit-tested algorithm development. It is not intended to be a byte-for-byte reproduction of the MATLAB/SuiteSparse experimental pipeline.

---

## Status

| Item | Status |
|---|---|
| Reference paper reproduction path | MATLAB/Octave implementation in the repository root |
| Python implementation status | Research prototype |
| Main IGG/SPO algorithmic logic | Implemented |
| IGG default gating mode | Implemented via `RunConfig.gating_mode="igg"` |
| Local SPO on low-information increments | Implemented |
| Block/node-level active-set pruning and expansion | Implemented |
| Normalized chi-squared metric | Implemented using scalar measurement denominator |
| ATE utilities | Implemented |
| CHOLMOD/scikit-sparse support | Optional hybrid backend when available |
| Dynamic CCOLAMD / dynamic elimination tree clone | Not claimed |
| Full reproduction of all paper tables | Not claimed |

---

## Important limitations

This Python implementation is intentionally conservative about its claims:

- Dynamic CCOLAMD and elimination-tree maintenance are approximated through available Python/SciPy/scikit-sparse interfaces.
- The CHOLMOD path is hybrid: fixed-state factor updates use update/downdate where possible, but adding new variables still rebuilds the factorization.
- The selective partial solve uses a reduced permuted normal subsystem / reach approximation, not a fully cached dynamic block-Schur engine.
- The full MATLAB/Octave paper experiment harness, sensitivity sweeps, FLOP accounting, and all benchmark outputs are not bundled here.

A precise statement is: this directory is a **paper-aligned Python prototype of the IGG-SPO algorithmic logic**, not an exact reproducibility archive for the paper tables.

---

## Installation

### Basic installation

The basic installation does **not** require CHOLMOD/scikit-sparse:

```bash
cd python
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

On Windows PowerShell:

```powershell
cd python
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -U pip
python -m pip install -e .
```

### Install with test tools

```bash
python -m pip install -e ".[test]"
```

or:

```bash
python -m pip install -r requirements-dev.txt
```

### Optional CHOLMOD/scikit-sparse backend

The CHOLMOD path is optional because `scikit-sparse` can be harder to install, especially on Windows.

On Linux, install SuiteSparse development headers first, then install the optional dependency:

```bash
sudo apt-get update
sudo apt-get install -y libsuitesparse-dev python3-dev build-essential
python -m pip install -e ".[cholmod]"
```

Alternatively:

```bash
python -m pip install -r requirements-cholmod.txt
```

If CHOLMOD is unavailable, use `linear_backend="legacy"` or keep the default `linear_backend="auto"`, which falls back when CHOLMOD cannot be used.

---

## Quick start

Default run, using CHOLMOD automatically when available and falling back otherwise:

```bash
python scripts/run_demo.py
```

Force the CHOLMOD hybrid backend:

```bash
python scripts/run_demo_cholmod.py
```

Force the legacy sparse/dense fallback backend:

```python
from incremental_slam.config import RunConfig
from incremental_slam.main import run_selective_demo

cfg = RunConfig(dataset="MIT", linear_backend="legacy")
g_final, hist = run_selective_demo(cfg)
```

Useful history fields:

```python
hist.information_eta
hist.information_delta_eta
hist.global_update_flags
hist.active_set_size
hist.global_error          # normalized chi-squared
hist.ate_rmse
```

---

## Backend modes

`RunConfig.linear_backend` supports:

- `"auto"` — default; use CHOLMOD when available, otherwise legacy
- `"cholmod"` — require the hybrid CHOLMOD backend
- `"legacy"` — force the older refactorization-based path

`RunConfig.gating_mode` supports:

- `"igg"` — default; paper-style information-gain gating
- `"lcg"` — loop-closure gating baseline
- `"always"` — global active set at every increment

Dataset threshold defaults in `config.py` follow the paper's dataset-specific operating points where available:

| Dataset key | `tau_d` | `tau_eta` |
|---|---:|---:|
| `MIT` | `1e-3` | `1.0` |
| `MITP` | `1e-3` | `1.0` |
| `FR079` | `1e-4` | `0.6` |
| `CSAIL` | `1e-5` | `0.95` |
| `Intel` | `1e-6` | `0.72` |
| `FRH` | `1e-7` | `0.45` |

---

## Directory structure

```text
python/
├── README.md
├── pyproject.toml
├── requirements.txt
├── requirements-cholmod.txt
├── requirements-dev.txt
├── incremental_slam/
│   ├── geometry/
│   ├── graph/
│   ├── io/
│   ├── linearization/
│   ├── metrics/
│   ├── solvers/
│   └── viz/
├── scripts/
└── tests/
```

Main modules:

- `incremental_slam/io`: graph readers
- `incremental_slam/geometry`: SE(2) helpers
- `incremental_slam/graph`: graph update, reorder, loop closure, affected-set logic
- `incremental_slam/linearization`: factor linearization and Jacobian assembly
- `incremental_slam/metrics`: ATE, normalized chi-squared, and log-det information metrics
- `incremental_slam/solvers`: full/selective solvers and CHOLMOD manager
- `incremental_slam/main.py`: end-to-end drivers

---

## Tests

Run the test suite with:

```bash
python -m pytest
```

The tests cover:

- paper formula for \(\eta_t\) and \(\Delta\eta_t\)
- IGG gating based on information gain rather than loop-closure/objective-change logic
- SPO prune-expand active-set behavior
- normalized chi-squared denominator
- selective-solver smoke tests
- CHOLMOD availability flag behavior

At the time this upload package was prepared, the included tests passed in the sandbox environment.

---

## Recommended top-level README wording

If this directory is added to the main repository, the top-level README can describe it as:

```markdown
## Python prototype

A Python prototype is provided in `python/`. It implements the main IGG/SPO algorithmic logic, including information-gain gating, selective active-set refinement, normalized chi-squared, and ATE metrics. This implementation is useful for readability, experimentation, and unit-tested algorithm development.

The MATLAB/Octave implementation remains the reference implementation for reproducing the published paper results. The C++ implementation in `cpp/` is the experimental compiled backend. The Python implementation is not intended as a byte-for-byte reproduction of the MATLAB/SuiteSparse pipeline.
```

---

## Citation

If you use this Python prototype, please cite the published paper:

```bibtex
@article{Arablouei2026IncrementalSLAM,
  author  = {Arablouei, Reza},
  title   = {Efficient Incremental SLAM via Information-Guided Gating and Selective Partial Optimization},
  journal = {Robotics},
  year    = {2026},
  volume  = {15},
  number  = {5},
  article-number = {87},
  doi     = {10.3390/robotics15050087}
}
```
