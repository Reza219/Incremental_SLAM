# C++ Experimental Implementation for Incremental SLAM

This directory contains an **experimental C++ implementation** of the incremental SLAM backend developed around IGG/SPO-style update logic, online graph/state expansion, paper-style metric reporting, and sparse numerical backend experiments.

The **MATLAB/Octave code in the repository root remains the reference implementation for reproducing the paper results**. The C++ code is provided to support further development, testing, and reuse, but it has **not yet been fully validated against all paper tables on all benchmark datasets**.

This implementation should be treated as a research prototype, not as a drop-in replacement for mature sparse linear-algebra packages such as CHOLMOD, nor as a fully validated production SLAM backend.

---

## Quick status

| Item | Status |
|---|---|
| Reference paper reproduction path | MATLAB/Octave implementation in the repository root |
| C++ implementation status | Experimental research prototype |
| Online graph/state expansion | Supported via `FactorManager::ensure_state_size()` |
| Paper-style metrics | Supported: normalized chi-squared, ATE, final/mean summaries |
| Paper-style FLOP accounting | Supported: algorithmic update/solve FLOPs, including full-solve comparison fields |
| Dense Eigen backend | Supported for small/debug cases |
| CHOLMOD backend | Supported when SuiteSparse/CHOLMOD is available |
| CHOLMOD rebuild-free state-dimension expansion | Not supported; the CHOLMOD factor is marked stale and refactorized at the grown size |
| Custom sparse expanding Cholesky backend | Experimental backend for research and diagnostics |
| Exact dynamic CCOLAMD clone | Partial/research symbolic prototype; not a byte-for-byte SuiteSparse CCOLAMD replacement |
| CHOLMOD-class dynamic supernodal/multifrontal solver | Not claimed |
| Full validation against all paper tables | Not yet completed |

---

## Directory layout

Recommended layout inside this `cpp/` directory:

```text
cpp/
├── CMakeLists.txt
├── README.md
├── include/
├── src/
├── tests/
└── scripts/
    ├── run_paper_pipeline_template.sh
    ├── run_backend_compare_template.sh
    ├── summarize_paper_results.py
    └── plot_paper_curves.py
```

The `scripts/` directory should contain only user-facing benchmark, summary, plotting, and release-check helpers.

---

## Numerical backend boundaries

The repo separates three numerical paths.

### 1. CHOLMOD backend

The CHOLMOD backend is the mature sparse baseline when SuiteSparse is available.

It supports same-dimension sparse rank update/downdate operations. However, public CHOLMOD interfaces do not provide the fully rebuild-free, dimension-changing online expansion path targeted by the custom research backend. When the state dimension grows, the CHOLMOD factor is therefore marked stale and refactorized at the grown size.

### 2. Dense Eigen backend

The dense Eigen backend is intended for small systems, debugging, sanity checks, and environments where SuiteSparse is unavailable.

It is not intended for large sparse SLAM benchmarks.

### 3. Custom sparse expanding Cholesky backend

The custom backend is a repo-owned experimental sparse backend. It includes research mechanisms such as online growth, sparse suffix refactorization, affected-closure updates, dependency-cache traversal, structural-pattern classification, certification, and trace/replay diagnostics.

This backend is useful for exploring rebuild-free expansion ideas, but it is **not yet a mature CHOLMOD-class sparse solver**. Its performance, fallback behavior, and numerical robustness should be measured per dataset.

---

## Dependencies

Minimum no-SuiteSparse build:

- CMake
- C++17 compiler
- Eigen 3.4 or compatible Eigen 3.x headers

Optional sparse backend build:

- SuiteSparse/CHOLMOD

Optional development/testing tools:

- Python 3
- `pandas` and `matplotlib` for summary and plotting scripts

---

## Build without SuiteSparse

A no-SuiteSparse build should work with only Eigen available:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DISLAM_EIGEN_INCLUDE_DIR=/path/to/eigen-3.4.0 \
  -DISLAM_WITH_SUITESPARSE=OFF \
  -DISLAM_WITH_GTSAM=OFF

cmake --build build -j1
ctest --test-dir build --output-on-failure
```

Use this mode first if you only want to check that the C++ prototype builds and runs the dense/debug path.

---

## Build with SuiteSparse/CHOLMOD

A SuiteSparse/CHOLMOD build can be tested separately:

```bash
cmake -S . -B build_cholmod \
  -DCMAKE_BUILD_TYPE=Release \
  -DISLAM_EIGEN_INCLUDE_DIR=/path/to/eigen-3.4.0 \
  -DISLAM_WITH_SUITESPARSE=ON \
  -DISLAM_WITH_GTSAM=OFF

cmake --build build_cholmod -j1
ctest --test-dir build_cholmod --output-on-failure
```

If SuiteSparse is unavailable or CMake cannot locate it, use the no-SuiteSparse build above.

---

## Backend-only validation

The custom sparse backend can be built and tested without running the full SLAM benchmark pipeline:

```bash
cmake --build build --target islam_test_sparse_expanding_cholesky_randomized -j1
cmake --build build --target islam_test_sparse_expanding_cholesky_trace_replay -j1
ctest --test-dir build -R sparse_expanding --output-on-failure
```

The trace/replay test supports reproducible debugging:

```bash
./build/islam_test_sparse_expanding_cholesky_trace_replay --trace-out /tmp/islam_trace.txt
./build/islam_test_sparse_expanding_cholesky_trace_replay --replay /tmp/islam_trace.txt
./build/islam_test_sparse_expanding_cholesky_trace_replay --replay /tmp/islam_trace.txt --stop-after 50
```

These tests are useful for backend development, but they are not a substitute for full paper-pipeline validation on benchmark datasets.

---

## Demo

```bash
./build/islam_demo path/to/graph.g2o 50 10 0.3 5
```

Arguments:

```text
graph_path [max_batches] [max_gn_iter] [selective_alpha] [lc_gap]
```

---

## Benchmark executable

Real graph positional arguments:

```text
islam_run_benchmark <graph.g2o|graph.graph> <csv_out> <max_batches> <max_gn_iter> \
                    <selective_alpha> <lc_gap> <eta_threshold> <gating_mode> \
                    <use_spo> <backend> <dx_threshold>
```

Examples:

```bash
./build/islam_run_benchmark path/to/graph.g2o benchmark_dense.csv 200 10 0.3 5 1.0 IGG 1 dense-eigen 1e-3

./build/islam_run_benchmark path/to/graph.g2o benchmark_cholmod.csv 200 10 0.3 5 1.0 IGG 1 cholmod 1e-3

./build/islam_run_benchmark path/to/graph.g2o benchmark_custom.csv 200 10 0.3 5 1.0 IGG 1 custom-sparse-expanding 1e-3
```

Synthetic graph positional arguments:

```text
islam_run_benchmark synthetic <num_poses> <csv_out> <loop_period> <max_batches> \
                    <max_gn_iter> <selective_alpha> <lc_gap> <eta_threshold> \
                    <gating_mode> <use_spo> <backend> <dx_threshold>
```

The generated CSV includes backend-selection fields such as:

```text
backend_requested
backend_actual
factor_cholmod_growth_refactorizations
factor_custom_sparse_growth_updates
```

It also includes paper-style metric/FLOP fields such as:

```text
scalar_measurements
normalized_chi2
ate
solve_flops
full_solve_flops
update_flops
cumulative_solve_flops
cumulative_full_solve_flops
cumulative_update_flops
```

---

## Backend comparison helper

Use the helper script to run the three backend selections with the same algorithmic settings:

```bash
scripts/run_backend_compare_template.sh path/to/graph.g2o benchmark_outputs
```

For release validation, compare at least these configurations:

```text
1. Dense Eigen fallback/backend sanity run
2. CHOLMOD backend, when SuiteSparse is available
3. Custom sparse expanding Cholesky backend
```

The backend comparison is intended to check consistency and fallback behavior. It is not, by itself, the full paper reproduction pipeline.

---

## Paper-style reproduction pipeline

To reproduce paper-style experimental outputs, run the seven configured approaches and summarize final/mean normalized chi-squared, final/mean ATE, and algorithmic solve/update FLOPs:

```bash
BUILD_DIR=build \
BACKEND=custom-sparse-expanding \
DX_THRESHOLD=1e-3 \
ETA_THRESHOLD=1.0 \
scripts/run_paper_pipeline_template.sh path/to/graph.g2o outputs/MIT
```

The summary script writes:

```text
outputs/MIT/paper_summary.csv
outputs/MIT/paper_summary.md
```

The intended seven approaches are:

```text
GNi-SPO-IGG
GNi-SPO-LCG
GNi-SPO
GNi-IGG
GNi-LCG
GNi
GN1
```

The paper-style metrics include:

```text
final normalized chi-squared
mean normalized chi-squared
final ATE
mean ATE
mean solve FLOPs
mean update FLOPs
mean full-solve FLOPs for SPO variants
```

See:

```text
docs/PAPER_REPRODUCTION.md
```

for dataset-specific threshold settings, plotting commands, and the current validation status.

Important: this C++ pipeline is intended to reproduce the **reporting structure** of the paper. It has not yet been fully validated against every numeric entry in the paper tables on all datasets. Use the MATLAB/Octave implementation in the repository root as the reference reproduction path until that validation is complete.

---

## FLOP accounting scope

The FLOP counters are algorithmic workload estimates. They are intended to compare backend workload across methods, especially update and solve work, rather than to represent wall-clock runtime on a particular CPU.

The key reported counters are:

```text
solve_flops
full_solve_flops
update_flops
```

For SPO variants, `full_solve_flops` estimates the solve cost that would have been incurred if the same update/relinearization path had used a full solve instead of a partial solve. This is useful for isolating the benefit of the partial-solve strategy.

---

## MIT-P note

The metric pipeline can report MIT-P results if a MIT-P graph is supplied.

However, this C++ directory does not currently claim to internally synthesize the MIT-P dataset from MIT by adding reference-based position priors every 50th pose. If MIT-P reproduction is needed, prepare the corresponding graph first or use the reference MATLAB/Octave workflow.

---

## Release preflight

Run the static preflight check before local build/release validation:

```bash
scripts/preflight_release_check.sh
```

Then follow:

```text
docs/RELEASE_CHECKLIST.md
```

Before presenting the C++ implementation as validated, run at least:

```text
1. No-SuiteSparse build
2. SuiteSparse/CHOLMOD build, if available
3. CTest suite
4. Backend comparison on at least one small graph
5. Paper-pipeline run on at least one benchmark graph
6. Summary script check
7. Plot script check, if plots are to be included
```

---

## Recommended public-repo wording

For the top-level repository README, the C++ implementation can be introduced with wording such as:

```markdown
## C++ prototype

An experimental C++ implementation is provided in `cpp/`. It includes online graph expansion, backend comparison hooks, paper-style metric reporting, and algorithmic FLOP accounting. The MATLAB/Octave implementation remains the reference implementation for reproducing the paper results. The C++ implementation is provided for development and reuse, but has not yet been fully validated against all paper tables on all benchmark datasets.
```

This wording is intentionally conservative. It lets users access the C++ work without confusing it with the validated paper artifact.

---

## Citation and reference implementation

If you use this repository for reproducing the paper results, use the MATLAB/Octave code in the repository root as the reference path unless the C++ implementation has been explicitly validated for the relevant dataset and configuration.

If you use the C++ prototype, please cite the paper and clearly state that the C++ implementation is an experimental implementation.
