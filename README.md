# Incremental_SLAM

Code and data for the paper:

**Efficient Incremental SLAM via Information-Guided Gating and Selective Partial Optimization**  
Reza Arablouei  
*Robotics*, 2026, 15(5), 87  
Published article: https://www.mdpi.com/2218-6581/15/5/87  
DOI: https://doi.org/10.3390/robotics15050087

This repository contains the MATLAB/Octave implementation used to reproduce the main experimental results in the paper. The method combines **information-guided gating (IGG)** and **selective partial optimization (SPO)** for efficient incremental SLAM back-end optimization.

For reproducing the published paper results, please treat the MATLAB/Octave implementation in the repository root as the reference implementation.

## Repository structure

- The repository root contains the MATLAB/Octave implementation and data used for the paper experiments.
- `cpp/` contains an experimental C++ implementation/prototype. It includes online graph/state expansion, backend comparison hooks, paper-style metric reporting, and algorithmic FLOP accounting. The C++ implementation is provided for development and reuse, but it has not yet been fully validated against all paper tables on all benchmark datasets.
- `python/` contains a Python prototype of the IGG/SPO algorithmic workflow. It is intended for readability, unit-tested experimentation, and algorithm development. It is not intended as a byte-for-byte reproduction of the MATLAB/Octave implementation or as the primary paper-reproduction path.

## Python prototype status

The Python implementation in `python/` is useful for inspecting and extending the main algorithmic ideas, including information-gain gating, selective active-set refinement, normalized chi-squared evaluation, and ATE-style trajectory assessment. It is provided as a readable research prototype and development aid.

Use the MATLAB/Octave implementation in the repository root for reproducing the published paper results. See `python/README.md` for Python installation instructions, testing commands, and implementation notes.

## C++ prototype status

The C++ implementation in `cpp/` is a research prototype rather than the primary reproduction artifact. It is useful for exploring compiled implementations of the proposed ideas, including sparse numerical backends and online graph expansion. However, it should not yet be interpreted as a fully validated replacement for the MATLAB/Octave pipeline.

See `cpp/README.md` for C++ build instructions, backend notes, and benchmark commands.

## Animated examples

- https://www.youtube.com/watch?v=J5eVcUjuUBw
- https://www.youtube.com/watch?v=59vPl3Zh7F8
- https://www.youtube.com/watch?v=fLFyXeX8Vbc
- https://www.youtube.com/watch?v=dFoqH3nJXW8
- https://www.youtube.com/watch?v=sFVPo_RIFBA
- https://www.youtube.com/watch?v=xeePpAjU87I

## Citation

If you use this code or data, please cite:

```bibtex
@article{Arablouei2026IncrementalSLAM,
  author  = {Arablouei, Reza},
  title   = {Efficient Incremental SLAM via Information-Guided Gating and Selective Partial Optimization},
  journal = {Robotics},
  year    = {2026},
  volume  = {15},
  number  = {5},
  article = {87},
  doi     = {10.3390/robotics15050087},
  url     = {https://www.mdpi.com/2218-6581/15/5/87}
}
```
