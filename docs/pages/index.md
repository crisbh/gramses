---
title: User guide
ordered_subpage: building-and-running.md
ordered_subpage: gr-solver.md
ordered_subpage: time-stepping.md
---

GRAMSES is a modified version of [RAMSES](https://github.com/ramses-organisation/ramses) that performs numerical relativity simulations in cosmology, using an N-body description for dark matter. It is based on a fully constrained, non-linear 3+1 formulation of general relativity, and extends the RAMSES adaptive mesh refinement (AMR) and multigrid infrastructure to solve the Einstein field equations while evolving particle geodesics.

Most of the code is upstream RAMSES. The GR-specific parts are a small set of `*_gr*` files in `poisson/` and `pm/`, plus GR branches in [[amr_step]]. These pages describe those parts; the API reference covers everything compiled into `bin/ramses3d`.

## Papers

- [GRAMSES Paper I](http://arxiv.org/abs/1905.08890): methodology and code description.
- [GRAMSES Paper II](http://arxiv.org/abs/2001.07968): gauge and generation of initial conditions.

Applications:

- [Gravitomagnetic effects in a shearing-dust universe](https://arxiv.org/abs/2003.08014)
- [Gravitomagnetic field in cosmology](https://arxiv.org/abs/2010.08257)
- [Gravitomagnetic distortion in light propagation](https://arxiv.org/abs/2109.02632)

## Citation

If you use GRAMSES in a paper, please cite Barrera-Hinojosa & Li, JCAP 01 (2020) 007, [doi:10.1088/1475-7516/2020/01/007](https://doi.org/10.1088/1475-7516/2020/01/007).
