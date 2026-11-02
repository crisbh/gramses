---
title: GR field solver
---

## Fields and sources

The GR fields and their matter sources are stored per cell in [[gr_commons]].

`gr_pot(:,igr)` holds the 10 GR fields:

| `igr` | Field                         | Equation type |
|-------|-------------------------------|---------------|
| 1–3   | $V_i$                         | linear        |
| 4     | $U$                           | linear        |
| 5     | $\alpha$ (solved as $\alpha\psi$) | non-linear |
| 6     | $\psi$                        | non-linear    |
| 7–9   | $B^i$                         | linear        |
| 10    | $b$                           | linear        |

`gr_mat(:,1:4)` holds the matter source terms $s_i$ and $s$; $s_0$ is stored in the RAMSES `rho` array. The same slots are reused for intermediate sources such as $\nabla\cdot V$, $\nabla\cdot B$ and $A_{ij}A^{ij}$.

## Multigrid driver

[[multigrid_fine_gr]] solves one field `igr` on one AMR level. Before relaxing, it builds the right-hand side that depends on previously solved fields:

| `igr` | Extra source computed                                       |
|-------|-------------------------------------------------------------|
| 1     | Fine-level mask (first call per level)                      |
| 4     | $\nabla\cdot V$ ([[source_fine_gr_scalar]])                 |
| 5     | $A_{ij}A^{ij}$ ([[source_fine_gr_aij_aij]])                 |
| 7     | $\partial_j A^{ij}$ vector sources ([[comp_gr_aij]], [[source_fine_gr_vector]]) |
| 10    | $\nabla\cdot B$ ([[source_fine_gr_scalar]])                 |

Fields 5 and 6 use the non-linear relaxation in `multigrid_fine_fine_gr_nl.f90` and `multigrid_fine_coarse_gr_nl.f90`; the rest use the linear relaxation in `multigrid_fine_fine_gr_ln.f90`. Each relaxation sweep is a red–black Gauss–Seidel step with ghost-zone exchange in between, and the non-linear solver subtracts the mean of the source for regularisation.

Storage used during the solve, at the finest AMR level and on the coarse multigrid levels:

| Quantity                 | Finest AMR level | Coarse MG levels                |
|--------------------------|------------------|---------------------------------|
| GR field                 | `gr_pot(:,igr)`  | `active_mg(myid,ilevel)%u(:,1)` |
| Matter RHS               | `gr_mat(:,igrm)` | `active_mg(myid,ilevel)%u(:,2)` |
| Residual                 | `f(:,1)`         | `active_mg(myid,ilevel)%u(:,3)` |
| BC-modified RHS          | `f(:,2)`         | n/a                             |
| Mask                     | `f(:,3)`         | `active_mg(myid,ilevel)%u(:,4)` |
| Restricted NL GR field   | n/a              | `active_mg(myid,ilevel)%u(:,5)` |
| Restricted density       | n/a              | `active_mg(myid,ilevel)%u(:,6)` |

The equations themselves are derived in [Paper I](http://arxiv.org/abs/1905.08890).
