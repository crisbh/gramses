---
title: GR time stepping
---

The GR branch of [[amr_step]] runs on each AMR level when `gr=.true.`:

1. [[rho_fine]] deposits the GR matter sources into `rho` and `gr_mat` (flag `gr2`).
2. [[multigrid_fine_gr]] is called for `igr=1,…,10` in order. On the initial step (`nstep_coarse==-1`) only fields 1–4 are solved, particles are synchronised with field 4 only, and the step returns; this requires `levelmin` to be the only level in the initial conditions.
3. Particles are kicked by [[force_fine_gr]] and [[synchro_fine_gr]], one field at a time, in the order given by `gr_ord1 = (9,8,7,5,6)` (smallest to largest contribution).
4. Later in the step, [[move_fine_gr]] is called for each field in `gr_ord2 = (6,5,7,8,9,10)` (largest to smallest); `igrp==10` performs the position update.

GR fields (`gr_pot`) and sources (`gr_mat`) are included in the ghost-zone exchanges done by `make_virtual_fine_dp` and are carried through refinement and load balancing.

With `gr_newtonian=.true.`, step 1 changes: the Newtonian potential is solved and particles are synchronised with it ([[force_fine]], [[synchro_fine]]), the GR matter sources are deposited with these Newtonian-synchronised velocities, and the Newtonian synchronisation is then reverted. Steps 2–4 are unchanged.
