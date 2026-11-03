---
title: Building and running
---

## Build

```bash
cd bin
make -j        # produces bin/ramses3d
make clean     # removes *.o and *.mod
```

Compile-time options are variables at the top of `bin/Makefile` and can be overridden on the command line:

| Variable  | Default | Meaning                                          |
|-----------|---------|--------------------------------------------------|
| `NDIM`    | 3       | Number of dimensions                             |
| `NVECTOR` | 64      | Vector sweep size                                |
| `NPRE`    | 8       | Real precision in bytes                          |
| `NVAR`    | 5       | Number of hydro variables                        |
| `NENER`   | 0       | Number of non-thermal energies                   |
| `SOLVER`  | hydro   | Hydro solver directory                           |
| `PATCH`   | (empty) | Patch directory searched before the core sources |

The compiler line is selected by uncommenting one of the blocks in `bin/Makefile`. The default is `mpif90` with gfortran flags; there is also a gfortran debug block with bounds checking and floating-point traps.

Source files are found through `VPATH`, searched in order: patch directories, `../$(SOLVER)`, `../aton`, `../hydro`, `../pm`, `../poisson`, `../amr`. A file in an earlier directory shadows one with the same name in a later directory. New source files must be added to the object lists in `bin/Makefile`.

## Run

```bash
mpirun -np 8 ./bin/ramses3d namelist/gr_sim.nml
```

`namelist/gr_sim.nml` is a cosmological GR run from Gadget initial conditions; `namelist/gr_test.nml` is a smaller test setup.

## GR switches in `&RUN_PARAMS`

| Parameter      | Default   | Effect                                                                                  |
|----------------|-----------|-----------------------------------------------------------------------------------------|
| `gr`           | `.false.` | Solve the GR field equations and move particles along geodesics                        |
| `gr_newtonian` | `.false.` | With `gr`, compute the GR matter sources using a temporary Newtonian synchronisation of the particles |

Both are declared in [[amr_parameters]] and read in [[read_params]]. Solver tuning (convergence criterion, Gauss–Seidel smoothing counts) is in [[gr_parameters]].
