# GRAMSES — Numerical Relativity N-body Simulations in Cosmology

`GRAMSES` is a modified version of [`RAMSES`](https://github.com/ramses-organisation/ramses) to perform **numerical relativity simulations** in cosmology using a particle (N-body) description for dark matter. 
It is based on a fully constrained and non-linear 3+1 formulation of general relativity. For details of the mathematical formalism and code implementation, see:  
- [GRAMSES Paper 1](http://arxiv.org/abs/1905.08890) – methodology and code description  
- [GRAMSES Paper 2](http://arxiv.org/abs/2001.07968) – gauge and generation of initial conditions


## Overview

`GRAMSES` extends `RAMSES`' adaptive mesh refinement (AMR) and multigrid infrastructure to solve **Einstein's field equations** in a cosmological context, while evolving the geodesic equations for matter in an N-body particle scheme.

**Applications**

- [Gravitomagnetic effects in a shearing-dust universe](https://arxiv.org/abs/2003.08014) - comparison with other numerical relativity codes
- [Gravitomagnetic field in cosmology](https://arxiv.org/abs/2010.08257) - vector modes in cosmology
- [Gravitomagnetic distortion in light propagation](https://arxiv.org/abs/2109.02632) - weak-lensing x kSZ effect signal


## Key Features

- **Full AMR hierarchy** for efficient simulations of multi-scale dynamics.
- **Multigrid solvers** for the 3+1 Einstein field equations.
- **Kick-Drift-Kick** time integration for particle geodesics.
- **MPI (Message Passing Interface) parallelisation** for distributed-memory HPC systems.


## Installation and usage
```bash
git clone https://github.com/crisbh/gramses.git
cd gramses
make -j
```

Then, you can run as
```bash
mpirun -np 8 ./ramses3d namelist/relativistic_cosmo.nml
```

---

## Citation

If you use `GRAMSES` as part of a paper, please cite our publication:

```
@ARTICLE{Barrera_Li:2020a,
       author = {{Barrera-Hinojosa}, Cristian and {Li}, Baojiu},
        title = "{GRAMSES: a new route to general relativistic N-body simulations in cosmology. Part I. Methodology and code description}",
      journal = {\jcap},
     keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics, General Relativity and Quantum Cosmology},
         year = 2020,
        month = jan,
       volume = {2020},
       number = {1},
          eid = {007},
        pages = {007},
          doi = {10.1088/1475-7516/2020/01/007}
}
```


