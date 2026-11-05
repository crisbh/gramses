module poisson_parameters
  !! Newtonian Poisson solver parameters.
  use amr_parameters

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4  !! Convergence criterion for the Poisson solvers

  ! Type of force computation
  integer ::gravity_type=0  !! Type of force computation

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0.0  !! Gravity parameters

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0  !! Maximum level for CIC dark matter interpolation

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999  !! Levels below use fine multigrid, levels at or above use conjugate gradient

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2  !! Gauss-Seidel sweeps on fine multigrid levels
  integer, parameter :: ngs_coarse = 2  !! Gauss-Seidel sweeps on coarse multigrid levels

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1  !! Multigrid cycles on coarse levels in safe mode

end module poisson_parameters
