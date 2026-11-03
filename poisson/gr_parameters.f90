module gr_parameters
   !! Tuning parameters of the GR multigrid solver.
   use amr_parameters

   ! Source mean value for solver regularisation
   !> Mean source of the linear equations, for regularisation
   real(dp) :: src_mean_ln =0.0d0
   !> Mean source of the non-linear equations at the finest level
   real(dp) :: src_mean_nl =0.0d0
   !> Mean source of the non-linear equations on coarse multigrid levels
   real(dp) :: src_mean_nl_coarse =0.0d0 

   !> Multigrid convergence criterion
   real(dp) :: epsilon_gr=1.0D-4

   !> Speed of light in m/s
   real(dp) :: sol=299792458.0D0

   ! Relaxation parameters for linear/non-linear GR equations
   integer :: ngs_fine_gr_ln_pre = 2   !! Linear pre-smoothing sweeps (fine levels)
   integer :: ngs_fine_gr_ln_pst = 2   !! Linear post-smoothing sweeps (fine levels)
   integer :: ngs_fine_gr_nl_pre = 2   !! Non-linear pre-smoothing sweeps (fine levels)
   integer :: ngs_fine_gr_nl_pst = 2   !! Non-linear post-smoothing sweeps (fine levels)
   
   integer :: ngs_coarse_gr_ln_pre = 2   !! Linear pre-smoothing sweeps (coarse levels)
   integer :: ngs_coarse_gr_ln_pst = 2   !! Linear post-smoothing sweeps (coarse levels)
   integer :: ngs_coarse_gr_nl_pre = 2   !! Non-linear pre-smoothing sweeps (coarse levels)
   integer :: ngs_coarse_gr_nl_pst = 2   !! Non-linear post-smoothing sweeps (coarse levels)

end module gr_parameters
