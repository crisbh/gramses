module gr_parameters
   use amr_parameters

   ! Source mean value for solver regularisation
   real(dp) :: src_mean_ln =0.0d0
   real(dp) :: src_mean_nl =0.0d0
   real(dp) :: src_mean_nl_coarse =0.0d0 

   ! convergence criterion
   real(dp) :: epsilon_gr=1.0D-4

   ! c in standard units
   real(dp) :: sol=299792458.0D0

   ! Relaxation parameters for linear/non-linear GR equations
   integer :: ngs_fine_gr_ln_pre = 2   ! Linear pre-smoothing factor
   integer :: ngs_fine_gr_ln_pst = 2   ! Linear post-smoothing factor 
   integer :: ngs_fine_gr_nl_pre = 2   ! Non-linear pre-smoothing factor 
   integer :: ngs_fine_gr_nl_pst = 2   ! Non-linear post-smoothing factor 
   
   integer :: ngs_coarse_gr_ln_pre = 2   ! Linear pre-smoothing factor
   integer :: ngs_coarse_gr_ln_pst = 2   ! Linear post-smoothing factor 
   integer :: ngs_coarse_gr_nl_pre = 2   ! Non-linear pre-smoothing factor 
   integer :: ngs_coarse_gr_nl_pst = 2   ! Non-linear post-smoothing factor 

end module gr_parameters
