module gr_parameters
    use amr_parameters

    ! convergence criterion
    real(dp) :: epsilon_gr=1.0D-5

    ! c in standard units
    real(dp) :: sol=299792458.0D0

    ! Relaxation parameters for linear/non-linear GR equations
    integer :: ngs_fine_gr_ln_pre = 5   ! Linear pre-smoothing factor
    integer :: ngs_fine_gr_ln_pst = 5   ! Linear post-smoothing factor 
    integer :: ngs_fine_gr_nl_pre = 5   ! Non-linear pre-smoothing factor 
    integer :: ngs_fine_gr_nl_pst = 5   ! Non-linear post-smoothing factor 
    
    integer :: ngs_coarse_gr_ln_pre = 5   ! Linear pre-smoothing factor
    integer :: ngs_coarse_gr_ln_pst = 5   ! Linear post-smoothing factor 
    integer :: ngs_coarse_gr_nl_pre = 5   ! Non-linear pre-smoothing factor 
    integer :: ngs_coarse_gr_nl_pst = 5   ! Non-linear post-smoothing factor 

end module gr_parameters
