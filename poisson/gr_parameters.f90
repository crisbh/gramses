module gr_parameters
    use amr_parameters

    ! convergence criterion
    real(dp) :: epsilon_gr=1.0D-5

    ! c in standard units
    real(dp) :: sol=299792458.0D0

    ! Relaxation parameters for linear/non-linear GR equations
    integer :: ngs_fine_gr_ln = 5   ! Linear 
    integer :: ngs_fine_gr_nl = 5   ! Non-linear 


end module gr_parameters
