module gr_parameters
    use amr_parameters

    ! convergence criterion
    real(dp) :: epsilon_gr=1.0D-5

    ! c in standard units
    real(dp) :: sol=299792458.0D0

end module gr_parameters
