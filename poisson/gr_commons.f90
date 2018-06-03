module gr_commons
    use amr_parameters

! we need to declare 3 scalars (alpha, psi, s), a 4-vector for (W, beta) and a 4-vector for Si.

    real(dp),allocatable,dimension(:)   :: gr_alpha  ! lapse function
    real(dp),allocatable,dimension(:)   :: gr_psi    ! conformal factor
    real(dp),allocatable,dimension(:,:) :: gr_wbeta  ! (U,dV^i/dx^i) & (beta, dB^i/dx)
    real(dp),allocatable,dimension(:,:) :: gr_si     ! source term s_i
    real(dp),allocatable,dimension(:)   :: gr_strace ! trace of s_{ij}

end module gr_commons
