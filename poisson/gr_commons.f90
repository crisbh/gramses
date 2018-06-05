module gr_commons
    use amr_parameters

! gr_pot holds the 10 geometric variables and gr_mat matter source terms.

    real(dp),allocatable,dimension(:,:) :: gr_pot  ! geometric fields 
    real(dp),allocatable,dimension(:,:) :: gr_mat  ! source terms

end module gr_commons
