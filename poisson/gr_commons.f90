module gr_commons
   use amr_parameters

! gr_pot is used to store the 10 geometric variables (Vi,U,alpha,psi,Bi,b) and gr_mat the matter source terms (si,s). 
! Notice that s0 is stored in the original rho array.

   real(dp),allocatable,dimension(:,:) :: gr_pot  ! geometric fields 
   real(dp),allocatable,dimension(:,:) :: gr_mat  ! source terms

end module gr_commons
