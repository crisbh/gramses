module gr_commons
   use amr_parameters

! gr_pot is used to store the 10 geometric variables (Vi,U,alpha,psi,Bi,b)
! gr_mat is used to store the 4 matter source terms (si,s) as well as (div(V),div(B),A_ijA^ij). 
! Notice that s0 is stored in the original rho array.

   real(dp),allocatable,dimension(:,:) :: gr_pot  ! geometric fields 
   real(dp),allocatable,dimension(:,:) :: gr_mat  ! source terms

end module gr_commons
