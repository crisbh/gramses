module gr_commons
   use amr_parameters

! gr_pot is used to store the 10 GR fields (V_i,U,alpha,psi,B^i,b)
! gr_mat is used to store the 4 matter source terms (si,s) as well as (div(V),div(B),A_ijA^ij). 
! Notice that s0 is stored in the original rho array.

   real(dp),allocatable,dimension(:,:) :: gr_pot  ! GR fields 
   real(dp),allocatable,dimension(:,:) :: gr_mat  ! GR source terms
   real(dp),allocatable,dimension(:,:) :: gr_mat2 !GR source terms for output

end module gr_commons
