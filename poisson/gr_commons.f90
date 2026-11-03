module gr_commons
   !! Cell-based storage for the GR fields and their matter sources.
   !!
   !! `gr_pot(:,igr)`, `igr`=1..10: \(V_i\) (1-3), \(U\) (4), \(\alpha\) solved as \(\alpha\psi\) (5),
   !! \(\psi\) (6), \(B^i\) (7-9), \(b\) (10).
   !!
   !! `gr_mat(:,1:4)`: matter sources \(s_i\) and \(s\); \(s_0\) is stored in `rho`.
   !! The same slots also hold \(\nabla\cdot V\), \(\nabla\cdot B\) and \(A_{ij}A^{ij}\) during the solve.
   use amr_parameters

! gr_pot is used to store the 10 GR fields (V_i,U,alpha,psi,B^i,b)
! gr_mat is used to store the 4 matter source terms (si,s) as well as (div(V),div(B),A_ijA^ij). 
! Notice that s0 is stored in the original rho array.

   real(dp),allocatable,dimension(:,:) :: gr_pot  !! GR fields, indexed by `igr`
   real(dp),allocatable,dimension(:,:) :: gr_mat  !! GR matter source terms
   real(dp),allocatable,dimension(:,:) :: gr_mat2 !! GR source terms for output

end module gr_commons
