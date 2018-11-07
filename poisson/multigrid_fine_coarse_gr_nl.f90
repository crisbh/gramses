! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-coarse-level related routines
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     GR pot fields        gr_pot(:,igr)  active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!     restricted GR pot    N/A            active_mg(myid,ilevel)%u(:,5)
!     restricted coeff     N/A            active_mg(myid,ilevel)%u(:,6)
!
! ------------------------------------------------------------------------

! --------------------------------------------------------------------------------
! Compute residual for MG levels, and stores it into active_mg(myid,ilevel)%u(:,3)
! --------------------------------------------------------------------------------

subroutine cmp_residual_mg_coarse_gr_nl(ilevel,igr)
   use amr_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none
   integer, intent(in) :: ilevel,igr
   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, dx2, phi_c, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: icell_mg, iskip_mg, igrid_nbor_mg, icell_nbor_mg
   integer  :: igrid_amr, iskip_amr, cpu_nbor_amr
   integer  :: igshift, igrid_nbor_amr

   real(dp) :: dtwondim = (twondim)
   real(dp) :: ctilde,twoac2,aomega
   real(dp) :: potc,gr_a,gr_b,op

   dx  = 0.5d0**ilevel

   ! Set constants
   dx2    = dx**2                     ! Cell size squared
   ctilde = sol/boxlen_ini/100000.0d0 ! Speed of light in code units
   twoac2 = 2.0D0*(aexp*ctilde)**2    ! 2a^2c^2 factor 
   aomega = 1.5D0*aexp*omega_m        ! Numerical coeff for S_0 in psi Eq.

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active_mg(myid,ilevel)%ngrid
   
   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in cmp_residual_mg_coarse_gr_nl. Please check.'
      call clean_stop
   end if

   ! Loop over cells myid
   do ind=1,twotondim
      iskip_mg  = (ind-1)*ngrid
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids myid
      do igrid_mg=1,ngrid
         igrid_amr = active_mg(myid,ilevel)%igrid(igrid_mg)
         icell_mg = igrid_mg + iskip_mg
        
         gr_b = active_mg(myid,ilevel)%u(icell_mg,6) ! Restricted coeff in coarse cell
         potc = active_mg(myid,ilevel)%u(icell_mg,1) ! Value of GR field on center cell
         nb_sum=0.0d0                                ! Sum of GR field on neighbors

         ! SCAN FLAG TEST
         if(.not. btest(active_mg(myid,ilevel)%f(icell_mg,1),0)) then ! NO SCAN
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if
                  igrid_nbor_mg = lookup_mg(igrid_nbor_amr)
                  icell_nbor_mg = igrid_nbor_mg + &
                      (jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid
                  nb_sum = nb_sum + &
                              active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
               end do
            end do
          
            if(igr==5) then 
               op = (nb_sum-6.0D0*potc)*(1.0D0-potc/twoac2) + &
                    dx2*(aomega*((1.0D0-potc/twoac2)**6-1.0D0)-0.25D0*gr_b/(1.0D0-potc/twoac2)**6)
            else
               op = nb_sum - 6.0D0*potc - dx2*potc*gr_b
            end if
         end if ! END SCAN TEST

         ! Store ***MINUS THE RESIDUAL***
         active_mg(myid,ilevel)%u(icell_mg,3) = -op/dx2 + active_mg(myid,ilevel)%u(icell_mg,2)
      end do
   end do

end subroutine cmp_residual_mg_coarse_gr_nl

subroutine cmp_uvar_norm2_coarse_gr_nl(ivar,ilevel,norm2,n_cell_c)
   use amr_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none

   integer,  intent(in)  :: ilevel,ivar
   real(dp), intent(out) :: norm2,n_cell_c

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg,icell_mg,iskip_mg

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active_mg(myid,ilevel)%ngrid

   norm2    = 0.0d0
   n_cell_c = 0.0d0
   do ind=1,twotondim
      iskip_mg = (ind-1)*ngrid
      do igrid_mg=1,ngrid
         icell_mg = iskip_mg+igrid_mg
         if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0 .and. ivar/=4) cycle
         norm2    = norm2+(active_mg(myid,ilevel)%u(icell_mg,ivar))**2.0
         n_cell_c = n_cell_c+1.0d0
      end do
   end do
!  norm2 = dx2*norm2
end subroutine cmp_uvar_norm2_coarse_gr_nl   

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_coarse_gr_nl(ilevel,safe,redstep,igr)
   use amr_commons
   use pm_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none
   integer, intent(in) :: ilevel, igr
   logical, intent(in) :: safe
   logical, intent(in) :: redstep

   integer, dimension(1:3,1:2,1:8) :: iii, jjj
   integer, dimension(1:3,1:4)     :: ired, iblack

   real(dp) :: dx2, nb_sum, weight
   integer  :: ngrid
   integer  :: ind, ind0, igrid_mg, idim, inbor
   integer  :: igrid_amr, cpu_nbor_amr
   integer  :: iskip_mg, igrid_nbor_mg, icell_mg, icell_nbor_mg
   integer  :: igshift, igrid_nbor_amr
   real(dp) :: dtwondim = (twondim)

   real(dp) :: ctilde,twoac2,aomega
   real(dp) :: potc,gr_a,gr_b,op,dop

   ! Set constants
   dx2    = (0.5d0**ilevel)**2        ! Cell size squared
   ctilde = sol/boxlen_ini/100000.0d0 ! Speed of light in code units
   twoac2 = 2.0D0*(aexp*ctilde)**2    ! 2a^2c^2 factor 
   aomega = 1.5D0*aexp*omega_m        ! Numerical coeff for S_0 in psi Eq.

   ired  (1,1:4)=(/1,0,0,0/)
   iblack(1,1:4)=(/2,0,0,0/)
   ired  (2,1:4)=(/1,4,0,0/)
   iblack(2,1:4)=(/2,3,0,0/)
   ired  (3,1:4)=(/1,4,6,7/)
   iblack(3,1:4)=(/2,3,5,8/)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active_mg(myid,ilevel)%ngrid
   
   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in gauss_seidel_mg_coarse_gr_nl. Please check.'
      call clean_stop
   end if

   ! Loop over cells, with red/black ordering
   do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep
      if(redstep) then
         ind = ired  (ndim,ind0)
      else
         ind = iblack(ndim,ind0)
      end if

      iskip_mg  = (ind-1)*ngrid

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active_mg(myid,ilevel)%igrid(igrid_mg)
         icell_mg  = iskip_mg  + igrid_mg

         ! Read scan flag
         if(.not. btest(active_mg(myid,ilevel)%f(icell_mg,1),0)) then
         
            gr_b = active_mg(myid,ilevel)%u(icell_mg,6) ! Restricted coeff in coarse cell
            potc = active_mg(myid,ilevel)%u(icell_mg,1) ! Value of GR field on center cell
            nb_sum=0.0d0                                ! Sum of GR field on neighbors

            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
            ! Those cells are active, have all their neighbors active
            ! and all neighbors are in the AMR+MG trees
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid shift
                  igshift = iii(idim,inbor,ind)
                  ! Get neighbor grid
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if
                  ! Get neighbor cpu
                  igrid_nbor_mg  = lookup_mg(igrid_nbor_amr)
                  icell_nbor_mg  = igrid_nbor_mg + &
                      (jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid
                  nb_sum = nb_sum + &
                              active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
               end do
            end do
            
            if(igr==5) then 
               op = (nb_sum-6.0D0*potc)*(1.0D0-potc/twoac2) + &
                    dx2*(aomega*((1.0D0-potc/twoac2)**6-1.0D0)-0.25D0*gr_b/(1.0D0-potc/twoac2)**6)
               dop= -nb_sum/twoac2 - 6.0D0 + 12.0D0*potc/twoac2 &
                    -dx2*(6.0D0*aomega*(1.0D0-potc/twoac2)**5/twoac2 + 1.5D0*gr_b/(1.0D0-potc/twoac2)**7/twoac2)
            else
               op = nb_sum - 6.0D0*potc - dx2*potc*gr_b
               dop= -6.0D0 - dx2*gr_b 
            end if
            
            op = op - active_mg(myid,ilevel)%u(icell_mg,2)*dx2
            ! Update the GR potential, solving for potential on icell_amr
            active_mg(myid,ilevel)%u(icell_mg,1)=active_mg(myid,ilevel)%u(icell_mg,1) -&
                                                 op/dop
         end if
      end do
   end do
end subroutine gauss_seidel_mg_coarse_gr_nl

! ------------------------------------------------------------------------
! Compute the physical rhs of the coarse level, including the corrections
! from the finer level, such as the residual etc.
! The physical rhs is stored in active_mg(myid,ilevel)%u(:,2) 
! ------------------------------------------------------------------------

subroutine make_physical_rhs_coarse_gr_nl(ilevel,igr)
   use amr_commons
   use pm_commons
   use poisson_commons
   use gr_commons
   use gr_parameters
   
   implicit none

   integer, intent(in) :: ilevel, igr
   integer, dimension(1:3,1:2,1:8) :: iii,jjj
   
   integer  :: ngrid
   integer  :: ind
   integer  :: igshift,inbor,idim
   integer  :: iskip_mg,igrid_mg,icell_mg
   integer  :: iskip_amr,igrid_amr,icell_amr
   integer  :: icell_nbor_mg,igrid_nbor_mg
   integer  :: icell_nbor_amr,igrid_nbor_amr
   integer  :: cpu_nbor_amr

   real(dp) :: dx,dx2,nb_sum
   real(dp) :: dtwondim = (twondim)
   real(dp) :: ctilde,twoac2,aomega
   real(dp) :: potc,gr_a,gr_b,op,dop

   ! Set constants
   dx2    = (0.5d0**ilevel)**2        ! Cell size squared
   ctilde = sol/boxlen_ini/100000.0d0 ! Speed of light in code units
   twoac2 = 2.0D0*(aexp*ctilde)**2    ! 2a^2c^2 factor 
   aomega = 1.5D0*aexp*omega_m        ! Numerical coeff for S_0 in psi Eq.

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active_mg(myid,ilevel)%ngrid

   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in make_physical_rhs_coarse_gr_nl. Please check.'
      call clean_stop
   end if

   ! Loop over cells
   do ind=1,twotondim
      iskip_mg  = (ind-1)*ngrid
      iskip_amr = (ind-1)*ngridmax+ncoarse

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active_mg(myid,ilevel)%igrid(igrid_mg)     ! amr grid index
         icell_amr = iskip_amr+igrid_amr                        ! amr cell index
         icell_mg  = iskip_mg+igrid_mg                          ! mg cell index

         if(.not. btest(active_mg(myid,ilevel)%f(icell_mg,1),0)) then
            gr_b = active_mg(myid,ilevel)%u(icell_mg,6)  ! Restricted coeff in coarse cell
            potc = active_mg(myid,ilevel)%u(icell_mg,5)  ! Value of GR field on center cell
            nb_sum=0.0d0                                 ! Sum of GR field on neighbors

            do inbor=1,2
               do idim=1,ndim
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if
                  igrid_nbor_mg = lookup_mg(igrid_nbor_amr)
                  icell_nbor_mg = igrid_nbor_mg+(jjj(idim,inbor,ind)-1)*active_mg(cpu_nbor_amr,ilevel)%ngrid
                  nb_sum = nb_sum + active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,5)
               end do
            end do

            if(igr==5) then 
               op = (nb_sum-6.0D0*potc)*(1.0D0-potc/twoac2) + &
                    dx2*(aomega*((1.0D0-potc/twoac2)**6-1.0D0)-0.25D0*gr_b/(1.0D0-potc/twoac2)**6)
            else
               op = nb_sum - 6.0D0*potc - dx2*potc*gr_b
            end if
         end if
         
         active_mg(myid,ilevel)%u(icell_mg,2) = active_mg(myid,ilevel)%u(icell_mg,2)+op/dx2
      end do
   end do

end subroutine make_physical_rhs_coarse_gr_nl

!-------------------------------------------------------------------
! Restriction of the gr_pot
!------------------------------------------------------------------

subroutine restrict_gr_pot_coarse_reverse_gr_nl(ifinelevel)

   use amr_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: sf1
   real(dp) :: dtwotondim = (twotondim)
   integer,dimension(:),allocatable::n_masked

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   allocate(n_masked(1:active_mg(myid,ifinelevel)%ngrid))
   n_masked(1:active_mg(myid,ifinelevel)%ngrid)=0

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid
      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         ! Count the number of masked fine cells in each fine grid
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) then
            n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
         end if
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg                                  ! mg fine cell index

         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)          ! amr fine grid index
         icell_c_amr=father(igrid_f_amr)                                   ! amr coarse cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1                     ! amr coarse cell position
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1                     ! amr coarse cell position
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax           ! amr coarse grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                              ! coarse cell cpu index

         ! Convert to MG index, get MG coarse cell id                      ! mg coarse grid index
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                                  ! mg coarse cell index

         ! If coarse cell is masked, it is boundary and R\tilde{u} is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0.or.n_masked(igrid_f_mg)==8) cycle

         ! Restriction to compute the sf value in the coarse cell 
         sf1=active_mg(myid,ifinelevel)%u(icell_f_mg,1)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)+sf1
      end do
   end do

   deallocate(n_masked)   

end subroutine restrict_gr_pot_coarse_reverse_gr_nl

subroutine restrict_coeff_coarse_reverse_gr_nl(ifinelevel)
   use amr_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: rho1
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   integer,dimension(:),allocatable::n_masked

   icoarselevel=ifinelevel-1

   allocate(n_masked(1:active_mg(myid,ifinelevel)%ngrid))
   n_masked(1:active_mg(myid,ifinelevel)%ngrid)=0

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid
      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         ! Count the number of masked fine cells in each fine grid
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) then
            n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
         end if
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg                                  ! mg fine cell index
         ! Is fine cell masked?
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)          ! amr fine grid index
         icell_c_amr=father(igrid_f_amr)                                   ! amr coarse cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1                     ! amr coarse cell position
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1                     ! amr coarse cell position
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax           ! amr coarse grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                              ! coarse cell cpu index

         ! Convert to MG index, get MG coarse cell id                      ! mg coarse grid index
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                                  ! mg coarse cell index

         ! If coarse cell is masked, it's outside boundary and R\rho is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0.or.n_masked(igrid_f_mg)==8) cycle

         ! Restriction to compute the sf value in the coarse cell 
         rho1=active_mg(myid,ifinelevel)%u(icell_f_mg,6)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)+rho1
      end do
   end do

   deallocate(n_masked)

end subroutine restrict_coeff_coarse_reverse_gr_nl


! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_coarse_gr_nl(ifinelevel,igr)
   use amr_commons
   use poisson_commons
   use gr_commons   
   implicit none
   integer, intent(in) :: ifinelevel, igr

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr, ngrid_f, istart, nbatch
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg, iskip_f_mg, icell_f_mg
   integer  :: icoarselevel, ind_c, cpu_c_amr

   real(dp) :: a, b, c, d, coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector), save                :: igrid_f_amr, icell_amr, cpu_amr
   integer,  dimension(1:nvector,1:threetondim), save  :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim), save    :: nbors_father_grids
   real(dp), dimension(1:nvector), save                :: corr

   ! Local constants
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.D0*a
   icoarselevel=ifinelevel-1

   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Loop over fine grids by vector sweeps
   ngrid_f=active_mg(myid,ifinelevel)%ngrid
   do istart=1,ngrid_f,nvector

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active_mg(myid,ifinelevel)%igrid(istart+i-1)
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))
         cpu_amr(i)  =cpu_map(icell_amr(i))
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather(icell_amr,nbors_father_cells,nbors_father_grids,nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax
         iskip_f_mg  = (ind_f-1)*ngrid_f

         do i=1,nbatch
            ! Compute fine cell indices
            icell_amr(i) = iskip_f_amr + igrid_f_amr(i)
         end do
         corr=0.0d0

         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            do i=1,nbatch
               icell_f_mg  = iskip_f_mg + istart+i-1
               if(active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,4)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
               ind_c       = (icell_c_amr-ncoarse)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               cpu_c_amr   = cpu_map(father(igrid_c_amr))
               if(igrid_c_mg<=0) cycle

               icell_c_mg  = (ind_c-1)*active_mg(cpu_c_amr,icoarselevel)%ngrid + igrid_c_mg
               ! only unmasked coarse cells contribute to the fine-cell correction
               if(active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,4)<=0.0) cycle
               corr(i)=corr(i)+coeff*active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,1) &
                      &       -coeff*active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,5)
            end do
         end do

         ! Correct potential
         do i=1,nbatch
            icell_f_mg  = iskip_f_mg + istart+i-1
            active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1) = active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1) + corr(i)
         end do

      end do
      ! End loop over cells

   end do
   ! End loop over grids
end subroutine interpolate_and_correct_coarse_gr_nl
