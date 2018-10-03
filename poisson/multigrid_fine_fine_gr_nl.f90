! ------------------------------------------------------------------------
! Multigrid GR fields solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-fine-level related routines
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     GR field             gr_pot(:,igr)  active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_fine_gr_nl(ilevel,igr)
   ! Computes the residual the fine (AMR) level, and stores it into f(:,1)
   use amr_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none
   integer, intent(in) :: ilevel,igr

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2, dx2, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)

   real(dp) :: ctilde,2ac2,aomega
   real(dp) :: potc,gr_a,gr_b,op
   
   ! Set constants
   dx2     = (0.5d0**ilevel)**2        ! Cell size squared
   ctilde  = sol/boxlen_ini/100000.0d0 ! Speed of light in code units
   2ac2    = 2.0D0*(aexp*ctilde)**2    ! 2a^2c^2 factor 
   aomega  = 1.5D0*aexp*omega_m        ! Numerical coeff for S_0 in psi Eq.

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid

   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in cmp_residual_mg_fine_gr_nl. Please check.'          
      call clean_stop
   end if

   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)     ! amr grid index
         icell_amr = iskip_amr + igrid_amr              ! amr cell index

         ! Calculate space-dependent coefficients
         if(igr==5) then 
            gr_a = aomega*rho   (icell_amr  )
            gr_b = 0.25D0*gr_mat(icell_amr,1)
         else
            gr_a = aomega*(rho(icell_amr)+      gr_mat(icell_amr,4)-      (1.0D0-gr_pot(icell_amr,5)/2ac2)**6) + &
                          gr_mat(icell_amr,1)/(1.0D0-gr_pot(icell_amr,5)/2ac2)**6
            gr_b = aomega*(rho(icell_amr)+2.0D0*gr_mat(icell_amr,4)+5.0D0*(1.0D0-gr_pot(icell_amr,5)/2ac2)**6) + &
                   1.75D0*gr_mat(icell_amr,1)/(1.0D0-gr_pot(icell_amr,5)/2ac2)**6
            gr_a = gr_a*dx2/(1.0D0-gr_pot(icell_amr,5)/2ac2)
            gr_b = gr_b*dx2/(1.0D0-gr_pot(icell_amr,5)/2ac2)**2/2ac2
         end if

         ! SCAN FLAG TEST
         if(flag2(icell_amr)/ngridmax==0) then
            potc  = gr_pot(icell_amr,igr)  ! Value of GR field on center cell
            nb_sum=0.0d0                   ! Sum of GR field on neighbors

            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if
                  icell_nbor_amr = igrid_nbor_amr + &
                      (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  nb_sum = nb_sum + gr_pot(icell_nbor_amr,igr)
               end do
            end do
         end if ! END SCAN TEST

         if(igr==5) then 
            op = (nb_sum-6.0D0*potc)*(1.0D0-potc/2ac2) - dx2*(gr_a-aomega*(1.0D0-potc/2ac2)**6+gr_b/(1.0D0-potc/2ac2)**6)
         else
            op =  nb_sum-6.0D0*potc - potc*gr_b - gr_a
         end if
         f(icell_amr,1) = -op*oneoverdx2
      end do
   end do

end subroutine cmp_residual_mg_fine_gr_nl

! ##################################################################
! ##################################################################

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_fine_gr_nl(ilevel,redstep,igr)
   use amr_commons
   use pm_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none
   integer, intent(in) :: ilevel,igr
   logical, intent(in) :: redstep

   integer, dimension(1:3,1:2,1:8) :: iii, jjj
   integer, dimension(1:3,1:4)     :: ired, iblack

   real(dp) :: dx2, nb_sum, weight
   integer  :: ngrid
   integer  :: ind, ind0, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)
   real(dp) :: ctilde,2ac2,aomega
   real(dp) :: potc,gr_a,gr_b,op,dop

   ! Set constants
   dx2     = (0.5d0**ilevel)**2        ! Cell size squared
   ctilde  = sol/boxlen_ini/100000.0d0 ! Speed of light in code units
   2ac2    = 2.0D0*(aexp*ctilde)**2    ! 2a^2c^2 factor 
   aomega  = 1.5D0*aexp*omega_m        ! Numerical coeff for S_0 in psi Eq.

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

   ngrid=active(ilevel)%ngrid
   
   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in gauss_seidel_mg_fine_gr_nl. Please check.'  
      call clean_stop
   end if

   ! Loop over cells, with red/black ordering
   do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep
      if(redstep) then
         ind = ired  (ndim,ind0)
      else
         ind = iblack(ndim,ind0)
      end if

      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         
         ! Calculate space-dependent coefficients
         if(igr==5) then 
            gr_a = aomega*rho   (icell_amr  )
            gr_b = 0.25D0*gr_mat(icell_amr,1)
         else
            gr_a = aomega*(rho(icell_amr)+      gr_mat(icell_amr,4)-      (1.0D0-gr_pot(icell_amr,5)/2ac2)**6) + &
                          gr_mat(icell_amr,1)/(1.0D0-gr_pot(icell_amr,5)/2ac2)**6
            gr_b = aomega*(rho(icell_amr)+2.0D0*gr_mat(icell_amr,4)+5.0D0*(1.0D0-gr_pot(icell_amr,5)/2ac2)**6) + &
                   1.75D0*gr_mat(icell_amr,1)/(1.0D0-gr_pot(icell_amr,5)/2ac2)**6
            gr_a = gr_a*dx2/(1.0D0-gr_pot(icell_amr,5)/2ac2)
            gr_b = gr_b*dx2/(1.0D0-gr_pot(icell_amr,5)/2ac2)**2/2ac2
         end if

         ! Read scan flag
         if(flag2(icell_amr)/ngridmax==0) then
            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
            ! Those cells are active, have all their neighbors active
            ! and all neighbors are in the AMR+MG trees
            potc = gr_pot(icell_amr,igr)    ! Value of GR field on center cell
            nb_sum=0.0d0                    ! Sum of GR field on neighbors

            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid shift
                  igshift = iii(idim,inbor,ind)
                  ! Get neighbor grid
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if
                  icell_nbor_amr = igrid_nbor_amr + &
                      (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  nb_sum = nb_sum + gr_pot(icell_nbor_amr,igr)
               end do
            end do

            if(igr==5) then 
               op = (nb_sum-6.0D0*potc)*(1.0D0-potc/2ac2) - dx2*(gr_a-aomega*(1.0D0-potc/2ac2)**6 + gr_b/(1.0D0-potc/2ac2)**6)
               dop= -nb_sum/2ac2 - 6.0D0 + 12.0D0*potc/2ac2 &
                    -dx2*(6.0D0*aomega*(1.0D0-potc/2ac2)**5/2ac2 + 6.0D0*gr_b/(1.0D0-potc/2ac2)**7/2ac2)
            else
               op =  nb_sum - 6.0D0*potc - potc*gr_b - gr_a
               dop= -6.0D0 - gr_b 
            end if

            ! Update the GR potential, solving for potential on icell_amr
            gr_pot(icell_amr,igr) = gr_pot(icell_amr,igr) - op/dop
         end if
      end do
   end do
end subroutine gauss_seidel_mg_fine_gr_nl

! ------------------------------------------------------------------------
! GR pot restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_gr_pot_fine_reverse_gr_nl(ifinelevel,igr)
   use amr_commons
   use poisson_commons
   use gr_commons
   implicit none
   integer, intent(in) :: ifinelevel,igr

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: sf1
   real(dp) :: dtwotondim = (twotondim)
   integer,dimension(:),allocatable::n_masked

   integer :: icoarselevel
   icoarselevel=ifinelevel-1
   
   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in restrict_gr_pot_fine_reverse_gr_nl. Please check.'  
      call clean_stop
   end if

   allocate(n_masked(1:active(ifinelevel)%ngrid))
   n_masked(1:active(ifinelevel)%ngrid)=0
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)               ! amr fine-grid index
         icell_f_amr=igrid_f_amr+iskip_f_amr                            ! amr fine-cell index
         if(f(icell_f_amr,3)<=0d0) n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)               ! amr fine-grid index
         icell_f_amr=igrid_f_amr+iskip_f_amr                            ! amr fine-cell index
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)                                ! amr coarse-cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax        ! amr coarse-grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                           ! cpu for coarse cell

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)                              ! mg coarse-grid index
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                               ! mg coarse-cell index

         ! If coarse cell masked, it is boundary and R\tilde{u} is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Restriction to compute the sf value in the coarse cell
         sf1=gr_pot(icell_f_amr,igr)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)+sf1
      end do
   end do

   deallocate(n_masked)

end subroutine restrict_gr_pot_fine_reverse_gr_nl

! ------------------------------------------------------------------------
! Restrict Coefficients of non-linear differential operators
! ------------------------------------------------------------------------

subroutine restrict_coeff_fine_reverse_gr_nl(ifinelevel,igr)

   use amr_commons
   use poisson_commons
   use gr_commons
   use gr_parameters

   implicit none
   integer, intent(in) :: ifinelevel,igr

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: rho1
   real(dp) :: dtwotondim = (twotondim)
   real(dp) :: ctilde,2ac2,aomega
   real(dp) :: gr_b

   integer  :: icoarselevel
   integer  :: ix,iy,iz
   real(dp) :: ctilde,xc,dx
   integer,dimension(:),allocatable::n_masked

   icoarselevel=ifinelevel-1
   
   ! Set constants
   dx2     = (0.5d0**ilevel)**2        ! Cell size squared
   ctilde  = sol/boxlen_ini/100000.0d0 ! Speed of light in code units
   2ac2    = 2.0D0*(aexp*ctilde)**2    ! 2a^2c^2 factor 
   aomega  = 1.5D0*aexp*omega_m        ! Numerical coeff for S_0 in psi Eq.

   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in restrict_coeff_fine_reverse_gr_nl. Please check.'  
      call clean_stop
   end if

   allocate(n_masked(1:active(ifinelevel)%ngrid))
   n_masked(1:active(ifinelevel)%ngrid)=0
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr 
         if(f(icell_f_amr,3)<=0d0) n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)               ! amr fine-grid index
         icell_f_amr=igrid_f_amr+iskip_f_amr                            ! amr fine-cell index

         ! Calculate space-dependent coefficients
         if(igr==6) then 
            gr_b = aomega*(rho(icell_f_amr)+2.0D0*gr_mat(icell_f_amr,4)+5.0D0*(1.0D0-gr_pot(icell_f_amr,5)/2ac2)**6) + &
                   1.75D0*gr_mat(icell_amr,1)/(1.0D0-gr_pot(icell_amr,5)/2ac2)**6
            gr_b = gr_b/(1.0D0-gr_pot(icell_f_amr,5)/2ac2)**2/2ac2
         end if

         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)                                ! amr coarse-cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax        ! amr coarse-grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                           ! cpu for coarse cell

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)                              ! mg coarse-grid index
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                               ! mg coarse-cell index

         ! If coarse cell masked, it is boundary and R rho is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Restriction to compute the rho value in the coarse cell
         if(igr==5) then
            rho1=gr_mat(icell_f_amr,1)  ! This is to restrict A_ijA^ij in \Psi Eq.
         else
            rho1=gr_b                   ! This is to restrict \xi coeff in \xi Eq.    
         end if    
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)+rho1/(dtwotondim-dble(n_masked(igrid_f_mg)))
      end do
   end do
  
   deallocate(n_masked)   

end subroutine restrict_coeff_fine_reverse_gr_nl

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_fine_gr_nl(ifinelevel,igr)
   use amr_commons
   use poisson_commons
   use gr_commons
   implicit none
   integer, intent(in) :: ifinelevel,igr

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr
   integer  :: ngrid_f, istart, nbatch
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg
   integer  :: icoarselevel, ind_c, cpu_amr

   real(dp) :: a, b, c, d, coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector), save               :: igrid_f_amr, icell_amr
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_father_grids
   real(dp), dimension(1:nvector), save               :: corr

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

   ! Sanity check for non-linear GR cases
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range in interpolate_and_correct_fine_gr_nl. Please check.'  
      call clean_stop
   end if

   ! Loop over fine grids by vector sweeps
   ngrid_f=active(ifinelevel)%ngrid
   do istart=1,ngrid_f,nvector

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active(ifinelevel)%igrid(istart+i-1)
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather(icell_amr,nbors_father_cells,nbors_father_grids, &
              nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax

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
               if(f(icell_amr(i),3)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
               ind_c       = (icell_c_amr-ncoarse-1)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               cpu_amr     = cpu_map(father(igrid_c_amr))
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               if(igrid_c_mg<=0) cycle

               icell_c_mg=(ind_c-1)*active_mg(cpu_amr,icoarselevel)%ngrid+igrid_c_mg
               ! only unmasked coarse cells contribute to the fine-cell correction
               if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0.0) cycle
               corr(i)=corr(i)+coeff*active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,1) &
                      &       -coeff*active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)
            end do
         end do

         ! Correct GR potential
         do i=1,nbatch
            gr_pot(icell_amr(i),igr) = gr_pot(icell_amr(i),igr)+corr(i)
         end do

      end do
      ! End loop over cells

   end do
   ! End loop over grids
end subroutine interpolate_and_correct_fine_gr_nl
