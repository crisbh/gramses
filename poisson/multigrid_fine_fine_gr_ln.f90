! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-fine-level related routines
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     GR geometric field   gr_pot(:,igr)  active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_fine_gr_ln(ilevel,igr)
   ! Computes the residual the fine (AMR) level, and stores it into f(:,1)
   use amr_commons
   use poisson_commons
   use gr_commons
   implicit none
   integer, intent(in) :: ilevel,igr
   integer :: igrp

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2, phi_c, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)

   ! Set constants
   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid

   ! Set field index
   igrp = igr
   if(igr>6) igrp = igr-6

   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         phi_c = gr_pot(icell_amr,igrp)   ! Value of GR field on center cell
         nb_sum=0.0d0                     ! Sum of phi on neighbors

         ! SCAN FLAG TEST
         if(flag2(icell_amr)/ngridmax==0) then
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
                  nb_sum = nb_sum + gr_pot(icell_nbor_amr,igrp)
               end do
            end do
         else ! PERFORM SCAN
            if(f(icell_amr,3)<=0.0) then
               f(icell_amr,1)=0.0d0
               cycle
            end if
            do idim=1,ndim
               do inbor=1,2
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell !
                     ! Virtual phi value on unrefined neighbor cell :
                     ! -phi_c/mask_c
                     ! (simulates mask=-1.0 for the nonexistent refined cell)
                     nb_sum = nb_sum - phi_c/f(icell_amr,3)
                  else
                     ! Fetch neighbor cell
                     icell_nbor_amr = igrid_nbor_amr + &
                         (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                     if(f(icell_nbor_amr,3)<=0.0) then
                        ! Neighbor cell is masked :
                        ! compute its virtual phi with the mask
                        nb_sum = nb_sum + &
                            phi_c*(f(icell_nbor_amr,3)/f(icell_amr,3))
                     else
                        ! Neighbor cell is active, use its true potential
                        nb_sum = nb_sum + gr_pot(icell_nbor_amr,igrp)
                     end if
                  end if
               end do
            end do
         end if ! END SCAN TEST

         ! Store ***MINUS THE RESIDUAL*** in f(:,1), using BC-modified RHS
         f(icell_amr,1) = -oneoverdx2*( nb_sum - dtwondim*phi_c )+f(icell_amr,2)
      end do
   end do

end subroutine cmp_residual_mg_fine_gr_ln

! ##################################################################
! ##################################################################

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_fine_gr_ln(ilevel,redstep,igr)
   use amr_commons
   use pm_commons
   use poisson_commons
   use gr_commons
   implicit none
   integer, intent(in) :: ilevel,igr
   integer :: igrp
   logical, intent(in) :: redstep

   integer, dimension(1:3,1:2,1:8) :: iii, jjj
   integer, dimension(1:3,1:4)     :: ired, iblack

   real(dp) :: dx2, nb_sum, weight
   integer  :: ngrid
   integer  :: ind, ind0, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)

   ! Set constants
   dx2  = (0.5d0**ilevel)**2

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
   
   ! Set field index
   igrp = igr
   if(igr>6) igrp = igr-6

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

         nb_sum=0.0d0                       ! Sum of GR field on neighbors

         ! Read scan flag
         if(flag2(icell_amr)/ngridmax==0) then
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
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if
                  icell_nbor_amr = igrid_nbor_amr + &
                      (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  nb_sum = nb_sum + gr_pot(icell_nbor_amr,igrp)
               end do
            end do
            ! Update the potential, solving for potential on icell_amr
            gr_pot(icell_amr,igrp) = (nb_sum - dx2*f(icell_amr,2)) / dtwondim
         else
            ! Use the finer "solve" Gauss-Seidel near boundaries,
            ! with all necessary checks
            if (f(icell_amr,3)<=0.0) cycle
            if (safe_mode(ilevel) .and. f(icell_amr,3)<1.0) cycle

            weight=0.0d0 ! Central weight for "Solve G-S"
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

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell,
                     ! set mask=-1 on nonexistent neighbor cell
                     weight = weight - 1.0d0/f(icell_amr,3)
                  else
                     ! Fetch neighbor cell
                     icell_nbor_amr = igrid_nbor_amr + &
                         (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                     if(f(icell_nbor_amr,3)<=0.0) then
                        ! Neighbor cell is masked
                        weight = weight + f(icell_nbor_amr,3)/f(icell_amr,3)
                     else
                        ! Neighbor cell is active, increment neighbor sum
                        nb_sum = nb_sum + gr_pot(icell_nbor_amr,igrp)
                     end if
                  end if
               end do
            end do
            ! Update the potential, solving for potential on icell_amr
            gr_pot(icell_amr,igrp) = (nb_sum - dx2*f(icell_amr,2)) / (dtwondim - weight)
         end if
      end do
   end do
end subroutine gauss_seidel_mg_fine_gr_ln

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_fine_gr_ln(ifinelevel,igr)
   use amr_commons
   use poisson_commons
   use gr_commons
   implicit none
   integer, intent(in) :: ifinelevel,igr
   integer :: igrp

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

   ! Set field index
   igrp = igr
   if(igr>6) igrp = igr-6

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
               corr(i)=corr(i)+coeff*active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,1)
            end do
         end do

         ! Correct potential
         do i=1,nbatch
            gr_pot(icell_amr(i),igrp) = gr_pot(icell_amr(i),igrp)+corr(i)
         end do

      end do
      ! End loop over cells

   end do
   ! End loop over grids
end subroutine interpolate_and_correct_fine_gr_ln


