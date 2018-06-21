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
   integer :: igrp

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2, dx2, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)

   real(dp) :: ctilde,ctilde2,omoverac2,a2K2,onevera4
   real(dp) :: potc,gr_a,gr_b,gr_c,op

   ! Set constants
   dx        = 0.5d0**ilevel
   oneoverdx2= 1.0d0/(dx*dx)
   dx2       = (0.5d0**ilevel)**2                ! Cell size squared
   ctilde    = sol/boxlen_ini/100000.0d0         ! Speed of light in code units
   ctilde2   = ctilde**2                         ! Speed of light squared
   omoverac2 = 0.75D0*omega_m/(a*ctilde2)        ! Numerical coeff for S_0
   a2K2      = a**2*K**2/12.0D0                  ! Background factor a^2K^2/12
   oneovera4 = 0.125D0/a**4                      ! Numerical coeff for A_ij^2   

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

   ! Calculate constant GR coefficients
   if(igrp==5) then 
      gr_c=a2K2
   else
      !TO DO
   end if

   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)     ! amr grid index
         icell_amr = iskip_amr + igrid_amr              ! amr cell index

         ! Calculate space-dependent GR coefficients
         if(igrp==5) then 
            gr_a =  a2K2 - omoverac2*rho(icell_amr)
            gr_b = -oneovera4*gr_mat(icell_amr,4)
         else
            !TO DO
         end if   

         ! SCAN FLAG TEST
         if(flag2(icell_amr)/ngridmax==0) then

            potc  = gr_pot(icell_amr,igrp)   ! Value of GR field on center cell
            nb_sum=0.0d0                     ! Sum of GR field on eighbors

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
         end if ! END SCAN TEST

         if(igrp==5) then 
            op = (nb_sum-6.0D0*potc)*(1.0D0+potc) &
                    - dx2*(gr_a + gr_b/(1.0D0+potc)**6 + gr_c*((1.0D0+potc)**6-1.0D0))
                   
            f(icell_amr,1) = op
         else
            op = !TO DO alpha eq

            f(icell_amr,1) = op
         end if

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
   real(dp) :: ctilde,ctilde2,omoverac2,a2K2,onevera4
   real(dp) :: potc,gr_a,gr_b,gr_c,op,dop

   ! Set constants
   dx2       = (0.5d0**ilevel)**2                ! Cell size squared
   ctilde    = sol/boxlen_ini/100000.0d0         ! Speed of light in code units
   ctilde2   = ctilde**2                         ! Speed of light squared
   omoverac2 = 0.75D0*omega_m/(a*ctilde2)        ! Numerical coeff for S_0
   a2K2      = a**2*K**2/12.0D0                  ! Background factor a^2K^2/12
   oneovera4 = 0.125D0/a**4                      ! Numerical coeff for A_ij^2   

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
   if(igr>6.or.igr<5) then
      print '(A)','igr out of range. Check.'             
      call clean_stop
   end if
   
   ! Calculate constant GR coefficients
   if(igrp==5) then 
      gr_c=a2K2
   else
      !TO DO 
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
         if(igrp==5) then 
            gr_a =  a2K2 - omoverac2*rho(icell_amr)
            gr_b = -oneovera4*gr_mat(icell_amr,4)
         else
            !TO DO
         end if

         ! Read scan flag
         if(flag2(icell_amr)/ngridmax==0) then
            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
            ! Those cells are active, have all their neighbors active
            ! and all neighbors are in the AMR+MG trees
            potc = gr_pot(icell_amr,igrp)
            nb_sum=0.0d0                       ! Sum of GR field on neighbors

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
            if(igrp==5) then 
               op = (nb_sum-6.0D0*potc)*(1.0D0+potc) &
                    - dx2*(gr_a + gr_b/(1.0D0+potc)**6 + gr_c*((1.0D0+potc)**6-1.0D0))
           
               dop = - 6.0D0 - 12.0D0*potc + 6.0D0*dx2*gr_b/(1.0D0+potc)**7 &
                     - 6.0D0*dx2*gr_c*(1.0D0+potc)**5
            else
               op = !TO DO alpha eq
            end if

            ! Update the GR potential, solving for potential on icell_amr
            gr_pot(icell_amr,igrp) = gr_pot(icell_amr,igrp) - op/dop
         end if
      end do
   end do
end subroutine gauss_seidel_mg_fine_gr_nl

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_fine_gr_nl(ifinelevel,igr)
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
end subroutine interpolate_and_correct_fine_gr_nl
