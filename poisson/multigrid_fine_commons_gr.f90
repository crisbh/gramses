! ------------------------------------------------------------------------
! Multigrid GR solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all generic fine multigrid routines, such as
!   * multigrid iterations @ fine and coarse MG levels
!   * communicator building
!   * MPI routines
!   * helper functions
!
! Used variables:
!                           finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     GR field                gr_pot(:,igr )  active_mg(myid,ilevel)%u(:,1)
!     Matter RHS for GR       gr_mat(:,igrm)  active_mg(myid,ilevel)%u(:,2)
!     residual                f(:,1)          active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS         f(:,2)                   N/A
!     mask                    f(:,3)          active_mg(myid,ilevel)%u(:,4)
!     restricted NL GR field   N/A            active_mg(myid,ilevel)%u(:,5)
!     restricted dens          N/A            active_mg(myid,ilevel)%u(:,6)
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
! Main multigrid routine for GR fields, called by amr_step
! ------------------------------------------------------------------------

subroutine multigrid_fine_gr(ilevel,icount,igr)
   use amr_commons
   use poisson_commons
   use poisson_parameters
   use gr_commons
   use gr_parameters

   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ilevel,icount,igr
   integer :: ivect
   integer, parameter  :: MAXITER  = 50
   real(dp), parameter :: SAFE_FACTOR = 0.5

   integer :: ifine, i, iter, icpu
   logical :: allmasked, gr_lin
   real(kind=8) :: err, last_err
   real(kind=8) :: res_norm2, i_res_norm2, residual, residual_old
   real(dp) :: n_cell_f,n_cell_f_tot,n_cell_c,n_cell_c_tot
   real(dp) :: trunc_norm2,trunc_norm2_tot,trunc_err

#ifndef WITHOUTMPI
   logical :: allmasked_tot
   real(kind=8) :: res_norm2_tot, i_res_norm2_tot
   integer :: info
#endif

   if(gravity_type>0)return
   if(numbtot(1,ilevel)==0)return

   gr_lin = .true.
   if(igr==5.or.igr==6) gr_lin = .false.

   if(verbose) print '(A,I2,A,I2)','Entering fine multigrid GR at level' ,ilevel,' for igr=',igr

   residual_old = 1.0d0
   ! ---------------------------------------------------------------------
   ! Prepare first guess, mask and BCs at finest level
   ! ---------------------------------------------------------------------

   if(ilevel>levelmin)then
      call make_initial_gr(ilevel,icount,igr)            ! Interpolate GR field
   endif

   call make_virtual_fine_dp(gr_pot(:,igr),ilevel)       ! Update boundaries
   ! call make_boundary_gr(ilevel)                       ! Update physical boundaries
   
   if(igr==1     ) then
      call make_fine_mask    (ilevel)                    ! Fill the fine mask
      call make_virtual_fine_dp(f(:,3),ilevel)           ! Communicate mask
      call make_boundary_mask(ilevel)                    ! Set mask to -1 in phys bounds
   end if
  
   if(igr==4     ) then
      call source_fine_gr_scalar (ilevel,icount,igr)     ! Calculate div(V)
      call make_virtual_fine_dp(f(:,2),ilevel)           ! Update boundaries
   else if(igr==5) then
      call source_fine_gr_aij_aij(ilevel,icount    )     ! Calculate A_ij*A^ij
      call make_virtual_fine_dp(gr_mat(:,1),ilevel)      ! Update boundaries
   else if(igr==7) then
      do ivect=1,6
         call comp_gr_aij          (ilevel,icount,ivect) ! Calculate argument for div(A_ij) vector sources
         call make_virtual_fine_dp(f(:,2),ilevel)        ! Update boundaries
         call source_fine_gr_vector(ilevel,icount,ivect) ! Calculate div(A_ij) vector sources
      end do
      call make_virtual_fine_dp(gr_mat(:,2),ilevel)      ! Update boundaries
      call make_virtual_fine_dp(gr_mat(:,3),ilevel)      ! Update boundaries
      call make_virtual_fine_dp(gr_mat(:,4),ilevel)      ! Update boundaries
   else if(igr==10) then
      call source_fine_gr_scalar(ilevel,icount,igr)      ! Calculate div(B)
      call make_virtual_fine_dp(f(:,2),ilevel)           ! Update boundaries
   end if

   if(gr_lin) then
      call make_fine_bc_rhs_gr_ln(ilevel,icount,igr)     ! Fill BC-modified RHS for linear
   end if

   ! ---------------------------------------------------------------------
   ! Build communicators up
   ! ---------------------------------------------------------------------

   if(igr==1) then
      ! @ finer level
      call build_parent_comms_mg(active(ilevel),ilevel)
      ! @ coarser levels
      do ifine=(ilevel-1),2,-1
          call build_parent_comms_mg(active_mg(myid,ifine),ifine)
      end do

      ! ---------------------------------------------------------------------
      ! Restrict mask up, then set scan flag
      ! ---------------------------------------------------------------------
      ! @ finer level

      if(ilevel>1) then
         ! Restrict and communicate mask
         call restrict_mask_fine_reverse(ilevel)
         call make_reverse_mg_dp(4,ilevel-1)
         call make_virtual_mg_dp(4,ilevel-1)

         ! Convert volume fraction to mask value
         do icpu=1,ncpu
            if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
            active_mg(icpu,ilevel-1)%u(:,4)=2d0*active_mg(icpu,ilevel-1)%u(:,4)-1d0
         end do

         ! Check active mask state
         if(active_mg(myid,ilevel-1)%ngrid>0) then
            allmasked=(maxval(active_mg(myid,ilevel-1)%u(:,4))<=0d0)
         else
            allmasked=.true.
         end if

         ! Allreduce on mask state
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(allmasked, allmasked_tot, 1, MPI_LOGICAL, &
              & MPI_LAND, MPI_COMM_WORLD, info)
         allmasked=allmasked_tot
#endif
      else
         allmasked=.true.
      endif

      ! @ coarser levels
      ! Restrict mask and compute levelmin_mg in the process
      if (.not. allmasked) then
         levelmin_mg=1
         do ifine=(ilevel-1),2,-1

            ! Restrict and communicate mask
            call restrict_mask_coarse_reverse(ifine)
            call make_reverse_mg_dp(4,ifine-1)
            call make_virtual_mg_dp(4,ifine-1)

            ! Convert volume fraction to mask value
            do icpu=1,ncpu
               if(active_mg(icpu,ifine-1)%ngrid==0) cycle
               active_mg(icpu,ifine-1)%u(:,4)=2d0*active_mg(icpu,ifine-1)%u(:,4)-1d0
            end do

            ! Check active mask state
            if(active_mg(myid,ifine-1)%ngrid>0) then
               allmasked=(maxval(active_mg(myid,ifine-1)%u(:,4))<=0d0)
            else
               allmasked=.true.
            end if

            ! Allreduce on mask state
#ifndef WITHOUTMPI
            call MPI_ALLREDUCE(allmasked,allmasked_tot,1,MPI_LOGICAL, &
                    & MPI_LAND,MPI_COMM_WORLD,info)
            allmasked=allmasked_tot
#endif

            if(allmasked) then ! Coarser level is fully masked: stop here
               levelmin_mg=ifine
               exit
            end if
         end do
      else
         levelmin_mg=ilevel
      end if
      if(nboundary>0)levelmin_mg=max(levelmin_mg,2)

      ! Update flag with scan flag
      call set_scan_flag_fine(ilevel)
      do ifine=levelmin_mg,ilevel-1
         call set_scan_flag_coarse(ifine)
      end do
   end if ! (igr==1)   

   ! ---------------------------------------------------------------------
   ! Initiate solve at fine level
   ! ---------------------------------------------------------------------

   iter = 0
   err = 1.0d0
   main_iteration_loop: do
      iter=iter+1
      ! Pre-smoothing
      if(gr_lin) then
         do i=1,ngs_fine_gr_ln_pre
            if (i.eq.1 .and. iter.eq.1) then
               src_mean = 0.0d0
            end if
            call gauss_seidel_mg_fine_gr_ln(ilevel,.true. ,igr) ! Red step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field
            call gauss_seidel_mg_fine_gr_ln(ilevel,.false.,igr) ! Black step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field 
         end do
      else
         do i=1,ngs_fine_gr_nl_pre
            if (i.eq.1 .and. iter.eq.1) then
               src_mean = 0.0d0
            end if
            call gauss_seidel_mg_fine_gr_nl(ilevel,.true. ,igr) ! Red step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field
            call gauss_seidel_mg_fine_gr_nl(ilevel,.false.,igr) ! Black step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field 
            call cmp_source_mean_gr_nl(ilevel,igr)              ! Mean of rhs for regularisation
         end do
      end if

      ! Compute residual and restrict into upper level RHS
      if(gr_lin) then
         call cmp_residual_mg_fine_gr_ln(ilevel,igr)
         call make_virtual_fine_dp(f(1,1),ilevel) ! Communicate residual
      else
         call cmp_residual_mg_fine_gr_nl(ilevel,igr)
         call make_virtual_fine_dp(f(1,1),ilevel) ! Communicate residual
      end if   

      if(iter==1.and.gr_lin) then
         call cmp_residual_norm2_fine(ilevel,i_res_norm2)
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1, &
                 & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
         i_res_norm2=i_res_norm2_tot
#endif
      end if

      ! First clear the rhs in coarser reception comms
      do icpu=1,ncpu
         if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
         active_mg(icpu,ilevel-1)%u(:,2)=0.0d0
         if(.not.gr_lin) then
            active_mg(icpu,ilevel-1)%u(:,5)=0.0d0
         end if
      end do

      ! Restrict and do reverse-comm
      call restrict_residual_fine_reverse(ilevel)
      call make_reverse_mg_dp(2,ilevel-1) ! communicate rhs
      call make_virtual_mg_dp(2,ilevel-1) ! communicate rhs
       
      if(.not.gr_lin)then
         call restrict_gr_pot_fine_reverse_gr_nl(ilevel,igr)
         call make_reverse_mg_dp(5,ilevel-1) ! communicate gr_pot
         call make_virtual_mg_dp(5,ilevel-1) ! communicate gr_pot 
         call make_physical_rhs_coarse_gr_nl(ilevel-1,igr) 
      end if

      if(ilevel>1) then
         ! Reset correction at upper level before solve
         do icpu=1,ncpu
            if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
            if(gr_lin) then
               active_mg(icpu,ilevel-1)%u(:,1)=0.0d0
            else
               active_mg(icpu,ilevel-1)%u(:,1)=active_mg(icpu,ilevel-1)%u(:,5)
            end if   
         end do

         ! Multigrid-solve the upper level
         if(gr_lin) then
            call recursive_multigrid_coarse   (ilevel-1, safe_mode(ilevel))
         else
            call recursive_multigrid_coarse_gr(ilevel-1, safe_mode(ilevel),igr)
         end if

         ! Interpolate coarse solution and correct fine solution
         if(gr_lin) then
            call interpolate_and_correct_fine_gr_ln(ilevel,igr)
         else
            call interpolate_and_correct_fine_gr_nl(ilevel,igr)
         end if
         call make_virtual_fine_dp(gr_pot(1,igr),ilevel)  ! Communicate GR field
      end if

      ! Post-smoothing
      if(gr_lin) then
         do i=1,ngs_fine_gr_ln_pst
            call gauss_seidel_mg_fine_gr_ln(ilevel,.true. ,igr) ! Red step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field
            call gauss_seidel_mg_fine_gr_ln(ilevel,.false.,igr) ! Black step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field
         end do
      else
         do i=1,ngs_fine_gr_nl_pst
            call gauss_seidel_mg_fine_gr_nl(ilevel,.true. ,igr) ! Red step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field
            call gauss_seidel_mg_fine_gr_nl(ilevel,.false.,igr) ! Black step
            call make_virtual_fine_dp(gr_pot(1,igr),ilevel)     ! Communicate GR field
            call cmp_source_mean_gr_nl(ilevel,igr)              ! Calculate rhs mean for regularisation 
         end do
      end if    

      ! Update fine residual
      if(gr_lin) then
         call cmp_residual_mg_fine_gr_ln(ilevel,igr)
      else
         call cmp_residual_mg_fine_gr_nl(ilevel,igr)
      end if
      call make_virtual_fine_dp(f(1,1),ilevel) ! Communicate residual

      if(gr_lin) then
         call cmp_residual_norm2_fine(ilevel,res_norm2)
      else
         call cmp_residual_norm2_fine_gr_nl(ilevel,res_norm2,n_cell_f)
      end if

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(res_norm2,res_norm2_tot,1, &
              & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      res_norm2=res_norm2_tot

      if(.not.gr_lin) then   
         call MPI_ALLREDUCE(n_cell_f ,n_cell_f_tot ,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
         n_cell_f  = n_cell_f_tot
      end if
#endif

      if(gr_lin) then
         last_err = err
         err = dsqrt(res_norm2/(i_res_norm2+1d-20))
      else
         residual = dsqrt(res_norm2)/n_cell_f
         call cmp_uvar_norm2_coarse_gr_nl(2,ilevel-1,trunc_norm2,n_cell_c)
      
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(trunc_norm2,trunc_norm2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
         call MPI_ALLREDUCE(n_cell_c   ,n_cell_c_tot   ,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
         trunc_norm2 = trunc_norm2_tot
         n_cell_c    = n_cell_c_tot
#endif

      trunc_err = dsqrt(trunc_norm2)/n_cell_c
      end if    
      
      if(gr_lin) then
         ! Verbosity
         if(verbose) print '(A,I5,A,I5,A,1pE10.3)','   ==> Step=', &
            & iter,'GR_pot:',igr,' Error=',err
        
         ! Converged? 
         if(err<epsilon_gr .or. iter>=MAXITER) exit

         ! Not converged, check error and possibly enable safe mode for the level
         if(err > last_err*SAFE_FACTOR .and. (.not. safe_mode(ilevel))) then
            if(verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
            safe_mode(ilevel) = .true.
         end if
      else
         ! Verbosity
         if(verbose) print '(A,I4,A,1e13.6,A,1e13.6,A,1e13.6,A,1e13.6)','   ==> Step=',iter,', Residual =',residual,', Truncation error= ',trunc_err, ', src_mean =', src_mean

         ! Converged? 
         if(residual<1.0d-8 .or. residual<0.001d0*trunc_err .or. dabs(residual-residual_old)<1.0d-10) exit

         if(iter>=MAXITER) then
            write(*,*) 'iter has exceeded MAXITER; the code is not converging!'
            call clean_stop
         end if
 
         residual_old = residual

      end if    

   end do main_iteration_loop

   if(gr_lin) then
      if(myid==1) print '(A,I5,A,I5,A,1pE10.3)','   ==> Level=',ilevel, &
             ' igr=',igr, ' Error=',err
      if(myid==1 .and. iter==MAXITER) print *,'WARN: Fine multigrid GR &
      &equation failed to converge...'
   else
      if(myid==1) print '(A,I5,A,I5,A,1E15.6, A, 1e15.6)','   ==> Level=',ilevel,' igr=',&
               igr,' Residual=',residual,' Truncation error=',trunc_err
   end if
   
   ! ---------------------------------------------------------------------
   ! Do NOT cleanup MG levels until all 10 GR fields are solved.
   ! Instead simply reset some containers to use them for the next GR field.
   ! Cleanup MG levels only after solve is complete.
   ! ---------------------------------------------------------------------
   if(igr<10) then
      do ifine=1,ilevel-1
         call restore_mg_level_gr(ifine,igr)
      end do
      call restore_amr_level_gr(ilevel)
   else
      do ifine=1,ilevel-1
         call cleanup_mg_level(ifine)
      end do
   end if

end subroutine multigrid_fine_gr

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------
recursive subroutine recursive_multigrid_coarse_gr(ifinelevel, safe, igr)
   use amr_commons
   use poisson_commons
   use gr_commons
   use gr_parameters
   
   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ifinelevel, igr
   logical, intent(in) :: safe

   integer :: i, icpu, icycle, ncycle
   logical :: gr_lin

   gr_lin = .true.
   if(igr==5 .or. igr==6) gr_lin = .false.

   if(verbose) then
      if(ifinelevel>levelmin_mg) then
         print '(A,I2,A,I2)','V-cycle: entering coarse multigrid at level ',ifinelevel, ' for gr_pot ',igr
      else
         print '(A,I2,A,I2)','V-cycle: entering coarsest multigrid level ' ,ifinelevel, ' for gr_pot ',igr
      end if
   end if
   
   if(ifinelevel<=levelmin_mg) then
      ! Solve 'directly' :
      if(gr_lin) then
         do i=1,ngs_coarse_gr_ln_pst
            call gauss_seidel_mg_coarse(ifinelevel,safe,.true. ) ! Red step
            call make_virtual_mg_dp(1,ifinelevel)                ! Communicate 
            call gauss_seidel_mg_coarse(ifinelevel,safe,.false.) ! Black step
            call make_virtual_mg_dp(1,ifinelevel)                ! Communicate 
         end do
      else
         do i=1,ngs_coarse_gr_nl_pst
            call gauss_seidel_mg_coarse_gr_nl(ifinelevel,safe,.true. ,igr)  ! Red step
            call make_virtual_mg_dp(1,ifinelevel)                           ! Communicate
            call gauss_seidel_mg_coarse_gr_nl(ifinelevel,safe,.false.,igr)  ! Black step
            call make_virtual_mg_dp(1,ifinelevel)                           ! Communicate
         end do
      end if
      return
   end if

   if(safe) then
      ncycle=ncycles_coarse_safe
   else
      ncycle=1
   endif

   do icycle=1,ncycle

      ! Pre-smoothing
      if(gr_lin) then
         do i=1,ngs_coarse_gr_ln_pre
            call gauss_seidel_mg_coarse(ifinelevel,safe,.true. ) ! Red step
            call make_virtual_mg_dp(1,ifinelevel)                ! Communicate 
            call gauss_seidel_mg_coarse(ifinelevel,safe,.false.) ! Black step
            call make_virtual_mg_dp(1,ifinelevel)                ! Communicate 
         end do
      else
         do i=1,ngs_coarse_gr_nl_pre
            call gauss_seidel_mg_coarse_gr_nl(ifinelevel,safe,.true. ,igr)  ! Red step
            call make_virtual_mg_dp(1,ifinelevel)                           ! Communicate
            call gauss_seidel_mg_coarse_gr_nl(ifinelevel,safe,.false.,igr)  ! Black step
            call make_virtual_mg_dp(1,ifinelevel)                           ! Communicate
         end do
      end if

      ! Compute residual and restrict into upper level RHS
      if(gr_lin) then
         call cmp_residual_mg_coarse      (ifinelevel)
      else
         call cmp_residual_mg_coarse_gr_nl(ifinelevel,igr)
      end if
      call make_virtual_mg_dp(3,ifinelevel)  ! Communicate residual

      ! First clear the rhs in coarser reception comms
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         active_mg(icpu,ifinelevel-1)%u(:,2)=0.0d0
         if(.not.gr_lin) then
            active_mg(icpu,ifinelevel-1)%u(:,5)=0.0d0
         end if
      end do

      ! Restrict and do reverse-comm
      call restrict_residual_coarse_reverse(ifinelevel)
      call make_reverse_mg_dp(2,ifinelevel-1) ! communicate rhs
      call make_virtual_mg_dp(2,ifinelevel-1)   
      if(.not.gr_lin) then
         call restrict_gr_pot_coarse_reverse_gr_nl(ifinelevel)
         call make_reverse_mg_dp(5,ifinelevel-1) ! communicate rhs
         call make_virtual_mg_dp(5,ifinelevel-1)
         call make_physical_rhs_coarse_gr_nl(ifinelevel-1,igr)
      end if

      ! Reset correction from upper level before solve
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         if(gr_lin) then 
            active_mg(icpu,ifinelevel-1)%u(:,1)=0.0d0
         else
            active_mg(icpu,ifinelevel-1)%u(:,1)=active_mg(icpu,ifinelevel-1)%u(:,5)
         end if
      end do

      ! Multigrid-solve the upper level
      call recursive_multigrid_coarse_gr(ifinelevel-1, safe, igr)

      ! Interpolate coarse solution and correct back into fine solution
      if(gr_lin) then
         call interpolate_and_correct_coarse      (ifinelevel)
      else
         call interpolate_and_correct_coarse_gr_nl(ifinelevel)
      end if
      call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution

      ! Post-smoothing
      if(gr_lin) then
         do i=1,ngs_coarse_gr_ln_pst
            call gauss_seidel_mg_coarse(ifinelevel,safe,.true. ) ! Red step
            call make_virtual_mg_dp(1,ifinelevel)                ! Communicate 
            call gauss_seidel_mg_coarse(ifinelevel,safe,.false.) ! Black step
            call make_virtual_mg_dp(1,ifinelevel)                ! Communicate 
         end do
      else
         do i=1,ngs_coarse_gr_nl_pst
            call gauss_seidel_mg_coarse_gr_nl(ifinelevel,safe,.true. ,igr) ! Red step
            call make_virtual_mg_dp(1,ifinelevel)                          ! Communicate
            call gauss_seidel_mg_coarse_gr_nl(ifinelevel,safe,.false.,igr) ! Black step
            call make_virtual_mg_dp(1,ifinelevel)                          ! Communicate
         end do
      end if
      
   end do

   if(verbose) print '(A,I2,A,I2)','V-cycle: leaving coarse multigrid at level ',ifinelevel, ' for gr_pot ',igr

end subroutine recursive_multigrid_coarse_gr

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Preprocess the fine (AMR) level RHS to account for boundary conditions
!
!  _____#_____
! |     #     |      Cell I is INSIDE active domain (mask > 0)
! |  I  #  O  |      Cell O is OUTSIDE (mask <= 0 or nonexistent cell)
! |_____#_____|      # is the boundary
!       #
!
! phi(I) and phi(O) must BOTH be set at call time, if applicable
! phi(#) is computed from phi(I), phi(O) and the mask values
! If AMR cell O does not exist, phi(O) is computed by interpolation
!
! Sets BC-modified RHS    into f(:,2)
!
! ------------------------------------------------------------------------
subroutine make_fine_bc_rhs_gr_ln(ilevel,icount,igr)

   use amr_commons
   use pm_commons
   use poisson_commons
   use gr_commons 
   use gr_parameters

   implicit none
   integer, intent(in) :: ilevel,icount,igr

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2, phi_b, nb_mask, nb_phi, w

   ! Arrays for vectorized interpol_phi
   real(dp), dimension(1:nvector,1:twotondim) :: phi_int
   integer,  dimension(1:nvector) :: ind_cell

   integer  :: ngrid, igrm
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr
   integer  :: ifathercell_nbor_amr

   integer  :: nx_loc
   real(dp) :: scale, fourpi

   ! Set constants
   nx_loc = icoarse_max-icoarse_min+1
   scale  = boxlen/dble(nx_loc)
   fourpi = 4.D0*ACOS(-1.0D0)*scale
   if(cosmo) fourpi = 1.5D0*omega_m*aexp*scale

   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid

   ! Pick correct gr_mat component where source term is stored 
   igrm=0
   if(igr<4           ) igrm=igr
   if(igr<10.and.igr>6) igrm=igr-5
   
   ! Calculate source mean value for regularisation
   call cmp_source_mean_gr_ln(ilevel,igrm)

   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         ! Init BC-modified RHS :
         if(igrm>0) then
            f(icell_amr,2) = gr_mat(icell_amr,igrm) - src_mean
         end if

         if(f(icell_amr,3)<=0.0) cycle ! Do not process masked cells

         ! Separate directions
         do idim=1,ndim
            ! Loop over the 2 neighbors
            do inbor=1,2
               ! Get neighbor grid shift
               igshift = iii(idim,inbor,ind)

               ! Get neighbor grid and its parent cell
               if(igshift==0) then
                  igrid_nbor_amr = igrid_amr
                  ifathercell_nbor_amr = father(igrid_nbor_amr)
               else
                  igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  ifathercell_nbor_amr = nbor(igrid_amr,igshift)
               end if

               if(igrid_nbor_amr==0) then
                  ! No neighbor: set mask to -1 and interp. phi
                  nb_mask = -1.0d0

                  ! Interpolate from upper level
                  ind_cell(1)=ifathercell_nbor_amr
                  call interpol_gr_pot(ind_cell,phi_int,1,ilevel,icount,igr)
                  nb_phi = phi_int(1,jjj(idim,inbor,ind))
               else
                  ! Fetch neighbor cell id
                  icell_nbor_amr = igrid_nbor_amr + (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  ! Check neighbor cell mask
                  nb_mask = f(icell_nbor_amr,3)
                  if(nb_mask>0) cycle ! Neighbor cell is active too: cycle
                  nb_phi  = gr_pot(icell_nbor_amr,igr)
               end if
               ! phi(#) interpolated with mask:
               w = nb_mask/(nb_mask-f(icell_amr,3)) ! Linear parameter
               phi_b = ((1.0d0-w)*nb_phi + w*gr_pot(icell_amr,igr))

               ! Increment correction for current cell
               f(icell_amr,2) = f(icell_amr,2) - 2.0d0*oneoverdx2*phi_b
            end do
         end do
      end do
   end do

end subroutine make_fine_bc_rhs_gr_ln

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
! ------------------------------------------------------------------------
! Multigrid level restoration for the Poisson equation solver next
! ------------------------------------------------------------------------

!subroutine restore_mg_level_gr(ilevel)
subroutine restore_mg_level_gr(ilevel,igr)
   use amr_commons
   use pm_commons
   use poisson_commons
   use gr_commons

   implicit none

   integer, intent(in) :: ilevel

   integer :: igrid, icpu, cur_grid, cur_cpu
   integer :: igr

   do icpu=1,ncpu
      if(active_mg(icpu,ilevel)%ngrid>0)then
         active_mg(icpu,ilevel)%u(:,1)=0.0d0
         active_mg(icpu,ilevel)%u(:,2)=0.0d0
         active_mg(icpu,ilevel)%u(:,3)=0.0d0
         if(igr==5.or.igr==6) then
            active_mg(icpu,ilevel)%u(:,5)=0.0d0  ! We only use these for the Non-linear equations
            active_mg(icpu,ilevel)%u(:,6)=0.0d0
         end if
      endif

      if(emission_mg(icpu,ilevel)%ngrid>0)then
         emission_mg(icpu,ilevel)%u(:,1)=0.0d0
         emission_mg(icpu,ilevel)%u(:,2)=0.0d0
         emission_mg(icpu,ilevel)%u(:,3)=0.0d0
         if(igr==5.or.igr==6) then
            emission_mg(icpu,ilevel)%u(:,5)=0.0d0  ! We only use these for the Non-linear equations
            emission_mg(icpu,ilevel)%u(:,6)=0.0d0
         end if
      endif

   end do

end subroutine restore_mg_level_gr

subroutine restore_amr_level_gr(ilevel)

   use amr_commons
   use pm_commons
   use poisson_commons
   use gr_commons 

   implicit none

   integer,intent(in) :: ilevel

   integer :: ngrid,igrid,ind,ind_cell,ind_grid,iskip,i

   ngrid = active(ilevel)%ngrid

   do ind=1,twotondim
      iskip = ncoarse+(ind-1)*ngridmax
      do igrid=1,ngrid
           ind_grid = active(ilevel)%igrid(igrid)
           ind_cell = ind_grid+iskip
           do i=1,ndim-1
              f(ind_cell,i) = 0.0d0
           end do
      end do
   end do

end subroutine restore_amr_level_gr
