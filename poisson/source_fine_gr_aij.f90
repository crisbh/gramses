! ------------------------------------------------------------------------
! Multigrid extradof solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-fine-level related routines
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     --------------------------------------------------------------------
!     potential              sf              active_mg(myid,ilevel)%u(:,1)
!     physical RHS          rho              active_mg(myid,ilevel)%u(:,2)
!     truncation error      N/A              active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)            active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)            N/A
!     mask                 f(:,3)            active_mg(myid,ilevel)%u(:,4)
!     restricted sf         N/A              active_mg(myid,ilevel)%u(:,5)
!     restricted dens       N/A              active_mg(myid,ilevel)%u(:,6)
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
! Computate residual of the fine (AMR level) and store into f(:,1)
! ------------------------------------------------------------------------

subroutine cmp_source_fine_gr_aij(ilevel)

   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer, intent(in) :: ilevel

   real(dp) :: dx,dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg,idim,inbor
   integer  :: iskip_amr
   integer  :: igshift,igrid_nbor_amr,icell_nbor_amr
   real(dp) :: dtwondim = (twondim)

   real(dp) :: alpha,Rc,eta,sfc,ctilde,op
   real(dp), dimension(1:3) :: nb_sum_sfi
   integer,  dimension(1:nvector), save               :: igrid_amr,icell_amr
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_grids
   integer  :: nbatch,i,jnbor
   logical  :: bdy

   dx  = 0.5d0**ilevel
   dx2 = dx*dx
   ctilde = sol/boxlen_ini/100000.0d0
   Rc     = 0.5d0/sqrt(param_o)
   alpha  = 3.0d0*beta_bar*aexp**4/Rc**2

   ngrid=active(ilevel)%ngrid
  
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      do igrid_mg=1,ngrid,nvector
         nbatch=MIN(nvector,ngrid-igrid_mg+1)
         do i=1,nbatch
            igrid_amr(i) = active(ilevel)%igrid(i+igrid_mg-1)
            icell_amr(i) = iskip_amr+igrid_amr(i)
         end do
         call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
         do i=1,nbatch
            sfc = sf(icell_amr(i))
            if(flag2(icell_amr(i))/ngridmax==0) then
               bdy = .false.
               do jnbor=1,27
                  if(flag2(nbors_cells(i,jnbor))/ngridmax.ne.0) then
                     bdy = .true.
                  end if
               end do
               if(bdy==.true.) cycle
               
               ! Diagonal terms (A_11,A_22,A_33)
               if(ivect==1)then
                  nb_sum = sf(nbors_cells(i,13))+sf(nbors_cells(i,15))
               end if        
               if(ivect==4)then
                  nb_sum = sf(nbors_cells(i,11))+sf(nbors_cells(i,17))
               end if        
               if(ivect==6)then
                  nb_sum = sf(nbors_cells(i,5 ))+sf(nbors_cells(i,23)) 
               end if        
               
               ! Non-diagonal terms (A_12,A_12,A_23)
               if(ivect==2)then
                  uij = sf(nbors_cells(i,18))+sf(nbors_cells(i,10))-sf(nbors_cells(i, 16))-sf(nbors_cells(i,12)) !D_xD_yU
               end if
               if(ivect==3)then
                  uij = sf(nbors_cells(i,24))+sf(nbors_cells(i, 4))-sf(nbors_cells(i, 22))-sf(nbors_cells(i, 6)) !D_xD_zU
               end if
               if(ivect==5)then
                  uij = sf(nbors_cells(i,26))+sf(nbors_cells(i, 2))-sf(nbors_cells(i, 20))-sf(nbors_cells(i, 8)) !D_yD_zU
               end if

               if(ivect==1.or.ivect==4.or.ivect==6)then
                  uij = nb_sum-2.0*sfc 
               else
                  uij = uij/4.0D0     
               end if        

               f(icell_amr(i),1) = uij/dx2
            else
               f(icell_amr(i),1) = 0.0d0
            end if
         end do
      end do
   end do

end subroutine cmp_source_fine_gr_aij
