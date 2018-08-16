!#########################################################
!#########################################################
subroutine comp_gr_mat5(ilevel,icount,ivect)
  use amr_commons
  use pm_commons
  use poisson_commons
  use gr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::ilevel,icount,ivect
  !------------------------------------------------------------
  ! This subroutine calls comp_gr_mat5 to calculate the correct
  ! Aij component and store it into gr_mat(5) for later use
  !------------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  if(gravity_type>0)then

     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim

           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather cell centre positions
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do

           ! Call analytical gravity routine
           call gravana(xx,ff,dx_loc,ngrid)

           ! Scatter variables
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=ff(i,idim)
              end do
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids

     ! Update boundaries
     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  !------------------------------
  ! Compute gradient of potential
  !------------------------------
  else
     ! Update physical boundaries
     call make_boundary_phi(ilevel)

     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Compute gr_mat(5) components
        call gr_mat5_components(ind_grid,ngrid,ilevel,icount,ivect)
     end do
     ! End loop over grids

#if NDIM==3
     if (sink)then
        call f_gas_sink(ilevel)
     end if
#endif
     ! Update boundaries
!     do idim=1,ndim
!        call make_virtual_fine_dp(f(1,idim),ilevel)
!     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  rho_loc =0.0; rho_all =0.0
  epot_loc=0.0; epot_all=0.0
  fourpi=4.0D0*ACOS(-1.0D0)
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx_loc**ndim/fourpi/2.0D0

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Loop over dimensions
        do idim=1,ndim
           do i=1,ngrid
              if(son(ind_cell(i))==0)then
                 epot_loc=epot_loc+fact*f(ind_cell(i),idim)**2
              end if
           end do
        end do
        ! End loop over dimensions
        do i=1,ngrid
           rho_loc=MAX(rho_loc,dble(abs(rho(ind_cell(i)))))
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(epot_loc,epot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(rho_loc ,rho_all ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     epot_loc=epot_all
     rho_loc =rho_all
#endif
     epot_tot=epot_tot+epot_loc
     rho_max(ilevel)=rho_loc

111 format('   Entering comp_gr_mat5 for level ',I2)

end subroutine comp_gr_mat5
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gr_mat5_components(ind_grid,ngrid,ilevel,icount,ivect)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use gr_commons
  use gr_parameters
  implicit none

  integer::ngrid,ilevel,icount,ivect
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the A_ij components from (U,V^i) 
  ! in grids ind_grid(:) at level ilevel, using a
  ! 3 nodes kernel (3 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  real(dp)::dx,dx2
  real(dp)::scale,dx_loc

  integer                        :: id1,id2
  integer                        :: ig1,ig2
  integer                        :: ih1,ih2
  integer,dimension(1:9,1:2,1:8) :: ggg,hhh

  integer, dimension(1:nvector                          ),save :: icelln
  integer, dimension(1:nvector                          ),save :: ind_cell
  integer, dimension(1:nvector,            1:threetondim),save :: igridn
  integer, dimension(1:nvector,            1:threetondim),save :: nbors_cells
  integer, dimension(1:nvector,1:twotondim              ),save :: nbors_grids
 
  real(dp),dimension(1:nvector                          ),save :: dv
  real(dp),dimension(1:nvector                          ),save :: pot1,pot2
  real(dp),dimension(1:nvector,1:twotondim              ),save :: aij
  logical, dimension(1:nvector,1:twotondim              ),save :: bdy
  
  integer :: sgn,igrp             
 
  real(dp):: ctilde,ctilde2,ac2

  ctilde   = sol/boxlen_ini/100000.0d0          ! Speed of light in code units
  ctilde2  = ctilde**2                          ! Speed of light squared
  ac2      = aexp**2*ctilde2                    ! (ac)^2 factor

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  dx2=dx**2

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! 27 neighbors in the 3-point kernel FDA
  !   |direction
  !   | |node
  !   | | |cell
  !   v v v
  ! Parallel directions
  ggg(1,1,1:8)=(/13,14,13,14,13,14,13,14/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/14,15,14,15,14,15,14,15/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/11,11,14,14,11,11,14,14/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/14,14,17,17,14,14,17,17/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/ 5, 5, 5, 5,14,14,14,14/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/14,14,14,14,23,23,23,23/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ! Diagonal directions
  ! xy plane
  ggg(4,1,1:8)=(/10,11,13,14,10,11,13,14/); hhh(4,1,1:8)=(/4,3,2,1,8,7,6,5/)
  ggg(4,2,1:8)=(/14,15,17,18,14,15,17,18/); hhh(4,2,1:8)=(/4,3,2,1,8,7,6,5/)
  ggg(5,1,1:8)=(/11,12,14,15,11,12,14,15/); hhh(5,1,1:8)=(/4,3,2,1,8,7,6,5/)
  ggg(5,2,1:8)=(/13,14,16,17,13,14,16,17/); hhh(5,2,1:8)=(/4,3,2,1,8,7,6,5/)
  ! zx plane
  ggg(6,1,1:8)=(/ 4, 5, 4, 5,13,14,13,14/); hhh(6,1,1:8)=(/6,5,8,7,2,1,4,3/)
  ggg(6,2,1:8)=(/14,15,14,15,23,24,23,24/); hhh(6,2,1:8)=(/6,5,8,7,2,1,4,3/)
  ggg(7,1,1:8)=(/ 5, 6, 5, 6,14,15,14,15/); hhh(7,1,1:8)=(/6,5,8,7,2,1,4,3/)
  ggg(7,2,1:8)=(/13,14,13,14,22,23,22,23/); hhh(7,2,1:8)=(/6,5,8,7,2,1,4,3/)
  ! zy plane
  ggg(8,1,1:8)=(/ 2, 2, 5, 5,11,11,14,14/); hhh(8,1,1:8)=(/7,8,5,6,3,4,1,2/)
  ggg(8,2,1:8)=(/14,14,17,17,23,23,26,26/); hhh(8,2,1:8)=(/7,8,5,6,3,4,1,2/)
  ggg(9,1,1:8)=(/ 5, 5, 8, 8,14,14,17,17/); hhh(9,1,1:8)=(/7,8,5,6,3,4,1,2/)
  ggg(9,2,1:8)=(/11,11,14,14,20,20,23,23/); hhh(9,2,1:8)=(/7,8,5,6,3,4,1,2/)

  ! Gather father cells of the central grids
  do i=1,ngrid
     icelln(i)=father(ind_grid(i))
  end do
 
  ! Gather neighbouring father cells and store in nbors_cells 
  call get3cubefather(icelln,nbors_cells,nbors_grids,ngrid,ilevel)
  
  ! Gather neighboring grids
  do inbor=1,27
     do i=1,ngrid
        igridn(i,inbor)=son(nbors_cells(i,inbor))
     end do
  end do

  ! Sanity check
  do i=1,ngrid
     if(igridn(i,14).ne.ind_grid(i)) then
        write(*,*) 'Neighbouring grids incorrectly labelled in aij_fine_gr'
        call clean_stop
     end if
  end do

  ! Loop over the four GR potentials V1, V2, V3 and U
  do igrp=1,4

     ! Skip igrp cases that are not needed for a given ivect
     if(ivect==2.and.igrp==3) cycle
     if(ivect==3.and.igrp==2) cycle
     if(ivect==5.and.igrp==1) cycle

     bdy(1:nvector,1:twotondim)=.false.
     ! Loop over fine cells
     do ind=1,twotondim
        aij(1:nvector,1:6)=0.0D0
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Loop over 9 directions
        do idim=1,9
           ! Only do idim=4-9 for igrp=4
           if(idim>3.and.igrp<4) cycle

           ! Skip cases that are not needed for a given ivect
           if(ivect==1.or.ivect==4.or.ivect==6)then
              if(igrp<4.and.igrp.ne.idim                        ) cycle
              if(igrp=4.and.idim.ne.FLOOR(dble((ivect+0.5D0)/2))) cycle
              if(idim>3) cycle
           end if
           if(ivect==2)then
              if(idim==3.or. idim>5    ) cycle 
              if(igrp=<3.and.igrp==idim) cycle
              if(igrp==4.and.idim<4    ) cycle
           end if
           if(ivect==3)then
              if(idim.ne.1.or.idim.ne.3.or.idim.ne.6.or.idim.ne.7) cycle 
              if(igrp=<3.and.igrp==idim      ) cycle
              if(igrp==4.and.idim<6.or.idim>7) cycle
           end if
           if(ivect==5)then
              if(idim.ne.2.or.idim.ne.3.or.idim.ne.8.or.idim.ne.9) cycle 
              if(igrp=<3.and.igrp==idim) cycle
              if(igrp==4.and.idim<8    ) cycle
           end if

           ! Loop over nodes
           id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
           id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
           ! Gather GR potentials V1, V2, V3, U for igrp=1,2,3,4
           ! Node 1
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 pot1(i)=gr_pot(igridn(i,ig1)+ih1,igrp)
              else
                 bdy(i,1:8)=.true.      
              end if
           end do
           ! Node 2    
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 pot2(i)=gr_pot(igridn(i,ig2)+ih2,igrp)
              else
                 bdy(i,1:8)=.true.      
              end if
           end do

           if(idim>3) sgn=(-1)**(idim)

           ! For V1, V2 & V3, 3-point FDA to find gradient in direction idim
           if(igrp<4) then
              do i=1,ngrid
                 dv(i)=(pot2(i)-pot1(i))/(2.0D0*dx)
              end do
           ! For U, 3-point FDA to find second derivatives in 9 directions
           else
              ! For idim=1-3, 3-point FDA to find second derivative in direction idim
              if(idim<4) then 
                 do i=1,ngrid
                    dv(i)=(pot2(i)-2.0D0*gr_pot(ind_cell(i),4)+pot1(i))/dx2
                 end do
              ! For idim=4-9, 3-point FDA to find second derivative in diagonal directions
              else
                 do i=1,ngrid
                    dv(i)=(pot2(i)+pot1(i))/(4.0D0*dx2)*sgn
                 end do
              end if
           end if
 
           ! Calculate each A_ij component 
           if(igrp==1) then
              SELECT CASE (idim)
              CASE (1)
                 do i=1,ngrid
                    aij(i,1)=aij(i,1)+1.5D0*dv(i)
                    aij(i,4)=aij(i,4)-0.5D0*dv(i)
                    aij(i,6)=aij(i,6)-0.5D0*dv(i)
                 end do
              CASE (2)
                 do i=1,ngrid
                    aij(i,2)=aij(i,2)+      dv(i)
                 end do
              CASE (3)
                 do i=1,ngrid
                    aij(i,3)=aij(i,3)+      dv(i)
                 end do
              CASE DEFAULT
                 write(*,*) 'unsupported input value of idim!'
                 call clean_stop
              END SELECT
           end if

           if(igrp==2) then
              SELECT CASE (idim)
              CASE (1)
                 do i=1,ngrid
                    aij(i,2)=aij(i,2)+      dv(i)
                 end do
              CASE (2)
                 do i=1,ngrid
                    aij(i,1)=aij(i,1)-0.5D0*dv(i)
                    aij(i,4)=aij(i,4)+1.5D0*dv(i)
                    aij(i,6)=aij(i,6)-0.5D0*dv(i)
                 end do
              CASE (3)
                 do i=1,ngrid
                    aij(i,5)=aij(i,5)+      dv(i)
                 end do
              CASE DEFAULT
                 write(*,*) 'unsupported input value of idim!'
                 call clean_stop
              END SELECT
           end if

           if(igrp==3) then
              SELECT CASE (idim)
              CASE (1)
                 do i=1,ngrid
                    aij(i,3)=aij(i,3)+      dv(i)
                 end do
              CASE (2)
                 do i=1,ngrid
                    aij(i,5)=aij(i,5)+      dv(i)
                 end do
              CASE (3)
                 do i=1,ngrid
                    aij(i,1)=aij(i,1)-0.5D0*dv(i)
                    aij(i,4)=aij(i,4)-0.5D0*dv(i)
                    aij(i,6)=aij(i,6)+1.5D0*dv(i)
                 end do
              CASE DEFAULT
                 write(*,*) 'unsupported input value of idim!'
                 call clean_stop
              END SELECT
           end if

           if(igrp==4) then
              SELECT CASE (idim)
              CASE (1)
                 do i=1,ngrid
                    aij(i,1)=aij(i,1)+2.0D0*dv(i)
                 end do
              CASE (2)
                 do i=1,ngrid
                    aij(i,4)=aij(i,4)+2.0D0*dv(i)
                 end do
              CASE (3)
                 do i=1,ngrid
                    aij(i,6)=aij(i,6)+2.0D0*dv(i)
                 end do
              CASE (4,5)
                 do i=1,ngrid
                    aij(i,2)=aij(i,2)+2.0D0*dv(i)
                 end do
              CASE (6,7)
                 do i=1,ngrid
                    aij(i,3)=aij(i,3)+2.0D0*dv(i)
                 end do
              CASE (8,9)
                 do i=1,ngrid
                    aij(i,5)=aij(i,5)+2.0D0*dv(i)
                 end do
              CASE DEFAULT
                 write(*,*) 'unsupported input value of idim!'
                 call clean_stop
              END SELECT
           end if
        end do ! End loop over idim
        
     end do    ! End loop over fine cells
  end do       ! End loop over igrp

  ! Store one of the Aij component with proper factors for div(A^ij)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     do i=1,ngrid
        if(.not.bdy(i,ind))then
           gr_mat(ind_cell(i),5)=aij(i,ivect)
           gr_mat(ind_cell(i),5)=gr_mat(ind_cell(i),5)*(1.0D0+gr_pot(ind_cell(i),6)/ac2)/(1.0D0+gr_pot(ind_cell(i),5)/ac2)**6          
        end if
     end do
  end do

end subroutine gr_mat5_components
