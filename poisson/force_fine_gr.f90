!#########################################################
!#########################################################
subroutine force_fine_gr(ilevel,icount,igrp)
  use amr_commons
  use pm_commons
  use poisson_commons
  use gr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::ilevel,icount,igrp
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
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
        ! Compute gradient of gr_pot
        call gradient_gr_pot(ind_grid,ngrid,ilevel,icount,igrp)
     end do
     ! End loop over grids

#if NDIM==3
     if (sink)then
        call f_gas_sink(ilevel)
     end if
#endif
     ! Update boundaries
     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
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

111 format('   Entering force_fine_gr for level ',I2)

end subroutine force_fine_gr
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gradient_gr_pot(ind_grid,ngrid,ilevel,icount,igrp)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use gr_commons
  use gr_parameters
  implicit none

  integer::ngrid,ilevel,icount,igrp
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4,id5,id6
  integer::ig1,ig2,ig3,ig4,ig5,ig6
  integer::ih1,ih2,ih3,ih4,ih5,ih6
  integer::ix,iy,iz,inbor
  real(dp)::dx,a,b,scale,dx_loc,dx2
  integer,dimension(1:3,1:4,1:8)::ggg,hhh
  integer,dimension(1:3,1:2,1:8)::kkk,lll
  integer,dimension(1:3,1:4,1:8)::iii,jjj

  integer ,dimension(1:nvector                    ),save::ind_cell
  integer ,dimension(1:nvector,             1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim          ),save::igridn
  integer ,dimension(1:nvector,1:threetondim      ),save::igridn2
  real(dp),dimension(1:nvector                    ),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim, 1:ndim),save::phi_left,phi_right
  real(dp),dimension(1:nvector,1:twotondim, 1:27  ),save::pot_sons
  integer, dimension(1:nvector                    ),save::icelln     
  integer, dimension(1:nvector,1:threetondim      ),save::nbors_cells
  integer, dimension(1:nvector,1:twotondim        ),save::nbors_grids

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  dx2=dx**2

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_gr_pot(ind_left (1,idim),phi_left  (1,1,idim),ngrid,ilevel,icount,igrp)
        call interpol_gr_pot(ind_right(1,idim),phi_right (1,1,idim),ngrid,ilevel,icount,igrp)
     end do
  end if

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather GR potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=gr_pot(igridn(i,ig1)+ih1,igrp)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=gr_pot(igridn(i,ig2)+ih2,igrp)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=gr_pot(igridn(i,ig3)+ih3,igrp)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=gr_pot(igridn(i,ig4)+ih4,igrp)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do

        ! Calculate the 'force' contribution (-grad) from a given gr_pot
        do i=1,ngrid
           f(ind_cell(i),idim)=a*(phi1(i)-phi2(i))-b*(phi3(i)-phi4(i))
        end do

     end do ! End loop over idim
  end do    ! End loop over cells

  !---------------------------------------------------------------------------
  ! This block is used to calculate the force contribution \partial_i(beta^j).
  ! This is needed to update velocities.
  ! Recall that beta^j = B^j + \partial^j(b)
  !---------------------------------------------------------------------------

  if(igrp<7.or.igrp==10) return

  ! 27 neighbors in the 3-point kernel FDA.
  !   |direction
  !   | |node
  !   | | |cell
  !   v v v
  ! Parallel directions
  kkk(1,1,1:8)=(/13,14,13,14,13,14,13,14/); lll(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  kkk(1,2,1:8)=(/14,15,14,15,14,15,14,15/); lll(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  kkk(2,1,1:8)=(/11,11,14,14,11,11,14,14/); lll(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  kkk(2,2,1:8)=(/14,14,17,17,14,14,17,17/); lll(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  kkk(3,1,1:8)=(/ 5, 5, 5, 5,14,14,14,14/); lll(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  kkk(3,2,1:8)=(/14,14,14,14,23,23,23,23/); lll(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ! diagonal directions
  ! xy plane
  iii(1,1,1:8)=(/10,11,13,14,10,11,13,14/); jjj(1,1,1:8)=(/4,3,2,1,8,7,6,5/)
  iii(1,2,1:8)=(/11,12,14,15,11,12,14,15/); jjj(1,2,1:8)=(/4,3,2,1,8,7,6,5/)
  iii(1,3,1:8)=(/13,14,16,17,13,14,16,17/); jjj(1,3,1:8)=(/4,3,2,1,8,7,6,5/)
  iii(1,4,1:8)=(/14,15,17,18,14,15,17,18/); jjj(1,4,1:8)=(/4,3,2,1,8,7,6,5/)
  ! zx plane
  iii(2,1,1:8)=(/ 4, 5, 4, 5,13,14,13,14/); jjj(2,1,1:8)=(/6,5,8,7,2,1,4,3/)
  iii(2,2,1:8)=(/ 5, 6, 5, 6,14,15,14,15/); jjj(2,2,1:8)=(/6,5,8,7,2,1,4,3/)
  iii(2,3,1:8)=(/13,14,13,14,22,23,22,23/); jjj(2,3,1:8)=(/6,5,8,7,2,1,4,3/)
  iii(2,4,1:8)=(/14,15,14,15,23,24,23,24/); jjj(2,4,1:8)=(/6,5,8,7,2,1,4,3/)
  ! zy plane
  iii(3,1,1:8)=(/ 2, 2, 5, 5,11,11,14,14/); jjj(3,1,1:8)=(/7,8,5,6,3,4,1,2/)
  iii(3,2,1:8)=(/ 5, 5, 8, 8,14,14,17,17/); jjj(3,2,1:8)=(/7,8,5,6,3,4,1,2/)
  iii(3,3,1:8)=(/11,11,14,14,20,20,23,23/); jjj(3,3,1:8)=(/7,8,5,6,3,4,1,2/)
  iii(3,4,1:8)=(/14,14,17,17,23,23,26,26/); jjj(3,4,1:8)=(/7,8,5,6,3,4,1,2/)

  ! Gather father cells of the central grids
  do i=1,ngrid
     icelln(i)=father(ind_grid(i))
  end do
 
  ! Gather neighbouring father cells and store in nbors_cells 
  call get3cubefather(icelln,nbors_cells,nbors_grids,ngrid,ilevel)
  
  ! Gather neighboring grids
  do inbor=1,27
     do i=1,ngrid
        igridn2(i,inbor)=son(nbors_cells(i,inbor))
     end do
  end do

  ! Sanity check
  do i=1,ngrid
     if(igridn2(i,14).ne.ind_grid(i)) then
        write(*,*) 'Neighbouring grids incorrectly labelled in force_fine_gr'
        call clean_stop
     end if
  end do

  ! Interpolate b GR potential from upper level
  if(ilevel>levelmin) then
     do inbor=1,27
        if(inbor==14) cycle
  
        iz=(inbor          -1)/9
        iy=(inbor-9*iz     -1)/3
        ix=(inbor-9*iz-3*iy-1)

        if(ix.ne.1.and.iy.ne.1.and.iz.ne.1) cycle

        if(igrp==7.and.ix==1.and.iy.mod.2==0.and.iz.mod.2==0) cycle
        if(igrp==8.and.iy==1.and.iz.mod.2==0.and.ix.mod.2==0) cycle
        if(igrp==9.and.iz==1.and.ix.mod.2==0.and.iy.mod.2==0) cycle

        call interpol_gr_pot(nbors_cells(1,inbor),pot_sons(1,1,inbor),ngrid,ilevel,icount,10)
     end do
  end if

  ! Loop over fine cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
    
     do idim=1,ndim
        ! Loop over nodes
        ! For direct neighbours
        id5=lll(idim,1,ind); ig5=kkk(idim,1,ind); ih5=ncoarse+(id5-1)*ngridmax
        id6=lll(idim,2,ind); ig6=kkk(idim,2,ind); ih6=ncoarse+(id6-1)*ngridmax
        ! For diagonal neighbours
        id1=jjj(idim,1,ind); ig1=iii(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=jjj(idim,2,ind); ig2=iii(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=jjj(idim,3,ind); ig3=iii(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=jjj(idim,4,ind); ig4=iii(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather GR potential
        if(idim==igrp-6) then                                ! case of 2 direct neighbours
           ! Node 1
           do i=1,ngrid
              if(igridn2(i,ig5)>0) then
                 phi1(i)=gr_pot(igridn2(i,ig5)+ih5,10)
              else
                 phi1(i)=pot_sons(i,id5,ig5) 
              end if
           end do
           ! Node 2    
           do i=1,ngrid
              if(igridn2(i,ig6)>0) then
                 phi2(i)=gr_pot(igridn2(i,ig6)+ih6,10)
              else
                 phi2(i)=pot_sons(i,id6,ig6) 
              end if
           end do
        else                                                 ! case of 4 diagonal neighbours
           ! Node 1
           do i=1,ngrid
              if(igridn2(i,ig1)>0) then
                 phi1(i)=gr_pot(igridn2(i,ig1)+ih1,10)
              else
                 phi1(i)=pot_sons(i,id1,ig1) 
              end if
           end do
           ! Node 2    
           do i=1,ngrid
              if(igridn2(i,ig2)>0) then
                 phi2(i)=gr_pot(igridn2(i,ig2)+ih2,10)
              else
                 phi2(i)=pot_sons(i,id2,ig2) 
              end if
           end do
           ! Node 3
           do i=1,ngrid
              if(igridn2(i,ig3)>0) then
                 phi3(i)=gr_pot(igridn2(i,ig3)+ih3,10)
              else
                 phi3(i)=pot_sons(i,id3,ig3) 
              end if
           end do
           ! Node 4    
           do i=1,ngrid
              if(igridn2(i,ig4)>0) then
                 phi4(i)=gr_pot(igridn2(i,ig4)+ih4,10)
              else
                 phi4(i)=pot_sons(i,id4,ig4) 
              end if
           end do
        end if

        if(idim==igrp-6) then
           ! Laplacian of b     
           do i=1,ngrid
              f(ind_cell(i),idim)=f(ind_cell(i),idim)-(phi1(i)+phi2(i)-2.0D0*gr_pot(ind_cell(i),10))/dx2
           end do
        else
           ! Mixed derivatives of b (\partial_i\partial_j{b})     
           do i=1,ngrid
              f(ind_cell(i),idim)=f(ind_cell(i),idim)-(phi1(i)-phi2(i)-phi3(i)+phi4(i))/(4.0D0*dx2)
           end do
        end if

     end do ! End loop over idim
  end do    ! End loop over fine cells 

end subroutine gradient_gr_pot
