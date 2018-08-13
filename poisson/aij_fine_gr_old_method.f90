!#########################################################
!#########################################################
subroutine comp_aij_gr(ilevel,icount,ivect)
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
        call aij_fine_gr(ind_grid,ngrid,ilevel,icount,ivect)
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

end subroutine comp_aij_gr
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine aij_fine_gr(ind_grid,ngrid,ilevel,icount,ivect)
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
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:9,1:2,1:8)::ggg,hhh,ccc

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,1:27),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  integer,dimension(1:nvector              ),save :: icelln
  integer,dimension(1:nvector,1:threetondim),save :: nbors_cells
  integer,dimension(1:nvector,1:twotondim  ),save :: nbors_grids
  real(dp):: ctilde,ctilde2,ac2
 
  real(dp),dimension(1:nvector,1:6),save :: aij
  real(dp),dimension(1:nvector    ),save :: dv,v1_l,v2_l,v3_l,uu_l,v1_r,v2_r,v3_r,uu_r
  integer :: idv

  ctilde   = sol/boxlen_ini/100000.0d0          ! Speed of light in code units
  ctilde2  = ctilde**2                          ! Speed of light squared
  ac2      = aexp**2*ctilde2                    ! (ac)^2 factor

  aij(1:nvector,1:6)=0.0D0

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! 27 neighbors in the 3-point kernel FDA
  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
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
  ccc(1,1,1:8)=(/ 2, 1, 5,11, 4,10,14,13/)
  ccc(1,2,1:8)=(//)
  ccc(2,1,1:8)=(//)
  ccc(2,2,1:8)=(//)
  ccc(3,1,1:8)=(//)
  ccc(3,2,1:8)=(//)
  ! Diagonal directions
  ! xy plane
  ggg(4,1,1:8)=(/10,11,13,14,10,11,13,14/); hhh(4,1,1:8)=(/4,3,2,1,8,7,6,5/)
  ggg(4,2,1:8)=(/14,15,17,18,14,15,17,18/); hhh(4,2,1:8)=(/4,3,2,1,8,7,6,5/)
  ggg(5,1,1:8)=(/11,12,14,15,11,12,14,15/); hhh(5,1,1:8)=(/4,3,2,1,8,7,6,5/)
  ggg(5,2,1:8)=(/13,14,16,17,13,14,16,17/); hhh(5,2,1:8)=(/4,3,2,1,8,7,6,5/)
  ccc(4,1,1:8)=(//)
  ccc(4,2,1:8)=(//)
  ccc(5,1,1:8)=(//)
  ccc(5,2,1:8)=(//)
  ! zx plane
  ggg(6,1,1:8)=(/ 4, 5, 4, 5,13,14,13,14/); hhh(6,1,1:8)=(/6,5,8,7,2,1,4,3/)
  ggg(6,2,1:8)=(/14,15,14,15,23,24,23,24/); hhh(6,2,1:8)=(/6,5,8,7,2,1,4,3/)
  ggg(7,1,1:8)=(/ 5, 6, 5, 6,14,15,14,15/); hhh(7,1,1:8)=(/6,5,8,7,2,1,4,3/)
  ggg(7,2,1:8)=(/13,14,13,14,22,23,22,23/); hhh(7,2,1:8)=(/6,5,8,7,2,1,4,3/)
  ccc(6,1,1:8)=(//)
  ccc(6,2,1:8)=(//)
  ccc(7,1,1:8)=(//)
  ccc(7,2,1:8)=(//)
  ! zy plane
  ggg(8,1,1:8)=(/ 2, 2, 5, 5,11,11,14,14/); hhh(8,1,1:8)=(/7,8,5,6,3,4,1,2/)
  ggg(8,2,1:8)=(/14,14,17,17,23,23,26,26/); hhh(8,2,1:8)=(/7,8,5,6,3,4,1,2/)
  ggg(9,1,1:8)=(/ 5, 5, 8, 8,14,14,17,17/); hhh(9,1,1:8)=(/7,8,5,6,3,4,1,2/)
  ggg(9,2,1:8)=(/11,11,14,14,20,20,23,23/); hhh(9,2,1:8)=(/7,8,5,6,3,4,1,2/)
  ccc(8,1,1:8)=(//)
  ccc(8,2,1:8)=(//)
  ccc(9,1,1:8)=(//)
  ccc(9,2,1:8)=(//)

  ! Gather neighboring grids
  do i=1,ngrid
!    igridn(i,14)=ind_grid(i)
     icelln(i   )=father(ind_grid(i))
  end do
  
  call get3cubefather(icelln,nbors_cells,nbors_grids,ngrid,ilevel)
  
  ! Gather neighboring grids??
  do inbor=1,27
     do i=1,ngrid
        igridn(i,inbor)=son(nbors_cells(i,inbor))
     end do
  end do

  do i=1,ngrid
     if(igridn(i,14).ne.ind_grid(i)) then
        write(*,*) 'something wrong in aij_fine_gr'
        call clean_stop
     end if
  end do

  do ind=1,twotondim
     iksip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=ind_grid(i)+iskip 
     end do
     ! idim=1 x derivatives
     do idim=1,3
        ig1=ggg(idim,1,ind)
        ig2=ggg(idim,2,ind)
        id1=hhh(idim,1,ind)
        id2=hhh(idim,2,ind)
        ih1=ncoarse+(id1-1)*ngridmax
        ih2=ncoarse+(id2-1)*ngridmax
        do i=1,ngrid
           if(igridn(i,ig1)>0) then
              v1_l(i)=gr_pot(igridn(i,ig1)+ih1,1)
              v2_l(i)=gr_pot(igridn(i,ig1)+ih1,2)
              v3_l(i)=gr_pot(igridn(i,ig1)+ih1,3)
              uu_l(i)=gr_pot(igridn(i,ig1)+ih1,4)
           else
              v1_l(i)=0.0D0
              v2_l(i)=0.0D0
              v3_l(i)=0.0D0
              uu_l(i)=0.0D0
              do ind_average=1,twotondim
                 icell=nbors_cells(i,ccc(ind_average,ind))
                 v1_l(i)=v1_l(i)+bbbb(ind_average)*gr_pot(icell,1)
                 v2_l(i)=v2_l(i)+bbbb(ind_average)*gr_pot(icell,2)
                 v3_l(i)=v3_l(i)+bbbb(ind_average)*gr_pot(icell,3)
                 uu_l(i)=uu_l(i)+bbbb(ind_average)*gr_pot(icell,4)
              end do
           end if
           if(igridn(i,ig2)>0) then
              v1_r(i)=gr_pot(igridn(i,ig2)+ih2,1)
              v2_r(i)=gr_pot(igridn(i,ig2)+ih2,2)
              v3_r(i)=gr_pot(igridn(i,ig2)+ih2,3)
              uu_r(i)=gr_pot(igridn(i,ig2)+ih2,4)
           else
              v1_r(i)=0.0D0
              v2_r(i)=0.0D0
              v3_r(i)=0.0D0
              uu_r(i)=0.0D0
              do ind_average=1,twotondim
                 icell=father(igridn(i,ccc(ind_average,ind)))
                 v1_r(i)=v1_r(i)+bbbb(ind_average)*gr_pot(icell,1)
                 v2_r(i)=v2_r(i)+bbbb(ind_average)*gr_pot(icell,2)
                 v3_r(i)=v3_r(i)+bbbb(ind_average)*gr_pot(icell,3)
                 uu_r(i)=uu_r(i)+bbbb(ind_average)*gr_pot(icell,4)
              end do
           end if
        end do
        if(idim==1) then
           do i=1,ngrid
              dv(i)=(v1_r(i)-v1_l(i))/(2.0D0*dx)
              aij(i,1)=aij(i,1)+1.5D0*dv(i)
              aij(i,4)=aij(i,4)-0.5D0*dv(i)
              aij(i,6)=aij(i,6)-0.5D0*dv(i)
              dv(i)=(v2_r(i)-v2_l(i))/(2.0D0*dx)
              aij(i,2)=aij(i,2)+      dv(i)
              dv(i)=(v3_r(i)-v3_l(i))/(2.0D0*dx)
              aij(i,3)=aij(i,3)+      dv(i)
              dv(i)=(uu_r(i)+uu_l(i)-2.0D0*gr_pot(ind_cell(i),4))/dx**2
              aij(i,1)=aij(i,1)+2.0D0*dv(i)
           end do
        end if
        if(idim==2) then
           do i=1,ngrid
              dv(i)=(v1_r(i)-v1_l(i))/(2.0D0*dx)
              aij(i,2)=aij(i,2)+      dv(i)
              dv(i)=(v2_r(i)-v2_l(i))/(2.0D0*dx)
              aij(i,4)=aij(i,4)+1.5D0*dv(i)
              aij(i,1)=aij(i,1)-0.5D0*dv(i)
              aij(i,6)=aij(i,6)-0.5D0*dv(i)
              dv(i)=(v3_r(i)-v3_l(i))/(2.0D0*dx)
              aij(i,5)=aij(i,5)+      dv(i)
              dv(i)=(uu_r(i)+uu_l(i)-2.0D0*gr_pot(ind_cell(i),4))/dx**2
              aij(i,4)=aij(i,4)+2.0D0*dv(i)
           end do
        end if
        if(idim==3) then
           do i=1,ngrid
              dv(i)=(v1_r(i)-v1_l(i))/(2.0D0*dx)
              aij(i,3)=aij(i,3)+      dv(i)
              dv(i)=(v2_r(i)-v2_l(i))/(2.0D0*dx)
              aij(i,5)=aij(i,5)+      dv(i)
              dv(i)=(v3_r(i)-v3_l(i))/(2.0D0*dx)
              aij(i,6)=aij(i,6)+1.5D0*dv(i)
              aij(i,1)=aij(i,1)-0.5D0*dv(i)
              aij(i,4)=aij(i,4)-0.5D0*dv(i)
              dv(i)=(uu_r(i)+uu_l(i)-2.0D0*gr_pot(ind_cell(i),4))/dx**2
              aij(i,6)=aij(i,6)+2.0D0*dv(i)
           end do
        end if
     end do
     ! idim=4 y derivatives
     idim=4
     ig1=ggg(idim,1,ind)
     ig2=ggg(idim,2,ind)
     id1=hhh(idim,1,ind)
     id2=hhh(idim,2,ind)
     ih1=ncoarse+(id1-1)*ngridmax
     ih2=ncoarse+(id2-1)*ngridmax
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
              phi2(i)=phi_left(i,id2,idim)
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
              phi4(i)=phi_left(i,id4,idim)
           end if
        end do

        ! Calculate each A_ij component 
        do i=1,ngrid
           f(ind_cell(i),idim)=a*(phi1(i)-phi2(i))-b*(phi3(i)-phi4(i))
        end do
     end do
  end do

end subroutine aij_fine_gr
