subroutine source_fine_gr_scalar(ilevel,icount,igr)
  use amr_commons
  use pm_commons
  use poisson_commons
  use gr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::ilevel,icount,igr
  !------------------------------------------------------------------------------------------
  ! This subroutine calculates the source terms div(V) or div(B) for some of the GR equations.
  !------------------------------------------------------------------------------------------
  integer::igrid,ngrid,ncache,i
  integer ,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Compute div( ) source term from gr_pot
     call source_from_gr_pot_scalar(ind_grid,ngrid,ilevel,icount,igr)
  end do
  ! End loop over grids

  ! Update boundaries
  ! call make_virtual_fine_dp(gr_mat(1,igrm),ilevel)

111 format('   Entering source_fine_gr_scalar for level ',I2)

end subroutine source_fine_gr_scalar
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine source_from_gr_pot_scalar(ind_grid,ngrid,ilevel,icount,igr)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use gr_commons
  use gr_parameters
  implicit none

  integer::ngrid,ilevel,icount,igr
  integer::igrp
  integer,dimension(1:nvector)::ind_grid
  !-----------------------------------------------------
  ! This routine computes div(V) or div(B) for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 3 nodes kernel (3 points FDA).
  !-----------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2
  integer::ig1,ig2
  integer::ih1,ih2
  real(dp)::dx,scale,dx_loc
  integer,dimension(1:3,1:2,1:8)::ggg,hhh

  integer ,dimension(1:nvector                   ),save :: ind_cell
  integer ,dimension(1:nvector,0:twondim         ),save :: igridn
  integer ,dimension(1:nvector,            1:ndim),save :: ind_left,ind_right
  real(dp),dimension(1:nvector,1:twotondim       ),save :: phi_left,phi_right
  real(dp),dimension(1:nvector                   ),save :: phi1,phi2
  real(dp),dimension(1:nvector                   ),save :: dvb

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  if(igr.ne.4.and.igr.ne.10) then
     write(*,*) 'igr out of range in source_from_gr_pot_scalar. Please check.'     
     call clean_stop
  end if

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

  ! Loop over cells
  do ind=1,twotondim
     dvb(1:nvector)=0.0D0  
     phi1(1:nvector)=0.0D0  
     phi2(1:nvector)=0.0D0  
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     do i=1,ngrid
        f(ind_cell(i),2)=0.0D0
     end do

     ! Loop over dimensions/igrp
     do idim=1,ndim
        ! Pick correct component to calculate derivative along a given direction
        if(igr==4)  igrp=idim
        if(igr==10) igrp=idim+6 
        
        ! Interpolate potential from upper level
        if (ilevel>levelmin)then
           call interpol_gr_pot(ind_left (1,idim),phi_left (1,1),ngrid,ilevel,icount,igrp)
           call interpol_gr_pot(ind_right(1,idim),phi_right(1,1),ngrid,ilevel,icount,igrp)
        end if

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax

        ! Gather GR potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=gr_pot(igridn(i,ig1)+ih1,igrp)
           else
              phi1(i)=phi_left (i,id1)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=gr_pot(igridn(i,ig2)+ih2,igrp)
           else
              phi2(i)=phi_right(i,id2)
           end if
        end do

        ! Accummulate the derivative along a given direction for div( )
        do i=1,ngrid
           dvb(i)=dvb(i)+(phi2(i)-phi1(i))/(2.0D0*dx)
        end do

     end do ! End loop over idim
    
     ! Store div(V) or div(B) directly into f(2) as a source term
     do i=1,ngrid
        f(ind_cell(i),2)=dvb(i)
     end do

  end do    ! End loop over cells

end subroutine source_from_gr_pot_scalar
