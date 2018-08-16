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
  integer::igrm
  !------------------------------------------------------------------------------------------
  ! This subroutine calculates the source terms div(V) or div(B) for some of the GR equations.
  !------------------------------------------------------------------------------------------
  integer::igrid,ngrid,ncache,i
  integer ,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(igr==4)then
     igrm=4
  else if(igr==10) then
     igrm=5        
  else      
     write(*,*) 'igr out of range in source_fine_gr_scalar. Please check.'     
     call clean_stop
  end if

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

111 format('   Entering source_fine_gr for level ',I2)

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
  integer::igrm,igrp
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector            ),save :: icelln
  integer ,dimension(1:nvector            ),save :: ind_cell
  integer ,dimension(1:nvector,0:twondim  ),save :: igridn
  real(dp),dimension(1:nvector            ),save :: phi1,phi2,phi3,phi4

  logical ,dimension(1:nvector,1:twotondim),save :: bdy
  real(dp),dimension(1:nvector,1:twotondim),save :: div_b_sons

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  
  if(igr==4)then
     igrm=4
  else if(igr==10)then
     igrm=5
  else      
     write(*,*) 'igr out of range in source_from_gr_pot_scalar. Please check.'     
     call clean_stop
  end if

  a=0.50D0*4.0D0/3.0D0/dx
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

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
 
  ! Gather father cells of central grids
  do i=1,ngrid
     icelln(i)=father(ind_grid(i))
  end do

  ! Interpolate div(B) from coarse level.
  ! Notice that div(V) was already interpolated to this level when in the coarse level. 
  if (ilevel>levelmin)then
     if(igr==10)then
        div_b_sons(1:nvector,1:twotondim)=0.0D0
        call interpol_source_gr_mat(icelln(1),div_b_sons(1,1),ngrid,ilevel,icount,5)
     else
        write(*,*) 'igr out of range for interpolation in source_from_gr_pot_scalar. Please check.'     
        call clean_stop
     end if
  end if

  bdy(1:nvector,1:twotondim)=.false.
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Set correct component to calculate div(V,B)
        igrp=idim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax

        ! Gather GR potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=gr_pot(igridn(i,ig1)+ih1,igrp)
           else
              bdy(i,ind)=.true.      
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=gr_pot(igridn(i,ig2)+ih2,igrp)
           else
              bdy(i,ind)=.true.      
           end if
        end do

        ! Calculate the div( ) contributions from the GR potential 
        if(.not.bdy(i,ind))then
           if(idim==1)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),igrm)=-a*(phi1(i)-phi2(i))
              end do
           else
              do i=1,ngrid
                 gr_mat(ind_cell(i),igrm)=gr_mat(ind_cell(i),igrm)-a*(phi1(i)-phi2(i))
              end do
           end if        
        end if
     end do ! End loop over idim
     
     if(bdy(i,ind).and.igr==10)then
        do i=1,ngrid
           gr_mat(ind_cell(i),igrm)=div_b_sons(i,ind)
        end do
     end if

  end do    ! End loop over cells

end subroutine source_from_gr_pot_scalar
