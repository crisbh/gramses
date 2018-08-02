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
  !----------------------------------------------------------
  !
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i
  integer ,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(igr==4)then
     igrm=4 
  else if(igr==7)then
     igrm=1 
  else if(igr==8)then
     igrm=2
  else if(igr==9)then
     igrm=3
  else if(igr==10)then
     igrm=4
  else      
     write(*,*) 'igr out of range in source_fine_gr. Please check.'     
     call clean_stop
  end if

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Compute source term from gr_pot
     call source_from_gr_pot_scalar(ind_grid,ngrid,ilevel,icount,igr)
  end do
  ! End loop over grids

  ! Update boundaries
  call make_virtual_fine_dp(gr_mat(1,igrm),ilevel)

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

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right
  real(dp):: ctilde,ctilde2,ac2

  ctilde   = sol/boxlen_ini/100000.0d0          ! Speed of light in code units
  ctilde2  = ctilde**2                          ! Speed of light squared
  ac2      = aexp**2*ctilde2                    ! (ac)^2 factor

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  
  if(igr==4)then
     igrm=4 
  else if(igr==7)then
     igrm=1 
  else if(igr==8)then
     igrm=2
  else if(igr==9)then
     igrm=3
  else if(igr==10)then
     igrm=4
  else      
     write(*,*) 'igr out of range in source_from_gr_pot. Please check.'     
     call clean_stop
  end if

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
        igrp=idim
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

        igrp=idim

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

        ! Calculate the force contribution from each gr_pot
        if(idim==1)then
           do i=1,ngrid
              gr_mat(ind_cell(i),igrm)=-a*(phi1(i)-phi2(i))+b*(phi3(i)-phi4(i))
           end do
        else
           do i=1,ngrid
              gr_mat(ind_cell(i),igrm)=gr_mat(ind_cell(i),igrm)-a*(phi1(i)-phi2(i))+b*(phi3(i)-phi4(i))
           end do
        end if        
     end do
  end do

end subroutine source_from_gr_pot_scalar
