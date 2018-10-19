!#########################################################
!#########################################################
subroutine comp_gr_aij(ilevel,icount,ivect)
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
  ! This subroutine is a wrapper for gr_mat_aij_components
  ! to calculate the correct Aij component from a given ivect
  ! The result is stored in f(2) for later use
  !------------------------------------------------------------
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
     ! Compute a given A^ij component
     call gr_aij_components(ind_grid,ngrid,ilevel,icount,ivect)
  end do
  ! End loop over grids

111 format('   Entering comp_gr_aij for level ',I2)

end subroutine comp_gr_aij
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gr_aij_components(ind_grid,ngrid,ilevel,icount,ivect)
  use amr_commons
  use amr_parameters
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
  ! The result is stored in f(2) for later use
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc,inbor
  real(dp)::dx,dx2
  real(dp)::scale,dx_loc

  integer                        :: id1,id2
  integer                        :: ig1,ig2
  integer                        :: ih1,ih2
  integer, dimension(1:9,1:2,1:8) :: ggg,hhh
  integer, dimension(1:9,1:6    ) :: xxx,yyy
  real(dp),dimension(1:9,1:6    ) :: ppp,qqq

  integer, dimension(1:nvector                          ),save :: icelln
  integer, dimension(1:nvector                          ),save :: ind_cell
  integer, dimension(1:nvector,            1:threetondim),save :: igridn
  integer, dimension(1:nvector,            1:threetondim),save :: nbors_cells
  integer, dimension(1:nvector,1:twotondim              ),save :: nbors_grids
 
  real(dp),dimension(1:nvector                          ),save :: dv
  real(dp),dimension(1:nvector                          ),save :: pot1,pot2
  real(dp),dimension(1:nvector                          ),save :: aij
  logical, dimension(1:nvector                          ),save :: bdy
  
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


  xxx(1:9,1)=(/1,0,0,0,1,0,0,0,1/); yyy(1:9,1)=(/1,0,0,0,0,0,0,0,0/)
  xxx(1:9,2)=(/0,1,0,1,0,0,0,0,0/); yyy(1:9,2)=(/0,0,0,1,1,0,0,0,0/)
  xxx(1:9,3)=(/0,0,1,0,0,0,1,0,0/); yyy(1:9,3)=(/0,0,0,0,0,1,1,0,0/)
  xxx(1:9,4)=(/1,0,0,0,1,0,0,0,1/); yyy(1:9,4)=(/0,1,0,0,0,0,0,0,0/)
  xxx(1:9,5)=(/0,0,0,0,0,1,0,1,0/); yyy(1:9,5)=(/0,0,0,0,0,0,0,1,1/)
  xxx(1:9,6)=(/1,0,0,0,1,0,0,0,1/); yyy(1:9,6)=(/0,0,1,0,0,0,0,0,0/)

  ppp(1:9,1)=(/ 1.5D0, 0.0D0, 0.0D0, 0.0D0,-0.5D0, 0.0D0, 0.0D0, 0.0D0,-0.5D0/)
  ppp(1:9,2)=(/ 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
  ppp(1:9,3)=(/ 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0/)
  ppp(1:9,4)=(/-0.5D0, 0.0D0, 0.0D0, 0.0D0, 1.5D0, 0.0D0, 0.0D0, 0.0D0,-0.5D0/)
  ppp(1:9,5)=(/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0/)
  ppp(1:9,6)=(/-0.5D0, 0.0D0, 0.0D0, 0.0D0,-0.5D0, 0.0D0, 0.0D0, 0.0D0, 1.5D0/)
 
  !U derivates
  qqq(1:9,1)=(/ 2.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
  qqq(1:9,2)=(/ 0.0D0, 0.0D0, 0.0D0, 2.0D0, 2.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
  qqq(1:9,3)=(/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 2.0D0, 2.0D0, 0.0D0, 0.0D0/)
  qqq(1:9,4)=(/ 0.0D0, 2.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
  qqq(1:9,5)=(/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 2.0D0, 2.0D0/)
  qqq(1:9,6)=(/ 0.0D0, 0.0D0, 2.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
  
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

  ! Loop over fine cells
  do ind=1,twotondim
     bdy(1:nvector)=.false.
     aij(1:nvector)=0.0D0
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do 
     
     ! Loop over the four GR potentials V1, V2, V3 and U
     do igrp=1,4
        ! Skip igrp cases that are not needed for a given ivect
        if(ivect==2.and.igrp==3) cycle
        if(ivect==3.and.igrp==2) cycle
        if(ivect==5.and.igrp==1) cycle

        ! Loop over 9 directions
        do idim=1,9
           ! Only do idim=4-9 for igrp=4
           if(idim>3.and.igrp<4) cycle
            
           if(igrp<=3) then
              if(xxx((idim-1)*3+igrp,ivect).eq.0) cycle
           else
              if(yyy( idim          ,ivect).eq.0) cycle
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
                 bdy(i)=.true.      
              end if
           end do
           ! Node 2    
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 pot2(i)=gr_pot(igridn(i,ig2)+ih2,igrp)
              else
                 bdy(i)=.true.      
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

           if(igrp<=3) then
              do i=1,ngrid
                 aij(i)=aij(i)+ppp((idim-1)*3+igrp,ivect)*dv(i)
              end do
           else  
              do i=1,ngrid
                 aij(i)=aij(i)+qqq(idim,           ivect)*dv(i)
              end do
           end if
        end do ! End loop over idim
        
     end do    ! End loop over igrp 
     
     do i=1,ngrid
        if(.not.bdy(i)) then
           f(ind_cell(i),2)=2.0D0*aij(i)*(1.0D0+gr_pot(ind_cell(i),6)/ac2/(1.0D0-0.5D0*gr_pot(ind_cell(i),5)))/(1.0D0-0.5D0*gr_pot(ind_cell(i),5)/ac2)**6          
        end if
     end do

  end do       ! End loop over fine cells

end subroutine gr_aij_components
