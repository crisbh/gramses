subroutine source_fine_gr_vector(ilevel,icount,ivect)
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
  !---------------------------------------------------------------------------------
  ! This subroutine calls the next one to calculate the div(A^ij) source terms 
  !---------------------------------------------------------------------------------
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
     ! Compute source term from gr_pot
     call source_from_gr_pot_vector(ind_grid,ngrid,ilevel,icount,ivect)
  end do
  ! End loop over grids

  ! Update boundaries when vector source terms are completed
  !if(ivect==6)then
  !   call make_virtual_fine_dp(gr_mat(1,1),ilevel)
  !   call make_virtual_fine_dp(gr_mat(1,2),ilevel)
  !   call make_virtual_fine_dp(gr_mat(1,3),ilevel)
  !end if

111 format('   Entering source_fine_gr_vector for level ',I2)

end subroutine source_fine_gr_vector
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine source_from_gr_pot_vector(ind_grid,ngrid,ilevel,icount,ivect)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use gr_commons
  use gr_parameters
  implicit none

  integer::ngrid,ilevel,icount,ivect
  integer,dimension(1:nvector)::ind_grid
  !---------------------------------------------------------------------------------
  ! This subroutine calculates div(A^ij) source terms 
  !---------------------------------------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  real(dp)::dx,dx2
  real(dp)::scale,dx_loc

  integer                        :: id1,id2
  integer                        :: ig1,ig2
  integer                        :: ih1,ih2
  integer,dimension(1:9,1:2,1:8) :: ggg,hhh

  integer, dimension(1:nvector                          ),save :: icelln
  integer ,dimension(1:nvector                          ),save :: ind_cell
  integer ,dimension(1:nvector,0:twondim                ),save :: igridn
  integer, dimension(1:nvector,            1:threetondim),save :: nbors_cells
  integer, dimension(1:nvector,1:twotondim              ),save :: nbors_grids
  real(dp),dimension(1:nvector                          ),save :: phi1,phi2

  logical, dimension(1:nvector,1:twotondim              ),save :: bdy
  logical, dimension(1:nvector,1:twotondim              ),save :: div_aij_sons

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  
  ! 27 neighbors in the 3-point kernel FDA.
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
        write(*,*) 'Neighbouring grids incorrectly labelled in source_fine_gr_vector'
        call clean_stop
     end if
  end do

  ! Decide which gr_mat component to interpolate
  interpol_div_aij=.false.
  if     (ivect==3)then
     interpol_div_aij=.true.
     igrm=1
  else if(ivect==5)then
     interpol_div_aij=.true.
     igrm=2
  else if(ivect==6)then
     interpol_div_aij=.true.
     igrm=3
  end if
  ! Interpolate potential from upper level
  if(ilevel>levelmin)then
     if(interpol_div_aij)then
        div_aij_sons(1:nvector,1:twotondim)=0.0D0
        do inbor=1,27 
           call interpol_source_gr_mat(icelln(1),div_aij_sons(1,1),ngrid,ilevel,icount,igrm)
        end do
     end if
  end if

  bdy(1:nvector,1:twotondim)=.false.
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over 9 directions
     do idim=1,9

        ! Skip cases that are not needed for a given ivect
        if(ivect==1.or.ivect==4.or.ivect==6)then
           if(idim>3) cycle
        end if
        if(ivect==2)then
           if(idim==3.or. idim>5    ) cycle 
        end if
        if(ivect==3)then
           if(idim.ne.1.or.idim.ne.3.or.idim.ne.6.or.idim.ne.7) cycle 
        end if
        if(ivect==5)then
           if(idim.ne.2.or.idim.ne.3.or.idim.ne.8.or.idim.ne.9) cycle 
        end if

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax

        ! Gather A_ij
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              if(idim<4) phi1(i)=gr_mat(igridn(i,ig1)+ih1,5)
           else
              bdy(i,1:8)=.true.      
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              if(idim<4) phi2(i)=gr_mat(igridn(i,ig2)+ih2,5)
           else
              bdy(i,1:8)=.true.      
           end if
        end do
        
        ! Calculate source term contribution from each one of the 6 A_ij components.
        ! This is only done when not on boundary cells nor next to them. 
        if(.not.bdy(i,ind))then
           if(ivect==1)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),1)=-(phi1(i)-phi2(i))/(2.0D0*dx)
              end do
           end if
           if(ivect==2.and.idim==2)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),1)=-(phi1(i)-phi2(i))/(2.0D0*dx)+gr_mat(ind_cell(i),1)
              end do
           end if
           if(ivect==2.and.idim==1)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),2)=-(phi1(i)-phi2(i))/(2.0D0*dx)
              end do
           end if
           if(ivect==3.and.idim==3)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),1)=-(phi1(i)-phi2(i))/(2.0D0*dx)+gr_mat(ind_cell(i),1)
              end do
           end if
           if(ivect==3.and.idim==1)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),3)=-(phi1(i)-phi2(i))/(2.0D0*dx)
              end do
           end if
           if(ivect==4)then        
              do i=1,ngrid
                 gr_mat(ind_cell(i),2)=-(phi1(i)-phi2(i))/(2.0D0*dx)+gr_mat(ind_cell(i),2)
              end do
           end if
           if(ivect==5.and.idim==3)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),2)=-(phi1(i)-phi2(i))/(2.0D0*dx)+gr_mat(ind_cell(i),2)
              end do
           end if
           if(ivect==5.and.idim==2)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),3)=-(phi1(i)-phi2(i))/(2.0D0*dx)+gr_mat(ind_cell(i),3)
              end do
           end if
           if(ivect==6)then
              do i=1,ngrid
                 gr_mat(ind_cell(i),3)=-(phi1(i)-phi2(i))/(2.0D0*dx)+gr_mat(ind_cell(i),3)
              end do
           end if
        end if

     end do ! End loop over idim
  
     if(interpolate_div_aij.and.bdy(i,ind)) gr_mat(ind_cell(i),igrm)=div_aij_sons(i,ind)

  end do    ! End loop over fine cells

end subroutine source_from_gr_pot_vector
