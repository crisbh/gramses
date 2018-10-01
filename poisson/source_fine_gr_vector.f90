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
  ! This subroutine is a wrapper for the next one to calculate div(A^ij) source terms.
  ! The A^ij components are taken from the f(2) array.
  ! The result is stored into gr_mat(2-4).
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
     ! Compute source term from gr_mat
     call source_from_gr_mat_vector(ind_grid,ngrid,ilevel,icount,ivect)
  end do
  ! End loop over grids

  ! Update boundaries when vector source terms are completed
  !if(ivect==6)then
  !   call make_virtual_fine_dp(gr_mat(1,2),ilevel)
  !   call make_virtual_fine_dp(gr_mat(1,3),ilevel)
  !   call make_virtual_fine_dp(gr_mat(1,4),ilevel)
  !end if

111 format('   Entering source_fine_gr_vector for level ',I2)

end subroutine source_fine_gr_vector

!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine source_from_gr_mat_vector(ind_grid,ngrid,ilevel,icount,ivect)
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
  ! This subroutine calculates div(A^ij) source terms.
  ! The A^ij components are taken from the f(2) array.
  !---------------------------------------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc,igrm
  real(dp)::dx,dx2
  real(dp)::scale,dx_loc

  integer                        :: id1,id2
  integer                        :: ig1,ig2
  integer                        :: ih1,ih2
  integer,dimension(1:9,1:2,1:8) :: ggg,hhh
  integer,dimension(1:3,1:6    ) :: ok_idim

  real(dp) :: accl      

  integer, dimension(1:nvector                          ),save :: icelln
  integer ,dimension(1:nvector                          ),save :: ind_cell
  integer ,dimension(1:nvector,            1:threetondim),save :: igridn
  integer, dimension(1:nvector,            1:threetondim),save :: nbors_cells
  integer, dimension(1:nvector,1:twotondim              ),save :: nbors_grids
  real(dp),dimension(1:nvector                          ),save :: phi1,phi2

  logical, dimension(1:nvector                          ),save :: bdy
  logical, dimension(1:nvector,1:twotondim              ),save :: div_aij_sons
  logical, dimension(1:nvector,1:3                      ),save :: div_aij

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

  ! Useful derivatives of A_ij for a given ivect
  ok_idim(1:3,1)=(/2,0,0/)
  ok_idim(1:3,2)=(/3,2,0/)
  ok_idim(1:3,3)=(/4,0,2/)
  ok_idim(1:3,4)=(/0,3,0/)
  ok_idim(1:3,5)=(/0,4,3/)
  ok_idim(1:3,6)=(/0,0,4/)

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

  if(ilevel>levelmin.and.ivect==1) then
     do igrm=2,4 
        call interpol_gr_mat(icelln(1),div_aij_sons(1,1),ngrid,ilevel,icount,igrm)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
              gr_mat(ind_cell(i),igrm)=div_aij_sons(i,ind)
           end do
        end do
     end do
  end if

  !!!!!DO COMMUNICATIONS!!!!!!!

  bdy    (1:nvector    )=.false.
  div_aij(1:nvector,1:3)=0.0D0

  ! Loop over cells
  do ind=1,twotondim

     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over 3 directions
     do idim=1,3
        
        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax

        ! Initialise phi1, phi2 to 0
        phi1(1:nvector) = 0.0D0
        phi2(1:nvector) = 0.0D0

        ! Gather A_ij
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=f(igridn(i,ig1)+ih1,2)
           else
              bdy(i)=.true.      
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=f(igridn(i,ig2)+ih2,2)
           else
              bdy(i)=.true.      
           end if
        end do
        
        if(ok_idim(idim,ivect)>0) then
           do i=1,ngrid
              div_aij(i,idim)=(phi2(i)-phi1(i))/(2.0D0*dx)
           end do
        end if

     end do ! End loop over idim
     
     do idim=1,3
        if(ok_idim(idim,ivect)==0) cycle

        if(idim.eq.1) accl=0.0D0
        if(idim.ne.1) accl=1.0D0
        do i=1,ngrid
           if(.not.bdy(i)) then
              gr_mat(ind_cell(i),ok_idim(idim,ivect))=gr_mat(ind_cell(i),ok_idim(idim,ivect))*accl+div_aij(i,idim)
           end if
        end do
  
     end do 

  end do    ! End loop over fine cells

  !igrm(1:6,1)=(/2,3,4,0,0,0/)
  !igrm(1:6,2)=(/0,2,0,3,4,0/)
  !igrm(1:6,3)=(/0,0,2,0,3,4/)

end subroutine source_from_gr_mat_vector
