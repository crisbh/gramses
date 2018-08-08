subroutine interpol_source_gr_mat(ind_cell,phi_int,ncell,ilevel,icount,igrm)
  use amr_commons
  use gr_commons, only:gr_mat
  implicit none
  integer::ncell,ilevel,icount,igrm
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! fixing boundary values of gr_mat on AMR boundary cellls.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order gr_mat interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0d0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*gr_mat(ind_cell(i),igrm)
           else
              add=coeff*gr_mat(indice,igrm)
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

end subroutine interpol_source_gr_mat
