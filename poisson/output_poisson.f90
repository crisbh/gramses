subroutine backup_poisson(filename)
  use amr_commons
  use poisson_commons
  use gr_commons    

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer::ix,iy,iz                          ! for the coordinates
  real(dp)::dx                               ! for the coordinates
  real(dp),dimension(1:twotondim,1:3)::xc    ! cell coordinates relative to grids
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

#ifndef WITHOUTMPI
  integer,parameter::tag=1123
  integer::dummy_io,info2
#endif

  if(verbose)write(*,*)'Entering backup_poisson'

  ilun=ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

 ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)ndim+1
  write(ilun)nlevelmax
  write(ilun)nboundary
  do ilevel=1,nlevelmax ! CBH 17-02-20 Default RAMSES line
!  do ilevel=levelmin,nlevelmax  ! CBH 11-01-19
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax

              dx=0.5d0**ilevel
              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              if(ndim>0)xc(ind,1)=(dble(ix)-0.5d0)*dx
              if(ndim>1)xc(ind,2)=(dble(iy)-0.5d0)*dx
              if(ndim>2)xc(ind,3)=(dble(iz)-0.5d0)*dx

              ! Write Newtonian potential
              do i=1,ncache
                 xdp(i)=phi(ind_grid(i)+iskip)
              end do
              write(ilun)xdp

              ! Write force
              do ivar=1,ndim
                 do i=1,ncache
                    xdp(i)=f(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do

              ! Output GR-only arrays 17-02-20
              ! Notice that this block needs to match the READ block in init_poisson_gr.f90 in order for GRAMSES to restart from a given snapshot
              if(gr) then
                 ! Write GR potentials  
                 do ivar=1,10   
                    do i=1,ncache
                       xdp(i)=gr_pot(ind_grid(i)+iskip,ivar)
                    end do
                    write(ilun)xdp
                 end do

                 ! Write s0 
                 do i=1,ncache
                    xdp(i)=rho(ind_grid(i)+iskip)
                 end do
                 write(ilun)xdp
     
                 ! Write s_i and s=Tr(s_ij)  
                 do ivar=1,4    
                    do i=1,ncache
                       xdp(i)=gr_mat2(ind_grid(i)+iskip,ivar)
                    end do
                    write(ilun)xdp
                 end do

                 ! Write positions
                 do ivar=1,ndim
                    do i=1,ncache
                        xdp(i)=xg(ind_grid(i),ivar)+xc(ind,ivar)
                    end do
                    write(ilun)xdp
                 end do
              end if

           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

end subroutine backup_poisson
