!-----------------------------------------------------------!
! Gravity Lightcone 
! Original implementation by Rasera and Reverdy
! Based on RAMSES particle LC implementations
!-----------------------------------------------------------!
!
!------------------ MODIF V. REVERDY 2011 ------------------!
!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! Select cells which are inside the lightcone and write data into files.
!----------------------------------------------------------------------!   
subroutine output_conegrav(is_fullsky,filedir,filename,cone_id,observer_x,observer_y,observer_z, &
                         & observer_redshift,cone_zmax,cone_mode,input_levelmin,input_levelmax,  &
                         & input_thetay,input_thetaz,input_theta,input_phi,future)
  use amr_commons
  use poisson_commons
  use gr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none

  ! Misc
  integer,parameter::tag=2800
  integer::dummy_io
  integer::info
  integer::ilun
  logical::is_opened

  ! Input
  logical::is_fullsky                                                   ! option to specify is the cone is fullsky or not
  character(LEN=*)::filedir                                             ! name of the output dir
  character(LEN=*)::filename                                            ! name of the output file
  integer::cone_id                                                      ! cone identifier
  real(dp)::observer_x                                                  ! x coordinate of the observer (usually 0.5)
  real(dp)::observer_y                                                  ! y coordinate of the observer (usually 0.5)
  real(dp)::observer_z                                                  ! z coordinate of the observer (usually 0.5)
  real(dp)::observer_redshift                                           ! redshift of the observer
  real(dp)::cone_zmax                                                   ! zmax of the cone
  integer::cone_mode                                                    ! mode = 0 -> only leaves, 1 -> only coarse, 2 -> all grid hierarchy
  integer::input_levelmin                                               ! levelmax in case of all grid hierarchy
  integer::input_levelmax                                               ! levelmax in case of all grid hierarchy
  real(dp)::input_thetay                                                ! angle if non-fullksy
  real(dp)::input_thetaz                                                ! angle if non-fullksy
  real(dp)::input_theta                                                 ! angle if non-fullksy
  real(dp)::input_phi                                                   ! angle if non-fullksy
  
  ! Level
  integer::cone_levelmin
  integer::cone_levelmax
  integer::cone_nlevel
  
  ! Variables
  real(kind=8)::z2                                                      ! old redshift
  real(kind=8)::z1                                                      ! recent redshift
  real(kind=8)::dist2                                                   ! distance max
  real(kind=8)::dist1                                                   ! distance min
  real(kind=8)::om0in                                                   ! cosmological param
  real(kind=8)::omLin                                                   ! cosmological param
  real(kind=8)::hubin                                                   ! cosmological param
  real(kind=8)::Lbox                                                    ! cosmological param
  real(kind=8)::observer(3)                                             ! position of the observer
  real(kind=8)::dx
  real(kind=8)::dxcell
  real(kind=8),parameter::disttol=0.D0
  
  ! Variables : filename
  character(LEN=200)::filename_info
  character(LEN=200)::filename_dat
  character(LEN=200)::filename_hdr
  character(LEN=5)::nchar_cpu
  character(LEN=5)::nchar_coarse
  character(LEN=5)::nchar_cone
  
  ! Variables : index
  integer::ilevel
  integer::ncache
  integer::idim
  integer::iskip
  
  ! Variables : cells centers
  integer::ind
  integer::ix
  integer::iy
  integer::iz
  real(kind=8),dimension(1:twotondim,1:3)::xc                           ! centers of cells
  
  ! Variables : common for arrays
! integer::nvar=3                                                       ! fx, fy, fz
! integer,parameter::nvar=5                                             ! fx, fy, fz, phi, rho
                                                                        ! 1    3          1    3          3       9 
! integer,parameter::nvar=20                                            ! Phi, grad(Phi), Psi, grad(Psi), beta^i, grad(beta^i) ! CBH_LC
  integer::nvar                                                         ! Baojiu: this is no longer a "parameter" but a normal variable
  integer,parameter::nbool=1                                            ! son grid index
  integer,parameter::nvectorgrid = floor(nvector/(2.D0**ndim))          ! nvector/8
  integer,dimension(1:nvectorgrid)::ind_grid                            ! indices of grids
  integer,dimension(1:nvector)::ind_cell                                ! indices of cells
  
  ! Variables : input arrays
  real(kind=8),dimension(1:nvector,1:ndim)::xc_in
! real(kind=8),dimension(1:nvector,1:nvar)::var_in                      ! Baojiu: commtend out
  real(kind=8),dimension(:,:),allocatable::var_in                       ! Baojiu: change to allocatable
  integer,dimension(1:nvector,1:nbool)::bool_in

  ! Variables : selection arrays
  real(kind=8),dimension(1:27*nvector,1:ndim)::xc_out
! real(kind=8),dimension(1:27*nvector,1:nvar)::var_out                  ! Baojiu: commented out
  real(kind=8),dimension(:,:),allocatable::var_out                      ! Baojiu: change to allocatable
  integer,dimension(1:27*nvector,1:nbool)::bool_out
  
  ! Variables : write arrays
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim)::xc_wr            ! x, y, z centers of cells for output
! real(kind=4),dimension(1:nstride+27*nvector,1:nvar)::var_wr           ! Baojiu: commtend out
  real(kind=4),dimension(:,:),allocatable::var_wr                       ! Baojiu: change to allocatable 
  integer,dimension(1:nstride+27*nvector,1:nbool)::bool_wr              ! son
 
  ! Hdr arrays
  integer::ihdr
  integer::hdr_arraysize
  integer,dimension(:),allocatable::hdr_ilevel
  real(kind=8),dimension(:),allocatable::hdr_dx
  integer,dimension(:),allocatable::hdr_ncellperlevel

  ! Variables : counters
  integer::ncell_level,ncell_total
  integer::igrid,ngrid,i
  integer::ivector,istride
  integer::icell_in,icell_out,icell_wr
  integer::ncell_in,ncell_out,ncell_wr
  
  ! Write variables
  integer(kind=8),dimension(1:2)::reduce_infolocal
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  integer::ierr
  
  ! Iogroupsize variables
  real(kind=8)::iovolume
  
  ! aexp for info file
  real(kind=8)::aendconem2_info,aendconem1_info,aendcone_info !expansion factor for end of the shell
  !the end of previous shell and the end of previous previous shell (which is the begining of the current shell)
  !when activating overlapping buffer option.
  real(kind=8)::zendconem2_info,zendconem1_info,zendcone_info
  real(kind=8)::dendconem2_info,dendconem1_info,dendcone_info

  real(kind=8)::Omega0,OmegaL,OmegaR,coverH0
  real(kind=8)::coord_distance
  integer::future
  integer::ibool                                                        ! Baojiu: added for looping over bool variables

  ! Baojiu: set nvar here depending on which version 
  if(.not.gr) then
     nvar=5                                                             ! 5  field quantities to output for Newtonian   simulations
  then
     nvar=20                                                            ! 20 field quantities to output for relativisic simulations
  endif
  ! There does not to be anything special about 27: it just accounts
  ! for the possibility that some cells may be chosen multiple times
  ! due to the replica of the simulation boxes to form the lightcone
  allocate(var_in (1:           nvector,1:nvar))                        ! allocate array to store cells to be sent for selection
  allocate(var_out(1:        27*nvector,1:nvar))                        ! allocate array to store cells that are selected
  allocate(var_wr (1:nstride+27*nvector,1:nvar))                        ! allocate array to store cells to be written
  ! Baojiu 

  ! Correct quantities
  cone_levelmin = input_levelmin                                        ! minimum level on which to output fields
  cone_levelmax = input_levelmax                                        ! maximum level on which to output fields
  if((cone_levelmin.eq.0).or.(cone_levelmin.gt.nlevelmax)) then
     cone_levelmin = levelmin                                           ! shouldn't be lower than the domain level
  endif
  if((cone_levelmax.eq.0).or.(cone_levelmax.gt.nlevelmax)) then
     cone_levelmax = nlevelmax                                          ! shouldn't be higher than the most-refined level
  endif
  cone_nlevel = (cone_levelmax-cone_levelmin+1)                         ! number of levels on which to output fields

  ! Filename
  ilun=30*ncpu+myid+10                                                  ! file handler
  call title(myid,nchar_cpu)                                            ! create file name
  call title(nstep_coarse,nchar_coarse)                                 ! create file name
  write(nchar_cone,'(I5.5)') cone_id
  !
  if(future.eq.0) then
     if(is_fullsky) then
        filename_info = TRIM(filedir)//'info_cone_grav_fullsky_pastnfut_'//nchar_cone//'_ncoarse_'//TRIM(nchar_coarse)//'.txt'
     else
        filename_info = TRIM(filedir)//'info_cone_grav_narrow_pastnfut_'//nchar_cone//'_ncoarse_'//TRIM(nchar_coarse)//'.txt'
     endif
  else
     print*,'cone grav not implemented future=',future
     stop
  endif
  
  filename_dat = TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.dat' ! create file name
  filename_hdr = TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.hdr' ! create file name

  ! Special cases
  if(nstep_coarse.lt.4) return
  
  ! Compute quantities
  z2=1.0D0/aexp_old-1.0D0                                               ! redshift of the further side of the shell
  z1=1.0D0/aexp    -1.0D0                                               ! redshift of the closer  side of the shell
  !
  observer(1)=observer_x                                                ! observer x coordinate
  observer(2)=observer_y                                                ! observer y coordinate
  observer(3)=observer_z                                                ! observer z coordinate
  !
  om0in=omega_m                                                         ! Omega_m
  omLin=omega_l                                                         ! Omega_Lambda
  hubin=h0/100.D0                                                       ! h
  Lbox=boxlen_ini/hubin                                                 ! box size in Mpc
  !
  if((use_aexp_restart).and.(nstep_coarse_after_restart.eq.2)) then
     z2=1.D0/aexp_restart_light_cone-1.D0
  endif
  !
  if(conegrav_overlap) then
     if((aendconem1.lt.aendconem2).or.(aendcone.lt.aendconem1))print*,'WARNING CONEGRAV SHELL RANGE NOT WELL ORDERED, AEXP ',aendconem2,aendconem1,aendcone 
     z2=1.D0/aendconem2-1.D0
     z1=1.D0/aendcone  -1.D0
     if (myid.eq.1) then
        print*,''
        print*,'TEST CONE LIMIT'
        print*,'nstepcoarse,z2,z1',nstep_coarse,z2,z1
        print*,'aendconem2,aendconem1,aendcone',aendconem2,aendconem1,aendcone
        print*,''
     endif
  endif

  if(z2.gt.cone_zmax) return
  if(dabs(z2-z1)<1.0D-6) return
  if((z1.lt.observer_redshift).and.(z2.lt.observer_redshift)) return

  ! Hdr arrays
  hdr_arraysize = cone_nlevel
  allocate(hdr_ilevel(1:cone_nlevel))
  allocate(hdr_dx(1:cone_nlevel))
  allocate(hdr_ncellperlevel(1:cone_nlevel))
  hdr_ilevel=0
  hdr_dx=0
  hdr_ncellperlevel=0

  ! Compute iogroupsize
  call perform_my_selection_simple_conegrav_ramsesunits(.true.,                         &
       &                                .false.,disttol,z1,z2,dist1,dist2,              &
       &                                om0in,omLin,hubin,Lbox,                         &
       &                                observer,observer_redshift,dxcell,nvar,nbool,   &
       &                                xc_in,var_in,bool_in,ncell_in,                  &
       &                                xc_out,var_out,bool_out,ncell_out,.false.)
  !
  iovolume = ((2.0D0*dist2)**3)-(((2.D0*dist1)/dsqrt(real(ndim,kind=8)))**3)
  if(.not.is_fullsky) iovolume = iovolume*(input_thetay*input_thetaz)/(41263.D0)
  if(adaptive_iogroupsize) IOGROUPSIZECONE = nint(iovolume*real(IOGROUPSIZE, kind=8))
  if(adaptive_iogroupsize.and.(.not.is_fullsky)) IOGROUPSIZECONE = IOGROUPSIZECONE/8
  if(IOGROUPSIZECONE.le.1) IOGROUPSIZECONE = 1
  if(IOGROUPSIZECONE.ge.IOGROUPSIZE) IOGROUPSIZECONE = IOGROUPSIZE
  if(verbose) write(*,*) cone_id, 'IOGROUPSIZECONE = ',IOGROUPSIZECONE, 'iovolume = ', iovolume

  ! Wait for the token
  if(conetoken) then 
     if(mod(myid-1,IOGROUPSIZECONE).ne.0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  endif
  
  ! Loop over levels
  reduce_infolocal = 0
  reduce_infoglobal = 0
  ncell_total = 0                                                       ! total number of cells on all relevant levels
  ihdr = 0
  is_opened = .FALSE.

  do ilevel=cone_levelmin,cone_levelmax                                 ! loop over all refinement levels of interest
     ! Initialization of counters
     ncell_level = 0                                                    ! ncells written for the given level
     ivector     = 0                                                    ! initialise local cell index
     istride     = 0
     icell_in    = 0                                                    ! for-selection-cells array pointer (unused)
     icell_out   = 0                                                    ! selected-cells array pointer (unused)
     icell_wr    = 0                                                    ! write-cell array pointer
     ncell_in    = 0                                                    ! number of cells sent for selection
     ncell_out   = 0                                                    ! number of cells finally selected
     ncell_wr    = 0                                                    ! number of cells written to file
    
     ! Initialization of arrays
     xc_in    = 0.D0                                                    ! array of cell coordinates to be sent for selection
     xc_out   = 0.D0                                                    ! array of cell coordinates finally selected
     xc_wr    = 0.D0                                                    ! array of cell coordinates to be written to file
     var_in   = 0.D0                                                    ! array of cell floating-point variables to be sent for selection
     var_out  = 0.D0                                                    ! array of cell floating-point variables finally selected
     var_wr   = 0.D0                                                    ! array of cell floating-point variables to be written to file
     !
     bool_in  = 0                                                       ! array of cell integer variables to be sent for selection
     bool_out = 0                                                       ! array of cell integer variables finally selected
     bool_wr  = 0                                                       ! array of cell integer variables to be written to file
    
     ! Index
     ihdr = ihdr+1
    
     ! Level constants
     dx=0.5D0**ilevel                                                   ! cell size on level ilevel (code unit)
     dxcell=dx                                                          ! need to check (maybe dxcell = (1/2.)*dx

     ! Set position of cell centers relative to grid center
     do ind=1,twotondim
        iz=(ind-1          )/4
        iy=(ind-1     -4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0) xc(ind,1)=(dble(ix)-0.5D0)*dx                        ! x distance of cell centre to grid centre (code unit)
        if(ndim>1) xc(ind,2)=(dble(iy)-0.5D0)*dx                        ! y distance of cell centre to grid centre (code unit)
        if(ndim>2) xc(ind,3)=(dble(iz)-0.5D0)*dx                        ! z distance of cell centre to grid centre (code unit)
     end do
     !
    
     ! Loop over myid grids by vectorgrid sweeps
     ncache=active(ilevel)%ngrid                                        ! total number of active grid at ilevels on my CPU
     do igrid=1,ncache,nvectorgrid
        ngrid=MIN(nvectorgrid,ncache-igrid+1)                           ! batch size for grid vectorisation
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)                  ! global index of grid
!          ind_grid(i)=active(ilevel)%pcomm%igrid(igrid+i-1)            ! RAMSES-LC line CBH_LC 16-02-2021
        end do
      
        ! Loop over cells
        ivector=0                                                       ! initialise local cell index
        do i=1,ngrid
           do ind=1,twotondim
              !
              ivector=ivector+1                                         ! increment local cell index
              !
              ! Gather cell indices
              iskip=ncoarse+(ind-1)*ngridmax
              ind_cell(ivector)=iskip+ind_grid(i)                       ! global index of cell
              ! 
              ! Compute cell centre coordinates
              do idim=1,ndim
                 xc_in(ivector,idim)=xg(ind_grid(i),idim)+xc(ind,idim)  ! gather cell coordinates in code units
              end do
              !
              ! Collect floating-point variables
              ! Baojiu: modified to be application to both Newtonian and GR simulations
              if(.not.gr) then
                 var_in(ivector,1)=f(ind_cell(ivector),1)               ! gather x component of force in cell
                 var_in(ivector,2)=f(ind_cell(ivector),2)               ! gather y component of force in cell
                 var_in(ivector,3)=f(ind_cell(ivector),3)               ! gather z component of force in cell
                 var_in(ivector,4)=phi(ind_cell(ivector))               ! gather potential in cell
                 var_in(ivector,5)=rho(ind_cell(ivector))               ! gather density in cell
              else
              ! CBH_LC
              ! Notice that gradients need to be copied from f() at the correct igrp sequence in move or synchro
              !At igrp = 10, f() = - grad_i(b)
                 var_in(ivector,1 )=gr_pot(ind_cell(ivector),5)         ! Phi
                 var_in(ivector,2 )=gr_pot(ind_cell(ivector),6)         ! Xi 
                 var_in(ivector,3 )=gr_pot(ind_cell(ivector),7)-f(ind_cell(ivector),1)  ! \beta^i=B^i+grad_i(b)
                 var_in(ivector,4 )=gr_pot(ind_cell(ivector),8)-f(ind_cell(ivector),2)
                 var_in(ivector,5 )=gr_pot(ind_cell(ivector),9)-f(ind_cell(ivector),3)
                 var_in(ivector,6 )=gr_pot_grad(ind_cell(ivector),1 )   ! \grad_i\Phi
                 var_in(ivector,7 )=gr_pot_grad(ind_cell(ivector),2 )
                 var_in(ivector,8 )=gr_pot_grad(ind_cell(ivector),3 )
                 var_in(ivector,9 )=gr_pot_grad(ind_cell(ivector),4 )   ! \grad_i\Xi
                 var_in(ivector,10)=gr_pot_grad(ind_cell(ivector),5 )
                 var_in(ivector,11)=gr_pot_grad(ind_cell(ivector),6 )
                 var_in(ivector,12)=gr_pot_grad(ind_cell(ivector),7 )   ! \grad_i\beta^j
                 var_in(ivector,13)=gr_pot_grad(ind_cell(ivector),8 )
                 var_in(ivector,14)=gr_pot_grad(ind_cell(ivector),9 )
                 var_in(ivector,15)=gr_pot_grad(ind_cell(ivector),10)
                 var_in(ivector,16)=gr_pot_grad(ind_cell(ivector),11)
                 var_in(ivector,17)=gr_pot_grad(ind_cell(ivector),12)
                 var_in(ivector,18)=gr_pot_grad(ind_cell(ivector),13)
                 var_in(ivector,19)=gr_pot_grad(ind_cell(ivector),14)
                 var_in(ivector,20)=gr_pot_grad(ind_cell(ivector),15)
              endif
              ! END CBH_LC
            
              ! Collect integer variables
              ! Baojiu: modified to work for general nbool
              do ibool=1,nbool
                 bool_in(ivector,ibool)=son(ind_cell(ivector))          ! gather index of the cells' son grid
              end do
              ! 
           end do 
        end do                                                          ! do i=1,ngrid 
        !
        ncell_in = ivector                                              ! number of cells to be sent for selection
        !
        ! Perfom selection on cells
        if(is_fullsky) then
           call perform_my_selection_simple_conegrav_ramsesunits(.false.,                 &
           &                                .false.,disttol,z1,z2,dist1,dist2,            &
           &                                om0in,omLin,hubin,Lbox,                       &
           &                                observer,observer_redshift,dxcell,nvar,nbool, &
           &                                xc_in,var_in,bool_in,ncell_in,                &
           &                                xc_out,var_out,bool_out,ncell_out,.false.)
        else
           call perform_my_selection_conegrav_ramsesunits(.false.,                        &
           &                                .false.,disttol,z1,z2,dist1,dist2,            &
           &                                om0in,omLin,hubin,Lbox,                       &
           &                                observer,observer_redshift,dxcell,nvar,nbool, &
           &                                xc_in,var_in,bool_in,ncell_in,                &
           &                                xc_out,var_out,bool_out,ncell_out,.false.,    &
           &                                input_thetay,input_thetaz,input_theta,input_phi)      
        endif
        !
        ! Transfer out_arrays of size ncell_out in write_arrays of size icell_wr
        ! Note that icell_wr is not reset between consecutive calls to this func
        if(ncell_out.gt.0) then
           call transfer_grav(nvar,nbool,ncell_out,icell_wr,ncell_level,ncell_total, &
                            & xc_out,var_out,bool_out,xc_wr,var_wr,bool_wr)
        end if
        !
        ! If nstride has been reached/exceeded, write arrays of size nstride and
        ! reset the write-array pointer istride
        if(icell_wr.ge.nstride) then
           call write_grav(.false.,is_opened,ilun,filename_dat,   &
                         & xc_wr,nvar,var_wr,nbool,bool_wr,       &
                         & icell_wr,cone_nlevel,reduce_infolocal)
        end if
        ! 
     end do                                                             ! do igrid=1,ncache,nvectorgrid (loop over grids)
    
     ! If some cells remain in write-array, write them
     if(icell_wr.gt.0) then
        call write_grav(.true.,is_opened,ilun,filename_dat,       &
                      & xc_wr,nvar,var_wr,nbool,bool_wr,          &
                      & icell_wr,cone_nlevel,reduce_infolocal)
     end if

     hdr_ilevel(ihdr)=ilevel
     hdr_dx(ihdr)=dx
     hdr_ncellperlevel(ihdr)=ncell_level
     !
  end do                                                                ! end of loop over levels
  !
  if(is_opened) close(ilun)                                             ! close file after writing complete
  !
  ! Write hdr files
  call write_hdr(ncell_total,filename_hdr,hdr_arraysize,hdr_ilevel,hdr_dx,hdr_ncellperlevel)
  
  ! Send the token
  if(conetoken) then 
     if(mod(myid,IOGROUPSIZECONE).ne.0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag,MPI_COMM_WORLD,info)
    end if
  endif

  ! Reduce
  reduce_infolocal(2) = ncell_total
  call MPI_REDUCE(reduce_infolocal,reduce_infoglobal,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD, ierr)
  
  ! Write info file
  if(myid.eq.1) then
     if(conegrav_overlap) then
        if(future.eq.0) then
           aendconem2_info=aendconem2
           aendconem1_info=aendconem1
           aendcone_info  =aendcone
           
           zendconem2_info=1.D0/aendconem2_info-1.D0
           zendconem1_info=1.D0/aendconem1_info-1.D0
           zendcone_info  =1.D0/aendcone_info  -1.D0
        else
           print*,'cone grav not implemented future=',future
           stop
        endif

        call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
        if (observer_redshift.gt.0) then
           dendconem2_info=((coord_distance(zendconem2_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
           dendconem1_info=((coord_distance(zendconem1_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
           dendcone_info  =((coord_distance(zendcone_info  ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
        else 
           dendconem2_info=(coord_distance(zendconem2_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
           dendconem1_info=(coord_distance(zendconem1_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
           dendcone_info  =(coord_distance(zendcone_info   ,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
        end if
     else
        aendconem2_info=0.0D0
        aendconem1_info=0.0D0
        aendcone_info  =0.0D0
        zendconem2_info=0.0D0
        zendconem1_info=0.0D0
        zendcone_info  =0.0D0
        dendconem2_info=0.0D0
        dendconem1_info=0.0D0
        dendcone_info  =0.0D0
     endif
     !
     call write_infoconegrav(ncell_total,cone_id,cone_zmax,observer_x,observer_y,observer_z,observer_redshift, &
                  &          cone_nlevel,cone_levelmin,cone_levelmax,                                          &
                  &          1.0D0/(z2+1.0D0),1.0D0/(z1+1.0D0),z2,z1,dist2,dist1,disttol,                      &
                  &          reduce_infoglobal,filename_info,                                                  &
                  &          is_fullsky,input_thetay,input_thetaz,input_theta,input_phi,                       &
                  &          aendconem2_info,aendconem1_info,aendcone_info,                                    &
                  &          zendconem2_info,zendconem1_info,zendcone_info,                                    &
                  &          dendconem2_info,dendconem1_info,dendcone_info,future)
  end if
  
  ! Finalize
  deallocate(hdr_ilevel)
  deallocate(hdr_dx)
  deallocate(hdr_ncellperlevel)
 
  ! Baojiu: deallcate these arrays
  deallocate(var_in)
  deallocate(var_out)
  deallocate(var_wr)
  ! Baojiu
 
end subroutine output_conegrav

!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine performs selections of cells for a full-sky lightcone.
! It first finds the minimum cube that encloses the lightcone section in 
! between z1 and z2, by tiling replicas of the simulation box if needed.
!----------------------------------------------------------------------!
subroutine perform_my_selection_simple_conegrav_ramsesunits(onlydist,               &
     &                                justcount,disttol,z1,z2,dist1,dist2,          &
     &                                om0in,omLin,hubin,Lbox,                       &
     &                                observer,observer_redshift,dxcell,nvar,nbool, &
     &                                xc_in,var_in,bool_in,ncell_in,                &
     &                                xc_out,var_out,bool_out,ncell_out,verbose)
  !------------------------------------------------------------------------!  
  ! WARNING: disttol, observer, dxcell and poscell are in RAMSES UNITS [0,1]
  !------------------------------------------------------------------------!  
  use amr_parameters,ONLY:nvector,ndim
  implicit none
  
  ! Variables from input
  logical::onlydist
  logical::justcount
  real(kind=8)::disttol
  real(kind=8)::z1
  real(kind=8)::z2
  real(kind=8)::dist1
  real(kind=8)::dist2
  real(kind=8)::om0in
  real(kind=8)::omLin
  real(kind=8)::hubin
  real(kind=8)::Lbox
  real(kind=8),dimension(1:3)::observer
  real(kind=8)::observer_redshift
  real(kind=8)::dxcell
  integer::nvar
  integer::nbool
  real(kind=8),dimension(1:nvector,1:ndim)::xc_in
  real(kind=8),dimension(1:nvector,1:nvar)::var_in
  integer,dimension(1:nvector,1:nbool)::bool_in
  integer::ncell_in
  real(kind=8),dimension(1:27*nvector,1:ndim)::xc_out
  real(kind=8),dimension(1:27*nvector,1:nvar)::var_out
  integer,dimension(1:27*nvector,1:nbool)::bool_out
  integer::ncell_out
  logical::verbose
  
  ! Variables local
  real(kind=8)::Omega0,OmegaL,OmegaR,coverH0
  real(kind=8)::coord_distance
  real(kind=8)::distmin,distmax
  real(kind=8)::xcoord,ycoord,zcoord
  real(kind=8)::dist,dxtest1,dxtest2,facnorm
  integer::nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
  integer::i,j,k
  integer::myint
  integer::icell_in
  integer::icell_out
  integer::ivar
  integer::ibool
  
  ! Variables for conegrav
  real(kind=8)::demidiag                                                ! demi great diagonal of a cube
  real(kind=8)::demidx                                                  ! dx/2
  
  ! Beginning of function
  if(verbose) write(*,*) 'Enter perform_my_selection_simple'
  
  ! Initialize cosmological parameters
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)

  if(verbose) write(*,*) 'After init_cosmo',Omega0,OmegaL,OmegaR,coverH0
  
  ! Compute comoving distance of the photon planes from the observer
  ! dist1,dist2=integral of c.dt/a between zero and z1,z2
  if(observer_redshift.gt.0.0D0) then
     dist1=((coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
     dist2=((coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)+disttol ! in [0, 1]
  else 
     dist1=(coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
     dist2=(coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)/Lbox)+disttol ! in [0, 1]
  end if
  
  if(verbose) write(*,*) 'After coord',dist1,dist2
  
  if(.not.onlydist) then
     ! Compute the demi diagonal of the cube
     demidx=dxcell/2.0D0                                                ! half cell size
     demidiag=dsqrt(3.0D0)*demidx                                       ! half of the length of cell diagonal

     ! Compute how many replica are needed
     nrepxm=myint(observer(1)-dist2)                                    ! index of leftmost  replica (smallest x)
     nrepxp=myint(observer(1)+dist2)                                    ! index of rightmost replica (largest  x)
     nrepym=myint(observer(2)-dist2)                                    ! index of rear      replica (smallest y)
     nrepyp=myint(observer(2)+dist2)                                    ! index of front     replica (largest  y)
     nrepzm=myint(observer(3)-dist2)                                    ! index of bottom    replica (smallest z)
     nrepzp=myint(observer(3)+dist2)                                    ! index of top       replica (largest  z)

     facnorm=1.0d0/(dist2-dist1)
  
     ncell_out=0                                                        ! total number of selected cells
     icell_out=0                                                        ! initialise cell pointer
     if(verbose) write(*,*) 'before loop',facnorm
     !
     ! loop on all the replica of potential interest
     do k=nrepzm,nrepzp,1
        do j=nrepym,nrepyp,1
           do i=nrepxm,nrepxp,1
              do icell_in=1,ncell_in
                 ! 
                 ! compute relative coordinate of cell centre to observer
                 xcoord=xc_in(icell_in,1)+dble(i)-observer(1)           ! relative x coordinate (in unit of boxsize)
                 ycoord=xc_in(icell_in,2)+dble(j)-observer(2)           ! relative y coordinate (in unit of boxsize)
                 zcoord=xc_in(icell_in,3)+dble(k)-observer(3)           ! relative z coordinate (in unit of boxsize)
                 !
                 dist=dsqrt(xcoord**2+ycoord**2+zcoord**2)              ! distance from cell centre to observer (code unit, or in unit of boxsize)
                 !
                 distmin=dist-demidiag                                  ! closest distance of cell's circumsphere to observer
                 distmax=dist+demidiag                                  ! biggest distance of cell's circumsphere to observer
                 !
                 !-----------------------------------------------------------------------
                 ! select a cell in two situations:
                 ! (1) if the cell's circumsphere intersects the shell enclosed by dist1 and dist2 
                 ! (2) if the cell's centre lies in the shell enclosed by dist1 and dist2
                 !-----------------------------------------------------------------------
                 if((distmax>dist1.and.distmin<=dist2).or.(dist>dist1.and.dist<=dist2)) then
                    ! 
                    icell_out=icell_out+1                               ! increment cell pointer
                    if(icell_out.ge.27*nvector) then 
                       write(*,*) 'PROBLEM : icell_out = ',icell_out, '27*nvector = ', 27*nvector, ' ncell_in = ',ncell_in
                       write(*,*) nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
                    endif
                    !
                    xc_out(icell_out,1)=xcoord                          ! gather cell centre's x coordinate
                    xc_out(icell_out,2)=ycoord                          ! gather cell centre's y coordinate
                    xc_out(icell_out,3)=zcoord                          ! gather cell centre's z coordinate
                    !
                    ! copy other variables
                    do ivar=1,nvar
                       var_out(icell_out,ivar)=var_in(icell_in,ivar)    ! gather cell floating-point quantities
                    end do
                    do ibool=1,nbool
                       bool_out(icell_out,ibool)=bool_in(icell_in,ibool)! gather cell integer quantities
                    end do

                 endif
              enddo
           enddo
        enddo
     enddo
     !
     ncell_out = icell_out                                              ! total number of cells selected
  endif
  ! 
  if (verbose) write(*,*) 'End of perform_my_selection_simple',ncell_out

end subroutine perform_my_selection_simple_conegrav_ramsesunits


!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! Perform the selection of cells for narrow cone described by an opening
! field of view (input_thetay,input_thetaz) and a direction with respect
! to the x-axis (input_theta,input_phi). All cells that either fall into
! the linecone section entirely, or intersect with the lightcone section
! will be selected. Note that all vectors are rotated so that the centre
! line of the field-of-view coincides with the x-axis after the rotation
!----------------------------------------------------------------------!   
subroutine perform_my_selection_conegrav_ramsesunits(onlydist,                      &
     &                                justcount,disttol,z1,z2,dist1,dist2,          &
     &                                om0in,omLin,hubin,Lbox,                       &
     &                                observer,observer_redshift,dxcell,nvar,nbool, &
     &                                xc_in,var_in,bool_in,ncell_in,                &
     &                                xc_out,var_out,bool_out,ncell_out,verbose,    &
     &                                input_thetay,input_thetaz,input_theta,input_phi)
  !===========================================================================
  ! WARNING : disttol, observer, dxcell and poscell are in RAMSES UNITS [0,1]
  !===========================================================================
  use amr_parameters, ONLY: nvector, ndim
  implicit none
  
  ! Variables from input
  logical::onlydist
  logical::justcount
  real(kind=8)::disttol
  real(kind=8)::z1
  real(kind=8)::z2
  real(kind=8)::dist1
  real(kind=8)::dist2
  real(kind=8)::om0in
  real(kind=8)::omLin
  real(kind=8)::hubin
  real(kind=8)::Lbox
  real(kind=8),dimension(1:3)::observer,observerLbox
  real(kind=8)::observer_redshift
  real(kind=8)::dxcell
  integer::nvar
  integer::nbool
  real(kind=8),dimension(1:nvector,1:ndim)::xc_in
  real(kind=8),dimension(1:nvector,1:nvar)::var_in
  integer,dimension(1:nvector,1:nbool)::bool_in
  integer::ncell_in
  real(kind=8),dimension(1:27*nvector,1:ndim)::xc_out
  real(kind=8),dimension(1:27*nvector,1:nvar)::var_out
  integer,dimension(1:27*nvector,1:nbool)::bool_out
  integer::ncell_out
  logical::verbose
  real(kind=8)::input_thetay
  real(kind=8)::input_thetaz
  real(kind=8)::input_theta
  real(kind=8)::input_phi
  real(kind=8)::small
  
  ! Variables local
  real(kind=8)::Omega0,OmegaL,OmegaR,coverH0
  real(kind=8)::coord_distance
  real(kind=8)::distmin,distmax
  real(kind=8)::dist
  integer::nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
  integer::i,j,k
  integer::myint
  integer::icell_in
  integer::icell_out
  integer::ivar
  integer::ibool
  real(kind=8)::pi,thetarad,phirad,thetayrad,thetazrad,tanybound,tanzbound
  real(kind=8)::dist1Lbox,dist2Lbox,tany,tanz
  real(kind=8),dimension(1:9)::xcoordfr,ycoordfr,zcoordfr,xcoord,ycoord,zcoord,dxcoord,dycoord,dzcoord
  real(kind=8),dimension(1:3, 1:3)::rot,rotm1
  integer::ivertex
  logical::ok
  
  ! Variables for conegrav
  real(kind=8)::demidiag                                                ! demi great diagonal of a cube
  real(kind=8)::demidx                                                  ! dx/2
  
  ! Beginning of function
  small=1.0D-5
  if(verbose) write(*,*) 'Enter perform_my_selection_conegrav_ramsesunits'
  
  ! Initialize cosmological parameters
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  
  ! Deal with angles
  pi = dacos(-1.0D0)
  thetarad  = input_theta *pi/180.0D0
  phirad    = input_phi   *pi/180.0D0
  thetayrad = input_thetay*pi/180.0D0
  thetazrad = input_thetaz*pi/180.0D0
  !
  call compute_rotation_matrix_grav(thetarad,phirad,rot,rotm1)
  
  if(verbose) write(*,*) 'After init_cosmo',Omega0,OmegaL,OmegaR,coverH0
  
  ! Compute comoving distance of the photon planes from the observer : dist1,dist2=integral of c.dt/a between zero and z1,z2
  if(observer_redshift.gt.0) then
     dist1=((coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
     dist2=((coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)+disttol ! in [0, 1]
  else 
     dist1=(coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
     dist2=(coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)/Lbox)+disttol ! in [0, 1]
  end if
  !
  dist1Lbox = dist1*Lbox                                                ! distance to the inner side of the shell in Mpc
  dist2Lbox = dist2*Lbox                                                ! distance to the outer side of the shell in Mpc
  observerLbox=observer*Lbox                                            ! observer coordinates in Mpc

  if(verbose) write(*,*) 'After coord',dist1,dist2
  
  ! Actual selection of cells 
  if(.not.onlydist) then
     ! Compute the min/max indices of the set of replicas to be considered
     call compute_replica_grav(thetayrad,thetazrad,dist1Lbox,dist2Lbox,observerLbox,Lbox,rot, &
             &                 nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp)    
      
     tanybound=dtan(thetayrad)
     tanzbound=dtan(thetazrad)
  
     ! Compute specific points of the cell
     demidx=dxcell/2.D0                                                 ! half of cell side size
     demidiag=dsqrt(3.D0)*demidx                                        ! half of the size of cell diagonal
     ! Distances of special points of the cell to cell centre, with the 
     ! following order: first the cell centre itself, and then followed
     ! by (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1)
     ! and (0,1,1), where "0" and "1" respectively mean the smaller and 
     ! larger of of the x, y and z coordinates.
     dxcoord=(/0.0D0,-demidx, demidx,-demidx, demidx,-demidx, demidx,-demidx,demidx/) ! x direction
     dycoord=(/0.0D0,-demidx,-demidx, demidx, demidx,-demidx,-demidx, demidx,demidx/) ! y direction 
     dzcoord=(/0.0D0,-demidx,-demidx,-demidx,-demidx, demidx, demidx, demidx,demidx/) ! z direction
     
     ! Loop over cells
     ncell_out=0   
     icell_out=0
     if(verbose) write(*,*) 'before loop'   
     do k=nrepzm,nrepzp,1
        do j=nrepym,nrepyp,1
           do i=nrepxm,nrepxp,1
              do icell_in=1,ncell_in
                 ! Distances of special points of the cell to observer (code unit)
                 xcoordfr=xc_in(icell_in,1)+dble(i)-observer(1)+dxcoord ! x distance
                 ycoordfr=xc_in(icell_in,2)+dble(j)-observer(2)+dycoord ! y distance
                 zcoordfr=xc_in(icell_in,3)+dble(k)-observer(3)+dzcoord ! z distance
                 !
                 ! Rotation to get in the framework of the photon plane
                 xcoord=xcoordfr*rotm1(1,1)+ &
                      & ycoordfr*rotm1(2,1)+ &
                      & zcoordfr*rotm1(3,1)                             ! x distance after rotation
                 ycoord=xcoordfr*rotm1(1,2)+ & 
                      & ycoordfr*rotm1(2,2)+ &
                      & zcoordfr*rotm1(3,2)                             ! y distance after rotation
                 zcoord=xcoordfr*rotm1(1,3)+ &
                      & ycoordfr*rotm1(2,3)+ &
                      & zcoordfr*rotm1(3,3)                             ! z distance after rotation
                 !     
                 ! Loop over special points (vertexes)
                 ok=.false.
                 ivertex=1
                 do while((.not.ok).and.(ivertex.le.9))                 ! loop over 9 special points of cell
                    if(xcoord(ivertex).gt.small) then                   ! only select cells not too close to observer
                       tany=dabs(ycoord(ivertex)/xcoord(ivertex))
                       tanz=dabs(zcoord(ivertex)/xcoord(ivertex))
                       !
                       dist=dsqrt(xcoord(1)**2+ &                       ! distance of cell centre to observer
                           &      ycoord(1)**2+ &
                           &      zcoord(1)**2)                         ! RY made it compatible with full cone
                       !
                       distmin=dist-demidiag                            ! min distance of cell's circumsphere to observer
                       distmax=dist+demidiag                            ! max distance of cell's circumsphere to observer
                       ! Select a cell if it satisfies the following conditions:
                       ! (1) its centre lies in the lightcone field of view, and
                       ! (2) its circumsphere intersects the lightcone shell, or
                       ! (3) its cell centre lies within the lightcone shell.
                       ! Note that (2) and (3) are not necessarily exclusive.
                       ok = (((tany.le.tanybound).and.(tanz.le.tanzbound)).and. &
                          &  ((distmax>dist1.and.distmin<=dist2).or.(dist>dist1.and.dist<=dist2)))
                    end if
                    ivertex=ivertex+1                                   ! increment pointer to special points
                 end do
 
                 ! If one vertex is ok, add it to the list
                 if(ok) then
                    icell_out=icell_out+1                               ! increment pointer to var_out
                    if(icell_out.ge.27*nvector) then 
                       write(*,*) 'PROBLEM : icell_out = ',icell_out, '27*nvector = ', 27*nvector, ' ncell_in = ',ncell_in
                       write(*,*) nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
                    endif
                    !
                    xc_out(icell_out,1)=xcoord(1)                       ! gather x coordinate of cell centre
                    xc_out(icell_out,2)=ycoord(1)                       ! gather y coordinate of cell centre
                    xc_out(icell_out,3)=zcoord(1)                       ! gather z coordinate of cell centre
                    !
                    ! Baojiu: modify the lines below for Newtonian vs GR runs
                    if(.not.gr) then 
                    ! Rotate and copy the force
                       var_out(icell_out,1)=var_in(icell_in,1)*rotm1(1,1)+ &
                            &               var_in(icell_in,2)*rotm1(2,1)+ &
                            &               var_in(icell_in,3)*rotm1(3,1)
                       var_out(icell_out,2)=var_in(icell_in,1)*rotm1(1,2)+ &
                            &               var_in(icell_in,2)*rotm1(2,2)+ &
                            &               var_in(icell_in,3)*rotm1(3,2)
                       var_out(icell_out,3)=var_in(icell_in,1)*rotm1(1,3)+ &
                            &               var_in(icell_in,2)*rotm1(2,3)+ &
                            &               var_in(icell_in,3)*rotm1(3,3)
   
                       ! Copy other variables
                       do ivar=4,nvar
                          var_out(icell_out,ivar)=var_in(icell_in,ivar)
                       end do
                    else 
                       ! Baojiu: these need to complete. Notice that scalar
                       ! quantities can be copied directly; however, vector
                       ! quantities need to be rotated before being copied.
                       var_out(icell_out,1 )=
                       var_out(icell_out,2 )=
                       var_out(icell_out,3 )=
                       var_out(icell_out,4 )=
                       var_out(icell_out,5 )=
                       var_out(icell_out,6 )=
                       var_out(icell_out,7 )=
                       var_out(icell_out,8 )=
                       var_out(icell_out,9 )=
                       var_out(icell_out,10)=
                       var_out(icell_out,11)=
                       var_out(icell_out,12)=
                       var_out(icell_out,13)=
                       var_out(icell_out,14)=
                       var_out(icell_out,15)=
                       var_out(icell_out,16)=
                       var_out(icell_out,17)=
                       var_out(icell_out,18)=
                       var_out(icell_out,19)=
                       var_out(icell_out,20)=
                    endif
                    !
                    do ibool=1,nbool
                       bool_out(icell_out,ibool)=bool_in(icell_in,ibool)
                    end do
                 endif
              enddo
           enddo
        enddo
     enddo
     ncell_out=icell_out                                                ! total number of cells selected in this call
  endif
  
  if(verbose) write(*,*) 'End of perform_my_selection',ncell_out

end subroutine perform_my_selection_conegrav_ramsesunits



!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine select cells for simple sample -- it is not a full-sky
! or a narrow lightcone, but a cubic region from the simulaton box whose
! minimum and maximum x, y and z coordinates are specificed as arguement
! of this subroutine. 
!----------------------------------------------------------------------!   
subroutine extract_samplegrav(filedir,filename,xmin,xmax,ymin,ymax,zmin,zmax, &
            &                 nsample,input_levelmin,input_levelmax)
  use amr_commons
  use poisson_commons
  use gr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none

  ! Misc
  integer,parameter::tag=2900
  integer::dummy_io
  integer::info
  integer::ilun
  logical::is_opened

  ! Input
  character(LEN=*)::filedir                                             ! name of the output dir
  character(LEN=*)::filename                                            ! name of the output file
  real(dp)::xmin
  real(dp)::xmax
  real(dp)::ymin
  real(dp)::ymax
  real(dp)::zmin
  real(dp)::zmax
  integer::nsample
  integer::input_levelmin
  integer::input_levelmax
  
  ! Variables : level
  integer::sample_levelmin
  integer::sample_levelmax
  integer::sample_nlevel

  ! Variables
  real(kind=8)::dx
  real(kind=8)::dxcell
  real(kind=8)::dxon2
  
  ! Variables : filename
  character(LEN=200)::filename_info
  character(LEN=200)::filename_dat
  character(LEN=200)::filename_hdr
  character(LEN=5)::nchar_cpu
  character(LEN=5)::nchar_coarse
  
  ! Variables : index
  integer::ilevel
  integer::ncache
  integer::idim
  integer::iskip
  
  ! Variables : cells centers
  integer::ind
  integer::ix
  integer::iy
  integer::iz
  real(kind=8),dimension(1:twotondim,1:3)::xc
  
  ! Variables : common for arrays
! integer::nvar = 3                                                     ! fx, fy, fz
! integer,parameter::nvar=5                                             ! fx, fy, fz, phi, rho
                                                                        ! 1    3          1    3          3       9 
! integer,parameter::nvar=20                                            ! Phi, grad(Phi), Psi, grad(Psi), beta^i, grad(beta^i) ! CBH_LC
  integer::nvar                                                         ! Baojiu: changed nvar to a normal variable to be given value later
  integer,parameter::nbool=1 ! son
  integer,parameter::nvectorgrid = floor(nvector/(2D0.**ndim))          ! nvector/8
  integer,dimension(1:nvectorgrid)::ind_grid                            ! indices of grids
  integer,dimension(1:nvector)::ind_cell                                ! indices of cells
  integer::ivar
  integer::ibool
  
  ! Variables : input arrays
  real(kind=8),dimension(1:nvector,1:ndim)::xc_in
! real(kind=8),dimension(1:nvector,1:nvar)::var_in
  real(kind=8),dimension(:,:),allocatable::var_in                       ! Baojiu: changed to allocatable
  integer,dimension(1:nvector,1:nbool)::bool_in

  ! Variables : selection arrays
  real(kind=8),dimension(1:27*nvector,1:ndim)::xc_out
! real(kind=8),dimension(1:27*nvector,1:nvar)::var_out
  real(kind=8),dimension(:,:),allocatable::var_out                      ! Baojiu: changed to allocatable
  integer,dimension(1:27*nvector,1:nbool)::bool_out
  
  ! Variables : write arrays
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim)::xc_wr            ! x, y, z centers of cells for output
! real(kind=4),dimension(1:nstride+27*nvector,1:nvar)::var_wr           ! gravity vars
  real(kind=4),dimension(:,:),allocatable::var_wr                       ! Baojiu: changed to allocatable
  integer,dimension(1:nstride+27*nvector,1:nbool)::bool_wr              ! son
 
  ! Hdr arrays
  integer::ihdr
  integer::hdr_arraysize
  integer,dimension(:),allocatable::hdr_ilevel
  real(kind=8),dimension(:),allocatable::hdr_dx
  integer,dimension(:),allocatable::hdr_ncellperlevel

  ! Variables : counters
  integer::ncell_level, ncell_total
  integer::igrid, ngrid, i
  integer::ivector, istride
  integer::icell_in, icell_out, icell_wr
  integer::ncell_in, ncell_out, ncell_wr

  ! Write variables
  integer(kind=8),dimension(1:2)::reduce_infolocal
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  integer::ierr

  ! Iogroupsize variables
  real(kind=8)::iovolume

  ierr=0

  ! Baojiu: allocate arrays
  if(.not.gr) then
     nvar=5
  else
     nvar=20
  endif
  allocate(var_in (1:           nvector,1:nvar))
  allocate(var_out(1:        27*nvector,1:nvar))
  allocate(var_wr (1:nstride+27*nvector,1:nvar))
  ! Baojiu: end of changes
  
  ! Correct quantities
  sample_levelmin=input_levelmin                                        ! minimum level to output cell fields
  sample_levelmax=input_levelmax                                        ! maximum level to output cell fields
  if((sample_levelmin.eq.0).or.(sample_levelmin.gt.nlevelmax)) then
     sample_levelmin = levelmin                                         ! should't be lower than the domain level
  endif
  if((sample_levelmax.eq.0).or.(sample_levelmax.gt.nlevelmax)) then
     sample_levelmax = nlevelmax                                        ! shouldn't be higher than finest AMR level
  endif
  sample_nlevel=sample_levelmax-sample_levelmin+1                       ! number of levels to outut cell fields

  ! Filename
  ilun=30*ncpu+myid+10                                                  ! file handler
  call title(myid,nchar_cpu)                                            ! create file name
  call title(nstep_coarse,nchar_coarse)                                 ! create file name
  filename_info=TRIM(filedir)//'info_sample_grav_ncoarse_'//TRIM(nchar_coarse)//'.txt'
  filename_dat =TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.dat'
  filename_hdr =TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.hdr'

  ! Hdr arrays
  hdr_arraysize=sample_nlevel
  allocate(hdr_ilevel(1:hdr_arraysize))
  allocate(hdr_dx(1:hdr_arraysize))
  allocate(hdr_ncellperlevel(1:hdr_arraysize))
  hdr_ilevel=0
  hdr_dx=0
  hdr_ncellperlevel=0

  iovolume=(dabs(xmax-xmin)*dabs(ymax-ymin)*dabs(zmax-zmin))
  if(adaptive_iogroupsize) IOGROUPSIZESAMPLE=nint(iovolume*real(IOGROUPSIZE,kind=8))
  if(IOGROUPSIZESAMPLE.le.1) IOGROUPSIZESAMPLE=1
  if(IOGROUPSIZESAMPLE.ge.IOGROUPSIZE) IOGROUPSIZESAMPLE=IOGROUPSIZE
  if(verbose) write(*,*) 'IOGROUPSIZESAMPLE = ',IOGROUPSIZESAMPLE, 'iovolume = ', iovolume

  ! Wait for the token
  if(sampletoken) then 
     if(mod(myid-1,IOGROUPSIZESAMPLE).ne.0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  endif
  
  ! Loop over levels
  ncell_total=0                                                         ! total number cells selected
  ihdr=0
  reduce_infolocal=0
  reduce_infoglobal=0
  is_opened=.FALSE.

  do ilevel=sample_levelmin,sample_levelmax
     ! Initialization of counters
     ncell_level=0                                                      ! total number of cells on a given leve
     ivector=0
     istride=0
     icell_in=0                                                         ! pointer to array for cells sent for selection
     icell_out=0                                                        ! pointer to array for cells selected
     icell_wr=0                                                         ! pointer to array for cells to be written
     ncell_in=0                                                         ! number of cells sent for selection
     ncell_out=0                                                        ! number of cells selected
     ncell_wr=0                                                         ! number of cells to be written
    
     ! Initialization of arrays
     xc_in  =0.D0                                                       ! array for coordinates of cells sent for selection
     xc_out =0.D0                                                       ! array for coordinates of cells selected
     xc_wr  =0.D0                                                       ! array for coordinates of cells to be written
     var_in =0.D0                                                       ! array for field-quantities in cells sent for selection
     var_out=0.D0                                                       ! array for field-quantities in cells selected
     var_wr =0.D0                                                       ! array for field quantities in cells to be written
     bool_in =0                                                         ! array for integer variables in cells send for selection
     bool_out=0                                                         ! array for integer variables in cells selected
     bool_wr =0                                                         ! array for integer variables in cells to be written
    
     ! Index
     ihdr=ihdr+1
    
     ! Level constants
     dx=0.5D0**ilevel                                                   ! cell size in code unit
     dxcell=dx                                                          ! need to check (maybe dxcell = (1/2.)*dx
     dxon2=0.5D0*dx                                                     ! half cell size

     ! Set position of cell centers relative to grid center
     do ind=1,twotondim
        iz=(ind-1          )/4
        iy=(ind-1     -4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0) xc(ind,1)=(dble(ix)-0.5D0)*dx                        ! relative x coordinate of cell centre vs greid centre
        if(ndim>1) xc(ind,2)=(dble(iy)-0.5D0)*dx                        ! relative y coordinate of cell centre vs greid centre
        if(ndim>2) xc(ind,3)=(dble(iz)-0.5D0)*dx                        ! relative z coordinate of cell centre vs greid centre
     end do
    
     ! Loop over myid grids by vectorgrid sweeps
     ncache=active(ilevel)%ngrid                                        ! total number of grids on myid
     do igrid=1,ncache,nvectorgrid
        ngrid=MIN(nvectorgrid,ncache-igrid+1)                           ! batch size of grid in current sweep
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)                  ! global grid indices
!          ind_grid(i)=active(ilevel)%pcomm%igrid(igrid+i-1)            ! RAMSES-LC line, commenbted out CBH_LC 16-02-2021
        end do
      
        ! Loop over cells
        ivector = 0                                                     ! initialise pointer to array for cells to be sent for selection
        do i=1,ngrid
           do ind=1,twotondim
              ! 
              ivector=ivector+1                                         ! increment ivector (or counter)
              !
              ! Gather cell indices
              iskip=ncoarse+(ind-1)*ngridmax
              ind_cell(ivector)=iskip+ind_grid(i)                       ! cell indices
              ! 
              ! Collect cell centre coordinates
              do idim=1,ndim
                 xc_in(ivector,idim)=xg(ind_grid(i),idim)+xc(ind,idim)  ! cell centre coordinates
              end do
              !
              ! Collect field quantities
              ! Baojiu: modified to deal with both Newtonian and GR runs
              if(.not.gr) then
                 var_in(ivector,1)=f(ind_cell(ivector),1)               ! gather x component of force
                 var_in(ivector,2)=f(ind_cell(ivector),2)               ! gather y component of force
                 var_in(ivector,3)=f(ind_cell(ivector),3)               ! gather z component of force
                 var_in(ivector,4)=phi(ind_cell(ivector))               ! gather potential
                 var_in(ivector,5)=rho(ind_cell(ivector))               ! gather density
              else
                 ! CBH_LC
                 ! Notice that gradients need to be copied from f() at the correct igrp sequence in move or synchro
                 ! At igrp = 10, f() = - grad_i(b)
                 var_in(ivector,1 )=gr_pot(ind_cell(ivector),5)         ! Phi
                 var_in(ivector,2 )=gr_pot(ind_cell(ivector),6)         ! Xi 
                 var_in(ivector,3 )=gr_pot(ind_cell(ivector),7)-f(ind_cell(ivector),1)  ! \beta^i = B^i + grad_i(b)
                 var_in(ivector,4 )=gr_pot(ind_cell(ivector),8)-f(ind_cell(ivector),2)
                 var_in(ivector,5 )=gr_pot(ind_cell(ivector),9)-f(ind_cell(ivector),3)
                 var_in(ivector,6 )=gr_pot_grad(ind_cell(ivector),1 )   ! grad_i Phi
                 var_in(ivector,7 )=gr_pot_grad(ind_cell(ivector),2 )
                 var_in(ivector,8 )=gr_pot_grad(ind_cell(ivector),3 )
                 var_in(ivector,9 )=gr_pot_grad(ind_cell(ivector),4 )   ! grad_i Xi
                 var_in(ivector,10)=gr_pot_grad(ind_cell(ivector),5 )
                 var_in(ivector,11)=gr_pot_grad(ind_cell(ivector),6 )
                 var_in(ivector,12)=gr_pot_grad(ind_cell(ivector),7 )   ! grad_i \beta^j
                 var_in(ivector,13)=gr_pot_grad(ind_cell(ivector),8 )
                 var_in(ivector,14)=gr_pot_grad(ind_cell(ivector),9 )
                 var_in(ivector,15)=gr_pot_grad(ind_cell(ivector),10)
                 var_in(ivector,16)=gr_pot_grad(ind_cell(ivector),11)
                 var_in(ivector,17)=gr_pot_grad(ind_cell(ivector),12)
                 var_in(ivector,18)=gr_pot_grad(ind_cell(ivector),13)
                 var_in(ivector,19)=gr_pot_grad(ind_cell(ivector),14)
                 var_in(ivector,20)=gr_pot_grad(ind_cell(ivector),15)
                 ! END CBH_LC
              endif
              ! Baojiu: end of modification
              !
              ! Collect integer quantities
              bool_in(ivector,1) = son(ind_cell(ivector))               ! gather son grid indices
              !
           end do 
        end do                                                          ! end of loop over cells (do i=1,ngrid)
                                                                        ! ivector cells have been gathered      
        ! Perfom selection on cells
        ncell_in=ivector                                                ! number of cells sent for selection
        icell_out=0                                                     ! initialise pointer to array for selected cells
        do icell_in=1,ncell_in
           ! This is a simple sample selection, so only cells that
           ! intersect with the desired cubic region are selected. 
           if((((xc_in(icell_in,1)+dxon2).gt.xmin).and.((xc_in(icell_in,1)-dxon2).lt.xmax)) .and. &
            & (((xc_in(icell_in,2)+dxon2).gt.ymin).and.((xc_in(icell_in,2)-dxon2).lt.ymax)) .and. &
            & (((xc_in(icell_in,3)+dxon2).gt.zmin).and.((xc_in(icell_in,3)-dxon2).lt.zmax))) then
              !
              icell_out=icell_out+1                                     ! increase pointer/counter for selected cells
              !
              do idim=1,ndim
                 xc_out(icell_out,idim)=xc_in(icell_in,idim)            ! gather cell centre coordinates
              end do
              !
              do ivar=1,nvar
                 var_out(icell_out,ivar)=var_in(icell_in,ivar)          ! gather field quantities in cells
              end do
              !
              do ibool=1,nbool
                 bool_out(icell_out,ibool)=bool_in(icell_in,ibool)      ! gather son grid indices
              end do
           end if
        end do
        !
        ncell_out = icell_out                                           ! total number of selected cells
        !   
        ! Transfer out arrays of size ncell_out in write arrays of size icell_wr
        ! Note that icell_wr is not reset to 0 between consecutive calls to this
        ! function
        if(ncell_out.gt.0) then
           call transfer_grav(nvar,nbool,ncell_out,icell_wr,ncell_level,ncell_total, &
                            & xc_out,var_out,bool_out,xc_wr,var_wr,bool_wr)
        end if
        !
        ! If nstride has been reached, write arrays of size icell_wr and update size
        if(icell_wr.ge.nstride) then
           call write_grav(.false.,is_opened,ilun,filename_dat,       &
                            & xc_wr,nvar,var_wr,nbool,bool_wr,        &
                            & icell_wr,sample_nlevel,reduce_infolocal)
        end if
        !
     end do ! loop over grids
    
     ! If some cells remain in write vector, write them
     if(icell_wr.gt.0) then
        call write_grav(.true.,is_opened,ilun,filename_dat,           &
                      & xc_wr,nvar,var_wr,nbool,bool_wr,              &
                      & icell_wr,sample_nlevel,reduce_infolocal)
     end if

     hdr_ilevel(ihdr)=ilevel
     hdr_dx(ihdr)=dx
     hdr_ncellperlevel(ihdr)=ncell_level

  end do ! loop over levels

  if(is_opened) close(ilun)                                             ! close file after writing
  
  ! Write hdr files
  call write_hdr(ncell_total,filename_hdr,hdr_arraysize,hdr_ilevel,hdr_dx,hdr_ncellperlevel)
  
  ! Send the token
  if(sampletoken) then 
     if(mod(myid,IOGROUPSIZESAMPLE).ne.0.and.(myid.lt.ncpu)) then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag,MPI_COMM_WORLD,info)
     end if
  endif

  ! Reduce
  reduce_infolocal(2)=ncell_total
  call MPI_REDUCE(reduce_infolocal,reduce_infoglobal,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
  ! Write info file
  if(myid.eq.1) then
     call write_infosamplegrav(ncell_total,aexp,sample_nlevel,sample_levelmin,sample_levelmax, &
                             & xmin,xmax,ymin,ymax,zmin,zmax,reduce_infoglobal,filename_info)
  end if
  
  ! Finalize
  deallocate(hdr_ilevel)
  deallocate(hdr_dx)
  deallocate(hdr_ncellperlevel)
   
  ! Baojiu: deallocate arrays
  deallocate(var_in ) 
  deallocate(var_out)
  deallocate(var_wr )
  ! Baojiu: end of modification 
  
end subroutine extract_samplegrav



!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine transfer data for selected cells from the *_out arrays
! to the corresponding *_wr arrays. The former stores data of cells that
! have got selected, while the latter stores data to be written. This is
! needed because we do not write the data immediately when some cell has
! been selected, but the data is accumulated until its size exceeds some
! predefined number (nstride) -- this means that the data must be stored
! somewhere temporarily until nstride is reached, so that arrays such as
! xc_out, var_out and bool_out can be released for selecting more cells.
!----------------------------------------------------------------------!   
subroutine transfer_grav(nvar,nbool,ncell_out,icell_wr,ncell_level,ncell_total, &
                       & xc_out,var_out,bool_out,xc_wr,var_wr,bool_wr)

  ! Initialization
  use amr_commons, ONLY: ndim, nstride, nvector, ncpu
  implicit none
  
  ! Variables from input
  integer::nvar
  integer::nbool
  integer::ncell_out
  integer::icell_wr
  integer::ncell_level
  integer::ncell_total
  !
  real(kind=8),dimension(1:27*nvector,1:ndim) ::xc_out
  real(kind=8),dimension(1:27*nvector,1:nvar) ::var_out
  integer     ,dimension(1:27*nvector,1:nbool)::bool_out
  !
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim) ::xc_wr
  real(kind=4),dimension(1:nstride+27*nvector,1:nvar) ::var_wr
  integer     ,dimension(1:nstride+27*nvector,1:nbool)::bool_wr
  
  ! local variables
  integer::icell
  integer::idim
  integer::ivar
  integer::ibool
  
  ! Copy loops
  do icell=1,ncell_out
    do idim=1,ndim
      xc_wr(icell_wr+icell,idim)=xc_out(icell,idim)                     ! cell centre coordinates
    end do
    do ivar=1,nvar
      var_wr(icell_wr+icell,ivar)=var_out(icell,ivar)                   ! cell floating-point variables
    end do
    do ibool=1,nbool
      bool_wr(icell_wr+icell,ibool)=bool_out(icell,ibool)               ! cell integer variables
    end do
  end do
  !
  icell_wr=icell_wr+ncell_out                                           ! update cell pointer
  ncell_level=ncell_level+ncell_out                                     ! increment total cells written for the given level
  ncell_total=ncell_total+ncell_out                                     ! increment total cells written
                        
end subroutine transfer_grav



!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine writes desired information of cells into a binary file
! after the size of the arrays xc_wr, var_wr and bool_wr exceed nstride.
! Note that it only writes a data size of nstride and the remaining data
! is not written but  moved to the beginning of the corresponding array.
! If the code is told that this is the last batch (if is_last.eq..true.) 
! of data to be written, then it writes the data even the array size has
! not reached nstride.
!----------------------------------------------------------------------!   
subroutine write_grav(is_last,is_opened,ilun,filename,          &
                    & xc_wr,nvar,var_wr,nbool,bool_wr,          &
                    & istride,cone_nlevel,reduce_infolocal)
                        
  use amr_commons,ONLY:ndim,nstride,nvector,ncpu
  implicit none
  
  ! Variables from input
  logical::is_last                                                      ! last writing
  logical::is_opened                                                    ! if the file is already opened
  integer::ilun                                                         ! unit for writing
  character(LEN=200)::filename                                          ! file location
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim)::xc_wr            ! center of cells array
  integer::nvar                                                         ! number of variables in var array
  real(kind=4),dimension(1:nstride+27*nvector,1:nvar)::var_wr           ! variable array
  integer::nbool                                                        ! number of variables in bool array
  integer,dimension(1:nstride+27*nvector,1:nbool)::bool_wr              ! boolean array
  integer::istride
  integer::cone_nlevel
  integer(kind=8),dimension(1:2)::reduce_infolocal
  
  ! Other variables
  integer::idim=0
  integer::ivar=0
  integer::ibool=0
  integer::wrsize=nstride
  integer::i=0
  
  ! Case for writing size
  if(is_last) then
    wrsize = istride                                                    ! last time for writing, nstride may not be reached
  else 
    wrsize = nstride                                                    ! not last time for writing, nstride must be reached
  endif
  
  ! Write the beginning of the file
  if(.not.is_opened) then                                               ! open file for writing if not yet
    open(ilun,file=TRIM(filename),form='unformatted')
    rewind(ilun)  
    write(ilun) ncpu
    write(ilun) nstride
    write(ilun) cone_nlevel
    is_opened=.TRUE.
  endif
 
  ! Baojiu: I find the order of writing pretty weird; also the lines 
  ! below are modified to work for both Newtonian and GR simulations
  !
  ! Write content
  if(.not.gr) then
     write(ilun) xc_wr  (1:wrsize,1)                                    ! write cell centre x coordinates
     write(ilun) var_wr (1:wrsize,1)                                    ! write cell centre x forces
     write(ilun) xc_wr  (1:wrsize,2)                                    ! write cell centre y coordinates
     write(ilun) var_wr (1:wrsize,2)                                    ! write cell centre y forces
     write(ilun) xc_wr  (1:wrsize,3)                                    ! write cell centre z coordinates
     write(ilun) var_wr (1:wrsize,3)                                    ! write cell centre z forces
     write(ilun) var_wr (1:wrsize,4)                                    ! write cell centre desities
     write(ilun) var_wr (1:wrsize,5)                                    ! write cell centre potentials
     write(ilun) bool_wr(1:wrsize,1)                                    ! write cell son grid indcies
  else
     do idim=1,ndim
        write(ilun) xc_wr(1:wrsize,idim)                                ! write cell centre x coordinates
     end do
     !
     do ivar=1,nvar
        write(ilun) var_wr(1:wrsize,ivar)                               ! write cell centre         
     end do
     !
     do ibool=1,nbool
        write(ilun) bool_wr(1:wrsize,ibool)                             ! write cell centre         
     end do
  endif

  ! Copy content
  if(.not.is_last) then
     do i=1,istride-nstride
        do idim=1,ndim
           xc_wr(i,idim)=xc_wr(i+nstride,idim)                          ! move unwritten quantities to the start of array
        end do
        do ivar=1,nvar
           var_wr(i,ivar)=var_wr(i+nstride,ivar)                        ! move unwritten quantities to the start of array
        end do
        do ibool=1,nbool
           bool_wr(i,ibool)=bool_wr(i+nstride,ibool)                    ! move unwritten quantities to the start of array
        end do
     end do
     istride=istride-nstride                                            ! reset write-array pointer
  endif
  
  if(reduce_infolocal(1).eq.0) then
     reduce_infolocal(1) = 1
  endif

end subroutine write_grav


!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine writes useful information of header into a file.
!----------------------------------------------------------------------!   
subroutine write_hdr(ncell_total,filename,hdr_arraysize,hdr_ilevel,hdr_dx,hdr_ncellperlevel)
  use amr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  
  ! Variables from input
  integer::ncell_total
  integer::hdr_arraysize
  integer,dimension(1:hdr_arraysize)::hdr_ilevel
  real(kind=8),dimension(1:hdr_arraysize)::hdr_dx
  integer,dimension(1:hdr_arraysize)::hdr_ncellperlevel
  character(LEN=*)::filename
  
  ! Local variables
  integer::ilun
  character(LEN=200)::filename_tmp
  
  ! Write hdr
  ilun=40*ncpu+myid+10
  if(ncell_total>0) then
     filename_tmp=filename
     open(ilun,file=TRIM(filename_tmp),form='unformatted')
     rewind(ilun)
     write(ilun) ncpu
     write(ilun) nstride
     write(ilun) ncell_total
     write(ilun) hdr_ncellperlevel
     close(ilun)
     if(conegrav_formattedhdr) then
        filename_tmp=TRIM(filename)//'_formatted'
        open(ilun,file=TRIM(filename_tmp),form='formatted')
        rewind(ilun)
        write(ilun,*) ncpu
        write(ilun,*) nstride
        write(ilun,*) ncell_total
        write(ilun,*) hdr_ncellperlevel
        close(ilun)
     endif
  endif
  
end subroutine write_hdr



!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine writes useful information of cone section into a file.
!----------------------------------------------------------------------!   
subroutine write_infoconegrav(ncell_tot,conegravid,conegravzlim,observer_x,observer_y,observer_z,observer_redshift, &
                            & nlevel,nlevel_min,nlevel_max,                                                         &
                            & amax,amin,zmax,zmin,dmax,dmin,dtol,                                                   &
                            & reduce_infoglobal,filename,                                                           &
                            & is_fullsky,input_thetay,input_thetaz,input_theta,input_phi,                           &
                            & aendconem2_info,aendconem1_info,aendcone_info,                                        &
                            & zendconem2_info,zendconem1_info,zendcone_info,                                        &
                            & dendconem2_info,dendconem1_info,dendcone_info,future)
  use amr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  
  ! Variables from input
  integer::ncell_tot
  integer::conegravid
  real(kind=8)::conegravzlim
  real(kind=8)::observer_x
  real(kind=8)::observer_y
  real(kind=8)::observer_z
  real(kind=8)::observer_redshift
  real(kind=8)::amax
  real(kind=8)::amin
  real(kind=8)::zmax
  real(kind=8)::zmin
  real(kind=8)::dmax
  real(kind=8)::dmin
  real(kind=8)::dtol
  integer::nlevel
  integer::nlevel_min
  integer::nlevel_max
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  character(LEN=200)::filename
  logical::is_fullsky
  real(kind=8)::input_thetay
  real(kind=8)::input_thetaz
  real(kind=8)::input_theta
  real(kind=8)::input_phi
  real(kind=8)::aendconem2_info,aendconem1_info,aendcone_info
  real(kind=8)::zendconem2_info,zendconem1_info,zendcone_info
  real(kind=8)::dendconem2_info,dendconem1_info,dendcone_info

  ! Local variables
  integer::ilun
  integer::isfullsky
  integer:: future

  ! Fullsky
  isfullsky=0
  if(is_fullsky) then
     isfullsky=1
  else
     isfullsky=0
  endif
    
  ! Write
  ilun=40*ncpu+myid+10
  open(ilun,file=TRIM(filename),form='formatted')
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("nstride     =",I11)')nstride
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,'("aexp        =",E23.15)')aexp
  write(ilun,'("nlevel      =",I11)')nlevel
  write(ilun,'("nlevel_min  =",I11)')nlevel_min
  write(ilun,'("nlevel_max  =",I11)')nlevel_max
  write(ilun,'("observer_x  =",E23.15)')observer_x
  write(ilun,'("observer_y  =",E23.15)')observer_y
  write(ilun,'("observer_z  =",E23.15)')observer_z
  write(ilun,'("observer_rds=",E23.15)')observer_redshift
  write(ilun,'("cone_id     =",I11)')conegravid
  write(ilun,'("cone_zlim   =",E23.15)')conegravzlim
  write(ilun,'("amax        =",E23.15)')amax
  write(ilun,'("amin        =",E23.15)')amin
  write(ilun,'("zmax        =",E23.15)')zmax
  write(ilun,'("zmin        =",E23.15)')zmin
  write(ilun,'("dmax        =",E23.15)')dmax
  write(ilun,'("dmin        =",E23.15)')dmin
  write(ilun,'("dtol        =",E23.15)')dtol
  write(ilun,'("nglobalfile =",I20)')reduce_infoglobal(1)
  write(ilun,'("nglobalcell =",I20)')reduce_infoglobal(2)
  write(ilun,'("isfullsky   =",I11)')isfullsky
  write(ilun,'("thetay      =",E23.15)')input_thetay
  write(ilun,'("thetaz      =",E23.15)')input_thetaz
  write(ilun,'("theta       =",E23.15)')input_theta
  write(ilun,'("phi         =",E23.15)')input_phi
  write(ilun,'("aendconem2  =",E23.15)')aendconem2_info
  write(ilun,'("aendconem1  =",E23.15)')aendconem1_info
  write(ilun,'("aendcone    =",E23.15)')aendcone_info
  write(ilun,'("zendconem2  =",E23.15)')zendconem2_info
  write(ilun,'("zendconem1  =",E23.15)')zendconem1_info
  write(ilun,'("zendcone    =",E23.15)')zendcone_info
  write(ilun,'("dendconem2  =",E23.15)')dendconem2_info
  write(ilun,'("dendconem1  =",E23.15)')dendconem1_info
  write(ilun,'("dendcone    =",E23.15)')dendcone_info
  write(ilun,'("future      =",I11)')future
  close(ilun)
  
end subroutine write_infoconegrav



!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine writes useful information of simple smaple into a file
!----------------------------------------------------------------------!   
subroutine write_infosamplegrav(ncell_tot,aexp_sample,nlevel,nlevel_min,nlevel_max,      &
                              & xmin,xmax,ymin,ymax,zmin,zmax,reduce_infoglobal,filename)
  use amr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  
  ! Variables from input
  integer::ncell_tot
  real(kind=8)::aexp_sample
  integer::nlevel
  integer::nlevel_min
  integer::nlevel_max
  real(kind=8)::xmin
  real(kind=8)::xmax
  real(kind=8)::ymin
  real(kind=8)::ymax
  real(kind=8)::zmin
  real(kind=8)::zmax
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  character(LEN=200)::filename
    
  ! Local variables
  integer::ilun
  
  ! Write
  ilun=40*ncpu+myid+10
  open(ilun,file=TRIM(filename),form='formatted')
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("nstride     =",I11)')nstride
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,'("aexp        =",E23.15)')aexp_sample
  write(ilun,'("nlevel      =",I11)')nlevel
  write(ilun,'("nlevel_min  =",I11)')nlevel_min
  write(ilun,'("nlevel_max  =",I11)')nlevel_max
  write(ilun,'("xmin        =",E23.15)')xmin
  write(ilun,'("xmax        =",E23.15)')xmax
  write(ilun,'("ymin        =",E23.15)')ymin
  write(ilun,'("ymax        =",E23.15)')ymax
  write(ilun,'("zmin        =",E23.15)')zmin
  write(ilun,'("zmax        =",E23.15)')zmax
  write(ilun,'("nglobalfile =",I20)')reduce_infoglobal(1)
  write(ilun,'("nglobalcell =",I20)')reduce_infoglobal(2)
  close(ilun)
  
end subroutine write_infosamplegrav


!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine calculates the rotation matrix and its transpose. Note
! that the transpose of a rotation matrix is its own inverse matrix, and
! that thetashiftrad and phishiftrad are respectively the rotation angle
! with respect to the y-axis (i.e., rotation in the x-z plane) and the z
! axis (i.e., rotation in the x-y plane). phishiftrad is the same as the
! usual phi in a spherical polar coordinate system, and thetashiftrad is
! the negative of the usual theta in a spheridal polar coordinate system
!----------------------------------------------------------------------!   
subroutine compute_rotation_matrix_grav(thetashiftrad,phishiftrad,rot,rotm1)
  implicit none
  real(kind=8) :: thetashiftrad,phishiftrad
  real(kind=8) :: rot(3,3),rotm1(3,3)
  integer :: i,j
  
  rot(1,1) =  dcos(thetashiftrad)*dcos(phishiftrad)
  rot(1,2) =  dcos(thetashiftrad)*dsin(phishiftrad)
  rot(1,3) = -dsin(thetashiftrad)
  rot(2,1) = -dsin(phishiftrad)
  rot(2,2) =  dcos(phishiftrad)
  rot(2,3) = 0.0D0
  rot(3,1) =  dcos(phishiftrad)*dsin(thetashiftrad)
  rot(3,2) =  dsin(phishiftrad)*dsin(thetashiftrad)
  rot(3,3) =  dcos(thetashiftrad)

  do j=1,3
     do i=1,3
        rotm1(i,j)=rot(j,i)
     enddo
  enddo

end subroutine compute_rotation_matrix_grav



!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine determines the minimum polygon that contains the whole
! lightcone section between z1 and z2. Note that currently the lightcone
! aligns with the x-axis, and this makes things easier: 2 of the 6 faces
! are perpendicular to the x-axis. The face closer to observer has its 4
! vertices coincident with the 4 vertices of the "inner" surface of cone 
! section, and the face farther apart from observer is tangential to the
! "outer" surface of the cone section.
!----------------------------------------------------------------------!   
subroutine compute_minimum_polygon_grav(x1,x2,thetayrad,thetazrad,sl)
  implicit none
  real(kind=8) :: x1,x2,thetayrad,thetazrad,sl(3,8)
  real(kind=8) :: r(3),axis(3)

  sl(1,1:4)=x1/dsqrt(1.0D0+dtan(thetayrad)**2+dtan(thetazrad)**2)       ! x coordinate of 4 inner vertices
  sl(2,1  )=-sl(1,1)*dtan(thetayrad)                                    ! y coordinate of inner vertex 1
  sl(3,1  )=-sl(1,1)*dtan(thetazrad)                                    ! z coordinate of inner vertex 1
  sl(2,2  )= sl(1,1)*dtan(thetayrad)                                    ! y coordinate of inner vertex 2
  sl(3,2  )=-sl(1,1)*dtan(thetazrad)                                    ! z coordinate of inner vertex 2
  sl(2,3  )=-sl(1,1)*dtan(thetayrad)                                    ! y coordinate of inner vertex 3
  sl(3,3  )= sl(1,1)*dtan(thetazrad)                                    ! z coordinate of inner vertex 3
  sl(2,4  )= sl(1,1)*dtan(thetayrad)                                    ! y coordinate of inner vertex 4
  sl(3,4  )= sl(1,1)*dtan(thetazrad)                                    ! z coordinate of inner vertex 4
  !
  sl(1,5:8)= x2                                                         ! x coordinate of 4 outer vertices
  sl(2,5  )=-x2*dtan(thetayrad)                                         ! y coordinate of outer vertex 1
  sl(3,5  )=-x2*dtan(thetazrad)                                         ! z coordinate of outer vertex 1 
  sl(2,6  )= x2*dtan(thetayrad)                                         ! y coordinate of outer vertex 2
  sl(3,6  )=-x2*dtan(thetazrad)                                         ! z coordinate of outer vertex 2
  sl(2,7  )=-x2*dtan(thetayrad)                                         ! y coordinate of outer vertex 3
  sl(3,7  )= x2*dtan(thetazrad)                                         ! z coordinate of outer vertex 3
  sl(2,8  )= x2*dtan(thetayrad)                                         ! y coordinate of outer vertex 4
  sl(3,8  )= x2*dtan(thetazrad)                                         ! z coordinate of outer vertex 4

end subroutine compute_minimum_polygon_grav


!----------------------------------------------------------------------!   
!----------------------------------------------------------------------!   
! This subroutine determines how many replicas are needed to enclose the
! lightcone section betwene z1 and z2, and finds their indices. It first
! finds a (minimal) polygon that encloses the cone section, and computes
! the coordinates of the 8 vertices of this polygon. By find the minimum
! and maximum of these coordinates, it gets how many replicas are needed
! to enclose the cone section. The replica that contains the observer is
! given index (0,0,0), and other replicas are given indices accordingly,
! e.g., (-1,-1,2) indicates the replica that is first to the left in the 
! x-axis, first to the left in the y-axis and second to the right in the 
! z-axis, with respect to the observer's replica. 
!----------------------------------------------------------------------!   
subroutine compute_replica_grav(thetayrad,thetazrad,dist1,dist2,observer,Lbox,rot, &
     &                          nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp)
  !
  implicit none
  !
  real(kind=8)::thetayrad,thetazrad,observer(3),Lbox,rot(3,3),dist1,dist2
  integer::nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
  integer::myint
  real(kind=8)::sl(3,8),slfr(3)
  real(kind=8)::xplmin,xplmax,yplmin,yplmax,zplmin,zplmax
  integer::i,j

  ! Compute the coordinates of the vertices of the polygon
  call compute_minimum_polygon_grav(dist1,dist2,thetayrad,thetazrad,sl)

  do j=1,8
     do i=1,3
        slfr(i)=sl(1,j)*rot(1,i) &
             & +sl(2,j)*rot(2,i) &
             & +sl(3,j)*rot(3,i)                                        ! rotate coordinates of the polygon
     enddo
     if(j.eq.1) then
        xplmin=slfr(1)
        yplmin=slfr(2)
        zplmin=slfr(3)
        xplmax=xplmin
        yplmax=yplmin
        zplmax=zplmin
     else
        xplmin=dmin1(xplmin,slfr(1))                                    ! smallest x coordinate of all vertices
        xplmax=dmax1(xplmax,slfr(1))                                    ! largest  x coordinate of all vertices
        yplmin=dmin1(yplmin,slfr(2))                                    ! smallest y coordinate of all vertices
        yplmax=dmax1(yplmax,slfr(2))                                    ! largest  y coordinate of all vertices
        zplmin=dmin1(zplmin,slfr(3))                                    ! smallest z coordinate of all vertices
        zplmax=dmax1(zplmax,slfr(3))                                    ! largest  z coordinate of all vertices
     endif
  enddo
  !
  nrepxm=myint((xplmin+observer(1))/Lbox)                               ! smallest x index of all replicas
  nrepxp=myint((xplmax+observer(1))/Lbox)                               ! largest  x index of all replicas
  nrepym=myint((yplmin+observer(2))/Lbox)                               ! smallest y index of all replicas
  nrepyp=myint((yplmax+observer(2))/Lbox)                               ! largest  y index of all replicas
  nrepzm=myint((zplmin+observer(3))/Lbox)                               ! smallest z index of all replicas
  nrepzp=myint((zplmax+observer(3))/Lbox)                               ! largest  z index of all replicas  

end subroutine compute_replica_grav
