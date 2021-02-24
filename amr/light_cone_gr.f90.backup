!-----------------------------------------------------------!
! Gravity Lightcone 
! Original implementation by Rasera and Reverdy
! Based on RAMSES particle LC implementations
!-----------------------------------------------------------!
!
!
!
!
!------------------ MODIF V. REVERDY 2011 ------------------!
!===========================================================================
! FULL SKY GRAVITY CONE
!===========================================================================
subroutine output_conegrav(is_fullsky, filedir, filename, cone_id, observer_x, observer_y, observer_z, observer_redshift, cone_zmax, cone_mode, input_levelmin, input_levelmax, input_thetay, input_thetaz, input_theta, input_phi,future)
  
  !=========================================================!
  
  ! Initialization
  use amr_commons
  use poisson_commons
  use gr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none

  !=========================================================!
  
  ! Misc
  integer,parameter::tag=2800
  integer::dummy_io
  integer::info
  integer::ilun
  logical::is_opened

  ! Input
  logical::is_fullsky         ! option to specify is the cone is fullsky or not
  character(LEN=*)::filedir   ! name of the output dir
  character(LEN=*)::filename  ! name of the output file
  integer::cone_id            ! cone identifier
  real(dp)::observer_x        ! x coordinate of the observer (usually 0.5)
  real(dp)::observer_y        ! y coordinate of the observer (usually 0.5)
  real(dp)::observer_z        ! z coordinate of the observer (usually 0.5)
  real(dp)::observer_redshift ! redshift of the observer
  real(dp)::cone_zmax         ! zmax of the cone
  integer::cone_mode          ! mode = 0 -> only leaves, 1 -> only coarse, 2 -> all grid hierarchy
  integer::input_levelmin     ! levelmax in case of all grid hierarchy
  integer::input_levelmax     ! levelmax in case of all grid hierarchy
  real(dp)::input_thetay      ! angle if non-fullksy
  real(dp)::input_thetaz      ! angle if non-fullksy
  real(dp)::input_theta       ! angle if non-fullksy
  real(dp)::input_phi         ! angle if non-fullksy
  
  ! Level
  integer::cone_levelmin
  integer::cone_levelmax
  integer::cone_nlevel
  
  ! Variables
  real(kind=8)::z2                ! old redshift
  real(kind=8)::z1                ! recent redshift
  real(kind=8)::dist2             ! distance max
  real(kind=8)::dist1             ! distance min
  real(kind=8)::om0in             ! cosmological param
  real(kind=8)::omLin             ! cosmological param
  real(kind=8)::hubin             ! cosmological param
  real(kind=8)::Lbox              ! cosmological param
  real(kind=8)::observer(3)       ! position of the observer
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
  real(kind=8),dimension(1:twotondim,1:3)::xc     ! centers of cells
  
  ! Variables : common for arrays
  !integer::nvar = 3 ! fx, fy, fz
  integer,parameter::nvar=5 ! fx, fy, fz, phi, rho
                             !  1       3       1       3        3          9 
  !integer,parameter::nvar=20 ! Phi, grad(Phi), Psi, grad(Psi), beta^i, grad(beta^i) ! CBH_LC
  integer,parameter::nbool=1 ! son
  integer,parameter::nvectorgrid = floor(nvector/(2.**ndim))  ! nvector/8
  integer,dimension(1:nvectorgrid)::ind_grid  ! indices of grids
  integer,dimension(1:nvector)::ind_cell      ! indices of cells
  
  ! Variables : input arrays
  real(kind=8),dimension(1:nvector,1:ndim)::xc_in
  real(kind=8),dimension(1:nvector,1:nvar)::var_in
  integer,dimension(1:nvector,1:nbool)::bool_in

  ! Variables : selection arrays
  real(kind=8),dimension(1:27*nvector,1:ndim)::xc_out
  real(kind=8),dimension(1:27*nvector,1:nvar)::var_out
  integer,dimension(1:27*nvector,1:nbool)::bool_out
  
  ! Variables : write arrays
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim)::xc_wr      ! x, y, z centers of cells for output
  real(kind=4),dimension(1:nstride+27*nvector,1:nvar)::var_wr     ! gravity vars
  integer,dimension(1:nstride+27*nvector,1:nbool)::bool_wr        ! son
  
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
  
  ! aexp for info file
  real(kind=8)::aendconem2_info,aendconem1_info,aendcone_info !expansion factor for end of the shell
  !the end of previous shell and the end of previous previous shell (which is the begining of the current shell)
  !when activating overlapping buffer option.
  real(kind=8)::zendconem2_info,zendconem1_info,zendcone_info
  real(kind=8)::dendconem2_info,dendconem1_info,dendcone_info

  real(kind=8)::Omega0,OmegaL,OmegaR,coverH0
  real(kind=8)::coord_distance
  integer:: future
  !=========================================================!
  
  ! Correct quantities
  cone_levelmin = input_levelmin
  cone_levelmax = input_levelmax
  if ((cone_levelmin .EQ. 0).OR.(cone_levelmin.GT.nlevelmax)) cone_levelmin = levelmin
  if ((cone_levelmax .EQ. 0).OR.(cone_levelmax.GT.nlevelmax)) cone_levelmax = nlevelmax
  cone_nlevel = (cone_levelmax-cone_levelmin+1)

  ! Filename
  ilun=30*ncpu+myid+10
  call title(myid, nchar_cpu)
  call title(nstep_coarse, nchar_coarse)
  write(nchar_cone,'(I5.5)') cone_id
  if(future==0) then
     if (is_fullsky) then
        filename_info = TRIM(filedir)//'info_cone_grav_fullsky_pastnfut_'//nchar_cone//'_ncoarse_'//TRIM(nchar_coarse)//'.txt'
     else
        filename_info = TRIM(filedir)//'info_cone_grav_narrow_pastnfut_'//nchar_cone//'_ncoarse_'//TRIM(nchar_coarse)//'.txt'
     endif
  else
     print*,'cone grav not implemented future=',future
     stop
  endif
  
  filename_dat = TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.dat'
  filename_hdr = TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.hdr'

  ! Special cases
  if(nstep_coarse .LT. 4) return
  
  ! Compute quantities
  z2=1./aexp_old-1.
  z1=1./aexp-1.
  observer(1)=observer_x
  observer(2)=observer_y
  observer(3)=observer_z
  om0in=omega_m
  omLin=omega_l
  hubin=h0/100.D0
  Lbox=boxlen_ini/hubin
  if((use_aexp_restart).AND.(nstep_coarse_after_restart==2)) z2=1./aexp_restart_light_cone-1.
  if(conegrav_overlap) then
     if((aendconem1.lt.aendconem2).or.(aendcone.lt.aendconem1))print*,'WARNING CONEGRAV SHELL RANGE NOT WELL ORDERED, AEXP ',aendconem2,aendconem1,aendcone 
     z2=1./aendconem2-1.
     z1=1./aendcone-1.
     if (myid==1) then
        print*,''
        print*,'TEST CONE LIMIT'
        print*,'nstepcoarse,z2,z1',nstep_coarse,z2,z1
        print*,'aendconem2,aendconem1,aendcone',aendconem2,aendconem1,aendcone
        print*,''
     endif
  endif
  


  !zmax_cone_full=0.25
  !if(myid.EQ.1) write(*,*)'Z2 = ',z2, 'zmax_cone_full = ',zmax_cone_full,'z2-z1 = ',(z2-z1)
  !if(z2.gt.zmax_cone_full)return
  if(z2.gt.cone_zmax)return
  if(abs(z2-z1)<1d-6)return
  if((z1.LT.observer_redshift).AND.(z2.LT.observer_redshift))return

  ! Hdr arrays
  hdr_arraysize = cone_nlevel
  allocate(hdr_ilevel(1:cone_nlevel))
  allocate(hdr_dx(1:cone_nlevel))
  allocate(hdr_ncellperlevel(1:cone_nlevel))
  hdr_ilevel = 0
  hdr_dx = 0
  hdr_ncellperlevel = 0

  ! Compute iogroupsize
  call perform_my_selection_simple_conegrav_ramsesunits(.true., &
       &                                .false.,disttol,z1,z2,dist1,dist2, &
       &                                om0in,omLin,hubin,Lbox, &
       &                                observer,observer_redshift,dxcell,nvar,nbool, &
       &                                xc_in,var_in,bool_in,ncell_in, &
       &                                xc_out,var_out,bool_out,ncell_out,.false.)
  iovolume = ((2.0D0*dist2)**3.0D0)-(((2.D0*dist1)/SQRT(real(ndim,kind=8)))**3.0D0)
  if (.NOT.is_fullsky) iovolume = iovolume*(input_thetay*input_thetaz)/(41263.D0)
  if (adaptive_iogroupsize) IOGROUPSIZECONE = nint(iovolume*real(IOGROUPSIZE, kind=8))
  if (adaptive_iogroupsize.AND.(.NOT.is_fullsky)) IOGROUPSIZECONE = IOGROUPSIZECONE/8
  if (IOGROUPSIZECONE .LE. 1) IOGROUPSIZECONE = 1
  if (IOGROUPSIZECONE .GE. IOGROUPSIZE) IOGROUPSIZECONE = IOGROUPSIZE
  if (verbose) write(*,*) cone_id, 'IOGROUPSIZECONE = ',IOGROUPSIZECONE, 'iovolume = ', iovolume

  ! Wait for the token
  if(conetoken) then 
    if (mod(myid-1,IOGROUPSIZECONE)/=0) then
      call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                  & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
    end if
  endif

  !=========================================================!
  
  ! Loop over levels
  reduce_infolocal = 0
  reduce_infoglobal = 0
  ncell_total = 0
  ihdr = 0
  is_opened = .FALSE.

  do ilevel=cone_levelmin,cone_levelmax

    !*******************************************************!
    ! Initialization of counters
    ncell_level = 0
    ivector = 0
    istride = 0
    icell_in = 0 
    icell_out = 0 
    icell_wr = 0
    ncell_in = 0
    ncell_out = 0 
    ncell_wr = 0
    
    ! Initialization of arrays
    xc_in = 0.D0
    var_in = 0.D0
    bool_in = 0
    xc_out = 0.D0
    var_out = 0.D0
    bool_out = 0
    xc_wr = 0.D0
    var_wr = 0.D0
    bool_wr = 0
    
    ! Index
    ihdr = ihdr+1

    
    !*******************************************************!
    ! Level constants
    dx=0.5D0**ilevel ! 
    dxcell=dx ! need to check (maybe dxcell = (1/2.)*dx

    ! Set position of cell centers relative to grid center
    do ind=1,twotondim
      iz=(ind-1)/4
      iy=(ind-1-4*iz)/2
      ix=(ind-1-2*iy-4*iz)
      if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
      if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
      if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
    end do
    
    !*******************************************************!
    
    ! Loop over myid grids by vectorgrid sweeps
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvectorgrid
      ngrid=MIN(nvectorgrid,ncache-igrid+1)
      do i=1,ngrid
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         !ind_grid(i)=active(ilevel)%pcomm%igrid(igrid+i-1) RAMSES-LC line CBH_LC 16-02-2021
      end do
      
      ! Loop over cells
      ivector = 0 ! ivector = icell_in
      do i=1,ngrid
        do ind=1,twotondim
        
          ! Increment ivector
          ivector=ivector+1

          ! Gather cell indices
          iskip=ncoarse+(ind-1)*ngridmax
          ind_cell(ivector)=iskip+ind_grid(i)
          
          ! Compute cell center
          do idim=1,ndim
            xc_in(ivector,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
          end do
          
          ! Compute var :
          var_in(ivector,1)=f(ind_cell(ivector),1)
          var_in(ivector,2)=f(ind_cell(ivector),2)
          var_in(ivector,3)=f(ind_cell(ivector),3)
          var_in(ivector,4)=phi(ind_cell(ivector))
          var_in(ivector,5)=rho(ind_cell(ivector))

          ! CBH_LC
          ! Notice that gradients need to be copied from f() at the correct igrp sequence in move or synchro
          ! At igrp = 10, f() = - grad_i(b)
          !var_in(ivector,1 )=gr_pot(ind_cell(ivector),5) ! Phi
          !var_in(ivector,2 )=gr_pot(ind_cell(ivector),6) ! Xi 
          !var_in(ivector,3 )=gr_pot(ind_cell(ivector),7) - f(ind_cell(ivector),1)  ! \beta^i = B^i + grad_i(b)
          !var_in(ivector,4 )=gr_pot(ind_cell(ivector),8) - f(ind_cell(ivector),2)
          !var_in(ivector,5 )=gr_pot(ind_cell(ivector),9) - f(ind_cell(ivector),3)
          !var_in(ivector,6 )=gr_pot_grad(ind_cell(ivector),1 )    ! grad_i Phi
          !var_in(ivector,7 )=gr_pot_grad(ind_cell(ivector),2 )
          !var_in(ivector,8 )=gr_pot_grad(ind_cell(ivector),3 )
          !var_in(ivector,9 )=gr_pot_grad(ind_cell(ivector),4 )    ! grad_i Xi
          !var_in(ivector,10)=gr_pot_grad(ind_cell(ivector),5 )
          !var_in(ivector,11)=gr_pot_grad(ind_cell(ivector),6 )
          !var_in(ivector,12)=gr_pot_grad(ind_cell(ivector),7 )    ! grad_i \beta^j
          !var_in(ivector,13)=gr_pot_grad(ind_cell(ivector),8 )
          !var_in(ivector,14)=gr_pot_grad(ind_cell(ivector),9 )
          !var_in(ivector,15)=gr_pot_grad(ind_cell(ivector),10)
          !var_in(ivector,16)=gr_pot_grad(ind_cell(ivector),11)
          !var_in(ivector,17)=gr_pot_grad(ind_cell(ivector),12)
          !var_in(ivector,18)=gr_pot_grad(ind_cell(ivector),13)
          !var_in(ivector,19)=gr_pot_grad(ind_cell(ivector),14)
          !var_in(ivector,20)=gr_pot_grad(ind_cell(ivector),15)
          
          ! END CBH_LC
          
          ! Compute bool
          bool_in(ivector,1) = son(ind_cell(ivector))
          
        end do 
      end do ! end of loop over cells : arrays of size ivector<=nvector of data have been generated
      
      !*******************************************************!
      
      ! Perfom selection on cells
      ncell_in = ivector
      if (is_fullsky) then
          call perform_my_selection_simple_conegrav_ramsesunits(.false., &
          &                                .false.,disttol,z1,z2,dist1,dist2, &
          &                                om0in,omLin,hubin,Lbox, &
          &                                observer,observer_redshift,dxcell,nvar,nbool, &
          &                                xc_in,var_in,bool_in,ncell_in, &
          &                                xc_out,var_out,bool_out,ncell_out,.false.)
      else
          call perform_my_selection_conegrav_ramsesunits(.false., &
          &                                .false.,disttol,z1,z2,dist1,dist2, &
          &                                om0in,omLin,hubin,Lbox, &
          &                                observer,observer_redshift,dxcell,nvar,nbool, &
          &                                xc_in,var_in,bool_in,ncell_in, &
          &                                xc_out,var_out,bool_out,ncell_out,.false., &
          &                                input_thetay, input_thetaz, input_theta, input_phi)      
      endif
      
      !*******************************************************!
      
      ! Transfer out arrays of size ncell_out in write arrays of size icell_wr
      if(ncell_out.GT.0) then
        call transfer_grav(nvar,nbool,ncell_out,icell_wr,ncell_level,ncell_total, &
                                & xc_out,var_out,bool_out, &
                                & xc_wr,var_wr,bool_wr)
      end if
      
      !*******************************************************!
      
      ! If nstride has been reached, write arrays of size icell_wr and update size
      ! icell_wr = istride
      if (icell_wr.GE.nstride) then
        call write_grav(.false., is_opened, ilun, filename_dat, &
                            & xc_wr, nvar, var_wr, nbool, bool_wr, &
                            & icell_wr, cone_nlevel, reduce_infolocal)
      end if
      
      !*******************************************************!
      
    end do ! loop over grids
    !*******************************************************!
    
    ! If some cells lasts in write vector, write them
    
    if (icell_wr.GT.0) then
      call write_grav(.true., is_opened, ilun, filename_dat, &
                         & xc_wr, nvar, var_wr, nbool, bool_wr, &
                         & icell_wr, cone_nlevel, reduce_infolocal)
    end if

  hdr_ilevel(ihdr) = ilevel
  hdr_dx(ihdr) = dx
  hdr_ncellperlevel(ihdr) = ncell_level

  end do ! loop over levels
  if (is_opened) close(ilun)
  
  !=========================================================!

  ! Write hdr files
  call write_hdr(ncell_total, filename_hdr, hdr_arraysize, hdr_ilevel, hdr_dx, hdr_ncellperlevel)
  
  !=========================================================!

  ! Send the token
  if(conetoken) then 
    if(mod(myid,IOGROUPSIZECONE)/=0 .and.(myid.lt.ncpu))then
      dummy_io=1
      call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
           & MPI_COMM_WORLD,info)
    end if
  endif

  !=========================================================!
  
  ! Reduce
  reduce_infolocal(2) = ncell_total
  call MPI_REDUCE(reduce_infolocal, reduce_infoglobal, 2, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  ! Write info file
  if (myid.EQ.1) then
     if(conegrav_overlap)then
        if(future==0) then
           aendconem2_info=aendconem2
           aendconem1_info=aendconem1
           aendcone_info=aendcone
           
           zendconem2_info=1./aendconem2_info-1.
           zendconem1_info=1./aendconem1_info-1.
           zendcone_info  =1./aendcone_info-1.
        else
           print*,'cone grav not implemented future=',future
           stop
        endif

        call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
        if (observer_redshift.GT.0) then
           dendconem2_info=((coord_distance(zendconem2_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
           dendconem1_info=((coord_distance(zendconem1_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
           dendcone_info  =((coord_distance(zendcone_info  ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
        else 
           dendconem2_info=(coord_distance(zendconem2_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
           dendconem1_info=(coord_distance(zendconem1_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
           dendcone_info  =(coord_distance(zendcone_info   ,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
        end if
     else
        aendconem2_info=0.
        aendconem1_info=0.
        aendcone_info=0.
        zendconem2_info=0.
        zendconem1_info=0.
        zendcone_info=0.
        dendconem2_info=0.
        dendconem1_info=0.
        dendcone_info=0.
     endif
  call write_infoconegrav(ncell_total, cone_id, cone_zmax, observer_x, observer_y, observer_z, observer_redshift, &
               & cone_nlevel, cone_levelmin, cone_levelmax, &
               & 1./(z2+1.),1./(z1+1.), z2, z1, dist2, dist1, disttol, &
               & reduce_infoglobal, filename_info, &
               & is_fullsky, input_thetay, input_thetaz, input_theta, input_phi,&
               aendconem2_info,aendconem1_info,aendcone_info,&
               zendconem2_info,zendconem1_info,zendcone_info,&
               dendconem2_info,dendconem1_info,dendcone_info,future)
  end if
  
  !=========================================================!
  
  ! Finalize
  deallocate(hdr_ilevel)
  deallocate(hdr_dx)
  deallocate(hdr_ncellperlevel)
  
  !=========================================================!
  
end subroutine output_conegrav

!===========================================================================
! PERFORM MY SELECTION : FULLSKY VERSION (NO ROTATION)
!===========================================================================
subroutine perform_my_selection_simple_conegrav_ramsesunits(onlydist, &
     &                                justcount,disttol,z1,z2,dist1,dist2, &
     &                                om0in,omLin,hubin,Lbox, &
     &                                observer,observer_redshift,dxcell,nvar,nbool, &
     &                                xc_in,var_in,bool_in,ncell_in, &
     &                                xc_out,var_out,bool_out,ncell_out,verbose)
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
  real(kind=8)::demidiag ! demi great diagonal of a cube
  real(kind=8)::demidx   ! dx/2
  
  ! Beginning of function
  if (verbose) write(*,*) 'Enter perform_my_selection_simple'
  
  ! Initialize cosmological parameters
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)

  if (verbose) write(*,*) 'After init_cosmo',Omega0,OmegaL,OmegaR,coverH0

  
  ! Compute comoving distance of the photon planes from the observer
  ! dist1,dist2=integral of c.dt/a between zero and z1,z2
  if (observer_redshift.GT.0) then
    dist1=((coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
    dist2=((coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)+disttol ! in [0, 1]
  else 
    dist1=(coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
    dist2=(coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)/Lbox)+disttol ! in [0, 1]
  end if
  
  if (verbose) write(*,*) 'After coord',dist1,dist2
  
  if (.NOT.onlydist) then
     ! Compute the demi diagonal of the cube
     demidx=dxcell/2.D0
     demidiag=sqrt(3.D0)*demidx

     ! Compute how many replica are needed
     nrepxm=myint(observer(1)-dist2) ! in [0, 1]
     nrepxp=myint(observer(1)+dist2) ! in [0, 1]
     nrepym=myint(observer(2)-dist2) ! in [0, 1]
     nrepyp=myint(observer(2)+dist2) ! in [0, 1]
     nrepzm=myint(observer(3)-dist2) ! in [0, 1]
     nrepzp=myint(observer(3)+dist2) ! in [0, 1]

     facnorm=1.0d0/(dist2-dist1)
  
     ncell_out=0   
     icell_out=0
     if (verbose) write(*,*) 'before loop',facnorm
     ! loop on all the replica of potential interest
     do k=nrepzm,nrepzp,1
        do j=nrepym,nrepyp,1
           do i=nrepxm,nrepxp,1
              do icell_in=1,ncell_in
           
                 ! Compute center (replica)
                 xcoord=xc_in(icell_in,1)+dble(i)-observer(1) ! in [0, n*replica]
                 ycoord=xc_in(icell_in,2)+dble(j)-observer(2) ! in [0, n*replica]
                 zcoord=xc_in(icell_in,3)+dble(k)-observer(3) ! in [0, n*replica]
                 
                 ! Compute distance to the center
                 dist=sqrt(xcoord**2+ycoord**2+zcoord**2)
                 distmin=dist-demidiag
                 distmax=dist+demidiag
                 
!RY modif
!!!!                 if ((distmin > dist1 .AND. distmax <= dist2) .OR. (dist > dist1 .AND. dist <= dist2)) then
                 if ((distmax > dist1 .AND. distmin <= dist2) .OR. (dist > dist1 .AND. dist <= dist2)) then
!end RY modif

!!!!!!WARNING
                 !if (dist > dist1 .AND. dist <= dist2) then
!!!!!!WARNING 
                    ! This particle is good, we can add it to the list
            
                    icell_out=icell_out+1
                    if (icell_out.GE.27*nvector) then 
                       write(*,*) 'PROBLEM : icell_out = ',icell_out, '27*nvector = ', 27*nvector, ' ncell_in = ',ncell_in
                       write(*,*) nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
                    endif
                    xc_out(icell_out,1)=xcoord
                    xc_out(icell_out,2)=ycoord
                    xc_out(icell_out,3)=zcoord
                 
                    ! Compute the redshift of the particle using linear
                    dxtest1=dist-dist1
                    dxtest2=dist2-dist

                    ! Copy other variables
                    do ivar=1,nvar
                       var_out(icell_out,ivar)=var_in(icell_in,ivar)
                    end do
                    do ibool=1,nbool
                       bool_out(icell_out,ibool)=bool_in(icell_in,ibool)
                    end do

                 endif
              enddo
           enddo
        enddo
     enddo
     ncell_out = icell_out
  endif
  
  if (verbose) write(*,*) 'End of perform_my_selection_simple',ncell_out
end subroutine perform_my_selection_simple_conegrav_ramsesunits
!===========================================================================



!===========================================================================
! PERFORM MY SELECTION : NARROW VERSION (ROTATION)
!===========================================================================
subroutine perform_my_selection_conegrav_ramsesunits(onlydist, &
     &                                justcount,disttol,z1,z2,dist1,dist2, &
     &                                om0in,omLin,hubin,Lbox, &
     &                                observer,observer_redshift,dxcell,nvar,nbool, &
     &                                xc_in,var_in,bool_in,ncell_in, &
     &                                xc_out,var_out,bool_out,ncell_out,verbose,&
     &                                input_thetay, input_thetaz, input_theta, input_phi)
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
  real(kind=8)::pi, thetarad, phirad, thetayrad, thetazrad, tanybound, tanzbound
  real(kind=8)::dist1Lbox, dist2Lbox, tany, tanz
  real(kind=8),dimension(1:9)::xcoordfr, ycoordfr, zcoordfr, xcoord, ycoord, zcoord, dxcoord, dycoord, dzcoord
  real(kind=8),dimension(1:3, 1:3)::rot, rotm1
  integer::ivertex
  logical::ok
  
  ! Variables for conegrav
  real(kind=8)::demidiag ! demi great diagonal of a cube
  real(kind=8)::demidx   ! dx/2
  
  ! Beginning of function
  small = 1.D-5
  if (verbose) write(*,*) 'Enter perform_my_selection_simple'
  
  ! Initialize cosmological parameters
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  
  ! Deal with angles
  pi = acos(-1.0d0)
  thetarad = input_theta*pi/180.0d0
  phirad = input_phi*pi/180.0d0
  thetayrad = input_thetay*pi/180.0d0
  thetazrad = input_thetaz*pi/180.0d0
  call compute_rotation_matrix_grav(thetarad,phirad,rot,rotm1)
  
  ! Verbose
  if (verbose) write(*,*) 'After init_cosmo',Omega0,OmegaL,OmegaR,coverH0
  
  ! Compute comoving distance of the photon planes from the observer : dist1,dist2=integral of c.dt/a between zero and z1,z2
  if (observer_redshift.GT.0) then
    dist1=((coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)-disttol ! in [0, 1]
    dist2=((coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox)+disttol ! in [0, 1]
  else 
    dist1=(coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)/Lbox)-disttol ! in [0, 1]
    dist2=(coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)/Lbox)+disttol ! in [0, 1]
  end if
  dist1Lbox = dist1*Lbox
  dist2Lbox = dist2*Lbox
  observerLbox=observer*Lbox

  ! Verbose
  if (verbose) write(*,*) 'After coord',dist1,dist2
  
  ! Actual computation
  if (.NOT.onlydist) then

    ! Compute the set of replica to be considered
    call compute_replica_grav(thetayrad,thetazrad,dist1Lbox,dist2Lbox,observerLbox,Lbox,rot, &
    &                       nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp)    
      
    tanybound=tan(thetayrad)
    tanzbound=tan(thetazrad)
  
    ! Compute specific points of the cube
    demidx=dxcell/2.D0
    demidiag=sqrt(3.D0)*demidx
    dxcoord = (/ 0.0D0, -demidx, demidx, -demidx, demidx, -demidx, demidx, -demidx, demidx /)
    dycoord = (/ 0.0D0, -demidx, -demidx, demidx, demidx, -demidx, -demidx, demidx, demidx /)
    dzcoord = (/ 0.0D0, -demidx, -demidx, -demidx, -demidx, demidx, demidx, demidx, demidx /)
     
    ! Loop over cells
    ncell_out=0   
    icell_out=0
    if (verbose) write(*,*) 'before loop'   
    do k=nrepzm,nrepzp,1
      do j=nrepym,nrepyp,1
        do i=nrepxm,nrepxp,1
          do icell_in=1,ncell_in
           
            ! Compute center (replica)
            xcoordfr=xc_in(icell_in,1)+dble(i)-observer(1)+dxcoord ! in [0, n*replica]
            ycoordfr=xc_in(icell_in,2)+dble(j)-observer(2)+dycoord ! in [0, n*replica]
            zcoordfr=xc_in(icell_in,3)+dble(k)-observer(3)+dzcoord ! in [0, n*replica]

            ! Rotation to get in the framework of the photon plane
            xcoord=xcoordfr*rotm1(1,1)+ &
                  & ycoordfr*rotm1(2,1)+ &
                  & zcoordfr*rotm1(3,1)
            ycoord=xcoordfr*rotm1(1,2)+ &
                  & ycoordfr*rotm1(2,2)+ &
                  & zcoordfr*rotm1(3,2)
            zcoord=xcoordfr*rotm1(1,3)+ &
                  & ycoordfr*rotm1(2,3)+ &
                  & zcoordfr*rotm1(3,3)
                     
            ! Loop over vertexes
            ok = .false.
            ivertex = 1
            do while ((.not.ok) .and. (ivertex.le.9))
              if (xcoord(ivertex) .gt. small) then
                tany=abs(ycoord(ivertex)/xcoord(ivertex))
                tanz=abs(zcoord(ivertex)/xcoord(ivertex))
                !Begin RY make it compatible with full sky
                !!!dist=sqrt(xcoord(ivertex)**2+ycoord(ivertex)**2+zcoord(ivertex)**2) 
                !!!ok = ((tany .le. tanybound) .and. (tanz .le. tanzbound) .and. (dist .ge. dist1) .and. (dist .le. dist2))
                dist=sqrt(xcoord(1)**2+ycoord(1)**2+zcoord(1)**2) !RY make it compatible with full cone
                distmin=dist-demidiag
                distmax=dist+demidiag
                ok = ((tany .le. tanybound) .and. (tanz .le. tanzbound) .and. ((distmax > dist1 .AND. distmin <= dist2) .OR. (dist > dist1 .AND. dist <= dist2)))
                !End RY make it compatible with full sky
                

                
              end if
              ivertex = ivertex+1
            end do

            ! If one vertex is ok, add it to the list
            if (ok) then
              icell_out=icell_out+1
              if (icell_out.GE.27*nvector) then 
                write(*,*) 'PROBLEM : icell_out = ',icell_out, '27*nvector = ', 27*nvector, ' ncell_in = ',ncell_in
                write(*,*) nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
              endif
              xc_out(icell_out,1)=xcoord(1)
              xc_out(icell_out,2)=ycoord(1)
              xc_out(icell_out,3)=zcoord(1)

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
              do ibool=1,nbool
                bool_out(icell_out,ibool)=bool_in(icell_in,ibool)
              end do
            endif
          
          enddo
        enddo
      enddo
    enddo
    ncell_out = icell_out
  endif
  
  if (verbose) write(*,*) 'End of perform_my_selection',ncell_out
end subroutine perform_my_selection_conegrav_ramsesunits
!===========================================================================



!===========================================================================
! GRAVITY SAMPLE (same as conegrav)
!===========================================================================
subroutine extract_samplegrav(filedir, filename, xmin, xmax, ymin, ymax, zmin, zmax, nsample, input_levelmin, input_levelmax)
  
  !=========================================================!
  
  ! Initialization
  use amr_commons
  use poisson_commons
  use gr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none

  !=========================================================!
  
  ! Misc
  integer,parameter::tag=2900
  integer::dummy_io
  integer::info
  integer::ilun
  logical::is_opened

  ! Input
  character(LEN=*)::filedir   ! name of the output dir
  character(LEN=*)::filename  ! name of the output file
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
  !integer::nvar = 3 ! fx, fy, fz
  integer,parameter::nvar=5 ! fx, fy, fz, phi, rho
                             !  1       3       1       3        3          9 
  !integer,parameter::nvar=20 ! Phi, grad(Phi), Psi, grad(Psi), beta^i, grad(beta^i) ! CBH_LC
  integer,parameter::nbool=1 ! son
  integer,parameter::nvectorgrid = floor(nvector/(2.**ndim))  ! nvector/8
  integer,dimension(1:nvectorgrid)::ind_grid  ! indices of grids
  integer,dimension(1:nvector)::ind_cell      ! indices of cells
  integer::ivar
  integer::ibool
  
  ! Variables : input arrays
  real(kind=8),dimension(1:nvector,1:ndim)::xc_in
  real(kind=8),dimension(1:nvector,1:nvar)::var_in
  integer,dimension(1:nvector,1:nbool)::bool_in

  ! Variables : selection arrays
  real(kind=8),dimension(1:27*nvector,1:ndim)::xc_out
  real(kind=8),dimension(1:27*nvector,1:nvar)::var_out
  integer,dimension(1:27*nvector,1:nbool)::bool_out
  
  ! Variables : write arrays
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim)::xc_wr   ! x, y, z centers of cells for output
  real(kind=4),dimension(1:nstride+27*nvector,1:nvar)::var_wr  ! gravity vars
  integer,dimension(1:nstride+27*nvector,1:nbool)::bool_wr     ! son
  
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

  ierr = 0

  !=========================================================!
  
  ! Correct quantities
  sample_levelmin = input_levelmin
  sample_levelmax = input_levelmax
  if ((sample_levelmin .EQ. 0).OR.(sample_levelmin.GT.nlevelmax)) sample_levelmin = levelmin
  if ((sample_levelmax .EQ. 0).OR.(sample_levelmax.GT.nlevelmax)) sample_levelmax = nlevelmax
  sample_nlevel = (sample_levelmax-sample_levelmin+1)

  ! Filename
  ilun=30*ncpu+myid+10
  call title(myid, nchar_cpu)
  call title(nstep_coarse, nchar_coarse)
  filename_info = TRIM(filedir)//'info_sample_grav_ncoarse_'//TRIM(nchar_coarse)//'.txt'
  filename_dat = TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.dat'
  filename_hdr = TRIM(filedir)//TRIM(filename)//TRIM(nchar_cpu)//'.hdr'

  ! Hdr arrays
  hdr_arraysize = sample_nlevel
  allocate(hdr_ilevel(1:hdr_arraysize))
  allocate(hdr_dx(1:hdr_arraysize))
  allocate(hdr_ncellperlevel(1:hdr_arraysize))
  hdr_ilevel = 0
  hdr_dx = 0
  hdr_ncellperlevel = 0

  iovolume = (ABS(xmax-xmin)*ABS(ymax-ymin)*ABS(zmax-zmin))
  if (adaptive_iogroupsize) IOGROUPSIZESAMPLE = nint(iovolume*real(IOGROUPSIZE, kind=8))
  if (IOGROUPSIZESAMPLE .LE. 1) IOGROUPSIZESAMPLE = 1
  if (IOGROUPSIZESAMPLE .GE. IOGROUPSIZE) IOGROUPSIZESAMPLE = IOGROUPSIZE
  if (verbose) write(*,*) 'IOGROUPSIZESAMPLE = ',IOGROUPSIZESAMPLE, 'iovolume = ', iovolume

  ! Wait for the token
  if(sampletoken) then 
    if (mod(myid-1,IOGROUPSIZESAMPLE)/=0) then
      call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                  & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
    end if
  endif
  
  !=========================================================!
  
  ! Loop over levels
  ncell_total = 0
  ihdr = 0
  reduce_infolocal = 0
  reduce_infoglobal = 0
  is_opened = .FALSE.

  do ilevel=sample_levelmin,sample_levelmax

    !*******************************************************!
    ! Initialization of counters
    ncell_level = 0
    ivector = 0
    istride = 0
    icell_in = 0 
    icell_out = 0 
    icell_wr = 0
    ncell_in = 0
    ncell_out = 0 
    ncell_wr = 0
    
    ! Initialization of arrays
    xc_in = 0.D0
    var_in = 0.D0
    bool_in = 0
    xc_out = 0.D0
    var_out = 0.D0
    bool_out = 0
    xc_wr = 0.D0
    var_wr = 0.D0
    bool_wr = 0
    
    ! Index
    ihdr = ihdr+1

    
    !*******************************************************!
    ! Level constants
    dx=0.5D0**ilevel ! 
    dxcell=dx ! need to check (maybe dxcell = (1/2.)*dx
    dxon2 = 0.5D0*dx

    ! Set position of cell centers relative to grid center
    do ind=1,twotondim
      iz=(ind-1)/4
      iy=(ind-1-4*iz)/2
      ix=(ind-1-2*iy-4*iz)
      if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
      if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
      if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
    end do
    
    !*******************************************************!
    
    ! Loop over myid grids by vectorgrid sweeps
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvectorgrid
      ngrid=MIN(nvectorgrid,ncache-igrid+1)
      do i=1,ngrid
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         !ind_grid(i)=active(ilevel)%pcomm%igrid(igrid+i-1) RAMSES-LC line CBH_LC 16-02-2021
      end do
      
      ! Loop over cells
      ivector = 0 ! ivector = icell_in
      do i=1,ngrid
        do ind=1,twotondim
        
          ! Increment ivector
          ivector=ivector+1

          ! Gather cell indices
          iskip=ncoarse+(ind-1)*ngridmax
          ind_cell(ivector)=iskip+ind_grid(i)
          
          ! Compute cell center
          do idim=1,ndim
            xc_in(ivector,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
          end do
          
          ! Compute var :
          var_in(ivector,1)=f(ind_cell(ivector),1)
          var_in(ivector,2)=f(ind_cell(ivector),2)
          var_in(ivector,3)=f(ind_cell(ivector),3)
          var_in(ivector,4)=phi(ind_cell(ivector))
          var_in(ivector,5)=rho(ind_cell(ivector))

          ! CBH_LC
          ! Notice that gradients need to be copied from f() at the correct igrp sequence in move or synchro
          ! At igrp = 10, f() = - grad_i(b)
          !var_in(ivector,1 )=gr_pot(ind_cell(ivector),5) ! Phi
          !var_in(ivector,2 )=gr_pot(ind_cell(ivector),6) ! Xi 
          !var_in(ivector,3 )=gr_pot(ind_cell(ivector),7) - f(ind_cell(ivector),1)  ! \beta^i = B^i + grad_i(b)
          !var_in(ivector,4 )=gr_pot(ind_cell(ivector),8) - f(ind_cell(ivector),2)
          !var_in(ivector,5 )=gr_pot(ind_cell(ivector),9) - f(ind_cell(ivector),3)
          !var_in(ivector,6 )=gr_pot_grad(ind_cell(ivector),1 )    ! grad_i Phi
          !var_in(ivector,7 )=gr_pot_grad(ind_cell(ivector),2 )
          !var_in(ivector,8 )=gr_pot_grad(ind_cell(ivector),3 )
          !var_in(ivector,9 )=gr_pot_grad(ind_cell(ivector),4 )    ! grad_i Xi
          !var_in(ivector,10)=gr_pot_grad(ind_cell(ivector),5 )
          !var_in(ivector,11)=gr_pot_grad(ind_cell(ivector),6 )
          !var_in(ivector,12)=gr_pot_grad(ind_cell(ivector),7 )    ! grad_i \beta^j
          !var_in(ivector,13)=gr_pot_grad(ind_cell(ivector),8 )
          !var_in(ivector,14)=gr_pot_grad(ind_cell(ivector),9 )
          !var_in(ivector,15)=gr_pot_grad(ind_cell(ivector),10)
          !var_in(ivector,16)=gr_pot_grad(ind_cell(ivector),11)
          !var_in(ivector,17)=gr_pot_grad(ind_cell(ivector),12)
          !var_in(ivector,18)=gr_pot_grad(ind_cell(ivector),13)
          !var_in(ivector,19)=gr_pot_grad(ind_cell(ivector),14)
          !var_in(ivector,20)=gr_pot_grad(ind_cell(ivector),15)
          
          ! END CBH_LC
  
          ! Compute bool
          bool_in(ivector,1) = son(ind_cell(ivector))
          
        end do 
      end do ! end of loop over cells : arrays of size ivector<=nvector of data have been generated
      
      !*******************************************************!
      
      ! Perfom selection on cells
      ncell_in = ivector
      icell_out = 0
      do icell_in=1,ncell_in
        if ((((xc_in(icell_in, 1)+dxon2) .GT. xmin) .AND. ((xc_in(icell_in, 1)-dxon2) .LT. xmax)) .AND. &
            (((xc_in(icell_in, 2)+dxon2) .GT. ymin) .AND. ((xc_in(icell_in, 2)-dxon2) .LT. ymax)) .AND. &
            (((xc_in(icell_in, 3)+dxon2) .GT. zmin) .AND. ((xc_in(icell_in, 3)-dxon2) .LT. zmax))) then

          icell_out=icell_out+1
          do idim=1,ndim
            xc_out(icell_out,idim)=xc_in(icell_in,idim)
          end do
          do ivar=1,nvar
            var_out(icell_out,ivar)=var_in(icell_in,ivar)
          end do
          do ibool=1,nbool
            bool_out(icell_out,ibool)=bool_in(icell_in,ibool)
          end do
        end if
      end do
      ncell_out = icell_out
      
      !*******************************************************!
      
      ! Transfer out arrays of size ncell_out in write arrays of size icell_wr
      if(ncell_out.GT.0) then
        call transfer_grav(nvar,nbool,ncell_out,icell_wr,ncell_level,ncell_total, &
                                & xc_out,var_out,bool_out, &
                                & xc_wr,var_wr,bool_wr)
      end if
      
      !*******************************************************!
      
      ! If nstride has been reached, write arrays of size icell_wr and update size
      if (icell_wr.GE.nstride) then
        call write_grav(.false., is_opened, ilun, filename_dat, &
                            & xc_wr, nvar, var_wr, nbool, bool_wr, &
                            & icell_wr, sample_nlevel, reduce_infolocal)
      end if
      
      !*******************************************************!
      
    end do ! loop over grids
    !*******************************************************!
    
    ! If some cells lasts in write vector, write them
    
    if (icell_wr.GT.0) then
      call write_grav(.true., is_opened, ilun, filename_dat, &
                         & xc_wr, nvar, var_wr, nbool, bool_wr, &
                         & icell_wr, sample_nlevel, reduce_infolocal)
    end if

  hdr_ilevel(ihdr) = ilevel
  hdr_dx(ihdr) = dx
  hdr_ncellperlevel(ihdr) = ncell_level

  end do ! loop over levels
  if (is_opened) close(ilun)
  
  !=========================================================!

  ! Write hdr files
  call write_hdr(ncell_total, filename_hdr, hdr_arraysize, hdr_ilevel, hdr_dx, hdr_ncellperlevel)
  
  !=========================================================!

  ! Send the token
  if(sampletoken) then 
    if(mod(myid,IOGROUPSIZESAMPLE)/=0 .and.(myid.lt.ncpu))then
      dummy_io=1
      call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
           & MPI_COMM_WORLD,info)
    end if
  endif

  !=========================================================!
  
  ! Reduce
  reduce_infolocal(2) = ncell_total
  call MPI_REDUCE(reduce_infolocal, reduce_infoglobal, 2, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  ! Write info file
  if (myid.EQ.1) then
    call write_infosamplegrav(ncell_total, aexp, sample_nlevel, sample_levelmin, sample_levelmax, &
                            & xmin, xmax, ymin, ymax, zmin, zmax, &
                            & reduce_infoglobal, filename_info)
  end if
  
  !=========================================================!
  
  ! Finalize
  deallocate(hdr_ilevel)
  deallocate(hdr_dx)
  deallocate(hdr_ncellperlevel)
  
  !=========================================================!
  
end subroutine extract_samplegrav
!===========================================================================



!===========================================================================
! TRANSFER TO WRITE ARRAYS
!===========================================================================
subroutine transfer_grav(nvar,nbool,ncell_out,icell_wr,ncell_level,ncell_total, &
                              & xc_out,var_out,bool_out, &
                              & xc_wr,var_wr,bool_wr)

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
  real(kind=8),dimension(1:27*nvector,1:ndim)::xc_out
  real(kind=8),dimension(1:27*nvector,1:nvar)::var_out
  integer,dimension(1:27*nvector,1:nbool)::bool_out
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim)::xc_wr
  real(kind=4),dimension(1:nstride+27*nvector,1:nvar)::var_wr
  integer,dimension(1:nstride+27*nvector,1:nbool)::bool_wr
  
  ! Variables locales
  integer::icell
  integer::idim
  integer::ivar
  integer::ibool
  
  ! Copy loops
  do icell=1,ncell_out
    do idim=1,ndim
      xc_wr(icell_wr+icell,idim)=xc_out(icell,idim)
    end do
    do ivar=1,nvar
      var_wr(icell_wr+icell,ivar)=var_out(icell,ivar)
    end do
    do ibool=1,nbool
      bool_wr(icell_wr+icell,ibool)=bool_out(icell,ibool)
    end do
  end do
  icell_wr=icell_wr+ncell_out
  ncell_level=ncell_level+ncell_out
  ncell_total=ncell_total+ncell_out
                        
end subroutine transfer_grav
!===========================================================================



!===========================================================================
! WRITE FILE
!===========================================================================
subroutine write_grav(is_last, is_opened, ilun, filename, &
                        & xc_wr, nvar, var_wr, nbool, bool_wr, &
                        & istride, cone_nlevel, reduce_infolocal)
                        
  ! Initialization
  use amr_commons, ONLY: ndim, nstride, nvector, ncpu
  implicit none
  
  ! Variables from input
  logical::is_last                                            ! last writing
  logical::is_opened                                          ! if the file is already opened
  integer::ilun                                               ! unit for writing
  character(LEN=200)::filename                                 ! file location
  real(kind=4),dimension(1:nstride+27*nvector,1:ndim)::xc_wr  ! center of cells array
  integer::nvar                                               ! number of variables in var array
  real(kind=4),dimension(1:nstride+27*nvector,1:nvar)::var_wr ! variable array
  integer::nbool                                              ! number of variables in bool array
  integer,dimension(1:nstride+27*nvector,1:nbool)::bool_wr    ! boolean array
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
    wrsize = istride
  else
    wrsize = nstride
  endif
  
  ! Write the beginning of the file
  if(.NOT.is_opened) then
    open(ilun,file=TRIM(filename),form='unformatted')
    rewind(ilun)  
    write(ilun)ncpu
    write(ilun)nstride
    write(ilun)cone_nlevel
    is_opened=.TRUE.
  endif
  
  ! Write content
  write(ilun)xc_wr(1:wrsize,1)
  write(ilun)var_wr(1:wrsize,1)
  write(ilun)xc_wr(1:wrsize,2)
  write(ilun)var_wr(1:wrsize,2)
  write(ilun)xc_wr(1:wrsize,3)
  write(ilun)var_wr(1:wrsize,3)
  write(ilun)var_wr(1:wrsize,4)
  write(ilun)var_wr(1:wrsize,5)
  write(ilun)bool_wr(1:wrsize,1)

  ! Copy content
  if(.NOT.is_last) then
    do i=1,istride-nstride
      do idim=1,ndim
        xc_wr(i,idim)=xc_wr(i+nstride,idim)
      end do
      do ivar=1,nvar
        var_wr(i,ivar)=var_wr(i+nstride,ivar)
      end do
      do ibool=1,nbool
        bool_wr(i,ibool)=bool_wr(i+nstride,ibool)
      end do
    end do
    istride=istride-nstride
  endif
  
  if(reduce_infolocal(1).EQ.0) then
    reduce_infolocal(1) = 1
  endif

end subroutine write_grav
!===========================================================================



!===========================================================================
! WRITE HDR FILE
!===========================================================================
subroutine write_hdr(ncell_total, filename, hdr_arraysize, hdr_ilevel, hdr_dx, hdr_ncellperlevel)
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
    filename_tmp = filename
    open(ilun,file=TRIM(filename_tmp),form='unformatted')
    rewind(ilun)
    write(ilun)ncpu
    write(ilun)nstride
    write(ilun)ncell_total
    write(ilun)hdr_ncellperlevel
    close(ilun)
    if (conegrav_formattedhdr) then
      filename_tmp = TRIM(filename)//'_formatted'
      open(ilun,file=TRIM(filename_tmp),form='formatted')
      rewind(ilun)
      write(ilun, *) ncpu
      write(ilun, *) nstride
      write(ilun, *) ncell_total
      write(ilun, *) hdr_ncellperlevel
      close(ilun)
    endif
  endif
  
end subroutine write_hdr
!===========================================================================



!===========================================================================
! WRITE INFO
!===========================================================================
subroutine write_infoconegrav(ncell_tot, conegravid, conegravzlim, observer_x, observer_y, observer_z, observer_redshift, &
                            & nlevel, nlevel_min, nlevel_max, &
                            & amax, amin, zmax, zmin, dmax, dmin, dtol, &
                            & reduce_infoglobal, filename, &
                            & is_fullsky, input_thetay, input_thetaz, input_theta, input_phi,&
                            aendconem2_info,aendconem1_info,aendcone_info,&
                            zendconem2_info,zendconem1_info,zendcone_info, &
                            dendconem2_info,dendconem1_info,dendcone_info,future)
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
  isfullsky = 0
  if (is_fullsky) then
    isfullsky = 1
  else
    isfullsky = 0
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
!===========================================================================



!===========================================================================
! WRITE INFO
!===========================================================================
subroutine write_infosamplegrav(ncell_tot, aexp_sample, nlevel, nlevel_min, nlevel_max, &
                        & xmin, xmax, ymin, ymax, zmin, zmax, &
                        & reduce_infoglobal, filename)
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
!===========================================================================



!===========================================================================
! ROTATION MATRIX (CONE NARROW)
!===========================================================================
subroutine compute_rotation_matrix_grav(thetashiftrad,phishiftrad,rot,rotm1)
  implicit none
  real(kind=8) :: thetashiftrad,phishiftrad
  real(kind=8) :: rot(3,3),rotm1(3,3)
  integer :: i,j
  rot(1,1) = cos(thetashiftrad)*cos(phishiftrad)
  rot(1,2) = cos(thetashiftrad)*sin(phishiftrad)
  rot(1,3) = -sin(thetashiftrad)
  rot(2,1) = -sin(phishiftrad)
  rot(2,2) = cos(phishiftrad)
  rot(2,3) = 0.0d0
  rot(3,1) = cos(phishiftrad)*sin(thetashiftrad)
  rot(3,2) = sin(phishiftrad)*sin(thetashiftrad)
  rot(3,3) = cos(thetashiftrad)
  do j=1,3
     do i=1,3
        rotm1(i,j)=rot(j,i)
     enddo
  enddo
end subroutine compute_rotation_matrix_grav
!===========================================================================



!===========================================================================
! MINIMUM POLYGON (CONE NARROW)
!===========================================================================
subroutine compute_minimum_polygon_grav(x1,x2,thetayrad,thetazrad,sl)
  implicit none
  real(kind=8) :: x1,x2,thetayrad,thetazrad,sl(3,8)
  real(kind=8) :: r(3),axis(3)
  sl(1,1:4)=x1/sqrt(1.0d0+tan(thetayrad)**2+tan(thetazrad)**2)
  sl(2,1)=-sl(1,1)*tan(thetayrad)
  sl(3,1)=-sl(1,1)*tan(thetazrad)
  sl(2,2)= sl(1,1)*tan(thetayrad)
  sl(3,2)=-sl(1,1)*tan(thetazrad)
  sl(2,3)=-sl(1,1)*tan(thetayrad)
  sl(3,3)= sl(1,1)*tan(thetazrad)
  sl(2,4)= sl(1,1)*tan(thetayrad)
  sl(3,4)= sl(1,1)*tan(thetazrad)
  sl(1,5:8)=x2
  sl(2,5)=-x2*tan(thetayrad)
  sl(3,5)=-x2*tan(thetazrad)
  sl(2,6)= x2*tan(thetayrad)
  sl(3,6)=-x2*tan(thetazrad)
  sl(2,7)=-x2*tan(thetayrad)
  sl(3,7)= x2*tan(thetazrad)
  sl(2,8)= x2*tan(thetayrad)
  sl(3,8)= x2*tan(thetazrad)
end subroutine compute_minimum_polygon_grav
!===========================================================================



!===========================================================================
! REPLICA (CONE NARROW)
!===========================================================================
subroutine compute_replica_grav(thetayrad,thetazrad,dist1,dist2,observer,Lbox,rot, &
     &                           nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp)
  implicit none
  real(kind=8) :: thetayrad,thetazrad,observer(3),Lbox,rot(3,3),dist1,dist2
  integer :: nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
  integer :: myint
  real(kind=8) :: sl(3,8),slfr(3)
  real(kind=8) :: xplmin,xplmax,yplmin,yplmax,zplmin,zplmax
  integer :: i,j
  call compute_minimum_polygon_grav(dist1,dist2,thetayrad,thetazrad,sl)
  do j=1,8
     do i=1,3
        slfr(i)=sl(1,j)*rot(1,i) &
             & +sl(2,j)*rot(2,i) &
             & +sl(3,j)*rot(3,i)
     enddo
     if (j.eq.1) then
        xplmin=slfr(1)
        xplmax=xplmin
        yplmin=slfr(2)
        yplmax=yplmin
        zplmin=slfr(3)
        zplmax=zplmin
     else
        xplmin=min(xplmin,slfr(1))
        xplmax=max(xplmax,slfr(1))
        yplmin=min(yplmin,slfr(2))
        yplmax=max(yplmax,slfr(2))
        zplmin=min(zplmin,slfr(3))
        zplmax=max(zplmax,slfr(3))
     endif
  enddo
  nrepxm=myint((xplmin+observer(1))/Lbox)
  nrepxp=myint((xplmax+observer(1))/Lbox)
  nrepym=myint((yplmin+observer(2))/Lbox)
  nrepyp=myint((yplmax+observer(2))/Lbox)
  nrepzm=myint((zplmin+observer(3))/Lbox)
  nrepzp=myint((zplmax+observer(3))/Lbox)   
end subroutine compute_replica_grav
!===========================================================================
!
