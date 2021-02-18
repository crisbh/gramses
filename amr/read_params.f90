! CBH_LC 17-02-2021
! Following RAMSES-LC, there are two different read_params subroutines.
! STANDARD VERSION
#ifndef WITHBCASTPARAMS
subroutine read_params
  use amr_commons
  use pm_parameters
  use poisson_parameters
  use gr_parameters
  use hydro_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,levelmax
  character(LEN=80)::infile, info_file
  character(LEN=80)::cmdarg
  character(LEN=5)::nchar
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok, info_ok
  integer,parameter::tag=1134
#ifndef WITHOUTMPI
  integer::dummy_io,ierr,info2
#endif

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/run_params/clumpfind,cosmo,pic,sink,lightcone,poisson,gr,gr_newtonian,hydro,rt,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering &
       & ,bisec_tol,static,geom,overload,cost_weighting,aton,nrestart_quad,restart_remap &
       & ,static_dm,static_gas,static_stars,convert_birth_times,use_proper_time,remap_pscalar
  namelist/output_params/noutput,foutput,fbackup,aout,tout &
       & ,tend,delta_tout,aend,delta_aout,gadget_output,walltime_hrs,minutes_dump
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen,nlevel_collapse
  namelist/poisson_params/epsilon,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax
  namelist/lightcone_params/thetay_cone,thetaz_cone,zmax_cone
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame,movie_vars &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only_frame &
       & ,imovout,imov,tstartmov,astartmov,tendmov,aendmov,proj_axis,movie_vars_txt &
       & ,theta_camera,phi_camera,dtheta_camera,dphi_camera,focal_camera,dist_camera,ddist_camera &
       & ,perspective_camera,smooth_frame,shader_frame,tstart_theta_camera,tstart_phi_camera &
       & ,tend_theta_camera,tend_phi_camera,method_frame,varmin_frame,varmax_frame


  ! CBH_LC 16-02-2021
  namelist/utils_params/withiocoarse,map_full,fmap_full,proj_map_full,nx_map_full,ny_map_full,map_zoom,fmap_zoom,proj_map,nx_map,ny_map,&
       &xmin_map,xmax_map,ymin_map,ymax_map,zmin_map,zmax_map,sample_full,fsample_full,nsample_sample_full,sample_zoom,samplegrav_zoom,fsample_zoom,&
       &nsample_sample,xmin_sample,xmax_sample,ymin_sample,ymax_sample,zmin_sample,zmax_sample,&
       &cone_full,cone_narrow,conegrav_full,conegrav_narrow,&
       &conefull_number,conefull_ok,conefull_id,conefull_observer_x,conefull_observer_y,conefull_observer_z,conefull_observer_redshift,conefull_zmax,&
       &conenarrow_number,conenarrow_ok,conenarrow_id,conenarrow_observer_x,conenarrow_observer_y,conenarrow_observer_z,conenarrow_observer_redshift,conenarrow_zmax,&
       &conenarrow_thetay,conenarrow_thetaz,conenarrow_theta,conenarrow_phi 


  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Careful with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif
  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  if(myid==1)then
  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 3.0                            '
  write(*,*)'       written by Romain Teyssier (University of Zurich)       '
  write(*,*)'               (c) CEA 1999-2007, UZH 2008-2014                '
  write(*,*)' '
  write(*,'(" Working with nproc = ",I4," for ndim = ",I1)')ncpu,ndim
  ! Check nvar is not too small
#ifdef SOLVERhydro
  write(*,'(" Using solver = hydro with nvar = ",I2)')nvar
  if(nvar<ndim+2)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')ndim+2
     call clean_stop
  endif
#endif
#ifdef SOLVERmhd
  write(*,'(" Using solver = mhd with nvar = ",I2)')nvar
  if(nvar<8)then
     write(*,*)'You should have: nvar>=8'
     write(*,'(" Please recompile with -DNVAR=8")')
     call clean_stop
  endif
#endif

  !Write I/O group size information
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '
  if(IOGROUPSIZE>0) write(*,*)'IOGROUPSIZE=',IOGROUPSIZE
  if(IOGROUPSIZECONE>0) write(*,*)'IOGROUPSIZECONE=',IOGROUPSIZECONE
  if(IOGROUPSIZEREP>0) write(*,*)'IOGROUPSIZEREP=',IOGROUPSIZEREP
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = iargc()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call clean_stop
  END IF
  CALL getarg(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------

  ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)then
        write(*,*)'File '//TRIM(infile)//' does not exist'
     endif
     call clean_stop
  end if

  !-------------------------------------------------
  ! Default passive scalar map
  !-------------------------------------------------
#if NVAR>NDIM+2
  allocate(remap_pscalar(1:nvar-(ndim+2)))
  do i=1,nvar-(ndim+2)
     remap_pscalar(i) = i+ndim+2
  enddo
#endif

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=lightcone_params,END=83)
83 continue
  rewind(1)
  read(1,NML=movie_params,END=82)
82 continue
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue

  ! CBH_LC 16-02-2021
  rewind(1)
  read(1,NML=utils_params,END=80)
80 continue

  if(debugreadparams) then
     if(myid==5) then
        write(42,*) cosmo,pic,sink,poisson,hydro,verbose,debug &
             & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering,bisec_tol,static,geom,overload,noutput,foutput,fbackup,aout,tout,&!output_mode,
             & levelmin,levelmax,ngridmax,ngridtot &
             & ,npartmax,nparttot, & !nsinkmax,nsinktot,
             & nexpand,boxlen, epsilon,gravity_type,gravity_params,cg_levelmin,cic_levelmax, & !quint,file_quint,exact_aend,& vfact_grafic,&
             & withiocoarse,map_full,fmap_full,proj_map_full,nx_map_full,ny_map_full,map_zoom,fmap_zoom,proj_map,nx_map,ny_map,xmin_map,xmax_map,ymin_map,ymax_map,zmin_map,&
             &zmax_map,sample_full,fsample_full,nsample_sample_full,sample_zoom,samplegrav_zoom,fsample_zoom,nsample_sample,xmin_sample,xmax_sample,ymin_sample,ymax_sample,zmin_sample,&
             &zmax_sample,cone_full,cone_narrow,conegrav_full,conegrav_narrow,conefull_number,conefull_ok,conefull_id,conefull_observer_x,conefull_observer_y,&
             &conefull_observer_z,conefull_observer_redshift,conefull_zmax,conenarrow_number,conenarrow_ok,conenarrow_id,conenarrow_observer_x,&
             &conenarrow_observer_y,conenarrow_observer_z,conenarrow_observer_redshift,conenarrow_zmax,conenarrow_thetay,conenarrow_thetaz,conenarrow_theta,conenarrow_phi
     endif
  endif

  ! END CBH_LC 16-02-2021


  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------
  if (myid==1 .and. narg == 2) then
     CALL getarg(2,cmdarg)
     read(cmdarg,*) nrestart
  endif

  if (myid==1 .and. nrestart .gt. 0) then
     call title(nrestart,nchar)
     info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
     inquire(file=info_file, exist=info_ok)
     do while(.not. info_ok .and. nrestart .gt. 1)
        nrestart = nrestart - 1
        call title(nrestart,nchar)
        info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
        inquire(file=info_file, exist=info_ok)
     enddo
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(info_ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

  if (nrestart .gt. 0 .and. .not. info_ok) then
     if (myid==1) then
         write(*,*) "Error: Could not find restart file"
     endif
     call clean_stop
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=MIN(int(tend/delta_tout),MAXOUT)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=MIN(int(aend/delta_aout),MAXOUT)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  noutput=MIN(noutput,MAXOUT)
  if(imovout>0) then
     allocate(tmovout(0:imovout))
     allocate(amovout(0:imovout))
     tmovout=1d100
     amovout=1d100
     if(tendmov>0)then
        do i=0,imovout
           tmovout(i)=(tendmov-tstartmov)*dble(i)/dble(imovout)+tstartmov
        enddo
     endif
     if(aendmov>0)then
        do i=0,imovout
           amovout(i)=(aendmov-astartmov)*dble(i)/dble(imovout)+astartmov
        enddo
     endif
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nml_ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=int(ngridtot/int(ncpu,kind=8),kind=4)
     endif
  end if
  if(npartmax==0)then
     npartmax=int(nparttot/int(ncpu,kind=8),kind=4)
  endif
  if(myid>1)verbose=.false.
  if(sink.and.(.not.pic))then
     pic=.true.
  endif
  !if(clumpfind.and.(.not.pic))then
  !   pic=.true.
  !endif
  !if(pic.and.(.not.poisson))then
  !   poisson=.true.
  !endif

  call read_hydro_params(nml_ok)
#ifdef RT
  call rt_read_hydro_params()
#endif
#if NDIM==3
  if (sink)call read_sink_params
  if (clumpfind .or. sink)call read_clumpfind_params
#endif
  if (movie)call set_movie_vars

  close(1)

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

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
  end do

  if(.not.cosmo)then
     use_proper_time=.false.
     convert_birth_times=.false.
  endif

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params

! CBH_LC 17-02-2021
! BCAST VERSION - might require some checks
#else
subroutine read_params
  use amr_commons
  use pm_parameters
  use poisson_parameters
  use gr_parameters
  use hydro_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,levelmax
  character(LEN=80)::infile, info_file
  character(LEN=80)::cmdarg
  character(LEN=5)::nchar
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok, info_ok
  integer,parameter::tag=1134
#ifndef WITHOUTMPI
  integer::dummy_io,ierr,info2
#endif

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/run_params/clumpfind,cosmo,pic,sink,lightcone,poisson,gr,gr_newtonian,hydro,rt,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering &
       & ,bisec_tol,static,geom,overload,cost_weighting,aton,nrestart_quad,restart_remap &
       & ,static_dm,static_gas,static_stars,convert_birth_times,use_proper_time,remap_pscalar
  namelist/output_params/noutput,foutput,fbackup,aout,tout &
       & ,tend,delta_tout,aend,delta_aout,gadget_output,walltime_hrs,minutes_dump
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen,nlevel_collapse
  namelist/poisson_params/epsilon,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax
  namelist/lightcone_params/thetay_cone,thetaz_cone,zmax_cone
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame,movie_vars &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only_frame &
       & ,imovout,imov,tstartmov,astartmov,tendmov,aendmov,proj_axis,movie_vars_txt &
       & ,theta_camera,phi_camera,dtheta_camera,dphi_camera,focal_camera,dist_camera,ddist_camera &
       & ,perspective_camera,smooth_frame,shader_frame,tstart_theta_camera,tstart_phi_camera &
       & ,tend_theta_camera,tend_phi_camera,method_frame,varmin_frame,varmax_frame


  ! CBH_LC 16-02-2021
  namelist/utils_params/withiocoarse,map_full,fmap_full,proj_map_full,nx_map_full,ny_map_full,map_zoom,fmap_zoom,proj_map,nx_map,ny_map,&
       &xmin_map,xmax_map,ymin_map,ymax_map,zmin_map,zmax_map,sample_full,fsample_full,nsample_sample_full,sample_zoom,samplegrav_zoom,fsample_zoom,&
       &nsample_sample,xmin_sample,xmax_sample,ymin_sample,ymax_sample,zmin_sample,zmax_sample,&
       &cone_full,cone_narrow,conegrav_full,conegrav_narrow,&
       &conefull_number,conefull_ok,conefull_id,conefull_observer_x,conefull_observer_y,conefull_observer_z,conefull_observer_redshift,conefull_zmax,&
       &conenarrow_number,conenarrow_ok,conenarrow_id,conenarrow_observer_x,conenarrow_observer_y,conenarrow_observer_z,conenarrow_observer_redshift,conenarrow_zmax,&
       &conenarrow_thetay,conenarrow_thetaz,conenarrow_theta,conenarrow_phi 


  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Careful with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif
  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  if(myid==1)then
  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 3.0                            '
  write(*,*)'       written by Romain Teyssier (University of Zurich)       '
  write(*,*)'               (c) CEA 1999-2007, UZH 2008-2014                '
  write(*,*)' '
  write(*,'(" Working with nproc = ",I4," for ndim = ",I1)')ncpu,ndim
  ! Check nvar is not too small
#ifdef SOLVERhydro
  write(*,'(" Using solver = hydro with nvar = ",I2)')nvar
  if(nvar<ndim+2)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')ndim+2
     call clean_stop
  endif
#endif
#ifdef SOLVERmhd
  write(*,'(" Using solver = mhd with nvar = ",I2)')nvar
  if(nvar<8)then
     write(*,*)'You should have: nvar>=8'
     write(*,'(" Please recompile with -DNVAR=8")')
     call clean_stop
  endif
#endif

  !Write I/O group size information
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '
  if(IOGROUPSIZE>0) write(*,*)'IOGROUPSIZE=',IOGROUPSIZE
  if(IOGROUPSIZECONE>0) write(*,*)'IOGROUPSIZECONE=',IOGROUPSIZECONE
  if(IOGROUPSIZEREP>0) write(*,*)'IOGROUPSIZEREP=',IOGROUPSIZEREP
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = iargc()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call clean_stop
  END IF
  CALL getarg(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

!-------------------------------------------------!
!!!!!!!!!!!!!! CBH_LC 17-02-2021 !!!!!!!!!!!!!!!!!!
!-------------------------------------------------!

  ! Memory allocation for the input parameters broadcast buffer
  ! 16*4                   : 16x         Logical
  !                           4x         Character           
  ! 3*bit_size(ngridtot)/8 :  3x         Integer(kind=8) 
  ! 29*bit_size(nrestart)/8: 29x         Integer(kind=4)
  ! 1*bit_size(nsubcycle)/8:  1x table x Integer(kind=4)
  ! 30*kind(vfact_grafic)  : 30x         Real(kind=dp) 
  !  1*kind(aout)          :  1x table x Real(kind=dp)
  !  1*kind(aout)          :  1x table x Real(kind=dp)


  ! Memory allocation for the input parameters broadcast buffer
  ! 20+maxconefull+maxconenarrow*4                   : 20+maxconefull+maxconenarrowx         Logical
  !                           4x         Character           
  ! 3*bit_size(ngridtot)/8 :  3x         Integer(kind=8) 
  ! 31*bit_size(nrestart)/8: 31x         Integer(kind=4)
  ! size(nsubcycle)*bit_size(nrestart)/8:  size x Integer(kind=4)
  ! MAXCONEFULL*bit_size(nrestart)/8: MAXCONEFULL*Integer(kind=4)
  ! MAXCONENARROW*bit_size(nrestart)/8: MAXCONENARROW*Integer(kind=4)
  !
  ! 26*kind(vfact_grafic)  : 26x         Real(kind=dp) 
  !  1*kind(aout)          :  1x table x Real(kind=dp)
  !  1*kind(aout)          :  1x table x Real(kind=dp)
  !5*MAXCONEFULL*kind(vfact_grafic)
  !9*MAXCONENARROW*kind(vfact_grafic)

!  h_length = (20+MAXCONEFULL+MAXCONENARROW)*4 + 258 + 3*bit_size(ngridtot)/8 + 31*bit_size(nrestart)/8 + size(nsubcycle)*bit_size(nrestart)/8+MAXCONEFULL*bit_size(nrestart)/8 +MAXCONENARROW*bit_size(nrestart)/8&
!       & + 26*kind(vfact_grafic) + (size(aout)+size(tout))*kind(aout)+(5*MAXCONEFULL+9*MAXCONENARROW)*kind(vfact_grafic)
  
  ! CBH_LC 17-02-2021
  ! CBH: Removed vfact_grafic terms as not under use. Not sure what should replace them.
  h_length = (20+MAXCONEFULL+MAXCONENARROW)*4 + 258 + 3*bit_size(ngridtot)/8 + 31*bit_size(nrestart)/8 + size(nsubcycle)*bit_size(nrestart)/8+MAXCONEFULL*bit_size(nrestart)/8 +MAXCONENARROW*bit_size(nrestart)/8 + (size(aout)+size(tout))*kind(aout)

  allocate (VBheader(0:h_length-1))
  h_pos = 0

!-------------------------------------------------!
!!!!!!!!!!!!!! CBH_LC 17-02-2021 !!!!!!!!!!!!!!!!!!
!-------------------------------------------------!

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------

  ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

!-------------------------------------------------!
!!!!!!!!!!!!!! CBH_LC 17-02-2021 !!!!!!!!!!!!!!!!!!
!-------------------------------------------------!

  IF(myid==1) THEN

     namelist_file=TRIM(infile)
     INQUIRE(file=infile,exist=nml_ok)
     if(.not. nml_ok)then
        if(myid==1)then
           write(*,*)'File '//TRIM(infile)//' does not exist'
        endif
        call clean_stop
     end if

     !-------------------------------------------------
     ! Default passive scalar map
     !-------------------------------------------------
#if NVAR>NDIM+2
     allocate(remap_pscalar(1:nvar-(ndim+2)))
     do i=1,nvar-(ndim+2)
        remap_pscalar(i) = i+ndim+2
     enddo
#endif

     open(1,file=infile)
     rewind(1)
     read(1,NML=run_params)
     rewind(1)
     read(1,NML=output_params)
     rewind(1)
     read(1,NML=amr_params)
     rewind(1)
     read(1,NML=lightcone_params,END=83)
83    continue
     rewind(1)
     read(1,NML=movie_params,END=82)
82    continue
     rewind(1)
     read(1,NML=poisson_params,END=81)
81    continue

     ! CBH_LC 16-02-2021
     rewind(1)
     read(1,NML=utils_params,END=80)
80    continue


     ! RUN_PARAMS
     Call Mpi_Pack(    cosmo,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(      pic,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(     sink,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(  poisson,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(    hydro,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(  verbose,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(    debug,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( nrestart,  1,             Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( ncontrol,  1,             Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( nstepmax,  1,             Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(nsubcycle,size(nsubcycle), Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   nremap,  1,             Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( ordering,128,           Mpi_Character, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(bisec_tol,  1,    Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   static,  1,             Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(     geom,  1,             Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( overload,  1,             Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)

     ! OUTPUT_PARAMS
     Call Mpi_Pack(    noutput,  1,                  Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(    foutput,  1,                  Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(    fbackup,  1,                  Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(       aout, size(aout), Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(       tout, size(tout), Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
!     Call Mpi_Pack(output_mode,  1,                  Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)

     ! AMR_PARAMS
     Call Mpi_Pack( levelmin,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( levelmax,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( ngridmax,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( ngridtot,  1,         Mpi_Integer8, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( npartmax,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack( nparttot,  1,         Mpi_Integer8, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
!     Call Mpi_Pack( nsinkmax,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
!     Call Mpi_Pack( nsinktot,  1,         Mpi_Integer8, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(  nexpand,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   boxlen,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)

     ! POISSON_PARAMS
     Call Mpi_Pack(       epsilon,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(  gravity_type,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(gravity_params, 10, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   cg_levelmin,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(  cic_levelmax,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)

     ! COSMO_PARAMS
!     Call Mpi_Pack(       quint,  1,          Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
!     Call Mpi_Pack(  file_quint,128,        Mpi_Character, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
!     Call Mpi_Pack(  exact_aend,  1,          Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
!     Call Mpi_Pack(vfact_grafic,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)

     ! UTILS_PARAMS
     Call Mpi_Pack(          withiocoarse,  1, Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           map_full,  1,          Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(          fmap_full,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(      proj_map_full,  1, Mpi_Character, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        nx_map_full,  1, Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        ny_map_full,  1, Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           map_zoom,  1, Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(          fmap_zoom,  1, Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           proj_map,  1, Mpi_Character, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(             nx_map,  1, Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(             ny_map,  1, Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           xmin_map,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           xmax_map,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           ymin_map,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           ymax_map,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           zmin_map,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(           zmax_map,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        sample_full,  1,          Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(       fsample_full,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr) 
     Call Mpi_Pack(nsample_sample_full,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        sample_zoom,  1,          Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(    samplegrav_zoom,  1,          Mpi_Logical, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(       fsample_zoom,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(     nsample_sample,  1,          Mpi_Integer, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        xmin_sample,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        xmax_sample,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        ymin_sample,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        ymax_sample,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        zmin_sample,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        zmax_sample,  1, Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        cone_full              ,  1            , Mpi_Logical         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(        cone_narrow            ,  1            , Mpi_Logical         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(    conegrav_full              ,  1            , Mpi_Logical         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(  conegrav_narrow              ,  1            , Mpi_Logical         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   conefull_number             ,  1            , Mpi_Integer         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   conefull_ok                 ,MAXCONEFULL    , Mpi_Logical         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   conefull_id                 ,MAXCONEFULL    , Mpi_Integer         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conefull_observer_x            ,MAXCONEFULL    , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conefull_observer_y            ,MAXCONEFULL    , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conefull_observer_z            ,MAXCONEFULL    , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conefull_observer_redshift     ,MAXCONEFULL    , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conefull_zmax                  ,MAXCONEFULL    , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   conenarrow_number           ,  1            , Mpi_Integer         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   conenarrow_ok               ,MAXCONENARROW  , Mpi_Logical         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(   conenarrow_id               ,MAXCONENARROW  , Mpi_Integer         , VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_observer_x          ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_observer_y          ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_observer_z          ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_observer_redshift   ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_zmax       ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_thetay     ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_thetaz     ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_theta      ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     Call Mpi_Pack(conenarrow_phi        ,MAXCONENARROW  , Mpi_Double_Precision, VBheader, h_length, h_pos, Mpi_Comm_World, ierr)
     
  ENDIF

  ! Broadcast of input parameters
  Call Mpi_Bcast(VBheader,h_length,Mpi_Packed,0,Mpi_Comm_World,ierr)

  ! Processes with ID != 1 unpack the input parameters
  IF(myid/=1) THEN

     ! RUN_PARAMS
     Call Mpi_Unpack(VBheader, h_length, h_pos,    cosmo,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,      pic,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,     sink,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,  poisson,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,    hydro,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,  verbose,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,    debug,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, nrestart,  1,             Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, ncontrol,  1,             Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, nstepmax,  1,             Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,nsubcycle,size(nsubcycle), Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,   nremap,  1,             Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, ordering,128,           Mpi_Character, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,bisec_tol,  1,    Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,   static,  1,             Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,     geom,  1,             Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, overload,  1,             Mpi_Integer, Mpi_Comm_World, ierr)

     ! OUTPUT_PARAMS
     Call Mpi_Unpack(VBheader, h_length, h_pos,    noutput,  1,                  Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,    foutput,  1,                  Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,    fbackup,  1,                  Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,       aout, size(aout), Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,       tout, size(tout), Mpi_Double_Precision, Mpi_Comm_World, ierr)
!     Call Mpi_Unpack(VBheader, h_length, h_pos,output_mode,  1,                  Mpi_Integer, Mpi_Comm_World, ierr)

     ! AMR_PARAMS
     Call Mpi_Unpack(VBheader, h_length, h_pos, levelmin,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, levelmax,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, ngridmax,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, ngridtot,  1,         Mpi_Integer8, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, npartmax,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos, nparttot,  1,         Mpi_Integer8, Mpi_Comm_World, ierr)
!     Call Mpi_Unpack(VBheader, h_length, h_pos, nsinkmax,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
!     Call Mpi_Unpack(VBheader, h_length, h_pos, nsinktot,  1,         Mpi_Integer8, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,  nexpand,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,   boxlen,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)

     ! POISSON_PARAMS
     Call Mpi_Unpack(VBheader, h_length, h_pos,       epsilon,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,  gravity_type,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,gravity_params, 10, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,   cg_levelmin,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,  cic_levelmax,  1,          Mpi_Integer, Mpi_Comm_World, ierr)

!     ! COSMO_PARAMS
!     Call Mpi_Unpack(VBheader, h_length, h_pos,       quint,  1,          Mpi_Logical, Mpi_Comm_World, ierr)
!     Call Mpi_Unpack(VBheader, h_length, h_pos,  file_quint,128,        Mpi_Character, Mpi_Comm_World, ierr)
!     Call Mpi_Unpack(VBheader, h_length, h_pos,  exact_aend,  1,          Mpi_Logical, Mpi_Comm_World, ierr)
!     Call Mpi_Unpack(VBheader, h_length, h_pos,vfact_grafic,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)

     ! UTILS_PARAMS
     Call Mpi_Unpack(VBheader, h_length, h_pos,          withiocoarse,  1, Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           map_full,  1,          Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,          fmap_full,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,      proj_map_full,  1, Mpi_Character, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        nx_map_full,  1, Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        ny_map_full,  1, Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           map_zoom,  1, Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,          fmap_zoom,  1, Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           proj_map,  1, Mpi_Character, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,             nx_map,  1, Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,             ny_map,  1, Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           xmin_map,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           xmax_map,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           ymin_map,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           ymax_map,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           zmin_map,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,           zmax_map,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        sample_full,  1,          Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,       fsample_full,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,nsample_sample_full,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        sample_zoom,  1,          Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,    samplegrav_zoom,  1,          Mpi_Logical, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,       fsample_zoom,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,     nsample_sample,  1,          Mpi_Integer, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        xmin_sample,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        xmax_sample,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        ymin_sample,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        ymax_sample,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        zmin_sample,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        zmax_sample,  1, Mpi_Double_Precision, Mpi_Comm_World, ierr)     
     Call Mpi_Unpack(VBheader, h_length, h_pos,        cone_full                   ,  1             , Mpi_Logical         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        cone_narrow                 ,  1             , Mpi_Logical         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conegrav_full               ,  1             , Mpi_Logical         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conegrav_narrow             ,  1             , Mpi_Logical         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_number             ,  1             , Mpi_Integer         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_ok                 ,MAXCONEFULL     , Mpi_Logical         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_id                 ,MAXCONEFULL     , Mpi_Integer         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_observer_x         ,MAXCONEFULL     , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_observer_y         ,MAXCONEFULL     , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_observer_z         ,MAXCONEFULL     , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_observer_redshift  ,MAXCONEFULL     , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conefull_zmax               ,MAXCONEFULL     , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_number           ,  1             , Mpi_Integer         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_ok               ,MAXCONENARROW   , Mpi_Logical         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_id               ,MAXCONENARROW   , Mpi_Integer         , Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_observer_x       ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_observer_y       ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_observer_z       ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_observer_redshift,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_zmax             ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_thetay           ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_thetaz           ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_theta            ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)
     Call Mpi_Unpack(VBheader, h_length, h_pos,        conenarrow_phi              ,MAXCONENARROW   , Mpi_Double_Precision, Mpi_Comm_World, ierr)

  ENDIF

  deallocate(VBheader)

  if(debugreadparams) then
     if(myid==5) then
        write(42,*) cosmo,pic,sink,poisson,hydro,verbose,debug &
             & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering,bisec_tol,static,geom,overload,noutput,foutput,fbackup,aout,tout,&!output_mode,
             & levelmin,levelmax,ngridmax,ngridtot &
             & ,npartmax,nparttot, & !nsinkmax,nsinktot,
             & nexpand,boxlen, epsilon,gravity_type,gravity_params,cg_levelmin,cic_levelmax, & !quint,file_quint,exact_aend,& vfact_grafic,&
             & withiocoarse,map_full,fmap_full,proj_map_full,nx_map_full,ny_map_full,map_zoom,fmap_zoom,proj_map,nx_map,ny_map,xmin_map,xmax_map,ymin_map,ymax_map,zmin_map,&
             &zmax_map,sample_full,fsample_full,nsample_sample_full,sample_zoom,samplegrav_zoom,fsample_zoom,nsample_sample,xmin_sample,xmax_sample,ymin_sample,ymax_sample,zmin_sample,&
             &zmax_sample,cone_full,cone_narrow,conegrav_full,conegrav_narrow,conefull_number,conefull_ok,conefull_id,conefull_observer_x,conefull_observer_y,&
             &conefull_observer_z,conefull_observer_redshift,conefull_zmax,conenarrow_number,conenarrow_ok,conenarrow_id,conenarrow_observer_x,&
             &conenarrow_observer_y,conenarrow_observer_z,conenarrow_observer_redshift,conenarrow_zmax,conenarrow_thetay,conenarrow_thetaz,conenarrow_theta,conenarrow_phi
     endif
  endif

  ! END CBH_LC 16-02-2021


  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------
  if (myid==1 .and. narg == 2) then
     CALL getarg(2,cmdarg)
     read(cmdarg,*) nrestart
  endif

  if (myid==1 .and. nrestart .gt. 0) then
     call title(nrestart,nchar)
     info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
     inquire(file=info_file, exist=info_ok)
     do while(.not. info_ok .and. nrestart .gt. 1)
        nrestart = nrestart - 1
        call title(nrestart,nchar)
        info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
        inquire(file=info_file, exist=info_ok)
     enddo
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(info_ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

  if (nrestart .gt. 0 .and. .not. info_ok) then
     if (myid==1) then
         write(*,*) "Error: Could not find restart file"
     endif
     call clean_stop
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=MIN(int(tend/delta_tout),MAXOUT)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=MIN(int(aend/delta_aout),MAXOUT)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  noutput=MIN(noutput,MAXOUT)
  if(imovout>0) then
     allocate(tmovout(0:imovout))
     allocate(amovout(0:imovout))
     tmovout=1d100
     amovout=1d100
     if(tendmov>0)then
        do i=0,imovout
           tmovout(i)=(tendmov-tstartmov)*dble(i)/dble(imovout)+tstartmov
        enddo
     endif
     if(aendmov>0)then
        do i=0,imovout
           amovout(i)=(aendmov-astartmov)*dble(i)/dble(imovout)+astartmov
        enddo
     endif
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nml_ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=int(ngridtot/int(ncpu,kind=8),kind=4)
     endif
  end if
  if(npartmax==0)then
     npartmax=int(nparttot/int(ncpu,kind=8),kind=4)
  endif
  if(myid>1)verbose=.false.
  if(sink.and.(.not.pic))then
     pic=.true.
  endif
  !if(clumpfind.and.(.not.pic))then
  !   pic=.true.
  !endif
  !if(pic.and.(.not.poisson))then
  !   poisson=.true.
  !endif

  call read_hydro_params(nml_ok)
#ifdef RT
  call rt_read_hydro_params()
#endif
#if NDIM==3
  if (sink)call read_sink_params
  if (clumpfind .or. sink)call read_clumpfind_params
#endif
  if (movie)call set_movie_vars

  close(1)

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

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
  end do

  if(.not.cosmo)then
     use_proper_time=.false.
     convert_birth_times=.false.
  endif

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params
#endif
