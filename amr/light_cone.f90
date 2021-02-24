!-----------------------------------------------------------! 
!--------------------- CBH_LC 09-02-2021 -------------------!
!-----------------------------------------------------------! 
subroutine output_cone_fullsky(filename,cone_id,observer_x,observer_y,observer_z,observer_redshift,cone_zmax,future)
  use amr_commons
  use pm_commons
  use mpi
  implicit none
  character(LEN=200)::filename

  integer::ilun,nx_loc,ipout,npout,npart_out
  integer::dummy_io,info
  integer,parameter::tag=1800
  character(LEN=200)::fileloc
  character(LEN=5)::nchar
  real(kind=8),dimension(1:ndim,1:nvector)::pos,vel
    
  integer(kind=8),dimension(1:nvector)::idtab

  real(sp),dimension(1:nvector)::pottab
  real(sp),dimension(1:27*nvector)::potout
  real(sp),dimension(1:nstride+27*nvector)::pot_out
  
  real(sp),dimension(1:ndim,1:nvector)::fparttab
  real(sp),dimension(1:ndim,1:27*nvector)::fpartout
  real(sp),dimension(1:ndim,1:nstride+27*nvector)::fpart_out

  real(kind=8),dimension(1:ndim,1:27*nvector)::posout,velout
  integer(kind=8),dimension(1:27*nvector)::idout
  real(kind=8),dimension(1:27*nvector)::zout
  
  real(sp),dimension(1:nstride+27*nvector,1:ndim)::xp_out,vp_out
  real(sp),dimension(1:nstride+27*nvector)::zp_out
  integer(kind=8),dimension(1:nstride+27*nvector)::id_out

  real(kind=8)::z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8)::observer(3)
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel
  integer::i,ig,ip,npart1
  integer,dimension(1:nvector)::ind_part
  logical::opened
  integer::future
!------------------ MODIF V. REVERDY 2012 ------------------!
  integer::cone_id
  integer::ierr
  real(kind=8)::cone_zmax
  
  real(kind=8)::observer_x
  real(kind=8)::observer_y
  real(kind=8)::observer_z
  real(kind=8)::observer_redshift
  
  real(kind=8) :: dist1cone, dist2cone, disttolcone

  integer(kind=8),dimension(1:2)::reduce_infolocal
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  real(kind=8)::iovolume

  logical::is_fullsky
  real(kind=8)::input_thetay
  real(kind=8)::input_thetaz
  real(kind=8)::input_theta
  real(kind=8)::input_phi

  ! aexp for info file
  real(kind=8)::aendconem2_info,aendconem1_info,aendcone_info,aexpold_info,aexp_info
  !aendconem* -> expansion factor for end of the shell
  !the end of previous shell and the end of previous previous shell (which is the begining of the current shell)
  !when activating overlapping buffer option.
  real(kind=8)::zendconem2_info,zendconem1_info,zendcone_info,zexpold_info,zexp_info
  real(kind=8)::dendconem2_info,dendconem1_info,dendcone_info,dexpold_info,dexp_info

  real(kind=8)::Omega0,OmegaL,OmegaR,coverH0
  real(kind=8)::coord_distance
  !
  !
  !
  !
  !
  reduce_infolocal = 0
  reduce_infoglobal = 0
  iovolume = 1.D0
  
  dist1cone = 0.D0
  dist2cone = 0.D0
  disttolcone = 0.D0

  is_fullsky = .true.
  input_thetay = 0.D0
  input_thetaz = 0.D0
  input_theta = 0.D0
  input_phi = 0.D0

!-----------------------------------------------------------!

  opened=.false.

  if(nstep_coarse.lt.4) return
  if(verbose) write(*,*)'Entering output_part_lightcone_fullsky'
  
  z2=1.D0/aexp_old-1.D0                                                 ! redshift of farther side of shell
  z1=1.D0/aexp    -1.D0                                                 ! redshift of closer  side of shell

!------------------ MODIF V. REVERDY 2011 ------------------!
  if(use_aexp_restart.AND.(nstep_coarse_after_restart.eq.2)) then
     z2=1.0D0/aexp_restart_light_cone-1.0D0
  end if

  if(myid.eq.1) then
     if(verbose) write(*,*) '*************************************************************'
     if(verbose) write(*,*) 'use_aexp_restart=',use_aexp_restart
     if(verbose) write(*,*) 'nstep_coarse_after_restart=',nstep_coarse_after_restart
     if(verbose) write(*,*) 'filename=',filename
     if(use_aexp_restart.AND.(nstep_coarse_after_restart.eq.2)) then
        if(verbose) write(*,*) 'aexp_restart_light_cone=',aexp_restart_light_cone,'    aexp=',aexp
     else
        if(verbose) write(*,*) 'aexp_old=',aexp_old,'    aexp=',aexp
     end if
     if(verbose) write(*,*) 'z2=',z2,'    z1=',z1
     if(verbose) write(*,*) '-------------------------------------------------------------'
  end if
!-----------------------------------------------------------!

  if(cone_overlap) then
     if((aendconem1.lt.aendconem2).or.(aendcone.lt.aendconem1)) print*,'WARNING CONE SHELL RANGE NOT WELL ORDERED, AEXP ',aendconem2,aendconem1,aendcone 

     if(future.eq.1) then
        z2=1.0D0/aendconem1-1.0D0
        z1=1.0D0/aendcone-1.0D0
     else if(future.eq.-1) then
        z2=1.0D0/aendconem2-1.0D0
        z1=1.0D0/aendconem1-1.0D0
     else
        print*,'not implemented cone part future=',future
        stop
     endif

     if(myid.eq.1) then
        print*, ''
        print*, 'TEST CONE PART LIMIT'
        print*, 'nstepcoarse,z2,z1',nstep_coarse,z2,z1
        print*, 'aendconem2,aendconem1,aendcone',aendconem2,aendconem1,aendcone
        print*, ''
     endif
  endif

  !zmax_cone_full=0.25
  !if(z2.gt.zmax_cone_full)return
  if(z2.gt.cone_zmax) return                                            ! shell outside intended redshift range
  if(dabs(z2-z1)<1d-6) return                                           ! shell too thin
  if((z1.lt.observer_redshift).and.(z2.lt.observer_redshift)) return    ! shell completely in the future of observer (VR)
  
  om0in=omega_m                                                         ! Omega_m
  omLin=omega_l                                                         ! Omega_Lambda
  hubin=h0/100.D0                                                       ! h
  Lbox=boxlen_ini/hubin                                                 ! box size in Mpc (not Mpc/h)
  observer(1) = Lbox*observer_x                                         ! observer x coordinate in Mpc (VR)
  observer(2) = Lbox*observer_y                                         ! observer y coordinate in Mpc (VR)
  observer(3) = Lbox*observer_z                                         ! observer z coordinate in Mpc (VR)

  ! Compute iogroupsize
  call perform_my_selection_simple(.true.,.false.,z1,z2,                                   &
       &                           om0in,omLin,hubin,Lbox,                                 &
       &                           observer,observer_redshift,                             &
       &                           pos,vel,idtab,pottab,fparttab,ip,                       &
       &                           posout,velout,idout,zout,potout,fpartout,npout,.false., &
       &                           dist1cone,dist2cone,disttolcone)
  !
  iovolume = ((2.0D0*dist2cone)**3)-(((2.D0*dist1cone)/dsqrt(real(ndim,kind=8)))**3)
  if(adaptive_iogroupsize) IOGROUPSIZECONE = nint(iovolume*real(IOGROUPSIZE,kind=8))
  if(IOGROUPSIZECONE.le.1) IOGROUPSIZECONE = 1
  if(IOGROUPSIZECONE.ge.IOGROUPSIZE) IOGROUPSIZECONE = IOGROUPSIZE
  if(verbose) write(*,*) cone_id,'IOGROUPSIZECONE = ',IOGROUPSIZECONE,'iovolume = ',iovolume

  ! Wait for the token
  if(conetoken) then
     if(mod(myid-1,IOGROUPSIZECONE).ne.0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  endif
  
  ilun=3*ncpu+myid+10                                                   ! file handler
  call title(myid,nchar)                                                ! create file name
  fileloc=TRIM(filename)//TRIM(nchar)//'.dat'                           ! create file name
    
  npart_out=0                                                           ! total number of selected particles
  npout=0                                                               ! number of selected particles from a batch of size nvector
  ipout=0 
  ilevel=levelmin
  !
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)                                           ! first grid on given CPU
     ip=0                                                               ! initialise particle counter
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)                                             ! number of particles in the grid
        if(npart1>0)then        
           ipart=headp(igrid)                                           ! firt particle on grid
           ! Loop over particles
           do jpart=1,npart1
              ip=ip+1                                                   ! increment particle counter
              ind_part(ip)=ipart                                        ! find global particle ID and store in ind_part
              if(ip.eq.nvector)then                                     ! collect data in batches of size nvector
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do i=1,ip
                       pos(idim,i)=xp(ind_part(i),idim)*Lbox            ! particle coordinates
                       vel(idim,i)=vp(ind_part(i),idim)                 ! particle velocities
#ifdef WITHPARTFORCE
                       fparttab(idim,i)=fpart(ind_part(i),idim)         ! particle force components
#endif
                    end do
                 end do
                 do i=1,ip
                    idtab(i)=idp(ind_part(i))                           ! particle ID
#ifdef OUTPUT_PARTICLE_POTENTIAL
                    pottab(i)=ptcl_phi(ind_part(i))                     ! particle potential value
#endif                    
                 end do
                 !
                 !===========================================================================
                 call perform_my_selection_simple(.false.,.false.,z1,z2,                                  &
                      &                           om0in,omLin,hubin,Lbox,                                 &
                      &                           observer,observer_redshift,                             &
                      &                           pos,vel,idtab,pottab,fparttab,ip,                       &
                      &                           posout,velout,idout,zout,potout,fpartout,npout,.false., &
                      &                           dist1cone,dist2cone,disttolcone)
                 !===========================================================================
                 ! 
                 if(npout>0)then                                        ! if some particles get selected
                    do idim=1,ndim
                       do i=1,npout
                          xp_out(ipout+i,idim)=posout(idim,i)/Lbox      ! collect coordinate of selected particles
                          vp_out(ipout+i,idim)=velout(idim,i)           ! collect velocities of selected particles
#ifdef WITHPARTFORCE
                          fpart_out(ipout+i,idim)=fpartout(idim,i)      ! collect the forces on selected particles
#endif

                       end do
                    end do
                    do i=1,npout
                       id_out(ipout+i)=idout(i)                         ! collect IDs of selected particles
                    end do
                    do i=1,npout
                       zp_out(ipout+i)=zout(i)                          ! collect redshifts of selected particles
#ifdef OUTPUT_PARTICLE_POTENTIAL
                       pot_out(ipout+i)=potout(i)                       ! collect potentials on selected particles
#endif                    
                    end do
                    ipout=ipout+npout
                    npart_out=npart_out+npout                           ! increment total number of selected particles
                 endif                                                  ! if(npout>0)
                 ip=0                                                   ! reset particle counter
              end if                                                    ! if(ip.eq.nvector)
              !
              if(ipout.ge.nstride)then                                  ! reached size for writing into file (nstride)
                 if(.not.opened) then                                   ! open file for writing if not already done so
                    open(ilun,file=TRIM(fileloc),form='unformatted')
                    rewind(ilun)  
                    write(ilun) ncpu
                    write(ilun) nstride
                    write(ilun) npart
                    opened=.true.
                 endif
                 !
                 do idim=1,ndim
                    write(ilun)xp_out(1:nstride,idim)                   ! write particle coordinates
                    write(ilun)vp_out(1:nstride,idim)                   ! write particle velocities
                 end do
                 write(ilun)zp_out(1:nstride)                           ! write particle redshifts
                 if(write_cone_idp) write(ilun) id_out(1:nstride)       ! write particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
                 write(ilun)pot_out(1:nstride)                          ! write particle potentials
#endif                    
#ifdef WITHPARTFORCE
                 do idim=1,ndim
                    write(ilun)fpart_out(1:nstride,idim)                ! write particle forces
                 end do
#endif
                 ! there can be more particles than nstride; these are
                 ! not wirtten right now and will be written later
                 do idim=1,ndim                                         ! move unwritten particles forward in arrays
                    do i=1,ipout-nstride
                       xp_out(i,idim)=xp_out(i+nstride,idim)            ! coordinates
                       vp_out(i,idim)=vp_out(i+nstride,idim)            ! velocities
#ifdef WITHPARTFORCE
                       fpart_out(i,idim)=fpart_out(i+nstride,idim)      ! forces
#endif
                    end do
                 end do
                 do i=1,ipout-nstride
                    id_out(i)=id_out(i+nstride)                         ! particle IDs
                 end do
                 do i=1,ipout-nstride
                    zp_out(i) =zp_out(i+nstride)                        ! particle redshifts
#ifdef OUTPUT_PARTICLE_POTENTIAL
                    pot_out(i)=pot_out(i+nstride)                       ! particle potentials
#endif                    
                 end do
                 !
                 ipout=ipout-nstride                                    ! reset quantity used to count particles
              endif                                                     ! if(ipout.ge.nstride)
              ipart=nextp(ipart)  ! Go to next particle
           end do
           ! End loop over particles           
        end if
        igrid=next(igrid)                                               ! go to next grid
     end do
     ! End loop over grids

     if(ip>0)then                                                       ! remaining particles on grids
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ip
              pos(idim,i)=xp(ind_part(i),idim)*Lbox                     ! particle coordinates
              vel(idim,i)=vp(ind_part(i),idim)                          ! particle velocities
#ifdef WITHPARTFORCE
              fparttab(idim,i)=fpart(ind_part(i),idim)                  ! particle force components
#endif
              
           end do
        end do
        do i=1,ip
           idtab(i) =idp(ind_part(i))                                   ! particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
           pottab(i)=ptcl_phi(ind_part(i))                              ! particle potential values
#endif                    
        end do
        !===========================================================================
        call perform_my_selection_simple(.false.,.false.,z1,z2,                                  &
             &                           om0in,omLin,hubin,Lbox,                                 &
             &                           observer,observer_redshift,                             &
             &                           pos,vel,idtab,pottab,fparttab,ip,                       &
             &                           posout,velout,idout,zout,potout,fpartout,npout,.false., &
             &                           dist1cone, dist2cone, disttolcone)
        !===========================================================================
        if(npout>0)then                                                 ! some of these particle get selected
           do idim=1,ndim
              do i=1,npout
                 xp_out(ipout+i,idim)=posout(idim,i)/Lbox               ! collect coordinate of selected particles
                 vp_out(ipout+i,idim)=velout(idim,i)                    ! collect velocities of selected particles
#ifdef WITHPARTFORCE
                 fpart_out(ipout+i,idim)=fpartout(idim,i)               ! collect the forces on selected particles
#endif
              end do
           end do
           do i=1,npout
              id_out(ipout+i)=idout(i)                                  ! collect IDs of selected particles
           end do
           do i=1,npout
              zp_out(ipout+i)=zout(i)                                   ! collect redshifts of selected particles
#ifdef OUTPUT_PARTICLE_POTENTIAL
              pot_out(ipout+i)=potout(i)                                ! collect potentials on selected particles
#endif                    
           end do
           ipout=ipout+npout
           npart_out=npart_out+npout                                    ! increment total number of selected particles 
        endif                                                           ! if(npout>0)
     endif                                                              ! if(ip>0)
     !
     if(ipout>=nstride)then                                             ! reached size for writing into file (nstride)
        if(.not.opened) then                                            ! open file for writing if not already done so
           open(ilun,file=TRIM(fileloc),form='unformatted')
           rewind(ilun)  
           write(ilun) ncpu
           write(ilun) nstride
           write(ilun) npart
           opened=.true.
        endif
        do idim=1,ndim
           write(ilun)xp_out(1:nstride,idim)                            ! write particle coordinates
           write(ilun)vp_out(1:nstride,idim)                            ! write particle velocities
        end do        
        write(ilun)zp_out(1:nstride)                                    ! write particle redshifts
        if (write_cone_idp) write(ilun)id_out(1:nstride)                ! write particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
        write(ilun)pot_out(1:nstride)                                   ! write particle potentials
#endif 
#ifdef WITHPARTFORCE
        do idim=1,ndim
           write(ilun)fpart_out(1:nstride,idim)                         ! write particle forces
        end do
#endif
        !
        do idim=1,ndim                                                  ! move unwritten particles forward in arrays
           do i=1,ipout-nstride
              xp_out(i,idim)=xp_out(i+nstride,idim)                     ! coordinates
              vp_out(i,idim)=vp_out(i+nstride,idim)                     ! velocities
#ifdef WITHPARTFORCE 
              fpart_out(i,idim)=fpart_out(i+nstride,idim)               ! forces
#endif
           end do
        end do
        do i=1,ipout-nstride
           id_out(i)=id_out(i+nstride)                                  ! IDs
        end do
        do i=1,ipout-nstride
           zp_out(i) =zp_out(i+nstride)                                 ! redshifts
#ifdef OUTPUT_PARTICLE_POTENTIAL
           pot_out(i)=pot_out(i+nstride)                                ! potential
#endif                    

        end do
        ipout=ipout-nstride                                             ! reset quantity used to count particles
     endif                                                              ! if(ipout>=nstride)
  end do
  ! End loop over cpus

  if(ipout>0)then                                                       ! write additional particles which are not yet written
     if(.not.opened) then
        open(ilun,file=TRIM(fileloc),form='unformatted')
        rewind(ilun)  
        write(ilun) ncpu
        write(ilun) nstride
        write(ilun) npart
        opened=.true.
     endif
     do idim=1,ndim
        write(ilun)xp_out(1:ipout,idim)                                 ! coordinates
        write(ilun)vp_out(1:ipout,idim)                                 ! velocities
     end do
     write(ilun)zp_out(1:ipout)                                         ! redshifts
     if(write_cone_idp) write(ilun)id_out(1:ipout)                      ! IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
     write(ilun)pot_out(1:ipout)                                        ! potentials
#endif                    
#ifdef WITHPARTFORCE
     do idim=1,ndim
        write(ilun)fpart_out(1:ipout,idim)                              ! forces
     end do
#endif 
  endif                                                                 ! if(ipout>0)

  if(opened) close(ilun)                                                ! close file after writing complete
  if(verbose) write(*,*)'cone_fullsky output=',myid,npart_out

  if(npart_out>0) then                                                  ! write header
     fileloc=TRIM(filename)//TRIM(nchar)//'.hdr'
     open(ilun,file=TRIM(fileloc),form='unformatted')
     rewind(ilun)
     write(ilun) ncpu
     write(ilun) nstride
     write(ilun) npart_out
     close(ilun)
     if(reduce_infolocal(1).eq.0) reduce_infolocal(1) = 1               ! VR 2012
  endif
  !
  if((opened.and.(npart_out.eq.0)).or.((.not.opened).and.(npart_out>0))) then
     write(*,*)'Error in output_cone_fullsky'
     write(*,*)'npart_out=',npart_out,'opened=',opened
     stop
  endif

  ! Send the token
  if(conetoken) then
     if(mod(myid,IOGROUPSIZECONE).ne.0.and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag,MPI_COMM_WORLD,info)
     end if
  endif

!------------------ MODIF V. REVERDY 2011 ------------------!
  ! Reduce
  reduce_infolocal(2) = npart_out
  call MPI_REDUCE(reduce_infolocal,reduce_infoglobal,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
  ! Write info file
  if (myid.eq.1) then
     if(cone_overlap)then
        aendconem2_info=aendconem2
        aendconem1_info=aendconem1
        aendcone_info=aendcone
        aexpold_info=aexp_old
        !
        if((use_aexp_restart).and.(nstep_coarse_after_restart.eq.2)) aexpold_info=aexp_restart_light_cone ! RY 14/02/2017
        aexp_info=aexp
        !
        zendconem2_info=1.D0/aendconem2_info-1.D0
        zendconem1_info=1.D0/aendconem1_info-1.D0
        zendcone_info  =1.D0/aendcone_info-1.D0
        zexpold_info=1.D0/aexpold_info-1.D0
        zexp_info=1.D0/aexp-1.D0
        !
        call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
        !
        if(observer_redshift.gt.0) then
           dendconem2_info=((coord_distance(zendconem2_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dendconem1_info=((coord_distance(zendconem1_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dendcone_info  =((coord_distance(zendcone_info  ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dexpold_info   =((coord_distance(zexpold_info   ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dexp_info      =((coord_distance(zexp_info      ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
        else 
           dendconem2_info=(coord_distance(zendconem2_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dendconem1_info=(coord_distance(zendconem1_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dendcone_info  =(coord_distance(zendcone_info   ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dexpold_info   =(coord_distance(zexpold_info    ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dexp_info      =(coord_distance(zexp_info       ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
        end if
     else
        aendconem2_info=0.D0
        aendconem1_info=0.D0
        aendcone_info=0.D0
        aexpold_info=0.D0
        aexp_info=0.D0
        zendconem2_info=0.D0
        zendconem1_info=0.D0
        zendcone_info=0.D0
        zexpold_info=0.D0
        zexp_info=0.D0
        dendconem2_info=0.D0
        dendconem1_info=0.D0
        dendcone_info=0.D0
        dexpold_info=0.D0
        dexp_info=0.D0
     endif

     call write_infoconepart(npart_out,cone_id,cone_zmax,observer_x,observer_y,observer_z,observer_redshift, &
               & 1.D0/(z2+1.D0),1.D0/(z1+1.D0),z2,z1,dist2cone,dist1cone,disttolcone,                        &
               & reduce_infoglobal,                                                                          &
               & is_fullsky,input_thetay,input_thetaz,input_theta,input_phi,                                 &
               & aendconem2_info,aendconem1_info,aendcone_info,aexpold_info,aexp_info,                       &
               & zendconem2_info,zendconem1_info,zendcone_info,zexpold_info,zexp_info,                       &
               & dendconem2_info,dendconem1_info,dendcone_info,dexpold_info,dexp_info,                       &
               & future)
  end if
!-----------------------------------------------------------!

end subroutine output_cone_fullsky


!===========================================================================
subroutine perform_my_selection_simple(onlydist,justcount,z1,z2,                                  &
     &                                 om0in,omLin,hubin,Lbox,                                    &
     &                                 observer,observer_redshift,                                &
     &                                 pos,vel,idtab,pottab,fparttab,npart,                       &
     &                                 posout,velout,idout,zout,potout,fpartout,npartout,verbose, &
     &                                 dist1cone,dist2cone,disttolcone)
  !===========================================================================
  ! All the quantities below are real*8 except
  !      juscount : logical
  !      npart,npartout : integer*4
  !
  ! juscount : .true. to just count the particles to be selected. Needed
  !            to allocate appropriately the arrays posout,velout and zout
  !            The parameter npartout becomes then an output given the number
  !            of particles selected.
  !            .false. to perform the actual selection: npartout is an input
  !            then  posout, velout, zout are appropriate outputs.
  !
  ! z1,z2    : the lightcone part of interest is in between z1 and z2, 
  !            with z1 < z2. If we consider a redshift z(t) where all the 
  !            particles are synchrone, and if coarse timestep is
  !            a fixed dt, it is most appropriate to choose z1 and z2 such
  !            that z1=z(t+dt/2) and z2=z(t-dt/2) to have best accuracy.
  !
  ! om0in    : the value of the cosmological parameter omega0 (typically 0.3)
  !
  ! omLin    : the value of the cosmological constant lambda (typically 0.7)
  !
  ! hubin    : the value of H0/100, where H0 is the present Hubble constant
  !            in km/s/Mpc (typically 0.7)
  !
  ! Lbox     : the comoving size of the simulation box in Mpc (NOT in Mpc/h)
  !
  ! observer(3) : the observer position in the box in Mpc, assuming that
  !            coordinates are in [0,Lbox[
  !
  ! pos(3,npart) : comoving positions of the input particles in Mpc, assumed to be
  !            in [0,Lbox[.
  !
  ! vel(3,npart) : velocities of the input particles (in any unit, it does not
  !            matter)
  !
  ! npart    : number of input particles to be treated
  !
  ! posout(3,npartout) : output comoving positions of selected particles in Mpc.
  !
  ! velout(3,npartout) : output velocities of selected particles
  !
  ! zout(npartout) : output redshift of selected particles
  !
  ! npartout : number of selected particles. To be computed appropriately,
  !            this routine must be called with juscount=.true., which will give
  !            npartout as an output. Then this routine must be called with 
  !            juscount=.false. with the correct value of npartout, after having
  !            allocated correctly arrays posout,velout,zout.
  !===========================================================================
  use amr_parameters, ONLY: nvector,sp
  implicit none
  logical::onlydist
  logical :: justcount,verbose
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
  real(kind=8) :: observer(3)
  integer :: npart,npartout
  real(kind=8) :: pos(3,nvector),vel(3,nvector)
  integer(kind=8) :: idtab(nvector)

  real(kind=8) :: posout(3,27*nvector),velout(3,27*nvector),zout(27*nvector)
  integer(kind=8) :: idout(27*nvector)
  real(sp) :: pottab(nvector)
  real(sp) :: potout(27*nvector)
  real(sp) :: fparttab(3,nvector)
  real(sp) :: fpartout(3,27*nvector)
  real(kind=8) :: coord_distance  
  real(kind=8) :: dist1,dist2
  real(kind=8) :: dist1cone,dist2cone,disttolcone ! modif V. REVERDY 2012
  real(kind=8) :: xcoord,ycoord,zcoord
  real(kind=8) :: dist,dxtest1,dxtest2,facnorm
  integer :: nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
  integer :: i,j,k,np,npartcount
  integer :: myint
  real(kind=8) :: pi
  real(kind=8)::observer_redshift ! V. REVERDY
  if (verbose) write(*,*) 'Enter perform_my_selection_simple'
  
  ! initialize cosmological parameters
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  
  if(verbose) write(*,*) 'After init_cosmo',Omega0,OmegaL,OmegaR,coverH0
  
  ! compute comoving distance of the photon planes from the observer
  ! dist1,dist2=integral of c.dt/a between zero and z1,z2
  ! dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0) ! VR
  ! dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0) ! VR
  if(observer_redshift.gt.0) then
     dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0) ! VR
     dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0) ! VR
  else 
     dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0) ! VR
     dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0) ! VR
  end if
  !
  dist1cone=dist1/Lbox                                                  ! convert to code unit (dimensionless)
  dist2cone=dist2/Lbox                                                  ! convert to code unit (dimensionless)
  
  if(verbose) write(*,*) 'After coord',dist1,dist2

  if(.not.onlydist) then
     ! compute how many replica are needed by finding a 
     ! cube containing both the whole shell and oberver
     nrepxm=myint((observer(1)-dist2)/Lbox)                             ! number of cubes needed from oberver to left
     nrepxp=myint((observer(1)+dist2)/Lbox)                             ! number of cubes needed from oberver to right
     nrepym=myint((observer(2)-dist2)/Lbox)                             ! number of cubes needed from oberver to front
     nrepyp=myint((observer(2)+dist2)/Lbox)                             ! number of cubes needed from oberver to back
     nrepzm=myint((observer(3)-dist2)/Lbox)                             ! number of cubes needed from oberver to top
     nrepzp=myint((observer(3)+dist2)/Lbox)                             ! number of cubes needed from oberver to bottom
     !
     facnorm=1.0d0/(dist2-dist1)                                        ! inverse thickness of shell

     npartcount=0                                                       ! initialise count of selected particles
     !   
     if(verbose) write(*,*) 'before loop',facnorm
     ! loop on all the replica of potential interest
     do k=nrepzm,nrepzp,1
        do j=nrepym,nrepyp,1
           do i=nrepxm,nrepxp,1
              do np=1,npart
                 xcoord=pos(1,np)+Lbox*dble(i)-observer(1)              ! x distance between particles and observer in Mpc
                 ycoord=pos(2,np)+Lbox*dble(j)-observer(2)              ! y distance between particles and observer in Mpc
                 zcoord=pos(3,np)+Lbox*dble(k)-observer(3)              ! z distance between particles and observer in Mpc
                 ! 
                 dist=dsqrt(xcoord**2+ycoord**2+zcoord**2)              ! distance between particles and observer in Mpc
                 if(dist>dist1.and.dist<=dist2) then                    ! select partice if in shell
                    ! this particle is good, we can add it to the list
                    npartcount=npartcount+1                             ! incremennt particle count
                    !              
                    posout(1,npartcount)=xcoord                         ! particle x coordinate in Mpc
                    posout(2,npartcount)=ycoord                         ! particle y coordinate in Mpc
                    posout(3,npartcount)=zcoord                         ! particle z coordinate in Mpc
                    velout(1,npartcount)=vel(1,np)                      ! particle x velocity
                    velout(2,npartcount)=vel(2,np)                      ! particle y velocity
                    velout(3,npartcount)=vel(3,np)                      ! particle z velocity
                    idout(npartcount)=idtab(np)                         ! particle ID
#ifdef OUTPUT_PARTICLE_POTENTIAL
                    potout(npartcount)=pottab(np)                       ! particle potential
#endif        
#ifdef WITHPARTFORCE
                    fpartout(1,npartcount)=fparttab(1,np)               ! particle x force
                    fpartout(2,npartcount)=fparttab(2,np)               ! particle y force
                    fpartout(3,npartcount)=fparttab(3,np)               ! particle z force     
#endif    
                    ! compute redshift of particle by linear interpolation
                    dxtest1=dist-dist1
                    dxtest2=dist2-dist
                    zout(npartcount)=(dxtest1*z2+dxtest2*z1)*facnorm    ! particle redshift
                    
                 endif
              enddo
           enddo                                                        ! do np=1,npart
        enddo
     enddo
     npartout=npartcount
     !
  endif
  !
  if(verbose) write(*,*) 'End of perform_my_selection_simple',npartout

end subroutine perform_my_selection_simple

!=======================================================================
function myint(x)
  !=======================================================================
  ! The REAL int function
  !=======================================================================
  real(kind=8) :: x
  integer :: myint

  if(x.ge.0.0D0) then
     myint=int(x)
  else
     myint=int(x)-1
  endif
end function myint

!===============================================!
!MAP OF PROJECTED DM DENSITY USING CIC SMOOTHING!
!Adapted from part2map.f90 by Yann Rasera       !
!Date: 11/09/08      
!Remarks: Not tested for serial                 !
!
!===============================================!
subroutine part2map_ramses(proj,nx,ny,xmin,xmax,ymin,ymax,zmin,zmax,periodic,filename)
  use amr_commons,ONLY:myid,next,headl,numbl,levelmin,ncpu,dp,sp,verbose
  use pm_commons,ONLY:xp,mp,nextp,headp,numbp
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  character(LEN=1)::proj
  integer::nx,ny
  real(dp)::xmin,xmax,ymin,ymax,zmin,zmax
  logical::periodic
  character(LEN=200)::filename
  integer::ilun
  
  real(dp),dimension(:,:),allocatable::map
  real(dp),dimension(:,:),allocatable::map_sum
  real(dp)::xxmin,xxmax,yymin,yymax,dx,dy,ddx,ddy,dex,dey
  logical::ok_part
  integer::npart_tot,npart1
  real(dp)::mtot
  integer::idim,jdim,icpu,igrid,jgrid,ilevel,ipart,jpart
  integer::ix,iy,ixp1,iyp1
  integer::nx_final,ny_final
  integer::info

  
#ifndef WITHOUTMPI
  if(myid==1) then
#endif  
     if (verbose) then
        write(*,*)'Entering part2map'
        write(*,*)'Working map =',nx,ny
        write(*,*)'projection  =',proj
        write(*,*)'xrange      =',xmin,xmax
        write(*,*)'yrange      =',ymin,ymax
        write(*,*)'zrange      =',zmin,zmax
        write(*,*)'periodic    =',periodic
        write(*,*)'filename    =',trim(filename)
     endif
#ifndef WITHOUTMPI   
  endif
#endif   
    
!Allocate map and set projection
  allocate(map(0:nx,0:ny))
  map=0.
  
  if (proj=='x')then
     idim=2
     jdim=3
     xxmin=ymin ; xxmax=ymax
     yymin=zmin ; yymax=zmax
  else if (proj=='y') then
     idim=1
     jdim=3
     xxmin=xmin ; xxmax=xmax
     yymin=zmin ; yymax=zmax
  else
     idim=1
     jdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
  end if
  dx=(xxmax-xxmin)/nx
  dy=(yymax-yymin)/ny
  
  
!Select particles and compute map
  npart_tot=0
  mtot=0.
  ilevel=levelmin
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)    
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           npart_tot=npart_tot+npart1
           ipart=headp(igrid)
           do jpart=1,npart1
              
              !Projection particle ipart
              if(periodic)then
                 ok_part=(xp(ipart,1)>=xmin.and.xp(ipart,1)<=xmax.and. &
                      &   xp(ipart,2)>=ymin.and.xp(ipart,2)<=ymax.and. &
                      &   xp(ipart,3)>=zmin.and.xp(ipart,3)<=zmax)
                 if(ok_part)then
                    ddx=(xp(ipart,idim)-xxmin)/dx
                    ddy=(xp(ipart,jdim)-yymin)/dy
                    ix=ddx
                    iy=ddy
                    ddx=ddx-ix
                    ddy=ddy-iy
                    dex=1.0-ddx
                    dey=1.0-ddy
                    if(ix<0)ix=ix+nx
                    if(ix>=nx)ix=ix-nx
                    if(iy<0)iy=iy+ny
                    if(iy>=ny)iy=iy-ny
                    ixp1=ix+1
                    iyp1=iy+1
                    if(ixp1<0)ixp1=ixp1+nx
                    if(ixp1>=nx)ixp1=ixp1-nx
                    if(iyp1<0)iyp1=iyp1+ny
                    if(iyp1>=ny)iyp1=iyp1-ny
                    map(ix  ,iy  )=map(ix  ,iy  )+mp(ipart)*dex*dey
                    map(ix  ,iyp1)=map(ix  ,iyp1)+mp(ipart)*dex*ddy
                    map(ixp1,iy  )=map(ixp1,iy  )+mp(ipart)*ddx*dey
                    map(ixp1,iyp1)=map(ixp1,iyp1)+mp(ipart)*ddx*ddy
                    mtot=mtot+mp(ipart)
                 end if
              else
                 ok_part=(xp(ipart,1)>=xmin.and.xp(ipart,1)<=xmax.and. &
                      &   xp(ipart,2)>=ymin.and.xp(ipart,2)<=ymax.and. &
                      &   xp(ipart,3)>=zmin.and.xp(ipart,3)<=zmax)
                 if(ok_part)then
                    ddx=(xp(ipart,idim)-xxmin)/dx
                    ddy=(xp(ipart,jdim)-yymin)/dy
                    ix=ddx
                    iy=ddy
                    ddx=ddx-ix
                    ddy=ddy-iy
                    dex=1.0-ddx
                    dey=1.0-ddy
                    if(ix>=0.and.ix<nx.and.iy>=0.and.iy<ny)then
                       map(ix  ,iy  )=map(ix  ,iy  )+mp(ipart)*dex*dey
                       map(ix+1,iy  )=map(ix+1,iy  )+mp(ipart)*ddx*dey
                       map(ix  ,iy+1)=map(ix  ,iy+1)+mp(ipart)*dex*ddy
                       map(ix+1,iy+1)=map(ix+1,iy+1)+mp(ipart)*ddx*ddy
                       mtot=mtot+mp(ipart)
                    endif
                 end if
              endif
              !End projection particle ipart
              
              ipart=nextp(ipart) ! Go to next particle
           end do
           ! End loop over particles    
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
  end do
  ! End loop over cpus
  
#ifndef WITHOUTMPI
  if(myid==1) then
#endif  
     if (verbose)then
        write(*,*)'Proc',myid
        write(*,*)'Found',npart_tot,' particles' 
        write(*,*)'Mass=',mtot
     endif
#ifndef WITHOUTMPI   
  endif
#endif   

  !Gather and write file
  if(periodic)then
     nx_final=nx
     ny_final=ny
  else
     nx_final=nx+1
     ny_final=ny+1
  endif
  
  allocate(map_sum(nx_final,ny_final))
  
#ifndef WITHOUTMPI 
  if (ncpu>1)then
     call MPI_REDUCE(map(0:nx_final-1,0:ny_final-1),map_sum,nx_final*ny_final,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
  else
     map_sum=map(0:nx_final-1,0:ny_final-1)
  endif
#else
  map_sum=map(0:nx_final-1,0:ny_final-1)  
#endif
  
  if (verbose)then
     write(*,*)'Total mass=',sum(map_sum)
  endif
  
  ilun=3*ncpu+myid+10
#ifndef WITHOUTMPI
  if(myid==1) then
#endif
     if(verbose) write(*,*)'Write file '//TRIM(filename)
     open(unit=ilun,file=trim(filename),form='unformatted')
     write(ilun)nx_final,ny_final
     write(ilun)real(map_sum,kind=sp)
     close(ilun)    
#ifndef WITHOUTMPI   
  endif
#endif   
  deallocate(map,map_sum)

#ifndef WITHOUTMPI
  if(myid==1) then 
#endif 
  if(verbose)write(*,*)'Part2map completed'
#ifndef WITHOUTMPI   
  endif
#endif   
 
end subroutine part2map_ramses

!---------------------------------------------------------!
!---------------------------------------------------------!
! This subroutine samples particles in a specified  cube by
! user by providing minimum and maximum x, y, z coordinates
!
! It writes a particle into file as soon as it is selected,
! rather than waiting for reaching some certain batch size.
!---------------------------------------------------------!
subroutine extract_sample(xmin,xmax,ymin,ymax,zmin,zmax,nsample,filename)
  use amr_commons,ONLY : myid,next,headl,numbl,levelmin,ncpu,dp,sp,verbose, &
                       & nstride,ndim,aexp,i8b,i4b,iogroupsize,sampletoken,
                       & IOGROUPSIZESAMPLE,write_sample_idp,adaptive_iogroupsize
  use pm_commons, ONLY : xp,vp,mp,idp,nextp,headp,numbp,npart
  use mpi

  implicit none
  real(dp)::xmin,xmax,ymin,ymax,zmin,zmax
  integer::nsample !Warning: type should be compatible with idp to compute mod(idp(ipart),nsample)!!!
  character(LEN=200)::filename

  integer::dummy_io,info
  integer,parameter::tag=1900
  character(LEN=5)::nchar
  character(LEN=200)::fileloc
  real(sp),dimension(1:nstride,1:ndim)::xp_out,vp_out
  integer,dimension(1:nstride)::ind_part
  integer::ilun
  logical::ok_part
  integer::npart1,npart_out
  integer::icpu,igrid,jgrid,ilevel,ipart,jpart,idim,jdim,ip,i
  integer::ierr
  logical::opened
 !------------------------- RV 2012 ------------------------! 
  integer(kind=8),dimension(1:nstride)::idp_out
  integer(kind=8),dimension(1:2)::reduce_infolocal
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  real(kind=8)::iovolume
  !
  opened=.false.
  reduce_infolocal = 0
  reduce_infoglobal = 0
!-------------------------- RV 2012 ------------------------! 

#ifndef WITHOUTMPI
  if(myid.eq.1) then
#endif 
  if(verbose) then
     write(*,*) 'Entering extract_sample with parameters'
     write(*,*) 'nsample     =',nsample
     write(*,*) 'xrange      =',xmin,xmax
     write(*,*) 'yrange      =',ymin,ymax
     write(*,*) 'zrange      =',zmin,zmax
  endif
#ifndef WITHOUTMPI   
  endif
#endif   

  iovolume=(dabs(xmax-xmin)*dabs(ymax-ymin)*dabs(zmax-zmin))
  if(adaptive_iogroupsize) IOGROUPSIZESAMPLE = nint(iovolume*real(IOGROUPSIZE, KIND=8))
  if(IOGROUPSIZESAMPLE.le.1) IOGROUPSIZESAMPLE=1
  if(IOGROUPSIZESAMPLE.ge.IOGROUPSIZE) IOGROUPSIZESAMPLE=IOGROUPSIZE
  if(verbose) write(*,*) 'IOGROUPSIZESAMPLE = ',IOGROUPSIZESAMPLE,'iovolume = ',iovolume

  ! Wait for the token
  if(sampletoken) then
    if(mod(myid-1,IOGROUPSIZESAMPLE).ne.0) then
      call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
    end if
  endif
   
  !!!!NOTE: Do I have to declare ipart as a integer(8) for big run??? 

  ! select particles and write them in file
  ilun=3*ncpu+myid+10                                                   ! file handler
  call title(myid,nchar)                                                ! create file name
  fileloc=TRIM(filename)//TRIM(nchar)//'.dat'                           ! create file name
  
  npart_out=0                                                           ! total number of selected particles
  ilevel=levelmin
  ip=0
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)                                           ! first grid on CPU  
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)                                             ! number of particles in the grid
        if(npart1>0)then
           ipart=headp(igrid)                                           ! first particle in the grid
           do jpart=1,npart1
              !Selection particle ipart
              ok_part = &
                  &(xp(ipart,1).ge.xmin.and.xp(ipart,1).le.xmax.and. &  ! particle in desired x range
                  & xp(ipart,2).ge.ymin.and.xp(ipart,2).le.ymax.and. &  ! particle in desired y range
                  & xp(ipart,3).ge.zmin.and.xp(ipart,3).le.zmax.and. &  ! particle in desired z range
                  &  mod(idp(ipart),int(nsample,kind=i8b))==0)          ! downsample particles      
                                                                        ! NOTE: a better sampling method would be with a random generator
              if(ok_part) then
                 ip=ip+1                                                ! increment number of selected particles 
                 ind_part(ip)=ipart                                     ! gather global particle ID and store in ind_part
                 if(ip.eq.nstride)then
                    do idim=1,ndim
                       do i=1,ip
                          xp_out(i,idim)=real(xp(ind_part(i),idim),kind=sp) ! particle coordinates
                          vp_out(i,idim)=real(vp(ind_part(i),idim),kind=sp) ! particle velocities
                       end do
                    end do
                    do i=1,ip
                       idp_out(i)=idp(ind_part(i))                      ! particle ID
                    end do
                    if(.not.opened) then                                ! open file for writing if not already
                       open(ilun,file=TRIM(fileloc),form='unformatted')
                       rewind(ilun)  
                       write(ilun) ncpu
                       write(ilun) nstride
                       write(ilun) npart
                       opened=.true.
                    endif
                    do idim=1,ndim
                       write(ilun) xp_out(1:nstride,idim)               ! write particle coordinates
                       write(ilun) vp_out(1:nstride,idim)               ! write particle velocities
                    end do
                    if (write_sample_idp) write(ilun)idp_out(1:nstride) ! write particle IDs
                    ip=0
                 endif
                 !
                 npart_out=npart_out+1                                  ! increment total number of selected particles
              endif                                                     ! if(ok_part)
              !
              ipart=nextp(ipart)                                        ! go to next particle
           end do                                                       ! end loop over particles    
           !
        end if
        !
        igrid=next(igrid)                                               ! go to next grid
     end do                                                             ! end loop over grids
     !
  end do                                                                ! end loop over cpus
  !
  if(ip>0) then                                                         ! gather and write additional particles
     do idim=1,ndim
        do i=1,ip
           xp_out(i,idim)=real(xp(ind_part(i),idim),kind=sp)            ! particle coordinates
           vp_out(i,idim)=real(vp(ind_part(i),idim),kind=sp)            ! particle velocities
        end do
     end do
     do i=1,ip
        idp_out(i)=idp(ind_part(i))                                     ! particle IDs
     end do
     if(.not.opened) then                                               ! open file for writing if not yet
        open(ilun,file=TRIM(fileloc),form='unformatted')
        rewind(ilun)  
        write(ilun) ncpu
        write(ilun) nstride
        write(ilun) npart
        opened=.true.
     endif
     do idim=1,ndim
        write(ilun) xp_out(1:ip,idim)                                   ! write particle coorindates
        write(ilun) vp_out(1:ip,idim)                                   ! write particle velocities
     end do 
     if(write_sample_idp) write(ilun)idp_out(1:ip)                      ! write particle IDs
  endif                                                                 ! if(ip>0)
  
  if(opened) close(ilun)                                                ! close file after writing
  !
  if(npart_out>0) then
     fileloc=TRIM(filename)//TRIM(nchar)//'.hdr'
     open(ilun,file=TRIM(fileloc),form='unformatted')
     rewind(ilun)
     write(ilun) ncpu
     write(ilun) nstride
     write(ilun) npart_out
     write(ilun) aexp
     close(ilun)
     if(reduce_infolocal(1).eq.0) reduce_infolocal(1)=1                 ! VR 2012
  endif
  
#ifndef WITHOUTMPI
  if(myid.eq.1) then
#endif 
     if(verbose) write(*,*)'Extract_sample completed'
#ifndef WITHOUTMPI   
  endif
#endif   

  if((opened.and.(npart_out.eq.0)).or.((.not.opened).and.(npart_out>0))) then
     write(*,*) 'Error in extract_sample'
     write(*,*) 'npart_out=',npart_out,'opened=',opened
     stop
  endif

  ! send the token
  if(sampletoken) then
     if(mod(myid,IOGROUPSIZESAMPLE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag,MPI_COMM_WORLD,info)
     end if
  endif

!------------------------- RV 2012 -------------------------! 
  ! write info file
  reduce_infolocal(2)=npart_out
  call MPI_REDUCE(reduce_infolocal,reduce_infoglobal,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! write info file
  if (myid.eq.1) then
    call write_infosamplepart(npart_out,aexp,                &
                            & xmin,xmax,ymin,ymax,zmin,zmax, &
                            & reduce_infoglobal)
  end if
!------------------------- RV 2012 -------------------------! 
  
end subroutine extract_sample

!---------------------------------------------------------!
!---------------------------------------------------------!
!
!---------------------------------------------------------!
subroutine output_cone(filename,cone_id,future) ! NOT DEFAULT RAMSES' LC!
  use amr_commons
  use pm_commons
  use mpi
  
  implicit none
  character(LEN=200)::filename
  integer::dummy_io,info
  integer,parameter::tag=2000
  
  integer::ilun,nx_loc,ipout,npout,npart_out
  character(LEN=200)::fileloc
  character(LEN=5)::nchar
  real(kind=8),dimension(1:ndim,1:nvector)::pos,vel
  integer(kind=8),dimension(1:nvector)::idtab
  real(kind=8),dimension(1:ndim,1:27*nvector)::posout,velout
  integer(kind=8),dimension(1:27*nvector)::idout
  real(kind=8),dimension(1:27*nvector)::zout
  real(sp),dimension(1:nstride+27*nvector,1:ndim)::xp_out,vp_out
  real(sp),dimension(1:nstride+27*nvector)::zp_out
  integer(kind=8),dimension(1:nstride+27*nvector)::id_out

  real(sp),dimension(1:nvector)::pottab
  real(sp),dimension(1:27*nvector)::potout
  real(sp),dimension(1:nstride+27*nvector)::pot_out

  real(sp),dimension(1:ndim,1:nvector)::fparttab
  real(sp),dimension(1:ndim,1:27*nvector)::fpartout
  real(sp),dimension(1:ndim,1:nstride+27*nvector)::fpart_out


  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: observer(3),thetay,thetaz,theta,phi
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel
  integer::i,ig,ip,npart1
  integer(kind=8),dimension(1:2)::reduce_infolocal
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  integer,dimension(1:nvector)::ind_part
  logical::opened
  integer:: future
  real(kind=8)::dist1cone, dist2cone,disttolcone
  integer::ierr
  ! aexp for info file
  real(kind=8)::aendconem2_info,aendconem1_info,aendcone_info,aexpold_info,aexp_info
  !aendconem* -> expansion factor for end of the shell
  !the end of previous shell and the end of previous previous shell (which is the begining of the current shell)
  !when activating overlapping buffer option.
  real(kind=8)::zendconem2_info,zendconem1_info,zendcone_info,zexpold_info,zexp_info
  real(kind=8)::dendconem2_info,dendconem1_info,dendcone_info,dexpold_info,dexp_info

  real(kind=8)::Omega0,OmegaL,OmegaR,coverH0
  real(kind=8)::coord_distance
  real(kind=8)::observer_redshift

  integer::cone_id
  real(kind=8)::iovolume
  
  iovolume = 1.D0

  opened=.false.
  reduce_infolocal=0
  reduce_infoglobal=0
  dist1cone=0.D0
  dist2cone=0.D0
  ierr=0
  
  if(nstep_coarse.lt.4)return
  if(verbose)write(*,*)'Entering output_part'
  
  z2=1.0D0/aexp_old-1.0D0
  z1=1.0D0/aexp    -1.0D0
  if(use_aexp_restart.and.(nstep_coarse_after_restart.eq.2)) z2=1.D0/aexp_restart_light_cone-1.D0 ! RY put it at this position 14/02/2017 (previously it was after)

  if(cone_overlap) then
     if((aendconem1.lt.aendconem2).or.(aendcone.lt.aendconem1)) print*,'WARNING CONE SHELL RANGE NOT WELL ORDERED, AEXP ',aendconem2,aendconem1,aendcone 
     if(future.eq.1) then
        z2=1.0D0/aendconem1-1.0D0
        z1=1.0D0/aendcone  -1.0D0
     else if(future.eq.-1) then
        z2=1.0D0/aendconem2-1.0D0
        z1=1.0D0/aendconem1-1.0D0
     else
        print*,'not implemented cone part future=',future
        stop
     endif
     if(myid.eq.1) then
        print*,''
        print*,'TEST CONE PART LIMIT'
        print*,'nstepcoarse,z2,z1',nstep_coarse,z2,z1
        print*,'aendconem2,aendconem1,aendcone',aendconem2,aendconem1,aendcone
        print*,''
     endif
  endif

  !zmax_cone=7.
  if(z2.gt.conenarrow_zmax(cone_id)) return
  if(dabs(z2-z1)<1d-6) return

  theta =conenarrow_theta (cone_id)                                     ! angular coordinate of cone centre
  phi   =conenarrow_phi   (cone_id)                                     ! angular coordinate of cone cnetre
  thetay=conenarrow_thetay(cone_id)                                     ! opening angle of cone
  thetaz=conenarrow_thetaz(cone_id)                                     ! opening angle

  om0in=omega_m                                                         ! Omega_m
  omLin=omega_l                                                         ! Omega_Lambda
  hubin=h0/100.0D0                                                      ! h
  Lbox=boxlen_ini/hubin                                                 ! box size in Mpc
  !
  observer(1) = Lbox*conenarrow_observer_x(cone_id)                     ! observer x coordinate (VR)
  observer(2) = Lbox*conenarrow_observer_y(cone_id)                     ! observer y coordinate (VR)
  observer(3) = Lbox*conenarrow_observer_z(cone_id)                     ! observer z coordinate (VR)

  ! compute iogroupsize
  call perform_my_selection_simple(.true.,.false.,z1,z2,                                   &
       &                           om0in,omLin,hubin,Lbox,                                 &
       &                           observer,conenarrow_observer_redshift(cone_id),         &
       &                           pos,vel,idtab,pottab,fparttab,ip,                       &
       &                           posout,velout,idout,zout,potout,fpartout,npout,.false., &
       &                           dist1cone,dist2cone,disttolcone)
  iovolume = ((2.0D0*dist2cone)**3)-(((2.D0*dist1cone)/dsqrt(real(ndim,kind=8)))**3)
  iovolume = iovolume*(conenarrow_thetay(cone_id)*conenarrow_thetaz(cone_id))/(41263.D0)
  if(adaptive_iogroupsize) IOGROUPSIZECONE=nint(iovolume*real(IOGROUPSIZE,kind=8))
  if(IOGROUPSIZECONE.le.1) IOGROUPSIZECONE=1
  if(IOGROUPSIZECONE.ge.IOGROUPSIZE) IOGROUPSIZECONE=IOGROUPSIZE
  if(verbose) write(*,*) cone_id, 'IOGROUPSIZECONE = ',IOGROUPSIZECONE, 'iovolume = ', iovolume
  
  ! wait for the token
  if(mod(myid-1,IOGROUPSIZE).ne.0) then
     call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
  end if

  ilun=3*ncpu+myid+10                                                   ! file handler
  call title(myid,nchar)                                                ! create file name
  fileloc=TRIM(filename)//TRIM(nchar)//'.dat'                           ! create file name

  npart_out=0
  ipout=0
  npout=0
  ilevel=levelmin
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)                                           ! first grid on CPU
     ip=0                                                               ! initialise particle counter
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)                                             ! number of particles in the grid
        if(npart1>0) then        
           ipart=headp(igrid)                                           ! first particle in the grid
           ! Loop over particles
           do jpart=1,npart1
              ip=ip+1                                                   ! increment particle counter
              ind_part(ip)=ipart                                        ! gather particle global ID and store in ind_part
              !
              if(ip.eq.nvector) then                                    ! number of particles reaches nvector
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do i=1,ip
                       pos(idim,i)=xp(ind_part(i),idim)*Lbox            ! gather particle coordinates in Mpc
                       vel(idim,i)=vp(ind_part(i),idim)                 ! gather particle velocities
#ifdef WITHPARTFORCE
                       fparttab(idim,i)=fpart(ind_part(i),idim)         ! gather particle forces
#endif
                    end do
                 end do
                 do i=1,ip
                    idtab(i)=idp(ind_part(i))                           ! gather particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
                    pottab(i)=ptcl_phi(ind_part(i))                     ! gather particle potentials
#endif 
                 end do
                 !===========================================================================
                 call perform_my_selection(.false.,z1,z2,                                                                 &
                      &                           om0in,omLin,hubin,Lbox,                                                 &
                      &                           observer,conenarrow_observer_redshift(cone_id),thetay,thetaz,theta,phi, &
                      &                           pos,vel,idtab,pottab,fparttab,ip,                                       &
                      &                           posout,velout,idout,zout,potout,fpartout,npout,.false.,                 &
                      &                           dist1cone,dist2cone)
                 !===========================================================================
                 if(npout>0) then
                    do idim=1,ndim
                       do i=1,npout
                          xp_out(ipout+i,idim)=posout(idim,i)/Lbox      ! collect particle coordinates for writing (dimensionless)
                          vp_out(ipout+i,idim)=velout(idim,i)           ! collect particle velocity for writing
#ifdef WITHPARTFORCE
                          fpart_out(ipout+i,idim)=fpartout(idim,i)      ! collect particle forces for writing
#endif
                       end do
                    end do
                    do i=1,npout
                       id_out(ipout+i)=idout(i)                         ! collect particle IDs for writing
                    end do
                    do i=1,npout
                       zp_out(ipout+i)=zout(i)                          ! collect particle redshifts for writing
#ifdef OUTPUT_PARTICLE_POTENTIAL
                       pot_out(ipout+i)=potout(i)                       ! collect particle potentials for writing
#endif     
                    end do
                    ipout=ipout+npout
                    npart_out=npart_out+npout                           ! increment total number of selected particles
                 endif
                 ip=0                                                   ! reset particle counter
              end if                                                    ! if(ip.eq.nvector)
              !
              if(ipout.ge.nstride) then                                 ! write into file when array size exceeds nstride
                 if(.not.opened) then
                    open(ilun,file=TRIM(fileloc),form='unformatted')    ! open file for write if not yet
                    rewind(ilun)  
                    write(ilun) ncpu
                    write(ilun) nstride
                    write(ilun) npart
                    opened=.true.
                 endif
                 do idim=1,ndim
                    write(ilun) xp_out(1:nstride,idim)                  ! write particle coordinates
                    write(ilun) vp_out(1:nstride,idim)                  ! write particle velocities
                 end do
                 write(ilun) zp_out(1:nstride)                          ! write particle redshifts
                 if(write_cone_idp) write(ilun) id_out(1:nstride)       ! write particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
                 write(ilun) pot_out(1:nstride)                         ! write particle potentials
#endif 
#ifdef WITHPARTFORCE
                 do idim=1,ndim
                    write(ilun)fpart_out(1:nstride,idim)                ! write particle forces
                 end do
#endif
                 ! for particles beyond nstride, move to the start of array and
                 ! write them later when particle number exceeds nstride again.
                 do idim=1,ndim
                    do i=1,ipout-nstride
                       xp_out(i,idim)=xp_out(i+nstride,idim)            ! coordinates
                       vp_out(i,idim)=vp_out(i+nstride,idim)            ! velocities
#ifdef WITHPARTFORCE
                       fpart_out(i,idim)=fpart_out(i+nstride,idim)      ! forcea
#endif
                    end do
                 end do
                 do i=1,ipout-nstride
                    id_out(i)=id_out(i+nstride)                         ! IDs
                 end do
                 do i=1,ipout-nstride
                    zp_out(i)=zp_out(i+nstride)                         ! redshifts
#ifdef OUTPUT_PARTICLE_POTENTIAL
                    pot_out(i)=pot_out(i+nstride)                       ! potentials
#endif   
                 end do
                 ipout=ipout-nstride                                    ! reset pointer to first particle
              endif                                                     ! if(ipout.ge.nstride)
              !
              ipart=nextp(ipart)                                        ! go to next particle
           end do
           ! End loop over particles           
        end if
        igrid=next(igrid)                                               ! go to next grid
     end do
     ! End loop over grids

     if(ip>0) then                                                      ! deal with remaining particles
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ip
              pos(idim,i)=xp(ind_part(i),idim)*Lbox                     ! gather particle coordinates in Mpc
              vel(idim,i)=vp(ind_part(i),idim)                          ! gather particle velocities
#ifdef WITHPARTFORCE
              fparttab(idim,i)=fpart(ind_part(i),idim)                  ! gather particle forces
#endif
           end do
        end do
        do i=1,ip
           idtab(i)=idp(ind_part(i))                                    ! gather particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
           pottab(i)=ptcl_phi(ind_part(i))                              ! gather particle potentials
#endif  
        end do
        !===========================================================================
        call perform_my_selection(.false.,z1,z2,                                                                 &
             &                           om0in,omLin,hubin,Lbox,                                                 &
             &                           observer,conenarrow_observer_redshift(cone_id),thetay,thetaz,theta,phi, &
             &                           pos,vel,idtab,pottab,fparttab,ip,                                       &
             &                           posout,velout,idout,zout,potout,fpartout,npout,.false.,                 &
             &                           dist1cone,dist2cone)
        !===========================================================================
        !
        if(npout>0) then
           do idim=1,ndim
              do i=1,npout
                 xp_out(ipout+i,idim)=posout(idim,i)/Lbox               ! collect particle coordinates for writing (dimensionless)
                 vp_out(ipout+i,idim)=velout(idim,i)                    ! collect particle velocities for writing
#ifdef WITHPARTFORCE
                 fpart_out(ipout+i,idim)=fpartout(idim,i)               ! collect particle forces for writing
#endif
              end do
           end do
           do i=1,npout
              id_out(ipout+i)=idout(i)                                  ! collect particle IDs for writing
           end do
           do i=1,npout
              zp_out(ipout+i)=zout(i)                                   ! collect particle redshifts for writing
#ifdef OUTPUT_PARTICLE_POTENTIAL
              pot_out(ipout+i)=potout(i)                                ! collect particle potential for writing
#endif  
           end do
           ipout=ipout+npout                                            ! move pointer to particle
           npart_out=npart_out+npout                                    ! increment total number of selected particles 
        endif                                                           ! if(npout>0)
        !
     endif                                                              ! if(ip>0)
     !
     if(ipout.ge.nstride) then                                          ! write into file when array size exceeds nstride
        if(.not.opened) then                                            ! open file for writing if not yet
           open(ilun,file=TRIM(fileloc),form='unformatted')
           rewind(ilun)  
           write(ilun) ncpu
           write(ilun) nstride
           write(ilun) npart
           opened=.true.
        endif
        do idim=1,ndim
           write(ilun) xp_out(1:nstride,idim)                           ! write particle coordinates
           write(ilun) vp_out(1:nstride,idim)                           ! write partilce velocities
        end do
        write(ilun) zp_out(1:nstride)                                   ! write particle redshifts
        if(write_cone_idp) write(ilun)id_out(1:nstride)                 ! write particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
        write(ilun) pot_out(1:nstride)                                  ! write particle potentials
#endif 
#ifdef WITHPARTFORCE
        do idim=1,ndim
           write(ilun) fpart_out(1:nstride,idim)                        ! write particle forces
        end do
#endif
        ! for particles beyond nstride, move to the start of array and
        ! write them later when particle number exceeds nstride again.
        do idim=1,ndim
           do i=1,ipout-nstride
              xp_out(i,idim)=xp_out(i+nstride,idim)                     ! coordinates
              vp_out(i,idim)=vp_out(i+nstride,idim)                     ! velocities
#ifdef WITHPARTFORCE
              fpart_out(i,idim)=fpart_out(i+nstride,idim)               ! forces
#endif
           end do
        end do
        do i=1,ipout-nstride
           id_out(i)=id_out(i+nstride)                                  ! IDs
        end do
        do i=1,ipout-nstride
           zp_out(i)=zp_out(i+nstride)                                  ! redshifts
#ifdef OUTPUT_PARTICLE_POTENTIAL 
           pot_out(i)=pot_out(i+nstride)                                ! potential
#endif  
        end do
        ipout=ipout-nstride                                             ! reset particle pointer
     endif                                                              ! if(ipout.ge.nstride)
  end do
  ! End loop over cpus

  if(ipout>0) then                                                      ! write the remaining particles
     if(.not.opened) then                                               ! open file for writing if not yet
        open(ilun,file=TRIM(fileloc),form='unformatted')
        rewind(ilun)  
        write(ilun) ncpu
        write(ilun) nstride
        write(ilun) npart
        opened=.true.
     endif
     do idim=1,ndim
        write(ilun) xp_out(1:ipout,idim)                                ! write particle coordinates
        write(ilun) vp_out(1:ipout,idim)                                ! write particle velocities
     end do
     write(ilun) zp_out(1:ipout)                                        ! write particle redshifts
     if (write_cone_idp) write(ilun) id_out(1:ipout)                    ! write particle IDs
#ifdef OUTPUT_PARTICLE_POTENTIAL
     write(ilun) pot_out(1:ipout)                                       ! write particle potentials
#endif  
#ifdef WITHPARTFORCE
     do idim=1,ndim
        write(ilun) fpart_out(1:ipout,idim)                             ! write particle forces
     end do
#endif 
  endif

  if(opened) close(ilun)                                                ! close file after writing
  
  if (verbose)write(*,*)'cone output=',myid,npart_out

  if(npart_out>0) then
     fileloc=TRIM(filename)//TRIM(nchar)//'.hdr'
     open(ilun,file=TRIM(fileloc),form='unformatted')                   ! open header file for writing
     rewind(ilun)
     write(ilun) ncpu
     write(ilun) nstride
     write(ilun) npart_out
     close(ilun)
     if(reduce_infolocal(1).eq.0) reduce_infolocal(1) = 1 ! modif V. REVERDY 2012 
  endif
  if((opened.and.(npart_out.eq.0)).or.((.not.opened).and.(npart_out>0))) then
     write(*,*) 'Error in output_cone'
     write(*,*) 'npart_out=',npart_out,'opened=',opened
     stop
  endif

  ! Send the token
  if(mod(myid,IOGROUPSIZE).ne.0.and.(myid.lt.ncpu))then
     dummy_io=1
     call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag,MPI_COMM_WORLD,info)
  end if

!-------------------------- VR 2012 ------------------------!           ! write info file 
  reduce_infolocal(2) = npart_out
  call MPI_REDUCE(reduce_infolocal,reduce_infoglobal,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! Write info file                               
  if (myid.eq.1) then
     if(cone_overlap)then
        aendconem2_info=aendconem2
        aendconem1_info=aendconem1
        aendcone_info=aendcone
        aexpold_info=aexp_old
        if((use_aexp_restart).and.(nstep_coarse_after_restart.eq.2))aexpold_info=aexp_restart_light_cone ! RY 14/02/2017
        aexp_info=aexp

        zendconem2_info=1./aendconem2_info-1.
        zendconem1_info=1./aendconem1_info-1.
        zendcone_info  =1./aendcone_info-1.
        zexpold_info=1./aexpold_info-1.
        zexp_info=1./aexp-1.

        call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
        if (conenarrow_observer_redshift(cone_id).gt.0) then
           dendconem2_info=((coord_distance(zendconem2_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dendconem1_info=((coord_distance(zendconem1_info,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dendcone_info  =((coord_distance(zendcone_info  ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dexpold_info   =((coord_distance(zexpold_info   ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
           dexp_info      =((coord_distance(zexp_info      ,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0))/Lbox) ! in [0, 1]
        else 
           dendconem2_info=(coord_distance(zendconem2_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dendconem1_info=(coord_distance(zendconem1_info ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dendcone_info  =(coord_distance(zendcone_info   ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dexpold_info   =(coord_distance(zexpold_info    ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
           dexp_info      =(coord_distance(zexp_info       ,Omega0,OmegaL,OmegaR,coverH0)/Lbox) ! in [0, 1]
        end if
     else
        aendconem2_info=0.D0
        aendconem1_info=0.D0
        aendcone_info  =0.D0
        aexpold_info   =0.D0
        aexp_info      =0.D0
        zendconem2_info=0.D0
        zendconem1_info=0.D0
        zendcone_info  =0.D0
        zexpold_info   =0.D0
        zexp_info      =0.D0
        dendconem2_info=0.D0
        dendconem1_info=0.D0
        dendcone_info  =0.D0
        dexpold_info   =0.D0
        dexp_info      =0.D0
     endif

     call write_infoconepart(npart_out,conenarrow_id(cone_id),conenarrow_zmax(cone_id),conenarrow_observer_x(cone_id),  &
                           & conenarrow_observer_y(cone_id),conenarrow_observer_z(cone_id),                             &
                           & conenarrow_observer_redshift(cone_id),                                                     &
                           & 1.D0/(z2+1.D0),1.D0/(z1+1.D0),z2,z1,dist2cone,dist1cone,0.D0,                              &
                           & reduce_infoglobal,                                                                         &
                           & .false.,thetay,thetaz,theta,phi,                                                           &
                           & aendconem2_info,aendconem1_info,aendcone_info,aexpold_info,aexp_info,                      &
                           & zendconem2_info,zendconem1_info,zendcone_info,zexpold_info,zexp_info,                      &
                           & dendconem2_info,dendconem1_info,dendcone_info,dexpold_info,dexp_info,future)
  end if
!-------------------------- VR 2012 ------------------------! 

end subroutine output_cone


!---------------------------------------------------------!
!---------------------------------------------------------!
!
!---------------------------------------------------------!
subroutine perform_my_selection(justcount,z1,z2,                                           &
     &                          om0in,omLin,hubin,Lbox,                                    &
     &                          observer,observer_redshift,thetay,thetaz,theta,phi,        &
     &                          pos,vel,idtab,pottab,fparttab,npart,                       &
     &                          posout,velout,idout,zout,potout,fpartout,npartout,verbose, &
     &                          dist1cone,dist2cone)
  !===========================================================================
  ! All the quantities below are real*8 except
  !      juscount : logical
  !      npart,npartout : integer*4
  !
  ! juscount : .true. to just count the particles to be selected. Needed
  !            to allocate appropriately the arrays posout,velout and zout
  !            The parameter npartout becomes then an output given the number
  !            of particles selected.
  !            .false. to perform the actual selection: npartout is an input
  !            then  posout, velout, zout are appropriate outputs.
  !
  ! z1,z2    : the lightcone part of interest is in between z1 and z2, 
  !            with z1 < z2. If we consider a redshift z(t) where all the 
  !            particles are synchrone, and if coarse timestep is
  !            a fixed dt, it is most appropriate to choose z1 and z2 such
  !            that z1=z(t+dt/2) and z2=z(t-dt/2) to have best accuracy.
  !
  ! om0in    : the value of the cosmological parameter omega0 (typically 0.3)
  !
  ! omLin    : the value of the cosmological constant lambda (typically 0.7)
  !
  ! hubin    : the value of H0/100, where H0 is the present Hubble constant
  !            in km/s/Mpc (typically 0.7)
  !
  ! Lbox     : the comoving size of the simulation box in Mpc (NOT in Mpc/h)
  !
  ! observer(3) : the observer position in the box in Mpc, assuming that
  !            coordinates are in [0,Lbox[
  !
  ! thetay   : half the opening angle in degrees of the lightcone along y direction
  !            (it should be obviously smaller than 90 degrees to avoid catastrophic 
  !            behavior). The lightcone is assume to be aligned with x axis (after
  !            appropriates rotations given by angles theta and phi)
  !
  ! thetaz   : half the opening angle in degrees of the lightcone along z direction
  !            Given thetay and thetaz, the area of the survey is thus 4.thetay.thetaz
  !
  ! theta, phi : 2 angles in degrees defining a rotation to avoid alignement of
  !            the lightcone with the major axes of the simulation box.
  !            Example : theta=21, phi=17.
  !
  ! pos(3,npart) : comoving positions of the input particles in Mpc, assumed to be
  !            in [0,Lbox[.
  !
  ! vel(3,npart) : velocities of the input particles (in any unit, it does not
  !            matter)
  !
  ! npart    : number of input particles to be treated
  !
  ! posout(3,npartout) : output comoving positions of selected particles in Mpc.
  !
  ! velout(3,npartout) : output velocities of selected particles
  !
  ! zout(npartout) : output redshift of selected particles
  !
  ! npartout : number of selected particles. To be computed appropriately,
  !            this routine must be called with juscount=.true., which will give
  !            npartout as an output. Then this routine must be called with 
  !            juscount=.false. with the correct value of npartout, after having
  !            allocated correctly arrays posout,velout,zout.
  !===========================================================================
  use amr_parameters, ONLY: nvector,sp
  implicit none
  logical :: justcount,verbose
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
  real(kind=8) :: observer(3),thetay,thetaz,theta,phi
  real(kind=8) :: pos(3,nvector),vel(3,nvector)
  integer(kind=8) :: idtab(nvector)
  real(kind=8) :: posout(3,27*nvector),velout(3,27*nvector),zout(27*nvector)
  integer(kind=8) :: idout(27*nvector)
  real(sp) :: pottab(nvector)
  real(sp) :: potout(27*nvector)
  real(sp) :: fparttab(3,nvector)
  real(sp) :: fpartout(3,27*nvector)
  real(sp) :: fxfr,fyfr,fzfr
  real(kind=8) :: coord_distance
  real(kind=8) :: thetarad,phirad,thetayrad,thetazrad,tanybound,tanzbound
  real(kind=8) :: rot(3,3),rotm1(3,3),dist1,dist2,cosy,cosz
  real(kind=8) :: xcoordfr,ycoordfr,zcoordfr,xcoord,ycoord,zcoord
  real(kind=8) :: tany,tanz,dist,vxfr,vyfr,vzfr,dxtest1,dxtest2,facnorm
  real(kind=8) :: pi
  real(kind=8) :: small=1d-5
  real(kind=8) :: dist1cone, dist2cone,observer_redshift
  
  integer :: nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
  integer :: i,j,k,np,npartcount
  integer :: npart,npartout
  
  if(verbose) write(*,*) 'Entering perform_my_selection'
  
  ! pi=3.14159
  pi=dacos(-1.0d0)
  
  ! Initialize cosmological parameters
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  
  ! Convert angles in radians
  thetarad=theta*pi/180.0d0
  phirad  =phi  *pi/180.0d0
  
  ! Compute the rotation matrix and its inverse to be in the appropriate frame
  call compute_rotation_matrix(thetarad,phirad,rot,rotm1)
  
  ! Compute comoving distance of the photon planes from the observer
  ! dist1,dist2=integral of c.dt/a between zero and z1,z2
  if(observer_redshift.gt.0) then
     dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0) ! VR
     dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)-coord_distance(observer_redshift,Omega0,OmegaL,OmegaR,coverH0) ! VR
  else 
     dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0) ! VR
     dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0) ! VR
  end if
  
  dist1cone = dist1/Lbox                                                ! distance to the closer  side of the shell (dimensionleass)
  dist2cone = dist2/Lbox                                                ! distance to the farther side of the shell (dimensionless)
  
  ! Convert angles in radians
  thetayrad=thetay*pi/180.0d0
  thetazrad=thetaz*pi/180.0d0
  
  ! Compute the set of replica to be considered
  call compute_replica(thetayrad,thetazrad,dist1,dist2,observer,Lbox,rot, &
       &                       nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp)    
  
  facnorm=1.0d0/(dist2-dist1)
  tanybound=dtan(thetayrad)
  tanzbound=dtan(thetazrad)
  
  npartcount=0    
  ! loop on all the replica of potential interest
  do k=nrepzm,nrepzp,1
     do j=nrepym,nrepyp,1
        do i=nrepxm,nrepxp,1
           do np=1,npart
              xcoordfr=pos(1,np)+Lbox*dble(i)-observer(1)               ! x distance to observer
              ycoordfr=pos(2,np)+Lbox*dble(j)-observer(2)               ! y distance to observer
              zcoordfr=pos(3,np)+Lbox*dble(k)-observer(3)               ! z distance to observer
              
              ! Rotation to get in the framework of the photon plane
              xcoord=xcoordfr*rotm1(1,1)+ &
                   & ycoordfr*rotm1(2,1)+ &
                   & zcoordfr*rotm1(3,1)                                ! rotated x coordinate relative to observer
              ycoord=xcoordfr*rotm1(1,2)+ &
                   & ycoordfr*rotm1(2,2)+ &
                   & zcoordfr*rotm1(3,2)                                ! rotated y coordinate relative to observer
              zcoord=xcoordfr*rotm1(1,3)+ &
                   & ycoordfr*rotm1(2,3)+ &
                   & zcoordfr*rotm1(3,3)                                ! rotated z coordinate relative to observer
              
              if(xcoord>small) then                                     ! to avoid divergences near the origin
                 tany=dabs(ycoord/xcoord)
                 tanz=dabs(zcoord/xcoord)
                 dist=dsqrt(xcoord**2+ycoord**2+zcoord**2)
                 if((tany.le.tanybound).and.(tanz.le.tanzbound) &
                      &  .and.(dist.gt.dist1).and.(dist.le.dist2)) then ! this particle is good, we can add it to the list
                    !
                    npartcount=npartcount+1                             ! increment total number of selected particles
                    !
                    posout(1,npartcount)=xcoord                         ! x coordinate relative to cone vertex
                    posout(2,npartcount)=ycoord                         ! y coordinate relative to cone vertex
                    posout(3,npartcount)=zcoord                         ! z coordinate relative to cone vector
                    idout (npartcount)=idtab(np)
#ifdef OUTPUT_PARTICLE_POTENTIAL
                    potout(npartcount)=pottab(np)                       ! potential
#endif  
                    ! Velocities are rotated
                    vxfr=vel(1,np)                                      ! x component of unrotated velocity
                    vyfr=vel(2,np)                                      ! y component of unrotated velocity
                    vzfr=vel(3,np)                                      ! z component of unrotated velocity
                    velout(1,npartcount)=vxfr*rotm1(1,1)+ &
                         &               vyfr*rotm1(2,1)+ &
                         &               vzfr*rotm1(3,1)                ! x component of rotated velocity
                    velout(2,npartcount)=vxfr*rotm1(1,2)+ &
                         &               vyfr*rotm1(2,2)+ &
                         &               vzfr*rotm1(3,2)                ! y component of rotated velocity
                    velout(3,npartcount)=vxfr*rotm1(1,3)+ &
                         &               vyfr*rotm1(2,3)+ &
                         &               vzfr*rotm1(3,3)                ! z component of rotated velocity
                    ! Forces are rotated
#ifdef WITHPARTFORCE
                    fxfr=fparttab(1,np)                                 ! x component of unrotated force
                    fyfr=fparttab(2,np)                                 ! y component of unrotated force
                    fzfr=fparttab(3,np)                                 ! z component of unrotated force
                    fpartout(1,npartcount)=fxfr*rotm1(1,1)+ &
                         &                 fyfr*rotm1(2,1)+ &
                         &                 fzfr*rotm1(3,1)              ! x component of rotated force
                    fpartout(2,npartcount)=fxfr*rotm1(1,2)+ &
                         &                 fyfr*rotm1(2,2)+ &
                         &                 fzfr*rotm1(3,2)              ! y component of rotated force
                    fpartout(3,npartcount)=fxfr*rotm1(1,3)+ &
                         &                 fyfr*rotm1(2,3)+ &
                         &                 fzfr*rotm1(3,3)              ! z component of rotated force
                    
#endif
                    ! Compute redshift of particle by linear interpolation
                    dxtest1=dist-dist1
                    dxtest2=dist2-dist
                    zout(npartcount)=(dxtest1*z2+dxtest2*z1)*facnorm
                    !
                 endif
              endif
           enddo
        enddo
     enddo
  enddo
  !
  npartout=npartcount
  !
  if (verbose) write(*,*) 'End of perform_my_selection'

end subroutine perform_my_selection

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
! Compute the rotation matrix rot and its transpose matrix rotm1. rot brings
! a vector along the x-axis to a new vector described by (theta, phi), while
! rotm1 does the opposite. Here phi is the same as usually used in spherical
! polar coordinate system, while theta is pi/2 minus the theta angle usually
! used in spherical polar coordinates. Both angles are expressed in radians.
!--------------------------------------------------------------------------!
subroutine compute_rotation_matrix(thetashiftrad,phishiftrad,rot,rotm1)
  implicit none
  real(kind=8) :: thetashiftrad,phishiftrad
  real(kind=8) :: rot(3,3),rotm1(3,3)

  integer :: i,j

  rot(1,1) =  dcos(thetashiftrad)*dcos(phishiftrad)
  rot(1,2) =  dcos(thetashiftrad)*dsin(phishiftrad)
  rot(1,3) = -dsin(thetashiftrad)
  rot(2,1) = -dsin(phishiftrad)
  rot(2,2) =  dcos(phishiftrad)
  rot(2,3) = 0.0d0
  rot(3,1) = dcos(phishiftrad)*dsin(thetashiftrad)
  rot(3,2) = dsin(phishiftrad)*dsin(thetashiftrad)
  rot(3,3) = dcos(thetashiftrad)

  do j=1,3
     do i=1,3
        rotm1(i,j)=rot(j,i)
     enddo
  enddo

end subroutine compute_rotation_matrix

!---------------------------------------------------------------------------!
!---------------------------------------------------------------------------!
! A slice of photons between redshifts z1 and z2 corresponding to coordinates
! x1 and x2 at its center and of opening angle (thetay, thetaz) is considered.
!
! We obtain the coordinates of the 8 vertices of the mimimum (simple) polygon
! containing it. The centre of the lighcone is taken to be alighed to x-axis.
!---------------------------------------------------------------------------!
subroutine compute_minimum_polygon(x1,x2,thetayrad,thetazrad,sl)
  implicit none
  real(kind=8)::x1,x2,thetayrad,thetazrad,sl(3,8)

  ! Part of the polygon close to the observer
  sl(1,1:4)=x1/dsqrt(1.0D0+dtan(thetayrad)**2+dtan(thetazrad)**2)       ! x coordinate
  sl(2,1)=-sl(1,1)*dtan(thetayrad)                                      ! y coordinate of vertex 1
  sl(3,1)=-sl(1,1)*dtan(thetazrad)                                      ! z coordinate of vertex 1
  sl(2,2)= sl(1,1)*dtan(thetayrad)                                      ! y coordinate of vertex 2
  sl(3,2)=-sl(1,1)*dtan(thetazrad)                                      ! z coordinate of vertex 2
  sl(2,3)=-sl(1,1)*dtan(thetayrad)                                      ! y coordinate of vertex 3
  sl(3,3)= sl(1,1)*dtan(thetazrad)                                      ! z coordinate of vertex 3
  sl(2,4)= sl(1,1)*dtan(thetayrad)                                      ! y coordinate of vertex 4
  sl(3,4)= sl(1,1)*dtan(thetazrad)                                      ! z coordinate of vertex 4

  ! Part of the polygon far away from the observer
  sl(1,5:8)= x2                                                         ! x coordinate
  sl(2,5)  =-x2*dtan(thetayrad)                                         ! y coordinate of vertex 5
  sl(3,5)  =-x2*dtan(thetazrad)                                         ! z coordinate of vertex 5
  sl(2,6)  = x2*dtan(thetayrad)                                         ! y coordinate of vertex 6
  sl(3,6)  =-x2*dtan(thetazrad)                                         ! z coordinate of vertex 6
  sl(2,7)  =-x2*dtan(thetayrad)                                         ! y coordinate of vertex 7
  sl(3,7)  = x2*dtan(thetazrad)                                         ! z coordinate of vertex 7
  sl(2,8)  = x2*dtan(thetayrad)                                         ! y coordinate of vertex 8
  sl(3,8)  = x2*dtan(thetazrad)                                         ! z coordinate of vertex 8

end subroutine compute_minimum_polygon

!--------------------------------------------------------------------------! 
!--------------------------------------------------------------------------! 
! Calculate the numbers of replics needed to contain the entire photon shell
!
! two*theta1, two*theta2 are the opening angles of the lightcone in degrees.
! The observer's position is expressed in comoving Mpc, so is the simulation
! boxsize (Lbox). Moreover, the positions of particles inside the simulation
! box are supposed to be in [0,Lbox], and z1 and z2 are the redshifts of the 
! successive photon planes, with z1 < z2
!--------------------------------------------------------------------------!
subroutine compute_replica(thetayrad,thetazrad,dist1,dist2,observer,Lbox,rot, &
     &                     nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp)
  implicit none
  real(kind=8) :: thetayrad,thetazrad,observer(3),Lbox,rot(3,3),dist1,dist2
  integer :: nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
  integer :: myint
  real(kind=8) :: sl(3,8),slfr(3)
  real(kind=8) :: xplmin=0.0d0,xplmax=0.0d0,yplmin=0.0d0,yplmax=0.0d0,zplmin=0.0d0,zplmax=0.0d0
  integer :: i,j

  ! compute the minimum polygon containing the two shells of photons
  call compute_minimum_polygon(dist1,dist2,thetayrad,thetazrad,sl)

  ! rotate the minimum polygon in the reference frame of the simulation
  do j=1,8
     do i=1,3
        slfr(i) = sl(1,j)*rot(1,i)+ &
                & sl(2,j)*rot(2,i)+ &
                & sl(3,j)*rot(3,i)                                      ! rotate to the physical angles
     enddo
     !
     if(j.eq.1) then
        xplmin=slfr(1)
        yplmin=slfr(2)
        zplmin=slfr(3)
        xplmax=xplmin
        yplmax=yplmin
        zplmax=zplmin
     else
        xplmin=dmin1(xplmin,slfr(1))                                    ! find minimum x coordinate of all 8 vertices
        xplmax=dmax1(xplmax,slfr(1))                                    ! find maximum x coordinate of all 8 vertices
        yplmin=dmin1(yplmin,slfr(2))                                    ! find minimum y coordinate of all 8 vertices
        yplmax=dmax1(yplmax,slfr(2))                                    ! find maximum y coordinate of all 8 vertices
        zplmin=dmin1(zplmin,slfr(3))                                    ! find minimum z coordinate of all 8 vertices
        zplmax=dmax1(zplmax,slfr(3))                                    ! find maximum z coordinate of all 8 vertices
     endif
  enddo

  ! Uses the fact that a cube will contain the minimum polygon if and only
  ! if all its edges are contained in the cube to compute relevant replica
  nrepxm=myint((xplmin+observer(1))/Lbox)                               ! index of leftmost  replica
  nrepxp=myint((xplmax+observer(1))/Lbox)                               ! index of rightmost replica
  nrepym=myint((yplmin+observer(2))/Lbox)                               ! index of front     replica
  nrepyp=myint((yplmax+observer(2))/Lbox)                               ! index of back      replica
  nrepzm=myint((zplmin+observer(3))/Lbox)                               ! index of top       replica
  nrepzp=myint((zplmax+observer(3))/Lbox)                               ! index of bottom    replica

end subroutine compute_replica

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
! Subroutine to write useful information
!--------------------------------------------------------------------------!
subroutine write_infoconepart(npart_tot,conepartid,conepartzlim,observer_x,observer_y,observer_z,observer_redshift, &
                            & amax,amin,zmax,zmin,dmax,dmin,dtol,                                                   &
                            & reduce_infoglobal,                                                                    &
                            & isfullsky,input_thetay,input_thetaz,input_theta,input_phi,                            &
                            & aendconem2_info,aendconem1_info,aendcone_info,aexpold_info,aexp_info,                 &
                            & zendconem2_info,zendconem1_info,zendcone_info,zexpold_info,zexp_info,                 &
                            & dendconem2_info,dendconem1_info,dendcone_info,dexpold_info,dexp_info,future)
  use amr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  
  ! Variables from input
  integer::npart_tot
  integer::conepartid
  real(kind=8)::conepartzlim
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
  integer(kind=8),dimension(1:2)::reduce_infoglobal
  logical::isfullsky
  integer::isfullskyint
  real(kind=8)::input_thetay
  real(kind=8)::input_thetaz
  real(kind=8)::input_theta
  real(kind=8)::input_phi

  ! Local variables
  integer::ilun
  character(LEN=5)::nchar
  character(LEN=200)::filename
  character(LEN=200)::filedir
  character(LEN=5)::nchar_cone
  
  real(kind=8)::aendconem2_info,aendconem1_info,aendcone_info,aexpold_info,aexp_info
  !aendconem* -> expansion factor for end of the shell
  !the end of previous shell and the end of previous previous shell (which is the begining of the current shell)
  !when activating overlapping buffer option.
  real(kind=8)::zendconem2_info,zendconem1_info,zendcone_info,zexpold_info,zexp_info
  real(kind=8)::dendconem2_info,dendconem1_info,dendcone_info,dexpold_info,dexp_info

  integer::future

  ! Write
  ilun=40*ncpu+myid+10
  call title(nstep_coarse,nchar)
  write(nchar_cone,'(I5.5)') conepartid
  filedir='output_ncoarse_'//TRIM(nchar)//'/'
  if (isfullsky) then
     isfullskyint=1
     if(future==1) then 
        filename=TRIM(filedir)//'info_cone_part_fullsky_future_'//nchar_cone//'_ncoarse_'//TRIM(nchar)//'.txt'
     else if (future==-1)then
        filename=TRIM(filedir)//'info_cone_part_fullsky_past_'//nchar_cone//'_ncoarse_'//TRIM(nchar)//'.txt'
     else
        print*,'not implemented cone part future=',future
        stop
     endif
  else
     isfullskyint=0
     if(future==1)then
        filename=TRIM(filedir)//'info_cone_part_narrow_future_'//nchar_cone//'_ncoarse_'//TRIM(nchar)//'.txt'
     else if (future==-1)then
        filename=TRIM(filedir)//'info_cone_part_narrow_past_'//nchar_cone//'_ncoarse_'//TRIM(nchar)//'.txt'
     else
        print*,'not implemented cone part future=',future
        stop
     endif
  endif
  open(ilun,file=TRIM(filename),form='formatted')
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("nstride     =",I11)')nstride
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,'("aexp        =",E23.15)')aexp
  write(ilun,'("observer_x  =",E23.15)')observer_x
  write(ilun,'("observer_y  =",E23.15)')observer_y
  write(ilun,'("observer_z  =",E23.15)')observer_z
  write(ilun,'("observer_rds=",E23.15)')observer_redshift
  write(ilun,'("cone_id     =",I11)')conepartid
  write(ilun,'("cone_zlim   =",E23.15)')conepartzlim
  write(ilun,'("amax        =",E23.15)')amax
  write(ilun,'("amin        =",E23.15)')amin
  write(ilun,'("zmax        =",E23.15)')zmax
  write(ilun,'("zmin        =",E23.15)')zmin
  write(ilun,'("dmax        =",E23.15)')dmax
  write(ilun,'("dmin        =",E23.15)')dmin
  write(ilun,'("dtol        =",E23.15)')dtol
  write(ilun,'("nglobalfile =",I20)')reduce_infoglobal(1)
  write(ilun,'("nglobalcell =",I20)')reduce_infoglobal(2)
  write(ilun,'("isfullsky   =",I11)')isfullskyint
  write(ilun,'("thetay      =",E23.15)')input_thetay
  write(ilun,'("thetaz      =",E23.15)')input_thetaz
  write(ilun,'("theta       =",E23.15)')input_theta
  write(ilun,'("phi         =",E23.15)')input_phi 
  write(ilun,'("aendconem2  =",E23.15)')aendconem2_info
  write(ilun,'("aendconem1  =",E23.15)')aendconem1_info
  write(ilun,'("aendcone    =",E23.15)')aendcone_info
  write(ilun,'("aexpold     =",E23.15)')aexpold_info
  write(ilun,'("aexp        =",E23.15)')aexp_info
  write(ilun,'("zendconem2  =",E23.15)')zendconem2_info
  write(ilun,'("zendconem1  =",E23.15)')zendconem1_info
  write(ilun,'("zendcone    =",E23.15)')zendcone_info
  write(ilun,'("zexpold     =",E23.15)')zexpold_info
  write(ilun,'("zexp        =",E23.15)')zexp_info
  write(ilun,'("dendconem2  =",E23.15)')dendconem2_info
  write(ilun,'("dendconem1  =",E23.15)')dendconem1_info
  write(ilun,'("dendcone    =",E23.15)')dendcone_info
  write(ilun,'("dexpold     =",E23.15)')dexpold_info
  write(ilun,'("dexp        =",E23.15)')dexp_info
  write(ilun,'("future      =",I11)')future
  close(ilun)
  
end subroutine write_infoconepart
!===========================================================================


!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
! Subroutine to write useful information
!--------------------------------------------------------------------------!
subroutine write_infosamplepart(npart_tot, aexp_sample, &
                        & xmin, xmax, ymin, ymax, zmin, zmax, &
                        & reduce_infoglobal)
  use amr_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  
  ! Variables from input
  integer::npart_tot
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
    
  ! Local variables
  integer::ilun
  character(LEN=5)::nchar
  character(LEN=200)::filename
  character(LEN=200)::filedir
  character::nchar_sample
  
  ! Write
  ilun=40*ncpu+myid+10
  call title(nstep_coarse,nchar)
  filedir='output_ncoarse_'//TRIM(nchar)//'/'
  filename=TRIM(filedir)//'info_sample_part_ncoarse_'//TRIM(nchar)//'.txt'
  open(ilun,file=TRIM(filename),form='formatted')
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("nstride     =",I11)')nstride
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,'("aexp        =",E23.15)')aexp_sample
  write(ilun,'("xmin        =",E23.15)')xmin
  write(ilun,'("xmax        =",E23.15)')xmax
  write(ilun,'("ymin        =",E23.15)')ymin
  write(ilun,'("ymax        =",E23.15)')ymax
  write(ilun,'("zmin        =",E23.15)')zmin
  write(ilun,'("zmax        =",E23.15)')zmax
  write(ilun,'("nglobalfile =",I20)')reduce_infoglobal(1)
  write(ilun,'("nglobalcell =",I20)')reduce_infoglobal(2)
  close(ilun)
  
end subroutine write_infosamplepart



!------------------ MODIF V. REVERDY 2012 ------------------! 

!-----------------------------------------------------------! 
!--------------------- CBH_LC 09-02-2021 -------------------!
!-----------------------------------------------------------! 

!===================
!cone cosmo routines
!===================

!===========================================================================
subroutine init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  !=========================================================================
  ! om0in : the value of omega0
  ! omLin : the value of Lambda
  !         We MUST have omega0+Lambda=1.0d0
  ! hubin : the value of H0/100 where H0 is the present Hubble constant
  !         in km/s/Mpc
  !=========================================================================
  implicit none
  real(kind=8) :: om0in,omLin,hubin
  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
  real(kind=8) :: verysmall=1d-6
  !
  omega0=om0in
  omegaL=omLin
  omegaR=1.0d0-omega0-omegaL
  !
  if(dabs(omegaR)>verysmall) then
     write(*,*) 'ERROR in propagate_photons, init_cosmo.'
     write(*,*) 'This routine works only for flat universes, omega0+Lambda=1.'
     STOP
  endif
  !
  coverH0=299792.5d0/(100.0d0*hubin)
  !
end subroutine init_cosmo_cone

!===========================================================================
function coord_distance(zz,Omega0,OmegaL,OmegaR,coverH0)
  !=========================================================================
  use amr_commons
  implicit none
  !real(kind=8) :: z,res,coord_distance,zz
  real(kind=8) :: z,res,coord_distance,zz,atmp
  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
  integer::i

  z=dabs(zz)
  atmp=1.D0/(z+1.D0)

  i=1
  do while(aexp_frw(i)>atmp.and.i<n_frw)
     i=i+1
  end do

  ! Interpolate proper time
  res = tprop_frw(i  )*(atmp-aexp_frw(i-1))/(aexp_frw(i  )-aexp_frw(i-1))+ &
      & tprop_frw(i-1)*(atmp-aexp_frw(i  ))/(aexp_frw(i-1)-aexp_frw(i  ))

!  call qromb(0.d0,z,res,omega0,omegaL,OmegaR)
!  coord_distance=coverH0*res
  coord_distance=-coverH0*res
  !
  if(zz.lt.0.0D0) coord_distance=-coord_distance
end function coord_distance

!===========================================================================
function funcE(z,Omega0,OmegaL,OmegaR)
  !=========================================================================
  implicit none
  real(kind=8) :: funcE,z,HsurH0
  real(kind=8) :: omega0,omegaL,OmegaR

  funcE=1.D0/HsurH0(z,Omega0,OmegaL,OmegaR)

end function funcE

!===========================================================================
function HsurH0(z,omega0,omegaL,OmegaR)
  !=========================================================================
  implicit none
  real(kind=8) :: z,omega0,omegaL,OmegaR,HsurH0

  HsurH0=dsqrt(Omega0*(1.D0+z)**3+OmegaR*(1.D0+z)**2+OmegaL)

end function HsurH0

!===========================================================================
SUBROUTINE qromb(a,b,ss,omega0,omegaL,OmegaR)
  !===========================================================================
  implicit none
  INTEGER :: JMAX,JMAXP,K,KM
  REAL(kind=8) :: a,b,ss,EPS,omega0,omegaL,OmegaR
  PARAMETER (EPS=1.d-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
  !  USES polint,trapzd
  INTEGER :: j
  REAL(kind=8) :: dss,h(JMAXP),s(JMAXP)

  h(1)=1.D0
  do j=1,JMAX
     call trapzd(a,b,s(j),j,omega0,omegaL,OmegaR)
     if (j.ge.K) then
        call polint(h(j-KM),s(j-KM),K,0.D0,ss,dss)
        if(dabs(dss).le.EPS*dabs(ss)) return
     endif
     s(j+1)=s(j)
     h(j+1)=0.25D0*h(j)
  enddo

  print *, 'too many steps in qromb'
END SUBROUTINE qromb

!===========================================================================
SUBROUTINE polint(xa,ya,n,x,y,dy)
  !===========================================================================
  implicit none
  INTEGER :: n,NMAX
  REAL(kind=8) :: dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER :: i,m,ns
  REAL(kind=8) ::den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  !
  ns=1
  dif=dabs(x-xa(1))
  do i=1,n
     dift=dabs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  !
  y=ya(ns)
  ns=ns-1
  !
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.D0) print *, 'failure in polint'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  !
  return

END SUBROUTINE polint

!===========================================================================
SUBROUTINE trapzd(a,b,s,n,omega0,omegaL,OmegaR)
  !===========================================================================
  implicit none
  INTEGER :: n
  REAL(kind=8) :: a,b,s,funcE,omega0,omegaL,OmegaR
  INTEGER :: it,j
  REAL(kind=8) :: del,sum,tnm,x

  if (n.eq.1) then
     s=0.5D0*(b-a)*(funcE(a,omega0,omegaL,OmegaR)+funcE(b,omega0,omegaL,OmegaR))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5D0*del
     sum=0.D0
     do j=1,it
        sum=sum+funcE(x,omega0,omegaL,OmegaR)
        x=x+del
     enddo
     s=0.5D0-*(s+(b-a)*sum/tnm)
  endif
  return
END SUBROUTINE trapzd

!=======================================================================
!=======================================================================
! DEFAULT RAMSES LIGHTCONE SUBROUTINES ! CBH_LC
!=======================================================================
!=======================================================================

!subroutine output_cone()
!  use amr_commons
!  use pm_commons
!  implicit none
!
!#ifndef WITHOUTMPI
!#include "mpif.h"
!  integer::info,info2,dummy_io
!#endif
!
!  integer,parameter::tag=1118
!
!  character(len=5) :: istep_str
!  character(len=100) :: conedir, conecmd, conefile
!
!  integer::ilun,ipout,npout,npart_out
!  character(LEN=80)::fileloc
!  character(LEN=5)::nchar
!  real(kind=8),dimension(1:3,1:nvector),save::pos,vel
!  real(kind=8),dimension(:,:),allocatable::posout,velout
!  real(kind=8),dimension(:),allocatable::zout
!  real(kind=8),dimension(:,:),allocatable::tmparr
!  real(sp),dimension(:,:),allocatable::xp_out,vp_out
!  real(sp),dimension(:),allocatable::zp_out
!  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
!  real(kind=8) :: observer(3),thetay,thetaz,theta,phi
!  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel
!  integer::i,ip,npart1
!  integer::nalloc1,nalloc2
!
!  integer,dimension(1:nvector),save::ind_part
!  logical::opened
!  opened=.false.
!
!  if(nstep_coarse.lt.2)return
!
!  z2=1./aexp_old-1.
!  z1=1./aexp-1.
!
!  if(z2.gt.zmax_cone)return
!  if(abs(z2-z1)<1d-6)return
!
!  theta=25.
!  phi=17.
!  thetay=thetay_cone
!  thetaz=thetaz_cone
!  om0in=omega_m
!  omLin=omega_l
!  hubin=h0/100.
!  Lbox=boxlen_ini/hubin
!  observer=(/Lbox/2.0,Lbox/2.0,Lbox/2.0/)
!
!  ilun=3*ncpu+myid+10
!
!  ! Determine the filename, dir, etc
!  if(myid==1)write(*,*)'Computing and dumping lightcone'
!
!  call title(nstep_coarse, istep_str)
!  conedir = "cone_" // trim(istep_str) // "/"
!  conecmd = "mkdir -p " // trim(conedir)
!  if(.not.withoutmkdir) then
!     if (myid==1) call system(conecmd)
!  endif
!
!#ifndef WITHOUTMPI
!  call MPI_BARRIER(MPI_COMM_WORLD, info)
!#endif
!
!  conefile = trim(conedir)//'cone_'//trim(istep_str)//'.out'
!  call title(myid,nchar)
!  fileloc=TRIM(conefile)//TRIM(nchar)
!
!  npart_out=0
!  ipout=0
!  npout=0
!
!  ! Pre-allocate arrays for particle selection -----
!  nalloc1=nvector
!  allocate(posout(1:3, 1:nalloc1))
!  allocate(velout(1:3, 1:nalloc1))
!  allocate(zout(1:nalloc1))
!
!  nalloc2=nvector+nstride
!  allocate(xp_out(1:nalloc2,1:3))
!  allocate(vp_out(1:nalloc2,1:3))
!  allocate(zp_out(1:nalloc2))
!
!  allocate(tmparr(1:3, 1:nalloc2))
!  ! ------------------------------------------------
!
!  ! Wait for the token
!#ifndef WITHOUTMPI
!  if(IOGROUPSIZECONE>0) then
!     if (mod(myid-1,IOGROUPSIZECONE)/=0) then
!        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
!             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
!     end if
!  endif
!#endif
!
!
!  ilevel=levelmin
!  ! Loop over cpus
!  do icpu=1,ncpu
!     ! Loop over grids
!     igrid=headl(icpu,ilevel)
!     ip=0
!     do jgrid=1,numbl(icpu,ilevel)
!        npart1=numbp(igrid)  ! Number of particles in the grid
!        if(npart1>0)then
!           ipart=headp(igrid)
!
!           ! Loop over particles
!           do jpart=1,npart1
!              ip=ip+1
!              ind_part(ip)=ipart
!              if(ip==nvector)then
!                 ! Lower left corner of 3x3x3 grid-cube
!                 do idim=1,ndim
!                    do i=1,ip
!                       pos(idim,i)=xp(ind_part(i),idim)*Lbox
!                       vel(idim,i)=vp(ind_part(i),idim)
!                    end do
!                 end do
!                 !===========================================================================
!                 ! Count selection particles
!                 call perform_my_selection(.true.,z1,z2, &
!                      &                           om0in,omLin,hubin,Lbox, &
!                      &                           observer,thetay,thetaz,theta,phi, &
!                      &                           pos,vel,ip, &
!                      &                           posout,velout,zout,npout,.false.)
!
!
!                 call extend_arrays_if_needed()
!
!
!                 ! Perform actual selection
!                 call perform_my_selection(.false.,z1,z2, &
!                      &                           om0in,omLin,hubin,Lbox, &
!                      &                           observer,thetay,thetaz,theta,phi, &
!                      &                           pos,vel,ip, &
!                      &                           posout,velout,zout,npout,.false.)
!                 !===========================================================================
!                 if(npout>0)then
!                    do idim=1,ndim
!                       do i=1,npout
!                          xp_out(ipout+i,idim)=real(posout(idim,i)/Lbox,kind=sp)
!                          vp_out(ipout+i,idim)=real(velout(idim,i),kind=sp)
!                       end do
!                    end do
!                    do i=1,npout
!                       zp_out(ipout+i)=real(zout(i),kind=sp)
!                    end do
!                    ipout=ipout+npout
!                    npart_out=npart_out+npout
!                 endif
!
!                 ip=0
!              end if
!              if(ipout>=nstride)then
!                 if(.not.opened) then
!                    open(ilun,file=TRIM(fileloc),form='unformatted')
!                    rewind(ilun)
!                    write(ilun)ncpu
!                    write(ilun)nstride
!                    write(ilun)npart
!                    opened=.true.
!                 endif
!                 do idim=1,ndim
!                    write(ilun)xp_out(1:nstride,idim)
!                    write(ilun)vp_out(1:nstride,idim)
!                 end do
!                 write(ilun)zp_out(1:nstride)
!                 do idim=1,ndim
!                    do i=1,ipout-nstride
!                       xp_out(i,idim)=xp_out(i+nstride,idim)
!                       vp_out(i,idim)=vp_out(i+nstride,idim)
!                    end do
!                 end do
!                 do i=1,ipout-nstride
!                    zp_out(i)=zp_out(i+nstride)
!                 end do
!                 ipout=ipout-nstride
!              endif
!              ipart=nextp(ipart)  ! Go to next particle
!           end do
!           ! End loop over particles
!        end if
!        igrid=next(igrid)   ! Go to next grid
!     end do
!     ! End loop over grids
!
!     if(ip>0)then
!        ! Lower left corner of 3x3x3 grid-cube
!        do idim=1,ndim
!           do i=1,ip
!              pos(idim,i)=xp(ind_part(i),idim)*Lbox
!              vel(idim,i)=vp(ind_part(i),idim)
!           end do
!        end do
!        !===========================================================================
!        ! Count selection particles
!        call perform_my_selection(.true.,z1,z2, &
!             &                           om0in,omLin,hubin,Lbox, &
!             &                           observer,thetay,thetaz,theta,phi, &
!             &                           pos,vel,ip, &
!             &                           posout,velout,zout,npout,.false.)
!
!        call extend_arrays_if_needed()
!
!        ! Perform actual selection
!        call perform_my_selection(.false.,z1,z2, &
!             &                           om0in,omLin,hubin,Lbox, &
!             &                           observer,thetay,thetaz,theta,phi, &
!             &                           pos,vel,ip, &
!             &                           posout,velout,zout,npout,.false.)
!        !===========================================================================
!        if(npout>0)then
!           do idim=1,ndim
!              do i=1,npout
!                 xp_out(ipout+i,idim)=real(posout(idim,i)/Lbox,kind=sp)
!                 vp_out(ipout+i,idim)=real(velout(idim,i),kind=sp)
!              end do
!           end do
!           do i=1,npout
!              zp_out(ipout+i)=real(zout(i),kind=sp)
!           end do
!           ipout=ipout+npout
!           npart_out=npart_out+npout
!        endif
!     endif
!     if(ipout>=nstride)then
!        if(.not.opened) then
!           open(ilun,file=TRIM(fileloc),form='unformatted')
!           rewind(ilun)
!           write(ilun)ncpu
!           write(ilun)nstride
!           write(ilun)npart
!           opened=.true.
!        endif
!        do idim=1,ndim
!           write(ilun)xp_out(1:nstride,idim)
!           write(ilun)vp_out(1:nstride,idim)
!        end do
!        write(ilun)zp_out(1:nstride)
!        do idim=1,ndim
!           do i=1,ipout-nstride
!              xp_out(i,idim)=xp_out(i+nstride,idim)
!              vp_out(i,idim)=vp_out(i+nstride,idim)
!           end do
!        end do
!        do i=1,ipout-nstride
!           zp_out(i)=zp_out(i+nstride)
!        end do
!        ipout=ipout-nstride
!     endif
!  end do
!  ! End loop over cpus
!
!  if(ipout>0)then
!     if(.not.opened) then
!        open(ilun,file=TRIM(fileloc),form='unformatted')
!        rewind(ilun)
!        write(ilun)ncpu
!        write(ilun)nstride
!        write(ilun)npart
!        opened=.true.
!     endif
!     do idim=1,ndim
!        write(ilun)xp_out(1:ipout,idim)
!        write(ilun)vp_out(1:ipout,idim)
!     end do
!     write(ilun)zp_out(1:ipout)
!  endif
!
!  if(opened)close(ilun)
!
!  if (verbose)write(*,*)'cone output=',myid,npart_out
!
!  if(npart_out>0) then
!     open(ilun,file=TRIM(fileloc)//".txt",form='formatted')
!     rewind(ilun)
!     write(ilun,*) ncpu
!     write(ilun,*) nstride
!     write(ilun,*) npart_out
!     write(ilun,*) aexp_old
!     write(ilun,*) aexp
!     close(ilun)
!  endif
!
!     ! Send the token
!#ifndef WITHOUTMPI
!  if(IOGROUPSIZECONE>0) then
!     if(mod(myid,IOGROUPSIZECONE)/=0 .and.(myid.lt.ncpu))then
!        dummy_io=1
!        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
!             & MPI_COMM_WORLD,info2)
!     end if
!  endif
!#endif
!
!
!   if((opened.and.(npart_out==0)).or.((.not.opened).and.(npart_out>0))) then
!     write(*,*)'Error in output_cone'
!     write(*,*)'npart_out=',npart_out,'opened=',opened
!     stop
!  endif
!
!
!contains
!
!    ! Extends (deallocates and reallocates) the arrays
!    ! posout, velout, zout, xp_out, vp_out and zp_out
!    ! after npout has been updated, so they can hold enough particles
!    !
!    ! Reallocation is done in chunks of size alloc_chunk_size, to avoid
!    ! reallocating too frequently.
!
!    subroutine extend_arrays_if_needed()
!
!        ! Allocation chunk size
!        integer, parameter :: alloc_chunk_size = 100
!        integer :: new_nalloc1, new_nalloc2
!        integer :: nchunks1, nchunks2
!
!        if (nalloc1 >= npout .and. nalloc2 >= npout+nstride) return
!
!
!        ! Compute new array sizes
!        nchunks1 = npout / alloc_chunk_size
!        if (mod(npout, alloc_chunk_size) > 0) nchunks1=nchunks1+1
!
!        nchunks2 = (npout+nstride) / alloc_chunk_size
!        if (mod(npout+nstride, alloc_chunk_size) > 0) nchunks2=nchunks2+1
!
!        new_nalloc1 = nchunks1 * alloc_chunk_size
!        new_nalloc2 = nchunks2 * alloc_chunk_size
!
!        ! Resize temp array
!        deallocate(tmparr)
!        allocate(tmparr(1:3,1:max(new_nalloc1,new_nalloc2)))
!
!
!        ! Resize xp_out, vp_out, zp_out
!        do idim=1,ndim
!            tmparr(idim,1:nalloc2)=xp_out(1:nalloc2,idim)
!        end do
!        deallocate(xp_out); allocate(xp_out(1:new_nalloc2,1:3))
!        do idim=1,ndim
!            xp_out(1:nalloc2,idim)=real(tmparr(idim,1:nalloc2),kind=sp)
!        end do
!
!        do idim=1,ndim
!            tmparr(idim,1:nalloc2)=vp_out(1:nalloc2,idim)
!        end do
!        deallocate(vp_out); allocate(vp_out(1:new_nalloc2,1:3))
!        do idim=1,ndim
!            vp_out(1:nalloc2,idim)=real(tmparr(idim,1:nalloc2),kind=sp)
!        end do
!
!        tmparr(1,1:nalloc2)=zp_out(1:nalloc2)
!        deallocate(zp_out); allocate(zp_out(1:new_nalloc2))
!        zp_out(1:nalloc2)=real(tmparr(1,1:nalloc2),kind=sp)
!
!        nalloc2 = new_nalloc2
!
!
!        ! Resize posout, velout, zout
!        do idim=1,ndim
!            tmparr(idim,1:nalloc1)=posout(idim,1:nalloc1)
!        deallocate(posout); allocate(posout(1:3,1:new_nalloc1))
!        end do
!        do idim=1,ndim
!            posout(idim,1:nalloc1)=tmparr(idim,1:nalloc1)
!        end do
!
!        do idim=1,ndim
!            tmparr(idim,1:nalloc1)=velout(idim,1:nalloc1)
!        end do
!        deallocate(velout); allocate(velout(1:3,1:new_nalloc1))
!        do idim=1,ndim
!            velout(idim,1:nalloc1)=tmparr(idim,1:nalloc1)
!        end do
!
!        tmparr(1,1:nalloc1)=zout(1:nalloc1)
!        deallocate(zout); allocate(zout(1:new_nalloc1))
!        zout(1:nalloc1)=tmparr(1,1:nalloc1)
!
!        nalloc1 = new_nalloc1
!
!    end subroutine extend_arrays_if_needed
!end subroutine output_cone

!subroutine perform_my_selection(justcount,z1,z2, &
!     &                          om0in,omLin,hubin,Lbox, &
!     &                          observer,thetay,thetaz,theta,phi, &
!     &                          pos,vel,npart, &
!     &                          posout,velout,zout,npartout,verbose)
!  !===========================================================================
!  ! All the quantities below are real*8 except
!  !      juscount : logical
!  !      npart,npartout : integer*4
!  !
!  ! juscount : .true. to just count the particles to be selected. Needed
!  !            to allocate appropriately the arrays posout,velout and zout
!  !            The parameter npartout becomes then an output given the number
!  !            of particles selected.
!  !            .false. to perform the actual selection: npartout is an input
!  !            then  posout, velout, zout are appropriate outputs.
!  !
!  ! z1,z2    : the lightcone part of interest is in between z1 and z2,
!  !            with z1 < z2. If we consider a redshift z(t) where all the
!  !            particles are synchrone, and if coarse timestep is
!  !            a fixed dt, it is most appropriate to choose z1 and z2 such
!  !            that z1=z(t+dt/2) and z2=z(t-dt/2) to have best accuracy.
!  !
!  ! om0in    : the value of the cosmological parameter omega0 (typically 0.3)
!  !
!  ! omLin    : the value of the cosmological constant lambda (typically 0.7)
!  !
!  ! hubin    : the value of H0/100, where H0 is the present Hubble constant
!  !            in km/s/Mpc (typically 0.7)
!  !
!  ! Lbox     : the comoving size of the simulation box in Mpc (NOT in Mpc/h)
!  !
!  ! observer(3) : the observer position in the box in Mpc, assuming that
!  !            coordinates are in [0,Lbox[
!  !
!  ! thetay   : half the opening angle in degrees of the lightcone along y direction
!  !            (it should be obviously smaller than 90 degrees to avoid catastrophic
!  !            behavior). The lightcone is assume to be aligned with x axis (after
!  !            appropriates rotations given by angles theta and phi)
!  !
!  ! thetaz   : half the opening angle in degrees of the lightcone along z direction
!  !            Given thetay and thetaz, the area of the survey is thus 4.thetay.thetaz
!  !
!  ! theta, phi : 2 angles in degrees defining a rotation to avoid alignement of
!  !            the lightcone with the major axes of the simulation box.
!  !            Example : theta=21, phi=17.
!  !
!  ! pos(3,npart) : comoving positions of the input particles in Mpc, assumed to be
!  !            in [0,Lbox[.
!  !
!  ! vel(3,npart) : velocities of the input particles (in any unit, it does not
!  !            matter)
!  !
!  ! npart    : number of input particles to be treated
!  !
!  ! posout(3,npartout) : output comoving positions of selected particles in Mpc.
!  !
!  ! velout(3,npartout) : output velocities of selected particles
!  !
!  ! zout(npartout) : output redshift of selected particles
!  !
!  ! npartout : number of selected particles. To be computed appropriately,
!  !            this routine must be called with juscount=.true., which will give
!  !            npartout as an output. Then this routine must be called with
!  !            juscount=.false. with the correct value of npartout, after having
!  !            allocated correctly arrays posout,velout,zout.
!  !===========================================================================
!  use amr_parameters, ONLY: nvector
!  implicit none
!  logical :: justcount,verbose
!  integer :: npart,npartout
!  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
!  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
!  real(kind=8) :: observer(3),thetay,thetaz,theta,phi
!  real(kind=8) :: pos(1:3,1:nvector),vel(1:3,1:nvector)
!  real(kind=8) :: posout(3,npartout),velout(3,npartout),zout(npartout)
!  real(kind=8) :: coord_distance
!  real(kind=8) :: thetarad,phirad,thetayrad,thetazrad,tanybound,tanzbound
!  real(kind=8) :: rot(3,3),rotm1(3,3),dist1,dist2
!  real(kind=8) :: xcoordfr,ycoordfr,zcoordfr,xcoord,ycoord,zcoord
!  real(kind=8) :: tany,tanz,dist,vxfr,vyfr,vzfr,dxtest1,dxtest2,facnorm
!  real(kind=8) :: pi
!  real(kind=8) :: small=1d-5
!
!  integer :: nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp
!  integer :: i,j,k,np,npartcount
!
!  if (verbose) write(*,*) 'Entering perform_my_selection'
!
!  ! pi=3.14159
!  pi=acos(-1.0d0)
!
!  ! Initialize cosmological parameters
!  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
!
!  ! Convert angles in radians
!  thetarad=theta*pi/180.0d0
!  phirad=phi*pi/180.0d0
!
!  ! Compute the rotation matrix and its inverse to be in the appropriate frame
!  call compute_rotation_matrix(thetarad,phirad,rot,rotm1)
!
!  ! Compute comoving distance of the photon planes from the observer
!  ! dist1,dist2=integral of c.dt/a between zero and z1,z2
!  dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)
!  dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)
!
!  ! Convert angles in radians
!  thetayrad=thetay*pi/180.0d0
!  thetazrad=thetaz*pi/180.0d0
!
!  ! Compute the set of replica to be considered
!  call compute_replica(thetayrad,thetazrad,dist1,dist2,observer,Lbox,rot, &
!       &                       nrepxm,nrepxp,nrepym,nrepyp,nrepzm,nrepzp)
!
!  facnorm=1.0d0/(dist2-dist1)
!  tanybound=tan(thetayrad)
!  tanzbound=tan(thetazrad)
!
!  npartcount=0
!  ! loop on all the replica of potential interest
!  do k=nrepzm,nrepzp,1
!     do j=nrepym,nrepyp,1
!        do i=nrepxm,nrepxp,1
!           do np=1,npart
!              xcoordfr=pos(1,np)+Lbox*dble(i)-observer(1)
!              ycoordfr=pos(2,np)+Lbox*dble(j)-observer(2)
!              zcoordfr=pos(3,np)+Lbox*dble(k)-observer(3)
!
!              ! Rotation to get in the framework of the photon plane
!              xcoord=xcoordfr*rotm1(1,1)+ &
!                   & ycoordfr*rotm1(2,1)+ &
!                   & zcoordfr*rotm1(3,1)
!              ycoord=xcoordfr*rotm1(1,2)+ &
!                   & ycoordfr*rotm1(2,2)+ &
!                   & zcoordfr*rotm1(3,2)
!              zcoord=xcoordfr*rotm1(1,3)+ &
!                   & ycoordfr*rotm1(2,3)+ &
!                   & zcoordfr*rotm1(3,3)
!
!              if (xcoord > small) then ! To avoid divergences near the origin
!                 tany=abs(ycoord/xcoord)
!                 tanz=abs(zcoord/xcoord)
!                 dist=sqrt(xcoord**2+ycoord**2+zcoord**2)
!                 if (tany <= tanybound .and. tanz <= tanzbound &
!                      &  .and. dist > dist1 .and. dist <= dist2) then
!                    ! This particle is good, we can add it to the list
!                    npartcount=npartcount+1
!
!                    if (.not. justcount) then
!                        posout(1,npartcount)=xcoord
!                        posout(2,npartcount)=ycoord
!                        posout(3,npartcount)=zcoord
!
!                        ! Velocities are rotated
!                        vxfr=vel(1,np)
!                        vyfr=vel(2,np)
!                        vzfr=vel(3,np)
!                        velout(1,npartcount)=vxfr*rotm1(1,1)+ &
!                            &               vyfr*rotm1(2,1)+ &
!                            &               vzfr*rotm1(3,1)
!                        velout(2,npartcount)=vxfr*rotm1(1,2)+ &
!                            &               vyfr*rotm1(2,2)+ &
!                            &               vzfr*rotm1(3,2)
!                        velout(3,npartcount)=vxfr*rotm1(1,3)+ &
!                            &               vyfr*rotm1(2,3)+ &
!                            &               vzfr*rotm1(3,3)
!
!                        ! Compute the redshift of the particle using linear
!                        ! interpolation
!                        dxtest1=dist-dist1
!                        dxtest2=dist2-dist
!                        zout(npartcount)=(dxtest1*z2+dxtest2*z1)*facnorm
!                    endif
!                 endif
!              endif
!           enddo
!        enddo
!     enddo
!  enddo
!  npartout=npartcount
!  if (verbose) write(*,*) 'End of perform_my_selection'
!end subroutine perform_my_selection

!=======================================================================
!=======================================================================
! DEFAULT RAMSES LIGHTCONE SUBROUTINES ! CBH_LC
!=======================================================================
!=======================================================================

