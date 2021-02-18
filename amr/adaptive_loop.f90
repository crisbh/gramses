subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use gr_commons
  use gr_parameters
  use cooling_module
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer(kind=8)::n_step
  integer::info,tot_pt
  real(kind=8)::tt1,tt2,muspt,muspt_this_step,wallsec,dumpsec
  real(kind=4)::real_mem,real_mem_tot
  real(kind=8),save::tstart=0.0
#endif
  integer::ilevel,idim,ivar
  integer::igrp,igrm

!  ! CBH_LC_mem
!  real(kind=4)::real_mem_mean,real_mem_min
!  character(LEN=5)::nchar
!  character(LEN=80)::filename_memstat
!  character(LEN=80)::filename_memall
!  character(LEN=80)::filename_timings
!  real::tt
!  tt = 0.
!  ! END CBH_LC_mem

#ifndef WITHOUTMPI
  tt1=MPI_WTIME()
  ! for calculating total run time
  if (tstart.eq.0.0) then
     tstart = MPI_WTIME(info)
  end if
#endif

  ! CBH_LC
  if(nrestart.EQ.0) then
    writencoarse = .true.
  else
    writencoarse = .false.
  endif
  ! END CBH_LC

!  ! CBH_LC_mem
!  call title(nrestart,nchar)
!  filename_memstat = TRIM(filemem1)//TRIM(nchar)//TRIM(filememtimeext)
!  filename_memall = TRIM(filemem2)//TRIM(nchar)//TRIM(filememtimeext)
!  filename_timings = TRIM(filetime)//TRIM(nchar)//TRIM(filememtimeext)
!  if((displaymem.EQ.1).AND.(myid.EQ.1)) then
!    if (okmemstat) then
!      open(ilun_mem1,file=TRIM(filename_memstat),form='formatted',status='unknown')
!      rewind ilun_mem1
!      close(ilun_mem1)
!    endif
!    if (okmemall) then
!      open(ilun_mem2,file=TRIM(filename_memall),form='formatted',status='unknown')
!      rewind ilun_mem2
!      close(ilun_mem2)
!    endif
!  endif
!!  if(dotiming.AND.(myid.EQ.1)) open(ilun_time1,file=TRIM(filename_timings),form='formatted',status='unknown') ! CBH_LC - not using dotiming
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_amr ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_amr ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem

  call init_amr                      ! Initialize AMR variables

!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_time ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_time ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem

  call init_time                     ! Initialize time variables

!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_hydro ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_hydro ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       & call rt_init_hydro          ! Initialize radiation variables
#endif
!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_poisson ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_poisson ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem
  if(poisson)call init_poisson_gr    ! Initialize poisson and GR variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_refine ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_refine ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem
  if(nrestart==0)call init_refine    ! Build initial AMR grid

#ifdef grackle
  if(use_grackle==0)then
     if(cooling.and..not.neq_chem) &
        call set_table(dble(aexp))    ! Initialize cooling look up table
  endif
#else
!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_table ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_table ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_part ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_part ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem

  if(pic)call init_part              ! Initialize particle variables

!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_tree ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_tree ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem

  if(pic)call init_tree              ! Initialize particle tree
!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' Before init_refine2 ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' Before init_refine2 ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem

  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again
!  ! BEGIN CBH_LC_mem
!  if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' After init ', TRIM(filename_memstat), ilun_mem1, tt, tt)
!  if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' After init ', TRIM(filename_memall), ilun_mem2, tt, tt)
!  ! END CBH_LC_mem


#ifndef WITHOUTMPI
  muspt=0.
  tot_pt=-1
  tt2=MPI_WTIME()
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration'

  ! CBH_LC 09-02-2021
  !------------------ MODIF V. REVERDY 2011 ------------------!
  if(use_aexp_restart) nstep_coarse_after_restart=nstep_coarse_after_restart+1
  if(myid==1) then
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)'====================== BEGIN TIME STEP ======================'
  end if
  if(nstep_coarse_after_restart==1) aexp_restart_light_cone=aexp
  !-----------------------------------------------------------!    
  ! END CBH_LC 09-02-2021


  do ! Main time loop
                               call timer('coarse levels','start')

                               
#ifndef WITHOUTMPI
     tt1=MPI_WTIME()
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(momentum_feedback)call make_virtual_fine_dp(pstarold(1),ilevel)
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
#ifdef RT
           if(rt)then
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           endif
#endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if

           if(gr)then
              do igrp=1,10
                 call make_virtual_fine_dp(gr_pot(1,igrp),ilevel)
              end do
              do igrm=1,4
                 call make_virtual_fine_dp(gr_mat (1,igrm),ilevel)
                 call make_virtual_fine_dp(gr_mat2(1,igrm),ilevel)
              end do
           end if

           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)
                               call timer('coarse levels','start')

     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro)then
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(momentum_feedback)call make_virtual_fine_dp(pstarold(1),ilevel)
              if(simple_boundary)call make_boundary_hydro(ilevel)
           end if
#ifdef RT
           ! Radiation book-keeping
           if(rt)then
              call rt_upload_fine(ilevel)
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           end if
#endif
           ! Gravity book-keeping
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if

           if(gr)then
              do igrp=1,10
                 call make_virtual_fine_dp(gr_pot(1,igrp),ilevel)
              end do
              do igrm=1,4
                 call make_virtual_fine_dp(gr_mat (1,igrm),ilevel)
                 call make_virtual_fine_dp(gr_mat2(1,igrm),ilevel)
              end do
           end if
        end do

        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
!     nstep_coarse=nstep_coarse+1 CBH_LC - this is done directly in amr_step for lightcone

!! BEGIN CBH_LC_mem
!!------------------ MODIF V. REVERDY 2011 ------------------!
!    if((okmemstat).AND.(displaymem.EQ.1)) call writememfilestep(1, 0, 2, nstep, ' End of loop', TRIM(filename_memstat), ilun_mem1, tt, tt)
!    if((okmemall).AND.(displaymem.EQ.1)) call writememfilestep(3, 0, 2, nstep, ' End of loop ', TRIM(filename_memall), ilun_mem2, tt, tt)
!    ! if((okmemstat).AND.(displaymem.EQ.1)) call writememfile(1, 0, 2, ' With memory map ', 'memstat.txt', ilun_mem1, 0, 0)
!    ! if((okmemall).AND.(displaymem.EQ.1)) call writememfile(3, 0, 2, ' Memory map ', 'memmap.txt', 975, 0, 0)
!    ! if((displaymem.EQ.1)) call writememdir(nstep, 2, 0, ' Test', 1, 1, 1, 1)
!!-----------------------------------------------------------! 
!! END CBH_LC_mem
     
#ifndef WITHOUTMPI
     tt2=MPI_WTIME()
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot ,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!        ! BEGIN CBH_LC_mem
!        call MPI_ALLREDUCE(real_mem,real_mem_mean,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,info)
!        call MPI_ALLREDUCE(real_mem,real_mem_min ,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
!        real_mem_mean=real_mem_mean/ncpu
!        ! END CBH_LC_mem
        if(myid==1)then
           if (tot_pt==0) muspt=0. ! dont count first timestep
           n_step = int(numbtot(1,levelmin),kind=8)*twotondim
           do ilevel=levelmin+1,nlevelmax
             n_step = n_step + int(numbtot(1,ilevel),kind=8)*product(nsubcycle(levelmin:ilevel-1))*(twotondim-1)
           enddo
           muspt_this_step = (tt2-tt1)*1e6/n_step*ncpu
           muspt = muspt + muspt_this_step
           tot_pt = tot_pt + 1
           write(*,'(a,f8.2,a,f12.2,a,f12.2,a)')' Time elapsed since last coarse step:',tt2-tt1 &
          ,' s',muspt_this_step,' mus/pt'  &
          ,muspt / max(tot_pt,1), ' mus/pt (av)'
!           ! CBH_LC_mem
!           write(*,*)'Min'
!           call writemem(real_mem_min)
!           write(*,*)'Mean'
!           call writemem(real_mem_mean)
!           write(*,*)'Max'
!           ! END CBH_LC_mem
           call writemem(real_mem_tot)
           write(*,*)'Total running time:', NINT((tt2-tstart)*100.0)*0.01,'s'
        endif
        if(walltime_hrs.gt.0d0) then
           wallsec = walltime_hrs*3600.     ! Convert from hours to seconds
           dumpsec = minutes_dump*60.       ! Convert minutes before end to seconds
           if(wallsec-dumpsec.lt.tt2-tstart) then
              output_now=.true.
              if(myid==1) write(*,*) 'Dumping snapshot before walltime runs out'
              ! Now set walltime to a negative number so we don't keep printing outputs
              walltime_hrs = -1d0
           endif
        endif
     endif
#endif

!------------------ MODIF V. REVERDY 2011 ------------------!
  if(myid==1) then
    if(verbose)write(*,*)'======================= END TIME STEP ======================='
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
    if(verbose)write(*,*)' '
  end if
!-----------------------------------------------------------!   

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
