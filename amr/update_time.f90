!=======================================================================
real(kind=8) function wallclock()
  implicit none
#ifdef WITHOUTMPI
  integer,      save :: tstart
  integer            :: tcur
  integer            :: count_rate
#else
  real(kind=8), save :: tstart
  real(kind=8)       :: tcur
#endif
  logical,      save :: first_call=.true.
  real(kind=8), save :: norm, offset=0.
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !---------------------------------------------------------------------
  if (first_call) then
#ifdef WITHOUTMPI
     call system_clock(count=tstart, count_rate=count_rate)
     norm=1d0/count_rate
#else
     norm = 1d0
     tstart = MPI_Wtime()
#endif
     first_call=.false.
  end if
#ifdef WITHOUTMPI
  call system_clock(count=tcur)
#else
  tcur = MPI_Wtime()
#endif
  wallclock = (tcur-tstart)*norm + offset
  if (wallclock < 0.) then
     offset = offset + 24d0*3600d0
     wallclock = wallclock + 24d0*3600d0
  end if
end function wallclock
!=======================================================================
module timer_m
  implicit none
  integer,            parameter         :: mtimer=200                    ! max nr of timers
  real(kind=8),       dimension(mtimer) :: start, time
  integer                               :: ntimer=0, itimer
  character(len=72), dimension(mtimer)  :: labels
contains
!-----------------------------------------------------------------------
subroutine findit (label)
  implicit none
  character(len=*) label
  do itimer=1,ntimer
     if (trim(label) == trim(labels(itimer))) return
  end do
  ntimer = ntimer+1
  itimer = ntimer
  labels(itimer) = label
  time(itimer) = 0.
end subroutine
end module
!=======================================================================
subroutine timer (label, cmd)
  use timer_m
  implicit none
  character(len=*)::label,cmd
  real(kind=8)::wallclock,current
!-----------------------------------------------------------------------
  current = wallclock()                                                 ! current time
  if (itimer > 0) then                                                  ! if timer is active ..
     time(itimer) = time(itimer) + current - start(itimer)              ! add to it
  end if
  call findit (label)                                                   ! locate timer slot
  if (cmd == 'start') then                                              ! start command
     start(itimer) = current                                            ! register start time
  else if (cmd == 'stop') then                                          ! stop command
     itimer = 0                                                         ! turn off timer
  end if
end subroutine
!=======================================================================
subroutine output_timer(write_file, filename)
  use amr_parameters
  use amr_commons
  use timer_m
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  real(kind=8) :: gtotal, avtime, rmstime
  real(kind=8), dimension(ncpu) :: vtime
  integer,      dimension(ncpu) :: all_ntimer
  logical,      dimension(ncpu) :: gprint_timer
  integer      :: imn, imx, mpi_err, icpu
  logical      :: print_timer
#endif
  real(kind=8) :: total
  integer      :: i
  logical      :: id_is_one, write_file
  integer      :: ilun=11
  character(LEN=80)::filename, fileloc !Optional for writing timing info
!-----------------------------------------------------------------------
  id_is_one = myid == 1
  total = 1e-9
  if (.not. write_file) ilun=6 ! 6 = std output
  if (id_is_one .and. write_file) then
     fileloc=TRIM(filename) ! Open file for timing info
     open(unit=ilun,file=fileloc,form='formatted')
  endif

  if (id_is_one .and. ncpu==1) write (ilun,'(/a,i7,a)') '     seconds         %    STEP (rank=',myid,')'
  do i = 1,ntimer
     total = total + time(i)
  end do
  if (ncpu==1) then
     do i = 1,ntimer
        if (id_is_one .and. time(i)/total > 0.001) write (ilun,'(f12.3,4x,f6.1,4x,a24)') &
          time(i), 100.*time(i)/total,labels(i)
     end do
     if (id_is_one) write (ilun,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  end if
#ifndef WITHOUTMPI
  if (ncpu > 1) then
     ! Check that timers are consistent across ranks
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_GATHER(ntimer,1,MPI_INTEGER,all_ntimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
     if (id_is_one) then
        if (maxval(all_ntimer) .ne. minval(all_ntimer)) then
           write (ilun,*)
           write (ilun,*) '--------------------------------------------------------------------'
           write (ilun,*) 'Error: Inconsistent number of timers on each rank. Min, max nr:', minval(all_ntimer), maxval(all_ntimer)
           write (ilun,*) 'Timing summary below can be misleading'
           write (ilun,*) 'Labels of timer on rank==1 :'
           write (ilun,*) '--------------------------------------------------------------------'
           do i=1,ntimer
              write(ilun,'(i3,1x,a)') i, labels(i)
           enddo
        endif
        ! Find first occurence of a rank with a different number of timers -- if it exists
        gprint_timer=.false.
        do icpu=1,ncpu
           if (all_ntimer(icpu) .ne. ntimer) then
              gprint_timer(icpu) = .true.
              exit
           endif
        enddo
        if (any(gprint_timer)) call sleep(1) ! Make sure that master rank finished, before we print from other rank.
     endif
     call MPI_SCATTER(gprint_timer,1,MPI_LOGICAL,print_timer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_err)
     if (print_timer) then
        write (ilun,*)
        write (ilun,*) 'Labels of timer on rank==',myid
        write (ilun,*) '--------------------------------------------------------------------'
        do i=1,ntimer
           write(ilun,'(i3,1x,a)') i, labels(i)
        enddo
        write (ilun,*)
     endif

     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(total,gtotal,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
     gtotal = gtotal / ncpu

     if (id_is_one) write (ilun,*) '--------------------------------------------------------------------'
     if (id_is_one) write (ilun,'(/a)') '     minimum       average       maximum' // &
                  '  standard dev        std/av       %   rmn   rmx  TIMER'
     do i = 1,ntimer
        call MPI_GATHER(real(time(i),kind=8),1,MPI_REAL8,vtime,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
        if (id_is_one) then
           if (maxval(vtime)/gtotal > 0.001) then
              avtime  = sum(vtime) / ncpu ! average time used
              imn     = minloc(vtime,1)
              imx     = maxloc(vtime,1)
              rmstime = sqrt(sum((vtime - avtime)**2)/ncpu)
              write (ilun,'(5(f12.3,2x),f6.1,2x,2i4,4x,a24)') &
                 vtime(imn), avtime, vtime(imx), rmstime, rmstime/avtime, 100.*avtime/gtotal, imn, imx, labels(i)
           endif
        endif
     end do
     if (id_is_one) write (ilun,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  endif
#endif
  if (id_is_one) close(ilun)
end subroutine
!=======================================================================
subroutine reset_timer
   use timer_m
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
!-----------------------------------------------------------------------
   do itimer = 1,ntimer
      time(itimer)=0.0
   end do
end subroutine
!=======================================================================
subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  real(kind=8)::ttend
  real(kind=8),save::ttstart=0
#endif
  integer::ilevel

  real(dp)::dt,econs,mcons
  integer::i,itest

  ! Local constants
  dt=dtnew(ilevel)
  itest=0

#ifndef WITHOUTMPI
  if(myid==1)then
     if(ttstart.eq.0.0)ttstart=MPI_WTIME()
  endif
#endif

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(nstep_coarse .ne. nstep_coarse_old)then

     !--------------------------
     ! Check mass conservation
     !--------------------------
     if(mass_tot_0==0.0D0)then
        mass_tot_0=mass_tot
        mcons=0.0D0
     else
        mcons=(mass_tot-mass_tot_0)/mass_tot_0
     end if

     !----------------------------
     ! Check energy conservation
     !----------------------------
     if(epot_tot_old.ne.0)then
        epot_tot_int=epot_tot_int + &
             & 0.5D0*(epot_tot_old+epot_tot)*log(aexp/aexp_old)
     end if
     epot_tot_old=epot_tot
     aexp_old=aexp
     if(einit==0.0D0)then
        einit=epot_tot+ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(ekin_tot+epot_tot-epot_tot_int-einit) / &
             &(-(epot_tot-epot_tot_int-einit)+ekin_tot)
     end if

     if(mod(nstep_coarse,ncontrol)==0.or.output_done)then
        if(myid==1)then

           !-------------------------------
           ! Output AMR structure to screen
           !-------------------------------
           write(*,*)'Mesh structure'
           do i=1,nlevelmax
              if(numbtot(1,i)>0)write(*,999)i,numbtot(1:4,i)
           end do

           !----------------------------------------------
           ! Output mass and energy conservation to screen
           !----------------------------------------------
           if(cooling.or.pressure_fix)then
              write(*,778)nstep_coarse,mcons,econs,epot_tot,ekin_tot,eint_tot
           else
              write(*,777)nstep_coarse,mcons,econs,epot_tot,ekin_tot
           end if
#ifdef SOLVERmhd
           write(*,'(" emag=",ES9.2)') emag_tot
#endif
           if(pic)then
              write(*,888)nstep,t,dt,aexp,&
                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
                   & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
           else
              write(*,888)nstep,t,dt,aexp,&
                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
           endif
           itest=1
        end if
        output_done=.false.
     end if

     !---------------
     ! Exit program
     !---------------
     if(t>=tout(noutput).or.aexp>=aout(noutput).or. &
          & nstep_coarse>=nstepmax)then
        if(myid==1)then
           write(*,*)'Run completed'
#ifndef WITHOUTMPI
           ttend=MPI_WTIME()
           write(*,*)'Total elapsed time:',ttend-ttstart
#endif
        endif
        call clean_stop
     end if

  end if
  nstep_coarse_old=nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(nstep,ncontrol)==0)then
     if(myid==1.and.itest==0)then
        if(pic)then
           write(*,888)nstep,t,dt,aexp,&
                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
                & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
        else
           write(*,888)nstep,t,dt,aexp,&
                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
        endif
     end if
  end if

  !------------------------
  ! Update time variables
  !------------------------
  t=t+dt
  nstep=nstep+1
  if(cosmo)then

     ! Find neighboring times
     i=1
     do while(tau_frw(i)>t.and.i<n_frw)
        i=i+1
     end do
     t=max(t,tau_frw(n_frw)) ! CBH_LC
     ! Interpolate expansion factor
     aexp = aexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & aexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     hexp = hexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & hexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     texp =    t_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          &    t_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
  else
     aexp = 1.0
     hexp = 0.0
     texp = t
  end if

777 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2)
778 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
888 format(' Fine step=',i7,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine update_time

subroutine clean_stop
  use amr_commons
  use poisson_commons
  use gr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer :: ilevel
  character(LEN=80)::str

  call output_timer(.false., str)

#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif

  ! allocations in read_params.f90
  if(allocated(remap_pscalar)) deallocate(remap_pscalar)


  ! allocations in init_amr.f90
  if(allocated(bound_key)) deallocate(bound_key)
  if(allocated(bound_key2)) deallocate(bound_key2)
  if(allocated(headl)) deallocate(headl)
  if(allocated(taill)) deallocate(taill)
  if(allocated(numbl)) deallocate(numbl)
  if(allocated(numbtot)) deallocate(numbtot)
  if(allocated(headb)) deallocate(headb)
  if(allocated(tailb)) deallocate(tailb)
  if(allocated(numbb)) deallocate(numbb)
  if(allocated(boundary)) deallocate(boundary)

  ! communicators
  if(allocated(active))then
     do ilevel=1,nlevelmax
        ! virtual_boundaries.f90
        ! TODO: cleaner solution (fortran 2003): s/pointer/allocatable/
        if(active(ilevel)%ngrid>0) deallocate(active(ilevel)%igrid)
     enddo
     deallocate(active)
  endif
  if(allocated(emission)) deallocate(emission)
  if(allocated(reception)) deallocate(reception)
  !
  if(allocated(lookup_mg)) deallocate(lookup_mg)
  !
  if(allocated(father)) deallocate(father)
  if(allocated(nbor)) deallocate(nbor)
  if(allocated(next)) deallocate(next)
  if(allocated(prev)) deallocate(prev)
  if(pic)then
     if(allocated(headp)) deallocate(headp)
     if(allocated(tailp)) deallocate(tailp)
     if(allocated(numbp)) deallocate(numbp)
  endif
  if(allocated(xg)) deallocate(xg)
  ! amr cell-based arrays
  if(allocated(flag1)) deallocate(flag1)
  if(allocated(flag2)) deallocate(flag2)
  if(allocated(son)) deallocate(son)
  ! mpi cell-based arrays
  if(allocated(cpu_map)) deallocate(cpu_map)
  if(allocated(cpu_map2)) deallocate(cpu_map2)
  if(allocated(hilbert_key)) deallocate(hilbert_key)


  ! allocations in init_poisson.f90
  if(allocated(safe_mode)) deallocate(safe_mode)
  if(allocated(active_mg)) deallocate(active_mg)
  if(allocated(emission_mg)) deallocate(emission_mg)
  ! cell-centred variables
  if(allocated(rho)) deallocate(rho)
  if(allocated(phi)) deallocate(phi)
  
  if(gr)then
     if(allocated(gr_pot)) deallocate(gr_pot)
     if(allocated(gr_mat)) deallocate(gr_mat)
     if(allocated(gr_mat2)) deallocate(gr_mat2)
  end if

  if(allocated(phi_old)) deallocate(phi_old)
  if(allocated(f)) deallocate(f)


  ! allocations in init_time.f90
  if(allocated(aexp_frw)) deallocate(aexp_frw)
  if(allocated(hexp_frw)) deallocate(hexp_frw)
  if(allocated(tau_frw)) deallocate(tau_frw)
  if(allocated(t_frw)) deallocate(t_frw)


  ! init_part.f90 - in general, BIG deallocations
  if(allocated(idp)) deallocate(idp)
  if(allocated(nextp)) deallocate(nextp)
  if(allocated(prevp)) deallocate(prevp)
  if(allocated(levelp)) deallocate(levelp)
  if(allocated(mp)) deallocate(mp)
  if(allocated(vp)) deallocate(vp)
  if(allocated(xp)) deallocate(xp)

  !if(gr) then
  !   if(allocated(fp_gr)) deallocate(fp_gr)
  ! end if

  stop
end subroutine clean_stop

!-------------------------------------------------!
!------------------- CBH_LC ----------------------!
!-------------------------------------------------!

!------------------ MODIF V. REVERDY 2011 ------------------!
subroutine writememdir(memnstep, memformat, openmode, comment, writecomment, writetitles, writememstat, writememmap)

  ! Options
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
#ifndef WITHOUTMPI
  use mpi
#endif 
  implicit none

  ! Input variables
  integer::memnstep                ! ncoarse step
  integer::memformat               ! -1 = auto, 0 = o, 1 = Ko, 2 = Mo, 3 = Go
  integer::openmode                ! 0 = append, 1 = old, 2 = new
  character(LEN=*)::comment        ! title of the current written part
  integer::writecomment            ! 1 = write the comment, 0 = dont write the comment
  integer::writetitles             ! 1 = write min, mean, max, tot, 0 = dont write
  integer::writememstat            ! 1 = write mem stat, 0 = dont write
  integer::writememmap             ! 1 = write mem map, 0 = dont write

  ! Other variables
  real::real_mem = 0.                               ! used memory for the current cpu
  real(KIND=4)::real_mem_4 = 0.                     ! used memory for the current cpu in kind = 4
  real::real_mem_min = 0.                           ! min of mem
  real::real_mem_mean = 0.                          ! mean of mem
  real::real_mem_max = 0.                           ! max of mem
  real::real_mem_tot = 0.                           ! total of memory
  real::formattedmem                                ! function call
  integer(KIND=4),allocatable::real_mem_table_cpu(:)! table for cpu indexes
  real(KIND=4),allocatable::real_mem_table_mem(:)   ! table for mem
  integer::info                                     ! info MPI
  integer::i                                        ! loop variable

  ! Strings
  integer::out_unit_stat = ilun_mem1                ! unit for min, mean, max
  integer::out_unit_map = ilun_mem2                 ! unit for map
  character(LEN=128)::filename_dir                  ! directory
  character(LEN=128)::filename_stat                 ! stat files
  character(LEN=128)::filename_map                  ! map file 
  character(LEN=5)::filename_nchar                  ! title nchar

#ifndef WITHOUTMPI

  ! Initialize
  allocate(real_mem_table_cpu(1:ncpu))
  allocate(real_mem_table_mem(1:ncpu))

  ! Get memory of current cpu
  call getmem(real_mem)
  real_mem_4 = real(real_mem, KIND=4)

  ! Compute min, mean, max and map
  call MPI_ALLREDUCE(real_mem,real_mem_max ,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(real_mem,real_mem_min ,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLGATHER(myid,1,MPI_REAL,real_mem_table_cpu,1,MPI_INTEGER,MPI_COMM_WORLD,info)
  call MPI_ALLGATHER(real_mem_4,1,MPI_REAL,real_mem_table_mem,1,MPI_REAL,MPI_COMM_WORLD,info)
  real_mem_mean=real_mem_tot/ncpu

  ! Write to files
  if (myid.EQ.1) then

     ! File names
     call title(memnstep, filename_nchar)
     filename_dir = TRIM('')
     filename_stat = TRIM(filename_dir)//TRIM('memstat_')//filename_nchar
     filename_map = TRIM(filename_dir)//TRIM('memmap_')//filename_nchar
     filename_stat = TRIM(filename_stat)
     filename_map = TRIM (filename_map)

     ! Mem stats
     if (writememstat .NE. 0) then
      !write(*,*) 'TEST_STAT = ', filename_stat
        open(unit=out_unit_stat, file=filename_stat, action="write", status="unknown", access="append")
           if(writecomment .NE. 0) write(out_unit_stat,*)'# ',trim(comment)
           if(writetitles .NE. 0) then
              write(out_unit_stat,*)'Min = ', formattedmem(real_mem_min, memformat)
              write(out_unit_stat,*)'Mean = ', formattedmem(real_mem_mean, memformat)
              write(out_unit_stat,*)'Max = ', formattedmem(real_mem_max, memformat)
              write(out_unit_stat,*)'Tot = ', formattedmem(real_mem_tot, memformat)
           else
              write(out_unit_stat,*)formattedmem(real_mem_min, memformat)
              write(out_unit_stat,*)formattedmem(real_mem_mean, memformat)
              write(out_unit_stat,*)formattedmem(real_mem_max, memformat)
              write(out_unit_stat,*)formattedmem(real_mem_tot, memformat)
           endif
           write(out_unit_stat,*)' '
        close(out_unit_stat)
     endif

     ! Mem map
     if (writememmap .NE. 0) then
        !write(*,*) 'TEST_MAP = ', filename_map
        open(unit=out_unit_map, file=filename_map, action="write", status="unknown", access="append")
           if(writecomment .NE. 0) write(out_unit_map,*)'# ',trim(comment)
           if(writetitles .NE. 0) then
              write(out_unit_map,*)'Mem map for ncpu = ', ncpu
           else
              write(out_unit_map,*)ncpu
           endif
           i = 0
           do while(i<=ncpu)
              write(out_unit_map,*)real_mem_table_cpu(i),' ',formattedmem(real_mem_table_mem(i),memformat)
              i = i+1
           end do
           write(out_unit_map,*)' '
        close(out_unit_map)
     endif

  ! End of writing
  endif
  
  ! Dealloc
  deallocate(real_mem_table_cpu)
  deallocate(real_mem_table_mem)
 
#endif
end subroutine writememdir

subroutine writememfilestep(writememmode, linuxmode, memformat, mem_nstep, comment, filename, fileunit, tt1, tt2)

  ! Options
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
#ifndef WITHOUTMPI
  use mpi
#endif 
  implicit none

  ! Input variables
  integer::writememmode            ! 0 = do nothing, 1 = print min, mean, max, 2 = print stat for each cpu, 3 = print all
  integer::linuxmode               ! 0 = first method (should work on most linux system), 1 = second method for titane
  integer::memformat               ! -1 = auto, 0 = o, 1 = Ko, 2 = Mo, 3 = Go
  character(LEN=*)::filename       ! filename
  integer::fileunit                ! unit
  character(LEN=*)::comment        ! title of the current written part
  real::tt1                        ! previous time
  real::tt2                        ! current time (if tt2 < 0 and tt1 < 0 no time is written in the file)
  integer::mem_nstep               ! simulation step

  ! Other variables
  real::real_mem = 0.                               ! used memory for the current cpu
  real(KIND=4)::real_mem_4 = 0.                     ! used memory for the current cpu in kind = 4
  real::real_mem_min = 0.                           ! min of mem
  real::real_mem_mean = 0.                          ! mean of mem
  real::real_mem_max = 0.                           ! max of mem
  real::real_mem_tot = 0.                           ! total of memory
  real::formattedmem                                ! function call
  integer::int_mem_tmp = 0                          ! used memory for linuxmode
  integer::out_unit = 891                           ! unit
  integer(KIND=4),allocatable::real_mem_table_cpu(:)! table for cpu indexes
  real(KIND=4),allocatable::real_mem_table_mem(:)   ! table for mem
  integer::info                                     ! info MPI
  integer::i                                        ! loop variable
  character(LEN=5)::mem_nchar                       ! step characters
  integer::mem_icpumax                              ! icpu of maximum mem
  integer::mem_icpumin                              ! icup of minimum mem
  integer::reduce_mem_icpumax                       ! icpu of maximum mem
  integer::reduce_mem_icpumin                       ! icup of minimum mem
  mem_icpumax = -1
  mem_icpumin = -1
#ifndef WITHOUTMPI
  ! First method (should work on most linux system)
  if (((linuxmode .EQ. 0).OR.(linuxmode .EQ. 1)).AND.((writememmode.EQ.1).OR.(writememmode.EQ.2).OR.(writememmode.EQ.3))) then

     ! Initialization
     out_unit = fileunit
     allocate(real_mem_table_cpu(1:ncpu))
     allocate(real_mem_table_mem(1:ncpu))

     ! Get mem
     if (linuxmode .EQ. 0) then
        call getmem(real_mem)         ! general case
     else if (linuxmode .EQ. 1) then
        !call cccmcur(int_mem_tmp)    ! uncomment for Titane
        real_mem=1024.*dble(int_mem_tmp)
     endif

     real_mem_4 = real(real_mem, KIND=4)

     ! Nothing
     if (writememmode .EQ. 0) then
        real_mem_min = real_mem
        real_mem_mean = real_mem
        real_mem_max = real_mem
        real_mem_tot = real_mem

        ! Min/Mean/Max
     else if (writememmode .EQ. 1) then
        call MPI_ALLREDUCE(real_mem,real_mem_max ,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_min ,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
        real_mem_mean=real_mem_tot/ncpu
        if (real_mem.GE.real_mem_max) mem_icpumax = myid-1
        if (real_mem.LE.real_mem_min) mem_icpumin = myid-1
        call MPI_REDUCE(mem_icpumax,reduce_mem_icpumax,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,info)
        call MPI_REDUCE(mem_icpumin,reduce_mem_icpumin,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,info)

        if (myid==1) then
           call title(mem_nstep, mem_nchar)
           open(unit=out_unit,file=filename,action="write",status="unknown",access="append")
           write(out_unit,*)'# ',mem_nchar,' : ',trim(comment)
           if ((tt1.GT.0).OR.(tt2.GT.0)) write(out_unit,*)'Coarse delta_t =', tt2-tt1
           write(out_unit,*)'Min = ', formattedmem(real_mem_min, memformat), reduce_mem_icpumin
           write(out_unit,*)'Mean = ', formattedmem(real_mem_mean, memformat)
           write(out_unit,*)'Max = ', formattedmem(real_mem_max, memformat), reduce_mem_icpumax
           write(out_unit,*)'Tot = ', formattedmem(real_mem_tot, memformat)
           write(out_unit,*)' '
           close(out_unit)
        endif

        ! Stats for each CPU
     else if (writememmode .EQ. 2) then
        call MPI_ALLGATHER(myid,1,MPI_REAL,real_mem_table_cpu,1,MPI_INTEGER,MPI_COMM_WORLD,info)
        call MPI_ALLGATHER(real_mem_4,1,MPI_REAL,real_mem_table_mem,1,MPI_REAL,MPI_COMM_WORLD,info)

        if (myid .EQ. 1) then
           call title(mem_nstep, mem_nchar)
           open(unit=out_unit,file=filename,action="write",status="unknown",access="append")
           write(out_unit,*)'# ',mem_nchar,' : ',trim(comment)
           if ((tt1.GT.0).OR.(tt2.GT.0)) write(out_unit,*)'Coarse delta_t =', tt2-tt1
           write(out_unit,*)'Mem map for ncpu = ',ncpu
           i = 1
           do while(i<=ncpu)
              write(out_unit,*)real_mem_table_cpu(i),' ',formattedmem(real_mem_table_mem(i),memformat)
              i = i+1
           end do
           write(out_unit,*)' '
           close(out_unit)
        endif

        ! Min/Mean/Max + Stats for each CPU
     else if (writememmode .EQ. 3) then
        call MPI_ALLREDUCE(real_mem,real_mem_max ,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_min ,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
        call MPI_ALLGATHER(myid,1,MPI_REAL,real_mem_table_cpu,1,MPI_INTEGER,MPI_COMM_WORLD,info)
        call MPI_ALLGATHER(real_mem_4,1,MPI_REAL,real_mem_table_mem,1,MPI_REAL,MPI_COMM_WORLD,info)
        real_mem_mean=real_mem_tot/ncpu
        if (real_mem.GE.real_mem_max) mem_icpumax = myid-1
        if (real_mem.LE.real_mem_min) mem_icpumin = myid-1
        call MPI_REDUCE(mem_icpumax,reduce_mem_icpumax,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,info)
        call MPI_REDUCE(mem_icpumin,reduce_mem_icpumin,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,info)

        if (myid==1) then
           call title(mem_nstep, mem_nchar)
           open(unit=out_unit,file=filename,action="write",status="unknown",access="append")
           write(out_unit,*)'# ',mem_nchar,' : ',trim(comment)
           if ((tt1.GT.0).OR.(tt2.GT.0)) write(out_unit,*)'Coarse delta_t =', tt2-tt1
           write(out_unit,*)'Min = ', formattedmem(real_mem_min, memformat), reduce_mem_icpumin
           write(out_unit,*)'Mean = ', formattedmem(real_mem_mean, memformat)
           write(out_unit,*)'Max = ', formattedmem(real_mem_max, memformat), reduce_mem_icpumax
           write(out_unit,*)'Tot = ', formattedmem(real_mem_tot, memformat)
           write(out_unit,*)'Mem map for ncpu = ',ncpu
           i = 1
           do while(i<=ncpu)
              write(out_unit,*)real_mem_table_cpu(i),' ',formattedmem(real_mem_table_mem(i),memformat)
              i = i+1
           end do
           write(out_unit,*)' '
           close(out_unit)
        endif

     endif

  ! Dealloc
  deallocate(real_mem_table_cpu)
  deallocate(real_mem_table_mem)

  endif

#endif
end subroutine writememfilestep

subroutine writememfile(writememmode, linuxmode, memformat, comment, filename, fileunit, tt1, tt2)

  ! Options
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
#ifndef WITHOUTMPI
  use mpi
#endif 
  implicit none

  ! Input variables
  integer::writememmode            ! 0 = do nothing, 1 = print min, mean, max, 2 = print stat for each cpu, 3 = print all
  integer::linuxmode               ! 0 = first method (should work on most linux system), 1 = second method for titane
  integer::memformat               ! -1 = auto, 0 = o, 1 = Ko, 2 = Mo, 3 = Go
  character(LEN=*)::filename       ! filename
  integer::fileunit                ! unit
  character(LEN=*)::comment        ! title of the current written part
  real::tt1                        ! previous time
  real::tt2                        ! current time (if tt2 < 0 and tt1 < 0 no time is written in the file)

  ! Other variables
  real::real_mem = 0.                               ! used memory for the current cpu
  real(KIND=4)::real_mem_4 = 0.                     ! used memory for the current cpu in kind = 4
  real::real_mem_min = 0.                           ! min of mem
  real::real_mem_mean = 0.                          ! mean of mem
  real::real_mem_max = 0.                           ! max of mem
  real::real_mem_tot = 0.                           ! total of memory
  real::formattedmem                                ! function call
  integer::int_mem_tmp = 0                          ! used memory for linuxmode
  integer::out_unit = 891                           ! unit
  integer(KIND=4),allocatable::real_mem_table_cpu(:)! table for cpu indexes
  real(KIND=4),allocatable::real_mem_table_mem(:)   ! table for mem
  integer::info                                     ! info MPI
  integer::i                                        ! loop variable

#ifndef WITHOUTMPI
  ! First method (should work on most linux system)
  if (((linuxmode .EQ. 0).OR.(linuxmode .EQ. 1)).AND.((writememmode.EQ.1).OR.(writememmode.EQ.2).OR.(writememmode.EQ.3))) then

     ! Initialization
     out_unit = fileunit
     allocate(real_mem_table_cpu(1:ncpu))
     allocate(real_mem_table_mem(1:ncpu))

     ! Get mem
     if (linuxmode .EQ. 0) then
        call getmem(real_mem)         ! general case
     else if (linuxmode .EQ. 1) then
        !call cccmcur(int_mem_tmp)    ! uncomment for Titane
        real_mem=1024.*dble(int_mem_tmp)
     endif

     real_mem_4 = real(real_mem, KIND=4)

     ! Nothing
     if (writememmode .EQ. 0) then
        real_mem_min = real_mem
        real_mem_mean = real_mem
        real_mem_max = real_mem
        real_mem_tot = real_mem

        ! Min/Mean/Max
     else if (writememmode .EQ. 1) then
        call MPI_ALLREDUCE(real_mem,real_mem_max ,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_min ,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
        real_mem_mean=real_mem_tot/ncpu

        if (myid==1) then
           open(unit=out_unit,file=filename,action="write",status="unknown",access="append")
           write(out_unit,*)'#',trim(comment)
           if ((tt1.GT.0).OR.(tt2.GT.0)) write(out_unit,*)'Coarse delta_t =', tt2-tt1
           write(out_unit,*)'Min = ', formattedmem(real_mem_min, memformat)
           write(out_unit,*)'Mean = ', formattedmem(real_mem_mean, memformat)
           write(out_unit,*)'Max = ', formattedmem(real_mem_max, memformat)
           write(out_unit,*)'Tot = ', formattedmem(real_mem_tot, memformat)
           write(out_unit,*)' '
           close(out_unit)
        endif

        ! Stats for each CPU
     else if (writememmode .EQ. 2) then
        call MPI_ALLGATHER(myid,1,MPI_REAL,real_mem_table_cpu,1,MPI_INTEGER,MPI_COMM_WORLD,info)
        call MPI_ALLGATHER(real_mem_4,1,MPI_REAL,real_mem_table_mem,1,MPI_REAL,MPI_COMM_WORLD,info)

        if (myid .EQ. 1) then
           open(unit=out_unit,file=filename,action="write",status="unknown",access="append")
           write(out_unit,*)'#',trim(comment)
           if ((tt1.GT.0).OR.(tt2.GT.0)) write(out_unit,*)'Coarse delta_t =', tt2-tt1
           write(out_unit,*)'Mem map for ncpu = ',ncpu
           i = 1
           do while(i<=ncpu)
              write(out_unit,*)real_mem_table_cpu(i),' ',formattedmem(real_mem_table_mem(i),memformat)
              i = i+1
           end do
           write(out_unit,*)' '
           close(out_unit)
        endif

        ! Min/Mean/Max + Stats for each CPU
     else if (writememmode .EQ. 3) then
        call MPI_ALLREDUCE(real_mem,real_mem_max ,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_min ,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
        call MPI_ALLGATHER(myid,1,MPI_REAL,real_mem_table_cpu,1,MPI_INTEGER,MPI_COMM_WORLD,info)
        call MPI_ALLGATHER(real_mem_4,1,MPI_REAL,real_mem_table_mem,1,MPI_REAL,MPI_COMM_WORLD,info)
        real_mem_mean=real_mem_tot/ncpu

        if (myid==1) then
           open(unit=out_unit,file=filename,action="write",status="unknown",access="append")
           write(out_unit,*)'# ',trim(comment)
           if ((tt1.GT.0).OR.(tt2.GT.0)) write(out_unit,*)'Coarse delta_t =', tt2-tt1
           write(out_unit,*)'Min = ', formattedmem(real_mem_min, memformat)
           write(out_unit,*)'Mean = ', formattedmem(real_mem_mean, memformat)
           write(out_unit,*)'Max = ', formattedmem(real_mem_max, memformat)
           write(out_unit,*)'Tot = ', formattedmem(real_mem_tot, memformat)
           write(out_unit,*)'Mem map for ncpu = ',ncpu
           i = 1
           do while(i<=ncpu)
              write(out_unit,*)real_mem_table_cpu(i),' ',formattedmem(real_mem_table_mem(i),memformat)
              i = i+1
           end do
           write(out_unit,*)' '
           close(out_unit)
        endif

     endif

  ! Dealloc
  deallocate(real_mem_table_cpu)
  deallocate(real_mem_table_mem)

  endif

#endif
end subroutine writememfile


! Returns the used memory in a given format (octet, Ko, Mo, Go)
real function formattedmem(usedmem, memformat)
  implicit none
  ! Input variables
  real::usedmem
  integer::memformat

  ! Local variables
  integer::getpagesize  = 0
  integer::ipagesize = 0
  !real::formattedmem

  ! Page size and octet memory
#ifdef NOSYSTEM
  !  call PXFSYSCONF(_SC_PAGESIZE,ipagesize,ierror)
  ipagesize=4096
#else
  !  ipagesize = getpagesize()
  ipagesize=4096
#endif 

  usedmem=dble(usedmem)*dble(ipagesize)

  ! Compute result
  if (memformat .EQ. 0) then
     formattedmem = usedmem
  else if (memformat .EQ. 1) then
     formattedmem = usedmem/1024
  else if (memformat .EQ. 2) then
     formattedmem = usedmem/1024.**2.
  else if (memformat .EQ. 3) then
     formattedmem = usedmem/1024.**3.
  else
     if(usedmem>1024.**3.)then
        formattedmem = usedmem/1024.**3.
     else if (usedmem>1024.**2.) then
        formattedmem = usedmem/1024.**2
     else if (usedmem>1024.) then
        formattedmem = usedmem/1024.
     endif
  endif

  return
end function formattedmem
!------------------ MODIF V. REVERDY 2011 ------------------!

!-------------------------------------------------!
!------------------ END CBH_LC  ------------------!
!-------------------------------------------------!


subroutine writemem(usedmem)
  real(kind=4)::usedmem

  usedmem=real(usedmem)*4096

  if(usedmem>1024.**4.)then
     write(*,999)usedmem/1024.**4.
  else if (usedmem>1024.**3.) then
     write(*,998)usedmem/1024.**3.
  else if (usedmem>1024.**2.) then
     write(*,997)usedmem/1024.**2.
  else if (usedmem>1024.) then
     write(*,996)usedmem/1024.
  endif

996 format(' Used memory:',F9.1,' kB')
997 format(' Used memory:',F9.1,' MB')
998 format(' Used memory:',F9.3,' GB')
999 format(' Used memory:',F9.3,' TB')

end subroutine writemem

subroutine getmem(outmem)
  use amr_commons,only:myid
#ifndef WITHOUTMPI
  use amr_commons,only:IOGROUPSIZE
  use amr_commons,only:ncpu
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::dummy_io,info2
#endif
  real(kind=4)::outmem
  character(len=300) :: dir, dir2, file
  integer::read_status
  integer,parameter::tag=1134
  integer::nmem,ind,j
  logical::file_exists

  file='/proc/self/stat'
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  inquire(file=file, exist=file_exists)
  if (file_exists) then
     open(unit=12,file=file,form='formatted')
     read(12,'(A300)',IOSTAT=read_status)dir
     close(12)
  else
     read_status=-1000
  endif

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


  if (read_status < 0)then
     outmem=0.
     if (myid==1 .and. read_status .ne. -1000)write(*,*)'Problem in checking free memory'
  else
     ind=300
     j=0
     do while (j<23)
        ind=index(dir,' ')
        dir2=dir(ind+1:300)
        j=j+1
        dir=dir2
     end do
     ind=index(dir,' ')
     dir2=dir(1:ind)
     read(dir2,'(I12)')nmem
     outmem=real(nmem,kind=4)
  end if

end subroutine getmem
!------------------------------------------------------------------------
SUBROUTINE getProperTime(tau,tproper)
! Calculate proper time tproper corresponding to conformal time tau (both
! in code units).
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::tau, tproper
  integer::i
  if(.not. cosmo .or. tau .eq. 0.d0) then ! this might happen quite often
     tproper = tau
     return
  endif
  i = 1
  do while( tau_frw(i) > tau .and. i < n_frw )
     i = i+1
  end do
  tproper = t_frw(i  )*(tau-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & t_frw(i-1)*(tau-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
END SUBROUTINE getProperTime
!------------------------------------------------------------------------
SUBROUTINE getAgeGyr(t_birth_proper, age)
! Calculate proper time passed, in Gyrs, since proper time t_birth_proper
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  implicit none
  real(dp):: t_birth_proper, age
  real(dp), parameter:: yr = 3.15569d+07
  real(dp),save:: scale_t_Gyr
  logical,save::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then
     ! The timescale has not been initialized
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_t_Gyr = (scale_t/aexp**2)/yr/1.e9
     scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_Gyr
END SUBROUTINE getAgeGyr
!------------------------------------------------------------------------
SUBROUTINE getAgeSec(t_birth_proper, age)
! Calculate proper time passed, in sec, since proper time t_birth_proper
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  implicit none
  real(dp):: t_birth_proper, age
  real(dp),save:: scale_t_sec
  logical::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then
     ! The timescale has not been initialized
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_t_sec = (scale_t/aexp**2)
     scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_sec
END SUBROUTINE getAgeSec
!------------------------------------------------------------------------
