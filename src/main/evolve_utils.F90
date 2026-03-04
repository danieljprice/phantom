!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module evolve_utils
!
! Utility routines for the evolve module
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: analysis, checkconserved, dim, dynamic_dtmax, energies,
!   evwrite, externalforces, fileutils, forcing, injection, io, io_control,
!   io_summary, mf_write, mpiutils, options, part, ptmass, readwrite_dumps,
!   readwrite_infile, subgroup, substepping, timestep, timestep_ind,
!   timing, utils_apr
!
 implicit none
 public :: update_time_and_dt
 public :: ptmass_create_and_update_forces
 public :: write_ev_files
 public :: check_and_write_dump
 public :: print_log

 private

contains

!----------------------------------------------------------------
!+
!  wrapper routine to update the time and get the dt for the
!  next timestep. Also print a log every step as required.
!+
!----------------------------------------------------------------
subroutine update_time_and_dt(nsteps,time,dtmax,dtmaxold,rhomaxnow,tlast,tcheck,tprint,dt,tall,tstep,tcpustep,&
                              istepfrac,nbinmaxprev,ntot,nalivetot,nmovedtot,at_dump_time)
 use dim,           only:ind_timesteps
 use io,            only:id,master,nprocs,iverbose,iprint,warning,fatal
 use io_control,    only:at_simulation_end
 use dynamic_dtmax, only:dtmax_ifactor
 use mpiutils,      only:bcast_mpi,reduceall_mpi
 use timestep,      only:dtrad,dtforce,dtinject,dtcourant,dterr,print_dtlog,tmax
 use timestep_ind,  only:print_dtind_efficiency,update_time_per_bin,print_dtlog_ind,change_nbinmax,&
                         nactivetot,nbinmax
 integer,         intent(inout) :: nsteps
 real,            intent(inout) :: time,tcheck,tprint,dt
 real,            intent(in)    :: dtmax,dtmaxold,rhomaxnow,tlast
 real(kind=4),    intent(inout) :: tall
 real(kind=4),    intent(in)    :: tstep,tcpustep
 integer,         intent(inout) :: istepfrac
 integer(kind=1), intent(in)    :: nbinmaxprev
 integer(kind=8), intent(in)    :: ntot,nalivetot
 integer(kind=8), intent(inout) :: nmovedtot
 logical,         intent(out)   :: at_dump_time
 real :: dtprint,timecheck
 integer :: inbin

 nsteps = nsteps + 1

 if (ind_timesteps) then
    tcheck = tcheck + dt

    !--update time in way that is free of round-off errors
    time = tlast + istepfrac/real(2**nbinmaxprev)*dtmaxold

    !--print efficiency of partial timestep
    if (id==master) call print_dtind_efficiency(iverbose,nalivetot,nactivetot,tall,tstep,1)

    call update_time_per_bin(tcpustep,istepfrac,nbinmaxprev,inbin)
    nmovedtot = nmovedtot + nactivetot

    !--check that time is as it should be, may indicate error in individual timestep routines
    if (abs(tcheck-time) > 1.e-4) call warning('evolve','time out of sync',var='error',val=abs(tcheck-time))

    if (id==master .and. (iverbose >= 1 .or. inbin <= 3)) &
       call print_dtlog_ind(iprint,istepfrac,2**nbinmaxprev,time,dt,nactivetot,tcpustep,ntot)

    !--if total number of bins has changed, adjust istepfrac and dt accordingly
    !  (ie., decrease or increase the timestep)
    if (nbinmax /= nbinmaxprev .or. dtmax_ifactor /= 0) then
       call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
    endif

 else
    ! advance time on master thread only
    if (id==master) time = time + dt
    call bcast_mpi(time)
!
!--set new timestep from Courant/forces condition
!
    ! constraint from time to next printout, must reach this exactly
    ! Following redefinitions are to avoid crashing if dtprint = 0 & to reach next output while avoiding round-off errors
    dtprint = min(tprint,tmax) - time + epsilon(dtmax)
    if (dtprint <= epsilon(dtmax) .or. dtprint >= (1.0-1e-8)*dtmax ) dtprint = dtmax + epsilon(dtmax)
    dt = min(dtforce,dtcourant,dterr,dtmax+epsilon(dtmax),dtprint,dtinject,dtrad)
!
!--write log every step (NB: must print after dt has been set in order to identify timestep constraint)
!
    if (id==master) call print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,dtrad,&
                                    dtprint,dtinject,ntot)
 endif

 !   check that MPI threads are synchronised in time
 timecheck = reduceall_mpi('+',time)
 if (abs(timecheck/nprocs - time) > 1.e-13) then
    call fatal('evolve','time differs between MPI threads',var='time',val=timecheck/nprocs)
 endif
!
!--Determine if this is the correct time to write to the data file, because:
!  1) we are at the end of an individual timestep (dtmax, so istepfrac==2**nbinmax),
!  2) we reached the time of the next dump (tprint)
!  3) we reached the end of the simulation (at_simulation_end)
!
 at_dump_time = at_simulation_end(time,nsteps,rhomaxnow)

 if (ind_timesteps) then
    if (istepfrac==2**nbinmax) at_dump_time = .true.
 else
    if (time >= tprint) at_dump_time = .true.
 endif

end subroutine update_time_and_dt

!----------------------------------------------------------------
!+
!  wrapper routine to create sinks and update forces
!+
!----------------------------------------------------------------
subroutine ptmass_create_and_update_forces(time,dtextforce)
 use dim,         only:gravity,use_apr
 use part,        only:nptmass,npart,xyzmh_ptmass,vxyz_ptmass,fxyzu,fext,divcurlv,poten,massoftype,dptmass,&
                       xyzh,vxyzu,fxyz_ptmass,fxyz_ptmass_tree,dsdt_ptmass,group_info,bin_info,nmatrix,&
                       n_group,n_ingroup,n_sing,fxyz_ptmass_sinksink
 use ptmass,      only:ptmass_create_all,use_regnbody,icreate_sinks
 use subgroup,    only:subgroup_search
 use substepping, only:get_force
 real, intent(in)    :: time
 real, intent(inout) :: dtextforce
 integer :: nptmass_old,dummy

 nptmass_old = nptmass
 if (gravity .and. icreate_sinks > 0 .and. .not.use_apr) then
    call ptmass_create_all(nptmass,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,poten,massoftype,&
                           xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink,dptmass,time)
 endif
 ! Need to recompute the force when sink or stars are created
 if (nptmass > nptmass_old) then
    if (use_regnbody) then
       call subgroup_search(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,nmatrix,&
                           new_ptmass=.true.,dtext=dtextforce)
    endif

    dummy = 0
    call get_force(nptmass,npart,0,1,time,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass,vxyz_ptmass,&
                   fxyz_ptmass,fxyz_ptmass_tree,dsdt_ptmass,0.,0.,dummy,.false.,bin_info,&
                   group_info,nmatrix)
 endif

end subroutine ptmass_create_and_update_forces

!----------------------------------------------------------------
!+
!  wrapper routine for writing the .ev files (main .ev file
!  and sink .ev files)
!+
!----------------------------------------------------------------
subroutine write_ev_files(ntot,time,dt,nskip,nskipped,nskipped_sink,at_dump_time)
 use io,             only:iverbose
 use checkconserved, only:check_conservation_errors
 use energies,       only:totmom,angtot,etot,mdust,mtot,hdivBonB_ave,hdivBonB_max,np_e_eq_0,np_cs_eq_0
 use evwrite,        only:write_evfile
 use externalforces, only:iext_spiral
 use options,        only:write_files,iexternalforce
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink
 use ptmass,         only:pt_write_sinkev
 use timing,         only:get_timings,increment_timer,itimer_ev
#ifdef MFLOW
 use mf_write,       only:mflow_write
#endif
#ifdef VMFLOW
 use mf_write,       only:vmflow_write
#endif
#ifdef BINPOS
 use mf_write,       only:binpos_write
#endif
 integer(kind=8), intent(in)    :: ntot
 real,            intent(in)    :: time,dt
 integer,         intent(inout) :: nskip,nskipped,nskipped_sink
 logical,         intent(in)    :: at_dump_time
 real(kind=4)    :: t1,t2,tcpu1,tcpu2
 integer(kind=8) :: nevwrite_threshold,nsinkwrite_threshold

 !
 !  For individual timesteps, we do not want to do this every step, but we want
 !  to do this as often as possible without a performance hit. The criteria
 !  here is that it is done once > 200% of particles (cumulatively) have been evolved.
 !  That is, either 2 steps with all particles, or e.g. 200 steps with 1% of particles.
 !
 nevwrite_threshold    = int(1.99*ntot,kind=8) ! every 2 full steps
 nsinkwrite_threshold  = int(0.99*ntot,kind=8)
 ! do this less often for external forces that are expensive to evaluate
 if (iexternalforce==iext_spiral) nevwrite_threshold = int(4.99*ntot,kind=8) ! every 5 full steps

 nskipped = nskipped + nskip
 if (nskipped >= nevwrite_threshold .or. at_dump_time .or. iverbose==5) then
    nskipped = 0
    call get_timings(t1,tcpu1)
    call write_evfile(time,dt) ! the write_files option is checked inside the routine
    call check_conservation_errors(totmom,angtot,etot,mdust,mtot,hdivBonB_ave,&
                                   hdivBonB_max,np_e_eq_0,np_cs_eq_0)

    !--write with the same ev file frequency also mass flux and binary position
#ifdef MFLOW
    call mflow_write(time,dt)
#endif
#ifdef VMFLOW
    call vmflow_write(time,dt)
#endif
#ifdef BINPOS
    call binpos_write(time,dt)
#endif
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_ev,t2-t1,tcpu2-tcpu1)  ! time taken for write_ev operation
 endif
!
!--print out the sink particle properties to the sink.ev files
!
 nskipped_sink = nskipped_sink + nskip
 if (nskipped_sink >= nsinkwrite_threshold .or. at_dump_time) then
    nskipped_sink = 0
    if (write_files) call pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass,&
                                          fxyz_ptmass,fxyz_ptmass_sinksink)
 endif

end subroutine write_ev_files

!----------------------------------------------------------------
!+
!  wrapper routine for writing the dump file. Mainly determines
!  whether to write a full or small dump, and whether the
!  time between dumps should be reduced if the calculation
!  is taking too long.
!+
!----------------------------------------------------------------
subroutine check_and_write_dump(time,tstart,tcpustart,rhomaxnow,nsteps,&
                                nout,noutput,noutput_dtmax,ncount_fulldumps,&
                                dumpfile,infile,evfile,logfile,abortrun)
 use dim,              only:ind_timesteps,inject_parts,driving,idumpfile,use_apr
 use io,               only:iprint,iwritein,id,master,flush_warnings
 use io_control,       only:at_simulation_end,check_for_full_dump
 use fileutils,        only:getnextfilename
 use forcing,          only:write_forcingdump
 use utils_apr,        only:write_aprtrack
 use options,          only:write_files
 use injection,        only:rkill
 use part,             only:ideadhead,shuffle_part,npart,nptmass,xyzmh_ptmass,accrete_particles_outside_sphere
 use ptmass,           only:calculate_mdot
 use readwrite_infile, only:write_infile
 use readwrite_dumps,  only:write_fulldump,write_smalldump
 use dynamic_dtmax,    only:check_for_restart_dump,idtmax_frac
 use timing,           only:timers,itimer_io,itimer_lastdump,itimer_fromstart,increment_timer,get_timings
#ifdef LIVE_ANALYSIS
 use analysis,         only:do_analysis
 use fileutils,        only:numfromfile
 use io,               only:ianalysis
 use part,             only:massoftype,igas,vxyzu,xyzh
#endif
 real,             intent(in)    :: time,rhomaxnow
 real(kind=4),     intent(inout) :: tstart,tcpustart
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: dumpfile,evfile,logfile
 integer,          intent(in)    :: nsteps,nout
 integer,          intent(inout) :: ncount_fulldumps,noutput,noutput_dtmax
 logical,          intent(out)   :: abortrun
 logical :: fulldump,writedump
 character(len=120) :: dumpfile_orig
 real(kind=4) :: t1,t2,tcpu1,tcpu2

 dumpfile_orig = trim(dumpfile)
!
! setting nout to > 1 will only write the dump file
! every nout steps. We always dump if nout <= 0
!
 writedump = ((nout <= 0) .or. (mod(noutput,nout)==0))
 if (writedump) then
    !--modify evfile and logfile names with new number
    if (noutput==1) then
       evfile  = getnextfilename(evfile)
       logfile = getnextfilename(logfile)
    endif
    ! Update values for restart dumps
    call check_for_restart_dump()

    if (idtmax_frac==0) then
       dumpfile = getnextfilename(dumpfile,idumpfile)
       dumpfile_orig = trim(dumpfile)
    else
       write(dumpfile,'(2a)') dumpfile(:index(dumpfile,'_')-1),'.restart'
    endif
 endif

 !--kill particles before writing dump, also calculate average mdot on sink particles
 if (rkill > 0) call accrete_particles_outside_sphere(rkill)
 if (.not.inject_parts) call calculate_mdot(nptmass,time,xyzmh_ptmass)

 !--do not dump dead particles into dump files
 if (ideadhead > 0) call shuffle_part(npart)
 !
 !--get timings since last dump and overall code scaling
 !  (get these before writing the dump so we can check whether or not we
 !   need to write a full dump based on the wall time)
 !
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_fromstart,t2-tstart,tcpu2-tcpustart)

 fulldump = check_for_full_dump(time,rhomaxnow,nsteps,noutput,timers(itimer_fromstart)%wall,&
                                timers(itimer_lastdump)%wall,abortrun)
!
!--flush any buffered warnings to the log file
!
 if (id==master) call flush_warnings()
!
!--write dump file
!
 call get_timings(t1,tcpu1)
 if (writedump .and. write_files) then
    if (fulldump) then
       call write_fulldump(time,dumpfile)
       if (id==master) then
          call write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
          if (driving) call write_forcingdump(time,dumpfile)
       endif
       ncount_fulldumps = ncount_fulldumps + 1
    else
       call write_smalldump(time,dumpfile)
    endif
    if (use_apr) call write_aprtrack(time,dumpfile)
 endif
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_io,t2-t1,tcpu2-tcpu1)

#ifdef LIVE_ANALYSIS
 if (id==master .and. idtmax_frac==0) then
    call do_analysis(dumpfile,numfromfile(dumpfile),xyzh,vxyzu, &
                     massoftype(igas),npart,time,ianalysis)
 endif
#endif

 dumpfile = trim(dumpfile_orig)

end subroutine check_and_write_dump

!----------------------------------------------------------------
!+
!  print out the log information after a dump
!+
!----------------------------------------------------------------
subroutine print_log(nsteps,nalivetot,nmovedtot,nsteplast,dtmax,tall,tcpulast,twalllast)
 use dim,          only:ind_timesteps
 use io,           only:iprint,id,master,iverbose
 use evwrite,      only:write_evlog
 use part,         only:npart,nptmass,iphase,ibin,xyzh
 use io_summary,   only:summary_printout
 use timing,       only:reduce_timers,timers,itimer_lastdump,reset_timers,get_timings
 use timestep_ind, only:reset_time_per_bin,print_dtind_efficiency,nbinmax,write_binsummary
 integer,         intent(in)    :: nsteps,nsteplast
 integer(kind=8), intent(in)    :: nalivetot,nmovedtot
 real,            intent(in)    :: dtmax
 real(kind=4),    intent(inout) :: tall
 real(kind=4),    intent(out)   :: tcpulast,twalllast

 call reduce_timers
 if (id==master) then
    call print_timinginfo(iprint,nsteps,nsteplast)
    !--Write out summary to log file
    call summary_printout(iprint,nptmass)
 endif

 if (ind_timesteps) then
    !--print summary of timestep bins
    if (iverbose >= 0) then
       call write_binsummary(npart,nbinmax,dtmax,iphase,ibin,xyzh)
       call reset_time_per_bin() ! reset bin timers
       if (id==master) call print_dtind_efficiency(iverbose,nalivetot,nmovedtot,tall,timers(itimer_lastdump)%wall,2)
    endif
 endif

 !--print summary of energies and other useful values to the log file
 if (id==master) call write_evlog(iprint)

 call get_timings(twalllast,tcpulast)
 call reset_timers

end subroutine print_log

!----------------------------------------------------------------
!+
!  routine to print out the timing information at each full dump
!+
!----------------------------------------------------------------
subroutine print_timinginfo(iprint,nsteps,nsteplast)
 use io,     only:formatreal
 use timing, only:timer,timers,print_timer,itimer_fromstart,itimer_lastdump,&
                  itimer_step,itimer_balance,itimer_dens,&
                  itimer_force,itimer_ev,itimer_io,ntimers
 integer,      intent(in) :: iprint,nsteps,nsteplast
 real                     :: dfrac,fracinstep
 real(kind=4)             :: time_fullstep
 character(len=20)        :: string,string1,string2,string3
 integer                  :: itimer

 write(string,"(i12)") nsteps
 call formatreal(real(timers(itimer_fromstart)%wall),string1)
 call formatreal(real(timers(itimer_fromstart)%cpu),string2)
 call formatreal(real(timers(itimer_fromstart)%cpu/(timers(itimer_fromstart)%wall+epsilon(0._4))),string3)
 write(iprint,"(1x,'Since code start: ',a,' timesteps, wall: ',a,'s cpu: ',a,'s cpu/wall: ',a)") &
       trim(adjustl(string)),trim(string1),trim(string2),trim(string3)

 write(string,"(i12)") nsteps-nsteplast
 call formatreal(real(timers(itimer_lastdump)%wall),string1)
 call formatreal(real(timers(itimer_lastdump)%cpu),string2)
 call formatreal(real(timers(itimer_lastdump)%cpu/(timers(itimer_lastdump)%wall+epsilon(0._4))),string3)
 write(iprint,"(1x,'Since last dump : ',a,' timesteps, wall: ',a,'s cpu: ',a,'s cpu/wall: ',a)") &
       trim(adjustl(string)),trim(string1),trim(string2),trim(string3)

 time_fullstep = timers(itimer_lastdump)%wall + timers(itimer_ev)%wall + timers(itimer_io)%wall
 write(iprint,"(/,25x,a)") '  wall         cpu  cpu/wall  load bal      frac'

! skip the first 2 timers
! 1: from start
! 2: from last dump
! 3: step
 do itimer = 3, ntimers
    call print_timer(iprint,itimer,time_fullstep)
 enddo

 dfrac = 1./(timers(itimer_lastdump)%wall + epsilon(0._4))
 fracinstep = timers(itimer_step)%wall*dfrac
 if (fracinstep < 0.99) then
    write(iprint,"(1x,a,f6.2,a)") 'WARNING: ',100.*(1.-fracinstep),'% of time was in unusual routines (not dens/force/tree)'
 endif

end subroutine print_timinginfo

end module evolve_utils
