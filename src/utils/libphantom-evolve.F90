!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: evolvesplit
!
!  DESCRIPTION:
!  evolves the simulation through all timesteps
!  this subroutine contains the main timestepping loop and calls
!  the output routines at the appropriate times
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, energies, evwrite, fileutils, forcing,
!    inject, io, mpiutils, options, part, quitdump, readwrite_dumps,
!    readwrite_infile, sort_particles, step_lf_global, timestep,
!    timestep_ind, timing
!+
!--------------------------------------------------------------------------
module evolvesplit
#if IND_TIMESTEPS
 use timestep_ind, only:maxbins
#endif
 implicit none
 public :: evol, evol_init, init_step, finalize_step, dtlast

 private
 integer         :: noutput,nsteplast,ncount_fulldumps
 real            :: dtlast
 real            :: tprint,tzero
 real(kind=4)    :: t1,t2,tcpu1,tcpu2,tcpu,tcputot,tstart,tcpustart
 real(kind=4)    :: twalllast,tcpulast,tstep,tcpustep,twallperdump,twallused
 real(kind=4)    :: tev,tcpuev
#ifdef IND_TIMESTEPS
 integer(kind=1) :: nbinmaxprev
 integer(kind=8) :: nmovedtot,ntot
 real            :: tlast,tcheck
 real(kind=4)    :: tall
 real(kind=4) :: timeperbin(0:maxbins)
#endif
 integer :: nskip,nskipped,nevwrite_threshold

contains

subroutine evol_init()
 use timestep, only:time,dt,dtmax,nsteps
 use part, only:npart
 use options, only:iexternalforce
 use timing, only:get_timings
 use io, only:iprint
#ifdef IND_TIMESTEPS
 use part,             only:ibin
 use timestep_ind,     only:istepfrac,nbinmax
 use timestep,         only:restartonshortest
#endif

#ifdef IND_TIMESTEPS
 integer :: i
#endif

 tprint    = 0.
 nsteps    = 0
 nsteplast = 0
 tzero     = time
 dtlast    = 0.

 noutput          = 1
 ncount_fulldumps = 0
 tprint = tzero + dtmax
#ifdef IND_TIMESTEPS
 istepfrac = 0
 tlast = tzero
 dt = dtmax/2**nbinmax
 nmovedtot = 0
 tall = 0.
 tcheck = time
 timeperbin(:) = 0.
!
! first time through, move all particles on shortest timestep
! then allow them to gradually adjust levels
!
 if (time < tiny(time) .or. restartonshortest) then
    !$omp parallel do schedule(static) private(i)
    do i=1,npart
       ibin(i) = nbinmax
    enddo
 endif

#else
 nskip = npart
 if (dt >= (tprint-time)) dt = tprint-time   ! reach tprint exactly
#endif
!
! threshold for writing to .ev file, to avoid repeatedly computing energies
! for all the particles which would add significantly to the cpu time
!
 nskipped = 0
 if (iexternalforce==8) then
    nevwrite_threshold = int(4.99*npart) ! every 5 full steps
 else
    nevwrite_threshold = int(1.99*npart) ! every 3 full steps
 endif
!
! timing between dumps
!
 call get_timings(twalllast,tcpulast)
 tstart = twalllast
 tcpustart = tcpulast
 tstep = 0.
 tcpustep = 0.
 tev = 0.
 tcpuev = 0.
 call flush(iprint)
end subroutine evol_init

subroutine init_step()
 use timing,           only:get_timings
 use io,               only:iprint,fatal,iverbose
 use timestep,         only:nsteps,dtmax
 use part,             only:npart,xyzh
 use mpiutils,         only:reduce_mpi
#ifdef IND_TIMESTEPS
 use part,             only:ibin,iphase
 use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,nactive,nactivetot,write_binsummary
#endif

#ifdef IND_TIMESTEPS
 integer :: nalive

 istepfrac = istepfrac + 1
 nbinmaxprev = nbinmax

 !--sanity check on istepfrac...
 if (istepfrac > 2**nbinmax) then
    write(iprint,*) 'ERROR: istepfrac = ',istepfrac,' / ',2**nbinmax
    call fatal('evolve','error in individual timesteps')
 endif

 !--flag particles as active or not for this timestep
 call set_active_particles(npart,nactive,nalive,iphase,ibin,xyzh)
 nactivetot = reduce_mpi('+',nactive)
 ntot = reduce_mpi('+',nalive)
 nskip = nactivetot

 !--print summary of timestep bins
 if (iverbose >= 2) call write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)
#endif

 nsteps = nsteps + 1
!
!--evolve data for one timestep
!  for individual timesteps this is the shortest timestep
!
 call get_timings(t1,tcpu1)
end subroutine init_step

subroutine finalize_step(infile, logfile, evfile, dumpfile)
 use timing,           only:get_timings,print_time
 use timestep,         only:time,tmax,dt,dtmax,nsteps,nmax,nout
 use io,               only:id,master,iverbose,iprint,warning,iwritein,flush_warnings
 use fileutils,        only:getnextfilename
 use quitdump,         only:quit
 use readwrite_dumps,  only:write_smalldump,write_fulldump
 use evwrite,          only:write_evfile,write_evlog
 use energies,         only:etot,totmom
 use mpiutils,         only:reduce_mpi
 use options,          only:nfulldump,twallmax,nmaxdumps
 use readwrite_infile, only:write_infile
#ifdef IND_TIMESTEPS
 use part,             only:npart,iphase,ibin,xyzh
 use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,update_time_per_bin,&
                            write_binsummary,change_nbinmax,nactivetot,maxbins
 use io,               only:fatal,warning
#else
 use timestep,         only:dtforce,dtcourant,dterr
#endif
#ifdef DRIVING
 use forcing,          only:write_forcingdump
#endif
#ifdef SORT
 use sort_particles,   only:sort_part
#endif
#ifdef CORRECT_BULK_MOTION
 use centreofmass,     only:correct_bulk_motion
#endif
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: logfile,evfile,dumpfile

#ifdef IND_TIMESTEPS
 integer :: inbin
 real :: fracactive, speedup
#endif
 logical :: fulldump, abortrun, at_dump_time

 !--timings for step call
 call get_timings(t2,tcpu2)
 tstep = tstep + t2-t1
 tcpustep = tcpustep + tcpu2-tcpu1

#ifdef IND_TIMESTEPS
 tcheck = tcheck + dt

 !--update time in way that is free of round-off errors
 time = tlast + istepfrac/real(2**nbinmaxprev)*dtmax

 !--print efficiency of partial timestep
 if (id==master .and. iverbose >= 0 .and. ntot > 0) then
    if (nactivetot==ntot) then
       tall = t2-t1
    elseif (tall > 0.) then
       fracactive = nactivetot/real(ntot)
       speedup = (t2-t1)/tall
       if (iverbose >= 2) &
             write(iprint,"(1x,'(',3(a,f6.2,'%'),')')") &
                  'moved ',100.*fracactive,' of particles in ',100.*speedup, &
                  ' of time, efficiency = ',100.*fracactive/speedup
    endif
 endif
 call update_time_per_bin(tcpu2-tcpu1,istepfrac,nbinmaxprev,timeperbin,inbin)
 nmovedtot = nmovedtot + nactivetot

 !--check that time is as it should be, may indicate error in individual timestep routines
 if (abs(tcheck-time) > 1.e-4) call warning('evolve','time out of sync',var='error',val=abs(tcheck-time))

 if (id==master .and. (iverbose >= 1 .or. inbin <= 3)) &
       write(iprint,5) istepfrac,2**nbinmaxprev,time,dt,nactivetot,tcpu2-tcpu1
5 format('> step ',i6,' /',i6,2x,'t = ',f10.5,2x,'dt = ',es10.3,' moved ',i10,' in ',f8.2,' cpu-s <')

 !--if total number of bins has changed, adjust istepfrac and dt accordingly
 !  (ie., decrease or increase the timestep)
 if (nbinmax /= nbinmaxprev) call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)

#else
 time = time + dt
 dtlast = dt
!
!--set new timestep from Courant/forces condition
!
 dt = min(dtforce,dtcourant,dterr,dtmax+epsilon(dtmax)) ! use dtmax + eps to avoid round-off problems

!
!--write log every step (NB: must print after dt has been set in order to identify timestep constraint)
!
 if (id==master) call print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax)

#endif

 if (abs(dt) < 1e-8*dtmax) then
    write(iprint,*) 'main loop: timestep too small, dt = ',dt
    call quit   ! also writes dump file, safe here because called by all threads
 endif
!
!--write to data file if time is right
!
 at_dump_time = .false.
#ifdef IND_TIMESTEPS
 if ( (istepfrac==2**nbinmax) &
#else
    if ( (time >= tprint)             &
#endif
     .or.(time >= tmax)            &
     .or.((mod(nsteps,nout)==0).and.(nout > 0)) &
     .or.(nsteps >= nmax)) then
 !--modify evfile and logfile names with new number
 if (noutput==1) then
    evfile = getnextfilename(evfile)
    logfile = getnextfilename(logfile)
 endif
 dumpfile = getnextfilename(dumpfile)

#ifdef MPI
 !--do not dump dead particles into dump files
 if (ideadhead > 0) call shuffle_part(npart)
#endif

!
!--get timings since last dump and overall code scaling
!  (get these before writing the dump so we can check whether or not we
!   need to write a full dump based on the wall time)
!
 call get_timings(t2,tcpu2)
 tcpu     = tcpu2-tcpulast
 tcpu     = reduce_mpi('+',tcpu)
 tcputot  = tcpu2 - tcpustart
 tcputot  = reduce_mpi('+',tcputot)
 tcpustep = reduce_mpi('+',tcpustep)
 tcpuev   = reduce_mpi('+',tcpuev)

 fulldump = (mod(noutput,nfulldump)==0)
!
!--if max wall time is set (> 1 sec) stop the run at the last full dump
!  that will fit into the walltime constraint, based on the wall time between
!  the last two dumps added to the current total walltime used.
!  If we are about to write a small dump but it looks like we won't make the next dump,
!  dump a full dump instead and stop the run
!
 abortrun = .false.
 if (twallmax > 1.) then
    twallused = t2-tstart
    twallperdump = t2-twalllast
    if (fulldump) then
       if ((twallused + abs(nfulldump)*twallperdump) > twallmax) then
          abortrun = .true.
       endif
    else
       if ((twallused + twallperdump) > twallmax) then
          fulldump = .true.
          abortrun = .true.
          write(iprint,"(1x,a)") '>> PROMOTING DUMP TO FULL DUMP BASED ON WALL TIME CONSTRAINTS... '
          nfulldump = 1  !  also set all future dumps to be full dumps (otherwise gets confusing)
       endif
    endif
 endif
!
!--flush any buffered warnings to the log file
!
 if (id==master) call flush_warnings()
!
!--write dump file
!
 at_dump_time = .true.
 if (fulldump) then
    call write_fulldump(time,dumpfile)
    if (id==master) then
       call write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
#ifdef DRIVING
       call write_forcingdump(time,dumpfile)
#endif
    endif
    ncount_fulldumps = ncount_fulldumps + 1

#ifndef IND_TIMESTEPS
#ifdef SORT
    if (time < tmax) call sort_part()
#endif
#endif
 else
    call write_smalldump(time,dumpfile)
 endif

 if (id==master) call print_timinginfo(iprint,nsteps,nsteplast,t2,tstart,&
                              tcpu,tcpustep,tcputot,twalllast,tstep,tev,tcpuev)
#ifdef IND_TIMESTEPS
 !--print summary of timestep bins
 if (iverbose >= 0 .and. id==master .and. abs(tall) > tiny(tall) .and. ntot > 0) then
    fracactive = nmovedtot/real(ntot)
    speedup = (t2-twalllast)/(tall + tiny(tall))
    write(iprint,"(/,a,f6.2,'%')") ' IND TIMESTEPS efficiency = ',100.*fracactive/speedup
    if (iverbose >= 1) then
       write(iprint,"(a,1pe14.2,'s')") '  wall time per particle (last full step) : ',tall/real(ntot)
       write(iprint,"(a,1pe14.2,'s')") '  wall time per particle (ave. all steps) : ',(t2-twalllast)/real(nmovedtot)
    endif
 endif
 if (iverbose >= 0) then
    call write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)
    timeperbin(:) = 0.
 endif
 tlast = tprint
 istepfrac = 0
 nmovedtot = 0
#endif
 !
 !--if twallmax > 1s stop the run at the last full dump that will fit into the walltime constraint,
 !  based on the wall time between the last two dumps added to the current total walltime used.
 !
 if (abortrun) then
    call print_time(t2-tstart,'>> WALL TIME = ',iprint)
    call print_time(twallmax,'>> NEXT DUMP WILL TRIP OVER MAX WALL TIME: ',iprint)
    write(iprint,"(1x,a)") '>> ABORTING... '
    return
 endif

 if (nmaxdumps > 0 .and. ncount_fulldumps >= nmaxdumps) then
    write(iprint,"(a)") '>> reached maximum number of full dumps as specified in input file, stopping...'
    return
 endif

 twalllast = t2
 tcpulast = tcpu2
 tstep = 0.
 tcpustep = 0.
 tev = 0.
 tcpuev = 0.

 noutput = noutput + 1
 tprint = tzero + noutput*dtmax
 nsteplast = nsteps
endif
!
!--Calculate total energy etc and write to ev file
!  For individual timesteps, we do not want to do this every step, but we want
!  to do this as often as possible without a performance hit. The criteria
!  here is that it is done once >10% of particles (cumulatively) have been evolved.
!  That is, either >10% are being stepped, or e.g. 1% have moved 10 steps.
!
nskipped = nskipped + nskip
if (nskipped >= nevwrite_threshold .or. at_dump_time) then
   nskipped = 0
   call get_timings(t1,tcpu1)
   call write_evfile(time,dt)
   !--timings for step call
   call get_timings(t2,tcpu2)
   tev= tev + t2-t1
   tcpuev = tcpuev + tcpu2-tcpu1
   if (at_dump_time .and. id==master) call write_evlog(iprint)
   if (id==master .and. iverbose >= 1) write(iprint,*) ' etot = ',etot,' momtot = ',totmom
endif

#ifdef CORRECT_BULK_MOTION
call correct_bulk_motion()
#endif

#ifndef IND_TIMESTEPS
 !--reach tprint exactly. must take this out for integrator to be symplectic
if (dt >= (tprint-time)) dt = tprint-time + epsilon(0.)
#endif

if (iverbose >= 1 .and. id==master) write(iprint,*)
call flush(iprint)
end subroutine finalize_step

subroutine evol(infile,logfile,evfile,dumpfile)
 use timestep,         only:time,tmax,nmax,nsteps,dt,dtextforce
#ifdef INJECT_PARTICLES
 use inject,           only:inject_particles
#endif
 use step_lf_global,   only:step
 use part, only:npart
#ifdef IND_TIMESTEPS
 use timestep_ind, only:nactive
#endif
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: logfile,evfile,dumpfile
#ifndef IND_TIMESTEPS
 integer :: nactive
#endif
 real :: dtnew

 call evol_init()

 timestepping: do while ((time < tmax).and.(nsteps < nmax))

    call init_step()

#ifdef INJECT_PARTICLES
    !
    !--injection of new particles into simulation
    !
    call inject_particles(time,dtlast)
#endif

#ifndef IND_TIMESTEPS
    nactive = npart
#endif
    call step(npart,nactive,time,dt,dtextforce,dtnew)


    call finalize_step(infile,logfile,evfile,dumpfile)

 enddo timestepping

!------------------------------------------------------------------------

 return
end subroutine evol

!----------------------------------------------------------------
!+
!  routine to print out the timing information at each full dump
!+
!----------------------------------------------------------------
subroutine print_timinginfo(iprint,nsteps,nsteplast,t2,tstart,tcpu,tcpustep,&
                            tcputot,twalllast,tstep,tev,tcpuev)
 use io, only:formatreal
 integer,      intent(in) :: iprint,nsteps,nsteplast
 real(kind=4), intent(in) :: t2,tstart,tcpu,tcpustep,tcputot,twalllast,tstep,tev,tcpuev
 real :: dfrac,fracinstep
 character(len=20) :: string,string1,string2,string3

 !write(iprint,*) 'Timings:'
 write(string,"(i12)") nsteps
 call formatreal(real(t2-tstart),string1)
 call formatreal(real(tcputot),string2)
 call formatreal(real(tcputot/(t2-tstart+epsilon(t2))),string3)
 write(iprint,"(1x,'Since code start: ',a,' timesteps, wall: ',a,'s cpu: ',a,'s cpu/wall: ',a)") &
       trim(adjustl(string)),trim(string1),trim(string2),trim(string3)

 write(string,"(i12)") nsteps-nsteplast
 call formatreal(real(t2-twalllast),string1)
 call formatreal(real(tcpu),string2)
 call formatreal(real(tcpu/(t2-twalllast+epsilon(twalllast))),string3)
 write(iprint,"(1x,'Since last dump : ',a,' timesteps, wall: ',a,'s cpu: ',a,'s cpu/wall: ',a)") &
       trim(adjustl(string)),trim(string1),trim(string2),trim(string3)


 dfrac = 1./(t2-twalllast + epsilon(t2))
 fracinstep = tstep*dfrac
 if (fracinstep < 0.99) then
    write(iprint,20) 100.*fracinstep,'step     ',tstep,tcpustep/(tstep + epsilon(tstep))
    write(iprint,20) 100.*tev*dfrac, 'write_ev ',tev,tcpuev/(tev + epsilon(tev))
    write(iprint,"(1x,a,f6.2,a)") 'WARNING: spending ',100.*(1.-fracinstep),'% of time in routines EXTERNAL to step'
 endif
20 format (1x,'|',f6.2,'% ',a,f8.2,'s, cpu/wall = ',f6.2)
 write(iprint,*)

end subroutine print_timinginfo

!-----------------------------------------------------------------
!+
!  routine to print out the timestep information to the log file
!+
!-----------------------------------------------------------------
subroutine print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax)
 integer, intent(in) :: iprint
 real,    intent(in) :: time,dt,dtforce,dtcourant,dterr,dtmax

 if (abs(dt-dtforce) < tiny(dt)) then
    write(iprint,10) time,dt,'(force)'
 elseif (abs(dt-dtcourant) < tiny(dt)) then
    write(iprint,10) time,dt,'(courant)'
 elseif (abs(dt-dterr) < tiny(dt)) then
    write(iprint,10) time,dt,'(error)'
 elseif (abs(dt-dtmax) < tiny(dt)) then
    write(iprint,10) time,dt,'(dtmax)'
 else
    !print*,dt,dtforce,dtcourant,dterr,dtmax
    write(iprint,10) time,dt,'(unknown)'
 endif
10 format(' t = ',f9.4,' dt = ',es10.3,1x,a)

end subroutine print_dtlog

end module evolvesplit
