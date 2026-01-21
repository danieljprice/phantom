!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module evolve
!
! Evolves the simulation through all timesteps
!   This subroutine contains the main timestepping loop and calls
!   the output routines at the appropriate times
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: HIIRegion, apr, boundary_dyn, centreofmass,
!   checkconserved, dim, dynamic_dtmax, easter_egg, energies, evolve_utils,
!   forcing, inject, io, io_control, io_summary, mpiutils, part,
!   partinject, ptmass, radiation_utils, step_lf_global, timestep,
!   timestep_ind, timing
!
 use dim, only:ind_timesteps
 implicit none
 public :: evol,evol_init,evol_prestep,evol_poststep

 private
 logical      :: initialized = .false.

 ! global counters and state variables
 integer         :: nsteplast
 integer         :: nskip,nskipped,nskipped_sink
 integer         :: noutput,noutput_dtmax,ncount_fulldumps
 real(kind=4)    :: tcpustart,tstart,tall
 real(kind=4)    :: twalllast,tcpulast
 real            :: tprint,tcheck,tzero,tlast,dtlast,rhomaxold,dtmaxold
 logical         :: abortrun_bdy,use_global_dt
 integer(kind=8) :: nmovedtot,nalivetot
 integer(kind=1) :: nbinmaxprev

contains

!----------------------------------------------------------------
!+
!  initialise counters and global variables used
!  in the evolve subroutine. Note that in general the
!  routine initialise() should be called prior to calling evol
!+
!----------------------------------------------------------------
subroutine evol_init(time,dtmax,rhomaxnow,dt)
 use io,             only:iprint
 use checkconserved, only:init_conservation_checks
 use part,           only:npart
 use ptmass,         only:set_integration_precision
 use step_lf_global, only:init_step
 use timestep,       only:dtinject,dtrad
 use timestep_ind,   only:reset_time_per_bin,nbinmax
 use timing,         only:get_timings,setup_timers
 real, intent(in)    :: time,dtmax,rhomaxnow
 real, intent(inout) :: dt

 dtlast    = 0.
 dtinject  = huge(dtinject)
 dtrad     = huge(dtrad)
 abortrun_bdy = .false.
 rhomaxold = rhomaxnow

 call init_conservation_checks()
 call init_counters(time,dtmax) ! set istepfrac to 0 as well as tprint and noutput
 !
 ! Set substepping integration precision depending on the system (default is FSI)
 !
 call set_integration_precision

 if (ind_timesteps) then
    use_global_dt = .false.
    dt = dtmax/2.**nbinmax  ! use 2.0 here to allow for step too small
    call reset_time_per_bin() ! initialise bin timers
    call init_step(npart,time,dtmax)
 else
    use_global_dt = .true.
    nbinmax   = 0 ! dummy value
    if (dt >= (tprint-time)) dt = tprint-time   ! reach tprint exactly
 endif
!
! code timings
!
 call get_timings(twalllast,tcpulast)
 tstart    = twalllast
 tcpustart = tcpulast

 call setup_timers
 call flush(iprint)

 initialized = .true.

end subroutine evol_init

!----------------------------------------------------------------
!+
!  evolve the simulation over all timesteps
!+
!----------------------------------------------------------------
subroutine evol(infile,logfile,evfile,dumpfile,flag)
 use dim,              only:do_radiation
 use io_control,       only:at_simulation_end
 use part,             only:npart,xyzh,fxyzu,vxyzu,rad,radprop
 use radiation_utils,  only:update_radenergy,exchange_radiation_energy,implicit_radiation
 use step_lf_global,   only:step
 use timestep,         only:time,dt,dtmax,nsteps,dtextforce,rhomaxnow
 use timestep_ind,     only:nactive
 integer, optional, intent(in)   :: flag
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: logfile,evfile,dumpfile
 real            :: dtnew
 real(kind=4)    :: t1,tcpu1
 logical         :: do_radiation_update,abortrun

 ! the following isrequired because evol is called multiple times in AMUSE... -SR
 if (.not. initialized) call evol_init(time,dtmax,rhomaxnow,dt)

 ! logical checks
 do_radiation_update = do_radiation .and. exchange_radiation_energy .and. .not.implicit_radiation

 !
 ! main timestepping loop
 !
 timestepping: do while (.not. at_simulation_end(time,nsteps,rhomaxnow))

    call evol_prestep(time,dtmax,dt,t1,tcpu1,nactive,present(flag))
    !
    ! Strang splitting: implicit update for half step
    !
    if (do_radiation_update) call update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,0.5*dt)
    !
    !--evolve data for one timestep
    !  for individual timesteps this is the shortest timestep
    !
    call step(npart,nactive,time,dt,dtextforce,dtnew)
    !
    ! Strang splitting: implicit update for another half step
    !
    if (do_radiation_update) call update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,0.5*dt)

    call evol_poststep(infile,logfile,evfile,dumpfile,&
                       time,t1,tcpu1,dt,dtmax,nactive,abortrun)
    if (abortrun) exit

 enddo timestepping

end subroutine evol

!----------------------------------------------------------------
!+
!  wrapper routine for everything done before the step call
!+
!----------------------------------------------------------------
subroutine evol_prestep(time,dtmax,dt,t1,tcpu1,nactive,inject_flag_present)
 use dim,          only:inject_parts,use_apr,ind_timesteps
 use io,           only:fatal,id,master,iprint,iverbose
 use apr,          only:update_apr
 use boundary_dyn, only:dynamic_bdy,update_boundaries
 use dynamic_dtmax,only:check_dtmax_for_decrease
 use evolve_utils, only:ptmass_create_and_update_forces
 use inject,       only:inject_particles
 use HIIRegion,    only:HII_feedback,iH2R,HIIuprate,nHIIsources
 use io_control,   only:nfulldump
 use mpiutils,     only:reduceall_mpi
 use part,         only:npart,npartoftype,nptmass,xyzh,vxyzu,fxyzu,apr_level,&
                        xyzmh_ptmass,vxyz_ptmass,gravity,iboundary,ntot,ibin,iphase,&
                        isionised
 use partinject,   only:update_injected_particles
 use ptmass,       only:icreate_sinks,ipart_createstars
 use timestep,     only:dtextforce,dtinject,rhomaxnow
 use timestep_ind, only:istepfrac,nbinmax,set_active_particles,write_binsummary,nactivetot,maxbins
 use timing,       only:get_timings,timers,itimer_lastdump
 real,            intent(in)    :: time
 real,            intent(inout) :: dt,dtmax
 real(kind=4),    intent(out)   :: t1,tcpu1
 integer,         intent(out)   :: nactive
 logical,         intent(in)    :: inject_flag_present
 real(kind=4) :: twallperdump
 integer :: npart_old,nalive,istepHII
 !
 ! injection of new particles into simulation
 !
 if (inject_parts .and. .not. inject_flag_present) then
    npart_old = npart
    call inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject)
    dtlast = dt
 endif

 if (use_apr) call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level) ! split or merge as required

 dtmaxold = dtmax
 nbinmaxprev = nbinmax  ! for global timesteps this is ignored

 if (ind_timesteps) then
    istepfrac   = istepfrac + 1
    if (nbinmax > maxbins) call fatal('evolve','timestep too small: try decreasing dtmax?')

    !--determine if dt needs to be decreased; if so, then this will be done
    !  in step the next time it is called;
    !  for global timestepping, this is called in the block where at_dump_time==.true.
    if (istepfrac == 2**nbinmax) then
       twallperdump = reduceall_mpi('max',timers(itimer_lastdump)%wall)
       call check_dtmax_for_decrease(iprint,time,dtmax,twallperdump,&
                                     rhomaxold,rhomaxnow,nfulldump,use_global_dt)
    endif

    !--sanity check on istepfrac...
    if (istepfrac > 2**nbinmax) then
       write(iprint,*) 'ERROR: istepfrac = ',istepfrac,' / ',2**nbinmax
       call fatal('evolve','error in individual timesteps')
    endif

    !--flag particles as active or not for this timestep
    call set_active_particles(npart,nactive,nalive,nactivetot,nalivetot,iphase,ibin,xyzh)
    nskip = int(nactivetot)

    !--print summary of timestep bins
    if (iverbose >= 2) call write_binsummary(npart,nbinmax,dtmax,iphase,ibin,xyzh)

    !--Implement dynamic boundaries (for individual-timestepping) once per dump
    if (dynamic_bdy .and. nactive==nalive .and. istepfrac==2**nbinmax) then
       call update_boundaries(nactive,nalive,npart,abortrun_bdy)
    endif
 else
    !--for global timestep, set nskip to total number of particles across all nodes
    nskip = int(ntot)
    nactive = npart
 endif

 if (gravity .and. icreate_sinks > 0) call ptmass_create_and_update_forces(time,dtextforce)

 if (iH2R > 0 .and. nHIIsources > 0 .and. id==master) then
    istepHII = 1
    if (ind_timesteps) then
       istepHII = 2**nbinmax/HIIuprate
       if (istepHII==0) istepHII = 1
    endif
    if (mod(istepfrac,istepHII) == 0 .or. istepfrac == 1 .or. (icreate_sinks == 2 .and. ipart_createstars /= 0)) then
       call HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,isionised)
    endif
 endif

 call get_timings(t1,tcpu1)

end subroutine evol_prestep

!----------------------------------------------------------------
!+
!  wrapper routine for everything done after the step call
!+
!----------------------------------------------------------------
subroutine evol_poststep(infile,logfile,evfile,dumpfile,time,t1,tcpu1,dt,dtmax,nactive,abortrun)
 use dim,           only:driving
 use io,            only:id,master,iprint
 use boundary_dyn,  only:dynamic_bdy,update_boundaries
 use centreofmass,  only:correct_bulk_motions
 use dynamic_dtmax, only:check_dtmax_for_decrease
 use easter_egg,    only:egged,bring_the_egg
 use evolve_utils,  only:update_time_and_dt,write_ev_files,check_and_write_dump,print_log
 use forcing,       only:correct_bulk_motion
 use io_control,    only:nout,nmaxdumps,twallmax,nfulldump
 use io_summary,    only:iosum_nreal,summary_counter,summary_printout,summary_printnow
 use part,          only:npart,nptmass,ntot
 use timing,        only:get_timings,increment_timer,itimer_step,itimer_lastdump,timers,print_time
 use timestep,      only:nsteps,rhomaxnow
 use timestep_ind,  only:istepfrac
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: logfile,evfile,dumpfile
 real,             intent(inout) :: time
 real(kind=4),     intent(in)    :: t1,tcpu1
 real,             intent(inout) :: dt,dtmax
 integer,          intent(inout) :: nactive
 logical,          intent(out)   :: abortrun
 real(kind=4) :: t2,tcpu2,twallperdump
 logical      :: at_dump_time,iexist

 abortrun =.false.

 !--timings for step call
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_step,t2-t1,tcpu2-tcpu1)
 call summary_counter(iosum_nreal,t2-t1)
 !
 !--update time = time + dt, update step counters, and get the dt for the next timestep
 !
 call update_time_and_dt(nsteps,time,dtmax,dtmaxold,rhomaxnow,tlast,tcheck,tprint,dt,tall,t2-t1,tcpu2-tcpu1,&
                         istepfrac,nbinmaxprev,ntot,nalivetot,nmovedtot,at_dump_time)
 !
 !--Update timer from last dump to see if dtmax needs to be reduced
 !
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_lastdump,t2-t1,tcpu2-tcpu1)
 !
 !--Calculate total energy etc and write to ev file. Do this before writing dumps
 !   so that values calculated in energies are correctly included in the dumpfiles
 !
 call write_ev_files(ntot,time,dt,nskip,nskipped,nskipped_sink,at_dump_time)
 !
 !--write to data file if time is right
 !
 if (at_dump_time) then
    !
    !--Global timesteps: Decrease dtmax if requested (done in step for individual timesteps)
    !
    if (.not. ind_timesteps) then
       twallperdump = timers(itimer_lastdump)%wall
       call check_dtmax_for_decrease(iprint,time,dtmax,twallperdump,rhomaxold,rhomaxnow,nfulldump,use_global_dt)
       dt = min(dt,dtmax) ! required if decreasing dtmax to ensure that the physically motivated timestep is not too long
    endif

    call check_and_write_dump(time,tstart,tcpustart,rhomaxnow,nsteps,&
                              nout,noutput,noutput_dtmax,ncount_fulldumps,&
                              dumpfile,infile,evfile,logfile,abortrun)

    ! reset counters for when the next dump should be written
    call print_log(nsteps,nalivetot,nmovedtot,nsteplast,dtmax,tall,tcpulast,twalllast)
    call update_dump_counters(nsteps,dtmax)

    !--Implement dynamic boundaries (for global timestepping)
    if (.not. ind_timesteps .and. dynamic_bdy) call update_boundaries(nactive,nactive,npart,abortrun_bdy)
    if (abortrun_bdy) return
    !
    !--if twallmax > 1s stop the run at the last full dump that will fit into the walltime constraint,
    !  based on the wall time between the last two dumps added to the current total walltime used.
    !
    if (abortrun) then
       if (id==master) then
          call print_time(t2-tstart,'>> WALL TIME = ',iprint)
          call print_time(twallmax,'>> NEXT DUMP WILL TRIP OVER MAX WALL TIME: ',iprint)
          write(iprint,"(1x,a)") '>> ABORTING... '
       endif
       return
    endif

    if (nmaxdumps > 0 .and. ncount_fulldumps >= nmaxdumps) then
       if (id==master) write(iprint,"(a)") '>> reached maximum number of full dumps as specified in input file, stopping...'
       abortrun = .true.
       return
    endif

 endif

 if (driving .and. correct_bulk_motion) call correct_bulk_motions()

 call flush(iprint)
 !--Write out log file prematurely (if requested based upon nstep, walltime)
 if ( summary_printnow() ) call summary_printout(iprint,nptmass)

 !--???
 inquire(file='egg.txt',exist=iexist)
 if (iexist .and. .not.egged) then
    call bring_the_egg
    egged = .true.
 endif
 if (.not.iexist) egged = .false.

end subroutine evol_poststep

!----------------------------------------------------------------
!+
!  initialize various counters
!+
!----------------------------------------------------------------
subroutine init_counters(time,dtmax)
 use energies,     only:np_cs_eq_0,np_e_eq_0
 use timestep,     only:nsteps
 use timestep_ind, only:istepfrac
 real, intent(in) :: time,dtmax

! total number of steps and number of steps since last dump
 nsteps    = 0
 nsteplast = 0
 tzero = time  ! starting time
 tlast = time  ! time of last dump

! number of snapshots that have been written
 noutput          = 1
 noutput_dtmax    = 1
 ncount_fulldumps = 0

! time to print next
 tprint = time + dtmax
 dtmaxold = dtmax

! counters for individual timesteps
 istepfrac = 0
 if (ind_timesteps) then
    nmovedtot = 0
    tall      = 0.
    tcheck    = time
 endif

! warning counters
 np_cs_eq_0 = 0
 np_e_eq_0  = 0

!
! threshold for writing to .ev file, to avoid repeatedly computing energies
! for all the particles which would add significantly to the cpu time
!
 nskipped = 0
 nskipped_sink = 0
 nskip = 0

end subroutine init_counters

!----------------------------------------------------------------
!+
!  update the dump counters
!+
!----------------------------------------------------------------
subroutine update_dump_counters(nsteps,dtmax)
 use dynamic_dtmax, only:idtmax_n_next,idtmax_frac,idtmax_frac_next,&
                         dtmax_ifactorWT,dtmax_ifactor,idtmax_n
 use timestep_ind,  only:istepfrac
 integer, intent(in)    :: nsteps
 real,    intent(in)    :: dtmax

 ! number of steps since last dump
 nsteplast = nsteps
 ! time of last dump
 tlast = tprint

 ! reset counters for individual timesteps
 if (ind_timesteps) then
    istepfrac = 0
    nmovedtot = 0
 endif

 if (idtmax_frac==0) noutput = noutput + 1  ! required to determine frequency of full dumps
 noutput_dtmax = noutput_dtmax + 1     ! required to adjust tprint; will account for varying dtmax

! factors associated with dynamically changing dtmax
 idtmax_n    = idtmax_n_next
 idtmax_frac = idtmax_frac_next

! time to print next dump
 tprint = tzero + noutput_dtmax*dtmaxold

! adjustments to tprint and other counters if dtmax has changed
 if (dtmax_ifactor /= 0) then
    tzero           = tprint - dtmaxold
    tprint          = tzero  + dtmax
    noutput_dtmax   = 1
    dtmax_ifactor   = 0
    dtmax_ifactorWT = 0
 endif

end subroutine update_dump_counters

end module evolve
