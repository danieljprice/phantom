!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Dependencies: HIIRegion, analysis, apr, boundary_dyn, centreofmass,
!   checkconserved, dim, easter_egg, energies, evwrite, externalforces,
!   fileutils, forcing, inject, io, io_summary, mf_write, mpiutils,
!   options, part, partinject, ptmass, quitdump, radiation_utils,
!   readwrite_dumps, readwrite_infile, step_lf_global, subgroup,
!   substepping, supertimestep, timestep, timestep_ind, timestep_sts,
!   timing
!
 use dim, only:ind_timesteps
 implicit none
 public :: evol

 private
 logical      :: initialized = .false.
 integer      :: nsteplast
 integer      :: nskip,nskipped,nskipped_sink
 integer      :: noutput,noutput_dtmax,ncount_fulldumps
 integer      :: istepfrac
 real(kind=4) :: tcpustart,tstart,tall
 real(kind=4) :: twalllast,tcpulast
 real         :: tprint,tcheck,tlast,dtau,dtlast,rhomaxold
 logical      :: abortrun_bdy,use_global_dt
 integer(kind=8) :: nmovedtot

contains

!----------------------------------------------------------------
!+
!  initialise counters and global variables used
!  in the evolve subroutine. Note that in general the
!  routine initialise() should be called prior to calling evol
!+
!----------------------------------------------------------------
subroutine evol_init(tzero,time,dtmax,rhomaxnow,dt)
 use io,             only:iprint
 use checkconserved, only:init_conservation_checks
 use energies,       only:np_cs_eq_0,np_e_eq_0
 use part,           only:npart,ntot
 use ptmass,         only:set_integration_precision
 use step_lf_global, only:init_step
 use timestep,       only:dtinject,dtrad,dtdiff,nsteps
 use timestep_ind,   only:reset_time_per_bin,nbinmax,nactive
 use timestep_sts,   only:sts_get_dtau_next,sts_init_step,use_sts
 use timing,         only:get_timings,setup_timers
 real, intent(in)    :: tzero,time,dtmax,rhomaxnow
 real, intent(inout) :: dt

 nsteps    = 0
 nsteplast = 0
 dtlast    = 0.
 dtinject  = huge(dtinject)
 dtrad     = huge(dtrad)
 np_cs_eq_0 = 0
 np_e_eq_0  = 0
 abortrun_bdy = .false.

 call init_conservation_checks()

 noutput          = 1
 noutput_dtmax    = 1
 ncount_fulldumps = 0
 tprint           = tzero + dtmax
 rhomaxold        = rhomaxnow

 !
 ! Set substepping integration precision depending on the system (default is FSI)
 !
 call set_integration_precision

 if (ind_timesteps) then
    use_global_dt = .false.
    istepfrac     = 0
    tlast         = tzero
    dt            = dtmax/2.**nbinmax  ! use 2.0 here to allow for step too small
    nmovedtot     = 0
    tall          = 0.
    tcheck        = time
    if (ind_timesteps) call reset_time_per_bin() ! initialise bin timers
    call init_step(npart,time,dtmax)
    if (use_sts) then
       call sts_get_dtau_next(dtau,dt,dtmax,dtdiff,nbinmax)
       call sts_init_step(npart,time,dtmax,dtau)  ! overwrite twas for particles requiring super-timestepping
    endif
 else
    use_global_dt = .true.
    nskip   = int(ntot)
    nactive = npart
    istepfrac = 0 ! dummy values
    nbinmax   = 0
    if (dt >= (tprint-time)) dt = tprint-time   ! reach tprint exactly
 endif
!
! threshold for writing to .ev file, to avoid repeatedly computing energies
! for all the particles which would add significantly to the cpu time
!
 nskipped = 0
 nskipped_sink = 0
!
! code timings
!
 call get_timings(twalllast,tcpulast)
 tstart    = twalllast
 tcpustart = tcpulast

 call setup_timers

 call flush(iprint)

end subroutine evol_init

!----------------------------------------------------------------
!+
!  evolve the simulation over all timesteps
!+
!----------------------------------------------------------------
subroutine evol(infile,logfile,evfile,dumpfile,flag)
 use dim,              only:maxvxyzu,mhd,periodic,use_apr,ind_timesteps,driving,inject_parts
 use io,               only:iprint,id,master,iverbose,fatal,warning
 use timestep,         only:time,tmax,dt,dtmax,nmax,nout,nsteps,dtinject,&
                            dtextforce,rhomaxnow,check_dtmax_for_decrease
 use easter_egg,       only:egged,bring_the_egg
 use options,          only:rhofinal1,nfulldump,nmaxdumps,twallmax
 use step_lf_global,   only:step
 use timing,           only:get_timings,print_time,timer,reset_timer,increment_timer,&
                            timers,ntimers,itimer_fromstart,itimer_lastdump,itimer_step,itimer_ev
 use mpiutils,         only:reduce_mpi,reduceall_mpi,barrier_mpi
 use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,&
                            write_binsummary,change_nbinmax,nactive,nactivetot,maxbins,&
                            get_newbin,print_dtind_efficiency,reset_time_per_bin
 use timestep_sts,     only:sts_get_dtau_next,sts_init_step
 use step_lf_global,   only:init_step
 use timestep_sts,     only:use_sts
 use supertimestep,    only:step_sts
 use forcing,          only:correct_bulk_motion
 use centreofmass,     only:correct_bulk_motions
 use inject,           only:inject_particles
 use partinject,       only:update_injected_particles
 use dim,              only:do_radiation
 use options,          only:exchange_radiation_energy,implicit_radiation
 use radiation_utils,  only:update_radenergy
 use apr,              only:update_apr
 use part,             only:npart,npartoftype,nptmass,xyzh,vxyzu,fxyzu,apr_level,&
                            xyzmh_ptmass,vxyz_ptmass,gravity,iboundary,ntot,ibin,iphase,&
                            isionised,rad,radprop,igas
 use quitdump,         only:quit
 use ptmass,           only:icreate_sinks,ipart_createstars
 use io_summary,       only:iosum_nreal,summary_counter,summary_printout,summary_printnow
 use boundary_dyn,     only:dynamic_bdy,update_boundaries
 use HIIRegion,        only:HII_feedback,iH2R,HIIuprate
 integer, optional, intent(in)   :: flag
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: logfile,evfile,dumpfile
 integer         :: nalive,npart_old,nskip,istepHII
 integer(kind=1) :: nbinmaxprev
 integer(kind=8) :: nmovedtot,nalivetot
 real            :: dtnew,tzero,dtmaxold
 real(kind=4)    :: t1,t2,tcpu1,tcpu2
 real(kind=4)    :: twallperdump
 real(kind=4)    :: tall
 logical         :: abortrun,abortrun_bdy,at_dump_time
 logical         :: use_global_dt,do_radiation_update,iexist

 do_radiation_update = do_radiation .and. exchange_radiation_energy .and. .not.implicit_radiation

 tzero = time
 if (.not. initialized) then  ! required because evol is called multiple times in AMUSE... -SR
    call evol_init(tzero,time,dtmax,rhomaxnow,dt)
    initialized = .true.
 endif
!
! --------------------- main loop ----------------------------------------
!
 timestepping: do while ((time < tmax).and.((nsteps < nmax) .or.  (nmax < 0)).and.(rhomaxnow*rhofinal1 < 1.0))

    !
    ! injection of new particles into simulation
    !
    if (inject_parts .and. .not. present(flag)) then
       npart_old = npart
       call inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
       call update_injected_particles(npart_old,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject)
       dtlast = dt
    endif

    if (use_apr) call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level) ! split or merge as required

    dtmaxold = dtmax
    if (ind_timesteps) then
       istepfrac   = istepfrac + 1
       nbinmaxprev = nbinmax
       if (nbinmax > maxbins) call fatal('evolve','timestep too small: try decreasing dtmax?')

       !--determine if dt needs to be decreased; if so, then this will be done
       !  in step the next time it is called;
       !  for global timestepping, this is called in the block where at_dump_time==.true.
       if (istepfrac == 2**nbinmax) then
          twallperdump = reduceall_mpi('max',timers(itimer_lastdump)%wall)
          call check_dtmax_for_decrease(iprint,dtmax,twallperdump,&
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
    endif

    if (gravity .and. icreate_sinks > 0) call ptmass_create_and_update_forces(time,dtextforce)

    if (iH2R > 0 .and. id==master) then
       istepHII = 1
       if (ind_timesteps) then
          istepHII = 2**nbinmax/HIIuprate
          if (istepHII==0) istepHII = 1
       endif
       if (mod(istepfrac,istepHII) == 0 .or. istepfrac == 1 .or. (icreate_sinks == 2 .and. ipart_createstars /= 0)) then
          call HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,isionised)
       endif
    endif

    nsteps = nsteps + 1

    call get_timings(t1,tcpu1)
    !
    ! Strang splitting: implicit update for half step
    !
    if (do_radiation_update) call update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,0.5*dt)
    !
    !--evolve data for one timestep
    !  for individual timesteps this is the shortest timestep
    !
    if ( use_sts ) then
       call step_sts(npart,nactive,time,dt,dtextforce,dtnew,iprint)
    else
       call step(npart,nactive,time,dt,dtextforce,dtnew)
    endif
    !
    ! Strang splitting: implicit update for another half step
    !
    if (do_radiation_update) call update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,0.5*dt)

    !--timings for step call
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_step,t2-t1,tcpu2-tcpu1)
    call summary_counter(iosum_nreal,t2-t1)
!
!--update time = time + dt, and get the dt for the next timestep
!
   at_dump_time = (time >= tmax) &   ! force dump if reached the end of the simulation
                   .or.((nsteps >= nmax).and.(nmax >= 0)).or.(rhomaxnow*rhofinal1 >= 1.0)

   call update_time_and_dt(time,dtmax,dtmaxold,tlast,tcheck,tmax,dt,tall,t2-t1,tcpu2-tcpu1,&
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
    call write_ev_files(ntot,time,dt,at_dump_time)
!
!--write to data file if time is right
!
    if (at_dump_time) then
       call check_and_write_dump(time,tmax,dtmax,dtmaxold,dt,tzero,tprint,rhomaxold,rhomaxnow,nsteps,&
                                 nmax,nout,noutput,noutput_dtmax,ncount_fulldumps,&
                                 dumpfile,infile,evfile,logfile,use_global_dt,abortrun)

       call print_log_and_reset_counters(nsteps,nsteplast,tlast,tprint,tall,&
                                         istepfrac,nalivetot,nmovedtot)

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

 enddo timestepping

end subroutine evol

!----------------------------------------------------------------
!+
!  wrapper routine to update the time and get the dt for the
!  next timestep. Also print a log every step as required.
!+
!----------------------------------------------------------------
subroutine update_time_and_dt(time,dtmax,dtmaxold,tlast,tcheck,tmax,dt,tall,tstep,tcpustep,&
                              istepfrac,nbinmaxprev,ntot,nalivetot,nmovedtot,at_dump_time)
 use dim,          only:ind_timesteps
 use io,           only:id,master,nprocs,iverbose,iprint,warning,fatal
 use mpiutils,     only:bcast_mpi,reduceall_mpi
 use timestep,     only:dtrad,dtforce,dtinject,dtcourant,dterr,print_dtlog,dtmax_ifactor
 use timestep_ind, only:print_dtind_efficiency,update_time_per_bin,print_dtlog_ind,change_nbinmax,&
                        nactivetot,nbinmax
 real,            intent(inout) :: time,tlast,tcheck,tmax,dt
 real,            intent(in)    :: dtmax,dtmaxold
 real(kind=4),    intent(inout) :: tall
 real(kind=4),    intent(in)    :: tstep,tcpustep
 integer,         intent(inout) :: istepfrac
 integer(kind=1), intent(in)    :: nbinmaxprev
 integer(kind=8), intent(in)    :: ntot,nalivetot
 integer(kind=8), intent(inout) :: nmovedtot
 logical,         intent(inout) :: at_dump_time
 real :: dtprint,timecheck
 integer :: inbin

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
!--Determine if this is the correct time to write to the data file
!
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
 use subgroup,    only:group_identify
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
       call group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,nmatrix,&
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
subroutine write_ev_files(ntot,time,dt,at_dump_time)
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
 integer(kind=8), intent(in) :: ntot
 real,            intent(in) :: time,dt
 logical,         intent(in) :: at_dump_time
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
subroutine check_and_write_dump(time,tmax,dtmax,dtmaxold,dt,tzero,tprint,rhomaxold,rhomaxnow,nsteps,&
                                nmax,nout,noutput,noutput_dtmax,ncount_fulldumps,&
                                dumpfile,infile,evfile,logfile,use_global_dt,abortrun)
 use dim,              only:ind_timesteps,inject_parts,driving,idumpfile
 use io,               only:iprint,iwritein,id,master,flush_warnings
 use fileutils,        only:getnextfilename
 use forcing,          only:write_forcingdump
 use options,          only:nfulldump,twallmax,write_files,rkill
 use part,             only:ideadhead,shuffle_part,npart,nptmass,xyzmh_ptmass,accrete_particles_outside_sphere
 use ptmass,           only:calculate_mdot
 use readwrite_infile, only:write_infile
 use readwrite_dumps,  only:write_fulldump,write_smalldump
 use timestep,         only:check_dtmax_for_decrease,check_for_restart_dump,idtmax_n_next,idtmax_frac,idtmax_frac_next,&
                            dtmax_ifactorWT,dtmax_ifactor,idtmax_n
 use timing,           only:timers,itimer_io,itimer_lastdump,itimer_fromstart,increment_timer,get_timings
#ifdef LIVE_ANALYSIS
 use analysis,         only:do_analysis
 use fileutils,        only:numfromfile
 use io,               only:ianalysis
#endif
 real,             intent(in)    :: time,tmax,rhomaxnow
 real,             intent(inout) :: rhomaxold
 real,             intent(inout) :: dtmax,dtmaxold,dt,tzero,tprint
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: dumpfile,evfile,logfile
 integer,          intent(in)    :: nsteps,nmax,nout
 integer,          intent(inout) :: ncount_fulldumps,noutput,noutput_dtmax
 logical,          intent(in)    :: use_global_dt
 logical,          intent(out)   :: abortrun
 logical :: fulldump,writedump
 character(len=120) :: dumpfile_orig
 real(kind=4) :: twallperdump,twallused,t1,t2,tcpu1,tcpu2
 !
 !--Global timesteps: Decrease dtmax if requested (done in step for individual timesteps)
 !
 if (.not. ind_timesteps) then
    twallperdump = timers(itimer_lastdump)%wall
    call check_dtmax_for_decrease(iprint,dtmax,twallperdump,rhomaxold,rhomaxnow,nfulldump,use_global_dt)
    dt = min(dt,dtmax) ! required if decreasing dtmax to ensure that the physically motivated timestep is not too long
 endif

 dumpfile_orig = trim(dumpfile)
 if ((nout <= 0) .or. (mod(noutput,nout)==0)) then
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
    writedump = .true.
 else
    writedump = .false.
 endif

 !--do not dump dead particles into dump files
 if (ideadhead > 0) call shuffle_part(npart)
!
!--get timings since last dump and overall code scaling
!  (get these before writing the dump so we can check whether or not we
!   need to write a full dump based on the wall time;
!   move timer_lastdump outside at_dump_time block so that dtmax can
!   be reduced it too long between dumps)
!
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_fromstart,t2-tstart,tcpu2-tcpustart)

 fulldump = (nout <= 0 .and. mod(noutput,nfulldump)==0) .or. (mod(noutput,nout*nfulldump)==0)
!
!--if max wall time is set (> 1 sec) stop the run at the last full dump
!  that will fit into the walltime constraint, based on the wall time between
!  the last two dumps added to the current total walltime used.  The factor of three for
!  changing to full dumps is to account for the possibility that the next step will take longer.
!  If we are about to write a small dump but it looks like we won't make the next dump,
!  write a full dump instead and stop the run
!
 abortrun = .false.
 if (twallmax > 1.) then
    twallused    = timers(itimer_fromstart)%wall
    twallperdump = timers(itimer_lastdump)%wall
    if (fulldump) then
       if ((twallused + abs(nfulldump)*twallperdump) > twallmax) then
          abortrun = .true.
       endif
    else
       if ((twallused + 3.0*twallperdump) > twallmax) then
          fulldump = .true.
          if (id==master) write(iprint,"(1x,a)") '>> PROMOTING DUMP TO FULL DUMP BASED ON WALL TIME CONSTRAINTS... '
          nfulldump = 1  !  also set all future dumps to be full dumps (otherwise gets confusing)
          if ((twallused + twallperdump) > twallmax) abortrun = .true.
       endif
    endif
 endif
!
!--Promote to full dump if this is the final dump
!
 if ( (time >= tmax) .or. ( (nmax > 0) .and. (nsteps >= nmax) ) ) fulldump = .true.
!
!--flush any buffered warnings to the log file
!
 if (id==master) call flush_warnings()
!
!--write dump file
!
 if (rkill > 0) call accrete_particles_outside_sphere(rkill)
 if (.not.inject_parts) call calculate_mdot(nptmass,time,xyzmh_ptmass)

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
 endif
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_io,t2-t1,tcpu2-tcpu1)

#ifdef LIVE_ANALYSIS
 if (id==master .and. idtmax_frac==0) then
    call do_analysis(dumpfile,numfromfile(dumpfile),xyzh,vxyzu, &
                     massoftype(igas),npart,time,ianalysis)
 endif
#endif

 ! reset counters for when the next dump should be written
 if (idtmax_frac==0) noutput = noutput + 1           ! required to determine frequency of full dumps

 noutput_dtmax = noutput_dtmax + 1     ! required to adjust tprint; will account for varying dtmax
 idtmax_n      = idtmax_n_next
 idtmax_frac   = idtmax_frac_next
 tprint        = tzero + noutput_dtmax*dtmaxold
 dumpfile      = trim(dumpfile_orig)
 if (dtmax_ifactor /= 0) then
    tzero           = tprint - dtmaxold
    tprint          = tzero  + dtmax
    noutput_dtmax   = 1
    dtmax_ifactor   = 0
    dtmax_ifactorWT = 0
 endif

end subroutine check_and_write_dump

!----------------------------------------------------------------
!+
!  routine to print out the log information after a dump
!  and reset istepfrac and other counters
!+
!----------------------------------------------------------------
subroutine print_log_and_reset_counters(nsteps,nsteplast,tlast,tprint,tall,istepfrac,nalivetot,nmovedtot)
 use dim,          only:ind_timesteps
 use io,           only:iprint,id,master,iverbose
 use evwrite,      only:write_evlog
 use part,         only:npart,nptmass,iphase,ibin,xyzh
 use io_summary,   only:summary_printout
 use timing,       only:reduce_timers,timers,itimer_lastdump,reset_timer,ntimers,get_timings
 use timestep,     only:dtmax
 use timestep_ind, only:reset_time_per_bin,print_dtind_efficiency,nbinmax,write_binsummary
 integer,         intent(in)    :: nsteps
 integer,         intent(inout) :: istepfrac,nsteplast
 real,            intent(inout) :: tlast
 real,            intent(in)    :: tprint
 real(kind=4),    intent(inout) :: tall
 integer(kind=8), intent(in)    :: nalivetot
 integer(kind=8), intent(inout) :: nmovedtot
 integer :: i

 call reduce_timers
 if (id==master) then
    call print_timinginfo(iprint,nsteps,nsteplast)
    !--Write out summary to log file
    call summary_printout(iprint,nptmass)
 endif
 nsteplast = nsteps

 if (ind_timesteps) then
    !--print summary of timestep bins
    if (iverbose >= 0) then
       call write_binsummary(npart,nbinmax,dtmax,iphase,ibin,xyzh)
       if (ind_timesteps) call reset_time_per_bin() ! reset bin timers
       if (id==master) call print_dtind_efficiency(iverbose,nalivetot,nmovedtot,tall,timers(itimer_lastdump)%wall,2)
    endif
    tlast = tprint
    istepfrac = 0
    nmovedtot = 0
 endif

 !--print summary of energies and other useful values to the log file
 if (id==master) call write_evlog(iprint)

 call get_timings(twalllast,tcpulast)
 do i = 1,ntimers
    call reset_timer(i)
 enddo

end subroutine print_log_and_reset_counters

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

end module evolve
