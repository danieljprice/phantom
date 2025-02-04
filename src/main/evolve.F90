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
!   checkconserved, dim, energies, evwrite, externalforces, fileutils,
!   forcing, inject, io, io_summary, mf_write, mpiutils, options, part,
!   partinject, ptmass, quitdump, radiation_utils, readwrite_dumps,
!   readwrite_infile, step_lf_global, subgroup, substepping, supertimestep,
!   timestep, timestep_ind, timestep_sts, timing
!
 implicit none
 public :: evol

 private

contains

subroutine evol(infile,logfile,evfile,dumpfile,flag)
 use io,               only:iprint,iwritein,id,master,iverbose,&
                            flush_warnings,nprocs,fatal,warning
 use timestep,         only:time,tmax,dt,dtmax,nmax,nout,nsteps,dtextforce,rhomaxnow,&
                            dtmax_ifactor,dtmax_ifactorWT,dtmax_dratio,check_dtmax_for_decrease,&
                            idtmax_n,idtmax_frac,idtmax_n_next,idtmax_frac_next
 use evwrite,          only:write_evfile,write_evlog
 use energies,         only:etot,totmom,angtot,mdust,np_cs_eq_0,np_e_eq_0,hdivBonB_ave,&
                            hdivBonB_max,mtot
 use checkconserved,   only:etot_in,angtot_in,totmom_in,mdust_in,&
                            init_conservation_checks,check_conservation_error,&
                            check_magnetic_stability,mtot_in
 use dim,              only:maxvxyzu,mhd,periodic,idumpfile,use_apr,ind_timesteps
 use fileutils,        only:getnextfilename
 use options,          only:nfulldump,twallmax,nmaxdumps,rhofinal1,iexternalforce,rkill
 use readwrite_infile, only:write_infile
 use readwrite_dumps,  only:write_smalldump,write_fulldump
 use step_lf_global,   only:step
 use timing,           only:get_timings,print_time,timer,reset_timer,increment_timer,&
                            setup_timers,timers,reduce_timers,ntimers,&
                            itimer_fromstart,itimer_lastdump,itimer_step,itimer_io,itimer_ev
 use mpiutils,         only:reduce_mpi,reduceall_mpi,barrier_mpi,bcast_mpi
#ifdef IND_TIMESTEPS
 use part,             only:ibin,iphase
 use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,update_time_per_bin,&
                            write_binsummary,change_nbinmax,nactive,nactivetot,maxbins,&
                            print_dtlog_ind,get_newbin,print_dtind_efficiency
 use timestep,         only:dtdiff
 use timestep_sts,     only:sts_get_dtau_next,sts_init_step
 use step_lf_global,   only:init_step
#else
 use timestep,         only:dtforce,dtcourant,dterr,print_dtlog
#endif
 use timestep_sts,     only:use_sts
 use supertimestep,    only:step_sts
#ifdef DRIVING
 use forcing,          only:write_forcingdump
#endif
#ifdef CORRECT_BULK_MOTION
 use centreofmass,     only:correct_bulk_motion
#endif
 use part,             only:ideadhead,shuffle_part
#ifdef INJECT_PARTICLES
 use inject,           only:inject_particles
 use part,             only:npartoftype
 use partinject,       only:update_injected_particles
#endif
 use dim,              only:do_radiation
 use options,          only:exchange_radiation_energy,implicit_radiation
 use part,             only:rad,radprop,igas
 use radiation_utils,  only:update_radenergy
 use timestep,         only:dtrad
#ifdef LIVE_ANALYSIS
 use analysis,         only:do_analysis
 use part,             only:igas
 use fileutils,        only:numfromfile
 use io,               only:ianalysis
#endif
 use apr,              only:update_apr,create_or_update_apr_clump
 use part,             only:npart,nptmass,xyzh,vxyzu,fxyzu,fext,divcurlv,massoftype, &
                            xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dptmass,gravity,iboundary, &
                            fxyz_ptmass_sinksink,ntot,poten,ndustsmall,&
                            accrete_particles_outside_sphere,apr_level,aprmassoftype,&
                            sf_ptmass,isionised,dsdt_ptmass,isdead_or_accreted
 use part,             only:n_group,n_ingroup,n_sing,group_info,bin_info,nmatrix
 use quitdump,         only:quit
 use ptmass,           only:icreate_sinks,ptmass_create,ipart_rhomax,pt_write_sinkev,calculate_mdot, &
                            set_integration_precision,ptmass_create_stars,use_regnbody,ptmass_create_seeds,&
                            ipart_createseeds,ipart_createstars
 use io_summary,       only:iosum_nreal,summary_counter,summary_printout,summary_printnow
 use externalforces,   only:iext_spiral
 use boundary_dyn,     only:dynamic_bdy,update_boundaries
 use HIIRegion,        only:HII_feedback,iH2R,HIIuprate
 use subgroup,         only:group_identify
 use substepping,      only:get_force
#ifdef MFLOW
 use mf_write,         only:mflow_write
#endif
#ifdef VMFLOW
 use mf_write,         only:vmflow_write
#endif
#ifdef BINPOS
 use mf_write,         only:binpos_write
#endif

 integer, optional, intent(in)   :: flag
 character(len=*), intent(in)    :: infile
 character(len=*), intent(inout) :: logfile,evfile,dumpfile
 integer         :: i,noutput,noutput_dtmax,nsteplast,ncount_fulldumps
 real            :: dtnew,dtlast,timecheck,rhomaxold,dtmax_log_dratio
 real            :: tprint,tzero,dtmaxold,dtinject
 real(kind=4)    :: t1,t2,tcpu1,tcpu2,tstart,tcpustart
 real(kind=4)    :: twalllast,tcpulast,twallperdump,twallused
#ifdef IND_TIMESTEPS
 integer         :: nalive,inbin
 integer(kind=1) :: nbinmaxprev
 integer(kind=8) :: nmovedtot,nalivetot
 real            :: tlast,tcheck,dtau
 real(kind=4)    :: tall
 real(kind=4)    :: timeperbin(0:maxbins)
 logical         :: dt_changed
#else
 real            :: dtprint
 integer         :: nactive,istepfrac
 integer(kind=1) :: nbinmax
 logical, parameter :: dt_changed = .false.
#endif
#ifdef INJECT_PARTICLES
 integer         :: npart_old
#endif
 logical         :: fulldump,abortrun,abortrun_bdy,at_dump_time,writedump
 logical         :: should_conserve_energy,should_conserve_momentum,should_conserve_angmom
 logical         :: should_conserve_dustmass,should_conserve_aprmass
 logical         :: use_global_dt
 integer         :: j,nskip,nskipped,nevwrite_threshold,nskipped_sink,nsinkwrite_threshold
 character(len=120) :: dumpfile_orig
 integer         :: dummy,istepHII,nptmass_old

 dummy = 0

 tprint    = 0.
 nsteps    = 0
 nsteplast = 0
 tzero     = time
 dtlast    = 0.
 dtinject  = huge(dtinject)
 dtrad     = huge(dtrad)
 np_cs_eq_0 = 0
 np_e_eq_0  = 0
 abortrun_bdy = .false.
 dumpfile_orig = trim(dumpfile)

 call init_conservation_checks(should_conserve_energy,should_conserve_momentum,&
                               should_conserve_angmom,should_conserve_dustmass,&
                               should_conserve_aprmass)

 noutput          = 1
 noutput_dtmax    = 1
 ncount_fulldumps = 0
 tprint           = tzero + dtmax
 rhomaxold        = rhomaxnow
 if (dtmax_dratio > 0.) then
    dtmax_log_dratio = log10(dtmax_dratio)
 else
    dtmax_log_dratio = 0.0
 endif

 !
 ! Set substepping integration precision depending on the system (default is FSI)
 !
 call set_integration_precision

#ifdef IND_TIMESTEPS
 use_global_dt = .false.
 istepfrac     = 0
 tlast         = tzero
 dt            = dtmax/2.**nbinmax  ! use 2.0 here to allow for step too small
 nmovedtot     = 0
 tall          = 0.
 tcheck        = time
 timeperbin(:) = 0.
 dt_changed    = .false.
 call init_step(npart,time,dtmax)
 if (use_sts) then
    call sts_get_dtau_next(dtau,dt,dtmax,dtdiff,nbinmax)
    call sts_init_step(npart,time,dtmax,dtau)  ! overwrite twas for particles requiring super-timestepping
 endif
#else
 use_global_dt = .true.
 nskip   = int(ntot)
 nactive = npart
 istepfrac = 0 ! dummy values
 nbinmax   = 0
 if (dt >= (tprint-time)) dt = tprint-time   ! reach tprint exactly
#endif
!
! threshold for writing to .ev file, to avoid repeatedly computing energies
! for all the particles which would add significantly to the cpu time
!

 nskipped = 0
 if (iexternalforce==iext_spiral) then
    nevwrite_threshold = int(4.99*ntot) ! every 5 full steps
 else
    nevwrite_threshold = int(1.99*ntot) ! every 2 full steps
 endif
 nskipped_sink = 0
 nsinkwrite_threshold  = int(0.99*ntot)
!
! code timings
!
 call get_timings(twalllast,tcpulast)
 tstart    = twalllast
 tcpustart = tcpulast

 call setup_timers

 call flush(iprint)

!
! --------------------- main loop ----------------------------------------
!
 timestepping: do while ((time < tmax).and.((nsteps < nmax) .or.  (nmax < 0)).and.(rhomaxnow*rhofinal1 < 1.0))

#ifdef INJECT_PARTICLES
    !
    ! injection of new particles into simulation
    !
    if (.not. present(flag)) then
       npart_old=npart
       call inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
       call update_injected_particles(npart_old,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject)
    endif
#endif

    if (use_apr) then
       ! split or merge as required
       call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
    endif

    dtmaxold    = dtmax
#ifdef IND_TIMESTEPS
    istepfrac   = istepfrac + 1
    nbinmaxprev = nbinmax
    if (nbinmax > maxbins) call fatal('evolve','timestep too small: try decreasing dtmax?')

    !--determine if dt needs to be decreased; if so, then this will be done
    !  in step the next time it is called;
    !  for global timestepping, this is called in the block where at_dump_time==.true.
    if (istepfrac==2**nbinmax) then
       twallperdump = reduceall_mpi('max', timers(itimer_lastdump)%wall)
       call check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_log_dratio,&
                                     rhomaxold,rhomaxnow,nfulldump,use_global_dt)
    endif

    !--sanity check on istepfrac...
    if (istepfrac > 2**nbinmax) then
       write(iprint,*) 'ERROR: istepfrac = ',istepfrac,' / ',2**nbinmax
       call fatal('evolve','error in individual timesteps')
    endif

    !--flag particles as active or not for this timestep
    call set_active_particles(npart,nactive,nalive,iphase,ibin,xyzh)
    nactivetot = reduceall_mpi('+', nactive)
    nalivetot = reduceall_mpi('+', nalive)
    nskip = int(nactivetot)

    !--print summary of timestep bins
    if (iverbose >= 2) call write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)

    !--Implement dynamic boundaries (for individual-timestepping) once per dump
    if (dynamic_bdy .and. nactive==nalive .and. istepfrac==2**nbinmax) then
       call update_boundaries(nactive,nalive,npart,abortrun_bdy)
    endif
#else
    !--If not using individual timestepping, set nskip to the total number of particles
    !  across all nodes
    nskip = int(ntot)
#endif
    nptmass_old = nptmass
    if (gravity .and. icreate_sinks > 0 .and. ipart_rhomax /= 0) then
       !
       ! creation of new sink particles
       !
       if (use_apr) then
          call create_or_update_apr_clump(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,aprmassoftype)
       else
          call ptmass_create(nptmass,npart,ipart_rhomax,xyzh,vxyzu,fxyzu,fext,divcurlv,&
                          poten,massoftype,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink,sf_ptmass,dptmass,time)
       endif
    endif

    if (icreate_sinks == 2) then
       !
       ! creation of new seeds into evolved sinks
       !
       if (ipart_createseeds /= 0) then
          call ptmass_create_seeds(nptmass,ipart_createseeds,sf_ptmass,time)
       endif
       !
       ! creation of new stars from sinks (cores)
       !
       if (ipart_createstars /= 0) then
          call ptmass_create_stars(nptmass,ipart_createstars,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink, &
                               sf_ptmass,time)
       endif
    endif

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

    ! Need to recompute the force when sink or stars are created
    if (nptmass > nptmass_old .or. ipart_createseeds /= 0 .or. ipart_createstars /= 0) then
       if (use_regnbody) then
          call group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,nmatrix,&
                              new_ptmass=.true.,dtext=dtextforce)
       endif
       call get_force(nptmass,npart,0,1,time,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass,vxyz_ptmass,&
                      fxyz_ptmass,dsdt_ptmass,0.,0.,dummy,.false.,sf_ptmass,bin_info,group_info,&
                      nmatrix)
       if (ipart_createseeds /= 0) ipart_createseeds = 0 ! reset pointer to zero
       if (ipart_createstars /= 0) ipart_createstars = 0 ! reset pointer to zero
       dummy = 0
    endif
    !
    ! Strang splitting: implicit update for half step
    !
    if (do_radiation  .and. exchange_radiation_energy  .and. .not.implicit_radiation) then
       call update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,0.5*dt)
    endif

    nsteps = nsteps + 1
!
!--evolve data for one timestep
!  for individual timesteps this is the shortest timestep
!
    call get_timings(t1,tcpu1)
    if ( use_sts ) then
       call step_sts(npart,nactive,time,dt,dtextforce,dtnew,iprint)
    else
       call step(npart,nactive,time,dt,dtextforce,dtnew)
    endif
    !
    ! Strang splitting: implicit update for another half step
    !
    if (do_radiation .and. exchange_radiation_energy .and. .not.implicit_radiation) then
       call update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,0.5*dt)
    endif

    dtlast = dt

    !--timings for step call
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_step,t2-t1,tcpu2-tcpu1)
    call summary_counter(iosum_nreal,t2-t1)

#ifdef IND_TIMESTEPS
    tcheck = tcheck + dt

    !--update time in way that is free of round-off errors
    time = tlast + istepfrac/real(2**nbinmaxprev)*dtmaxold

    !--print efficiency of partial timestep
    if (id==master) call print_dtind_efficiency(iverbose,nalivetot,nactivetot,tall,t2-t1,1)

    call update_time_per_bin(tcpu2-tcpu1,istepfrac,nbinmaxprev,timeperbin,inbin)
    nmovedtot = nmovedtot + nactivetot

    !--check that time is as it should be, may indicate error in individual timestep routines
    if (abs(tcheck-time) > 1.e-4) call warning('evolve','time out of sync',var='error',val=abs(tcheck-time))

    if (id==master .and. (iverbose >= 1 .or. inbin <= 3)) &
       call print_dtlog_ind(iprint,istepfrac,2**nbinmaxprev,time,dt,nactivetot,tcpu2-tcpu1,ntot)

    !--if total number of bins has changed, adjust istepfrac and dt accordingly
    !  (ie., decrease or increase the timestep)
    if (nbinmax /= nbinmaxprev .or. dtmax_ifactor /= 0) then
       call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
       dt_changed = .true.
    endif

#else

    ! advance time on master thread only
    if (id == master) time = time + dt
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
    if (id==master) call print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,dtrad,dtprint,dtinject,ntot)
#endif

!   check that MPI threads are synchronised in time
    timecheck = reduceall_mpi('+',time)
    if (abs(timecheck/nprocs - time) > 1.e-13) then
       call fatal('evolve','time differs between MPI threads',var='time',val=timecheck/nprocs)
    endif
!
!--Update timer from last dump to see if dtmax needs to be reduced
!
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_lastdump,t2-t1,tcpu2-tcpu1)
!
!--Determine if this is the correct time to write to the data file
!
    at_dump_time = (time >= tmax) &
                   .or.((nsteps >= nmax).and.(nmax >= 0)).or.(rhomaxnow*rhofinal1 >= 1.0)
#ifdef IND_TIMESTEPS
    if (istepfrac==2**nbinmax) at_dump_time = .true.
#else
    if (time >= tprint) at_dump_time = .true.
#endif
!
!--Calculate total energy etc and write to ev file
!  For individual timesteps, we do not want to do this every step, but we want
!  to do this as often as possible without a performance hit. The criteria
!  here is that it is done once >10% of particles (cumulatively) have been evolved.
!  That is, either >10% are being stepped, or e.g. 1% have moved 10 steps.
!  Perform this prior to writing the dump files so that diagnostic values calculated
!  in energies can be correctly included in the dumpfiles
!
    nskipped = nskipped + nskip
    if (nskipped >= nevwrite_threshold .or. at_dump_time .or. dt_changed .or. iverbose==5) then
       nskipped = 0
       call get_timings(t1,tcpu1)
       call write_evfile(time,dt)
       if (should_conserve_momentum) call check_conservation_error(totmom,totmom_in,1.e-1,'linear momentum')
       if (should_conserve_angmom)   call check_conservation_error(angtot,angtot_in,1.e-1,'angular momentum')
       if (should_conserve_energy)   call check_conservation_error(etot,etot_in,1.e-1,'energy')
       if (should_conserve_dustmass) then
          do j = 1,ndustsmall
             call check_conservation_error(mdust(j),mdust_in(j),1.e-1,'dust mass',decrease=.true.)
          enddo
       endif
       if (mhd) call check_magnetic_stability(hdivBonB_ave,hdivBonB_max)
       if (should_conserve_aprmass) call check_conservation_error(mtot,mtot_in,massoftype(igas),'total mass')
       if (id==master) then
          if (np_e_eq_0  > 0) call warning('evolve','N gas particles with energy = 0',var='N',ival=int(np_e_eq_0,kind=4))
          if (np_cs_eq_0 > 0) call warning('evolve','N gas particles with sound speed = 0',var='N',ival=int(np_cs_eq_0,kind=4))
       endif

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
!-- Print out the sink particle properties & reset dt_changed.
!-- Added total force on sink particles and sink-sink forces to write statement (fxyz_ptmass,fxyz_ptmass_sinksink)
    nskipped_sink = nskipped_sink + nskip
    if (nskipped_sink >= nsinkwrite_threshold .or. at_dump_time .or. dt_changed) then
       nskipped_sink = 0
       call pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink)
#ifdef IND_TIMESTEPS
       dt_changed = .false.
#endif
    endif
!
!--write to data file if time is right
!
    if (at_dump_time) then
#ifndef IND_TIMESTEPS
!
!--Global timesteps: Decrease dtmax if requested (done in step for individual timesteps)
       twallperdump = timers(itimer_lastdump)%wall
       call check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_log_dratio,&
                                     rhomaxold,rhomaxnow,nfulldump,use_global_dt)
       dt = min(dt,dtmax) ! required if decreasing dtmax to ensure that the physically motivated timestep is not too long
#endif

       !--modify evfile and logfile names with new number
       if ((nout <= 0) .or. (mod(noutput,nout)==0)) then
          if (noutput==1) then
             evfile  = getnextfilename(evfile)
             logfile = getnextfilename(logfile)
          endif
!         Update values for restart dumps
          if (dtmax_ifactorWT ==0) then
             idtmax_n_next    =  idtmax_n
             idtmax_frac_next =  idtmax_frac
          elseif (dtmax_ifactorWT > 0) then
             idtmax_n_next    =  idtmax_n   *dtmax_ifactorWT
             idtmax_frac_next =  idtmax_frac*dtmax_ifactorWT
          elseif (dtmax_ifactorWT < 0) then
             idtmax_n_next    = -idtmax_n   /dtmax_ifactorWT
             idtmax_frac_next = -idtmax_frac/dtmax_ifactorWT
          endif
          idtmax_frac_next = idtmax_frac_next + 1
          idtmax_frac_next = mod(idtmax_frac_next,idtmax_n_next)
          dumpfile_orig = trim(dumpfile)
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
#ifndef INJECT_PARTICLES
       call calculate_mdot(nptmass,time,xyzmh_ptmass)
#endif
       call get_timings(t1,tcpu1)
       if (writedump) then
          if (fulldump) then
             call write_fulldump(time,dumpfile)
             if (id==master) then
                call write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
#ifdef DRIVING
                call write_forcingdump(time,dumpfile)
#endif
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
       call reduce_timers
       if (id==master) then
          call print_timinginfo(iprint,nsteps,nsteplast)
          !--Write out summary to log file
          call summary_printout(iprint,nptmass)
       endif
#ifdef IND_TIMESTEPS
       !--print summary of timestep bins
       if (iverbose >= 0) then
          call write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)
          timeperbin(:) = 0.
          if (id==master) call print_dtind_efficiency(iverbose,nalivetot,nmovedtot,tall,timers(itimer_lastdump)%wall,2)
       endif
       tlast = tprint
       istepfrac = 0
       nmovedtot = 0
#endif
       !--print summary of energies and other useful values to the log file
       if (id==master) call write_evlog(iprint)

#ifndef IND_TIMESTEPS
       !--Implement dynamic boundaries (for global timestepping)
       if (dynamic_bdy) call update_boundaries(nactive,nactive,npart,abortrun_bdy)
#endif

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

       if (abortrun_bdy) then
          write(iprint,"(1x,a)") 'Will likely surpass maxp_hard next time we need to add particles.'
          write(iprint,"(1x,a)") 'Recompile with larger maxp_hard.'
          write(iprint,"(1x,a)") '>> ABORTING... '
          return
       endif

       if (nmaxdumps > 0 .and. ncount_fulldumps >= nmaxdumps) then
          if (id==master) write(iprint,"(a)") '>> reached maximum number of full dumps as specified in input file, stopping...'
          return
       endif

       twalllast = t2
       tcpulast = tcpu2
       do i = 1,ntimers
          call reset_timer(i)
       enddo

       if (idtmax_frac==0) then
          noutput    = noutput + 1           ! required to determine frequency of full dumps
       endif
       noutput_dtmax = noutput_dtmax + 1     ! required to adjust tprint; will account for varying dtmax
       idtmax_n      = idtmax_n_next
       idtmax_frac   = idtmax_frac_next
       tprint        = tzero + noutput_dtmax*dtmaxold
       nsteplast     = nsteps
       dumpfile      = trim(dumpfile_orig)
       if (dtmax_ifactor/=0) then
          tzero           = tprint - dtmaxold
          tprint          = tzero  + dtmax
          noutput_dtmax   = 1
          dtmax_ifactor   = 0
          dtmax_ifactorWT = 0
       endif
    endif

#ifdef CORRECT_BULK_MOTION
    call correct_bulk_motion()
#endif

    if (iverbose >= 1 .and. id==master) write(iprint,*)
    call flush(iprint)
    !--Write out log file prematurely (if requested based upon nstep, walltime)
    if ( summary_printnow() ) call summary_printout(iprint,nptmass)

 enddo timestepping

end subroutine evol

!----------------------------------------------------------------
!+
!  routine to print out the timing information at each full dump
!+
!----------------------------------------------------------------
subroutine print_timinginfo(iprint,nsteps,nsteplast)
 use io,     only:formatreal
 use timing, only:timer,timers,print_timer,itimer_fromstart,itimer_lastdump,&
                  itimer_step,itimer_link,itimer_balance,itimer_dens,&
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
    write(iprint,"(1x,a,f6.2,a)") 'WARNING: ',100.*(1.-fracinstep),'% of time was in unusual routines (not dens/force/link)'
 endif

end subroutine print_timinginfo

end module evolve
