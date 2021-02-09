! AMUSE interface library for Phantom
! (c) 2019 - 2021 Steven Rieder

!
! Initialize Phantom and set default parameters
!
subroutine amuse_initialize_code()
    use dim, only:maxp,maxp_hard,maxvxyzu
    use memory, only:allocate_memory
    use units, only:set_units,utime,umass,udist,unit_density
    use physcon, only:gram,seconds,solarm,pc
    use timestep, only:dtmax,dtextforce
    use initial, only:initialise
#ifdef IND_TIMESTEPS
    use timestep_ind,     only:istepfrac,ibinnow
    use part,             only:ibin,ibin_old,ibin_wake
#else
    use timestep,         only:dtcourant,dtforce
#endif
    implicit none
    call allocate_memory(maxp_hard)
    call amuse_set_defaults()
    call amuse_set_polyk(0.)
    !print*, "maxvxyzu: ", maxvxyzu
    !maxvxyzu = 4
    !call set_units(dist=50.*pc,mass=4600.*solarm,G=1.)
    !call set_units(dist=1.d20,mass=1.d40,G=1.)
    !call set_units(dist=1.*pc,mass=1.*solarm,G=1.)
    ! call set_units(dist=1.,mass=1.,time=1.)
    call initialise()

    umass = 1.98892d33 * gram  ! AMUSE MSun
    utime = 60 * 60 * 24 * 365.25 * 1d6 * seconds  ! 10 Julian Kyr
    !call set_units(time=utime,mass=umass,G=1.)
    call set_units(dist=0.1*pc,mass=1.*solarm,G=1.)
    !print*, "utime/mass/dist/udens: ", utime, umass, udist, unit_density
    !dtmax = 0.01 * utime
    dtextforce = 1.e-6

#ifdef IND_TIMESTEPS
    ibin(:)       = 0
    ibin_old(:)   = 0
    ibin_wake(:)  = 0
    istepfrac     = 0
    ibinnow       = 0
#else
    dtcourant = huge(dtcourant)
    dtforce   = huge(dtforce)
#endif

end subroutine amuse_initialize_code

subroutine amuse_commit_particles()
    !use eos, only:ieos,init_eos
    use deriv, only:derivs
    use part, only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
            dustprop,ddustprop,dustfrac,ddustevol,dens,drad,radprop,dustevol,&
            eos_vars,metrics,pxyzu,rad
    use timestep, only:time,dtmax
    use units, only:udist,utime,umass
    !use timestep, only:dtmax
    implicit none
    integer :: ierr
    integer :: nactive
    double precision :: dt, dtnew_first

    dtnew_first = dtmax
    dt = 0
    nactive = npart
    !double precision :: dtmax_tmp
    !dtmax_tmp = dtmax
    !call code_init()
    !call initialise()
    !call init_eos(ieos,ierr)
    call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
            Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
            dustevol,ddustevol,dustfrac,eos_vars,time,dt,dtnew_first,pxyzu,dens,metrics)
    ! call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
    !         Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtnew_first)
    !dtmax = dtmax_tmp !code_init sets dtmax to 1 and we can't have that.
    print*, "COMMIT_PARTICLES, calling amuse_init_evol"
    call amuse_init_evol()
    !*print*, "udist utime umass:", udist, utime, umass
end subroutine amuse_commit_particles

subroutine amuse_recommit_particles()
    use deriv, only:derivs
    use part, only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
            dustprop,ddustprop,dustfrac,ddustevol,dens,drad,radprop,dustevol,&
            eos_vars,metrics,pxyzu,rad
    use timestep, only:time,dtmax
    implicit none
    integer :: ierr
    integer :: nactive
    double precision :: dt, dtnew_first
    dtnew_first = dtmax
    dt = 0
    nactive = npart
    call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
            Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
            dustevol,ddustevol,dustfrac,eos_vars,time,dt,dtnew_first,pxyzu,dens,metrics)
    print*, "RECOMMIT_PARTICLES"
end subroutine amuse_recommit_particles

subroutine amuse_cleanup_code()
    use initial, only:endrun
    implicit none

    call endrun()
end subroutine amuse_cleanup_code

! Get npart
subroutine amuse_get_npart(npart_out, nodisabled)
    use part, only:npart,xyzh
    implicit none
    integer, intent(out) :: npart_out
    logical, intent(in)  :: nodisabled

    integer :: i

    if (nodisabled) then
        npart_out = 0
        do i=1,npart
            if (xyzh(4,i)  >  0.) then
                npart_out = npart_out + 1
            endif
        enddo
    else
        npart_out = npart
    endif
end subroutine amuse_get_npart

! Set default parameters
subroutine amuse_set_defaults()
    use options, only: set_default_options,iexternalforce
    implicit none

   call set_default_options()
end subroutine amuse_set_defaults

! This initialises things. This really should only be called once, before the first step.
subroutine amuse_init_evol()
 use io,               only:iprint,iwritein,id,master,iverbose,flush_warnings,nprocs,fatal,warning
 use timestep,         only:time,tmax,dt,dtmax,nmax,nout,nsteps,rhomaxnow,&
                            dtmax_ifactor,dtmax_dratio,check_dtmax_for_decrease
 use energies,         only:etot,totmom,angtot,mdust,np_cs_eq_0,np_e_eq_0
 use dim,              only:maxvxyzu,mhd,periodic
 use fileutils,        only:getnextfilename
 use options,          only:nfulldump,twallmax,nmaxdumps,rhofinal1,use_dustfrac,iexternalforce,&
                            icooling,ieos,ipdv_heating,ishock_heating,iresistive_heating
 use step_lf_global,   only:step
 use timing,           only:get_timings,print_time,timer,reset_timer,increment_timer
 use mpiutils,         only:reduce_mpi,reduceall_mpi,barrier_mpi,bcast_mpi
#ifdef SORT
 use sort_particles,   only:sort_part
#endif
#ifdef IND_TIMESTEPS
 use dim,              only:maxp
 use part,             only:maxphase,ibin,iphase
 use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,update_time_per_bin,&
                            write_binsummary,change_nbinmax,nactive,nactivetot,maxbins,&
                            print_dtlog_ind,get_newbin
 use timestep,         only:dtdiff
 use timestep_sts,     only:sts_get_dtau_next,sts_init_step
 use step_lf_global,   only:init_step
#else
 use timestep,         only:dtforce,dtcourant,dterr,print_dtlog
#endif
 use timestep_sts,     only: use_sts
 use supertimestep,    only: step_sts
#ifdef DRIVING
 use forcing,          only:write_forcingdump
#endif
#ifdef CORRECT_BULK_MOTION
 use centreofmass,     only:correct_bulk_motion
#endif
#ifdef MPI
 use part,             only:ideadhead,shuffle_part
#endif
#ifdef IND_TIMESTEPS
 use part,             only:twas
 use timestep_ind,     only:get_dt
#endif
 use part,             only:npart,nptmass,xyzh,vxyzu,fxyzu,fext,divcurlv,massoftype, &
                            xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gravity,iboundary,npartoftype, &
                            fxyz_ptmass_sinksink,ntot,poten,ndustsmall
 use quitdump,         only:quit
 use ptmass,           only:icreate_sinks,ptmass_create,ipart_rhomax,pt_write_sinkev
 use io_summary,       only:iosum_nreal,summary_counter,summary_printout,summary_printnow
 use externalforces,   only:iext_spiral
#ifdef MFLOW
 use mf_write,         only:mflow_write
#endif
#ifdef VMFLOW
 use mf_write,         only:vmflow_write
#endif
#ifdef BINPOS
 use mf_write,         only:binpos_write
#endif

 real            :: dtnew,dtlast,timecheck,rhomaxold,dtmax_log_dratio
 real            :: tprint,tzero,dtmaxold,dtinject
 real(kind=4)    :: t1,t2,tcpu1,tcpu2,tstart,tcpustart
 real(kind=4)    :: twalllast,tcpulast,twallperdump,twallused
#ifdef IND_TIMESTEPS
 integer         :: i,nalive,inbin,iamtypei
 integer(kind=1) :: nbinmaxprev
 integer(kind=8) :: nmovedtot,nalivetot
 real            :: tlast,fracactive,speedup,tcheck,dtau,efficiency
 real(kind=4)    :: tall
 real(kind=4)    :: timeperbin(0:maxbins)
 logical         :: dt_changed
 integer         :: iloop,npart_old
#else
 integer         :: nactive
 logical, parameter :: dt_changed = .false.
#endif
 logical         :: fulldump,abortrun,at_dump_time
 logical         :: use_global_dt
 integer         :: j,nskip,nskipped,nevwrite_threshold,nskipped_sink,nsinkwrite_threshold
 type(timer)     :: timer_fromstart,timer_lastdump,timer_step,timer_ev,timer_io

 nsteps    = 0
 tzero     = time
 dtlast    = 0.
 dtinject  = huge(dtinject)
 np_cs_eq_0 = 0
 np_e_eq_0  = 0

 rhomaxold        = rhomaxnow
 if (dtmax_dratio > 0.) then
    dtmax_log_dratio = log10(dtmax_dratio)
 else
    dtmax_log_dratio = 0.0
 endif

#ifdef IND_TIMESTEPS
 use_global_dt = .false.
 istepfrac     = 0
 tlast         = tzero
 dt            = dtmax/2**nbinmax
 nmovedtot     = 0
 tall          = 0.
 tcheck        = time
 timeperbin(:) = 0.
 dt_changed    = .false.
 nbinmax       = 0
!
! first time through, move all particles on shortest timestep
! then allow them to gradually adjust levels.
! Keep boundary particles on level 0 since forces are never calculated
! and to prevent boundaries from limiting the timestep
!
 if (time < tiny(time)) then
    !$omp parallel do schedule(static) private(i,iamtypei)
    do i=1,npart
       ibin(i) = nbinmax
       if (maxphase==maxp) then
          if (abs(iphase(i))==iboundary) ibin(i) = 0
       endif
    enddo
 endif
 !call init_step(npart,time,dtmax)
 !if (use_sts) then
 !   call sts_get_dtau_next(dtau,dt,dtmax,dtdiff,nbinmax)
 !   call sts_init_step(npart,time,dtmax,dtau)  ! overwrite twas for particles requiring super-timestepping
 !endif
#else
 use_global_dt = .true.
 nskip   = npart
 nactive = npart
#endif

 npart_old = 0
end subroutine amuse_init_evol

subroutine amuse_new_step(tlast)
 use io,               only:iprint,iwritein,id,master,iverbose,flush_warnings,nprocs,fatal,warning
 use timestep,         only:time,tmax,dt,dtmax,nmax,nout,nsteps,dtextforce,rhomaxnow,&
                            dtmax_ifactor,dtmax_dratio,check_dtmax_for_decrease
 use energies,         only:etot,totmom,angtot,mdust,np_cs_eq_0,np_e_eq_0
 use dim,              only:maxvxyzu,mhd,periodic
 use fileutils,        only:getnextfilename
 use options,          only:nfulldump,twallmax,nmaxdumps,rhofinal1,use_dustfrac,iexternalforce,&
                            icooling,ieos,ipdv_heating,ishock_heating,iresistive_heating
 use step_lf_global,   only:step
 use timing,           only:get_timings,print_time,timer,reset_timer,increment_timer
 use mpiutils,         only:reduce_mpi,reduceall_mpi,barrier_mpi,bcast_mpi
#ifdef SORT
 use sort_particles,   only:sort_part
#endif
#ifdef IND_TIMESTEPS
 use dim,              only:maxp
 use part,             only:maxphase,ibin,iphase
 use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,update_time_per_bin,&
                            write_binsummary,change_nbinmax,nactive,nactivetot,maxbins,&
                            print_dtlog_ind,get_newbin
 use timestep,         only:dtdiff,C_cool
 use timestep_sts,     only:sts_get_dtau_next,sts_init_step
 use step_lf_global,   only:init_step
#else
 use timestep,         only:dtforce,dtcourant,dterr,print_dtlog
#endif
 use timestep_sts,     only: use_sts
 use supertimestep,    only: step_sts
#ifdef DRIVING
 use forcing,          only:write_forcingdump
#endif
#ifdef CORRECT_BULK_MOTION
 use centreofmass,     only:correct_bulk_motion
#endif
#ifdef MPI
 use part,             only:ideadhead,shuffle_part
#endif
#ifdef IND_TIMESTEPS
 use part,             only:twas
 use timestep_ind,     only:get_dt
#endif
 use part,             only:npart,nptmass,xyzh,vxyzu,fxyzu,fext,divcurlv,massoftype, &
                            xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gravity,iboundary,npartoftype, &
                            fxyz_ptmass_sinksink,ntot,poten,ndustsmall
 use quitdump,         only:quit
 use ptmass,           only:icreate_sinks,ptmass_create,ipart_rhomax,pt_write_sinkev
 use io_summary,       only:iosum_nreal,summary_counter,summary_printout,summary_printnow
 use externalforces,   only:iext_spiral
#ifdef MFLOW
 use mf_write,         only:mflow_write
#endif
#ifdef VMFLOW
 use mf_write,         only:vmflow_write
#endif
#ifdef BINPOS
 use mf_write,         only:binpos_write
#endif

 real            :: dtnew,dtlast,timecheck,rhomaxold,dtmax_log_dratio
 real            :: tprint,dtmaxold,dtinject
 real(kind=4)    :: t1,t2,tcpu1,tcpu2,tstart,tcpustart
 real(kind=4)    :: twalllast,tcpulast,twallperdump,twallused
#ifdef IND_TIMESTEPS
 integer         :: i,nalive,inbin,iamtypei
 integer(kind=1) :: nbinmaxprev
 integer(kind=8) :: nmovedtot,nalivetot
 real            :: fracactive,speedup,tcheck,dtau,efficiency,tbegin
 real, intent(in) :: tlast
 real(kind=4)    :: tall
 real(kind=4)    :: timeperbin(0:maxbins)
 logical         :: dt_changed
 integer         :: iloop,npart_old
#else
 real            :: dtprint
 integer         :: nactive
 logical, parameter :: dt_changed = .false.
#endif
 logical         :: use_global_dt
 integer         :: j,nskip,nskipped,nevwrite_threshold,nskipped_sink,nsinkwrite_threshold
 type(timer)     :: timer_fromstart,timer_lastdump,timer_step,timer_ev,timer_io
!
! --------------------- main loop ----------------------------------------
!
 tbegin = time
 tcheck = time
 npart_old = npart
 !timestepping: do while ((time < tmax).and.((nsteps < nmax) .or.  (nmax < 0)).and.(rhomaxnow*rhofinal1 < 1.0))

 if (istepfrac==0) then
    twallperdump = reduceall_mpi('max', timer_lastdump%wall)
    call check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_ifactor,dtmax_log_dratio,&
                                  rhomaxold,rhomaxnow,nfulldump,use_global_dt)
 endif
 if (dt > dtmax) dt = dtmax
 print*, "time, tcheck, dt, dtmax, tlast+dtmax, npart_old, npart: ", time, tcheck, dt, dtmax, (tlast+dtmax), npart_old, npart

 timesubstepping: do while (istepfrac < 2**nbinmax)
    dtmaxold    = dtmax
#ifdef IND_TIMESTEPS
    use_global_dt = .false.

    !if (nbinmax==0) nbinmax = 6
    istepfrac   = istepfrac + 1
    nbinmaxprev = nbinmax
    !--determine if dt needs to be decreased; if so, then this will be done
    !  in step the next time it is called;
    !  for global timestepping, this is called in the block where at_dump_time==.true.
    ! if (istepfrac==2**nbinmax) then
    !    twallperdump = reduceall_mpi('max', timer_lastdump%wall)
    !    call check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_ifactor,dtmax_log_dratio,&
    !                                  rhomaxold,rhomaxnow,nfulldump,use_global_dt)
    ! endif

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
#endif

    nsteps = nsteps + 1
!
!--evolve data for one timestep
!  for individual timesteps this is the shortest timestep
!
    !C_cool = 0.1
    dt = dtmax/2**nbinmax !FIXME
    call get_timings(t1,tcpu1)
    if ( use_sts ) then
       call step_sts(npart,nactive,time,dt,dtextforce,dtnew,iprint)
    else
       call step(npart,nactive,time,dt,dtextforce,dtnew)
    endif
    dtlast = dt

    !--timings for step call

    call get_timings(t2,tcpu2)
    call increment_timer(timer_step,t2-t1,tcpu2-tcpu1)
    call summary_counter(iosum_nreal,t2-t1)

#ifdef IND_TIMESTEPS
    tcheck = tcheck + dt

    !--update time in way that is free of round-off errors
    !print*, "CHECK tlast, time: ", tlast, time
    time = tbegin + istepfrac/real(2**nbinmaxprev)*dtmaxold
    !print*, "CHECK new time (istepfrac, nbinmaxprev, dtmaxold): ", time, istepfrac, nbinmaxprev, dtmaxold

    !--print efficiency of partial timestep
    if (id==master .and. iverbose >= 0 .and. nalivetot > 0) then
       if (nactivetot==nalivetot) then
          tall = t2-t1
       elseif (tall > 0.) then
          fracactive = nactivetot/real(nalivetot)
          speedup = (t2-t1)/tall
          if (iverbose >= 2) then
             if (speedup > 0) then
                efficiency = 100.*fracactive/speedup
             else
                efficiency = 0.
             endif
             write(iprint,"(1x,'(',3(a,f6.2,'%'),')')") &
                  'moved ',100.*fracactive,' of particles in ',100.*speedup, &
                  ' of time, efficiency = ',efficiency
          endif
       endif
    endif
    call update_time_per_bin(tcpu2-tcpu1,istepfrac,nbinmaxprev,timeperbin,inbin)
    nmovedtot = nmovedtot + nactivetot

    !--check that time is as it should be, may indicate error in individual timestep routines
    if (abs(tcheck-time) > 1.e-4) call warning('evolve','time out of sync',var='error',val=abs(tcheck-time))
    if (abs(tcheck-time) > 1.e-4) print*, "time, tcheck: ", time, tcheck

    if (id==master .and. (iverbose >= 1 .or. inbin <= 3)) &
       call print_dtlog_ind(iprint,istepfrac,2**nbinmaxprev,time,dt,nactivetot,tcpu2-tcpu1,npart)

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
    dt = min(dtforce,dtcourant,dterr,dtmax+epsilon(dtmax))
!
!--write log every step (NB: must print after dt has been set in order to identify timestep constraint)
!
    if (id==master) call print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,dtprint,dtinject,npart)
#endif

!    if (abs(dt) < 1e-8*dtmax) then
!       write(iprint,*) 'main loop: timestep too small, dt = ',dt
!       call quit   ! also writes dump file, safe here because called by all threads
!    endif

!   check that MPI threads are synchronised in time
    if (nprocs > 0) then
        timecheck = reduceall_mpi('+',time)
        if (abs(timecheck/nprocs - time) > 1.e-13) then
           call fatal('evolve','time differs between MPI threads',var='time',val=timecheck/nprocs)
        endif
    endif
!
!--Update timer from last dump to see if dtmax needs to be reduced
!
    call get_timings(t2,tcpu2)
    call increment_timer(timer_lastdump,t2-t1,tcpu2-tcpu1)

#ifdef CORRECT_BULK_MOTION
    call correct_bulk_motion()
#endif

    if (iverbose >= 1 .and. id==master) write(iprint,*)
    call flush(iprint)
    !--Write out log file prematurely (if requested based upon nstep, walltime)
    if ( summary_printnow() ) call summary_printout(iprint,nptmass)


 !enddo timestepping
 enddo timesubstepping
 !tlast = time
 !dtlast = dt
end subroutine amuse_new_step


! New particles
subroutine amuse_new_sph_particle(i, mass, x, y, z, vx, vy, vz, u, h)
    use part, only:igas,npart,npartoftype,xyzh,vxyzu,massoftype
    use part, only:abundance,iHI,ih2ratio
    use partinject, only:add_or_update_particle
#ifdef IND_TIMESTEPS
    use part, only:twas,ibin
    use timestep_ind, only:istepfrac,get_dt,nbinmax,change_nbinmax,get_newbin
    use timestep, only:dt,time,dtmax
    use timestep, only:C_cour,rhomaxnow
    use eos, only:gamma
#endif
    use timestep, only:dtextforce
    use units, only:umass,udist,utime
    implicit none
    integer :: n, i, itype
    double precision :: mass, x, y, z, vx, vy, vz, u, h
    double precision :: position(3), velocity(3)
#ifdef IND_TIMESTEPS
    integer(kind=1) :: nbinmaxprev
    real :: dtinject
    dtinject = huge(dtinject)
    ! dtmax = 0.01 !TODO This is arbitrarily set. Probably bad.
#endif
    ! Adding a particle of unknown temperature -> use a small cooling step
    dtextforce = 1.e-6

    itype = igas
    i = npart + 1
    ! print*, "**************   X position of particle ", i, x, "   ***************"
    position(1) = x
    position(2) = y
    position(3) = z
    velocity(1) = vx
    velocity(2) = vy
    velocity(3) = vz

    if (npartoftype(itype) == 0) then
        massoftype(itype) = mass
    endif
    call add_or_update_particle(itype,position,velocity,h, &
        u,i,npart,npartoftype,xyzh,vxyzu)
    abundance(:,i) = 0.
    abundance(ih2ratio,i) = 0.5
    abundance(iHI,i) = 1.  ! assume all gas is atomic hydrogen initially
    if (i == 1) then
        print*, "xyz vxyz u mass = ", x, y, z, vx, vy, vz, u, mass
        print*, "udist, utime, umass = ", udist, utime, umass
        print*, "x vx u mass in cgs = ", x*udist, vx*udist/utime, u*udist*udist/utime/utime, mass*umass
    endif

#ifdef IND_TIMESTEPS
    dtinject = C_cour * h / (gamma*(gamma-1)*u)**0.5
    nbinmaxprev = nbinmax
    call get_newbin(dtinject,dtmax,nbinmax,allow_decrease=.false.)
    ! not doing clever stuff: all particles go in the shortest possible bin.
    ! FIXME rethink this later...
    nbinmax = 30
    if (nbinmax > nbinmaxprev) then ! update number of bins if needed
       call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
       print*, "nbinmax (prev), time: ", nbinmax, nbinmaxprev, time
       print*, "npart:", npart
    endif
    ! put all injected particles on shortest bin
    ibin(i) = nbinmax
    twas(i) = time + 0.5*get_dt(dtmax,ibin(i))
    if (i == 1) then
        print*, x, y, z, vx, vy, vz, mass
    endif
#endif

end subroutine amuse_new_sph_particle

subroutine amuse_new_dm_particle(i, mass, x, y, z, vx, vy, vz)
    use part, only:idarkmatter,npart,npartoftype,xyzh,vxyzu,massoftype
    use partinject, only:add_or_update_particle
    implicit none
    integer :: n, i, itype
    double precision :: mass, x, y, z, vx, vy, vz, h_smooth, u
    double precision :: position(3), velocity(3)
  
    u = 0
    itype = idarkmatter
    i = npart + 1
    position(1) = x
    position(2) = y
    position(3) = z
    velocity(1) = vx
    velocity(2) = vy
    velocity(3) = vz
    h_smooth = 0.1 ! TODO set this to some default
    if (npartoftype(itype) == 0) then
        massoftype(itype) = mass
    endif

    call add_or_update_particle(itype,position,velocity,h_smooth, &
        u,i,npart,npartoftype,xyzh,vxyzu)
end subroutine

subroutine amuse_new_sink_particle(j, mass, x, y, z, vx, vy, vz, &
        radius, h_smooth)
    use io, only:fatal
    use part, only:nptmass,maxptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
    implicit none
    integer :: i, j
    double precision :: mass, x, y, z, vx, vy, vz, radius, h_smooth
    double precision :: position(3), velocity(3)
  
    nptmass = nptmass + 1
    ! Replace this with something AMUSE can handle
    if (nptmass > maxptmass) call fatal('creating new sink', 'nptmass > maxptmass')
    i = nptmass
    j = -i
    xyzmh_ptmass(:,i) = 0.
    xyzmh_ptmass(1,i) = x
    xyzmh_ptmass(2,i) = y
    xyzmh_ptmass(3,i) = z
    xyzmh_ptmass(4,i) = mass
    xyzmh_ptmass(ihacc,i) = radius
    xyzmh_ptmass(ihsoft,i) = h_smooth
    vxyz_ptmass(1,i) = vx
    vxyz_ptmass(2,i) = vy
    vxyz_ptmass(3,i) = vz
end subroutine


subroutine amuse_delete_particle(i)
    use part, only:kill_particle,xyzmh_ptmass
    integer, intent(in) :: i
    integer :: j
    if (i == abs(i)) then
        call kill_particle(i)
    else
        j = -i
        ! Sink particles can't be killed - so we just set its mass to zero
        xyzmh_ptmass(4,j) = 0
    endif
end subroutine

subroutine amuse_get_unit_length(unit_length_out)
    use units, only: udist
    implicit none
    double precision, intent(out) :: unit_length_out
    unit_length_out = udist
end subroutine amuse_get_unit_length

subroutine amuse_set_unit_length(unit_length_in)
    use units, only: udist,utime,umass,set_units
    implicit none
    double precision, intent(in) :: unit_length_in
    udist = unit_length_in
    !call set_units(dist=udist,time=utime,G=1.)
    print*, "set_unit_length called: utime/mass/dist: ", utime, umass, udist
end subroutine amuse_set_unit_length

subroutine amuse_get_unit_mass(unit_mass_out)
    use units, only: umass
    implicit none
    double precision, intent(out) :: unit_mass_out
    unit_mass_out = umass
end subroutine amuse_get_unit_mass

subroutine amuse_set_unit_mass(unit_mass_in)
    use units, only: umass,utime,udist,set_units
    implicit none
    double precision, intent(in) :: unit_mass_in
    umass = unit_mass_in
    call set_units(mass=umass,time=utime,G=1.)
    print*, "set_unit_mass called: utime/mass/dist: ", utime, umass, udist
end subroutine amuse_set_unit_mass

subroutine amuse_get_unit_time(unit_time_out)
    use units, only: utime
    implicit none
    double precision, intent(out) :: unit_time_out
    unit_time_out = utime
end subroutine amuse_get_unit_time

subroutine amuse_set_unit_time(unit_time_in)
    use units, only: utime,umass,udist,set_units
    implicit none
    double precision, intent(out) :: unit_time_in
    utime = unit_time_in
    call set_units(time=utime,mass=umass,G=1.)
    print*, "set_unit_time called: utime/mass/dist: ", utime, umass, udist
end subroutine amuse_set_unit_time

subroutine amuse_get_constant_solarm(solarm_out)
    use physcon, only:solarm
    implicit none
    double precision, intent(out) :: solarm_out
    solarm_out = solarm
end subroutine amuse_get_constant_solarm

subroutine amuse_get_constant_pc(pc_out)
    use physcon, only:pc
    implicit none
    double precision, intent(out) :: pc_out
    pc_out = pc
end subroutine amuse_get_constant_pc

subroutine amuse_get_constant_planckh(planckh_out)
    use physcon, only:planckh
    implicit none
    double precision, intent(out) :: planckh_out
    planckh_out = planckh
end subroutine amuse_get_constant_planckh

subroutine amuse_get_potential_energy(epot_out)
    use energies, only:epot
    implicit none
    double precision, intent(out) :: epot_out
    epot_out = epot
end subroutine amuse_get_potential_energy

subroutine amuse_get_kinetic_energy(ekin_out)
    use energies, only:ekin
    implicit none
    double precision, intent(out) :: ekin_out
    ekin_out = ekin
end subroutine amuse_get_kinetic_energy

subroutine amuse_get_thermal_energy(etherm_out)
    use energies, only:etherm
    implicit none
    double precision, intent(out) :: etherm_out
    etherm_out = etherm
end subroutine amuse_get_thermal_energy

subroutine amuse_get_time_step(dt_out)
    use units, only:utime
    use timestep, only:dtmax
    implicit none
    double precision, intent(out) :: dt_out
    dt_out = dtmax
end subroutine amuse_get_time_step

subroutine amuse_get_number_of_sph_particles(n)
    use part, only:npartoftype,igas
    implicit none
    integer, intent(out) :: n
    logical :: nodisabled
    nodisabled = .true.
    n = npartoftype(igas)
end subroutine amuse_get_number_of_sph_particles

subroutine amuse_get_number_of_particles(n)
    use part, only:npart
    implicit none
    integer, intent(out) :: n
    logical :: nodisabled
    nodisabled = .true.
    call amuse_get_npart(n, nodisabled)
end subroutine amuse_get_number_of_particles

subroutine amuse_get_time(time_out)
    use timestep, only:time
    implicit none
    double precision, intent(out) :: time_out
    time_out = time
end subroutine amuse_get_time

subroutine amuse_get_density(i, rho)
    use part, only:rhoh,iphase,massoftype,xyzh
    implicit none
    integer :: i
    double precision :: pmassi
    double precision, intent(out) :: rho
    pmassi = massoftype(abs(iphase(i)))
    rho = rhoh(xyzh(4,i), pmassi)
end subroutine amuse_get_density

subroutine amuse_get_pressure(i, p)
    use part, only:rhoh,iphase,massoftype,xyzh
    use eos, only:ieos,equationofstate
    implicit none
    integer :: i, eos_type
    double precision :: pmassi, ponrho, rho, spsound, x, y, z
    double precision, intent(out) :: p
    eos_type = ieos
    pmassi = massoftype(abs(iphase(i)))
    call amuse_get_density(i, rho)
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    call equationofstate(eos_type,ponrho,spsound,rho,x,y,z)
    p = ponrho * rho
end subroutine amuse_get_pressure

subroutine amuse_get_mass(i, part_mass)
    use part, only:iphase,massoftype,xyzmh_ptmass
    implicit none
    double precision, intent(out) :: part_mass
    integer, intent(in) :: i
    integer :: j
    if (i == abs(i)) then
        part_mass = massoftype(abs(iphase(i)))
    else
        j = -i
        part_mass = xyzmh_ptmass(4,j)
    endif
end subroutine amuse_get_mass

subroutine amuse_get_state_gas(i, mass, x, y, z, vx, vy, vz, u, h)
    implicit none
    integer :: i
    double precision, intent(inout) :: mass, x, y, z, vx, vy, vz, u, h
    call amuse_get_mass(i, mass)
    call amuse_get_position(i, x, y, z)
    call amuse_get_velocity(i, vx, vy, vz)
    call amuse_get_internal_energy(i, u)
    call amuse_get_smoothing_length(i, h)
end subroutine amuse_get_state_gas

subroutine amuse_get_state_dm(i, mass, x, y, z, vx, vy, vz)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz
    call amuse_get_mass(i, mass)
    call amuse_get_position(i, x, y, z)
    call amuse_get_velocity(i, vx, vy, vz)
end subroutine amuse_get_state_dm

subroutine amuse_get_state_sink(j, mass, x, y, z, vx, vy, vz, radius)
    implicit none
    integer :: j
    double precision :: mass, x, y, z, vx, vy, vz, radius
    call amuse_get_mass(j, mass)
    call amuse_get_position(j, x, y, z)
    call amuse_get_velocity(j, vx, vy, vz)
    call amuse_get_sink_radius(j, radius)
end subroutine amuse_get_state_sink

subroutine amuse_get_sink_radius(j, radius)
    use part, only:xyzmh_ptmass, ihacc
    implicit none
    integer :: j
    double precision :: radius
    radius = xyzmh_ptmass(ihacc, j)
end subroutine amuse_get_sink_radius

subroutine amuse_get_position(i, x, y, z)
    use part, only:xyzh,xyzmh_ptmass
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(out) :: x, y, z
    if (i == abs(i)) then
        x = xyzh(1, i)
        y = xyzh(2, i)
        z = xyzh(3, i)
    else
        j = -i
        x = xyzmh_ptmass(1, j)
        y = xyzmh_ptmass(2, j)
        z = xyzmh_ptmass(3, j)
    endif
end subroutine amuse_get_position

subroutine amuse_get_velocity(i, vx, vy, vz)
    use part, only:vxyzu,vxyz_ptmass
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(out) :: vx, vy, vz
    if (i == abs(i)) then
        vx = vxyzu(1, i)
        vy = vxyzu(2, i)
        vz = vxyzu(3, i)
    else
        j = -i
        vx = vxyz_ptmass(1, j)
        vy = vxyz_ptmass(2, j)
        vz = vxyz_ptmass(3, j)
    endif
end subroutine amuse_get_velocity

subroutine amuse_get_smoothing_length(i, h)
    use part, only:xyzh,xyzmh_ptmass,ihsoft
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(out) :: h
    if (i == abs(i)) then
        h = xyzh(4, i)
    else
        j = -i
        h = xyzmh_ptmass(ihsoft,j)
    endif
end subroutine amuse_get_smoothing_length

subroutine amuse_get_radius(i, radius)
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(out) :: radius
    if (i == abs(i)) then
        call amuse_get_smoothing_length(i, radius)
    else
        j = -i
        call amuse_get_sink_radius(j, radius)
    endif
end subroutine amuse_get_radius

subroutine amuse_get_internal_energy(i, u)
    use dim, only:maxvxyzu
    use part, only:vxyzu
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: u

    if (maxvxyzu >= 4) then
        u = vxyzu(4, i)
    else
        u = 0
    endif
end subroutine

subroutine amuse_set_hi_abundance(i, hi_abundance)
    use part, only:abundance,iHI
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: hi_abundance

    abundance(iHI, i) = hi_abundance
end subroutine

subroutine amuse_get_hi_abundance(i, hi_abundance)
    use part, only:abundance,iHI
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: hi_abundance

    hi_abundance = abundance(iHI, i)
end subroutine

subroutine amuse_set_proton_abundance(i, proton_abundance)
    use part, only:abundance,iproton
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: proton_abundance

    abundance(iproton, i) = proton_abundance
end subroutine

subroutine amuse_get_proton_abundance(i, proton_abundance)
    use part, only:abundance,iproton
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: proton_abundance

    proton_abundance = abundance(iproton, i)
end subroutine

subroutine amuse_set_electron_abundance(i, electron_abundance)
    use part, only:abundance,ielectron
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: electron_abundance

    abundance(ielectron, i) = electron_abundance
end subroutine

subroutine amuse_get_electron_abundance(i, electron_abundance)
    use part, only:abundance,ielectron
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: electron_abundance

    electron_abundance = abundance(ielectron, i)
end subroutine

subroutine amuse_set_co_abundance(i, co_abundance)
    use part, only:abundance,ico
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: co_abundance

    abundance(ico, i) = co_abundance
end subroutine

subroutine amuse_get_co_abundance(i, co_abundance)
    use part, only:abundance,ico
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: co_abundance

    co_abundance = abundance(ico, i)
end subroutine

subroutine amuse_set_h2ratio(i, h2ratio)
    use part, only:abundance,iHI,ih2ratio
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: h2ratio

    abundance(ih2ratio, i) = h2ratio
end subroutine

subroutine amuse_get_h2ratio(i, h2ratio)
    use part, only:abundance,iHI,ih2ratio
    implicit none
    integer, intent(in) :: i
    double precision, intent(out) :: h2ratio

    h2ratio = abundance(ih2ratio, i)
end subroutine

subroutine amuse_set_time_step(dt_in)
    use units, only:utime
    use timestep, only:dtmax
    implicit none
    double precision, intent(in) :: dt_in
    dtmax = dt_in
    print*, "dtmax: ", dtmax
end subroutine

subroutine amuse_set_mass(i, part_mass)
    use part, only:iphase,massoftype,xyzmh_ptmass
    implicit none
    double precision, intent(in) :: part_mass
    integer, intent(in) :: i
    integer :: j
    if (i == abs(i)) then
        massoftype(abs(iphase(i))) = part_mass
    else
        j = -i
        xyzmh_ptmass(4,j) = part_mass
    endif
end subroutine

subroutine amuse_set_state_gas(i, mass, x, y, z, vx, vy, vz, u, h)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz, u, h
    call amuse_set_mass(i, mass)
    call amuse_set_position(i, x, y, z)
    call amuse_set_velocity(i, vx, vy, vz)
    call amuse_set_internal_energy(i, u)
    call amuse_set_smoothing_length(i, h)
end subroutine

subroutine amuse_set_state_dm(i, mass, x, y, z, vx, vy, vz)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz
    call amuse_set_mass(i, mass)
    call amuse_set_position(i, x, y, z)
    call amuse_set_velocity(i, vx, vy, vz)
end subroutine

subroutine amuse_set_state_sink(j, mass, x, y, z, vx, vy, vz, radius)
    implicit none
    integer :: j
    double precision :: mass, x, y, z, vx, vy, vz, radius
    call amuse_set_mass(j, mass)
    call amuse_set_position(j, x, y, z)
    call amuse_set_velocity(j, vx, vy, vz)
    call amuse_set_sink_radius(j, radius)
end subroutine

subroutine amuse_set_sink_radius(j, radius)
    use part, only:xyzmh_ptmass, ihacc
    implicit none
    integer :: j
    double precision :: radius
    xyzmh_ptmass(ihacc, j) = radius
end subroutine

subroutine amuse_set_position(i, x, y, z)
    use part, only:xyzh,xyzmh_ptmass
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(in) :: x, y, z
    if (i == abs(i)) then
        xyzh(1, i) = x
        xyzh(2, i) = y
        xyzh(3, i) = z
    else
        j = -i
        xyzmh_ptmass(1, j) = x
        xyzmh_ptmass(2, j) = y
        xyzmh_ptmass(3, j) = z
    endif
end subroutine

subroutine amuse_set_velocity(i, vx, vy, vz)
    use part, only:vxyzu,vxyz_ptmass
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(in) :: vx, vy, vz
    if (i == abs(i)) then
        vxyzu(1, i) = vx
        vxyzu(2, i) = vy
        vxyzu(3, i) = vz
    else
        j = -i
        vxyz_ptmass(1, j) = vx
        vxyz_ptmass(2, j) = vy
        vxyz_ptmass(3, j) = vz
    endif
end subroutine

subroutine amuse_set_smoothing_length(i, h)
    use part, only:xyzh,xyzmh_ptmass,ihsoft
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(in) :: h
    if (i == abs(i)) then
        xyzh(4, i) = h
    else
        j = -i
        xyzmh_ptmass(ihsoft, j) = h
    endif
end subroutine

subroutine amuse_set_radius(i, radius)
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision :: radius
    if (i == abs(i)) then
        call amuse_set_smoothing_length(i, radius)
    else
        j = -i
        call amuse_set_sink_radius(j, radius)
    endif
end subroutine

subroutine amuse_set_internal_energy(i, u)
    use dim, only:maxvxyzu
    use part, only:vxyzu
    use timestep, only:dtextforce
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: u
    if (maxvxyzu >= 4) then
        vxyzu(4, i) = u
    endif
    ! Changing temperature -> better use a small cooling step
    dtextforce = 1.e-8
end subroutine

subroutine amuse_evolve_model(tmax_in)
    use timestep, only:tmax, time, dt, dtmax, rhomaxnow
    ! use evolvesplit, only:init_step, finalize_step
#ifdef IND_TIMESTEPS
    use timestep_ind, only:istepfrac
#endif
    use options, only:rhofinal1
    use ptmass, only:rho_crit
    use part, only:npart
    use step_lf_global, only:init_step
    implicit none
    double precision, intent(in) :: tmax_in
    logical :: maximum_density_reached
    real :: tlast

    tmax = tmax_in + epsilon(tmax_in)
    !dtmax = (tmax - time)
    
    tlast = time
    timestepping: do while (time < tmax)
#ifdef IND_TIMESTEPS
        istepfrac = 0
        print*, "*****init timestep"
        call init_step(npart, time, dtmax)
#endif
        print*, "*****new_step - Time; tmax: ", time, tmax
        call amuse_new_step(tlast)
        print*, "*****new_step done - Time; tmax: ", time, tmax
    enddo timestepping
end subroutine

!
! Setters and getters for parameters
!

! Setters

subroutine amuse_set_c_courant(C_cour_in)
    use timestep, only:C_cour
    implicit none
    double precision, intent(in) :: C_cour_in
    C_cour = C_cour_in
end subroutine

subroutine amuse_set_c_force(C_force_in)
    use timestep, only:C_force
    implicit none
    double precision, intent(in) :: C_force_in
    C_force = C_force_in
end subroutine

subroutine amuse_set_c_cool(C_cool_in)
    use timestep, only:C_cool
    implicit none
    double precision, intent(in) :: C_cool_in
    C_cool = C_cool_in
end subroutine

subroutine amuse_set_tolv(tolv_in)
    use timestep, only:tolv
    implicit none
    double precision, intent(in) :: tolv_in
    tolv = tolv_in
end subroutine

subroutine amuse_set_hfact(hfact_in)
    use part, only:hfact
    implicit none
    double precision, intent(in) :: hfact_in
    hfact = hfact_in
end subroutine

subroutine amuse_set_tolh(tolh_in)
    use options, only:tolh
    implicit none
    double precision, intent(in) :: tolh_in
    tolh = tolh_in
end subroutine

subroutine amuse_set_tree_accuracy(tree_accuracy_in)
    use kdtree, only:tree_accuracy
    implicit none
    double precision, intent(in) :: tree_accuracy_in
    tree_accuracy = tree_accuracy_in
end subroutine

subroutine amuse_set_alpha(alpha_in)
    use options, only:alpha
    implicit none
    double precision, intent(in) :: alpha_in
    alpha = alpha_in
end subroutine

subroutine amuse_set_alphamax(alphamax_in)
    use options, only:alphamax
    implicit none
    double precision, intent(in) :: alphamax_in
    alphamax = alphamax_in
end subroutine

subroutine amuse_set_beta(beta_in)
    use options, only:beta
    implicit none
    double precision, intent(in) :: beta_in
    beta = beta_in
end subroutine

subroutine amuse_set_avdecayconst(avdecayconst_in)
    use options, only:avdecayconst
    implicit none
    double precision, intent(in) :: avdecayconst_in
    avdecayconst = avdecayconst_in
end subroutine

subroutine amuse_set_idamp(idamp_in)
    use options, only:idamp
    implicit none
    integer, intent(in) :: idamp_in
    idamp = idamp_in
end subroutine

subroutine amuse_set_ieos(ieos_in)
    use eos, only:ieos
    implicit none
    integer, intent(in) :: ieos_in
    ieos = ieos_in
end subroutine

subroutine amuse_set_icooling(icooling_in)
    use io, only:id,master,iprint
    use options, only:icooling,iexternalforce
    use part, only:h2chemistry
    use chem, only:init_chem
    use cooling, only:init_cooling,init_cooling_type
    use h2cooling, only:init_h2cooling
    implicit none
    integer :: ierr
    integer, intent(in) :: icooling_in
    icooling = icooling_in
    if (icooling > 0) then
        if (h2chemistry) then
            if (id==master) write(iprint,*) 'initialising cooling function...'
            call init_chem()
            call init_h2cooling()
        else
            call init_cooling(ierr)
            if (ierr /= 0) call fatal('initial','error initialising cooling')
        endif
        call init_cooling_type(h2chemistry)
    endif
end subroutine

subroutine amuse_set_polyk(polyk_in)
    use eos, only:polyk
    implicit none
    double precision, intent(in) :: polyk_in
    polyk = polyk_in
end subroutine

subroutine amuse_set_mu(mu_in)
    use eos, only:gmw
    implicit none
    double precision, intent(in) :: mu_in
    gmw = mu_in
end subroutine

subroutine amuse_set_rhofinal(rhofinal_in)
    use options, only:rhofinal_cgs, rhofinal1
    use units, only:unit_density
    implicit none
    double precision, intent(in) :: rhofinal_in
    rhofinal_cgs = rhofinal_in * unit_density
    if (rhofinal_cgs > 0.) then
        rhofinal1 = unit_density/rhofinal_cgs
    else
        rhofinal1 = 0.0
    endif
end subroutine

subroutine amuse_set_rho_crit(rho_crit_in)
    use units, only:unit_density
    use ptmass, only:rho_crit_cgs
    implicit none
    double precision, intent(in) :: rho_crit_in
    rho_crit_cgs = rho_crit_in * unit_density
end subroutine

subroutine amuse_set_r_crit(r_crit_in)
    use ptmass, only:r_crit
    implicit none
    double precision, intent(in) :: r_crit_in
    r_crit = r_crit_in
end subroutine

subroutine amuse_set_h_acc(h_acc_in)
    use ptmass, only:h_acc
    implicit none
    double precision, intent(in) :: h_acc_in
    h_acc = h_acc_in
end subroutine

subroutine amuse_set_h_soft_sinkgas(h_soft_sinkgas_in)
    use ptmass, only:h_soft_sinkgas
    implicit none
    double precision, intent(in) :: h_soft_sinkgas_in
    h_soft_sinkgas = h_soft_sinkgas_in
end subroutine

subroutine amuse_set_h_soft_sinksink(h_soft_sinksink_in)
    use ptmass, only:h_soft_sinksink
    implicit none
    double precision, intent(in) :: h_soft_sinksink_in
    h_soft_sinksink = h_soft_sinksink_in
end subroutine

subroutine amuse_set_f_acc(f_acc_in)
    use ptmass, only:f_acc
    implicit none
    double precision, intent(in) :: f_acc_in
    f_acc = f_acc_in
end subroutine

subroutine amuse_set_iexternalforce(iexternalforce_in)
    use options, only:iexternalforce
    implicit none
    integer, intent(in) :: iexternalforce_in
    iexternalforce = iexternalforce_in
end subroutine

subroutine amuse_set_irealvisc(irealvisc_in)
    use viscosity, only:irealvisc
    implicit none
    integer, intent(in) :: irealvisc_in
    irealvisc = irealvisc_in
end subroutine

subroutine amuse_set_shearparam(shearparam_in)
    use viscosity, only:shearparam
    implicit none
    double precision, intent(in) :: shearparam_in
    shearparam = shearparam_in
end subroutine

subroutine amuse_set_bulkvisc(bulkvisc_in)
    use viscosity, only:bulkvisc
    implicit none
    double precision, intent(in) :: bulkvisc_in
    bulkvisc = bulkvisc_in
end subroutine

subroutine amuse_set_gamma(gamma_in)
    use eos, only:gamma
    implicit none
    double precision, intent(in) :: gamma_in
    gamma = gamma_in
end subroutine

subroutine amuse_set_umass(umass_in)
    use units, only:umass
    implicit none
    double precision, intent(in) :: umass_in
    umass = umass_in
end subroutine

subroutine amuse_set_udist(udist_in)
    use units, only:udist
    implicit none
    double precision, intent(in) :: udist_in
    udist = udist_in
end subroutine

subroutine amuse_set_utime(utime_in)
    use units, only:utime
    implicit none
    double precision, intent(in) :: utime_in
    utime = utime_in
end subroutine

! End of Setters

! Getters

subroutine amuse_get_c_courant(C_cour_out)
    use timestep, only:C_cour
    implicit none
    double precision, intent(out) :: C_cour_out
    C_cour_out = C_cour
end subroutine

subroutine amuse_get_c_force(C_force_out)
    use timestep, only:C_force
    implicit none
    double precision, intent(out) :: C_force_out
    C_force_out = C_force
end subroutine

subroutine amuse_get_c_cool(C_cool_out)
    use timestep, only:C_cool
    implicit none
    double precision, intent(out) :: C_cool_out
    C_cool_out = C_cool
end subroutine

subroutine amuse_get_tolv(tolv_out)
    use timestep, only:tolv
    implicit none
    double precision, intent(out) :: tolv_out
    tolv_out = tolv
end subroutine

subroutine amuse_get_hfact(hfact_out)
    use part, only:hfact
    implicit none
    double precision, intent(out) :: hfact_out
    hfact_out = hfact
end subroutine

subroutine amuse_get_tolh(tolh_out)
    use options, only:tolh
    implicit none
    double precision, intent(out) :: tolh_out
    tolh_out = tolh
end subroutine

subroutine amuse_get_tree_accuracy(tree_accuracy_out)
    use kdtree, only:tree_accuracy
    implicit none
    double precision, intent(out) :: tree_accuracy_out
    tree_accuracy_out = tree_accuracy
end subroutine

subroutine amuse_get_alpha(alpha_out)
    use options, only:alpha
    implicit none
    double precision, intent(out) :: alpha_out
    alpha_out = alpha
end subroutine

subroutine amuse_get_alphamax(alphamax_out)
    use options, only:alphamax
    implicit none
    double precision, intent(out) :: alphamax_out
    alphamax_out = alphamax
end subroutine

subroutine amuse_get_beta(beta_out)
    use options, only:beta
    implicit none
    double precision, intent(out) :: beta_out
    beta_out = beta
end subroutine

subroutine amuse_get_avdecayconst(avdecayconst_out)
    use options, only:avdecayconst
    implicit none
    double precision, intent(out) :: avdecayconst_out
    avdecayconst_out = avdecayconst
end subroutine

subroutine amuse_get_idamp(idamp_out)
    use options, only:idamp
    implicit none
    integer, intent(out) :: idamp_out
    idamp_out = idamp
end subroutine

subroutine amuse_get_ieos(ieos_out)
    use eos, only:ieos
    implicit none
    integer, intent(out) :: ieos_out
    ieos_out = ieos
end subroutine

subroutine amuse_get_icooling(icooling_out)
    use options, only:icooling
    implicit none
    integer, intent(out) :: icooling_out
    icooling_out = icooling
end subroutine

subroutine amuse_get_polyk(polyk_out)
    use eos, only:polyk
    implicit none
    double precision, intent(out) :: polyk_out
    polyk_out = polyk
end subroutine

subroutine amuse_get_mu(mu_out)
    use eos, only:gmw
    implicit none
    double precision, intent(out) :: mu_out
    mu_out = gmw
end subroutine

subroutine amuse_get_rhofinal(rhofinal_out)
    use options, only:rhofinal_cgs
    use units, only:unit_density
    implicit none
    double precision, intent(out) :: rhofinal_out
    rhofinal_out = rhofinal_cgs / unit_density
end subroutine

subroutine amuse_get_rho_crit(rho_crit_out)
    use units, only:unit_density
    use ptmass, only:rho_crit_cgs
    implicit none
    double precision, intent(out) :: rho_crit_out
    rho_crit_out = rho_crit_cgs / unit_density
end subroutine

subroutine amuse_get_r_crit(r_crit_out)
    use ptmass, only:r_crit
    implicit none
    double precision, intent(out) :: r_crit_out
    r_crit_out = r_crit
end subroutine

subroutine amuse_get_h_acc(h_acc_out)
    use ptmass, only:h_acc
    implicit none
    double precision, intent(out) :: h_acc_out
    h_acc_out = h_acc
end subroutine

subroutine amuse_get_h_soft_sinkgas(h_soft_sinkgas_out)
    use ptmass, only:h_soft_sinkgas
    implicit none
    double precision, intent(out) :: h_soft_sinkgas_out
    h_soft_sinkgas_out = h_soft_sinkgas
end subroutine

subroutine amuse_get_h_soft_sinksink(h_soft_sinksink_out)
    use ptmass, only:h_soft_sinksink
    implicit none
    double precision, intent(out) :: h_soft_sinksink_out
    h_soft_sinksink_out = h_soft_sinksink
end subroutine

subroutine amuse_get_f_acc(f_acc_out)
    use ptmass, only:f_acc
    implicit none
    double precision, intent(out) :: f_acc_out
    f_acc_out = f_acc
end subroutine

subroutine amuse_get_iexternalforce(iexternalforce_out)
    use options, only:iexternalforce
    implicit none
    integer, intent(out) :: iexternalforce_out
    iexternalforce_out = iexternalforce
end subroutine

subroutine amuse_get_irealvisc(irealvisc_out)
    use viscosity, only:irealvisc
    implicit none
    integer, intent(out) :: irealvisc_out
    irealvisc_out = irealvisc
end subroutine

subroutine amuse_get_shearparam(shearparam_out)
    use viscosity, only:shearparam
    implicit none
    double precision, intent(out) :: shearparam_out
    shearparam_out = shearparam
end subroutine

subroutine amuse_get_bulkvisc(bulkvisc_out)
    use viscosity, only:bulkvisc
    implicit none
    double precision, intent(out) :: bulkvisc_out
    bulkvisc_out = bulkvisc
end subroutine

subroutine amuse_get_gamma(gamma_out)
    use eos, only:gamma
    implicit none
    double precision, intent(out) :: gamma_out
    gamma_out = gamma
end subroutine

subroutine amuse_get_umass(umass_out)
    use units, only:umass
    implicit none
    double precision, intent(out) :: umass_out
    umass_out = umass
end subroutine

subroutine amuse_get_utime(utime_out)
    use units, only:utime
    implicit none
    double precision, intent(out) :: utime_out
    utime_out = utime
end subroutine

subroutine amuse_get_udist(udist_out)
    use units, only:udist
    implicit none
    double precision, intent(out) :: udist_out
    udist_out = udist
end subroutine

! End of Getters
