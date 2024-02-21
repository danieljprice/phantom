! AMUSE interface library for Phantom
! (c) 2019 - 2024 Steven Rieder

!
! Initialize Phantom and set default parameters
!

subroutine amuse_initialize_code_new()
 use dim,             only:tagline
 use dim, only: maxp_hard
 use memory, only:allocate_memory
 use mpiutils,        only:init_mpi,finalise_mpi
 use initial,         only:initialise,finalise,startrun,endrun
 use io,              only:id,master,nprocs,set_io_unit_numbers,die
 use evolve,          only:evol
 use units, only:set_units,utime,umass,udist,unit_density
 use physcon, only:gram,seconds,solarm,pc,au
 implicit none
 integer            :: nargs
 character(len=120) :: infile,logfile,evfile,dumpfile

 id = 0

 call init_mpi(id,nprocs)
 call allocate_memory(int(maxp_hard,kind=8))
 call set_io_unit_numbers()
 call initialise()
 call amuse_set_defaults() ! replaces reading infile
 call set_units(dist=1 * au,mass=1.*solarm,G=1.)
 
end subroutine amuse_initialize_code_new

subroutine amuse_initialize_code()
    use dim, only:maxp,maxp_hard,maxvxyzu
    use io, only:id,nprocs,iverbose
    use mpiutils, only:init_mpi
    use memory, only:allocate_memory
    use units, only:set_units,utime,umass,udist,unit_density
    use physcon, only:gram,seconds,solarm,pc,au
    use timestep, only:dtmax,dtextforce
    use initial, only:initialise
#ifdef IND_TIMESTEPS
    use timestep_ind,     only:istepfrac,ibinnow
    use part,             only:ibin,ibin_old,ibin_wake
#else
    use timestep,         only:dtcourant,dtforce
#endif
    implicit none
    iverbose=5
    call init_mpi(id,nprocs)
    call allocate_memory(int(maxp_hard,kind=8))
    call initialise()
    call amuse_set_defaults()
    call amuse_set_polyk(0.)
    !print*, "maxvxyzu: ", maxvxyzu
    !maxvxyzu = 4
    !call set_units(dist=50.*pc,mass=4600.*solarm,G=1.)
    !call set_units(dist=1.d20,mass=1.d40,G=1.)
    !call set_units(dist=1.*pc,mass=1.*solarm,G=1.)
    ! call set_units(dist=1.,mass=1.,time=1.)

    !umass = 1.98892d33 * gram  ! AMUSE MSun
    !utime = 60 * 60 * 24 * 365.25 * 1d6 * seconds  ! 10 Julian Kyr
    !call set_units(time=utime,mass=umass,G=1.)
    !call set_units(dist=0.1*pc,mass=1.*solarm,G=1.)
    call set_units(dist=1 * au,mass=1.*solarm,G=1.)
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

    !call amuse_initialize_wind()
end subroutine amuse_initialize_code

subroutine amuse_set_phantom_option(name, valstring, imatch)
    ! This subroutine is meant to be a replacement for read_infile
    use inject,           only:read_options_inject
    use dust,             only:read_options_dust
    use growth,           only:read_options_growth
    use metric,           only:read_options_metric
    use nicil_sup,        only:read_options_nicil
    use dust_formation,   only:read_options_dust_formation
    use ptmass_radiation, only:read_options_ptmass_radiation
    use eos,              only:read_options_eos
    use cooling,          only:read_options_cooling
    use ptmass,           only:read_options_ptmass
    use damping,          only:read_options_damping
    use gravwaveutils,    only:read_options_gravitationalwaves
    use boundary_dyn,     only:read_options_boundary
    implicit none
    character(*), intent(inout):: name, valstring
    logical:: imatch, igotall
    integer:: ierr
    imatch = .false.
    if (.not.imatch) call read_options_inject(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_dust_formation(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_ptmass_radiation(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_dust(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_growth(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_metric(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_nicil(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_eos(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_cooling(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_damping(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_ptmass(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_gravitationalwaves(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) call read_options_boundary(name,valstring,imatch,igotall,ierr)
    if (.not.imatch) write(*,*) "Could not set option ", name
end subroutine amuse_set_phantom_option

subroutine amuse_initialize_wind()
    ! instead of reading a wind setup, set values here
    use inject, only:read_options_inject
    use dust_formation, only: read_options_dust_formation
    use ptmass_radiation, only: read_options_ptmass_radiation
    logical:: imatch, igotall
    integer:: ierr

    call read_options_inject("sonic_type", "0", imatch, igotall, ierr)
    call read_options_inject("wind_velocity", "20.", imatch, igotall, ierr)
    call read_options_inject("wind_inject_radius", "2.000", imatch, igotall, ierr)    ! wind injection radius (au, if 0 takes Rstar)
    call read_options_inject("wind_mass_rate", "1.000E-05", imatch, igotall, ierr)    ! wind mass loss rate (Msun/yr)
    call read_options_inject("wind_temperature", "2500.", imatch, igotall, ierr)    ! wind temperature at injection radius (K, if 0 takes Teff)
    call read_options_inject("iwind_resolution", "5", imatch, igotall, ierr)    ! if<>0 set number of particles on the sphere, reset particle mass
    call read_options_inject("nfill_domain", "0", imatch, igotall, ierr)    ! number of spheres used to set the background density profile
    call read_options_inject("wind_shell_spacing", "1.000", imatch, igotall, ierr)    ! desired ratio of sphere spacing to particle spacing
    call read_options_inject("iboundary_spheres", "5", imatch, igotall, ierr)    ! number of boundary spheres (integer)
    call read_options_inject("outer_boundary", "30.", imatch, igotall, ierr)    ! delete gas particles outside this radius (au)
    call read_options_inject("rkill", "-1.000", imatch, igotall, ierr)    ! deactivate particles outside this radius (<0 is off)

    !# options controlling dust
    call read_options_dust_formation("idust_opacity", "0", imatch, igotall, ierr)    ! compute dust opacity (0=off,1 (bowen), 2 (nucleation))

    !# options controlling radiation pressure from sink particles
    call read_options_ptmass_radiation("isink_radiation", "1", imatch, igotall, ierr)    ! sink radiation pressure method (0=off,1=alpha,2=dust,3=alpha+dust)
    call read_options_ptmass_radiation("alpha_rad", "1.000", imatch, igotall, ierr)

end subroutine amuse_initialize_wind

subroutine amuse_commit_parameters()
end subroutine amuse_commit_parameters

subroutine amuse_commit_particles()
    !use eos, only:ieos,init_eos
    !use deriv, only:derivs
    !use part, only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
    !        dustprop,ddustprop,dustfrac,ddustevol,dens,drad,radprop,dustevol,&
    !        eos_vars,metrics,pxyzu,rad
    !use timestep, only:time,dtmax
    !use units, only:udist,utime,umass
    use initial, only:startrun
    use energies, only: mtot
    !use timestep, only:dtmax
    implicit none
    !integer :: ierr
    !integer :: nactive
    !double precision :: dt, dtnew_first
    character(len=120) :: infile,logfile,evfile,dumpfile
    call startrun(infile,logfile,evfile,dumpfile, .true.)
    print*, "total mass in code unit: ", mtot

    !dtnew_first = dtmax
    !dt = 0
    !nactive = npart
    print*, "COMMIT_PARTICLES, calling amuse_init_evol"
    !call amuse_init_evol()
    call amuse_evol(.true.)
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
   use options, only: set_default_options,iexternalforce,write_files
   implicit none

   call set_default_options()
   ! A few changes to Phantom's defaults
   call amuse_set_gamma(1.)
   call amuse_set_ieos(1)
   write_files=.false.

end subroutine amuse_set_defaults

! This initialises things. This really should only be called once, before the first step.
subroutine amuse_evol(amuse_initialise)
 use io,               only:iprint,iwritein,id,master,iverbose,&
                            flush_warnings,nprocs,fatal,warning
 use timestep,         only:time,tmax,dt,dtmax,nmax,nout,nsteps,dtextforce,rhomaxnow,&
                            dtmax_ifactor,dtmax_ifactorWT,dtmax_dratio,check_dtmax_for_decrease,&
                            idtmax_n,idtmax_frac,idtmax_n_next,idtmax_frac_next
 use evwrite,          only:write_evfile,write_evlog
 use energies,         only:etot,totmom,angtot,mdust,np_cs_eq_0,np_e_eq_0,hdivBB_xa,&
                            compute_energies
 use checkconserved,   only:etot_in,angtot_in,totmom_in,mdust_in,&
                            init_conservation_checks,check_conservation_error,&
                            check_magnetic_stability
 use dim,              only:maxvxyzu,mhd,periodic,idumpfile
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
 use part,             only:rad,radprop
 use radiation_utils,  only:update_radenergy
 use timestep,         only:dtrad
#ifdef LIVE_ANALYSIS
 use analysis,         only:do_analysis
 use part,             only:igas
 use fileutils,        only:numfromfile
 use io,               only:ianalysis
#endif
 use part,             only:npart,nptmass,xyzh,vxyzu,fxyzu,fext,divcurlv,massoftype, &
                            xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gravity,iboundary, &
                            fxyz_ptmass_sinksink,ntot,poten,ndustsmall,accrete_particles_outside_sphere
 use quitdump,         only:quit
 use ptmass,           only:icreate_sinks,ptmass_create,ipart_rhomax,pt_write_sinkev,calculate_mdot
 use io_summary,       only:iosum_nreal,summary_counter,summary_printout,summary_printnow
 use externalforces,   only:iext_spiral
 use boundary_dyn,     only:dynamic_bdy,update_boundaries
#ifdef MFLOW
 use mf_write,         only:mflow_write
#endif
#ifdef VMFLOW
 use mf_write,         only:vmflow_write
#endif
#ifdef BINPOS
 use mf_write,         only:binpos_write
#endif
 implicit none
 logical, intent(in) :: amuse_initialise

 integer   :: flag
 character(len=256)    :: infile
 character(len=256)    :: logfile,evfile,dumpfile
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
 integer         :: npart_old
#ifdef INJECT_PARTICLES
#endif
 logical         :: fulldump,abortrun,abortrun_bdy,at_dump_time,writedump
 logical         :: should_conserve_energy,should_conserve_momentum,should_conserve_angmom
 logical         :: should_conserve_dustmass
 logical         :: use_global_dt
 integer         :: j,nskip,nskipped,nevwrite_threshold,nskipped_sink,nsinkwrite_threshold
 character(len=120) :: dumpfile_orig

 tzero     = time
 tlast         = time
 if (amuse_initialise) then
 tprint    = 0.
 nsteps    = 0
 nsteplast = 0
 dtlast    = 0.
 dtinject  = huge(dtinject)
 dtrad     = huge(dtrad)
 np_cs_eq_0 = 0
 np_e_eq_0  = 0
 abortrun_bdy = .false.

 call init_conservation_checks(should_conserve_energy,should_conserve_momentum,&
                               should_conserve_angmom,should_conserve_dustmass)
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

#ifdef IND_TIMESTEPS
 use_global_dt = .false.
 istepfrac     = 0
 tlast         = time
 write(*,*) "\n\n\n***********tlast: ", tlast, "\n\n\n"
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

 endif !(amuse_initialise)
!end subroutine amuse_init_evol

 if (.not.amuse_initialise) then
!subroutine amuse_new_step(tlast)
! use io,               only:iprint,iwritein,id,master,iverbose,flush_warnings,nprocs,fatal,warning
! use timestep,         only:time,tmax,dt,dtmax,nmax,nout,nsteps,dtextforce,rhomaxnow,&
!                            dtmax_ifactor,dtmax_dratio,check_dtmax_for_decrease
! use energies,         only:etot,totmom,angtot,mdust,np_cs_eq_0,np_e_eq_0
! use dim,              only:maxvxyzu,mhd,periodic
! use fileutils,        only:getnextfilename
! use options,          only:nfulldump,twallmax,nmaxdumps,rhofinal1,use_dustfrac,iexternalforce,&
!                            icooling,ieos,ipdv_heating,ishock_heating,iresistive_heating
! use step_lf_global,   only:step
! use timing,           only:get_timings,print_time,timer,reset_timer,increment_timer,&
!                            setup_timers,timers,reduce_timers,ntimers,&
!                            itimer_fromstart,itimer_lastdump,itimer_step,itimer_io,itimer_ev
! use mpiutils,         only:reduce_mpi,reduceall_mpi,barrier_mpi,bcast_mpi
!#ifdef SORT
! use sort_particles,   only:sort_part
!#endif
!#ifdef IND_TIMESTEPS
! use dim,              only:maxp
! use part,             only:maxphase,ibin,iphase
! use timestep_ind,     only:istepfrac,nbinmax,set_active_particles,update_time_per_bin,&
!                            write_binsummary,change_nbinmax,nactive,nactivetot,maxbins,&
!                            print_dtlog_ind,get_newbin
! use timestep,         only:dtdiff,C_cool
! use timestep_sts,     only:sts_get_dtau_next,sts_init_step
! use step_lf_global,   only:init_step
!#else
! use timestep,         only:dtforce,dtcourant,dterr,print_dtlog
!#endif
! use timestep_sts,     only: use_sts
! use supertimestep,    only: step_sts
!#ifdef DRIVING
! use forcing,          only:write_forcingdump
!#endif
!#ifdef CORRECT_BULK_MOTION
! use centreofmass,     only:correct_bulk_motion
!#endif
!#ifdef MPI
! use part,             only:ideadhead,shuffle_part
!#endif
!#ifdef IND_TIMESTEPS
! use part,             only:twas
! use timestep_ind,     only:get_dt
!#endif
! use part,             only:npart,nptmass,xyzh,vxyzu,fxyzu,fext,divcurlv,massoftype, &
!                            xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gravity,iboundary,npartoftype, &
!                            fxyz_ptmass_sinksink,ntot,poten,ndustsmall
! use quitdump,         only:quit
! use ptmass,           only:icreate_sinks,ptmass_create,ipart_rhomax,pt_write_sinkev
! use io_summary,       only:iosum_nreal,summary_counter,summary_printout,summary_printnow
! use externalforces,   only:iext_spiral
!#ifdef MFLOW
! use mf_write,         only:mflow_write
!#endif
!#ifdef VMFLOW
! use mf_write,         only:vmflow_write
!#endif
!#ifdef BINPOS
! use mf_write,         only:binpos_write
!#endif
!
! implicit none
!
! real            :: dtnew,dtlast,timecheck,rhomaxold,dtmax_log_dratio
! real            :: tprint,dtmaxold,dtinject
! real(kind=4)    :: t1,t2,tcpu1,tcpu2,tstart,tcpustart
! real(kind=4)    :: twalllast,tcpulast,twallperdump,twallused
!#ifdef IND_TIMESTEPS
! integer         :: i,nalive,inbin,iamtypei
! integer(kind=1) :: nbinmaxprev
! integer(kind=8) :: nmovedtot,nalivetot
! real            :: fracactive,speedup,tcheck,dtau,efficiency,tbegin
! real, intent(in) :: tlast
! real(kind=4)    :: tall
! real(kind=4)    :: timeperbin(0:maxbins)
! logical         :: dt_changed
! integer         :: iloop,npart_old
!#else
! real            :: dtprint
! integer         :: nactive
! logical, parameter :: dt_changed = .false.
!#endif
! logical         :: use_global_dt
! integer         :: j,nskip,nskipped,nevwrite_threshold,nskipped_sink,nsinkwrite_threshold
! logical         :: fulldump,abortrun,abortrun_bdy,at_dump_time,writedump
!
!
! --------------------- main loop ----------------------------------------
!
 !tbegin = time
 tcheck = time
 npart_old = npart

 !timestepping: do while ((time < tmax).and.((nsteps < nmax) .or.  (nmax < 0)).and.(rhomaxnow*rhofinal1 < 1.0))
 !write(*,*) "is istepfrac (",istepfrac,") smaller than 2**nbinmax (",2**nbinmax, ")?"
 timesubstepping: do while (istepfrac < 2**nbinmax)
 !write(*,*) "istepfrac (",istepfrac,") is smaller than 2**nbinmax (",2**nbinmax, "), continuing"
 
#ifdef INJECT_PARTICLES
    !
    ! injection of new particles into simulation
    !
    !if (.not. present(flag)) then
       npart_old=npart
       !write(*,*) "INJECTING"
       call inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
       call update_injected_particles(npart_old,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject)
    !endif
#endif

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

    if (gravity .and. icreate_sinks > 0 .and. ipart_rhomax /= 0) then
       !
       ! creation of new sink particles
       !
       call ptmass_create(nptmass,npart,ipart_rhomax,xyzh,vxyzu,fxyzu,fext,divcurlv,&
                          poten,massoftype,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,time)
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
    write(*,*) "dtlast: ", dtlast

    !--timings for step call
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_step,t2-t1,tcpu2-tcpu1)
    call summary_counter(iosum_nreal,t2-t1)

#ifdef IND_TIMESTEPS
    tcheck = tcheck + dt

    !--update time in way that is free of round-off errors
    time = tlast + istepfrac/real(2**nbinmaxprev)*dtmaxold
    write(*,*) "new time: ", time, tlast, istepfrac, nbinmaxprev, dtmaxold

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
!    at_dump_time = (time >= tmax) &
!                   .or.((nsteps >= nmax).and.(nmax >= 0)).or.(rhomaxnow*rhofinal1 >= 1.0)
!#ifdef IND_TIMESTEPS
!    if (istepfrac==2**nbinmax) at_dump_time = .true.
!#else
!    if (time >= tprint) at_dump_time = .true.
!#endif
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
       ! We don't want to write the evfile, but we do want to calculate the energies
       !call write_evfile(time,dt)
       call compute_energies(time)
       if (should_conserve_momentum) call check_conservation_error(totmom,totmom_in,1.e-1,'linear momentum')
       if (should_conserve_angmom)   call check_conservation_error(angtot,angtot_in,1.e-1,'angular momentum')
       if (should_conserve_energy)   call check_conservation_error(etot,etot_in,1.e-1,'energy')
       if (should_conserve_dustmass) then
          do j = 1,ndustsmall
             call check_conservation_error(mdust(j),mdust_in(j),1.e-1,'dust mass',decrease=.true.)
          enddo
       endif
       if (mhd) call check_magnetic_stability(hdivBB_xa)
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
       !call pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink)
#ifdef IND_TIMESTEPS
       dt_changed = .false.
#endif
 write(*,*) "\n\n\n***********tlast: ", tlast, "\n\n\n"
    endif
!
!--write to data file if time is right
!
!    if (at_dump_time) then
!#ifndef IND_TIMESTEPS
!!
!!--Global timesteps: Decrease dtmax if requested (done in step for individual timesteps)
!       twallperdump = timers(itimer_lastdump)%wall
!       call check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_log_dratio,&
!                                     rhomaxold,rhomaxnow,nfulldump,use_global_dt)
!       dt = min(dt,dtmax) ! required if decreasing dtmax to ensure that the physically motivated timestep is not too long
!#endif
!
!       !--modify evfile and logfile names with new number
!       if ((nout <= 0) .or. (mod(noutput,nout)==0)) then
!          if (noutput==1) then
!             evfile  = getnextfilename(evfile)
!             logfile = getnextfilename(logfile)
!          endif
!!         Update values for restart dumps
!          if (dtmax_ifactorWT ==0) then
!             idtmax_n_next    =  idtmax_n
!             idtmax_frac_next =  idtmax_frac
!          elseif (dtmax_ifactorWT > 0) then
!             idtmax_n_next    =  idtmax_n   *dtmax_ifactorWT
!             idtmax_frac_next =  idtmax_frac*dtmax_ifactorWT
!          elseif (dtmax_ifactorWT < 0) then
!             idtmax_n_next    = -idtmax_n   /dtmax_ifactorWT
!             idtmax_frac_next = -idtmax_frac/dtmax_ifactorWT
!          endif
!          idtmax_frac_next = idtmax_frac_next + 1
!          idtmax_frac_next = mod(idtmax_frac_next,idtmax_n_next)
!          dumpfile_orig = trim(dumpfile)
!          if (idtmax_frac==0) then
!             dumpfile = getnextfilename(dumpfile,idumpfile)
!             dumpfile_orig = trim(dumpfile)
!          else
!             write(dumpfile,'(2a)') dumpfile(:index(dumpfile,'_')-1),'.restart'
!          endif
!          writedump = .true.
!       else
!          writedump = .false.
!       endif
!
!       !--do not dump dead particles into dump files
!       if (ideadhead > 0) call shuffle_part(npart)
!!
!!--get timings since last dump and overall code scaling
!!  (get these before writing the dump so we can check whether or not we
!!   need to write a full dump based on the wall time;
!!   move timer_lastdump outside at_dump_time block so that dtmax can
!!   be reduced it too long between dumps)
!!
!       call increment_timer(itimer_fromstart,t2-tstart,tcpu2-tcpustart)
!
!       fulldump = (nout <= 0 .and. mod(noutput,nfulldump)==0) .or. (mod(noutput,nout*nfulldump)==0)
!!
!!--if max wall time is set (> 1 sec) stop the run at the last full dump
!!  that will fit into the walltime constraint, based on the wall time between
!!  the last two dumps added to the current total walltime used.  The factor of three for
!!  changing to full dumps is to account for the possibility that the next step will take longer.
!!  If we are about to write a small dump but it looks like we won't make the next dump,
!!  write a full dump instead and stop the run
!!
!       abortrun = .false.
!       if (twallmax > 1.) then
!          twallused    = timers(itimer_fromstart)%wall
!          twallperdump = timers(itimer_lastdump)%wall
!          if (fulldump) then
!             if ((twallused + abs(nfulldump)*twallperdump) > twallmax) then
!                abortrun = .true.
!             endif
!          else
!             if ((twallused + 3.0*twallperdump) > twallmax) then
!                fulldump = .true.
!                if (id==master) write(iprint,"(1x,a)") '>> PROMOTING DUMP TO FULL DUMP BASED ON WALL TIME CONSTRAINTS... '
!                nfulldump = 1  !  also set all future dumps to be full dumps (otherwise gets confusing)
!                if ((twallused + twallperdump) > twallmax) abortrun = .true.
!             endif
!          endif
!       endif
!!
!!--Promote to full dump if this is the final dump
!!
!       if ( (time >= tmax) .or. ( (nmax > 0) .and. (nsteps >= nmax) ) ) fulldump = .true.
!!
!!--flush any buffered warnings to the log file
!!
!       if (id==master) call flush_warnings()
!!
!!--write dump file
!!
!       if (rkill > 0) call accrete_particles_outside_sphere(rkill)
!#ifndef INJECT_PARTICLES
!       call calculate_mdot(nptmass,time,xyzmh_ptmass)
!#endif
!       call get_timings(t1,tcpu1)
!       if (writedump) then
!          if (fulldump) then
!             call write_fulldump(time,dumpfile)
!             if (id==master) then
!                call write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
!#ifdef DRIVING
!                call write_forcingdump(time,dumpfile)
!#endif
!             endif
!             ncount_fulldumps = ncount_fulldumps + 1
!          else
!             call write_smalldump(time,dumpfile)
!          endif
!       endif
!       call get_timings(t2,tcpu2)
!       call increment_timer(itimer_io,t2-t1,tcpu2-tcpu1)
!
!#ifdef LIVE_ANALYSIS
!       if (id==master .and. idtmax_frac==0) then
!          call do_analysis(dumpfile,numfromfile(dumpfile),xyzh,vxyzu, &
!                           massoftype(igas),npart,time,ianalysis)
!       endif
!#endif
!       call reduce_timers
!       if (id==master) then
!          call print_timinginfo(iprint,nsteps,nsteplast)
!          !--Write out summary to log file
!          call summary_printout(iprint,nptmass)
!       endif
!#ifdef IND_TIMESTEPS
!       !--print summary of timestep bins
!       if (iverbose >= 0) then
!          call write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)
!          timeperbin(:) = 0.
!          if (id==master) call print_dtind_efficiency(iverbose,nalivetot,nmovedtot,tall,timers(itimer_lastdump)%wall,2)
!       endif
!       tlast = tprint
!       istepfrac = 0
!       nmovedtot = 0
!#endif
!       !--print summary of energies and other useful values to the log file
!       if (id==master) call write_evlog(iprint)
!
!#ifndef IND_TIMESTEPS
!       !--Implement dynamic boundaries (for global timestepping)
!       if (dynamic_bdy) call update_boundaries(nactive,nactive,npart,abortrun_bdy)
!#endif
!
!       !
!       !--if twallmax > 1s stop the run at the last full dump that will fit into the walltime constraint,
!       !  based on the wall time between the last two dumps added to the current total walltime used.
!       !
!       if (abortrun) then
!          if (id==master) then
!             call print_time(t2-tstart,'>> WALL TIME = ',iprint)
!             call print_time(twallmax,'>> NEXT DUMP WILL TRIP OVER MAX WALL TIME: ',iprint)
!             write(iprint,"(1x,a)") '>> ABORTING... '
!          endif
!          return
!       endif
!
!       if (abortrun_bdy) then
!          write(iprint,"(1x,a)") 'Will likely surpass maxp_hard next time we need to add particles.'
!          write(iprint,"(1x,a)") 'Recompile with larger maxp_hard.'
!          write(iprint,"(1x,a)") '>> ABORTING... '
!          return
!       endif
!
!       if (nmaxdumps > 0 .and. ncount_fulldumps >= nmaxdumps) then
!          if (id==master) write(iprint,"(a)") '>> reached maximum number of full dumps as specified in input file, stopping...'
!          return
!       endif
!
!       twalllast = t2
!       tcpulast = tcpu2
!       do i = 1,ntimers
!          call reset_timer(i)
!       enddo
!
!       if (idtmax_frac==0) then
!          noutput    = noutput + 1           ! required to determine frequency of full dumps
!       endif
!       noutput_dtmax = noutput_dtmax + 1     ! required to adjust tprint; will account for varying dtmax
!       idtmax_n      = idtmax_n_next
!       idtmax_frac   = idtmax_frac_next
!       tprint        = tzero + noutput_dtmax*dtmaxold
!       nsteplast     = nsteps
!       dumpfile      = trim(dumpfile_orig)
!       if (dtmax_ifactor/=0) then
!          tzero           = tprint - dtmaxold
!          tprint          = tzero  + dtmax
!          noutput_dtmax   = 1
!          dtmax_ifactor   = 0
!          dtmax_ifactorWT = 0
!       endif
!    endif

#ifdef CORRECT_BULK_MOTION
    call correct_bulk_motion()
#endif

    if (iverbose >= 1 .and. id==master) write(iprint,*)
    call flush(iprint)
    !--Write out log file prematurely (if requested based upon nstep, walltime)
    if ( summary_printnow() ) call summary_printout(iprint,nptmass)

 enddo timesubstepping

endif
!end subroutine amuse_new_step
end subroutine amuse_evol


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
    !abundance(:,i) = 0.
    !abundance(ih2ratio,i) = 0.5
    !abundance(iHI,i) = 1.  ! assume all gas is atomic hydrogen initially
    if (i == 1) then
        print*, "xyz vxyz u mass = ", x, y, z, vx, vy, vz, u, mass
        print*, "udist, utime, umass = ", udist, utime, umass
        print*, "x vx u mass in cgs = ", x*udist, vx*udist/utime, u*udist*udist/utime/utime, mass*umass
    endif

#ifdef IND_TIMESTEPS
    dtinject = C_cour * h / (gamma*(gamma-1)*u)**0.5
    nbinmaxprev = nbinmax
    call get_newbin(dtinject,dtmax,nbinmax,allow_decrease=.false.)
    !! not doing clever stuff: all particles go in the shortest possible bin.
    !! FIXME rethink this later...
    nbinmax = 3
    !if (nbinmax > nbinmaxprev) then ! update number of bins if needed
    !   call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
    !   print*, "nbinmax (prev), time: ", nbinmax, nbinmaxprev, time
    !   print*, "npart:", npart
    !endif
    ! put all injected particles on shortest bin
    ibin(i) = nbinmax
    twas(i) = time + 0.5*get_dt(dtmax,ibin(i))
    if (i == 5) then
        print*, "particle ", i, ": "
        print*, x, y, z, vx, vy, vz, mass, ibin(i), twas(i)
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
    !udist = unit_length_in
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
    !umass = unit_mass_in
    !call set_units(mass=umass,time=utime,G=1.)
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
    !utime = unit_time_in
    !call set_units(time=utime,mass=umass,G=1.)
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
    !nodisabled = .true.
    nodisabled = .false.
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
    real :: tempi
    eos_type = ieos
    pmassi = massoftype(abs(iphase(i)))
    call amuse_get_density(i, rho)
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    call equationofstate(eos_type,ponrho,spsound,rho,x,y,z,tempi)
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
    write(*,*) 'getting dm ', i
end subroutine amuse_get_state_dm

subroutine amuse_get_state_sink(i, mass, x, y, z, vx, vy, vz, radius)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz, radius
    call amuse_get_mass(i, mass)
    call amuse_get_position(i, x, y, z)
    call amuse_get_velocity(i, vx, vy, vz)
    call amuse_get_sink_radius(i, radius)
    write(*,*) 'getting sink ', i, ': radius is ', radius
end subroutine amuse_get_state_sink

subroutine amuse_get_sink_radius(j, radius)
    use part, only:xyzmh_ptmass, ihacc
    implicit none
    integer :: j
    double precision :: radius
    radius = xyzmh_ptmass(ihacc, -j)
end subroutine amuse_get_sink_radius

subroutine amuse_get_sink_temperature(j, temperature)
    use part, only:xyzmh_ptmass, iTeff
    implicit none
    integer, intent(in) :: j
    double precision :: temperature
    temperature = xyzmh_ptmass(iTeff, -j)
end subroutine

subroutine amuse_get_sink_luminosity(j, luminosity)
    use part, only:xyzmh_ptmass, iLum
    implicit none
    integer, intent(in) :: j
    double precision :: luminosity
    luminosity = xyzmh_ptmass(iLum, -j)
end subroutine

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

subroutine amuse_get_acceleration(i, fx, fy, fz)
    use part, only:fxyzu,fxyz_ptmass
    implicit none
    integer, intent(in) :: i
    integer :: j
    double precision, intent(out) :: fx, fy, fz
    if (i == abs(i)) then
        fx = fxyzu(1, i)
        fy = fxyzu(2, i)
        fz = fxyzu(3, i)
    else
        j = -i
        fx = fxyz_ptmass(1, j)
        fy = fxyz_ptmass(2, j)
        fz = fxyz_ptmass(3, j)
    endif
end subroutine amuse_get_acceleration

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
        write(*,*) "not a sink"
        call amuse_get_smoothing_length(i, radius)
    else
        write(*,*) "a sink"
        j = -i
        call amuse_get_sink_radius(i, radius)
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

subroutine amuse_set_state_sink(i, mass, x, y, z, vx, vy, vz, radius)
    implicit none
    integer :: i
    double precision :: mass, x, y, z, vx, vy, vz, radius
    call amuse_set_mass(i, mass)
    call amuse_set_position(i, x, y, z)
    call amuse_set_velocity(i, vx, vy, vz)
    call amuse_set_sink_radius(i, radius)
end subroutine

subroutine amuse_set_sink_radius(j, radius)
    use part, only:xyzmh_ptmass, ihacc
    implicit none
    integer :: j
    double precision :: radius
    xyzmh_ptmass(ihacc, -j) = radius
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

subroutine amuse_set_sink_temperature(j, temperature)
    use part, only:xyzmh_ptmass, iTeff
    implicit none
    integer, intent(in) :: j
    double precision :: temperature
    xyzmh_ptmass(iTeff, -j) = temperature
end subroutine

subroutine amuse_set_sink_luminosity(j, luminosity)
    use part, only:xyzmh_ptmass, iLum
    implicit none
    integer, intent(in) :: j
    double precision :: luminosity
    xyzmh_ptmass(iLum, -j) = luminosity
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
#ifdef INJECT_PARTICLES
    use inject,           only:inject_particles
    use part,             only:npartoftype
    use partinject,       only:update_injected_particles
#endif
    use options, only:rhofinal1
    use ptmass, only:rho_crit
    use part, only:npart
    use step_lf_global, only:init_step

    use part, only: xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass
    implicit none
    logical :: amuse_initialise
    double precision, intent(in) :: tmax_in
    logical :: maximum_density_reached
    real :: tlast
    real :: dtinject,dtlast
    integer(kind=1) :: nbinmax
#ifdef INJECT_PARTICLES
    integer         :: npart_old
#endif
#ifndef IND_TIMESTEPS
    integer :: istepfrac
    istepfrac = 0 ! dummy values
#endif
    dtinject  = huge(dtinject)
    dtlast = 0
    nbinmax = 0
    
    tmax = tmax_in ! - epsilon(tmax_in)
    !dtmax = (tmax - time)
    
    tlast = time
    write(*,*) "TIMESTEPPING: evolve from ", time, " to ", tmax
    timestepping: do while (time < tmax)

#ifdef IND_TIMESTEPS
        istepfrac = 0
        !print*, "*****init timestep"
        !call init_step(npart, time, dtmax)
#endif
        !print*, "*****new_step - Time; tmax: ", time, tmax
        !call amuse_new_step(tlast)
        amuse_initialise = .false.
        call amuse_evol(amuse_initialise)
        !print*, "*****new_step done - Time; tmax: ", time, tmax
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
    use dim, only:h2chemistry
    use chem, only:init_chem
    use cooling, only:init_cooling,Tfloor
    !use cooling_ism, only:init_cooling
    implicit none
    integer :: ierr
    integer, intent(in) :: icooling_in
    icooling = icooling_in
    if (icooling > 0) then
        Tfloor = 1  ! K
        call init_cooling(id,master,iprint,ierr)
        !if (h2chemistry) then
        !    if (id==master) write(iprint,*) 'initialising cooling function...'
        !    call init_chem()
        !    call init_h2cooling()
        !else
        !    call init_cooling(ierr)
        !    if (ierr /= 0) call fatal('initial','error initialising cooling')
        !endif
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

subroutine amuse_inject()
    use timestep, only:time
    use part, only:xyzh,vxyzu,npart,npartoftype
    use ptmass, only:xyzmh_ptmass,vxyz_ptmass,
    use inject, only:inject_particles
    implicit none
    real :: dtinject, dtlast
    dtlast = 0
    dtinject = huge(dtinject)
    
    ! if npart > 0, this will also delete 'too far away particles'
    ! also used to determine number of already released shells
    ! time and dtlast are used together to determine 'time passed'
    ! 
    call inject_particles(&
        time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,dtinject&
    )
end subroutine

subroutine amuse_update_injected()
    use partinject, only:update_injected_particles
    use part, only:npart
    use timestep_ind, only:istepfrac,nbinmax
    use timestep, only:time,dtmax,dt
    implicit none
    integer :: npart_old
    call update_injected_particles(&
        npart_old,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject&
    )
end subroutine
