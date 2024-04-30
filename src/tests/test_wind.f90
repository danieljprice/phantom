!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testwind
!
! Unit tests of the wind injection module
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checksetup, dim, eos, inject, io, options, part,
!   partinject, physcon, step_lf_global, testutils, timestep, timestep_ind,
!   units, wind
!
 implicit none
 public :: test_wind

 private

contains
!----------------------------------------------------------
!+
!  Unit tests of timestepping and boundary crossing
!+
!----------------------------------------------------------
subroutine test_wind(ntests,npass)
 use io,        only:iprint,id,master,iverbose!,iwritein
 use boundary,  only:set_boundary
 use options,   only:ieos!,icooling
 use physcon,   only:au,solarm,solarl
 use units,     only:umass,set_units,utime,unit_energ,udist
 use inject,    only:init_inject,inject_particles,set_default_options_inject,inject_type
 use eos,       only:gmw,ieos,init_eos,gamma,polyk
 use part,      only:npart,init_part,nptmass,xyzmh_ptmass,vxyz_ptmass,xyzh,vxyzu,&
                     nptmass,npartoftype,igas,iTeff,iLum,iReff,massoftype,ntot
 use timestep,  only:time,tmax,dt,dtmax,nsteps,dtrad,dtforce,dtcourant,dterr,print_dtlog
 use step_lf_global, only:step,init_step
 use testutils,      only:checkval,update_test_scores
 use dim,            only:isothermal,inject_parts,mpi
 use partinject,     only:update_injected_particles
 use timestep_ind,   only:nbinmax
 use wind,           only:trvurho_1D
 use checksetup,     only:check_setup
 !use readwrite_infile, only:read_infile,write_infile

 integer, intent(inout) :: ntests,npass

 integer :: i,ierr,nerror,istepfrac,npart_old,nfailed(9),nwarn
 real :: dtinject,dtlast,t,default_particle_mass,dtext,dtnew,dtprint,dtmaxold,tprint

 if (mpi) then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING WIND TEST (currently not working with MPI)'
    return
 elseif (inject_type /= 'wind') then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING WIND TEST (need to compile with wind injection module)'
    return
 else
    if (id==master) write(*,"(/,a,/)") '--> TESTING WIND MODULE'
 endif

 call set_units(dist=au,mass=solarm,G=1.d0)
 call set_boundary(-50.,50.,-50.,50.,-50.,50.)

 call init_part()

 ! set properties of mass-losing sink particle
 nptmass = 1
 xyzmh_ptmass(4,1)     = 1.2*solarm/umass
 xyzmh_ptmass(5,1)     = au/udist
 xyzmh_ptmass(iTeff,1) = 50000.
 xyzmh_ptmass(iReff,1) = au/udist
 xyzmh_ptmass(iLum,1)  = 2e4 *solarl * utime / unit_energ

 !
 ! for binary wind simulations the particle mass is IRRELEVANT
 ! since it will be over-written on the first call to init_inject
 !
 default_particle_mass = 1.e-11
 massoftype(igas) = default_particle_mass * (solarm / umass)

 gmw = 2.381
 if (isothermal) then
    gamma = 1.
    ieos = 1
 else
    gamma = 1.4
    ieos = 2
 endif
 polyk = 0.

 call init_eos(ieos,ierr)

 iverbose = 0
 !icooling = 0
 dtmax    = 1.
 tmax     = 8.
 tprint   = tmax
 dt       = 0.
 dtinject = huge(dtinject)
 dtrad    = huge(dtrad)
 t        = 0.
 dtnew    = 0.

 ! trans-sonic wind
 call set_default_options_inject(1)
 call check_setup(nerror,nwarn)

 istepfrac  = 0
 nfailed(:) = 0
 call init_inject(nerror)

 !debug if (id==master) call write_infile('w.in','w.log','w.ev','w_00000',iwritein,iprint)

 ! check particle's mass
 call checkval(massoftype(igas),1.490822861042279E-09,epsilon(0.),&
      nfailed(1),'no errors in setting particle mass')
 npart_old = npart
 call inject_particles(t,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                       npart,npart_old,npartoftype,dtinject)
 call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)

 ! check 1D wind profile
 i = size(trvurho_1D(1,:))
 !print '((5(1x,es22.15)))',trvurho_1D(:,i),massoftype(igas)
 call checkval(trvurho_1D(2,i),7.058624412798283E+13,epsilon(0.),nfailed(2),'outer wind radius')
 call checkval(trvurho_1D(3,i),1.112160584479353E+06,epsilon(0.),nfailed(3),'outer wind velocity')
 call checkval(trvurho_1D(4,i),2.031820842001706E+12,epsilon(0.),nfailed(4),'outer wind internal energy')
 call checkval(trvurho_1D(5,i),8.878887149408118E-15,epsilon(0.),nfailed(5),'outer wind density')

 dt = dtinject
 dtlast = 0.
 time = 0.

 call init_step(npart_old,time,dtmax)

 do while (t < tmax)

    dtext = dt
    !
    ! injection of new particles into simulation
    !
    npart_old=npart
    call inject_particles(t,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)
    dtmaxold = dtmax
    nsteps = nsteps+1

    call step(npart,npart,t,dt,dtext,dtnew)

    t = t + dt
    time = t

    dtlast = dt
    dtprint = min(tprint,tmax) - t + epsilon(dtmax)
    if (dtprint <= epsilon(dtmax) .or. dtprint >= (1.0-1e-8)*dtmax ) dtprint = dtmax + epsilon(dtmax)
    dt = min(dtforce,dtcourant,dterr,dtmax+epsilon(dtmax),dtprint,dtinject,dtrad)
    if (id==master) call print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,dtrad,dtprint,dtinject,ntot)

 enddo

 !print '((3(1x,es22.15),i8))',xyzmh_ptmass(4,1),xyzmh_ptmass(7,1),xyzmh_ptmass(15,1),npart
 call checkval(xyzmh_ptmass(4,1),1.199987894518367E+00,epsilon(0.),nfailed(6),'sink particle mass')
 call checkval(xyzmh_ptmass(7,1),0.,epsilon(0.),nfailed(7),'mass accreted')
 call checkval(npart,12180,0,nfailed(8),'number of ejected particles')
 call checkval(xyzmh_ptmass(15,1),1.591640703559762E-06,epsilon(0.),nfailed(9),'wind mass loss rate')
 call update_test_scores(ntests,nfailed,npass)

 if (id==master) write(*,"(/,a)") '<-- WIND TEST COMPLETE'

end subroutine test_wind

end module testwind
