!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Dependencies: allocutils, boundary, checksetup, dim, dust_formation,
!   eos, inject, io, options, part, partinject, physcon, ptmass,
!   ptmass_radiation, readwrite_infile, step_lf_global, testutils,
!   timestep, timestep_ind, units, wind
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
 use io,         only:id,master
 use inject,     only:inject_type
 use boundary,   only:set_boundary
 use physcon,    only:au,solarm,solarl
 use units,      only:set_units,udist
 use part,       only:npart,xyzmh_ptmass,xyzh,vxyzu,dust_temp,iReff
 use dim,        only:mpi,maxTdust,maxp,sink_radiation,nucleation,ind_timesteps,&
                      disc_viscosity,nalpha,use_dust,isothermal
 use allocutils, only:allocate_array
 use options,    only:alpha
 use timestep,   only:tmax,tolv
 use readwrite_dumps, only:write_fulldump
 integer, intent(inout) :: ntests,npass
 integer :: npart_old,istepfrac,npart_prefill
 real    :: dtinject,rmax,mstar0
 logical :: testcyl,use_shock_switch

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

 use_shock_switch = (nalpha >= 0)
 testcyl = .not.sink_radiation .and. .not.nucleation .and. use_shock_switch .and. ind_timesteps
 disc_viscosity = .false.
 if (testcyl) then
    alpha = 1.
    disc_viscosity = .true.
 endif
 tolv = 1.e2
 rmax = 10.*au/udist

 ! test trans-sonic wind - no radiation, no dust
 if (.not.isothermal) then
    if (id==master) write(*,"(/,a,/)") '--> testing transonic wind'
    call init_testwind(1,ntests,npass,npart_old,istepfrac,dtinject,npart_prefill,mstar0,rmax)

    call integrate_wind(npart_old,istepfrac,dtinject)

    ! check injected mass
    call test_injected_mass(ntests,npass,npart,npart_prefill,1,xyzmh_ptmass,mstar0,tmax)

    ! check wind against 1D solution
    call test_against_1D_profile(ntests,npass,npart,xyzh,vxyzu,&
                                 1,xyzmh_ptmass,2.*xyzmh_ptmass(iReff,1),rmax)

    ! uncomment the following line to write a dump file for debugging
    !call write_fulldump(tmax,'testwind_00000',int(npart,kind=8))
 else
    if (id==master) write(*,"(/,a,/)") '--> skipping transonic wind (isothermal)'
 endif

 if (sink_radiation) then

! test wind with bowen dust + radiative acceleration
   if (id==master) write(*,"(/,a,/)") '--> testing bowen dust + radiative acceleration'

    maxTdust = maxp
    if (allocated(dust_temp)) deallocate(dust_temp)
    call allocate_array('dust_temp',dust_temp,maxTdust)

    call init_testwind(2,ntests,npass,npart_old,istepfrac,dtinject,npart_prefill,mstar0,rmax)

    call integrate_wind(npart_old,istepfrac,dtinject)

    ! check injected mass
    call test_injected_mass(ntests,npass,npart,npart_prefill,1,xyzmh_ptmass,mstar0,tmax)

    ! check wind against 1D solution
    !call test_against_1D_profile(ntests,npass,npart,xyzh,vxyzu,&
    !                             1,xyzmh_ptmass,2.*xyzmh_ptmass(iReff,1),rmax)

    ! uncomment the following line to write a dump file for debugging
    !call write_fulldump(tmax,'testbowen_00000',int(npart,kind=8))
 else
    if (id==master) write(*,"(/,a,/)") '--> skipping sink radiation test'
 endif

 if (id==master) write(*,"(/,a)") '<-- WIND TEST COMPLETE'

end subroutine test_wind

!-----------------------------------------------------------------------
!+
!  Initialises the wind test, and perform checks on the wind profile
!+
!-----------------------------------------------------------------------
subroutine init_testwind(icase,ntests,npass,npart_old,istepfrac,dtinject,npart_prefill,mstar0,rmax)
 use io,         only:iverbose,error,id,master
 use inject,     only:init_inject,inject_particles,set_default_options_inject
 use units,      only:umass,udist,unit_mdot,unit_velocity,unit_luminosity,utime
 use physcon,    only:au,solarm,solarl,km,seconds,years,pi,gg
 use eos,        only:gmw,ieos,init_eos,gamma,polyk
 use part,       only:npart,init_part,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,xyzh,vxyzu,&
                     npartoftype,igas,iTeff,iLum,iReff,massoftype,iTwind,ivwind,imloss
 use timestep,   only:tmax,dt,dtmax,dtrad
 use dim,        only:isothermal
 use wind,       only:trvurho_1D,rfill_domain_au
 use ptmass,     only:set_integration_precision
 use testutils,  only:checkval,checkvalbuf,checkvalbuf_end,update_test_scores
 use checksetup, only:check_setup
 use partinject, only:update_injected_particles
 use timestep_ind,     only:nbinmax
 use ptmass_radiation, only:alpha_rad,isink_radiation
 use dust_formation,   only:idust_opacity

 integer, intent(in) :: icase
 integer, intent(inout) :: ntests,npass
 integer, intent(out) :: npart_old,istepfrac,npart_prefill
 real,    intent(out) :: dtinject,mstar0
 real,    intent(in) :: rmax

 integer :: i,ierr,nerror,nwarn,ngrid
 integer :: nfailed(2),ncheck(2)
 real :: errmax(2)
 real :: t,default_particle_mass,dtnew
 real :: mdot0,mdot,mstar,r,v,rho,u,e,e0
 real, parameter :: tol_e = 2.e-4, tol_mdot = 5.e-16

 call init_part()
 call set_integration_precision()

 ! set properties of mass-losing sink particle
 nptmass = 1
 xyzmh_ptmass(:,:) = 0.
 xyzmh_ptmass(4,1)  = 1.2*solarm/umass
 xyzmh_ptmass(5,1)  = au/udist
 if (icase == 1) then      !trans-sonic wind
    xyzmh_ptmass(iTeff,1) = 50000.
    xyzmh_ptmass(ivwind,1) = 0.
    xyzmh_ptmass(iTwind,1) = 50000.
    !xyzmh_ptmass(imloss,1) = 1.d-5*(solarm/years)/unit_mdot
    xyzmh_ptmass(imloss,1) = 1.d-5*(solarm/umass)/(years/utime)
 elseif (icase == 2) then      !super sonic-wind
    xyzmh_ptmass(iTeff,1) = 3000.
    xyzmh_ptmass(ivwind,1) = 20.*km/seconds/unit_velocity
    xyzmh_ptmass(iTwind,1) = 2500.
    xyzmh_ptmass(imloss,1) = 1.d-5*(solarm/years)/unit_mdot
 endif
 xyzmh_ptmass(iReff,1) = 2.*au/udist
 xyzmh_ptmass(iLum,1)  = 2e4 *solarl / unit_luminosity
 vxyz_ptmass = 0.
 fxyz_ptmass = 0.
 !
 ! for binary wind simulations the particle mass is IRRELEVANT
 ! since it will be over-written on the first call to init_inject
 !
 default_particle_mass = 1.e-9
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
 dtmax = 1.
 tmax  = 12.
 !wind + bowen dust + radiation force
 if (icase == 1) then
    alpha_rad = 0.
    isink_radiation = 0 !radiation + alpha_rad
    idust_opacity = 0   !bowen opacity
 elseif (icase == 2) then
    alpha_rad = 1.
    isink_radiation = 3 !radiation + alpha_rad
    idust_opacity = 1   !bowen opacity
 else
    call error('test_wind','unknown test case ',ival=icase)
    return
 endif
 dt       = 0.
 dtinject = huge(dtinject)
 dtrad    = huge(dtrad)
 t        = 0.
 dtnew    = 0.
 mstar0   = xyzmh_ptmass(4,1)

 ! trans-sonic wind
 call set_default_options_inject(icase)
 call check_setup(nerror,nwarn)

 ! set how much of the domain to pre-fill
 rfill_domain_au = rmax*udist/au  ! pre-fill out to rmax

 nfailed(:) = 0; ncheck(:) = 0; errmax(:) = 0.
 istepfrac  = 0
 call init_inject(nerror)

 npart_old = npart
 ! the tests below work for both trans-sonic wind and radiation-driven wind
 if (icase == 1 .or. icase == 2) then
    ! pre-fill domain with boundary + fill shells (like setup_bondiinject): dtlast=0 triggers init path
    call inject_particles(t,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                          npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)

    npart_prefill = npart
    xyzmh_ptmass(4,1) = mstar0  ! reset star mass as we pre-fill twice in order to get npart_prefill correct

    ! check 1D wind profile satisfies Mdot = 4 pi r^2 rho v = constant
    ngrid = size(trvurho_1D,dim=2)
    if (id==master) write(*,"(/,a)") '--> checking 1D wind profile'
    mdot0 = xyzmh_ptmass(imloss,1)
    mstar = xyzmh_ptmass(4,1)*umass
    do i=1,ngrid
       r = trvurho_1D(2,i)
       v = trvurho_1D(3,i)
       rho = trvurho_1D(5,i)
       u = trvurho_1D(4,i)

       ! check that the mass flux is constant at every radius in the 1D wind profile
       mdot = 4.*pi*r**2*rho*v/unit_mdot
       call checkvalbuf(mdot,mdot0,tol_mdot,'dM/dt = 4 pi r^2 rho v',nfailed(1),ncheck(1),errmax(1))

       ! check Bernoulli energy is constant at every radius in the 1D wind profile
       if (isothermal) then
          e = 0.5*v**2 - (1.-alpha_rad)*gg*mstar/r + polyk*log(rho)
       else
          e = 0.5*v**2 - (1.-alpha_rad)*gg*mstar/r + gamma*u
       endif
       e = 0.5*v**2 - (1.-alpha_rad)*gg*mstar/r + gamma*u
       if (i == 1) e0 = e
       call checkvalbuf(e,e0,tol_e,'Bernoulli constant',nfailed(2),ncheck(2),errmax(2))
    enddo
    call checkvalbuf_end('dM/dt = 4 pi r^2 rho v',ncheck(1),nfailed(1),errmax(1),tol_mdot)
    call checkvalbuf_end('Bernoulli constant',ncheck(2),nfailed(2),errmax(2),tol_e)
    call update_test_scores(ntests,nfailed,npass)
 endif

end subroutine init_testwind

!-----------------------------------------------------------------------
!+
!  checks the wind profile against the 1D solution
!+
!-----------------------------------------------------------------------
subroutine test_against_1D_profile(ntests,npass,npart,xyzh,vxyzu,isink,xyzmh_ptmass,rmin,rmax)
 use io,        only:id,master
 use wind,      only:interp_wind_profile_at_r
 use part,      only:rhoh,massoftype,isdead_or_accreted,igas
 use units,     only:udist,unit_density,unit_ergg,unit_velocity
 use physcon,   only:au
 use testutils, only:checkvalbuf,checkvalbuf_end,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: npart,isink
 real, intent(in) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),rmin,rmax
 integer :: i,nfailed(3),ncheck(3)
 real :: dx(3),r,rhoi,ui,vi,rho,u,v,errmax(3)
 real, parameter :: tol_v = 1.3e-1, tol_u = 9.e-2, tol_rho = 5.e-16

 if (id==master) write(*,"(/,a,2(f7.2,a))") &
    '--> checking wind profile against 1D for r between ',rmin*udist/au,' and ',rmax*udist/au,' au'

 nfailed(:) = 0; ncheck(:) = 0; errmax(:) = 0.
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       dx = xyzh(1:3,i)-xyzmh_ptmass(1:3,isink)
       r = sqrt(dot_product(dx,dx))
       if (r > rmin .and. r < rmax) then
          rhoi = rhoh(xyzh(4,i),massoftype(igas))*unit_density
          vi = dot_product(vxyzu(1:3,i),dx/r)*unit_velocity
          ui = vxyzu(4,i)*unit_ergg
          r = r*udist
          call interp_wind_profile_at_r(r,v,u,rho,isink)
          call checkvalbuf(v,vi,tol_v,'v',nfailed(1),ncheck(1),errmax(1))
          call checkvalbuf(u,ui,tol_u,'u',nfailed(2),ncheck(2),errmax(2))
          call checkvalbuf(rho,rhoi,tol_rho,'rho',nfailed(3),ncheck(3),errmax(3))
       endif
    endif
 enddo

 call checkvalbuf_end('v against 1D wind profile',ncheck(1),nfailed(1),errmax(1),tol_v)
 call checkvalbuf_end('u against 1D wind profile',ncheck(2),nfailed(2),errmax(2),tol_u)
 call checkvalbuf_end('rho against 1D wind profile',ncheck(3),nfailed(3),errmax(3),tol_rho)
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_against_1D_profile

!-----------------------------------------------------------------------
!+
!  checks the injected mass
!+
!-----------------------------------------------------------------------
subroutine test_injected_mass(ntests,npass,npart,npart_prefill,isink,xyzmh_ptmass,mass0,tmax)
 use part,      only:massoftype,igas,ieject,imloss,imacc
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: npart,npart_prefill,isink
 real,    intent(in)    :: xyzmh_ptmass(:,:),mass0,tmax
 integer :: neject,npart_per_sphere
 integer :: nfailed(5)
 real :: mstar,minject,mprefill,tol_mass

 ! injected mass is the total mass of particles minus the boundary shells
 mstar = mass0
 minject = xyzmh_ptmass(imloss,isink)*tmax
 mprefill = npart_prefill*massoftype(igas)
 neject  = nint(minject/massoftype(igas))
 nfailed(:) = 0
 npart_per_sphere = nint(xyzmh_ptmass(ieject,isink))
 tol_mass = npart_per_sphere*massoftype(igas)/minject
 call checkval(xyzmh_ptmass(4,isink),mstar-minject-mprefill,4.e-6,nfailed(1),'sink particle mass')
 call checkval(xyzmh_ptmass(imacc,isink),0.,epsilon(0.),nfailed(2),'mass accreted')
 call checkval(minject,(npart-npart_prefill)*massoftype(igas),tol_mass,nfailed(3),'mass injected')
 call checkval(npart-npart_prefill,neject,npart_per_sphere,nfailed(4),'number of ejected particles')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_injected_mass

!-----------------------------------------------------------------------
!+
!  integrates the wind test
!+
!-----------------------------------------------------------------------
subroutine integrate_wind(npart_old,istepfrac,dtinject)
 use io,        only:id,iprint,master
 use timestep,  only:time,tmax,dt,dtmax,nsteps,dtrad,dtforce,dtcourant,dterr,print_dtlog
 use part,      only:npart,init_part,xyzmh_ptmass,vxyz_ptmass,xyzh,vxyzu,npartoftype,ntot
 use timestep_ind, only:nbinmax
 use step_lf_global, only:step,init_step
 use partinject,     only:update_injected_particles
 use inject,     only:inject_particles

 integer, intent(inout) :: istepfrac,npart_old
 real, intent(inout) :: dtinject

 real :: dtlast,t,dtext,dtnew,dtprint,dtmaxold,tprint

 dt     = dtinject
 dtlast = 0.
 time   = 0.
 tprint = tmax
 t      = 0.

 if (id==master) write(*,"(/,a)") '--> integrating wind'
 call init_step(npart_old,time,dtmax)

 do while (t < tmax)

    dtext = dt
    !
    ! injection of new particles into simulation
    !
    npart_old = npart
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

end subroutine integrate_wind

end module testwind
