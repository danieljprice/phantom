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
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, boundary, checksetup, dim, dust_formation,
!   eos, inject, io, mathfunc, mpiutils, options, part, partinject,
!   physcon, ptmass, ptmass_radiation, step_lf_global, testutils, timestep,
!   timestep_ind, units, wind
!
 use io,        only:id,master,iverbose,error
 use testutils, only:checkval,checkvalbuf,checkvalbuf_end,update_test_scores
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
 use inject,     only:inject_type
 use boundary,   only:set_boundary
 use physcon,    only:au,solarm
 use units,      only:set_units,udist
 use part,       only:npart,xyzmh_ptmass,xyzh,vxyzu,dust_temp,iReff
 use dim,        only:maxTdust,maxp,sink_radiation,nucleation,ind_timesteps,disc_viscosity,&
                      nalpha,isothermal
 use ptmass,     only:set_integration_precision
 use allocutils, only:allocate_array
 use options,    only:alpha
 use timestep,   only:tmax,tolv
 integer, intent(inout) :: ntests,npass
 integer :: npart_old,istepfrac,npart_prefill
 real    :: dtinject,rmax,mstar0
 integer :: i
 integer, parameter :: nshell_test(3) = (/12,252,1000/)
 logical :: testcyl,use_shock_switch

 if (inject_type /= 'wind') then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING WIND TEST (need to compile with wind injection module)'
    return
 else
    if (id==master) write(*,"(/,a,/)") '--> TESTING WIND MODULE'
 endif

 call set_units(dist=au,mass=solarm,G=1.d0)
 call set_boundary(-50.,50.,-50.,50.,-50.,50.)
 call set_integration_precision

 use_shock_switch = (nalpha >= 0)
 testcyl = .not.sink_radiation .and. .not.nucleation .and. use_shock_switch .and. ind_timesteps
 disc_viscosity = .false.
 if (testcyl) then
    alpha = 1.
    disc_viscosity = .true.
 endif
 tolv = 1.e2
 rmax = 10.*au/udist

 ! test the 1D wind solver against the exact solution for an isothermal wind
 call test_parker_wind(ntests,npass,rmax)

 ! winds from several sinks, independent of the order of the particles
 call test_multiple_wind_sinks(ntests,npass)

 ! test trans-sonic wind - no radiation, no dust
 if (.not.isothermal) then
    ! particle mass set from the number of particles per shell must give exactly that number
    do i=1,size(nshell_test)
       if (id==master) write(*,"(/,a,i0,a,/)") '--> testing wind with ',nshell_test(i),' particles per shell'
       call init_testwind(1,ntests,npass,npart_old,istepfrac,dtinject,npart_prefill,mstar0,rmax,nshell_test(i))
    enddo

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
subroutine init_testwind(icase,ntests,npass,npart_old,istepfrac,dtinject,npart_prefill,mstar0,rmax,nshell)
 use inject,     only:init_inject,inject_particles,set_default_options_inject
 use units,      only:umass,udist,unit_mdot,unit_velocity,unit_luminosity,utime
 use physcon,    only:au,solarm,solarl,km,seconds,years,pi,gg
 use eos,        only:gmw,ieos,init_eos,gamma,polyk
 use part,       only:npart,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,xyzh,vxyzu,init_part,&
                      npartoftype,igas,iTeff,iLum,iReff,massoftype,iTwind,ivwind,imloss,iwalpha,rho,&
                      sync_aprmassoftype,ieject
 use timestep,   only:tmax,dt,dtmax,dtrad
 use dim,        only:isothermal
 use wind,       only:wind_profiles,rfill_domain_au,get_wind_particle_mass
 use checksetup, only:check_setup
 use partinject, only:update_injected_particles
 use timestep_ind,     only:nbinmax
 use ptmass_radiation, only:isink_radiation,get_Teff
 use dust_formation,   only:idust_opacity

 integer, intent(in)    :: icase
 integer, intent(inout) :: ntests,npass
 integer, intent(out)   :: npart_old,istepfrac,npart_prefill
 real,    intent(out)   :: dtinject,mstar0
 real,    intent(in)    :: rmax
 integer, intent(in), optional :: nshell

 integer :: i,ierr,nerror,nwarn,ngrid,wind_type
 integer :: nfailed(2),ncheck(2)
 real :: errmax(2)
 real :: t,default_particle_mass,dtnew,mdot_in,pmass_in
 real :: mdot0,mdot,mstar,r,v,rho1d,u,e,e0
 real, parameter :: tol_e = 2.e-4, tol_mdot = 5.e-16

 call init_part

 ! set properties of mass-losing sink particle
 nptmass = 1
 xyzmh_ptmass = 0.
 xyzmh_ptmass(4,1)  = 1.2*solarm/umass
 xyzmh_ptmass(5,1)  = au/udist
 if (icase == 1) then      !trans-sonic wind
    xyzmh_ptmass(iTeff,1) = 50000.
    xyzmh_ptmass(ivwind,1) = 0.
    xyzmh_ptmass(iTwind,1) = 50000.
    !xyzmh_ptmass(imloss,1) = 1.d-5*(solarm/years)/unit_mdot
    xyzmh_ptmass(imloss,1) = 1.d-5*(solarm/umass)/(years/utime)
 elseif (icase == 2) then      !super sonic-wind
    ! stellar radius of 1 au as in Siess et al. (2022), inside the injection radius
    xyzmh_ptmass(iTeff,1) = get_Teff(2.e4,1.)
    xyzmh_ptmass(ivwind,1) = 20.*km/seconds/unit_velocity
    xyzmh_ptmass(iTwind,1) = 2500.
    xyzmh_ptmass(imloss,1) = 1.d-5*(solarm/years)/unit_mdot
 endif
 xyzmh_ptmass(iReff,1) = 2.*au/udist
 xyzmh_ptmass(iLum,1)  = 2e4 *solarl / unit_luminosity
 vxyz_ptmass = 0.
 fxyz_ptmass = 0.
 default_particle_mass = 1.e-9

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
    xyzmh_ptmass(iwalpha,1) = 0.
    isink_radiation = 0 !radiation + alpha_rad
    idust_opacity = 0   !bowen opacity
 elseif (icase == 2) then
    xyzmh_ptmass(iwalpha,1) = 1.
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
 ! wind type set by set_default_options_inject: trans-sonic for icase 1
 wind_type = 0
 if (icase == 1 .and. .not.isothermal) wind_type = 1
 call check_setup(nerror,nwarn)

 ! set how much of the domain to pre-fill (not needed to check particles per shell)
 rfill_domain_au = rmax*udist/au  ! pre-fill out to rmax
 if (present(nshell)) rfill_domain_au = 0.

 ! set the particle mass as phantomsetup does: from the number of particles per shell,
 ! by resolving the sonic point (trans-sonic wind) or directly
 if (present(nshell)) then
    massoftype(igas) = get_wind_particle_mass(xyzmh_ptmass(:,1),wind_type,nshell,0)
 elseif (icase == 1) then
    massoftype(igas) = get_wind_particle_mass(xyzmh_ptmass(:,1),wind_type,0,8)
 else
    massoftype(igas) = default_particle_mass * (solarm / umass)
 endif

 ! check
 nfailed = 0; ncheck = 0; errmax = 0.
 istepfrac  = 0
 mdot_in  = xyzmh_ptmass(imloss,1)
 pmass_in = massoftype(igas)
 call init_inject(nerror)
 call sync_aprmassoftype()

 ! the mass loss rate and particle mass are inputs: init_inject must not change them
 call checkval(xyzmh_ptmass(imloss,1),mdot_in,epsilon(0.),nfailed(1),'Mdot unchanged by init_inject')
 call checkval(massoftype(igas),pmass_in,epsilon(0.),nfailed(2),'particle mass unchanged by init_inject')
 call update_test_scores(ntests,nfailed,npass)
 nfailed = 0

 npart_old = npart
 npart_prefill = npart
 if (present(nshell)) then
    call checkval(nint(xyzmh_ptmass(ieject,1)),nshell,0,nfailed(1),'particles per shell')
    call update_test_scores(ntests,nfailed(1:1),npass)
    return
 endif

 ! the tests below work for both trans-sonic wind and radiation-driven wind
 if (icase == 1 .or. icase == 2) then
    ! pre-fill domain with boundary + fill shells (like setup_bondiinject): dtlast=0 triggers init path
    call inject_particles(t,0.,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,&
                          npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)

    npart_prefill = npart

    ! check 1D wind profile satisfies Mdot = 4 pi r^2 rho v = constant
    ngrid = size(wind_profiles(1)%trvurho,dim=2)
    if (id==master) write(*,"(/,a)") '--> checking 1D wind profile'
    mdot0 = xyzmh_ptmass(imloss,1)
    ! use the mass of the star the profile was computed with (before the prefill was emitted)
    mstar = mstar0*umass
    do i=1,ngrid
       r = wind_profiles(1)%trvurho(2,i)
       v = wind_profiles(1)%trvurho(3,i)
       rho1d = wind_profiles(1)%trvurho(5,i)
       u = wind_profiles(1)%trvurho(4,i)

       ! check that the mass flux is constant at every radius in the 1D wind profile
       mdot = 4.*pi*r**2*rho1d*v/unit_mdot
       call checkvalbuf(mdot,mdot0,tol_mdot,'dM/dt = 4 pi r^2 rho v',nfailed(1),ncheck(1),errmax(1))

       ! check Bernoulli energy is constant at every radius in the 1D wind profile
       ! only possible if the Eddington factor (alpha) is constant!
       if (icase == 1) then
          if (isothermal) then
             e = 0.5*v**2 - (1.-xyzmh_ptmass(iwalpha,1))*gg*mstar/r + polyk*log(rho1d)
          else
             e = 0.5*v**2 - (1.-xyzmh_ptmass(iwalpha,1))*gg*mstar/r + gamma*u
          endif
          e = 0.5*v**2 - (1.-xyzmh_ptmass(iwalpha,1))*gg*mstar/r + gamma*u
          if (i == 1) e0 = e
          call checkvalbuf(e,e0,tol_e,'Bernoulli constant',nfailed(2),ncheck(2),errmax(2))
       endif
    enddo
    call checkvalbuf_end('dM/dt = 4 pi r^2 rho v',ncheck(1),nfailed(1),errmax(1),tol_mdot)
    if (icase == 1) call checkvalbuf_end('Bernoulli constant',ncheck(2),nfailed(2),errmax(2),tol_e)
    call update_test_scores(ntests,nfailed,npass)
 endif

end subroutine init_testwind

!-----------------------------------------------------------------------
!+
!  checks the wind profile against the 1D solution
!+
!-----------------------------------------------------------------------
subroutine test_against_1D_profile(ntests,npass,npart,xyzh,vxyzu,isink,xyzmh_ptmass,rmin,rmax)
 use wind,      only:interp_wind_profile_at_r
 use part,      only:rho,isdead_or_accreted
 use units,     only:udist,unit_density,unit_ergg,unit_velocity
 use physcon,   only:au
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: npart,isink
 real,    intent(in)    :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),rmin,rmax
 integer :: i,nfailed(3),ncheck(3)
 real :: dx(3),r,rhoi,ui,vi,rho1d,u,v,errmax(3)
 real, parameter :: tol_v = 1.6e-1, tol_u = 1.2e-1, tol_rho = 9.e-16

 if (id==master) write(*,"(/,a,2(f7.2,a))") &
    '--> checking wind profile against 1D for r between ',rmin*udist/au,' and ',rmax*udist/au,' au'

 nfailed = 0
 ncheck = 0
 errmax = 0.
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       dx = xyzh(1:3,i)-xyzmh_ptmass(1:3,isink)
       r = sqrt(dot_product(dx,dx))
       if (r > rmin .and. r < rmax) then
          rhoi = rho(i)*unit_density
          vi = dot_product(vxyzu(1:3,i),dx/r)*unit_velocity
          ui = vxyzu(4,i)*unit_ergg
          r = r*udist
          call interp_wind_profile_at_r(r,v,u,rho1d,isink)
          call checkvalbuf(v,vi,tol_v,'v',nfailed(1),ncheck(1),errmax(1))
          call checkvalbuf(u,ui,tol_u,'u',nfailed(2),ncheck(2),errmax(2))
          call checkvalbuf(rho1d,rhoi,tol_rho,'rho',nfailed(3),ncheck(3),errmax(3))
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
!  checks the 1D wind solver against the exact solution for an
!  isothermal (Parker) wind with constant Eddington factor alpha,
!  written in terms of the Lambert W function (Cranmer 2004)
!+
!-----------------------------------------------------------------------
subroutine test_parker_wind(ntests,npass,rmax)
 use physcon, only:au,solarm,years,gg,Rg
 use units,   only:umass,udist,utime
 use part,    only:nptmass,xyzmh_ptmass,nsinkproperties,iTeff,iTwind,ivwind,imloss,iReff,iwalpha
 use eos,     only:gmw,gamma,ieos,init_eos,polyk
 use options, only:icooling,write_files
 use wind,    only:wind_params,get_wind_solution,save_windprofile,wind_profiles
 use ptmass_radiation, only:isink_radiation
 use dust_formation,   only:idust_opacity
 integer, intent(inout) :: ntests,npass
 real,    intent(in)    :: rmax
 integer, parameter :: nalpha = 2
 real,    parameter :: alpha_test(nalpha) = (/0.,0.2/)
 real,    parameter :: tol_rs = 1.1e-2, tol_v0 = 3.5e-4, tol_v = 5.e-3
 type(wind_params) :: params
 integer :: i,j,ierr,nfailed(3),ncheck,icooling_prev,isink_radiation_prev
 real    :: xyzmh(nsinkproperties),rsonic,cs,rc,r,v,vexact,errmax,tcross,tfill
 logical :: write_files_prev

 if (id==master) write(*,"(/,a)") '--> testing 1D wind solver against exact isothermal (Parker) wind'

 ! mass-losing sink particle with unknown wind velocity, to be found by the 1D solver
 nptmass = 1
 xyzmh_ptmass = 0.
 xyzmh_ptmass(4,1)      = 1.2*solarm/umass
 xyzmh_ptmass(5,1)      = au/udist
 xyzmh_ptmass(iTeff,1)  = 50000.
 xyzmh_ptmass(iTwind,1) = 50000.
 xyzmh_ptmass(ivwind,1) = 0.
 xyzmh_ptmass(imloss,1) = 1.d-5*(solarm/umass)/(years/utime)
 xyzmh_ptmass(iReff,1)  = 2.*au/udist

 ! isothermal gas, no cooling and no dust
 gmw   = 2.381
 gamma = 1.
 ieos  = 1
 polyk = 0.
 call init_eos(ieos,ierr)
 icooling_prev        = icooling
 isink_radiation_prev = isink_radiation
 write_files_prev     = write_files
 icooling        = 0
 idust_opacity   = 0
 isink_radiation = 1   ! constant Eddington factor, alpha_rad
 write_files     = .false.

 cs = sqrt(Rg*xyzmh_ptmass(iTwind,1)/gmw)
 do j=1,nalpha
    if (id==master) write(*,"(/,a,f4.2)") '--> Parker wind with alpha = ',alpha_test(j)
    xyzmh_ptmass(iwalpha,1) = alpha_test(j)
    nfailed = 0; ncheck = 0; errmax = 0.

    ! trans-sonic solution: sonic radius and velocity at the injection radius
    call get_wind_solution(xyzmh_ptmass(:,1),1,xyzmh,params,rsonic)
    rc = gg*params%Mstar*(1.-alpha_test(j))/(2.*cs**2)
    call checkval(rsonic,rc,tol_rs,nfailed(1),'sonic radius')
    vexact = parker_wind_velocity(params%Rinject,rc,cs)
    call checkval(params%vwind,vexact,tol_v0,nfailed(2),'injection velocity')

    ! velocity profile through and beyond the sonic point
    call save_windprofile(params,real(rmax*udist),real(rmax*udist),0.,tcross,tfill,'parker_profile.dat',1)
    do i=1,size(wind_profiles(1)%trvurho,dim=2)
       r = wind_profiles(1)%trvurho(2,i)
       v = wind_profiles(1)%trvurho(3,i)
       call checkvalbuf(v,parker_wind_velocity(r,rc,cs),tol_v,'v(r)',nfailed(3),ncheck,errmax)
    enddo
    call checkvalbuf_end('v(r) against Lambert W solution',ncheck,nfailed(3),errmax,tol_v)
    call update_test_scores(ntests,nfailed,npass)
 enddo

 ! one 1D profile per sink, for several sinks with different winds
 call test_profiles_per_sink(ntests,npass,rmax)

 icooling        = icooling_prev
 isink_radiation = isink_radiation_prev
 write_files     = write_files_prev

end subroutine test_parker_wind

!-----------------------------------------------------------------------
!+
!  checks that each sink has its own 1D wind profile: three sinks where
!  sink 2 has no wind and sinks 1 and 3 have different mass loss rates.
!  Each profile must carry the mass loss rate of its own sink, and the
!  age of the wind must be recovered from its radius
!+
!-----------------------------------------------------------------------
subroutine test_profiles_per_sink(ntests,npass,rmax)
 use physcon, only:pi,solarm,years
 use units,   only:umass,udist,utime
 use part,    only:nptmass,xyzmh_ptmass,nsinkproperties,imloss,iwalpha
 use wind,    only:wind_params,get_wind_solution,save_windprofile,wind_profiles,&
                   has_wind_profile,interp_wind_profile,wind_age_at_r
 integer, intent(inout) :: ntests,npass
 real,    intent(in)    :: rmax
 integer, parameter :: nsinks = 3
 real,    parameter :: tol_mdot = 4.e-16, tol_age = 1.e-14
 real,    parameter :: mdot_msun_yr(nsinks) = (/1.e-5,0.,3.e-6/)
 type(wind_params) :: params
 integer :: i,k,nfailed(5),ncheck(2)
 real    :: xyzmh(nsinkproperties),rsonic,tcross,tfill,mdot,errmax(2),age,r,v,u,rho

 if (id==master) write(*,"(/,a)") '--> checking 1D wind profiles for several sinks'

 nptmass = nsinks
 do k=1,nsinks
    xyzmh_ptmass(:,k) = xyzmh_ptmass(:,1)
    xyzmh_ptmass(iwalpha,k) = 0.
    xyzmh_ptmass(imloss,k)  = mdot_msun_yr(k)*(solarm/umass)/(years/utime)
 enddo
 if (allocated(wind_profiles)) deallocate(wind_profiles)
 do k=1,nsinks
    if (xyzmh_ptmass(imloss,k) <= 0.) cycle
    call get_wind_solution(xyzmh_ptmass(:,k),1,xyzmh,params,rsonic)
    call save_windprofile(params,real(rmax*udist),real(rmax*udist),0.,tcross,tfill,'profile.dat',k)
 enddo

 nfailed = 0; ncheck = 0; errmax = 0.
 call checkval(has_wind_profile(1),.true.,nfailed(1),'sink 1 has a wind profile')
 call checkval(has_wind_profile(2),.false.,nfailed(2),'sink 2 has no wind profile')
 call checkval(has_wind_profile(3),.true.,nfailed(3),'sink 3 has a wind profile')
 do k=1,nsinks,2
    ! the mass flux 4 pi r^2 rho v at every radius is the mass loss rate of that sink
    do i=1,size(wind_profiles(k)%trvurho,dim=2)
       mdot = 4.*pi*wind_profiles(k)%trvurho(2,i)**2*wind_profiles(k)%trvurho(5,i)*&
              wind_profiles(k)%trvurho(3,i)/(mdot_msun_yr(k)*solarm/years)
       call checkvalbuf(mdot,1.,tol_mdot,'dM/dt of own sink',nfailed(4),ncheck(1),errmax(1))
    enddo
    ! r -> age is the inverse of age -> r (from row 2, since row 1 is at age zero)
    do i=2,size(wind_profiles(k)%trvurho,dim=2),10
       age = 0.7*wind_profiles(k)%trvurho(1,i)/utime
       call interp_wind_profile(k,age,r,v,u,rho)
       call checkvalbuf(wind_age_at_r(k,r)/age,1.,tol_age,'age from radius',nfailed(5),ncheck(2),errmax(2))
    enddo
 enddo
 call checkvalbuf_end('dM/dt = 4 pi r^2 rho v of own sink',ncheck(1),nfailed(4),errmax(1),tol_mdot)
 call checkvalbuf_end('age from radius',ncheck(2),nfailed(5),errmax(2),tol_age)
 call update_test_scores(ntests,nfailed,npass)
 nptmass = 1

end subroutine test_profiles_per_sink

!-----------------------------------------------------------------------
!+
!  winds from several sinks: sink 2 has no wind, sinks 1 and 3 have
!  different mass loss rates (and hence different numbers of particles
!  per shell). Particles are moved ballistically between injections.
!  Checks the particles and mass emitted by each sink, that boundary
!  shells are complete and keep a common radius, and that reversing the
!  order of the particles and injecting again as on a restart changes
!  nothing. All checks are global, so the test also works with MPI
!+
!-----------------------------------------------------------------------
subroutine test_multiple_wind_sinks(ntests,npass)
 use inject,       only:init_inject,inject_particles,set_default_options_inject,iboundary_spheres
 use units,        only:umass,udist,utime,unit_velocity,unit_luminosity
 use physcon,      only:au,solarm,solarl,km,years,kboltz,mass_proton_cgs
 use eos,          only:gmw,ieos,init_eos,gamma,polyk
 use part,         only:npart,nptmass,xyzmh_ptmass,vxyz_ptmass,xyzh,vxyzu,init_part,npartoftype,&
                        igas,iboundary,iTeff,iLum,iReff,massoftype,iTwind,ivwind,imloss,iwalpha,&
                        rho,ieject,iseed_sink,iphase,iamtype,sync_aprmassoftype
 use timestep,     only:tmax,dtmax
 use timestep_ind, only:nbinmax
 use wind,         only:rfill_domain_au,get_wind_particle_mass
 use partinject,   only:update_injected_particles
 use mpiutils,     only:reduceall_mpi
 use ptmass_radiation, only:isink_radiation,get_Teff
 use dust_formation,   only:idust_opacity
 use dim,              only:isothermal
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nsinks = 3, nsteps = 40
 real,    parameter :: mdot_msun_yr(nsinks) = (/1.e-5,0.,3.e-6/), xsink(nsinks) = (/0.,-20.,20./)
 real,    parameter :: tol_sum = 2.e-14, tol_m = 1.e-15
 integer :: i,k,ierr,istep,istepfrac,npart_old,nfailed(9),nshell(nsinks),nps(nsinks),nb(nsinks)
 integer(kind=8) :: nk(nsinks)
 real    :: t,dt,dtinject,mtot0,mtot,m0(nsinks),dt_shell(nsinks),sums(2,nsinks),sums_new(2,nsinks)
 logical :: one_radius(nsinks)

 if (id==master) write(*,"(/,a)") '--> testing winds from several sinks'
 call init_part
 nptmass = nsinks
 xyzmh_ptmass = 0.
 vxyz_ptmass  = 0.
 do k=1,nsinks
    xyzmh_ptmass(1,k)      = xsink(k)*au/udist
    xyzmh_ptmass(4,k)      = 1.2*solarm/umass
    xyzmh_ptmass(5,k)      = au/udist
    xyzmh_ptmass(iReff,k)  = 2.*au/udist
    xyzmh_ptmass(iTeff,k)  = get_Teff(2.e4,1.)
    xyzmh_ptmass(iLum,k)   = 2.e4*solarl/unit_luminosity
    xyzmh_ptmass(iTwind,k) = 2500.
    xyzmh_ptmass(iwalpha,k) = 0.
    xyzmh_ptmass(imloss,k) = mdot_msun_yr(k)*(solarm/umass)/(years/utime)
    if (mdot_msun_yr(k) > 0.) xyzmh_ptmass(ivwind,k) = 40.*km/unit_velocity
    m0(k) = xyzmh_ptmass(4,k)
 enddo
 gmw = 2.381
 if (isothermal) then
    ! isothermal wind at the wind temperature, set by polyk = cs^2
    gamma = 1.
    ieos  = 1
    polyk = kboltz*xyzmh_ptmass(iTwind,1)/(mass_proton_cgs*gmw*unit_velocity**2)
 else
    gamma = 1.4
    ieos  = 2
    polyk = 0.
 endif
 call init_eos(ieos,ierr)
 isink_radiation = 0
 idust_opacity   = 0
 dtmax = 1.
 tmax  = 12.
 rfill_domain_au = 0.
 call set_default_options_inject(0)
 iboundary_spheres = 5
 massoftype(igas) = get_wind_particle_mass(xyzmh_ptmass(:,1),0,100,0)
 call init_inject(ierr)
 call sync_aprmassoftype()
 do k=1,nsinks
    nps(k) = nint(xyzmh_ptmass(ieject,k))
    if (mdot_msun_yr(k) > 0.) dt_shell(k) = massoftype(igas)*nps(k)/xyzmh_ptmass(imloss,k)
 enddo
 mtot0 = sum(xyzmh_ptmass(4,1:nsinks))

 ! emit the boundary shells, then move the particles ballistically between injections
 t = 0.
 dt = 0.
 istepfrac = 0
 npart_old = npart
 call inject_particles(t,0.,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)
 do istep=1,nsteps
    dt = 0.2*minval(dt_shell,mask=mdot_msun_yr > 0.)
    do i=1,npart
       xyzh(1:3,i) = xyzh(1:3,i) + dt*vxyzu(1:3,i)
    enddo
    t = t + dt
    npart_old = npart
    call inject_particles(t,dt,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)
 enddo

 ! particles, boundary particles and mass emitted by each sink
 nfailed = 0
 do k=1,nsinks
    nk(k) = reduceall_mpi('+',int(count(iseed_sink(1:npart) == k),kind=8))
    nb(k) = int(reduceall_mpi('+',count(iseed_sink(1:npart) == k .and. iamtype(iphase(1:npart)) == iboundary)))
    nshell(k) = 0
    if (mdot_msun_yr(k) > 0.) nshell(k) = nint(nb(k)/real(nps(k))) + floor(t/dt_shell(k))
 enddo
 call checkval(nb(1),iboundary_spheres*nps(1),0,nfailed(1),'boundary shells complete (sink 1)')
 call checkval(nb(3),iboundary_spheres*nps(3),0,nfailed(1),'boundary shells complete (sink 3)')
 call checkval(nk(2),0_8,0,nfailed(2),'no particles from sink 2')
 call checkval(xyzmh_ptmass(4,2),m0(2),0.,nfailed(3),'mass of sink 2 unchanged')
 call checkval(nk(1),int(nps(1)*nshell(1),kind=8),0,nfailed(4),'particles emitted by sink 1')
 call checkval(nk(3),int(nps(3)*nshell(3),kind=8),0,nfailed(5),'particles emitted by sink 3')
 call checkval(xyzmh_ptmass(4,1),m0(1)-nk(1)*massoftype(igas),tol_m,nfailed(6),'mass lost by sink 1')
 call checkval(xyzmh_ptmass(4,3),m0(3)-nk(3)*massoftype(igas),tol_m,nfailed(7),'mass lost by sink 3')
 mtot = sum(xyzmh_ptmass(4,1:nsinks)) + reduceall_mpi('+',npart)*massoftype(igas)
 call checkval(mtot,mtot0,tol_m,nfailed(8),'total mass conserved')
 call get_boundary_sums(nsinks,npart,xyzh,vxyzu,xyzmh_ptmass,sums,one_radius)
 call checkval(all(one_radius(1:nsinks:2)),.true.,nfailed(9),'boundary shells have one radius')
 call update_test_scores(ntests,nfailed,npass)

 ! reverse the order of the particles and inject again, as on a restart
 nfailed = 0
 xyzh(:,1:npart)    = xyzh(:,npart:1:-1)
 vxyzu(:,1:npart)   = vxyzu(:,npart:1:-1)
 iphase(1:npart)    = iphase(npart:1:-1)
 iseed_sink(1:npart) = iseed_sink(npart:1:-1)
 npart_old = npart
 mtot0 = sum(xyzmh_ptmass(4,1:nsinks))
 call inject_particles(t,0.,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)
 call get_boundary_sums(nsinks,npart,xyzh,vxyzu,xyzmh_ptmass,sums_new,one_radius)
 call checkval(npart,npart_old,0,nfailed(1),'no particles emitted on restart')
 call checkval(sum(xyzmh_ptmass(4,1:nsinks)),mtot0,0.,nfailed(2),'sink masses unchanged on restart')
 call checkval(sums_new(1,1),sums(1,1),tol_sum,nfailed(3),'sum of boundary radii, sink 1')
 call checkval(sums_new(2,1),sums(2,1),tol_sum,nfailed(4),'sum of boundary speeds, sink 1')
 call checkval(sums_new(1,3),sums(1,3),tol_sum,nfailed(5),'sum of boundary radii, sink 3')
 call checkval(sums_new(2,3),sums(2,3),tol_sum,nfailed(6),'sum of boundary speeds, sink 3')
 call update_test_scores(ntests,nfailed(1:6),npass)

end subroutine test_multiple_wind_sinks

!-----------------------------------------------------------------------
!+
!  sums over all MPI tasks of the distance to the sink and the speed
!  relative to the sink of the boundary particles of each sink, and
!  whether the boundary particles of each shell share one radius
!  (i.e. there are as many distinct radii as boundary shells)
!+
!-----------------------------------------------------------------------
subroutine get_boundary_sums(nsinks,npart,xyzh,vxyzu,xyzmh_ptmass,sums,one_radius)
 use part,     only:iboundary,iseed_sink,iphase,iamtype,ieject
 use mpiutils, only:reduceall_mpi
 integer, intent(in)  :: nsinks,npart
 real,    intent(in)  :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:)
 real,    intent(out) :: sums(2,nsinks)
 logical, intent(out) :: one_radius(nsinks)
 integer, parameter :: maxradii = 100
 real,    parameter :: tol_r = 1.e-10
 integer :: i,k,nradii,nb
 real    :: r,radii(maxradii)

 sums = 0.
 do k=1,nsinks
    nradii = 0
    nb = 0
    do i=1,npart
       if (iseed_sink(i) /= k .or. iamtype(iphase(i)) /= iboundary) cycle
       nb = nb + 1
       r = sqrt(dot_product(xyzh(1:3,i)-xyzmh_ptmass(1:3,k),xyzh(1:3,i)-xyzmh_ptmass(1:3,k)))
       sums(1,k) = sums(1,k) + r
       sums(2,k) = sums(2,k) + sqrt(dot_product(vxyzu(1:3,i),vxyzu(1:3,i)))
       if (all(abs(radii(1:nradii)-r) > tol_r*r) .and. nradii < maxradii) then
          nradii = nradii + 1
          radii(nradii) = r
       endif
    enddo
    ! each MPI task holds whole shells, so the shells on this task have distinct radii
    one_radius(k) = (nradii*nint(xyzmh_ptmass(ieject,k)) == nb)
 enddo
 sums = reduceall_mpi('+',sums)
 do k=1,nsinks
    one_radius(k) = (reduceall_mpi('+',merge(0,1,one_radius(k))) == 0)
 enddo

end subroutine get_boundary_sums

!-----------------------------------------------------------------------
!+
!  exact velocity of an isothermal trans-sonic wind with sound speed cs
!  and sonic radius rc = GM(1-alpha)/(2 cs^2) (Cranmer 2004):
!
!    v^2/cs^2 = -W[-D(r)],  D(r) = (rc/r)^4 exp(3 - 4 rc/r)
!
!  where W is the Lambert W function, using the principal branch for
!  the subsonic flow (r < rc) and the -1 branch for supersonic flow
!+
!-----------------------------------------------------------------------
real function parker_wind_velocity(r,rc,cs) result(v)
 use mathfunc, only:lambertw
 real, intent(in) :: r,rc,cs
 real :: d

 d = (rc/r)**4*exp(3. - 4.*rc/r)
 if (r < rc) then
    v = cs*sqrt(-lambertw(-d,0))
 else
    v = cs*sqrt(-lambertw(-d,-1))
 endif

end function parker_wind_velocity

!-----------------------------------------------------------------------
!+
!  checks the injected mass
!+
!-----------------------------------------------------------------------
subroutine test_injected_mass(ntests,npass,npart,npart_prefill,isink,xyzmh_ptmass,mass0,tmax)
 use part,      only:massoftype,igas,ieject,imloss,imacc
 use mpiutils,  only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: npart,npart_prefill,isink
 real,    intent(in)    :: xyzmh_ptmass(:,:),mass0,tmax
 integer :: neject,npart_per_sphere,nall,nprefill
 integer :: nfailed(5)
 real :: mstar,minject,mprefill,tol_mass

 ! numbers of particles summed over all MPI tasks
 nall     = int(reduceall_mpi('+',npart))
 nprefill = int(reduceall_mpi('+',npart_prefill))

 ! injected mass is the total mass of particles minus the boundary shells
 mstar = mass0
 minject = xyzmh_ptmass(imloss,isink)*tmax
 mprefill = nprefill*massoftype(igas)
 neject  = nint(minject/massoftype(igas))
 nfailed = 0
 npart_per_sphere = nint(1.2*xyzmh_ptmass(ieject,isink))  ! 1.2 is a safety factor
 tol_mass = npart_per_sphere*massoftype(igas)/minject
 call checkval(xyzmh_ptmass(4,isink),mstar-minject-mprefill,8e-6,nfailed(1),'sink particle mass')
 call checkval(xyzmh_ptmass(imacc,isink),0.,epsilon(0.),nfailed(2),'mass accreted')
 call checkval(minject,(nall-nprefill)*massoftype(igas),tol_mass,nfailed(3),'mass injected')
 call checkval(nall-nprefill,neject,npart_per_sphere,nfailed(4),'number of ejected particles')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_injected_mass

!-----------------------------------------------------------------------
!+
!  integrates the wind test
!+
!-----------------------------------------------------------------------
subroutine integrate_wind(npart_old,istepfrac,dtinject)
 use io,             only:iprint
 use timestep,       only:time,tmax,dt,dtmax,nsteps,dtrad,dtforce,dtcourant,dterr,print_dtlog
 use part,           only:npart,xyzmh_ptmass,vxyz_ptmass,xyzh,vxyzu,npartoftype,ntot,rho
 use timestep_ind,   only:nbinmax
 use step_lf_global, only:step,init_step
 use partinject,     only:update_injected_particles
 use inject,         only:inject_particles

 integer, intent(inout) :: istepfrac,npart_old
 real,    intent(inout) :: dtinject

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
    call inject_particles(t,dtlast,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
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
