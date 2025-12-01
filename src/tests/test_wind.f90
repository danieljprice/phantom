!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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

 logical :: vb = .true.
 real, parameter :: eps_sum = 4.e-14

contains
!----------------------------------------------------------
!+
!  Unit tests of timestepping and boundary crossing
!+
!----------------------------------------------------------
subroutine test_wind(ntests,npass)
 use io,         only:id,master,iprint,iwritein
 use inject,     only:inject_type
 use boundary,   only:set_boundary
 use physcon,    only:au,solarm,solarl
 use units,      only:set_units
 use part,       only:npart,xyzmh_ptmass,vxyzu,dust_temp,igas,massoftype,imloss
 use testutils,  only:checkval,update_test_scores
 use dim,        only:mpi,maxTdust,maxp,sink_radiation,nucleation,ind_timesteps,disc_viscosity,nalpha
 use allocutils, only:allocate_array
 use options,    only:alpha
 use timestep,   only:tmax
 use readwrite_infile, only:read_infile,write_infile

 integer, intent(inout) :: ntests,npass

 integer, parameter :: npart_per_shell = 812, nboundary = 5
 integer :: npart_old,nfailed(5),istepfrac,neject
 real :: dtinject,eint,ekin,mstar,minject
 logical :: testkd,testcyl,test2,use_shock_switch

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
 testkd  = sink_radiation .and. nucleation .and. use_shock_switch .and. ind_timesteps
 test2   = .not.sink_radiation .and. .not.nucleation .and. .not.use_shock_switch .and. .not.ind_timesteps
 testcyl = .not.sink_radiation .and. .not.nucleation .and. use_shock_switch .and. ind_timesteps
!transonic, testcyl=F, testkd=F, test2=T  1.199987894792037E+00  0.000000000000000E+00  1.591640703559762E-06  3.366824949389495E+03  5.525582106704587E+01   12180
 disc_viscosity = .false.
 if (testcyl) then
    alpha = 1.
    disc_viscosity = .true.
 endif

! test trans-sonic wind - no radiation, no dust

 call init_testwind(1,ntests,npass,npart_old,istepfrac,dtinject)
 mstar = xyzmh_ptmass(4,1)
 if (id==master) call write_infile('w.in','w.log','w.ev','w_00000',iwritein,iprint)
 call integrate_wind(npart_old,istepfrac,dtinject)
! remove particles from the boundary shells
 minject = (npart-nboundary*npart_per_shell)*massoftype(igas)
 nfailed(:) = 0
 eint = sum(vxyzu(4,1:npart))
 ekin = sqrt(sum(vxyzu(1,1:npart)**2+vxyzu(2,1:npart)**2+vxyzu(3,1:npart)**2))
 if (vb) print '("transonic, testcyl=",l1,", testkd=",l1,", test2=",l1,5(1x,es22.15),i8)',&
      testcyl,testkd,test2,xyzmh_ptmass(4,1),xyzmh_ptmass(7,1),xyzmh_ptmass(15,1),eint,ekin,npart
 call checkval(xyzmh_ptmass(4,1),mstar-minject,1e-3*massoftype(igas),nfailed(1),'sink particle mass')
 call checkval(xyzmh_ptmass(7,1),0.,epsilon(0.),nfailed(2),'mass accreted')
 neject = nint(xyzmh_ptmass(imloss,1)*tmax/massoftype(igas))
 call checkval(npart-nboundary*npart_per_shell,neject,npart_per_shell,nfailed(3),'number of ejected particles')
 if (testcyl) then  ! alpha is constant and equal to 1, disc_viscosity=T, no nucleation or sink radiation
    call checkval(eint,3.067302718051912E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,5.401484078064862E+01,eps_sum,nfailed(5),'total kinetic energy')
 elseif (testkd) then ! sink radiation, nucleation, ind_timesteps=T, disc_viscosity=F
    call checkval(eint,2.887208554583773E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,5.879138880015775E+01,eps_sum,nfailed(5),'total kinetic energy')
 elseif (test2) then ! no sink radiation, no nucleation, alpha=1.0, ind_timesteps=F, disc_viscosity=F
    call checkval(eint,3.366824949389491E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,5.525582106704594E+01,eps_sum,nfailed(5),'total kinetic energy')
 else
    call checkval(eint,2.909690458803881E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,5.743680106088664E+01,eps_sum,nfailed(5),'total kinetic energy')
 endif
 call update_test_scores(ntests,nfailed,npass)

 if (sink_radiation) then

! test wind with bowen dust + radiative acceleration

    maxTdust = maxp
    if (allocated(dust_temp)) deallocate(dust_temp)
    call allocate_array('dust_temp',dust_temp,maxTdust)

    call init_testwind(2,ntests,npass,npart_old,istepfrac,dtinject)
    !if (id==master) call write_infile('w2.in','w2.log','w2.ev','w2_00000',iwritein,iprint)
    call integrate_wind(npart_old,istepfrac,dtinject)
    minject = (npart-nboundary*npart_per_shell)*massoftype(igas)
    nfailed(:) = 0
    eint = sum(vxyzu(4,1:npart))
    ekin = sqrt(sum(vxyzu(1,1:npart)**2+vxyzu(2,1:npart)**2+vxyzu(3,1:npart)**2))
    if (vb) print '("sink_rad, testkd=",l1,5(1x,es22.15),i8)',testkd,&
         xyzmh_ptmass(4,1),xyzmh_ptmass(7,1),xyzmh_ptmass(15,1),eint,ekin,npart
    call checkval(xyzmh_ptmass(4,1),mstar-minject,1e-3*massoftype(igas),nfailed(1),'sink particle mass')
    call checkval(xyzmh_ptmass(7,1),0.,epsilon(0.),nfailed(2),'mass accreted')
    neject = nint(xyzmh_ptmass(imloss,1)*tmax/massoftype(igas))
    call checkval(npart-nboundary*npart_per_shell,neject,npart_per_shell,nfailed(3),'number of ejected particles')
    if (testkd) then
       call checkval(eint,2.160437358204532E+02,eps_sum,nfailed(4),'total internal energy')
       call checkval(ekin,1.658273173102032E+02,eps_sum,nfailed(5),'total kinetic energy')
    else
       call checkval(eint,2.184418217766686E+02,eps_sum,nfailed(4),'total internal energy')
       call checkval(ekin,1.658858256327306E+02,eps_sum,nfailed(5),'total kinetic energy')
    endif
 else
    if (id==master) write(*,"(/,a,/)") '    SKIPPING SINK RADIATION TEST'
 endif
 call update_test_scores(ntests,nfailed,npass)

 if (id==master) write(*,"(/,a)") '<-- WIND TEST COMPLETE'

end subroutine test_wind

!-----------------------------------------------------------------------
!
subroutine init_testwind(icase,ntests,npass,npart_old,istepfrac,dtinject)
!
!-----------------------------------------------------------------------

 use io,         only:iverbose
 use inject,     only:init_inject,inject_particles,set_default_options_inject
 use units,      only:umass,udist,unit_mdot,unit_velocity,unit_luminosity,utime
 use physcon,    only:au,solarm,solarl,km,seconds,years
 use eos,        only:gmw,ieos,init_eos,gamma,polyk
 use part,       only:npart,init_part,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,xyzh,vxyzu,&
                     npartoftype,igas,iTeff,iLum,iReff,massoftype,iTwind,ivwind,imloss
 use timestep,   only:tmax,dt,dtmax,dtrad
 use dim,        only:isothermal
 use wind,       only:trvurho_1D
 use ptmass,     only:set_integration_precision
 use testutils,  only:checkval,update_test_scores
 use checksetup, only:check_setup
 use partinject, only:update_injected_particles
 use timestep_ind,     only:nbinmax
 use ptmass_radiation, only:alpha_rad,isink_radiation
 use dust_formation,   only:idust_opacity

 integer, intent(in) :: icase
 integer, intent(inout) :: ntests,npass
 integer, intent(out) :: npart_old,istepfrac
 real,    intent(out) :: dtinject

 integer :: i,ierr,nerror,nwarn,nfailed(5)
 real :: t,default_particle_mass,dtnew

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
 dtmax = 1.
 tmax  = 8.
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
    stop '[test wind] unknown test '
 endif
 dt       = 0.
 dtinject = huge(dtinject)
 dtrad    = huge(dtrad)
 t        = 0.
 dtnew    = 0.

 ! trans-sonic wind
 call set_default_options_inject(icase)
 call check_setup(nerror,nwarn)

 nfailed(:) = 0
 istepfrac  = 0
 call init_inject(nerror)
 npart_old = npart

!trans-sonic wind - no radiation
 if (icase == 1) then
    ! check particle's mass
    call inject_particles(t,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
        npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)

    ! check 1D wind profile
    i = size(trvurho_1D(1,:))
    if (vb) print '("trans-sonic1",(6(1x,es22.15)))',massoftype(igas),trvurho_1D(:,i)
    call checkval(massoftype(igas),1.587420476277492E-09,eps_sum,nfailed(1),'setting particle mass')
    call checkval(trvurho_1D(2,i),7.099825176736505E+13, eps_sum,nfailed(2),'1D wind terminal radius')
    call checkval(trvurho_1D(3,i),1.113551988490835E+06, eps_sum,nfailed(3),'1D wind terminal velocity')
    call checkval(trvurho_1D(4,i),2.021389819449251E+12, eps_sum,nfailed(4),'1D wind internal energy')
    call checkval(trvurho_1D(5,i),8.936063664906353E-15, eps_sum,nfailed(5),'1D wind terminal density')
    call update_test_scores(ntests,nfailed,npass)
 endif

 !wind + radiation
 if (icase == 2) then
    ! check particle's mass
    call inject_particles(t,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)

    ! check 1D wind profile
    i = size(trvurho_1D(1,:))
    if (vb) print '("wind+rad",(6(1x,es22.15)))',massoftype(igas),trvurho_1D(:,i)
    call checkval(massoftype(igas),7.262861965649780E-10,eps_sum,nfailed(1),'setting particle mass')
    call checkval(trvurho_1D(2,i), 1.123761571188968E+14,eps_sum,nfailed(2),'1D wind terminal radius')
    call checkval(trvurho_1D(3,i), 2.098365055449723E+06,eps_sum,nfailed(3),'1D wind terminal velocity')
    call checkval(trvurho_1D(4,i), 7.427528334581337E+10,eps_sum,nfailed(4),'1D wind internal energy')
    call checkval(trvurho_1D(5,i), 1.892878173438733E-15,eps_sum,nfailed(5),'1D wind terminal density')
    call update_test_scores(ntests,nfailed,npass)
 endif

end subroutine init_testwind

!-----------------------------------------------------------------------
!
subroutine integrate_wind(npart_old,istepfrac,dtinject)
!
!-----------------------------------------------------------------------

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
