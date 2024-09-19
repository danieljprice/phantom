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
! :Dependencies: allocutils, boundary, checksetup, dim, dust_formation,
!   eos, inject, io, options, part, partinject, physcon, ptmass,
!   ptmass_radiation, readwrite_infile, step_lf_global, testutils,
!   timestep, timestep_ind, units, wind
!
 implicit none
 public :: test_wind

 private

 logical :: vb = .false.

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
 use part,       only:npart,xyzmh_ptmass,vxyzu,dust_temp
 use testutils,  only:checkval,update_test_scores
 use dim,        only:mpi,maxTdust,maxp,sink_radiation,nucleation,ind_timesteps
 use allocutils, only:allocate_array
 use options,    only:alphamax
 use readwrite_infile, only:read_infile,write_infile

 integer, intent(inout) :: ntests,npass

 real, parameter :: eps_sum = 1e-14
 integer :: npart_old,nfailed(5),istepfrac
 real :: dtinject,eint,ekin
 logical :: testkd,testcyl,test2

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

 testkd  = sink_radiation .and. nucleation .and. alphamax == 1. .and. ind_timesteps
 test2   = .not.sink_radiation .and. .not.nucleation .and. alphamax == 1. .and. .not.ind_timesteps
 testcyl = .not.sink_radiation .and. .not.nucleation .and. alphamax == 1. .and. ind_timesteps

! test trans-sonic wind - no radiation, no dust

 call init_testwind(1,ntests,npass,npart_old,istepfrac,dtinject)
 if (id==master) call write_infile('w.in','w.log','w.ev','w_00000',iwritein,iprint)
 call integrate_wind(npart_old,istepfrac,dtinject)
 nfailed(:) = 0
 eint = sum(vxyzu(4,1:npart))
 ekin = sqrt(sum(vxyzu(1,1:npart)**2+vxyzu(2,1:npart)**2+vxyzu(3,1:npart)**2))
 if (vb) print '(5(1x,es22.15),i8)',xyzmh_ptmass(4,1),xyzmh_ptmass(7,1),xyzmh_ptmass(15,1),eint,ekin,npart
 call checkval(xyzmh_ptmass(4,1),1.199987894518367E+00,epsilon(0.),nfailed(1),'sink particle mass')
 call checkval(xyzmh_ptmass(7,1),0.,epsilon(0.),nfailed(2),'mass accreted')
 call checkval(npart,12180,0,nfailed(3),'number of ejected particles')
 if (testcyl) then
    call checkval(eint,3.360686893182378E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,5.605632523862468E+01,eps_sum,nfailed(5),'total kinetic energy')
 elseif (testkd) then
    call checkval(eint,3.164153170427767E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,6.101010545772693E+01,eps_sum,nfailed(5),'total kinetic energy')
 elseif (test2) then
    call checkval(eint,3.367417540822784E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,5.524867074648306E+01,eps_sum,nfailed(5),'total kinetic energy')
 else
    call checkval(eint,3.179016341424608E+03,eps_sum,nfailed(4),'total internal energy')
    call checkval(ekin,6.005124961952793E+01,eps_sum,nfailed(5),'total kinetic energy')
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
    nfailed(:) = 0
    eint = sum(vxyzu(4,1:npart))
    ekin = sqrt(sum(vxyzu(1,1:npart)**2+vxyzu(2,1:npart)**2+vxyzu(3,1:npart)**2))
    if (vb) print '(5(1x,es22.15),i8)',xyzmh_ptmass(4,1),xyzmh_ptmass(7,1),xyzmh_ptmass(15,1),eint,ekin,npart
    call checkval(xyzmh_ptmass(4,1),1.199987815414834E+00,epsilon(0.),nfailed(1),'sink particle mass')
    call checkval(xyzmh_ptmass(7,1),0.,epsilon(0.),nfailed(2),'mass accreted')
    call checkval(npart,21924,0,nfailed(3),'number of ejected particles')
    if (testkd) then
       call checkval(eint,2.187465510809545E+02,eps_sum,nfailed(4),'total internal energy')
       call checkval(ekin,1.709063901093157E+02,eps_sum,nfailed(5),'total kinetic energy')
    else
       call checkval(eint,2.218461223513102E+02,eps_sum,nfailed(4),'total internal energy')
       call checkval(ekin,1.709669096834302E+02,eps_sum,nfailed(5),'total kinetic energy')
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

 use io,        only:iverbose
 use inject,    only:init_inject,inject_particles,set_default_options_inject
 use units,     only:umass,utime,unit_energ,udist
 use physcon,   only:au,solarm,solarl
 use eos,       only:gmw,ieos,init_eos,gamma,polyk
 use part,      only:npart,init_part,nptmass,xyzmh_ptmass,vxyz_ptmass,xyzh,vxyzu,&
                     npartoftype,igas,iTeff,iLum,iReff,massoftype
 use timestep,  only:tmax,dt,dtmax,dtrad
 use wind,           only:trvurho_1D
 use timestep_ind,   only:nbinmax
 use dim,            only:isothermal
 use checksetup,     only:check_setup
 use partinject,     only:update_injected_particles
 use testutils,      only:checkval,update_test_scores
 use ptmass,         only:set_integration_precision
 use ptmass_radiation, only:alpha_rad,isink_radiation
 use dust_formation,   only:idust_opacity

 integer, intent(in) :: icase
 integer, intent(inout) :: ntests,npass
 integer, intent(out) :: npart_old,istepfrac
 real, intent(out) :: dtinject

 integer :: i,ierr,nerror,nwarn,nfailed(5)
 real :: t,default_particle_mass,dtnew

 call init_part()
 call set_integration_precision()

 ! set properties of mass-losing sink particle
 nptmass = 1
 xyzmh_ptmass(4,1)  = 1.2*solarm/umass
 xyzmh_ptmass(5,1)  = au/udist
 if (icase == 1) then
    xyzmh_ptmass(iTeff,1) = 50000.
 elseif (icase == 2) then
    xyzmh_ptmass(iTeff,1) = 3000.
 endif
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
    if (vb) print '((6(1x,es22.15)))',trvurho_1D(:,i),massoftype(igas)
    call checkval(massoftype(igas),1.490822861042279E-9,epsilon(0.),nfailed(1),'setting particle mass')
    call checkval(trvurho_1D(2,i),7.058624412798283E+13,epsilon(0.),nfailed(2),'1D wind terminal radius')
    call checkval(trvurho_1D(3,i),1.112160584479353E+06,epsilon(0.),nfailed(3),'1D wind terminal velocity')
    call checkval(trvurho_1D(4,i),2.031820842001706E+12,epsilon(0.),nfailed(4),'1D wind internal energy')
    call checkval(trvurho_1D(5,i),8.878887149408118E-15,epsilon(0.),nfailed(5),'1D wind terminal density')
    call update_test_scores(ntests,nfailed,npass)
 endif

 !wind + radiation
 if (icase == 2) then
    ! check particle's mass
    call inject_particles(t,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
    call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)

    ! check 1D wind profile
    i = size(trvurho_1D(1,:))
    if (vb) print '((6(1x,es22.15)))',trvurho_1D(:,i),massoftype(igas)
    call checkval(massoftype(igas),6.820748526700016E-10,epsilon(0.),nfailed(1),'setting particle mass')
    call checkval(trvurho_1D(2,i), 1.546371444697654E+14,epsilon(0.),nfailed(2),'1D wind terminal radius')
    call checkval(trvurho_1D(3,i), 4.298693548460183E+06,epsilon(0.),nfailed(3),'1D wind terminal velocity')
    call checkval(trvurho_1D(4,i), 4.318674031561777E+10,epsilon(0.),nfailed(4),'1D wind internal energy')
    call checkval(trvurho_1D(5,i), 4.879641694552266E-16,epsilon(0.),nfailed(5),'1D wind terminal density')
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

 dt = dtinject
 dtlast = 0.
 time = 0.
 tprint   = tmax
 t = 0.

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
