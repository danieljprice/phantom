!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module testwind
!
! Unit tests of the step module
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checksetup, deriv, dim, eos, io, mpidomain,
!   mpiutils, options, part, physcon, step_lf_global, testutils, timestep,
!   timestep_ind, timing, unifdis, viscosity
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
 use io,        only:iprint,id,master,iverbose
 use boundary,   only:set_boundary
 use options,   only:icooling,ieos,tolh,alpha,alphau,alphaB,avdecayconst,beta
 use physcon,   only:au,solarm,mass_proton_cgs,kboltz,solarl
 use units,     only:umass,set_units,utime,unit_energ,udist!,unit_velocity
 use inject,    only:init_inject,inject_particles
 use eos,       only:gmw,ieos,init_eos,gamma!,polyk
 use part,      only:npart,init_part,nptmass,xyzmh_ptmass,vxyz_ptmass,xyzh,vxyzu,nptmass,&
      npartoftype,igas,iTeff,iLum,iReff,massoftype,ntot,hfact
 use physcon,   only:au,solarm,mass_proton_cgs,kboltz, solarl
 use timestep,  only:time,tmax,dt,dtmax,nsteps,dtrad,dtforce,dtcourant,dterr,print_dtlog,&
      tolv,C_cour,C_force
 use testutils,      only:checkval,update_test_scores
 use dim,            only:isothermal
 use step_lf_global, only:step,init_step
 use partinject,     only:update_injected_particles
 use timestep_ind,   only:nbinmax
 use wind,           only:trvurho_1D
 use damping,        only:idamp
 use checksetup,     only:check_setup

 integer, intent(inout) :: ntests,npass

 integer :: i,ierr,nerror,istepfrac,npart_old,nfailed(9),nwarn
 real :: dtinject,dtlast,t,default_particle_mass,dtext,dtnew,dtprint,dtmaxold,tprint

 if (id==master) write(*,"(/,a,/)") '--> TESTING WIND MODULE'

 call set_units(dist=au,mass=solarm,G=1.)
 call set_boundary(-50.,50.,-50.,50.,-50.,50.)

 call init_part()

 !set properties of mass loosing sink particle
 nptmass = 1
 xyzmh_ptmass(4,1)     = 1.2*solarm/umass
 xyzmh_ptmass(5,1)     = 0.6*au/udist
 xyzmh_ptmass(iTeff,1) = 2500.
 xyzmh_ptmass(iReff,1) = au/udist
 xyzmh_ptmass(iLum,1)  = 9150.47 *solarl * utime / unit_energ

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
    gamma = 1.66666667
    ieos = 2
 endif
 !polyk = kboltz*T_wind/(mass_proton_cgs * gmw * unit_velocity**2)

 call init_eos(ieos,ierr)

 iverbose = 1
 hfact    = 1.2
 tolv     = 1.e-2
 tolh     = 1.e-4
 C_cour   = 0.3
 C_force  = 0.25
 nsteps   = 0
 alpha    = 0.
 alphau   = 1.
 alphaB   = 0.
 beta     = 2.
 idamp    = 0
 avdecayconst = 0.1
 icooling = 0
 dtmax    = 1.
 tmax     = 4.
 tprint   = tmax
 dt       = 0.
 dtinject = huge(dtinject)
 dtrad    = huge(dtrad)
 t        = 0.
 dtnew    = 0.


 call check_setup(nerror,nwarn)

 istepfrac  = 0
 nfailed(:) = 0
 call init_inject(nerror)
! check particle's mass
 call checkval(massoftype(igas),6.820748526700016e-10,epsilon(0.),&
      nfailed(1),'no errors in setting particle mass')
 npart_old = npart
 call inject_particles(t,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                       npart,npartoftype,dtinject)
 call update_injected_particles(npart_old,npart,istepfrac,nbinmax,t,dtmax,dt,dtinject)

! check 1D wind profile
 i = size(trvurho_1D(1,:))
 call checkval(trvurho_1D(2,i),4.847727063707113E+13,epsilon(0.),nfailed(2),'outer wind radius')
 call checkval(trvurho_1D(3,i),4.201251511248741E+05,epsilon(0.),nfailed(3),'outer wind velocity')
 call checkval(trvurho_1D(4,i),1.772793525117932E+11,epsilon(0.),nfailed(4),'outer wind temperature')
 call checkval(trvurho_1D(5,i),5.080390084918874E-14,epsilon(0.),nfailed(5),'outer wind density')

 dt = dtinject
 dtlast = 0.
 time = 0.

 do while (t < tmax)

    dtext = dt
    if (id==master) write(*,*) ' t = ',t,' dt = ',dt
    !
    ! injection of new particles into simulation
    !
    npart_old=npart
    call inject_particles(t,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,dtinject)
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

 call checkval(xyzmh_ptmass(4,1),1.199993907707417d0,epsilon(0.),nfailed(6),'sink particle mass')
 call checkval(xyzmh_ptmass(7,1),0.,epsilon(0.),nfailed(7),'mass accreted')
 call checkval(npart,12992,0,nfailed(8),'number of ejected particles')
 call checkval(xyzmh_ptmass(15,1),1.591640703559762d-6,epsilon(0.),nfailed(9),'wind mass loss rate')

 call update_test_scores(ntests,nfailed,npass)

 if (id==master) write(*,"(/,a)") '<-- STEP TEST COMPLETE'

end subroutine test_wind

end module testwind
