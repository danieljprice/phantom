!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: teststep
!
!  DESCRIPTION:
!  Unit tests of the step module
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, eos, io, mpiutils, options, part, physcon,
!    step_lf_global, testutils, timestep, timestep_ind, timing, unifdis,
!    viscosity
!+
!--------------------------------------------------------------------------
module teststep
 implicit none
 public :: test_step

 private

contains

subroutine test_step(ntests,npass)
 use io,       only:id,master
#ifdef PERIODIC
 use io,       only:iverbose
 use dim,      only:maxp,maxvxyzu,maxalpha
 use boundary, only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use eos,      only:polyk,gamma,use_entropy
 use mpiutils, only:reduceall_mpi
 use options,  only:tolh,alpha,alphau,alphaB,ieos
 use part,     only:npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu,divcurlv, &
                    Bevol,dBevol,Bextx,Bexty,Bextz,alphaind,fext, &
                    maxphase,mhd,maxBevol,igas
 use unifdis,  only:set_unifdis
 use physcon,  only:pi
 use timing,   only:getused
 use step_lf_global,  only:step,init_step
 use timestep,        only:dtmax
 use viscosity,       only:irealvisc
 use part,            only:iphase,isetphase,igas
 use timestep,        only:dtmax
 use testutils,       only:checkval,checkvalf
#endif
#ifdef IND_TIMESTEPS
 use timestep_ind, only: nbinmax
 use part,         only: ibin
#endif
 integer, intent(inout) :: ntests,npass
#ifdef PERIODIC
 real                   :: psep,hzero,totmass,dt,t,dtext,dtnew_dum
 real                   :: rhozero
 integer                :: i,nsteps
 integer :: nfailed(9)

 if (id==master) write(*,"(/,a,/)") '--> TESTING STEP MODULE / boundary crossing'

 npart = 0
 psep = dxbound/50.
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,npart,xyzh)

 npartoftype(:) = 0
 npartoftype(1) = npart
 !print*,' thread ',id,' npart = ',npart
 iverbose = 0

 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)

 rhozero = 7.5
 totmass = rhozero/(dxbound*dybound*dzbound)
 massoftype(igas) = totmass/real(reduceall_mpi('+',npart))
 hzero = hfact*(massoftype(igas)/rhozero)**(1./3.)
!
!--set constant velocity (in all components)
!
 vxyzu(1:3,:) = 1.
!
!--set everything else to zero
!
 if (maxvxyzu >= 4) vxyzu(4,:) = 0.
 fxyzu(:,:) = 0.
 fext(:,:)  = 0.
 Bevol(:,:) = 0.
 Bextx = 0.
 Bexty = 0.
 Bextz = 0.
 dBevol(:,:) = 0.
 divcurlv(:,:) = 0.
 polyk = 0.
!
!--make sure AV is off
!
 alpha = 0.
 alphau = 0.
 alphaB = 0.
 tolh = 1.e-5
 if (maxalpha==maxp) alphaind = 0.
 irealvisc = 0
!
!--use isothermal or adiabatic equation of state
!
 if (maxvxyzu==4) then
    ieos = 2
    gamma = 5./3.
 else
    ieos = 1
    gamma = 1.0
 endif

 nsteps  = 10
 dt      = 2.0/(nsteps)
 dtmax   = dt
 t = 0.

 ! If using individual timesteps, ibin may be uninitialised
#ifdef IND_TIMESTEPS
 do i = 1, npart
    ibin(i) = nbinmax
 enddo
#endif

 call init_step(npart,t,dtmax)

 nfailed(:) = 0
 do i=1,nsteps
    t = t + dt
    dtext = dt
    if (id==master) write(*,*) ' t = ',t,' dt = ',dt
    call step(npart,npart,t,dt,dtext,dtnew_dum)
    nfailed(:) = 0
    call checkval(npart,xyzh(4,:),hzero,3.e-4,nfailed(1),'h (density)')
    call checkval(npart,fxyzu(1,:),0.,tiny(fxyzu),nfailed(2),'fx')
    call checkval(npart,fxyzu(2,:),0.,tiny(fxyzu),nfailed(3),'fy')
    call checkval(npart,fxyzu(3,:),0.,tiny(fxyzu),nfailed(4),'fz')
    if (maxvxyzu >= 4) call checkval(npart,fxyzu(4,:),0.,tiny(fxyzu),nfailed(5),'du/dt')
    if (mhd) then
       call checkval(npart,dBevol(1,:),0.,tiny(0.),nfailed(6),'dBevolx/dt')
       call checkval(npart,dBevol(2,:),0.,tiny(0.),nfailed(7),'dBevoly/dt')
       call checkval(npart,dBevol(3,:),0.,tiny(0.),nfailed(8),'dBevolz/dt')
       if (maxBevol==4) call checkval(npart,dBevol(4,:),0.,tiny(0.),nfailed(9),'dpsi/dt')
    endif
    ntests = ntests + 1
    if (all(nfailed(:)==0)) npass = npass + 1
 enddo

 if (id==master) write(*,"(/,a)") '<-- STEP TEST COMPLETE'

#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF STEP MODULE (need -DPERIODIC)'

#endif

end subroutine test_step

end module teststep
