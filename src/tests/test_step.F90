!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module teststep
!
! Unit tests of the step module
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checksetup, deriv, dim, eos, io, mpidomain,
!   mpiutils, options, part, physcon, step_lf_global, testutils, timestep,
!   timestep_ind, timing, unifdis, viscosity
!
 implicit none
 public :: test_step

 private

contains
!----------------------------------------------------------
!+
!  Unit tests of timestepping and boundary crossing
!+
!----------------------------------------------------------
subroutine test_step(ntests,npass)
 use io,       only:id,master
#ifdef PERIODIC
 use io,       only:iverbose
 use dim,      only:maxp,maxvxyzu,maxalpha,periodic
 use boundary, only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use eos,      only:polyk,gamma,init_eos
 use mpiutils, only:reduceall_mpi
 use options,  only:tolh,alpha,alphau,alphaB,ieos
 use part,     only:init_part,npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu, &
                    fext,dBevol,alphaind,maxphase,mhd,igas
 use unifdis,  only:set_unifdis
 use physcon,  only:pi
 use timing,   only:getused
 use step_lf_global,  only:step,init_step
 use timestep,        only:dtmax
 use viscosity,       only:irealvisc
 use part,            only:iphase,isetphase,igas
 use timestep,        only:dtmax
 use testutils,       only:checkval,checkvalf,update_test_scores
 use mpidomain,       only:i_belong
 use checksetup,      only:check_setup
 use deriv,           only:get_derivs_global
#ifdef IND_TIMESTEPS
 use part,            only:ibin
 use timestep_ind,    only:nbinmax
#endif
#endif
 integer, intent(inout) :: ntests,npass
#ifdef PERIODIC
 real                   :: psep,hzero,totmass,dt,t,dtext,dtnew_dum
 real                   :: rhozero
 integer                :: i,nsteps,nerror,nwarn,ierr
 integer :: nfailed(9)

 if (id==master) write(*,"(/,a,/)") '--> TESTING STEP MODULE / boundary crossing'

 call init_part()
 npart = 0
 psep = dxbound/50.
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                 psep,hfact,npart,xyzh,periodic,mask=i_belong)
 npartoftype(:) = 0
 npartoftype(1) = npart
 !print*,' thread ',id,' npart = ',npart
 iverbose = 0
 fxyzu(:,:) = 0.
 fext(:,:)  = 0.
 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)

 rhozero = 7.5
 totmass = rhozero/(dxbound*dybound*dzbound)
 massoftype(igas) = totmass/real(reduceall_mpi('+',npart))
 hzero = hfact*(massoftype(igas)/rhozero)**(1./3.)
!
!--set constant velocity (in all components)
!
 vxyzu(1:3,:) = 1.
 if (maxvxyzu>=4) vxyzu(4,:) = 0.
!
!--make sure AV is off
!
 polyk = 0.
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
 call init_eos(ieos,ierr)

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
 nfailed(:) = 0
 call check_setup(nerror,nwarn)
 call checkval(nerror,0,0,nfailed(1),'no errors in setup')
 call update_test_scores(ntests,nfailed,npass)

 fxyzu = 0.
 fext = 0.
 call get_derivs_global()
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
       call checkval(npart,dBevol(4,:),0.,tiny(0.),nfailed(9),'dpsi/dt')
    endif
    call update_test_scores(ntests,nfailed,npass)
 enddo

 if (id==master) write(*,"(/,a)") '<-- STEP TEST COMPLETE'

#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF STEP MODULE (need -DPERIODIC)'

#endif

end subroutine test_step

end module teststep
