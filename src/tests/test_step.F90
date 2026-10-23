!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module teststep
!
! Unit tests of the step module
!
! :References: None
!
! :Owner: Not Committed Yet
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checksetup, damping, deriv, dim, eos,
!   extern_staticsine, externalforces, io, mpidomain, mpiutils, options,
!   part, physcon, ptmass, step_lf_global, testutils, timestep,
!   timestep_ind, timing, unifdis, units, viscosity
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
subroutine test_step(ntests,npass,string)
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
 character(len=*), intent(in), optional :: string
#ifdef PERIODIC
 real                   :: psep,hzero,totmass,dt,t,dtext,dtnew_dum
 real                   :: rhozero
 integer                :: i,nsteps,nerror,nwarn,ierr
 integer :: nfailed(9)

 if (present(string)) then
    if (string=='mass-atmosphere') then
       call test_mass_atmosphere(ntests,npass)
       return
    endif
 endif

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

!-----------------------------------------------------------------------
! Hydrostatic periodic atmosphere: rho=exp(-Phi/cs^2), light particles above
! heavy ones, with a small three-dimensional position perturbation.
!-----------------------------------------------------------------------
subroutine test_mass_atmosphere(ntests,npass)
 use dim, only:periodic,gravity,mpi,use_apr,use_dust,mhd,gr,do_radiation,ind_timesteps,maxvxyzu
 use part, only:init_part,npart,npartoftype,massoftype,particle_mass,variable_mass, &
               xyzh,vxyzu,fxyzu,fext,hfact,igas,set_particle_type,alphaind
 use boundary, only:set_boundary
 use extern_staticsine, only:amplitude,wavek,phase,inclination
 use externalforces, only:iext_staticsine,externalforce
 use ptmass, only:set_integration_precision
 use options, only:ieos,iexternalforce,alpha,alphau,tolh,idamp
 use damping, only:damp
 use eos, only:polyk,gamma,init_eos
 use deriv, only:get_derivs_global
 use step_lf_global, only:init_step,step
 use timestep, only:dtmax
 use testutils, only:checkval,update_test_scores
 use io, only:id,master,iverbose
 use units, only:set_units
 integer, intent(inout) :: ntests,npass
 integer :: level,nx,nt,ix,iy,iz,i,ierr,nfailed(6),verbose_saved,nx_coarse,penetrated,maxpenetrated
 real :: dx,t,dt,dtext,dtnew,mass_initial,maxmach,displacement,delta(3),initial_peak,phi,support
 real :: saved_potential(4)
 real, allocatable :: initial(:,:),masses(:)

 if (.not.periodic.or.gravity.or.mpi.or.use_apr.or.use_dust.or.mhd.or.gr.or.do_radiation.or.ind_timesteps) then
    if (id==master) print*, 'SKIPPING mass-atmosphere: use PERIODIC=yes, GRAVITY=no and gas-only shared timesteps'
    return
 endif
 saved_potential = [amplitude,wavek,phase,inclination]
 verbose_saved = iverbose
 iverbose = 0
 call set_units(1.d0,1.d0,1.d0)
 ieos = 1
 gamma = 1.
 polyk = 1.
 call init_eos(ieos,ierr)
 call set_integration_precision()
 iexternalforce = iext_staticsine
 amplitude = 0.5*log(1000.)
 wavek = 2.*acos(-1.)
 phase = 0.
 inclination = 0.
 hfact = 1.2
 alpha = 1.
 alphau = 0.
 tolh = 1.e-8
 do level=1,3
    call init_part()
    nx = 32*min(level,2)
    if (level==3) nx = 32
    nt = nx/4
    dx = 1./nx
    amplitude = 0.5*log(1000.)
    if (level==3) amplitude = 0.
    call set_boundary(0.,1.,0.,0.25,0.,0.25)
    npart = 0
    nx_coarse = nx
    if (level==3) nx_coarse = nx/2
    do ix=1,nx_coarse
       do iy=1,nt
          do iz=1,nt
             npart = npart + 1
             delta = 1.e-3*dx*sin(real(npart)*[13.,17.,23.])
             xyzh(:,npart) = [(ix-0.5)*dx+delta(1),(iy-0.5)*dx+delta(2), &
                              (iz-0.5)*dx+delta(3),hfact*dx]
             particle_mass(npart) = exp(-amplitude*cos(wavek*xyzh(1,npart)))*dx**3
             call set_particle_type(npart,igas)
          enddo
       enddo
    enddo
    if (level==3) then
       ! Equal physical density across a sharp 8:1 particle-mass interface.
       do ix=1,nx
          do iy=1,2*nt
             do iz=1,2*nt
                npart = npart + 1
                delta = 5.e-4*dx*sin(real(npart)*[13.,17.,23.])
                xyzh(:,npart) = [0.5+(ix-0.5)*dx/2.+delta(1),(iy-0.5)*dx/2.+delta(2), &
                                 (iz-0.5)*dx/2.+delta(3),hfact*dx/2.]
                particle_mass(npart) = (dx/2.)**3
                call set_particle_type(npart,igas)
             enddo
          enddo
       enddo
    endif
    npartoftype = 0
    npartoftype(igas) = npart
    variable_mass = .true.
    mass_initial = sum(particle_mass(1:npart))
    massoftype(igas) = mass_initial/npart
    vxyzu = 0.
    if (maxvxyzu>=4) vxyzu(4,1:npart) = 1.5
    alphaind = 1.
    initial = xyzh(1:3,1:npart)
    masses = particle_mass(1:npart)
    fext = 0.
    fxyzu = 0.
    call get_derivs_global()
    dtmax = 0.1*dx
    if (level==3) dtmax = 0.05*dx
    t = 0.
    call init_step(npart,t,dtmax)
    maxmach = 0.
    initial_peak = 0.
    maxpenetrated = 0
    idamp = 1
    damp = 5.
    do while (t < 2.5)
       dt = min(dtmax,2.5-t)
       if (t < 0.5) dt = min(dt,0.5-t)
       dtext = dt
       call step(npart,npart,t,dt,dtext,dtnew)
       t = t + dt
       if (idamp/=0) then
          initial_peak = max(initial_peak,maxval(sqrt(sum(vxyzu(1:3,1:npart)**2,dim=1))))
          if (t>=0.5) then
             ! Turn the existing drag off entirely and seed non-axisymmetric motion.
             idamp = 0
             damp = 0.
             do i=1,npart
                vxyzu(1:3,i) = vxyzu(1:3,i) + 1.e-3*sin(real(i)*[13.,17.,23.])
                call externalforce(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),t, &
                                    fext(1,i),fext(2,i),fext(3,i),phi)
             enddo
             call get_derivs_global()
          endif
       else
          maxmach = max(maxmach,maxval(sqrt(sum(vxyzu(1:3,1:npart)**2,dim=1))))
       endif
       if (level==3) then
          penetrated = 0
          do i=1,npart
             support = 2.*hfact*dx
             if (i>nx_coarse*nt*nt) support = support/2.
             if (i<=nx_coarse*nt*nt) then
                if (xyzh(1,i)>0.5+support.and.xyzh(1,i)<1.-support) penetrated = penetrated+1
             else
                if (xyzh(1,i)>support.and.xyzh(1,i)<0.5-support) penetrated = penetrated+1
             endif
          enddo
          maxpenetrated = max(maxpenetrated,penetrated)
       endif
    enddo
    displacement = 0.
    do i=1,npart
       delta = xyzh(1:3,i)-initial(:,i)
       delta(1) = modulo(delta(1)+0.5,1.)-0.5
       delta(2:3) = modulo(delta(2:3)+0.125,0.25)-0.125
       displacement = displacement + sum(delta**2)
    enddo
    displacement = sqrt(displacement/npart)
    if (id==master) print*, 'Atmosphere N, max Mach, RMS displacement/spacing: ',npart,maxmach,displacement/dx
    if (id==master) print*, 'Equilibration max Mach, deep interface crossings: ',initial_peak,maxpenetrated
    nfailed = 0
    call checkval(ierr,0,0,nfailed(1),'atmosphere EOS initialises')
    call checkval(maxmach < 0.1,.true.,nfailed(2),'hydrostatic unequal-mass atmosphere remains subsonic')
    call checkval(displacement < 0.5*dx,.true.,nfailed(3),'sub-cell RMS displacement in static atmosphere')
    call checkval(all(particle_mass(1:npart)==masses),.true.,nfailed(4),'atmosphere particle masses stay fixed')
    call checkval(npart,size(masses),0,nfailed(5),'atmosphere particle count stays fixed')
    if (level==3) call checkval(maxpenetrated,0,0,nfailed(6),'no penetration beyond the initial kernel support')
    call update_test_scores(ntests,nfailed,npass)
    deallocate(initial,masses)
 enddo
 amplitude = saved_potential(1)
 wavek = saved_potential(2)
 phase = saved_potential(3)
 inclination = saved_potential(4)
 iverbose = verbose_saved
 call init_part()

end subroutine test_mass_atmosphere

end module teststep
