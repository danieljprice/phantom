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
! :Owner: Francesco Radica
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checksetup, damping, deriv, dim, eos,
!   extern_staticsine, externalforces, io, kernel, mpidomain, mpiutils,
!   options, part, physcon, ptmass, shock_capturing, step_lf_global, testutils,
!   timestep, timestep_ind, timing, unifdis, units, viscosity
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
 use kernel,          only:kernelname
 use mpidomain,       only:i_belong
 use checksetup,      only:check_setup
 use deriv,           only:get_derivs_global
#ifdef IND_TIMESTEPS
 use part,            only:ibin
 use timestep_ind,    only:nbinmax
#endif
 real                   :: psep,hzero,totmass,dt,t,dtext,dtnew_dum
 real                   :: rhozero,tol_h
 integer                :: i,nsteps,nerror,nwarn,ierr
 integer :: nfailed(9)
#endif
 integer, intent(inout) :: ntests,npass
 character(len=*), intent(in), optional :: string

 if (present(string)) then
    if (index(string,'mass-atmosphere')==1) then
       call test_mass_atmosphere(ntests,npass,string)
       return
    endif
 endif

#ifdef PERIODIC
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
 ! Wendland kernels (used with APR) have larger errors on a cubic lattice
 tol_h = 3.e-4
 if (index(kernelname,'Wendland') /= 0) tol_h = 4.e-3
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
    call checkval(npart,xyzh(4,:),hzero,tol_h,nfailed(1),'h (density)')
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
subroutine test_mass_atmosphere(ntests,npass,string)
 use dim, only:periodic,gravity,mpi,use_apr,use_dust,mhd,gr,do_radiation,ind_timesteps,maxvxyzu
 use part, only:init_part,npart,npartoftype,massoftype,particle_mass,variable_mass, &
               xyzh,vxyzu,fxyzu,fext,hfact,igas,set_particle_type,alphaind,rho,eos_vars,ics
 use boundary, only:set_boundary
 use extern_staticsine, only:amplitude,wavek,phase,inclination
 use externalforces, only:iext_staticsine,externalforce
 use ptmass, only:set_integration_precision
 use options, only:ieos,iexternalforce,alpha,alphau,tolh,idamp
 use damping, only:damp
 use eos, only:polyk,gamma,init_eos
 use deriv, only:get_derivs_global
 use step_lf_global, only:init_step,step
 use timestep, only:dtmax,dtcourant,dtforce
 use testutils, only:checkval,update_test_scores
 use io, only:id,master,iverbose
 use units, only:set_units
 use unifdis, only:set_unifdis
 use kernel, only:kernelname,hfact_default,radkern
 integer, intent(inout) :: ntests,npass
 character(len=*), intent(in) :: string
 integer :: level,nx,nt,ix,iy,iz,i,ierr,nfailed(8),verbose_saved,nx_coarse,penetrated,maxpenetrated,first_case,last_case
 real :: dx,t,dt,dtext,dtnew,mass_initial,maxmach,displacement,delta(3),initial_peak,phi,support
 real :: saved_potential(4)
 real :: lengths(3),displacement_xyz(3),release_mach,release_force,release_t,initial_displacement,profile_error,err
 real :: bulk_mach,bulk_force,profile_relax,kinetic_ratio
 real, allocatable :: initial(:,:),masses(:),released(:,:)
 logical :: equilibrated
 real, parameter :: max_relax_time = 30.

 if (.not.periodic .or. gravity .or. mpi .or. use_apr .or. use_dust .or. mhd .or. gr .or. &
     do_radiation .or. ind_timesteps .or. kernelname/='Wendland 2/3D C^4') then
    if (id==master) print*, 'ERROR: mass-atmosphere requires PERIODIC=yes, GRAVITY=no, WendlandC4 and gas-only shared timesteps'
    ntests = ntests + 1
    return
 endif
 if (index(string,'lattice')>0) then
    call test_mass_lattice(ntests,npass)
    return
 endif
 saved_potential = [amplitude,wavek,phase,inclination]
 verbose_saved = iverbose
 iverbose = 0
 call set_units(1.d0,1.d0,1.d0)
 ieos = 1
 gamma = 1.
 if (index(string,'isothermal')==0.and.maxvxyzu>=4) then
    ieos = 2
    gamma = 5./3.
 endif
 polyk = 1.
 call init_eos(ieos,ierr)
 if (id==master) print*, 'Atmosphere EOS, gamma: ',ieos,gamma
 call set_integration_precision()
 iexternalforce = iext_staticsine
 amplitude = 0.5*log(1000.)
 wavek = 2.*acos(-1.)
 phase = 0.
 inclination = 0.
 hfact = hfact_default
 alpha = 1.
 alphau = 0.
 tolh = 1.e-8
 first_case = 1
 last_case = 7
 if (index(string,'small')>0.or.index(string,'resolved')>0.or.index(string,'convergence')>0) last_case = 1
 if (index(string,'control')>0) then
    first_case = 6
    last_case = 7
 endif
 if (index(string,'interface')>0) then
    first_case = 3
    last_case = 3
 endif
 do level=first_case,last_case
    call init_part()
    nx = 32*min(level,2)
    if (level==3) nx = 32
    if (level==4.or.level==5) nx = 32*(level-3)
    if (index(string,'convergence')>0) nx = 64
    if (index(string,'resolved')>0) nx = 128
    if (index(string,'convergence')>0) hfact = hfact_default + 0.3*log(real(nx)/32.)/log(2.)
    nt = nx/4
    dx = 1./nx
    amplitude = 0.5*log(1000.)
    if (level==3.or.level>=6) amplitude = 0.
    lengths = [1.,0.25,0.25]
    if (level==4.or.level==5) lengths = [1.,sqrt(3.)/8.,sqrt(2./3.)/4.]
    call set_boundary(0.,lengths(1),0.,lengths(2),0.,lengths(3))
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
    if (level==4.or.level==5) then
       npart = 0
       call set_unifdis('hcp',id,master,0.,lengths(1),0.,lengths(2),0.,lengths(3), &
                        dx,hfact,npart,xyzh,.true.,verbose=.false.)
       do i=1,npart
          xyzh(1:3,i) = xyzh(1:3,i) + 1.e-3*dx*sin(real(i)*[13.,17.,23.])
          particle_mass(i) = exp(-amplitude*cos(wavek*xyzh(1,i)))*product(lengths)/npart
          call set_particle_type(i,igas)
       enddo
    endif
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
    variable_mass = level/=6
    mass_initial = sum(particle_mass(1:npart))
    massoftype(igas) = mass_initial/npart
    vxyzu = 0.
    if (maxvxyzu>=4) vxyzu(4,1:npart) = 1.5
    alphaind = 1.
    initial = xyzh(1:3,1:npart)
    released = initial
    masses = particle_mass(1:npart)
    fext = 0.
    fxyzu = 0.
    call get_derivs_global()
    if (ieos==2) then
       ! Assign the analytic entropy, not an entropy polluted by the initial
       ! density-sum error, before adiabatic relaxation preserves it.
       do i=1,npart
          phi = exp(-amplitude*cos(wavek*xyzh(1,i)))
          vxyzu(4,i) = (rho(i)/phi)**(gamma-1.)/(gamma-1.)
       enddo
       call get_derivs_global()
    endif
    dtmax = 0.1*dx
    if (level==3) dtmax = 0.05*dx
    if (id==master) print*, 'Atmosphere case, N, hfact, dtmax: ',level,npart,hfact,dtmax
    t = 0.
    call init_step(npart,t,dtmax)
    maxmach = 0.
    initial_peak = 0.
    maxpenetrated = 0
    release_mach = 0.
    release_force = 0.
    bulk_mach = 0.
    bulk_force = 0.
    profile_relax = 0.
    kinetic_ratio = huge(kinetic_ratio)
    profile_error = 0.
    equilibrated = .false.
    release_t = max_relax_time
    idamp = 1
    damp = 1.
    do while (t < release_t+2.)
       dt = min(dtmax,dtcourant,dtforce,release_t+2.-t)
       dtext = dt
       call step(npart,npart,t,dt,dtext,dtnew)
       t = t + dt
       if (.not.equilibrated) then
          release_mach = maxval(sqrt(sum(vxyzu(1:3,1:npart)**2,dim=1))/eos_vars(ics,1:npart))
          initial_peak = max(initial_peak,release_mach)
          release_force = 0.
          bulk_force = 0.
          err = 0.
          do i=1,npart
             call externalforce(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),t, &
                                delta(1),delta(2),delta(3),phi)
             delta = fxyzu(1:3,i)+delta
             release_force = max(release_force,sqrt(sum(delta**2)))
             bulk_force = bulk_force + particle_mass(i)*sum(delta**2)
             err = err + particle_mass(i)*(rho(i)*exp(phi)-1.)**2
          enddo
          bulk_mach = sqrt(sum(particle_mass(1:npart)*sum(vxyzu(1:3,1:npart)**2,dim=1) &
                           /eos_vars(ics,1:npart)**2)/mass_initial)
          bulk_force = sqrt(bulk_force/mass_initial)
          profile_relax = sqrt(err/mass_initial)
          kinetic_ratio = 0.5*sum(particle_mass(1:npart)*sum(vxyzu(1:3,1:npart)**2,dim=1)) &
                          /sum(particle_mass(1:npart)*eos_vars(ics,1:npart)**2)
          ! Use a mass-weighted kinetic criterion, as in stellar relaxation,
          ! while retaining pointwise velocity and force bounds for light tracers.
          if (t>=0.5.and.kinetic_ratio<1.e-7.and.release_mach<1.e-2.and.release_force<1.e-2) then
             equilibrated = .true.
             release_t = t
             released = xyzh(1:3,1:npart)
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
          maxmach = max(maxmach,maxval(sqrt(sum(vxyzu(1:3,1:npart)**2,dim=1))/eos_vars(ics,1:npart)))
          err = 0.
          do i=1,npart
             phi = exp(-amplitude*cos(wavek*xyzh(1,i)))
             err = err + particle_mass(i)*(rho(i)/phi-1.)**2
          enddo
          profile_error = max(profile_error,sqrt(err/mass_initial))
       endif
       penetrated = 0
       do i=1,npart
          support = radkern*hfact*dx
          if (level==3) then
             if (i>nx_coarse*nt*nt) support = support/2.
             if (i<=nx_coarse*nt*nt) then
                if (xyzh(1,i)>0.5+support.and.xyzh(1,i)<1.-support) penetrated = penetrated+1
             else
                if (xyzh(1,i)>support.and.xyzh(1,i)<0.5-support) penetrated = penetrated+1
             endif
          else
             phi = abs(modulo(xyzh(1,i)-initial(1,i)+0.5,1.)-0.5)
             if (phi>support) penetrated = penetrated+1
          endif
       enddo
       maxpenetrated = max(maxpenetrated,penetrated)
       if (int(t)>int(t-dt).and.id==master) print*, 'Atmosphere case, t, relaxed, Mach, force, crossings: ', &
                                                  level,t,equilibrated,release_mach,release_force,maxpenetrated
       if (int(t)>int(t-dt).and.id==master) print*, 'Bulk Mach, bulk force, density RMS error: ', &
                                                  bulk_mach,bulk_force,profile_relax
       if (maxpenetrated>0) exit
       if (.not.equilibrated.and.t>=max_relax_time) exit
    enddo
    idamp = 0
    damp = 0.
    displacement = 0.
    displacement_xyz = 0.
    initial_displacement = 0.
    do i=1,npart
       delta = xyzh(1:3,i)-initial(:,i)
       delta = modulo(delta+0.5*lengths,lengths)-0.5*lengths
       initial_displacement = initial_displacement + sum(delta**2)
       delta = xyzh(1:3,i)-released(:,i)
       delta = modulo(delta+0.5*lengths,lengths)-0.5*lengths
       displacement = displacement + sum(delta**2)
       displacement_xyz = displacement_xyz + delta**2
    enddo
    displacement = sqrt(displacement/npart)
    if (id==master) print*, 'Atmosphere N, max Mach, RMS displacement/spacing: ',npart,maxmach,displacement/dx
    if (id==master) print*, 'Equilibration max Mach, deep interface crossings: ',initial_peak,maxpenetrated
    if (id==master) print*, 'Case, RMS displacement xyz/spacing: ',level,sqrt(displacement_xyz/npart)/dx
    if (id==master) print*, 'Release max Mach, max residual acceleration: ',release_mach,release_force
    if (id==master) print*, 'Bulk Mach, bulk force, relaxed density RMS error: ',bulk_mach,bulk_force,profile_relax
    if (id==master) print*, 'Released kinetic/pressure energy ratio: ',kinetic_ratio
    if (id==master) print*, 'Equilibrated, release time, density error, total shift/spacing: ', &
                           equilibrated,release_t,profile_error,sqrt(initial_displacement/npart)/dx
    nfailed = 0
    call checkval(ierr,0,0,nfailed(1),'atmosphere EOS initialises')
    call checkval(equilibrated.and.maxmach < 0.1,.true.,nfailed(2),'released atmosphere stays below Mach 0.1')
    call checkval(equilibrated.and.displacement < 0.5*dx,.true.,nfailed(3),'sub-cell RMS displacement after equilibrium release')
    call checkval(all(particle_mass(1:npart)==masses),.true.,nfailed(4),'atmosphere particle masses stay fixed')
    call checkval(npart,size(masses),0,nfailed(5),'atmosphere particle count stays fixed')
    call checkval(maxpenetrated,0,0,nfailed(6),'no penetration beyond nominal initial support during relaxation or release')
    call checkval(equilibrated,.true.,nfailed(7),'atmosphere reaches discrete force equilibrium before release')
    call checkval(equilibrated.and.profile_error < 0.05,.true.,nfailed(8),'released mass-density profile stays within five percent')
    call update_test_scores(ntests,nfailed,npass)
    deallocate(initial,masses,released)
 enddo
 amplitude = saved_potential(1)
 wavek = saved_potential(2)
 phase = saved_potential(3)
 inclination = saved_potential(4)
 hfact = hfact_default
 iverbose = verbose_saved
 call init_part()

end subroutine test_mass_atmosphere

subroutine test_mass_lattice(ntests,npass)
 use part, only:init_part,npart,npartoftype,massoftype,particle_mass,variable_mass, &
               xyzh,vxyzu,fxyzu,fext,hfact,igas,set_particle_type,alphaind,rho
 use boundary, only:set_boundary
 use options, only:ieos,iexternalforce,alpha,alphau,tolh
 use shock_capturing, only:beta
 use eos, only:polyk,gamma,init_eos
 use deriv, only:get_derivs_global
 use kernel, only:kernelname,hfact_default
 use testutils, only:checkval,update_test_scores
 use physcon, only:pi
 use units, only:set_units
 integer, intent(inout) :: ntests,npass
 integer :: ix,iy,iz,i,iside,ilevel,ierr,nfailed(4),eos_saved,external_saved
 real :: saved(7),base(4,512),mode(3,512),accel(3,512,2),energy(2),energy0,norm,curvature,energy_curve
 real :: delta,force_curve(2),thermal_curve(2)

 if (kernelname/='Wendland 2/3D C^4') then
    ntests = ntests + 1
    return
 endif
 saved = [hfact,gamma,polyk,alpha,alphau,tolh,beta]
 eos_saved = ieos
 external_saved = iexternalforce
 call init_part()
 call set_units(1.d0,1.d0,1.d0)
 call set_boundary(0.,8.,0.,8.,0.,8.)
 hfact = hfact_default
 ieos = 1
 iexternalforce = 0
 gamma = 1.
 polyk = 1.
 alpha = 0.
 alphau = 0.
 beta = 0.
 tolh = 1.e-13
 call init_eos(ieos,ierr)
 npart = 0
 do ix=1,8
    do iy=1,8
       do iz=1,8
          npart = npart+1
          xyzh(:,npart) = [ix-0.5,iy-0.5,iz-0.5,hfact]
          particle_mass(npart) = 1.
          call set_particle_type(npart,igas)
       enddo
    enddo
 enddo
 npartoftype = 0
 npartoftype(igas) = npart
 massoftype(igas) = 1.
 variable_mass = .true.
 vxyzu = 0.
 fext = 0.
 alphaind = 0.
 call get_derivs_global()
 base = xyzh(:,1:npart)
 energy0 = sum(log(rho(1:npart)))
 nfailed = 0
 call checkval(ierr,0,0,nfailed(1),'lattice isothermal EOS initialises')
 call checkval(maxval(abs(fxyzu(1:3,1:npart))),0.,1.e-9,nfailed(2),'uniform periodic lattice has zero force')
 do i=1,npart
    mode(:,i) = cos(0.25*pi*dot_product(base(1:3,i),[1.,1.,2.]))*[1.,-1.,0.]/sqrt(2.)
 enddo
 norm = sum(mode**2)
 ! Remove the leading finite-difference error before comparing energy and force.
 do ilevel=1,2
    delta = 1.e-3/ilevel
    do iside=1,2
       xyzh(:,1:npart) = base
       xyzh(1:3,1:npart) = base(1:3,:) + (2*iside-3)*delta*mode
       call get_derivs_global()
       accel(:,:,iside) = fxyzu(1:3,1:npart)
       energy(iside) = sum(log(rho(1:npart)))
    enddo
    force_curve(ilevel) = -sum(mode*(accel(:,:,2)-accel(:,:,1)))/(2.*delta*norm)
    thermal_curve(ilevel) = (sum(energy)-2.*energy0)/(delta**2*norm)
 enddo
 curvature = (4.*force_curve(2)-force_curve(1))/3.
 energy_curve = (4.*thermal_curve(2)-thermal_curve(1))/3.
 print*, 'Transverse lattice curvature, energy curvature: ',curvature,energy_curve
 call checkval(energy_curve,curvature,1.e-5,nfailed(3),'periodic force differentiates isothermal energy')
 call checkval(curvature>=-1.e-7,.true.,nfailed(4),'oblique transverse lattice mode is not unstable')
 call update_test_scores(ntests,nfailed,npass)
 hfact = saved(1)
 gamma = saved(2)
 polyk = saved(3)
 alpha = saved(4)
 alphau = saved(5)
 tolh = saved(6)
 beta = saved(7)
 ieos = eos_saved
 iexternalforce = external_saved
 call init_part()
end subroutine test_mass_lattice

end module teststep
