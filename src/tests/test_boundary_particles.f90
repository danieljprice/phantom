!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testboundary
!
! Unit tests of the boundary particles module
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: boundary, boundarypart, deriv, dim, energies, eos, io,
!   kdtree, options, part, physcon, ptmass, readwrite_dumps, spherical,
!   step_lf_global, stretchmap, testutils, timestep, units, vectorutils
!
 use testutils, only:checkval,update_test_scores
 implicit none
 public :: test_boundary

 private

 ! n=1 polytrope with M = R = G = 1
 real, parameter :: rstar = 1.
 real, parameter :: mstar = 1.

 ! tolerances for rigid-body motion and momentum conservation of the core.
 real, parameter :: tol_vrigid = 1.e-11
 real, parameter :: tol_xrigid = 1.e-12
 real, parameter :: tol_angmom = 1.e-12
 real, parameter :: tol_omega  = 1.e-12
 real, parameter :: tol_return = 1.e-12

 ! indices into the error array returned by run_free_spin
 integer, parameter :: ierr_angmom = 1, ierr_mom = 2, ierr_shape1 = 3, ierr_shape2 = 4, &
                       ierr_vrigid = 5, ierr_omega = 6, ierr_return = 7, ierr_evomega = 8, &
                       nerr_spin = 8

 ! indices into the error array returned by run_known_torque
 integer, parameter :: ierrt_angmom = 1, ierrt_omega = 2, ierrt_vcm = 3, ierrt_xcm = 4, &
                       ierrt_shape = 5, ierrt_vrigid = 6, ierrt_xfinal = 7, ierrt_linear = 8, &
                       nerr_torque = 8

 logical :: write_dumps = .true.
 integer, parameter :: ndumpfreq = 10

contains

subroutine test_boundary(ntests,npass)
 use io, only:id,master
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING BOUNDARY PARTICLES'

 call test_core_rigid_drift(ntests,npass)
 call test_core_star_sink_orbit(ntests,npass)
 call test_boundary_star_gas_star_orbit(ntests,npass)
 call test_core_free_spin(ntests,npass)
 call test_core_known_torque(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- BOUNDARY PARTICLE TEST COMPLETE'

end subroutine test_boundary

!-----------------------------------------------------------------------
!+
!  n=1 polytrope whose inner half (r < R/2) is a rigid core of boundary
!  particles, on a circular Keplerian orbit around a sink particle.
!
!  Checks that
!   - the core stays rigid: velocities are those of a rigid body
!     (v = v_cm + omega x r) and pairwise distances are preserved.
!   - total linear and angular momentum are conserved
!   - total energy is conserved once the thermal energy discarded by
!     freezing u in the core (ecore_lost) is added back
!+
!-----------------------------------------------------------------------
subroutine test_core_star_sink_orbit(ntests,npass)
 use dim,            only:gravity,maxvxyzu,mpi,periodic,gr
 use io,             only:id,master,iverbose
 use physcon,        only:pi
 use part,           only:init_part,npart,npartoftype,massoftype,xyzh,vxyzu,fxyzu,fext,hfact,&
                          igas,iboundary,set_particle_type,iphase,iamtype,iamboundary,&
                          isdead_or_accreted,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,&
                          dsdt_ptmass,ihacc,ihsoft,update_npartoftypetot
 use eos,            only:ieos,gamma,polyk
 use options,        only:alpha
 use eos,            only:ipdv_heating,ishock_heating
 use units,          only:set_units
 use boundary,       only:set_boundary
 use spherical,      only:set_sphere
 use stretchmap,     only:rho_func
 use kdtree,         only:tree_accuracy
 use ptmass,         only:get_accel_sink_gas,h_soft_sinkgas,set_integration_precision
 use deriv,          only:get_derivs_global
 use step_lf_global, only:init_step,step
 use timestep,       only:dtmax,tolv
 use energies,       only:compute_energies,etot
 use boundarypart,   only:get_boundary_particle_forces,calculate_boundary_core_CM
 use checkconserved, only:ecore_lost
 integer, intent(inout) :: ntests,npass
 integer :: i,istep,nsteps,nfailed(6),ncore,nalive
 integer(kind=8) :: nptot
 real    :: psep,rcore,msink,a,omega,xstar(3),vstar(3),xsink(3),vsink(3)
 real    :: t,dt,dtext,dtnew,tmax,phii,dum,dum2,pmassi
 real    :: mom0(3),angmom0(3),mom(3),angmom(3),etot0,pscale,lscale
 real    :: xcore(3),vcore(3),dvmax,dxmax,sep0,sep
 real, allocatable :: xyz0(:,:)
 procedure(rho_func), pointer :: density_func

 if (.not.gravity .or. maxvxyzu < 4 .or. mpi .or. gr) then
    if (id==master) write(*,"(/,a)") &
       ' skipping core-star-sink orbit test (requires GRAVITY=yes, non-isothermal, no MPI, no GR)'
    return
 endif
 if (id==master) write(*,"(/,a)") '--> testing polytrope with boundary-particle core orbiting a sink'

 call init_part()
 call set_units(mass=1.d0,dist=1.d0,G=1.d0)
 iverbose       = 0
 tree_accuracy  = 0.     ! direct summation, so gravity conserves momentum to round-off
 tolv           = 1.e-2
 alpha          = 1.
 ipdv_heating   = 1
 ishock_heating = 1
 if (periodic) call set_boundary(-10.,10.,-10.,10.,-10.,10.)
 !
 ! n=1 polytrope: rho = rhoc sin(pi r/R)/(pi r/R), P = K rho^2 with R = sqrt(pi K/2G)
 !
 ieos  = 2
 gamma = 2.
 polyk = 2.*rstar**2/pi
 !
 ! set up the star centred on the origin
 !
 psep = 0.15
 npartoftype(:) = 0
 npart = 0
 nptot = 0
 density_func => rho_n1_polytrope
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,nptot,rhofunc=density_func)
 massoftype(igas)      = mstar/real(npart)
 massoftype(iboundary) = massoftype(igas)
 !
 ! particles inside half the stellar radius become a boundary-particle core
 !
 rcore = 0.5*rstar
 do i=1,npart
    if (sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i))) < rcore) then
       call set_particle_type(i,iboundary)
    else
       call set_particle_type(i,igas)
    endif
 enddo
 ncore = count(iamboundary(iamtype(iphase(1:npart))))
 npartoftype(iboundary) = ncore
 npartoftype(igas)      = npart - ncore
 call update_npartoftypetot()  ! step only distinguishes particle types if ntypes > 1
 if (id==master) print "(2(a,i6))",' npart = ',npart,' ncore = ',ncore
 !
 ! circular Keplerian orbit of the star around a sink, centre of mass at rest at the origin
 !
 msink = 1.
 a     = 3.  ! Roche lobe of star ~0.38a > R
 omega = sqrt((mstar + msink)/a**3)
 xstar = [-msink/(mstar + msink)*a, 0., 0.]
 xsink = [ mstar/(mstar + msink)*a, 0., 0.]
 vstar = [0., xstar(1)*omega, 0.]
 vsink = [0., xsink(1)*omega, 0.]

 do i=1,npart
    vxyzu(4,i)   = polyk*rho_n1_polytrope(sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i))))/(gamma - 1.)
    xyzh(1:3,i)  = xyzh(1:3,i) + xstar
    vxyzu(1:3,i) = vstar
 enddo

 nptmass = 1
 xyzmh_ptmass(:,1)    = 0.
 xyzmh_ptmass(1:3,1)  = xsink
 xyzmh_ptmass(4,1)    = msink
 xyzmh_ptmass(ihacc,1)  = 0.05
 xyzmh_ptmass(ihsoft,1) = 0.
 h_soft_sinkgas       = 0.
 vxyz_ptmass(:,1)     = vsink
 !
 ! initial sink-gas forces, with the core's sink force averaged as in step_extern
 !
 fext(:,:) = 0.
 fxyz_ptmass(:,1:nptmass) = 0.
 dsdt_ptmass(:,1:nptmass) = 0.
 do i=1,npart
    pmassi = massoftype(iamtype(iphase(i)))
    call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                            fext(1,i),fext(2,i),fext(3,i),phii,pmassi,fxyz_ptmass,dsdt_ptmass,dum,dum2)
 enddo
 call get_boundary_particle_forces(npart,iphase,xyzh(1:3,:),fext)
 !
 ! SPH forces; fixed timestep well below the Courant limit (cs ~ 1, h ~ 0.2)
 !
 t      = 0.
 dt     = 5.e-3
 dtmax  = dt
 tmax   = 1.     ! ~ one stellar dynamical time, ~1/20 of the orbit
 nsteps = nint(tmax/dt)
 fxyzu(:,:) = 0.
 call get_derivs_global()
 !
 ! initial conserved quantities
 !
 call get_momenta(mom0,angmom0)
 call compute_energies(t)
 etot0  = etot
 ecore_lost = 0.   ! accumulated by step
 pscale = mstar*sqrt(dot_product(vstar,vstar))   ! orbital momentum of the star
 lscale = sqrt(dot_product(angmom0,angmom0))
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
 sep0 = sqrt(dot_product(xcore - xyzmh_ptmass(1:3,1),xcore - xyzmh_ptmass(1:3,1)))
 allocate(xyz0(3,npart))
 xyz0 = xyzh(1:3,1:npart)
 !
 ! evolve (set the sink-gas substepping coefficients, as done in initial)
 !
 call set_integration_precision
 call init_step(npart,t,dtmax)
 if (write_dumps) call write_test_dump('boundary_test_',t,0)
 do istep=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    if (write_dumps .and. mod(istep,ndumpfreq)==0) call write_test_dump('boundary_test_',t,istep/ndumpfreq)
 enddo
 !
 ! rigidity of the core
 !
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
 nalive = 0
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) nalive = nalive + 1
 enddo
 call get_core_rigidity(xyz0,dvmax,dxmax)
 deallocate(xyz0)
 sep = sqrt(dot_product(xcore - xyzmh_ptmass(1:3,1),xcore - xyzmh_ptmass(1:3,1)))
 !
 ! conservation
 !
 call get_momenta(mom,angmom)
 call compute_energies(t)

 nfailed(:) = 0
 call checkval(nalive,npart,0,nfailed(1),'no particles accreted')
 call checkval(dvmax/sqrt(dot_product(vcore,vcore)),0.,tol_vrigid,nfailed(2),'core velocities rigid')
 call checkval(dxmax/rcore,0.,tol_xrigid,nfailed(3),'core pairwise distances preserved')
 call checkval(sqrt(dot_product(mom - mom0,mom - mom0))/pscale,0.,1.e-12,nfailed(4),'linear momentum')
 call checkval(sqrt(dot_product(angmom - angmom0,angmom - angmom0))/lscale,0.,tol_angmom,nfailed(5),'angular momentum')
 call checkval(etot + ecore_lost,etot0,1.e-5,nfailed(6),'total energy (incl. core energy lost)')
 call update_test_scores(ntests,nfailed,npass)

 nfailed(:) = 0
 call checkval(sep,sep0,1.e-3,nfailed(1),'core-sink separation (circular orbit)')
 call update_test_scores(ntests,nfailed(1:1),npass)

 ! reset
 nptmass = 0
 iverbose = 0

end subroutine test_core_star_sink_orbit

!-----------------------------------------------------------------------
!+
!  Two identical n=1 polytropes on a circular Keplerian orbit, one made
!  entirely of boundary particles and the other entirely of gas, with
!  no sink particles. The stars are far enough apart that they share no
!  SPH neighbours, so they interact only via self-gravity.
!
!  The boundary star can only feel the gas star through the averaged
!  SPH/self-gravity acceleration (fxyzu) computed in derivs, so this
!  checks that this acceleration is actually applied. With direct
!  summation gravity is pairwise symmetric, so
!   - total linear momentum is conserved to round-off
!   - the separation between the two stars' centres of mass stays fixed
!  We also check the boundary star stays rigid.
!+
!-----------------------------------------------------------------------
subroutine test_boundary_star_gas_star_orbit(ntests,npass)
 use dim,            only:gravity,maxvxyzu,mpi,periodic,gr
 use io,             only:id,master,iverbose
 use part,           only:init_part,npart,npartoftype,massoftype,xyzh,vxyzu,fxyzu,fext,hfact,&
                          igas,iboundary,set_particle_type,iphase,iamtype,iamboundary,&
                          isdead_or_accreted,nptmass,update_npartoftypetot
 use eos,            only:ieos,gamma,polyk
 use options,        only:alpha
 use eos,            only:ipdv_heating,ishock_heating
 use units,          only:set_units
 use boundary,       only:set_boundary
 use physcon,        only:pi
 use kdtree,         only:tree_accuracy
 use deriv,          only:get_derivs_global
 use step_lf_global, only:init_step,step
 use timestep,       only:dtmax,tolv
 use energies,       only:compute_energies,etot
 use boundarypart,   only:calculate_boundary_core_CM
 integer, intent(inout) :: ntests,npass
 integer :: i,istep,nsteps,nfailed(6),nbound,ngas,nalive
 real    :: a,omega,x1(3),x2(3),v1(3),v2(3),t,dt,dtext,dtnew,tmax
 real    :: mom0(3),angmom0(3),mom(3),angmom(3),etot0,pscale,lscale
 real    :: xbound(3),vbound(3),xgas(3),vgas(3),sep0,sep,dvmax,dxmax
 real, allocatable :: xyz0(:,:)

 if (.not.gravity .or. maxvxyzu < 4 .or. mpi .or. gr) then
    if (id==master) write(*,"(/,a)") &
       ' skipping boundary star-gas star orbit test (requires GRAVITY=yes, non-isothermal, no MPI, no GR)'
    return
 endif
 if (id==master) write(*,"(/,a)") '--> testing boundary-particle star orbiting a gas star'

 call init_part()
 call set_units(mass=1.d0,dist=1.d0,G=1.d0)
 iverbose       = 0
 tree_accuracy  = 0.     ! direct summation, so gravity conserves momentum to round-off
 tolv           = 1.e-2
 alpha          = 1.
 ipdv_heating   = 1
 ishock_heating = 1
 if (periodic) call set_boundary(-10.,10.,-10.,10.,-10.,10.)
 ieos  = 2
 gamma = 2.
 polyk = 2.*rstar**2/pi
 !
 ! circular orbit of two equal-mass stars about the origin; separation is
 ! large enough that each star is well inside its Roche lobe (~0.38a) and
 ! the stars share no SPH neighbours
 !
 a     = 3.
 omega = sqrt(2.*mstar/a**3)
 x1    = [-0.5*a, 0., 0.]
 x2    = [ 0.5*a, 0., 0.]
 v1    = [0., x1(1)*omega, 0.]
 v2    = [0., x2(1)*omega, 0.]
 !
 ! star 1 is made of boundary particles, star 2 of gas
 !
 npart = 0
 npartoftype(:) = 0
 call add_polytrope(x1,v1,iboundary)
 nbound = npart
 call add_polytrope(x2,v2,igas)
 ngas   = npart - nbound
 npartoftype(iboundary) = nbound
 npartoftype(igas)      = ngas
 massoftype(iboundary)  = mstar/real(nbound)
 massoftype(igas)       = mstar/real(ngas)
 call update_npartoftypetot()
 if (id==master) print "(2(a,i6))",' nboundary = ',nbound,' ngas = ',ngas

 nptmass   = 0
 fext(:,:) = 0.
 !
 ! SPH forces; fixed timestep well below the Courant limit (cs ~ 1, h ~ 0.2)
 !
 t      = 0.
 dt     = 5.e-3
 dtmax  = dt
 tmax   = 1.     ! ~ one stellar dynamical time, ~1/20 of the orbit
 nsteps = nint(tmax/dt)
 fxyzu(:,:) = 0.
 call get_derivs_global()
 !
 ! initial conserved quantities and orbit
 !
 call get_momenta(mom0,angmom0)
 call compute_energies(t)
 etot0  = etot
 pscale = mstar*sqrt(dot_product(v1,v1))   ! orbital momentum of one star
 lscale = sqrt(dot_product(angmom0,angmom0))
 call get_star_CMs(xbound,vbound,xgas,vgas)
 sep0 = sqrt(dot_product(xbound - xgas,xbound - xgas))
 allocate(xyz0(3,npart))
 xyz0 = xyzh(1:3,1:npart)
 !
 ! evolve
 !
 call init_step(npart,t,dtmax)
 if (write_dumps) call write_test_dump('boundary_twostar_',t,0)
 do istep=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    if (write_dumps .and. mod(istep,ndumpfreq)==0) call write_test_dump('boundary_twostar_',t,istep/ndumpfreq)
 enddo
 !
 ! rigidity of the boundary star
 !
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xbound,vbound)
 nalive = 0
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) nalive = nalive + 1
 enddo
 call get_core_rigidity(xyz0,dvmax,dxmax)
 deallocate(xyz0)
 !
 ! orbit and conservation
 !
 call get_star_CMs(xbound,vbound,xgas,vgas)
 sep = sqrt(dot_product(xbound - xgas,xbound - xgas))
 call get_momenta(mom,angmom)
 call compute_energies(t)

 nfailed(:) = 0
 call checkval(nalive,npart,0,nfailed(1),'no particles lost')
 call checkval(dvmax/sqrt(dot_product(vbound,vbound)),0.,tol_vrigid,nfailed(2),'boundary star velocities rigid')
 call checkval(dxmax/rstar,0.,tol_xrigid,nfailed(3),'boundary star pairwise distances preserved')
 call checkval(sqrt(dot_product(mom - mom0,mom - mom0))/pscale,0.,1.e-12,nfailed(4),'linear momentum')
 call checkval(sqrt(dot_product(angmom - angmom0,angmom - angmom0))/lscale,0.,tol_angmom,nfailed(5),'angular momentum')
 call checkval(etot,etot0,1.e-2,nfailed(6),'total energy')
 call update_test_scores(ntests,nfailed,npass)

 nfailed(:) = 0
 call checkval(sep,sep0,1.e-3,nfailed(1),'boundary star-gas star separation (circular orbit)')
 call update_test_scores(ntests,nfailed(1:1),npass)

 iverbose = 0

end subroutine test_boundary_star_gas_star_orbit

!-----------------------------------------------------------------------
!+
!  A boundary-particle sphere spinning about a principal axis, evolved
!  with the full step (kicks from derivs, rigid drift of the core).
!
!  Checks that
!   - spin angular momentum and linear momentum are conserved
!   - the core stays rigid (pairwise distances after one and two periods,
!     rigid velocities) and omega stays constant
!   - after one period particles return to their initial positions
!   - all of these hold to round-off at dt/2 as well, i.e. the error
!     does not depend on dt
!   - adding a distant sink of negligible mass, which forces the
!     sink-gas substepping path (kickdrift with fext), gives the same
!     result. This checks the core is drifted exactly once per substep
!   - the core angular velocity written to the .ev file (omegacore_*)
!     matches the known spin
!+
!-----------------------------------------------------------------------
subroutine test_core_free_spin(ntests,npass)
 use dim,     only:mpi,gr
 use io,      only:id,master
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nper = 157     ! steps per rotation period at the base timestep
 real,    parameter :: omega0 = 0.5   ! in units of the break-up spin sqrt(GM/R^3) = 1
 real,    parameter :: pi = acos(-1.)
 integer :: nfailed(10)
 real    :: period,dt,err1(nerr_spin),err2(nerr_spin),errs(nerr_spin),dxsink
 real, allocatable :: xfinal1(:,:),xfinals(:,:)

 if (mpi .or. gr) then
    if (id==master) write(*,"(/,a)") ' skipping free-spin test (requires no MPI, no GR)'
    return
 endif
 if (id==master) write(*,"(/,a)") '--> testing free spin of an isolated boundary-particle core'

 period = 2.*pi/omega0
 dt     = period/real(nper)
 if (id==master) print "(a,f6.3)",' omega*dt = ',omega0*dt

 call run_free_spin(omega0,dt,nper,.false.,write_dumps,err1,xfinal1)
 call run_free_spin(omega0,0.5*dt,2*nper,.false.,.false.,err2)
 call run_free_spin(omega0,dt,nper,.true.,.false.,errs,xfinals)
 dxsink = maxval(abs(xfinals - xfinal1))/rstar

 nfailed(:) = 0
 call checkval(err1(ierr_angmom),0.,tol_angmom,nfailed(1),'spin angular momentum')
 call checkval(err1(ierr_mom),0.,tol_angmom,nfailed(2),'linear momentum')
 call checkval(max(err1(ierr_shape1),err1(ierr_shape2)),0.,tol_xrigid,nfailed(3),&
               'core pairwise distances preserved')
 call checkval(err1(ierr_vrigid),0.,tol_vrigid,nfailed(4),'core velocities rigid')
 call checkval(err1(ierr_omega),0.,tol_omega,nfailed(5),'omega constant')
 call checkval(err1(ierr_return),0.,tol_return,nfailed(6),'return to start after one period')
 call checkval(err1(ierr_evomega),0.,tol_omega,nfailed(10),'core omega in .ev file')
 call checkval(max(err2(ierr_shape2)/tol_xrigid,err2(ierr_vrigid)/tol_vrigid,&
                   err2(ierr_omega)/tol_omega,err2(ierr_return)/tol_return),0.,1.,nfailed(7),&
               'round-off at dt/2 too (relative to tolerances)')
 call checkval(max(err2(ierr_angmom),err2(ierr_mom)),0.,tol_angmom,nfailed(8),'momenta at dt/2')
 ! the sink's (real) tidal pull on the core is ~1e-10 and sets this tolerance
 call checkval(dxsink,0.,1.e-10,nfailed(9),'same result with sink substepping')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_core_free_spin

!-----------------------------------------------------------------------
!+
!  Unit test of the rigid-rotation drift (drift_boundary_core) on its
!  own, with no forces: the core is a free rigid rotor, translating at
!  constant velocity.
!
!  Case 1, spherical core (cubic lattice, isotropic inertia), spinning
!  about an oblique axis: the drift is exact, so after one period every
!  particle is back at its initial position (shifted by V*t), and
!  pairwise distances, rigid velocities, linear and angular momentum
!  are all preserved to round-off, independent of dt.
!
!  Case 2, prolate symmetric top (lattice stretched along z), with
!  omega not along a principal axis: the symmetry axis e3 should
!  precess about the fixed L at the rate |L|/I1, so that
!  I(t) = I1 1 + (I3 - I1) e3(t) e3(t)^T. Shape, rigid velocities and
!  momenta are still exact, and the error in I(t) should converge at
!  second order in dt (explicit midpoint rotor).
!+
!-----------------------------------------------------------------------
subroutine test_core_rigid_drift(ntests,npass)
 use dim,     only:mpi
 use io,      only:id,master
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(9)
 real    :: errs(5),errt1(5),errt2(5)

 if (mpi) then
    if (id==master) write(*,"(/,a)") ' skipping rigid drift test (requires no MPI)'
    return
 endif
 if (id==master) write(*,"(/,a)") '--> testing rigid-rotation drift of a free boundary-particle core'

 call run_rigid_drift(1.,100,errs)
 call run_rigid_drift(1.5,100,errt1)
 call run_rigid_drift(1.5,200,errt2)

 nfailed(:) = 0
 call checkval(errs(1),0.,tol_angmom,nfailed(1),'sphere: angular momentum')
 call checkval(errs(2),0.,tol_angmom,nfailed(2),'sphere: linear momentum')
 call checkval(errs(3),0.,1.e-13,nfailed(3),'sphere: pairwise distances preserved')
 call checkval(errs(4),0.,1.e-13,nfailed(4),'sphere: velocities rigid')
 call checkval(errs(5),0.,1.e-12,nfailed(5),'sphere: exact return after one period')
 call checkval(max(errt1(1),errt1(2)),0.,tol_angmom,nfailed(6),'top: linear and angular momentum')
 call checkval(max(errt1(3),errt1(4)),0.,1.e-13,nfailed(7),'top: pairwise distances, rigid velocities')
 call checkval(errt1(5),0.,1.e-3,nfailed(8),'top: precession matches analytic')
 call checkval(errt1(5)/errt2(5),4.,0.5,nfailed(9),'top: 2nd order convergence of precession')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_core_rigid_drift

!-----------------------------------------------------------------------
!+
!  set up a free core of boundary particles, with the lattice stretched
!  by zstretch along z, moving at constant velocity and spinning, and
!  apply nstep rigid drifts covering one rotation period (sphere) or
!  one precession period (top). Returns errors in
!  (angular momentum, linear momentum, pairwise distances,
!   velocity rigidity, final positions or inertia tensor)
!+
!-----------------------------------------------------------------------
subroutine run_rigid_drift(zstretch,nstep,err)
 use io,           only:id,master
 use part,         only:init_part,npart,npartoftype,massoftype,xyzh,vxyzu,iboundary,&
                        iphase,update_npartoftypetot
 use boundarypart, only:calculate_boundary_core_CM,drift_boundary_core
 use vectorutils,  only:cross_product3D
 real,    intent(in)  :: zstretch
 integer, intent(in)  :: nstep
 real,    intent(out) :: err(5)
 real, parameter :: omega0 = 0.5
 ! physcon's pi has 11 digits, which would make the reference period
 ! overshoot one turn by 2e-11 rad; the drift itself is exact to round-off
 real, parameter :: pi = acos(-1.)
 real    :: vcen(3),omegavec0(3),xcore(3),vcore(3),spin(3),inertia0(3,3),inertia(3,3)
 real    :: mom0(3),angmom0(3),mom(3),angmom(3),e3(3),lhat(3),lxe3(3),inertia_an(3,3)
 real    :: tend,dt,dvmax,dxmax,dxret,theta,smag
 integer :: i,istep,k
 real, allocatable :: xyz0(:,:)

 call init_part()
 npart = 0
 npartoftype(:) = 0
 vcen = [0.1,-0.2,0.3]
 call add_polytrope([0.,0.,0.],vcen,iboundary,psep_in=0.2)
 npartoftype(iboundary) = npart
 massoftype(iboundary)  = mstar/real(npart)
 call update_npartoftypetot()
 xyzh(3,1:npart) = zstretch*xyzh(3,1:npart)

 if (zstretch == 1.) then
    omegavec0 = omega0*[1.,2.,2.]/3.               ! oblique axis, any axis is principal
 else
    omegavec0 = omega0*[sin(pi/6.),0.,cos(pi/6.)]  ! 30 degrees off the symmetry axis
 endif
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
 do i=1,npart
    call cross_product3D(omegavec0,xyzh(1:3,i) - xcore,vxyzu(1:3,i))
    vxyzu(1:3,i) = vxyzu(1:3,i) + vcen
 enddo
 call get_free_core_moments(spin,inertia0)
 smag = sqrt(dot_product(spin,spin))
 if (zstretch == 1.) then
    tend = 2.*pi/omega0
 else
    tend = 2.*pi*inertia0(1,1)/smag                ! precession period of e3 about L
 endif
 dt = tend/real(nstep)

 call get_momenta(mom0,angmom0)
 allocate(xyz0(3,npart))
 xyz0 = xyzh(1:3,1:npart)
 do istep=1,nstep
    call drift_boundary_core(dt,npart,iphase,xyzh,vxyzu)
 enddo
 call get_momenta(mom,angmom)
 call get_core_rigidity(xyz0,dvmax,dxmax)

 err(1) = sqrt(dot_product(angmom - angmom0,angmom - angmom0))/sqrt(dot_product(angmom0,angmom0))
 err(2) = sqrt(dot_product(mom - mom0,mom - mom0))/sqrt(dot_product(mom0,mom0))
 err(3) = dxmax/rstar
 err(4) = dvmax/(omega0*rstar)
 if (zstretch == 1.) then
    ! exact rotor: every particle returns to its start, shifted by V*t
    dxret = 0.
    do i=1,npart
       dxret = max(dxret,sqrt(dot_product(xyzh(1:3,i) - xyz0(:,i) - vcen*tend,&
                                          xyzh(1:3,i) - xyz0(:,i) - vcen*tend)))
    enddo
    err(5) = dxret/rstar
 else
    ! analytic symmetric top: e3 = z axis initially, rotated about L by |L|t/I1
    call get_free_core_moments(spin,inertia)
    lhat  = spin/smag
    theta = smag*tend/inertia0(1,1)
    e3    = [0.,0.,1.]
    call cross_product3D(lhat,e3,lxe3)
    e3 = e3*cos(theta) + lxe3*sin(theta) + lhat*dot_product(lhat,e3)*(1. - cos(theta))
    do k=1,3
       inertia_an(:,k) = (inertia0(3,3) - inertia0(1,1))*e3(:)*e3(k)
       inertia_an(k,k) = inertia_an(k,k) + inertia0(1,1)
    enddo
    err(5) = maxval(abs(inertia - inertia_an))/inertia0(1,1)
 endif
 if (id==master) print "(a,f4.2,a,i4,a,5es10.3)",' zstretch = ',zstretch,' nstep = ',nstep,&
                       ' errors (L,P,shape,vrigid,final): ',err
 deallocate(xyz0)

end subroutine run_rigid_drift

!-----------------------------------------------------------------------
!+
!  mean specific spin L = <r x v> and inertia tensor I = <r^2 1 - r r^T>
!  of the boundary particles about their centre of mass
!+
!-----------------------------------------------------------------------
subroutine get_free_core_moments(spin,inertia)
 use part,         only:npart,xyzh,vxyzu,iphase,iamtype,iamboundary
 use boundarypart, only:calculate_boundary_core_CM
 use vectorutils,  only:cross_product3D
 real, intent(out) :: spin(3),inertia(3,3)
 real    :: xcore(3),vcore(3),ri(3),li(3)
 integer :: i,k,nb

 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
 spin    = 0.
 inertia = 0.
 nb      = 0
 do i=1,npart
    if (.not.iamboundary(iamtype(iphase(i)))) cycle
    nb = nb + 1
    ri = xyzh(1:3,i) - xcore
    call cross_product3D(ri,vxyzu(1:3,i) - vcore,li)
    spin = spin + li
    do k=1,3
       inertia(:,k) = inertia(:,k) - ri(:)*ri(k)
       inertia(k,k) = inertia(k,k) + dot_product(ri,ri)
    enddo
 enddo
 spin    = spin/real(nb)
 inertia = inertia/real(nb)

end subroutine get_free_core_moments

!-----------------------------------------------------------------------
!+
!  Spin-up of a free core by a known applied force and torque
!
!  The input accelerations are a_i = a0 + b x r_i + S r_i, with r_i
!  relative to the current core centre of mass and S symmetric and
!  traceless. On the cubic lattice I = Ibar 1 exactly, so this gives a
!  constant net acceleration a0 and constant specific torque
!  tau = I b = Ibar b. The S r term is a non-rigid (tidal-like) stretch
!  with zero net force and torque, which the kick must discard.
!
!  Leapfrog with an exact rotational drift is exact for constant linear
!  and angular acceleration, so the following hold to round-off:
!   - L(t) = L0 + tau t,  omega(t) = omega0 + b t
!   - V(t) = V0 + a0 t,   X(t) = X0 + V0 t + a0 t^2/2
!   - the core stays rigid (pairwise distances, rigid velocities)
!   - for b parallel to omega0, every particle is at
!     x_i = X(t) + R_z(theta) r_i(0), with theta = omega0 t + b t^2/2
!  and for b perpendicular to omega0 (tilting the spin axis), all but
!  the last. Also checks that the kick is linear in its input, so that
!  applying it separately to two force contributions (as for fxyzu in
!  derivs and fext in substepping) gives the same result as applying
!  it to their sum.
!+
!-----------------------------------------------------------------------
subroutine test_core_known_torque(ntests,npass)
 use dim, only:mpi
 use io,  only:id,master
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(13)
 real    :: erra(nerr_torque),errb(nerr_torque)

 if (mpi) then
    if (id==master) write(*,"(/,a)") ' skipping known-torque test (requires no MPI)'
    return
 endif
 if (id==master) write(*,"(/,a)") '--> testing spin-up of a boundary-particle core by a known torque'

 call run_known_torque([0.,0.,0.02],erra)   ! spin up about the spin axis
 call run_known_torque([0.02,0.,0.],errb)   ! tilt the spin axis

 nfailed(:) = 0
 call checkval(erra(ierrt_angmom),0.,tol_angmom,nfailed(1),'spin-up: L = L0 + tau t')
 call checkval(erra(ierrt_omega),0.,tol_omega,nfailed(2),'spin-up: omega = omega0 + b t')
 call checkval(max(erra(ierrt_vcm),erra(ierrt_xcm)),0.,tol_vrigid,nfailed(3),'spin-up: centre of mass motion')
 call checkval(erra(ierrt_shape),0.,tol_xrigid,nfailed(4),'spin-up: pairwise distances preserved')
 call checkval(erra(ierrt_vrigid),0.,tol_vrigid,nfailed(5),'spin-up: velocities rigid')
 call checkval(erra(ierrt_xfinal),0.,tol_return,nfailed(6),'spin-up: positions match analytic rotation')
 call checkval(erra(ierrt_linear),0.,tol_angmom,nfailed(7),'kick is linear in its input')
 call checkval(errb(ierrt_angmom),0.,tol_angmom,nfailed(8),'tilt: L = L0 + tau t')
 call checkval(errb(ierrt_omega),0.,tol_omega,nfailed(9),'tilt: omega = omega0 + b t')
 call checkval(max(errb(ierrt_vcm),errb(ierrt_xcm)),0.,tol_vrigid,nfailed(10),'tilt: centre of mass motion')
 call checkval(errb(ierrt_shape),0.,tol_xrigid,nfailed(11),'tilt: pairwise distances preserved')
 call checkval(errb(ierrt_vrigid),0.,tol_vrigid,nfailed(12),'tilt: velocities rigid')
 call checkval(errb(ierrt_linear),0.,tol_angmom,nfailed(13),'kick is linear in its input (tilt)')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_core_known_torque

!-----------------------------------------------------------------------
!+
!  evolve a spinning, translating core of boundary particles for one
!  initial rotation period, over which the spin changes by |b| t ~ 50%,
!  under the applied accelerations
!  a_i = a0 + b x r_i + S r_i (see test_core_known_torque), and return
!  the errors listed by the ierrt_* indices
!+
!-----------------------------------------------------------------------
subroutine run_known_torque(bvec,err)
 use io,           only:id,master
 use part,         only:init_part,npart,npartoftype,massoftype,xyzh,vxyzu,iboundary,&
                        iphase,update_npartoftypetot
 use boundarypart, only:get_boundary_particle_forces,drift_boundary_core,&
                        calculate_boundary_core_CM,calculate_boundary_core_omega
 use vectorutils,  only:cross_product3D
 real, intent(in)  :: bvec(3)
 real, intent(out) :: err(nerr_torque)
 real, parameter :: pi = acos(-1.)   ! physcon's pi has 11 digits
 real, parameter :: omega0 = 0.5
 integer, parameter :: nstep = 200
 real :: vcen(3),a0(3),smat(3,3),omegavec0(3),xcm0(3),vcm0(3),xcm(3),vcm(3)
 real :: spin0(3),spin(3),inertia0(3,3),inertia(3,3),ibar,omega(3),omega_exp(3)
 real :: tend,dt,dvmax,dxmax,theta,ct,st,rot(3,3),xexp(3),dxret
 integer :: i,istep
 real, allocatable :: acc(:,:),acc1(:,:),acc2(:,:),xyz0(:,:)

 call init_part()
 npart = 0
 npartoftype(:) = 0
 vcen = [0.1,-0.2,0.3]
 call add_polytrope([0.,0.,0.],vcen,iboundary,psep_in=0.2)
 npartoftype(iboundary) = npart
 massoftype(iboundary)  = mstar/real(npart)
 call update_npartoftypetot()

 omegavec0 = [0.,0.,omega0]
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcm0,vcm0)
 do i=1,npart
    call cross_product3D(omegavec0,xyzh(1:3,i) - xcm0,vxyzu(1:3,i))
    vxyzu(1:3,i) = vxyzu(1:3,i) + vcen
 enddo
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcm0,vcm0)
 call get_free_core_moments(spin0,inertia0)
 ibar = (inertia0(1,1) + inertia0(2,2) + inertia0(3,3))/3.   ! isotropic on the cubic lattice

 a0   = [0.01,0.02,-0.03]
 smat = reshape([0.1, 0.05, 0., &
                 0.05,-0.05,0.02,&
                 0.,  0.02,-0.05],[3,3])
 tend = 2.*pi/omega0
 dt   = tend/real(nstep)

 allocate(acc(3,npart),acc1(3,npart),acc2(3,npart),xyz0(3,npart))
 xyz0 = xyzh(1:3,1:npart)
 !
 ! linearity of the kick: converting two contributions separately and
 ! summing gives the same as converting their sum
 !
 call get_applied_accel(a0,[0.,0.,0.],smat,acc1)
 call get_applied_accel([0.,0.,0.],bvec,0.*smat,acc2)
 acc = acc1 + acc2
 call get_boundary_particle_forces(npart,iphase,xyzh(1:3,:),acc1)
 call get_boundary_particle_forces(npart,iphase,xyzh(1:3,:),acc2)
 call get_boundary_particle_forces(npart,iphase,xyzh(1:3,:),acc)
 err(ierrt_linear) = maxval(abs(acc1 + acc2 - acc))/maxval(abs(acc))
 !
 ! kick-drift-kick, as in kickdrift
 !
 do istep=1,nstep
    do i=1,npart
       vxyzu(1:3,i) = vxyzu(1:3,i) + 0.5*dt*acc(:,i)
    enddo
    call drift_boundary_core(dt,npart,iphase,xyzh,vxyzu)
    call get_applied_accel(a0,bvec,smat,acc)
    call get_boundary_particle_forces(npart,iphase,xyzh(1:3,:),acc)
    do i=1,npart
       vxyzu(1:3,i) = vxyzu(1:3,i) + 0.5*dt*acc(:,i)
    enddo
 enddo
 !
 ! compare with the analytic solution
 !
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcm,vcm)
 call get_free_core_moments(spin,inertia)
 call calculate_boundary_core_omega(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),omega)
 call get_core_rigidity(xyz0,dvmax,dxmax)
 omega_exp = omegavec0 + bvec*tend
 err(ierrt_angmom) = norm2(spin - (spin0 + ibar*bvec*tend))/norm2(spin0)
 err(ierrt_omega)  = norm2(omega - omega_exp)/norm2(omega_exp)
 ! centre of mass errors relative to the final speed and the distance
 ! travelled. Each drift rebuilds v_i = V + omega x r_i and V is then
 ! re-averaged, so round-off in V grows by ~eps*omega*R per step and is
 ! checked with the velocity tolerance (measured up to ~3e-13)
 err(ierrt_vcm)    = norm2(vcm - (vcm0 + a0*tend))/norm2(vcm0 + a0*tend)
 err(ierrt_xcm)    = norm2(xcm - (xcm0 + vcm0*tend + 0.5*a0*tend**2))/norm2(vcm0*tend + 0.5*a0*tend**2)
 err(ierrt_shape)  = dxmax/rstar
 err(ierrt_vrigid) = dvmax/(norm2(omega_exp)*rstar)
 err(ierrt_xfinal) = 0.
 if (abs(bvec(1)) + abs(bvec(2)) <= 0.) then
    ! spin-up about z: rotation by theta = omega0 t + b t^2/2 about z
    theta = omega0*tend + 0.5*bvec(3)*tend**2
    ct = cos(theta)
    st = sin(theta)
    rot = reshape([ct,st,0., -st,ct,0., 0.,0.,1.],[3,3])
    dxret = 0.
    do i=1,npart
       xexp  = xcm + matmul(rot,xyz0(:,i) - xcm0)
       dxret = max(dxret,norm2(xyzh(1:3,i) - xexp))
    enddo
    err(ierrt_xfinal) = dxret/rstar
 endif
 if (id==master) print "(a,3f6.3,a,/,4x,8es10.3)",' b = (',bvec,&
    ') errors (L,omega,V,X,shape,vrigid,xfinal,linear):',err(1:8)
 deallocate(acc,acc1,acc2,xyz0)

end subroutine run_known_torque

!-----------------------------------------------------------------------
!+
!  applied accelerations a_i = a0 + b x r_i + S r_i on the boundary
!  particles, with r_i relative to the current core centre of mass
!+
!-----------------------------------------------------------------------
subroutine get_applied_accel(a0,bvec,smat,acc)
 use part,         only:npart,xyzh,vxyzu,iphase
 use boundarypart, only:calculate_boundary_core_CM
 use vectorutils,  only:cross_product3D
 real, intent(in)  :: a0(3),bvec(3),smat(3,3)
 real, intent(out) :: acc(:,:)
 real    :: xcm(3),vcm(3),ri(3),bxr(3)
 integer :: i

 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcm,vcm)
 do i=1,npart
    ri = xyzh(1:3,i) - xcm
    call cross_product3D(bvec,ri,bxr)
    acc(:,i) = a0 + bxr + matmul(smat,ri)
 enddo

end subroutine get_applied_accel

!-----------------------------------------------------------------------
!+
!  set up an isolated spinning boundary-particle core and evolve it for
!  two rotation periods (2*nper steps of dt). If use_sink, a sink of
!  negligible mass is placed far away so that step uses the sink-gas
!  substepping. Returns the errors listed by the ierr_* indices, and
!  optionally the final particle positions
!+
!-----------------------------------------------------------------------
subroutine run_free_spin(omega0,dt,nper,use_sink,dumps,err,xfinal)
 use io,             only:id,master,iverbose
 use part,           only:init_part,npart,npartoftype,massoftype,xyzh,vxyzu,fxyzu,fext,&
                          iboundary,iphase,iamtype,nptmass,xyzmh_ptmass,vxyz_ptmass,&
                          fxyz_ptmass,dsdt_ptmass,ihacc,ihsoft,update_npartoftypetot
 use eos,            only:ieos,gamma,polyk
 use options,        only:iexternalforce
 use units,          only:set_units
 use boundary,       only:set_boundary
 use dim,            only:periodic
 use physcon,        only:pi
 use kdtree,         only:tree_accuracy
 use ptmass,         only:get_accel_sink_gas,h_soft_sinkgas,set_integration_precision
 use deriv,          only:get_derivs_global
 use step_lf_global, only:init_step,step
 use timestep,       only:dtmax,tolv
 use boundarypart,   only:get_boundary_particle_forces,calculate_boundary_core_CM,&
                          calculate_boundary_core_omega
 use vectorutils,    only:cross_product3D
 use evwrite,        only:init_evfile
 use energies,       only:compute_energies,ev_data,iev_sum,iev_omegacore
 real,    intent(in)  :: omega0,dt
 integer, intent(in)  :: nper
 logical, intent(in)  :: use_sink,dumps
 real,    intent(out) :: err(nerr_spin)
 real, allocatable, intent(out), optional :: xfinal(:,:)
 integer :: i,istep
 real    :: t,dtext,dtnew,phii,dum,dum2,pmassi,lscale,dxret
 real    :: omegavec0(3),omega(3),xcore(3),vcore(3),mom0(3),angmom0(3),mom(3),angmom(3)
 real    :: dvmax,dxmax
 real, allocatable :: xyz0(:,:)

 call init_part()
 call set_units(mass=1.d0,dist=1.d0,G=1.d0)
 iverbose       = 0
 iexternalforce = 0
 tree_accuracy  = 0.
 tolv           = 1.e-2
 if (periodic) call set_boundary(-10.,10.,-10.,10.,-10.,10.)
 ieos  = 2
 gamma = 2.
 polyk = 2.*rstar**2/pi
 !
 ! isolated core of boundary particles centred on the origin, spinning
 ! rigidly about z with zero centre-of-mass velocity
 !
 npart = 0
 npartoftype(:) = 0
 call add_polytrope([0.,0.,0.],[0.,0.,0.],iboundary,psep_in=0.2)
 npartoftype(iboundary) = npart
 massoftype(iboundary)  = mstar/real(npart)
 call update_npartoftypetot()
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
 omegavec0 = [0.,0.,omega0]
 do i=1,npart
    call cross_product3D(omegavec0,xyzh(1:3,i) - xcore,vxyzu(1:3,i))
 enddo
 if (id==master) print "(a,i6,a,l1,a,f6.4)",' ncore = ',npart,' sink = ',use_sink,' dt = ',dt
 !
 ! optional distant sink of negligible mass: forces the substepping path
 !
 nptmass   = 0
 fext(:,:) = 0.
 if (use_sink) then
    nptmass = 1
    xyzmh_ptmass(:,1)      = 0.
    xyzmh_ptmass(1:3,1)    = [8.,0.,0.]
    xyzmh_ptmass(4,1)      = 1.e-10
    xyzmh_ptmass(ihacc,1)  = 0.01
    xyzmh_ptmass(ihsoft,1) = 0.
    h_soft_sinkgas         = 0.
    vxyz_ptmass(:,1)       = 0.
    fxyz_ptmass(:,1)       = 0.
    dsdt_ptmass(:,1)       = 0.
    do i=1,npart
       pmassi = massoftype(iamtype(iphase(i)))
       call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                               fext(1,i),fext(2,i),fext(3,i),phii,pmassi,fxyz_ptmass,dsdt_ptmass,dum,dum2)
    enddo
    call get_boundary_particle_forces(npart,iphase,xyzh(1:3,:),fext)
 endif
 !
 ! initial forces and conserved quantities
 !
 t     = 0.
 dtmax = dt
 fxyzu(:,:) = 0.
 call get_derivs_global()
 call get_momenta(mom0,angmom0)
 lscale = sqrt(dot_product(angmom0,angmom0))
 allocate(xyz0(3,npart))
 xyz0 = xyzh(1:3,1:npart)
 !
 ! evolve for two rotation periods
 !
 call set_integration_precision
 call init_step(npart,t,dtmax)
 if (dumps) call write_test_dump('boundary_spin_',t,0)
 err = 0.
 do istep=1,2*nper
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    if (dumps .and. mod(istep,ndumpfreq)==0) call write_test_dump('boundary_spin_',t,istep/ndumpfreq)
    if (istep==nper .or. istep==2*nper) then
       call get_core_rigidity(xyz0,dvmax,dxmax)
       call calculate_boundary_core_omega(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),omega)
       err(ierr_vrigid) = max(err(ierr_vrigid),dvmax/(omega0*rstar))
       err(ierr_omega)  = max(err(ierr_omega),sqrt(dot_product(omega - omegavec0,omega - omegavec0))/omega0)
       if (istep==nper) then
          err(ierr_shape1) = dxmax/rstar
          dxret = 0.
          do i=1,npart
             dxret = max(dxret,sqrt(dot_product(xyzh(1:3,i) - xyz0(:,i),xyzh(1:3,i) - xyz0(:,i))))
          enddo
          err(ierr_return) = dxret/rstar
       else
          err(ierr_shape2) = dxmax/rstar
       endif
    endif
 enddo
 call get_momenta(mom,angmom)
 err(ierr_angmom) = sqrt(dot_product(angmom - angmom0,angmom - angmom0))/lscale
 err(ierr_mom)    = sqrt(dot_product(mom,mom))/(mstar*omega0*rstar)
 !
 ! core angular velocity written to the .ev file (no file is opened)
 !
 call init_evfile(1,'boundary_spin.ev',.false.)
 call compute_energies(t)
 err(ierr_evomega) = max(sqrt(dot_product(ev_data(iev_sum,iev_omegacore(1:3)) - omegavec0,&
                                          ev_data(iev_sum,iev_omegacore(1:3)) - omegavec0)),&
                         abs(ev_data(iev_sum,iev_omegacore(4)) - omega0))/omega0
 if (id==master) print "(a,8es10.3)",' errors (L,P,shape1,shape2,vrigid,omega,return,ev omega): ',err
 if (present(xfinal)) then
    ! positions relative to the core centre of mass, so that the (real)
    ! pull of the sink on the core as a whole does not count as a difference
    call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
    allocate(xfinal(3,npart))
    do i=1,npart
       xfinal(:,i) = xyzh(1:3,i) - xcore
    enddo
 endif
 deallocate(xyz0)

 ! reset
 nptmass = 0
 iverbose = 0

end subroutine run_free_spin

!-----------------------------------------------------------------------
!+
!  measure how far the boundary particles are from rigid-body motion:
!   dvmax = max |v_i - v_cm - omega x (x_i - x_cm)|, with omega = I^-1 L
!   dxmax = max | |x_i - x_j| - |x0_i - x0_j| | over all pairs of
!           boundary particles, where x0 are the initial positions.
!  Pairwise distances are preserved by any rigid motion, so dxmax needs
!  no fit for the rotation of the core
!+
!-----------------------------------------------------------------------
subroutine get_core_rigidity(xyz0,dvmax,dxmax)
 use io,           only:id,master
 use part,         only:npart,xyzh,vxyzu,iphase,iamtype,iamboundary,isdead_or_accreted
 use boundarypart, only:calculate_boundary_core_CM,calculate_boundary_core_omega
 use vectorutils,  only:cross_product3D
 real, intent(in)  :: xyz0(:,:)
 real, intent(out) :: dvmax,dxmax
 integer, allocatable :: ilist(:)
 integer :: i,j,ii,jj,ncore
 real    :: xcore(3),vcore(3),omega(3),omega_x_r(3),dv(3),dij,dij0

 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
 call calculate_boundary_core_omega(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),omega)
 if (id==master) print "(a,es10.3)",' core |omega| = ',sqrt(dot_product(omega,omega))

 allocate(ilist(npart))
 ncore = 0
 dvmax = 0.
 do i=1,npart
    if (iamboundary(iamtype(iphase(i))) .and. .not.isdead_or_accreted(xyzh(4,i))) then
       ncore = ncore + 1
       ilist(ncore) = i
       call cross_product3D(omega,xyzh(1:3,i) - xcore,omega_x_r)
       dv = vxyzu(1:3,i) - vcore - omega_x_r
       dvmax = max(dvmax,sqrt(dot_product(dv,dv)))
    endif
 enddo

 dxmax = 0.
 !$omp parallel do default(none) &
 !$omp shared(ncore,ilist,xyzh,xyz0) &
 !$omp private(ii,jj,i,j,dij,dij0) &
 !$omp reduction(max:dxmax) &
 !$omp schedule(dynamic,16)
 do ii=1,ncore
    i = ilist(ii)
    do jj=ii+1,ncore
       j = ilist(jj)
       dij  = sqrt(dot_product(xyzh(1:3,i) - xyzh(1:3,j),xyzh(1:3,i) - xyzh(1:3,j)))
       dij0 = sqrt(dot_product(xyz0(:,i) - xyz0(:,j),xyz0(:,i) - xyz0(:,j)))
       dxmax = max(dxmax,abs(dij - dij0))
    enddo
 enddo
 !$omp end parallel do
 deallocate(ilist)

end subroutine get_core_rigidity

!-----------------------------------------------------------------------
!+
!  append an n=1 polytrope (M = R = G = 1) centred on xcen, moving with
!  velocity vcen, made entirely of particles of type itype.
!  Particle masses must be set by the caller.
!+
!-----------------------------------------------------------------------
subroutine add_polytrope(xcen,vcen,itype,psep_in)
 use io,         only:id,master
 use part,       only:npart,xyzh,vxyzu,hfact,set_particle_type
 use eos,        only:gamma,polyk
 use spherical,  only:set_sphere
 use stretchmap, only:rho_func
 real,    intent(in) :: xcen(3),vcen(3)
 integer, intent(in) :: itype
 real,    intent(in), optional :: psep_in
 procedure(rho_func), pointer :: density_func
 integer :: i,npin
 integer(kind=8) :: nptot
 real    :: psep,r

 psep = 0.15
 if (present(psep_in)) psep = psep_in
 npin = npart
 density_func => rho_n1_polytrope
 nptot = npart
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,nptot,&
                 rhofunc=density_func,xyz_origin=xcen)
 do i=npin+1,npart
    call set_particle_type(i,itype)
    r = sqrt(dot_product(xyzh(1:3,i) - xcen,xyzh(1:3,i) - xcen))
    vxyzu(1:3,i) = vcen
    if (size(vxyzu,1) >= 4) vxyzu(4,i) = polyk*rho_n1_polytrope(r)/(gamma - 1.)
 enddo

end subroutine add_polytrope

!-----------------------------------------------------------------------
!+
!  centre of mass position and velocity of the boundary particles and
!  of the gas particles (all particles of a given type have equal mass)
!+
!-----------------------------------------------------------------------
subroutine get_star_CMs(xbound,vbound,xgas,vgas)
 use part, only:npart,xyzh,vxyzu,iphase,iamtype,igas,iboundary,isdead_or_accreted
 real, intent(out) :: xbound(3),vbound(3),xgas(3),vgas(3)
 integer :: i,itype,nb,ng

 nb = 0
 ng = 0
 xbound = 0.
 vbound = 0.
 xgas   = 0.
 vgas   = 0.
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    itype = iamtype(iphase(i))
    if (itype==iboundary) then
       nb = nb + 1
       xbound = xbound + xyzh(1:3,i)
       vbound = vbound + vxyzu(1:3,i)
    elseif (itype==igas) then
       ng = ng + 1
       xgas = xgas + xyzh(1:3,i)
       vgas = vgas + vxyzu(1:3,i)
    endif
 enddo
 if (nb > 0) then
    xbound = xbound/real(nb)
    vbound = vbound/real(nb)
 endif
 if (ng > 0) then
    xgas = xgas/real(ng)
    vgas = vgas/real(ng)
 endif

end subroutine get_star_CMs

!-----------------------------------------------------------------------
!+
!  write a full dump for debugging
!+
!-----------------------------------------------------------------------
subroutine write_test_dump(prefix,t,idump)
 use io,              only:id,master
 use readwrite_dumps, only:write_fulldump
 character(len=*), intent(in) :: prefix
 real,             intent(in) :: t
 integer,          intent(in) :: idump
 character(len=40)            :: filename

 write(filename,"(a,i5.5)") trim(prefix),idump
 if (id==master) print "(a,es10.3)",' writing '//trim(filename)//' at t = ',t
 call write_fulldump(t,filename)

end subroutine write_test_dump

!-----------------------------------------------------------------------
!+
!  total linear and angular momentum of particles and sinks
!+
!-----------------------------------------------------------------------
subroutine get_momenta(mom,angmom)
 use part,        only:npart,xyzh,vxyzu,massoftype,iphase,iamtype,isdead_or_accreted,&
                       nptmass,xyzmh_ptmass,vxyz_ptmass,ispinx,ispinz
 use vectorutils, only:cross_product3D
 real, intent(out) :: mom(3),angmom(3)
 real    :: pmassi,li(3)
 integer :: i

 mom    = 0.
 angmom = 0.
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    pmassi = massoftype(iamtype(iphase(i)))
    mom = mom + pmassi*vxyzu(1:3,i)
    call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),li)
    angmom = angmom + pmassi*li
 enddo
 do i=1,nptmass
    pmassi = xyzmh_ptmass(4,i)
    if (pmassi < 0.) cycle
    mom = mom + pmassi*vxyz_ptmass(1:3,i)
    call cross_product3D(xyzmh_ptmass(1:3,i),vxyz_ptmass(1:3,i),li)
    angmom = angmom + pmassi*li + xyzmh_ptmass(ispinx:ispinz,i)
 enddo

end subroutine get_momenta

!-----------------------------------------------------------------------
!+
!  density profile of an n=1 polytrope with M = R = G = 1
!+
!-----------------------------------------------------------------------
real function rho_n1_polytrope(r)
 use physcon, only:pi
 real, intent(in) :: r
 real :: x

 x = pi*r/rstar
 if (x < 1.e-6) then
    rho_n1_polytrope = pi*mstar/(4.*rstar**3)
 else
    rho_n1_polytrope = pi*mstar/(4.*rstar**3)*sin(x)/x
 endif

end function rho_n1_polytrope

end module testboundary
