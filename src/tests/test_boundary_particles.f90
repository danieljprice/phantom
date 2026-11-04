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
!   step_lf_global, stretchmap, testutils, timestep, units
!
 use testutils, only:checkval,update_test_scores
 implicit none
 public :: test_boundary

 private

 ! n=1 polytrope with M = R = G = 1
 real, parameter :: rstar = 1.
 real, parameter :: mstar = 1.

 logical :: write_dumps = .true.
 integer, parameter :: ndumpfreq = 10

contains

subroutine test_boundary(ntests,npass)
 use io, only:id,master
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING BOUNDARY PARTICLES'

 call test_core_star_sink_orbit(ntests,npass)
 call test_boundary_star_gas_star_orbit(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- BOUNDARY PARTICLE TEST COMPLETE'

end subroutine test_boundary

!-----------------------------------------------------------------------
!+
!  n=1 polytrope whose inner half (r < R/2) is a rigid core of boundary
!  particles, on a circular Keplerian orbit around a sink particle.
!
!  Checks that
!   - the core stays rigid (same velocity, fixed shape)
!   - total linear momentum is conserved, which requires the core to
!     feel the reaction to the forces it exerts on the envelope
!   - total angular momentum and energy are approximately conserved
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
 integer, intent(inout) :: ntests,npass
 integer :: i,istep,nsteps,nfailed(6),ncore,nalive
 integer(kind=8) :: nptot
 real    :: psep,rcore,msink,a,omega,xstar(3),vstar(3),xsink(3),vsink(3)
 real    :: t,dt,dtext,dtnew,tmax,phii,dum,dum2,pmassi
 real    :: mom0(3),angmom0(3),mom(3),angmom(3),etot0,pscale,lscale
 real    :: xcore(3),vcore(3),dvmax,dxmax,sep0,sep
 real, allocatable :: dxcore0(:,:)
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
 pscale = mstar*sqrt(dot_product(vstar,vstar))   ! orbital momentum of the star
 lscale = sqrt(dot_product(angmom0,angmom0))
 call calculate_boundary_core_CM(npart,iphase,xyzh(1:3,:),vxyzu(1:3,:),xcore,vcore)
 sep0 = sqrt(dot_product(xcore - xyzmh_ptmass(1:3,1),xcore - xyzmh_ptmass(1:3,1)))
 allocate(dxcore0(3,npart))
 do i=1,npart
    dxcore0(:,i) = xyzh(1:3,i) - xcore
 enddo
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
 dvmax = 0.
 dxmax = 0.
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) nalive = nalive + 1
    if (iamboundary(iamtype(iphase(i))) .and. .not.isdead_or_accreted(xyzh(4,i))) then
       dvmax = max(dvmax,sqrt(dot_product(vxyzu(1:3,i) - vcore,vxyzu(1:3,i) - vcore)))
       dxmax = max(dxmax,sqrt(dot_product(xyzh(1:3,i) - xcore - dxcore0(:,i),&
                                          xyzh(1:3,i) - xcore - dxcore0(:,i))))
    endif
 enddo
 deallocate(dxcore0)
 sep = sqrt(dot_product(xcore - xyzmh_ptmass(1:3,1),xcore - xyzmh_ptmass(1:3,1)))
 !
 ! conservation
 !
 call get_momenta(mom,angmom)
 call compute_energies(t)

 nfailed(:) = 0
 call checkval(nalive,npart,0,nfailed(1),'no particles accreted')
 call checkval(dvmax/sqrt(dot_product(vcore,vcore)),0.,1.e-12,nfailed(2),'core velocities equal')
 call checkval(dxmax/rcore,0.,1.e-12,nfailed(3),'core shape preserved')
 call checkval(sqrt(dot_product(mom - mom0,mom - mom0))/pscale,0.,1.e-12,nfailed(4),'linear momentum')
 call checkval(sqrt(dot_product(angmom - angmom0,angmom - angmom0))/lscale,0.,1.e-4,nfailed(5),'angular momentum')
 call checkval(etot,etot0,1.e-2,nfailed(6),'total energy')
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
 real, allocatable :: dxbound0(:,:)

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
 allocate(dxbound0(3,npart))
 do i=1,npart
    dxbound0(:,i) = xyzh(1:3,i) - xbound
 enddo
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
 dvmax  = 0.
 dxmax  = 0.
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    nalive = nalive + 1
    if (iamboundary(iamtype(iphase(i)))) then
       dvmax = max(dvmax,sqrt(dot_product(vxyzu(1:3,i) - vbound,vxyzu(1:3,i) - vbound)))
       dxmax = max(dxmax,sqrt(dot_product(xyzh(1:3,i) - xbound - dxbound0(:,i),&
                                          xyzh(1:3,i) - xbound - dxbound0(:,i))))
    endif
 enddo
 deallocate(dxbound0)
 !
 ! orbit and conservation
 !
 call get_star_CMs(xbound,vbound,xgas,vgas)
 sep = sqrt(dot_product(xbound - xgas,xbound - xgas))
 call get_momenta(mom,angmom)
 call compute_energies(t)

 nfailed(:) = 0
 call checkval(nalive,npart,0,nfailed(1),'no particles lost')
 call checkval(dvmax/sqrt(dot_product(vbound,vbound)),0.,1.e-12,nfailed(2),'boundary star velocities equal')
 call checkval(dxmax/rstar,0.,1.e-12,nfailed(3),'boundary star shape preserved')
 call checkval(sqrt(dot_product(mom - mom0,mom - mom0))/pscale,0.,1.e-12,nfailed(4),'linear momentum')
 call checkval(sqrt(dot_product(angmom - angmom0,angmom - angmom0))/lscale,0.,1.e-4,nfailed(5),'angular momentum')
 call checkval(etot,etot0,1.e-2,nfailed(6),'total energy')
 call update_test_scores(ntests,nfailed,npass)

 nfailed(:) = 0
 call checkval(sep,sep0,1.e-3,nfailed(1),'boundary star-gas star separation (circular orbit)')
 call update_test_scores(ntests,nfailed(1:1),npass)

 iverbose = 0

end subroutine test_boundary_star_gas_star_orbit

!-----------------------------------------------------------------------
!+
!  append an n=1 polytrope (M = R = G = 1) centred on xcen, moving with
!  velocity vcen, made entirely of particles of type itype.
!  Particle masses must be set by the caller.
!+
!-----------------------------------------------------------------------
subroutine add_polytrope(xcen,vcen,itype)
 use io,         only:id,master
 use part,       only:npart,xyzh,vxyzu,hfact,set_particle_type
 use eos,        only:gamma,polyk
 use spherical,  only:set_sphere
 use stretchmap, only:rho_func
 real,    intent(in) :: xcen(3),vcen(3)
 integer, intent(in) :: itype
 procedure(rho_func), pointer :: density_func
 integer :: i,npin
 integer(kind=8) :: nptot
 real    :: psep,r

 psep = 0.15
 npin = npart
 density_func => rho_n1_polytrope
 nptot = npart
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,nptot,&
                 rhofunc=density_func,xyz_origin=xcen)
 do i=npin+1,npart
    call set_particle_type(i,itype)
    r = sqrt(dot_product(xyzh(1:3,i) - xcen,xyzh(1:3,i) - xcen))
    vxyzu(1:3,i) = vcen
    vxyzu(4,i)   = polyk*rho_n1_polytrope(r)/(gamma - 1.)
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
