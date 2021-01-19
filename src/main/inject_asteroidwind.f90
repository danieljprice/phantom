!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! None
!
! :References: None
!
! :Owner: Bec Nealon
!
! :Runtime parameters:
!   - dndt_type     : *injection rate (0=const, 1=cos(t), 2=r^(-2))*
!   - mdot          : *mass injection rate in grams/second*
!   - npartperorbit : *particle injection rate in particles/binary orbit*
!   - vlag          : *percentage lag in velocity of wind*
!
! :Dependencies: binaryutils, externalforces, infile_utils, io, options,
!   part, partinject, physcon, random, units
!
 use io, only:error
 use physcon, only:pi
 implicit none
 character(len=*), parameter, public :: inject_type = 'asteroidwind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 private

 real         :: mdot          = 5.e8     ! mass injection rate in grams/second
 real         :: npartperorbit = 100.     ! particle injection rate in particles per orbit
 real         :: vlag          = 0.1      ! percentage lag in velocity of wind
 integer      :: dndt_type     = 0        ! injection rate (0=const, 1=cos(t), 2=r^(-2))
 real,save    :: dndt_scaling             ! scaling to get ninject correct
 logical,save :: scaling_set              ! has the scaling been set (initially false)

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(inout) :: ierr

 scaling_set = .false.

 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Inject particles
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,            only:fatal
 use part,          only:nptmass,massoftype,igas,hfact,ihsoft
 use partinject,    only:add_or_update_particle
 use physcon,       only:twopi,gg,kboltz,mass_proton_cgs
 use random,        only:ran2
 use units,         only:udist, umass, utime
 use options,       only:iexternalforce
 use externalforces,only:mass1
 use binaryutils,   only:get_orbit_bits
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real,    dimension(3)  :: xyz,vxyz,r1,r2,v2,vhat,v1
 integer :: i,ipart,npinject,seed,pt,test_integral
 real    :: dmdt,dndt,rasteroid,h,u,speed,inject_this_step
 real    :: m1,m2,mu,period,r,q
 real    :: phi,theta,mod_time,dt,func_now, phi_facing

 real, save :: have_injected,func_old,t_old
 real, save :: semia, ra, rp, ecc

 if (nptmass < 2 .and. iexternalforce == 0) call fatal('inject_asteroidwind','not enough point masses for asteroid wind injection')
 if (nptmass > 2) call fatal('inject_asteroidwind','too many point masses for asteroid wind injection')

 if (nptmass == 2) then
    pt = 2
    r1 = xyzmh_ptmass(1:3,1)
    m1 = xyzmh_ptmass(4,1)
    v1 = vxyz_ptmass(1:3,1)
 else
    pt = 1
    r1 = 0.
    m1 = mass1
    v1 = 0.
 endif

 r2        = xyzmh_ptmass(1:3,pt)
 rasteroid = xyzmh_ptmass(ihsoft,pt)
 m2        = xyzmh_ptmass(4,pt)
 v2        = vxyz_ptmass(1:3,pt)

 speed     = sqrt(dot_product(v2,v2))
 vhat      = v2/speed

 r         = sqrt(dot_product(r1-r2,r1-r2))
 q         = m2/m1
 mu        = 1./(1 + q)

 ! Calculate semi major axis and eccentricity from energy and radius
 if (.not.scaling_set) call get_orbit_bits(v2-v1,r2-r1,m1,iexternalforce,semia,ecc,ra,rp)

 period    = twopi*sqrt((semia*udist)**3/(gg*(m1+m2)*umass))   ! period of orbit
 period    = period/utime                                      ! in code units

 dmdt      = mdot/(umass/utime)                                ! convert grams/sec to code units
 dndt      = npartperorbit/period                              ! convert particles per orbit into code units
 mod_time  = mod(time,period)/period

!
!-- Mass of gas particles is set by mass accretion rate and particle injection rate
!
 massoftype(igas) = dmdt/dndt

! If it hasn't yet been calculated, find the scaling for dndt to give the correct ninject
 if (.not.scaling_set) then
    ! First guess (best to be big)
    dndt_scaling = 1000.

    ! Integrate across one orbit
    if (dndt_type < 2) then
       call integrate_it(0.,1.,ra,rp,ecc,test_integral)
    else
       call integrate_it_with_r(0.,1.,ra,rp,ecc,semia,test_integral)
    endif

    ! Calculate scaling for the rest of the simulation
    dndt_scaling = real(npartperorbit/(test_integral/dndt_scaling))

    t_old = mod(time*0.99,period)/period
    have_injected = npartoftype(igas)
    if (time < tiny(time)) have_injected = 0.
    func_old = dndt_func(t_old,r,ra,rp,ecc)

    ! Save these values for future
    scaling_set = .true.
 endif

!
!-- How many particles do we need to inject?
!   (Seems to need at least eight gas particles to not crash) <-- This statement may or may not be true...
!
 if (npartoftype(igas)<8) then
    npinject = 8-npartoftype(igas)
 else
    ! Calculate how many extra particles from previous step to now
    dt = mod_time - t_old
    if (dt < 0.) dt = dt + 1.0 !if it's just ticked over to the next orbit

    ! Just trapezoidal rule between the previous step and this one
    func_now = dndt_func(mod_time,r,ra,rp,ecc)
    inject_this_step = 0.5*dt*(func_old + func_now)

    npinject = max(0, int(0.5 + have_injected + inject_this_step - npartoftype(igas) ))

    ! Save for next step (faster than integrating the whole thing each time)
    t_old = mod_time
    have_injected = have_injected + inject_this_step
    func_old = func_now
 endif

!
!-- Randomly inject particles around the asteroids outer 'radius'
!-- Only inject them on the side that is facing the central sink
!
 do i=1,npinject
    phi_facing = atan2(xyzmh_ptmass(2,pt),xyzmh_ptmass(1,pt)) + 0.5*pi
    phi       = phi_facing + ran2(seed)*pi
    theta     = ran2(seed)*pi
    xyz       = r2 + (/rasteroid*cos(phi)*sin(theta),rasteroid*sin(phi)*sin(theta),rasteroid*cos(theta)/)
    vxyz      = (1.-vlag/100)*speed*vhat
    u         = 0. ! setup is isothermal so utherm is not stored
    h         = hfact*(rasteroid/2.)
    ipart     = npart + 1
    call add_or_update_particle(igas,xyz,vxyz,h,u,ipart,npart,npartoftype,xyzh,vxyzu)
 enddo

 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Returns dndt(t) depending on which function is chosen
!  Note that time in this function is strictly the fraction
!  of the orbit, not absolute time
!+
!-----------------------------------------------------------------------

function dndt_func(t,r,ra,rp,ecc) result(dndt)
 real, intent(in) :: t,ra,rp,ecc,r

 real :: dndt

 select case (dndt_type)
 case (0)
    dndt = 1.0
 case (1)
    dndt = 0.5*(cos(2.*pi*(t-0.5)) + 1.)
 case (2)
    dndt = (ra*rp/(r**2)) - ((1.-ecc)/(1.+ecc))
 case default
    call error('inject_asteroid','dndt choice does not exist, setting to zero')
    dndt = 0.
 end select

 dndt = dndt*dndt_scaling

end function dndt_func

!-----------------------------------------------------------------------
!+
!  Cheap dirty integration using trapezoidal rule
!  Takes into account that the integral must be int
!+
!-----------------------------------------------------------------------

subroutine integrate_it(t_start,t_end,ra,rp,ecc,integral)
 real, intent(in)    :: t_start, t_end,ra,rp,ecc
 integer, intent(out)   :: integral
 integer :: ii, nint
 real :: ya,yb,t_int,dt,re_integral

 re_integral = 0.
 nint = 2*npartperorbit
 dt = real((t_end - t_start)/nint)
 t_int = t_start

 do ii = 1,nint
    ya = dndt_func(t_int,1.,ra,rp,ecc)
    yb = dndt_func(t_int+dt,1.,ra,rp,ecc)

    re_integral = re_integral + dt*0.5*(ya + yb)

    t_int = t_int + dt
 enddo

 integral = int(0.5 + re_integral)

end subroutine integrate_it

!-----------------------------------------------------------------------
!+
!  Cheap dirty integration using trapezoidal rule
!  This is the more complicated version because r is not a
!  simple function of t, so needs to be treated differently
!+
!-----------------------------------------------------------------------

subroutine integrate_it_with_r(t_start,t_end,ra,rp,ecc,semia,integral)
 use binaryutils, only:get_E
 real, intent(in)    :: t_start, t_end,ra,rp,ecc,semia
 integer, intent(out)   :: integral
 integer :: ii, nint
 real :: ya,yb,t_int,dt,re_integral
 real :: r_new,theta,E,denom,numer

 re_integral = 0.
 nint = 2*npartperorbit
 dt = real((t_end - t_start)/nint)
 t_int = t_start

 ! set up for first step
 call get_E(1.0,ecc,0.0,E)
 theta = atan2(sqrt(1.-ecc**2)*sin(E),(cos(E) - ecc))
 r_new = semia*(1. - ecc**2)/(1. + ecc*cos(theta))
 ya = dndt_func(0.,r_new,ra,rp,ecc)

 do ii = 1,nint

    call get_E(1.0,ecc,t_int,E)

    numer = sqrt(1. - ecc**2)*sin(E)
    denom = cos(E) - ecc
    theta = atan2(numer,denom)

    r_new = semia*(1. - ecc**2)/(1. + ecc*cos(theta))

    yb = dndt_func(t_int+dt,r_new,ra,rp,ecc)

    re_integral = re_integral + dt*0.5*(ya + yb)

    t_int = t_int + dt
    ya = yb
 enddo

 integral = int(0.5 + re_integral)

end subroutine integrate_it_with_r

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(mdot         ,'mdot'         ,'mass injection rate in grams/second'              ,iunit)
 call write_inopt(npartperorbit,'npartperorbit','particle injection rate in particles/binary orbit',iunit)
 call write_inopt(vlag         ,'vlag'         ,'percentage lag in velocity of wind'               ,iunit)
 call write_inopt(dndt_type    ,'dndt_type'    ,'injection rate (0=const, 1=cos(t), 2=r^(-2))'     ,iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 select case(trim(name))
 case('mdot')
    read(valstring,*,iostat=ierr) mdot
    ngot = ngot + 1
    if (mdot  <  0.) call fatal(label,'mdot < 0 in input options')
 case('npartperorbit')
    read(valstring,*,iostat=ierr) npartperorbit
    ngot = ngot + 1
    if (npartperorbit < 0.) call fatal(label,'npartperorbit < 0 in input options')
 case('vlag')
    read(valstring,*,iostat=ierr) vlag
    ngot = ngot + 1
 case('dndt_type')
    read(valstring,*,iostat=ierr) dndt_type
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
