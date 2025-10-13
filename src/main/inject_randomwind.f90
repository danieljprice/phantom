!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection module for wind from an orbiting body, as used
! in Trevascus et al. (2021)
!
! :References:
!   Trevascus et al. (2021), MNRAS 505, L21-L25
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - delta_theta : *standard deviation for the gaussion distribution*
!   - inject_pt   : *the particle that excites wind (when wind_type=1)*
!   - mdot        : *mass injection rate with unit, e.g. 1e8*g/s, 1e-7M_s/yr*
!   - mdot_type   : *injection rate (0=const, 2=r^(-2))*
!   - phi         : *the tilt orientation of the star, (random_type=1)*
!   - r_inject    : *inject radius with units, e.g. 1*AU, 1e8m, (when wind_type=1)*
!   - r_ref       : *radius at whieh Mdot=mdot for 1/r^2 injection type*
!   - random_type : *random position on the surface, 0 for random, 1 for gaussian*
!   - theta       : *the tilt inclination of the star or planet (random_type=1)*
!   - vlag        : *percentage lag in velocity of wind*
!   - wind_type   : *wind setup (0=asteroidwind, 1=randomwind, 2=boil-off)*
!
! :Dependencies: evolveplanet, externalforces, infile_utils, io, options,
!   part, partinject, physcon, random, units, vectorutils
!
 use io, only:error
 use physcon, only:pi
 use random,        only:get_random_pos_on_sphere, get_gaussian_pos_on_sphere, ran2
 implicit none
 character(len=*), parameter, public :: inject_type = 'randomwind'
 character(len=20), public :: mdot_str = "5.e8*g/s"
 character(len=20), public :: r_inject_str = "0.5*au"
 real, public :: mdot = 1.e8     ! mass injection rate in grams/second

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par

 real, public :: have_injected = 0.
 real, public :: inject_this_step = 0.

 private

 real         :: npartperorbit = 1000.     ! particle injection rate in particles per orbit
 integer      :: wind_type     = 0        ! wind setup (0=asteroidwind, 1=randomwind)
 real         :: vlag          = 0.0      ! percentage lag in velocity of wind
 integer      :: mdot_type     = 0        ! injection rate (0=const, 1=cos(t), 2=r^(-2))
 integer      :: random_type   = 0        ! random position on the surface, 0 for random, 1 for gaussian
 real         :: delta_theta   = 0.5      ! standard deviation for the gaussion distribution (random_type=1)
 real         :: t_old         = 0.
 real         :: r_ref         = 1.       ! reference radius for mdot_type=2
 real         :: theta         = 0.       ! the inclination of the star or planet
 real         :: phi           = 0.       ! the orientation of the star
 integer      :: inject_pt     = 2        ! the partical that produce wind (when wind_type=1)
 real         :: wind_speed    = 1.0      ! wind speed in code unit (when wind_type=1)
 real         :: wind_speed_factor = 1.2  ! factor to scale the wind speed based on the Keplerian speed at rinject
 !real         :: rinject       = 1.0

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use options, only:need_pressure_on_sinks
 integer, intent(inout) :: ierr

 ierr = 0
 if (wind_type==2) need_pressure_on_sinks= .true.

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Inject particles
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,            only:fatal
 use part,          only:nptmass,massoftype,igas,hfact,ihsoft,ipbondi,irbondi
 use partinject,    only:add_or_update_particle
 use physcon,       only:twopi,gg,kboltz,mass_proton_cgs
 use random,        only:get_random_pos_on_sphere, get_gaussian_pos_on_sphere
 use units,         only:in_code_units
 use vectorutils,   only:cross_product3D, rotatevec
 use options,       only:iexternalforce
 use externalforces,only:mass1
 use evolveplanet,  only:evolve_planet
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: ierr
 real,    dimension(3)  :: xyz,vxyz,r1,r2,v2,vhat,v1
 integer :: i,ipart,npinject,seed,pt
 real    :: dmdt,rinject,h,u,speed,m1,m2,r,dt
 real    :: dx(3), vecz(3), veczprime(3), rotaxis(3), cs
 real    :: theta_rad,phi_rad,cost,sint,cosp,sinp,mdotacc,mdotwind
 real    :: Minject, frac_extra

 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)
 if (dtlast <= 0) return

 ! initialise some parameter to avoid warning...
 pt = 1
 rinject = 1.0
 r = 1.0

 ! calculate the wind velocity and other quantities for different wind type
 select case (wind_type)
 case(1,2) ! set up random wind
    if (inject_pt > nptmass) call fatal('inject_randomwind', 'not enough point masses for inject target, check inject_pt')
    if (wind_type==2) then
       rinject = xyzmh_ptmass(irbondi,inject_pt)
       cs = sqrt(xyzmh_ptmass(4,inject_pt)/rinject)
       wind_speed = cs
       call evolve_planet(xyzmh_ptmass(ipbondi,inject_pt),xyzmh_ptmass(irbondi,inject_pt),mdotacc,mdotwind)
       if (rinject<=xyzmh_ptmass(ihsoft,inject_pt)) call fatal('inject_randomwind', 'Bondi radius not set for wind_type=2')
    else
       rinject = in_code_units(r_inject_str, ierr)
    endif
    r2 = xyzmh_ptmass(1:3,inject_pt)
    v2        = vxyz_ptmass(1:3,inject_pt)
    wind_speed = wind_speed_factor*sqrt(xyzmh_ptmass(4, inject_pt)/rinject)
    u         = 0. ! setup is isothermal so utherm is not stored
    h         = hfact
 case default ! set up asteroid wind
    if (nptmass < 1 .and. iexternalforce == 0) &
      call fatal('inject_asteroidwind','not enough point masses for asteroid wind injection')
    if (nptmass > 2) &
      call fatal('inject_asteroidwind','too many point masses for asteroid wind injection')
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
    rinject   = xyzmh_ptmass(ihsoft,pt)
    m2        = xyzmh_ptmass(4,pt)
    v2        = vxyz_ptmass(1:3,pt)
    speed     = sqrt(dot_product(v2,v2))
    vhat      = v2/speed
    r         = sqrt(dot_product(r1-r2,r1-r2))
    wind_speed  = (1.-vlag/100)*speed
    u         = 0. ! setup is isothermal so utherm is not stored
    h         = hfact*(rinject/2.)
 end select

 !
 ! Add any dependency on radius to mass injection rate (and convert to code units)
 !
 if (wind_type==2) then
    dmdt = mdotwind
 else
    mdot = in_code_units(mdot_str,ierr)
    dmdt = mdot*mdot_func(r,r_ref) ! r_ref is the radius for which mdot_fund = mdot
 endif
 !
 !-- How many particles do we need to inject?
 !   (Seems to need at least eight gas particles to not crash) <-- This statement may or may not be true...
 !
 if (npartoftype(igas) < 8) then
    npinject = 8-npartoftype(igas)
 else
    ! Calculate how many particles to inject based on the mass injection rate
    dt = time - t_old
    Minject = dtlast*dmdt
    npinject = int(Minject/massoftype(igas))
    ! Add one particle with probability equal to the fractional part
    frac_extra = Minject/massoftype(igas) - npinject
    if (ran2(seed) < frac_extra) npinject = npinject + 1
    ! Save time for next step
    t_old = time
 endif
 !
 !-- set up the tilt of the star, and vectors for rotation
 !
 theta_rad = theta*pi/180.
 phi_rad = phi*pi/180.
 cost = cos(theta_rad)
 sint = sin(theta_rad)
 cosp = cos(phi_rad)
 sinp = sin(phi_rad)
 vecz = (/0.,0.,1./)
 veczprime = (/sint*cosp,sint*sinp,cost/)
 if (abs(theta_rad-0)<1e-6) then
    rotaxis = (/0.,0.,1./)
 else
    call cross_product3D(vecz, veczprime, rotaxis)
 endif
 !
 !-- Randomly inject particles around the body's outer 'radius'.
 !
 do i=1,npinject
    select case (wind_type)
    case (1)
       dx = get_pos_on_sphere(seed, delta_theta)
       call rotatevec(dx, rotaxis, theta_rad)
       call cross_product3D(veczprime, dx, vhat)
       vxyz      = v2 + wind_speed*vhat
    case default
       ! Get random position on sphere
       dx = get_pos_on_sphere(seed, delta_theta)
       ! Position is planet position + Bondi radius * random direction
       xyz = r2 + rinject*dx
       ! Velocity is planet velocity + sound speed * normal direction
       vxyz = v2 + wind_speed*dx
    end select
    ipart     = npart + 1
    call add_or_update_particle(igas,xyz,vxyz,h,u,ipart,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Updates the injected particles
!+
!-----------------------------------------------------------------------
subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Returns dndt(t) depending on which function is chosen
!  Note that time in this function is strictly the fraction
!  of the orbit, not absolute time
!+
!-----------------------------------------------------------------------
real function mdot_func(r,r_ref)
 real, intent(in) :: r,r_ref

 select case (mdot_type)
 case (2)
    mdot_func = (r_ref/r)**2
 case default
    mdot_func = 1.0
 end select

end function mdot_func

!-----------------------------------------------------------------------
!+
!  Returns a random location on a sperical surface
!+
!-----------------------------------------------------------------------
function get_pos_on_sphere(iseed, delta_theta) result(dx)
 use random,        only:get_random_pos_on_sphere, get_gaussian_pos_on_sphere
 integer, intent(inout) :: iseed
 real, intent(inout) :: delta_theta
 real  :: dx(3)

 select case (random_type)
 case(1)
    dx = get_gaussian_pos_on_sphere(iseed, delta_theta)
 case(0)
    dx = get_random_pos_on_sphere(iseed)
 end select

end function get_pos_on_sphere

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_type, 'wind_type', 'wind setup (0=asteroidwind, 1=randomwind, 2=boil-off)', iunit)
 if (wind_type /= 2) then
    call write_inopt(mdot_str,'mdot','mass injection rate with unit, e.g. 1e8*g/s, 1e-7M_s/yr',iunit)
    call write_inopt(npartperorbit,'npartperorbit',&
                  'particle injection rate in particles/binary orbit',iunit)
    call write_inopt(vlag,'vlag','percentage lag in velocity of wind',iunit)
    call write_inopt(mdot_type,'mdot_type','injection rate (0=const, 2=r^(-2))',iunit)
    if (mdot_type==2) then
       call write_inopt(r_ref,'r_ref','radius at whieh Mdot=mdot for 1/r^2 injection type',iunit)
    endif
 endif
 call write_inopt(random_type, 'random_type', 'random position on the surface, 0 for random, 1 for gaussian', iunit)
 if (random_type==1) then
    call write_inopt(delta_theta, 'delta_theta', 'standard deviation for the gaussion distribution', iunit)
    call write_inopt(theta, 'theta', 'the tilt inclination of the star or planet (random_type=1)', iunit)
    call write_inopt(phi, 'phi', 'the tilt orientation of the star, (random_type=1)', iunit)
 endif
 if (wind_type==1) then
    call write_inopt(inject_pt, 'inject_pt', 'the particle that excites wind (when wind_type=1)', iunit)
    call write_inopt(r_inject_str, 'r_inject', 'inject radius with units, e.g. 1*AU, 1e8m, (when wind_type=1)', iunit)
 endif
 call write_inopt(wind_speed_factor, &
 & 'wind_speed_factor', 'factor to scale the wind speed based on the Keplerian speed at rinject', iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(wind_type,'wind_type',db,errcount=nerr,min=0,max=2)
 if (wind_type /= 2) then
    call read_inopt(mdot_str,'mdot',db,errcount=nerr)
    call read_inopt(npartperorbit,'npartperorbit',db,errcount=nerr,min=0.)
    call read_inopt(vlag,'vlag',db,errcount=nerr)
    call read_inopt(mdot_type,'mdot_type',db,errcount=nerr,min=0,max=2)
    if (mdot_type==2) then
       call read_inopt(r_ref,'r_ref',db,errcount=nerr,min=0.)
    endif
 endif
 call read_inopt(random_type,'random_type',db,errcount=nerr,min=0,max=1)
 if (random_type==1) then
    call read_inopt(delta_theta,'delta_theta',db,errcount=nerr)
    call read_inopt(theta,'theta',db,errcount=nerr)
    call read_inopt(phi,'phi',db,errcount=nerr)
 endif
 if (wind_type==1) then
    call read_inopt(inject_pt,'inject_pt',db,errcount=nerr)
    call read_inopt(r_inject_str,'r_inject',db,errcount=nerr)
 endif
 call read_inopt(wind_speed_factor,'wind_speed_factor',db,errcount=nerr)

end subroutine read_options_inject

!-----------------------------------------------------------------------
!+
!  Sets default options for the injection module
!+
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

end subroutine set_default_options_inject

end module inject
