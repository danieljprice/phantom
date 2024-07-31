!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Owner: David Liptai
!
! :Runtime parameters:
!   - delta_theta : *standard deviation for the gaussion distribution (random_type=1)*
!   - mdot        : *mass injection rate with unit, e.g. 1e8*g/s, 1e-7M_s/yr*
!   - mdot_type   : *injection rate (0=const, 1=cos(t), 2=r^(-2))*
!   - random_type : *random position on the surface, 0 for random, 1 for gaussian*
!   - vlag        : *percentage lag in velocity of wind*
!
! :Dependencies: binaryutils, externalforces, infile_utils, io, options,
!   part, partinject, physcon, random, units
!
 use io, only:error
 use physcon, only:pi
 implicit none
 character(len=*), parameter, public :: inject_type = 'randomwind'
 character(len=20), public :: mdot_str = "5.e8*g/s"
 real, public :: mdot = 1.e8     ! mass injection rate in grams/second

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par

 private

 real         :: npartperorbit = 1000.     ! particle injection rate in particles per orbit
 real         :: vlag          = 0.0      ! percentage lag in velocity of wind
 integer      :: mdot_type     = 0        ! injection rate (0=const, 1=cos(t), 2=r^(-2))
 integer      :: random_type   = 0        ! random position on the surface, 0 for random, 1 for gaussian
 real         :: delta_theta   = 0.5      ! standard deviation for the gaussion distribution (random_type=1)
 real         :: have_injected = 0.
 real         :: t_old = 0.

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(inout) :: ierr

 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Inject particles
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,            only:fatal
 use part,          only:nptmass,massoftype,igas,hfact,ihsoft
 use partinject,    only:add_or_update_particle
 use physcon,       only:twopi,gg,kboltz,mass_proton_cgs
 use random,        only:get_random_pos_on_sphere, get_gaussian_pos_on_sphere
 use units,         only:in_code_units
 use options,       only:iexternalforce
 use externalforces,only:mass1
 use binaryutils,   only:get_orbit_bits
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: ierr
 real,    dimension(3)  :: xyz,vxyz,r1,r2,v2,vhat,v1
 integer :: i,ipart,npinject,seed,pt
 real    :: dmdt,rbody,h,u,speed,inject_this_step
 real    :: m1,m2,r
 real    :: dt

 if (nptmass < 1 .and. iexternalforce == 0) &
    call fatal('inject_randomwind','not enough point masses for random wind injection')
 if (nptmass > 2) &
    call fatal('inject_randomwind','too many point masses for random wind injection')

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
 rbody     = xyzmh_ptmass(ihsoft,pt)
 m2        = xyzmh_ptmass(4,pt)
 v2        = vxyz_ptmass(1:3,pt)

 speed     = sqrt(dot_product(v2,v2))
 vhat      = v2/speed

 r         = sqrt(dot_product(r1-r2,r1-r2))

 !
 ! Add any dependency on radius to mass injection rate (and convert to code units)
 !
 mdot = in_code_units(mdot_str,ierr)
 dmdt = mdot*mdot_func(r,rbody) ! Use rbody as r_ref, currently the softening length of the body

 !
 !-- How many particles do we need to inject?
 !   (Seems to need at least eight gas particles to not crash) <-- This statement may or may not be true...
 !
 if (npartoftype(igas) < 8) then
    npinject = 8-npartoftype(igas)
 else
    ! Calculate how many extra particles from previous step to now
    dt = time - t_old
    inject_this_step = dt*dmdt/massoftype(igas)

    npinject = max(0, int(0.5 + have_injected + inject_this_step - npartoftype(igas) ))

    ! Save for next step (faster than integrating the whole thing each time)
    t_old = time
    have_injected = have_injected + inject_this_step
 endif

 !
 !-- Randomly inject particles around the body's outer 'radius'.
 !
 do i=1,npinject
    xyz       = r2 + rbody*get_pos_on_sphere(seed, delta_theta)
    vxyz      = (1.-vlag/100)*speed*vhat
    u         = 0. ! setup is isothermal so utherm is not stored
    h         = hfact*(rbody/2.)
    ipart     = npart + 1
    call add_or_update_particle(igas,xyz,vxyz,h,u,ipart,npart,npartoftype,xyzh,vxyzu)
 enddo

 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)

end subroutine inject_particles

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

 call write_inopt(mdot_str,'mdot','mass injection rate with unit, e.g. 1e8*g/s, 1e-7M_s/yr',iunit)
 call write_inopt(npartperorbit,'npartperorbit',&
                  'particle injection rate in particles/binary orbit',iunit)
 call write_inopt(vlag,'vlag','percentage lag in velocity of wind',iunit)
 call write_inopt(mdot_type,'mdot_type','injection rate (0=const, 1=cos(t), 2=r^(-2))',iunit)
 call write_inopt(random_type, 'random_type', 'random position on the surface, 0 for random, 1 for gaussian', iunit)
 call write_inopt(delta_theta, 'delta_theta', 'standard deviation for the gaussion distribution (random_type=1)', iunit)

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
    read(valstring,'(A)',iostat=ierr) mdot_str
    ngot = ngot + 1
    ! if (mdot  <  0.) call fatal(label,'mdot < 0 in input options')
 case('npartperorbit')
    read(valstring,*,iostat=ierr) npartperorbit
    ngot = ngot + 1
    if (npartperorbit < 0.) call fatal(label,'npartperorbit < 0 in input options')
 case('vlag')
    read(valstring,*,iostat=ierr) vlag
    ngot = ngot + 1
 case('mdot_type')
    read(valstring,*,iostat=ierr) mdot_type
    ngot = ngot + 1
 case('random_type')
    read(valstring,*,iostat=ierr) random_type
    ngot = ngot + 1
 case('delta_theta')
    read(valstring,*,iostat=ierr) delta_theta
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

end subroutine set_default_options_inject

end module inject
