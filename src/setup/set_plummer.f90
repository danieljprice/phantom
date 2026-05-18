!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setplummer
!
! Utility routines for sampling Plummer and Hernquist density profiles
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none

 integer, parameter, public :: iprofile_plummer    = 1
 integer, parameter, public :: iprofile_hernquist  = 2
 character(len=*), parameter, public :: profile_label(2) = &
     (/'Plummer sphere  ','Hernquist sphere'/)

 real,    parameter :: default_cut_mass = 0.999

 public :: get_accel_profile,density_profile,radius_from_mass

 private

contains

!-----------------------------------------------------------------------
!+
!  Density profile for Plummer and Hernquist spheres
!+
!-----------------------------------------------------------------------
real function density_profile(iprofile,r,rsoft,mass)
 use physcon, only:pi
 integer, intent(in) :: iprofile
 real,    intent(in) :: r,rsoft,mass

 select case(iprofile)
 case(iprofile_plummer)
    density_profile = 3.*mass*rsoft**2/(4.*pi*(rsoft**2 + r**2)**2.5)
 case(iprofile_hernquist)
    density_profile = mass*rsoft/(2.*pi*r*(rsoft + r)**3)
 case default
    density_profile = 0.
 end select

end function density_profile

!-----------------------------------------------------------------------
!+
!  Analytic gravitational acceleration for supported profiles
!+
!-----------------------------------------------------------------------
subroutine get_accel_profile(iprofile,pos,rsoft,mass,acc)
 integer, intent(in)  :: iprofile
 real,    intent(in)  :: pos(3),rsoft,mass
 real,    intent(out) :: acc(3)
 real :: r,denom

 r = sqrt(dot_product(pos,pos))
 if (r <= epsilon(0.)) then
    acc = 0.
    return
 endif

 select case(iprofile)
 case(iprofile_plummer)
    denom = (r*r + rsoft*rsoft)**1.5
    acc   = -mass*pos/denom
 case(iprofile_hernquist)
    denom = (r + rsoft)**2*r
    acc   = -mass*pos/denom
 case default
    acc = 0.
 end select

end subroutine get_accel_profile

!-----------------------------------------------------------------------
!+
!  radius corresponding to supplied mass fraction
!+
!-----------------------------------------------------------------------
pure real function radius_from_mass(iprofile,rmass,rsoft)
 integer, intent(in) :: iprofile
 real,    intent(in) :: rmass,rsoft
 real :: rmass23
 real :: fraction
 real, parameter :: eps = epsilon(0.)

 fraction = min(max(rmass,0.),1.)

 select case(iprofile)
 case(iprofile_plummer)
    rmass23 = fraction**(2./3.)
    radius_from_mass = rsoft*sqrt(rmass23/max(eps,1. - rmass23))
 case(iprofile_hernquist)
    radius_from_mass = rsoft*(fraction + sqrt(fraction))/max(eps,1. - fraction)
 case default
    radius_from_mass = 0.
 end select

end function radius_from_mass

end module setplummer

