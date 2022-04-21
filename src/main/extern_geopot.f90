!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module extern_geopot
!
! Implementation of external forces from geopotential model
!
! Currently only implements J2, i.e. effect of oblateness
! but could be extended to deal with higher order terms
!
! Spin vector direction is arbitrary
!
! :References: https://en.wikipedia.org/wiki/Geopotential_model
!              Hong et al. (2021), ApJ 920, 151
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - J2    : *J2 parameter*
!
! :Dependencies: infile_utils, io, kernel, physcon
!
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, public    :: J2 = 0.
 real, public    :: tilt_angle = 0.
 real, private   :: sin_angle = 0.
 real, private   :: cos_angle = 1.

 public :: get_geopot_force
 public :: write_options_geopot, read_options_geopot
 private

contains

!------------------------------------------------
!+
!  Compute higher order terms in the acceleration
!  namely the J2 term caused by oblateness
!+
!------------------------------------------------
subroutine get_geopot_force(xi,yi,zi,dr,mdr3,Rp,fextxi,fextyi,fextzi,phi)
 real, intent(in)    :: xi,yi,zi
 real, intent(in)    :: dr    !  1/r
 real, intent(in)    :: mdr3  !  GM/r^3
 real, intent(in)    :: Rp    !  radius of bodys
 real, intent(inout) :: fextxi,fextyi,fextzi
 real, intent(inout) :: phi
 real :: spinvec(3),r_dot_s,term,term1,term2

 call get_spinvec(spinvec)

 ! Equation 1 of Hong et al. (2021)
 r_dot_s = (xi*spinvec(1) + yi*spinvec(2) + zi*spinvec(3))*dr
 term = 1.5*J2*(Rp*dr)**2*mdr3
 term1 = term*(5.*r_dot_s**2 - 1.)
 term2 = term*(-2.*r_dot_s)/dr

 fextxi = fextxi + term1*xi + term2*spinvec(1)
 fextyi = fextyi + term1*yi + term2*spinvec(2)
 fextzi = fextzi + term1*zi + term2*spinvec(3)

 ! potential is as given in wikipedia except we replace z/r with r_dot_s
 phi    = phi + 0.5*J2*(Rp**2)*mdr3*(3.*r_dot_s**2 - 1.)

end subroutine get_geopot_force

!---------------------------------------------------------------
!+
!  Define speed and direction of rotation
!  At present direction is hard-wired to rotation in x-y plane
!+
!---------------------------------------------------------------
subroutine get_spinvec(spinvec)
 real, intent(out) :: spinvec(3)

 spinvec = (/sin_angle,0.,cos_angle/)

end subroutine get_spinvec

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_geopot(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(J2,'J2','J2 value in code units',iunit)
 call write_inopt(tilt_angle,'tilt_angle','tilt angle (obliquity) in degrees',iunit)

end subroutine write_options_geopot

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_geopot(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use physcon, only:deg_to_rad
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_geopot'

 igotall = .false.
 imatch = .true.
 select case(trim(name))
 case('J2')
    read(valstring,*,iostat=ierr) J2
    ngot = ngot + 1
 case('tilt_angle')
    read(valstring,*,iostat=ierr) tilt_angle
    sin_angle = sin(tilt_angle*deg_to_rad)
    cos_angle = cos(tilt_angle*deg_to_rad)
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_geopot

end module extern_geopot
