!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_geopot
!
! Implementation of external forces from geopotential model
!
! Currently only implements J2, i.e. effect of oblateness
! but could be extended to deal with higher order terms
!
! Spin vector direction is specified by tilt_angle
!
! :References: https://en.wikipedia.org/wiki/Geopotential_model
!              Hong et al. (2021), ApJ 920, 151
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - J2         : *J2 value in code units*
!   - tilt_angle : *tilt angle (obliquity) in degrees*
!
! :Dependencies: infile_utils, io, physcon
!
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, public  :: J2 = 0.
 real, public  :: spinvec(3) = 0.
 real, private :: tilt_angle = 0.

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
subroutine get_geopot_force(xi,yi,zi,dr,mdr3,Rp,J2i,si,fxi,fyi,fzi,phi,&
                            dsx,dsy,dsz,fxj,fyj,fzj)
 real, intent(in)    :: xi,yi,zi
 real, intent(in)    :: dr    !  1/r
 real, intent(in)    :: mdr3  !  GM/r^3
 real, intent(in)    :: Rp    !  radius of body
 real, intent(in)    :: J2i   !  J2
 real, intent(in)    :: si(3) ! unit spin vector
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(inout) :: phi
 real, intent(inout), optional :: dsx,dsy,dsz,fxj,fyj,fzj
 real :: r_dot_s,term,term1,term2

 ! Equation 1 of Hong et al. (2021)
 r_dot_s = (xi*si(1) + yi*si(2) + zi*si(3))*dr
 term = 1.5*J2i*(Rp*dr)**2*mdr3
 term1 = term*(5.*r_dot_s**2 - 1.)
 term2 = term*(-2.*r_dot_s)/dr

 fxi = fxi + term1*xi + term2*si(1)
 fyi = fyi + term1*yi + term2*si(2)
 fzi = fzi + term1*zi + term2*si(3)
 if (present(dsx)) then
    ! reaction torque on extended body (time derivative of spin, r x F)
    dsx = dsx - term2*(yi*si(3) - zi*si(2))
    dsy = dsy - term2*(zi*si(1) - xi*si(3))
    dsz = dsz - term2*(xi*si(2) - yi*si(1))
 endif
 if (present(fxj)) then
    ! acceleration on j due to i, needs to be multiplied by mi/mj later
    fxj = fxj - term1*xi - term2*si(1)
    fyj = fyj - term1*yi - term2*si(2)
    fzj = fzj - term1*zi - term2*si(3)
 endif

 ! potential is as given in wikipedia except we replace z/r with r_dot_s
 phi = phi + 0.5*J2*(Rp**2)*mdr3*(3.*r_dot_s**2 - 1.)

end subroutine get_geopot_force

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
    spinvec = (/sin(tilt_angle*deg_to_rad),0.,cos(tilt_angle*deg_to_rad)/)
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_geopot

end module extern_geopot
