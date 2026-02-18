!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric
!
! Kerr metric in Kerr-Schild coordinates
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dump_utils, infile_utils
!

!------ The Kerr metric in Kerr-Schild coordinates

 implicit none
 character(len=*), parameter :: metric_type = 'Kerr-Schild'
 integer,          parameter :: imetric     = 4

 real, public  :: mass1 = 1.       ! mass of central object
 real, public  :: a     = 0.       ! spin of central object

contains

!----------------------------------------------------------------
!+
!  Kerr-Schild metric tensor in CARTESIAN-like form
!+
!----------------------------------------------------------------
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: x,y,z,x2,y2,z2,a2
 real :: r2spherical,r2,r
 real :: rho2,r2a2,term
 real :: rs

 rs = 2.*mass1

 x  = position(1)
 y  = position(2)
 z  = position(3)
 x2 = x**2
 y2 = y**2
 z2 = z**2
 a2 = a**2
 r2spherical = x2+y2+z2
 r2        = 0.5*(r2spherical-a2+sqrt((r2spherical-a2))**2 + 4.*a2*z2)
 r         = sqrt(r2)
 rho2      = r2 + a2*(z2/r2)

 r2a2 = r2+a2
 term = rs/rho2

 gcov(0,0) = -1. + term*r                                 !gtt
 gcov(0,1) = term*r*(r*x+a*y)/r2a2                        !gtx
 gcov(0,2) = term*r*(r*y-a*x)/r2a2                        !gty
 gcov(0,3) = term*z                                       !gtz
 gcov(1,0) = gcov(0,1)                                    !gxt
 gcov(1,1) = 1.+term*r*((r*x+a*y)/r2a2)**2                !gxx
 gcov(1,2) = term*r*(r*x+a*y)*(r*y-a*x)/r2a2**2           !gxy
 gcov(1,3) = term*z*(r*x+a*y)/r2a2                        !gxz
 gcov(2,0) = gcov(0,2)                                    !gyt
 gcov(2,1) = gcov(1,2)                                    !gyx
 gcov(2,2) = 1.+term*r*((r*y-a*x)/r2a2)**2                !gyy
 gcov(2,3) = term*z*(r*y-a*x)/r2a2                        !gyz
 gcov(3,0) = gcov(0,3)                                    !gzt
 gcov(3,1) = gcov(1,3)                                    !gzx
 gcov(3,2) = gcov(2,3)                                    !gzy
 gcov(3,3) = 1.+term*z2/r                                 !gzz

 if (present(gcon)) then
    gcon(0,0) = -1. - term*r                                 !gtt
    gcon(0,1) = term*r*(r*x+a*y)/r2a2                        !gtx
    gcon(0,2) = term*r*(r*y-a*x)/r2a2                        !gty
    gcon(0,3) = term*z                                       !gtz
    gcon(1,0) = gcon(0,1)                                    !gxt
    gcon(1,1) = 1.-term*r*((r*x+a*y)/r2a2)**2                !gxx
    gcon(1,2) = -term*r*(r*x+a*y)*(r*y-a*x)/r2a2**2          !gxy
    gcon(1,3) = -term*z*(r*x+a*y)/r2a2                       !gxz
    gcon(2,0) = gcon(0,2)                                    !gyt
    gcon(2,1) = gcon(1,2)                                    !gyx
    gcon(2,2) = 1.-term*r*((r*y-a*x)/r2a2)**2                !gyy
    gcon(2,3) = -term*z*(r*y-a*x)/r2a2                       !gyz
    gcon(3,0) = gcon(0,3)                                    !gzt
    gcon(3,1) = gcon(1,3)                                    !gzx
    gcon(3,2) = gcon(2,3)                                    !gzy
    gcon(3,3) = 1.-term*z2/r                                 !gzz
 endif

end subroutine get_metric_cartesian

!-----------------------------------------------------------------------
!+
!  dummy routine to get the metric in spherical coordinates (not used)
!+
!-----------------------------------------------------------------------
pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg

 gcov = 0.
 if (present(gcon)) gcon = 0.
 if (present(sqrtg)) sqrtg = 1.

end subroutine get_metric_spherical

!-----------------------------------------------------------------------
!+
!  cartesian metric derivatives (not used, must be done numerically)
!+
!-----------------------------------------------------------------------
pure subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
 real,    intent(in)  :: position(3)
 real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)

 dgcovdx = 0.
 dgcovdy = 0.
 dgcovdz = 0.

end subroutine metric_cartesian_derivatives

!-----------------------------------------------------------------------
!+
!  dummy routine for spherical metric derivatives, not used
!+
!-----------------------------------------------------------------------
pure subroutine metric_spherical_derivatives(position,dgcovdr,dgcovdtheta,dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi

 dgcovdr = 0.
 dgcovdtheta = 0.
 dgcovdphi = 0.

end subroutine metric_spherical_derivatives

!-----------------------------------------------------------------------
!+
!  dummy routine to convert cartesian to spherical coordinates
!+
!-----------------------------------------------------------------------
pure subroutine cartesian2spherical(xcart,xspher)
 real, intent(in)  :: xcart(3)
 real, intent(out) :: xspher(3)

 xspher = xcart

end subroutine cartesian2spherical

!-------------------------------------------------------------------------------
!+
!  Subroutine to update the metric inputs if time dependent
!+
!-------------------------------------------------------------------------------
subroutine update_metric(time)
 real, intent(in) :: time

end subroutine update_metric

!-----------------------------------------------------------------------
!+
!  Check if a particle should be accreted by the black hole
!+
!-----------------------------------------------------------------------
subroutine accrete_particles_metric(xi,yi,zi,mi,ti,accradius,accreted)
 real,    intent(in)  :: xi,yi,zi,mi,ti,accradius
 logical, intent(out) :: accreted

 accreted = .false.

end subroutine accrete_particles_metric

!-----------------------------------------------------------------------
!+
!  writes relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_metric(hdr,time,accradius,ierr)
 use dump_utils, only:dump_h
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time,accradius
 integer,      intent(out)   :: ierr

 ierr = 0

end subroutine write_headeropts_metric

!-----------------------------------------------------------------------
!+
!  reads relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_metric(hdr,ierr)
 use dump_utils, only:dump_h
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr  = 0

end subroutine read_headeropts_metric

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 !use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 !write(iunit,"(/,a)") '# There are no options relating to the '//trim(metric_type)//' metric'

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(db,nerr)
 use infile_utils, only:inopts!,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 !call read_inopt(metric_file,'metric_file',db,errcount=nerr)

end subroutine read_options_metric

end module metric
