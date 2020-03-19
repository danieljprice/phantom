module metric
 implicit none
 character(len=*), parameter :: metric_type = 'Schwarzschild'
 integer,          parameter :: imetric     = 2

 real, public :: mass1 = 1.       ! mass of central object

contains

!----------------------------------------------------------------
!+
!  Metric tensors
!+
!----------------------------------------------------------------

!--- The metric tensor in 'CARTESIAN-like form'
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
 real, intent(in) :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: r,r2,r3,rs_on_r3,coeff,x,y,z,x2,y2,z2,term
 real :: rs
 rs = 2.*mass1

 r2 = dot_product(position,position)
 r  = sqrt(r2)
 r3 = r*r2
 rs_on_r3 = rs/r3
 x  = position(1)
 y  = position(2)
 z  = position(3)
 x2 = x**2
 y2 = y**2
 z2 = z**2


 !--- The Schwarzschild metric tensor in CARTESIAN-like form
 if (present(sqrtg)) sqrtg = 1.

 term  = (1.-rs/r)
 coeff = 1./term

 gcov = 0.
 gcov(0,0) = -term

 gcov(1,1) = coeff*(1.-rs_on_r3*(y2+z2))
 gcov(2,1) = coeff*x*y*rs_on_r3
 gcov(3,1) = coeff*x*z*rs_on_r3

 gcov(1,2) = gcov(2,1)
 gcov(2,2) = coeff*(1.-rs_on_r3*(x2+z2))
 gcov(3,2) = coeff*y*z*rs_on_r3

 gcov(1,3) = gcov(3,1)
 gcov(2,3) = gcov(3,2)
 gcov(3,3) = coeff*(1.-rs_on_r3*(x2+y2))

 if (present(gcon)) then
    gcon      = 0.
    gcon(0,0) = -1./term

    gcon(1,1) = 1.-rs_on_r3*x2
    gcon(2,1) = -rs_on_r3*x*y
    gcon(3,1) = -rs_on_r3*x*z

    gcon(1,2) = gcon(2,1)
    gcon(2,2) = 1.-rs_on_r3*y2
    gcon(3,2) = -rs_on_r3*y*z

    gcon(1,3) = gcon(3,1)
    gcon(2,3) = gcon(3,2)
    gcon(3,3) = 1.-rs_on_r3*z2
 endif

end subroutine get_metric_cartesian

pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: r,theta,sintheta,r2
 real :: rs
 rs = 2.*mass1

 r = position(1)
 r2 = r**2
 theta = position(2)
 sintheta = sin(theta)

 gcov(0,0) = -1. + rs/r
 gcov(1,0) = 0.
 gcov(2,0) = 0.
 gcov(3,0) = 0.
 gcov(0,1) = 0.
 gcov(1,1) = 1./(1. - rs/r)
 gcov(2,1) = 0.
 gcov(3,1) = 0.
 gcov(0,2) = 0.
 gcov(1,2) = 0.
 gcov(2,2) = r2
 gcov(3,2) = 0.
 gcov(0,3) = 0.
 gcov(1,3) = 0.
 gcov(2,3) = 0.
 gcov(3,3) = r2*sintheta**2

 if (present(gcon)) then
    gcon(0,0) = -(r/(r - rs))
    gcon(1,0) = 0.
    gcon(2,0) = 0.
    gcon(3,0) = 0.
    gcon(0,1) = 0.
    gcon(1,1) = 1. - rs/r
    gcon(2,1) = 0.
    gcon(3,1) = 0.
    gcon(0,2) = 0.
    gcon(1,2) = 0.
    gcon(2,2) = 1./r2
    gcon(3,2) = 0.
    gcon(0,3) = 0.
    gcon(1,3) = 0.
    gcon(2,3) = 0.
    gcon(3,3) = 1./(r2*sintheta**2)
 endif

 if (present(sqrtg)) sqrtg = r2*sintheta

end subroutine get_metric_spherical

!----------------------------------------------------------------
!+
!  Metric tensors derivatives
!+
!----------------------------------------------------------------

!--- Derivatives of the covariant 'CARTEISAN' metric
pure subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
 real,    intent(in)  :: position(3)
 real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 real :: x,y,z,r,r2,r3,r4,r5,rs_on_r3,x2,y2,z2,rs2
 real :: rs
 rs = 2.*mass1

 dgcovdx = 0.
 dgcovdy = 0.
 dgcovdz = 0.
 x = position(1)
 y = position(2)
 z = position(3)
 x2= x**2
 y2= y**2
 z2= z**2
 r2 = dot_product(position,position)
 r  = sqrt(r2)
 r3 = r*r2
 r4 = r2*r2
 r5 = r*r4

 rs_on_r3 = rs/r3
 rs2 = rs**2

 ! dx
 dgcovdx(0,0) = -((rs*x)/r3)
 dgcovdx(1,0) = 0.
 dgcovdx(2,0) = 0.
 dgcovdx(3,0) = 0.

 dgcovdx(0,1) = 0.
 dgcovdx(1,1) = (x*(2*r4 - 2*r3*rs - 2*rs2*(y2 + z2) - 2*r4 + r*rs*(x2 + 4*(y2 + z2))))/(r4*(r - rs)**2)
 dgcovdx(2,1) = (rs*(r3 - r2*rs - 3*r*x2 + 2*rs*x2)*y)/(r4*(r - rs)**2)
 dgcovdx(3,1) = (rs*(r3 - r2*rs - 3*r*x2 + 2*rs*x2)*z)/(r4*(r - rs)**2)

 dgcovdx(0,2) = 0.
 dgcovdx(1,2) = dgcovdx(2,1)
 dgcovdx(2,2) = (x*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) - 2*rs2*(x2 + z2) + r*rs*(4*x2 + y2 + 4*z2)))/(r4*(r - rs)**2)
 dgcovdx(3,2) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)

 dgcovdx(0,3) = 0.
 dgcovdx(1,3) = dgcovdx(3,1)
 dgcovdx(2,3) = dgcovdx(3,2)
 dgcovdx(3,3) = (-2*(r - rs)**2*x*(-z2) + r*(-2*r + rs)*x*z2)/(r4*(r - rs)**2)

 ! dy
 dgcovdy(0,0) = -((rs*y)/r3)
 dgcovdy(1,0) = 0.
 dgcovdy(2,0) = 0.
 dgcovdy(3,0) = 0.

 dgcovdy(0,1) = 0.
 dgcovdy(1,1) = (y*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) - 2*rs2*(y2 + z2) + r*rs*(x2 + 4*(y2 + z2))))/(r4*(r - rs)**2)
 dgcovdy(2,1) = (rs*x*(r3 - r2*rs - 3*r*y2 + 2*rs*y2))/(r4*(r - rs)**2)
 dgcovdy(3,1) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)

 dgcovdy(0,2) = 0.
 dgcovdy(1,2) = dgcovdy(2,1)
 dgcovdy(2,2) = (y*(2*r4 - 2*r3*rs - 2*rs2*(x2 + z2) - 2*r4 + r*rs*(4*x2 + y2 + 4*z2)))/(r4*(r - rs)**2)
 dgcovdy(3,2) = (rs*(r3 - r2*rs - 3*r*y2 + 2*rs*y2)*z)/(r4*(r - rs)**2)

 dgcovdy(0,3) = 0.
 dgcovdy(1,3) = dgcovdy(3,1)
 dgcovdy(2,3) = dgcovdy(3,2)
 dgcovdy(3,3) = (-2*(r - rs)**2*y*(-z2) + r*(-2*r + rs)*y*z2)/(r4*(r - rs)**2)

 ! dz
 dgcovdz(0,0) = -((rs*z)/r3)
 dgcovdz(1,0) = 0.
 dgcovdz(2,0) = 0.
 dgcovdz(3,0) = 0.

 dgcovdz(0,1) = 0.
 dgcovdz(1,1) = (z*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) - 2*rs2*(y2 + z2) + r*rs*(x2 + 4*(y2 + z2))))/(r4*(r - rs)**2)
 dgcovdz(2,1) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)
 dgcovdz(3,1) = (rs*x*(r3 - r2*rs - 3*r*z2 + 2*rs*z2))/(r4*(r - rs)**2)

 dgcovdz(0,2) = 0.
 dgcovdz(1,2) = dgcovdz(2,1)
 dgcovdz(2,2) = (z*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) -2*rs2*(x2 + z2) + r*rs*(4*x2 + y2 + 4*z2)))/(r4*(r - rs)**2)
 dgcovdz(3,2) = (rs*y*(r3 - r2*rs - 3*r*z2 + 2*rs*z2))/(r4*(r - rs)**2)

 dgcovdz(0,3) = 0.
 dgcovdz(1,3) = dgcovdz(3,1)
 dgcovdz(2,3) = dgcovdz(3,2)
 dgcovdz(3,3) = (z*(2*(r - rs)*(r3 - r*(x2 + y2) + rs*(x2 + y2)) + r*(-2*r + rs)*z2))/(r4*(r - rs)**2)

end subroutine metric_cartesian_derivatives

!--- Derivatives of the covariant 'SPHERICAL' metric
pure subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi
 real :: r, theta
 real :: rs
 rs = 2.*mass1

 r = position(1)
 theta = position(2)

 dgcovdr = 0.
 dgcovdtheta = 0.
 dgcovdphi = 0.

 dgcovdr(0,0) = -rs/r**2
 dgcovdr(1,1) = -rs/(r-rs)**2
 dgcovdr(2,2) = 2.*r
 dgcovdr(3,3) = 2*r*sin(theta)**2

 dgcovdtheta(3,3) = 2.*r**2*cos(theta)*sin(theta)

end subroutine metric_spherical_derivatives

!----------------------------------------------------------------
!+
!  Coordinate transformations
!+
!----------------------------------------------------------------

!--- (Jacobian tensor) Derivatives of Schwarzschild 'Spherical' with respect to 'Cartesian' coordinates
pure subroutine get_jacobian(position,dxdx)
 real, intent(in), dimension(3) :: position
 real, intent(out), dimension(0:3,0:3) :: dxdx
 real, dimension(3) :: dSPHERICALdx,dSPHERICALdy,dSPHERICALdz
 real :: drdx,drdy,drdz
 real :: dthetadx,dthetady,dthetadz
 real :: dphidx,dphidy,dphidz
 real :: x,y,z,x2,y2,z2,r2,r,rcyl2,rcyl

 x  = position(1)
 y  = position(2)
 z  = position(3)
 x2 = x**2
 y2 = y**2
 z2 = z**2
 r2 = x2+y2+z2
 r  = sqrt(r2)
 rcyl2 = x2+y2
 rcyl  = sqrt(x2+y2)

 drdx = x/r
 drdy = y/r
 drdz = z/r

 dthetadx = x*z/(r2*rcyl)
 dthetady = y*z/(r2*rcyl)
 dthetadz = -rcyl/r2

 dphidx = -y/(x2+y2)
 dphidy = x/(x2+y2)
 dphidz = 0.

 dSPHERICALdx=(/drdx,dthetadx,dphidx/)
 dSPHERICALdy=(/drdy,dthetady,dphidy/)
 dSPHERICALdz=(/drdz,dthetadz,dphidz/)

 dxdx        = 0.
 dxdx(0,0)   = 1.
 dxdx(1:3,1) = dSPHERICALdx
 dxdx(1:3,2) = dSPHERICALdy
 dxdx(1:3,3) = dSPHERICALdz

end subroutine get_jacobian

pure subroutine cartesian2spherical(xcart,xspher)
 real, intent(in) :: xcart(3)
 real, intent(out) ::xspher(3)
 real :: x,y,z
 real :: r,theta,phi

 x  = xcart(1)
 y  = xcart(2)
 z  = xcart(3)

 r  = sqrt(x**2+y**2+z**2)
 theta = acos(z/r)
 phi   = atan2(y,x)

 xspher   = (/r,theta,phi/)

end subroutine cartesian2spherical

pure subroutine spherical2cartesian(xspher,xcart)
 real, intent(in) :: xspher(3)
 real, intent(out) :: xcart(3)
 real :: x,y,z,r,theta,phi

 r     = xspher(1)
 theta = xspher(2)
 phi   = xspher(3)
 x = r*sin(theta)*cos(phi)
 y = r*sin(theta)*sin(phi)
 z = r*cos(theta)

 xcart = (/x,y,z/)

end subroutine spherical2cartesian

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to the '//trim(metric_type)//' metric'

 call write_inopt(mass1,'mass1','black hole mass in code units',iunit)

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,warn
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 character(len=*), parameter :: tag = 'metric'
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('mass1')
    read(valstring,*,iostat=ierr) mass1
    if (mass1 < 0.)  call fatal(tag,'black hole mass: mass1 < 0')
    if (mass1 == 0.) call warn(tag,'black hole mass: mass1 = 0')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_metric

end module metric
