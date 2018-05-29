module metric

!------ The Kerr metric in Boyer-Lindquist coordinates

 implicit none
 character(len=*), parameter :: metric_type = 'Kerr'
 integer,          parameter :: imetric     = 3

 real, public  :: mass1 = 1.       ! mass of central object
 real, public  :: a     = 1.       ! spin of central object

contains

!----------------------------------------------------------------
!+
!  Metric tensors
!+
!----------------------------------------------------------------

!--- The metric tensor in 'CARTESIAN-like form'
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3), sqrtg
 real :: x,y,z,x2,y2,z2,a2
 real :: r2spherical,r2,r
 real :: rho2,delta,sintheta2
 real :: gphiphi,gtphi,gtt
 real :: rs
 rs = 2.*mass1

 sqrtg = 1.

 x  = position(1)
 y  = position(2)
 z  = position(3)
 x2 = x**2
 y2 = y**2
 z2 = z**2
 a2 = a**2
 r2spherical = x2+y2+z2
 r2        = 0.5*(r2spherical-a2)+0.5*sqrt( (r2spherical-a2)**2 + 4.*a2*z2 )
 r         = sqrt(r2)
 rho2      = r2 + a2*(z2/r2)

 !--- The Boyer-Lindquist metric tensor in CARTESIAN-like form

 delta     = r2 - rs*r + a2
 sintheta2 = 1. - z2/r2
 gtt       = -1. + (r*rs)/rho2
 gphiphi   = sintheta2*(a2 + r2 + (a2*r*rs*sintheta2)/rho2)
 gtphi     = -((a*r*rs*sintheta2)/rho2)

 gcov(0,0) = -1. + (r*rs)/rho2
 gcov(1,0) = -((gtphi*y)/(x2 + y2))
 gcov(2,0) = (gtphi*x)/(x2 + y2)
 gcov(3,0) = 0.
 gcov(0,1) = -((gtphi*y)/(x2 + y2))
 gcov(1,1) = (r2*x2)/(delta*rho2) + (gphiphi*y2)/(x2 + y2)**2 + (x2*z2)/(rho2*(r2 - z2))
 gcov(2,1) = (r2*x*y)/(delta*rho2) - (gphiphi*x*y)/(x2 + y2)**2 + (x*y*z2)/(rho2*(r2 - z2))
 gcov(3,1) = ((a2 + r2)*x*z)/(delta*rho2) - (x*z*(1. - ((a2 + r2)*z2)/(r2*rho2)))/(r2 - z2)
 gcov(0,2) = (gtphi*x)/(x2 + y2)
 gcov(1,2) = (r2*x*y)/(delta*rho2) - (gphiphi*x*y)/(x2 + y2)**2 + (x*y*z2)/(rho2*(r2 - z2))
 gcov(2,2) = (r2*y2)/(delta*rho2) + (gphiphi*x2)/(x2 + y2)**2 + (y2*z2)/(rho2*(r2 - z2))
 gcov(3,2) = ((a2 + r2)*y*z)/(delta*rho2) - (y*z*(1. - ((a2 + r2)*z2)/(r2*rho2)))/(r2 - z2)
 gcov(0,3) = 0.
 gcov(1,3) = ((a2 + r2)*x*z)/(delta*rho2) - (x*z*(1. - ((a2 + r2)*z2)/(r2*rho2)))/(r2 - z2)
 gcov(2,3) = ((a2 + r2)*y*z)/(delta*rho2) - (y*z*(1. - ((a2 + r2)*z2)/(r2*rho2)))/(r2 - z2)
 gcov(3,3) = ((a2 + r2)**2*z2)/(delta*r2*rho2) + (rho2*(1. - ((a2 + r2)*z2)/(r2*rho2))**2)/(r2 - z2)

 gcon(0,0) = gphiphi/(-gtphi**2 + gphiphi*gtt)
 gcon(1,0) = -((gtphi*y)/(gtphi**2 - gphiphi*gtt))
 gcon(2,0) = (gtphi*x)/(gtphi**2 - gphiphi*gtt)
 gcon(3,0) = 0.
 gcon(0,1) = -((gtphi*y)/(gtphi**2 - gphiphi*gtt))
 gcon(1,1) = (delta*r2*x2)/((a2 + r2)**2*rho2) + (gtt*y2)/(-gtphi**2 + gphiphi*gtt) + (x2*z2)/(rho2*(r2 - z2))
 gcon(2,1) = -((gtt*x*y)/(-gtphi**2 + gphiphi*gtt)) + (delta*r2*x*y)/((a2 + r2)**2*rho2) + (x*y*z2)/(rho2*(r2 - z2))
 gcon(3,1) = -((x*z)/rho2) + (delta*x*z)/((a2 + r2)*rho2)
 gcon(0,2) = (gtphi*x)/(gtphi**2 - gphiphi*gtt)
 gcon(1,2) = -((gtt*x*y)/(-gtphi**2 + gphiphi*gtt)) + (delta*r2*x*y)/((a2 + r2)**2*rho2) + (x*y*z2)/(rho2*(r2 - z2))
 gcon(2,2) = (gtt*x2)/(-gtphi**2 + gphiphi*gtt) + (delta*r2*y2)/((a2 + r2)**2*rho2) + (y2*z2)/(rho2*(r2 - z2))
 gcon(3,2) = -((y*z)/rho2) + (delta*y*z)/((a2 + r2)*rho2)
 gcon(0,3) = 0.
 gcon(1,3) = -((x*z)/rho2) + (delta*x*z)/((a2 + r2)*rho2)
 gcon(2,3) = -((y*z)/rho2) + (delta*y*z)/((a2 + r2)*rho2)
 gcon(3,3) = (r2 - z2)/rho2 + (delta*z2)/(r2*rho2)

end subroutine get_metric_cartesian

!--- The metric tensor in SPHERICAL-like form
subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real,    intent(in)  :: position(3)
 real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 real :: a2,r2,r,rho2,delta,sintheta2
 real :: gtt,grr,gthetatheta,gtphi,gphiphi
 real :: phi,theta
 real :: rs
 rs = 2.*mass1

 a2 = a**2
 r     = position(1)
 theta = position(2)
 phi   = position(3)
 r2    = r**2
 rho2  = r2 + a2*cos(theta)**2
 delta = r2 - rs*r + a2
 sintheta2 = sin(theta)**2

 gcov(0,0) = -1. + (r*rs)/rho2
 gcov(1,0) = 0.
 gcov(2,0) = 0.
 gcov(3,0) = -((a*r*rs*sintheta2)/rho2)
 gcov(0,1) = 0.
 gcov(1,1) = rho2/delta
 gcov(2,1) = 0.
 gcov(3,1) = 0.
 gcov(0,2) = 0.
 gcov(1,2) = 0.
 gcov(2,2) = rho2
 gcov(3,2) = 0.
 gcov(0,3) = -((a*r*rs*sintheta2)/rho2)
 gcov(1,3) = 0.
 gcov(2,3) = 0.
 gcov(3,3) = sintheta2*(a2 + r2 + (a2*r*rs*sintheta2)/rho2)

 gtt         = gcov(0,0)
 grr         = gcov(1,1)
 gthetatheta = gcov(2,2)
 gtphi       = gcov(0,3)
 gphiphi     = gcov(3,3)

 gcon(0,0) = gphiphi/(-gtphi**2 + gphiphi*gtt)
 gcon(1,0) = 0.
 gcon(2,0) = 0.
 gcon(3,0) = gtphi/(gtphi**2 - gphiphi*gtt)
 gcon(0,1) = 0.
 gcon(1,1) = delta/rho2
 gcon(2,1) = 0.
 gcon(3,1) = 0.
 gcon(0,2) = 0.
 gcon(1,2) = 0.
 gcon(2,2) = 1./rho2
 gcon(3,2) = 0.
 gcon(0,3) = gtphi/(gtphi**2 - gphiphi*gtt)
 gcon(1,3) = 0.
 gcon(2,3) = 0.
 gcon(3,3) = gtt/(-gtphi**2 + gphiphi*gtt)

 sqrtg = grr*gthetatheta*(-gtphi**2+gphiphi*gtt)

end subroutine get_metric_spherical

!----------------------------------------------------------------
!+
!  Metric tensors derivatives
!+
!----------------------------------------------------------------

!--- Derivatives of the covariant 'CARTEISAN' metric
subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdx,dgcovdy,dgcovdz
 dgcovdx=0.
 dgcovdy=0.
 dgcovdz=0.
 STOP 'No analytic metric derivatives implemented'
end subroutine metric_cartesian_derivatives

!--- Derivatives of the covariant 'SPHERICAL' metric
subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi
 real :: r, theta, sintheta, costheta, rho, delta
 real :: rs
 rs = 2.*mass1

 r = position(1)
 theta = position(2)
 sintheta = sin(theta)
 costheta = cos(theta)
 rho = sqrt(r**2+a**2*costheta**2)
 delta = r**2-r*rs+a**2

 dgcovdr = 0.
 dgcovdtheta = 0.
 dgcovdphi = 0.

 dgcovdr(0,0) = (-2.*r**2*rs)/rho**4 + rs/rho**2
 dgcovdr(3,0) = (2.*a*r**2*rs*sintheta**2)/rho**4 - (a*rs*sintheta**2)/rho**2
 dgcovdr(1,1) = (2.*r)/delta - (rho**2*(2.*r - rs))/delta**2
 dgcovdr(2,2) = 2.*r
 dgcovdr(0,3) = (2.*a*r**2*rs*sintheta**2)/rho**4 - (a*rs*sintheta**2)/rho**2
 dgcovdr(3,3) = sintheta**2*(2.*r - (2.*a**2*r**2*rs*sintheta**2)/rho**4 + (a**2*rs*sintheta**2)/rho**2)

 dgcovdtheta(0,0) = (2.*a**2*costheta*r*rs*sintheta)/rho**4
 dgcovdtheta(3,0) = (-2.*a*costheta*r*rs*sintheta)/rho**2 - (2.*a**3*costheta*r*rs*sintheta**3)/rho**4
 dgcovdtheta(1,1) = (-2.*a**2*costheta*sintheta)/delta
 dgcovdtheta(2,2) = -2.*a**2*costheta*sintheta
 dgcovdtheta(0,3) = (-2.*a*costheta*r*rs*sintheta)/rho**2 - (2.*a**3*costheta*r*rs*sintheta**3)/rho**4
 dgcovdtheta(3,3) = 2.*costheta*sintheta*(a**2 + r**2 + (a**2*r*rs*sintheta**2)/rho**2)                                     &
 &                  + sintheta**2*((2.*a**2*costheta*r*rs*sintheta)/rho**2 + (2.*a**4*costheta*r*rs*sintheta**3)/rho**4)

end subroutine metric_spherical_derivatives


!----------------------------------------------------------------
!+
!  Coordinate transformations
!+
!----------------------------------------------------------------

!--- (Jacobian tensor) Derivatives of Boyer-Lindquist 'Spherical' with respect to 'Cartesian' coordinates
subroutine get_jacobian(position,dxdx)
 real, intent(in), dimension(3) :: position
 real, intent(out), dimension(0:3,0:3) :: dxdx
 real, dimension(3) :: dBLdx,dBLdy,dBLdz
 real :: drdx,drdy,drdz
 real :: dthetadx,dthetady,dthetadz
 real :: dphidx,dphidy,dphidz
 real :: x,y,z,x2,y2,z2
 real :: a2,r2spherical,r2,r,rho2,delta
 real :: sintheta
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
 r2        = 0.5*(r2spherical-a2)+0.5*sqrt( (r2spherical-a2)**2 + 4.*a2*z2 )
 r         = sqrt(r2)
 rho2      = r2 + a2*(z2/r2)
 delta     = r2 - rs*r + a2
 sintheta = sqrt(1.-z2/r2)

 drdx = r*x/rho2
 drdy = r*y/rho2
 drdz = (a2+r2)*z/(r*rho2)

 dthetadx = x*z/(r*rho2*sintheta)
 dthetady = y*z/(r*rho2*sintheta)
 dthetadz = -1./(r*sintheta)*(1.-(a2+r2)*z2/(r2*rho2))

 dphidx = -y/(x2+y2)
 dphidy = x/(x2+y2)
 dphidz = 0.

 dBLdx=(/drdx,dthetadx,dphidx/)
 dBLdy=(/drdy,dthetady,dphidy/)
 dBLdz=(/drdz,dthetadz,dphidz/)

 ! Temporal coordinates are the same.
 dxdx        = 0.
 dxdx(0,0)   = 1.
 dxdx(1:3,1) = dBLdx
 dxdx(1:3,2) = dBLdy
 dxdx(1:3,3) = dBLdz

end subroutine get_jacobian

!--- Boyer-Lindquist coordinate transformations from CARTEISAN to SPHERICAL
subroutine cartesian2spherical(xcart,xspher)
 real, intent(in) :: xcart(3)
 real, intent(out) ::xspher(3)
 real :: x,y,z,x2,y2,z2,a2,r2spherical,r2,r
 real :: theta,phi

 x  = xcart(1)
 y  = xcart(2)
 z  = xcart(3)
 x2 = x**2
 y2 = y**2
 z2 = z**2
 a2 = a**2
 r2spherical = x2+y2+z2
 r2        = 0.5*(r2spherical-a2)+0.5*sqrt( (r2spherical-a2)**2 + 4.*a2*z2 )
 r         = sqrt(r2)
 theta     = acos(z/r)
 phi       = atan2(y,x)
 xspher   = (/r,theta,phi/)

end subroutine cartesian2spherical

!--- Boyer-Lindquist coordinate transformations from SPHERICAL to CARTEISAN
subroutine spherical2cartesian(xspher,xcart)
 real, intent(in) :: xspher(3)
 real, intent(out) :: xcart(3)
 real :: x,y,z,r,theta,phi,r2,a2

 a2 = a**2
 r  = xspher(1)
 r2 = r**2
 theta = xspher(2)
 phi   = xspher(3)
 x = sqrt(r2+a2)*sin(theta)*cos(phi)
 y = sqrt(r2+a2)*sin(theta)*sin(phi)
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
 call write_inopt(a,'a','spin parameter for Kerr metric',iunit)

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,warning
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
    if (mass1 == 0.) call warning(tag,'black hole mass: mass1 = 0')
    ngot = ngot + 1
 case('a')
    read(valstring,*,iostat=ierr) a
    if (abs(a) > 1.)  call fatal(tag,'black hole spin: |a| > 1')
    if (a == 0.) call warning(tag,'black hole spin: a = 0')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_metric

end module metric
