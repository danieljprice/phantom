module metric

!------ The Kerr metric in Boyer-Lindquist coordinates

 implicit none
 character(len=*), parameter :: metric_type = 'Kerr'
 integer,          parameter :: imetric     = 3

 real, public  :: mass1 = 1.       ! mass of central object
 real, public  :: a     = 0.9       ! spin of central object

contains

!----------------------------------------------------------------
!+
!  Metric tensors
!+
!----------------------------------------------------------------

!--- The metric tensor in 'CARTESIAN-like form'
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: x,y,z,x2,y2,z2,a2
 real :: r2spherical,r2,r
 real :: rho2,delta,sintheta2
 real :: gphiphi,gtphi,gtt
 real :: rs
 real :: a2pr2,term1,dx2py2,drho2delta,dr2,drho2,dr2mz2
 rs = 2.*mass1

 if (present(sqrtg)) sqrtg = 1.

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
 dr2       = 1./r2
 rho2      = r2 + a2*(z2*dr2)
 drho2     = 1./rho2
 a2pr2     = a2 + r2

 !--- The Boyer-Lindquist metric tensor in CARTESIAN-like form

 delta     = a2pr2 - rs*r
 sintheta2 = 1. - z2*dr2
 gtt       = -1. + (r*rs)*drho2
 gphiphi   = sintheta2*(a2pr2 + (a2*r*rs*sintheta2)*drho2)
 gtphi     = -((a*r*rs*sintheta2)*drho2)

 term1     = 1. - (a2pr2*z2)*dr2*drho2
 dx2py2    = 1./(x2 + y2)
 drho2delta = drho2/delta
 dr2mz2    = 1./(r2 - z2)

 gcov(0,0) = -1. + (r*rs)*drho2
 gcov(1,0) = -((gtphi*y)*dx2py2)
 gcov(2,0) = (gtphi*x)*dx2py2
 gcov(3,0) = 0.
 gcov(0,1) = -((gtphi*y)*dx2py2)
 gcov(1,1) = (r2*x2)*drho2delta + (gphiphi*y2)*dx2py2**2 + (x2*z2)*drho2*dr2mz2
 gcov(2,1) = (r2*x*y)*drho2delta - (gphiphi*x*y)*dx2py2**2 + (x*y*z2)*drho2*dr2mz2
 gcov(3,1) = (a2pr2*x*z)*drho2delta - (x*z*term1)*dr2mz2
 gcov(0,2) = (gtphi*x)*dx2py2
 gcov(1,2) = (r2*x*y)*drho2delta - (gphiphi*x*y)*dx2py2**2 + (x*y*z2)*drho2*dr2mz2
 gcov(2,2) = (r2*y2)*drho2delta + (gphiphi*x2)*dx2py2**2 + (y2*z2)*drho2*dr2mz2
 gcov(3,2) = (a2pr2*y*z)*drho2delta - (y*z*term1)*dr2mz2
 gcov(0,3) = 0.
 gcov(1,3) = (a2pr2*x*z)*drho2delta - (x*z*term1)*dr2mz2
 gcov(2,3) = (a2pr2*y*z)*drho2delta - (y*z*term1)*dr2mz2
 gcov(3,3) = (a2pr2**2*z2)*dr2*drho2delta + (rho2*term1**2)*dr2mz2

 if (present(gcon)) then
    gcon(0,0) = gphiphi/(-gtphi**2 + gphiphi*gtt)
    gcon(1,0) = -((gtphi*y)/(gtphi**2 - gphiphi*gtt))
    gcon(2,0) = (gtphi*x)/(gtphi**2 - gphiphi*gtt)
    gcon(3,0) = 0.
    gcon(0,1) = -((gtphi*y)/(gtphi**2 - gphiphi*gtt))
    gcon(1,1) = (delta*r2*x2)/(a2pr2**2*rho2) + (gtt*y2)/(-gtphi**2 + gphiphi*gtt) + (x2*z2)/(rho2*(r2 - z2))
    gcon(2,1) = -((gtt*x*y)/(-gtphi**2 + gphiphi*gtt)) + (delta*r2*x*y)/(a2pr2**2*rho2) + (x*y*z2)/(rho2*(r2 - z2))
    gcon(3,1) = -((x*z)*drho2) + (delta*x*z)/(a2pr2*rho2)
    gcon(0,2) = (gtphi*x)/(gtphi**2 - gphiphi*gtt)
    gcon(1,2) = -((gtt*x*y)/(-gtphi**2 + gphiphi*gtt)) + (delta*r2*x*y)/(a2pr2**2*rho2) + (x*y*z2)/(rho2*(r2 - z2))
    gcon(2,2) = (gtt*x2)/(-gtphi**2 + gphiphi*gtt) + (delta*r2*y2)/(a2pr2**2*rho2) + (y2*z2)/(rho2*(r2 - z2))
    gcon(3,2) = -((y*z)*drho2) + (delta*y*z)/(a2pr2*rho2)
    gcon(0,3) = 0.
    gcon(1,3) = -((x*z)*drho2) + (delta*x*z)/(a2pr2*rho2)
    gcon(2,3) = -((y*z)*drho2) + (delta*y*z)/(a2pr2*rho2)
    gcon(3,3) = (r2 - z2)*drho2 + (delta*z2)/(r2*rho2)
 endif

end subroutine get_metric_cartesian

!--- The metric tensor in SPHERICAL-like form
pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
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

 if (present(gcon)) then
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
 endif

 if (present(sqrtg)) sqrtg = grr*gthetatheta*(-gtphi**2+gphiphi*gtt)

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
 real :: rs,x,y,z,x2,y2,z2,a2,r2spherical,r2,r,rho2
 real :: a2pr2,delta,sintheta2,gtt,gphiphi,gtphi
 real :: r21,rho21,r2mz2,r2mz21,x2py2,x2py21,delta1
 real :: term1,term2
 real :: dr2dx,dr2dy,dr2dz
 real :: drho2dx,drho2dy,drho2dz
 real :: ddeltadx,ddeltady,ddeltadz
 real :: dgtphidx,dgtphidy,dgtphidz
 real :: dsintheta2dx,dsintheta2dy,dsintheta2dz
 real :: dgphiphidx,dgphiphidy,dgphiphidz

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

 a2pr2     = a2+r2
 delta     = a2pr2 - rs*r
 delta1    = 1./delta
 sintheta2 = 1. - z2/r2
 gtt       = -1. + (r*rs)/rho2
 gphiphi   = sintheta2*(a2 + r2 + (a2*r*rs*sintheta2)/rho2)
 gtphi     = -((a*r*rs*sintheta2)/rho2)

 r21    = 1./r2
 rho21  = 1./rho2
 r2mz2  = r2 - z2
 r2mz21 = 1./r2mz2
 x2py2  = x2 + y2
 x2py21 = 1./x2py2

 term1 = a2pr2*rho21
 term2 = (2*a2pr2**2*rho21**2*z**3)/r2**2
 dr2dx = (2*r**2*x)/rho2
 dr2dy = (2*r**2*y)/rho2
 dr2dz = (2*a2pr2*z)/rho2
 drho2dx = 2*r21*rho21*x*(r2**2 - a2*z2)
 drho2dy = 2*r21*rho21*y*(r2**2 - a2*z2)
 drho2dz = 2*r21**2*z*(a2*r2 + a2pr2*rho21*(r2**2 - a2*z2))
 ddeltadx = r*rho21*(2*r - rs)*x
 ddeltady = r*rho21*(2*r - rs)*y
 ddeltadz = (a2pr2*rho21*(2*r - rs)*z)/r
 dgtphidx = (a*rho21**2*rs*(Drho2dx*r2mz2 - x*(r2 + z2)))/r
 dgtphidy = (a*rho21**2*rs*(Drho2dy*r2mz2 - y*(r2 + z2)))/r
 dgtphidz = a*Drho2dz*r*rho21**2*rs*sintheta2 - &
            (a*a2pr2*rho21**2*rs*sintheta2*z)/r - a*r*rho21*rs*(-2*r21*z + 2*a2pr2*r21**2*rho21*z**3)
 dsintheta2dx = 2*r21*rho21*x*z2
 dsintheta2dy = 2*r21*rho21*y*z2
 dsintheta2dz = -2*r21*z + 2*a2pr2*r21**2*rho21*z**3
 dgphiphidx = dsintheta2dx*r2 + Dr2dx*sintheta2 + a2*rho21**2*(dsintheta2dx*rho2**2 - drho2dx*r*rs*sintheta2**2 + &
            r*rs*sintheta2*(2*dsintheta2dx*rho2 + sintheta2*x))
 dgphiphidy = dsintheta2dy*r2 + Dr2dy*sintheta2 + a2*rho21**2*(dsintheta2dy*rho2**2 - drho2dy*r*rs*sintheta2**2 + &
            r*rs*sintheta2*(2*dsintheta2dy*rho2 + sintheta2*y))
 dgphiphidz = dsintheta2dz*r2 + Dr2dz*sintheta2 + (a2**2*rho21**2*rs*sintheta2**2*z)/r + a2*rho21**2*(dsintheta2dz*rho2**2 - &
            drho2dz*r*rs*sintheta2**2 + r*rs*sintheta2*(2*dsintheta2dz*rho2 + sintheta2*z))

 dgcovdx=0.
 dgcovdy=0.
 dgcovdz=0.

 dgcovdx(0,0) = r*rho21**2*rs*(-Drho2dx + x)
 dgcovdx(0,1) = x2py21*(-Dgtphidx + 2*gtphi*x*x2py21)*y
 dgcovdx(1,1) = -(Ddeltadx*delta1**2*r2*rho21*x2) + delta1*r2*rho21*(2*x + 2*rho21*x**3 - Drho2dx*rho21*x2) + &
               x2py21**2*(Dgphiphidx - 4*gphiphi*x*x2py21)*y2 - r2mz21*rho21*(-2*x + 2*r2*r2mz21*rho21*x**3 + Drho2dx*rho21*x2)*z2
 dgcovdx(0,2) = x2py21*(gtphi + Dgtphidx*x - 2*gtphi*x2*x2py21)
 dgcovdx(1,2) = -(y*(Ddeltadx*delta1**2*r2*rho21*x + delta1*r2*rho21*(-1 + Drho2dx*rho21*x - 2*rho21*x2) + x2py21**2*(gphiphi + &
               Dgphiphidx*x - 4*gphiphi*x2*x2py21) + r2mz21*rho21*(-1 + Drho2dx*rho21*x + 2*r2*r2mz21*rho21*x2)*z2))
 dgcovdx(2,2) = Dgphiphidx*x2*x2py21**2 + 2*gphiphi*x*x2py21**2*(1 - 2*x**2*x2py21) - rho21*y2*(delta1*r2*(Ddeltadx*delta1 + &
               rho21*(Drho2dx - 2*x)) + r2mz21*rho21*(Drho2dx + 2*r2*r2mz21*x)*z2)
 dgcovdx(0,3) = 0
 dgcovdx(1,3) = -(z*(r2mz21 + delta1*term1*(-1 + Ddeltadx*delta1*x + Drho2dx*rho21*x) - 2*delta1*r2*rho21**2*x2 + &
               r2mz21*(-2*rho21**2*x**2 + r21*term1*(-1 + rho21*x*(Drho2dx + 2*x)))*z2 + 2*r2*r2mz21**2*rho21*x2*(-1 + &
               r21*term1*z2)))
 dgcovdx(2,3) = -(y*z*(Ddeltadx*delta1**2*term1 + delta1*rho21*(Drho2dx*term1 - 2*r2*rho21*x) + r2mz21*rho21*(-2*r2*r2mz21*x - &
               2*rho21*x*z2 + r21*term1*(Drho2dx + 2*(x + r2*r2mz21*x))*z2)))
 dgcovdx(3,3) = -(a2pr2**2*delta1*Drho2dx*r21*rho21**2*z2) - a2pr2**2*delta1*r21*rho21*(Ddeltadx*delta1 + 2*rho21*x)*z2 - &
               2*r2*r2mz21**2*x*(-1 + r21*term1*z2)**2 - Drho2dx*r2mz21*(-1 + r21*term1*z2)*(1 + r21*(-1 + &
               2*rho2*rho21)*term1*z2) + 4*rho21*x*z2*(delta1*term1 + r2mz21*rho2*(rho21 - r21*term1)*(-1 + r21*term1*z2))

 dgcovdy(0,0) = r*rho21**2*rs*(-Drho2dy + y)
 dgcovdy(0,1) = x2py21*(-(Dgtphidy*y) + gtphi*(-1 + 2*x2py21*y2))
 dgcovdy(1,1) = -(Ddeltady*delta1**2*r2*rho21*x2) - delta1*r2*rho21**2*x2*(Drho2dy - 2*y) + x2py21**2*(2*gphiphi*y - &
               4*gphiphi*x2py21*y**3 + Dgphiphidy*y2) - r2mz21*rho21**2*x2*(Drho2dy + 2*r2*r2mz21*y)*z2
 dgcovdy(0,2) = x*x2py21*(Dgtphidy - 2*gtphi*x2py21*y)
 dgcovdy(1,2) = -(x*(Ddeltady*delta1**2*r2*rho21*y + delta1*r2*rho21*(-1 + Drho2dy*rho21*y - 2*rho21*y2) + x2py21**2*(gphiphi + &
               Dgphiphidy*y - 4*gphiphi*x2py21*y2) + r2mz21*rho21*(-1 + Drho2dy*rho21*y + 2*r2*r2mz21*rho21*y2)*z2))
 dgcovdy(2,2) = Dgphiphidy*x2*x2py21**2 + 2*delta1*r2*rho21*y - 4*gphiphi*x2*x2py21**3*y - delta1*r2*rho21*(-2*rho21*y**3 + &
               Ddeltady*delta1*y2 + Drho2dy*rho21*y2) - r2mz21*rho21*(-2*y + 2*r2*r2mz21*rho21*y**3 + Drho2dy*rho21*y2)*z2
 dgcovdy(0,3) = 0
 dgcovdy(1,3) = -(x*z*(Ddeltady*delta1**2*term1 + delta1*rho21*(Drho2dy*term1 - 2*r2*rho21*y) + r2mz21*rho21*(-2*r2*r2mz21*y - &
               2*rho21*y*z2 + r21*term1*(Drho2dy + 2*(y + r2*r2mz21*y))*z2)))
 dgcovdy(2,3) = -(z*(r2mz21 + delta1*term1*(-1 + Ddeltady*delta1*y + Drho2dy*rho21*y) - 2*delta1*r2*rho21**2*y2 + &
               r2mz21*(-2*rho21**2*y**2 + r21*term1*(-1 + rho21*y*(Drho2dy + 2*y)))*z2 + 2*r2*r2mz21**2*rho21*y2*(-1 + &
               r21*term1*z2)))
 dgcovdy(3,3) = -(a2pr2**2*delta1*Drho2dy*r21*rho21**2*z2) - a2pr2**2*delta1*r21*rho21*(Ddeltady*delta1 + 2*rho21*y)*z2 - &
               2*r2*r2mz21**2*y*(-1 + r21*term1*z2)**2 - Drho2dy*r2mz21*(-1 + r21*term1*z2)*(1 + r21*(-1 + &
               2*rho2*rho21)*term1*z2) + 4*rho21*y*z2*(delta1*term1 + r2mz21*rho2*(rho21 - r21*term1)*(-1 + r21*term1*z2))

 dgcovdz(0,0) = -(Drho2dz*r*rho21**2*rs) + (rho21*rs*term1*z)/r
 dgcovdz(0,1) = -(Dgtphidz*x2py21*y)
 dgcovdz(1,1) = -(Ddeltadz*delta1**2*r2*rho21*x2) - delta1*Drho2dz*r2*rho21**2*x2 + Dgphiphidz*x2py21**2*y2 + &
               2*r2mz21*rho21*x2*z + 2*delta1*rho21*term1*x2*z - Drho2dz*r2mz21*rho21**2*x2*z2 - &
               r2mz21**2*rho21*x2*(-2*z + 2*term1*z)*z2
 dgcovdz(0,2) = Dgtphidz*x*x2py21
 dgcovdz(1,2) = -(Ddeltadz*delta1**2*r2*rho21*x*y) - delta1*Drho2dz*r2*rho21**2*x*y - Dgphiphidz*x*x2py21**2*y + &
               2*r2mz21*rho21*x*y*z + 2*delta1*rho21*term1*x*y*z - Drho2dz*r2mz21*rho21**2*x*y*z2 - &
               r2mz21**2*rho21*x*y*(-2*z + 2*term1*z)*z2
 dgcovdz(2,2) = Dgphiphidz*x2*x2py21**2 - Ddeltadz*delta1**2*r2*rho21*y2 - delta1*Drho2dz*r2*rho21**2*y2 + 2*r2mz21*rho21*y2*z + &
               2*delta1*rho21*term1*y2*z - Drho2dz*r2mz21*rho21**2*y2*z2 - r2mz21**2*rho21*y2*(-2*z + 2*term1*z)*z2
 dgcovdz(0,3) = 0
 dgcovdz(1,3) = delta1*term1*x - Ddeltadz*delta1**2*term1*x*z - delta1*Drho2dz*rho21*term1*x*z + 2*delta1*rho21*term1*x*z2 - &
               r2mz21*x*(1 - r21*term1*z2) + r2mz21**2*x*z*(-2*z + 2*term1*z)*(1 - r21*term1*z2) - r2mz21*x*z*(term2 - &
               2*r21*term1*z - 2*r21*rho21*term1*z**3 + Drho2dz*r21*rho21*term1*z2)
 dgcovdz(2,3) = delta1*term1*y - Ddeltadz*delta1**2*term1*y*z - delta1*Drho2dz*rho21*term1*y*z + 2*delta1*rho21*term1*y*z2 - &
               r2mz21*y*(1 - r21*term1*z2) + r2mz21**2*y*z*(-2*z + 2*term1*z)*(1 - r21*term1*z2) - r2mz21*y*z*(term2 - &
               2*r21*term1*z - 2*r21*rho21*term1*z**3 + Drho2dz*r21*rho21*term1*z2)
 dgcovdz(3,3) = 2*a2pr2**2*delta1*r21*rho21*z - (2*a2pr2**3*delta1*rho21**2*z**3)/r2**2 + 4*a2pr2**2*delta1*r21*rho21**2*z**3 - &
               a2pr2**2*Ddeltadz*delta1**2*r21*rho21*z2 - a2pr2**2*delta1*Drho2dz*r21*rho21**2*z2 + &
               Drho2dz*r2mz21*(1 - r21*term1*z2)**2 - r2mz21**2*rho2*(-2*z + 2*term1*z)*(1 - r21*term1*z2)**2 + &
               2*r2mz21*rho2*(1 - r21*term1*z2)*(term2 - 2*r21*term1*z - 2*r21*rho21*term1*z**3 + Drho2dz*r21*rho21*term1*z2)

 dgcovdx(1,0) = dgcovdx(0,1)
 dgcovdx(2,0) = dgcovdx(0,2)
 dgcovdx(2,1) = dgcovdx(1,2)
 dgcovdx(3,0) = dgcovdx(0,3)
 dgcovdx(3,1) = dgcovdx(1,3)
 dgcovdx(3,2) = dgcovdx(2,3)

 dgcovdy(1,0) = dgcovdy(0,1)
 dgcovdy(2,0) = dgcovdy(0,2)
 dgcovdy(2,1) = dgcovdy(1,2)
 dgcovdy(3,0) = dgcovdy(0,3)
 dgcovdy(3,1) = dgcovdy(1,3)
 dgcovdy(3,2) = dgcovdy(2,3)

 dgcovdz(1,0) = dgcovdz(0,1)
 dgcovdz(2,0) = dgcovdz(0,2)
 dgcovdz(2,1) = dgcovdz(1,2)
 dgcovdz(3,0) = dgcovdz(0,3)
 dgcovdz(3,1) = dgcovdz(1,3)
 dgcovdz(3,2) = dgcovdz(2,3)

 ! STOP 'No analytic metric derivatives implemented'
end subroutine metric_cartesian_derivatives

!--- Derivatives of the covariant 'SPHERICAL' metric
pure subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
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
pure subroutine get_jacobian(position,dxdx)
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
pure subroutine cartesian2spherical(xcart,xspher)
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
pure subroutine spherical2cartesian(xspher,xcart)
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
