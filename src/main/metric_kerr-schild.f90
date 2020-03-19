module metric

!------ The Kerr metric in Kerr-Schild coordinates

 implicit none
 character(len=*), parameter :: metric_type = 'Kerr-Schild'
 integer,          parameter :: imetric     = 4

 real, public  :: mass1 = 1.       ! mass of central object
 real, public  :: a     = 0.       ! spin of central object

contains

!----------------------------------------------------------------
!+
!  Metric tensors
!+
!----------------------------------------------------------------

!--- The Kerr-Schild metric tensor in CARTESIAN-like form
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: x,y,z,x2,y2,z2,a2
 real :: r2spherical,r2,r
 real :: rho2,delta,r2a2,term,sintheta2
 real :: gphiphi,gtphi,gtt
 real :: rs
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

end module metric
