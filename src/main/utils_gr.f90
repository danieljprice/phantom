module utils_gr
 implicit none

 public :: dot_product_gr, get_u0, get_bigv, rho2dens, h2dens

 private

contains

!----------------------------------------------------------------
!+
!  Function to perform a dot product in general relativity
!  i.e. g_\mu\nu v^\mu \v^nu
!+
!----------------------------------------------------------------
pure real function dot_product_gr(vec1,vec2,gcov)
 real, intent(in) :: vec1(:)
 real, intent(in) :: vec2(size(vec1))
 real, intent(in) :: gcov(size(vec1),size(vec2))
 integer :: i

 dot_product_gr = 0.
 do i=1,size(vec1)
    dot_product_gr = dot_product_gr + dot_product(gcov(:,i),vec1(i)*vec2(:))
 enddo

 return
end function dot_product_gr

!-------------------------------------------------------------------------------

subroutine get_u0(gcov,v,U0)
 real, intent(in)  :: gcov(0:3,0:3), v(1:3)
 real, intent(out) :: U0
 real :: v4(0:3)

 v4(0) = 1.
 v4(1:3) = v(1:3)
 U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))

end subroutine get_u0

!-------------------------------------------------------------------------------

subroutine get_bigv(x,v,bigv,bigv2,alpha,lorentz)
 use metric_tools, only:get_metric3plus1
 real, intent(in)  :: x(1:3),v(1:3)
 real, intent(out) :: bigv(1:3),bigv2,alpha,lorentz
 real :: betadown(1:3),betaUP(1:3),gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)

 call get_metric3plus1(x,alpha,betadown,betaUP,gammaijdown,gammaijUP)
 bigv = (v + betaUP)/alpha
 bigv2 = dot_product_gr(bigv,bigv,gammaijdown)
 lorentz = 1./sqrt(1.-bigv2)

end subroutine get_bigv

!-------------------------------------------------------------------------------

subroutine h2dens(dens,xyzh,grpacki,v)
 use part, only: rhoh,massoftype,igas
 real, intent(in) :: xyzh(1:4),grpacki(:,:,:),v(1:3)
 real, intent(out):: dens
 real :: rho, h, xyz(1:3)

 xyz = xyzh(1:3)
 h   = xyzh(4)
 rho = rhoh(h,massoftype(igas))
 call rho2dens(dens,rho,xyz,grpacki,v)

end subroutine h2dens

subroutine rho2dens(dens,rho,position,grpacki,v)
 use metric_tools, only:unpack_grpacki
 real, intent(in) :: rho,position(1:3),grpacki(:,:,:),v(1:3)
 real, intent(out):: dens
 real :: gcov(0:3,0:3), sqrtg, U0

 ! Hard coded sqrtg=1 since phantom is always in cartesian coordinates
 sqrtg = 1.
 call unpack_grpacki(grpacki,gcov=gcov)
 call get_u0(gcov,v,U0)
 dens = rho/(sqrtg*U0)

end subroutine rho2dens

! This is not being used at the moment.
! subroutine dens2rho(rho,dens,position,v)
!  use metric_tools, only: get_metric
!  real, intent(in) :: dens,position(1:3),v(1:3)
!  real, intent(out):: rho
!  real :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg, U0
!
!  call get_metric(position,gcov,gcon,sqrtg)
!  call get_u0(gcov,v,U0)
!
!  rho = sqrtg*U0*dens
!
! end subroutine dens2rho

!-------------------------------------------------------------------------------

end module utils_gr
