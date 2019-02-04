module metric_tools
 use metric, only:imetric

 implicit none

!-------------------------------------------------------------------------------
!
! This module contains wrapper subroutines to get:
!      - The metric (covariant and contravariant)
!      - Derivatives of the covariant metric
! As well as some general tools that are not specfic to each metric:
!      - Numerical metric derivatives
!      - Tensor transformations
!
!-------------------------------------------------------------------------------

!--- List of coordinates
 integer, public, parameter :: &
    icoord_cartesian  = 1,     &    ! Cartesian coordinates
    icoord_spherical  = 2           ! Spherical coordinates

!--- List of metrics
 integer, public, parameter :: &
    imet_minkowski      = 1,   &    ! Minkowski metric
    imet_schwarzschild  = 2,   &    ! Schwarzschild metric
    imet_kerr           = 3         ! Kerr metric

!--- Choice of coordinate system
!    (When using this with PHANTOM, it should always be set to cartesian)
 integer, public, parameter :: icoordinate = icoord_cartesian

!--- Choice for contravariant metric
!    false  ->  use analytic contravariant metric
!    true   ->  invert the covariant metric
 logical, private, parameter :: useinv4x4 = .true.

!-------------------------------------------------------------------------------

 public :: get_metric, get_metric_derivs, get_metric3plus1, print_metricinfo, init_metric, pack_metric, unpack_metric
 public :: pack_metricderivs
 public :: imetric

 interface get_metric3plus1
  module procedure get_metric3plus1_only, get_metric3plus1_both
 end interface get_metric3plus1

 private

contains

!-------------------------------------------------------------------------------

!--- This is a wrapper subroutine to get the metric tensor in both covariant (gcov) and
!    contravariant (gcon) form.
pure subroutine get_metric(position,gcov,gcon,sqrtg)
 use metric,     only: get_metric_cartesian,get_metric_spherical,cartesian2spherical
 use inverse4x4, only: inv4x4
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 real :: det

 select case(icoordinate)
 case(icoord_cartesian)
    if (useinv4x4) then
       call get_metric_cartesian(position,gcov,sqrtg=sqrtg)
       call inv4x4(gcov,gcon,det)
    else
       call get_metric_cartesian(position,gcov,gcon=gcon,sqrtg=sqrtg)
    endif
 case(icoord_spherical)
    if (useinv4x4) then
       call get_metric_spherical(position,gcov,sqrtg=sqrtg)
       call inv4x4(gcov,gcon,det)
    else
       call get_metric_spherical(position,gcov,gcon=gcon,sqrtg=sqrtg)
    endif
 end select

end subroutine get_metric

! This is a wrapper subroutine to get the derivatives of the covariant metric tensor.
! The actual analytic metric derivaties are in the metric module, which are different for each type
! of metric.
subroutine get_metric_derivs(position,dgcovdx1, dgcovdx2, dgcovdx3)
 use metric, only: metric_cartesian_derivatives, metric_spherical_derivatives, imetric
 real, intent(in)  :: position(3)
 real, intent(out) :: dgcovdx1(0:3,0:3), dgcovdx2(0:3,0:3), dgcovdx3(0:3,0:3)

 select case(icoordinate)

 case(icoord_cartesian)
    call metric_cartesian_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
 case(icoord_spherical)
    call metric_spherical_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
 end select

end subroutine get_metric_derivs

!-------------------------------------------------------------------------------
! (not being used at the moment)
!--- The numerical derivatives of the covariant metric tensor
! pure subroutine numerical_metric_derivs(position,dgcovdx, dgcovdy, dgcovdz)
!  real, intent(in) :: position(3)
!  real, intent(out), dimension(0:3,0:3) :: dgcovdx,dgcovdy,dgcovdz
!  real :: gblah(0:3,0:3), temp(3), gplus(0:3,0:3),gminus(0:3,0:3),dx,dy,dz,di,sqrtgblag
!  di = 1.e-8
!  dx = di
!  dy = di
!  dz = di
!
!  temp      = position
!  temp(1)   = temp(1)+dx
!  call get_metric(temp,gplus,gblah,sqrtgblag)
!  temp      = position
!  temp(1)   = temp(1)-dx
!  call get_metric(temp,gminus,gblah,sqrtgblag)
!  dgcovdx = 0.5*(gplus-gminus)/dx
!
!  temp      = position
!  temp(2)   = temp(2)+dy
!  call get_metric(temp,gplus,gblah,sqrtgblag)
!  temp      = position
!  temp(2)   = temp(2)-dy
!  call get_metric(temp,gminus,gblah,sqrtgblag)
!  dgcovdy = 0.5*(gplus-gminus)/dy
!
!  temp      = position
!  temp(3)   = temp(3)+dz
!  call get_metric(temp,gplus,gblah,sqrtgblag)
!  temp      = position
!  temp(3)   = temp(3)-dz
!  call get_metric(temp,gminus,gblah,sqrtgblag)
!  dgcovdz = 0.5*(gplus-gminus)/dz
! end subroutine numerical_metric_derivs

!-------------------------------------------------------------------------------
pure subroutine get_metric3plus1_only(x,alpha,betadown,betaUP,gammaijdown,gammaijUP)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha,betadown(1:3),betaUP(1:3),gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real              :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 call metric3p1(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
end subroutine get_metric3plus1_only

pure subroutine get_metric3plus1_both(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha,betadown(1:3),betaUP(1:3),gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
 call metric3p1(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
end subroutine get_metric3plus1_both

pure subroutine metric3p1(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha,betadown(1:3),betaUP(1:3), gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
 integer :: i,j

 call get_metric(x,gcov,gcon,sqrtg)
 betadown    = gcov(0,1:3)
 gammaijdown = gcov(1:3,1:3)
 alpha       = sqrt(-1./gcon(0,0))
 betaUP      = gcon(0,1:3)*alpha**2
 gammaijUP   = 0.
 do i=1,3
    do j=1,3
       gammaijUP(i,j) = gcon(i,j) + betaUP(i)*betaUP(j)/alpha**2
    enddo
 enddo
end subroutine metric3p1
!-------------------------------------------------------------------------------

! This is not being used at the moment...
!-- Do a coordinate transformation of a 4x4 rank-2 tensor with both indices down
! subroutine tensortransform_dd(position,T_old,T_new)
!  use metric, only: get_jacobian
!  real, intent(in), dimension(3) :: position
!  real, intent(in), dimension(0:3,0:3) :: T_old
!  real, intent(out), dimension(0:3,0:3) :: T_new
!  real, dimension(0:3,0:3) :: dxdx
!  integer :: i,j,k,l
!  call get_jacobian(position,dxdx)
!  T_new = 0.
!  do i=0,3
!     do j=0,3
!        do k=0,3
!           do l=0,3
!              T_new(i,j) = T_new(i,j)+dxdx(k,i)*dxdx(l,j)*T_old(k,l)
!           enddo
!        enddo
!     enddo
!  enddo
! end subroutine tensortransform_dd

subroutine print_metricinfo(iprint)
 use metric, only:metric_type
 integer, intent(in) :: iprint

 write(iprint,*) 'Metric = ',trim(metric_type)

end subroutine print_metricinfo

subroutine init_metric(npart,xyzh,metrics,metricderivs)
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: metrics(:,:,:,:), metricderivs(:,:,:,:)
 integer :: i

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,metrics,metricderivs) &
 !$omp private(i)
 do i=1,npart
    call pack_metric(xyzh(1:3,i),metrics(:,:,:,i))
    call pack_metricderivs(xyzh(1:3,i),metricderivs(:,:,:,i))
 enddo
 !omp end parallel do

end subroutine init_metric

!
!--- Subroutine to pack the metric (cov and con) into a single array
!
pure subroutine pack_metric(xyz,metrici)
 real, intent(in)  :: xyz(3)
 real, intent(out) :: metrici(:,:,:)
 real :: sqrtg

 call get_metric(xyz,gcov=metrici(:,:,1),gcon=metrici(:,:,2),sqrtg=sqrtg)

end subroutine pack_metric

subroutine pack_metricderivs(xyzi,metricderivsi)
 real, intent(in)  :: xyzi(3)
 real, intent(out) :: metricderivsi(0:3,0:3,3)

 call get_metric_derivs(xyzi,metricderivsi(:,:,1),metricderivsi(:,:,2),metricderivsi(:,:,3))

end subroutine pack_metricderivs

!
!--- Subroutine to return metric/components from metrici array
!
pure subroutine unpack_metric(metrici,gcov,gcon,gammaijdown,gammaijUP,alpha,betadown,betaUP)
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 real, intent(in), dimension(0:3,0:3,2) :: metrici
 real, intent(out), dimension(0:3,0:3), optional :: gcov,gcon
 real, intent(out), dimension(1:3,1:3), optional :: gammaijdown,gammaijUP
 real, intent(out),                     optional :: alpha,betadown(3),betaUP(3)
 integer :: i,j

#ifdef FINVSQRT
 if (present(alpha)) alpha  = finvsqrt(-metrici(0,0,2))
#else
 if (present(alpha)) alpha  = sqrt(-1./metrici(0,0,2))
#endif

 if (present(betaUP)) betaUP = metrici(0,1:3,2) * (-1./metrici(0,0,2)) ! = gcon(0,1:3)*alpha**2

 if (present(gammaijUP)) then
    gammaijUP = 0.
    do j=1,3
      do i=1,3
         gammaijUP(i,j) = metrici(i,j,2) + metrici(0,i,2)*metrici(0,j,2)*(-1./metrici(0,0,2)) ! = gcon(i,j) + betaUP(i)*betaUP(j)/alpha**2
      enddo
    enddo
 endif

 if (present(gcov))        gcov        = metrici(0:3,0:3,1)
 if (present(gcon))        gcon        = metrici(0:3,0:3,2)
 if (present(gammaijdown)) gammaijdown = metrici(1:3,1:3,1)
 if (present(betadown))    betadown    = metrici(0,1:3,1)

end subroutine unpack_metric

end module metric_tools
