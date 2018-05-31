module metric_tools

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

 public :: get_metric, get_metric_derivs, get_metric3plus1, print_metricinfo, init_metric, get_grpacki, unpack_grpacki

 interface get_metric3plus1
  module procedure get_metric3plus1_only, get_metric3plus1_both
 end interface get_metric3plus1

 private

contains

!-------------------------------------------------------------------------------

!--- This is a wrapper subroutine to get the metric tensor in both covariant (gcov) and
!    contravariant (gcon) form.
subroutine get_metric(position,gcov,gcon,sqrtg)
 use metric,     only: get_metric_cartesian,get_metric_spherical,cartesian2spherical
 use inverse4x4, only: inv4x4
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 real :: det

 select case(icoordinate)
 case(icoord_cartesian)
    call get_metric_cartesian(position,gcov,gcon,sqrtg)
 case(icoord_spherical)
    call get_metric_spherical(position,gcov,gcon,sqrtg)
 end select

 if (useinv4x4) call inv4x4(gcov,gcon,det)

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
    if (imetric /= imet_kerr) then
       call metric_cartesian_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
    else if (imetric == imet_kerr) then
       call numerical_metric_derivs(position,dgcovdx1, dgcovdx2, dgcovdx3)
    else
       STOP 'No derivatives being used...'
    end if
 case(icoord_spherical)
    call metric_spherical_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
 end select

end subroutine get_metric_derivs

!-------------------------------------------------------------------------------

!--- The numerical derivatives of the covariant metric tensor
subroutine numerical_metric_derivs(position,dgcovdx, dgcovdy, dgcovdz)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdx,dgcovdy,dgcovdz
 real :: gblah(0:3,0:3), temp(3), gplus(0:3,0:3),gminus(0:3,0:3),dx,dy,dz,di,sqrtgblag
 di = 1.e-8
 dx = di
 dy = di
 dz = di

 temp      = position
 temp(1)   = temp(1)+dx
 call get_metric(temp,gplus,gblah,sqrtgblag)
 temp      = position
 temp(1)   = temp(1)-dx
 call get_metric(temp,gminus,gblah,sqrtgblag)
 dgcovdx = 0.5*(gplus-gminus)/dx

 temp      = position
 temp(2)   = temp(2)+dy
 call get_metric(temp,gplus,gblah,sqrtgblag)
 temp      = position
 temp(2)   = temp(2)-dy
 call get_metric(temp,gminus,gblah,sqrtgblag)
 dgcovdy = 0.5*(gplus-gminus)/dy

 temp      = position
 temp(3)   = temp(3)+dz
 call get_metric(temp,gplus,gblah,sqrtgblag)
 temp      = position
 temp(3)   = temp(3)-dz
 call get_metric(temp,gminus,gblah,sqrtgblag)
 dgcovdz = 0.5*(gplus-gminus)/dz
end subroutine numerical_metric_derivs

!-------------------------------------------------------------------------------
subroutine get_metric3plus1_only(x,alpha,betadown,betaUP,gammaijdown,gammaijUP)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha,betadown(1:3),betaUP(1:3),gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real              :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 call metric3p1(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
end subroutine get_metric3plus1_only

subroutine get_metric3plus1_both(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha,betadown(1:3),betaUP(1:3),gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
 call metric3p1(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
end subroutine get_metric3plus1_both

subroutine metric3p1(x,alpha,betadown,betaUP,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
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

subroutine init_metric(npart,xyzh,grpack,metricderivs)
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: grpack(:,:,:,:), metricderivs(:,:,:,:)
 integer :: i

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,grpack,metricderivs) &
 !$omp private(i)
 do i=1,npart
    call get_grpacki(xyzh(1:3,i),grpack(:,:,:,i))
    call pack_metricderivs(xyzh(1:3,i),metricderivs(:,:,:,i))
 enddo
 !omp end parallel do

end subroutine init_metric

!
!--- Subroutine to pack the metric (cov and con) + its derivatives into a single array
!
subroutine get_grpacki(xyz,grpacki)
 real, intent(in)  :: xyz(3)
 real, intent(out) :: grpacki(:,:,:)
 real :: sqrtg

 call get_metric(xyz,gcov=grpacki(:,:,1),gcon=grpacki(:,:,2),sqrtg=sqrtg)
 call get_metric_derivs(xyz,dgcovdx1=grpacki(:,:,3),dgcovdx2=grpacki(:,:,4),dgcovdx3=grpacki(:,:,5))

end subroutine get_grpacki

subroutine pack_metricderivs(xyzi,metricderivsi)
 real, intent(in)  :: xyzi(3)
 real, intent(out) :: metricderivsi(0:3,0:3,3)

 call get_metric_derivs(xyzi,metricderivsi(:,:,1),metricderivsi(:,:,2),metricderivsi(:,:,3))

end subroutine pack_metricderivs

!
!--- Subroutine to return metric/components from grpack
!
subroutine unpack_grpacki(grpacki,gcov,gcon,gammaijdown,gammaijUP,dg1,dg2,dg3,alpha,betadown,betaUP)
 real, intent(in),  dimension(0:3,0:3,5)         :: grpacki
 real, intent(out), dimension(0:3,0:3), optional :: gcov,gcon,dg1,dg2,dg3
 real, intent(out), dimension(1:3,1:3), optional :: gammaijdown,gammaijUP
 real, intent(out),                     optional :: alpha,betadown(3),betaUP(3)
 real, allocatable, dimension(:,:) :: gammaijUP_  ! <-- Only allocate when needed (expensive?)
 real :: alpha_,betaUP_(3)                        ! <-- These are cheap to allocate ?
 integer :: i,j

 if (present(alpha).or.present(betaUP).or.present(gammaijUP)) alpha_  = sqrt(-1./grpacki(0,0,2))
 if (present(betaUP).or.present(gammaijUP))                   betaUP_ = grpacki(0,1:3,2) * (alpha_**2)
 !                                                                     |^ gcon_(0,1:3) ^|
 if (present(gammaijUP)) then
    allocate(gammaijUP_(1:3,1:3))
    gammaijUP_   = 0.
    do j=1,3
      do i=1,3
         gammaijUP_(i,j) = grpacki(i,j,2) + betaUP_(i)*betaUP_(j)/alpha_**2
         !                |^ gcon_(i,j) ^|
      enddo
    enddo
    gammaijUP   = gammaijUP_
 endif

 if (present(gcov))        gcov        = grpacki(:,:,1)
 if (present(gcon))        gcon        = grpacki(:,:,2)
 if (present(gammaijdown)) gammaijdown = grpacki(1:3,1:3,1)
 if (present(alpha))       alpha       = alpha_
 if (present(betadown))    betadown    = grpacki(0,1:3,1)
 if (present(betaUP))      betaUP      = betaUP_
 if (present(dg1))         dg1         = grpacki(:,:,3)
 if (present(dg2))         dg2         = grpacki(:,:,4)
 if (present(dg3))         dg3         = grpacki(:,:,5)

end subroutine unpack_grpacki

end module metric_tools
