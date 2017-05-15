module metric_tools
implicit none
!
! This module contains wrapper subroutines to get:
!      - The metric (covariant and contravariant)
!      - Derivatives of the covariant metric
! As well as some general tools that are not specfic to each metric:
!      - Numerical metric derivatives
!      - Tensor transformations
!
character(len=*), parameter :: coordinate_sys = 'Cartesian'

contains

! This is a wrapper subroutine to get the metric tensor in both covariant (gcov) and
! contravariant (gcon) form.
subroutine get_metric(position,gcov,gcon,sqrtg)
   use metric, only: get_metric_cartesian,get_metric_spherical,cartesian2spherical
   use inverse4x4, only: inv4x4
   real,    intent(in)  :: position(3)
   real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
   real :: det
   ! real :: gcovs(0:3,0:3), gcons(0:3,0:3), sqrtgs, dxdx(0:3,0:3), det, xsphere(3)
   logical, parameter :: useinv4x4 = .true.

   select case(coordinate_sys)

   case('Cartesian')
      call get_metric_cartesian(position,gcov,gcon,sqrtg)
   case('Spherical')
      call get_metric_spherical(position,gcov,gcon,sqrtg)
   end select

   if (useinv4x4) call inv4x4(gcov,gcon,det)
   ! call cartesian2spherical(position,xsphere)
   ! sqrtg=1.
   ! call inv4x4(gcov,gcon,det)
   ! sqrtg = sqrt(-det)
end subroutine get_metric

! This is a wrapper subroutine to get the derivatives of the covariant metric tensor.
! The actual analytic metric derivaties are in the metric module, which are different for each type
! of metric.
subroutine get_metric_derivs(position,dgcovdx1, dgcovdx2, dgcovdx3)
   use metric, only: metric_cartesian_derivatives, metric_spherical_derivatives, metric_type
   real,    intent(in)  :: position(3)
   real,    intent(out) :: dgcovdx1(0:3,0:3), dgcovdx2(0:3,0:3), dgcovdx3(0:3,0:3)

   select case(coordinate_sys)

   case('Cartesian')
      if (.not. metric_type=='Kerr') then
         call metric_cartesian_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
      else if (metric_type=='Kerr') then
         call numerical_metric_derivs(position,dgcovdx1, dgcovdx2, dgcovdx3)
      else
         STOP 'No derivatives being used...'
      end if
   case('Spherical')
      call metric_spherical_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
   end select

end subroutine get_metric_derivs

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

!-- Do a coordinate transformation of a 4x4 rank-2 tensor with both indices down
subroutine tensortransform_dd(position,T_old,T_new)
   use metric, only: get_jacobian
   real, intent(in), dimension(3) :: position
   real, intent(in), dimension(0:3,0:3) :: T_old
   real, intent(out), dimension(0:3,0:3) :: T_new
   real, dimension(0:3,0:3) :: dxdx
   integer :: i,j,k,l
   call get_jacobian(position,dxdx)
   T_new = 0.
   do i=0,3
      do j=0,3
         do k=0,3
            do l=0,3
               T_new(i,j) = T_new(i,j)+dxdx(k,i)*dxdx(l,j)*T_old(k,l)
            enddo
         enddo
      enddo
   enddo
end subroutine tensortransform_dd

end module metric_tools
