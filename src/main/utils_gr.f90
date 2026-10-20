!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_gr
!
! Utility routines for the GR code
!
! :References: Liptai & Price (2019), MNRAS 485, 819-842
!
! :Owner: Spencer Magnall
!
! :Runtime parameters: None
!
! :Dependencies: io, metric, metric_tools
!
 implicit none

 public :: dot_product_gr, get_u0, get_bigv, rho2dens, get_geodesic_accel, get_sqrtg, get_sqrt_gamma
 public :: perturb_metric,get_u_los,get_coord_accel,get_b0,get_b_lower

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

end function dot_product_gr

!----------------------------------------------------------------
!+
!  Function to return U^0, the time component of the 4-velocity
!+
!----------------------------------------------------------------
pure subroutine get_u0(gcov,v,U0,ierr)
 real,    intent(in)  :: gcov(0:3,0:3), v(1:3)
 real,    intent(out) :: U0
 integer, intent(out) :: ierr
 real :: v4(0:3),vv

 ierr    = 0
 v4(0)   = 1.
 v4(1:3) = v(1:3)
 vv      = dot_product_gr(v4,v4,gcov)
 if (vv < 0.) then
    U0 = 1./sqrt(-vv)
 else
    U0   = huge(U0)  ! velocity faster than light, avoid floating point exception
    ierr = 1
 endif

end subroutine get_u0

!----------------------------------------------------------------
!+
!  Covariant 4-velocity projected along the unit vector rhat,
!  U_i rhat^i with U_i = U^0 (g_{0i} + g_{ij} v^j), for the
!  coordinate velocity v = projv*rhat. Used for the velocity
!  jump in the shock viscosity; for a fixed metric this is a
!  monotonic function of projv.
!
!  Near a horizon, v = projv*rhat can lie outside the light cone
!  even when the full velocity does not. U^0 is therefore limited
!  to u0max by bounding g_mu nu v^mu v^nu <= -1/u0max^2. This
!  keeps U_los finite and monotonic in projv, and does not change
!  states with U^0 < u0max
!+
!----------------------------------------------------------------
pure real function get_u_los(gcov,projv,rhat,u0max)
 real, intent(in) :: gcov(0:3,0:3),projv,rhat(3),u0max
 real :: v(3),v4(0:3),vv,U0

 v = projv*rhat
 v4(0)   = 1.
 v4(1:3) = v
 vv = min(dot_product_gr(v4,v4,gcov),-1./u0max**2)
 U0 = 1./sqrt(-vv)
 get_u_los = U0*dot_product(gcov(0,1:3) + matmul(gcov(1:3,1:3),v),rhat)

end function get_u_los

!----------------------------------------------------------------
!+
!  Coordinate acceleration a^i = dv^i/dt of a particle, given
!  the rate of change of its canonical momentum dp_i/dt, where
!  p_i = w U_i, w is the specific enthalpy and
!  U_i = U^0 (g_{0i} + g_{ij} v^j). By the chain rule
!
!   dp_i/dt = w U^0 (g_ik + U_i U_k) a^k + w v^k d_k U_i
!
!  neglecting dw/dt. The matrix is inverted analytically using
!  the Sherman-Morrison formula with the inverse spatial metric
!  gamma^ij. Returns a = 0 if v is superluminal
!+
!----------------------------------------------------------------
pure subroutine get_coord_accel(metrici,metricderivsi,v,enth,dpdt,a)
 use metric_tools, only:unpack_metric
 real, intent(in)  :: metrici(0:3,0:3,2),metricderivsi(0:3,0:3,3)
 real, intent(in)  :: v(3),enth,dpdt(3)
 real, intent(out) :: a(3)
 real    :: gcov(0:3,0:3),gammaijUP(3,3),dgdt(0:3,0:3),v4(0:3),dgv4(0:3)
 real    :: U0,Ulow(3),Uup(3),rhs(3)
 integer :: k,ierr

 a = 0.
 v4(0)   = 1.
 v4(1:3) = v
 call unpack_metric(metrici,gcov=gcov,gammaijUP=gammaijUP)
 call get_u0(gcov,v,U0,ierr)
 if (ierr /= 0) return
 Ulow = U0*matmul(gcov(1:3,:),v4)

 ! rate of change of the metric along the particle path, v^k d_k g_mu nu
 dgdt = 0.
 do k=1,3
    dgdt = dgdt + v(k)*metricderivsi(:,:,k)
 enddo
 dgv4 = matmul(dgdt,v4)

 ! subtract w v^k d_k U_i, using d_k U^0 = 0.5 (U^0)^3 d_k g_mu nu v^mu v^nu
 rhs = dpdt - enth*U0*(dgv4(1:3) + 0.5*U0*dot_product(v4,dgv4)*Ulow)

 ! solve w U^0 (g_ij + U_i U_j) a^j = rhs_i
 Uup = matmul(gammaijUP,Ulow)
 a = (matmul(gammaijUP,rhs) - Uup*dot_product(Uup,rhs)/(1. + dot_product(Ulow,Uup)))/(enth*U0)

end subroutine get_coord_accel

!----------------------------------------------------------------
!+
!  Function to return V^i, the velocity of an Eulerian observer
!+
!----------------------------------------------------------------
subroutine get_bigv(metrici,v,bigv,bigv2,alpha,lorentz)
 use metric_tools, only:unpack_metric
 use io,           only:fatal
 real, intent(in)  :: metrici(0:3,0:3,2),v(1:3)
 real, intent(out) :: bigv(1:3),bigv2,alpha,lorentz
 real :: betaUP(1:3),gammaijdown(1:3,1:3)

 call unpack_metric(metrici,betaUP=betaUP,alpha=alpha,gammaijdown=gammaijdown)
 bigv = (v + betaUP)/alpha
 bigv2 = dot_product_gr(bigv,bigv,gammaijdown)
 if (bigv2 > 1.) call fatal('get_bigv','velocity faster than speed of light -- bigv2',val=bigv2)
 lorentz = 1./sqrt(1.-bigv2)

end subroutine get_bigv

!----------------------------------------------------------------
!+
!  get primitive rest-frame density (dens) from conserved rho
!  (kernel-summed mass density in Phantom)
!+
!----------------------------------------------------------------
subroutine rho2dens(dens,rho,metrici,v)
 use metric_tools, only:unpack_metric
 use io,           only:error
 real, intent(in)  :: rho,metrici(:,:,:),v(1:3)
 real, intent(out) :: dens
 integer :: ierror
 real :: gcov(0:3,0:3), sqrtg, U0

 call unpack_metric(metrici,gcov=gcov)
 call get_sqrtg(gcov, sqrtg)
 call get_u0(gcov,v,U0,ierror)
 dens = rho/(sqrtg*U0)

 if (ierror > 0) call error('get_u0 in rho2dens','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')

end subroutine rho2dens

!----------------------------------------------------------------
!+
!  get terms required on the RHS of the geodesic equation
!  in the form dp_i/dt = a_i, as described in Liptai & Price (2019)
!+
!----------------------------------------------------------------
subroutine get_geodesic_accel(axyz,npart,vxyz,metrics,metricderivs)
 use metric_tools, only:unpack_metric
 integer, intent(in)  :: npart
 real,    intent(in)  :: vxyz(:,:), metrics(:,:,:,:), metricderivs(:,:,:,:)
 real,    intent(out) :: axyz(3,npart)
 real :: gcon(0:3,0:3), v(0:3), gderiv(0:3,0:3,0:3), a(3)
 integer :: i,lambda,mu,sigma

 !$omp parallel do default(none) &
 !$omp shared(npart,vxyz,axyz,metricderivs,metrics) &
 !$omp private(i,a,v,gderiv,gcon,lambda,mu,sigma)
 do i=1,npart
    v = (/1.,vxyz(:,i)/)
    gderiv = 0.
    gderiv(1:3,1:3,1:3) = metricderivs(:,:,:,i)
    call unpack_metric(metrics(:,:,:,i),gcon=gcon)
    a = 0.
    do lambda=0,3
       do mu=0,3
          do sigma=0,3
             a(:) = a(:) + (gcon(1:3,lambda) - v(1:3)*gcon(0,lambda)) * &
                        (gderiv(mu,lambda,sigma) - 0.5*gderiv(mu,sigma,lambda))*v(mu)*v(sigma)
          enddo
       enddo
    enddo
    axyz(1:3,i) = a
 enddo
 !omp end parallel do

end subroutine get_geodesic_accel

!----------------------------------------------------------------
!+
!  get determininant of the 4-metric
!+
!----------------------------------------------------------------
subroutine get_sqrtg(gcov, sqrtg)
 use metric_tools, only:imetric,imet_binarybh,imet_et
 real, intent(in)  :: gcov(0:3,0:3)
 real, intent(out) :: sqrtg
 real :: det
 real :: a11,a12,a13,a14
 real :: a21,a22,a23,a24
 real :: a31,a32,a33,a34
 real :: a41,a42,a43,a44

 if (imetric == imet_et .or. imetric == imet_binarybh) then

    a11 = gcov(0,0)
    a21 = gcov(1,0)
    a31 = gcov(2,0)
    a41 = gcov(3,0)
    a12 = gcov(0,1)
    a22 = gcov(1,1)
    a32 = gcov(2,1)
    a42 = gcov(3,1)
    a13 = gcov(0,2)
    a23 = gcov(1,2)
    a33 = gcov(2,2)
    a43 = gcov(3,2)
    a14 = gcov(0,3)
    a24 = gcov(1,3)
    a34 = gcov(2,3)
    a44 = gcov(3,3)

    ! Calculate the determinant
    det = a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22*a33*a41 + a12*a24*a33*a41 + &
       a13*a22*a34*a41 - a12*a23*a34*a41 - a14*a23*a31*a42 + a13*a24*a31*a42 + &
       a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34*a42 + a11*a23*a34*a42 + &
       a14*a22*a31*a43 - a12*a24*a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 + &
       a12*a21*a34*a43 - a11*a22*a34*a43 - a13*a22*a31*a44 + a12*a23*a31*a44 + &
       a13*a21*a32*a44 - a11*a23*a32*a44 - a12*a21*a33*a44 + a11*a22*a33*a44

    sqrtg = sqrt(-det)
 else
    ! If we are not using an evolving metric then
    ! Sqrtg = 1
    sqrtg = 1.
 endif

end subroutine get_sqrtg

!----------------------------------------------------------------
!+
!  get determininant of the 3-metric
!+
!----------------------------------------------------------------
subroutine get_sqrt_gamma(gcov,sqrt_gamma)
 use metric, only:metric_type
 real, intent(in)  :: gcov(0:3,0:3)
 real, intent(out) :: sqrt_gamma
 real :: a11,a12,a13
 real :: a21,a22,a23
 real :: a31,a32,a33
 !real :: a41,a42,a43
 real :: det

 if (metric_type == 'et') then
    ! Calculate the determinant of a 3x3 matrix
    ! Spatial metric is just the physical metric
    ! without the tt component

    a11 = gcov(1,1)
    a12 = gcov(1,2)
    a13 = gcov(1,3)
    a21 = gcov(2,1)
    a22 = gcov(2,2)
    a23 = gcov(2,3)
    a31 = gcov(3,1)
    a32 = gcov(3,2)
    a33 = gcov(3,3)

    det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32-a22*a31)
    sqrt_gamma = sqrt(det)

 else
    sqrt_gamma = 1.
 endif

end subroutine get_sqrt_gamma

!----------------------------------------------------------------
!+
!  add a Newtonian gravitational perturbation to the metric
!+
!----------------------------------------------------------------
subroutine perturb_metric(phi,gcovper,gcov)
 real, intent(in)  :: phi
 real, intent(out) :: gcovper(0:3,0:3)
 real, intent(in), optional :: gcov(0:3,0:3)

 if (present(gcov)) then
    gcovper = gcov
 else
    gcovper = 0.
    gcovper(0,0) = -1.
    gcovper(1,1) = 1.
    gcovper(2,2) = 1.
    gcovper(3,3) = 1.
 endif

 ! Set the perturbed metric based on the Bardeen formulation
 gcovper(0,0) = gcovper(0,0) - 2.*phi
 gcovper(1,1) = gcovper(1,1) - 2.*phi
 gcovper(2,2) = gcovper(2,2) - 2.*phi
 gcovper(3,3) = gcovper(3,3) - 2.*phi

end subroutine perturb_metric

!----------------------------------------------------------------
!+
!  time component of the comoving magnetic field four-vector, b^0,
!  given the spatial components b^i and the coordinate velocity v^i,
!  from the orthogonality condition b^mu U_mu = 0
!+
!----------------------------------------------------------------
pure subroutine get_b0(gcov,v,b,b0)
 real, intent(in)  :: gcov(0:3,0:3),v(1:3),b(1:3)
 real, intent(out) :: b0
 real :: vdown(0:3)
 integer :: mu

 ! v_mu = g_mu nu v^nu with v^0 = 1
 do mu=0,3
    vdown(mu) = gcov(mu,0) + dot_product(gcov(mu,1:3),v)
 enddo
 b0 = -dot_product(b,vdown(1:3))/vdown(0)

end subroutine get_b0

!----------------------------------------------------------------
!+
!  lower the comoving magnetic field four-vector, given
!  bxyz = (b^1,b^2,b^3,b^0), returning the spatial components
!  of b_mu and b^2 = b^mu b_mu
!+
!----------------------------------------------------------------
pure subroutine get_b_lower(gcov,bxyz,bdown,b2)
 real, intent(in)  :: gcov(0:3,0:3),bxyz(4)
 real, intent(out) :: bdown(3),b2
 real :: b4(0:3),b4down(0:3)
 integer :: mu

 b4(0)   = bxyz(4)
 b4(1:3) = bxyz(1:3)
 do mu=0,3
    b4down(mu) = dot_product(gcov(mu,:),b4)
 enddo
 b2    = dot_product(b4,b4down)
 bdown = b4down(1:3)

end subroutine get_b_lower

end module utils_gr
