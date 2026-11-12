!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module boundarypart
!
! This module contains subroutines related to boundary particles
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: None
!

 implicit none
 public :: get_boundary_particle_forces,set_boundary_particle_velocity,&
           calculate_boundary_core_CM,calculate_boundary_core_omega,&
           drift_boundary_core,get_core_residual_accel

 !
 ! sum of m du/dt over the core from the last call with du/dt present,
 ! i.e. the rate at which thermal energy is discarded by freezing u
 !
 real, public :: core_dudt = 0.

 !
 ! if .false., get_boundary_particle_forces leaves the forces on the core
 ! untouched. Only for analysis tools that need the raw forces on the core
 ! (e.g. to measure the acceleration discarded by the rigid constraint)
 !
 logical, public :: rigid_core_forces = .true.

 private

contains
!----------------------------------------------------------------
!+
!  Replace the accelerations of boundary particles with the kick of a
!  rigid body subject to the same net force and torque:
!
!    a_i = a_cm + alpha x r_i,    alpha = I^-1 tau
!
!  where r_i = x_i - x_cm, a_cm = <a> is the mean acceleration,
!  tau = <r x (a - a_cm)> the mean specific torque about the core centre
!  of mass and I is the moment of inertia tensor per unit mass.
!+
!---------------------------------------------------------------
subroutine get_boundary_particle_forces(npart,iphase,xyz,fxyzu,dBevol,drad,ddustprop,ddustevol)
 use part,   only:iamboundary,iamtype,iradxi,massoftype,iboundary
 use dim,    only:mhd,do_radiation,use_dust,use_dustgrowth,maxvxyzu
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:)
 real, intent(inout) :: fxyzu(:,:)
 real, intent(inout), optional :: dBevol(:,:),drad(:,:),ddustevol(:,:),ddustprop(:,:)
 integer             :: nboundary,i,ierr
 real, dimension(3)  :: avg_lin_accel,avg_torque,xyz_core_CM,torquei,ri
 real, dimension(3)  :: alpha,alpha_x_r
 real                :: inertia(3,3),inertia_inv(3,3),dudt_sum
 logical             :: has_dudt

 if (.not.rigid_core_forces) return
 has_dudt = (maxvxyzu >= 4 .and. size(fxyzu,1) >= 4)
 nboundary = 0
 avg_lin_accel = 0.
 xyz_core_CM = 0.
 dudt_sum = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,fxyzu,xyz,has_dudt) &
 !$omp private(i) &
 !$omp reduction(+:avg_lin_accel,nboundary,xyz_core_CM,dudt_sum)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       nboundary = nboundary + 1
       avg_lin_accel = avg_lin_accel + fxyzu(1:3,i)
       xyz_core_CM = xyz_core_CM + xyz(1:3,i)
       if (has_dudt) dudt_sum = dudt_sum + fxyzu(4,i)
    endif
 enddo
 !$omp end parallel do
 if (has_dudt) core_dudt = massoftype(iboundary)*dudt_sum

 if (nboundary > 0) then
    xyz_core_CM = xyz_core_CM / real(nboundary)
    avg_lin_accel = avg_lin_accel / real(nboundary)

    ! mean specific torque about the core centre of mass
    avg_torque = 0.
    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyz,fxyzu,avg_lin_accel,xyz_core_CM) &
    !$omp reduction(+:avg_torque) &
    !$omp private(i,torquei)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          call cross_product3D(xyz(1:3,i)-xyz_core_CM, fxyzu(1:3,i)-avg_lin_accel, torquei)
          avg_torque = avg_torque + torquei
       endif
    enddo
    !$omp end parallel do
    avg_torque = avg_torque / real(nboundary)

    ! angular acceleration of the rigid core
    call get_core_inertia(npart,iphase,xyz,xyz_core_CM,nboundary,inertia)
    call invert_core_inertia(inertia,inertia_inv,ierr)
    if (ierr /= 0) then
       alpha = 0.   ! degenerate core (single particle or collinear): no rotation
    else
       alpha = matmul(inertia_inv,avg_torque)
    endif

    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyz,fxyzu,avg_lin_accel,xyz_core_CM,alpha) &
    !$omp private(i,ri,alpha_x_r)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          ri = xyz(1:3,i) - xyz_core_CM
          call cross_product3D(alpha,ri,alpha_x_r)
          fxyzu(1:3,i) = avg_lin_accel + alpha_x_r
          if (size(fxyzu,1) >= 4) fxyzu(4,i) = 0.
       endif
    enddo
    !$omp end parallel do

    if (mhd .and. present(dBevol)) call zero_boundary_derivs(npart,iphase,dBevol)
    if (do_radiation .and. present(drad)) call zero_boundary_derivs(npart,iphase,drad(iradxi:iradxi,:))
    if (use_dust .and. present(ddustevol)) call zero_boundary_derivs(npart,iphase,ddustevol)
    if (use_dustgrowth .and. present(ddustprop)) call zero_boundary_derivs(npart,iphase,ddustprop(1:1,:))
 endif

end subroutine get_boundary_particle_forces

!----------------------------------------------------------------
!+
!  Mean specific inertia tensor of the boundary particles about
!  xyz_CM, I = <r^2 1 - r r^T>
!+
!---------------------------------------------------------------
subroutine get_core_inertia(npart,iphase,xyz,xyz_CM,nboundary,inertia)
 use part, only:iamboundary,iamtype
 integer, intent(in) :: npart,nboundary
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:),xyz_CM(3)
 real, intent(out)   :: inertia(3,3)
 integer :: i,k
 real    :: ri(3),rrt(3,3)

 rrt = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,xyz,xyz_CM) &
 !$omp private(i,k,ri) &
 !$omp reduction(+:rrt)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       ri = xyz(1:3,i) - xyz_CM
       do k=1,3
          rrt(:,k) = rrt(:,k) + ri(:)*ri(k)
       enddo
    endif
 enddo
 !$omp end parallel do
 rrt = rrt / real(nboundary)

 inertia = -rrt
 do k=1,3
    inertia(k,k) = inertia(k,k) + rrt(1,1) + rrt(2,2) + rrt(3,3)
 enddo

end subroutine get_core_inertia

!----------------------------------------------------------------
!+
!  Invert the moment-of-inertia tensor
!+
!---------------------------------------------------------------
subroutine invert_core_inertia(inertia,inertia_inv,ierr)
 use vectorutils, only:matrixinvert3D,cross_product3D
 real,    intent(in)  :: inertia(3,3)
 real,    intent(out) :: inertia_inv(3,3)
 integer, intent(out) :: ierr
 real, parameter :: tol = 1.e-12
 real :: scale,det,row23(3)

 inertia_inv = 0.
 ierr = 1
 scale = (inertia(1,1) + inertia(2,2) + inertia(3,3))/3.
 if (scale <= 0.) return
 call cross_product3D(inertia(2,:),inertia(3,:),row23)
 det = dot_product(inertia(1,:),row23)
 if (abs(det) < tol*scale**3) return
 call matrixinvert3D(inertia,inertia_inv,ierr)

end subroutine invert_core_inertia

!----------------------------------------------------------------
!+
!  Mean specific spin angular momentum of the boundary particles
!  about its centre of mass, L = Σ (r_i x (v_i - v_cm)) / N
!+
!---------------------------------------------------------------
subroutine get_core_spin(npart,iphase,xyz,vxyz,xyz_CM,vxyz_CM,nboundary,spin)
 use part, only:iamboundary,iamtype
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart,nboundary
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:),vxyz(:,:),xyz_CM(3),vxyz_CM(3)
 real, intent(out)   :: spin(3)
 integer :: i
 real    :: li(3)

 spin = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,xyz,vxyz,xyz_CM,vxyz_CM) &
 !$omp private(i,li) &
 !$omp reduction(+:spin)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       call cross_product3D(xyz(1:3,i)-xyz_CM,vxyz(1:3,i)-vxyz_CM,li)
       spin = spin + li
    endif
 enddo
 !$omp end parallel do
 spin = spin / real(nboundary)

end subroutine get_core_spin

!----------------------------------------------------------------
!+
!  Angular velocity of the boundary-particle core about its centre
!  of mass, omega = I^-1 L
!+
!---------------------------------------------------------------
subroutine calculate_boundary_core_omega(npart,iphase,xyz,vxyz,omega)
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:),vxyz(:,:)
 real, intent(out)   :: omega(3)
 integer :: nboundary,ierr
 real    :: xcm(3),vcm(3),spin(3),inertia(3,3),inertia_inv(3,3)

 omega = 0.
 call get_core_state(npart,iphase,xyz,vxyz,xcm,vcm,spin,inertia,nboundary)
 if (nboundary <= 0) return
 call invert_core_inertia(inertia,inertia_inv,ierr)
 if (ierr /= 0) return
 omega = matmul(inertia_inv,spin)

end subroutine calculate_boundary_core_omega

!----------------------------------------------------------------
!+
!  Return the boundary-particle core's centre of mass, centre-of-mass
!  velocity, moment of inertia tensor, and specific spin angular
!  momentum
!+
!---------------------------------------------------------------
subroutine get_core_state(npart,iphase,xyz,vxyz,xcm,vcm,spin,inertia,nboundary)
 integer, intent(in)  :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real,    intent(in)  :: xyz(:,:),vxyz(:,:)
 real,    intent(out) :: xcm(3),vcm(3),spin(3),inertia(3,3)
 integer, intent(out) :: nboundary

 spin = 0.
 inertia = 0.
 call calculate_boundary_core_CM(npart,iphase,xyz,vxyz,xcm,vcm,nboundary)
 if (nboundary <= 0) return
 call get_core_inertia(npart,iphase,xyz,xcm,nboundary,inertia)
 call get_core_spin(npart,iphase,xyz,vxyz,xcm,vcm,nboundary,spin)

end subroutine get_core_state

!----------------------------------------------------------------
!+
!  Rotation matrix about vector phi by angle |phi| using
!  Rodrigues' formula, R = I + sin(theta) K + (1 - cos(theta)) K^2,
!  where K is the cross-product matrix
!+
!---------------------------------------------------------------
pure subroutine get_rotation_matrix(phi,rotmat)
 real, intent(in)  :: phi(3)
 real, intent(out) :: rotmat(3,3)
 real :: theta,k(3),kmat(3,3)
 integer :: j

 rotmat = 0.
 do j=1,3
    rotmat(j,j) = 1.
 enddo
 theta = sqrt(dot_product(phi,phi))
 if (theta <= tiny(theta)) return
 k = phi/theta
 kmat = reshape((/ 0., k(3), -k(2), &
                  -k(3), 0., k(1), &
                   k(2), -k(1), 0. /),(/3,3/))
 rotmat = rotmat + sin(theta)*kmat + (1. - cos(theta))*matmul(kmat,kmat)

end subroutine get_rotation_matrix

!----------------------------------------------------------------
!+
!  Returns new angular velocity vector of a rigid-body core under
!  the action of spin angular momentum L over time dt. The explicit
!  midpoint method is used (second order in dt, exact for a
!  perfectly spherically symmetric body):
!
!    omega_0   = I^-1 L                 (initial angular velocity)
!    R_1/2     = R(omega_0 dt/2)        (rotation by dt/2)
!    omega_1/2 = R_1/2 I^-1 R_1/2^T L   (rotated MoI is R_1/2 * I * R_1/2^T)
!    R         = R(omega_1/2 dt)        (rotate by dt/2)
!+
!---------------------------------------------------------------
subroutine get_core_rotation(dt,spin,inertia,rotmat,omega_new,ierr)
 real,    intent(in)  :: dt,spin(3),inertia(3,3)
 real,    intent(out) :: rotmat(3,3),omega_new(3)
 integer, intent(out) :: ierr
 real :: inertia_inv(3,3),omega(3),rot_half(3,3)

 omega_new = 0.
 call get_rotation_matrix((/0.,0.,0./),rotmat)
 call invert_core_inertia(inertia,inertia_inv,ierr)
 if (ierr /= 0) return

 omega = matmul(inertia_inv,spin)
 call get_rotation_matrix(0.5*dt*omega,rot_half)
 omega = matmul(rot_half,matmul(inertia_inv,matmul(transpose(rot_half),spin)))
 call get_rotation_matrix(dt*omega,rotmat)
 omega_new = matmul(rotmat,matmul(inertia_inv,matmul(transpose(rotmat),spin)))

end subroutine get_core_rotation

!----------------------------------------------------------------
!+
!  Drift the boundary-particle core over a time dt as a rigid body:
!
!    x_i <- x_CM + V dt + R (x_i - x_CM)
!    v_i <- v_CM + omega_new x (x_i - x_CM_new)
!
!  where X, V are the core centre of mass and velocity, R the
!  free-rotor orientation change at fixed spin L (get_core_rotation)
!  and omega_new = I_new^-1 L. Pairwise distances are preserved to
!  round-off, and V and L are unchanged. Any non-rigid part of the
!  incoming core velocities is discarded.
!
!  This replaces the straight-line drift x_i <- x_i + v_i dt for
!  boundary particles, so callers must skip them in their own drift.
!+
!---------------------------------------------------------------
subroutine drift_boundary_core(dt,npart,iphase,xyzh,vxyzu)
 use part,        only:iamboundary,iamtype,npartoftype,iboundary
 use vectorutils, only:cross_product3D
 real,    intent(in)    :: dt
 integer, intent(in)    :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,nboundary,ierr
 real    :: xcm(3),vcm(3),spin(3),inertia(3,3),rotmat(3,3),omega_new(3)
 real    :: xcm_new(3),ri(3),vrot(3)

 if (npartoftype(iboundary) <= 0) return

 call get_core_state(npart,iphase,xyzh,vxyzu,xcm,vcm,spin,inertia,nboundary)
 if (nboundary <= 0) return
 call get_core_rotation(dt,spin,inertia,rotmat,omega_new,ierr)

 if (ierr /= 0) then  ! translate only
    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyzh,vxyzu,dt) &
    !$omp private(i)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) xyzh(1:3,i) = xyzh(1:3,i) + dt*vxyzu(1:3,i)
    enddo
    !$omp end parallel do
 else  ! translation + rotation
    xcm_new = xcm + dt*vcm
    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyzh,vxyzu,xcm,xcm_new,vcm,rotmat,omega_new) &
    !$omp private(i,ri,vrot)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          ri = matmul(rotmat,xyzh(1:3,i) - xcm)
          call cross_product3D(omega_new,ri,vrot)
          xyzh(1:3,i)  = xcm_new + ri
          vxyzu(1:3,i) = vcm + vrot
       endif
    enddo
    !$omp end parallel do
 endif

end subroutine drift_boundary_core

!----------------------------------------------------------------
!+
!  Residual acceleration of the boundary particles, i.e. the part of
!  the raw acceleration that the rigid-body constraint discards:
!
!    da_i = a_i - [ a_cm + alpha x r_i + omega x (omega x r_i) ]
!
!  where r_i = x_i - x_cm, omega = I^-1 L is the core angular velocity
!  and alpha = I^-1 (tau - omega x I omega) is given by Euler's equations,
!  with tau = <r x (a - a_cm)> the mean specific torque. With this choice
!  the residual has no net force and no net torque, <da> = <r x da> = 0,
!  so it contains only deformation (see boundary_core_validity notes).
!
!  accel must be the raw (unprojected) accelerations, e.g. from a derivs
!  call with rigid_core_forces = .false. plus any external/sink forces.
!  On output dela is zero for non-boundary particles.
!+
!---------------------------------------------------------------
subroutine get_core_residual_accel(npart,iphase,xyz,vxyz,accel,dela,xcm,omega,nboundary)
 use part,        only:iamboundary,iamtype
 use vectorutils, only:cross_product3D
 integer, intent(in)  :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real,    intent(in)  :: xyz(:,:),vxyz(:,:),accel(:,:)
 real,    intent(out) :: dela(:,:),xcm(3),omega(3)
 integer, intent(out) :: nboundary
 integer :: i,ierr
 real    :: vcm(3),spin(3),inertia(3,3),inertia_inv(3,3)
 real    :: acm(3),torque(3),torquei(3),ri(3),alpha(3),gyro(3)
 real    :: alpha_x_r(3),omega_x_r(3),centrip(3)

 dela(:,1:npart) = 0.
 omega = 0.
 call get_core_state(npart,iphase,xyz,vxyz,xcm,vcm,spin,inertia,nboundary)
 if (nboundary <= 0) return

 acm = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,accel) &
 !$omp private(i) &
 !$omp reduction(+:acm)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) acm = acm + accel(1:3,i)
 enddo
 !$omp end parallel do
 acm = acm / real(nboundary)

 torque = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,xyz,accel,acm,xcm) &
 !$omp private(i,torquei) &
 !$omp reduction(+:torque)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       call cross_product3D(xyz(1:3,i)-xcm,accel(1:3,i)-acm,torquei)
       torque = torque + torquei
    endif
 enddo
 !$omp end parallel do
 torque = torque / real(nboundary)

 call invert_core_inertia(inertia,inertia_inv,ierr)
 if (ierr /= 0) then
    alpha = 0.   ! degenerate core: translation only
 else
    omega = matmul(inertia_inv,spin)
    call cross_product3D(omega,matmul(inertia,omega),gyro)
    alpha = matmul(inertia_inv,torque - gyro)
 endif

 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,xyz,accel,dela,acm,xcm,alpha,omega) &
 !$omp private(i,ri,alpha_x_r,omega_x_r,centrip)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       ri = xyz(1:3,i) - xcm
       call cross_product3D(alpha,ri,alpha_x_r)
       call cross_product3D(omega,ri,omega_x_r)
       call cross_product3D(omega,omega_x_r,centrip)
       dela(1:3,i) = accel(1:3,i) - (acm + alpha_x_r + centrip)
    endif
 enddo
 !$omp end parallel do

end subroutine get_core_residual_accel

!----------------------------------------------------------------
!+
!  Zero a time-derivative array for boundary particles
!+
!---------------------------------------------------------------
subroutine zero_boundary_derivs(npart,iphase,darr)
 use part, only:iamboundary,iamtype
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(inout) :: darr(:,:)
 integer :: i

 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,darr) &
 !$omp private(i)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) darr(:,i) = 0.
 enddo
 !$omp end parallel do

end subroutine zero_boundary_derivs


!----------------------------------------------------------------
!+
!  Set velocity of boundary particles to be their centre of mass
!  velocity
!+
!---------------------------------------------------------------
subroutine set_boundary_particle_velocity(npart,iphase,xyz,vxyz,xyz_CM,vxyz_CM)
 use part, only:iamboundary,iamtype
 integer, intent(in) :: npart
 integer(kind=1), intent(in)   :: iphase(:)
 real, intent(in)              :: xyz(:,:)
 real, intent(inout)           :: vxyz(:,:)
 real, intent(out), optional   :: xyz_CM(3),vxyz_CM(3)
 integer                       :: i
 real                          :: xcm(3),vcm(3)

 call calculate_boundary_core_CM(npart,iphase,xyz,vxyz,xcm,vcm)

 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,vxyz,vcm) &
 !$omp private(i)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) vxyz(1:3,i) = vcm
 enddo
 !$omp end parallel do
 if (present(xyz_CM))  xyz_CM  = xcm
 if (present(vxyz_CM)) vxyz_CM = vcm

end subroutine set_boundary_particle_velocity


!----------------------------------------------------------------
!+
!  Calculates centre of mass and velocity centre of mass for
!  boundary particles
!+
!---------------------------------------------------------------
subroutine calculate_boundary_core_CM(npart,iphase,xyz,vxyz,xyz_CM,vxyz_CM,nboundary_out)
 use part, only:iamboundary,iamtype
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:),vxyz(:,:)
 real, intent(out)   :: xyz_CM(3),vxyz_CM(3)
 integer, intent(out), optional :: nboundary_out
 integer             :: nboundary,i

 nboundary = 0
 xyz_CM = 0.
 vxyz_CM = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,xyz,vxyz) &
 !$omp private(i) &
 !$omp reduction(+:nboundary,xyz_CM,vxyz_CM)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       nboundary = nboundary + 1
       xyz_CM = xyz_CM + xyz(1:3,i)
       vxyz_CM = vxyz_CM + vxyz(1:3,i)
    endif
 enddo
 !$omp end parallel do
 if (nboundary > 0) then
    xyz_CM = xyz_CM / real(nboundary)
    vxyz_CM = vxyz_CM / real(nboundary)
 endif
 if (present(nboundary_out)) nboundary_out = nboundary

end subroutine calculate_boundary_core_CM

end module boundarypart
