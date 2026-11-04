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
           calculate_boundary_core_CM

 private

contains
!----------------------------------------------------------------
!+
!  Reassign every boundary particle their average linear acceleration
!  and a rotational acceleration that preserves their collective
!  angular momentum while enforcing co-rotation
!+
!---------------------------------------------------------------
subroutine get_boundary_particle_forces(npart,iphase,xyz,fxyzu,dBevol,drad,ddustprop,ddustevol)
 use part,   only:iamboundary,iamtype,iradxi
 use dim,    only:mhd,do_radiation,use_dust,use_dustgrowth
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:)
 real, intent(inout) :: fxyzu(:,:)
 real, intent(inout), optional :: dBevol(:,:),drad(:,:),ddustevol(:,:),ddustprop(:,:)
 integer             :: nboundary,i
 real, dimension(3)  :: avg_lin_accel,rot_accel_i,avg_torque,xyz_core_CM,torquei

 nboundary = 0
 avg_lin_accel = 0.
 xyz_core_CM = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,fxyzu,xyz) &
 !$omp private(i) &
 !$omp reduction(+:avg_lin_accel,nboundary,xyz_core_CM)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       nboundary = nboundary + 1
       avg_lin_accel = avg_lin_accel + fxyzu(1:3,i)
       xyz_core_CM = xyz_core_CM + xyz(1:3,i)
    endif
 enddo
 !$omp end parallel do

 if (nboundary > 0) then
    xyz_core_CM = xyz_core_CM / real(nboundary)
    avg_lin_accel = avg_lin_accel / real(nboundary)

    avg_torque = 0.
    ! loop to calculate solid-body acceleration
    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyz,fxyzu,avg_lin_accel,xyz_core_CM) &
    !$omp reduction(+:avg_torque) &
    !$omp private(i,torquei)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          call cross_product3D(xyz(1:3,i)-xyz_core_CM, fxyzu(1:3,i)-avg_lin_accel, torquei)  ! specific torque
          avg_torque = avg_torque + torquei
       endif
    enddo
    !$omp end parallel do
    avg_torque = avg_torque / real(nboundary)

    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyz,fxyzu,avg_lin_accel,avg_torque,xyz_core_CM) &
    !$omp private(i,rot_accel_i)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          call cross_product3D(avg_torque, xyz(1:3,i)-xyz_core_CM, rot_accel_i)
          fxyzu(1:3,i) = avg_lin_accel !+ rot_accel_i
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
subroutine calculate_boundary_core_CM(npart,iphase,xyz,vxyz,xyz_CM,vxyz_CM)
 use part, only:iamboundary,iamtype
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:),vxyz(:,:)
 real, intent(out)   :: xyz_CM(3),vxyz_CM(3)
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

end subroutine calculate_boundary_core_CM

end module boundarypart
