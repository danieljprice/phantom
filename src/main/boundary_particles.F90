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
 public :: get_boundary_particle_forces

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
 use part, only:iamboundary,iamtype,iradxi
 use dim, only:mhd,do_radiation,use_dust,use_dustgrowth,maxvxyzu
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:)
 real, intent(inout) :: fxyzu(:,:)
 real, intent(out)   :: dBevol(:,:)
 real, intent(out)   :: drad(:,:)
 real, intent(out)   :: ddustevol(:,:),ddustprop(:,:)
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
    avg_lin_accel = avg_lin_accel / real(nboundary)

    avg_torque = 0.
    !loop to calculate solid-body acceleration
    !$omp parallel do default(none) &
    !$omp shared(npart,xyz,fxyzu,avg_lin_accel,xyz_core_CM) &
    !$omp reduction(+:avg_torque) &
    !$omp private(i,torquei)
    do i=1,npart
       call cross_product3D(xyz(1:3,i)-xyz_core_CM, fxyzu(1:3,i)-avg_lin_accel, torquei)  ! specific torque
       avg_torque = avg_torque + torquei
    enddo
    !$omp end parallel do
    avg_torque = avg_torque / real(nboundary)

    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyz,fxyzu,avg_lin_accel,avg_torque,xyz_core_CM) &
    !$omp shared(dBevol,drad,ddustevol,ddustprop) &
    !$omp private(i,rot_accel_i)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          call cross_product3D(avg_torque, xyz(1:3,i)-xyz_core_CM, rot_accel_i)
          fxyzu(1:3,i) = avg_lin_accel + rot_accel_i
          if (maxvxyzu==4) fxyzu(4,i) = 0.
          if (mhd) dBevol(:,i) = 0.
          if (do_radiation) drad(iradxi,i) = 0.
          if (use_dust) ddustevol(:,i) = 0.
          if (use_dustgrowth) ddustprop(1,i) = 0.
       endif
    enddo
    !$omp end parallel do
 endif

end subroutine get_boundary_particle_forces

end module boundarypart
