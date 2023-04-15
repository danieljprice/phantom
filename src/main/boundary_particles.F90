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
!  Reassign every boundary particle the average acceleration of
!  all boundary particles.
!+
!---------------------------------------------------------------
subroutine get_boundary_particle_forces(npart,iphase,fxyzu,dBevol,drad,ddustprop,ddustevol)
 use part, only:iamboundary,iamtype,iradxi
 use dim, only:mhd,do_radiation,use_dust,use_dustgrowth,maxvxyzu
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(inout) :: fxyzu(:,:)
 integer             :: nboundary,i
 real                :: avg_accel(1:3)
 real, intent(out)   :: dBevol(:,:)
 real, intent(out)   :: drad(:,:)
 real, intent(out)   :: ddustevol(:,:),ddustprop(:,:)

 nboundary = 0
 avg_accel = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,fxyzu) &
 !$omp private(i) &
 !$omp reduction(+:avg_accel,nboundary)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       nboundary = nboundary + 1
       avg_accel = avg_accel + fxyzu(1:3,i)
    endif
 enddo
 !$omp end parallel do

 if (nboundary > 0) then
    avg_accel = avg_accel / real(nboundary)
    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,fxyzu,avg_accel,dBevol,drad,ddustevol,ddustprop) &
    !$omp private(i)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          fxyzu(1:3,i) = avg_accel
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


!----------------------------------------------------------------
!+
!  Reassign every boundary particle a unique velocity that would
!  that preserves their collective angular momentum while enforcing
!  co-rotation
!+
!---------------------------------------------------------------
subroutine average_boundary_particle_rotation(npart,iphase,fxyzu,xyz)
 use part,        only:iamboundary,iamtype
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:)
 real, intent(inout) :: fxyzu(:,:)
 integer             :: i,nboundary
 real, dimension(3)  :: accel_i,torquei,avg_torque,xyz_core_CM
 real                :: sep2

 ! Calculate core centre of mass
 nboundary = 0
 xyz_core_CM = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,iphase,fxyzu,xyz) &
 !$omp private(i) &
 !$omp reduction(+:xyz_core_CM,nboundary)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       nboundary = nboundary + 1
       xyz_core_CM = xyz_core_CM + xyz(1:3,i)
    endif
 enddo
!$omp end parallel do

 if (nboundary > 0) then
    xyz_core_CM = xyz_core_CM / nboundary

    ! Calculate average torque w.r.t. centre of mass
    avg_torque = 0.
    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,fxyzu,xyz_core_CM,xyz) &
    !$omp private(i,torquei) &
    !$omp reduction(+:avg_torque)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          call cross_product3D(xyz(1:3,i) - xyz_core_CM, fxyzu(1:3,i), torquei)
          avg_torque = avg_torque + torquei
       endif
    enddo
    !$omp end parallel do
    avg_torque = avg_torque / nboundary

    !$omp parallel do default(none) &
    !$omp shared(npart,iphase,xyz,fxyzu,avg_torque,xyz_core_CM) &
    !$omp private(i,accel_i,sep2)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          sep2 = dot_product(xyz(1:3,i)-xyz_core_CM, xyz(1:3,i)-xyz_core_CM)
          call cross_product3D(avg_torque / sep2, xyz(1:3,i)-xyz_core_CM, accel_i)
          fxyzu(1:3,i) = fxyzu(1:3,i) + accel_i
       endif
    enddo
    !$omp end parallel do
 endif

end subroutine average_boundary_particle_rotation

end module boundarypart
