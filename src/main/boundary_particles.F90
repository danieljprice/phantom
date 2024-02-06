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
 public :: get_boundary_particle_forces,set_boundary_particle_velocity

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
 use dim,    only:mhd,do_radiation,use_dust,use_dustgrowth,maxvxyzu
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:)
 real, intent(inout) :: fxyzu(:,:)
 real, intent(out), optional :: dBevol(:,:),drad(:,:),ddustevol(:,:),ddustprop(:,:)
 integer             :: nboundary,i
 real, dimension(3)  :: avg_lin_accel,rot_accel_i,avg_torque,xyz_core_CM,torquei

 nboundary = 0
 avg_lin_accel = 0.
 xyz_core_CM = 0.
 !!$omp parallel do default(none) &
 !!$omp shared(npart,iphase,fxyzu,xyz) &
 !!$omp private(i) &
 !!$omp reduction(+:avg_lin_accel,nboundary,xyz_core_CM)
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       nboundary = nboundary + 1
       avg_lin_accel = avg_lin_accel + fxyzu(1:3,i)
       xyz_core_CM = xyz_core_CM + xyz(1:3,i)
    endif
 enddo
 !!$omp end parallel do

 if (nboundary > 0) then
    xyz_core_CM = xyz_core_CM / real(nboundary)
    avg_lin_accel = avg_lin_accel / real(nboundary)

    avg_torque = 0.
    !!loop to calculate solid-body acceleration
    !!$omp parallel do default(none) &
    !!$omp shared(npart,iphase,xyz,fxyzu,avg_lin_accel,xyz_core_CM) &
    !!$omp reduction(+:avg_torque) &
    !!$omp private(i,torquei)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          call cross_product3D(xyz(1:3,i)-xyz_core_CM, fxyzu(1:3,i)-avg_lin_accel, torquei)  ! specific torque
          avg_torque = avg_torque + torquei
       endif
    enddo
    !!$omp end parallel do
    avg_torque = avg_torque / real(nboundary)

    !!$omp parallel do default(none) &
    !!$omp shared(npart,iphase,xyz,fxyzu,avg_lin_accel,avg_torque,xyz_core_CM) &
    !!$omp shared(dBevol,drad,ddustevol,ddustprop) &
    !!$omp private(i,rot_accel_i)
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          call cross_product3D(avg_torque, xyz(1:3,i)-xyz_core_CM, rot_accel_i)
          fxyzu(1:3,i) = avg_lin_accel !+ rot_accel_i
          if (maxvxyzu==4) fxyzu(4,i) = 0.
          if (mhd) dBevol(:,i) = 0.
          if (do_radiation) drad(iradxi,i) = 0.
          if (use_dust) ddustevol(:,i) = 0.
          if (use_dustgrowth) ddustprop(1,i) = 0.
       endif
    enddo
    !!$omp end parallel do
 endif

end subroutine get_boundary_particle_forces


!----------------------------------------------------------------
!+
!  Calculates centre of mass and velocity centre of mass for
!  boundary particles
!+
!---------------------------------------------------------------
subroutine set_boundary_particle_velocity(npart,iphase,xyz,xyz_CM,vxyz_CM,vxyzu)
 use part, only:iamboundary,iamtype
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: iphase(:)
 real, intent(in)    :: xyz(:,:)
 real, intent(inout), optional :: vxyzu(:,:)
 real, intent(out)   :: xyz_CM(3),vxyz_CM(3)
 integer             :: nboundary,i

 nboundary = 0
 xyz_CM = 0.
 vxyz_CM = 0.
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       nboundary = nboundary + 1
       xyz_CM = xyz_CM + xyz(1:3,i)
       if (present(vxyzu)) vxyz_CM = vxyz_CM + vxyzu(1:3,i)
    endif
 enddo
 xyz_CM = xyz_CM / real(nboundary)
 vxyz_CM = vxyz_CM / real(nboundary)

 if (present(vxyzu)) then
    do i=1,npart
       if (iamboundary(iamtype(iphase(i)))) then
          vxyzu(1:3,i) = vxyz_CM
       endif
    enddo
 endif

end subroutine set_boundary_particle_velocity

! subroutine inspect_forces(npart,iphase,fxyzu)
!  use part, only:iamboundary,iamtype
!  integer, intent(in) :: npart
!  integer(kind=1), intent(in) :: iphase(:)
!  real, intent(inout) :: fxyzu(:,:)
!  integer :: i
!  real, dimension(3) :: minf,maxf

!  minf = huge(0.)
!  maxf = 0.
!  do i=1,npart
!     if (iamboundary(iamtype(iphase(i)))) then
!        minf(1) = min(minf(1),abs(fxyzu(1,i)))
!        minf(2) = min(minf(2),abs(fxyzu(2,i)))
!        minf(3) = min(minf(3),abs(fxyzu(3,i)))
!        maxf(1) = max(maxf(1),abs(fxyzu(1,i)))
!        maxf(2) = max(maxf(2),abs(fxyzu(2,i)))
!        maxf(3) = max(maxf(3),abs(fxyzu(3,i)))
!     endif
!  enddo

!  print*,'minf = ',minf
!  print*,'maxf = ',maxf
!  read*

! end subroutine inspect_forces

end module boundarypart
