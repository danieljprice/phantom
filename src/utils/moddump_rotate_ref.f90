!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Solid body rotation for all particles to a given mean inclination and mean position angle
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters: None
!
! :Dependencies: part, physcon, vectorutils
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use physcon,              only:pi
 use part,                 only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use vectorutils,          only:rotatevec
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: alpha,gamma
 real                   :: ref_incl,ref_Omega
 real                   :: incl,Omega
 real                   :: a,b,c,d,e,f,g,h,i
 real                   :: temp_x,temp_y,temp_z
 real                   :: temp_vx,temp_vy,temp_vz
 real                   :: temp(3),temp_v(3)
 integer                :: j

 ! Angles of the current disc in deg
 incl = 50.864
 Omega = 34.36 !38.5

 ! Reference angles in deg : the disc will be put into that plane
 ref_incl = 54.6
 ref_Omega = 53    ! Omega in code units

 ! Rotation angles & coeffs
 alpha = (ref_incl - incl) *pi/180    !about x
 gamma = (ref_Omega - Omega) *pi/180     !about z

 !
 a = cos(alpha)
 b = -sin(alpha)*cos(gamma)
 c = sin(alpha)*sin(gamma)
 d = sin(alpha)
 e = cos(alpha)*cos(gamma)
 f = -cos(alpha)*sin(gamma)
 g = 0
 h = sin(gamma)
 i = cos(gamma)

 ! New positions & velocities
 do j = 1,npart
    temp(1) = xyzh(1,j)
    temp(2) = xyzh(2,j)
    temp(3) = xyzh(3,j)
    call rotatevec(temp, (/1.0,0.,0./), alpha)
    call rotatevec(temp, (/0.,0.,1.0/), gamma)

    temp_v(1) = vxyzu(1,j)
    temp_v(2) = vxyzu(2,j)
    temp_v(3) = vxyzu(3,j)
    call rotatevec(temp_v, (/1.0,0.,0./), alpha)
    call rotatevec(temp_v, (/0.,0.,1.0/), gamma)

!    temp_x = a*xyzh(1,j) + b*xyzh(2,j) + c*xyzh(3,j)
!    temp_y = d*xyzh(1,j) + e*xyzh(2,j) + f*xyzh(3,j)
!    temp_z = g*xyzh(1,j) + h*xyzh(2,j) + i*xyzh(3,j)
!
!    temp_vx = a*vxyzu(1,j) + b*vxyzu(2,j) + c*vxyzu(3,j)
!    temp_vy = d*vxyzu(1,j) + e*vxyzu(2,j) + f*vxyzu(3,j)
!    temp_vz = g*vxyzu(1,j) + h*vxyzu(2,j) + i*vxyzu(3,j)

    temp_x = temp(1)
    temp_y = temp(2)
    temp_z = temp(3)
    xyzh(1,j) = temp_x
    xyzh(2,j) = temp_y
    xyzh(3,j) = temp_z

    temp_vx = temp_v(1)
    temp_vy = temp_v(2)
    temp_vz = temp_v(3)
    vxyzu(1,j) = temp_vx
    vxyzu(2,j) = temp_vy
    vxyzu(3,j) = temp_vz

 enddo

 do j = 1,nptmass

    temp(1) = xyzmh_ptmass(1,j)
    temp(2) = xyzmh_ptmass(2,j)
    temp(3) = xyzmh_ptmass(3,j)
    call rotatevec(temp, (/1.0,0.,0./), alpha)
    call rotatevec(temp, (/0.,0.,1.0/), gamma)

    temp_v(1) = vxyz_ptmass(1,j)
    temp_v(2) = vxyz_ptmass(2,j)
    temp_v(3) = vxyz_ptmass(3,j)
    call rotatevec(temp_v, (/1.0,0.,0./), alpha)
    call rotatevec(temp_v, (/0.,0.,1.0/), gamma)

    temp_x = temp(1)
    temp_y = temp(2)
    temp_z = temp(3)
    xyzmh_ptmass(1,j) = temp_x
    xyzmh_ptmass(2,j) = temp_y
    xyzmh_ptmass(3,j) = temp_z

    temp_vx = temp_v(1)
    temp_vy = temp_v(2)
    temp_vz = temp_v(3)
    vxyz_ptmass(1,j) = temp_vx
    vxyz_ptmass(2,j) = temp_vy
    vxyz_ptmass(3,j) = temp_vz
!    temp_x = a*xyzmh_ptmass(1,j) + b*xyzmh_ptmass(2,j) + c*xyzmh_ptmass(3,j)
!    temp_y = d*xyzmh_ptmass(1,j) + e*xyzmh_ptmass(2,j) + f*xyzmh_ptmass(3,j)
!    temp_z = g*xyzmh_ptmass(1,j) + h*xyzmh_ptmass(2,j) + i*xyzmh_ptmass(3,j)
!
!    temp_vx = a*vxyz_ptmass(1,j) + b*vxyz_ptmass(2,j) + c*vxyz_ptmass(3,j)
!    temp_vy = d*vxyz_ptmass(1,j) + e*vxyz_ptmass(2,j) + f*vxyz_ptmass(3,j)
!    temp_vz = g*vxyz_ptmass(1,j) + h*vxyz_ptmass(2,j) + i*vxyz_ptmass(3,j)
!
!
!    xyzmh_ptmass(1,j) = temp_x
!    xyzmh_ptmass(2,j) = temp_y
!    xyzmh_ptmass(3,j) = temp_z
!
!    vxyz_ptmass(1,j) = temp_vx
!    vxyz_ptmass(2,j) = temp_vy
!    vxyz_ptmass(3,j) = temp_vz
!
 enddo

end subroutine modify_dump

end module moddump

