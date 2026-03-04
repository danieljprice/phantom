!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module set_streamer
!
! This module is contains utilities for setting up particles on a parabolic orbit
! in order to be used as a streamer particle.
! Our conventions for angles are the same as in Xiang-Gruess (2016).
! Eccentricity is set to unity, i.e. for a parabolic orbit:
!
!
!   R_p            : pericentre distance (minimum approach)
!   R_in           : injection distance from the star
!   R_imp          : cylindrical impact radius on the disc (node, z=0)
!   incl_imp_deg   : inclination (degrees) of the orbit at impact (0 prograde)
!
! :References:
!   - Xiang-Gruess (2016), MNRAS 455, 3086-3100
!   - Cuello et al. (2019), MNRAS 483, 4114-4139
!
! :Owner: Cristiano Longarini
!
! :Runtime parameters: None
!
! :Dependencies: physcon, vectorutils
!
 use physcon,    only: pi
 use vectorutils,only: rotatevec
 implicit none
 private
 public :: set_streamer_particle

 integer, parameter :: ierr_ok=0, ierr_badinput=1, ierr_badgeom=2

contains

subroutine set_streamer_particle(mu, R_p, R_in, R_imp, incl_imp_deg, &
                                   x, v, ierr, phi_imp_deg, ingoing)

 real,    intent(in)  :: mu, R_p, R_in, R_imp, incl_imp_deg
 real,    intent(in), optional :: phi_imp_deg
 logical, intent(in), optional :: ingoing
 real,    intent(out) :: x(3), v(3)
 integer, intent(out) :: ierr
 real :: p, Omega_deg, Omega, i_deg, i_rad
 real :: cosw, omega_arg
 real :: ctheta, stheta, r, vp
 real :: ct, st
 real :: xp(3), vpv(3)
 logical :: preperi
 real :: rot_axis(3), incl_angle
 real :: z_axis(3)

 ierr = ierr_ok
 x = 0.0; v = 0.0

 if (mu <= 0.0 .or. R_p <= 0.0 .or. R_in <= 0.0 .or. R_imp <= 0.0) then
    ierr = ierr_badinput
    return
 endif

 ! Semilatus rectum for a parabola
 p = 2.0*R_p

 ! Geometry consistency: real node impact requires p < 2 R_imp
 if (p >= 2.0*R_imp) then
    ierr = ierr_badgeom
    return
 endif

 ! Default angles/flags
 Omega_deg = 0.0
 if (present(phi_imp_deg)) Omega_deg = phi_imp_deg
 preperi   = .true.
 if (present(ingoing)) preperi = ingoing

 Omega = Omega_deg * (pi/180.0)
 i_deg = incl_imp_deg
 i_rad = i_deg * (pi/180.0)

 ! Computing argument of pericentre from the desired R_imp
 cosw  = max(-1.0, min(1.0, p/R_imp - 1.0))
 omega_arg = acos(cosw)

 ! Computing true anomaly at injection R_in
 ctheta = max(-1.0, min(1.0, p/R_in - 1.0))
 stheta = sqrt(max(0.0, 1.0 - ctheta*ctheta))
 if (preperi) stheta = -stheta   ! pre-pericentre => vr < 0 => sin theta < 0

 ct = ctheta; st = stheta
 r  = p / (1.0 + ct)

 ! 2D state (x-y plane)
 vp = sqrt(mu/p)
 xp  = (/ r*ct,      r*st,        0.0 /)
 vpv = (/ -vp*st, vp*(1.0+ct),    0.0 /)

 ! 1) Rotate within the orbital plane by +omega_arg about +z to fix R_imp
 z_axis = (/0.0,0.0,1.0/)
 call rotatevec(xp,  z_axis, omega_arg)
 call rotatevec(vpv, z_axis, omega_arg)

 ! 2) Incline about ascending-node axis as in set_flyby
 incl_angle = i_rad
 rot_axis   = (/ sin(Omega), -cos(Omega), 0.0 /)
 call rotatevec(xp,  rot_axis, incl_angle)
 call rotatevec(vpv, rot_axis, incl_angle)

 x = xp
 v = vpv

end subroutine set_streamer_particle

end module set_streamer
