!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dem
!
! This module implements the soft sphere discrete element method (DEM)
! for sink-sink interactions.
!
! epsilon is a parameter [0,1] that controls the restitution of the collision.
! epsilon = 1 is perfectly elastic, epsilon = 0 is perfectly inelastic.
! epsilon = 0.5 is a reasonable default.
!
! :References: Schwartz+2012, Granular Matter 14, 363-380
!
! Optional tensile cohesion (linear spring when spheres are slightly separated)
! mimics van der Waals / regolith bond strength at sub-contact scale.
!
! :Owner: Daniel Price
!
 implicit none
 private

 public :: get_ssdem_force,dem_cohesion_summary

 real, public :: ct_dem = 0.1         ! Tangential damping coefficient
 real, public :: epsilon_n_dem = 0.5  ! Normal coefficient of restitution (user-settable)
 real, public :: kn_cgs = 1e7         ! Spring constant (e.g. 10^4 kg/s^2 = 10^7 g/s^2)
 real, public :: kt_cgs = 0.          ! Tensile spring constant (g/s^2 per cm gap); 0 = no cohesion
 real, public :: coh_gap_max_cgs = 0. ! Max surface gap (cm) for cohesive bond; 0 = use 1% of mean radius

contains

!----------------------------------------------------------------
!+
!  Soft-sphere DEM normal force (Hooke's law)
!  Implements Eq. (3) from Schwartz+2012 for overlapping spheres
!+
!----------------------------------------------------------------
subroutine get_ssdem_force(Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,fx,fy,fz,veli,velj,wi,wj,dtmin)
 use physcon,     only:pi
 use vectorutils, only:cross_product
 use units,       only:umass,utime,udist
 real, intent(in)    :: Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,veli(3),velj(3),wi(3),wj(3)
 real, intent(inout) :: fx,fy,fz,dtmin
 real :: r,overlap,gap,kn,kn_dem,kt_dem,coh_gap_max
 real :: cn,ct,reduced_mass,log_epsilon_n_dem,li,lj
 real :: nvec(3),vrel(3),n_cross_wi(3),n_cross_wj(3),u_dot_n,u_n(3),u_t(3)

 !----------------------------------------------------------------
 ! Normal force
 !----------------------------------------------------------------
 r = 1.0 / ddr
  ! Normal unit vector
 nvec(1) = dx * ddr
 nvec(2) = dy * ddr
 nvec(3) = dz * ddr

 overlap = Rsinki + Rsinkj - r
 gap = r - Rsinki - Rsinkj
 kn = 0.
 kn_dem = kn_cgs / (umass/utime**2)  ! convert to code units
 kt_dem = kt_cgs / (umass/utime**2)
 reduced_mass = mj * mi / (mj + mi)
 if (coh_gap_max_cgs > 0.) then
    coh_gap_max = coh_gap_max_cgs / udist
 else
    coh_gap_max = 0.01 * (Rsinki + Rsinkj)
 endif
 if (overlap > 0.0) then
    kn = kn_dem
    fx = fx + kn * overlap * nvec(1) / mj
    fy = fy + kn * overlap * nvec(2) / mj
    fz = fz + kn * overlap * nvec(3) / mj

    !----------------------------------------------------------------
    ! Damping force (only when overlapping; Schwartz+2012 Eqs. 8-15;
    ! contact-only by design, see docs/DEM_COHESION_MIGRATION.md #1.4)
    !----------------------------------------------------------------
    ! Cross products: n x omega, where omega is the spin vector of the sphere
    n_cross_wi = cross_product(nvec,wi)
    n_cross_wj = cross_product(nvec,wj)

    ! Eqs. (9) and (10)from Schwartz+2012
    li = (Rsinki**2 - Rsinkj**2 + r**2) / (2.0 * r)
    lj = (Rsinkj**2 - Rsinki**2 + r**2) / (2.0 * r)

    ! Relative velocity at contact point (Eq. 8 from Schwartz+2012)
    vrel = veli - velj + li * n_cross_wi - lj * n_cross_wj

    ! Normal and tangential components
    u_dot_n = dot_product(vrel, nvec)
    u_n = u_dot_n * nvec
    u_t = vrel - u_n

    ! Eqn (15) from Schwartz+2012
    log_epsilon_n_dem = log(epsilon_n_dem)
    cn = -2.0 * sqrt(reduced_mass * kn) * log_epsilon_n_dem / sqrt(pi**2 + log_epsilon_n_dem**2)
    ct = 0.

    ! Damping forces
    fx = fx - cn * u_n(1) / mj - ct * u_t(1)
    fy = fy - cn * u_n(2) / mj - ct * u_t(2)
    fz = fz - cn * u_n(3) / mj - ct * u_t(3)
 elseif (kt_dem > 0. .and. gap > 0. .and. gap < coh_gap_max) then
    ! tensile spring when spheres are bonded but not overlapping
    kn = kt_dem
    fx = fx - kn * gap * nvec(1) / mj
    fy = fy - kn * gap * nvec(2) / mj
    fz = fz - kn * gap * nvec(3) / mj
 endif

 ! spring timescale sets the stable timestep in both regimes (overlap or
 ! tensile bond), independent of whether damping is applied (contact-only)
 if (kn > 0.) dtmin = min(dtmin,sqrt(reduced_mass/kn))

end subroutine get_ssdem_force

!----------------------------------------------------------------
!+
!  Print DEM cohesion settings after reading the .in file
!+
!----------------------------------------------------------------
subroutine dem_cohesion_summary
 use io,      only:iprint
 use units,   only:umass,utime
 real :: kt_dem

 if (kt_cgs <= 0.) return
 kt_dem = kt_cgs / (umass/utime**2)
 write(iprint,"(/,a)") ' DEM tensile cohesion enabled'
 write(iprint,"(a,1pg12.4,a)") '   kt_cgs = ',kt_cgs,' g/s^2 per cm surface gap'
 if (coh_gap_max_cgs > 0.) then
    write(iprint,"(a,1pg12.4,a)") '   coh_gap_max = ',coh_gap_max_cgs,' cm'
 else
    write(iprint,"(a)") '   coh_gap_max = 1% of (R_i + R_j) per pair (default)'
 endif
 write(iprint,"(a,1pg12.4)") '   kt (code units) = ',kt_dem
end subroutine dem_cohesion_summary

end module dem
