!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: bowen_dust
!
!  DESCRIPTION:
!  Acceleration on gas due to dust grains at radiative equilibrium
!  The model is described in the following article:
!  G.H. Bowen - Dynamical modeling of long-period variable star atmospheres (1988)
!
!  REFERENCES: None
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, kernel, part, physcon, units
!+
!--------------------------------------------------------------------------
module bowen_dust
 implicit none

 public :: radiative_acceleration, init_bowen
 logical, parameter :: use_alpha_wind = .false.

 private
 integer, parameter :: N = 1024
 integer, parameter :: wind_emitting_sink = 1
 logical, parameter :: verbose = .false.
 real, parameter :: alpha_wind = 1.d0
 real :: kappa_gas, kmax, L, c_light, specific_energy_to_T_ratio, Cprime, Tcond, deltaT, Teff,&
      Reff, omega_osc, deltaR_osc, usteboltz,wind_inject, Rmin
 integer :: nwall_particles

contains


!-----------------------------------------------------------------------
!+
!  Computes acceleration due to radiation pressure on dust on all SPH particles
!  is is assumed that optical depth is roughly spherically symmetric around the central star
!+
!-----------------------------------------------------------------------
subroutine radiative_acceleration(npart, xyzh, vxyzu, dt, fext, fxyzu_shock, time)
 use part,     only:isdead_or_accreted,rhoh,xyzmh_ptmass,massoftype,igas,divcurlv
 use physcon,  only:pi,mass_proton_cgs
 use eos,      only:gamma
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(in)    :: dt,time
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(inout) :: fext(:,:),fxyzu_shock(:,:)

 real :: r(3, npart), d(npart), z_coord(npart), d2_axis(npart)
 real :: dmin, dmax, dr, part_mass
 integer :: nfound, found(npart)
 real :: rho_grid(N), r_grid(N), Teq(N), a_over_r(N), kappa_dust(N)!, T_grid(N), dlnkap_dlnT(N), tau_prime(N)
 real :: dQ_dt,Mstar,R_star,L_star,Rforce,alpha_surface,dist,frad,frepuls,Teq_part,xr(3),x_star(3)
 real :: force,fac,Trad,kap_dust,alpha_w,rho,T,Qcool,dQcool_dlnT,mass_per_H
 logical :: skip
 integer :: i,j,icooling


 if (npart == nwall_particles) return
 mass_per_H = 1.4 * mass_proton_cgs
 part_mass = massoftype(igas)
 x_star(1:3) = xyzmh_ptmass(1:3,wind_emitting_sink)
 Mstar = xyzmh_ptmass(4,wind_emitting_sink)
 !R_star = Reff + deltaR_osc*sin(omega_osc*time)
 R_star = wind_inject + deltaR_osc*sin(omega_osc*time)
 L_star = 4.*pi*usteboltz*R_star**2*Teff**4
 print *,'time ###  ##',time,R_star,Reff,npart,nwall_particles,L,L_star

 if (use_alpha_wind) then
    force = 1.8
    Rforce = force*R_star
!$omp parallel do default(none) &
!$omp shared(npart,xyzh,vxyzu,fext,divcurlv,x_star,Rforce,Mstar,R_star,gamma,dt,nwall_particles) &
!$omp private(xr,alpha_surface,dist)
   do i=nwall_particles+1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          xr(1:3) = xyzh(1:3,i)-x_star(1:3)
          dist = sqrt(xr(1)**2 + xr(2)**2 + xr(3)**2)
          ! repulsive force, Ayliffe & Bate 2010, MNRAS 408, 876
          if (dist < Rforce) then
             alpha_surface = ((Rforce-dist)/(Rforce-R_star))**4
          else
             alpha_surface = 0.
          endif
          fext(1:3,i) = fext(1:3,i) + Mstar*max(alpha_wind,alpha_surface)/dist**3*xr(1:3)
          !update internal energy: add pdV work
          vxyzu(4,i) = vxyzu(4,i) -(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
       endif
    enddo
!$omp end parallel do
 else
    icooling = 3
    force = 1.2
    Rforce = force*R_star
    fac = L_star/(4.*pi*c_light*Mstar)
    print *,'kmax =====>',fac*kmax
!!$omp parallel do default(none) &
!!$omp shared(npart,xyzh,vxyzu,fext,divcurlv,x_star,Rforce,fac,Teff,R_star,Mstar,gamma,kappa_gas,dt,nwall_particles) &
!!$omp private(xr,alpha_surface,alpha_w,kap_dust,Trad,dist)
   do i=nwall_particles+1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          xr(1:3) = xyzh(1:3,i)-x_star(1:3)
          dist = sqrt(xr(1)**2 + xr(2)**2 + xr(3)**2)
          ! repulsive force, Ayliffe & Bate 2010, MNRAS 408, 876
          if (dist < Rforce) then
             alpha_surface = ((Rforce-dist)/(Rforce-R_star))**4
          else
             alpha_surface = 0.
          endif
          Trad = Teff*sqrt(R_star/dist)
          call get_kappa_dust(Trad, kap_dust)
          alpha_w = fac*(kap_dust+kappa_gas)
          !print *,i,alpha_w,alpha_surface,kap_dust,kappa_gas,Trad,Teff
          fext(1:3,i) = fext(1:3,i) + Mstar*max(alpha_w,alpha_surface)/dist**3*xr(1:3)
          ! !update internal energy: add pdV work
          ! vxyzu(4,i) = vxyzu(4,i) -(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
          ! !cooling terms
          ! T = vxyzu(4,i)/specific_energy_to_T_ratio
          ! rho = rhoh(xyzh(4,i), part_mass)
          ! call calc_cooling(icooling,mass_per_H,T,Trad,rho,Qcool,dQcool_dlnT)
          ! !print *,i,vxyzu(4,i),Qcool*dt,Trad,T
          ! if (dt  >=  Cprime/rhoh(xyzh(4,i),part_mass)) then   ! assume thermal equilibrium
          !   vxyzu(4,i) = Trad*specific_energy_to_T_ratio
          !   fxyzu(4,i) = 0.d0
          ! else
          !   !Qcool = -(vxyzu(4,i) - Trad*specific_energy_to_T_ratio) &
          !   !     * (rhoh(xyzh(4,i),part_mass)/Cprime)
          !   vxyzu(4,i) = vxyzu(4,i) + Qcool*dt
          !   !fext(4,i) = fext(4,i) + Qcool
          ! endif
       endif
    enddo
!!$omp end parallel do
    call implicit_wind_cooling(icooling,npart,xyzh,vxyzu,fxyzu_shock,fext,R_star,x_star,dt)
 endif
end subroutine radiative_acceleration

!-----------------------------------------------------------------------
!+
!  Convert parameters into code units.
!+
!-----------------------------------------------------------------------
subroutine init_bowen(u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,wind_injection_radius,&
       bowen_Cprime, bowen_Tcond, bowen_delta,bowen_Teff,wind_osc_vamplitude,wind_osc_period,nwall)
 use physcon,     only: solarl,c,steboltz,pi,radconst
 use units,       only: udist, umass, utime
 integer, intent(in) :: nwall
 real, intent(in)  :: u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,wind_injection_radius,&
         bowen_Cprime, bowen_Tcond, bowen_delta, bowen_Teff, wind_osc_vamplitude,wind_osc_period

 kappa_gas = bowen_kappa / (udist**2/umass)
 kmax = bowen_kmax / (udist**2/umass)
 L = bowen_L /(umass*udist**2/utime**3)
 c_light = c / (udist/utime)
 specific_energy_to_T_ratio = u_to_temperature_ratio!/(udist/utime)**2
 Cprime = bowen_Cprime / (umass*utime/udist**3)
 Tcond = bowen_Tcond
 deltaT = bowen_delta
 Teff  = bowen_Teff
 Reff = sqrt(bowen_L/(4.*pi*steboltz*bowen_Teff**4))/udist
 omega_osc = 2.*pi/wind_osc_period
 deltaR_osc = wind_osc_vamplitude/omega_osc
 nwall_particles = nwall
 wind_inject = wind_injection_radius
 usteboltz = L/(4.*pi*Reff**2*Teff**4)
 Rmin = Reff - deltaR_osc

end subroutine init_bowen

!-----------------------------------------------------------------------
!+
!  Interpolates a quantity computed on the discretized line of sight for all SPH particles
!  (spherical symmetry assumed)
!+
!-----------------------------------------------------------------------
subroutine interpolate_on_particles(npart, d, N, dmax, quantity, output)
 integer, intent(in)  :: npart
 real,    intent(in)  :: d(npart)
 integer, intent(in)  :: N
 real,    intent(in)  :: dmax, quantity(N)
 real,    intent(out) :: output(npart)

 real :: r, dr
 integer :: i, j

 dr = dmax / N
 do i=nwall_particles+1,npart
    r = d(i)
    j = min(int(r/dr),N-1)
    output(i) = (r-dr*j)*(quantity(j+1)-quantity(j))/dr + quantity(j)
 enddo
end subroutine interpolate_on_particles

!-----------------------------------------------------------------------
!+
!  Calculates kd, the dust opacity along the half-line
!   Inputs:
!     * N: number of steps
!     * Teq(N): radiative equilibrium temperature along the half-line
!   Outputs:
!     * kd(N): dust mass opacity the half-line
!+
!-----------------------------------------------------------------------
subroutine calculate_kd(N, Teq, kd)
 integer, intent(in)  :: N
 real,    intent(in)  :: Teq(N)
 real,    intent(out) :: kd(N)

 kd = kmax/(1. + exp((Teq-Tcond)/deltaT))
end subroutine calculate_kd

subroutine get_kappa_dust(Teq, kappa_dust)
 real,    intent(in)  :: Teq
 real,    intent(out) :: kappa_dust

 kappa_dust = kmax/(1.d0 + exp((Teq-Tcond)/deltaT))

end subroutine get_kappa_dust

!-----------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature along the half-line
!+
!-----------------------------------------------------------------------
subroutine calculate_Teq(N, R_star, tau_prime, OR2, Teq)
 integer, intent(in)  :: N
 real,    intent(in)  :: R_star, tau_prime(N), OR2(N)
 real,    intent(out) :: Teq(N)

 where(OR2  <  R_star**2)
    Teq = Teff
 elsewhere
    Teq = Teff*(0.5*(1.-sqrt(1.-(R_star**2/OR2))) + 0.75*tau_prime)**(1./4.)
 end where
end subroutine calculate_Teq

!-----------------------------------------------------------------------
!+
!  Calculates tau_prime along the discretized line of sight
!+
!-----------------------------------------------------------------------
subroutine calculate_tau_prime(N, dmax, R_star, kappa, kap, rho_over_r2, OR, tau_prime)
 integer, intent(in)  :: N
 real,    intent(in)  :: dmax, R_star, kappa, kap(N), rho_over_r2(2*N+1)
 real,    intent(out) :: OR(N), tau_prime(N)

 real :: dr, fact(N)
 integer :: i

 tau_prime(N) = 0.
 dr = dmax/N
 fact = dr/4. * R_star**2 * (kappa+kap)
 do i=1,N-1
    OR(i) = i*dr
    tau_prime(N-i) = tau_prime(N-i+1) +&
                       fact(i)*(rho_over_r2(i)+rho_over_r2(i+1)+rho_over_r2(2*N-i+1)+rho_over_r2(2*N-i+2))
 enddo
 OR(N) = dmax
end subroutine calculate_tau_prime

!-----------------------------------------------------------------------
!+
!  compute radiative equilibrium structure (Lucy approximation)
!+
!-----------------------------------------------------------------------
subroutine density(npart, position, distance2, part_h, part_mass, nfound, found,&
                     rmin, rmax, r_star, N, rho, rho_over_r2)
 use kernel, only:cnormk,wkern
 integer, intent(in)  :: npart
 real,    intent(in)  :: position(npart), distance2(npart), part_h(npart), part_mass
 integer, intent(in)  :: nfound, found(npart)
 real,    intent(in)  :: rmin, rmax, r_star
 integer :: N
 real,    intent(out) :: rho(N), rho_over_r2(N)

 real :: OH, PH2, OR(N), OR2, HR, q2, q, fact0, fact, h, h2
 real :: delta_r, rmin_o, rmin_p, rmax_p, dr, r_star2
 integer :: i, np, j, j_min, j_max

 ! Discretization of the line :
 dr = (rmax-rmin)/(N-1.)
 !     OR(i) = dr*(i-1)+rmin
 !     i = (OR(i)-rmin)/dr + 1

 ! (Very) little optimization :
 rmin_o = rmin - dr
 !     OR(i) = dr*i+rmin_o
 !     i = (OR(i)-rmin_o)/dr

 rho(:) = 0.
 fact0 = part_mass*cnormk
 do i=1,N
    OR(i) = dr*i+rmin_o
 enddo
 do i=1,nfound
    np = found(i)
    OH = position(np)
    PH2 = distance2(np)
    h = part_h(np)
    h2 = h*h
    delta_r = sqrt(4.*h2 - PH2)
    ! rmin_p and rmax_p are the positions on the line of the two intersections between the line and the interaction sphere
    rmin_p = OH-delta_r
    rmax_p = OH+delta_r
    j_min = ceiling((rmin_p-rmin_o)/dr)
    j_max = floor((rmax_p-rmin_o)/dr)
    j_min = max(1, j_min)
    j_max = min(N, j_max)
    ! Adds the contribution of the particle i on the density at all the discretized locations in the interaction sphere
    fact = fact0/h**3
    do j=j_min, j_max
       HR = OR(j) - OH
       q2 = (PH2+HR**2)/h2
       q = sqrt(q2)
       rho(j) = rho(j) + fact*wkern(q2,q)
    enddo
 enddo

 ! rho_over_r2 = 0 inside the star so that we do not divide by zero!
 r_star2 = r_star**2
 do i=1,N
    OR2 = OR(i)**2
    if (OR2  <  r_star2) then
       rho_over_r2(i) = 0.
    else
       rho_over_r2(i) = rho(i)/OR2
    endif
 enddo

end subroutine density

!-----------------------------------------------------------------------------------
!+
!  Select interesting particles (whose interaction zone intersects the line of sight)
!+
!-----------------------------------------------------------------------------------
subroutine select_particles(npart, d2_axis, part_h, nfound, found)
 integer, intent(in)  :: npart
 real,    intent(in)  :: d2_axis(npart), part_h(npart)
 integer, intent(out) :: nfound, found(npart)

 integer :: i
 real :: h

 nfound = 0
 do i=nwall_particles+1,npart
    h = part_h(i)
    if (d2_axis(i)  <  4.*h**2) then
       nfound = nfound + 1
       found(nfound) = i
    endif
 enddo
end subroutine select_particles

!-----------------------------------------------------------------------
!+
!  Project SPH particles on the z-axis and calculate distance to axis
!+
!-----------------------------------------------------------------------
subroutine project_on_z(npart, r, z_coord, d2_axis)
 integer, intent(in)  :: npart
 real,    intent(in)  :: r(3,npart)
 real,    intent(out) :: z_coord(npart), d2_axis(npart)

 z_coord = r(3,:)
 d2_axis = r(1,:)**2+r(2,:)**2
end subroutine project_on_z

!-----------------------------------------------------------------------
!+
!  Project SPH particles on a given axis specified by the angles theta and phi
!  and calculate distance to axis
!+
!-----------------------------------------------------------------------
!  subroutine project_on_line(npart, r, theta, phi, z_coord, d2_axis)
!    integer, intent(in)  :: npart
!    real,    intent(in)  :: r(3,npart), theta, phi
!    real,    intent(out) :: z_coord(npart), d2_axis(npart)

!    integer :: i
!    real :: u(3), OP(3), HP(3), p

 ! Direction vector
!    u(1) = sin(theta)*cos(phi)
!    u(2) = sin(theta)*sin(phi)
!    u(3) = cos(theta)

!    do i=1,npart
!      OP = r(:,i)
!      p = dot_product(OP,u)
!      z_coord(i) = p
!      HP = OP - p*u
!      d2_axis(i) = dot_product(HP,HP)
!    enddo
!  end subroutine


!-----------------------------------------------------------------------
!+
! Computes the coordinates and radial distance of SPH particles in the
! reference frame attached to the wind emitting star
!+
!-----------------------------------------------------------------------
subroutine center_star(npart, xyzh, xzy_origin, r, d, dmin, dmax)
 use dim,      only:maxp
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(4,maxp), xzy_origin(3)
 real,    intent(out) :: r(3,npart), d(npart), dmin, dmax

 integer :: i

 do i=1,npart
    r(:,i) = xyzh(1:3,i)-xzy_origin
 enddo
 d = sqrt(r(1,:)**2 + r(2,:)**2 + r(3,:)**2)
 dmin = minval(d(nwall_particles+1:npart))
 dmax = maxval(d(nwall_particles+1:npart))
end subroutine center_star

!-----------------------------------------------------------------------
!+
!  wind cooling. Implicit method to get the internal energy
!+
!-----------------------------------------------------------------------
subroutine implicit_wind_cooling(icooling,npart,xyzh,vxyzu,fxyzu,fext,R_star,x_star,dt)
 use eos,      only:gamma
 use part,     only:isdead_or_accreted,rhoh,massoftype,igas,divcurlv
 use physcon,  only:mass_proton_cgs

 integer, intent(in)    :: npart,icooling
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(in)    :: dt,R_star,x_star(3)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(inout) :: fxyzu(:,:),fext(:,:)

 real, parameter :: tol = 1.d-4 ! to be adjusted
 real :: T,Trad,rho,term1,term2,term3,delta_e,mass_per_H,part_mass,Qcool,dQcool_dlnT
 real :: dist,xr(3),q1,q2,dq1,dq2
 integer, parameter :: iter_max = 200
 integer :: i,iter

 mass_per_H = 1.4 * mass_proton_cgs
 part_mass = massoftype(igas)
! !$omp parallel do default(none) &
! !$omp shared(icooling,npart,xyzh,vxyzu,fxyzu,divcurlv) &
! !$omp shared(gamma,part_mass,mass_per_H,dt,nwall_particles) &
! !$omp private(iter,term1,term2,term3,delta_e,Qcool,dQcool_dlnT,T,Tg,rho)
 do i=nwall_particles+1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       xr(1:3) = xyzh(1:3,i)-x_star(1:3)
       dist = sqrt(xr(1)**2 + xr(2)**2 + xr(3)**2)
       iter = 0
       delta_e = 1.d-3
       Trad = Teff*sqrt(R_star/dist)
       rho = rhoh(xyzh(4,i), part_mass)
       term1 = vxyzu(4,i)+fxyzu(4,i)*dt !includes shock heating term
       term2 = 1.+(gamma-1.)*dt*divcurlv(1,i)
       do while (abs(delta_e) > tol .and. iter < iter_max)
          T = vxyzu(4,i)/specific_energy_to_T_ratio
          call calc_cooling(1,mass_per_H,T,Trad,rho,Q1,dQ1)
          call calc_cooling(2,mass_per_H,T,Trad,rho,Q2,dQ2)
          Qcool = q1+q2
          dQcool_dlnT = dq1+dq2
          term3 = vxyzu(4,i)*term2-Qcool*dt
          delta_e = (term1-term3)/(term2-dQcool_dlnT*dt/vxyzu(4,i))
          !vxyzu(4,i) = vxyzu(4,i)*(1.+delta_e)
          vxyzu(4,i) = vxyzu(4,i)+delta_e
          iter = iter + 1
       enddo
    endif
!    fxyzu(4,i) = 0.
    if (vxyzu(4,i) < 0. .or. isnan(vxyzu(4,i))) then
       print *,i,vxyzu(4,i)
       stop ' u<0'
    endif
 enddo
! !$omp end parallel do

end subroutine implicit_wind_cooling

!----------------------------------------------------------------
!+
!  Driver for the cooling function
!+
!----------------------------------------------------------------
subroutine calc_cooling (icooling, mass_per_H, Tgas, Trad, rho, Qcool, dQcool_dlnT)

  integer, intent(in) :: icooling
  real, intent(in) :: mass_per_H
  real, intent(in) :: Tgas, Trad, rho
  real, intent(out) :: Qcool, dQcool_dlnT

  Qcool = 0.d0
  dQcool_dlnT = 0.d0
  if (icooling == 1 .or. mod(icooling,10) == 3 .or. icooling == 11) then
     call cooling_neutral_hydrogen(mass_per_H, Tgas, rho, Qcool, dQcool_dlnT)
  end if
  if (icooling == 2 .or. mod(icooling,10) == 3 .or. icooling == 12) then
     call cooling_Bowen_relaxation(Tgas, Trad, rho, Qcool, dQcool_dlnT)
  end if
  if (icooling >= 10) then
     call cooling_radiative_relaxation(Tgas, Trad, Qcool, dQcool_dlnT)
  end if
end subroutine calc_cooling

!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_Bowen_relaxation(T, Trad, rho, Qcool, dQcool_dlnT)
  real, intent(in) :: T, Trad, rho
  real, intent(inout) :: Qcool, dQcool_dlnT
  real :: fac

  fac = specific_energy_to_T_ratio*rho/ Cprime
  Qcool = Qcool + fac*(Trad-T)
  dQcool_dlnT = dQcool_dlnT - fac*T
end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  Woitke (2006 A&A) cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_radiative_relaxation(T, Trad, Qcool, dQcool_dlnT)
  real, intent(in) :: T, Trad
  real, intent(inout) :: Qcool, dQcool_dlnT
  real :: fac

  fac = 4.*kappa_gas*usteboltz
  Qcool = Qcool + fac*(Trad**4-T**4)
  dQcool_dlnT = dQcool_dlnT - 4.*fac*T**3
end subroutine cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to neutral H (Spitzer)
!+
!-----------------------------------------------------------------------
subroutine cooling_neutral_hydrogen( mass_per_H, T, rho, Qcool, dQcool_dlnT)
  use units, only: umass, utime,udist
  real, intent(in) :: mass_per_H
  real, intent(in) :: T, rho
  real, intent(inout) :: Qcool, dQcool_dlnT

  real, parameter :: f = 0.2d0
  real :: eps_e, Q_H0, fQ

  if (T > 3000.d0) then
     fQ = utime**2*udist/umass
     call calc_eps_e(T, eps_e)
     Q_H0  = -f*fQ*7.3d-19 * eps_e * exp(-118400.d0/T) *rho / (mass_per_H)**2
     Qcool = Qcool + Q_H0
     dQcool_dlnT = dQcool_dlnT + 118400.d0/T * Q_H0
  endif
end subroutine cooling_neutral_hydrogen

!-----------------------------------------------------------------------
!+
!  compute electron equilibrium abundance (Palla et al 1983)
!+
!-----------------------------------------------------------------------
subroutine calc_eps_e(T, eps_e)
! used for cooling_neutral_hydrogen_radiation
  real, intent(in) :: T
  real, intent(out) :: eps_e

  real :: k1, k2, k3, k8, k9, p, q

!  if (T > 3000.) then
     k1 = 1.88d-10 / T**6.44d-1
     k2 = 1.83d-18 * T
     k3 = 1.35d-9
     k8 = 5.80d-11 * sqrt(T) * exp(-1.58d5/T)
     k9 = 1.7d-4 * k8

     p = .5d0*k8/k9
     q = k1*(k2+k3)/(k3*k9)

     eps_e = (p + sqrt(q+p**2))/q
!  else
!     eps_e = 0.d0
!  endif
end subroutine calc_eps_e

end module bowen_dust
