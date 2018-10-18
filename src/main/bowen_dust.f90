!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
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
!  OWNER: Lionel
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

 public :: radiative_acceleration, bowen_init
 logical, parameter :: use_alpha_wind = .true.

 private
 integer, parameter :: N = 1024
 integer, parameter :: wind_emitting_sink = 1
 logical, parameter :: verbose = .false.
 real, parameter :: alpha_wind = 1.d0
 real :: kappa, kmax, L, c_light, specific_energy_to_T_ratio, Cprime, Tcond, delta, Teff,&
      Reff, omega_osc, deltaR_osc
 integer :: nwall_particles

contains

!-----------------------------------------------------------------------
!+
!  Computes the acceleration due to radiation pressure on dust, on all SPH particles
!  Here we assume that optical depth is roughly spherically symmetric around the central star
!  Execution time: a few milliseconds with one million particles
!   Inputs:
!     * npart: number of particles
!     * xyzh(4,maxp): big phantom's array of particles
!     * O(3): location of the central star
!     * part_mass: mass of the particles (assumed constant)
!   Outputs:
!     * rad_acceleration(3,npart): acceleration vectors
!+
!-----------------------------------------------------------------------
subroutine radiative_acceleration(npart, xyzh, vxyzu, dt, fxyzu, time)
 use io,       only:iprint
 use part,     only:isdead_or_accreted,rhoh,xyzmh_ptmass,massoftype,igas
 use physcon,  only:pi
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(in)    :: dt,time
 real,    intent(inout) :: fxyzu(:,:)

 real :: O(3), h(npart), r(3, npart), d(npart), dmin, dmax, part_mass
 real :: position(npart), distance2(npart)
 integer :: nfound, found(npart)
 real :: rho(2*N+1), rho_over_r2(2*N+1)
 real :: OR(N), OR2(N), tau_prime(N)
 real :: Teq(N), kd(N), a(N), a_over_OR(N), kap(N)
 real :: Teq_part(npart), a_over_d_part(npart)
 real :: dQ_dt,Mstar,Rforce,R_star,alpha_surface
 integer :: i

 O = xyzmh_ptmass(1:3,1)
 !R_star = xyzmh_ptmass(5,1) should be the actual stellar radius
 R_star = Reff + deltaR_osc*sin(omega_osc*time)
 Rforce = R_star + deltaR_osc
 kap = 0.
 part_mass = massoftype(igas)
 if (npart == 0) then
    dmin = 0.
    dmax = 0.
 else
    h = xyzh(4,1:npart)
    call center_star(npart, xyzh, O, r, d, dmin, dmax)
    if (use_alpha_wind) then
       Mstar = xyzmh_ptmass(4,wind_emitting_sink)
       do i=1,npart
          if (.not.isdead_or_accreted(xyzh(4,i))) then
             !repulsive force, Ayliffe & Bate 2010, 408, 876
             if (d(i) < 2.*Rforce) then
                alpha_surface = ((2.*Rforce-d(i))/Rforce)**4
             else
                alpha_surface = 0.
             endif
             if (i <= nwall_particles) then
                fxyzu(:,i) = 0.
             else
                fxyzu(1:3,i) = fxyzu(1:3,i) + Mstar*max(alpha_wind,alpha_surface)/d(i)**3*r(1:3,i)
             endif
           endif
        enddo
     else
       if (abs((dmin-dmax)/dmax)  <  1.0d-10) then
          rho_over_r2 = 0.
       else
          !call project_on_line(npart, r, 0.d0, 0.d0, position, distance2)
          call project_on_z(npart, r, position, distance2)
          call select_particles(npart, distance2, h, nfound, found)
          call density(npart, position, distance2, h, part_mass, nfound, found,&
                    -dmax, dmax, R_star, 2*N+1, rho, rho_over_r2)
       endif
       call calculate_tau_prime(N, dmax, R_star, kappa, kap, rho_over_r2, OR, tau_prime)
       OR2 = OR**2
       call calculate_Teq(N, R_star, tau_prime, OR2, Teq)
       call calculate_kd(N, Teq, kd)
       !do i=1,100
       !   print '(i4,8(1x,es12.4))',i,Teq(i),kmax,kappa,kd(i),OR(i),0.5*(1.-sqrt(1.-(R_star**2/OR2(i)))),0.75*tau_prime(i)
       !enddo
       a = L/(4.*pi*c_light) * kd/OR2
       a_over_OR = a/OR
       call interpolate_on_particles(npart, d, N, dmax, a_over_OR, a_over_d_part)
       call interpolate_on_particles(npart, d, N, dmax, Teq, Teq_part)
       do i=1,npart
          if (h(i)  >  0.) then
             fxyzu(1,i) = fxyzu(1,i) + a_over_d_part(i) * r(1,i)
             fxyzu(2,i) = fxyzu(2,i) + a_over_d_part(i) * r(2,i)
             fxyzu(3,i) = fxyzu(3,i) + a_over_d_part(i) * r(3,i)
             if (dt  >=  Cprime/rhoh(xyzh(4,i),part_mass)) then   ! assume thermal equilibrium
                vxyzu(4,i) = Teq_part(i)*specific_energy_to_T_ratio
                fxyzu(4,i) = 0.d0
             else
                dQ_dt = -(vxyzu(4,i) - Teq_part(i)*specific_energy_to_T_ratio) &
                * (rhoh(xyzh(4,i),part_mass)/Cprime)
                fxyzu(4,i) = fxyzu(4,i) + dQ_dt
             endif
          endif
          !if (i<20) print '("a",i4,18(1x,es12.4))',i,fxyzu(1:4,i),xyzh(1:3,i)
       enddo
    endif
 endif
 !stop

 if (verbose) then
    write(iprint,*) "##############"
    write(iprint,*) "# BOWEN DUST #"
    write(iprint,*) "##############"
    write(iprint,*) " * star radius: ", R_star
    write(iprint,*) " * min(Teq):", minval(Teq)
    write(iprint,*) " * max(Teq):", maxval(Teq)
    write(iprint,*) " * amin:", minval(a)
    write(iprint,*) " * amax:", maxval(a)
    write(iprint,*) ""
 endif
end subroutine radiative_acceleration


!-----------------------------------------------------------------------
!+
!  Convert parameters into simulation units.
!+
!-----------------------------------------------------------------------
subroutine bowen_init(u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,&
       bowen_Cprime, bowen_Tcond, bowen_delta,bowen_Teff,wind_osc_vamplitude,wind_osc_period,nwall)
 use physcon,     only: solarl,c,steboltz,pi
 use units,       only: udist, umass, utime
 integer, intent(in) :: nwall
 real, intent(in)  :: u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,&
         bowen_Cprime, bowen_Tcond, bowen_delta, bowen_Teff, wind_osc_vamplitude,wind_osc_period

 kappa = bowen_kappa / (udist**2/umass)
 kmax = bowen_kmax / (udist**2/umass)
 L = bowen_L /(umass*udist**2/utime**3)
 c_light = c / (udist/utime)
 specific_energy_to_T_ratio = u_to_temperature_ratio/(udist/utime)**2
 Cprime = bowen_Cprime / (umass*utime/udist**3)
 Tcond = bowen_Tcond
 delta = bowen_delta
 Teff  = bowen_Teff
 Reff = sqrt(bowen_L/(4.*pi*steboltz*bowen_Teff**4))/udist
 omega_osc = 2.*pi/wind_osc_period
 deltaR_osc = wind_osc_vamplitude/omega_osc
 nwall_particles = nwall

end subroutine bowen_init

!-----------------------------------------------------------------------
!+
!  Interpolates a quantity computed on the discretized half line on all SPH particles, assuming spherical symmetry
!  Method: basic linear interpolation
!   Inputs:
!     * npart: number of particles
!     * d(npart): distance of particles
!     * N: number of steps
!     * dmax: maximum distance
!     * quantity(N): the quantity to interpolate
!   Outputs:
!     * output(npart): interpolated values
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
 do i=1,npart
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

 kd = kmax/(1. + exp((Teq-Tcond)/delta))
end subroutine calculate_kd

!-----------------------------------------------------------------------
!+
!  Calculates Teq, the radiative equilibrium temperature along the half-line
!   Inputs:
!     * N: number of steps
!     * Teff: effective temperature of central star
!     * R_star: radius of central star
!     * tau_prime(N): optical depth
!     * OR2(N): squared distance
!   Outputs:
!     * Teq(N): radiative equilibrium temperature
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
!  Calculates tau_prime, a rough estimation of the optical depth along the half-line
!  Method: two-front trapezoidal rule, assuming that rho_over_r2 is discretized on a symmetric manner
!          the output is the mean value of the optical depth along the two half-lines
!   Inputs:
!     * N: number of  (on the half-line)
!     * dmax: maximum distance
!     * R_star: radius of central star
!     * kappa: gas opacity (assumed constant everywhere!)
!     * rho_over_r2(2N+1): density divided by squared distance (on the full line)
!   Outputs:
!     * OR(N): location of the steps on the half-line
!     * tau_prime(N): mean optical depth along the half-line
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
!  Computes density along the line
!   Inputs:
!     * npart: number of particles
!     * position(npart): location of the projection of the particles on the line (can be either positive or negative)
!     * distance2(npart): squared distance between the particles and the line
!     * part_h(npart): smoothing length of particles
!     * part_mass: mass of particles (assumed the same for all particles)
!     * nfound: number of interesting particles
!     * found(1:nfound): list of interesting particles
!     * rmin, rmax: boundaries of the line (can be either positive or negative but rmax>rmin)
!     * r_star: radius of a central star
!     * N: number of steps along the line
!   Outputs:
!     * rho(N): density along the line (at N locations between rmin and rmax)
!     * rho_over_r2(N): rho/r² outside of the central star
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

!-----------------------------------------------------------------------
!+
!  Select interesting particles (whose interaction zone intersects the line)
!   Inputs:
!     * npart: number of particles
!     * distance2(npart): squared distance to the line
!     * part_h(npart): smoothing length of particle's SPH kernel
!   Outputs:
!     * nfound: number of particles found
!     * found(1:nfound): index of interesting particles
!+
!-----------------------------------------------------------------------
subroutine select_particles(npart, distance2, part_h, nfound, found)
 integer, intent(in)  :: npart
 real,    intent(in)  :: distance2(npart), part_h(npart)
 integer, intent(out) :: nfound, found(npart)

 integer :: i
 real :: h

 nfound = 0
 do i=1,npart
    h = part_h(i)
    if (distance2(i)  <  4.*h**2) then
       nfound = nfound + 1
       found(nfound) = i
    endif
 enddo
end subroutine select_particles

!-----------------------------------------------------------------------
!+
!  Project SPH particles on a vertical line
!   Inputs:
!     * npart: number of particles
!     * r(3,npart): location of the particles
!   Outputs:
!     * position(npart): z coordinates of the particles
!     * distance2(npart): squared distance between the particles and z axis
!+
!-----------------------------------------------------------------------
subroutine project_on_z(npart, r, position, distance2)
 integer, intent(in)  :: npart
 real,    intent(in)  :: r(3,npart)
 real,    intent(out) :: position(npart), distance2(npart)

 position = r(3,:)
 distance2 = r(1,:)**2+r(2,:)**2
end subroutine project_on_z

!-----------------------------------------------------------------------
!+
!  Project SPH particles on a line
!   Inputs:
!     * npart: number of particles
!     * r(3,npart): location of the particles
!     * theta, phi: direction of the line (theta: colatitude; phi: longitude)
!   Outputs:
!     * position(npart): location of the projection of the particles on the line
!     * distance2(npart): squared distance between the particles and the line
!+
!-----------------------------------------------------------------------
!  subroutine project_on_line(npart, r, theta, phi, position, distance2)
!    integer, intent(in)  :: npart
!    real,    intent(in)  :: r(3,npart), theta, phi
!    real,    intent(out) :: position(npart), distance2(npart)

!    integer :: i
!    real :: u(3), OP(3), HP(3), p

 ! Director vector
!    u(1) = sin(theta)*cos(phi)
!    u(2) = sin(theta)*sin(phi)
!    u(3) = cos(theta)

!    do i=1,npart
!      OP = r(:,i)
!      p = dot_product(OP,u)
!      position(i) = p
!      HP = OP - p*u
!      distance2(i) = dot_product(HP,HP)
!    enddo
!  end subroutine


!-----------------------------------------------------------------------
!+
! Computes the location of SPH particles in a new frame
!  Inputs:
!    * npart: number of particles
!    * xyzh(4,maxp): big phantom's array of particles
!    * O(3): location of the new frame
!  Outputs:
!    * r(3,npart): location of SPH particles in the new frame
!    * d(npart): distance from O to SPH particles
!    * dmin: minimum of d
!    * dmax: maximum of d
!+
!-----------------------------------------------------------------------
subroutine center_star(npart, xyzh, O, r, d, dmin, dmax)
 use dim,      only:maxp
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(4,maxp), O(3)
 real,    intent(out) :: r(3,npart), d(npart), dmin, dmax

 integer :: i

 do i=1,npart
    r(:,i) = xyzh(1:3,i)-O
 enddo
 d = sqrt(r(1,:)**2 + r(2,:)**2 + r(3,:)**2)
 dmin = minval(d)
 dmax = maxval(d)
end subroutine center_star

end module bowen_dust
