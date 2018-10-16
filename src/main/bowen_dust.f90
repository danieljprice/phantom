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

 public :: radiative_acceleration, pulsating_bowen_wind_profile, bowen_init
 logical, parameter :: use_alpha_wind = .true.

 private
 integer, parameter :: N = 1024
 integer, parameter :: wind_emitting_sink = 1
 logical, parameter :: verbose = .false.
 real, parameter :: alpha_wind = 1.d0
 real :: kappa, kmax, L, c_light, specific_energy_to_T_ratio, Cprime, Tcond, delta, Teff,&
      Reff, omega_osc, deltaR_osc

contains

!-----------------------------------------------------------------------
!+
!  Compute the radius, velocity and temperature of a sphere as a function of its local time
!+
!-----------------------------------------------------------------------
subroutine pulsating_bowen_wind_profile(local_time, r, v, u, rho, e, sphere_number,&
       wind_mass_rate,wind_injection_radius,wind_velocity,wind_osc_vamplitude,&
       wind_osc_period,shift_spheres,central_star_mass,time_between_spheres,&
       wind_temperature,wind_gamma)
 use part,    only: nptmass, xyzmh_ptmass
 use physcon, only: pi
 integer, intent(in) :: sphere_number
 real, intent(in) :: local_time,wind_injection_radius,wind_velocity,wind_osc_vamplitude,&
       wind_mass_rate,wind_osc_period,shift_spheres,central_star_mass,time_between_spheres,&
       wind_temperature,wind_gamma
 real, intent(out) :: r, v, u, rho, e

 real :: GM
 real :: base_time,base_radius, base_velocity, base_temperature, acceleration, rad_acc

 base_time = (sphere_number-shift_spheres)*time_between_spheres
 base_velocity = wind_velocity + wind_osc_vamplitude* cos(2.*pi*base_time/wind_osc_period)
 base_radius = wind_injection_radius + wind_osc_vamplitude*wind_osc_period/(2.*pi)*sin(2.*pi*base_time/wind_osc_period)
 base_temperature = wind_temperature
 acceleration = -central_star_mass / base_radius**2

 if (use_alpha_wind) then
    rad_acc = -alpha_wind*acceleration
 else
    call guess_acceleration(base_radius, Reff, rad_acc)
 endif
 acceleration = acceleration + rad_acc

 r = base_radius + local_time*base_velocity + acceleration*local_time**2/2.
 v = base_velocity + acceleration*local_time
 u = base_temperature * specific_energy_to_T_ratio
 rho = wind_mass_rate / (4.*pi*r**2*v)
 if (nptmass == 0) then
    GM = 0.
 else
    GM = xyzmh_ptmass(4,wind_emitting_sink)
 endif
 e = .5*v**2 - GM/r + wind_gamma*u
 !if (sphere_number<3320) print '("$$",i4,i2,8(1x,es12.5))',sphere_number,shift_spheres,base_radius,base_velocity,local_time&
 !     ,time_between_spheres,acceleration,r,v

end subroutine

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
!     * R_star: radius of the central star
!   Outputs:
!     * rad_acceleration(3,npart): acceleration vectors
!+
!-----------------------------------------------------------------------
subroutine radiative_acceleration(npart, xyzh, vxyzu, dt, fxyzu, time)
 use io,       only:iprint
 use part,     only:rhoh,xyzmh_ptmass,massoftype,igas
 use physcon,  only:pi
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(in)    :: dt,time
 real,    intent(inout) :: fxyzu(:,:)

 real :: O(3), h(npart), r(3, npart), d(npart), dmin, dmax, R_star, part_mass
 real :: position(npart), distance2(npart)
 integer :: nfound, found(npart)
 real :: rho(2*N+1), rho_over_r2(2*N+1)
 real :: OR(N), OR2(N), tau_prime(N)
 real :: Teq(N), kd(N), a(N), a_over_OR(N), kap(N)
 real :: Teq_part(npart), a_over_d_part(npart)
 real :: dQ_dt,fgrav
 integer :: i

 O = xyzmh_ptmass(1:3,1)
 !R_star = xyzmh_ptmass(5,1) should be the actual stellar radius
 R_star = Reff + deltaR_osc*sin(omega_osc*time)
 kap = 0.
 part_mass = massoftype(igas)
 if (npart == 0) then
    dmin = 0.
    dmax = 0.
 else
    h = xyzh(4,1:npart)
    call center_star(npart, xyzh, O, r, d, dmin, dmax)
    if (use_alpha_wind) then
       fgrav = alpha_wind*xyzmh_ptmass(4,wind_emitting_sink)
       do i=1,npart
          if (xyzh(4,i)  >  0.) then
             !if (i<4) print '(6(1x,es12.4))',fxyzu(1:3,i),fgrav/d(i)**3*r(1:3,i)
             fxyzu(1:3,i) = fxyzu(1:3,i) + fgrav/d(i)**3*r(1:3,i)
             vxyzu(4,i) = Teff*specific_energy_to_T_ratio
             fxyzu(4,i) = 0.d0
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
!  Easy routine to estimate the acceleration due to
!  radiation pressure at a given radius.
!+
!-----------------------------------------------------------------------
subroutine guess_acceleration(rinject, r_star, acceleration)
 use io,       only:iprint
 use physcon,  only: pi
 real, intent(in)  :: rinject, r_star
 real, intent(out) :: acceleration

 real :: OR2(1), tau_prime(1), Teq(1), kd(1), a(1)
 integer, parameter :: N=1

 tau_prime = 0.
 OR2 = rinject**2
 call calculate_Teq(N, r_star, tau_prime, OR2, Teq)
 call calculate_kd(N, Teq, kd)
 !call calculate_acceleration(N, OR2, kd, a)
 !acceleration = a(1)
 acceleration = L/(4.*pi*c_light) * kd(1)/OR2(1)

 if (verbose) then
    write(iprint,*) " * Teq(rinject) : ", Teq(1),tau_prime
    write(iprint,*) " * kd(rinject)  : ", kd(1)
    write(iprint,*) " * rinject      : ", sqrt(OR2(1)),rinject
    write(iprint,*) " * a(rinject)   : ", acceleration
    write(iprint,*) ""
 endif

end subroutine guess_acceleration

!-----------------------------------------------------------------------
!+
!  Convert parameters into simulation units.
!+
!-----------------------------------------------------------------------
subroutine bowen_init(u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,&
       bowen_Cprime, bowen_Tcond, bowen_delta,bowen_Teff,wind_osc_vamplitude,wind_osc_period)
 use physcon,     only: solarl,c,steboltz,pi
 use units,       only: udist, umass, utime
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

end subroutine

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
end subroutine

!-----------------------------------------------------------------------
!+
!  Calculates the acceleration due to radiative pressure on dust grains, along the half-line
!   Inputs:
!     * N: number of steps
!     * L: luminosity of central star
!     * c: speed of light
!     * OR2(N): squared distance
!     * kd(N): density divided by squared distance
!   Outputs:
!     * a(N): acceleration along the half-line
!+
!-----------------------------------------------------------------------
 ! subroutine calculate_acceleration(N, L, c, OR2, kd, a)
 !   use physcon, only:pi
 !   integer, intent(in)  :: N
 !   real,    intent(in)  :: L, c, OR2(N), kd(N)
 !   real,    intent(out) :: a(N)

 !   a = L/(4.*pi*c) * kd/OR2

 ! end subroutine

!-----------------------------------------------------------------------
!+
!  Calculates kd, the dust opacity along the half-line
!   Inputs:
!     * N: number of steps
!     * kmax: maximum dust mass ratio (chosen so that radiative pressure is of the same order as gravity)
!     * Tcond: condensation temperature of dust grains
!     * delta: parameter for the width of the condensation temperature range
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
end subroutine

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
    !Teq = ((Teff**4/2.)*(1.-sqrt(1.-(R_star**2/OR2))) + (Teff**4*3./4.)*tau_prime)**(1./4.)
 end where
end subroutine

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
end subroutine

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

end subroutine

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
end subroutine

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
end subroutine

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
end subroutine

end module bowen_dust
