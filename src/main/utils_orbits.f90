!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module orbits
!
! Utility routines for orbital dynamics and converting between
! various orbital elements
!
! Includes helper routines for GR orbits
!
! :References:
!   Eggleton (1983) ApJ 268, 368-369 (ref:eggleton83)
!   https://en.wikipedia.org/wiki/Orbital_elements
!   https://en.wikipedia.org/wiki/Eccentricity_vector
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 real, parameter, public :: pi = 4.*atan(1.)
 real, parameter, public :: deg_to_rad = pi/180.
 real, parameter, public :: rad_to_deg = 180.0/pi

 public :: Rochelobe_estimate, L1_point

 ! routines to return orbital elements from position and velocity
 public :: get_semimajor_axis
 public :: get_orbital_period
 public :: get_eccentricity,get_eccentricity_vector
 public :: get_pericentre_distance
 public :: get_inclination
 public :: get_argument_of_periapsis
 public :: get_longitude_of_ascending_node
 public :: get_true_anomaly
 public :: get_orbital_elements
 public :: get_specific_energy
 public :: orbit_is_parabolic

 ! generic interface for semimajor axis
 interface get_semimajor_axis
  module procedure get_semimajor_axis_from_period,get_semimajor_axis_from_posvel,&
                   get_semimajor_axis_from_energy
 end interface get_semimajor_axis

 ! generic interface for orbital period
 interface get_orbital_period
  module procedure get_period_from_posvel,get_period_from_semimajor_axis
 end interface get_orbital_period

 ! general interface for eccentricity vector
 interface get_eccentricity_vector
  module procedure get_eccentricity_vector_posvel,get_eccentricity_vector_sinks
 end interface get_eccentricity_vector

 ! generic interface for scalar eccentricity
 interface get_eccentricity
  module procedure get_eccentricity_posvel,get_eccentricity_posvel_scalar
 end interface get_eccentricity

 ! specific energy, angular momentum, escape velocity
 interface get_specific_energy
  module procedure get_specific_energy,get_specific_energy_posvel
 end interface get_specific_energy

 public :: get_specific_energy_gr
 public :: get_angmom_vector,get_angmom_unit_vector,get_mean_angmom_vector
 public :: escape,get_escape_velocity_at_infinity

 ! conversion routines between true anomaly, mean anomaly,
 ! eccentric anomaly, time to pericentre
 public :: get_T_flyby_hyp, get_T_flyby_par
 public :: get_E,get_E_from_mean_anomaly,get_E_from_true_anomaly
 public :: get_time_to_separation,get_true_anomaly_from_separation
 public :: get_time_between_true_anomalies

 ! convert between different sets of elements
 public :: convert_flyby_to_elements,convert_posvel_to_flyby
 public :: get_dx_dv_ptmass

 interface get_orbital_elements
  module procedure get_orbital_elements,get_elements_mean_anomaly,get_orbparams
 end interface get_orbital_elements

 ! time derivative of semi-major axis and eccentric anomaly
 public :: get_a_dot

 ! routines for GR orbits
 public :: isco_kerr, refine_velocity

 private

contains
!------------------------------------
! Compute estimate of the Roche Lobe
! Eggleton (1983) ApJ 268, 368-369
!------------------------------------
real function Rochelobe_estimate(m1,m2,sep)
 real, intent(in) :: m1,m2,sep
 real :: q,q13,q23

 if (m1 > 0. .and. m2 > 0.) then
    q = m2/m1
    q13 = q**(1./3.)
    q23 = q13*q13
    Rochelobe_estimate = sep * 0.49*q23/(0.6*q23 + log(1. + q13))
 else
    Rochelobe_estimate = sep
 endif

end function Rochelobe_estimate

!---------------------------------------------
! Find first Lagrange point (L1)
! via Newton-Raphson solution of quintic
!
! INPUT: mass ratio of binary
! OUTPUT: L1 point, as distance from primary
!---------------------------------------------
real function L1_point(qinv)
 real, intent(in) :: qinv
 real :: fL, dfL, dL, L, q11

 q11 = 1./(1.+qinv)
 L = 0.5 + 0.2222222*log10(qinv)

 dL = 1.e7
 do while (abs(dL)>1.e-6)
    fL = qinv/L**2- 1./(1.-L)**2 - (1.+qinv)*L + 1.
    dfL=-2*qinv/L**3 - 2./(1.-L)**3 - (1.+qinv)
    dL = -fL/(dfL*L)
    L = L*(1.+dL)
 enddo

 L1_point = L

end function L1_point

!----------------------------------------------------------------
!+
!  This function calculates semimajor axis of the orbit
!  from the position and velocity vectors
!+
!----------------------------------------------------------------
real function get_semimajor_axis_from_posvel(mu,dx,dv) result(a)
 real, intent(in) :: mu
 real, intent(in) :: dx(3),dv(3)
 real :: e,rp

 e = get_eccentricity(mu,dx,dv)
 rp = get_pericentre_distance(mu,dx,dv)

 if (orbit_is_parabolic(e)) then
    a = rp  ! for parabolic orbit return pericentre distance
 else
    a = rp / (1.-e) ! return semi-major axis for elliptical or hyperbolic orbit
 endif

end function get_semimajor_axis_from_posvel

!-------------------------------------------------------------
! Function to determine the semi-major axis given the period
!-------------------------------------------------------------
real function get_semimajor_axis_from_period(mu,period) result(a)
 real, intent(in) :: mu,period

 a = (mu*(period/(2.*pi))**2)**(1./3.)

end function get_semimajor_axis_from_period

!----------------------------------------------------------------
!+
!  Compute the semi-major axis a = (r*mu)/(2*mu-r*v2)
!  from the separation, mass and square of velocity
!+
!----------------------------------------------------------------
real function get_semimajor_axis_from_energy(mu,r,v2) result(a)
 real, intent(in) :: mu,r,v2
 real :: denom

 ! a = (r*mu)/(2.*mu-r*v2) for elliptical or hyperbolic orbit
 denom = 2.*mu-r*v2
 if (abs(denom) > epsilon(denom)) then
    a = (r*mu)/denom
 else ! parabolic orbit
    a = huge(a)
 endif

end function get_semimajor_axis_from_energy

!----------------------------------------------------------------
!+
!  This function calculates period of the orbit
!  from the position and velocity vectors
!+
!----------------------------------------------------------------
real function get_period_from_posvel(mu,dx,dv) result(period)
 real, intent(in) :: mu
 real, intent(in) :: dx(3),dv(3)
 real :: a

 a = get_semimajor_axis(mu,dx,dv)
 period = get_orbital_period(mu,a)

end function get_period_from_posvel

!-------------------------------------------------------------
! Function to determine the period given the semi-major axis
!-------------------------------------------------------------
function get_period_from_semimajor_axis(mu,a) result(period)
 real, intent(in) :: mu,a
 real :: period

 ! Kepler's 3rd law
 period = 2.*pi*sqrt(abs(a)**3/mu)

end function get_period_from_semimajor_axis

!----------------------------------------------------------------
!+
!  Compute the pericentre distance from the position and velocity
!  in a way that works for all orbit types
!+
!----------------------------------------------------------------
real function get_pericentre_distance(mu,dx,dv) result(rperi)
 real, intent(in) :: mu,dx(3),dv(3)
 real :: h,e

 h = get_angmom(dx,dv)
 e = get_eccentricity(mu,dx,dv)
 ! recover rp from angular momentum: rp = h^2/(mu*(1+e))
 rperi = h**2/(mu*(1.+e))

end function get_pericentre_distance

!----------------------------------------------------------------
!+
!  Compute the specific energy e = 0.5*v2 - mu/r
!+
!----------------------------------------------------------------
real function get_specific_energy(mu,v2,r) result(energy)
 real, intent(in) :: mu,v2,r

 energy = 0.5*v2 - mu/r

end function get_specific_energy

!----------------------------------------------------------------
!+
!  Specific energy from relative position and velocity vectors
!+
!----------------------------------------------------------------
real function get_specific_energy_posvel(mu,dx,dv) result(energy)
 real, intent(in) :: mu,dx(3),dv(3)

 energy = 0.5*dot_product(dv,dv) - mu/sqrt(dot_product(dx,dx))

end function get_specific_energy_posvel

!----------------------------------------------------------------
!+
!  Compute the time derivative of the semi-major axis
!  adot = 2*(mu2*v+r2*v*acc)/((2*mu-r*v2)**2)
!+
!----------------------------------------------------------------
real function get_a_dot(mu,r2,r,v2,v,acc) result(adot)
 real, intent(in) :: mu,r2,r,v2,v,acc

 adot = 2.*(mu*mu*v+r2*v*acc)/((2.*mu-r*v2)**2)

end function get_a_dot

!----------------------------------------------------------------
!+
!  This function calculates eccentricity vector of the orbit
!  from the relative position and velocity
!
!  https://en.wikipedia.org/wiki/Eccentricity_vector
!+
!----------------------------------------------------------------
function get_eccentricity_vector_posvel(mu,dx,dv) result(e_vec)
 real, intent(in)   :: mu,dx(3),dv(3)
 real :: e_vec(3)
 real :: vcrossh_vec(3)
 real :: r

 vcrossh_vec = vcrossh(dx,dv)
 r = sqrt(dot_product(dx,dx))

 ! eccentricity vector = (rdot cross h )/(G(m1+m2)) - rhat
 e_vec(:) = vcrossh_vec(:)/mu - dx/r

end function get_eccentricity_vector_posvel

!----------------------------------------------------
! interface to above assuming two sink particles
!----------------------------------------------------
function get_eccentricity_vector_sinks(xyzmh_ptmass,vxyz_ptmass,i1,i2) result(e_vec)
 real,    intent(in) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: i1, i2
 real :: e_vec(3),mu,dx(3),dv(3)

 if (i1 > 0 .and. i2 > 0) then
    mu = xyzmh_ptmass(4,i1) + xyzmh_ptmass(4,i2)
    dx = xyzmh_ptmass(1:3,i2) - xyzmh_ptmass(1:3,i1)
    dv = vxyz_ptmass(1:3,i2) - vxyz_ptmass(1:3,i1)
    e_vec = get_eccentricity_vector(mu,dx,dv)
 else
    e_vec = 0.
 endif

end function get_eccentricity_vector_sinks

!----------------------------------------------------------------
!+
!  This function calculates eccentricity of the orbit
!  from the relative position and velocity
!+
!----------------------------------------------------------------
real function get_eccentricity_posvel(mu,dx,dv) result(e)
 real, intent(in)   :: mu,dx(3),dv(3)
 real :: e_vec(3)

 e_vec = get_eccentricity_vector(mu,dx,dv)
 e = sqrt(dot_product(e_vec,e_vec))

end function get_eccentricity_posvel

!----------------------------------------------------------------
!+
!  Compute the eccentricity e = sqrt(ex**2+ey**2+ez**2)
!  this version takes scalar positions and velocities
!  and a scalar separation to avoid recomputing square roots
!+
!----------------------------------------------------------------
real function get_eccentricity_posvel_scalar(mu,r,x,y,z,vx,vy,vz) result(e)
 real, intent(in) :: mu,r,x,y,z,vx,vy,vz
 real :: ex,ey,ez
 real :: hx,hy,hz

 hx = y*vz-z*vy
 hy = z*vx-x*vz
 hz = x*vy-y*vx

 ex = (vy*hz-vz*hy)/mu - x/r
 ey = (vz*hx-vx*hz)/mu - y/r
 ez = (vx*hy-hx*vy)/mu - z/r

 e = sqrt(ex**2+ey**2+ez**2)

end function get_eccentricity_posvel_scalar

!------------------------------------------------------------
!+
!  Whether orbit is parabolic (e=1 with some small tolerance)
!+
!------------------------------------------------------------
logical function orbit_is_parabolic(e)
 real, parameter :: tol_eccentricity = 1.e-12
 real, intent(in) :: e

 if (abs(e-1.0) < tol_eccentricity) then
    orbit_is_parabolic = .true.
 else
    orbit_is_parabolic = .false.
 endif

end function orbit_is_parabolic

!----------------------------------------------------------------
!+
!  Compute the specific angular momentum vector h = r x v
!+
!----------------------------------------------------------------
function get_angmom_vector(dx,dv) result(h_vec)
 real, intent(in) :: dx(3),dv(3)
 real :: h_vec(3)

 h_vec = cross_product(dx,dv)

end function get_angmom_vector

!----------------------------------------------------------------
!+
!  Compute the specific angular momentum unit vector h_hat = h/|h|
!+
!----------------------------------------------------------------
function get_angmom_unit_vector(dx,dv) result(h_vec)
 real, intent(in) :: dx(3),dv(3)
 real :: h_vec(3)

 h_vec = get_angmom_vector(dx,dv)
 h_vec = h_vec / sqrt(dot_product(h_vec,h_vec))

end function get_angmom_unit_vector

!----------------------------------------------------------------
!+
!  Magnitude of the specific angular momentum h = |h|
!+
!----------------------------------------------------------------
real function get_angmom(dx,dv) result(h)
 real, intent(in) :: dx(3),dv(3)
 real :: h_vec(3)

 h_vec = cross_product(dx,dv)
 h = sqrt(dot_product(h_vec,h_vec))

end function get_angmom

!----------------------------------------------------------------
!+
!  This function calculates v x h
!+
!----------------------------------------------------------------
function vcrossh(dx,dv)
 real, intent(in) :: dx(3),dv(3)
 real :: h(3),vcrossh(3)

 h = get_angmom_vector(dx,dv)
 vcrossh = cross_product(dv,h)

end function vcrossh

!----------------------------------------------------------------
!+
!  Compute the n vector = k cross h, vector along line of nodes
!+
!----------------------------------------------------------------
function get_line_of_nodes_vector(dx,dv) result(n_vec)
 real, intent(in) :: dx(3),dv(3)
 real, parameter :: k_vec(3) = (/0,0,1/)
 real :: h_vec(3),n_vec(3)

 h_vec = get_angmom_vector(dx,dv)

 ! n vector = k cross h, vector along line of nodes
 n_vec = cross_product(k_vec,h_vec)

end function get_line_of_nodes_vector

!-------------------------------------------------------------
! Function to find mean angular momentum vector from a list
! of positions and velocities
!-------------------------------------------------------------
function get_mean_angmom_vector(n,xyz,vxyz) result(l)
 integer, intent(in) :: n
 real,    intent(in) :: xyz(:,:),vxyz(:,:)
 real    :: l(3),li(3)
 integer :: i

 l = 0.
 do i=1,n
    li = cross_product(xyz(:,i),vxyz(:,i))
    l = l + li
 enddo
 l = l/real(n)

end function get_mean_angmom_vector

!---------------------------------------------------------------
!+
! Get eccentric anomaly given time since pericentre
!+
!---------------------------------------------------------------
subroutine get_E(period,ecc,deltat,E)
 real, intent(in)  :: period,ecc,deltat
 real, intent(out) :: E
 real :: mu,M_ref

 mu = 2.*pi/period
 M_ref = mu*deltat ! mean anomaly

 E = get_E_from_mean_anomaly(M_ref,ecc)

end subroutine get_E

!---------------------------------------------------------------
!+
! Get eccentric anomaly from mean anomaly (this function uses
! bisection to guarantee convergence, which is not guaranteed for
! small M or E)
!+
!---------------------------------------------------------------
real function get_E_from_mean_anomaly(M_ref,ecc) result(E)
 real, intent(in) :: M_ref,ecc
 real :: E_left,E_right,E_guess,M_guess
 real, parameter :: tol = 1.e-10

 ! first guess
 E_left = 0.
 E_right = 2.*pi
 E_guess = pi
 M_guess = M_ref - 2.*tol

 do while (abs(M_ref - M_guess) > tol)
    M_guess = get_mean_anomaly_from_E(E_guess,ecc)

    if (M_guess > M_ref) then
       E_right = E_guess
    else
       E_left = E_guess
    endif
    E_guess = 0.5*(E_left + E_right)
 enddo

 E = E_guess

end function get_E_from_mean_anomaly

!---------------------------------------------------------------
!+
!  Get mean anomaly from eccentric anomaly
!+
!---------------------------------------------------------------
real function get_mean_anomaly_from_E(E,ecc) result(M)
 real, intent(in) :: E,ecc

 if (orbit_is_parabolic(ecc)) then
    M = E + 1./3.*E**3
 elseif (ecc > 1.) then ! hyperbolic
    M = ecc*sinh(E) - E
 else    ! eccentric
    M = E - ecc*sin(E)
 endif

end function get_mean_anomaly_from_E

!---------------------------------------------------------------
!+
!  Get eccentric (or parabolic/hyperbolic) anomaly from true anomaly
!  https://space.stackexchange.com/questions/23128/design-of-an-elliptical-transfer-orbit/23130#23130
!+
!---------------------------------------------------------------
real function get_E_from_true_anomaly(nu,ecc) result(E)
 real, intent(in) :: nu  ! true anomaly in radians
 real, intent(in) :: ecc    ! eccentricity

 if (orbit_is_parabolic(ecc)) then
    E = tan(0.5*nu)
 elseif (ecc > 1.) then
    E = 2.*atanh(sqrt((ecc - 1.)/(ecc + 1.))*tan(0.5*nu))
 else ! ecc < 1.0
    E = atan2(sqrt(1. - ecc**2)*sin(nu),(ecc + cos(nu)))
 endif

end function get_E_from_true_anomaly

!------------------------------------------------------------
!+
!  Flyby time for hyperbolic trajectory:
!  returns the elapsed time symmetric about pericentre
!  i.e. twice the time to reach pericentre from the input
!  true anomaly, or the time to travel from -f to f
!+
!------------------------------------------------------------
function get_T_flyby_hyp(mu,ecc,f,a) result(T)
 real, intent(in) :: mu,ecc,f,a
 real :: T,h,M_h0,F0,nu

!--True anomaly in radians
 nu = f * deg_to_rad

!--Specific angular momentum
 h = sqrt(mu * a * (ecc**2 - 1.0))

!--Eccentric anomaly
 F0 = 2.0 * atanh(sqrt((ecc - 1.0) / (ecc + 1.0)) * tan(nu / 2.0))

!--Mean anomaly
 M_h0 = ecc * sinh(F0) - F0

!--Time of flight
 T = 2.0 * abs((h**3 / mu**2) * (1.0 / (ecc**2 - 1.0)**1.5) * M_h0)

end function get_T_flyby_hyp

!------------------------------------------------------------
!+
!  Flyby time for parabolic trajectory: see get_T_flyby_hyp
!+
!------------------------------------------------------------
function get_T_flyby_par(mu,dma,n0) result(T)
 real, intent(in) :: mu,dma,n0
 real :: T,nu,xi,yi,Di,Df,p

!--Semi-latus rectum
 p = 2.0 * dma

!--Initial position
 xi = -2.0 * sqrt(n0 - 1.0) * dma
 yi = dma * (1.0 - (xi / p)**2)

!--True anomaly
 nu = pi - atan(abs(xi / yi))

!--Barker's equation
 Di = tan(-nu / 2.0)
 Df = tan(nu / 2.0)

 T = 0.5 * sqrt(p**3 / mu) * (Df + (1.0 / 3.0) * Df**3 - Di - (1.0 / 3.0) * Di**3)

end function get_T_flyby_par

!------------------------------------------------------------
!+
!  Escape velocity at infinity for an unbound orbit
!+
!------------------------------------------------------------
real function get_escape_velocity_at_infinity(mu,dx,dv) result(v_esc)
 real, intent(in) :: mu,dx(3),dv(3)
 real :: energy

 energy = get_specific_energy(mu,dot_product(dv,dv),sqrt(dot_product(dx,dx)))
 if (energy > 0.) then
    v_esc = sqrt(2.*energy)
 else
    v_esc = -1.0
 endif

end function get_escape_velocity_at_infinity

!----------------------------------------------------------------
!+
!  Whether velocity is greater than escape velocity
!+
!----------------------------------------------------------------
logical function escape(vel,mu,r)
 real, intent(in) :: vel,mu,r
 real :: v_esc

 ! we require the position of object on orbit/ wrt central object, velocity on
 ! orbit is velocity of object wrt central obj, mu is G*mass of central object

 v_esc = sqrt(2.*mu/r)
 if (vel > v_esc) then
    print "(a,g0)",' Star has escaped: v/v_esc = ',vel/v_esc
    escape = .true.
 else
    escape = .false.
 endif

end function escape

!----------------------------------------------------------------
!+
!  Compute the inclination angle in degrees
!+
!----------------------------------------------------------------
real function get_inclination(dx,dv) result(inc)
 real, intent(in) :: dx(3),dv(3)
 real :: h(3)

 h = get_angmom_unit_vector(dx,dv)
 inc = acos(h(3))*rad_to_deg

end function get_inclination

!-------------------------------------------------------------
!+
!  Calculate true anomaly from separation for given orbital elements
!+
!-------------------------------------------------------------
real function get_true_anomaly_from_separation(a,e,r) result(f)
 real, intent(in) :: a, e, r
 real :: cos_f

 if (orbit_is_parabolic(e)) then
    ! Parabolic orbit
    cos_f = 2.0 * a / r - 1.0
 elseif (e > 0.) then    ! Elliptical orbit, e < 1.0 or hyperbolic orbit, e > 1.0 (a < 0)
    cos_f = (a * (1.0 - e*e) / r - 1.0) / e
 else ! for a circular orbit this routine is ambiguous, return f = 180.0
    cos_f = -1.0
 endif

 ! Ensure cos_f is within valid range
 cos_f = max(-1.0, min(1.0, cos_f))
 f = acos(cos_f)*rad_to_deg

end function get_true_anomaly_from_separation

!------------------------------------------------------------
!+
!  Get true anomaly (in degrees)from position and velocity
!+
!------------------------------------------------------------
real function get_true_anomaly(mu,dx,dv) result(f)
 real, intent(in) :: mu,dx(3),dv(3)
 real :: ecc_vec(3),h_hat(3)
 real :: cos_f,sin_f

 ecc_vec = get_eccentricity_vector(mu,dx,dv)
 h_hat = get_angmom_unit_vector(dx,dv)

 ! true anomaly: angle between eccentricity vector and position vector
 sin_f = dot_product(cross_product(ecc_vec,dx),h_hat)
 cos_f = dot_product(ecc_vec,dx)

 ! use atan2 for proper quadrant determination
 f = atan2(sin_f,cos_f)*rad_to_deg

end function get_true_anomaly

!------------------------------------------------------------
!+
!  Get longitude of ascending node
!+
!------------------------------------------------------------
real function get_longitude_of_ascending_node(mu,dx,dv) result(Omega)
 real, intent(in) :: mu,dx(3),dv(3)
 real :: n_vec(3), n_norm

 ! Longitude of ascending node: angle between line of nodes and x-axis
 ! -90.0 is because we define Omega as East of North
 n_vec = get_line_of_nodes_vector(dx,dv)
 n_norm = sqrt(dot_product(n_vec,n_vec))
 if (n_norm < tiny(1.0)) then
    Omega = 90.0
 else
    Omega = atan2(n_vec(2),n_vec(1))*rad_to_deg - 90.0
 endif

end function get_longitude_of_ascending_node

!------------------------------------------------------------
!+
!  Get argument of periapsis from position and velocity
!+
!------------------------------------------------------------
real function get_argument_of_periapsis(mu,dx,dv) result(w)
 real, intent(in) :: mu,dx(3),dv(3)
 real :: n_vec(3),ecc_vec(3),h_hat(3)
 real :: ecc_proj_x,ecc_proj_y,n_norm

 ecc_vec = get_eccentricity_vector(mu,dx,dv)
 h_hat = get_angmom_unit_vector(dx,dv)
 n_vec = get_line_of_nodes_vector(dx,dv)

 n_norm = sqrt(dot_product(n_vec,n_vec))
 if (n_norm < tiny(1.0)) then
    ! i ~ 0 or 180: define w from the sky-plane periapsis direction.
    w = atan2(ecc_vec(2),ecc_vec(1))*rad_to_deg
    ! handle singularity with Omega=90 (Node along x-axis)
    if (h_hat(3) < 0.0) then
       ! retrograde: w_sky = 180 - w_orbit
       w = 180.0 - w
    else
       ! prograde: w_sky = 180 + w_orbit
       w = w - 180.0
    endif
 else
    ! Argument of periapsis: angle between line of nodes and eccentricity vector
    ! Project eccentricity vector onto orbital plane and use atan2 for proper quadrant
    ecc_proj_x = dot_product(n_vec,ecc_vec)  ! component along line of nodes
    ecc_proj_y = dot_product(cross_product(h_hat,n_vec),ecc_vec)  ! component perpendicular to line of nodes
    w = atan2(ecc_proj_y,ecc_proj_x)*rad_to_deg
 endif

end function get_argument_of_periapsis

!----------------------------------------------------------------
!+
!  This routine outputs the orbital elements
!  semi-major axis (a), eccentricity (e), inclination angle (i),
!  argument of periastron (w), longitude of ascending node (O)
!  and true anomaly (f) from input position, velocity
!  and gravitational parameter mu = G(m1+m2)
!+
!----------------------------------------------------------------
subroutine get_orbital_elements(mu,dx,dv,a,e,inc,O,w,f)
 real, intent(in)  :: mu,dx(3),dv(3)
 real, intent(out) :: a,e,inc,O,w,f

 a = get_semimajor_axis(mu,dx,dv)
 e = get_eccentricity(mu,dx,dv)
 inc = get_inclination(dx,dv)
 O = get_longitude_of_ascending_node(mu,dx,dv)
 w = get_argument_of_periapsis(mu,dx,dv)
 f = get_true_anomaly(mu,dx,dv)

end subroutine get_orbital_elements

!----------------------------------------------------------
!+
! routine to extract main orbital parameters needed in
! the subgroup module
!+
!----------------------------------------------------------
subroutine get_orbparams(dr,dv,mu,r,v2,aij,eij,apoij,Tij)
 real, intent(in)    :: dr(3),dv(3),r,v2,mu
 real, intent(out)   :: aij,eij,apoij,Tij

 aij = get_semimajor_axis(mu,r,v2)

 eij = get_eccentricity(mu,r,dr(1),dr(2),dr(3),dv(1),dv(2),dv(3))

 Tij = get_orbital_period(mu,aij)

 apoij = aij*(1+eij)

end subroutine get_orbparams

!----------------------------------------------------------------
!+
!  Compute the Keplerian elements from the position and velocity:
!  returns mean anomaly in degrees instead of true anomaly
!+
!----------------------------------------------------------------
subroutine get_elements_mean_anomaly(mu,r,x,y,z,vx,vy,vz,a,ecc,inc,Omega,w,M)
 real, intent(in)  :: mu,r,x,y,z,vx,vy,vz
 real, intent(out) :: a,ecc,inc,Omega,w,M
 real :: E,nu,f
 real :: dx(3),dv(3)

 dx = [x,y,z]
 dv = [vx,vy,vz]
 call get_orbital_elements(mu,dx,dv,a,ecc,inc,Omega,w,f)

 ! convert true anomaly to mean anomaly
 nu = f*deg_to_rad ! radians
 E = get_E_from_true_anomaly(nu,ecc)
 M = get_mean_anomaly_from_E(E,ecc)*rad_to_deg

end subroutine get_elements_mean_anomaly

!-------------------------------------------------------------
!+
!  Convert rp,d,e to semi-major axis and initial true anomaly
!
!  For hyperbolic orbits (e > 1) we use a = rp/(1-e)
!
!  For parabolic orbits (e = 1), we set a = rp since this
!  is how the set_binary routine handles parabolic orbits
!+
!-------------------------------------------------------------
subroutine convert_flyby_to_elements(rp,d,e,a,f)
 real, intent(in)  :: rp,d,e
 real, intent(out) :: a,f

 if (orbit_is_parabolic(e)) then
    a = rp  ! for parabolic we return a = rp
 else
    a = rp / (1.-e)  ! semimajor axis is +ve or -ve for e < 1 and e > 1, respectively
 endif
 f = -abs(get_true_anomaly_from_separation(a,e,d))

end subroutine convert_flyby_to_elements

!-----------------------------------------------------------------
!+
!  Convert observed separation and velocity to rp,e,d,O,w,i
!
!  Inputs:
!   mu: G(m1+m2), the gravitational parameter
!   dx(3): observed separation vector
!   dv(3): observed velocity difference vector
!   initial_sep: initial separation
!
!  Outputs:
!   rp: pericentre distance
!   e: eccentricity
!   d: current separation
!   O: position angle of the ascending node
!   w: argument of periapsis
!   i: inclination
!   time_to_obs: time to reach observed separation from initial separation
!+
!-----------------------------------------------------------------
subroutine convert_posvel_to_flyby(mu,dx,dv,rp,d,e,O,w,inc,initial_sep,time_to_obs)
 real, intent(in)  :: mu,dx(3),dv(3),initial_sep
 real, intent(out) :: rp,d,e,O,w,inc
 real, intent(out) :: time_to_obs
 real :: a,f
 real :: current_sep

 ! Calculate current separation
 current_sep = sqrt(dot_product(dx,dx))
 d = current_sep

 ! Calculate orbital elements from position and velocity
 call get_orbital_elements(mu,dx,dv,a,e,inc,O,w,f)

 ! Calculate pericenter distance
 rp = get_pericentre_distance(mu,dx,dv)

 ! Calculate time to reach observed separation from initial separation (could be negative)
 time_to_obs = get_time_to_separation(mu,a,e,current_sep,initial_sep)

 if (time_to_obs < 0.0) then
    print*, 'Warning: time to reach observed separation is negative'
 endif

end subroutine convert_posvel_to_flyby

!-------------------------------------------------------------
!+
!  Calculate time to reach a specific separation for
!  elliptical, hyperbolic and parabolic orbits
!
!  Inputs:
!   a: semi-major axis
!   e: eccentricity
!   r1: initial separation
!   r2: final separation
!   mu: gravitational parameter
!
!  Outputs:
!   time: time to reach the final separation, can be negative
!+
!-------------------------------------------------------------
real function get_time_to_separation(mu,a,e,r1,r2) result(time)
 real, intent(in) :: mu,a,e,r1,r2
 real :: f1,f2

 ! Calculate true anomalies for both separations
 f1 = get_true_anomaly_from_separation(a,e,r1)
 f2 = get_true_anomaly_from_separation(a,e,r2)

 time = get_time_between_true_anomalies(mu,a,e,f1,f2)

end function get_time_to_separation

!-------------------------------------------------------------
!+
!  Calculate time between two true anomalies
!+
!-------------------------------------------------------------
real function get_time_between_true_anomalies(mu,a,e,f1,f2) result(time)
 real, intent(in) :: mu,a,e,f1,f2
 real :: nu1,nu2,E1,E2,M1,M2

 ! convert to eccentric anomalies
 nu1 = f1*deg_to_rad
 nu2 = f2*deg_to_rad
 E1 = get_E_from_true_anomaly(nu1,e)
 E2 = get_E_from_true_anomaly(nu2,e)
 ! convert to mean anomalies
 M1 = get_mean_anomaly_from_E(E1,e)
 M2 = get_mean_anomaly_from_E(E2,e)

 if (orbit_is_parabolic(e)) then
    ! parabolic case
    time = (M2 - M1) * sqrt(2.0 * abs(a)**3 / mu)
 else
    ! time from mean anomalies for hyperbolic and elliptic orbits
    time = (M2 - M1) * sqrt(abs(a)**3 / mu)  ! = (M2 - M1) * period / 2*pi for elliptic
 endif

end function get_time_between_true_anomalies

!----------------------------------------------------------------
!+
!  Calculate relative position and velocity vectors between
!  two point masses
!+
!----------------------------------------------------------------
subroutine get_dx_dv_ptmass(xyzmh_ptmass,vxyz_ptmass,dx,dv,i,j)
 real, intent(in)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(out) :: dx(3),dv(3)
 integer, intent(in), optional :: i,j
 integer :: i1,i2

 i1 = 1; i2 = 2
 if (present(i)) i1 = i
 if (present(j)) i2 = j
 dx = xyzmh_ptmass(1:3,i2) - xyzmh_ptmass(1:3,i1)
 dv = vxyz_ptmass(1:3,i2)  - vxyz_ptmass(1:3,i1)

end subroutine get_dx_dv_ptmass

!----------------------------------------------------------------
!+
!  This subroutine calculates the Kerr metric's Innermost stable
!  circular orbit. These formulae were obtained from
!  Bardeen & Teukolsky 1972 Astrophysical Journal, Vol. 178, pp. 347-370
!+
!----------------------------------------------------------------
subroutine isco_kerr(a,mass_bh,r_isco)
 real, intent(in) :: a,mass_bh
 real, intent(out) :: r_isco
 real  :: z1,z2

 ! Modified the formulas for Z1 and Z2 so that it uses the spin parameter
 ! used in the main code, which is a/M.

 z1 = 1 + (1-a**2.)**(1./3.)*((1+a)**(1./3.) + (1-a)**(1./3.))
 z2 = ((3*a**2.) + z1**2.)**(1./2.)

 !now we check if the value of a is +ve or -ve
 !+ve implies prograde while -ve implies retrograde

 if (a>=0.) then
    print*, "prograde rotation of BH wrt orbit"
    r_isco = mass_bh*(3. + z2 - sqrt((3.-z1)*(3.+z1+2.*z2)))
 else
    r_isco = mass_bh*(3. + z2 + sqrt((3.-z1)*(3.+z1+2.*z2)))
    print*, "retrograde rotation of BH wrt orbit"
 endif
 print*,"----------------------------------------------------------"
 print*, "ISCO of KERR metric with spin ",a," is: ",r_isco
 print*,"----------------------------------------------------------"

end subroutine isco_kerr

!--------------------------------------------------------------------------
!+
!  Obtain velocity in Boyer–Lindquist coordinates using gradient descent method. Iterations are
!  performed until the target specific energy epsilon is reached
!  References: Tejeda et al. (2017), MNRAS 469, 4483–4503
!  Owner: Mario Aguilar Faúndez
!+
!--------------------------------------------------------------------------
subroutine refine_velocity(x,y,z,vx,vy,vz,M_h,a,r,epsilon_target,alpha,delta_v,tol,max_iters)
 real,    intent(in)    :: x,y,z,M_h,a,r,epsilon_target,alpha,delta_v,tol
 integer, intent(in)    :: max_iters
 real,    intent(inout) :: vx,vy,vz
 real :: epsilon_0, epsilon_x, epsilon_y
 real :: d_eps_dx, d_eps_dy
 real :: temp_vx, temp_vy
 real :: sign_epsilon
 integer :: iter

 print*, 'Initial velocities: vx = ', vx, ', vy = ', vy, ', vz = ', vz, 'velocity magnitude = ', sqrt(vx**2 + vy**2 + vz**2)
 iter = 0
 do while (iter < max_iters)
    ! Compute epsilon at current velocity
    call get_specific_energy_gr(x,y,z,vx,vy,vz,M_h,a,r,epsilon_0)

    ! Check convergence
    if (abs(epsilon_0 - epsilon_target) < tol) exit

    ! Determine sign of epsilon - epsilon_target
    sign_epsilon = sign(1.0, epsilon_0 - epsilon_target)

    ! We only iterate in the x-y plane.
    ! Compute epsilon_x by changing vx by a small delta
    temp_vx = vx + delta_v
    call get_specific_energy_gr(x,y,z,temp_vx,vy,vz,M_h,a,r,epsilon_x)
    d_eps_dx = (epsilon_x - epsilon_0) / delta_v

    ! Compute epsilon_y by changing vy by a small delta
    temp_vy = vy + delta_v
    call get_specific_energy_gr(x,y,z,vx,temp_vy,vz,M_h,a,r,epsilon_y)
    d_eps_dy = (epsilon_y - epsilon_0) / delta_v

    ! Update velocities using gradient descent with sign adjustment
    vx = vx - alpha * sign_epsilon * d_eps_dx
    vy = vy - alpha * sign_epsilon * d_eps_dy

    iter = iter + 1
 enddo
 print*, 'Total number of iterations for refining velocity: ', iter
 print*, 'Updated velocities: vx = ', vx, ', vy = ', vy, ', vz = ', vz, 'velocity magnitude = ', sqrt(vx**2 + vy**2 + vz**2)
 call get_specific_energy_gr(x,y,z,vx,vy,vz,M_h,a,r,epsilon_0)
 if (iter == max_iters) then
    print *, 'Warning: Gradient descent did not converge after ', max_iters, ' iterations.'
 endif

end subroutine refine_velocity

!--------------------------------------------------------------------------
!+
!  Compute GR total specific energy
!  References: Tejeda et al. (2017), MNRAS 469, 4483–4503
!+
!--------------------------------------------------------------------------
subroutine get_specific_energy_gr(x,y,z,vx,vy,vz,M_h,a,r,epsilon_0)
 real, intent(in)  :: x, y, z, vx, vy, vz, M_h, a, r
 real, intent(out) :: epsilon_0
 real :: xy_term, r_dot_term, Gamma_flat, Gamma_curved, Gamma_0
 real :: rho_squared, Delta

 rho_squared = r**2 + (a**2 * z**2) / r**2
 Delta = r**2 - 2 * M_h * r + a**2
 xy_term = x * vy - y * vx
 r_dot_term = r**2 * (x * vx + y * vy) + (r**2 + a**2) * z * vz
 Gamma_flat = 1.0d0 - vx**2 - vy**2 - vz**2
 Gamma_curved = (2.0d0 * M_h * r / rho_squared) * ((1.0d0 - a * xy_term / (r**2 + a**2))**2 + &
                        r_dot_term**2 / (r**2 * Delta * (r**2 + a**2)))
 Gamma_0 = 1.0d0 / sqrt(Gamma_flat - Gamma_curved)

 epsilon_0 = Gamma_0 * (1.0d0 - (2.0d0 * M_h * r / rho_squared) * (1.0d0 - a * xy_term / (r**2 + a**2)))

end subroutine get_specific_energy_gr

!-------------------------------------------------------------
!
! cross product routine
!
!-------------------------------------------------------------
pure function cross_product(veca,vecb) result(vecc)
 real, intent(in)  :: veca(3),vecb(3)
 real :: vecc(3)

 vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
 vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
 vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)

end function cross_product

end module orbits
