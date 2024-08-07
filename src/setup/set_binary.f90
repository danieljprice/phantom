!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setbinary
!
! This module is contains utilities for setting up binaries
!   Our conventions for binary orbital parameters are consistent with
!   those produced by the imorbel code (Pearce, Wyatt & Kennedy 2015)
!   which can be used to produce orbits matching observed orbital
!   arcs of binary companions on the sky
!
! :References:
!   Eggleton (1983) ApJ 268, 368-369 (ref:eggleton83)
!   Lucy (2014), A&A 563, A126
!   Pearce, Wyatt & Kennedy (2015), MNRAS 448, 3679
!   https://en.wikipedia.org/wiki/Orbital_elements
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: binaryutils
!
 implicit none
 public :: set_binary,Rochelobe_estimate,L1_point,get_a_from_period
 public :: get_mean_angmom_vector,get_eccentricity_vector

 private
 interface get_eccentricity_vector
  module procedure get_eccentricity_vector,get_eccentricity_vector_sinks
 end interface get_eccentricity_vector

 real, parameter :: pi = 4.*atan(1.)
 real, parameter :: deg_to_rad = pi/180.
 integer, parameter :: &
   ierr_m1   = 1, &
   ierr_m2   = 2, &
   ierr_ecc  = 3, &
   ierr_semi = 4, &
   ierr_HIER1 = 5, &
   ierr_HIER2 = 6, &
   ierr_subststar = 7, &
   ierr_Omegasubst = 8, &
   ierr_missstar = 9

contains

!------------------------------------------------------------------------------
!+
!  setup for a binary orbit
!
!  INPUT:
!    m1 - mass of object 1
!    m2 - mass of object 2
!    semimajoraxis - semimajor axis (e/=1) or pericentre distance (e=1)
!    eccentricity - eccentricity
!    accretion_radius1 - accretion radius for point mass 1
!    accretion_radius2 - accretion radius for point mass 2
!    [optional] posang_ascnode - position angle of the ascending node (Omega, deg)
!    [optional] arg_peri - argument of periapsis (w, deg)
!    [optional] incl - orbital inclination (i, deg)
!    [optional] f - true anomaly (nu, deg)
!    [optional] mean_anomaly - mean anomaly (M, deg; replaces true anomaly)
!
!  OUTPUT: cartesian positions and velocities for both objects
!+
!------------------------------------------------------------------------------
subroutine set_binary(m1,m2,semimajoraxis,eccentricity, &
                      accretion_radius1,accretion_radius2, &
                      xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                      posang_ascnode,arg_peri,incl,f,mean_anomaly,verbose)
 use binaryutils, only:get_E,get_E_from_mean_anomaly,get_E_from_true_anomaly
 real,    intent(in)    :: m1,m2
 real,    intent(in)    :: semimajoraxis,eccentricity
 real,    intent(in)    :: accretion_radius1,accretion_radius2
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 integer, intent(out)   :: ierr
 real,    intent(in),  optional :: posang_ascnode,arg_peri,incl,f,mean_anomaly
 real,    intent(out), optional :: omega_corotate
 logical, intent(in),  optional :: verbose
 integer :: i1,i2,i
 real    :: mtot,dx(3),dv(3),Rochelobe1,Rochelobe2,period,bigM,rperi,rapo
 real    :: x1(3),x2(3),v1(3),v2(3),omega0,cosi,sini,xangle,reducedmass,angmbin
 real    :: a,E,E_dot,P(3),Q(3),omega,big_omega,inc,ecc,tperi
 real    :: term1,term2,term3,term4,theta,theta_max,energy
 logical :: do_verbose
 character(len=12) :: orbit_type

 ierr = 0
 do_verbose = .true.
 if (present(verbose)) do_verbose = verbose

 i1 = nptmass + 1
 i2 = nptmass + 2
 nptmass = nptmass + 2

 ! masses
 mtot = m1 + m2
 reducedmass = m1*m2/mtot

 ! check for stupid parameter choices
 if (m1 <= 0.) then
    print "(1x,a)",'ERROR: set_binary: primary mass <= 0'
    ierr = ierr_m1
 endif
 if (m2 < 0.) then
    print "(1x,a)",'ERROR: set_binary: secondary mass < 0'
    ierr = ierr_m2
 endif
 if (abs(semimajoraxis) <= tiny(0.)) then
    print "(1x,a)",'ERROR: set_binary: semi-major axis = 0'
    ierr = ierr_semi
 endif
 if (semimajoraxis < 0. .and. eccentricity <= 1.) then
    print "(1x,a)",'ERROR: set_binary: using a < 0 requires e > 1'
    ierr = ierr_semi
 endif
 if (eccentricity < 0.) then
    print "(1x,a)",'ERROR: set_binary: eccentricity must be positive'
    ierr = ierr_ecc
 endif
 if (eccentricity > 1. .and. present(f)) then
    theta = f*pi/180.
    theta_max = acos(-1./eccentricity)
    if (abs(theta) > theta_max) then
       print "(1x,2(a,f8.2))",'ERROR: max true anomaly for e = ',eccentricity, &
                              ' is |nu| < ',theta_max*180./pi
       ierr = ierr_ecc
    endif
 endif
 ! exit routine if cannot continue
 if (ierr /= 0) return

 ! set parameters that depend on the orbit type
 if (eccentricity < 1.) then
    a = abs(semimajoraxis)
    rperi = a*(1. - eccentricity)
    rapo  = semimajoraxis*(1. + eccentricity)
    period = sqrt(4.*pi**2*a**3/mtot)
    angmbin = reducedmass*sqrt(mtot*a*(1. - eccentricity**2))
    energy = -mtot/(2.*a)
 elseif (eccentricity > 1.) then
    a = -abs(semimajoraxis)
    rperi = a*(1. - eccentricity)
    rapo  = huge(rapo)
    period = huge(period)
    angmbin = reducedmass*sqrt(mtot*a*(1. - eccentricity**2))
    energy = -mtot/(2.*a)
 else
    a = huge(a)
    rperi = abs(semimajoraxis) ! for parabolic orbit we must give the pericentre distance
    rapo  = huge(rapo)
    period = huge(period)
    angmbin = reducedmass*sqrt(2.*mtot*rperi)
    energy = 0.
 endif

 Rochelobe1 = Rochelobe_estimate(m2,m1,rperi)
 Rochelobe2 = Rochelobe_estimate(m1,m2,rperi)

 if (do_verbose) then
    print "(/,2x,a)",'---------- binary parameters ----------- '
    print "(8(2x,a,1pg14.6,/),2x,a,1pg14.6)", &
        'primary mass     :',m1, &
        'secondary mass   :',m2, &
        'mass ratio m2/m1 :',m2/m1, &
        'reduced mass     :',reducedmass, &
        'semi-major axis  :',a, &
        'period           :',period, &
        'eccentricity     :',eccentricity, &
        'pericentre       :',rperi, &
        'apocentre        :',rapo
 endif
 if (accretion_radius1 > Rochelobe1) then
    print "(1x,a)",'WARNING: set_binary: accretion radius of primary > Roche lobe at periastron'
 endif
 if (accretion_radius2 > Rochelobe2) then
    print "(1x,a)",'WARNING: set_binary: accretion radius of secondary > Roche lobe at periastron'
 endif

 dx = 0.
 dv = 0.
 if (present(posang_ascnode) .and. present(arg_peri) .and. present(incl)) then
    ! Campbell elements
    ecc = eccentricity
    omega     = arg_peri*deg_to_rad
    ! our conventions here are Omega is measured East of North
    big_omega = posang_ascnode*deg_to_rad + 0.5*pi
    inc       = incl*deg_to_rad

    if (present(f)) then
       ! get eccentric, parabolic or hyperbolic anomaly from true anomaly
       ! (https://en.wikipedia.org/wiki/Eccentric_anomaly#From_the_true_anomaly)
       theta = f*deg_to_rad
       E = get_E_from_true_anomaly(theta,ecc)
    elseif (present(mean_anomaly)) then
       ! get eccentric anomaly from mean anomaly by solving Kepler equation
       bigM = mean_anomaly*deg_to_rad
       E = get_E_from_mean_anomaly(bigM,ecc)
    else
       ! set binary at apastron
       tperi = 0.5*period ! time since periastron: half period = apastron

       ! Solve Kepler equation for eccentric anomaly
       call get_E(period,eccentricity,tperi,E)
    endif

    ! Positions in plane (Thiele-Innes elements)
    P(1) = cos(omega)*cos(big_omega) - sin(omega)*cos(inc)*sin(big_omega)
    P(2) = cos(omega)*sin(big_omega) + sin(omega)*cos(inc)*cos(big_omega)
    P(3) = sin(omega)*sin(inc)
    Q(1) = -sin(omega)*cos(big_omega) - cos(omega)*cos(inc)*sin(big_omega)
    Q(2) = -sin(omega)*sin(big_omega) + cos(omega)*cos(inc)*cos(big_omega)
    Q(3) = sin(inc)*cos(omega)

    if (eccentricity < 1.) then ! eccentric
       orbit_type = 'Eccentric'
       term1 = a*(cos(E)-ecc)
       term2 = a*(sqrt(1. - ecc*ecc)*sin(E))
       E_dot = sqrt((m1 + m2)/(a**3))/(1.-ecc*cos(E))
       term3 = a*(-sin(E)*E_dot)
       term4 = a*(sqrt(1.- ecc*ecc)*cos(E)*E_dot)
    elseif (eccentricity > 1.) then ! hyperbolic
       orbit_type = 'Hyperbolic'
       term1 = a*(cosh(E)-ecc)
       term2 = -a*(sqrt(ecc*ecc - 1.)*sinh(E))
       E_dot = sqrt((m1 + m2)/(abs(a)**3))/(ecc*cosh(E)-1.)
       term3 = a*(sinh(E)*E_dot)
       term4 = -a*(sqrt(ecc*ecc - 1.)*cosh(E)*E_dot)
    else ! parabolic
       orbit_type = 'Parabolic'
       term1 = rperi*(1. - E*E)
       term2 = rperi*(2.*E)
       E_dot = sqrt(2.*(m1 + m2)/(rperi**3))/(1. + E*E)
       term3 = -E*(rperi*E_dot)
       term4 = rperi*E_dot
    endif

    if (do_verbose) then
       print "(4(2x,a,1pg14.6,/),2x,a,1pg14.6)", &
             trim(orbit_type)//' anomaly:',E, &
             'E_dot            :',E_dot, &
             'inclination     (i, deg):',incl, &
             'angle asc. node (O, deg):',posang_ascnode, &
             'arg. periapsis  (w, deg):',arg_peri
       if (present(f)) print "(2x,a,1pg14.6)", &
             'true anomaly    (f, deg):',f
       if (present(mean_anomaly)) print "(2x,a,1pg14.6)", &
             'mean anomaly    (M, deg):',mean_anomaly
    endif

    ! Rotating everything
    ! Set the positions for the primary and the central secondary
    dx(:) = term1*P(:) + term2*Q(:)

    ! Set the velocities
    dv(:) = term3*P(:) + term4*Q(:)

 else
    ! set binary at apastron
    dx = (/semimajoraxis*(1. + eccentricity),0.,0./)
    dv = (/0.,sqrt(semimajoraxis*(1.-eccentricity**2)*mtot)/dx(1),0./)
 endif

 ! positions of each star so centre of mass is at zero
 x1 = -dx*m2/mtot
 x2 =  dx*m1/mtot

 ! velocities
 v1 = -dv*m2/mtot
 v2 =  dv*m1/mtot

 omega0 = v2(2)/x2(1)

 ! print info about positions and velocities
 if (do_verbose) then
    print "(9(2x,a,1pg14.6,/),2x,a,1pg14.6)", &
        'energy (mtot/2a) :',energy,&
        'energy (KE+PE)   :',-mtot/sqrt(dot_product(dx,dx)) + 0.5*dot_product(dv,dv),&
        'angular momentum :',angmbin, &
        'mean ang. speed  :',omega0, &
        'Omega_0 (prim)   :',v2(2)/x2(1), &
        'Omega_0 (second) :',v2(2)/x2(1), &
        'R_accretion (1)  :',accretion_radius1, &
        'R_accretion (2)  :',accretion_radius2, &
        'Roche lobe  (1)  :',Rochelobe1, &
        'Roche lobe  (2)  :',Rochelobe2
 endif

 if (present(omega_corotate)) then
    if (do_verbose) print "(a)",' SETTING VELOCITIES FOR COROTATING FRAME: '
    omega_corotate = omega0
    v1(2) = v1(2) - omega0*x1(1)
    v2(2) = v2(2) - omega0*x2(1)
    if (do_verbose) print "(2(2x,a,1pg14.6,/))", &
     'Omega_0 (primary)     :',v1(2)/x1(1), &
     'Omega_0 (secondary)   :',v2(2)/x2(1)
 endif

 ! conclude printout
 if (do_verbose) print "(2x,40('-'),/)"

!
!--positions and accretion radii
!
 xyzmh_ptmass(:,i1:i2) = 0.
 xyzmh_ptmass(1:3,i1) = x1
 xyzmh_ptmass(1:3,i2) = x2
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(5,i1) = accretion_radius1
 xyzmh_ptmass(5,i2) = accretion_radius2
 xyzmh_ptmass(6,i1) = 0.0
 xyzmh_ptmass(6,i2) = 0.0
!
!--velocities
!
 vxyz_ptmass(:,i1) = v1
 vxyz_ptmass(:,i2) = v2
!
! rotate if inclination is non-zero
!
 if (present(incl) .and. .not.(present(arg_peri) .and. present(posang_ascnode))) then
    xangle = incl*deg_to_rad
    cosi = cos(xangle)
    sini = sin(xangle)
    do i=i1,i2
       call rotate(xyzmh_ptmass(1:3,i),cosi,sini)
       call rotate(vxyz_ptmass(1:3,i),cosi,sini)
    enddo
 endif

end subroutine set_binary

pure subroutine rotate(xyz,cosi,sini)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: cosi,sini
 real :: xi,yi,zi

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) =  xi*cosi + zi*sini
 xyz(2) =  yi
 xyz(3) = -xi*sini + zi*cosi

end subroutine rotate

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

!-------------------------------------------------------------
! Function to determine the semi-major axis given the period
!-------------------------------------------------------------
function get_a_from_period(m1,m2,period) result(a)
 real, intent(in) :: m1,m2,period
 real :: a

 a = ((m1 + m2)*(period/(2.*pi))**2)**(1./3.)

end function get_a_from_period

!----------------------------------------------------
! Eccentricity vector, for second body
! https://en.wikipedia.org/wiki/Eccentricity_vector
!----------------------------------------------------
function get_eccentricity_vector(m1,m2,x1,x2,v1,v2)
 real, intent(in) :: m1,m2
 real, intent(in) :: x1(3),x2(3),v1(3),v2(3)
 real :: x0(3),v0(3),r(3),v(3),dr,mu
 real :: get_eccentricity_vector(3)

 ! centre of mass position and velocity
 x0 = (m1*x1 + m2*x2)/(m1 + m2)
 v0 = (m1*v1 + m2*v2)/(m1 + m2)

 ! position and velocity vectors relative to each other
 r = x2 - x1
 v = v2 - v1

 ! intermediate quantities
 dr = 1./sqrt(dot_product(r,r))
 mu = m1 + m2  ! "standard gravitational parameter"

 ! formula for eccentricity vector
 get_eccentricity_vector = (dot_product(v,v)/mu - dr)*r - dot_product(r,v)/mu*v

end function get_eccentricity_vector

!----------------------------------------------------
! interface to above assuming two sink particles
!----------------------------------------------------
function get_eccentricity_vector_sinks(xyzmh_ptmass,vxyz_ptmass,i1,i2)
 real,    intent(in) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: i1, i2
 real :: get_eccentricity_vector_sinks(3)

 if (i1 > 0 .and. i2 > 0) then
    get_eccentricity_vector_sinks = get_eccentricity_vector(&
        xyzmh_ptmass(4,i1),xyzmh_ptmass(4,i2),&
        xyzmh_ptmass(1:3,i1),xyzmh_ptmass(1:3,i2),&
        vxyz_ptmass(1:3,i1),vxyz_ptmass(1:3,i2))
 else
    get_eccentricity_vector_sinks = 0.
 endif

end function get_eccentricity_vector_sinks

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
    call get_cross_product(xyz(:,i),vxyz(:,i),li)
    l = l + li
 enddo
 l = l/real(n)

end function get_mean_angmom_vector

!-------------------------------------------------------------
!
! cross product routine
!
!-------------------------------------------------------------
pure subroutine get_cross_product(veca,vecb,vecc)
 real, intent(in)  :: veca(3),vecb(3)
 real, intent(out) :: vecc(3)

 vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
 vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
 vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)

end subroutine get_cross_product

end module setbinary
