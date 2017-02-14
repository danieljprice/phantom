!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Rebecca Nealon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: inject, injectwind, io, part, physcon, setbinary,
!    setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,igas,maxp,maxphase,isetphase,iphase
 use setbinary,  only:set_binary
 use physcon,    only:au,solarm,years
 use units,      only:udist,umass,utime,set_units
 use inject,     only:inject_type,inj_wind
 use injectwind, only:wind_resolution,wind_star_radius,wind_mass_rate,wind_velocity, &
                      wind_star_mass,wind_star_radius,wind_osc_vamplitude, &
                      wind_sphdist, wind_injection_radius
 use injectwind, only:binary_mass_ratio,binary_semi_major_axis,binary_eccentricity,binary_companion_accretion_radius
 use setup_params, only:rhozero,npart_total
 use io,        only:master
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: m1,hacc1,massr,a,ecc,hacc2

 integer :: resolution, particles_per_sphere
 real :: phi, neighbour_distance, wind_velocity_u, wind_mass_rate_u
 real :: sphere_distance, wind_injection_radius_u

 real, parameter :: pi = 3.1415926536
 real :: rmin,rmax,psep,totvol,dr,t,dt,totmass,r
 integer :: i,nx,np,maxvxyzu
!
! units
!
 call set_units(mass=solarm,dist=au,G=1.d0)

 m1     = wind_star_mass * solarm / umass ! Depends on parameters in the input file
 hacc1  = wind_star_radius * au / udist ! Depends on parameters in the input file
 hacc2 = binary_companion_accretion_radius * au / udist ! Depends on parameters in the input file
 massr = binary_mass_ratio ! Depends on parameters in the input file
 a = binary_semi_major_axis * au / udist ! Depends on parameters in the input file
 ecc = binary_eccentricity ! Depends on parameters in the input file

 wind_velocity_u = max(wind_velocity, wind_osc_vamplitude) * (1.0d5/udist) / (1.0d0/utime) ! Input unit : kilometer per second
 wind_mass_rate_u = wind_mass_rate * (solarm/umass) / (years/utime) ! Input unit : solar mass per year
 wind_injection_radius_u = wind_injection_radius * (au/udist)
 resolution = int(wind_resolution)
 phi = (sqrt(5.)+1.)/2. ! Golden ratio
 neighbour_distance = 2./((2.*resolution-1.)*sqrt(sqrt(5.)*phi))
 particles_per_sphere = 20 * (2*resolution*(resolution-1)) + 12
 inject_type = inj_wind

!
!--general parameters
!
 time = 0.
 polyk = 0.
 gamma = 5./3.

!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 sphere_distance = wind_sphdist * neighbour_distance
 massoftype = sphere_distance*wind_injection_radius_u*wind_mass_rate_u &
     /(particles_per_sphere*wind_velocity_u) ! Depends on parameters in the input file inject_wind.f90

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0
 ! need to use E = 4.722: periapse at 2014.25
 call set_binary(m1,massr,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,&
                 posang_ascnode=1.4294,arg_peri=1.69646,incl=2.0612)

!*** INCLUDE ENVIRONMENT ***
 rmin = 0.1
 rmax = 8.10! 0.16

 np = 100 !size(xyzh(1,:))
 maxvxyzu = 4!size(vxyzu(:,1))
 totvol = 4./3.*pi*rmax**3
 nx = int(np**(1./3.))
 psep = totvol**(1./3.)/real(nx)
 print*,' Setup for medium around black hole '
 print*,' total volume = ',totvol,' particle separation = ',psep
 print*,totvol/psep**3

 npart = 0
 npart_total = 0

call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,xyzh,rhofunc=rhofunc,nptot=npart_total)

!
!--set particle properties
!
 totmass = 1.87d-3! 9.0d-6
 rhozero = totmass/totvol
 polyk   = 41249.51
 print*,' total mass = ',totmass,' mean density = ',rhozero,' polyk = ',polyk
 print*,' free fall time = ',sqrt(3.*pi/(32.*rhozero))

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 massoftype(1) = totmass/npart_total
 print*,' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1:3,i) = 0.
       !
       !--Note: Running the polytrope with u stored is not quite
       !  the same as using P = K rho^gamma because we really
       !  should use the actual rho, not rhozero.
       !
          r = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))
          vxyzu(4,i) = 4.3d6/(2*(gamma-1)*r) !polyk*rhozero**(gamma-1)/(gamma-1)
 enddo

end subroutine setpart

real function rhofunc(r)
 real, intent(in) :: r

 rhofunc = 1./r

end function rhofunc


end module setup
