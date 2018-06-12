!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: externalforces, io, part, physcon, prompting,
!    setup_params, spherical, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 real, parameter :: rho_crit = 10.0     ! g/cm^3
 real, parameter :: temp_crit = 1.0e7   ! K

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,maxvxyzu
 use spherical, only:set_sphere
 use io,        only:master
 use setup_params, only:rhozero,npart_total
 use externalforces, only:mass1
 use units,          only:set_units,udist,umass,utime
 use physcon,        only:pi,solarr,solarm,Rg,gg
 use prompting,      only: prompt
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: massr,m1,m2,a,hacc1,hacc2,ecc
 real :: P1,P2,P3,Q1,Q2,Q3,big_omega,omega,inc,E,v1,v2,E_dot
 real :: rmin,rmax,psep,totvol,dr,t,dt,total_mass,r
 integer :: i,nx,np,npmax
 real :: Rg_codeunits

 call set_units(dist=solarr, mass=solarm, G=1.0)

!
!--parameters for this setup
!
 ! cgs units
 a = sqrt((Rg * temp_crit) / (pi * gg * rho_crit))
 total_mass = sqrt(16.0 * pi / rho_crit) * (Rg * temp_crit / gg)**(1.5)

 ! convert to code units
 a = a / udist
 total_mass = total_mass / umass

!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 7./3.
 rmin = 0.0
 rmax = pi*a

!
!--setup
!
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.

 np    = min(10000000,int(2.0/3.0*size(xyzh(1,:)))) ! approx max number allowed in sphere given size(xyzh(1,:))
 npmax = np
 call prompt('Enter the approximate number of particles in the sphere ',np,0,npmax)

 totvol = 4./3.*pi*rmax**3
 nx = int((2.2/3.0*np)**(1./3.))
 psep = totvol**(1./3.)/real(nx)
 print *, ' rho_crit = ', rho_crit, ' in g/cm^3'
 print *, ' radius = ', rmax, ' in solar radius'
 print *, ' total volume = ',totvol,' particle separation = ',psep
 print *, ' totalvol / partsep**3 = ', totvol/psep**3
 print *, ''

 npart = 0
 npart_total = 0

 call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,xyzh,rhofunc=rhofunc,nptot=npart_total)

!
!--set particle properties
!

 npart_total = npart
 rhozero = total_mass / totvol
 print *, ' total mass = ',total_mass,' mean density = ',rhozero

 npartoftype(:) = 0
 npartoftype(1) = npart
 print *, ' npart = ',npart_total

 massoftype(1) = total_mass / npart_total
 print *, ' particle mass = ',massoftype(1)
 print *, ''

 ! convert Rg to code units
 Rg_codeunits = Rg * utime * utime / (udist * udist)

 do i=1,npart
    vxyzu(1:3,i) = 0.
    r = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))
    if (maxvxyzu >= 4) then
       vxyzu(4,i) = 1.5 * Rg_codeunits * temp_crit * a * sin(r/a) / r
    endif
 enddo


end subroutine setpart

real function rhofunc(r)
 use units,          only:udist,umass
 use physcon,        only:pi,Rg,gg,solarm,solarr
 real, intent(in) :: r
 real :: a

 ! cgs units
 a = sqrt((Rg * temp_crit) / (pi * gg * rho_crit))

 ! convert to code units
 a = a / udist

 ! cgs units
 rhofunc = rho_crit * a * sin(r/a) / r

 ! convert to code units
 rhofunc = rhofunc / (umass / udist**3)

end function rhofunc

end module setup
