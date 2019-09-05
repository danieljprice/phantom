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
!  DEPENDENCIES: dim, eos, kernel, part, physcon, units
!+
!--------------------------------------------------------------------------
module radiative_accel
 implicit none
 integer, public  :: irad_accel = 0
 real,    private :: alpha_rad = 0.0

 public :: get_rad_accel_from_sinks,read_options_radiative_accel,write_options_radiative_accel
contains

subroutine init_radiative_accel(ierr)
  integer, intent(out) :: ierr

  ierr = 0

end subroutine init_radiative_accel

subroutine get_rad_accel_from_sinks(nptmass,xi,yi,zi,xyzmh_ptmass,fextrad,ii)
 use dim,  only:maxptmass
 use part, only:nsinkproperties
 integer,           intent(in)    :: nptmass,ii
 real,              intent(in)    :: xi,yi,zi
 real,              intent(in)    :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real,    optional, intent(out)   :: fextrad(3)
 real                             :: dx,dy,dz,r,pmassj,plumj,pmlossj,ax,ay,az
 integer                          :: j

 fextrad = 0.
 do j=1,nptmass
   dx = xi - xyzmh_ptmass(1,j)
   dy = yi - xyzmh_ptmass(2,j)
   dz = zi - xyzmh_ptmass(3,j)
   r = sqrt(dx**2 + dy**2 + dz**2)
   pmassj  = xyzmh_ptmass(4,j)
   plumj   = xyzmh_ptmass(12,j)
   pmlossj = xyzmh_ptmass(13,j)
   call get_radiative_acceleration_from_star(r,dx,dy,dz,pmassj,plumj,pmlossj,ax,ay,az,ii)
   fextrad(1) = fextrad(1) + ax
   fextrad(2) = fextrad(2) + ay
   fextrad(3) = fextrad(3) + az
 end do

 end subroutine get_rad_accel_from_sinks

 subroutine get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar,Lstar,Mdot,ax,ay,az,ii)
#ifdef NUCLEATION
  use part,           only:nucleation
  use dust_formation, only:calc_alpha_dust
#elif BOWEN
  use dust_formation, only:calc_alpha_bowen
#endif
  integer, intent(in) :: ii
  real, intent(in)  :: r,dx,dy,dz,Mstar,Lstar,Mdot
  real, intent(out) :: ax,ay,az
  real :: fac,alpha_dust

#ifdef NUCLEATION
  call calc_alpha_dust(Mstar,Lstar,Mdot,nucleation(5,ii),alpha_dust)
  fac = alpha_dust*Mstar/r**3
#elif BOWEN
  if (alpha_rad > 0 ) then
     alpha_dust = alpha_rad
  else
     call calc_alpha_bowen(Mstar,Lstar,alpha_dust)
  endif
  fac = alpha_dust*Mstar/r**3
#else
  ! alpha wind
  fac = alpha_rad*Mstar/r**3
#endif
  ax = fac*dx
  ay = fac*dy
  az = fac*dz
end subroutine get_radiative_acceleration_from_star

subroutine write_options_radiative_accel(iunit)
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit
 call write_inopt(alpha_rad,'alpha_rad','fraction of the gravitational acceleration imparted to the gas',iunit)

 end subroutine write_options_radiative_accel

 subroutine read_options_radiative_accel(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_radiative_accel'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('alpha_rad')
    read(valstring,*,iostat=ierr) alpha_rad
    ngot = ngot + 1
    if (alpha_rad < 0.) call fatal(label,'invalid setting for alpha_rad (must be > 0)')
    if (alpha_rad > 0.) irad_accel = 1
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)
end subroutine read_options_radiative_accel

end module radiative_accel
