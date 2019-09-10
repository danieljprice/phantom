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

subroutine get_rad_accel_from_sinks(nptmass,npart,xyzh,xyzmh_ptmass,fext)
  use part,    only:isdead_or_accreted
#ifdef NUCLEATION
 use part,  only:nucleation
#endif
#ifdef SINKRADIATION
  use part,  only:dust_temp
#endif
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:)
 real,     intent(in)    :: xyzmh_ptmass(:,:)
 real,     intent(inout) :: fext(:,:)
 real                    :: dx,dy,dz,xa,ya,za,r,pmassj,plumj,pmlossj,ax,ay,az
 integer                 :: i,j

 do j=1,nptmass
    pmassj  = xyzmh_ptmass(4,j)
    plumj   = xyzmh_ptmass(12,j)
    pmlossj = xyzmh_ptmass(14,j)
    xa = xyzmh_ptmass(1,j)
    ya = xyzmh_ptmass(2,j)
    za = xyzmh_ptmass(3,j)
    !$omp parallel  do default(none) &
#ifdef NUCLEATION
    !$omp shared(nucleation)&
#endif
#ifdef SINKRADIATION
    !$omp shared(dust_temp)&
#endif
    !$omp shared(npart,xa,ya,za,pmassj,plumj,pmlossj,xyzh,fext) &
    !$omp private(i,dx,dy,dz,ax,ay,az,r)
    do i=1,npart
      if (.not.isdead_or_accreted(xyzh(4,i))) then
        dx = xyzh(1,i) - xa
        dy = xyzh(2,i) - ya
        dz = xyzh(3,i) - za
        r = sqrt(dx**2 + dy**2 + dz**2)
#ifdef NUCLEATION
        call get_radiative_acceleration_from_star(r,dx,dy,dz,pmassj,plumj,pmlossj,ax,ay,az,K3=nucleation(5:i))
#elif BOWEN
        call get_radiative_acceleration_from_star(r,dx,dy,dz,pmassj,plumj,pmlossj,ax,ay,az,Tdust=dust_temp(i))
#else
        call get_radiative_acceleration_from_star(r,dx,dy,dz,pmassj,plumj,pmlossj,ax,ay,az)
#endif
        fext(1,i) = fext(1,i) + ax
        fext(2,i) = fext(2,i) + ay
        fext(3,i) = fext(3,i) + az
      endif
    end do
    !$omp end parallel do
 enddo

 end subroutine get_rad_accel_from_sinks

 subroutine get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar,Lstar,Mdot,ax,ay,az,K3,Tdust)
#ifdef NUCLEATION
  use dust_formation, only:calc_alpha_dust
#elif BOWEN
  use dust_formation, only:calc_alpha_bowen
#endif
  real, intent(in)    :: r,dx,dy,dz,Mstar,Lstar,Mdot
  real, optional, intent(in)  :: K3,Tdust
  real, intent(out)   :: ax,ay,az
  real :: fac,alpha_dust

#ifdef NUCLEATION
  call calc_alpha_dust(Mstar,Lstar,Mdot,K3,alpha_dust)
  fac = alpha_dust*Mstar/r**3
#elif BOWEN
  if (alpha_rad > 0. ) then
     alpha_dust = alpha_rad
  else
     call calc_alpha_bowen(Mstar,Lstar,Tdust,alpha_dust)
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
    if (alpha_rad < 0.) call fatal(label,'invalid setting for alpha_rad (must be >= 0)')
    if (alpha_rad > 0.) irad_accel = 1
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)
#ifdef BOWEN
 irad_accel = 1
#endif
end subroutine read_options_radiative_accel

end module radiative_accel
