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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, io, options, part, physcon, setup_params,
!    spherical, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for bondi accretion
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master,fatal
 use spherical,    only:set_sphere
 use options,      only:ieos,iexternalforce
 use timestep,     only:tmax,dtmax
 use centreofmass, only:reset_centreofmass
 use units,        only:udist,umass,utime,set_units
 use physcon,      only:pc,solarm,gg,pi
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,igas
#ifdef GR
 use metric,      only:metric_type
#endif
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: psep,totmass,tff
 real    :: vol,rmax,rmin
 real    :: accretion_radius,gcode
 integer :: i,np,nx,ilattice,maxvxyzu,acc_rad
#ifdef GR
 real :: r,pos(3)
#endif

 call set_units(mass=solarm,G=1.d0,c=1.d0)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 5./3.

acc_rad = 2.1
#ifndef GR
 !--set sink properties
  nptmass = 1
  accretion_radius = acc_rad
  xyzmh_ptmass(:,nptmass) = 0.
  xyzmh_ptmass(4,nptmass) = 1.
  xyzmh_ptmass(ihacc,nptmass) = accretion_radius
  vxyz_ptmass(:,:) = 0.
#endif

 !--setup particles
 rmin  = 5.*acc_rad
 rmax  = 100.
 np = 1000
 vol = 4./3.*pi*rmax**3
 nx = int(np**(1./3.))
 psep = vol**(1./3.)/real(nx)
 print*,'Setup for gas: '
 print*,' volume = ',vol,' particle separation = ',psep
 print*,' maximum radius = ',rmax
 print*,' vol/psep**3 = ',vol/psep**3
 totmass = 10.*solarm/umass
 rhozero = totmass/vol
 tff = sqrt(3.*pi/(32.*rhozero))
 print*,' free fall time = ',tff

 tmax = tff*5
 dtmax = tmax/100.

 iexternalforce = 1
 ilattice = 2
 ieos = 11

 maxvxyzu = size(vxyzu(:,1))
 if (maxvxyzu < 4) then
    call fatal('setup','Setup requires thermal energy to be stored')
 endif
#ifdef GR
 if (.not. metric_type == 'Schwarzschild') then
    call fatal('setup','GR on, but metric is not Schwarzschild. Change the metric to Schwarzschild for gr_bondi setup.')
 endif
#endif

!--add stretched sphere
 npart = 0
 npart_total = 0
 call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,xyzh,rhofunc=rhofunc,nptot=npart_total)
 massoftype(igas) = totmass/npart_total

 gcode = gg*umass*utime**2/udist**3

 print *,'gcode = ',gcode
 print*,'npart = ',npart

 do i=1,npart
#ifdef GR
    pos = xyzh(1:3,i)
    r = sqrt(dot_product(pos,pos))
    vxyzu(1:3,i) = -vfunc(r)*pos/r
#endif

    vxyzu(4,i) = 0.
 enddo

!--reset centre of mass to the origin
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 npartoftype(:) = 0
 npartoftype(igas) = npart_total

end subroutine setpart


real function rhofunc(r)
#ifdef GR
use utils_gr, only:dot_product_gr,get_metric3plus1
real :: x(3),v(3),alpha,beta(3),gammaijdown(3,3),gammaijUP(3,3),gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
#endif
real, intent(in) :: r

#ifdef GR
 x = (/r,0.,0./)
 v = (/vfunc(r),0.,0./)
 call get_metric3plus1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 rhofunc = (sqrtg/alpha)*dfunc(r)
#else
 rhofunc = 1./r**2
#endif

end function rhofunc

real function dfunc(r)
 real, intent(in) :: r
 real :: d,m
 m = 1.
 d = 12.

 dfunc = d/(r**2*sqrt(2.*m/r*(1.- 2.*m/r)))

end function dfunc

#ifdef GR
real function gammafunc(r)
 use utils_gr, only:dot_product_gr,get_metric3plus1
 real, intent(in) :: r
 real :: m,x(3),v(3),v2,alpha,beta(3),gammaijdown(3,3),gammaijUP(3,3),gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
 m = 1.

 x = (/r,0.,0./)
 v = (/vfunc(r),0.,0./)

 call get_metric3plus1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)

 v2 = dot_product_gr(v,v,gammaijdown)

 gammafunc = 1./sqrt(1.-v2)

end function gammafunc

real function vfunc(r)
 real, intent(in) :: r
 real :: m
  m = 1.
  vfunc = sqrt(2.*m/r)*(1. - 2.*m/r)
end function vfunc
#endif

end module setup
