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
 use dim,            only:gr
 use physcon,        only:pi
 use externalforces, only:accradius1
#ifdef GR
 use metric,         only:mass1
#else
 use externalforces, only:mass1
#endif
 implicit none
 public :: setpart

 private

!-- Choice of constants for non-relativistic solution
 real, parameter :: rc = 5.
 real, parameter :: rhocrit = 1.

!-- Choice of constants for GR solution
 real, parameter :: den0 = 1.    !  12.
 real, parameter :: en0  = 1.e-9 !  0.000297118

 real, save  :: mdot
 real, save  :: gamma_eos

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
 use options,      only:ieos,iexternalforce,nfulldump
 use timestep,     only:tmax,dtmax
 use centreofmass, only:reset_centreofmass
 use units,        only:udist,umass,utime,set_units
 use physcon,      only:pc,solarm,gg
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,igas,set_particle_type,iboundary
 use stretchmap,   only:get_mass_r
 use kernel,       only:radkern
 use externalforces,only:accradius1_hard
#ifdef GR
 use metric,       only:imetric
 use metric_tools, only:imet_schwarzschild
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
 real    :: psep,tff
 real    :: vol,rmax,rmin
 real    :: gcode
 integer :: i,np,nx,ilattice,maxvxyzu,nbound
 real    :: r,pos(3),cs2
 real    :: totmass
 real    :: approx_m,approx_h

!-- Set code units
 call set_units(mass=solarm,G=1.d0,c=1.d0)
 gcode = gg*umass*utime**2/udist**3
 print*,' gcode = ',gcode

 maxvxyzu = size(vxyzu(:,1))

!--Set general parameters
 time = 0.
 if (maxvxyzu >=4) then
    gamma      = 5./3.
    ieos = 2
 else
    gamma = 1.
    ieos = 1
 endif
 gamma_eos = gamma
 if(.not. gr ) then
    iexternalforce = 1
 endif

!--- Setup particles
 rmax = 18.1
 np   = 100*1000
 vol  = 4./3.*pi*rmax**3
 nx   = int(np**(1./3.))
 psep = vol**(1./3.)/real(nx)
 if(.not.gr) then
    cs2  = mass1/(2.*rc)
    polyk = cs2
    mdot = rhocrit*4.*pi*rc**2*sqrt(cs2)
 endif

 accradius1 = 2.
 if(gr) accradius1 = 3.
 approx_m    = get_mass_r(rhofunc,rmax,accradius1)/np
 approx_h    = hfact*(approx_m/rhofunc(accradius1))**(1./3.)
 accradius1_hard = (accradius1 - radkern*approx_h)

 rmin  = accradius1_hard
 if(gr.and.rmin <=2.) STOP 'rmin is less than Schwarzschild radius!'
 totmass = get_mass_r(rhofunc,rmax,rmin)
 rhozero = totmass/vol
 tff     = sqrt(3.*pi/(32.*rhozero))

 print*,''
 print*,' Setup for gas: '
 print*,' min,max radius = ',rmin,rmax
 print*,' volume         = ',vol        ,' particle separation = ',psep
 print*,' vol/psep**3    = ',vol/psep**3,' totmass             = ',totmass
 print*,' free fall time = ',tff        ,' tmax                = ',tmax
 print*,''

 tmax = totmass/mdot
 if(gr) tmax = 10.*tff
 dtmax = tmax/500.

 ilattice       = 2
 nfulldump      = 1
 iexternalforce = 1

#ifdef GR
 if (imetric /= imet_schwarzschild) then
    call fatal('setup','GR on, but metric is not Schwarzschild. Change the metric to Schwarzschild for gr_bondi setup.')
 endif
#endif

!--- Add stretched sphere
 npart = 0
 npart_total = 0
 call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,xyzh,rhofunc=rhofunc,nptot=npart_total)
 massoftype(:) = totmass/npart
 print*,' npart = ',npart
 print*,''

 nbound = 0
 do i=1,npart
    pos = xyzh(1:3,i)
    r = sqrt(dot_product(pos,pos))
    if(gr) then
       if (maxvxyzu >= 4) vxyzu(4,i)   = efunc(r)/dfunc(r)
       vxyzu(1:3,i) = vfunc(r)*pos/r
    else
       if (maxvxyzu >= 4) vxyzu(4,i)   = cs2/(gamma-1.)
       vxyzu(1:3,i) = vfunc(r,rc)*pos/r
    endif
    if (r + radkern*xyzh(4,i)>rmax) then
       call set_particle_type(i,iboundary)
       nbound = nbound + 1
    else
       call set_particle_type(i,igas)
    endif
 enddo

!--- Reset centre of mass to the origin
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 npartoftype(:) = 0
 npartoftype(igas) = npart_total-nbound
 npartoftype(iboundary) = nbound

end subroutine setpart
!=====----------------------------------------------------------------=====

!--- Density rhostar(r)
real function rhofunc(r)
#ifdef GR
use metric_tools, only:get_metric3plus1
use utils_gr,     only:dot_product_gr
real :: x(3),v(3),alpha,beta(3),gammaijdown(3,3),gammaijUP(3,3),gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
#endif
real, intent(in) :: r

#ifdef GR
 x = (/r,0.,0./)
 v = (/vfunc(r),0.,0./)
 call get_metric3plus1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 rhofunc = (sqrtg/alpha)*dfunc(r)
#else
 rhofunc = mdot/(4.*pi*abs(vfunc(r,rc))*r**2)
#endif
end function rhofunc

!--- Radial velocity Vr(r)
real function vfunc(r,r0)
 real, intent(in) :: r
 real, intent(in), optional :: r0
 real :: cs2,m
 m = mass1
 if(gr) vfunc = -sqrt(2.*m/r)*(1. - 2.*m/r)
 if(gr.and.present(r0)) vfunc = (1. - 2.*m/r)/sqrt(1.-2.*m/r0)*sqrt(2.*m*(1./r - 1./r0))
 if(.not.gr) then
     cs2 = m/(2.*r0)
     if (r>=r0) then
        vfunc = -sqrt(-cs2*lambertw_0(-func(r,r0)))
     else
        vfunc = -sqrt(-cs2*lambertw_neg1(-func(r,r0)))
     endif
  endif
end function vfunc


!=====----------------------------------------------------------------------------------------------=====
!
!---Functions for GR solution
!
!--------------------------------------------------------------------------------------------------------
! See Hawley, Smarr & Wilson 1984, for D(r) & E(r)
!
!--- Density D(r)
real function dfunc(r)
 real, intent(in) :: r
 dfunc = den0/(r**2*sqrt(2.*mass1/r*(1.- 2.*mass1/r)))
end function dfunc

!--- Energy E(r)
real function efunc(r)
 real, intent(in) :: r
 efunc = en0/((sqrt(2.*mass1/r)*r**2)**gamma_eos * (1.- 2.*mass1/r)**((gamma_eos + 1.)/4.))
end function efunc

#ifdef GR
!--- Lorentz factor in GR
real function gammafunc(r)
 use metric_tools,   only:get_metric3plus1
 use utils_gr,       only:dot_product_gr
 real, intent(in) :: r
 real :: x(3),v(3),v2,alpha,beta(3),gammaijdown(3,3),gammaijUP(3,3),gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
 x = (/r,0.,0./)
 v = (/vfunc(r),0.,0./)
 call get_metric3plus1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 v2 = dot_product_gr(v,v,gammaijdown)
 gammafunc = 1./sqrt(1.-v2)
end function gammafunc
#endif


!=====----------------------------------------------------------------------------------------------=====
!
!---Functions for non-relativistic solution
!
!--------------------------------------------------------------------------------------------------------
! See Barry, Parlange & Li 2000, for the analytic approximations used for the Lambert W function.
!
!--- Lambert W function for the principal branch (k=0) approaching from the negative
real function lambertw_0(x)
 real, intent(in) :: x
 real, parameter  :: exp1 = exp(1.)
 real :: eta, N1, N2
 eta = 2. + 2.*exp1*x
 N2  = 6. + 3.*Sqrt(2.) - ((-5764. - 4108.*Sqrt(2.) + (2237. + 1457.*Sqrt(2.))*exp1)*eta)/&
      (-796. - 430.*Sqrt(2.) + (215. + 199.*Sqrt(2.))*exp1)
 N1  = (1. - 1./Sqrt(2.))*(Sqrt(2.) + N2)
 lambertw_0 = -1 + Sqrt(eta)/(1 + (Sqrt(eta)*N1)/(Sqrt(eta) + N2))
end function lambertw_0

!--- Lambert W function for the k=-1 branch
real function lambertw_neg1(x)
 real, intent(in) :: x
 real, parameter  :: M1 = 0.3361, M2 = -0.0042, M3 = -0.0201
 real :: sigma
 sigma = -1. - Log(-x)
 lambertw_neg1 = -1. - sigma - (2.*(1. - 1./(1. + (M1*Sqrt(sigma)*(1. + exp(M3*Sqrt(sigma))*M2*sigma))/Sqrt(2.))))/M1
end function lambertw_neg1

!--- Function used in the non-rel solution of velocity [D(r)] (See: Cranmer 2004)
real function func(r,rc)
 real, intent(in) :: r,rc
 func = (rc/r)**4 * exp(4.*(1.-rc/r) - 1.)
end function func

end module setup
