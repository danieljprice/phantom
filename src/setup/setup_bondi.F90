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
!
!
!  REFERENCES:
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!
!  DEPENDENCIES:
!+
!--------------------------------------------------------------------------
module setup
 use physcon,        only:pi
 use externalforces, only:accradius1,accradius1_hard
 use dim,            only:gr
#ifdef GR
 use metric,         only:mass1,imetric
 use metric_tools,   only:imet_schwarzschild
#else
 use externalforces, only:mass1
#endif
 use setup_params,   only:rhozero,npart_total
 use io,             only:master,fatal
 use spherical,      only:set_sphere
 use options,        only:ieos,iexternalforce,nfulldump
 use timestep,       only:tmax,dtmax
 use centreofmass,   only:reset_centreofmass
 use units,          only:udist,umass,utime,set_units
 use physcon,        only:pc,solarm,gg
 use part,           only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,igas,set_particle_type,iboundary
 use stretchmap,     only:get_mass_r
 use kernel,         only:radkern
 use prompting,      only:prompt
 use bondiexact,     only:get_bondi_solution,rcrit

 implicit none

 public :: setpart

 private

 real :: gamma_eos

contains

!----------------------------------------------------------------
!+
!  setup for bondi accretion
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: vol,rmax,rmin,psep,tff,gcode,rhor,vr,ur
 real    :: r,pos(3),cs2,totmass,approx_m,approx_h
 integer :: i,np,nx,maxvxyzu,nbound

!-- Set code units
 call set_units(G=1.d0,c=1.d0)
 gcode = gg*umass*utime**2/udist**3
 print*,' gcode = ',gcode

 maxvxyzu = size(vxyzu(:,1))

!--Set general parameters
 time           = 0.
 iexternalforce = 1

#ifdef GR
 if (imetric/=imet_schwarzschild) call fatal('setup_bondi','You are not using the Schwarzschild metric.')
#endif

 if (gr) then
    ieos  = 2
    gamma = 5./3.
    polyk = 1.
 else
    gamma = 1.
    ieos  = 1
    cs2   = mass1/(2.*rcrit)
    polyk = cs2
 endif

 gamma_eos      = gamma             ! Note, since non rel bondi is isothermal, solution doesn't depend on gamma

 accradius1      = 0.
 accradius1_hard = 0.

 rmin = 2.2!2.
 rmax = 20.!10.                    !20.!18.1
 if (gr) then
    rmin = rmin*mass1
    rmax = rmax*mass1
 endif
 np   = 50*1000
 call prompt(' Enter the desired number of particles: ',np,0)

 vol  = 4./3.*pi*(rmax**3 - rmin**3)
 nx   = int(np**(1./3.))
 psep = vol**(1./3.)/real(nx)

 totmass  = get_mass_r(rhofunc,rmax,rmin)
 approx_m = totmass/np
 approx_h = hfact*(approx_m/rhofunc(rmin))**(1./3.)
 rhozero  = totmass/vol

 tff   = sqrt(3.*pi/(32.*rhozero))
 tmax  = 10.*tff
 dtmax = tmax/150.

 print*,''
 print*,' Setup for gas: '
 print*,' min,max radius = ',rmin,rmax
 print*,' volume         = ',vol        ,' particle separation = ',psep
 print*,' vol/psep**3    = ',vol/psep**3,' totmass             = ',totmass
 print*,' free fall time = ',tff        ,' tmax                = ',tmax
 print*,''

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
    call get_bondi_solution(rhor,vr,ur,r,mass1,gamma)
    vxyzu(1:3,i) = vr*pos/r
    vxyzu(4,i)   = ur

    if (r + radkern*xyzh(4,i)>rmax .or. r - radkern*xyzh(4,i)<rmin) then
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


real function rhofunc(r)
 real, intent(in) :: r
 real :: rho,v,u
 call get_bondi_solution(rho,v,u,r,mass1,gamma_eos)
 rhofunc = rho
end function rhofunc

end module setup
