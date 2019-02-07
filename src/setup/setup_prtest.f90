!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!  DEPENDENCIES: externalforces, io, options, physcon, setup_params,
!    spherical, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for a spherical blob of noninteracting particles on circular orbits
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master,warning
 use options,      only:ieos,iexternalforce,alpha,beta
 use externalforces,only:mass1,iext_prdrag
 use spherical,    only:set_sphere
 use units,        only:set_units
 use physcon,      only:km
 real, parameter :: pi = 3.1415926536
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real    :: totmass,totvol,psep,rmax,rmin,xyz_orig(3)
 real    :: r, vphi, phi
 integer :: i,np,nx,maxvxyzu

 call set_units(dist=10.*km,c=1.)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.
 rmin  = 0.
 rmax  = 0.05

!
!--setup particles
!
 np = 100 !size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 totvol = 4./3.*pi*rmax**3
 nx = int(np**(1./3.))
 psep = totvol**(1./3.)/real(nx)

 print*,' total volume = ',totvol,' particle separation = ',psep
 print*,totvol/psep**3

 npart = 0
 npart_total = 0

 xyz_orig(:) = (/0.,100.,0./)

 call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,xyzh, &
                  nptot=npart_total,xyz_origin=xyz_orig)
!
!--set particle properties
!
 totmass = 1.
 rhozero = totmass/totvol
 polyk   = 0.0
 print*,' total mass = ',totmass,' mean density = ',rhozero,' polyk = ',polyk

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 massoftype(1) = totmass/npart_total
 print*,' particle mass = ',massoftype(1)

 do i=1,npart
    r = sqrt(dot_product(xyzh(1:3,i), xyzh(1:3,i)))
    vphi = sqrt(mass1/r)
    phi = atan2(xyzh(2,i), xyzh(1,i))

    vxyzu(1,i) = -vphi*sin(phi)
    vxyzu(2,i) =  vphi*cos(phi)
    vxyzu(3,i) = 0.
    if (maxvxyzu==4) then
       vxyzu(4,i) = 0.
       call warning('setup_prtest','maxvxyzu should not be 4, so set temp=0 for you')
    endif
 enddo

!
! --- set input defaults for nonviscous particles in prdrag
!
 ieos = 2
 iexternalforce = iext_prdrag
 alpha=0.
 beta=0.


end subroutine setpart

end module setup

