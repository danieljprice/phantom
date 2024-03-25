!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for simple MHD sine wave decay test
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, io, mpidomain, part, physcon, prompting,
!   setup_params, unifdis
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  Sets up Bx = sin(2*pi*x), zero magnetic and velocity field otherwise
!
!  The amplitude should decay as exp(-eta * 4 * pi^2 * t)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd,periodic,igas
 use io,           only:master
 use physcon,      only:pi
 use prompting,    only:prompt
 use mpidomain,    only:i_belong
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: i,maxvxyzu,npartx,maxp
 real    :: bzero,przero,uuzero,gam1
 real    :: deltax,totmass
!
!--boundaries
!
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.0625,0.0625)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
!
!--setup parameters
!
 bzero = 1.0
 rhozero = 2.0
 przero = 10.0
 if (maxvxyzu >= 4) then
    gamma = 5./3.
    gam1 = gamma - 1.
    uuzero = przero/(gam1*rhozero)
    polyk = 0.
 else
    gamma = 1.
    polyk = przero/rhozero
    uuzero = 1.5*polyk
 endif

 print "(/,a)",' Setup for MHD sine problem...'

 npartx = 64
 print "(a,f6.2,a)",' (NB: should be multiple of ',(xmax-xmin)/(zmax-zmin),' given box dimensions)'
 call prompt('Enter number of particles in x ',npartx,0)
 deltax = dxbound/npartx

 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                  hfact,npart,xyzh,periodic,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(igas)

 do i=1,npart
    vxyzu(1,i) = 0.
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(2,i) = Bzero*sin(2.*pi*(xyzh(1,i)-xmin))
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup
