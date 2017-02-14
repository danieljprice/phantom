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
!  OWNER: Terrence Tricco
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, io, part, physcon, setup_params, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  Sets up Bx = sin(2*pi*x), zero magnetic and velocity field otherwise
!
!  The amplitude should decay as exp(-eta * 4 * pi^2 * t)
!
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bevol,maxvecp,mhd,maxBevol
 use io,           only:master
 use physcon,      only:pi
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: deltax,totmass
 integer :: i,maxvxyzu,npartx,maxp
 real :: machzero,betazero,const,vzero,bzero,przero,uuzero,gam1
!
!--boundaries
!
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.0625,0.0625)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 5./3.
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

 gam1 = gamma - 1.
 uuzero = przero/(gam1*rhozero)

 print "(/,a)",' Setup for MHD sine problem...'

 polyk = 0.
 if (maxvxyzu < 4) stop 'need maxvxyzu=4 for orszag-tang setup'

 print*,'Enter number of particles in x (max = ',nint((maxp)**(1/3.)),')'
 print "(a,f6.2,a)",' (NB: should be multiple of ',(xmax-xmin)/(zmax-zmin),' given box dimensions)'
 read*,npartx
 deltax = dxbound/npartx

 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)

 npartoftype(:) = 0
 npartoftype(1) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1,i) = 0.
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    vxyzu(4,i) = uuzero
    if (mhd) then
       Bevol(:,i) = 0.
       if (maxvecp==maxp) then
          stop 'euler potentials + mhd sine not implemented'
       else
          Bevol(2,i) = Bzero*sin(2.*pi*(xyzh(1,i)-xmin))
       endif
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup

