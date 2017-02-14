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
!   Setup for the Orszag-Tang Vortex problem in 3D
!
!  REFERENCES:
!    Orszag S. A., Tang C.-M., 1979, J. Fluid Mech., 90, 129
!    Dahlburg R. B., Picone J. M., 1989, Physics of Fluids B, 1, 2153
!    Picone J. M., Dahlburg R. B., 1991, Physics of Fluids B, 3, 29
!    Price D. J., Monaghan J. J., 2005, MNRAS, 364, 384
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, io, mpiutils, part, physcon, prompting,
!    setup_params, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for 3D Orszag-Tang vortex problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bevol,maxvecp,mhd,maxBevol
 use io,           only:master,real4
 use prompting,    only:prompt
 use mpiutils,     only:bcast_mpi
 use physcon,      only:pi
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: deltax,totmass,dz
 integer :: i,maxvxyzu,nx,maxp
 real :: machzero,betazero,const,vzero,bzero,przero,uuzero,gam1
 logical :: use_closepacked
!
!--general parameters
!
 time = 0.
 gamma = 5./3.
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
!
!--setup parameters
!
 const = 4.*pi
 betazero = 10./3.
 machzero = 1.0
 vzero = 1.0
 bzero = 1.0/sqrt(const)
 przero = 0.5*bzero**2*betazero
 rhozero = gamma*przero*machzero

 gam1 = gamma - 1.
 uuzero = przero/(gam1*rhozero)

 print "(/,a)",' Setup for 3D Orszag-Tang vortex problem...'
 print 10,betazero,machzero,bzero,rhozero,przero
10 format(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)

 polyk = 0.
 if (maxvxyzu < 4) stop 'need maxvxyzu=4 for orszag-tang setup'

 use_closepacked = .true.

 nx = 128
 if (id==master) call prompt('Enter resolution (number of particles in x)',nx,8)
 call bcast_mpi(nx)
 deltax = dxbound/nx
!
!--boundaries
!
 if (use_closepacked) then
    dz = 2.*sqrt(6.)/nx
 else
    dz = 6.*deltax
 endif
 call set_boundary(-0.5,0.5,-0.5,0.5,-dz,dz)

 if (use_closepacked) then
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 else
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 endif
 npartoftype(:) = 0
 npartoftype(1) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1,i) = -vzero*sin(2.*pi*(xyzh(2,i)-ymin))
    vxyzu(2,i) = vzero*sin(2.*pi*(xyzh(1,i)-xmin))
    vxyzu(3,i) = 0.
    vxyzu(4,i) = uuzero
    if (mhd) then
       Bevol(:,i) = 0.
       if (maxvecp==maxp) then
          if (maxBevol /= 3) stop 'euler potentials + orstang not implemented'
          Bevol(3,i) = real4(0.5/pi*Bzero*(cos(2.*pi*(xyzh(2,i)-ymin)) &
                                  + 0.5*cos(4.*pi*(xyzh(1,i)-xmin))))
       else
          Bevol(1,i) = real4(-Bzero*sin(2.*pi*(xyzh(2,i)-ymin)))
          Bevol(2,i) = real4(Bzero*sin(4.*pi*(xyzh(1,i)-xmin)))

          !--non-zero divB test
          !Bevol(1,i) = 0.5/pi*sin(2.*pi*(xyzh(1,i)-xmin))
       endif
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup

