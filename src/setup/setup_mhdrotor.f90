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
!   Setup for MHD rotor problem
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, io, mpiutils, part, physcon, prompting,
!    setup_params, timestep, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxp,maxvxyzu,mhd
 use setup_params, only:npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:igas,Bxyz
 use prompting,    only:prompt
 use physcon,      only:pi
 use timestep,     only:tmax,dtmax
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: totvol,totmass,deltax,deltadisk,przero,vzero,denszero,densdisk
 real :: rdisk,radius,dx(3),Bzero(3),xorigin(3),const,densi
 integer :: i,nx
!
!--general parameters
!
 time = 0.
 gamma = 1.4
!
!--set particles
!
 if (id==master) then
    nx = 64
    call prompt('enter number of particles in x direction ',nx,1,int(sqrt(maxp/12.)))
 endif
 call bcast_mpi(nx)
!
!--boundary
!
 call set_boundary(-0.5,0.5,-0.5,0.5,-2.*sqrt(6.)/nx,2.*sqrt(6.)/nx)
 deltax = dxbound/nx
!
!--problem parameters
!
 const = sqrt(4.*pi)
 xorigin = 0.
 rdisk = 0.1             ! radius of the initial disk
 vzero = 2.0             ! rotation speed of initial disk
 Bzero(:) = (/5.0/const, 0., 0./)       ! uniform field in bx direction
 przero = 1.0              ! initial pressure
 denszero = 1.0            ! ambient density
 densdisk = 10.0           ! density of rotating disk

 if (id==master) then
    print "(/,a)",' MHD rotor problem '
    print 10,densdisk,rdisk,bzero(1),vzero,przero
 endif
10 format(/,' central density  = ',f10.3,', disk radius = ',f6.3,/, &
            ' initial Bx       = ',f10.3,',    rotation = ',f6.3,/, &
            ' pressure         = ',f10.3,/)

 npart = 0
 npart_total = 0
 deltadisk = deltax*(denszero/densdisk)**(1./3.)

 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,nptot=npart_total,rcylmin=rdisk)
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltadisk,hfact,npart,xyzh,nptot=npart_total,rcylmax=rdisk)

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 totvol = (dxbound*dybound - pi*rdisk**2)*dzbound
 totmass = (denszero*dxbound*dybound + densdisk*pi*rdisk**2)*dzbound
 massoftype(igas) = totmass/npart_total
 print*,' particle mass = ',massoftype(igas)

 do i=1,npart
    dx = xyzh(1:3,i) - xorigin(:)
    radius = sqrt(dot_product(dx(1:2),dx(1:2)))
    if (radius <= rdisk) then
       vxyzu(1,i) = -vzero*dx(2)/rdisk
       vxyzu(2,i) =  vzero*dx(1)/rdisk
       densi = densdisk
    else
       vxyzu(:,i) = 0.
       densi = denszero
    endif
    if (maxvxyzu >= 4) vxyzu(4,i) = przero/((gamma - 1.)*densi)
    if (mhd) Bxyz(1:3,i) = Bzero(:)
 enddo
 if (mhd) ihavesetupB = .true.
!
! input file parameters
!
 tmax = 0.15
 dtmax = 0.01

end subroutine setpart

end module setup

