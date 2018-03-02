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
!   Setup for the MHD Vortex problem. The centrifugal and centripetal
!   accelerations from the rotational velocity are balanced with
!   those from the magnetic tension. The thermal and magnetic pressure
!   are in balance. At t=10, the vortex should return to its original
!   position unchanged.
!
!  REFERENCES:
!    Balsara, 2004, APJS, 151, 149-184
!    Dumbser, Balsara, Toro, Munz, 2008, J. Comp. Phys., 227, 8209-8253.
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
!  setup for MHD vortex problem in 3D
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd
 use io,           only:master
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
 real :: deltax,totmass !,dz,rfact,expo
 integer :: i,maxvxyzu,nx,maxp
 real :: u, k, const, halfsqrt2,rsq1,PplusdeltaP
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
 k = 1.0
 u = 1.0

 print "(/,a)",' Setup for Balsara (2004) MHD vortex problem...'

 polyk = 0.
 if (maxvxyzu < 4) stop 'need maxvxyzu=4 for MHD vortex setup'
!
!--boundaries
!
 call set_boundary(-5.0,5.0,-5.0,5.0,-5.0,5.0)

 use_closepacked = .true.
 nx = 64
 if (id==master) call prompt('Enter resolution (number of particles in x)',nx,8)

 call bcast_mpi(nx)
 deltax = dxbound/nx

 if (use_closepacked) then
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 else
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 endif
 npartoftype(:) = 0
 npartoftype(1) = npart

 rhozero = 1.0
 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(1)

 const = 1.0 / (2.0 * pi)
 halfsqrt2 = 0.5 * sqrt(2.0)
 do i=1,npart
    !!
    !! I believe this is the 3D version from Dumbser et al (2008), which includes lots of rotation factors
    !!

!    rfact = 1.0 - halfsqrt2*xyzh(1,i) - halfsqrt2*xyzh(3,i)
!    rsq1 = 1.0 - rfact*rfact * (xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))
!    expo = exp(0.5 * rsq1)
!
!    vxyzu(1,i) = halfsqrt2 - k * const * expo * halfsqrt2 * rfact * xyzh(2,i)
!    vxyzu(2,i) = 1         + k * const * expo * halfsqrt2 * rfact * (xyzh(1,i) -xyzh(3,i))
!    vxyzu(3,i) = halfsqrt2 + k * const * expo * halfsqrt2 * rfact * xyzh(2,i)
!
!    Bxyz(1,i) = -u * const * expo * halfsqrt2 * rfact * xyzh(2,i)
!    Bxyz(2,i) =  u * const * expo * halfsqrt2 * rfact * (xyzh(1,i) - xyzh(3,i))
!    Bxyz(3,i) =  u * const * expo * halfsqrt2 * rfact * xyzh(2,i)
!
!    PplusdeltaP = 1.0 + exp(rsq1) / (32.0 * pi**3) * (u*u*rsq1 - 4*k*k*pi)
!    vxyzu(4,i) = 1.5 * PplusdeltaP


    !!
    !! This is the original Balsara (2004) version with an added v_z = 1 component, otherwise identical.
    !!

    rsq1 = 1.0 - xyzh(1,i)*xyzh(1,i) - xyzh(2,i)*xyzh(2,i)

    vxyzu(1,i) = 1.0 - xyzh(2,i) * k * const * exp(0.5*rsq1)
    vxyzu(2,i) = 1.0 + xyzh(1,i) * k * const * exp(0.5*rsq1)
    vxyzu(3,i) = 1.0

    Bxyz(1,i) = -xyzh(2,i) * u * const * exp(0.5*rsq1)
    Bxyz(2,i) =  xyzh(1,i) * u * const * exp(0.5*rsq1)
    Bxyz(3,i) = 0.0

    PplusdeltaP = 1.0 + exp(rsq1) / (32.0 * pi**3) * (u*u*rsq1 - 0.5*k*k)
    vxyzu(4,i) = 1.5 * PplusdeltaP

 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup
