!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the MHD Vortex problem. The centrifugal and centripetal
!   accelerations from the rotational velocity are balanced with
!   those from the magnetic tension. The thermal and magnetic pressure
!   are in balance. At t=10, the vortex should return to its original
!   position unchanged.
!
! :References:
!    Balsara, 2004, APJS, 151, 149-184
!    Dumbser, Balsara, Toro, Munz, 2008, J. Comp. Phys., 227, 8209-8253.
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, io, mpidomain, mpiutils, part, physcon,
!   prompting, setup_params, unifdis
!
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
 use setup_params, only:rhozero,ihavesetupB,npart_total
 use slab,         only:set_slab
 use boundary,     only:dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd,igas,maxvxyzu
 use io,           only:master
 use prompting,    only:prompt
 use mpiutils,     only:bcast_mpi
 use physcon,      only:pi
 use kernel,       only:hfact_default
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
 integer :: i,nx
 real :: u, k, const, halfsqrt2,rsq1,PplusdeltaP
!
!--general parameters
!
 time = 0.
 gamma = 5./3.
 hfact = hfact_default
!
!--setup parameters
!
 k = 1.0
 u = 1.0

 print "(/,a)",' Setup for Balsara (2004) MHD vortex problem...'

 polyk = 0.
 if (maxvxyzu < 4) stop 'need maxvxyzu=4 for MHD vortex setup'

 nx = 64
 if (id==master) call prompt('Enter resolution (number of particles in x)',nx,8)
 call bcast_mpi(nx)

 deltax = dxbound/nx
 call set_slab(id,master,nx,-5.,5.,-5.,5.,deltax,hfact,npart,npart_total,xyzh,'closepacked')

 npartoftype(:) = 0
 npartoftype(igas) = npart

 rhozero = 1.0
 totmass = rhozero*dxbound*dybound*dzbound
 massoftype(igas) = totmass/npart_total
 print*,'npart = ',npart,' particle mass = ',massoftype(igas)

 const = 1.0 / (2.0 * pi)
 halfsqrt2 = 0.5 * sqrt(2.0)
 do i=1,npart
    !
    ! This is the original Balsara (2004) version with 
    ! an added v_z = 1 component, otherwise identical
    !
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
