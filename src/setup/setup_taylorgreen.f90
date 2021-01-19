!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup for Taylor-Green Vortex test (used in Phantom paper)
!
! :References: Price et al. (2017)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, domain, io, mpiutils, part, physcon, prompting,
!   setup_params, unifdis
!
 implicit none
 public :: setpart

 real,    private :: polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for the Taylor-Green problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use physcon,      only:pi
 use prompting,    only:prompt
 use domain,       only:i_belong
 use part,         only:periodic
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: totmass,deltax,vzero,dz
 integer :: ipart,i,maxp,maxvxyzu,nx
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 nx = 128
 if (id==master) then
    print *, ''
    call prompt('Enter resolution (number of particles in x)',nx,8, nint((maxp)**(1/3.)))
 endif
 call bcast_mpi(nx)
 deltax = dxbound/nx

 dz = 2.*sqrt(6.)/nx
 call set_boundary(0.,1.,0.,1.,-dz,dz)

 rhozero = 1.
 if (maxvxyzu < 4) then
    polyk = 1.
    print*,' polyk = ',polyk
 else
    polyk = 0.
    polykset = 0.
 endif
 npart = 0
 npart_total = 0

 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                  hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',ipart,npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 print*,' particle mass = ',massoftype(1)

 vzero = 0.1

 do i=1,npart
    vxyzu(1:3,i) = 0.
    !--velocity field for the Taylor-Green problem
    vxyzu(1,i) =  vzero*sin(2.*pi*xyzh(1,i))*cos(2.*pi*xyzh(2,i))
    vxyzu(2,i) = -vzero*cos(2.*pi*xyzh(1,i))*sin(2.*pi*xyzh(2,i))
 enddo

end subroutine setpart

end module setup
