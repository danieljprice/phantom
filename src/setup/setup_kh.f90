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
!   Setup for Kelvin-Helmholtz instability from Robertson et al. (2010)
!
!  REFERENCES:
!   Robertson et al. (2010), MNRAS 401, 2463-2476
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, io, mpiutils, options, part, physcon, prompting,
!    setup_params, timestep, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private
 real, parameter :: rho1 = 1., rho2 = 2.

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:npart_total
 use io,           only:master
 use options,      only:nfulldump
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:igas
 use prompting,    only:prompt
 use physcon,      only:pi
 use timestep,     only:dtmax,tmax
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=26)                :: filename
 logical :: iexist
 integer :: i,maxp,maxvxyzu,npartx
 real :: totmass,deltax,przero,v1,v2
!
!--general parameters
!
 time = 0.
 gamma = 5./3
 filename= trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax      = 2.00
    dtmax     = 0.1
    nfulldump = 1
 endif
!
!--set particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 if (id==master) then
    npartx = 64
    call prompt('enter number of particles in x direction ',npartx,1,nint(sqrt(maxp/12.)))
 endif
 call bcast_mpi(npartx)
 deltax = dxbound/npartx
!
!--boundary
!
 call set_boundary(0.,1.,0.,1.,-2.*sqrt(6.)/npartx,2.*sqrt(6.)/npartx)

 npart = 0
 npart_total = 0
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,nptot=npart_total,rhofunc=rhofunc,dir=2)

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 totmass = 0.5*dxbound*dzbound*rho2 + 0.5*dxbound*dzbound*rho1
 massoftype(igas) = totmass/npart_total
 print*,' particle mass = ',massoftype(1)

 v1 = -0.5
 v2 = 0.5
 przero = 2.5
 do i=1,npart
    vxyzu(1,i) = v1 + Rfunc(xyzh(2,i))*(v2 - v1)
    vxyzu(2,i) = 0.1*sin(2.*pi*xyzh(1,i))
    vxyzu(3,i) = 0.
    if (maxvxyzu > 3) then
       vxyzu(4,i) = przero/((gamma - 1.)*rhofunc(xyzh(2,i)))
    endif
 enddo

end subroutine setpart

real function rhofunc(y)
 real, intent(in) :: y
 real, parameter :: rho1 = 1., rho2 = 2.

 rhofunc = rho1 + Rfunc(y)*(rho2 - rho1)

end function rhofunc

real function Rfunc(y)
 real, parameter  :: delta = 0.05
 real, intent(in) :: y
 real :: fac1,fac2

 fac1 = (1. - 1./(1. + exp(2.*(y-0.25)/delta)))
 fac2 = (1. - 1./(1. + exp(2.*(0.75-y)/delta)))
 Rfunc = fac1*fac2

end function Rfunc

end module setup

