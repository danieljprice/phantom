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
! Set up for exploding star problem
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, io, part, physcon, setup_params, timestep,
!    unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for exploding star problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use io,           only:master
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound
 use physcon,      only:pi
 use timestep,     only:tmax,dtmax
 use part,         only:mhd,Bevol
 use units,        only:set_units
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: deltax,totmass
 integer :: i,npartx
!
!--boundaries
!
 call set_boundary(-0.2,0.2,-0.2,0.2,-0.2,0.2)
!
!--units
!
 call set_units(G=1.d0)
!
! general parameters
!
 time = 0.
 hfact = 1.2
!
! setup particles
!
 npartx = 12
 deltax = dxbound/npartx
 rhozero = 1.0

 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)

 npartoftype(:) = 0
 npartoftype(1) = npart

 totmass = 1.
 massoftype = totmass/npart
 print*,' particle mass = ',massoftype(1)
!
! thermal energy
!
 polyk   = 0.0364
 gamma   = 1.
 do i=1,npart
    vxyzu(:,i) = 0.

    if (mhd) then
       Bevol(:,i) = 0.
       Bevol(3,i) = 2.38_4
    endif
 enddo
 ihavesetupB = .true.
!
!--set default runtime options for this setup
!
 tmax = 1.e-2
 dtmax = 1.e-4

end subroutine setpart

end module setup

