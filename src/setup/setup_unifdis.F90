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
!   Setup routine for uniform distribution
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, io, mpiutils, part, physcon,
!    prompting, setup_params, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:use_dustfrac
 use setup_params, only:rhozero,npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:dustfrac,ndusttypes
 use prompting,    only:prompt
 use physcon,      only:pi
 use dust,         only:set_dustfrac
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: totmass,deltax,dust_to_gas(ndusttypes)
 integer :: i,j,maxp,maxvxyzu
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
 if (id==master) then
    npartx = 50
    print*,' uniform setup... (max = ',nint((maxp)**(1/3.)),')'
    call prompt('enter number of particles in x direction ',npartx,1)
 endif
 call bcast_mpi(npartx)
 deltax = dxbound/npartx

 if (id==master) then
    rhozero = 1.
    call prompt(' enter density (gives particle mass)',rhozero,0.)
 endif
 call bcast_mpi(rhozero)

 if (maxvxyzu < 4) then
    if (id==master) then
       polykset = 1.
       call prompt(' enter sound speed in code units (sets polyk)',polykset,0.)
    endif
    call bcast_mpi(polykset)
    polyk = polykset**2
    print*,' polyk = ',polyk
 else
    polyk = 0.
    polykset = 0.
 endif

 if (use_dustfrac) then
    dust_to_gas(:) = 1.e-2
    do i = 1,ndusttypes
       if (id==master) call prompt(' enter dust-to-gas ratio ',dust_to_gas(i),0.)
       call bcast_mpi(dust_to_gas(i))
    enddo
 endif

 if (id==master) then
    ilattice = 1
    call prompt(' select lattice type (1=cubic, 2=closepacked)',ilattice,1)
 endif
 call bcast_mpi(ilattice)

 npart = 0
 npart_total = 0

 select case(ilattice)
 case(1)
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 case(2)
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 case default
    print*,' error: chosen lattice not available, using cubic'
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 end select

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 print*,' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1:3,i) = 0.
 enddo

 if (use_dustfrac) then
    do i=1,npart
       do j = 1,ndusttypes
          call set_dustfrac(dust_to_gas(j),dustfrac(j,i))
       enddo
    enddo
 endif

end subroutine setpart

end module setup

