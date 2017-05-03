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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, io, mpiutils, part, physcon, setup_params,
!    unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: deltax,rhozero,polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:idust,labeltype,set_particle_type
 use units,        only:udist,umass,unit_density,set_units
 use physcon,      only:pi,solarm,au
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: totmass,deltax
 real :: grainsize_ini, graindens_ini
 real :: grainsizecgs = 1.
 real :: graindenscgs = 1.
 integer :: ipart,i,maxp,maxvxyzu
 integer :: itype,ntypes
 integer :: npart_previous
!
! units
!
 call set_units(mass=solarm,dist=au,G=1.d0)
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
 npart = 0
 npart_total = 0
 npartoftype(:) = 0

 ntypes = 2
 overtypes: do itype=1,ntypes
    if (itype /= 2) then
       print*,'TEST',itype,ntypes
       if (id==master) then
          if (ntypes > 1) then
             print "(/,a,/)",'  >>> Setting up '//trim(labeltype(itype))//' particles <<<'
          endif
          print*,' uniform cubic setup...'
          print*,' enter number of particles in x (max = ',nint((maxp)**(1/3.))/real(ntypes),')'
          read*,npartx
       endif
       call bcast_mpi(npartx)
       deltax = dxbound/npartx

       if (id==master) then
          print*,' enter density (gives particle mass)'
          read*,rhozero
       endif
       call bcast_mpi(rhozero)

       if (itype==1) then
          if (maxvxyzu < 4) then
             if (id==master) then
                print*,' enter sound speed in code units (sets polyk)'
                read*,polykset
             endif
             call bcast_mpi(polykset)
             polyk = polykset**2
             print*,' polyk = ',polyk
          else
             polyk = 0.
             polykset = 0.
          endif
       endif

       if (id==master) then
          print*,' select lattice type (1=cubic, 2=closepacked)'
          read*,ilattice
       endif
       call bcast_mpi(ilattice)

       npart_previous = npart

       select case(ilattice)
       case(1)
          call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                             hfact,npart,xyzh,nptot=npart_total)
       case(2)
          call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                             hfact,npart,xyzh,nptot=npart_total)
       case default
          print*,' error: chosen lattice not available, using cubic'
          call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                             hfact,npart,xyzh,nptot=npart_total)
       end select

       !--set which type of particle it is
       do i=npart_previous+1,npart
          call set_particle_type(i,itype)
       enddo

       do i=npart_previous+1,npart
          vxyzu(1:3,i) = 0.
       enddo

       npartoftype(itype) = npart - npart_previous
       print*,' npart = ',ipart,npart,npart_total

       totmass = rhozero*dxbound*dybound*dzbound
       massoftype(itype) = totmass/npartoftype(itype)
       print*,' particle mass = ',massoftype(itype)
    endif
 enddo overtypes

 npart = npart + 1
 i     = npart
!--setup the quantities for one spherical dust particle

 if (id==master) then
    itype = idust
    call set_particle_type(i,itype)
    print*,'WARNING: be sure that the values of graindenscgs and grainsizecgs are consistent with dust.F90'
    grainsize_ini = grainsizecgs/udist
    graindens_ini = graindenscgs/unit_densisty
!     massoftype(itype) = 4./3.*pi*graindens_ini*grainsize_ini**3
    massoftype(itype) = tiny(0.)
    print*,' grain mass = ',massoftype(itype)
    xyzh(1:3,i) = 0.
    xyzh(4,i) = hfact*deltax
    print*,' enter grains velocity in code units'
    read*, vxyzu(1,i)
    vxyzu(2:3,i) = 0.
    npartoftype(itype) = 1
 endif
 print*,'npart',npart

 return
end subroutine setpart

end module setup

