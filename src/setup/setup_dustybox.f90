!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the dustybox problem in dust-gas mixtures
!
! :References: Laibe & Price (2011), MNRAS 418, 1491
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, io, mpidomain, mpiutils, part, physcon,
!   prompting, setup_params, unifdis, units
!
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: deltax,polykset
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
 use part,         only:labeltype,set_particle_type,igas,idust,periodic
 use physcon,      only:pi,solarm,au
 use units,        only:set_units
 use prompting,    only:prompt
 use mpidomain,    only:i_belong
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: totmass
 integer :: i,j,maxp,maxvxyzu
 integer :: itype,ntypes
 integer :: npart_previous
 logical, parameter :: ishift_box =.true.
!
! units (needed if physical drag is used)
!
 call set_units(mass=solarm,dist=au,G=1.d0)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.
 rhozero = 1.
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 npart = 0
 npart_total = 0
 npartoftype(:) = 0

 ntypes = 2
 overtypes: do j=1,ntypes
    if (j==1) then
       itype = igas
    else
       itype = idust + (j-2)
    endif
    if (id==master) then
       if (itype==1) npartx = 64
       if (ntypes > 1) then
          print "(/,a,/)",'  >>> Setting up '//trim(labeltype(itype))//' particles <<<'
       endif
       call prompt('enter number of particles in x direction ',npartx,1)
    endif
    call bcast_mpi(npartx)
    deltax = dxbound/npartx

    if (id==master) call prompt(' enter density (gives particle mass)',rhozero,0.)
    call bcast_mpi(rhozero)

    if (itype==1) then
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
    endif

    if (id==master) then
       ilattice = 1
       call prompt(' select lattice type (1=cubic, 2=closepacked)',ilattice,1)
    endif
    call bcast_mpi(ilattice)

    npart_previous = npart

    select case(ilattice)
    case(1)
       if (ishift_box .eqv. .false.) then
          call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                            hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
       else
          if (itype == igas) then
             call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                               hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
          else
             call set_unifdis('cubic',id,master,xmin+0.01*deltax,xmax+0.01*deltax,ymin+0.05*deltax, &
                              ymax+0.05*deltax,zmin+0.05*deltax,zmax+0.05*deltax,deltax, &
                              hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
             !call set_unifdis('cubic',id,master,xmin+0.5*deltax,xmax+0.5*deltax,ymin+0.5*deltax, &
             !                  ymax+0.5*deltax,zmin+0.5*deltax,zmax+0.5*deltax,deltax, &
             !                   hfact,npart,xyzh,nptot=npart_total)
             !--Use this previous setup to how bad the spline is
          endif
       endif
    case(2)
       call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                         hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
    case default
       print*,' error: chosen lattice not available, using cubic'
       call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                         hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
    end select

    !--set which type of particle it is
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
    enddo

    do i=npart_previous+1,npart
       if (itype==igas) then
          vxyzu(1,i)   = 0.
          vxyzu(2:3,i) = 0.
       else
          vxyzu(1,i)   = 1.
          vxyzu(2:3,i) = 0.
       endif
    enddo

    npartoftype(itype) = npart - npart_previous
    print*,' npart = ',npart,npart_total

    totmass = rhozero*dxbound*dybound*dzbound
    massoftype(itype) = totmass/npartoftype(itype)
    print*,' particle mass = ',massoftype(itype)

 enddo overtypes

end subroutine setpart

end module setup
