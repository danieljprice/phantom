!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Sedov blast wave problem with dust
!
! :References:
!   Laibe & Price (2012a), MNRAS 420, 2345
!   Laibe & Price (2012b), MNRAS 420, 2365
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

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use unifdis,      only:set_unifdis
 use io,           only:master
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:labeltype,set_particle_type,igas,idust,periodic
 use units,        only:umass,utime,unit_density,udist,set_units
 use physcon,      only:pc,solarm,pi
 use prompting,    only:prompt
 use mpidomain,    only:i_belong
 use mpiutils,     only:bcast_mpi
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: i,maxp,maxvxyzu,itype,npartx,npart_previous,ifluid,nfluid
 real    :: deltax,totmass,dust_to_gas_ratio
 real    :: rblast,prblast,enblast,gam1,spsoundzero
 real(8) :: uenergy

 call set_units(dist=pc,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.
 dust_to_gas_ratio = 1.d-2
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 npart = 0
 npart_total = 0
 npartoftype(:) = 0

 print*,' uniform cubic setup for dusty Sedov blast wave problem...'

 nfluid = 2
 overtypes: do ifluid=1,nfluid
    select case(ifluid)
    case(1)
       itype = igas
    case(2)
       itype = idust
    end select
    if (id==master) then
       if (ifluid==1) npartx = 64
       if (nfluid > 1) print "(/,a,/)",'  >>> Setting up '//trim(labeltype(itype))//' particles <<<'
       call prompt('enter number of particles in x direction ',npartx,1)
    endif
    call bcast_mpi(npartx)
    deltax = dxbound/npartx

    rhozero = 1.0
    if (itype==idust) then
       rhozero = rhozero*dust_to_gas_ratio
    endif
    print*,' density in code units = ',rhozero
    print*,' density in physical units = ',rhozero*unit_density,' g/cm^3'
    polyk = 0.
    if (maxvxyzu < 4) stop 'need maxvxyzu=4 for sedov setup'

    enblast = 1.0
    print*,' energy in code units = ',enblast
    uenergy = umass*(udist/utime)**2
    print*,' energy in physical units = ',enblast*uenergy

    print*,' time unit = ',utime,' s'
    rblast = 2.*hfact*deltax
    gamma = 5./3.
    gam1 = gamma - 1.
    prblast = gam1*enblast/(4./3.*pi*rblast**3)
    print*,' speed of sound in blast (code units) = ',sqrt(gamma*prblast/rhozero)

    spsoundzero = 2.d4/(udist/utime)
    print*,' speed of sound = ',2.d4,' cm/s, in code units = ',spsoundzero

    npart_previous = npart

    if (itype==igas) then
       call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                               hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
    else
       call set_unifdis('cubic',id,master,xmin+0.5*deltax,xmax+0.5*deltax,ymin+0.5*deltax, &
                              ymax+0.5*deltax,zmin+0.5*deltax,zmax+0.5*deltax,deltax, &
                              hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
    endif

    !--set which type of particle it is
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
    enddo

    do i=npart_previous+1,npart
       vxyzu(:,i) = 0.
       if (itype == igas) then
          if ((xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2) < rblast*rblast) then
             vxyzu(4,i) = prblast/(gam1*rhozero)
          else
             vxyzu(4,i) = spsoundzero/(gam1*gamma)
          endif
       endif
    enddo

    npartoftype(itype) = npart - npart_previous
    print*,' npart = ',npartoftype(itype),npart_total

    totmass = rhozero*dxbound*dybound*dzbound
    massoftype(itype) = totmass/npartoftype(itype)
    print*,' particle mass = ',massoftype(itype)

 enddo overtypes

end subroutine setpart

end module setup
