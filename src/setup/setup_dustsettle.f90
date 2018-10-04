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
!   Setup of dust settling problem from PL15
!
!  REFERENCES: Price & Laibe (2015), MNRAS 451, 5332
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, externalforces, io, mpiutils,
!    options, part, physcon, prompting, set_dust, setup_params,
!    table_utils, timestep, unifdis, units
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
 use setup_params,   only:npart_total,rhozero
 use io,             only:master
 use unifdis,        only:set_unifdis
 use boundary,       only:set_boundary,xmin,xmax,zmin,zmax,dxbound,dzbound
 use mpiutils,       only:bcast_mpi
 use part,           only:labeltype,set_particle_type,igas,dustfrac,&
                          ndusttypes,ndustsmall,grainsize,graindens
 use physcon,        only:pi,au,solarm
 use dim,            only:maxvxyzu,use_dust,maxp,maxdustsmall
 use prompting,      only:prompt
 use externalforces, only:mass1,Rdisc,iext_discgravity
 use options,        only:iexternalforce,use_dustfrac
 use timestep,       only:dtmax,tmax
 use units,          only:set_units,udist,umass
 use dust,           only:init_drag,idrag,grainsizecgs,graindenscgs,get_ts
 use set_dust,       only:set_dustfrac,set_dustbinfrac
 use table_utils,    only:logspace
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: totmass,deltax,dz,length
 integer :: i,iregime,ierr
 integer :: itype,npartx
 integer :: npart_previous
 real    :: H0,HonR,omega,ts(maxdustsmall),dustbinfrac(maxdustsmall)
 real    :: xmini,xmaxi,ymaxdisc,cs,t_orb
 real    :: dtg,smincgs,smaxcgs,sindex

!
! default options
!
 npartx = 32
 rhozero = 1.e-3
 dtg = 0.01
 grainsize = 0.
 graindens = 0.
 grainsizecgs = 0.1
 graindenscgs = 3.
 ndustsmall = 1
 smincgs = 1.e-5
 smaxcgs = 1.
 sindex = 3.5

 if (id==master) then
    itype = igas
    print "(/,a,/)",'  >>> Setting up particles for dust settling test <<<'
    call prompt(' enter number of '//trim(labeltype(itype))//' particles in x ',npartx,8,maxp/144)
 endif
 call bcast_mpi(npartx)
 if (id==master) call prompt('enter '//trim(labeltype(itype))//&
                      ' midplane density (gives particle mass)',rhozero,0.)
 call bcast_mpi(rhozero)
 call set_units(dist=10.*au,mass=solarm,G=1.)
 if (use_dust) then
    !--currently assume one fluid dust
    use_dustfrac = .true.
    if (id==master) then
       call prompt('Enter total dust to gas ratio',dtg,0.)
       call prompt('How many grain sizes do you want?',ndustsmall,1,maxdustsmall)
       ndusttypes = ndustsmall
       if (ndusttypes > 1) then
          !--grainsizes
          call prompt('Enter minimum grain size in cm',smincgs,0.)
          call prompt('Enter maximum grain size in cm',smaxcgs,0.)
          call logspace(grainsize(1:ndusttypes),smincgs,smaxcgs)
          grainsize(1:ndusttypes) = grainsize(1:ndusttypes)/udist
          !--mass distribution
          call prompt('Enter power-law index, e.g. MRN',sindex)
          call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes))
          !--grain density
          call prompt('Enter grain density in g/cm^3',graindens(1),0.)
          graindens(1:ndusttypes) = graindens(1)/umass*udist**3
       else
          call prompt('Enter grain size in cm',grainsizecgs,0.)
          call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
          grainsize(1) = grainsizecgs/udist
          graindens(1) = graindenscgs/umass*udist**3
       endif
    endif
    call bcast_mpi(dtg)
    call bcast_mpi(ndustsmall)
    call bcast_mpi(ndusttypes)
    call bcast_mpi(smincgs)
    call bcast_mpi(smaxcgs)
    call bcast_mpi(sindex)
    call bcast_mpi(grainsize)
    call bcast_mpi(graindens)
    call bcast_mpi(dustbinfrac)
    call bcast_mpi(grainsizecgs)
    call bcast_mpi(graindenscgs)
 endif
!
! general parameters
!
 HonR  = 0.05
 Rdisc = 5.
 mass1 = 1.
 H0    = HonR*Rdisc
 omega = sqrt(mass1/Rdisc**3)
 t_orb = 2.*pi/omega
 cs    = H0*omega
 time   = 0.
 iexternalforce = iext_discgravity
 dtmax = 0.1*t_orb
 tmax  = 15.*t_orb
!
! equation of state
!
 if (maxvxyzu >= 4) then
    gamma = 5./3.
 else
    gamma  = 1.
    polyk  = cs**2
 endif
!
! get stopping time information
!
 call init_drag(ierr)
 do i=1,ndustsmall
    call get_ts(idrag,grainsize(i),graindens(i),rhozero,0.0*rhozero,cs,0.,ts(i),iregime)
    print*,'s (cm) =',grainsize(i),'   ','St = ts * Omega =',ts(i)*omega
 enddo
!
! boundaries
!
 xmini = -0.25
 xmaxi =  0.25
 length = xmaxi - xmini
 deltax = length/npartx
 dz = 2.*sqrt(6.)/npartx
 !deltay = fac*deltax*sqrt(0.75)
 call set_boundary(xmini,xmaxi,-10.*H0,10.*H0,-dz,dz)

 npart = 0
 npart_total = 0
 npartoftype(:) = 0

 !--only works for one-fluid dust
 itype = igas

 !--get total mass from integration of density profile
 ymaxdisc = 3.*H0
 totmass = 2.*rhozero*sqrt(0.5*pi)*H0*erf(ymaxdisc/(sqrt(2.)*H0))*dxbound*dzbound

 npart_previous = npart

 call set_unifdis('closepacked',id,master,xmin,xmax,-ymaxdisc,ymaxdisc,zmin,zmax,deltax, &
                   hfact,npart,xyzh,nptot=npart_total,rhofunc=rhofunc,dir=2)

 !--set which type of particle it is
 do i=npart_previous+1,npart
    call set_particle_type(i,itype)

    vxyzu(:,i) = 0.

    !--set internal energy if necessary
    if (maxvxyzu >= 4) then
       if (gamma > 1.) then
          vxyzu(4,i) = cs**2/(gamma-1.)
       else
          vxyzu(4,i) = 1.5*cs**2
       endif
    endif

    !--one fluid dust: set dust fraction on gas particles
    if (use_dustfrac) then
       call set_dustfrac(dtg,dustfrac(:,i))
    else
       dustfrac(:,i) = 0.
    endif
 enddo

 npartoftype(itype) = npart - npart_previous
 if (id==master) print*,' npart = ',npart,npart_total

 massoftype(itype) = totmass/npartoftype(itype)*(1. + dtg)
 if (id==master) print*,' particle mass = ',massoftype(itype)

contains

real function rhofunc(x)
 real, intent(in) :: x

 rhofunc = exp(-0.5*(x/H0)**2)

end function rhofunc

end subroutine setpart

end module setup
