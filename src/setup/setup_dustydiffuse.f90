!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION: Setup for the Sedov blast wave problem
!
!  REFERENCES: None
!
!  OWNER: Mark Hutchison
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, eos, io, kernel, mpiutils, part,
!    prompting, setup_params, timestep, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------
!+
!  3D dust diffusion test from Price & Laibe (2015)
!+
!----------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero
 use unifdis,      only:set_unifdis
 use io,           only:master,fatal,iverbose
 use boundary,     only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax, &
                        set_boundary
 use timestep,     only:tmax,dtmax
 use prompting,    only:prompt
 use kernel,       only:hfact_default
 use part,         only:igas,dustfrac,Bevol,set_particle_type
 use mpiutils,     only:bcast_mpi,reduceall_mpi
 use dim,          only:maxp,maxvxyzu,maxtypes,mhd,ndusttypes
 use eos,          only:ieos
 use dust,         only:K_code,idrag
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: npartx,i
 integer :: eps_type
 real    :: deltax,totmass
 real    :: epstot,epsi(ndusttypes),rc,rc2,r2
 !
 !--general parameters
 !
 time = 0.
 hfact = hfact_default
 !
 ! setup uniform box
 !
 npartx = 32
 if (id==master) then
    print*,'Setup for DUSTYDIFFUSE problem...'
    call prompt(' Enter number of particles in x ',npartx,8,nint((maxp)**(1/3.)))
 endif
 call bcast_mpi(npartx)
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
 deltax = dxbound/npartx
 !
 ! runtime options
 !
 K_code = 0.1
 ieos = 1
 idrag = 3
 polyk = 1.
 gamma = 1.
 iverbose = 0
 rhozero = 3.

! if (maxvxyzu < 4) call fatal('setup','only evolve dustfrac so no need to bother with energy')

 npart = 0
 npartoftype(:) = 0
 iverbose = 2
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,verbose=.false.)
 npartoftype(igas) = npart
 npartoftypetot(igas) = reduceall_mpi('+',npartoftype(igas))
 totmass = rhozero*dxbound*dybound*dzbound
 massoftype(igas)  = totmass/npartoftypetot(igas)
 if (id==master) print*,' particle mass = ',massoftype(igas)

 vxyzu = 0.
 if (mhd) Bevol = 0.

 !
 ! setup dust fraction in the box
 !
 epstot = 0.1
 eps_type = 2
 select case(eps_type)
 case(1)
    !--Equal dust fractions
    epsi(:) = epstot/real(ndusttypes)
 case(2)
    !--Unequal dust fractions
    do i=1,ndusttypes
       epsi(i) = 1./real(i)
    enddo
    epsi = epstot/sum(epsi)*epsi
 case default
    stop 'eps_type not valid!'
 end select

 !--check that individual dust fractions add up to the total dust fraction
 if (abs(sum(epsi)-epstot)/epstot>1.e-14) then
    write(*,"(/,a)") 'ERROR! SUM(epsilon_k) /= epsilon'
    print*,'SUM(epsilon_k) = ',sum(epsi)
    print*,'       epsilon = ',epstot
 endif

 rc   = 0.25
 rc2  = rc**2
 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (r2 < rc2) then
       dustfrac(:,i) = epsi(:)*(1. - r2/rc2)
    else
       dustfrac(:,i) = 0.
    endif
    call set_particle_type(i,igas)
 enddo

 !
 ! evolve only the dust fraction with a simple predictor-corrector scheme
 !
 tmax  = 10.
 dtmax = 0.05

end subroutine setpart

end module setup


