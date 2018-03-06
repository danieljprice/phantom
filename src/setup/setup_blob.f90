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
!   Setup for the Agertz et al. (2007) evaporating blob problem
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, io, kernel, physcon, prompting, setup_params,
!    timestep, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for blob problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use unifdis,      only:set_unifdis
 use io,           only:master
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use physcon,      only:pi
 use kernel,       only:hfact_default
 use timestep,     only:dtmax,tmax
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: i,maxp,maxvxyzu,npartx,npartmed
 real    :: deltax,totmass,totvol
 real    :: rcloud,przero,denszero,denscloud,vzero,gam1,deltaxcloud
 real    :: velmachno,spsound,taukh,tcrush
 !
 ! general parameters
 !
 time = 0.
 hfact = hfact_default
 !
 ! problem-specific parameters
 !
 gamma = 1.66666666666667
 rcloud = 0.1
 przero = 1.0
 denszero  = 1.0
 denscloud = 10.0
 velmachno = 2.7
 gam1 = gamma - 1.
 spsound = sqrt(gamma*przero/denszero)
 vzero = velmachno*spsound
 rhozero = denszero
 polyk = 0. ! not necessary, we require maxvxyzu==4 for this test

 tcrush = 2.*rcloud*sqrt(denscloud/denszero)/vzero
 taukh = 1.6*tcrush
 !
 ! set box size
 !
 call set_boundary(-5.*rcloud,15.*rcloud,-5.*rcloud,5.*rcloud,-5.*rcloud,5.*rcloud)
 !
 ! setup particles
 !
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))

 write(*,*) 'Evaporating blob problem '
 if (maxvxyzu < 4) stop 'need maxvxyzu=4 for blob setup'
 write(*,"(/,'   cloud density = ',f10.3,',  R_cloud = ',f6.3,/,"// &
            "' ambient density = ',f10.3,', pressure = ',f6.3,/,"// &
            "'     mach number = ',f10.3,',      c_s = ',f6.3,/,"// &
            "'          tau_kh = ',f10.3,',      v_0 = ',f6.3,/)") &
       denscloud,rcloud,denszero,przero,vzero/spsound,spsound,taukh,vzero

 call prompt('enter number of particles in x dir',npartx,8)
 deltax = dxbound/npartx

 npart = 0
 deltaxcloud = deltax*(denszero/denscloud)**(1./3.)
 write(*,*) 'psep in cloud = ',deltaxcloud,' psep/Rcl = ',deltax/rcloud

 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                  hfact,npart,xyzh,rmin=rcloud,nptot=npart_total)
 npartmed = npart
 totvol   = dxbound*dybound*dzbound - 4./3.*pi*rcloud**3
 totmass  = denszero*totvol
 massoftype = totmass/npart_total
 print*,' particle mass = ',massoftype(1),totmass,totvol

 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltaxcloud, &
                  hfact,npart,xyzh,rmax=rcloud,nptot=npart_total)
 print*,'-------------------------------------------------------------'
 print*,' number of particles in surrounding medium = ',npartmed
 print*,'               number of particles in blob = ',npart-npartmed
 print*,'                 total number of particles = ',npart
 print*,'-------------------------------------------------------------'

 npartoftype(:) = 0
 npartoftype(1) = npart

 !
 ! set thermal energy and velocities
 !
 do i=1,npart
    vxyzu(:,i) = 0.
    if ((xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2) < rcloud**2) then
       vxyzu(4,i) = przero/(gam1*denscloud)
    else
       vxyzu(4,i) = przero/(gam1*denszero)
       vxyzu(1,i) = vzero
    endif
 enddo
 !
 ! input file options
 !
 dtmax = 0.05*taukh
 tmax  = 2.5*taukh

end subroutine setpart

end module setup

