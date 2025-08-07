!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for simple MHD wave propagation test
! as per section 5.1 of Iwasaki (2015)
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, io, kernel, options, part, physcon,
!   setup_params, slab, timestep
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup particles in uniform box with velocity perturbation
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB,npart_total
 use part,         only:Bxyz,mhd,igas,maxvxyzu
 use io,           only:master
 use timestep,     only:dtmax,tmax
 use options,      only:nfulldump
 use physcon,      only:pi
 use kernel,       only:hfact_default
 use infile_utils, only:get_options
 use slab,         only:set_slab,get_options_slab
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer                          :: i,ierr,nx
 real                             :: bzero,przero,uuzero,gam1,hzero,plasmabzero
!
!--general parameters
!
 time      = 0.
 tmax      = 0.6
 dtmax     = 0.06
 nfulldump = 1
 hfact     = hfact_default
!
!--setup parameters
!
 nx = 64
 plasmabzero = 3.
 rhozero     = 2.0
 przero      = 1.0
 bzero       = sqrt(2.0*przero/plasmabzero)
 if (maxvxyzu >= 4) then
    gamma = 5./3.
    gam1 = gamma - 1.
    uuzero = przero/(gam1*rhozero)
    polyk = przero/rhozero**gamma
 else
    gamma = 1.
    polyk = przero/rhozero
    uuzero = 1.5*polyk
 endif

 print "(/,a)",' Setup for MHD wave problem...'

 call get_options_slab(fileprefix,id,master,nx,rhozero,ierr,plasmab=plasmabzero)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 bzero  = sqrt(2.0*przero/plasmabzero)

 call set_slab(id,master,nx,-2.,2.,-1.,1.,hfact,npart,npart_total,xyzh,npartoftype,&
               rhozero,massoftype,igas)

 Bxyz = 0.0
 hzero = hfact*(massoftype(igas)/rhozero)**(1./3.)
 do i=1,npart
    vxyzu(1,i) = 0.01*exp(-(xyzh(1,i)/(3.0*hzero))**2)
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(1,i) = bzero
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup
