!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for simple MHD sine wave decay test
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, io, kernel, part, physcon, setup_params, slab
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  Sets up Bx = sin(2*pi*x), zero magnetic and velocity field otherwise
!
!  The amplitude should decay as exp(-eta * 4 * pi^2 * t)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB,npart_total
 use slab,         only:set_slab,get_options_slab
 use boundary,     only:xmin
 use part,         only:Bxyz,mhd,igas,maxvxyzu
 use io,           only:master
 use physcon,      only:pi
 use kernel,       only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: i,ierr,nx
 real    :: bzero,przero,uuzero,gam1
!
!--general parameters
!
 time = 0.
 hfact = hfact_default
!
!--setup parameters
!
 bzero = 1.0
 rhozero = 2.0
 przero = 10.0
 if (maxvxyzu >= 4) then
    gamma = 5./3.
    gam1 = gamma - 1.
    uuzero = przero/(gam1*rhozero)
    polyk = 0.
 else
    gamma = 1.
    polyk = przero/rhozero
    uuzero = 1.5*polyk
 endif

 print "(/,a)",' Setup for MHD sine problem...'

 nx = 64 ! default number of particles in x-direction
 call get_options_slab(trim(fileprefix),id,master,nx,rhozero,ierr)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 ! setup particles and boundaries for slab geometry
 call set_slab(id,master,nx,-0.5,0.5,-0.5,0.5,hfact,npart,npart_total,xyzh,&
               npartoftype,rhozero,massoftype,igas)

 do i=1,npart
    vxyzu(1,i) = 0.
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(2,i) = Bzero*sin(2.*pi*(xyzh(1,i)-xmin))
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup
