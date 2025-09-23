!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for Taylor-Green Vortex test (used in Phantom paper)
!
! :References: Price et al. (2017)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, kernel, part, physcon, setup_params, slab
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for the Taylor-Green problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master
 use slab,         only:set_slab,get_options_slab
 use physcon,      only:pi
 use part,         only:igas,maxvxyzu
 use kernel,       only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: vzero
 integer :: i,nx,ierr
!
!--general parameters
!
 time = 0.
 hfact = hfact_default
 gamma = 1.
!
!--setup particles
!
 nx = 128
 rhozero = 1.0
 call get_options_slab(fileprefix,id,master,nx,rhozero,ierr)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 polyk = 0.
 if (maxvxyzu < 4) polyk = 1.
 npart = 0; npartoftype(:) = 0; npart_total = 0

 call set_slab(id,master,nx,0.,1.,0.,1.,hfact,npart,npart_total,xyzh,&
               npartoftype,rhozero,massoftype,igas)

 vzero = 0.1

 do i=1,npart
    vxyzu(1:3,i) = 0.
    !--velocity field for the Taylor-Green problem
    vxyzu(1,i) =  vzero*sin(2.*pi*xyzh(1,i))*cos(2.*pi*xyzh(2,i))
    vxyzu(2,i) = -vzero*cos(2.*pi*xyzh(1,i))*sin(2.*pi*xyzh(2,i))
 enddo

end subroutine setpart

end module setup
