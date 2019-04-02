!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup for tokamak torus with equilibrium density profiles
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, extern_Bfield, externalforces, geometry, io, kernel,
!    mpiutils, options, part, physcon, random, setup_params, stretchmap
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for tokamak torus, originally by D. Price & C. Toniolo
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:periodic,maxvxyzu,maxp
 use setup_params, only:rhozero
 use io,           only:master,fatal
 use mpiutils,     only:bcast_mpi
 use physcon,      only:pi
 use random,       only:ran2
 use kernel,       only:hfact_default
 use options,      only:iexternalforce
 use part,         only:igas
 use extern_Bfield, only:Rtorus,a_on_Rtorus,currJ0
 use stretchmap,    only:set_density_profile
 use geometry,      only:coord_transform
 use externalforces, only:iext_externB
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: deltax,deltaz,deltaphi,randphi
 real :: xtorus,ztorus,atorus,rintorus2
 real :: phi,rcyl,ra2,pri,da2,gam1,hzero,totmass,totvol
 integer :: ipart,i,j,k,nx,nphi,nz,iseed

 if (periodic) then
    call fatal('setup','periodic boundaries not compatible with torus setup')
 endif
 !if (maxvxyzu < 4) then
 !   call fatal('setup','need ISOTHERMAL=no for this setup')
 !endif
!
!--general parameters
!
 time = 0.
 hfact = hfact_default
 gamma = 5./3.
 iseed = 1
!
!--setup particles
!
 atorus = Rtorus*a_on_Rtorus
 nx = 15
 nz = 15
 deltaz = 2.*atorus/nz
 deltax = 2.*atorus/nx
!
! particle and eos properties/constants
 rhozero = 1.
 hzero = 0.6*deltax
 da2 = 1./atorus**2
 ra2 = 0.
 pri = currJ0*currJ0*atorus*atorus* &
       (47. - 12.*ra2**5 + 75.*ra2**4 - 200.*ra2**3 + 270.*ra2**2 - 180.*ra2)/720.
 polyk = pri/rhozero**gamma
 gam1  = gamma - 1.

 ipart = 0
 do k=1,nz
    ztorus = (k-0.5)*deltaz - atorus
    do j=1,nx
       xtorus = (j-0.5)*deltax - atorus
       rintorus2 = xtorus**2 + ztorus**2

       if (rintorus2 < atorus**2) then
          rcyl     = xtorus + Rtorus
          deltaphi = deltaz*Rtorus/rcyl
          nphi     = int(2.*pi/deltaphi)
          deltaphi = 2.*pi/nphi
          randphi  = ran2(iseed)*deltaphi
          !--make ring of particles at r=rcyl
          do i=1,nphi
             ipart = ipart + 1
             if (ipart > maxp) stop 'setup_tokamak: ipart>maxp; recompile with MAXP=big number'
             phi = (i-1)*deltaphi
             xyzh(1,ipart) = rcyl*cos(phi+randphi)
             xyzh(2,ipart) = rcyl*sin(phi+randphi)
             xyzh(3,ipart) = ztorus
             xyzh(4,ipart) = 2.*hzero
             ra2 = rintorus2*da2
             pri = (currJ0*atorus)**2*(47.-12.*ra2**5+ &
                  75.*ra2**4-200.*ra2**3+270.*ra2**2-180.*ra2)/720.
             ! zero velocities
             vxyzu(1:3,i) = 0.
             if (maxvxyzu >= 4) vxyzu(4,ipart) = pri/gam1/rhozero
          enddo
       endif
    enddo
 enddo

 npart = ipart
 npartoftype(igas) = ipart

! choose an initial density
 totvol = (pi*atorus**2)*2.*pi*rtorus
 totmass = totvol*rhozero
 massoftype(igas) = totmass/npart

!
!--Stretching the spatial distribution to have the desired density distribution
!
 call set_density_profile(npart,xyzh,rhofunc=densfunc,min=0.,max=atorus,geom=4)
!
!--set input options
!
 iexternalforce = iext_externB

contains

real function densfunc(r)
 real, intent(in) :: r
 real :: ra2, pri

 !densfunc = 1. - (r/atorus)**2
 ra2 = (r/atorus)**2
 pri = (currJ0*atorus)**2*(47.-12.*ra2**5+ &
                  75.*ra2**4-200.*ra2**3+270.*ra2**2-180.*ra2)/720.
 densfunc = (pri/polyk)**(1./gamma)

end function densfunc

end subroutine setpart

end module setup

