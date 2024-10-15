!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
 !
 ! Setup for simulations of HII Region expansion
 ! :References:
 !
 ! :Owner:
 !
 ! :Runtime parameters:
 !
 ! :Dependencies: datafiles, dim, eos, infile_utils, io, part, physcon,
 !   prompting, spherical, timestep, units
 !
 implicit none
 public :: setpart

 private

contains

 !----------------------------------------------------------------
 !+
 !  setup for galactic centre simulation (no gas)
 !+
 !----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,isetphase,iphase
 use units,     only:set_units,umass,unit_velocity,udist
 use physcon,   only:solarm,pc,pi,au,kboltz,mass_proton_cgs
 use io,        only:fatal,master
 use eos,       only:gmw,ieos
 use dim,       only: maxp,maxphase
 use timestep,  only:dtmax
 use spherical, only:set_sphere
 use datafiles, only:find_phantom_datafile
 use HIIRegion, only:iH2R
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: i,np,nx
 real    :: psep,rmax,rmin,totmass,totvol,rho0,temp
 !
 ! units (mass = solarm, length = 1 pc)
 !
 call set_units(dist=pc,mass=1.e4*solarm,G=1.d0)
 !
 ! general parameters
 !
 time = 0.
 hfact = 1.2
 gamma = 1.
 gmw = 1.  ! completely ionized, solar abu; eventually needs to be WR abu
 dtmax = 0.01
 rmin  = 0.
 rmax  = 6*pc/udist
 ieos  = 21
 IH2R  = 1
 temp = 1000.
 totmass  = 6.4e4*solarm/umass
 totvol   = 4./3.*pi*rmax**3
 rho0 = totmass/totvol
 polyk = ((gamma*kboltz*temp)/(gmw*mass_proton_cgs))*(1./unit_velocity)**2


 !
 ! space available for injected gas particles
 !
 npart = 0
 npartoftype(:) = 0

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 !
 ! Set up the black hole at the Galactic centre
 !
 nptmass = 1
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:)  = 0.
 xyzmh_ptmass(4,1) = .004        ! M=1 in code units by definition
 xyzmh_ptmass(ihacc,1)  = 0.0001 ! accretion radius
 xyzmh_ptmass(ihsoft,1) = 0.1    ! no softening

 !
 ! setup initial sphere of particles
 !
 np       = 8000000
 nx       = int(np**(1./3.))
 psep     = totvol**(1./3.)/real(nx)
 npart    = 0
 ! only set up particles on master, otherwise we will end up with n duplicates
 if (id==master) then
    call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart,xyzh,np_requested=np)
 endif
 np       = npart

!
!--set particle properties
!

 npartoftype(:) = 0
 npartoftype(igas) = npart
 massoftype(:)  = 0.0
 massoftype(igas)  = totmass/npartoftype(igas)
 if (maxphase==maxp) then
    do i=1,npart
       iphase(i) = isetphase(igas,iactive=.true.)
    enddo
 endif


 if (nptmass == 0) call fatal('setup','no particles setup')

end subroutine setpart

end module setup

