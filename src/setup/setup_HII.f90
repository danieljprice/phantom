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
 integer :: nsrc = 1

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
 use io,        only:fatal,master,iprint
 use eos,       only:gmw,ieos
 use dim,       only: maxp,maxphase,maxptmass
 use random,    only: ran2
 use timestep,  only:dtmax
 use spherical, only:set_sphere
 use datafiles, only:find_phantom_datafile
 use HIIRegion, only:iH2R
 use utils_shuffleparticles, only:shuffleparticles

 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer(kind=8) :: nptot
 integer :: i,np,nx,seed
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
 rmax  = 2.91*pc/udist
 ieos  = 21
 IH2R  = 2
 temp = 1000.
 totmass  = 8.e3*solarm/umass
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
 nptmass = nsrc
 if (id==master) then
    call get_input_from_prompts()
 endif

 nptmass = nsrc

 if (nptmass > 1) then
    seed = 12
    do i=1,nptmass
       xyzmh_ptmass(1,i) = 0.5-ran2(seed)
       xyzmh_ptmass(2,i) = 0.5-ran2(seed)
       xyzmh_ptmass(3,i) = 0.5-ran2(seed)
       vxyz_ptmass(:,i)  = 0.
       xyzmh_ptmass(4,i) = .002*ran2(seed)+0.002
    enddo
    xyzmh_ptmass(ihacc,:)  = 0.0001 ! accretion radius
    xyzmh_ptmass(ihsoft,:) = 0.1    ! no softening
 else
    xyzmh_ptmass(1,1) = 0.
    xyzmh_ptmass(2,1) = 0.
    xyzmh_ptmass(3,1) = 0.
    vxyz_ptmass(:,1)  = 0.
    xyzmh_ptmass(4,1) = .004
    xyzmh_ptmass(ihacc,1)  = 0.0001
    xyzmh_ptmass(ihsoft,1) = 0.1
 endif

 !
 ! setup initial sphere of particles
 !
 np       = 2500000
 nx       = int(np**(1./3.))
 psep     = totvol**(1./3.)/real(nx)
 npart    = 0
 ! only set up particles on master, otherwise we will end up with n duplicates
 if (id==master) then
    call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart,xyzh,nptot,np_requested=np)
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

 call shuffleparticles(iprint,npart,xyzh,massoftype(1),rsphere=rmax,dsphere=rho0,dmedium=0.,&
 is_setup=.true.,prefix=trim(fileprefix))


 if (nptmass == 0) call fatal('setup','no particles setup')

end subroutine setpart


!----------------------------------------------------------------
!
!  Prompt user for inputs
!
!----------------------------------------------------------------
subroutine get_input_from_prompts()
 use prompting, only:prompt
 use dim,       only:maxptmass

 call prompt('Enter the number of ionizing sources in the sphere',nsrc,1,maxptmass)

end subroutine get_input_from_prompts

end module setup

