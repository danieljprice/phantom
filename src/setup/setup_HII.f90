!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! setup
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: HIIRegion, datafiles, dim, eos, io, part, physcon,
!   prompting, random, spherical, timestep, units, utils_shuffleparticles
!
 implicit none
 public :: setpart

 private
 integer :: nsrc = 1
 logical :: add_blob   = .false. ! dense blob next to the centre, to test the HII sectors
 real    :: rho_blob   = 100.    ! density contrast of the blob
 real    :: r_blob     = 0.1     ! radius of the blob [pc]
 real    :: d_blob     = 0.3     ! distance of the blob centre from the origin, along x [pc]

contains

 !----------------------------------------------------------------
 !+
 !  setup for testing HII region expansion in a uniform cloud
 !+
 !----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,isetphase,iphase
 use units,     only:set_units,umass,unit_velocity,udist
 use physcon,   only:solarm,pc,pi,au,kboltz,mass_proton_cgs
 use io,        only:fatal,master
 use eos,       only:gmw,ieos
 use dim,       only: maxp,maxphase,maxptmass,maxvxyzu
 use random,    only: ran2
 use timestep,  only:dtmax
 use spherical, only:set_sphere
 use datafiles, only:find_phantom_datafile
 use HIIRegion, only:iH2R
 use utils_shuffleparticles, only:shuffleparticles

 integer,          intent(in)    :: id
 integer,          intent(inout) :: npart
 integer,          intent(out)   :: npartoftype(:)
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: massoftype(:)
 real,             intent(out)   :: polyk,gamma,hfact
 real,             intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 real,             intent(out)   :: vxyzu(:,:)
 integer(kind=8) :: nptot
 integer :: i,np,nx,seed,nblob,nkeep
 real    :: psep,rmax,rmin,totmass,totvol,rho0,temp,pmass,rb,xblob(3)
 !
 ! units (mass = solarm, length = 1 pc)
 !
 call set_units(dist=pc,mass=1.e4*solarm,G=1.d0)
 !
 ! general parameters
 !
 time = 0.
 hfact = 1.2
 gamma = 5./3.
 gmw = 1.  ! completely ionized, solar abu; eventually needs to be WR abu
 dtmax = 0.01
 rmin  = 0.
 rmax  = 1.257*pc/udist
 ieos  = 21
 IH2R  = 1
 temp = 100.
 totmass  = 6.4e2*solarm/umass
 totvol   = 4./3.*pi*rmax**3
 rho0 = totmass/totvol
 !--isothermal c_s^2 = kT/mu of the cold gas, as the HII eos expects (Tcold)
 polyk = ((kboltz*temp)/(gmw*mass_proton_cgs))*(1./unit_velocity)**2

 npart = 0
 npartoftype(:) = 0

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.

 if (maxvxyzu >= 4) vxyzu(4,:) = polyk/(gamma-1.)

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
 np       = 1000000
 nx       = int(np**(1./3.))
 psep     = totvol**(1./3.)/real(nx)
 npart    = 0
 nptot    = 0
 ! only set up particles on master, otherwise we will end up with n duplicates
 if (id==master) then
    call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart,xyzh,nptot,np_requested=np)
 endif
 np       = npart
 pmass    = totmass/real(np)

 !
 ! dense blob: remove the cloud particles inside it and refill it with
 ! particles of the same mass at rho_blob times the cloud density
 !
 if (add_blob .and. id==master) then
    rb    = r_blob*pc/udist
    xblob = (/d_blob*pc/udist,0.,0./)
    nkeep = 0
    do i=1,npart
       if (sum((xyzh(1:3,i) - xblob)**2) >= rb**2) then
          nkeep = nkeep + 1
          xyzh(:,nkeep) = xyzh(:,i)
       endif
    enddo
    npart = nkeep
    nptot = npart
    nblob = nint(rho_blob*rho0*4./3.*pi*rb**3/pmass)
    call set_sphere('random',id,master,0.,rb,psep,hfact,npart,xyzh,nptot,&
                    np_requested=nblob,xyz_origin=xblob)
    write(*,"(1x,a,i0,a,i0,a)") 'blob: ',npart-nkeep,' particles, ',np-nkeep,' cloud particles removed'
 endif

!
!--set particle properties
!

 npartoftype(:) = 0
 npartoftype(igas) = npart
 massoftype(:)  = 0.0
 massoftype(igas)  = pmass
 if (maxphase==maxp) then
    do i=1,npart
       iphase(i) = isetphase(igas,iactive=.true.)
    enddo
 endif

 !call shuffleparticles(iprint,npart,xyzh,massoftype(1),rsphere=rmax,dsphere=rho0,dmedium=0.,&
 !is_setup=.true.,prefix=trim(fileprefix))

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
 call prompt('Add a dense blob next to the centre?',add_blob)
 if (add_blob) then
    call prompt('Enter the density contrast of the blob',rho_blob,1.)
    call prompt('Enter the radius of the blob (pc)',r_blob,0.)
    call prompt('Enter the distance of the blob from the centre along x (pc)',d_blob,0.)
 endif

end subroutine get_input_from_prompts

end module setup

