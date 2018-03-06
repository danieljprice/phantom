!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!   Moves two distant spheres into a binary orbit
!   Author: Bernard Field (supervisor: James Wurster)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, extern_gwinspiral, externalforces, io,
!    options, part, physcon, prompting, timestep, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

 logical, parameter :: use_defaults = .false.  ! if .true. will automatically use default values
 ! if .false., will ask user to prompt new values
 !--The default values
 real,    private   :: separation   = 50.
 logical, private   :: usegw        = .true.

contains
!-----------------------------------------------------------------------

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use io,             only: iprint,fatal
 use prompting,      only: prompt
 use options,        only: iexternalforce,nfulldump,damp
 use part,           only: igas
 use units,          only: unit_velocity
 use physcon,        only: c,pi
 use timestep,       only: tmax, dtmax
 use centreofmass,   only: get_centreofmass
 use externalforces, only: iext_gwinspiral
 use extern_gwinspiral, only: Nstar
 integer, intent(inout)    :: npart
 integer, intent(inout)    :: npartoftype(:)
 real,    intent(inout)    :: massoftype(:)
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i
 real                      :: com(3),com_star1(3),com_star2(3),vcom(3)
 real                      :: rad1,rad2,speed1,speed2,mstar1,mstar2,mtotal
 real                      :: c_code,tmax0

 !
 !--From the .setup file, get the number of particles per star
 if (Nstar(1) <= 0 .or. Nstar(2) <= 0) call fatal('moddump','Require particle numbers in both stars')
 !
 !--Request parameters (unless hardcoded to use defaults)
 ! now determine the parameters of their new orbit
 if (.not.use_defaults) then
    call prompt('Enter desired separation:',separation,0.)
    call prompt('Use Gravitational Wave inspiral?',usegw,.false.)
 endif
 !
 !--Locate mass centre of mass and mass of stars (assuming no particles added or lost)
 call get_centreofmass(com,vcom,npart,xyzh,vxyzu)
 call get_centreofmass(com_star1,vcom,nstar(1),xyzh(:,1:nstar(1)),vxyzu(:,1:nstar(1)))
 call get_centreofmass(com_star2,vcom,nstar(2),xyzh(:,nstar(1)+1:npart),vxyzu(:,nstar(1)+1:npart))
 mstar1 = nstar(1) * massoftype(igas)
 mstar2 = nstar(2) * massoftype(igas)
 mtotal = npart    * massoftype(igas)
 !
 !--Calcuate the new orbital parameters
 rad1   =  separation * mstar2/mtotal           ! distance of star 1 from the CoM
 rad2   =  separation * mstar1/mtotal           ! distance of star 1 from the CoM
 speed1 = -rad1 * sqrt(mtotal/separation**3)    ! speed of star 1
 speed2 =  rad2 * sqrt(mtotal/separation**3)    ! speed of star 2
 !
 !--Place stars on new orbits
 !  for simplicity, assume stars are on the x-axis
 vxyzu(1:3,:) = 0.0                             ! reset velocity
 xyzh (1:3,:) = 0.0                             ! reset location
 do i=1,nstar(1)
    xyzh(1,i)  =  rad1
    vxyzu(2,i) =  speed1
 enddo
 do i=nstar(1)+1,npart
    xyzh(1,i)  = -rad2
    vxyzu(2,i) =  speed2
 enddo
 !
 !--Set new runtime parameters
 tmax           = 1000.
 dtmax          =  100.
 damp           =    0.
 nfulldump      =   10
 iexternalforce =    0
 !
 !--Modify time and external forces, if including gravitational waves
 if (usegw) then
    iexternalforce = iext_gwinspiral
    c_code         = c/unit_velocity
    tmax0          = 5./256.*c_code**5*separation**4/(mstar1*mstar2*mtotal)
    write(iprint,'(1x,a,f16.8)') 'moddump_nsmerger: estimated time to collision = ', tmax0
 else
    tmax0 = 0.0
 endif
 tmax = max(tmax,2.0*tmax0)
 !
 write(iprint,'(1x,a)')       'moddump_nsmerger: This will move two distant objects into a binary system'
 write(iprint,'(1x,a,f16.8)') 'moddump_nsmerger: period of orbit = ', 2.*pi*sqrt(separation**3/mtotal)
 write(iprint,'(1x,a,f16.8)') 'moddump_nsmerger: distance of star 1 from origin = ', rad1
 write(iprint,'(1x,a,f16.8)') 'moddump_nsmerger: distance of star 2 from origin = ', rad2
 write(iprint,'(1x,a,f16.8)') 'moddump_nsmerger: speed of star 1 = ', speed1
 write(iprint,'(1x,a,f16.8)') 'moddump_nsmerger: speed of star 2 = ', speed2
 !
 return
end subroutine modify_dump
!-----------------------------------------------------------------------
end module moddump
