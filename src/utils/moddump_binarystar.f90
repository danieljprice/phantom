!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  Input is a relaxed star, output is two relaxed stars in binary orbit
!
!  REFERENCES: None
!
!  OWNER: Terrence Tricco
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, part, prompting, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type
 use units,        only: set_units,udist,umass,utime
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real :: rad

 rad = 10.0
 print *, 'Distance unit is: ', udist
 call prompt('Enter radial separation between stars (in multiple of distance unit)', rad, 0.)


 npart = npartoftype(igas)

 ! make any adjustments to first star here
 do i = 1, npart

 enddo


 ! duplicate relaxed star
 do i = npart+1, 2*npart
     ! place star a distance rad away
     xyzh(1,i) = xyzh(1,i-npart) + rad
     xyzh(2,i) = xyzh(2,i-npart) + rad
     xyzh(3,i) = xyzh(3,i-npart)
     xyzh(4,i) = xyzh(4,i-npart)

     vxyzu(1,i) = vxyzu(1,i-npart)
     vxyzu(2,i) = vxyzu(2,i-npart)
     vxyzu(3,i) = vxyzu(3,i-npart)

     call set_particle_type(i,igas)
 enddo

 npart = 2 * npart
 npartoftype(igas) = npart


 ! reset centre of mass of the binary system
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

end subroutine modify_dump

end module moddump

