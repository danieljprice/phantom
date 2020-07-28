!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: does not make any modifications
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, splitpart
!
 implicit none
 integer, parameter :: nchildren = 8

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart, only:split_particles
 use io,        only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr

 call split_particles(nchildren,npart,npartoftype,xyzh,massoftype,ierr)
 if (ierr /= 0) call fatal('moddump','could not split particles')
 print*,' got npart = ',npart

end subroutine modify_dump

end module moddump
