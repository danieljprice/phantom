!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  merges particles; input simulation with npart, get npart/nchild back
!
!  REFERENCES: Vacondio et al. 2013
!
!  OWNER: Rebecca Nealon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, splitpart
!+
!--------------------------------------------------------------------------
module moddump
 implicit none
 integer, parameter :: nchild = 10

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart, only:merge_all_particles
 use part,      only:igas,kill_particle,delete_dead_or_accreted_particles
 use part,      only:isdead_or_accreted,copy_particle
 use io,        only:fatal,error
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,nactive

 !-- how many active particles
 nactive = 0
 do i = 1,npart
   if (.not.isdead_or_accreted(xyzh(4,i))) then
     nactive = nactive + 1
   else
     call kill_particle(i,npartoftype)
   endif
 enddo

 if (nactive < npart) then
   call delete_dead_or_accreted_particles(npart,npartoftype)
   print*,' discarding inactive particles'
 endif

 ! Merge 'em!
 call merge_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                                nchild,nactive)

 print*,' new npart = ',npart

end subroutine modify_dump

end module moddump
