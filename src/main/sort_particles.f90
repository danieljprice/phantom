!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module sort_particles
!
! sorts the particles so neighbours are also close in memory
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, part, sortutils
!
 implicit none
 public :: sort_part_radius, sort_part_id

 private

contains

!----------------------------------------------------------------
!+
!  this version sorts the particles by radius (e.g. for a disc)
!+
!----------------------------------------------------------------
subroutine sort_part_radius(np)
 use io,        only:iprint,error
 use part,      only:xyzh,reorder_particles,npart,ll
 use sortutils, only:indexxfunc,r2func
 integer, intent(in) :: np
 real :: t1,t2

 call cpu_time(t1)
 write(iprint,*) '> sorting particles in radius...'
 if (np /= npart) call error('sort','np /= npart')

 call indexxfunc(npart,r2func,xyzh,ll)

 write(iprint,*) ' copying arrays...'
!
!--copy arrays into correct order
!
 call reorder_particles(ll,npart)

 call cpu_time(t2)
 write(iprint,*) '> sort completed in ',t2-t1,'s'

end subroutine sort_part_radius

!----------------------------------------------------------------
!+
!  this version sorts the particles by ID
!+
!----------------------------------------------------------------
subroutine sort_part_id
 use io,        only:iprint,error
 use part,      only:reorder_particles,npart,ll,iorig
 use sortutils, only:indexx
 real :: t1,t2

 call cpu_time(t1)
 write(iprint,*) '> sorting particles by ID...'

 call indexx(npart,iorig,ll)

 write(iprint,*) ' copying arrays...'
 !
 !--copy arrays into correct order
 !
 call reorder_particles(ll,npart)

 call cpu_time(t2)
 write(iprint,*) '> sort completed in ',t2-t1,'s'

end subroutine sort_part_id

end module sort_particles
