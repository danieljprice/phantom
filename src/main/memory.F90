!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: memory
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: allocutils, dim, io, kdtree, linklist, part, photoevap
!+
!--------------------------------------------------------------------------
module memory
 implicit none

contains

 !
 !--Allocate all allocatable arrays: mostly part arrays, and tree structures
 !
subroutine allocate_memory(n, part_only)
 use io, only:iprint,error,fatal,nprocs,id,master
 use dim, only:update_max_sizes,maxp
 use allocutils, only:nbytes_allocated,bytes2human
 use part, only:allocate_part
 use kdtree, only:allocate_kdtree
 use linklist, only:allocate_linklist
#ifdef PHOTO
 use photoevap, only:allocate_photoevap
#endif

 integer,           intent(in) :: n
 logical, optional, intent(in) :: part_only

 logical :: part_only_
 character(len=11) :: sizestring

 if (present(part_only)) then
    part_only_ = part_only
 else
    part_only_ = .false.
 endif

 if (nbytes_allocated > 0.0 .and. n <= maxp) then
    !print "(a)",' ARRAYS ALREADY ALLOCATED... SKIPPING'
    return ! just silently skip if arrays are already large enough
 endif
 call update_max_sizes(n)

 if (nprocs == 1) then
    write(iprint, *)
    if (part_only_) then
       write(iprint, '(a)') '--> ALLOCATING PART ARRAYS'
    else
       write(iprint, '(a)') '--> ALLOCATING ALL ARRAYS'
    endif
    write(iprint, '(a)') '---------------------------------------------------------'
 endif

 if (nbytes_allocated > 0.0) then
    call error('memory', 'Attempting to allocate memory, but memory is already allocated. &
    & Deallocating and then allocating again.')
    call deallocate_memory(part_only=part_only_)
    call update_max_sizes(n)
 endif

 call allocate_part
 if (.not. part_only_) then
    call allocate_kdtree
    call allocate_linklist
#ifdef PHOTO
    call allocate_photoevap
#endif
 endif

 call bytes2human(nbytes_allocated, sizestring)
 if (nprocs == 1) then
    write(iprint, '(a)') '---------------------------------------------------------'
    write(iprint, *) 'Total memory allocated to arrays: ', sizestring
    write(iprint, '(a)') '---------------------------------------------------------'
 else
    write(iprint, *) id, 'allocated ', sizestring
 endif

end subroutine allocate_memory

subroutine deallocate_memory(part_only)
 use dim, only:update_max_sizes
 use part, only:deallocate_part
 use kdtree, only:deallocate_kdtree
 use linklist, only:deallocate_linklist
#ifdef PHOTO
 use photoevap, only:deallocate_photoevap
#endif
 use allocutils, only:nbytes_allocated

 logical, optional, intent(in) :: part_only
 logical :: part_only_

 if (present(part_only)) then
    part_only_ = part_only
 else
    part_only_ = .false.
 endif

 call deallocate_part
 if (.not. part_only_) then
    call deallocate_kdtree
    call deallocate_linklist
#ifdef PHOTO
    call deallocate_photoevap
#endif
 endif

 nbytes_allocated = 0
 call update_max_sizes(0)

end subroutine deallocate_memory

end module memory
