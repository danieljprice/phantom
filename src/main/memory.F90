module memory
   use io,  only:fatal
 implicit none

 public :: allocate_memory

contains

subroutine allocate_memory
 use dim,   only:maxp
 use part,  only:xyzh

 integer :: allocstat

 allocate(xyzh(4, maxp), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','xyzh allocation error')

end subroutine allocate_memory

end module memory
