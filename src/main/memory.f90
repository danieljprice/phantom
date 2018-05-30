module memory
 implicit none

contains

 !
 !--Allocate all allocatable arrays: mostly part arrays, and tree structures
 !
 subroutine allocate_memory(n)
  use io, only:iprint,error
  use dim, only:update_max_sizes
  use allocutils, only:nbytes_allocated,bytes2human
  use part, only:allocate_part
  use kdtree, only:allocate_kdtree
  use linklist, only:allocate_linklist

   integer, intent(in) :: n

   character(len=11) :: sizestring

   write(iprint, *)
   write(iprint, '(a)') '--> ALLOCATING ARRAYS'
   write(iprint, '(a)') '--------------------------------------------------------'

   if (nbytes_allocated > 0.0) then
      call error('part', 'Attempting to allocate memory, but memory is already allocated. &
      & Deallocating and then allocating again.')
      call deallocate_memory
   endif

   call update_max_sizes(n)
   call allocate_part
   call allocate_kdtree
   call allocate_linklist

   call bytes2human(nbytes_allocated, sizestring)
   write(iprint, '(a)') '--------------------------------------------------------'
   write(iprint, *) 'Total memory allocated to arrays: ', sizestring
   write(iprint, '(a)') '--------------------------------------------------------'

 end subroutine allocate_memory

 subroutine deallocate_memory
    use part, only:deallocate_part
    use kdtree, only:deallocate_kdtree
    use linklist, only:deallocate_linklist

    call deallocate_part
    call deallocate_kdtree
    call deallocate_linklist
 end subroutine deallocate_memory

end module memory
