module memory
 implicit none

contains

 !
 !--Allocate all allocatable arrays: mostly part arrays, and tree structures
 !
 subroutine allocate_memory(n, part_only)
  use io, only:iprint,error,fatal
  use dim, only:update_max_sizes,maxp_hard
  use allocutils, only:nbytes_allocated,bytes2human
  use part, only:allocate_part
  use kdtree, only:allocate_kdtree
  use linklist, only:allocate_linklist

   integer,           intent(in) :: n
   logical, optional, intent(in) :: part_only

   character(len=11) :: sizestring

   if (n > maxp_hard) call fatal('memory', 'Trying to allocate above maxp hard limit.')

   call update_max_sizes(n)

   write(iprint, *)
   write(iprint, '(a)') '--> ALLOCATING ARRAYS'
   write(iprint, '(a)') '---------------------------------------------------------'

   if (nbytes_allocated > 0.0) then
      call error('part', 'Attempting to allocate memory, but memory is already allocated. &
      & Deallocating and then allocating again.')
      call deallocate_memory
   endif

   call allocate_part
   if (present(part_only)) then
     if (part_only) then
       call allocate_kdtree
       call allocate_linklist
     endif
   endif

   call bytes2human(nbytes_allocated, sizestring)
   write(iprint, '(a)') '---------------------------------------------------------'
   write(iprint, *) 'Total memory allocated to arrays: ', sizestring
   write(iprint, '(a)') '---------------------------------------------------------'

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
