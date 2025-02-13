!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module allocutils
!
! Utilities related to memory allocation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dtypekdtree, io
!
 use io,           only:fatal,error,iprint,id,master,iverbose
 use dtypekdtree,  only:kdnode,kdnode_bytes

 implicit none

 public :: allocate_array
 real, public :: nbytes_allocated = 0.0
 public :: bytes2human

 interface allocate_array
  module procedure &
      allocate_array_real8_1d, &
      allocate_array_real8_2d, &
      allocate_array_real8_3d, &
      allocate_array_real8_4d, &
      allocate_array_real4_1d, &
      allocate_array_real4_2d, &
      allocate_array_real4_3d, &
      allocate_array_real4_4d, &
      allocate_array_integer8_1d, &
      allocate_array_integer4_1d, &
      allocate_array_integer4_1d_long, &
      allocate_array_integer4_2d, &
      allocate_array_integer4_3d, &
      allocate_array_integer1_1d, &
      allocate_array_integer1_2d, &
      allocate_array_integer1_3d, &
      allocate_array_kdnode_1d, &
      allocate_array_kdnode_1d_long, &
      allocate_array_logical
 end interface allocate_array

 private

contains

subroutine allocate_array_real8_1d(name, x, n1)
 character(*),                intent(in)     :: name
 real(kind=8), allocatable,   intent(inout)  :: x(:)
 integer,                     intent(in)     :: n1
 integer                                     :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1/),kind=8), 'real(8)')

end subroutine allocate_array_real8_1d

subroutine allocate_array_real8_2d(name, x, n1, n2)
 character(len=*),            intent(in)     :: name
 real(kind=8), allocatable,   intent(inout)  :: x(:,:)
 integer,                     intent(in)     :: n1, n2
 integer                                     :: allocstat

 allocate(x(n1, n2), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2/),kind=8), 'real(8)')

end subroutine allocate_array_real8_2d

subroutine allocate_array_real8_3d(name, x, n1, n2, n3)
 character(len=*),            intent(in)     :: name
 real(kind=8), allocatable,   intent(inout)  :: x(:, :, :)
 integer,                     intent(in)     :: n1, n2, n3
 integer                                     :: allocstat

 allocate(x(n1, n2, n3), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2, n3/),kind=8), 'real(8)')

end subroutine allocate_array_real8_3d

subroutine allocate_array_real8_4d(name, x, n1, n2, n3, n4)
 character(len=*),            intent(in)     :: name
 real(kind=8), allocatable,   intent(inout)  :: x(:,:,:,:)
 integer,                     intent(in)     :: n1, n2, n3, n4
 integer                                     :: allocstat

 allocate(x(n1, n2, n3, n4), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2, n3, n4/),kind=8), 'real(8)')

end subroutine allocate_array_real8_4d

subroutine allocate_array_real4_1d(name, x, n1)
 character(len=*),            intent(in)     :: name
 real(kind=4), allocatable,   intent(inout)  :: x(:)
 integer,                     intent(in)     :: n1
 integer                                     :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1/),kind=8), 'real(4)')

end subroutine allocate_array_real4_1d

subroutine allocate_array_real4_2d(name, x, n1, n2)
 character(len=*),            intent(in)     :: name
 real(kind=4), allocatable,   intent(inout)  :: x(:,:)
 integer,                     intent(in)     :: n1, n2
 integer                                     :: allocstat

 allocate(x(n1, n2), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2/),kind=8), 'real(4)')

end subroutine allocate_array_real4_2d

subroutine allocate_array_real4_3d(name, x, n1, n2, n3)
 character(len=*),            intent(in)     :: name
 real(kind=4), allocatable,   intent(inout)  :: x(:, :, :)
 integer,                     intent(in)     :: n1, n2, n3
 integer                                     :: allocstat

 allocate(x(n1, n2, n3), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2, n3/),kind=8), 'real(4)')

end subroutine allocate_array_real4_3d

subroutine allocate_array_real4_4d(name, x, n1, n2, n3, n4)
 character(len=*),            intent(in)     :: name
 real(kind=4), allocatable,   intent(inout)  :: x(:,:,:,:)
 integer,                     intent(in)     :: n1, n2, n3, n4
 integer                                     :: allocstat

 allocate(x(n1, n2, n3, n4), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2, n3, n4/),kind=8), 'real(4)')

end subroutine allocate_array_real4_4d

subroutine allocate_array_integer8_1d(name, x, n1)
 character(len=*),               intent(in)     :: name
 integer(kind=8), allocatable,   intent(inout)  :: x(:)
 integer,                        intent(in)     :: n1
 integer                                        :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1/),kind=8), 'integer(8)')
end subroutine allocate_array_integer8_1d

subroutine allocate_array_integer4_1d(name, x, n1)
 character(len=*),               intent(in)     :: name
 integer(kind=4), allocatable,   intent(inout)  :: x(:)
 integer,                        intent(in)     :: n1
 integer                                        :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1/),kind=8), 'integer(4)')

end subroutine allocate_array_integer4_1d

subroutine allocate_array_integer4_1d_long(name, x, n1)
 character(len=*),               intent(in)     :: name
 integer(kind=4), allocatable,   intent(inout)  :: x(:)
 integer(kind=8),                intent(in)     :: n1
 integer                                        :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, (/n1/), 'integer(4)')

end subroutine allocate_array_integer4_1d_long

subroutine allocate_array_integer4_2d(name, x, n1, n2)
 character(len=*),               intent(in)     :: name
 integer(kind=4), allocatable,   intent(inout)  :: x(:,:)
 integer,                        intent(in)     :: n1, n2
 integer                                        :: allocstat

 allocate(x(n1, n2), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2/),kind=8), 'integer(4)')

end subroutine allocate_array_integer4_2d

subroutine allocate_array_integer4_3d(name, x, n1, n2, n3)
 character(len=*),               intent(in)     :: name
 integer(kind=4), allocatable,   intent(inout)  :: x(:, :, :)
 integer,                        intent(in)     :: n1, n2, n3
 integer                                        :: allocstat

 allocate(x(n1, n2, n3), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2, n3/),kind=8), 'integer(4)')

end subroutine allocate_array_integer4_3d

subroutine allocate_array_integer1_1d(name, x, n1)
 character(len=*),               intent(in)     :: name
 integer(kind=1), allocatable,   intent(inout)  :: x(:)
 integer,                        intent(in)     :: n1
 integer                                        :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1/),kind=8), 'integer(1)')

end subroutine allocate_array_integer1_1d

subroutine allocate_array_integer1_2d(name, x, n1, n2)
 character(len=*),               intent(in)     :: name
 integer(kind=1), allocatable,   intent(inout)  :: x(:,:)
 integer,                        intent(in)     :: n1, n2
 integer                                        :: allocstat

 allocate(x(n1, n2), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2/),kind=8), 'integer(1)')

end subroutine allocate_array_integer1_2d

subroutine allocate_array_integer1_3d(name, x, n1, n2, n3)
 character(len=*),               intent(in)     :: name
 integer(kind=1), allocatable,   intent(inout)  :: x(:, :, :)
 integer,                        intent(in)     :: n1, n2, n3
 integer                                        :: allocstat

 allocate(x(n1, n2, n3), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1, n2, n3/),kind=8), 'integer(1)')

end subroutine allocate_array_integer1_3d

subroutine allocate_array_kdnode_1d(name, x, n1)
 character(len=*),               intent(in)     :: name
 type(kdnode), allocatable,      intent(inout)  :: x(:)
 integer,                        intent(in)     :: n1
 integer                                        :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1/),kind=8), 'kdnode')

end subroutine allocate_array_kdnode_1d

subroutine allocate_array_kdnode_1d_long(name, x, n1)
 character(len=*),               intent(in)     :: name
 type(kdnode), allocatable,      intent(inout)  :: x(:)
 integer(kind=8),                intent(in)     :: n1
 integer                                        :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, (/n1/), 'kdnode')

end subroutine allocate_array_kdnode_1d_long

subroutine allocate_array_logical(name, x, n1)
 character(len=*),          intent(in)     :: name
 logical, allocatable,      intent(inout)  :: x(:)
 integer,                   intent(in)     :: n1
 integer                                   :: allocstat

 allocate(x(n1), stat = allocstat)
 call check_allocate(name, allocstat)
 call print_allocation_stats(name, int((/n1/),kind=8), 'integer(4)')

end subroutine allocate_array_logical

subroutine check_allocate(name, allocstat)
 character(len=*),   intent(in) :: name
 integer,            intent(in) :: allocstat

 if (allocstat /= 0) call fatal('memory', name // ' allocation error')

end subroutine check_allocate

subroutine print_allocation_stats(name, xdim, type)
 character(len=*),   intent(in) :: name
 integer(kind=8),    intent(in) :: xdim(:)
 character(len=*),   intent(in) :: type
 character(len=10)              :: number
 character(len=14)              :: dimstring
 character(len=11)              :: sizestring
 integer                        :: i
 real                           :: nbytes
 integer                        :: databytes

 databytes = 0
 if (type == 'real(8)') then
    databytes = 8
 elseif (type == 'real(4)') then
    databytes = 4
 elseif (type == 'integer(8)') then
    databytes = 8
 elseif (type == 'integer(4)') then
    databytes = 4
 elseif (type == 'integer(1)') then
    databytes = 1
 elseif (type == 'kdnode') then
    databytes = kdnode_bytes
 else
    call fatal('memory', 'invalid data type chosen for memory allocation')
 endif

 nbytes = real(databytes)

 dimstring = '('
 do i = 1, size(xdim)
    ! Calculate size of array
    nbytes = nbytes * real(xdim(i))

    ! Make pretty string
    write(number, '(i0)') xdim(i)
    dimstring = trim(dimstring) // number
    if (i < size(xdim)) then
       dimstring = trim(dimstring) // ','
    endif
 enddo
 dimstring = trim(dimstring) // ')'

 nbytes_allocated = nbytes_allocated + nbytes

 if (id==master .and. nbytes > 0 .and. iverbose >= 2) then
    call bytes2human(nbytes, sizestring)
    write(iprint, '(a10, a22, a14, a11)') type, name, dimstring, sizestring
 endif

end subroutine print_allocation_stats

subroutine bytes2human(bytes, sizestring)
 real,                intent(in)  :: bytes
 character(len=11),   intent(out) :: sizestring

 if (bytes > 1073741824.0) then
    write(sizestring, '(f8.3, a3)') bytes / 1073741824.0, ' GB'
 elseif (bytes > 1048576.0) then
    write(sizestring, '(f8.3, a3)') bytes / 1048576.0, ' MB'
 elseif (bytes > 1024.0) then
    write(sizestring, '(f8.3, a3)') bytes / 1024.0, ' KB'
 else
    write(sizestring, '(f8.3, a3)') bytes, ' B '
 endif

end subroutine bytes2human

end module allocutils
