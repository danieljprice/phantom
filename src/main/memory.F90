module memory
   use io,  only:fatal,iprint
 implicit none

 public :: allocate_memory

 real :: nbytes_allocated

 interface allocate_array
    module procedure &
      allocate_array_real8_1d, &
      allocate_array_real8_2d, &
      allocate_array_real4_1d, &
      allocate_array_real4_2d, &
      allocate_array_integer4_1d, &
      allocate_array_integer4_2d, &
      allocate_array_integer1_1d, &
      allocate_array_integer1_2d
 end interface

contains

 subroutine allocate_array_real8_1d(name, x, n1)
    character(*),                intent(in)     :: name
    real(kind=8), allocatable,   intent(inout)  :: x(:)
    integer,                     intent(in)     :: n1
    integer                                     :: allocstat

    allocate(x(n1), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1/), 'real(4)')
 end subroutine allocate_array_real8_1d

 subroutine allocate_array_real8_2d(name, x, n1, n2)
    character(len=*),            intent(in)     :: name
    real(kind=8), allocatable,   intent(inout)  :: x(:, :)
    integer,                     intent(in)     :: n1, n2
    integer                                     :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'real(8)')
 end subroutine allocate_array_real8_2d

 subroutine allocate_array_real4_1d(name, x, n1)
    character(len=*),            intent(in)     :: name
    real(kind=4), allocatable,   intent(inout)  :: x(:)
    integer,                     intent(in)     :: n1
    integer                                     :: allocstat

    allocate(x(n1), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1/), 'real(4)')
 end subroutine allocate_array_real4_1d

 subroutine allocate_array_real4_2d(name, x, n1, n2)
    character(len=*),            intent(in)     :: name
    real(kind=4), allocatable,   intent(inout)  :: x(:, :)
    integer,                     intent(in)     :: n1, n2
    integer                                     :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'real(4)')
 end subroutine allocate_array_real4_2d

 subroutine allocate_array_integer4_1d(name, x, n1)
    character(len=*),               intent(in)     :: name
    integer(kind=4), allocatable,   intent(inout)  :: x(:)
    integer,                        intent(in)     :: n1
    integer                                        :: allocstat

    allocate(x(n1), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1/), 'integer(4)')
 end subroutine allocate_array_integer4_1d

 subroutine allocate_array_integer4_2d(name, x, n1, n2)
    character(len=*),               intent(in)     :: name
    integer(kind=4), allocatable,   intent(inout)  :: x(:, :)
    integer,                        intent(in)     :: n1, n2
    integer                                        :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'integer(4)')
 end subroutine allocate_array_integer4_2d

 subroutine allocate_array_integer1_1d(name, x, n1)
    character(len=*),               intent(in)     :: name
    integer(kind=1), allocatable,   intent(inout)  :: x(:)
    integer,                        intent(in)     :: n1
    integer                                        :: allocstat

    allocate(x(n1), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1/), 'integer(1)')
 end subroutine allocate_array_integer1_1d

 subroutine allocate_array_integer1_2d(name, x, n1, n2)
    character(len=*),               intent(in)     :: name
    integer(kind=1), allocatable,   intent(inout)  :: x(:, :)
    integer,                        intent(in)     :: n1, n2
    integer                                        :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'integer(1)')
 end subroutine allocate_array_integer1_2d

 subroutine check_allocate(name, allocstat)
    character(len=*),   intent(in) :: name
    integer,            intent(in) :: allocstat

    if (allocstat /= 0) call fatal('memory', name // ' allocation error')
 end subroutine check_allocate

 subroutine print_allocation_stats(name, xdim, type)
    character(len=*),   intent(in) :: name
    integer,            intent(in) :: xdim(:)
    character(len=*),   intent(in) :: type
    character(len=10)              :: number
    character(len=14)              :: dimstring
    character(len=10)               :: sizestring
    integer                        :: i
    real                           :: nbytes
    integer                        :: databytes

    if (type == 'real(8)') then
       databytes = 8
    else if (type == 'real(4)') then
       databytes = 4
    else if (type == 'integer(4)') then
       databytes = 4
    else if (type == 'integer(1)') then
       databytes = 1
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
          dimstring = trim(dimstring) // ':'
       endif
    enddo
    dimstring = trim(dimstring) // ')'

    nbytes_allocated = nbytes_allocated + nbytes

    call bytes2human(nbytes, sizestring)

    write(iprint, '(a10, a22, a14, a10)') type, name, dimstring, sizestring
 end subroutine print_allocation_stats

subroutine bytes2human(bytes, sizestring)
   real,                intent(in)  :: bytes
   character(len=10),   intent(out) :: sizestring

   if (bytes > 1073741824.0) then
      write(sizestring, '(f7.3, a3)') bytes / 1073741824.0, ' GB'
   else if (bytes > 1048576.0) then
      write(sizestring, '(f7.3, a3)') bytes / 1048576.0, ' MB'
   else if (bytes > 1024.0) then
      write(sizestring, '(f7.3, a3)') bytes / 1024.0, ' KB'
   else
      write(sizestring, '(f7.3, a3)') bytes, ' B '
   endif
end subroutine bytes2human


subroutine allocate_memory
 use dim,   only:maxp

 use part,  only:xyzh
 use part,  only:xyzh_soa

 use dim,   only:maxvxyzu
 use part,  only:vxyzu

 use dim,   only:nalpha, maxalpha
 use part,  only:alphaind

 use dim,   only:ndivcurlv
 use part,  only:divcurlv

 use dim,   only:ndivcurlB
 use part,  only:divcurlB

 use dim,   only:maxBevol,maxmhd
 use part,  only:Bevol

 use part,  only:Bxyz

 use dim,   only:maxp_growth
 use part,  only:dustprop

 use part,  only:St

 use dim,   only:maxstrain
 use part,  only:straintensor

 use dim,   only:nabundances,maxp_h2
 use part,  only:abundance

 use dim,   only:maxtemp
 use part,  only:temperature

 use dim,   only:maxp_dustfrac
 use part,  only:dustfrac,dustevol,deltav

 use dim,   only:nsinkproperties,maxptmass
 use part,  only:xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink

 use dim,   only:maxgrav
 use part,  only:poten

 use dim,   only:maxmhdni
 use part,  only:n_R

 use dim,   only:maxne
 use part,  only:n_electronT

 use part,  only:eta_nimhd

 use dim,   only:maxlum
 use part,  only:luminosity

 use part,  only:fxyzu

 use part,  only:dBevol

 use part, only:divBsymm

 use part, only:fext

 use part, only:ddustfrac

 use part, only:ddustprop

 use part, only:vpred

 use part, only:dustpred

 use part, only:Bpred

 use part, only:dustproppred

#ifdef IND_TIMESTEPS
 use part, only:ibin

 use part, only:ibin_old

 use part, only:ibin_wake

 use part, only:dt_in

 use part, only:twas
#endif

 use part, only:iphase, iphase_soa

 use dim,  only:ngradh
 use part, only:gradh

 use part, only:tstop

 use part, only:ll

 !
 !--for analysis routines, do not allocate any more storage
 !  than is strictly necessary. This will eventually be deprecated by the
 !  memory manager.
 !
#ifdef ANALYSIS
 integer, parameter :: maxan = 0
 integer, parameter :: maxmhdan = 0
 integer, parameter :: maxdustan = 0
#else
 integer, parameter :: maxan = maxp
 integer, parameter :: maxmhdan = maxmhd
 integer, parameter :: maxdustan = maxp_dustfrac
#endif

 integer, parameter :: maxphase = maxan
 integer, parameter :: maxgradh = maxan

 character(len=10) :: sizestring

 write(iprint, *)
 write(iprint, '(a)') '--> ALLOCATING ARRAYS'
 write(iprint, '(a)') '--------------------------------------------------------'

 nbytes_allocated = 0.0

 call allocate_array('xyzh', xyzh, 4, maxp)
 call allocate_array('xyzh_soa', xyzh_soa, maxp, 4)
 call allocate_array('vxyzu', vxyzu, maxvxyzu, maxp)
 call allocate_array('alphaind', alphaind, nalpha, maxalpha)
 call allocate_array('divcurlv', divcurlv, ndivcurlv, maxp)
 call allocate_array('divcurlB', divcurlB, ndivcurlB, maxp)
 call allocate_array('Bevol', Bevol, maxBevol, maxmhd)
 call allocate_array('Bxyz', Bxyz, 3, maxmhd)
 call allocate_array('dustprop', dustprop, 5, maxp_growth)
 call allocate_array('St', St, maxp_growth)
 call allocate_array('straintensor', straintensor, 6, maxstrain)
 call allocate_array('abundance', abundance, nabundances, maxp_h2)
 call allocate_array('temperature', temperature, maxtemp)
 call allocate_array('dustfrac', dustfrac, maxp_dustfrac)
 call allocate_array('dustevol', dustevol, maxp_dustfrac)
 call allocate_array('deltav', deltav, 3, maxp_dustfrac)
 call allocate_array('xyzmh_ptmass', xyzmh_ptmass, nsinkproperties, maxptmass)
 call allocate_array('vxyz_ptmass', vxyz_ptmass, 3, maxptmass)
 call allocate_array('fxyz_ptmass', fxyz_ptmass, 4, maxptmass)
 call allocate_array('fxyz_ptmass_sinksink', fxyz_ptmass_sinksink, 4, maxptmass)
 call allocate_array('poten', poten, maxgrav)
 call allocate_array('n_R', n_R, 4, maxmhdni)
 call allocate_array('n_electronT', n_electronT, maxne)
 call allocate_array('eta_nimhd', eta_nimhd, 4, maxmhdni)
 call allocate_array('luminosity', luminosity, maxlum)
 call allocate_array('fxyzu', fxyzu, maxvxyzu, maxan)
 call allocate_array('dBevol', dBevol, maxBevol, maxmhdan)
 call allocate_array('divBsumm', divBsymm, maxmhdan)
 call allocate_array('fext', fext, 3, maxan)
 call allocate_array('ddustfrac', ddustfrac, maxdustan)
 call allocate_array('ddustprop', ddustprop, 5, maxp_growth)
 call allocate_array('vpred', vpred, maxvxyzu, maxan)
 call allocate_array('dustpred', dustpred, maxdustan)
 call allocate_array('Bpred', Bpred, maxBevol, maxmhdan)
 call allocate_array('dustproppred', dustproppred, 5, maxp_growth)
#ifdef IND_TIMESTEPS
 call allocate_array('ibin', ibin, maxan)
 call allocate_array('ibin_old', ibin_old, maxan)
 call allocate_array('ibin_wake', ibin_wake, maxan)
 call allocate_array('dt_in', dt_in, maxan)
 call allocate_array('twas', twas, maxan)
#endif
 call allocate_array('iphase', iphase, maxphase)
 call allocate_array('iphase_soa', iphase_soa, maxphase)
 call allocate_array('gradh', gradh, ngradh, maxgradh)
 call allocate_array('tstop', tstop, maxan)
 call allocate_array('ll', ll, maxan)

 call bytes2human(nbytes_allocated, sizestring)
 write(iprint, '(a)') '--------------------------------------------------------'
 write(iprint, *) 'Total memory allocated to arrays: ', sizestring
 write(iprint, '(a)') '--------------------------------------------------------'

end subroutine allocate_memory

end module memory
