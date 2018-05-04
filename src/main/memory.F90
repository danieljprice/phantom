module memory
use io,  only:fatal,error,iprint

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

use dim,   only:maxp_dustfrac,ndusttypes
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

 implicit none

 public :: allocate_memory
 public :: deallocate_memory

 real :: nbytes_allocated = 0.0

 interface allocate_array
    module procedure &
      allocate_array_real8_1d, &
      allocate_array_real8_2d, &
      allocate_array_real8_3d, &
      allocate_array_real4_1d, &
      allocate_array_real4_2d, &
      allocate_array_real4_3d, &
      allocate_array_integer4_1d, &
      allocate_array_integer4_2d, &
      allocate_array_integer4_3d, &
      allocate_array_integer1_1d, &
      allocate_array_integer1_2d, &
      allocate_array_integer1_3d
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
    real(kind=8), allocatable,   intent(inout)  :: x(:,:)
    integer,                     intent(in)     :: n1, n2
    integer                                     :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'real(8)')
 end subroutine allocate_array_real8_2d

 subroutine allocate_array_real8_3d(name, x, n1, n2, n3)
    character(len=*),            intent(in)     :: name
    real(kind=8), allocatable,   intent(inout)  :: x(:, :, :)
    integer,                     intent(in)     :: n1, n2, n3
    integer                                     :: allocstat

    allocate(x(n1, n2, n3), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2, n3/), 'real(8)')
 end subroutine allocate_array_real8_3d

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
    real(kind=4), allocatable,   intent(inout)  :: x(:,:)
    integer,                     intent(in)     :: n1, n2
    integer                                     :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'real(4)')
 end subroutine allocate_array_real4_2d

 subroutine allocate_array_real4_3d(name, x, n1, n2, n3)
    character(len=*),            intent(in)     :: name
    real(kind=4), allocatable,   intent(inout)  :: x(:, :, :)
    integer,                     intent(in)     :: n1, n2, n3
    integer                                     :: allocstat

    allocate(x(n1, n2, n3), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2, n3/), 'real(4)')
 end subroutine allocate_array_real4_3d

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
    integer(kind=4), allocatable,   intent(inout)  :: x(:,:)
    integer,                        intent(in)     :: n1, n2
    integer                                        :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'integer(4)')
 end subroutine allocate_array_integer4_2d

 subroutine allocate_array_integer4_3d(name, x, n1, n2, n3)
    character(len=*),               intent(in)     :: name
    integer(kind=4), allocatable,   intent(inout)  :: x(:, :, :)
    integer,                        intent(in)     :: n1, n2, n3
    integer                                        :: allocstat

    allocate(x(n1, n2, n3), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2, n3/), 'integer(4)')
 end subroutine allocate_array_integer4_3d

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
    integer(kind=1), allocatable,   intent(inout)  :: x(:,:)
    integer,                        intent(in)     :: n1, n2
    integer                                        :: allocstat

    allocate(x(n1, n2), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2/), 'integer(1)')
 end subroutine allocate_array_integer1_2d

 subroutine allocate_array_integer1_3d(name, x, n1, n2, n3)
    character(len=*),               intent(in)     :: name
    integer(kind=1), allocatable,   intent(inout)  :: x(:, :, :)
    integer,                        intent(in)     :: n1, n2, n3
    integer                                        :: allocstat

    allocate(x(n1, n2, n3), stat = allocstat)
    call check_allocate(name, allocstat)
    call print_allocation_stats(name, (/n1, n2, n3/), 'integer(1)')
 end subroutine allocate_array_integer1_3d

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

subroutine allocate_memory(maxp)
   integer, intent(in)              :: maxp
   logical :: do_run = .true.
   logical :: do_mhd = .false.
   logical :: do_mhdni = .false.
   logical :: do_Dust = .false.
   character(len=10) :: sizestring

#ifdef ANALYSIS
   do_run = .false.
#endif

#ifdef MHD
   do_mhd = .true.
#ifdef NONIDEALMHD
   do_mhdni = .true.
#endif
#endif

#ifdef DUST
   do_dust = .true.
#endif

 write(iprint, *)
 write(iprint, '(a)') '--> ALLOCATING ARRAYS'
 write(iprint, '(a)') '--------------------------------------------------------'

 if (nbytes_allocated > 0.0) then
    call error('memory', 'Attempting to allocate memory, but memory is already allocated. Deallocating and then allocating again.')
    call deallocate_memory
 endif

 call allocate_array('xyzh', xyzh, 4, maxp)
 call allocate_array('xyzh_soa', xyzh_soa, maxp, 4)
 call allocate_array('vxyzu', vxyzu, maxvxyzu, maxp)
 call allocate_array('alphaind', alphaind, nalpha, maxalpha)
 call allocate_array('divcurlv', divcurlv, ndivcurlv, maxp)
 call allocate_array('divcurlB', divcurlB, ndivcurlB, maxp)
 if (do_mhd) call allocate_array('Bevol', Bevol, maxBevol, maxp)
 if (do_mhd) call allocate_array('Bxyz', Bxyz, 3, maxp)
 call allocate_array('dustprop', dustprop, 5, maxp_growth)
 call allocate_array('St', St, maxp_growth)
 call allocate_array('straintensor', straintensor, 6, maxstrain)
 call allocate_array('abundance', abundance, nabundances, maxp_h2)
 call allocate_array('temperature', temperature, maxtemp)
 if (do_dust .and. do_run) call allocate_array('dustfrac', dustfrac, ndusttypes, maxp)
 if (do_dust .and. do_run) call allocate_array('dustevol', dustevol,ndusttypes, maxp)
 if (do_dust .and. do_run) call allocate_array('deltav', deltav, 3, ndusttypes, maxp)
 call allocate_array('xyzmh_ptmass', xyzmh_ptmass, nsinkproperties, maxptmass)
 call allocate_array('vxyz_ptmass', vxyz_ptmass, 3, maxptmass)
 call allocate_array('fxyz_ptmass', fxyz_ptmass, 4, maxptmass)
 call allocate_array('fxyz_ptmass_sinksink', fxyz_ptmass_sinksink, 4, maxptmass)
 call allocate_array('poten', poten, maxgrav)
 if (do_mhdni) call allocate_array('n_R', n_R, 4, maxp)
 call allocate_array('n_electronT', n_electronT, maxne)
 if (do_mhdni) call allocate_array('eta_nimhd', eta_nimhd, 4, maxp)
 call allocate_array('luminosity', luminosity, maxlum)
 if (do_run) call allocate_array('fxyzu', fxyzu, maxvxyzu, maxp)
 if (do_mhd .and. do_run) call allocate_array('dBevol', dBevol, maxBevol, maxp)
 if (do_mhd .and. do_run) call allocate_array('divBsumm', divBsymm, maxp)
 if (do_run) call allocate_array('fext', fext, 3, maxp)
 if (do_dust .and. do_run) call allocate_array('ddustfrac', ddustfrac, ndusttypes, maxp)
 call allocate_array('ddustprop', ddustprop, 5, maxp_growth)
 if (do_run) call allocate_array('vpred', vpred, maxvxyzu, maxp)
 if (do_dust .and. do_run) call allocate_array('dustpred', dustpred, ndusttypes, maxp)
 if (do_mhd .and. do_run) call allocate_array('Bpred', Bpred, maxBevol, maxp)
 call allocate_array('dustproppred', dustproppred, 5, maxp_growth)
#ifdef IND_TIMESTEPS
 if (do_run) call allocate_array('ibin', ibin, maxp)
 if (do_run) call allocate_array('ibin_old', ibin_old, maxp)
 if (do_run) call allocate_array('ibin_wake', ibin_wake, maxp)
 if (do_run) call allocate_array('dt_in', dt_in, maxp)
 if (do_run) call allocate_array('twas', twas, maxp)
#endif
 if (do_run) call allocate_array('iphase', iphase, maxp)
 if (do_run) call allocate_array('iphase_soa', iphase_soa, maxp)
 if (do_run) call allocate_array('gradh', gradh, ngradh, maxp)
 if (do_run) call allocate_array('tstop', tstop, ndusttypes, maxp)
 if (do_run) call allocate_array('ll', ll, maxp)

 call bytes2human(nbytes_allocated, sizestring)
 write(iprint, '(a)') '--------------------------------------------------------'
 write(iprint, *) 'Total memory allocated to arrays: ', sizestring
 write(iprint, '(a)') '--------------------------------------------------------'

end subroutine allocate_memory

subroutine deallocate_memory
deallocate(xyzh)
deallocate(xyzh_soa)
deallocate(vxyzu)
deallocate(alphaind)
deallocate(divcurlv)
deallocate(divcurlB)
deallocate(Bevol)
deallocate(Bxyz)
deallocate(dustprop)
deallocate(St)
deallocate(straintensor)
deallocate(abundance)
deallocate(temperature)
deallocate(dustfrac)
deallocate(dustevol)
deallocate(deltav)
deallocate(xyzmh_ptmass)
deallocate(vxyz_ptmass)
deallocate(fxyz_ptmass)
deallocate(fxyz_ptmass_sinksink)
deallocate(poten)
deallocate(n_R)
deallocate(n_electronT)
deallocate(eta_nimhd)
deallocate(luminosity)
deallocate(fxyzu)
deallocate(dBevol)
deallocate(divBsymm)
deallocate(fext)
deallocate(ddustfrac)
deallocate(ddustprop)
deallocate(vpred)
deallocate(dustpred)
deallocate(Bpred)
deallocate(dustproppred)
#ifdef IND_TIMESTEPS
deallocate(ibin)
deallocate(ibin_old)
deallocate(ibin_wake)
deallocate(dt_in)
deallocate(twas)
#endif
deallocate(iphase)
deallocate(iphase_soa)
deallocate(gradh)
deallocate(tstop)
deallocate(ll)

end subroutine deallocate_memory

end module memory
