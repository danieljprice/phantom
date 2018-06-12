!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: mpiutils
!
!  DESCRIPTION:
! This module contains MPI-related quantities and utilities
! though most can be called safely from non-MPI code
! (but obviously do nothing).
!
! In particular, we implement a number of generic interfaces
! to MPI calls that make the Fortran code a LOT nicer.
! These include:
!
! * reduce_mpi:
!
!   var = reduce_mpi('+',var)
!
!   where var can be int,int*8,real*4,real*8 or even arrays
!   and the reduction operators are '+','max' or 'min'.
!
! * reduceall_mpi:
!
!   As above but gets result on all processors.
!
! * bcast_mpi:
!
!   call bcast_mpi(var)
!
!   where var can be int,int*8,real*4,real*8.
!
! * send_recv:
!
!   call send_recv(arr,listpart,isendto,irecvfrom,itag,ioffset,nrecv)
!
!   sends and receives selected values from array of arbitrary type
!
! * barrier_mpi
!
!   calls MPI_BARRIER, no-op if called from non-MPI code
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, mpi
!+
!--------------------------------------------------------------------------
module mpiutils
#ifdef MPI
 use mpi
 implicit none
 integer, public :: mpierr
 integer, public :: status(MPI_STATUS_SIZE)
 integer, public :: MPI_DEFAULT_REAL

 integer, public :: comm_cellexchange, comm_cellcount, comm_balance, comm_balancecount

!
!--generic interface send_recv
!
 interface send_recv
  module procedure send_recv_arr2, send_recv_arr2_r4, send_recv_arr1_r4, &
                    send_recv_int, send_recv_arr1_int1, send_recv_arr2_buf, &
                    send_recv_arr1_buf_r4
 end interface

 public :: send_recv, cart_shift_diag
 logical, parameter, public :: use_mpi = .true.
#else
 implicit none
 logical, parameter, public :: use_mpi = .false.
#endif
!
!--generic interface reduce_mpi
!
 interface reduce_mpi
  module procedure reduce_mpi_real, reduce_mpi_real4, reduce_mpi_int, reduce_mpi_int8, &
                     reduce_mpi_int_arr, reduce_mpi_int8_arr, reduce_mpi_real4arr, reduce_mpi_real8arr
 end interface
!
!--generic interface reduceall_mpi
!
 interface reduceall_mpi
  module procedure reduceall_mpi_real, reduceall_mpi_real4, reduceall_mpi_int, reduceall_mpi_int8, reduceall_mpi_int1, &
                     reduceall_mpi_realarr, reduceall_mpi_real4arr, reduceall_mpi_int4arr
 end interface
 !
 !--generic interface reduceloc_mpi
 !
 interface reduceloc_mpi
  module procedure reduceloc_mpi_real4,reduceloc_mpi_real8,reduceloc_mpi_int
 end interface
!
!  generic interface reduce_in_place
!
 interface reduce_in_place_mpi
  module procedure reduce_in_place_mpi_real8arr2, reduce_in_place_mpi_real4arr2
 end interface
!
!--generic interface bcast_mpi
!
 interface bcast_mpi
  module procedure bcast_mpi_int1, bcast_mpi_int, bcast_mpi_int8, bcast_mpi_real4, bcast_mpi_real8, &
                   bcast_mpi_real8arr, bcast_mpi_real4arr, bcast_mpi_real8arr2, bcast_mpi_real4arr2
 end interface
!
!--generic interface fill_buffer
!
 interface fill_buffer
  module procedure fill_buffer_r8,fill_buffer_r4,fill_buffer_r8val,fill_buffer_r4val,fill_buffer_ival
 end interface
!
!--generic interface unfill_buf
!
 interface unfill_buf
  module procedure unfill_bufarr,unfill_buf1
 end interface

 public :: init_mpi, finalise_mpi
 public :: waitmyturn,endmyturn
 public :: reduce_mpi, reduceall_mpi, reduce_in_place_mpi
 public :: bcast_mpi
 public :: barrier_mpi
 public :: fill_buffer, unfill_buf
 public :: reduceloc_mpi

 private

contains
!--------------------------------------------------------------------
!+
!  initialisation of MPI, including initialisation of extra MPI types
!+
!--------------------------------------------------------------------
subroutine init_mpi(id,nprocs)
#ifdef MPI
 use io, only:fatal,master
#endif
 integer, intent(out) :: id,nprocs
#ifdef MPI
 real :: xtemp
!
!--initialise MPI - get number of threads and id of current processor
!
 call MPI_INIT(mpierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,mpierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,id,mpierr)
 if (mpierr /= 0) call fatal('init_mpi','error starting mpi')

 call MPI_COMM_DUP(MPI_COMM_WORLD,comm_cellexchange,mpierr)
 call MPI_COMM_DUP(MPI_COMM_WORLD,comm_cellcount,mpierr)
 call MPI_COMM_DUP(MPI_COMM_WORLD,comm_balance,mpierr)
 call MPI_COMM_DUP(MPI_COMM_WORLD,comm_balancecount,mpierr)

 if (id==master) print "(a,i0,a)",' running in MPI on ',nprocs,' threads'
!
!--also work out the size of default real for MPI communications
!
 select case(kind(xtemp))
 case(8)
    MPI_DEFAULT_REAL = MPI_REAL8
 case(4)
    MPI_DEFAULT_REAL = MPI_REAL4
 case default
    call fatal('init_mpi','cannot determine kind for default real')
 end select
#else
 id = 0
 nprocs = 1
 print "(a)",' MPI parallelisation is OFF'
#endif

end subroutine init_mpi

!---------------------------------------------
!+
!  finalisation of MPI, no-op if MPI not set
!+
!---------------------------------------------
subroutine finalise_mpi()
#ifdef MPI
 use io, only:fatal

 call MPI_COMM_FREE(comm_cellexchange,mpierr)
 call MPI_COMM_FREE(comm_cellcount,mpierr)
 call MPI_COMM_FREE(comm_balance,mpierr)
 call MPI_COMM_FREE(comm_balancecount,mpierr)

 call MPI_FINALIZE(mpierr)
 if (mpierr /= 0) call fatal('reduce','error in mpi_finalize call')
#endif

end subroutine finalise_mpi

!--------------------------------------------------------------------
!+
!  MPI utility to write in turn (no-op if called in non-mpi code)
!+
!--------------------------------------------------------------------
subroutine waitmyturn(myid)
 integer, intent(in) :: myid
#ifdef MPI
 integer :: nowgo

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

 nowgo = 0
 if (myid > 0) call MPI_RECV(nowgo,1,MPI_INTEGER,myid-1,99,MPI_COMM_WORLD,status,mpierr)

#endif

end subroutine waitmyturn

!--------------------------------------------------------------------
!+
! MPI utility to end writing in turn (no-op if called in non-mpi code)
!+
!--------------------------------------------------------------------
subroutine endmyturn(myid)
#ifdef MPI
 use io, only:nprocs
#endif
 integer, intent(in) :: myid
#ifdef MPI
 integer :: nowgo

 if (myid < nprocs-1) then
    nowgo = 1
    call MPI_SEND(nowgo,1,MPI_INTEGER,myid+1,99,MPI_COMM_WORLD,mpierr)
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

#endif

end subroutine endmyturn

!--------------------------------------------------------------------
!+
!  MPI barrier interface (no-op if called in non-mpi code)
!+
!--------------------------------------------------------------------
subroutine barrier_mpi()
#ifdef MPI
 use io, only:fatal

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('barrier_mpi','error in mpi_barrier call')

#endif

end subroutine barrier_mpi


!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on real numbers
!  can be called from non-MPI routines
!+
!--------------------------------------------------------------------------
real(kind=8) function reduce_mpi_real(string,xproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 real(kind=8),     intent(in) :: xproc
#ifdef MPI
 real(kind=8) :: xred,xsend

 xsend = xproc  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_REDUCE(xsend,xred,1,MPI_REAL8,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(xsend,xred,1,MPI_REAL8,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(xsend,xred,1,MPI_REAL8,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_real = xred
#else
 reduce_mpi_real = xproc
#endif

end function reduce_mpi_real

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on real*4 numbers
!  can be called from non-MPI routines
!+
!--------------------------------------------------------------------------
real(kind=4) function reduce_mpi_real4(string,xproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 real(kind=4),     intent(in) :: xproc
#ifdef MPI
 real(kind=4) :: xred,xsend

 xsend = xproc  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_REDUCE(xsend,xred,1,MPI_REAL4,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(xsend,xred,1,MPI_REAL4,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(xsend,xred,1,MPI_REAL4,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_real4 = xred
#else
 reduce_mpi_real4 = xproc
#endif

end function reduce_mpi_real4

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on array
!  of real*4 numbers. Can be called from non-MPI routines.
!+
!--------------------------------------------------------------------------
function reduce_mpi_real4arr(string,xproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 real(kind=4),     intent(in) :: xproc(:)
 real(kind=4) :: reduce_mpi_real4arr(size(xproc))
#ifdef MPI
 real(kind=4) :: xred(size(xproc)),xsend(size(xproc))

 xsend(:) = xproc(:)  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_REDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_real4arr(:) = xred(:)
#else
 reduce_mpi_real4arr(:) = xproc(:)
#endif

end function reduce_mpi_real4arr

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on array
!  of real*4 numbers. Can be called from non-MPI routines.
!+
!--------------------------------------------------------------------------
function reduce_mpi_real8arr(string,xproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 real(kind=8),     intent(in) :: xproc(:)
 real(kind=8) :: reduce_mpi_real8arr(size(xproc))
#ifdef MPI
 real(kind=8) :: xred(size(xproc)),xsend(size(xproc))

 xsend(:) = xproc(:)  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_REDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_real8arr(:) = xred(:)
#else
 reduce_mpi_real8arr(:) = xproc(:)
#endif

end function reduce_mpi_real8arr

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on integers
!  can be called from non-MPI routines
!  NB: returns INT*8
!+
!--------------------------------------------------------------------------
integer(kind=8) function reduce_mpi_int(string,iproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 integer,          intent(in) :: iproc
#ifdef MPI
 integer(kind=8) :: isend,ired

 isend = iproc  ! convert to int*8
 select case(trim(string))
 case('+')
    call MPI_REDUCE(isend,ired,1,MPI_INTEGER8,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(isend,ired,1,MPI_INTEGER8,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(isend,ired,1,MPI_INTEGER8,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_int = ired
#else
 reduce_mpi_int = iproc
#endif

end function reduce_mpi_int

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on
!  array of integers
!  can be called from non-MPI routines
!  NB: returns INT*8
!+
!--------------------------------------------------------------------------
function reduce_mpi_int_arr(string,iproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 integer,          intent(in) :: iproc(:)
 integer(kind=8) :: reduce_mpi_int_arr(size(iproc))
#ifdef MPI
 integer(kind=8) :: isend(size(iproc)),ired(size(iproc))

 isend(:) = iproc(:)  ! convert to int*8
 select case(trim(string))
 case('+')
    call MPI_REDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_int_arr(:) = ired(:)
#else
 reduce_mpi_int_arr(:) = iproc(:)
#endif

end function reduce_mpi_int_arr

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on int 8's
!  can be called from non-MPI routines
!+
!--------------------------------------------------------------------------
integer(kind=8) function reduce_mpi_int8(string,iproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 integer(kind=8),  intent(in) :: iproc
#ifdef MPI
 integer(kind=8) :: isend,ired

 isend = iproc  ! copy
 select case(trim(string))
 case('+')
    call MPI_REDUCE(isend,ired,1,MPI_INTEGER8,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(isend,ired,1,MPI_INTEGER8,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(isend,ired,1,MPI_INTEGER8,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_int8 = ired
#else
 reduce_mpi_int8 = iproc
#endif

end function reduce_mpi_int8

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on
!  array of integers
!  can be called from non-MPI routines
!  NB: returns INT*8
!+
!--------------------------------------------------------------------------
function reduce_mpi_int8_arr(string,iproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 integer(kind=8),  intent(in) :: iproc(:)
 integer(kind=8) :: reduce_mpi_int8_arr(size(iproc))
#ifdef MPI
 integer(kind=8) :: isend(size(iproc)),ired(size(iproc))

 isend(:) = iproc(:)  ! copy
 select case(trim(string))
 case('+')
    call MPI_REDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_REDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_MAX,master,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_REDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_MIN,master,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduce_mpi_int8_arr(:) = ired(:)
#else
 reduce_mpi_int8_arr(:) = iproc(:)
#endif

end function reduce_mpi_int8_arr

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on real numbers
!  can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
real(kind=8) function reduceall_mpi_real(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 real(kind=8),     intent(in) :: xproc
#ifdef MPI
 real(kind=8) :: xred,xsend

 xsend = xproc  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')

 reduceall_mpi_real = xred
#else
 reduceall_mpi_real = xproc
#endif

end function reduceall_mpi_real

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on real numbers
!  can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
real(kind=4) function reduceall_mpi_real4(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 real(kind=4),     intent(in) :: xproc
#ifdef MPI
 real(kind=4) :: xred,xsend

 xsend = xproc  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_REAL4,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_REAL4,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')

 reduceall_mpi_real4 = xred
#else
 reduceall_mpi_real4 = xproc
#endif

end function reduceall_mpi_real4

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on array
!  of real*8 numbers. Can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
function reduceall_mpi_realarr(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 real(kind=8),     intent(in) :: xproc(:)
 real(kind=8) :: reduceall_mpi_realarr(size(xproc))
#ifdef MPI
 real(kind=8) :: xred(size(xproc)),xsend(size(xproc))

 xsend(:) = xproc(:)  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')

 reduceall_mpi_realarr(:) = xred(:)
#else
 reduceall_mpi_realarr(:) = xproc(:)
#endif

end function reduceall_mpi_realarr

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on array
!  of real*4 numbers. Can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
function reduceall_mpi_real4arr(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 real(kind=4),     intent(in) :: xproc(:)
 real(kind=4) :: reduceall_mpi_real4arr(size(xproc))
#ifdef MPI
 real(kind=4) :: xred(size(xproc)),xsend(size(xproc))

 xsend(:) = xproc(:)  ! mpi calls don't like it if send and receive addresses are the same
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')

 reduceall_mpi_real4arr(:) = xred(:)
#else
 reduceall_mpi_real4arr(:) = xproc(:)
#endif

end function reduceall_mpi_real4arr

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on
!  array of integers
!  can be called from non-MPI routines
!  NB: returns INT*8
!+
!--------------------------------------------------------------------------
function reduceall_mpi_int4arr(string,iproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 integer(kind=4),  intent(in) :: iproc(:)
 integer(kind=8) :: reduceall_mpi_int4arr(size(iproc))
#ifdef MPI
 integer(kind=8) :: isend(size(iproc)),ired(size(iproc))

 isend(:) = iproc(:)  ! copy
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,size(isend),MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduceall_mpi_int4arr(:) = ired(:)
#else
 reduceall_mpi_int4arr(:) = iproc(:)
#endif

end function reduceall_mpi_int4arr


!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on int 8's
!  can be called from non-MPI routines
!+
!--------------------------------------------------------------------------
integer(kind=8) function reduceall_mpi_int(string,iproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 integer,          intent(in) :: iproc
#ifdef MPI
 integer(kind=8) :: isend,ired

 isend = iproc  ! copy
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')

 reduceall_mpi_int = ired
#else
 reduceall_mpi_int = iproc
#endif

end function reduceall_mpi_int

!--------------------------------------------------------------------------
!+
!  As above but for int8 arg
!+
!--------------------------------------------------------------------------
integer(kind=8) function reduceall_mpi_int8(string,iproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 integer(kind=8),  intent(in) :: iproc
#ifdef MPI
 integer(kind=8) :: isend,ired

 isend = iproc  ! copy
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')

 reduceall_mpi_int8 = ired
#else
 reduceall_mpi_int8 = iproc
#endif

end function reduceall_mpi_int8

!--------------------------------------------------------------------------
!+
!  function performing MPI reduction operations (+,max,min) on int 8's
!  can be called from non-MPI routines
!+
!--------------------------------------------------------------------------
integer(kind=1) function reduceall_mpi_int1(string,iproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 integer(kind=1),  intent(in) :: iproc
#ifdef MPI
 integer :: isend,ired

 isend = iproc  ! copy
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')

 if (ired > huge(0_1)) then
    call fatal('reduceall (mpi)',' overflow in integer*1 during reduction operation '//trim(string))
    reduceall_mpi_int1 = int(ired,kind=1) ! to prevent compiler warnings
 else
    reduceall_mpi_int1 = int(ired,kind=1)
 endif
#else
 reduceall_mpi_int1 = iproc
#endif

end function reduceall_mpi_int1


!--------------------------------------------------------------------------
!+
!  function performing in-place MPI reduction operations (+,max,min) on array
!  of real*8 numbers. Can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
subroutine reduce_in_place_mpi_real8arr2(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in)    :: string
 real(kind=8),     intent(inout) :: xproc(:,:)
#ifdef MPI
 real(kind=8) :: xsend(size(xproc(:,1)),size(xproc(1,:)))

 xsend = xproc
 select case(trim(string))
 case('+')
!    call MPI_ALLREDUCE(MPI_IN_PLACE,xproc,size(xproc),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(xsend,xproc,size(xsend),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
! case('max')
!    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
! case('min')
!    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce_in_place (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce_in_place','error in mpi_reduce call')

#endif

end subroutine reduce_in_place_mpi_real8arr2

!--------------------------------------------------------------------------
!+
!  function performing in-place MPI reduction operations (+,max,min) on array
!  of real*4 numbers. Can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
subroutine reduce_in_place_mpi_real4arr2(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in)    :: string
 real(kind=4),     intent(inout) :: xproc(:,:)
#ifdef MPI
 real(kind=4) :: xsend(size(xproc(:,1)),size(xproc(1,:)))

 xsend = xproc
 select case(trim(string))
 case('+')
!    call MPI_ALLREDUCE(MPI_IN_PLACE,xproc,size(xproc),MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(xsend,xproc,size(xsend),MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,mpierr)
! case('max')
!    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_MAX,MPI_COMM_WORLD,mpierr)
! case('min')
!    call MPI_ALLREDUCE(xsend,xred,size(xsend),MPI_REAL4,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce_in_place (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce_in_place','error in mpi_reduce call')

#endif

end subroutine reduce_in_place_mpi_real4arr2

!--------------------------------------------------------------------------
!+
!  min/max reduction identifying which proc has the maximum (real*8)
!+
!--------------------------------------------------------------------------
subroutine reduceloc_mpi_real8(string,xproc,loc)
 use io, only:fatal,id
 character(len=*), intent(in)    :: string
 real(kind=8),     intent(inout) :: xproc
 integer,          intent(out)   :: loc
#ifdef MPI
 real(kind=8) :: xred(2),xsend(2)

 xsend(1) = xproc
 xsend(2) = float(id)
 select case(trim(string))
 case('max')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')
 xproc = xred(1)
 loc = int(xred(2))
#else
 loc = id
#endif

end subroutine reduceloc_mpi_real8

!--------------------------------------------------------------------------
!+
!  min/max reduction identifying which proc has the maximum (real*4)
!+
!--------------------------------------------------------------------------
subroutine reduceloc_mpi_real4(string,xproc,loc)
 use io, only:fatal,id
 character(len=*), intent(in)    :: string
 real(kind=4),     intent(inout) :: xproc
 integer,          intent(out)   :: loc
#ifdef MPI
 real(kind=4) :: xred(2),xsend(2)

 xsend(1) = xproc
 xsend(2) = float(id)
 select case(trim(string))
 case('max')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_2REAL,MPI_MAXLOC,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_2REAL,MPI_MINLOC,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')
 xproc = xred(1)
 loc = int(xred(2))
#else
 loc = id
#endif

end subroutine reduceloc_mpi_real4

!--------------------------------------------------------------------------
!+
!  min/max reduction identifying which proc has the maximum (int*4)
!+
!--------------------------------------------------------------------------
subroutine reduceloc_mpi_int(string,xproc,loc)
 use io, only:fatal,id
 character(len=*), intent(in)    :: string
 integer,          intent(inout) :: xproc
 integer,          intent(out)   :: loc
#ifdef MPI
 integer :: xred(2),xsend(2)

 xsend(1) = xproc
 xsend(2) = id
 select case(trim(string))
 case('max')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_2INTEGER,MPI_MAXLOC,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xred,1,MPI_2INTEGER,MPI_MINLOC,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduceall','error in mpi_reduce call')
 xproc = xred(1)
 loc = xred(2)
#else
 loc = id
#endif

end subroutine reduceloc_mpi_int

#ifdef MPI
!----------------------------------------------------------------
!+
!  MPI utility to send and receive two-dimensional
!  indexed arrays between processors (real*8)
!+
!----------------------------------------------------------------
subroutine send_recv_arr2(xarr,listpart,isendto,irecvfrom,itag,ioffset,nrecv)
 real(kind=8), intent(in)    :: xarr(:,:)
 integer,      intent(in)    :: listpart(:)
 integer :: lblocklengths(size(listpart))
 integer,      intent(in)    :: isendto,irecvfrom,itag,ioffset
 integer,      intent(inout) :: nrecv
 integer :: itypearr,isendtype

 lblocklengths(:) = 1

! print*,id,' sending  ',size(listpart),' to   ',isendto

 CALL MPI_TYPE_CONTIGUOUS(size(xarr(:,1)),MPI_REAL8,itypearr,mpierr)
 CALL MPI_TYPE_COMMIT(itypearr,mpierr)

 if (size(listpart) > 0) then
    !CALL MPI_TYPE_CREATE_INDEXED_BLOCK(size(listpart),1,listpart,itypearr,isendtype,mpierr)
    CALL MPI_TYPE_INDEXED(size(listpart),lblocklengths,listpart,itypearr,isendtype,mpierr)
    CALL MPI_TYPE_COMMIT(isendtype,mpierr)

! print*,id,' sending  x= ',xarr(:,listpart(1)+1),listpart(1)+1,' to   ',isendto

    CALL MPI_SENDRECV(xarr,1,isendtype,isendto,itag,xarr(1,ioffset),nrecv,itypearr,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
    CALL MPI_TYPE_FREE(isendtype,mpierr)
 else
    CALL MPI_SENDRECV(xarr,0,itypearr,isendto,itag,xarr(1,ioffset),nrecv,itypearr,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
 endif

 CALL MPI_GET_COUNT(status,itypearr,nrecv,mpierr)
! print*,id,' received ',nrecv,' from ',status(MPI_SOURCE),' x= ',xarr(:,ioffset)
! print*,id,' received x= ',xarr(:,ioffset),ioffset,' from ',irecvfrom

 CALL MPI_TYPE_FREE(itypearr,mpierr)

 return
end subroutine send_recv_arr2

!----------------------------------------------------------------
!+
!  MPI utility to send and receive two dimensional
!  indexed arrays between processors (real*4)
!+
!----------------------------------------------------------------
subroutine send_recv_arr2_r4(xarr,listpart,isendto,irecvfrom,itag,ioffset,nrecv)
 real(kind=4), intent(in)    :: xarr(:,:)
 integer,      intent(in)    :: listpart(:)
 integer :: lblocklengths(size(listpart))
 integer,      intent(in)    :: isendto,irecvfrom,itag,ioffset
 integer,      intent(inout) :: nrecv
 integer :: itypearr,isendtype

 lblocklengths(:) = 1
 CALL MPI_TYPE_CONTIGUOUS(size(xarr(:,1)),MPI_REAL4,itypearr,mpierr)
 CALL MPI_TYPE_COMMIT(itypearr,mpierr)

 if (size(listpart) > 0) then
    CALL MPI_TYPE_INDEXED(size(listpart),lblocklengths,listpart,itypearr,isendtype,mpierr)
    CALL MPI_TYPE_COMMIT(isendtype,mpierr)

    ! print*,id,' sending  x= ',xarr(:,listpart(1)+1),listpart(1)+1,' to   ',isendto

    CALL MPI_SENDRECV(xarr,1,isendtype,isendto,itag,xarr(1,ioffset),nrecv,itypearr,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
    CALL MPI_TYPE_FREE(isendtype,mpierr)
 else
    CALL MPI_SENDRECV(xarr,0,itypearr,isendto,itag,xarr(1,ioffset),nrecv,itypearr,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
 endif

 CALL MPI_GET_COUNT(status,itypearr,nrecv,mpierr)

! print*,id,' received ',nrecv,' from ',status(MPI_SOURCE),' x= ',xarr(:,ioffset)
! print*,id,' received x= ',xarr(:,ioffset),ioffset,' from ',irecvfrom

 CALL MPI_TYPE_FREE(itypearr,mpierr)

 return
end subroutine send_recv_arr2_r4

!----------------------------------------------------------------
!+
!  MPI utility to send and receive one dimensional
!  indexed arrays between processors (real*4)
!+
!----------------------------------------------------------------
subroutine send_recv_arr1_r4(xarr,listpart,isendto,irecvfrom,itag,ioffset,nrecv)
 real(kind=4), intent(in)    :: xarr(:)
 integer,      intent(in)    :: listpart(:)
 integer :: lblocklengths(size(listpart))
 integer,      intent(in)    :: isendto,irecvfrom,itag,ioffset
 integer,      intent(inout) :: nrecv
 integer :: isendtype

 lblocklengths(:) = 1

 if (size(listpart) > 0) then
    CALL MPI_TYPE_INDEXED(size(listpart),lblocklengths,listpart,MPI_REAL4,isendtype,mpierr)
    CALL MPI_TYPE_COMMIT(isendtype,mpierr)

    ! print*,id,' sending  x= ',xarr(listpart(1)+1),listpart(1)+1,' to   ',isendto

    CALL MPI_SENDRECV(xarr,1,isendtype,isendto,itag,xarr(ioffset),nrecv,MPI_REAL4,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
    CALL MPI_TYPE_FREE(isendtype,mpierr)
 else
    CALL MPI_SENDRECV(xarr,0,MPI_REAL4,isendto,itag,xarr(ioffset),nrecv,MPI_REAL4,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
 endif

 CALL MPI_GET_COUNT(status,MPI_REAL4,nrecv,mpierr)

! print*,id,' received ',nrecv,' from ',status(MPI_SOURCE),' x= ',xarr(:,ioffset)
! print*,id,' received x= ',xarr(ioffset),ioffset,' from ',irecvfrom

 return
end subroutine send_recv_arr1_r4

!----------------------------------------------------------------
!+
!  MPI utility to send and receive one dimensional
!  indexed arrays between processors (int*1)
!+
!----------------------------------------------------------------
subroutine send_recv_arr1_int1(iarr,listpart,isendto,irecvfrom,itag,ioffset,nrecv)
 integer(kind=1), intent(inout) :: iarr(:)
 integer,         intent(in)    :: listpart(:)
 integer :: lblocklengths(size(listpart))
 integer,         intent(in)    :: isendto,irecvfrom,itag,ioffset
 integer,         intent(inout) :: nrecv
 integer :: isendtype

 lblocklengths(:) = 1

 if (size(listpart) > 0) then
    CALL MPI_TYPE_INDEXED(size(listpart),lblocklengths,listpart,MPI_INTEGER1,isendtype,mpierr)
    CALL MPI_TYPE_COMMIT(isendtype,mpierr)

    ! print*,id,' sending  x= ',xarr(listpart(1)+1),listpart(1)+1,' to   ',isendto

    CALL MPI_SENDRECV(iarr,1,isendtype,isendto,itag,iarr(ioffset),nrecv,MPI_INTEGER1,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
    CALL MPI_TYPE_FREE(isendtype,mpierr)
 else
    CALL MPI_SENDRECV(iarr,0,MPI_INTEGER1,isendto,itag,iarr(ioffset),nrecv,MPI_INTEGER1,irecvfrom,&
                      itag,MPI_COMM_WORLD,status,mpierr)
 endif

 CALL MPI_GET_COUNT(status,MPI_INTEGER1,nrecv,mpierr)

! print*,id,' received ',nrecv,' from ',status(MPI_SOURCE),' i= ',iarr(:,ioffset)
! print*,id,' received i= ',iarr(ioffset),ioffset,' from ',irecvfrom

 return
end subroutine send_recv_arr1_int1

!----------------------------------------------------------------
!+
!  MPI utility to send and receive a single integer
!  (part of generic send_recv interface)
!+
!----------------------------------------------------------------
subroutine send_recv_int(ival,ivalrecv,isendto,irecvfrom,itag)
 integer, intent(in)  :: ival
 integer, intent(out) :: ivalrecv
 integer, intent(in)  :: isendto,irecvfrom,itag

 CALL MPI_SENDRECV(ival,1,MPI_INTEGER,isendto,itag,ivalrecv,1,MPI_INTEGER,irecvfrom,&
                   itag,MPI_COMM_WORLD,status,mpierr)

 return
end subroutine send_recv_int

!----------------------------------------------------------------
!+
!  MPI utility to send and receive a two dimensional
!  contiguous array and return the result in a different array
!+
!----------------------------------------------------------------
subroutine send_recv_arr2_buf(xarr,xrecv,ioffset,nsend,isendto,irecvfrom,itag,nrecv)
 real,    intent(in)  :: xarr(:,:)
 real,    intent(out) :: xrecv(:,:)
 integer, intent(in)  :: ioffset,nsend,isendto,irecvfrom,itag
 integer, intent(out) :: nrecv
 integer :: itypearr

 CALL MPI_TYPE_CONTIGUOUS(size(xarr(:,1)),MPI_DEFAULT_REAL,itypearr,mpierr)
 CALL MPI_TYPE_COMMIT(itypearr,mpierr)

! print*,id,' sending  x= ',xarr(:,ioffset),ioffset,' to   ',isendto

 CALL MPI_SENDRECV(xarr(1,ioffset),nsend,itypearr,isendto,itag,xrecv,size(xrecv(1,:)),itypearr,irecvfrom,&
                   itag,MPI_COMM_WORLD,status,mpierr)

 CALL MPI_GET_COUNT(status,itypearr,nrecv,mpierr)

! print*,id,' received ',nrecv,' from ',status(MPI_SOURCE),' x= ',xarr(:,ioffset)
! print*,id,' received x= ',xrecv(:,1),1,' from ',irecvfrom

 CALL MPI_TYPE_FREE(itypearr,mpierr)

 return
end subroutine send_recv_arr2_buf

!-----------------------------------------------------------------------
!+
!  MPI utility to send and receive a one dimensional
!  contiguous array and return the result in a different array (real*4)
!+
!-----------------------------------------------------------------------
subroutine send_recv_arr1_buf_r4(xarr,xrecv,ioffset,nsend,isendto,irecvfrom,itag,nrecv)
 real(kind=4), intent(in)  :: xarr(:)
 real(kind=4), intent(out) :: xrecv(:)
 integer,      intent(in)  :: ioffset,nsend,isendto,irecvfrom,itag
 integer,      intent(out) :: nrecv

! print*,id,' sending  x= ',xarr(ioffset),ioffset,' to   ',isendto

 CALL MPI_SENDRECV(xarr(ioffset),nsend,MPI_REAL4,isendto,itag,xrecv,size(xrecv),MPI_REAL4, &
                   irecvfrom,itag,MPI_COMM_WORLD,status,mpierr)

 CALL MPI_GET_COUNT(status,MPI_REAL4,nrecv,mpierr)

! print*,id,' received x= ',xrecv(1),1,' from ',irecvfrom

 return
end subroutine send_recv_arr1_buf_r4

!-----------------------------------------------------------------------
!+
!  MPI utility to perform arbitrary cartesian shift
!+
!-----------------------------------------------------------------------
subroutine cart_shift_diag(comm_cart,icoords,ishift,irecvfrom,isendto)
 integer, intent(in)  :: comm_cart
 integer, intent(in)  :: icoords(:),ishift(:)
 integer, intent(out) :: isendto,irecvfrom
 integer :: icoordsshift(size(icoords))

 icoordsshift = icoords + ishift
 call MPI_CART_RANK(comm_cart,icoordsshift,isendto,mpierr)

 icoordsshift = icoords - ishift
 call MPI_CART_RANK(comm_cart,icoordsshift,irecvfrom,mpierr)

 return
end subroutine cart_shift_diag
#endif

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (integer*1)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_int1(ival,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 integer(kind=1), intent(inout) :: ival
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(ival,1,MPI_INTEGER1,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_int1

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (integer)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_int(ival,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 integer, intent(inout) :: ival
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(ival,1,MPI_INTEGER,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_int

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (integer*8)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_int8(ival,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 integer(kind=8), intent(inout) :: ival
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(ival,1,MPI_INTEGER8,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_int8

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (real*4)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_real4(rval,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 real(kind=4), intent(inout) :: rval
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(rval,1,MPI_REAL4,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_real4

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (real*8)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_real8(dval,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 real(kind=8), intent(inout) :: dval
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(dval,1,MPI_REAL8,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_real8

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (real*8 1d array)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_real8arr(dval,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 real(kind=8), intent(inout) :: dval(:)
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif
 call MPI_BCAST(dval,size(dval),MPI_REAL8,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_real8arr

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (real*4 1d array)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_real4arr(dval,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 real(kind=4), intent(inout) :: dval(:)
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(dval,size(dval),MPI_REAL4,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_real4arr

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (real*8 2d array)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_real8arr2(dval,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 real(kind=8), intent(inout) :: dval(:,:)
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif
 call MPI_BCAST(dval,size(dval),MPI_REAL8,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_real8arr2

!--------------------------------------------------------------------------
!+
!  function performing MPI BROADCAST (real*4 2d array)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_real4arr2(dval,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 real(kind=4), intent(inout) :: dval(:,:)
 integer, optional, intent(in) :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(dval,size(dval),MPI_REAL4,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_real4arr2

!--------------------------------------------------------------------------
!+
!  functions to fill a buffer to send to different processors
!+
!--------------------------------------------------------------------------
subroutine fill_buffer_r8(xbuffer,xarr,nbuf)
 real,         intent(inout) :: xbuffer(:)
 real(kind=8), intent(in)    :: xarr(:)
 integer,      intent(inout) :: nbuf
 integer :: i1

 i1 = nbuf + 1
 nbuf = i1 + size(xarr) - 1
 xbuffer(i1:nbuf) = real(xarr(:))

end subroutine fill_buffer_r8

subroutine fill_buffer_r4(xbuffer,xarr,nbuf)
 real,         intent(inout) :: xbuffer(:)
 real(kind=4), intent(in)    :: xarr(:)
 integer,      intent(inout) :: nbuf
 integer :: i1

 i1 = nbuf + 1
 nbuf = i1 + size(xarr) - 1
 xbuffer(i1:nbuf) = real(xarr(:))

end subroutine fill_buffer_r4

subroutine fill_buffer_r8val(xbuffer,xval,nbuf)
 real,         intent(inout) :: xbuffer(:)
 real(kind=8), intent(in)    :: xval
 integer,      intent(inout) :: nbuf

 nbuf = nbuf + 1
 xbuffer(nbuf) = real(xval)

end subroutine fill_buffer_r8val

subroutine fill_buffer_r4val(xbuffer,xval,nbuf)
 real,         intent(inout) :: xbuffer(:)
 real(kind=4), intent(in)    :: xval
 integer,      intent(inout) :: nbuf

 nbuf = nbuf + 1
 xbuffer(nbuf) = xval

end subroutine fill_buffer_r4val

subroutine fill_buffer_ival(xbuffer,ival,nbuf)
 real,            intent(inout) :: xbuffer(:)
 integer(kind=1), intent(in)    :: ival
 integer,         intent(inout) :: nbuf

 nbuf = nbuf + 1
 xbuffer(nbuf) = real(ival)

end subroutine fill_buffer_ival

!--------------------------------------------------------------------------
!+
!  functions to unfill a buffer received from different processors
!+
!--------------------------------------------------------------------------
function unfill_bufarr(xbuf,ioffset,len) result(xnew)
 integer, intent(inout) :: ioffset
 real,    intent(in)    :: xbuf(:)
 integer, intent(in)    :: len
 real :: xnew(len)
 integer :: i1,i2

 i1 = ioffset + 1
 i2 = ioffset + len
 xnew = xbuf(i1:i2)
 ioffset = i2

end function unfill_bufarr

function unfill_buf1(xbuf,ioffset) result(xnew)
 integer, intent(inout) :: ioffset
 real,    intent(in)    :: xbuf(:)
 real    :: xnew
 integer :: i1

 i1 = ioffset + 1
 xnew = xbuf(i1)
 ioffset = i1

end function unfill_buf1

end module mpiutils
