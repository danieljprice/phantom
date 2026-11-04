!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpiutils
!
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
! * reduce_in_place_mpi:
!
!   call reduce_in_place_mpi('+',var,n)
!
!   where var is a 1D array and n is the number of elements to reduce
!   (n may be smaller than the declared size). Also accepts a whole
!   2D array, with no length argument.
!
! * allgather_mpi:
!
!   call allgather_mpi(x,xall)
!
!   gathers one int x from each task into xall(1:nprocs)
!   on all tasks.
!
! * allgatherv_mpi:
!
!   call allgatherv_mpi(xsend,xrecv[,counts,displs])
!
!   concatenates the int or real*8 arrays xsend of all tasks, in task
!   order, into xrecv on all tasks. The number of entries from each task
!   and the offset of each block in xrecv may be passed in counts and
!   displs (e.g. from get_displs_mpi), otherwise they are gathered.
!
! * alltoall_mpi, alltoallv_mpi:
!
!   call alltoall_mpi(nsend,nrecv)
!   call alltoallv_mpi(xsend,nsend,xrecv,nrecv)
!
!   exchange integer counts, then blocks of real*8, between every pair
!   of tasks: block p of xsend (nsend(p) entries) goes to task p-1, and
!   block p of xrecv (nrecv(p) entries) comes from task p-1.
!
! * bcast_mpi:
!
!   call bcast_mpi(var)
!
!   where var can be int,int*8,real*4,real*8.
!
! * barrier_mpi
!
!   calls MPI_BARRIER, no-op if called from non-MPI code
!
! * isend_mpi, irecv_mpi, waitall_mpi
!
!   call isend_mpi(buf,iproc,itag,ireq), likewise irecv_mpi, then
!   call waitall_mpi(nreq,ireq): non-blocking messages
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, mpi, mpidens, mpiforce
!
#ifdef MPI
 use mpi
 implicit none
#else
 implicit none
 integer, parameter, private :: MPI_STATUS_SIZE = 0
#endif

 integer, public :: mpierr
 integer, public :: status(MPI_STATUS_SIZE)
 integer, public :: MPI_DEFAULT_REAL
 integer, public :: comm_cellexchange,comm_cellcount,comm_balance,comm_balancecount
#ifdef MPI
 integer, private :: fh_dump = -1   ! the dump file's MPI-IO handle (mpiio_open)
#endif

!
!--generic interface reduce_mpi
!
 interface reduce_mpi
  module procedure reduce_mpi_real,reduce_mpi_real4,reduce_mpi_int,reduce_mpi_int8, &
                     reduce_mpi_int_arr,reduce_mpi_int8_arr,reduce_mpi_real4arr,reduce_mpi_real8arr
 end interface reduce_mpi
!
!--generic interface reduceall_mpi
!
 interface reduceall_mpi
  module procedure reduceall_mpi_real,reduceall_mpi_real4,reduceall_mpi_int,reduceall_mpi_int8,reduceall_mpi_int1, &
                     reduceall_mpi_realarr,reduceall_mpi_realarr2,reduceall_mpi_real4arr,reduceall_mpi_real4arr2, &
                     reduceall_mpi_int4arr,reduceall_mpi_int1arr
 end interface reduceall_mpi
 !
 !--generic interface reduceloc_mpi
 !
 interface reduceloc_mpi
  module procedure reduceloc_mpi_real4,reduceloc_mpi_real8,reduceloc_mpi_int
 end interface reduceloc_mpi
!
!  generic interface reduce_in_place
!
 interface reduce_in_place_mpi
  module procedure reduce_in_place_mpi_real8arr,reduce_in_place_mpi_real4arr, &
                   reduce_in_place_mpi_real8arr2,reduce_in_place_mpi_real4arr2
 end interface reduce_in_place_mpi
!
!--generic interface bcast_mpi
!
 interface bcast_mpi
  module procedure bcast_mpi_int1,bcast_mpi_int,bcast_mpi_int8,bcast_mpi_real4,bcast_mpi_real8, &
                   bcast_mpi_real8arr,bcast_mpi_real4arr,bcast_mpi_real8arr2,bcast_mpi_real4arr2,&
                   bcast_mpi_logical,bcast_mpi_intarr2
 end interface bcast_mpi
!
!--generic interface fill_buffer
!
 interface fill_buffer
  module procedure fill_buffer_r8,fill_buffer_r4,fill_buffer_r8val,fill_buffer_r4val,fill_buffer_i1val,fill_buffer_i8val
 end interface fill_buffer
!
!--generic interface unfill_buf
!
 interface unfill_buf
  module procedure unfill_bufarr,unfill_buf1
 end interface unfill_buf
!
!--generic interface allgather_mpi
!
 interface allgather_mpi
  module procedure allgather_mpi_int
 end interface allgather_mpi
!
!--generic interface allgatherv_mpi
!
 interface allgatherv_mpi
  module procedure allgatherv_mpi_int,allgatherv_mpi_real4,allgatherv_mpi_real8
 end interface allgatherv_mpi
!
!--generic interface alltoallv_mpi
!
 interface alltoallv_mpi
  module procedure alltoallv_mpi_real4,alltoallv_mpi_real8
 end interface alltoallv_mpi
!
!--generic interfaces for non-blocking point-to-point messages
!
 interface isend_mpi
  module procedure isend_mpi_real4,isend_mpi_real8
 end interface isend_mpi

 interface irecv_mpi
  module procedure irecv_mpi_real4,irecv_mpi_real8
 end interface irecv_mpi

 public :: init_mpi, finalise_mpi
 public :: reduce_mpi, reduceall_mpi, reduce_in_place_mpi
 public :: bcast_mpi
 public :: allgather_mpi, allgatherv_mpi, get_displs_mpi
 public :: alltoall_mpi, alltoallv_mpi
 public :: isend_mpi, irecv_mpi, waitall_mpi
 public :: barrier_mpi
 public :: exscan_mpi_int8
 public :: mpiio_open,mpiio_write_at,mpiio_close,mpiio_error_string
 public :: fill_buffer, unfill_buf
 public :: reduceloc_mpi
 public :: waitmyturn,endmyturn
 public :: start_threadwrite,end_threadwrite

 private

contains
!--------------------------------------------------------------------
!+
!  initialisation of MPI, including initialisation of extra MPI types
!+
!--------------------------------------------------------------------
subroutine init_mpi(id,nprocs)
#ifdef MPI
 use io,       only:fatal,master
 use mpidens,  only:get_mpitype_of_celldens
 use mpiforce, only:get_mpitype_of_cellforce
#endif
 integer, intent(out) :: id,nprocs
#ifdef MPI
 real :: xtemp
 integer :: provided
!
!--initialise MPI - get number of threads and id of current processor
!
 call MPI_INIT_THread(MPI_THREAD_MULTIPLE,provided,mpierr)
 if (provided /= MPI_THREAD_MULTIPLE) call fatal('init_mpi','error initialising threaded MPI')

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
!  utility for initialising each thread
!+
!--------------------------------------------------------------------
subroutine start_threadwrite(id,iunit,filename)
#ifdef MPI
 use mpi
#endif
 use io, only:error,iverbose
 implicit none
 integer,          intent(in) :: id, iunit
 character(len=*), intent(in) :: filename
 integer :: nowgo,ierr

 if (iverbose >= 3) print *,id,' : starting write...'
 nowgo = 0
 if (id  >  0) then
#ifdef MPI
    call MPI_RECV(nowgo,1,MPI_INTEGER,id-1,99,MPI_COMM_WORLD,status,mpierr)
#endif
    open(unit=iunit,file=filename,status='old',form='unformatted',position='append',iostat=ierr)
    if (ierr /= 0) then
       call error('start_threadwrite','can''t append to dumpfile '//trim(filename))
    else
       if (iverbose >= 3) print*,'thread ',id,': opened file '//trim(filename)
    endif
 endif

end subroutine start_threadwrite

!--------------------------------------------------------------------
!+
!  utility for finalising each thread
!+
!--------------------------------------------------------------------
subroutine end_threadwrite(id)
 use io,       only:iverbose
#ifdef MPI
 use mpi
 use io,       only:nprocs
#endif
 implicit none
 integer, intent(in) :: id
#ifdef MPI
 integer :: nowgo
#endif

 if (iverbose >= 3) print *,' thread ',id,' : finished write.'
#ifdef MPI
 if (id  <  nprocs-1) then
    nowgo = 1
    call MPI_SEND(nowgo,1,MPI_INTEGER,id+1,99,MPI_COMM_WORLD,mpierr)
 endif
#endif

end subroutine end_threadwrite

!--------------------------------------------------------------------
!+
!  exclusive prefix sum over the MPI tasks: the sum of ival over the
!  tasks below this one (0 on the first, and without MPI)
!+
!--------------------------------------------------------------------
function exscan_mpi_int8(ival) result(isum)
#ifdef MPI
 use io, only:fatal,id,master
#endif
 integer(kind=8), intent(in) :: ival
 integer(kind=8) :: isum
#ifdef MPI
 isum = 0_8
 call MPI_EXSCAN(ival,isum,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('exscan','error in mpi_exscan')
 if (id == master) isum = 0_8   ! undefined there in MPI
#else
 isum = 0_8
#endif

end function exscan_mpi_int8

!--------------------------------------------------------------------
!+
!  MPI-IO on one file shared by all the MPI tasks: open an existing file
!  for writing with the ROMIO hint direct_write, and say whether direct I/O
!  is in effect (direct): the MPI library reports direct_write true on its
!  Lustre driver (ROMIO-based MPI-IO: Intel MPI, MPICH, Open MPI with
!  --mca io romio341; Open MPI's default OMPIO reports no hints, ROMIO's
!  generic driver does no direct I/O).  Write bytes at a 0-based offset;
!  close.  ierr is the MPI error code (0 = success); without MPI they
!  return ierr=1.
!+
!--------------------------------------------------------------------
subroutine mpiio_open(filename,ierr,direct)
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 logical,          intent(out) :: direct
#ifdef MPI
 integer :: info,infoused,ierrinfo
 logical :: got
 character(len=64) :: val

 direct = .false.
 call MPI_INFO_CREATE(info,ierrinfo)
 if (ierrinfo == MPI_SUCCESS) then
    call MPI_INFO_SET(info,'direct_write','true',ierrinfo)
 else
    info = MPI_INFO_NULL
 endif
 call MPI_FILE_open(MPI_COMM_WORLD,filename,MPI_MODE_WRONLY,info,fh_dump,ierr)
 if (info /= MPI_INFO_NULL) call MPI_INFO_FREE(info,ierrinfo)
 if (ierr /= MPI_SUCCESS) then
    fh_dump = -1
    return
 endif
 call MPI_FILE_GET_INFO(fh_dump,infoused,ierrinfo)
 if (ierrinfo /= MPI_SUCCESS) return
 ! got is set by the MPI library, whose .true. need not be the compiler's (Intel MPI's is 1,
 ! ifx's is -1, so that .not.got is true too): use it only to test, never in an expression
 call MPI_INFO_GET(infoused,'direct_write',len(val),val,got,ierrinfo)
 if (ierrinfo /= MPI_SUCCESS) got = .false.
 if (got) direct = (trim(val) == 'true')
 if (direct) then
    direct = .false.
    call MPI_INFO_GET(infoused,'romio_filesystem_type',len(val),val,got,ierrinfo)
    if (ierrinfo /= MPI_SUCCESS) got = .false.
    if (got) direct = (index(val,'LUSTRE') == 1)
 endif
 call MPI_INFO_FREE(infoused,ierrinfo)
#else
 ierr   = 1
 direct = .false.
#endif

end subroutine mpiio_open

subroutine mpiio_write_at(off,b,ierr)
 integer(kind=8), intent(in)  :: off
 integer(kind=1), intent(in), contiguous :: b(:)
 integer,         intent(out) :: ierr
#ifdef MPI
 integer(kind=MPI_OFFSET_KIND) :: moff
 integer(kind=8) :: n,done,c
 integer :: st(MPI_STATUS_SIZE),nwritten,ierrc

 ierr = 0
 n    = size(b,kind=8)
 done = 0
 do while (done < n .and. ierr == 0)
    c    = min(n - done,1073741824_8)
    moff = int(off + done,kind=MPI_OFFSET_KIND)
    call MPI_FILE_WRITE_AT(fh_dump,moff,b(done+1:done+c),int(c),MPI_BYTE,st,ierr)
    if (ierr == MPI_SUCCESS) then
       ! a write that ends short, with no error, is an error too
       call MPI_GET_COUNT(st,MPI_BYTE,nwritten,ierrc)
       if (ierrc /= MPI_SUCCESS .or. nwritten /= int(c)) ierr = MPI_ERR_IO
    endif
    done = done + c
 enddo
#else
 ierr = 1
#endif

end subroutine mpiio_write_at

subroutine mpiio_close(ierr)
 integer, intent(out) :: ierr
#ifdef MPI
 ierr = 0
 if (fh_dump /= -1) call MPI_FILE_close(fh_dump,ierr)
 fh_dump = -1
#else
 ierr = 0
#endif

end subroutine mpiio_close

! the MPI library's description of an MPI error code
function mpiio_error_string(ierr) result(string)
 integer, intent(in) :: ierr
 character(len=256)  :: string
#ifdef MPI
 integer :: lenstring,ierr2

 string = ' '
 call MPI_ERROR_STRING(ierr,string,lenstring,ierr2)
 if (ierr2 /= MPI_SUCCESS) write(string,'(a,i0)') 'MPI error code ',ierr
#else
 write(string,'(a,i0)') 'error code ',ierr
#endif

end function mpiio_error_string

!--------------------------------------------------------------------
!+
!  MPI barrier interface (no-op if called in non-mpi code)
!+
!--------------------------------------------------------------------
subroutine barrier_mpi()
#ifdef MPI
 use io, only:fatal
 integer :: mpierr_local

 ! call MPI_BARRIER on the master task only, but also barrier the OMP threads

 !$omp barrier
 !$omp masked
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr_local)
 if (mpierr_local /= 0) call fatal('barrier_mpi','error in mpi_barrier call')
 !$omp end masked
 !$omp barrier

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
!  function performing MPI reduction operations (+,max,min) on 2-d array
!  of real*8 numbers. Can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
function reduceall_mpi_realarr2(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 real(kind=8),     intent(in) :: xproc(:,:)
 real(kind=8) :: reduceall_mpi_realarr2(size(xproc,1),size(xproc,2))
#ifdef MPI
 real(kind=8) :: xred(size(xproc,1),size(xproc,2)),xsend(size(xproc,1),size(xproc,2))

 xsend(:,:) = xproc(:,:)  ! mpi calls don't like it if send and receive addresses are the same
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

 reduceall_mpi_realarr2(:,:) = xred(:,:)
#else
 reduceall_mpi_realarr2(:,:) = xproc(:,:)
#endif

end function reduceall_mpi_realarr2

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
!  function performing MPI reduction operations (+,max,min) on 2-d array
!  of real*4 numbers. Can be called from non-MPI routines.
!  Sends result to all threads.
!+
!--------------------------------------------------------------------------
function reduceall_mpi_real4arr2(string,xproc)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in) :: string
 real(kind=4),     intent(in) :: xproc(:,:)
 real(kind=4) :: reduceall_mpi_real4arr2(size(xproc,1),size(xproc,2))
#ifdef MPI
 real(kind=4) :: xred(size(xproc,1),size(xproc,2)),xsend(size(xproc,1),size(xproc,2))

 xsend(:,:) = xproc(:,:)  ! mpi calls don't like it if send and receive addresses are the same
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

 reduceall_mpi_real4arr2(:,:) = xred(:,:)
#else
 reduceall_mpi_real4arr2(:,:) = xproc(:,:)
#endif

end function reduceall_mpi_real4arr2

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
!  function performing MPI reduction operations (+,max,min) on
!  array of integers
!  can be called from non-MPI routines
!  NB: returns INT*1
!+
!--------------------------------------------------------------------------
function reduceall_mpi_int1arr(string,iproc)
#ifdef MPI
 use io, only:fatal,master
#endif
 character(len=*), intent(in) :: string
 integer(kind=1),  intent(in) :: iproc(:)
 integer(kind=1) :: reduceall_mpi_int1arr(size(iproc))
#ifdef MPI
 integer(kind=1) :: isend(size(iproc)),ired(size(iproc))

 isend(:) = iproc(:)  ! copy
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,size(isend),MPI_INTEGER1,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,size(isend),MPI_INTEGER1,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,size(isend),MPI_INTEGER1,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce','error in mpi_reduce call')

 reduceall_mpi_int1arr(:) = ired(:)
#else
 reduceall_mpi_int1arr(:) = iproc(:)
#endif

end function reduceall_mpi_int1arr

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
!  In-place MPI reduction (+, max or min) of the first n elements of a
!  real*8 array. n may be smaller than the actual argument. The result
!  is returned in those elements on every MPI rank. No-op without MPI.
!+
!--------------------------------------------------------------------------
subroutine reduce_in_place_mpi_real8arr(string,xproc,n)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in)    :: string
 integer,          intent(in)    :: n
 real(kind=8),     intent(inout) :: xproc(n)
#ifdef MPI
 real(kind=8) :: xsend(n)

 xsend = xproc
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(xsend,xproc,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(xsend,xproc,n,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xproc,n,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce_in_place (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce_in_place','error in mpi_reduce call')

#endif

end subroutine reduce_in_place_mpi_real8arr

!--------------------------------------------------------------------------
!+
!  In-place MPI reduction (+, max or min) of the first n elements of a
!  real*4 array. n may be smaller than the actual argument. The result
!  is returned in those elements on every MPI rank. No-op without MPI.
!+
!--------------------------------------------------------------------------
subroutine reduce_in_place_mpi_real4arr(string,xproc,n)
#ifdef MPI
 use io, only:fatal
#endif
 character(len=*), intent(in)    :: string
 integer,          intent(in)    :: n
 real(kind=4),     intent(inout) :: xproc(n)
#ifdef MPI
 real(kind=4) :: xsend(n)

 xsend = xproc
 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(xsend,xproc,n,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,mpierr)
 case('max')
    call MPI_ALLREDUCE(xsend,xproc,n,MPI_REAL4,MPI_MAX,MPI_COMM_WORLD,mpierr)
 case('min')
    call MPI_ALLREDUCE(xsend,xproc,n,MPI_REAL4,MPI_MIN,MPI_COMM_WORLD,mpierr)
 case default
    call fatal('reduce_in_place (mpi)','unknown reduction operation')
 end select
 if (mpierr /= 0) call fatal('reduce_in_place','error in mpi_reduce call')

#endif

end subroutine reduce_in_place_mpi_real4arr

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
 xsend(2) = real(id)
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
 xsend(2) = real(id,kind=4)
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
 integer,         intent(in), optional :: src
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
 integer, intent(in), optional :: src
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
!  function performing MPI BROADCAST (logical)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_logical(lval,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 logical, intent(inout) :: lval
 integer, intent(in), optional :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif

 call MPI_BCAST(lval,1,MPI_LOGICAL,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_logical

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
 integer,         intent(in), optional :: src
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
 integer,      intent(in), optional :: src
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
 integer,      intent(in), optional :: src
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
 integer,      intent(in), optional :: src
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
 integer,      intent(in), optional :: src
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
 integer,      intent(in), optional :: src
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
 integer,      intent(in), optional :: src
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
!  function performing MPI BROADCAST (int 2D array)
!+
!--------------------------------------------------------------------------
subroutine bcast_mpi_intarr2(ival,src)
#ifdef MPI
 use io, only:fatal,master
#endif
 integer, intent(inout) :: ival(:,:)
 integer, intent(in), optional :: src
#ifdef MPI
 integer :: sendsrc
 if (present(src)) then
    sendsrc = src
 else
    sendsrc = master
 endif
 call MPI_BCAST(ival,size(ival),MPI_INTEGER4,sendsrc,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('bcast','error in mpi_bcast')

#endif

end subroutine bcast_mpi_intarr2

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

subroutine fill_buffer_i1val(xbuffer,ival,nbuf)
 real,            intent(inout) :: xbuffer(:)
 integer(kind=1), intent(in)    :: ival
 integer,         intent(inout) :: nbuf

 nbuf = nbuf + 1
 xbuffer(nbuf) = real(ival)

end subroutine fill_buffer_i1val

subroutine fill_buffer_i8val(xbuffer,ival,nbuf)
 real,            intent(inout) :: xbuffer(:)
 integer(kind=8), intent(in)    :: ival
 integer,         intent(inout) :: nbuf

 nbuf = nbuf + 1
 xbuffer(nbuf) = real(ival)

end subroutine fill_buffer_i8val

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

!--------------------------------------------------------------------------
!+
!  gather one integer from each task into an array on all tasks
!+
!--------------------------------------------------------------------------
subroutine allgather_mpi_int(x,xall)
#ifdef MPI
 use io, only:fatal
#endif
 integer, intent(in)  :: x
 integer, intent(out) :: xall(:)

#ifdef MPI
 call MPI_ALLGATHER(x,1,MPI_INTEGER,xall,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('allgather_mpi','error in MPI_ALLGATHER call')
#else
 xall(:) = x
#endif

end subroutine allgather_mpi_int

!--------------------------------------------------------------------------
!+
!  offsets of the blocks of each task in a concatenated array,
!  given the number of entries from each task
!+
!--------------------------------------------------------------------------
pure function get_displs_mpi(counts) result(displs)
 integer, intent(in) :: counts(:)
 integer :: displs(size(counts))
 integer :: i

 if (size(counts) < 1) return

 ! each block starts after the previous one
 displs(1) = 0
 do i = 2,size(counts)
    displs(i) = displs(i-1) + counts(i-1)
 enddo

end function get_displs_mpi

!--------------------------------------------------------------------------
!+
!  get the number of entries each task sends, if not given,
!  and the offset of each task's block
!+
!--------------------------------------------------------------------------
subroutine get_counts_displs(nsend,counts_out,displs_out,counts,displs)
 integer, intent(in)  :: nsend
 integer, intent(out) :: counts_out(:),displs_out(:)
 integer, intent(in), optional :: counts(:),displs(:)

 if (present(counts)) then
    counts_out = counts
 else
    call allgather_mpi_int(nsend,counts_out)
 endif
 if (present(displs)) then
    displs_out = displs
 else
    displs_out = get_displs_mpi(counts_out)
 endif

end subroutine get_counts_displs

!--------------------------------------------------------------------------
!+
!  concatenate integer arrays from all tasks, in task order, on all tasks
!+
!--------------------------------------------------------------------------
subroutine allgatherv_mpi_int(xsend,xrecv,counts,displs)
 use io, only:nprocs,fatal
 integer, intent(in)  :: xsend(:)
 integer, intent(out) :: xrecv(:)
 integer, intent(in), optional :: counts(:),displs(:)
 integer :: icounts(nprocs),idispls(nprocs)

 call get_counts_displs(size(xsend),icounts,idispls,counts,displs)
 if (sum(icounts) > size(xrecv)) call fatal('allgatherv_mpi','receive array too small')
#ifdef MPI
 call MPI_ALLGATHERV(xsend,size(xsend),MPI_INTEGER,xrecv,icounts,idispls,MPI_INTEGER,&
                     MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('allgatherv_mpi','error in MPI_ALLGATHERV call')
#else
 xrecv(1:size(xsend)) = xsend
#endif

end subroutine allgatherv_mpi_int

!--------------------------------------------------------------------------
!+
!  concatenate real*4 arrays from all tasks, in task order, on all tasks
!+
!--------------------------------------------------------------------------
subroutine allgatherv_mpi_real4(xsend,xrecv,counts,displs)
 use io, only:nprocs,fatal
 real(kind=4), intent(in)  :: xsend(:)
 real(kind=4), intent(out) :: xrecv(:)
 integer,      intent(in), optional :: counts(:),displs(:)
 integer :: icounts(nprocs),idispls(nprocs)

 call get_counts_displs(size(xsend),icounts,idispls,counts,displs)
 if (sum(icounts) > size(xrecv)) call fatal('allgatherv_mpi','receive array too small')
#ifdef MPI
 call MPI_ALLGATHERV(xsend,size(xsend),MPI_REAL4,xrecv,icounts,idispls,MPI_REAL4,&
                     MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('allgatherv_mpi','error in MPI_ALLGATHERV call')
#else
 xrecv(1:size(xsend)) = xsend
#endif

end subroutine allgatherv_mpi_real4

!--------------------------------------------------------------------------
!+
!  concatenate real*8 arrays from all tasks, in task order, on all tasks
!+
!--------------------------------------------------------------------------
subroutine allgatherv_mpi_real8(xsend,xrecv,counts,displs)
 use io, only:nprocs,fatal
 real(kind=8), intent(in)  :: xsend(:)
 real(kind=8), intent(out) :: xrecv(:)
 integer,      intent(in), optional :: counts(:),displs(:)
 integer :: icounts(nprocs),idispls(nprocs)

 call get_counts_displs(size(xsend),icounts,idispls,counts,displs)
 if (sum(icounts) > size(xrecv)) call fatal('allgatherv_mpi','receive array too small')
#ifdef MPI
 call MPI_ALLGATHERV(xsend,size(xsend),MPI_REAL8,xrecv,icounts,idispls,MPI_REAL8,&
                     MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('allgatherv_mpi','error in MPI_ALLGATHERV call')
#else
 xrecv(1:size(xsend)) = xsend
#endif

end subroutine allgatherv_mpi_real8

!--------------------------------------------------------------------------
!+
!  send nsend(p) to task p-1 and receive nrecv(p) from it, for all p
!+
!--------------------------------------------------------------------------
subroutine alltoall_mpi(nsend,nrecv)
#ifdef MPI
 use io, only:fatal
#endif
 integer, intent(in)  :: nsend(:)
 integer, intent(out) :: nrecv(:)

#ifdef MPI
 call MPI_ALLTOALL(nsend,1,MPI_INTEGER,nrecv,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('alltoall_mpi','error in MPI_ALLTOALL call')
#else
 nrecv = nsend
#endif

end subroutine alltoall_mpi

!--------------------------------------------------------------------------
!+
!  exchange blocks of real*4 between all pairs of tasks: block p of
!  xsend (nsend(p) entries) goes to task p-1, block p of xrecv
!  (nrecv(p) entries, as given by alltoall_mpi) comes from task p-1
!+
!--------------------------------------------------------------------------
subroutine alltoallv_mpi_real4(xsend,nsend,xrecv,nrecv)
 use io, only:fatal
 real(kind=4), intent(in)  :: xsend(:)
 integer,      intent(in)  :: nsend(:),nrecv(:)
 real(kind=4), intent(out) :: xrecv(:)

 if (sum(nrecv) > size(xrecv)) call fatal('alltoallv_mpi','receive array too small')
#ifdef MPI
 call MPI_ALLTOALLV(xsend,nsend,get_displs_mpi(nsend),MPI_REAL4,&
                    xrecv,nrecv,get_displs_mpi(nrecv),MPI_REAL4,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('alltoallv_mpi','error in MPI_ALLTOALLV call')
#else
 xrecv(1:nsend(1)) = xsend(1:nsend(1))
#endif

end subroutine alltoallv_mpi_real4

!--------------------------------------------------------------------------
!+
!  exchange blocks of real*8 between all pairs of tasks: block p of
!  xsend (nsend(p) entries) goes to task p-1, block p of xrecv
!  (nrecv(p) entries, as given by alltoall_mpi) comes from task p-1
!+
!--------------------------------------------------------------------------
subroutine alltoallv_mpi_real8(xsend,nsend,xrecv,nrecv)
 use io, only:fatal
 real(kind=8), intent(in)  :: xsend(:)
 integer,      intent(in)  :: nsend(:),nrecv(:)
 real(kind=8), intent(out) :: xrecv(:)

 if (sum(nrecv) > size(xrecv)) call fatal('alltoallv_mpi','receive array too small')
#ifdef MPI
 call MPI_ALLTOALLV(xsend,nsend,get_displs_mpi(nsend),MPI_REAL8,&
                    xrecv,nrecv,get_displs_mpi(nrecv),MPI_REAL8,MPI_COMM_WORLD,mpierr)
 if (mpierr /= 0) call fatal('alltoallv_mpi','error in MPI_ALLTOALLV call')
#else
 xrecv(1:nsend(1)) = xsend(1:nsend(1))
#endif

end subroutine alltoallv_mpi_real8

!--------------------------------------------------------------------------
!+
!  start a non-blocking send of buf (real*4) to task iproc
!  (buf must not change until the request completes, see waitall_mpi)
!+
!--------------------------------------------------------------------------
subroutine isend_mpi_real4(buf,iproc,itag,ireq)
 use io, only:fatal
 integer,      intent(in)  :: iproc,itag
 real(kind=4), intent(in), asynchronous, contiguous :: buf(:)
 integer,      intent(out) :: ireq
#ifdef MPI
 integer :: ierr

 call MPI_ISEND(buf,size(buf),MPI_REAL4,iproc,itag,MPI_COMM_WORLD,ireq,ierr)
 if (ierr /= 0) call fatal('isend_mpi','error in MPI_ISEND call')
#else
 ireq = 0
 if (size(buf) > 0 .and. iproc /= 0) call fatal('isend_mpi','send to another task without MPI')
#endif

end subroutine isend_mpi_real4

!--------------------------------------------------------------------------
!+
!  start a non-blocking receive into buf (real*4) from task iproc
!  (buf is filled once the request completes, see waitall_mpi)
!+
!--------------------------------------------------------------------------
subroutine irecv_mpi_real4(buf,iproc,itag,ireq)
 use io, only:fatal
 integer,      intent(in)    :: iproc,itag
 real(kind=4), intent(inout), asynchronous, contiguous :: buf(:)
 integer,      intent(out)   :: ireq
#ifdef MPI
 integer :: ierr

 call MPI_IRECV(buf,size(buf),MPI_REAL4,iproc,itag,MPI_COMM_WORLD,ireq,ierr)
 if (ierr /= 0) call fatal('irecv_mpi','error in MPI_IRECV call')
#else
 ireq = 0
 if (size(buf) > 0 .and. iproc /= 0) call fatal('irecv_mpi','receive from another task without MPI')
#endif

end subroutine irecv_mpi_real4

!--------------------------------------------------------------------------
!+
!  start a non-blocking send of buf (real*8) to task iproc
!  (buf must not change until the request completes, see waitall_mpi)
!+
!--------------------------------------------------------------------------
subroutine isend_mpi_real8(buf,iproc,itag,ireq)
 use io, only:fatal
 integer,      intent(in)  :: iproc,itag
 real(kind=8), intent(in), asynchronous, contiguous :: buf(:)
 integer,      intent(out) :: ireq
#ifdef MPI
 integer :: ierr

 call MPI_ISEND(buf,size(buf),MPI_REAL8,iproc,itag,MPI_COMM_WORLD,ireq,ierr)
 if (ierr /= 0) call fatal('isend_mpi','error in MPI_ISEND call')
#else
 ireq = 0
 if (size(buf) > 0 .and. iproc /= 0) call fatal('isend_mpi','send to another task without MPI')
#endif

end subroutine isend_mpi_real8

!--------------------------------------------------------------------------
!+
!  start a non-blocking receive into buf (real*8) from task iproc
!  (buf is filled once the request completes, see waitall_mpi)
!+
!--------------------------------------------------------------------------
subroutine irecv_mpi_real8(buf,iproc,itag,ireq)
 use io, only:fatal
 integer,      intent(in)    :: iproc,itag
 real(kind=8), intent(inout), asynchronous, contiguous :: buf(:)
 integer,      intent(out)   :: ireq
#ifdef MPI
 integer :: ierr

 call MPI_IRECV(buf,size(buf),MPI_REAL8,iproc,itag,MPI_COMM_WORLD,ireq,ierr)
 if (ierr /= 0) call fatal('irecv_mpi','error in MPI_IRECV call')
#else
 ireq = 0
 if (size(buf) > 0 .and. iproc /= 0) call fatal('irecv_mpi','receive from another task without MPI')
#endif

end subroutine irecv_mpi_real8

!--------------------------------------------------------------------------
!+
!  wait for the non-blocking sends and receives in ireq(1:nreq)
!+
!--------------------------------------------------------------------------
subroutine waitall_mpi(nreq,ireq)
 use io, only:fatal
 integer, intent(in)    :: nreq
 integer, intent(inout) :: ireq(:)
#ifdef MPI
 integer :: ierr

 call MPI_WAITALL(nreq,ireq,MPI_STATUSES_IGNORE,ierr)
 if (ierr /= 0) call fatal('waitall_mpi','error in MPI_WAITALL call')
#endif

end subroutine waitall_mpi

end module mpiutils
