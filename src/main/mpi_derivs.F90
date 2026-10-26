!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpiderivs
!
! This module handles the MPI exchange of cells during the
!   density and force routines
!
! Each task posts nrecvslots receives for each other task, shared by the
! openMP threads. A thread polling for cells tests the receives (unless
! another thread is already doing so), then processes the cells that have
! arrived and re-posts their receives. Each thread sends from a ring of
! nsendbufs buffers, so it only waits for a send to complete when the ring
! wraps around.
!
! All MPI calls here are made inside the critical section mpi_cell_exchange,
! so that only one openMP thread at a time is inside the MPI library.
! Concurrent calls from many threads (allowed by MPI_THREAD_MULTIPLE)
! were found to lose messages and crash inside Open MPI.
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, dim, io, mpi, mpidens, mpiforce, mpimemory,
!   mpiutils, omputils, timing
!
#ifdef MPI
 use mpi
#endif
 use io,             only:id,nprocs
 use mpiutils,       only:comm_cellexchange,comm_cellcount
 use mpidens,        only:celldens
 use mpiforce,       only:cellforce

 implicit none
 interface init_cell_exchange
  module procedure init_celldens_exchange,init_cellforce_exchange
 end interface init_cell_exchange

 interface send_cell
  module procedure send_celldens,send_cellforce
 end interface send_cell

 interface recv_cells
  module procedure recv_celldens,recv_cellforce
 end interface recv_cells

 interface complete_sends
  module procedure complete_sends_dens,complete_sends_force
 end interface complete_sends

 interface combine_cells
  module procedure combine_celldens,combine_cellforce
 end interface combine_cells

 interface finish_cell_exchange
  module procedure finish_celldens_exchange,finish_cellforce_exchange
 end interface finish_cell_exchange

 interface recv_while_wait
  module procedure recv_while_wait_force,recv_while_wait_dens
 end interface recv_while_wait

 public :: allocate_cell_comms_arrays
 public :: deallocate_cell_comms_arrays

 public :: init_cell_exchange
 public :: send_cell
 public :: recv_cells
 public :: finish_cell_exchange
 public :: recv_while_wait
 public :: reset_cell_counters
 public :: combine_cells ! only to prevent compiler warning
 public :: print_cell_exchange_stats

 !
 !--the counters are module variables, but must be passed through as arguments
 !  in order to be threadsafe
 !
 integer, public, allocatable :: cell_counters(:,:)
 integer, parameter :: isent   = 1 ! counter for number of cells sent to i
 integer, parameter :: iexpect = 2 ! counter for number of cells expecting from i
 integer, parameter :: irecv   = 3 ! counter for number of cells received from i

 integer, parameter, public :: nrecvslots = 4 ! receives posted for each other task
 integer, parameter, public :: nsendbufs  = 8 ! send buffers for each openMP thread

 private

 integer, allocatable :: countrequest(:)     ! receives of the number of cells sent by each task
 integer, allocatable :: countsendrequest(:) ! sends of the number of cells sent to each task
 integer, allocatable :: countsent(:)        ! send buffer for the number of cells sent to each task
 integer :: dtype_dens  = -1                 ! MPI datatypes for celldens and cellforce,
 integer :: dtype_force = -1                 ! created once and reused
 integer(kind=8) :: ncells_sent = 0, ncells_recv = 0 ! cells sent and received since last summary

 ! receive buffers and requests, shared by all threads: (slot, source task)
 type(celldens),  allocatable :: recvbuf_dens(:,:)
 type(cellforce), allocatable :: recvbuf_force(:,:)
 integer,         allocatable :: recvreq(:,:)
 ! send buffers for each thread: (buffer, thread), and for each buffer the
 ! requests of the sends in progress: (request, buffer, thread) and their number
 type(celldens),  allocatable :: sendbuf_dens(:,:)
 type(cellforce), allocatable :: sendbuf_force(:,:)
 integer,         allocatable :: sendreq(:,:,:)
 integer,         allocatable :: nsendreq(:,:)
 integer,         allocatable :: lastbuf(:)  ! last send buffer used by each thread

 logical, allocatable :: thread_complete(:)  ! whether each thread has finished sending
 integer :: ncomplete_mpi = 0                ! number of tasks whose cells have all been received
 integer :: npolling = 0                     ! number of threads trying to poll the receives

contains

subroutine allocate_cell_comms_arrays
 use allocutils, only:allocate_array
 use dim,        only:mpi
 if (mpi) then
    call allocate_array('cell_counters', cell_counters, nprocs, 3)
    call allocate_array('countrequest',  countrequest,  nprocs)
    call allocate_array('countsendrequest', countsendrequest, nprocs)
    call allocate_array('countsent', countsent, nprocs)
    countsendrequest = null_request()
 else
    ! dummy cell counters that are required to prevent runtime errors
    ! in dens and force
    call allocate_array('cell_counters', cell_counters, 0, 0)
 endif
end subroutine allocate_cell_comms_arrays

subroutine deallocate_cell_comms_arrays
 if (allocated(cell_counters)) deallocate(cell_counters)
 if (allocated(countrequest )) deallocate(countrequest )
 if (allocated(countsendrequest)) deallocate(countsendrequest)
 if (allocated(countsent)) deallocate(countsent)
 if (allocated(recvbuf_dens)) deallocate(recvbuf_dens)
 if (allocated(recvbuf_force)) deallocate(recvbuf_force)
 if (allocated(recvreq)) deallocate(recvreq)
 if (allocated(sendbuf_dens)) deallocate(sendbuf_dens)
 if (allocated(sendbuf_force)) deallocate(sendbuf_force)
 if (allocated(sendreq)) deallocate(sendreq)
 if (allocated(nsendreq)) deallocate(nsendreq)
 if (allocated(lastbuf)) deallocate(lastbuf)
 if (allocated(thread_complete)) deallocate(thread_complete)
end subroutine deallocate_cell_comms_arrays

!-----------------------------------------------------------------------
!+
!  null request handle
!  (MPI_REQUEST_NULL is 0 in open mpi, but 0x2c000000 for MPICH-derived MPI)
!+
!-----------------------------------------------------------------------
integer function null_request()

#ifdef MPI
 null_request = MPI_REQUEST_NULL
#else
 null_request = 0
#endif

end function null_request

!-----------------------------------------------------------------------
!+
!  allocate the buffers and requests that do not depend on the cell type
!  (called by one thread)
!+
!-----------------------------------------------------------------------
subroutine allocate_exchange_requests
 use omputils, only:omp_num_threads

 if (.not.allocated(recvreq)) allocate(recvreq(nrecvslots,nprocs))
 if (allocated(sendreq)) then
    if (size(sendreq,3) /= omp_num_threads) deallocate(sendreq,nsendreq,lastbuf,thread_complete)
 endif
 if (.not.allocated(sendreq)) then
    allocate(sendreq(nprocs,nsendbufs,omp_num_threads),nsendreq(nsendbufs,omp_num_threads))
    allocate(lastbuf(omp_num_threads),thread_complete(omp_num_threads))
 endif
 sendreq = null_request()
 nsendreq = 0
 lastbuf = 0
 thread_complete = .false.
 ncomplete_mpi = 0
 npolling = 0

end subroutine allocate_exchange_requests

!-----------------------------------------------------------------------
!+
!  query whether all of the given send requests have completed
!+
!-----------------------------------------------------------------------
logical function sends_complete(n,ireq)
 integer, intent(in)    :: n
 integer, intent(inout) :: ireq(n)
#ifdef MPI
 integer :: mpierr

 !$omp critical (mpi_cell_exchange)
 call MPI_TESTALL(n,ireq,sends_complete,MPI_STATUSES_IGNORE,mpierr)
 !$omp end critical (mpi_cell_exchange)
#else
 sends_complete = .true.
#endif

end function sends_complete

!-----------------------------------------------------------------------
!+
!  claim the right to poll the receives: returns true for at most one
!  thread at a time, which must then call release_polling
!+
!-----------------------------------------------------------------------
logical function claim_polling()
 integer :: nprev

 !$omp atomic capture
 nprev = npolling
 npolling = npolling + 1
 !$omp end atomic
 claim_polling = (nprev == 0)

end function claim_polling

subroutine release_polling

 !$omp atomic
 npolling = npolling - 1

end subroutine release_polling

!----------------------------------------------------------------
!+
!  start the exchange of denss: post the receives
!  (called by all threads)
!+
!----------------------------------------------------------------
subroutine init_celldens_exchange(stack)
 use mpidens, only:stackdens,get_mpitype_of_celldens
 type(stackdens), intent(inout) :: stack
#ifdef MPI
 integer :: iproc,islot,mpierr
#endif

 !$omp single
 if (dtype_dens == -1) call get_mpitype_of_celldens(dtype_dens)
 call allocate_exchange_requests
 if (.not.allocated(recvbuf_dens)) allocate(recvbuf_dens(nrecvslots,nprocs))
 if (.not.allocated(sendbuf_dens)) allocate(sendbuf_dens(nsendbufs,size(sendreq,3)))
 if (size(sendbuf_dens,2) /= size(sendreq,3)) then
    deallocate(sendbuf_dens)
    allocate(sendbuf_dens(nsendbufs,size(sendreq,3)))
 endif
#ifdef MPI
 !$omp critical (mpi_cell_exchange)
 do iproc=1,nprocs
    do islot=1,nrecvslots
       call MPI_IRECV(recvbuf_dens(islot,iproc),1,dtype_dens,iproc-1,MPI_ANY_TAG,&
                      comm_cellexchange,recvreq(islot,iproc),mpierr)
    enddo
 enddo
 !$omp end critical (mpi_cell_exchange)
#endif
 !$omp end single

end subroutine init_celldens_exchange

!-----------------------------------------------------------------------
!+
!  send a cell to the tasks in targets, from the next buffer of this
!  thread, receiving cells into stack while waiting for that buffer
!+
!-----------------------------------------------------------------------
subroutine send_celldens(cell,targets,stack,counters)
 use io,       only:fatal
 use omputils, only:omp_thread_num
 use mpidens, only:stackdens
 type(celldens),  intent(in)    :: cell
 logical,         intent(in)    :: targets(nprocs)
 type(stackdens), intent(inout) :: stack
 integer,         intent(inout) :: counters(nprocs,3)
#ifdef MPI
 integer :: newproc,mpierr,ith,ibuf,n

 ith  = omp_thread_num() + 1
 ibuf = mod(lastbuf(ith),nsendbufs) + 1
 lastbuf(ith) = ibuf

 !--wait for the previous sends from this buffer, continuing to receive
 call complete_sends(ibuf,ith,stack,counters)

 sendbuf_dens(ibuf,ith) = cell
 n = 0
 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       n = n + 1
       call MPI_ISEND(sendbuf_dens(ibuf,ith),1,dtype_dens,newproc,1,comm_cellexchange,&
                      sendreq(n,ibuf,ith),mpierr)
       if (mpierr /= 0) call fatal('send_celldens','error in MPI_ISEND')
       !$omp atomic
       counters(newproc+1,isent) = counters(newproc+1,isent) + 1
    endif
 enddo
 !$omp end critical (mpi_cell_exchange)
 nsendreq(ibuf,ith) = n
#endif

end subroutine send_celldens

!-----------------------------------------------------------------------
!+
!  wait for the sends from one buffer of this thread to complete,
!  receiving cells into stack meanwhile
!+
!-----------------------------------------------------------------------
subroutine complete_sends_dens(ibuf,ith,stack,counters)
 use mpidens, only:stackdens
 integer,         intent(in)    :: ibuf,ith
 type(stackdens), intent(inout) :: stack
 integer,         intent(inout) :: counters(nprocs,3)

 do while (nsendreq(ibuf,ith) > 0)
    if (sends_complete(nsendreq(ibuf,ith),sendreq(:,ibuf,ith))) then
       nsendreq(ibuf,ith) = 0
    else
       call recv_cells(stack,counters)
    endif
 enddo

end subroutine complete_sends_dens

!------------------------------------------------
!+
!  receive cells: if no other thread is already
!  doing so, test the receives, then process the
!  cells that have arrived and re-post the receives
!  (other threads may test the receives meanwhile)
!+
!------------------------------------------------
subroutine recv_celldens(target_stack,counters)
 use io,        only:fatal
 use mpimemory, only:push_onto_stack,get_cell_ptr
 use mpidens, only:stackdens
 type(stackdens), intent(inout) :: target_stack
 integer,         intent(inout) :: counters(nprocs,3)
#ifdef MPI
 type(celldens),  pointer :: target_cell
 integer :: indices(nrecvslots*nprocs)
 integer :: nout,k,islot,iproc,mpierr

 nout = 0
 mpierr = 0
 if (claim_polling()) then
    !$omp critical (mpi_cell_exchange)
    call MPI_TESTSOME(nrecvslots*nprocs,recvreq,nout,indices,MPI_STATUSES_IGNORE,mpierr)
    !$omp end critical (mpi_cell_exchange)
 endif
 call release_polling
 if (mpierr /= 0) call fatal('recv_cells','error in MPI_TESTSOME')
 if (nout == MPI_UNDEFINED) nout = 0

 !--completed receives are inactive until re-posted, so no other thread will see them
 do k=1,nout
    islot = mod(indices(k)-1,nrecvslots) + 1
    iproc = (indices(k)-1)/nrecvslots + 1
    if (recvbuf_dens(islot,iproc)%owner == id) then
       target_cell => get_cell_ptr(target_stack,recvbuf_dens(islot,iproc)%waiting_index)
       call combine_cells(target_cell,recvbuf_dens(islot,iproc))
    else
       call push_onto_stack(target_stack,recvbuf_dens(islot,iproc))
    endif
    !$omp atomic
    counters(iproc,irecv) = counters(iproc,irecv) + 1
    !$omp critical (mpi_cell_exchange)
    call MPI_IRECV(recvbuf_dens(islot,iproc),1,dtype_dens,iproc-1,MPI_ANY_TAG,&
                   comm_cellexchange,recvreq(islot,iproc),mpierr)
    !$omp end critical (mpi_cell_exchange)
 enddo
#endif

end subroutine recv_celldens

!-----------------------------------------------------------------------
!+
!  wait until all of the cells sent by this thread have gone,
!  and until all of the cells sent to this task have been received
!  (called by all threads at the end of each phase of the exchange)
!+
!-----------------------------------------------------------------------
subroutine recv_while_wait_dens(stack,counters)
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_thread_num
 use timing,   only:get_timings,increment_timer,itimer_dens_wait
 use mpidens, only:stackdens
 type(stackdens), intent(inout) :: stack
 integer,         intent(inout) :: counters(nprocs,3)
#ifdef MPI
 real(kind=4) :: t1,t2,tcpu1,tcpu2
 integer      :: ith,ibuf

 !$omp masked
 call get_timings(t1,tcpu1)
 !$omp end masked

 !--complete the sends from this thread
 ith = omp_thread_num() + 1
 do ibuf=1,nsendbufs
    call complete_sends(ibuf,ith,stack,counters)
 enddo

 !--signal to other OMP threads that this thread has finished sending
 call set_thread_complete(.true.)

 !--continue receiving cells until all OMP threads have finished sending
 do while (.not. all_threads_complete())
    call recv_cells(stack,counters)
 enddo

 !--signal to other MPI tasks that this task has finished sending
 !$omp masked
 call send_cell_counts(counters)
 !$omp end masked

 !--continue receiving cells until all MPI tasks have finished sending
 do while (.not. all_tasks_complete())
    call recv_cells(stack,counters)
    !$omp masked
    call check_complete(counters)
    !$omp end masked
 enddo

 call barrier_mpi

 !$omp masked
 call finish_cell_counts
 ncomplete_mpi = 0
 ncells_recv = ncells_recv + sum(counters(:,irecv))
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_dens_wait,t2-t1,tcpu2-tcpu1)
 !$omp end masked
 call set_thread_complete(.false.)
#endif

end subroutine recv_while_wait_dens

!----------------------------------------------------------------
!+
!  finish the exchange: cancel the receives that are still posted
!  (called by all threads)
!+
!----------------------------------------------------------------
subroutine finish_celldens_exchange(stack)
 use io,       only:fatal
 use mpiutils, only:barrier_mpi
 use mpidens, only:stackdens
 type(stackdens), intent(inout) :: stack
#ifdef MPI
 integer :: iproc,islot,mpierr
 integer :: status(MPI_STATUS_SIZE)
 logical :: cancelled

 !$omp single
 !$omp critical (mpi_cell_exchange)
 do iproc=1,nprocs
    do islot=1,nrecvslots
       call MPI_CANCEL(recvreq(islot,iproc),mpierr)
       call MPI_WAIT(recvreq(islot,iproc),status,mpierr)
       call MPI_TEST_CANCELLED(status,cancelled,mpierr)
       if (.not.cancelled) call fatal('finish_cell_exchange','cell received after the exchange finished')
    enddo
 enddo
 !$omp end critical (mpi_cell_exchange)
 !$omp end single

 !--no task may start the next exchange until all receives are cancelled
 call barrier_mpi
#endif

end subroutine finish_celldens_exchange

!----------------------------------------------------------------
!+
!  start the exchange of forces: post the receives
!  (called by all threads)
!+
!----------------------------------------------------------------
subroutine init_cellforce_exchange(stack)
 use mpiforce, only:stackforce,get_mpitype_of_cellforce
 type(stackforce), intent(inout) :: stack
#ifdef MPI
 integer :: iproc,islot,mpierr
#endif

 !$omp single
 if (dtype_force == -1) call get_mpitype_of_cellforce(dtype_force)
 call allocate_exchange_requests
 if (.not.allocated(recvbuf_force)) allocate(recvbuf_force(nrecvslots,nprocs))
 if (.not.allocated(sendbuf_force)) allocate(sendbuf_force(nsendbufs,size(sendreq,3)))
 if (size(sendbuf_force,2) /= size(sendreq,3)) then
    deallocate(sendbuf_force)
    allocate(sendbuf_force(nsendbufs,size(sendreq,3)))
 endif
#ifdef MPI
 !$omp critical (mpi_cell_exchange)
 do iproc=1,nprocs
    do islot=1,nrecvslots
       call MPI_IRECV(recvbuf_force(islot,iproc),1,dtype_force,iproc-1,MPI_ANY_TAG,&
                      comm_cellexchange,recvreq(islot,iproc),mpierr)
    enddo
 enddo
 !$omp end critical (mpi_cell_exchange)
#endif
 !$omp end single

end subroutine init_cellforce_exchange

!-----------------------------------------------------------------------
!+
!  send a cell to the tasks in targets, from the next buffer of this
!  thread, receiving cells into stack while waiting for that buffer
!+
!-----------------------------------------------------------------------
subroutine send_cellforce(cell,targets,stack,counters)
 use io,       only:fatal
 use omputils, only:omp_thread_num
 use mpiforce, only:stackforce
 type(cellforce),  intent(in)    :: cell
 logical,         intent(in)    :: targets(nprocs)
 type(stackforce), intent(inout) :: stack
 integer,         intent(inout) :: counters(nprocs,3)
#ifdef MPI
 integer :: newproc,mpierr,ith,ibuf,n

 ith  = omp_thread_num() + 1
 ibuf = mod(lastbuf(ith),nsendbufs) + 1
 lastbuf(ith) = ibuf

 !--wait for the previous sends from this buffer, continuing to receive
 call complete_sends(ibuf,ith,stack,counters)

 sendbuf_force(ibuf,ith) = cell
 n = 0
 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       n = n + 1
       call MPI_ISEND(sendbuf_force(ibuf,ith),1,dtype_force,newproc,0,comm_cellexchange,&
                      sendreq(n,ibuf,ith),mpierr)
       if (mpierr /= 0) call fatal('send_cellforce','error in MPI_ISEND')
       !$omp atomic
       counters(newproc+1,isent) = counters(newproc+1,isent) + 1
    endif
 enddo
 !$omp end critical (mpi_cell_exchange)
 nsendreq(ibuf,ith) = n
#endif

end subroutine send_cellforce

!-----------------------------------------------------------------------
!+
!  wait for the sends from one buffer of this thread to complete,
!  receiving cells into stack meanwhile
!+
!-----------------------------------------------------------------------
subroutine complete_sends_force(ibuf,ith,stack,counters)
 use mpiforce, only:stackforce
 integer,         intent(in)    :: ibuf,ith
 type(stackforce), intent(inout) :: stack
 integer,         intent(inout) :: counters(nprocs,3)

 do while (nsendreq(ibuf,ith) > 0)
    if (sends_complete(nsendreq(ibuf,ith),sendreq(:,ibuf,ith))) then
       nsendreq(ibuf,ith) = 0
    else
       call recv_cells(stack,counters)
    endif
 enddo

end subroutine complete_sends_force

!------------------------------------------------
!+
!  receive cells: if no other thread is already
!  doing so, test the receives, then process the
!  cells that have arrived and re-post the receives
!  (other threads may test the receives meanwhile)
!+
!------------------------------------------------
subroutine recv_cellforce(target_stack,counters)
 use io,        only:fatal
 use mpimemory, only:push_onto_stack,get_cell_ptr
 use mpiforce, only:stackforce
 type(stackforce), intent(inout) :: target_stack
 integer,         intent(inout) :: counters(nprocs,3)
#ifdef MPI
 type(cellforce),  pointer :: target_cell
 integer :: indices(nrecvslots*nprocs)
 integer :: nout,k,islot,iproc,mpierr

 nout = 0
 mpierr = 0
 if (claim_polling()) then
    !$omp critical (mpi_cell_exchange)
    call MPI_TESTSOME(nrecvslots*nprocs,recvreq,nout,indices,MPI_STATUSES_IGNORE,mpierr)
    !$omp end critical (mpi_cell_exchange)
 endif
 call release_polling
 if (mpierr /= 0) call fatal('recv_cells','error in MPI_TESTSOME')
 if (nout == MPI_UNDEFINED) nout = 0

 !--completed receives are inactive until re-posted, so no other thread will see them
 do k=1,nout
    islot = mod(indices(k)-1,nrecvslots) + 1
    iproc = (indices(k)-1)/nrecvslots + 1
    if (recvbuf_force(islot,iproc)%owner == id) then
       target_cell => get_cell_ptr(target_stack,recvbuf_force(islot,iproc)%waiting_index)
       call combine_cells(target_cell,recvbuf_force(islot,iproc))
    else
       call push_onto_stack(target_stack,recvbuf_force(islot,iproc))
    endif
    !$omp atomic
    counters(iproc,irecv) = counters(iproc,irecv) + 1
    !$omp critical (mpi_cell_exchange)
    call MPI_IRECV(recvbuf_force(islot,iproc),1,dtype_force,iproc-1,MPI_ANY_TAG,&
                   comm_cellexchange,recvreq(islot,iproc),mpierr)
    !$omp end critical (mpi_cell_exchange)
 enddo
#endif

end subroutine recv_cellforce

!-----------------------------------------------------------------------
!+
!  wait until all of the cells sent by this thread have gone,
!  and until all of the cells sent to this task have been received
!  (called by all threads at the end of each phase of the exchange)
!+
!-----------------------------------------------------------------------
subroutine recv_while_wait_force(stack,counters)
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_thread_num
 use timing,   only:get_timings,increment_timer,itimer_force_wait
 use mpiforce, only:stackforce
 type(stackforce), intent(inout) :: stack
 integer,         intent(inout) :: counters(nprocs,3)
#ifdef MPI
 real(kind=4) :: t1,t2,tcpu1,tcpu2
 integer      :: ith,ibuf

 !$omp masked
 call get_timings(t1,tcpu1)
 !$omp end masked

 !--complete the sends from this thread
 ith = omp_thread_num() + 1
 do ibuf=1,nsendbufs
    call complete_sends(ibuf,ith,stack,counters)
 enddo

 !--signal to other OMP threads that this thread has finished sending
 call set_thread_complete(.true.)

 !--continue receiving cells until all OMP threads have finished sending
 do while (.not. all_threads_complete())
    call recv_cells(stack,counters)
 enddo

 !--signal to other MPI tasks that this task has finished sending
 !$omp masked
 call send_cell_counts(counters)
 !$omp end masked

 !--continue receiving cells until all MPI tasks have finished sending
 do while (.not. all_tasks_complete())
    call recv_cells(stack,counters)
    !$omp masked
    call check_complete(counters)
    !$omp end masked
 enddo

 call barrier_mpi

 !$omp masked
 call finish_cell_counts
 ncomplete_mpi = 0
 ncells_recv = ncells_recv + sum(counters(:,irecv))
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_force_wait,t2-t1,tcpu2-tcpu1)
 !$omp end masked
 call set_thread_complete(.false.)
#endif

end subroutine recv_while_wait_force

!----------------------------------------------------------------
!+
!  finish the exchange: cancel the receives that are still posted
!  (called by all threads)
!+
!----------------------------------------------------------------
subroutine finish_cellforce_exchange(stack)
 use io,       only:fatal
 use mpiutils, only:barrier_mpi
 use mpiforce, only:stackforce
 type(stackforce), intent(inout) :: stack
#ifdef MPI
 integer :: iproc,islot,mpierr
 integer :: status(MPI_STATUS_SIZE)
 logical :: cancelled

 !$omp single
 !$omp critical (mpi_cell_exchange)
 do iproc=1,nprocs
    do islot=1,nrecvslots
       call MPI_CANCEL(recvreq(islot,iproc),mpierr)
       call MPI_WAIT(recvreq(islot,iproc),status,mpierr)
       call MPI_TEST_CANCELLED(status,cancelled,mpierr)
       if (.not.cancelled) call fatal('finish_cell_exchange','cell received after the exchange finished')
    enddo
 enddo
 !$omp end critical (mpi_cell_exchange)
 !$omp end single

 !--no task may start the next exchange until all receives are cancelled
 call barrier_mpi
#endif

end subroutine finish_cellforce_exchange

!-----------------------------------------------------------------------
!+
!  set the flag that this OMP thread has finished sending
!  (flushed, as other threads spin on it)
!+
!-----------------------------------------------------------------------
subroutine set_thread_complete(val)
 use omputils, only:omp_thread_num
 logical, intent(in) :: val

 !$omp flush
 !$omp atomic write
 thread_complete(omp_thread_num()+1) = val
 !$omp flush

end subroutine set_thread_complete

!-----------------------------------------------------------------------
!+
!  query whether all OMP threads on this task have finished sending
!+
!-----------------------------------------------------------------------
logical function all_threads_complete()
 logical :: done
 integer :: i

 all_threads_complete = .true.
 !$omp flush
 do i=1,size(thread_complete)
    !$omp atomic read
    done = thread_complete(i)
    if (.not.done) then
       all_threads_complete = .false.
       return
    endif
 enddo
 !$omp flush

end function all_threads_complete

!-----------------------------------------------------------------------
!+
!  query whether all cells sent by all MPI tasks have been received
!  (ncomplete_mpi is updated by the master thread in check_complete)
!+
!-----------------------------------------------------------------------
logical function all_tasks_complete()
 integer :: ncomplete

 !$omp flush
 !$omp atomic read
 ncomplete = ncomplete_mpi
 all_tasks_complete = (ncomplete >= nprocs)
 !$omp flush

end function all_tasks_complete

!-----------------------------------------------------------------------
!+
!  send the number of cells sent to each task, so the receiver knows
!  when it has received them all. The values are copied into a separate
!  buffer, since the counters are reset before the sends complete
!+
!-----------------------------------------------------------------------
subroutine send_cell_counts(counters)
 use io, only:fatal
 integer, intent(in) :: counters(nprocs,3)
#ifdef MPI
 integer :: newproc,mpierr

 countsent(:) = counters(:,isent)
 ncells_sent = ncells_sent + sum(countsent)
 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    if (newproc /= id) then
       call MPI_ISEND(countsent(newproc+1),1,MPI_INTEGER4,newproc,0,comm_cellcount,&
                      countsendrequest(newproc+1),mpierr)
       if (mpierr /= 0) call fatal('send_cell_counts','error in MPI_ISEND')
    endif
 enddo
 !$omp end critical (mpi_cell_exchange)
#endif

end subroutine send_cell_counts

!-----------------------------------------------------------------------
!+
!  complete the sends of the cell counts (called once all tasks have
!  received them, so this does not block)
!+
!-----------------------------------------------------------------------
subroutine finish_cell_counts
#ifdef MPI
 integer :: mpierr

 !$omp critical (mpi_cell_exchange)
 call MPI_WAITALL(nprocs,countsendrequest,MPI_STATUSES_IGNORE,mpierr)
 !$omp end critical (mpi_cell_exchange)
#endif

end subroutine finish_cell_counts

!------------------------------------------------
!+
!  Adds new_cell onto target_cell
!+
!------------------------------------------------
subroutine combine_celldens(target_cell,new_cell)
 use mpidens,   only:celldens
 use dim,       only:maxrhosum

 type(celldens), intent(inout) :: target_cell
 type(celldens), intent(inout) :: new_cell

 integer :: k, l

 do k = 1,new_cell%npcell
    do l = 1,maxrhosum
       !$omp atomic
       target_cell%rhosums(l,k) = target_cell%rhosums(l,k) + new_cell%rhosums(l,k)
    enddo
    !$omp atomic
    target_cell%nneigh(k) = target_cell%nneigh(k) + new_cell%nneigh(k)
 enddo
 !$omp atomic
 target_cell%nneightry = target_cell%nneightry + new_cell%nneightry

end subroutine combine_celldens

subroutine combine_cellforce(target_cell,new_cell)
 use mpiforce,  only:cellforce
 use dim,       only:maxfsum

 type(cellforce), intent(inout) :: target_cell
 type(cellforce), intent(inout) :: new_cell

 integer :: k, l

 do k = 1,new_cell%npcell
    do l = 1,maxfsum
       !$omp atomic
       target_cell%fsums(l,k) =      target_cell%fsums(l,k) +  new_cell%fsums(l,k)
    enddo
    !$omp atomic
    target_cell%tsmin(k)     = min( target_cell%tsmin(k),     new_cell%tsmin(k)     )
    !$omp atomic
    target_cell%vsigmax(k)   = max( target_cell%vsigmax(k),   new_cell%vsigmax(k)   )
#ifdef IND_TIMESTEPS
    !$omp atomic
    target_cell%ibinneigh(k) = max( target_cell%ibinneigh(k), new_cell%ibinneigh(k) )
#endif
 enddo

#ifdef GRAVITY
 do k = 1,20
    !$omp atomic
    target_cell%fgrav(k) = target_cell%fgrav(k) + new_cell%fgrav(k)
 enddo
#endif

 !$omp atomic
 target_cell%ndrag   = target_cell%ndrag   + new_cell%ndrag
 !$omp atomic
 target_cell%nstokes = target_cell%nstokes + new_cell%nstokes
 !$omp atomic
 target_cell%nsuper  = target_cell%nsuper  + new_cell%nsuper

end subroutine combine_cellforce

!----------------------------------------------------------------
!+
!  check which tasks have sent all of their cells to this one
!+
!----------------------------------------------------------------

subroutine check_complete(counters)
 use io, only:fatal
 integer, intent(inout) :: counters(nprocs,3)
#ifdef MPI
 integer :: i,ncomplete
 logical :: countreceived
 integer :: mpierr
 integer :: status(MPI_STATUS_SIZE)

 ncomplete = 1 ! self
 do i=1,nprocs
    if (i /= id + 1) then
       !$omp critical (mpi_cell_exchange)
       call MPI_TEST(countrequest(i),countreceived,status,mpierr)
       !$omp end critical (mpi_cell_exchange)
       if (countreceived) then
          if (counters(i,irecv) == counters(i,iexpect)) then
             ncomplete = ncomplete + 1
          elseif (counters(i,irecv) > counters(i,iexpect)) then
             print*,'on',id,'from',i-1
             print*,'nrecv',counters(i,irecv)
             print*,'nexpect',counters(i,iexpect)
             call fatal('mpiderivs', 'received more cells than expected')
          endif
       endif
    endif
 enddo
 ! other threads spin on this value
 !$omp atomic write
 ncomplete_mpi = ncomplete
 !$omp flush
#else
 ncomplete_mpi = 1
#endif
end subroutine check_complete

!----------------------------------------------------------------
!+
!  reset counters for checking arrival of all cells
!+
!----------------------------------------------------------------
subroutine reset_cell_counters(counters)
#ifdef MPI
 use io, only:fatal
#endif
 integer, intent(inout) :: counters(:,:)
#ifdef MPI
 integer :: iproc
 integer :: mpierr

 !$omp masked
 counters(:,isent)   = 0
 counters(:,iexpect) = -1
 counters(:,irecv)   = 0

 do iproc=1,nprocs
    if (iproc /= id + 1) then
       !$omp critical (mpi_cell_exchange)
       call MPI_IRECV(counters(iproc,iexpect),1,MPI_INTEGER4,iproc-1, &
       MPI_ANY_TAG,comm_cellcount,countrequest(iproc),mpierr)
       !$omp end critical (mpi_cell_exchange)
       if (mpierr /= 0) call fatal('reset_cell_counters','error in MPI_IRECV')
    endif
 enddo
 !$omp end masked

#endif
end subroutine reset_cell_counters

!----------------------------------------------------------------
!+
!  print the number of cells sent and received by each task
!  since the last call, and the size of the stacks (must be
!  called on all tasks)
!+
!----------------------------------------------------------------
subroutine print_cell_exchange_stats(iunit)
 use io,        only:master
 use mpiutils,  only:reduceall_mpi
 use mpimemory, only:dens_stack_1,dens_stack_2,dens_stack_3,force_stack_1,force_stack_2
 integer, intent(in) :: iunit
 integer(kind=8) :: nsent(3),nrecv(3)
 integer         :: maxstack

 call get_min_mean_max(ncells_sent,nsent)
 call get_min_mean_max(ncells_recv,nrecv)
 maxstack = max(dens_stack_1%maxlength,dens_stack_2%maxlength,dens_stack_3%maxlength,&
                force_stack_1%maxlength,force_stack_2%maxlength)
 maxstack = int(reduceall_mpi('max',maxstack))
 if (id==master) then
    write(iunit,"(1x,a,3(1x,i12))") 'MPI cells sent per task     (min,mean,max):',nsent
    write(iunit,"(1x,a,3(1x,i12))") 'MPI cells received per task (min,mean,max):',nrecv
    write(iunit,"(1x,a,i12,/)")     'MPI largest stack size (cells):',maxstack
 endif
 ncells_sent = 0
 ncells_recv = 0

end subroutine print_cell_exchange_stats

!----------------------------------------------------------------
!+
!  minimum, mean and maximum of a value over all tasks
!+
!----------------------------------------------------------------
subroutine get_min_mean_max(val,stats)
 use mpiutils, only:reduceall_mpi
 integer(kind=8), intent(in)  :: val
 integer(kind=8), intent(out) :: stats(3)

 stats(1) = reduceall_mpi('min',val)
 stats(2) = reduceall_mpi('+',val)/nprocs
 stats(3) = reduceall_mpi('max',val)

end subroutine get_min_mean_max

end module mpiderivs
