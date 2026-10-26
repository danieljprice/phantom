!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpiderivs
!
! This module handles the MPI exchange of information during the
!   density and force routines
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
!   mpiutils, omputils
!
#ifdef MPI
 use mpi
#endif
 use io,             only:id,nprocs
 use mpiutils,       only:comm_cellexchange,comm_cellcount

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

 interface combine_cells
  module procedure combine_celldens,combine_cellforce
 end interface combine_cells

 interface finish_cell_exchange
  module procedure finish_celldens_exchange,finish_cellforce_exchange
 end interface finish_cell_exchange

 interface recv_while_wait
  module procedure recv_while_wait_force,recv_while_wait_dens
 end interface recv_while_wait

 interface post_recv
  module procedure post_recv_dens,post_recv_force
 end interface post_recv

 public :: allocate_cell_comms_arrays
 public :: deallocate_cell_comms_arrays

 public :: init_cell_exchange
 public :: init_send_requests
 public :: send_cell
 public :: recv_cells
 public :: check_send_finished
 public :: finish_cell_exchange
 public :: recv_while_wait
 public :: reset_cell_counters
 public :: check_complete
 public :: combine_cells ! only to prevent compiler warning

 !
 !--the counters are module variables, but must be passed through as arguments
 !  in order to be threadsafe
 !
 integer, public, allocatable :: cell_counters(:,:)
 integer, parameter :: isent   = 1 ! counter for number of cells sent to i
 integer, parameter :: iexpect = 2 ! counter for number of cells expecting from i
 integer, parameter :: irecv   = 3 ! counter for number of cells received from i

 private

 integer, allocatable :: countrequest(:)     ! receives of the number of cells sent by each task
 integer, allocatable :: countsendrequest(:) ! sends of the number of cells sent to each task
 integer, allocatable :: countsent(:)        ! send buffer for the number of cells sent to each task
 integer :: dtype_dens  = -1                 ! MPI datatypes for celldens and cellforce,
 integer :: dtype_force = -1                 ! created once and reused

contains

subroutine allocate_cell_comms_arrays
 use allocutils, only:allocate_array
 use dim,        only:mpi
 if (mpi) then
    call allocate_array('cell_counters', cell_counters, nprocs, 3)
    call allocate_array('countrequest',  countrequest,  nprocs)
    call allocate_array('countsendrequest', countsendrequest, nprocs)
    call allocate_array('countsent', countsent, nprocs)
    call init_send_requests(countsendrequest)
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
end subroutine deallocate_cell_comms_arrays

!----------------------------------------------------------------
!+
!  initialise the receive type for each thread
!+
!----------------------------------------------------------------
subroutine init_celldens_exchange(xbufrecv,ireq,thread_complete,ncomplete_mpi,dtype)
 use mpidens,  only:celldens,get_mpitype_of_celldens
 use omputils, only:omp_num_threads

 type(celldens), intent(inout) :: xbufrecv(nprocs)
 integer,        intent(out)   :: ireq(nprocs)
 logical,        intent(inout) :: thread_complete(omp_num_threads)
 integer,        intent(out)   :: ncomplete_mpi
 integer,        intent(out)   :: dtype
#ifdef MPI
 integer                           :: iproc
#endif

 !$omp single
 if (dtype_dens == -1) call get_mpitype_of_celldens(dtype_dens)
 !$omp end single
 dtype = dtype_dens

#ifdef MPI

!
!--post a receive for EACH processor, to match the number of sends
!  (each is re-posted in recv_cells after it completes)
!
 do iproc=1,nprocs
    call post_recv(xbufrecv(iproc),dtype,iproc,ireq(iproc))
 enddo

 !$omp masked
 ncomplete_mpi = 0
 !$omp end masked
 call set_thread_complete(thread_complete,.false.)
#else
 ncomplete_mpi = 0
 ireq = 0
#endif

end subroutine init_celldens_exchange

subroutine init_cellforce_exchange(xbufrecv,ireq,thread_complete,ncomplete_mpi,dtype)
 use mpiforce, only:cellforce,get_mpitype_of_cellforce
 use omputils, only:omp_num_threads

 type(cellforce), intent(inout) :: xbufrecv(nprocs)
 integer,         intent(out)   :: ireq(nprocs)
 logical,         intent(inout) :: thread_complete(omp_num_threads)
 integer,         intent(out)   :: ncomplete_mpi
 integer,         intent(out)   :: dtype
#ifdef MPI
 integer                           :: iproc
#endif

 !$omp single
 if (dtype_force == -1) call get_mpitype_of_cellforce(dtype_force)
 !$omp end single
 dtype = dtype_force

#ifdef MPI

!
!--post a receive for EACH processor, to match the number of sends
!  (each is re-posted in recv_cells after it completes)
!
 do iproc=1,nprocs
    call post_recv(xbufrecv(iproc),dtype,iproc,ireq(iproc))
 enddo

 !$omp masked
 ncomplete_mpi = 0
 !$omp end masked
 call set_thread_complete(thread_complete,.false.)
#else
 ncomplete_mpi = 0
 ireq = 0
#endif
end subroutine init_cellforce_exchange

!-----------------------------------------------------------------------
!+
!  post a non-blocking receive for a cell from task iproc-1
!  (non-persistent: persistent receives started concurrently from
!   several threads are not thread safe in some MPI libraries)
!+
!-----------------------------------------------------------------------
subroutine post_recv_dens(xbuf,dtype,iproc,ireq)
 use io,      only:fatal
 use mpidens, only:celldens
 type(celldens), intent(inout) :: xbuf
 integer,        intent(in)    :: dtype,iproc
 integer,        intent(out)   :: ireq
#ifdef MPI
 integer :: mpierr

 !$omp critical (mpi_cell_exchange)
 call MPI_IRECV(xbuf,1,dtype,iproc-1,MPI_ANY_TAG,comm_cellexchange,ireq,mpierr)
 !$omp end critical (mpi_cell_exchange)
 if (mpierr /= 0) call fatal('post_recv','error in MPI_IRECV')
#else
 ireq = 0
#endif

end subroutine post_recv_dens

subroutine post_recv_force(xbuf,dtype,iproc,ireq)
 use io,       only:fatal
 use mpiforce, only:cellforce
 type(cellforce), intent(inout) :: xbuf
 integer,         intent(in)    :: dtype,iproc
 integer,         intent(out)   :: ireq
#ifdef MPI
 integer :: mpierr

 !$omp critical (mpi_cell_exchange)
 call MPI_IRECV(xbuf,1,dtype,iproc-1,MPI_ANY_TAG,comm_cellexchange,ireq,mpierr)
 !$omp end critical (mpi_cell_exchange)
 if (mpierr /= 0) call fatal('post_recv','error in MPI_IRECV')
#else
 ireq = 0
#endif

end subroutine post_recv_force

!-----------------------------------------------------------------------
!+
!  Subroutine to initialize send request handles to null
!+
!-----------------------------------------------------------------------
subroutine init_send_requests(irequestsend)
 integer, intent(out) :: irequestsend(nprocs)

!--instead of simply assigning 0, assign MPI_REQUEST_NULL
!  it is infact 0 in open mpi, but 0x2c000000 for MPICH-derived MPI.

#ifdef MPI
 irequestsend = MPI_REQUEST_NULL
#else
 irequestsend = 0
#endif

end subroutine init_send_requests

!-----------------------------------------------------------------------
!+
!  Subroutine to broadcast particle buffer to a bunch of processors
!+
!-----------------------------------------------------------------------
subroutine send_celldens(cell,targets,irequestsend,xsendbuf,counters,dtype)
 use io,       only:fatal
 use mpidens,  only:celldens

 type(celldens), intent(in)    :: cell
 logical,        intent(in)    :: targets(nprocs)
 integer,        intent(inout) :: irequestsend(nprocs)
 type(celldens), intent(out)   :: xsendbuf
 integer,        intent(inout) :: counters(nprocs,3)
 integer,        intent(in)    :: dtype
#ifdef MPI
 integer                            :: newproc,mpierr

 xsendbuf = cell
 irequestsend = MPI_REQUEST_NULL

 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       call MPI_ISEND(xsendbuf,1,dtype,newproc,1,comm_cellexchange,irequestsend(newproc+1),mpierr)
       if (mpierr /= 0) call fatal('send_celldens','error in MPI_ISEND')
       !$omp atomic
       counters(newproc+1,isent) = counters(newproc+1,isent) + 1
    endif
 enddo
 !$omp end critical (mpi_cell_exchange)
#else
 xsendbuf = cell
#endif

end subroutine send_celldens

subroutine send_cellforce(cell,targets,irequestsend,xsendbuf,counters,dtype)
 use io,       only:fatal
 use mpiforce, only:cellforce

 type(cellforce), intent(in)    :: cell
 logical,         intent(in)    :: targets(nprocs)
 integer,         intent(inout) :: irequestsend(nprocs)
 type(cellforce), intent(out)   :: xsendbuf
 integer,         intent(inout) :: counters(nprocs,3)
 integer,         intent(in)    :: dtype
#ifdef MPI
 integer                            :: newproc,mpierr

 xsendbuf = cell
 irequestsend = MPI_REQUEST_NULL

 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       call MPI_ISEND(xsendbuf,1,dtype,newproc,0,comm_cellexchange,irequestsend(newproc+1),mpierr)
       if (mpierr /= 0) call fatal('send_cellforce','error in MPI_ISEND')
       !$omp atomic
       counters(newproc+1,isent) = counters(newproc+1,isent) + 1
    endif
 enddo
 !$omp end critical (mpi_cell_exchange)
#else
 xsendbuf = cell
#endif

end subroutine send_cellforce

!-----------------------------------------------------------------------
!+
!  Subroutine to check that non-blocking send has completed
!+
!-----------------------------------------------------------------------
subroutine check_send_finished(irequestsend,idone)
 integer, intent(inout) :: irequestsend(nprocs)
 logical, intent(out)   :: idone(nprocs)

#ifdef MPI
 integer :: newproc
 integer :: mpierr
 integer :: status(MPI_STATUS_SIZE)

 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    if (newproc /= id) call MPI_TEST(irequestsend(newproc+1),idone(newproc+1),status,mpierr)
 enddo
 !$omp end critical (mpi_cell_exchange)
 !--never test self; always set to true
 idone(id+1) = .true.
#else
 idone = .true.
#endif

end subroutine check_send_finished

subroutine recv_while_wait_dens(stack,xrecvbuf,irequestrecv,thread_complete,counters,ncomplete_mpi)
 use mpidens,  only:stackdens,celldens
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_num_threads
 type(stackdens), intent(inout) :: stack
 type(celldens),  intent(inout) :: xrecvbuf(nprocs)
 integer,         intent(inout) :: irequestrecv(nprocs)
 logical,         intent(inout) :: thread_complete(omp_num_threads)
 integer,         intent(inout) :: counters(nprocs,3)
 integer,         intent(inout) :: ncomplete_mpi
#ifdef MPI

 !--signal to other OMP threads that this thread has finished sending
 call set_thread_complete(thread_complete,.true.)

 !--continue receiving cells until all OMP threads have finished sending
 do while (.not. all_threads_complete(thread_complete))
    call recv_cells(stack,xrecvbuf,irequestrecv,counters)
 enddo

 !--signal to other MPI tasks that this task has finished sending
 !$omp masked
 call send_cell_counts(counters)
 !$omp end masked

 !--continue receiving cells until all MPI tasks have finished sending
 do while (.not. all_tasks_complete(ncomplete_mpi))
    call recv_cells(stack,xrecvbuf,irequestrecv,counters)
    !$omp masked
    call check_complete(counters,ncomplete_mpi)
    !$omp end masked
 enddo

 call barrier_mpi

 !$omp masked
 call finish_cell_counts
 ncomplete_mpi = 0
 !$omp end masked
 call set_thread_complete(thread_complete,.false.)

#endif

end subroutine recv_while_wait_dens

subroutine recv_while_wait_force(stack,xrecvbuf,irequestrecv,thread_complete,counters,ncomplete_mpi)
 use mpiforce, only:stackforce,cellforce
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_num_threads
 type(stackforce), intent(inout) :: stack
 type(cellforce),  intent(inout) :: xrecvbuf(nprocs)
 integer,          intent(inout) :: irequestrecv(nprocs)
 logical,          intent(inout) :: thread_complete(omp_num_threads)
 integer,          intent(inout) :: counters(nprocs,3)
 integer,          intent(inout) :: ncomplete_mpi
#ifdef MPI

 !--signal to other OMP threads that this thread has finished sending
 call set_thread_complete(thread_complete,.true.)

 !--continue receiving cells until all OMP threads have finished sending
 do while (.not. all_threads_complete(thread_complete))
    call recv_cells(stack,xrecvbuf,irequestrecv,counters)
 enddo

 !--signal to other MPI tasks that this task has finished sending
 !$omp masked
 call send_cell_counts(counters)
 !$omp end masked

 !--continue receiving cells until all MPI tasks have finished sending
 do while (.not. all_tasks_complete(ncomplete_mpi))
    call recv_cells(stack,xrecvbuf,irequestrecv,counters)
    !$omp masked
    call check_complete(counters,ncomplete_mpi)
    !$omp end masked
 enddo

 call barrier_mpi

 !$omp masked
 call finish_cell_counts
 ncomplete_mpi = 0
 !$omp end masked
 call set_thread_complete(thread_complete,.false.)

#endif

end subroutine recv_while_wait_force

!-----------------------------------------------------------------------
!+
!  set the flag that this OMP thread has finished sending
!  (flushed, as other threads spin on it)
!+
!-----------------------------------------------------------------------
subroutine set_thread_complete(thread_complete,val)
 use omputils, only:omp_num_threads,omp_thread_num
 logical, intent(inout) :: thread_complete(omp_num_threads)
 logical, intent(in)    :: val

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
logical function all_threads_complete(thread_complete)
 use omputils, only:omp_num_threads
 logical, intent(in) :: thread_complete(omp_num_threads)
 logical :: done
 integer :: i

 all_threads_complete = .true.
 !$omp flush
 do i=1,omp_num_threads
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
logical function all_tasks_complete(ncomplete_mpi)
 integer, intent(in) :: ncomplete_mpi
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
!  Receive cells
!+
!------------------------------------------------
subroutine recv_celldens(target_stack,xbuf,irequestrecv,counters)
 use io,        only:fatal
 use mpimemory, only:push_onto_stack,get_cell_ptr
 use mpidens,   only:stackdens,celldens

 type(celldens),  intent(inout) :: xbuf(:)  ! just need memory address
 type(stackdens), intent(inout) :: target_stack
 integer,         intent(inout) :: irequestrecv(nprocs)
 integer,         intent(inout) :: counters(nprocs,3)
#ifdef MPI
 type(celldens),  pointer           :: target_cell
 integer                            :: iproc,iwait
 logical                            :: igot
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)

 ! receive MPI broadcast
 do iproc=1,nprocs
    !$omp critical (mpi_cell_exchange)
    call MPI_TEST(irequestrecv(iproc),igot,status,mpierr)
    !$omp end critical (mpi_cell_exchange)
    if (mpierr /= 0) call fatal('recv_cell','error in MPI_TEST call')
    ! unpack results
    if (igot) then
       if (xbuf(iproc)%owner == id) then
          iwait = xbuf(iproc)%waiting_index
          target_cell => get_cell_ptr(target_stack,iwait)
          call combine_cells(target_cell,xbuf(iproc))
       else
          call push_onto_stack(target_stack, xbuf(iproc))
       endif
       !$omp atomic
       counters(iproc,irecv) = counters(iproc,irecv) + 1
       call post_recv(xbuf(iproc),dtype_dens,iproc,irequestrecv(iproc))
    endif
 enddo
#endif
end subroutine recv_celldens

subroutine recv_cellforce(target_stack,xbuf,irequestrecv,counters)
 use io,        only:fatal
 use mpimemory, only:push_onto_stack,get_cell_ptr
 use mpiforce,  only:stackforce,cellforce

 type(cellforce),  intent(inout) :: xbuf(:)  ! just need memory address
 type(stackforce), intent(inout) :: target_stack
 integer,          intent(inout) :: irequestrecv(nprocs)
 integer,          intent(inout) :: counters(nprocs,3)
#ifdef MPI
 type(cellforce),  pointer          :: target_cell
 integer                            :: iproc,iwait
 logical                            :: igot
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)

 ! receive MPI broadcast
 do iproc=1,nprocs
    !$omp critical (mpi_cell_exchange)
    call MPI_TEST(irequestrecv(iproc),igot,status,mpierr)
    !$omp end critical (mpi_cell_exchange)
    if (mpierr /= 0) call fatal('recv_cell','error in MPI_TEST call')
    ! unpack results
    if (igot) then
       if (xbuf(iproc)%owner == id) then
          iwait = xbuf(iproc)%waiting_index
          target_cell => get_cell_ptr(target_stack,iwait)
          call combine_cells(target_cell,xbuf(iproc))
       else
          call push_onto_stack(target_stack, xbuf(iproc))
       endif
       !$omp atomic
       counters(iproc,irecv) = counters(iproc,irecv) + 1
       call post_recv(xbuf(iproc),dtype_force,iproc,irequestrecv(iproc))
    endif
 enddo
#endif
end subroutine recv_cellforce

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
!  finish/clean up of load balancing process.
!+
!----------------------------------------------------------------
subroutine finish_celldens_exchange(irequestrecv,xsendbuf,dtype)
 use io,       only:fatal
 use mpidens,  only:celldens
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_thread_num
 integer,        intent(inout) :: irequestrecv(nprocs)
 type(celldens), intent(in)    :: xsendbuf
 integer,        intent(inout) :: dtype
#ifdef MPI
 integer                            :: newproc,iproc
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)

!
!--each thread sends a dummy message to every task to clear the last remaining receives
!
 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    call MPI_SEND(xsendbuf,1,dtype,newproc,1,comm_cellexchange,mpierr)
 enddo
 !$omp end critical (mpi_cell_exchange)

!
!--sync all threads here
!
 call barrier_mpi

!--complete the outstanding receives (this also frees them)
 !$omp critical (mpi_cell_exchange)
 do iproc=1,nprocs
    call MPI_WAIT(irequestrecv(iproc),status,mpierr)
 enddo
 !$omp end critical (mpi_cell_exchange)

#endif
end subroutine finish_celldens_exchange

subroutine finish_cellforce_exchange(irequestrecv,xsendbuf,dtype)
 use io,       only:fatal
 use mpiforce, only:cellforce
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_thread_num
 integer,         intent(inout) :: irequestrecv(nprocs)
 type(cellforce), intent(in)    :: xsendbuf
 integer,         intent(inout) :: dtype
#ifdef MPI
 integer                            :: newproc,iproc
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)

!
!--each thread sends a dummy message to every task to clear the last remaining receives
!
 !$omp critical (mpi_cell_exchange)
 do newproc=0,nprocs-1
    call MPI_SEND(xsendbuf,1,dtype,newproc,0,comm_cellexchange,mpierr)
 enddo
 !$omp end critical (mpi_cell_exchange)

!
!--sync all threads here
!
 call barrier_mpi

!--complete the outstanding receives (this also frees them)
 !$omp critical (mpi_cell_exchange)
 do iproc=1,nprocs
    call MPI_WAIT(irequestrecv(iproc),status,mpierr)
 enddo
 !$omp end critical (mpi_cell_exchange)

#endif
end subroutine finish_cellforce_exchange

!----------------------------------------------------------------
!+
!  check which threads have completed
!+
!----------------------------------------------------------------

subroutine check_complete(counters,ncomplete_mpi)
 use io, only:fatal
 integer, intent(inout) :: counters(nprocs,3)
 integer, intent(inout) :: ncomplete_mpi
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

end module mpiderivs
