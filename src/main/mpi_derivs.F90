!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module mpiderivs
!
! This module handles the MPI exchange of information during the
!   density and force routines
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, dim, dtypekdtree, io, mpi, mpidens, mpiforce,
!   mpimemory, mpiutils, omputils
!
#ifdef MPI
 use mpi
#endif
 use io,             only:id,nprocs
 use mpiutils,       only:comm_cellexchange,comm_cellcount
 use dtypekdtree,    only:kdnode,ndimtree

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

 interface finish_cell_exchange
  module procedure finish_celldens_exchange,finish_cellforce_exchange
 end interface finish_cell_exchange

 interface recv_while_wait
  module procedure recv_while_wait_force,recv_while_wait_dens
 end interface recv_while_wait

 interface reduce_group
  module procedure reduce_group_real, reduce_group_int
 end interface reduce_group

 public :: allocate_comms_arrays
 public :: deallocate_comms_arrays

 public :: init_cell_exchange
 public :: send_cell
 public :: recv_cells
 public :: check_send_finished
 public :: finish_cell_exchange
 public :: recv_while_wait
 public :: get_group_cofm
 public :: reduce_group

 public :: tree_sync
 public :: tree_bcast
 public :: finish_tree_comms
 public :: reset_cell_counters
 public :: check_complete

 !
 !--the counters are module variables, but must be passed through as arguments
 !  in order to be threadsafe
 !
 integer, public, allocatable :: cell_counters(:,:)
 integer, parameter :: isent   = 1 ! counter for number of cells sent to i
 integer, parameter :: iexpect = 2 ! counter for number of cells expecting from i
 integer, parameter :: irecv   = 3 ! counter for number of cells received from i

 private

#ifdef MPI
 integer :: ncomplete

 integer :: dtype_celldens

 integer :: globallevel
#endif
 integer,allocatable :: comm_cofm(:)  ! only comms up to globallevel are used
 integer,allocatable :: comm_owner(:) ! only comms up to globallevel are used

 integer,allocatable :: countrequest(:)

contains

subroutine allocate_comms_arrays
 use allocutils, only:allocate_array
 use dim,        only:mpi
 if (mpi) then
    call allocate_array('cell_counters', cell_counters, nprocs, 3)
    call allocate_array('countrequest',  countrequest,  nprocs)
    call allocate_array('comm_cofm',     comm_cofm,     nprocs)
    call allocate_array('comm_owner',    comm_owner,    nprocs)
    call init_tree_comms
 else
    ! dummy cell counters that are required to prevent runtime errors
    ! in dens and force
    call allocate_array('cell_counters', cell_counters, 0, 0)
 endif
end subroutine allocate_comms_arrays

subroutine deallocate_comms_arrays
 if (allocated(cell_counters)) deallocate(cell_counters)
 if (allocated(countrequest )) deallocate(countrequest )
 if (allocated(comm_cofm    )) deallocate(comm_cofm    )
 if (allocated(comm_owner   )) deallocate(comm_owner   )
end subroutine deallocate_comms_arrays

!----------------------------------------------------------------
!+
!  initialise the receive type for each thread
!+
!----------------------------------------------------------------
subroutine init_celldens_exchange(xbufrecv,ireq)
 use io,       only:fatal
 use mpidens,  only:get_mpitype_of_celldens,celldens

 type(celldens),     intent(inout) :: xbufrecv(nprocs)
 integer,            intent(out)   :: ireq(nprocs) !,nrecv

#ifdef MPI
 integer                           :: iproc
 integer                           :: mpierr
!
!--use persistent communication type for receives
!  cannot do same for sends as there are different destinations,
!  unless we make a request for each processor
!
!  We post a receive for EACH processor, to match the number of sends
!

 call get_mpitype_of_celldens(dtype_celldens)

 do iproc=1,nprocs
    call MPI_RECV_INIT(xbufrecv(iproc),1,dtype_celldens,iproc-1, &
                       MPI_ANY_TAG,comm_cellexchange,ireq(iproc),mpierr)
    if (mpierr /= 0) call fatal('init_cell_exchange','error in MPI_RECV_INIT')
!
!--start the persistent communication channel
!
    call MPI_START(ireq(iproc),mpierr)
    if (mpierr /= 0) call fatal('init_cell_exchange','error in MPI_START')
 enddo

 ncomplete = 0
#endif
end subroutine init_celldens_exchange

subroutine init_cellforce_exchange(xbufrecv,ireq,thread_complete,ncomplete_mpi,any_tag)
 use io,       only:fatal
 use mpiforce, only:dtype_cellforce,cellforce
 use omputils, only:omp_thread_num,omp_num_threads

 type(cellforce),    intent(inout) :: xbufrecv(nprocs)
 integer,            intent(out)   :: ireq(nprocs) !,nrecv
 logical,            intent(inout) :: thread_complete(omp_num_threads)
 integer,            intent(out)   :: ncomplete_mpi
 logical,            intent(in)    :: any_tag
#ifdef MPI
 integer :: iproc, mpierr, tag

!
!--use persistent communication type for receives
!  cannot do same for sends as there are different destinations,
!  unless we make a request for each processor
!
!  We post a receive for EACH processor, to match the number of sends
!
!  Optionally, we can also force each reciever to only accept messages with
!  tags that match the omp thread id. By default they accept MPI messages with
!  any tag.

 if (any_tag) then
    tag = MPI_ANY_TAG
 else
    tag = omp_thread_num()
 endif

 do iproc=1,nprocs
    call MPI_RECV_INIT(xbufrecv(iproc),1,dtype_cellforce,iproc-1, &
                       tag,comm_cellexchange,ireq(iproc),mpierr)
    if (mpierr /= 0) call fatal('init_cell_exchange','error in MPI_RECV_INIT')
!
!--start the persistent communication channel
!
    call MPI_START(ireq(iproc),mpierr)
    if (mpierr /= 0) call fatal('init_cell_exchange','error in MPI_START')
 enddo

 !$omp master
 ncomplete_mpi = 0
 !$omp end master
 thread_complete(omp_thread_num()+1) = .false.
#endif
end subroutine init_cellforce_exchange

!-----------------------------------------------------------------------
!+
!  Subroutine to broadcast particle buffer to a bunch of processors
!+
!-----------------------------------------------------------------------
subroutine send_celldens(cell,targets,irequestsend,xsendbuf,counters)
 use io,       only:fatal
 use mpidens,  only:celldens

 type(celldens),     intent(in)     :: cell
 logical,            intent(in)     :: targets(nprocs)
 integer,            intent(inout)  :: irequestsend(nprocs)
 type(celldens),     intent(out)    :: xsendbuf
 integer,            intent(inout)  :: counters(nprocs,3)
#ifdef MPI
 integer                            :: newproc
 integer                            :: mpierr

 xsendbuf = cell
 irequestsend = MPI_REQUEST_NULL

 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       call MPI_ISEND(xsendbuf,1,dtype_celldens,newproc,1,comm_cellexchange,irequestsend(newproc+1),mpierr)
       if (mpierr /= 0) call fatal('send_celldens','error in MPI_ISEND')
       !$omp atomic
       counters(newproc+1,isent) = counters(newproc+1,isent) + 1
    endif
 enddo
#endif

end subroutine send_celldens

subroutine send_cellforce(cell,targets,irequestsend,xsendbuf,counters)
 use io,       only:fatal
 use mpiforce, only:cellforce,dtype_cellforce

 type(cellforce),    intent(in)     :: cell
 logical,            intent(in)     :: targets(nprocs)
 integer,            intent(inout)  :: irequestsend(nprocs)
 type(cellforce),    intent(out)    :: xsendbuf
 integer,            intent(inout)  :: counters(nprocs,3)
#ifdef MPI
 integer                            :: newproc
 integer                            :: mpierr

 xsendbuf = cell
 irequestsend = MPI_REQUEST_NULL

 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       call MPI_ISEND(xsendbuf,1,dtype_cellforce,newproc,cell%owner_thread,comm_cellexchange,irequestsend(newproc+1),mpierr)
       if (mpierr /= 0) call fatal('send_cellforce','error in MPI_ISEND')
       !$omp atomic
       counters(newproc+1,isent) = counters(newproc+1,isent) + 1
    endif
 enddo
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

 do newproc=0,nprocs-1
    if (newproc /= id) call MPI_TEST(irequestsend(newproc+1),idone(newproc+1),status,mpierr)
 enddo
 !--never test self; always set to true
 idone(id+1) = .true.
#endif

end subroutine check_send_finished

subroutine recv_while_wait_dens(stack,xrecvbuf,irequestrecv,irequestsend)
 use mpidens,  only:stackdens,celldens
 type(stackdens),  intent(inout) :: stack
 type(celldens),   intent(inout) :: xrecvbuf(nprocs)
 integer,          intent(inout) :: irequestrecv(nprocs),irequestsend(nprocs)

#ifdef MPI
 integer             :: newproc
 integer             :: mpierr

 do newproc=0,nprocs-1
    if (newproc /= id) then
       !--tag=0 to signal done
       call MPI_ISEND(cell_counters(newproc+1,isent),1,MPI_INTEGER4,newproc,0,comm_cellcount,irequestsend(newproc+1),mpierr)
    endif
 enddo

 !--do not need to MPI_WAIT, because the following code requires the sends to go through
 do while (ncomplete < nprocs)
    call recv_celldens(stack,xrecvbuf,irequestrecv)
    call check_complete(cell_counters)
 enddo

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 !--reset counter for next round
 ncomplete = 0
#endif

end subroutine recv_while_wait_dens

subroutine recv_while_wait_force(stack,xrecvbuf,irequestrecv,irequestsend,thread_complete,counters,ncomplete_mpi)
 use mpiforce, only:stackforce,cellforce
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_num_threads,omp_thread_num
 type(stackforce), intent(inout) :: stack
 type(cellforce),  intent(inout) :: xrecvbuf(nprocs)
 integer,          intent(inout) :: irequestrecv(nprocs),irequestsend(nprocs)
 logical,          intent(inout) :: thread_complete(omp_num_threads)
 integer,          intent(inout) :: counters(nprocs,3)
 integer,          intent(inout) :: ncomplete_mpi
#ifdef MPI
 integer             :: newproc
 integer             :: mpierr

 !--signal to other OMP threads that this thread has finished sending
 thread_complete(omp_thread_num()+1) = .true.

!--continue receiving cells until all OMP threads have finished sending
 do while (.not. all(thread_complete))
    call recv_cellforce(stack,xrecvbuf,irequestrecv,counters)
 enddo

 !--signal to other MPI tasks that this task has finished sending
 !$omp master
 do newproc=0,nprocs-1
    if (newproc /= id) then
       call MPI_ISEND(counters(newproc+1,isent),1,MPI_INTEGER4,newproc,0,comm_cellcount,irequestsend(newproc+1),mpierr)
    endif
 enddo
 !$omp end master

 !--continue receiving cells until all MPI tasks have finished sending
 do while (ncomplete_mpi < nprocs)
    call recv_cellforce(stack,xrecvbuf,irequestrecv,counters)
    !$omp master
    call check_complete_force(counters,ncomplete_mpi)
    !$omp end master
 enddo

 call barrier_mpi

#endif

end subroutine recv_while_wait_force

!------------------------------------------------
!+
!  Receive cells
!+
!------------------------------------------------
subroutine recv_celldens(target_stack,xbuf,irequestrecv)
 use io,        only:fatal
 use mpimemory, only:push_onto_stack
 use mpidens,   only:stackdens,celldens

 type(celldens),     intent(inout)  :: xbuf(:)  ! just need memory address
 type(stackdens),    intent(inout)  :: target_stack
 integer,            intent(inout)  :: irequestrecv(nprocs)
#ifdef MPI
 integer                            :: iproc,k,iwait
 logical                            :: igot
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)

 ! receive MPI broadcast
 do iproc=1,nprocs
    call MPI_TEST(irequestrecv(iproc),igot,status,mpierr)
    if (mpierr /= 0) call fatal('recv_cell','error in MPI_TEST call')
    ! unpack results
    if (igot) then
       if (xbuf(iproc)%owner == id) then
          iwait = xbuf(iproc)%waiting_index
          do k = 1,xbuf(iproc)%npcell
             target_stack%cells(iwait)%rhosums(:,k) = target_stack%cells(iwait)%rhosums(:,k) + xbuf(iproc)%rhosums(:,k)
             target_stack%cells(iwait)%nneigh(k) = target_stack%cells(iwait)%nneigh(k) + xbuf(iproc)%nneigh(k)
          enddo
          target_stack%cells(iwait)%nneightry = target_stack%cells(iwait)%nneightry + xbuf(iproc)%nneightry
          cell_counters(iproc,irecv) = cell_counters(iproc,irecv) + 1
       else
          call push_onto_stack(target_stack, xbuf(iproc))
          cell_counters(iproc,irecv) = cell_counters(iproc,irecv) + 1
       endif
       call MPI_START(irequestrecv(iproc),mpierr)
    endif
 enddo
#endif
end subroutine recv_celldens

subroutine recv_cellforce(target_stack,xbuf,irequestrecv,counters)
 use io,        only:fatal
 use mpimemory, only:push_onto_stack
 use mpiforce,  only:stackforce,cellforce

 type(cellforce),    intent(inout)  :: xbuf(:)  ! just need memory address
 type(stackforce),   intent(inout)  :: target_stack
 integer,            intent(inout)  :: irequestrecv(nprocs)
 integer,            intent(inout)  :: counters(nprocs,3)
#ifdef MPI
 integer                            :: iproc,k,iwait
 logical                            :: igot
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)

 ! receive MPI broadcast
 do iproc=1,nprocs
    call MPI_TEST(irequestrecv(iproc),igot,status,mpierr)
    if (mpierr /= 0) call fatal('recv_cell','error in MPI_TEST call')
    ! unpack results
    if (igot) then
       if (xbuf(iproc)%owner == id) then
          iwait = xbuf(iproc)%waiting_index
          do k = 1,xbuf(iproc)%npcell
             target_stack%cells(iwait)%fsums(:,k) = target_stack%cells(iwait)%fsums(:,k) + xbuf(iproc)%fsums(:,k)
             target_stack%cells(iwait)%tsmin(k) = min(target_stack%cells(iwait)%tsmin(k), xbuf(iproc)%tsmin(k))
             target_stack%cells(iwait)%vsigmax(k) = max(target_stack%cells(iwait)%vsigmax(k), xbuf(iproc)%vsigmax(k))
#ifdef IND_TIMESTEPS
             target_stack%cells(iwait)%ibinneigh(k) = max(target_stack%cells(iwait)%ibinneigh(k), xbuf(iproc)%ibinneigh(k))
#endif
          enddo
#ifdef GRAVITY
          do k = 1,20
             target_stack%cells(iwait)%fgrav(k) = target_stack%cells(iwait)%fgrav(k) + xbuf(iproc)%fgrav(k)
          enddo
#endif
          target_stack%cells(iwait)%ndrag = target_stack%cells(iwait)%ndrag + xbuf(iproc)%ndrag
          target_stack%cells(iwait)%nstokes = target_stack%cells(iwait)%nstokes + xbuf(iproc)%nstokes
          target_stack%cells(iwait)%nsuper = target_stack%cells(iwait)%nsuper + xbuf(iproc)%nsuper
       else
          call push_onto_stack(target_stack, xbuf(iproc))
       endif
       !$omp atomic
       counters(iproc,irecv) = counters(iproc,irecv) + 1
       call MPI_START(irequestrecv(iproc),mpierr)
    endif
 enddo
#endif
end subroutine recv_cellforce

!----------------------------------------------------------------
!+
!  finish/clean up of load balancing process.
!+
!----------------------------------------------------------------
subroutine finish_celldens_exchange(irequestrecv,xsendbuf)
 use io,       only:fatal
 use mpidens,  only:celldens
 integer,            intent(inout)  :: irequestrecv(nprocs)
 type(celldens), intent(in)         :: xsendbuf
#ifdef MPI
 integer                            :: newproc,iproc
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)
!
!--each processor do a dummy send to next processor to clear the last remaining receive
!  (we know the receive has been posted for this, so use RSEND)
!
 do newproc=0,nprocs-1
    call MPI_RSEND(xsendbuf,1,dtype_celldens,newproc,1,comm_cellexchange,mpierr)
 enddo

!
!--sync all threads here
!
 call MPI_BARRIER(comm_cellexchange,mpierr)

!--free request handle
 do iproc=1,nprocs
    call MPI_WAIT(irequestrecv(iproc),status,mpierr)
    call MPI_REQUEST_FREE(irequestrecv(iproc),mpierr)
 enddo
#endif
end subroutine finish_celldens_exchange

subroutine finish_cellforce_exchange(irequestrecv,xsendbuf)
 use io,       only:fatal
 use mpiforce, only:cellforce,dtype_cellforce
 use mpiutils, only:barrier_mpi
 use omputils, only:omp_thread_num
 integer,            intent(inout)  :: irequestrecv(nprocs)
 type(cellforce), intent(in)        :: xsendbuf
#ifdef MPI
 integer                            :: newproc,iproc
 integer                            :: mpierr
 integer                            :: status(MPI_STATUS_SIZE)
!
!--each processor do a dummy send to next processor to clear the last remaining receive
!  (we know the receive has been posted for this, so use RSEND)
!
 do newproc=0,nprocs-1
    call MPI_RSEND(xsendbuf,1,dtype_cellforce,newproc,omp_thread_num(),comm_cellexchange,mpierr)
 enddo

!
!--sync all threads here
!
 call barrier_mpi

!--free request handle
 do iproc=1,nprocs
    call MPI_WAIT(irequestrecv(iproc),status,mpierr)
    call MPI_REQUEST_FREE(irequestrecv(iproc),mpierr)
 enddo
#endif
end subroutine finish_cellforce_exchange

!----------------------------------------------------------------
!+
!  initialise communicators for tree construction
!+
!----------------------------------------------------------------
subroutine init_tree_comms()
#ifdef MPI
 integer :: level,groupsize,color
 integer                            :: mpierr

 globallevel = int(ceiling(log(real(nprocs)) / log(2.0)))

 ! cofm group is all procs at level 0
 call MPI_COMM_DUP(MPI_COMM_WORLD, comm_cofm(1), mpierr)
 ! owner comm is never needed at level 0, left uninitialised

 do level = 1, globallevel
    groupsize = 2**(globallevel - level)

    ! cofm group
    color = id / groupsize
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, id, comm_cofm(level+1), mpierr)

    ! owner group
    if (groupsize == 1) then
       call MPI_COMM_DUP(MPI_COMM_WORLD, comm_owner(level+1), mpierr)
    else
       if (mod(id, groupsize) == 0) then
          color = 0
       else
          color = 1
       endif
       call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, id, comm_owner(level+1), mpierr)
    endif
 enddo
#endif
end subroutine init_tree_comms

subroutine finish_tree_comms()
#ifdef MPI
 integer :: level
 integer :: mpierr

 do level = 0, globallevel
    call MPI_COMM_FREE(comm_cofm(level+1), mpierr)
 enddo
 do level = 1, globallevel
    call MPI_COMM_FREE(comm_owner(level+1), mpierr)
 enddo
#endif
end subroutine finish_tree_comms

!----------------------------------------------------------------
!+
!  get the COFM of the group
!+
!----------------------------------------------------------------
subroutine get_group_cofm(xyzcofm,totmass_node,level,cofmsum,totmassg)
 real,      intent(in)        :: xyzcofm(3)
 real,      intent(in)        :: totmass_node
 integer,   intent(in)        :: level

 real,      intent(out)       :: cofmsum(3)
 real,      intent(out)       :: totmassg

#ifdef MPI
 real                         :: cofmpart(3)
 integer                      :: mpierr

 cofmpart = xyzcofm * totmass_node
 call MPI_ALLREDUCE(totmass_node,totmassg,1,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 call MPI_ALLREDUCE(cofmpart,cofmsum,3,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 cofmsum = cofmsum / totmassg
#endif

end subroutine get_group_cofm

!----------------------------------------------------------------
!+
!  tree group reductions
!+
!----------------------------------------------------------------

function reduce_group_real(x,string,level) result(xg)
 use io, only:fatal
 real,               intent(in)        :: x
 character(len=*),   intent(in)        :: string
 integer,            intent(in)        :: level
 real                                  :: xg

#ifdef MPI
 real                                  :: isend, ired
 integer                               :: mpierr

 isend = x

 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,1,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,1,MPI_REAL8,MPI_MAX,comm_cofm(level+1),mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,1,MPI_REAL8,MPI_MIN,comm_cofm(level+1),mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select

 xg = ired
#else
 xg = x
#endif

end function reduce_group_real

function reduce_group_int(x,string,level) result(xg)
 use io, only:fatal
 integer,            intent(in)        :: x
 character(len=*),   intent(in)        :: string
 integer,            intent(in)        :: level
 integer                               :: xg

#ifdef MPI
 integer                               :: isend, ired
 integer                               :: mpierr

 isend = x

 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_SUM,comm_cofm(level+1),mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_MAX,comm_cofm(level+1),mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_MIN,comm_cofm(level+1),mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select

 xg = ired
#else
 xg = x
#endif
end function reduce_group_int

!----------------------------------------------------------------
!+
!  synchronize the global tree, placing nodes in the correct position
!+
!----------------------------------------------------------------
subroutine tree_sync(node_in,n_in,node_synced,n_synced,ifirstingroup,level)
 use dtypekdtree, only:get_mpitype_of_kdnode

 integer, intent(in)         :: ifirstingroup,level
 integer, intent(in)         :: n_in      ! nodes sent per proc
 integer, intent(in)         :: n_synced  ! nodes in the synchronised array
 type(kdnode), intent(in)    :: node_in(n_in)
 type(kdnode), intent(inout) :: node_synced(n_synced)

#ifdef MPI
 integer                     :: dtype_kdnode
 integer                     :: mpierr

!  If there is only 1 owner, do a direct copy
 if (n_in == n_synced) then
    node_synced(:) = node_in(:)
 else
    call get_mpitype_of_kdnode(dtype_kdnode)
    ! skip if we are not an owner
    if (id == ifirstingroup) then
       ! perform node exchange
       call MPI_ALLGATHER(node_in,n_in,dtype_kdnode, &
                          node_synced,n_in,dtype_kdnode, &
                          comm_owner(level+1),mpierr)
    endif
 endif
#endif

end subroutine tree_sync

!----------------------------------------------------------------
!+
!  broadcast the tree from owners to non-owners
!+
!----------------------------------------------------------------
subroutine tree_bcast(node, nnode, level)
 use dtypekdtree, only:get_mpitype_of_kdnode

 integer,      intent(in)        :: nnode
 type(kdnode), intent(inout)     :: node(nnode)
 integer,      intent(in)        :: level

#ifdef MPI
 integer                         :: dtype_kdnode
 integer                         :: mpierr

 call get_mpitype_of_kdnode(dtype_kdnode)
 ! the bcast root is relative to the communicator (i.e. it needs to be 0, not ifirstingroup)
 call MPI_BCAST(node, nnode, dtype_kdnode, 0, comm_cofm(level+1), mpierr)
#endif

end subroutine tree_bcast

!----------------------------------------------------------------
!+
!  check which threads have completed
!+
!----------------------------------------------------------------
#ifdef MPI
subroutine check_complete(counters)
 use io, only:fatal
 integer, intent(inout) :: counters(nprocs,3)
 integer :: i
 logical :: countreceived
 integer :: mpierr
 integer :: status(MPI_STATUS_SIZE)

 ncomplete = 1 !self
 do i=1,nprocs
    if (i /= id + 1) then
       call MPI_TEST(countrequest(i),countreceived,status,mpierr)
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
end subroutine check_complete
#endif

#ifdef MPI
subroutine check_complete_force(counters,ncomplete_mpi)
 use io, only:fatal
 integer, intent(inout) :: counters(nprocs,3)
 integer, intent(out)   :: ncomplete_mpi
 integer :: i
 logical :: countreceived
 integer :: mpierr
 integer :: status(MPI_STATUS_SIZE)

 ncomplete_mpi = 1 !self
 do i=1,nprocs
    if (i /= id + 1) then
       call MPI_TEST(countrequest(i),countreceived,status,mpierr)
       if (countreceived) then
          if (counters(i,irecv) == counters(i,iexpect)) then
             ncomplete_mpi = ncomplete_mpi + 1
          elseif (counters(i,irecv) > counters(i,iexpect)) then
             print*,'on',id,'from',i-1
             print*,'nrecv',counters(i,irecv)
             print*,'nexpect',counters(i,iexpect)
             call fatal('mpiderivs', 'received more cells than expected')
          endif
       endif
    endif
 enddo
end subroutine check_complete_force
#endif

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

 counters(:,isent)   = 0
 counters(:,iexpect) = -1
 counters(:,irecv)   = 0

 !$omp master
 do iproc=1,nprocs
    if (iproc /= id + 1) then
       call MPI_IRECV(counters(iproc,iexpect),1,MPI_INTEGER4,iproc-1, &
       MPI_ANY_TAG,comm_cellcount,countrequest(iproc),mpierr)
       if (mpierr /= 0) call fatal('reset_cell_counters','error in MPI_IRECV')
    endif
 enddo
 !$omp end master
#endif
end subroutine reset_cell_counters

end module mpiderivs
