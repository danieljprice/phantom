!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: mpiderivs
!
!  DESCRIPTION:
!   This module handles the MPI exchange of information during the
!   density and force routines
!
!  REFERENCES: None
!
!  OWNER: Conrad Chan
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, dtypekdtree, io, mpi, mpidens, mpiforce, mpiutils,
!    stack
!+
!--------------------------------------------------------------------------
#ifdef MPI
module mpiderivs
 use mpi
 use io,             only:id,nprocs
 use dim,            only:maxprocs
 use mpiutils,       only:mpierr,status,MPI_DEFAULT_REAL, &
                          comm_cellexchange,comm_cellcount
 use dtypekdtree,    only:kdnode,ndimtree

 implicit none

 interface init_cell_exchange
  module procedure init_celldens_exchange,init_cellforce_exchange
 end interface

 interface send_cell
  module procedure send_celldens,send_cellforce
 end interface

 interface check_send_finished
  module procedure check_send_finished_dens,check_send_finished_force
 end interface

 interface recv_cells
  module procedure recv_celldens,recv_cellforce
 end interface

 interface finish_cell_exchange
  module procedure finish_celldens_exchange,finish_cellforce_exchange
 end interface

 interface recv_while_wait
  module procedure recv_while_wait_force,recv_while_wait_dens
 end interface

 private

 integer :: ncomplete

 integer :: dtype_celldens
 integer :: dtype_cellforce

 integer :: globallevel
 integer :: comm_cofm(maxprocs)  ! only comms up to globallevel are used
 integer :: comm_owner(maxprocs) ! only comms up to globallevel are used

 integer :: nsent(maxprocs)     ! counter for number of cells sent to i
 integer :: nexpect(maxprocs)   ! counter for number of cells expecting from i
 integer :: nrecv(maxprocs)     ! counter for number of cells received from i

 integer :: countrequest(maxprocs)

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
 public :: init_tree_comms
 public :: finish_tree_comms
 public :: reset_cell_counters

 interface reduce_group
  module procedure reduce_group_real, reduce_group_int
 end interface

contains

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
 integer                           :: iproc

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

end subroutine init_celldens_exchange

subroutine init_cellforce_exchange(xbufrecv,ireq)
 use io,       only:fatal
 use mpiforce, only:get_mpitype_of_cellforce,cellforce

 type(cellforce),    intent(inout) :: xbufrecv(nprocs)
 integer,            intent(out)   :: ireq(nprocs) !,nrecv
 integer                           :: iproc

!
!--use persistent communication type for receives
!  cannot do same for sends as there are different destinations,
!  unless we make a request for each processor
!
!  We post a receive for EACH processor, to match the number of sends
!

 call get_mpitype_of_cellforce(dtype_cellforce)

 do iproc=1,nprocs
    call MPI_RECV_INIT(xbufrecv(iproc),1,dtype_cellforce,iproc-1, &
                       MPI_ANY_TAG,comm_cellexchange,ireq(iproc),mpierr)
    if (mpierr /= 0) call fatal('init_cell_exchange','error in MPI_RECV_INIT')
!
!--start the persistent communication channel
!
    call MPI_START(ireq(iproc),mpierr)
    if (mpierr /= 0) call fatal('init_cell_exchange','error in MPI_START')
 enddo

 ncomplete = 0

end subroutine init_cellforce_exchange

!-----------------------------------------------------------------------
!+
!  Subroutine to broadcast particle buffer to a bunch of processors
!+
!-----------------------------------------------------------------------
subroutine send_celldens(cell,direction,irequestsend,xsendbuf)
 use mpidens,  only:celldens

 type(celldens),     intent(in)     :: cell
 integer,            intent(in)     :: direction
 integer,            intent(inout)  :: irequestsend(nprocs)
 type(celldens),     intent(out)    :: xsendbuf

 logical                            :: targets(nprocs)
 integer                            :: newproc

 xsendbuf = cell

 ! export
 if (direction == 0) then
    targets = cell%remote_export(1:nprocs)
    ! return
 elseif (direction == 1) then
    targets = .false.
    targets(cell%owner+1) = .true.
 endif

 irequestsend = MPI_REQUEST_NULL

 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       call MPI_ISEND(xsendbuf,1,dtype_celldens,newproc,1+direction,comm_cellexchange,irequestsend(newproc+1),mpierr)
       nsent(newproc+1) = nsent(newproc+1) + 1
    endif
 enddo

end subroutine send_celldens

subroutine send_cellforce(cell,direction,irequestsend,xsendbuf)
 use mpiforce, only:cellforce

 type(cellforce),    intent(in)     :: cell
 integer,            intent(in)     :: direction
 integer,            intent(inout)  :: irequestsend(nprocs)
 type(cellforce),    intent(out)    :: xsendbuf

 logical                            :: targets(nprocs)
 integer                            :: newproc

 xsendbuf = cell

 ! export
 if (direction == 0) then
    targets = cell%remote_export(1:nprocs)
    ! return
 elseif (direction == 1) then
    targets = .false.
    targets(cell%owner+1) = .true.
 endif

 irequestsend = MPI_REQUEST_NULL

 do newproc=0,nprocs-1
    if ((newproc /= id) .and. (targets(newproc+1))) then ! do not send to self
       call MPI_ISEND(xsendbuf,1,dtype_cellforce,newproc,1+direction,comm_cellexchange,irequestsend(newproc+1),mpierr)
       nsent(newproc+1) = nsent(newproc+1) + 1
    endif
 enddo

end subroutine send_cellforce

!-----------------------------------------------------------------------
!+
!  Subroutine to check that non-blocking send has completed
!+
!-----------------------------------------------------------------------
subroutine check_send_finished_dens(stack,irequestsend,irequestrecv,xrecvbuf)
 use mpidens,  only:stackdens,celldens
 type(stackdens),    intent(inout)  :: stack
 integer,            intent(inout)  :: irequestsend(nprocs),irequestrecv(nprocs)
 type(celldens),     intent(inout)  :: xrecvbuf(nprocs)

 logical :: idone(nprocs)
 integer :: newproc
 !
 !--wait for broadcast to complete, continue to receive whilst doing so
 !
 idone(:) = .false.
 idone(id+1) = .true.
 do while(.not.all(idone))
    do newproc=0,nprocs-1
       if (newproc /= id) call MPI_TEST(irequestsend(newproc+1),idone(newproc+1),status,mpierr)
    enddo
    !--post receives
    call recv_celldens(stack,xrecvbuf,irequestrecv)
 enddo

end subroutine check_send_finished_dens

subroutine check_send_finished_force(stack,irequestsend,irequestrecv,xrecvbuf)
 use mpiforce, only:stackforce,cellforce
 type(stackforce),   intent(inout)  :: stack
 integer,            intent(inout)  :: irequestsend(nprocs),irequestrecv(nprocs)
 type(cellforce),    intent(inout)  :: xrecvbuf(nprocs)

 logical :: idone(nprocs)
 integer :: newproc
 !
 !--wait for broadcast to complete, continue to receive whilst doing so
 !
 idone(:) = .false.
 idone(id+1) = .true.
 do while(.not.all(idone))
    do newproc=0,nprocs-1
       if (newproc /= id) call MPI_TEST(irequestsend(newproc+1),idone(newproc+1),status,mpierr)
    enddo
    !--post receives
    call recv_cellforce(stack,xrecvbuf,irequestrecv)
 enddo

end subroutine check_send_finished_force

subroutine recv_while_wait_dens(stack,xrecvbuf,irequestrecv,irequestsend)
 use mpidens,  only:stackdens,celldens
 type(stackdens),  intent(inout) :: stack
 type(celldens),   intent(inout) :: xrecvbuf(nprocs)
 integer,          intent(inout) :: irequestrecv(nprocs),irequestsend(nprocs)

 integer             :: newproc

 do newproc=0,nprocs-1
    if (newproc /= id) then
       !--tag=0 to signal done
       call MPI_ISEND(nsent(newproc+1),1,MPI_INTEGER4,newproc,0,comm_cellcount,irequestsend(newproc+1),mpierr)
    endif
 enddo

 !--do not need to MPI_WAIT, because the following code requires the sends to go through
 do while (ncomplete < nprocs)
    call recv_celldens(stack,xrecvbuf,irequestrecv)
    call check_complete
 enddo

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 !--reset counter for next round
 ncomplete = 0

end subroutine recv_while_wait_dens

subroutine recv_while_wait_force(stack,xrecvbuf,irequestrecv,irequestsend)
 use mpiforce, only:stackforce,cellforce
 type(stackforce), intent(inout) :: stack
 type(cellforce),  intent(inout) :: xrecvbuf(nprocs)
 integer,          intent(inout) :: irequestrecv(nprocs),irequestsend(nprocs)
 integer             :: newproc

 do newproc=0,nprocs-1
    if (newproc /= id) then
       !--tag=0 to signal done
       call MPI_ISEND(nsent(newproc+1),1,MPI_INTEGER4,newproc,0,comm_cellcount,irequestsend(newproc+1),mpierr)
    endif
 enddo

 !--do not need to MPI_WAIT, because the following code requires the sends to go through
 do while (ncomplete < nprocs)
    call recv_cellforce(stack,xrecvbuf,irequestrecv)
    call check_complete
 enddo

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 !--reset counter for next round
 ncomplete = 0

end subroutine recv_while_wait_force

!------------------------------------------------
!+
!  Receive cells
!+
!------------------------------------------------
subroutine recv_celldens(target_stack,xbuf,irequestrecv)
 use io,       only:fatal
 use stack,    only:push_onto_stack
 use mpidens,  only:stackdens,celldens

 type(celldens),     intent(inout)  :: xbuf(:)  ! just need memory address
 type(stackdens),    intent(inout)  :: target_stack
 integer,            intent(inout)  :: irequestrecv(nprocs)
 integer                            :: iproc,k,iwait
 logical                            :: igot

 ! receive MPI broadcast
 do iproc=1,nprocs
    call MPI_TEST(irequestrecv(iproc),igot,status,mpierr)
    if (mpierr /= 0) call fatal('recv_cell','error in MPI_TEST call')
    ! unpack results
    if (igot) then
       if (status(MPI_TAG) == 2) then
          iwait = xbuf(iproc)%waiting_index
          do k = 1,xbuf(iproc)%npcell
             target_stack%cells(iwait)%rhosums(:,k) = target_stack%cells(iwait)%rhosums(:,k) + xbuf(iproc)%rhosums(:,k)
             target_stack%cells(iwait)%nneigh(k) = target_stack%cells(iwait)%nneigh(k) + xbuf(iproc)%nneigh(k)
          enddo
          do k = 1,nprocs
             target_stack%cells(iwait)%remote_export(k) = target_stack%cells(iwait)%remote_export(k) &
                  .and. xbuf(iproc)%remote_export(k)
          enddo
          target_stack%cells(iwait)%nneightry = target_stack%cells(iwait)%nneightry + xbuf(iproc)%nneightry
          nrecv(iproc) = nrecv(iproc) + 1
       elseif (status(MPI_TAG) == 1) then
          call push_onto_stack(target_stack, xbuf(iproc))
          nrecv(iproc) = nrecv(iproc) + 1
       endif
       call MPI_START(irequestrecv(iproc),mpierr)
    endif
 enddo
end subroutine recv_celldens

subroutine recv_cellforce(target_stack,xbuf,irequestrecv)
 use io,       only:fatal
 use stack,    only:push_onto_stack
 use mpiforce, only:stackforce,cellforce

 type(cellforce),    intent(inout)  :: xbuf(:)  ! just need memory address
 type(stackforce),   intent(inout)  :: target_stack
 integer,            intent(inout)  :: irequestrecv(nprocs)
 integer                            :: iproc,k,iwait
 logical                            :: igot

 ! receive MPI broadcast
 do iproc=1,nprocs
    call MPI_TEST(irequestrecv(iproc),igot,status,mpierr)
    if (mpierr /= 0) call fatal('recv_cell','error in MPI_TEST call')

    ! unpack results
    if (igot) then
       if (status(MPI_TAG) == 2) then
          iwait = xbuf(iproc)%waiting_index
          do k = 1,xbuf(iproc)%npcell
             target_stack%cells(iwait)%fsums(:,k) = target_stack%cells(iwait)%fsums(:,k) + xbuf(iproc)%fsums(:,k)
             target_stack%cells(iwait)%dtdrag(k) = min(target_stack%cells(iwait)%dtdrag(k), xbuf(iproc)%dtdrag(k))
             target_stack%cells(iwait)%vsigmax(k) = max(target_stack%cells(iwait)%vsigmax(k), xbuf(iproc)%vsigmax(k))
          enddo
#ifdef GRAVITY
          do k = 1,20
             target_stack%cells(iwait)%fgrav(k) = target_stack%cells(iwait)%fgrav(k) + xbuf(iproc)%fgrav(k)
          enddo
#endif
          do k =1,nprocs
             target_stack%cells(iwait)%remote_export(k) = target_stack%cells(iwait)%remote_export(k) &
               .and. xbuf(iproc)%remote_export(k)
          enddo
          target_stack%cells(iwait)%ndrag = target_stack%cells(iwait)%ndrag + xbuf(iproc)%ndrag
          target_stack%cells(iwait)%nstokes = target_stack%cells(iwait)%nstokes + xbuf(iproc)%nstokes
          target_stack%cells(iwait)%nsuper = target_stack%cells(iwait)%nsuper + xbuf(iproc)%nsuper
          nrecv(iproc) = nrecv(iproc) + 1
       elseif (status(MPI_TAG) == 1) then
          call push_onto_stack(target_stack, xbuf(iproc))
          nrecv(iproc) = nrecv(iproc) + 1
       endif
       call MPI_START(irequestrecv(iproc),mpierr)
    endif
 enddo
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
 integer                            :: newproc,iproc
 type(celldens), intent(in)         :: xsendbuf

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

end subroutine finish_celldens_exchange

subroutine finish_cellforce_exchange(irequestrecv,xsendbuf)
 use io,       only:fatal
 use mpiforce, only:cellforce
 integer,            intent(inout)  :: irequestrecv(nprocs)
 integer                            :: newproc,iproc
 type(cellforce), intent(in)        :: xsendbuf

!
!--each processor do a dummy send to next processor to clear the last remaining receive
!  (we know the receive has been posted for this, so use RSEND)
!
 do newproc=0,nprocs-1
    call MPI_RSEND(xsendbuf,1,dtype_cellforce,newproc,1,comm_cellexchange,mpierr)
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

end subroutine finish_cellforce_exchange

!----------------------------------------------------------------
!+
!  initialise communicators for tree construction
!+
!----------------------------------------------------------------
subroutine init_tree_comms()
 integer :: level,groupsize,color

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

end subroutine init_tree_comms

subroutine finish_tree_comms()
 integer :: level

 do level = 0, globallevel
    call MPI_COMM_FREE(comm_cofm(level+1), mpierr)
 enddo
 do level = 1, globallevel
    call MPI_COMM_FREE(comm_owner(level+1), mpierr)
 enddo
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

 real                         :: cofmpart(3)

 cofmpart = xyzcofm * totmass_node
 call MPI_ALLREDUCE(totmass_node,totmassg,1,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 call MPI_ALLREDUCE(cofmpart,cofmsum,3,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 cofmsum = cofmsum / totmassg

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
 real                                  :: isend, ired
 real                                  :: xg

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

end function reduce_group_real

function reduce_group_int(x,string,level) result(xg)
 use io, only:fatal
 integer,            intent(in)        :: x
 character(len=*),   intent(in)        :: string
 integer,            intent(in)        :: level
 integer                               :: isend, ired
 integer                               :: xg

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

end function reduce_group_int

!----------------------------------------------------------------
!+
!  synchronize the global tree, placing nodes in the correct position
!+
!----------------------------------------------------------------
subroutine tree_sync(nodeentry, nnodes, node, ifirstingroup, groupsize, level)
 use dtypekdtree, only:get_mpitype_of_kdnode

 integer, intent(in)         :: ifirstingroup, groupsize, level
 integer, intent(in)         :: nnodes ! nodes sent per proc
 type(kdnode), intent(in)    :: nodeentry(nnodes)
 type(kdnode), intent(inout) :: node(nprocs/groupsize)

 integer                     :: dtype_kdnode

 integer                     :: nowners

 nowners = nprocs / groupsize

 ! only exchange if there is more than 1 owner (every level except top)
 if (nowners > 1) then
    call get_mpitype_of_kdnode(dtype_kdnode)
    ! skip if we are not an owner
    if (id == ifirstingroup) then
       ! perform node exchange
       call MPI_ALLGATHER(nodeentry,nnodes,dtype_kdnode,node,nnodes,dtype_kdnode,comm_owner(level+1),mpierr)
    endif
 else
    node = nodeentry(1)
 endif

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

 integer                         :: dtype_kdnode

 call get_mpitype_of_kdnode(dtype_kdnode)
 ! the bcast root is relative to the communicator (i.e. it needs to be 0, not ifirstingroup)
 call MPI_BCAST(node, nnode, dtype_kdnode, 0, comm_cofm(level+1), mpierr)

end subroutine tree_bcast

!----------------------------------------------------------------
!+
!  check which threads have completed
!+
!----------------------------------------------------------------
subroutine check_complete
 use io, only:fatal
 integer :: i
 logical :: countreceived

 ncomplete = 1 !self
 do i=1,nprocs
    if (i /= id + 1) then
       call MPI_TEST(countrequest(i),countreceived,status,mpierr)
       if (countreceived) then
          if (nrecv(i) == nexpect(i)) then
             ncomplete = ncomplete + 1
          elseif (nrecv(i) > nexpect(i)) then
             print*,'on',id,'from',i-1
             print*,'nrecv',nrecv(i)
             print*,'nexpect',nexpect(i)
             call fatal('mpiderivs', 'received more cells than expected')
          endif
       endif
    endif
 enddo
end subroutine check_complete

!----------------------------------------------------------------
!+
!  reset counters for checking arrival of all cells
!+
!----------------------------------------------------------------
subroutine reset_cell_counters
 use io, only:fatal
 integer :: iproc
 nsent(:) = 0
 nexpect(:) = -1
 nrecv(:) = 0

 do iproc=1,nprocs
    if (iproc /= id + 1) then
       call MPI_IRECV(nexpect(iproc),1,MPI_INTEGER4,iproc-1, &
       MPI_ANY_TAG,comm_cellcount,countrequest(iproc),mpierr)
       if (mpierr /= 0) call fatal('reset_cell_counters','error in MPI_IRECV')
    endif
 enddo
end subroutine reset_cell_counters

end module mpiderivs
#endif
