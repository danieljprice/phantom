!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpibalance
!
! This module moves the particles onto their correct processor
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, io, mpi, mpiutils, part, timing
!
#ifdef MPI
 use mpi
 use io,       only:id
 use mpiutils, only:mpierr,status,MPI_DEFAULT_REAL,reduceall_mpi, &
                    comm_balance,comm_balancecount,count_collective
 use part,     only:ipartbufsize
#endif
 use io,       only:nprocs

 implicit none

 public :: allocate_balance_arrays
 public :: deallocate_balance_arrays
 public :: balancedomains

 ! particles sent to other tasks, and calls to balancedomains (in total and skipped
 ! because no particles needed to move), for the MPI statistics
 integer(kind=8), public :: nmigrated = 0, nbalance_calls = 0, nbalance_skipped = 0

 private

#ifdef MPI
 ! particles are sent in batches of up to nbatch per message, with two buffers
 ! per destination so that one can be filled while the other is being sent
 integer, parameter   :: nbatch = 32
 real                 :: xsendbuf(ipartbufsize),xbuffer(ipartbufsize)
 real                 :: xrecvbatch(ipartbufsize*nbatch)  ! receives one message
 real,    allocatable :: xsendbatch(:,:,:)                ! (values, buffer, destination)
 integer, allocatable :: nbatched(:),icurbuf(:)           ! particles in, and index of, the current buffer
 integer, allocatable :: isendreq(:,:)                    ! (buffer, destination) send requests
 integer, allocatable :: countsendreq(:)                  ! sends of the particle counts
 integer              :: npartbuf                         ! reals per particle in the buffers
 integer :: irequestrecv(1)
 integer(kind=8)                           :: ntot_start
 integer                                   :: npartnew, ncomplete
#endif
 integer, allocatable :: nsent(:),nexpect(:),nrecv(:)
 integer, allocatable :: countrequest(:)

contains

!----------------------------------------------------------------
!+
!  allocate memory
!+
!----------------------------------------------------------------
subroutine allocate_balance_arrays()
 use allocutils, only:allocate_array

 call allocate_array('nsent',       nsent,       nprocs)
 call allocate_array('nexpect',     nexpect,     nprocs)
 call allocate_array('nrecv',       nrecv,       nprocs)
 call allocate_array('countrequest',countrequest,nprocs)

end subroutine allocate_balance_arrays

!----------------------------------------------------------------
!+
!  deallocate memory
!+
!----------------------------------------------------------------
subroutine deallocate_balance_arrays()

 if (allocated(nsent       )) deallocate(nsent       )
 if (allocated(nexpect     )) deallocate(nexpect     )
 if (allocated(nrecv       )) deallocate(nrecv       )
 if (allocated(countrequest)) deallocate(countrequest)
#ifdef MPI
 if (allocated(xsendbatch))   deallocate(xsendbatch,nbatched,icurbuf,isendreq,countsendreq)
#endif

end subroutine deallocate_balance_arrays

!----------------------------------------------------------------
!+
!  routine which moves particles onto their correct processor
!  (for any domain decomposition)
!+
!----------------------------------------------------------------
subroutine balancedomains(npart)
#ifndef MPI
 integer, intent(inout) :: npart
#else
 use io,     only:id,master,iverbose,fatal
 use part,   only:shuffle_part,count_dead_particles,ibelong,update_npartoftypetot,&
                   recount_npartoftype
 use timing, only:getused,printused
 use mpiutils, only:barrier_mpi
 integer, intent(inout) :: npart
 integer :: i,newproc
 integer(kind=8) :: ntot
 real(kind=4) :: tstart
 ! initialise to check for receives every 2 particles
 ! note that this is saved across calls, so it remembers the previous value
 integer :: check_interval=2
 integer :: recv_gap
 logical :: gotpart

 gotpart = .false.

 ! Balance is only necessary when there are more than 1 MPI tasks
 ! and some task has particles to send or dead particles to remove
 if (nprocs > 1) then
    nbalance_calls = nbalance_calls + 1
    if (reduceall_mpi('max',nparts_to_move(npart)) == 0) then
       nbalance_skipped = nbalance_skipped + 1
       return
    endif

    if (id==master .and. iverbose >= 3) call getused(tstart)
    if (id==master .and. iverbose >= 5) print*,'starting balance',npart

    call balance_init(npart)
    ntot_start = reduceall_mpi('+',npart - count_dead_particles())

    recv_gap = 0
    do i=1,npart
       !
       !--attempt to receive a particle
       !
       if (mod(i,check_interval) == 0) then
          call recv_part(gotpart)
          !
          !--if there is a long gap between receiving particles, then decrease
          !  the check frequency. if there is no gap, then increase frequency.
          !  the frequency should roughly match the rate that particles are
          !  available to be received.
          !
          if (gotpart) then
             ! Halve period if 2 consecutive checks both return a particle
             if (recv_gap == 0) check_interval = max(1, check_interval/2)
             recv_gap = 0
          else
             recv_gap = recv_gap + 1
             ! Double period if 4 consecutive checks do not return a particle
             if (recv_gap > 4) then
                check_interval = check_interval * 2
                ! Half the gap counter to correspond to current frequency
                recv_gap = recv_gap / 2
             endif
          endif
       endif
       !
       !--send particles which belong to other processors
       !
       newproc = ibelong(i)
       if (newproc /= id) call send_part(i,newproc)
    enddo
    call flush_all_batches()

    if (iverbose >= 5) then
       print*,id,' finished send, nsent = ',sum(nsent(1:nprocs)),' npart = ',npartnew
       print*,id,' received so far ',sum(nrecv(1:nprocs))
    endif
    call balance_finish(npart)
    nmigrated = nmigrated + sum(nsent(1:nprocs))
    ! ndead = count_dead_particles()
    ! print*,' thread ',id,' before shuffle, got ',npart,' dead ',ndead,' actual = ',npart - ndead,ideadhead
    call shuffle_part(npart)
    !--count types only once the particles sent away have been removed
    call recount_npartoftype
    call barrier_mpi()
    ! ndead = count_dead_particles()
    ! print*,' thread ',id,' after shuffle, got ',npart,' dead ',ndead,' actual = ',npart - ndead

    ntot = reduceall_mpi('+',npart)
    if (iverbose >= 4) print*,'>> shuffle: thread ',id,' got ',npart,' of ',ntot

    if (ntot /= ntot_start) then
       print*,id,'ntot',ntot
       print*,id,'ntot_start',ntot_start
       call fatal('balance','number of particles before and after balance not equal')
    endif

    !  Update particle types
    call update_npartoftypetot

    if (id==master .and. iverbose >= 3) call printused(tstart)

 endif
#endif
end subroutine balancedomains

#ifdef MPI
!----------------------------------------------------------------
!+
!  number of particles on this task that belong to another task,
!  plus one if there are dead particles to be removed
!+
!----------------------------------------------------------------
integer function nparts_to_move(npart)
 use io,   only:id
 use part, only:ibelong,ideadhead
 integer, intent(in) :: npart
 integer :: i,n

 n = 0
 !$omp parallel do default(none) shared(npart,ibelong,id) private(i) reduction(+:n)
 do i=1,npart
    if (ibelong(i) /= id) n = n + 1
 enddo
 !$omp end parallel do
 if (ideadhead /= 0) n = n + 1
 nparts_to_move = n

end function nparts_to_move
#endif

#ifdef MPI
!----------------------------------------------------------------
!+
!  initialisation of load balancing process,
!  declaration of types etc.
!+
!----------------------------------------------------------------
subroutine balance_init(npart)
 use part, only:fill_sendbuf
 integer, intent(in) :: npart
 integer :: i

 ! use a dummy call to fill_sendbuf to find out the number of reals per particle
 call fill_sendbuf(1,xbuffer,npartbuf) ! just fill for particle #1

 if (.not.allocated(xsendbatch)) then
    allocate(xsendbatch(npartbuf*nbatch,2,nprocs),nbatched(nprocs),icurbuf(nprocs))
    allocate(isendreq(2,nprocs),countsendreq(nprocs))
 endif
 nbatched = 0
 icurbuf  = 1
 isendreq = MPI_REQUEST_NULL
 countsendreq = MPI_REQUEST_NULL

!--use persistent communication type for receives, of up to nbatch particles
!
 call MPI_RECV_INIT(xrecvbatch,size(xrecvbatch),MPI_DEFAULT_REAL,MPI_ANY_SOURCE, &
                    MPI_ANY_TAG,comm_balance,irequestrecv(1),mpierr)
!
!--post a non-blocking receive so that we can receive particles
!
 call MPI_START(irequestrecv(1),mpierr)

 !
 !--count checking
 !
 nsent(:) = 0
 nexpect(:) = -1
 nrecv(:) = 0
 do i=1,nprocs
    if (id /= i-1) then
       call MPI_IRECV(nexpect(i),1,MPI_INTEGER4,i-1, &
                      MPI_ANY_TAG,comm_balancecount,countrequest(i),mpierr)
    endif
 enddo

 npartnew = npart

 !
 !--count number of proceses that have completed sending particles
 !
 ncomplete = 0

end subroutine balance_init

!-----------------------------------------------------------------------
!+
!  checks for a message of particles and, if one has arrived, places
!  them. With replace (the default) they fill the slots of particles
!  already sent (the dead particle list headed by ideadhead), otherwise
!  they are added to the end of the array.
!+
!-----------------------------------------------------------------------
subroutine recv_part(replace,gotpart)
 logical, intent(in),  optional :: replace
 logical, intent(out), optional :: gotpart
 logical :: igotpart,doreplace
 integer :: ncount,np,k

 doreplace = .true.
 if (present(replace)) doreplace = replace

 igotpart = .false.
 call MPI_TEST(irequestrecv(1),igotpart,status,mpierr)

 if (present(gotpart)) gotpart = igotpart

 if (igotpart) then
    call MPI_GET_COUNT(status,MPI_DEFAULT_REAL,ncount,mpierr)
    np = ncount/npartbuf
    do k=1,np
       xbuffer(1:npartbuf) = xrecvbatch((k-1)*npartbuf+1:k*npartbuf)
       call place_particle(doreplace)
    enddo
    nrecv(status(MPI_SOURCE)+1) = nrecv(status(MPI_SOURCE)+1) + np
    !
    !--post another receive ready for the next message
    !
    call MPI_START(irequestrecv(1),mpierr)
 endif

end subroutine recv_part

!-----------------------------------------------------------------------
!+
!  place the particle in xbuffer in the slot of a particle already sent,
!  or at the end of the particle arrays
!+
!-----------------------------------------------------------------------
subroutine place_particle(replace)
 use io,   only:fatal,id
 use part, only:isdead,unfill_buffer,maxp,ll,ideadhead,ibelong
 logical, intent(in) :: replace
 integer :: inew

 inew = 0
 if (replace) inew = ideadhead

 if (inew > 0 .and. inew <= maxp) then
    if (.not.isdead(inew)) call fatal('balance','replacing non-dead particle')
    !--replace a particle which has already been sent
    call unfill_buffer(inew,xbuffer)
    !--assume that this particle landed in the right place
    ibelong(inew) = id
    ideadhead = ll(inew)
 else
    if (inew /= 0) call fatal('balance','error in dead particle list',inew)
    !--make a new particle
    npartnew = npartnew + 1
    if (npartnew > maxp) call fatal('recv_part','npartnew > maxp',npartnew)
    call unfill_buffer(npartnew,xbuffer)
    ibelong(npartnew) = id
 endif

end subroutine place_particle

!-----------------------------------------------------------------------
!+
!  add a particle to the batch for its new processor, sending the batch
!  when it is full, and kill it on this processor
!+
!-----------------------------------------------------------------------
subroutine send_part(i,newproc,replace)
 use io,   only:fatal
 use part, only:fill_sendbuf,kill_particle
 integer, intent(in) :: i,newproc
 logical, intent(in), optional :: replace
 logical :: doreplace
 integer :: nbuf,ip,ioff

 doreplace = .true.
 if (present(replace)) doreplace = replace

 if (newproc < 0 .or. newproc > nprocs-1) then
    call fatal('balance','error in ibelong',ival=newproc,var='ibelong')
 endif
 ip = newproc + 1
 call fill_sendbuf(i,xsendbuf,nbuf)
 ioff = nbatched(ip)*npartbuf
 xsendbatch(ioff+1:ioff+nbuf,icurbuf(ip),ip) = xsendbuf(1:nbuf)
 nbatched(ip) = nbatched(ip) + 1

 !--kill particle on this processor
 call kill_particle(i)
 nsent(ip) = nsent(ip) + 1

 if (nbatched(ip) == nbatch) call flush_batch(ip,doreplace)

end subroutine send_part

!-----------------------------------------------------------------------
!+
!  send the current batch for processor ip-1, then make the other buffer
!  current, waiting (and receiving meanwhile) until its last send is done
!+
!-----------------------------------------------------------------------
subroutine flush_batch(ip,replace)
 integer, intent(in) :: ip
 logical, intent(in) :: replace
 integer :: ib

 if (nbatched(ip) == 0) return
 ib = icurbuf(ip)
 ! tag cannot be i, because some MPI implementations do not support large values for the tag
 call MPI_ISEND(xsendbatch(1,ib,ip),nbatched(ip)*npartbuf,MPI_DEFAULT_REAL,ip-1,0,&
                comm_balance,isendreq(ib,ip),mpierr)
 nbatched(ip) = 0
 icurbuf(ip)  = 3 - ib
 call wait_for_send(isendreq(icurbuf(ip),ip),replace)

end subroutine flush_batch

!-----------------------------------------------------------------------
!+
!  send the remaining partial batches and wait until all sends are done
!+
!-----------------------------------------------------------------------
subroutine flush_all_batches()
 integer :: ip,ib

 do ip=1,nprocs
    call flush_batch(ip,.true.)
 enddo
 do ip=1,nprocs
    do ib=1,2
       call wait_for_send(isendreq(ib,ip),.true.)
    enddo
 enddo

end subroutine flush_all_batches

!-----------------------------------------------------------------------
!+
!  wait for a send to complete, receiving particles meanwhile
!+
!-----------------------------------------------------------------------
subroutine wait_for_send(ireq,replace)
 integer, intent(inout) :: ireq
 logical, intent(in)    :: replace
 logical :: idone

 idone = .false.
 do while(.not.idone)
    call MPI_TEST(ireq,idone,status,mpierr)
    if (.not.idone) call recv_part(replace=replace)
 enddo

end subroutine wait_for_send

!----------------------------------------------------------------
!+
!  finish/clean up of load balancing process.
!+
!----------------------------------------------------------------
subroutine balance_finish(npart,replace)
 use io,    only:id,fatal,iverbose
 integer, intent(out) :: npart
 logical, intent(in), optional :: replace
 integer             :: newproc
 logical, parameter  :: iamcomplete = .true.
 logical             :: doreplace

!
!--send the complete signal to all other threads;
!  start receiving the complete signal from other threads
!
 do newproc=0,nprocs-1
    if (newproc /= id) then
       call MPI_ISEND(nsent(newproc+1),1,MPI_INTEGER4,newproc,0,comm_balancecount,countsendreq(newproc+1),mpierr)
    endif
 enddo

 if (present(replace)) then
    doreplace = replace
 else
    doreplace = .true.
 endif
!
!--continue to check for receive signals until all sends are complete
!
 do while (ncomplete < nprocs)
    call recv_part(replace=doreplace)
    call check_complete
 enddo
 npart = npartnew
 call MPI_WAITALL(nprocs,countsendreq,MPI_STATUSES_IGNORE,mpierr)

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 call count_collective
!
!--each processor do a dummy send to next processor to clear the last remaining receive
!  (we know the receive has been posted for this, so use RSEND)
!
 newproc = mod(id+1,nprocs)
 call MPI_RSEND(xsendbuf,0,MPI_DEFAULT_REAL,newproc,0,comm_balance,mpierr)

 if (iverbose >= 4 .or. (iverbose >= 3 .and. (sum(nsent(1:nprocs)) > 0 .or. sum(nrecv(1:nprocs)) > 0))) then
    print*,'>> balance: thread ',id,' sent:',sum(nsent(1:nprocs)),' received:',sum(nrecv(1:nprocs)),' npart =',npartnew
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 call count_collective

!--double check that all receives are complete and free request handle
 call MPI_WAIT(irequestrecv(1),status,mpierr)
 call MPI_REQUEST_FREE(irequestrecv(1),mpierr)

end subroutine balance_finish

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
             call fatal('mpibalance', 'received more particles than expected')
          endif
       endif
    endif
 enddo
end subroutine check_complete
#endif

end module mpibalance
