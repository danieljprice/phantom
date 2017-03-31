! $Id$
!----------------------------------------------------------------
! This module is part of the Phantom SPH code
! Use is by specific, written permission of the author
! (c) 2007 Daniel Price
!----------------------------------------------------------------
!+
!  This module moves the particles onto their correct processor
!+
!----------------------------------------------------------------
#ifdef MPI
module balance
 use mpi
 use mpiutils, only:mpierr,status,MPI_DEFAULT_REAL,reduceall_mpi
 use part,     only:ipartbufsize
 implicit none
 integer :: comm_done,nsent,nrecv,npartnew
 integer, dimension(1), private :: irequestrecv, irequestcomplete,&
                                   irequestsend, irequestend
 real, dimension(ipartbufsize) :: xsendbuf,xbuffer

 private
 public :: balancedomains
 public :: balance_init,balance_finish,send_part,recv_part
 integer(kind=8) :: ntot_start

contains

!----------------------------------------------------------------
!+
!  initialisation of load balancing process,
!  declaration of types etc.
!+
!----------------------------------------------------------------
subroutine balance_init(npart)
 implicit none
 integer, intent(in) :: npart
!
!--create a new communicator for the `completed' signal.
!  Could just use a different tag on those signals (in fact we do anyway),
!  but the particle exchange is done by sending the old label of the particle
!  on the previous processor as the tag, thus requiring MPI_ANY_TAG
!  for the receiving end. This way it is possible in future
!  to be able to use these tags to track particle identities.
!
 call MPI_COMM_DUP(MPI_COMM_WORLD,comm_done,mpierr)

!--use persistent communication type for receives
!  cannot do same for sends as there are different destination,
!  unless we make a request for each processor
!
 call MPI_RECV_INIT(xbuffer,size(xbuffer),MPI_DEFAULT_REAL,MPI_ANY_SOURCE, &
                    MPI_ANY_TAG,MPI_COMM_WORLD,irequestrecv(1),mpierr)
!
!--post a non-blocking receive so that we can receive particles
!
 call MPI_START(irequestrecv(1),mpierr)

 nsent = 0
 nrecv = 0
 npartnew = npart
 ntot_start = reduceall_mpi('+',npart)

 return
end subroutine balance_init

!----------------------------------------------------------------
!+
!  routine which moves particles onto their correct processor
!  (for any domain decomposition)
!+
!----------------------------------------------------------------
subroutine balancedomains(npart)
 use io,     only:id,master,iverbose,fatal
 use domain, only:ibelong
 use part,   only:shuffle_part,count_dead_particles
 use timing, only:getused,printused
 use mpiutils, only:barrier_mpi
 implicit none
 integer, intent(inout) :: npart
 integer :: i,newproc,ndead
 integer(kind=8) :: ntot
 real(kind=4) :: tstart

 if (id.eq.master .and. iverbose.ge.3) call getused(tstart)
 if (id.eq.master .and. iverbose.ge.5) print*,'starting balance',npart

 call balance_init(npart)

 do i=1,npart
!
!--attempt to receive particles
!
    call recv_part()
!
!--send particles which belong to other processors
!
    newproc = ibelong(i)
    if (newproc /= id) call send_part(i,newproc)

 enddo

 if (iverbose.ge.5) then
    print*,id,' finished send, nsent = ',nsent,' npart = ',npartnew
    print*,id,' received so far ',nrecv
 endif
 call balance_finish(npart)
 ndead = count_dead_particles()
 !print*,' thread ',id,' before shuffle, got ',npart,' dead ',ndead,' actual = ',npart - ndead,ideadhead
 call shuffle_part(npart)
 call barrier_mpi()
 ndead = count_dead_particles()
 !print*,' thread ',id,' after shuffle, got ',npart,' dead ',ndead,' actual = ',npart - ndead

 ntot = reduceall_mpi('+',npart)
 if (iverbose >= 3) print*,'>> shuffle: thread ',id,' got ',npart,' of ',ntot

 if (ntot /= ntot_start) call fatal('balance','number of particles before and after balance not equal')
 if (id.eq.master .and. iverbose.ge.3) call printused(tstart)

 return
end subroutine balancedomains

!-----------------------------------------------------------------------
!+
!  function which checks for particles to receive and receives them
!  if necessary. Non-zero ideadhead refers to start of dead particle
!  list in which to place particles. Sending in zero for ideadhead means
!  simply add new particles to the end of the array.
!+
!-----------------------------------------------------------------------
subroutine recv_part(replace)
 use io,      only:fatal,id
 use part,    only:isdead,unfill_buffer,maxp,ll,ideadhead
 use domain,  only:ibelong
 implicit none
 logical, intent(in), optional :: replace
 logical :: igotpart
 integer :: jpart,inew

 igotpart = .false.
 call MPI_TEST(irequestrecv(1),igotpart,status,mpierr)

 if (igotpart) then
    !print*,id,' got particle ',status(MPI_TAG),' from ',status(MPI_SOURCE)
    jpart = status(MPI_TAG)
    if (jpart.gt.maxp .or. jpart.le.0) call fatal('balance','error in receive tag',jpart)
!$omp critical
    nrecv = nrecv + 1
!$omp end critical
    if (present(replace)) then
       if (replace) then
          inew = ideadhead
       else
          inew = 0
       endif
    else
       inew = ideadhead
    endif

    if (inew.gt.0 .and. inew.le.maxp) then
       if (.not.isdead(inew)) &
          call fatal('balance','replacing non-dead particle')
    !
    !--replace a particle which has already been sent
    !
       call unfill_buffer(inew,xbuffer)
       !
       !--assume that this particle landed in the right place
       !
       ibelong(inew) = id
       !if (ibelong(inew).ne.id) &
       !   call fatal('balance','received particle does not belong')
!$omp critical
       ideadhead = ll(inew)
!$omp end critical
    else
       if (inew.ne.0) call fatal('balance','error in dead particle list',inew)
    !
    !--make a new particle
    !
!$omp critical
       npartnew = npartnew + 1
!$omp end critical
       if (npartnew.gt.maxp) call fatal('recv_part','npartnew > maxp',npartnew)
       call unfill_buffer(npartnew,xbuffer)
       ibelong(npartnew) = id

    endif
    !
    !--post another receive ready for next particle
    !
    call MPI_START(irequestrecv(1),mpierr)
 endif

 return
end subroutine recv_part

!-----------------------------------------------------------------------
!+
!  function which sends a particle onto a remote processor
!+
!-----------------------------------------------------------------------
subroutine send_part(i,newproc,replace)
 use io,   only:fatal,nprocs
 use part, only:fill_sendbuf,kill_particle
 implicit none
 integer, intent(in) :: i,newproc
 logical, intent(in), optional :: replace
 logical :: idone,doreplace

 if (present(replace)) then
    doreplace = replace
 else
    doreplace = .true.
 endif
 ! print*,id,' sending ',i,' to ',newproc
 !--copy the particle to the new processor
 if (newproc.lt.0 .or. newproc.gt.nprocs-1) then
    call fatal('balance','error in ibelong',ival=newproc,var='ibelong')
 else
    call fill_sendbuf(i,xsendbuf)
    call MPI_ISEND(xsendbuf,size(xsendbuf),MPI_DEFAULT_REAL,newproc,i,MPI_COMM_WORLD,irequestsend(1),mpierr)

    !--wait for send to complete, receive whilst doing so
    idone = .false.
    do while(.not.idone)
       call MPI_TEST(irequestsend(1),idone,status,mpierr)
       call recv_part(replace=doreplace)
    enddo
 endif
 !--kill particle on this processor
 call kill_particle(i)

 nsent = nsent + 1
 return
end subroutine send_part

!----------------------------------------------------------------
!+
!  finish/clean up of load balancing process.
!+
!----------------------------------------------------------------
subroutine balance_finish(npart,replace)
 use dim, only:maxp
 use io,  only:id,nprocs,fatal,iverbose
 implicit none
 integer, intent(out) :: npart
 logical, intent(in), optional :: replace
 integer :: ncomplete,newproc
 logical, parameter :: iamcomplete = .true.
 logical :: jiscomplete,idone,doreplace
!
!--send the complete signal to all other threads;
!  start receiving the complete signal from other threads
!
 !print*,' thread ',id,' sending complete signal'
 do newproc=0,nprocs-1
    call MPI_ISEND(iamcomplete,1,MPI_LOGICAL,newproc,maxp+1,comm_done,irequestend(1),mpierr)
 enddo
 jiscomplete = .false.
 call MPI_IRECV(jiscomplete,1,MPI_LOGICAL,MPI_ANY_SOURCE,maxp+1,comm_done,irequestcomplete(1),mpierr)

 if (present(replace)) then
    doreplace = replace
 else
    doreplace = .true.
 endif
!
!--continue to check for receive signals until all sends are complete
!
 ncomplete = 0
 do while (.not.allcomplete())
    call recv_part(replace=doreplace)
 enddo
 npart = npartnew

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!
!--each processor do a dummy send to next processor to clear the last remaining receive
!  (we know the receive has been posted for this, so use RSEND)
!
 newproc = mod(id+1,nprocs)
 call MPI_RSEND(xsendbuf,0,MPI_DEFAULT_REAL,newproc,0,MPI_COMM_WORLD,mpierr)
! call send_particle(0,newproc)

 if (iverbose.ge.3 .or. (iverbose.ge.1 .and. (nsent.gt.0 .or. nrecv.gt.0))) then
    print*,'>> balance: thread ',id,' sent:',nsent,' received:',nrecv,' npart =',npartnew
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

!--double check that all receives are complete and free request handle
 call MPI_TEST(irequestrecv(1),idone,status,mpierr)
 call MPI_REQUEST_FREE(irequestrecv(1),mpierr)
 if (.not.idone) call fatal('balance','receive has not completed')

!--deallocate communicator
 call MPI_COMM_FREE(comm_done,mpierr)

 return

contains
 !-----------------------------------------------------------------------
 !+
 !  query function to determine whether all threads have finished send
 !+
 !-----------------------------------------------------------------------
 logical function allcomplete()
  implicit none
  logical :: isdone

  allcomplete = .false.
 !
 !--test if the receive we have already posted has been completed
 !
  CALL MPI_TEST(irequestcomplete(1),isdone,status,mpierr)
 !
 !--if so, update the counter for the number of threads which have completed
 !  and, if there are threads remaining, post a non-blocking receive
 !  to get the signal from the next processor to complete
 !
  if (isdone) then
     ncomplete = ncomplete + 1
     !print*,id,' ncomplete = ',ncomplete,' got complete signal from ',status(MPI_SOURCE)
     if (ncomplete.lt.nprocs) then
        CALL MPI_IRECV(jiscomplete,1,MPI_LOGICAL,MPI_ANY_SOURCE,maxp+1,comm_done,irequestcomplete(1),mpierr)
     else
        allcomplete = .true.
     endif
  endif
 end function allcomplete

end subroutine balance_finish

end module balance
#endif
