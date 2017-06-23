!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: balance
!
!  DESCRIPTION:
!  This module moves the particles onto their correct processor
!
!  REFERENCES: None
!
!  OWNER: Conrad Chan
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, domain, io, mpi, mpiutils, part, timing
!+
!--------------------------------------------------------------------------
#ifdef MPI
module balance
 use mpi
 use mpiutils, only:mpierr,status,MPI_DEFAULT_REAL,reduceall_mpi
 use part,     only:ipartbufsize
 implicit none
 integer :: nsent,nrecv,npartnew,ncomplete
 integer, dimension(1), private :: irequestrecv,irequestsend
 real, dimension(ipartbufsize)  :: xsendbuf,xbuffer

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

 !
 !--count number of proceses that have completed sending particles
 !
 ncomplete = 0

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

 if (id==master .and. iverbose >= 3) call getused(tstart)
 if (id==master .and. iverbose >= 5) print*,'starting balance',npart

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

 if (iverbose >= 5) then
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
 if (iverbose >= 4) print*,'>> shuffle: thread ',id,' got ',npart,' of ',ntot

 if (ntot /= ntot_start) call fatal('balance','number of particles before and after balance not equal')
 if (id==master .and. iverbose >= 3) call printused(tstart)

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
    jpart = status(MPI_TAG)
    if (jpart == 0) then ! signal the end
       ncomplete = ncomplete + 1
    else
       if (jpart > maxp .or. jpart <= 0) call fatal('balance','error in receive tag',jpart)
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

       if (inew > 0 .and. inew <= maxp) then
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
!$omp critical
          ideadhead = ll(inew)
!$omp end critical
       else
          if (inew /= 0) call fatal('balance','error in dead particle list',inew)
          !
          !--make a new particle
          !
!$omp critical
          npartnew = npartnew + 1
!$omp end critical
          if (npartnew > maxp) call fatal('recv_part','npartnew > maxp',npartnew)
          call unfill_buffer(npartnew,xbuffer)
          ibelong(npartnew) = id
       endif
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

 !--copy the particle to the new processor
 if (newproc < 0 .or. newproc > nprocs-1) then
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
 integer :: newproc
 logical, parameter :: iamcomplete = .true.
 logical :: doreplace

!
!--send the complete signal to all other threads;
!  start receiving the complete signal from other threads
!
 do newproc=0,nprocs-1
    !
    !-- tag=0, rest is junk
    !
    call MPI_ISEND(xsendbuf,0,MPI_DEFAULT_REAL,newproc,0,MPI_COMM_WORLD,irequestsend(1),mpierr)
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
 enddo
 npart = npartnew

 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!
!--each processor do a dummy send to next processor to clear the last remaining receive
!  (we know the receive has been posted for this, so use RSEND)
!
 newproc = mod(id+1,nprocs)
 call MPI_RSEND(xsendbuf,0,MPI_DEFAULT_REAL,newproc,0,MPI_COMM_WORLD,mpierr)

 if (iverbose >= 4 .or. (iverbose >= 3 .and. (nsent > 0 .or. nrecv > 0))) then
    print*,'>> balance: thread ',id,' sent:',nsent,' received:',nrecv,' npart =',npartnew
 endif
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

!--double check that all receives are complete and free request handle
 call MPI_WAIT(irequestrecv(1),status,mpierr)
 call MPI_REQUEST_FREE(irequestrecv(1),mpierr)

end subroutine balance_finish

end module balance
#endif
