!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, mpi, mpiutils, part
!+
!--------------------------------------------------------------------------
#ifdef MPI
module mpiderivs
 use mpi
 use mpiutils, only:mpierr,status,MPI_DEFAULT_REAL
 implicit none
 private
 public :: init_results_exchange
 public :: send_results
 public :: recv_density_results
 public :: recv_force_results
 public :: check_send_finished
 public :: finish_results_exchange

contains

!----------------------------------------------------------------
!+
!  initialise the receive type for each thread
!+
!----------------------------------------------------------------
subroutine init_results_exchange(nprocs,xbufrecv,ireq)
 integer, intent(in)    :: nprocs
 real,    intent(inout) :: xbufrecv(:,:)
 integer, intent(out)   :: ireq(nprocs) !,nrecv
 integer :: iproc
!
!--use persistent communication type for receives
!  cannot do same for sends as there are different destinations,
!  unless we make a request for each processor
!
!  We post a receive for EACH processor, to match the number of sends
!
 do iproc=1,nprocs
    call MPI_RECV_INIT(xbufrecv(:,iproc),size(xbufrecv(:,iproc)),MPI_DEFAULT_REAL,MPI_ANY_SOURCE, &
                       MPI_ANY_TAG,MPI_COMM_WORLD,ireq(iproc),mpierr)
!
!--start the persistent communication channel
!
    call MPI_START(ireq(iproc),mpierr)
 enddo

end subroutine init_results_exchange

!-----------------------------------------------------------------------
!+
!  Subroutine to broadcast particle buffer to all remote processors
!+
!-----------------------------------------------------------------------
subroutine send_results(nprocs,xsendbuf,nsend,i,irequestsend)
 use io, only:id
 real,    intent(in) :: xsendbuf(:)
 integer, intent(in) :: nprocs,nsend,i
 integer, intent(inout) :: irequestsend(nprocs)
 integer :: newproc

 do newproc=0,nprocs-1
    ! send particle info to remote processor, using particle id as the tag
    if (newproc /= id) then ! do not send to self
       call MPI_ISEND(xsendbuf,nsend,MPI_DEFAULT_REAL,newproc,i,MPI_COMM_WORLD,irequestsend(newproc+1),mpierr)
    endif
 enddo

end subroutine send_results

!-----------------------------------------------------------------------
!+
!  Subroutine to check that non-blocking send has completed
!+
!-----------------------------------------------------------------------
subroutine check_send_finished(nprocs,irequestsend,irequestrecv,xrecvbuf,nrecv,getdv,getdB,realviscosity)
 use io, only:id
 integer, intent(in) :: nprocs
 integer, intent(inout) :: irequestsend(nprocs),irequestrecv(nprocs),nrecv
 real,    intent(inout) :: xrecvbuf(:,:)
 logical, intent(in), optional    :: getdv,getdB,realviscosity
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
    if (present(getdv) .and. present(getdB)) then
       call recv_density_results(nprocs,xrecvbuf,irequestrecv,nrecv,getdv,getdB,realviscosity)
    else
       call recv_force_results(nprocs,xrecvbuf,irequestrecv,nrecv)
    endif
 enddo

end subroutine check_send_finished

!-----------------------------------------------------------------------------
!+
!  Subroutine to receive broadcast of density results from remote processors
!+
!-----------------------------------------------------------------------------
subroutine recv_density_results(nprocs,xbuf,irequestrecv,nrecv,getdv,getdB,realviscosity)
 use part, only:xyzh,gradh,divcurlv,divcurlB,straintensor,maxgradh,npart,ngradh
 use dim,  only:ndivcurlv,ndivcurlB,maxp,maxstrain
 use io,   only:fatal
 real,    intent(inout) :: xbuf(:,:)  ! just need memory address
 integer, intent(in)    :: nprocs
 integer, intent(inout) :: irequestrecv(nprocs),nrecv
 logical, intent(in)    :: getdv,getdB,realviscosity
 integer :: n,j,iproc
 logical :: igotpart

 ! receive MPI broadcast
 do iproc=1,nprocs
    call MPI_TEST(irequestrecv(iproc),igotpart,status,mpierr)
    if (mpierr /= 0) call fatal('recv_density_results','error in MPI_TEST call')

    !
    ! unpack results
    ! there is no need for critical sections here as the density routine does not
    ! write to particles that are active on other processors, and particles are
    ! unique to each remote proc, so no danger of multiple threads writing
    ! to same particle
    !
    if (igotpart) then
       j = status(MPI_TAG)
       if (j <= 0 .or. j > npart) call fatal('recv_density_results','receive out of range',j)
       n = 1
       xyzh(4,j) = xbuf(1,iproc)
       if (maxgradh==maxp) then
          gradh(1,j) = real(xbuf(2,iproc),kind=kind(gradh))
#ifdef GRAVITY
          gradh(2,j) = real(xbuf(3,iproc),kind=kind(gradh))
#endif
          n = n + 2
       endif
       if (getdv) then
          divcurlv(:,j) = real(xbuf(n+1:n+ndivcurlv,iproc),kind=kind(divcurlv))
          n = n + ndivcurlv
       endif
       if (getdB) then
          divcurlB(:,j) = real(xbuf(n+1:n+ndivcurlB,iproc),kind=kind(divcurlB))
          n = n + ndivcurlB
       endif
       if (realviscosity .and. maxstrain==maxp) then
          straintensor(:,j) = real(xbuf(n+1:n+6,iproc),kind=kind(straintensor))
          n = n + 6
       endif
       nrecv = nrecv + 1
       !
       !--post another receive ready for next particle
       !
       call MPI_START(irequestrecv(iproc),mpierr)
    endif
 enddo

end subroutine recv_density_results

!------------------------------------------------
!+
!  As above but for results from force routine
!+
!------------------------------------------------
subroutine recv_force_results(nprocs,xbuf,irequestrecv,nrecv)
 use dim,  only:mhd,maxBevol,maxvxyzu,ndivcurlv,gravity
 use part, only:fxyzu,divcurlv,divBsymm,dBevol,poten,npart
 use io,   only:fatal
#ifdef IND_TIMESTEPS
 use part, only:ibin
#endif
 real,    intent(inout) :: xbuf(:,:)  ! just need memory address
 integer, intent(in)    :: nprocs
 integer, intent(inout) :: irequestrecv(nprocs),nrecv
 integer :: n,j,iproc
 logical :: igotpart

 ! receive MPI broadcast
 do iproc=1,nprocs
    call MPI_TEST(irequestrecv(iproc),igotpart,status,mpierr)
    if (mpierr /= 0) call fatal('recv_density_results','error in MPI_TEST call')

    ! unpack results
    if (igotpart) then
       j = status(MPI_TAG)
       if (j <= 0 .or. j > npart) call fatal('recv_force_results','receive out of range',j)
       fxyzu(1:maxvxyzu,j) = xbuf(1:maxvxyzu,iproc)
       n=maxvxyzu
#ifdef IND_TIMESTEPS
       ibin(j) = nint(xbuf(n+1,iproc),kind=kind(ibin))
       n = n + 1
#endif
       if (ndivcurlv >= 1) then
          divcurlv(1,j) = real(xbuf(n+1,iproc),kind=kind(divcurlv))
          n = n + 1
       endif
       if (mhd) then
          divBsymm(j) = real(xbuf(n+1,iproc),kind=kind(divBsymm))
          n = n + 1
          dBevol(:,j) = real(xbuf(n+1:n+maxBevol,iproc),kind=kind(dBevol))
          n = n + maxBevol
       endif
       if (gravity) then
          poten(j) = real(xbuf(n+1,iproc),kind=kind(poten))
          n = n + 1
       endif
       nrecv = nrecv + 1
       !
       !--post another receive ready for next particle
       !
       call MPI_START(irequestrecv(iproc),mpierr)
    endif
 enddo
end subroutine recv_force_results

!----------------------------------------------------------------
!+
!  finish/clean up of load balancing process.
!+
!----------------------------------------------------------------
subroutine finish_results_exchange(nprocs,xrecvbuf,irequestrecv,nrecv,nactive,getdv,getdB,realviscosity)
 use dim, only:maxp
 use io,  only:fatal
 integer, intent(in)    :: nprocs
 real,    intent(inout) :: xrecvbuf(:,:)
 integer, intent(inout) :: irequestrecv(nprocs),nrecv
 integer, intent(in)    :: nactive
 logical, intent(in), optional :: getdv,getdB,realviscosity
 integer :: newproc,iproc
 logical :: idone
 real    :: xsendbuf(1)

 do while (nrecv < nactive)
    !print*,id,' finishing ',nrecv
    if (present(getdv) .and. present(getdB)) then
       call recv_density_results(nprocs,xrecvbuf,irequestrecv,nrecv,getdv,getdB,realviscosity)
    else
       call recv_force_results(nprocs,xrecvbuf,irequestrecv,nrecv)
    endif
 enddo
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!
!--each processor do a dummy send to next processor to clear the last remaining receive
!  (we know the receive has been posted for this, so use RSEND)
!
 do newproc=0,nprocs-1
    call MPI_RSEND(xsendbuf,0,MPI_DEFAULT_REAL,newproc,0,MPI_COMM_WORLD,mpierr)
 enddo
!
!--sync all threads here
!
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

!--free request handle
 do iproc=1,nprocs
    call MPI_TEST(irequestrecv(iproc),idone,status,mpierr)
    call MPI_REQUEST_FREE(irequestrecv(iproc),mpierr)
    if (.not.idone) call fatal('finish_density_exchange_mpi','receive has not completed')
 enddo

 return
end subroutine finish_results_exchange

end module mpiderivs
#endif
