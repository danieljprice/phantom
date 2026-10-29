!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module omputils
!
! Utility subroutines specific to openMP parallelisation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!

!$ use dim, only:maxptmass
 implicit none
!$ integer, parameter :: nlocks = maxptmass
!$ integer(kind=8), dimension(0:nlocks) :: ipart_omp_lock

 integer :: omp_num_threads

 !--ids of the loops whose split between threads is tuned (see limits_omp)
 integer, parameter :: id1 = 1, id2 = 2, id3 = 3, id4 = 4, id5 = 5, id6 = 6
 integer, parameter, private :: nloops = 6
 real(kind=8), allocatable, private :: share(:,:)  ! cumulative share of each loop, up to each thread
 real(kind=8), allocatable, private :: tstart(:,:),tused(:,:)
 integer,      allocatable, private :: nused(:,:)
 integer, private :: ndone(nloops)

contains
!----------------------------------------------------------------
!+
!  info routine
!+
!----------------------------------------------------------------
subroutine info_omp
!$ integer, external :: omp_get_num_threads

!$omp parallel default(none)
!$omp single
!$ print "(a,i4,a)",' Running in openMP on',omp_get_num_threads(),' threads'
!$omp end single
!$omp end parallel

!$ if (.false.) then
    print "(a)",' openMP parallelisation is OFF'
!$ endif

end subroutine info_omp

!----------------------------------------------------------------
!+
!  initialisation routine necessary if locks are used
!+
!----------------------------------------------------------------
subroutine init_omp
 integer :: i
!$ external :: omp_init_lock
!$ integer, external :: omp_get_num_threads

 omp_num_threads = 1

!$ do i = 0, nlocks
!$  call omp_init_lock(ipart_omp_lock(i))
!$ enddo

!$omp parallel default(none) shared(omp_num_threads)
!$ omp_num_threads = omp_get_num_threads()
!$omp end parallel

 ! start with equal shares, as for schedule(static)
 if (allocated(share)) deallocate(share,tstart,tused,nused)
 allocate(share(0:omp_num_threads,nloops),tstart(0:omp_num_threads-1,nloops))
 allocate(tused(0:omp_num_threads-1,nloops),nused(0:omp_num_threads-1,nloops))
 do i = 0,omp_num_threads
    share(i,:) = i/real(omp_num_threads,kind=8)
 enddo
 tused = 0.
 nused = 0
 ndone = 0

end subroutine init_omp

!----------------------------------------------------------------
!+
!  interface to omp_get_num_threads that works even if
!  openMP is off
!+
!----------------------------------------------------------------
integer function omp_thread_num()
!$ integer, external :: omp_get_thread_num

 omp_thread_num = 0
!$ omp_thread_num = omp_get_thread_num()

end function omp_thread_num

!----------------------------------------------------------------
!+
!  range of the loop n1..n2 done by this thread, a contiguous range
!  as for schedule(static), so that each thread works on the same
!  particles every time, keeping them in its cache and memory.
!  For the loops given by iloop, the share of each thread is tuned
!  so that threads on slower cores, or with less memory bandwidth,
!  get fewer particles; call end_limits_omp after the loop
!+
!----------------------------------------------------------------
subroutine limits_omp(n1,n2,i1,i2,iloop)
 integer, intent(in)  :: n1,n2
 integer, intent(out) :: i1,i2
 integer, intent(in), optional :: iloop
 integer :: it,nt,n
!$ integer, external :: omp_get_num_threads
!$ real(kind=8), external :: omp_get_wtime

 it = omp_thread_num()
 nt = 1
!$ nt = omp_get_num_threads()
 n  = n2 - n1 + 1
 if (present(iloop) .and. nt > 1 .and. allocated(share)) then
    if (size(share,1) == nt+1) then
       i1 = n1 + int(n*share(it,iloop))
       i2 = n1 - 1 + int(n*share(it+1,iloop))
       if (it == nt-1) i2 = n2
       nused(it,iloop) = i2 - i1 + 1
!$     tstart(it,iloop) = omp_get_wtime()
       return
    endif
 endif
 i1 = n1 + int((int(n,kind=8)*it)/nt)
 i2 = n1 - 1 + int((int(n,kind=8)*(it+1))/nt)

end subroutine limits_omp

!----------------------------------------------------------------
!+
!  record the time this thread took for its share of a tuned loop;
!  the last thread to finish then updates the shares
!+
!----------------------------------------------------------------
subroutine end_limits_omp(iloop)
 integer, intent(in) :: iloop
 integer :: it,nt,nfinished
!$ integer, external :: omp_get_num_threads
!$ real(kind=8), external :: omp_get_wtime

 nt = 1
!$ nt = omp_get_num_threads()
 if (nt <= 1 .or. .not.allocated(share)) return
 if (size(share,1) /= nt+1) return

 it = omp_thread_num()
!$ tused(it,iloop) = omp_get_wtime() - tstart(it,iloop)
!$omp flush
!$omp atomic capture
 ndone(iloop) = ndone(iloop) + 1
 nfinished = ndone(iloop)
!$omp end atomic
 if (nfinished == nt) then
!$omp flush
    ndone(iloop) = 0
    call update_shares(iloop,nt)
 endif

end subroutine end_limits_omp

!----------------------------------------------------------------
!+
!  move each thread's share of a loop towards the rate at which
!  it got through its particles, averaged over many calls
!+
!----------------------------------------------------------------
subroutine update_shares(iloop,nt)
 integer, intent(in) :: iloop,nt
 real(kind=8), parameter :: tmin = 1.d-5, weight = 0.2
 real(kind=8) :: rate(0:nt-1),newshare(0:nt-1)
 integer :: it

 ! loops too short to time reliably keep their shares
 if (any(tused(0:nt-1,iloop) < tmin)) return

 rate     = max(nused(0:nt-1,iloop),1)/tused(0:nt-1,iloop)
 newshare = share(1:nt,iloop) - share(0:nt-1,iloop)
 newshare = (1. - weight)*newshare + weight*rate/sum(rate)
 do it = 1,nt-1
    share(it,iloop) = share(it-1,iloop) + newshare(it-1)/sum(newshare)
 enddo
 share(nt,iloop) = 1.

end subroutine update_shares

end module omputils
