!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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

contains
!----------------------------------------------------------------
!+
!  info routine
!+
!----------------------------------------------------------------
subroutine info_omp
#ifdef _OPENMP
 integer, external :: omp_get_num_threads

!$omp parallel
!$omp master
!$ print "(a,i4,a)",' Running in openMP on',omp_get_num_threads(),' threads'
!$omp end master
!$omp end parallel

#else

 print "(a)",' openMP parallelisation is OFF'

#endif

end subroutine info_omp

!----------------------------------------------------------------
!+
!  initialisation routine necessary if locks are used
!+
!----------------------------------------------------------------
subroutine init_omp
#ifdef _OPENMP
!$ integer :: i
!$ external :: omp_init_lock
 integer, external :: omp_get_num_threads

!$ do i = 0, nlocks
!$  call omp_init_lock(ipart_omp_lock(i))
!$ enddo

!$omp parallel
 omp_num_threads = omp_get_num_threads()
!$omp end parallel
#else
 omp_num_threads = 1
#endif

end subroutine init_omp

!----------------------------------------------------------------
!+
! routine to compute OpenMP loop limits
! required by some analysis routines
! originally by Aake Nordlund
!+
!----------------------------------------------------------------
subroutine limits_omp (n1,n2,i1,i2)
 integer, intent(in)  :: n1,n2
 integer, intent(out) :: i1,i2
#ifdef _OPENMP
 integer, external :: omp_get_num_threads, omp_get_thread_num
 logical, external :: omp_in_parallel

 if (omp_in_parallel()) then
    i1 = n1 + ((omp_get_thread_num()  )*n2)/omp_get_num_threads()
    i2 =      ((omp_get_thread_num()+1)*n2)/omp_get_num_threads()
 else
    i1 = max(1,n1)
    i2 = n2
 endif
#else
 i1 = max(1,n1)
 i2 = n2
#endif
end subroutine limits_omp

!----------------------------------------------------------------
!+
! routine to compute OpenMP loop limits
! when the work per iteration is known
!+
!----------------------------------------------------------------
subroutine limits_omp_work (n1,n2,i1,i2,work,mask,iskip)
 integer, intent(in)  :: n1,n2
 integer, intent(out) :: i1,i2,iskip
 real, intent(in) :: work(n2)
 integer, intent(in) :: mask(n2)

#ifdef _OPENMP
 integer, external :: omp_get_num_threads, omp_get_thread_num
 integer :: num_threads,id
 real :: chunk,my_chunk
 integer :: my_thread,i

 num_threads = omp_get_num_threads()
 id = omp_get_thread_num()
 iskip = 1

 chunk = sum(work,mask=(mask>0))/num_threads
 if (chunk < epsilon(0.)) then
    ! default to static scheduling
    !call limits_omp(n1,n2,i1,i2)
    i1 = 1 + id
    i2 = n2
    iskip = num_threads
    return
 endif
 i1 = 1
 i2 = 0
 my_chunk = 0.
 my_thread = -1
 do i=n1,n2
    if (mask(i) > 0) my_chunk = my_chunk + work(mask(i))
    if (my_chunk >= chunk) then
       my_thread = my_thread + 1
       if (my_thread == id) then
          i2 = i
          exit
       else
          i1 = i+1
          my_chunk = 0.
       endif
    endif
 enddo
 if (id==num_threads-1) i2 = n2
#else
 i1 = max(1,n1)
 i2 = n2
 iskip = 1
#endif
 !print*,'thread ',id,' limits  = ',i1,i2,my_chunk,' out of ',n1,n2,chunk*num_threads

end subroutine limits_omp_work

integer function omp_thread_num()
#ifdef _OPENMP
 integer, external :: omp_get_thread_num
 omp_thread_num = omp_get_thread_num()
#else
 omp_thread_num = 0
#endif
end function omp_thread_num

end module omputils
