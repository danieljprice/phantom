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

contains
!----------------------------------------------------------------
!+
!  info routine
!+
!----------------------------------------------------------------
subroutine info_omp
!$ integer, external :: omp_get_num_threads

!$omp parallel
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
!$ integer :: i
!$ external :: omp_init_lock
!$ integer, external :: omp_get_num_threads

 omp_num_threads = 1

!$ do i = 0, nlocks
!$  call omp_init_lock(ipart_omp_lock(i))
!$ enddo

!$omp parallel
!$ omp_num_threads = omp_get_num_threads()
!$omp end parallel

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

end module omputils
