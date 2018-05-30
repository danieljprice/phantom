!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: omputils
!
!  DESCRIPTION:
!  Utility subroutines specific to openMP parallelisation
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module omputils
!$ use dim, only:maxp_hard,maxptmass
!$ integer, parameter :: nlockgrp = 10
!$ integer, parameter :: nlocks = max(maxp_hard/nlockgrp,maxptmass)
!$ integer(kind=8), dimension(0:nlocks) :: ipart_omp_lock

contains

!----------------------------------------------------------------
!+
!  info routine
!+
!----------------------------------------------------------------
subroutine info_omp
#ifdef _OPENMP
 integer omp_get_num_threads

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
!$ integer :: i

!$ do i = 0, nlocks
!$    call omp_init_lock(ipart_omp_lock(i))
!$ enddo

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
 integer :: omp_get_num_threads, omp_get_thread_num
 logical :: omp_in_parallel

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

end module omputils
