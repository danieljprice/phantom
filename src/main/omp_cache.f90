!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module omp_cache

 implicit none

 integer, parameter :: maxcellcache = 50000

 real, allocatable :: ompcache(:,:,:)

 public :: allocate_cache
 public :: deallocate_cache
 public :: ompcache
 public :: maxcellcache

 private

contains

subroutine allocate_cache
 use allocutils, only:allocate_array
 integer :: OMP_GET_MAX_THREADS, nthread

 nthread = OMP_GET_MAX_THREADS()
 call allocate_array('ompcache', ompcache, maxcellcache, 4, nthread)

end subroutine allocate_cache

subroutine deallocate_cache

 if (allocated(ompcache))     deallocate(ompcache)

end subroutine deallocate_cache

end module omp_cache
