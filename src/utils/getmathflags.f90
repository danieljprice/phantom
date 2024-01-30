!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program getmathflags
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: getmathflags [no arguments]
!
! :Dependencies: mpiutils, testmath
!
 use testmath, only:test_math
 use mpiutils, only:init_mpi, finalise_mpi
 implicit none
 integer :: ntests,npass,id,nprocs
 logical :: usefsqrt,usefinvsqrt

 ntests = 0
 npass  = 0
 call init_mpi(id,nprocs)
 call test_math(ntests,npass,usefsqrt,usefinvsqrt)
 if (usefinvsqrt) then
    print "(a)",'yes'
 endif
 call finalise_mpi()

end program getmathflags
