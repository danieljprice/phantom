!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: getmathflags
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: getmathflags [no arguments]
!
!  DEPENDENCIES: mpiutils, testmath
!+
!--------------------------------------------------------------------------
program getmathflags
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
