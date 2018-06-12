!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testmath
!
!  DESCRIPTION:
!  This module performs unit tests of the fast sqrt routines
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: fastmath, io, mpiutils, random
!+
!--------------------------------------------------------------------------
module testmath
 implicit none
 public :: test_math

 private

contains

subroutine test_math(ntests,npass,usefsqrt,usefinvsqrt)
 use io,       only:id,master
 use fastmath, only:checksqrt,testsqrt,finvsqrt
 use random,   only:ran2
 use mpiutils, only:barrier_mpi
 integer, intent(inout) :: ntests,npass
 logical, intent(out)   :: usefsqrt,usefinvsqrt
 integer, parameter :: n = 1000000
 integer, parameter :: stderr = 0
 integer            :: ierr,iseed,nerr,i
 real, allocatable :: x(:),f(:),y(:)
 real :: t1,t2,errmax,tnative

 if (id==master) write(stderr,"(a,/)") '--> TESTING FAST SQRT ROUTINES'

 usefsqrt    = .true.
 usefinvsqrt = .true.
!
!--check for errors first
!
 call testsqrt(ierr,output=.false.)
 if (ierr /= 0) then
    ! report errors on any threads
    write(stderr, "(a)") ' *** ERROR with fast sqrt on this architecture ***'
    usefsqrt    = .false.
    usefinvsqrt = .false.
    write(stderr,"(/,a,/)") '<-- FAST SQRT TEST FAILED'
    return
 endif

 allocate(x(n),f(n),y(n),stat=ierr)
 if (ierr /= 0) return

 iseed = -5234
 do i=1,n
    x(i) = ran2(iseed)*1.e8
 enddo

 ntests = ntests + 1
 nerr = 0
 do i=1,n
    call checksqrt(x(i),5.e-7*x(i),ierr,.false.)
    nerr = max(ierr,nerr)
 enddo
 if (nerr > 0) then
    usefsqrt    = .false.
    usefinvsqrt = .false.
 else
    npass = npass + 1
 endif

!
!--check timings for inverse sqrt
!
 call cpu_time(t1)
 do i=1,n
    f(i) = 1./sqrt(x(i))
 enddo
 call cpu_time(t2)
 tnative = t2-t1
 if (id==master) write(stderr,"(a,es10.3,a)") ' native 1/sqrt done in ',tnative,' cpu-s'
 y = f

 call barrier_mpi

 call cpu_time(t1)
 do i=1,n
    f(i) = finvsqrt(x(i))
 enddo
 call cpu_time(t2)

 ! run tests on all threads, but only report detailed results on master thread
 if (id==master) write(stderr,"(a,es10.3,a)") '   fast 1/sqrt done in ',t2-t1,' cpu-s'

 if ((t2-t1) > tnative) then
    if (id==master) write(stderr,"(a,f4.1)") ' so finvsqrt(x) is SLOWER than 1/sqrt(x) by factor of ',&
                               (t2-t1)/tnative
    usefinvsqrt = .false.
 else
    if (id==master) write(stderr,"(a,f4.1)") ' so finvsqrt(x) is FASTER than 1/sqrt(x) by factor of ', &
                                tnative/(t2-t1)
 endif

 errmax = 0.
 do i=1,n
    errmax = max(errmax,abs(y(i) - f(i))/y(i))
 enddo
 if (id==master) write(stderr,"(1x,a,es10.3)") 'max relative error is ',errmax
 if (errmax > 1.e-7) usefinvsqrt = .false.

 if (id==master) write(stderr,*)
 call barrier_mpi
!
!--check timings for sqrt
!
 call cpu_time(t1)
 do i=1,n
    f(i) = sqrt(x(i))
 enddo
 call cpu_time(t2)
 tnative = t2-t1
 if (id==master) write(stderr,"(a,es10.3,a)") '   native sqrt done in ',tnative,' cpu-s'
 y = f
 call barrier_mpi

 call cpu_time(t1)
 do i=1,n
    f(i) = x(i)*finvsqrt(x(i))
 enddo
 call cpu_time(t2)
 if (id==master) write(stderr,"(a,es10.3,a)") ' x*finvsqrt(x) done in ',t2-t1,' cpu-s'

 if ((t2-t1) > tnative) then
    if (id==master) write(stderr,"(a,f4.1)") ' so x*finvsqrt(x) is SLOWER than sqrt(x) by factor of ',&
                             (t2-t1)/tnative
    usefsqrt = .false.
 else
    if (id==master) write(stderr,"(a,f4.1)") ' so x*finvsqrt(x) is FASTER than sqrt(x) by factor of ',tnative/(t2-t1)
 endif

 errmax = 0.
 do i=1,n
    errmax = max(errmax,abs(y(i) - f(i))/(y(i) + epsilon(y)))
 enddo
 if (id==master) write(stderr,"(1x,a,es10.3)") 'max relative error is ',errmax
 if (errmax > 1.e-7) usefinvsqrt = .false.

 if (allocated(x)) deallocate(x)
 if (allocated(f)) deallocate(f)
 if (allocated(y)) deallocate(y)

 if (id==master) write(stderr,"(/,a,/)") '<-- FAST SQRT TEST COMPLETE'
 call barrier_mpi

end subroutine test_math

end module testmath
