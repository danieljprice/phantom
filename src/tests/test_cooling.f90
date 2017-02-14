!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testcooling
!
!  DESCRIPTION:
!   Unit tests of the cooling module
!
!  REFERENCES:
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: cooling, io, testutils
!+
!--------------------------------------------------------------------------
module testcooling
 use testutils, only:checkval
 use io,        only:id,master
 implicit none
 public :: test_cooling

 private

contains

!--------------------------------------------
!+
!  Various tests of the cooling module
!+
!--------------------------------------------
subroutine test_cooling(ntests,npass)
 integer, intent(inout) :: ntests,npass
 !integer :: nfailed(10),ierr,iregime

 if (id==master) write(*,"(/,a)") '--> TESTING COOLING MODULE'

 !call set_units(mass=solarm,dist=au,G=1.d0)

 call test_cooling_rate(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- COOLING TEST COMPLETE'

end subroutine test_cooling

!--------------------------------------------
!+
!  Check that cooling function is continuous
!+
!--------------------------------------------
subroutine test_cooling_rate(ntests,npass)
 use cooling, only:cooling_rate_sd93
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nt = 100
 real :: logtmin,logtmax,logt,dlogt,t,crate
 integer :: i

 if (id==master) write(*,"(/,a)") '--> testing SD93 cooling rate'

 logtmax = log10(1.e8)
 logtmin = log10(1.e4)

 dlogt = (logtmax - logtmin)/real(nt)
 do i=1,nt
    logt = logtmin + (i-1)*dlogt
    t = 10**logt
    crate = cooling_rate_sd93(t)
    write(1,*) t,-crate
 enddo

end subroutine test_cooling_rate

end module testcooling
