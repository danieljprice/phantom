!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
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
!  DEPENDENCIES: coolfunc, cooling, io, testutils
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

 !call test_cooling_rate(ntests,npass)

 call test_coolfunc(ntests,npass)

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

!--------------------------------------------
!+
!  Check that cooling function is continuous
!+
!--------------------------------------------
subroutine test_coolfunc(ntests,npass)
 use coolfunc,  only:init_coolfunc,find_in_table
 use testutils, only:checkvalbuf,checkvalbuf_start
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nt = 100
 integer :: i,k,ndiff,ncheck
 real    :: table(nt),val
 logical :: my_test

 if (id==master) write(*,"(/,a)") '--> testing find_in_table routine'
 !
 ! set up table
 !
 do i=1,nt
    table(i) = i
 enddo

 ndiff = 0
 ncheck = 0
 ntests = ntests + 1
 call checkvalbuf_start('table(i) < val < table(i+1)')
 do i=1,nt-1
    val = i+0.5
    k = find_in_table(nt,table,val)
    !print*,k,table(k),val,table(k+1)
    if (k < nt) then
       my_test = (val > table(k) .and. val < table(k+1))
       call checkvalbuf(my_test,.true.,'table(i) < val < table(i+1)',ndiff,ncheck)
    endif
 enddo
 if (ndiff==0) npass = npass + 1

 !if (id==master) write(*,"(/,a)") '--> testing cooling tables'

end subroutine test_coolfunc

end module testcooling
