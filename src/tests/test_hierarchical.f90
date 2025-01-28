!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testsethier
!
! Unit tests of heirarchical setup module
!
! :References: None
!
! :Owner: Simone Ceppi
!
! :Runtime parameters: None
!
! :Dependencies: io
!
 implicit none
 public :: test_sethier

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of the set_heirarchical routine
!+
!-----------------------------------------------------------------------
subroutine test_sethier(ntests,npass)
 use io,         only:id,master
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING DISC SETUP'

 call test_heirarchical_string(ntests,npass)

! call test_heirarchical_setup(ntests,npass)

 !call test_readwrite_heirarchy(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- CHESS SETUP TESTS COMPLETE'


end subroutine test_sethier

subroutine test_heirarchical_string(ntests,npass)
 integer, intent(inout) :: ntests,npass

 print*,"test heirarchical string"
!  result = parse_string(¨111,112,113¨)
 ! call checkval(result,result_ref,"check 111,112,113")

end subroutine test_heirarchical_string

!subroutine test_heirarchical_setup

 ! call sethIERARCHICAL()

!
!--check that set_disc passes check_setup routine
!
!    call check_setup(nerr,nwarn)
!    call checkval(nerr,0,0,nfailed(1),'setup errors')
!    call update_test_scores(ntests,nfailed(1:1),npass)

!    call checkval(nwarn,0,0,nfailed(1),'setup warnings')
!    call update_test_scores(ntests,nfailed(1:1),npass)


 ! DO TESTS HERE
!end subroutine test_heirarchical_setup

end module testsethier
