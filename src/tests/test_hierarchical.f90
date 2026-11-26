!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testsethier
!
! Unit tests of heirarchical setup module
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, sethier_utils, testutils
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

 if (id==master) write(*,"(/,a,/)") '--> TESTING HIERARCHICAL SETUP'

 call test_heirarchical_string(ntests,npass)

! call test_heirarchical_setup(ntests,npass)

 !call test_readwrite_heirarchy(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- CHESS SETUP TESTS COMPLETE'

end subroutine test_sethier

!-----------------------------------------------------------------------
!+
!   Test conversion between hierarchies given with brackets, e.g. (o-o)-o,
!   and the numbered labels of the stars, e.g. 111,112,12, both ways
!+
!-----------------------------------------------------------------------
subroutine test_heirarchical_string(ntests,npass)
 use sethier_utils, only:hierarchy_from_brackets,hierarchy_to_brackets,lenhierstring,print_hierarchy
 use testutils,     only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ncases = 7, nbad = 4
 character(len=*), parameter :: brackets(ncases) = &
    (/'o                  ','o-o                ','(o-o)-o            ','o-(o-o)            ', &
      '(o-o)-(o-o)        ','((o-o)-o)-o        ','(o-(o-o))-((o-o)-o)'/)
 character(len=*), parameter :: labels(ncases) = &
    (/'1                           ','11,12                       ','111,112,12                  ', &
      '11,121,122                  ','111,112,121,122             ','1111,1112,112,12            ', &
      '111,1121,1122,1211,1212,122 '/)
 ! strings that are not valid hierarchies
 character(len=*), parameter :: bad_brackets(nbad) = (/'(o-o      ','o-o-o     ','(o)-o     ','o-x       '/)
 character(len=*), parameter :: bad_labels(nbad)   = (/'11,121    ','11,12,121 ','11,13     ','21,22     '/)
 character(len=lenhierstring) :: str
 integer :: i,ierr,nfailed(2)

 do i=1,ncases
    nfailed = 0
    call hierarchy_from_brackets(trim(brackets(i)),str,ierr)
    call checkval(ierr,0,0,nfailed(1),'convert '//trim(brackets(i)))
    call checkval(trim(str),trim(labels(i)),nfailed(2),trim(brackets(i))//' gives '//trim(labels(i)))
    call update_test_scores(ntests,nfailed,npass)

    nfailed = 0
    call hierarchy_to_brackets(trim(labels(i)),str,ierr)
    call checkval(ierr,0,0,nfailed(1),'convert '//trim(labels(i)))
    call checkval(trim(str),trim(brackets(i)),nfailed(2),trim(labels(i))//' gives '//trim(brackets(i)))
    call update_test_scores(ntests,nfailed,npass)
 enddo

 ! the brackets list the stars from left to right, whatever the order of the labels
 nfailed = 0
 call hierarchy_to_brackets('111,12,112',str,ierr)
 call checkval(ierr,0,0,nfailed(1),'convert 111,12,112')
 call checkval(trim(str),'(o-o)-o',nfailed(2),'111,12,112 gives (o-o)-o')
 call update_test_scores(ntests,nfailed,npass)

 ! print the tree of a quadruple system (output only, not checked)
 call print_hierarchy(trim(labels(ncases)))

 do i=1,nbad
    nfailed = 0
    call hierarchy_from_brackets(trim(bad_brackets(i)),str,ierr)
    call checkval(ierr /= 0,.true.,nfailed(1),'error for '//trim(bad_brackets(i)))
    call hierarchy_to_brackets(trim(bad_labels(i)),str,ierr)
    call checkval(ierr /= 0,.true.,nfailed(2),'error for '//trim(bad_labels(i)))
    call update_test_scores(ntests,nfailed,npass)
 enddo

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
