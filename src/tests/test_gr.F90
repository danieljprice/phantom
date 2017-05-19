!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testindtstep
!
!  DESCRIPTION:
!  test module for individual timestepping utilities
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!+
!--------------------------------------------------------------------------
module testgr
 implicit none
 public :: test_gr

 private

contains

subroutine test_gr(ntests,npass)
 use io,              only:id,master
 use testutils,       only:checkvalbuf,checkvalbuf_end,checkval
 use testmetric,      only:test_metric
 use testcons2prim,   only:test_cons2prim
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING GENERAL RELATIVITY'
 call test_metric(ntests,npass)
 call test_cons2prim(ntests,npass)
 if (id==master) write(*,"(/,a)") '<-- GR TESTS COMPLETE'

end subroutine test_gr

end module testgr
