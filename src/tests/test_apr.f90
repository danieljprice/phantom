!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testapr
!
! Unit test for adaptive particle refinement
!
! :References:
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: apr, apr_region, linklist
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_apr,setup_apr_region_for_test

 private

contains

!--------------------------------------------
!+
!  Various tests of the apr module
!+
!--------------------------------------------
subroutine test_apr(ntests,npass)
 use physcon, only:solarm,kpc
 use units,   only:set_units
 integer, intent(inout) :: ntests,npass
 !integer :: nfailed(10),ierr,iregime

 if (id==master) write(*,"(/,a)") '--> ADDING APR TO TEST'

 if (id==master) write(*,"(/,a)") '<-- APR TEST ADDITION COMPLETE'

end subroutine test_apr

!--------------------------------------------
!+
!  Set up an APR region that is used in other tests
!+
!--------------------------------------------
subroutine setup_apr_region_for_test()
 use apr,  only:init_apr,update_apr,apr_max_in,ref_dir
 use apr,  only:apr_type,apr_rad
 use part, only:npart,xyzh,vxyzu,fxyzu,apr_level
 use linklist, only:set_linklist
 !real :: ratesq(nrates)
 integer :: ierr

 if (id==master) write(*,"(/,a)") '--> adding an apr region'

 ! set parameters for the region
  apr_max_in  =   1    ! number of additional refinement levels (3 -> 2x resolution)
  ref_dir     =   1     ! increase (1) or decrease (-1) resolution
  apr_type    =  -1     ! choose this so you get the default option which is
                        ! reserved for the test suite
  apr_rad     =   0.25  ! radius of innermost region


 ! initialise
 call init_apr(apr_level,ierr)
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)

end subroutine setup_apr_region_for_test

end module testapr
