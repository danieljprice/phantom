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
