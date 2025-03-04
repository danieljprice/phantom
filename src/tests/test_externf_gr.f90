!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testexternf
!
! Unit tests of the externalforces module
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: io
!
 implicit none
 public :: test_externf

 private

contains
!----------------------------------------------------------
!+
!  unit tests of external forces
!+
!----------------------------------------------------------
subroutine test_externf(ntests,npass)
 use io,       only:id,master
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> SKIPPING EXTERNAL FORCES TESTS FOR GR'

end subroutine test_externf

end module testexternf
