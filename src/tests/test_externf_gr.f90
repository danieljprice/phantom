!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testexternf
!
!  DESCRIPTION:
!  Unit tests of the externalforces module
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io
!+
!--------------------------------------------------------------------------
module testexternf
 implicit none
 public :: test_externf

 private

contains

subroutine test_externf(ntests,npass)
! use dim,      only:maxp
 use io,       only:id,master
! use part,     only:npart,xyzh,hfact,massoftype,igas
! use testutils,only:checkval,checkvalf,checkvalbuf_start,checkvalbuf,checkvalbuf_end
! use unifdis,  only:set_unifdis
! use units,    only:set_units
! use physcon,  only:pc,solarm
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING EXTERNAL FORCES MODULE'

 if (id==master) write(*,"(/,a)") '<-- EXTERNAL FORCE TESTS COMPLETE'

end subroutine test_externf

end module testexternf
