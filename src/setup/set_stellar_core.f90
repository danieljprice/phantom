!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setstellarcore
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Mike Lau
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module setstellarcore
 implicit none
contains
!-----------------------------------------------------------------------
!+
!  Add a sink particle as a stellar core
!+
!-----------------------------------------------------------------------
subroutine set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,mcore,hsoft,ihsoft)
 integer, intent(out) :: nptmass
 integer, intent(in)  :: ihsoft
 real, intent(out)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(in)     :: mcore,hsoft
 integer              :: n
 nptmass                = 1
 n                      = nptmass
 xyzmh_ptmass(:,n)      = 0. ! zero all quantities by default
 xyzmh_ptmass(4,n)      = mcore
 xyzmh_ptmass(ihsoft,n) = hsoft
 vxyz_ptmass(:,n)       = 0.
end subroutine set_stellar_core
end module setstellarcore
