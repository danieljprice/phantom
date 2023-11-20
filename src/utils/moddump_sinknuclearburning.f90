!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! set up sink heating from nuclear burning
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: part, ptmass_heating, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only:xyzmh_ptmass,vxyz_ptmass,nptmass,ilum
 use ptmass_heating, only:Lnuc
 use units,          only:unit_energ,utime
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: massoftype(:),xyzh(:,:),vxyzu(:,:)
 real                   :: Lnuc_cgs
 integer                :: core_sink_id

 Lnuc_cgs = 1.672e38
 Lnuc = Lnuc_cgs / unit_energ * utime

 core_sink_id = 1
 xyzmh_ptmass(ilum,core_sink_id) = Lnuc

 return
end subroutine modify_dump

end module moddump
