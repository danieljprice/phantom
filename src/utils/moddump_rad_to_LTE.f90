!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Convert radiation dump to LTE dump (ieos=12)
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, io, part
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim, only:do_radiation
 use io,  only:fatal
 use eos, only:ieos,gmw
 use part,only:rad,iradxi
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i

 if (.not. do_radiation) call fatal("moddump_rad_to_LTE","Not compiled with radiation")
 do i=1,npart
    if (isnan(rad(iradxi,i))) call fatal("moddump_rad_to_LTE","rad array contains NaNs")
    vxyzu(4,i) = vxyzu(4,i) + rad(iradxi,i)
 enddo

 ieos = 12
 gmw = 0.6   ! CHANGE MU HERE for writing into infile
 print*,'mu has been changed to',gmw  ! mu should not change from what was assumed with radiation

end subroutine modify_dump

end module moddump

