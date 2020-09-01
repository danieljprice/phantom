!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! None
!
! :References: None
!
! :Owner: Enrico Ragusa
!
! :Runtime parameters: None
!
! :Dependencies: part
!

 use part, only:xyzmh_ptmass

 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 implicit none
 integer, intent(inout) :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout) :: massoftype
 real, dimension(:,:), intent(inout) :: xyzh,vxyzu
 !integer :: i
 real :: radius

 print*,' Phantommoddump: Change sink radii '

 radius=0.1
 !
 !--------Changing sink radius
 !
 print*,'Changing sink radii to',radius,' ...'
 !xyzmh_ptmass(5,1)=radius
 xyzmh_ptmass(5,2)=radius

 return
end subroutine modify_dump

end module moddump

