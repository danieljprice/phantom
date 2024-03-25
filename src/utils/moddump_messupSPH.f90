!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: does not make any modifications
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: mess_up_SPH, part, prompting, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use mess_up_SPH ! module from MCFOST
 use part,      only:xyzmh_ptmass,nptmass,kill_particle,shuffle_part
 use prompting, only:prompt
 use units,     only:udist
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical, allocatable :: mask(:)
 integer :: ioption,i

 print*,'udist=',udist

 ioption = 1
 print "(3(/,a))",' 1) randomize azimuth ',&
                  ' 2) randomize gap ',&
                  ' 3) delete Hill sphere'

 call prompt('enter option ',ioption)
 allocate(mask(npart))
 mask = .false.
 select case(ioption)
 case(3)
    call mask_Hill_sphere(npart,nptmass,xyzh,xyzmh_ptmass,udist,mask)
 case(2)
    print*,' randomizing gap...'
    call randomize_gap(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,udist)
 case(1)
    print*,' randomizing azimuths...'
    call randomize_azimuth(npart,xyzh,vxyzu)
 case default
    print*,' no modifications performed '
 end select
 print*,' done ',count(mask)
 if (any(mask)) then
    do i=1,npart
       if (mask(i)) call kill_particle(i,npartoftype)
    enddo
    call shuffle_part(npart)
    print*,' NEW NUMBER OF PARTICLES = ',npart
 endif
 deallocate(mask)

end subroutine modify_dump

end module moddump

