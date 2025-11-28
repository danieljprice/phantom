!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
 character(len=*), parameter, public :: moddump_flags = ''

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
 integer, allocatable :: mask(:)
 real(kind=8) :: factor
 integer :: ioption,i
 logical :: inside

 print*,'udist=',udist

 ioption = 1
 print "(3(/,a))",' 1) randomize azimuth ',&
                  ' 2) randomize gap ',&
                  ' 3) delete Hill sphere'

 call prompt('enter option ',ioption)
 allocate(mask(npart))
 mask = 0
 select case(ioption)
 case(3)
    call mask_Hill_sphere(npart,nptmass,xyzh,xyzmh_ptmass,udist,mask)
 case(2)
    print*,' randomizing gap...'
    factor = 1.d0
    inside = .false.
    call prompt('enter factor ',factor)
    call prompt('do you want to randomize inside or outside the Hill sphere? ',inside)
    call randomize_gap(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,udist,factor,inside)
 case(1)
    print*,' randomizing azimuths...'
    call randomize_azimuth(npart,xyzh,vxyzu)
 case default
    print*,' no modifications performed '
 end select
 print*,' done ',count(mask==1),' particles killed'
 if (any(mask==1)) then
    do i=1,npart
       if (mask(i) == 1) call kill_particle(i,npartoftype)
    enddo
    call shuffle_part(npart)
    print*,' NEW NUMBER OF PARTICLES = ',npart
 endif
 deallocate(mask)

end subroutine modify_dump

end module moddump

