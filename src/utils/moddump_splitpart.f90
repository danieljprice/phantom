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
! :Dependencies: injectutils, io, part, splitpart
!
 implicit none
 integer            :: nchild = 12
 integer, parameter :: lattice_type = 0 ! 0 for lattice, 1 for random
 integer, parameter :: ires = 1         ! use 12 particles per sphere

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,    only:split_all_particles
 use io,           only:fatal,error
 use injectutils,  only:get_parts_per_sphere
 use part,         only:delete_dead_or_accreted_particles
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr

 ierr = 0

 !-- if using the regular grid, set nchild to get desired resolution
 if (lattice_type == 0) then
    nchild = get_parts_per_sphere(ires) + 1
 endif

 !-- don't split accreted particles
 call delete_dead_or_accreted_particles(npart,npartoftype)

 ! Split 'em!
 call split_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                                nchild,lattice_type,ires)

 print*,' new npart = ',npart

end subroutine modify_dump

end module moddump
