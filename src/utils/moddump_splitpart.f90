!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
! :Dependencies: io, splitpart
!
 implicit none
 integer            :: nchild = 3
 integer, parameter :: lattice_type = 0 ! 0 for lattice, 1 for random
 integer, parameter :: ires = 1         ! use 12 particles per sphere

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,    only:split_particles
 use io,           only:fatal,error
 use injectutils,  only:get_parts_per_sphere
 use part,         only:igas,copy_particle
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr,ichild,iparent

 ierr = 0

 !-- if using the regular grid, set nchild to get desired resolution
 if (lattice_type == 0) then
   nchild = get_parts_per_sphere(ires) + 1
 endif

 !--check there is enough memory
 if (size(xyzh(1,:)) < npart*nchild) then
    call error('split_particles','not enough memory, increase MAXP and recompile')
    ierr = 1
    return
 endif

 !--update npartoftype
 npartoftype(:) = npartoftype*nchild

 !--find positions of the new particles
 
 ichild = npart !to keep track of the kids

 do iparent=1,npart
    ! send in the parent, children return
    ! (the parent acts as the first child, this routine generates nchild-1 new particles
    ! and adjusts the smoothing length on the parent)
    call split_particles(nchild,iparent,xyzh,vxyzu,lattice_type,ires,ichild)

    ! for next children
    ichild = ichild + nchild - 1
 enddo

 !-- new npart
 npart = npart * nchild

 !--new masses
 massoftype(:) = massoftype(:)/nchild


 if (ierr /= 0) call fatal('moddump','could not split particles')
 print*,' got npart = ',npart

end subroutine modify_dump

end module moddump
