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
 integer            :: nchild = 12
 integer, parameter :: lattice_type = 0 ! 0 for lattice, 1 for random
 integer, parameter :: ires = 1         ! use 12 particles per sphere

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,    only:split_particles
 use io,           only:fatal,error
 use injectutils,  only:get_parts_per_sphere
 use part,         only:igas,set_particle_type
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, allocatable, dimension(:,:) :: xyzh_child,vxyzu_child
 integer :: ierr,ichild,iold,inew,j

 ierr = 0

 !-- if using the regular grid, set nchild to get desired resolution
 if (lattice_type == 0) then
   nchild = get_parts_per_sphere(ires) + 1
 endif

 allocate(xyzh_child(4,nchild),vxyzu_child(3,nchild))

 !--check there is enough memory
 if (size(xyzh(1,:)) < npart*nchild) then
    call error('split_particles','not enough memory, increase MAXP and recompile')
    ierr = 1
    return
 endif

 !--update npartoftype
 npartoftype(:) = npartoftype*nchild

 !--find positions of the new particles
 inew = npart ! put new particles after original parts
 ichild = npart !to keep track of the kids

 do iold=1,npart
    ! send in the parent, children return
    call split_particles(nchild,iold,xyzh(1:4,iold),vxyzu(1:3,iold), &
                         lattice_type,ires,xyzh_child,vxyzu_child)

   ! copy children over, first replaces the parent
   xyzh(1:4,iold) = xyzh_child(1:4,1)
   vxyzu(1:3,iold) = vxyzu_child(1:3,1)
   do j = 2,nchild
     xyzh(1:4,ichild+j-1) = xyzh_child(1:4,j)
     vxyzu(1:4,ichild+j-1) = vxyzu_child(1:4,j)
     call set_particle_type(ichild+j-1,igas)
   enddo
   ichild = ichild + (nchild-1)
 enddo

 !-- new npart
 npart = npart * nchild

 !--new masses
 massoftype(:) = massoftype(:)/nchild

 !--tidy up
 deallocate(xyzh_child,vxyzu_child)

 if (ierr /= 0) call fatal('moddump','could not split particles')
 print*,' got npart = ',npart

end subroutine modify_dump

end module moddump
