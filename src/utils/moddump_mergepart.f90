!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  merges particles; input simulation with npart, get npart/nchild back
!
!  REFERENCES: Vacondio et al. 2013
!
!  OWNER: Rebecca Nealon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, splitpart
!+
!--------------------------------------------------------------------------
module moddump
 implicit none
 integer, parameter :: nchild = 4

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart, only:merge_particles
 use part,      only:igas,kill_particle,delete_dead_or_accreted_particles
 use io,        only:fatal,error
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr,nactive,nparent,remainder, i,j,k
 integer :: on_list(npart), children_list(nchild)
 integer :: ichild,iparent,child_found
 logical :: iactive(npart)
 real, allocatable, dimension(:,:) :: xyzh_parent,vxyzu_parent,xyzh_in,vxyzu_in
 real    :: mchild,rik(3),rik2,rikmax,mparent

ierr = 0

 !-- how many active particles (this routine currently ignores inactive particles)
 nactive = 0
 iactive = .false.
 do i=1,npart
   if (xyzh(4,i) > 0.) then
     nactive = nactive + 1
     iactive(i) = .true.
   endif
 enddo
 if (nactive < npart) print*,'Discarding inactive particles'

 !-- check number of parent particles
 mchild = massoftype(igas)/nchild
 nparent = floor(real(nactive)/real(nchild))
 remainder = mod(nactive,nchild)
 if (remainder/nactive > 0.01) then
   ! discarding a couple of particles is ok, just don't want it to be too many
    call error('merge_particles','need to merge evenly, make sure npart(active)/nchild ~ integer')
    ierr = 1
    return
  else
    ! ignore a couple
    do i = 1,remainder
      if (xyzh(4,npart-i) > 0.) iactive(npart-i) = .false.
    enddo
    print*,'Ignoring',remainder,'particles to get an even split'
 endif

 !-- arrays for new parent particles
 allocate(xyzh_parent(4,nparent),vxyzu_parent(3,nparent))
 allocate(xyzh_in(4,nchild),vxyzu_in(3,nchild))

 iparent = 0
 ichild = 0
 on_list = -1

over_parent: do i=1,npart
  !already on the list or inactive
  if (on_list(i) > 0 .or. .not.iactive(i)) cycle over_parent
  children_list = 0
  ichild = 0

  ! first child
  iparent = iparent + 1
  ichild = ichild + 1
  on_list(i) = 1
  children_list(ichild) = i

  ! find nearby children
  over_child: do j=1,nchild-1
    rikmax = huge(rikmax)
    !-- choose the next closest particle as child
    !-- (probably a much faster way to do this bit)
    child_found = -1
    over_neighbours: do k=1,npart
      if (on_list(k) > 0 .or. .not.iactive(k)) cycle over_neighbours
      rik = xyzh(1:3,k) - xyzh(1:3,i)
      rik2 = dot_product(rik,rik)
      if (rik2 < rikmax) then
        rikmax = rik2
        child_found = k
      endif
    enddo over_neighbours

    if (child_found < 0) then
      ! no children found for the parent particle
      call error('mergepart','no child found for parent particle')
      print*,'parent particle is',i
      ierr = 1
    endif

    ! if child found, save the child to the list
    ichild = ichild + 1
    children_list(ichild) = child_found
    on_list(child_found) = 1

  enddo over_child

! send in children, parent returns
call merge_particles(nchild,children_list,mchild, npart, &
     xyzh,vxyzu,xyzh_parent(1:4,iparent),vxyzu_parent(1:3,iparent),mparent)

enddo over_parent

!-- store results
  xyzh(1:4,1:nparent) = xyzh_parent(:,:)
  vxyzu(1:3,1:nparent) = vxyzu_parent(:,:)

!-- kill all the useless children
  do i=nparent+1,npart
    call kill_particle(i,npartoftype)
  enddo

!--tidy up
  call delete_dead_or_accreted_particles(npart,npartoftype)
  deallocate(xyzh_parent,vxyzu_parent,xyzh_in,vxyzu_in)

!--update npartoftype
  npartoftype(igas) = nparent
  npart = nparent
  massoftype(:) = mparent

  if (ierr /= 0) call fatal('moddump','could not merge particles')

end subroutine modify_dump

end module moddump
