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
 integer, parameter :: nchild = 10

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart, only:merge_particles
 use part,      only:igas,kill_particle,delete_dead_or_accreted_particles
 use part,      only:isdead_or_accreted,copy_particle
 use io,        only:fatal,error
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr,nparent,remainder, i,j,k
 integer :: on_list(npart), children_list(nchild)
 integer :: ichild,iparent,child_found,nactive
 real    :: rik(3),rik2,rikmax

ierr = 0

 !-- how many active particles (this routine currently ignores inactive particles)
 nactive = 0
 do i = 1,npart
   if (.not.isdead_or_accreted(xyzh(4,i))) then
     nactive = nactive + 1
   else
     call kill_particle(i,npartoftype)
   endif
 enddo

 print*,'nactive',nactive,'npart',npart

 if (nactive < npart) then
   call delete_dead_or_accreted_particles(npart,npartoftype)
   print*,'Discarding inactive particles'
 endif

 print*,'nactive',nactive,'npart',npart

 !-- check number of parent particles
 nparent = floor(real(nactive)/real(nchild))
 remainder = mod(nactive,nchild)
 if (remainder/nactive > 0.01) then
   ! discarding a couple of particles is ok, just don't want it to be too many
    call error('merge_particles','need to merge evenly, make sure npart(active)/nchild ~ integer')
    ierr = 1
    return
  elseif (remainder > 0) then
    ! we just forget about these later on
    print*,'Ignoring',remainder,'particles to get an even split'
 endif

 !--check there is enough memory
 if (size(xyzh(1,:)) < nparent + npart) then
    call error('merge_particles','not enough memory, increase MAXP and recompile')
    ierr = 1
    return
 endif

 iparent = 0
 ichild = 0
 on_list = -1
 children_list = 0
 i = 0

over_parent: do while (iparent < nparent) ! until all the children are found
  i = i + 1
  !already on the list
  if (on_list(i) > 0) cycle over_parent

  ! first child
  iparent = iparent + 1
  ichild = 1
  on_list(i) = i
  children_list(ichild) = i

  ! find nearby children
  over_child: do j=1,nchild-1
    !-- choose the next closest particle as child
    !-- (there *must* be a more accurate way to group them)
    child_found = -1
    rikmax = huge(rikmax)
    over_neighbours: do k=1,npart
      if (on_list(k) > 0) cycle over_neighbours
      rik = xyzh(1:3,k) - xyzh(1:3,i)
      rik2 = dot_product(rik,rik)
      if (rik2 < rikmax) then
        rikmax = rik2
        child_found = k
      endif
    enddo over_neighbours

    if (child_found > 0) then
      ! if child found, save the child to the list
      ichild = ichild + 1
      children_list(ichild) = child_found
      on_list(child_found) = k
      !print*,'parent',iparent,'paired with',child_found
    else
      ! no children found for the parent particle
      call error('mergepart','no child found for parent particle')
      print*,'parent particle is',iparent
      ierr = 1
    endif

  enddo over_child
  !read*

! send in children, parent returns
! parents temporarily stored after all the children
call merge_particles(nchild,children_list,massoftype(igas), npart, &
     xyzh,vxyzu,npart+iparent)

enddo over_parent

!-- move the new parents
do i = 1,nparent
  call copy_particle(npart+i,i)
enddo

!-- kill all the useless children
do i=nparent+1,npart
  call kill_particle(i)
enddo

!--tidy up
call delete_dead_or_accreted_particles(npart,npartoftype)

!--update npartoftype
  npartoftype(igas) = nparent
  npart = nparent
  massoftype(:) = massoftype(:) * nchild

  if (ierr /= 0) call fatal('moddump','could not merge particles')

end subroutine modify_dump

end module moddump
