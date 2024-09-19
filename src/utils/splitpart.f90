!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module splitpart
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: getneighbours, io, part, splitmergeutils, timestep_ind
!

 use splitmergeutils, only:split_a_particle,fancy_merge_into_a_particle

 implicit none

contains

subroutine split_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                                 nchild,lattice_type,ires)
 use io,    only:fatal,error
 use part,  only:igas,copy_particle
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in)    :: lattice_type,ires,nchild
 integer :: ierr,ichild,iparent

 ierr = 0

 !--check there is enough memory
 if (size(xyzh(1,:)) < npart*nchild) then
    call error('split_all_particles','not enough memory, increase MAXP and recompile')
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
    call split_a_particle(nchild,iparent,xyzh,vxyzu,lattice_type,ires,ichild)

    ! for next children
    ichild = ichild + nchild - 1
 enddo

 !-- new npart
 npart = npart * nchild

 !--new masses
 massoftype(:) = massoftype(:)/nchild

 if (ierr /= 0) call fatal('splitpart','could not split particles')

end subroutine split_all_particles

subroutine merge_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                                nchild,nactive_here,fancy_merging)
 use part,          only:igas,kill_particle,delete_dead_or_accreted_particles
 use part,          only:isdead_or_accreted,copy_particle
 use timestep_ind,  only:nactive
 use io,            only:fatal,error
 use getneighbours, only:generate_neighbour_lists,neighb,neighcount,neighmax
 integer, intent(inout) :: npart,npartoftype(:)
 integer, intent(in)    :: nchild
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical, optional, intent(in) :: fancy_merging
 integer, optional, intent(in) :: nactive_here
 integer :: ierr,nparent,remainder, i,k
 integer :: on_list(npart), children_list(nchild)
 integer :: neighbours(neighmax),neigh_count
 integer :: ichild,iparent,child_found,m
 real    :: rik(3),rik2,rikmax,r2_neighbours(neighmax)
 logical :: merge_stochastically,found

 ierr = 0

 !-- how should the particles be merged? stochastic is fastest but not *quite*
 !   as accurate as averaging properties over the children
 merge_stochastically = .true.
 if (present(fancy_merging)) then
    merge_stochastically = .false.
    print*,'Doing fancy merging, this takes a bit longer'
 endif

 !-- how many active particles? If called from moddump, it must be provided
 if (present(nactive_here)) nactive = nactive_here

 if (nactive < npart) then
    call delete_dead_or_accreted_particles(npart,npartoftype)
    print*,'Discarding inactive particles'
 endif

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
    print*,' ignoring',remainder,'particles to get an even split'
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
 neighbours = 0
 neigh_count = 0

 if (merge_stochastically) then
    !-- quick stochastic merging
    do i = 1,npart,nchild
       iparent = iparent + 1
       call copy_particle(i,npart+iparent,.true.)
       xyzh(4,npart+iparent) = xyzh(4,i) * (nchild)**(1./3.)
    enddo
 else
    !-- slower merging that averages properties of children to find new parent
    !-- find all the neighbours first
    call generate_neighbour_lists(xyzh,vxyzu,npart,'dummy',write_neighbour_list=.false.)

    !-- now find the children
    over_parent: do while (iparent < nparent) ! until all the children are found
       i = i + 1
       ! already on the list
       if (on_list(i) > 0) cycle over_parent

       ! first child
       iparent = iparent + 1
       ichild = 1
       on_list(i) = i
       children_list(ichild) = i
       neigh_count = neighcount(i)
       neighbours(:) = neighb(i,:)

       ! calculate distance to neighbours
       do k = 1, neigh_count
          m = neighbours(k)
          rik = xyzh(1:3,m) - xyzh(1:3,i)
          r2_neighbours(k) = dot_product(rik,rik)
       enddo

       finding_children: do while (ichild < nchild)
          found = .false.
          child_found = -1
          rikmax = huge(rikmax)

          !first check over the neighbours
          if (.not.found) then
             over_neighbours: do k=1, neigh_count
                m = neighbours(k)
                if (on_list(m) > 0) cycle over_neighbours
                if (r2_neighbours(k) < rikmax) then
                   rikmax = r2_neighbours(k)
                   child_found = m
                endif
             enddo over_neighbours
             if (child_found > 0) found = .true.
          endif

          ! otherwise, check everyone else
          if (.not.found) then
             over_all: do k = 1,npart
                if (on_list(k) > 0) cycle over_all
                rik = xyzh(1:3,k) - xyzh(1:3,i)
                rik2 = dot_product(rik,rik)
                if (rik2 < rikmax) then
                   rikmax = rik2
                   child_found = k
                endif
             enddo over_all
             if (child_found > 0) found = .true.
          endif

          ! error if no children found for the parent particle
          if (.not.found) then
             call error('mergepart','no child found for parent particle')
             ierr = 1
          endif

          ! save the child to the list
          ichild = ichild + 1
          children_list(ichild) = child_found
          on_list(child_found)  = child_found
       enddo finding_children


       ! send in children, parent returns
       ! parents temporarily stored after all the children
       call fancy_merge_into_a_particle(nchild,children_list, massoftype(igas), &
                               npart,xyzh,vxyzu,npart+iparent)

    enddo over_parent
 endif

 !-- move the new parents
 do i = 1,nparent
    call copy_particle(npart+i,i,.true.)
 enddo

 !-- kill all the useless children
 do i=nparent+1,npart
    call kill_particle(i,npartoftype)
 enddo

 !--tidy up
 call delete_dead_or_accreted_particles(npart,npartoftype)

 !--update npartoftype
 npartoftype(igas) = nparent
 npart = nparent
 massoftype(:) = massoftype(:) * nchild

 if (ierr /= 0) call fatal('moddump','could not merge particles')

end subroutine merge_all_particles

end module splitpart
