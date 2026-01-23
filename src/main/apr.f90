!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr
!
! Everything needed for live adaptive particle refinement
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: apr_region, cons2primsolver, dim, eos, extern_gr,
!   externalforces, get_apr_level, io, io_summary, kdtree, metric_tools,
!   mpiforce, neighkdtree, options, part, physcon, quitdump, random,
!   relaxem, timestep_ind, utils_apr, vectorutils
!
 use dim, only:gr,use_apr
 use apr_region
 use utils_apr

 implicit none

 public :: init_apr,update_apr
 public :: use_apr

 private
 real    :: sep_factor = 0.2
 logical :: apr_verbose = .false.
 logical :: do_relax = .false.
 logical :: adjusted_split = .true.

contains

!-----------------------------------------------------------------------
!+
!  Initialising all the apr arrays and properties
!+
!-----------------------------------------------------------------------
subroutine init_apr(apr_level,ierr)
 use part,          only:npart,massoftype,aprmassoftype
 use apr_region,    only:set_apr_centre,set_apr_regions
 use utils_apr,     only:ntrack_max
 use get_apr_level, only:set_get_apr
 use io_summary,    only:print_apr,iosum_apr
 use io,            only:warning,fatal
 use dim,           only:maxvxyzu
 integer,         intent(inout) :: ierr
 integer(kind=1), intent(inout) :: apr_level(:)
 logical :: previously_set
 integer :: i

 ! the resolution levels are in addition to the base resolution
 apr_max = apr_max_in + 1
 if (split_dir == 2) do_relax = .true.

 ! if we're reading in a file that already has the levels set,
 ! don't override these
 previously_set = .false.
 if (sum(int(apr_level(1:npart))) > npart) then
    previously_set = .true.
    if (split_dir /= 2) do_relax = .false.
 endif

 if (.not.previously_set) then
    ! initialise the base resolution level
    if (ref_dir == 1) then
       apr_level(1:npart) = int(1,kind=1)
    else
       apr_level(1:npart) = int(apr_max,kind=1)
    endif

    ! also set the massoftype array
    ! if we are derefining we make sure that
    ! massoftype(igas) is associated with the
    ! largest particle (don't do it twice accidentally!)
    if (ref_dir == -1) then
       massoftype(:) = massoftype(:) * 2.**(apr_max -1)
       top_level = 1
    else
       top_level = apr_max
    endif
 endif

 ! now set the aprmassoftype array, this stores all the masses for the different resolution levels
 do i = 1,apr_max
    aprmassoftype(:,i) = massoftype(:)/(2.**(i-1))
 enddo

 ! how many regions do we need
 if (apr_type == 3) then
    ntrack_max = 999
    ntrack = 0 ! to start with
 elseif (apr_type == -1) then
    ntrack_max = 2
 else
    ntrack_max = 1
 endif

 if ((ntrack_max > 1) .and. (split_dir /= 3)) then
    split_dir = 3 ! no directional splitting for creating/multiple regions
    call warning('init_apr','resetting split_dir=3 because using multiple regions')
 endif

 allocate(apr_centre(3,ntrack_max),track_part(ntrack_max))
 apr_centre(:,:) = 0.

 ! initialise the shape of the region
 call set_get_apr()

 ! initiliase the regions
 call set_apr_centre(apr_type,apr_centre,ntrack,track_part)
 if (.not.allocated(apr_regions)) allocate(apr_regions(apr_max),npart_regions(apr_max))
 call set_apr_regions(ref_dir,apr_max,apr_regions,apr_rad,apr_drad)
 npart_regions = 0
 icentre = 1 ! to initialise

 ! certain splitdir need certain things
 if (maxvxyzu < 4 .and. split_dir == 2) then
   call fatal('init_apr','split_dir == 2 not compatible with choice of eos')
 endif

 ierr = 0

 ! print summary please
 print_apr = .true.
 iosum_apr(1) = ntrack
 iosum_apr(2) = apr_max

 if (apr_verbose) print*,'initialised apr'

end subroutine init_apr

!-----------------------------------------------------------------------
!+
!  Subroutine to check if particles need to be split or merged
!+
!-----------------------------------------------------------------------
subroutine update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
!$ use omp_lib
 use dim,        only:maxp,ind_timesteps,maxvxyzu
 use part,       only:ntot,isdead_or_accreted,igas,aprmassoftype,&
                    shuffle_part,iphase,iactive,maxp,npartoftype
 use quitdump,   only:quit
 use relaxem,    only:relax_particles
 use utils_apr,  only:find_closest_region,icentre
 use apr_region, only:set_apr_centre
 use io,         only:fatal
 use get_apr_level, only:get_apr,create_or_update_apr_clump
 use io_summary, only:iosum_apr,print_apr
 real,    intent(inout)         :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:)
 integer, intent(inout)         :: npart
 integer(kind=1), intent(inout) :: apr_level(:)
 integer :: ii,jj,kk,npartnew,nsplit_total,apri,npartold,ll,idx_len,j,apr_last
 integer :: n_ref,nrelax,nmerge,nkilled,nmerge_total,mm,n_to_split,should_split(maxp)
 real, allocatable :: xyzh_ref(:,:),force_ref(:,:),pmass_ref(:)
 real, allocatable :: xyzh_merge(:,:),vxyzu_merge(:,:), rneighs(:)
 integer, allocatable :: relaxlist(:),mergelist(:),iclosest
 integer, allocatable :: idx_merge(:),should_merge(:),scan_array(:),idx_split(:)
 real :: get_apr_in(3),xi,yi,zi,dx,dy,dz,rmin_local
 logical :: relax_in_loop

 ! if this routine doesn't need to be used, just skip it
 if (apr_max == 1) return

 if (npart >= 0.9*maxp) then
    call fatal('apr','maxp is not large enough; set --maxp on the command line to something larger than ',var='maxp',ival=maxp)
 endif

 ! if the centre of the region can move, update it
 call set_apr_centre(apr_type,apr_centre,ntrack,track_part)

 ! if we don't have any regions, skip routine
 if (ntrack == 0) return

 ! Just a metric
 if (apr_verbose) print*,'original npart is',npart

 ! Before adjusting the particles, if we're going to
 ! relax them then let's save the reference particles
 if (do_relax) then
    allocate(xyzh_ref(4,maxp),force_ref(3,maxp),pmass_ref(maxp),relaxlist(maxp))
    relaxlist = -1

    n_ref = 0
    xyzh_ref = 0.
    force_ref = 0.
    pmass_ref = 0.

    do ii = 1,npart
       if (.not.isdead_or_accreted(xyzh(4,ii))) then ! ignore dead particles
          n_ref = n_ref + 1
          xyzh_ref(1:4,n_ref) = xyzh(1:4,ii)
          pmass_ref(n_ref) = aprmassoftype(igas,apr_level(ii))
          force_ref(1:3,n_ref) = fxyzu(1:3,ii)*pmass_ref(n_ref)
       endif
    enddo
 else
    allocate(relaxlist(1))  ! it is passed but not used in merge
 endif

 ! Do any particles need to be split?
 npartnew = npart
 npartold = npart
 nsplit_total = 0
 nrelax = 0
 apri = 0 ! to avoid compiler errors
 apr_last = 0
 ! generally a safe guess, gets checked later
 allocate(scan_array(npart*apr_max),rneighs(npart*apr_max),idx_split(npart*apr_max))

 if (apr_verbose) print*,'started splitting'

 do jj = 1,apr_max-1
    do ll = 1,ntrack ! for multiple regions
       icentre = ll
       npartold = npartnew ! to account for new particles as they are being made
       should_split(:) = 0 ! reset
       rneighs(:) = 0.
       idx_split(:) = 0
       n_to_split = 0

       !$omp parallel default(none) &
       !$omp shared(npartold,iphase,apr_level,xyzh,should_split,get_apr,icentre) &
       !$omp shared(idx_split) &
       !$omp private(ii,get_apr_in,apri) &
       !$omp reduction(+:nsplit_total,n_to_split) reduction(max:apr_last)
       !$omp do
       split_over_active: do ii = 1,npartold
          ! only do this on active particles
          if (ind_timesteps) then
             if (.not.iactive(iphase(ii))) cycle split_over_active
          endif

          get_apr_in(1:3) = xyzh(1:3,ii)
          ! this is the refinement level it *should* have based
          ! on it's current position
          call get_apr(get_apr_in,icentre,apri)
          ! if the level it should have is greater than the
          ! level it does have, increment it up one
          if (apri > apr_level(ii)) then
             should_split(ii) = 1 ! record that this should be split
             nsplit_total = nsplit_total + 1
             n_to_split = n_to_split + 1
             apr_last = apri
          endif
       enddo split_over_active
       !$omp end do
       !$omp end parallel

       ! reallocate if required; if this happens even once just use the biggest possible
       if (n_to_split > size(scan_array)) then
         deallocate(scan_array,rneighs,idx_split)
         allocate(scan_array(maxp),rneighs(maxp),idx_split(maxp))
       endif

       ! create the scan array - this loop should *not* be parallelised
       scan_array(:) = 0
       do ii = 2,npartold
         scan_array(ii) = scan_array(ii-1) + should_split(ii-1)
       enddo

       ! make the particle list
       idx_len = n_to_split
       npartnew = npartnew + idx_len ! total number of particles (for now)
       npartoftype(igas) = npartoftype(igas) + n_to_split ! add to npartoftype
       npart = npartnew ! for splitpart

       ! exit here if there's nothing more to do
       if (n_to_split == 0) cycle

       !$omp parallel default(none) &
       !$omp shared(npartold,should_split,idx_split,scan_array,idx_len) &
       !$omp shared(rneighs,xyzh,adjusted_split) &
       !$omp private(ii,mm,rmin_local,j,xi,yi,zi,dx,dy,dz)
       !$omp do
       do ii = 1,npartold
         if (should_split(ii) == 1) then
            idx_split(scan_array(ii) + 1) = ii
         endif
       enddo
       !$omp end do

       if (adjusted_split) then
       !$omp do schedule(dynamic)
          do ii = 1,idx_len
             mm = idx_split(ii) ! original particle that should be split          
             xi = xyzh(1,mm)
             yi = xyzh(2,mm)
             zi = xyzh(3,mm)

             rmin_local = huge(1.0)

             do j = 1,npartold
                if (j == mm) cycle
                dx = xi - xyzh(1,j)
                dy = yi - xyzh(2,j)
                dz = zi - xyzh(3,j)
                rmin_local = min(rmin_local,dx*dx + dy*dy + dz*dz)
             enddo
             rneighs(ii) = sqrt(rmin_local)
          enddo
       !$omp end do
       endif
       !$omp end parallel

       ! if relaxing, make some adjustments here:
       ! just use the first particle that has been marked to split
       ! to establish if we should be relaxing at all
       relax_in_loop = (do_relax .and. (gr .or. apr_last == top_level))

       ! now go through and actually split them - this should *probably* not be parallelised
       ! due to the content of the nested functions, idx_len probably isn't that long either
       if (adjusted_split) then
         do ii = 1,idx_len
            mm = idx_split(ii) ! original particle that should be split
            kk = npartold + ii ! location in array for new particle
            call splitpart(mm,kk,rneigh=rneighs(ii))
            if (relax_in_loop) then
                relaxlist(nrelax + ii) = mm
                relaxlist(nrelax + n_to_split + ii) = kk
            endif   
         enddo     
       else
         do ii = 1,idx_len
            mm = idx_split(ii) ! original particle that should be split
            kk = npartold + ii ! location in array for new particle
            call splitpart(mm,kk)
            if (relax_in_loop) then
                relaxlist(nrelax + ii) = mm
                relaxlist(nrelax + n_to_split + ii) = kk
            endif   
         enddo
       endif

       ! if relaxing, update the total number that will be relaxed
       if (relax_in_loop) nrelax = nrelax + 2*n_to_split
    enddo
 enddo

 ! Take into account all the added particles
 npart = npartnew
 ntot = npartnew
 if (apr_verbose) then
    print*,'split: ',nsplit_total
    print*,'npart: ',npart
 endif

 ! Do any particles need to be merged?
 deallocate(scan_array)
 allocate(mergelist(npart),xyzh_merge(4,npart),vxyzu_merge(maxvxyzu,npart))
 allocate(idx_merge(npart),should_merge(npart),scan_array(npart))
 npart_regions = 0
 nmerge_total = 0
 iclosest = 1
 do jj = 1,apr_max-1
    do ll = 1, ntrack
       icentre = ll
       kk = apr_max - jj + 1             ! to go from apr_max -> 2
       mergelist = -1 ! initialise
       nmerge = 0
       nkilled = 0
       xyzh_merge = 0.
       vxyzu_merge = 0.

       should_merge(:) = 0
       scan_array(:) = 0
       idx_merge(:) = 0

       ! identify what should be merged
       !$omp parallel do default(none) &
       !$omp shared(npart,apr_level,kk,xyzh,vxyzu,ntrack,ll,should_merge,iphase,apr_centre) &
       !$omp private(ii,iclosest) &
       !$omp reduction(+:nmerge)
       merge_over_active: do ii = 1,npart
          if ((apr_level(ii) == kk) .and. (.not.isdead_or_accreted(xyzh(4,ii)))) then ! avoid already dead particles
             if (ind_timesteps) then
                if (.not.iactive(iphase(ii))) cycle merge_over_active
             endif
             if (ntrack > 1) call find_closest_region(xyzh(1:3,ii),ntrack,apr_centre,iclosest)

             if ((ntrack == 1) .or. (iclosest == ll)) then
                should_merge(ii) = 1
                nmerge = nmerge + 1
             endif
          endif
       enddo merge_over_active
       !$omp end parallel do

       ! create the scan array - this loop should *not* be parallelised
       scan_array(:) = 0
       do ii = 2,npart
         scan_array(ii) = scan_array(ii-1) + should_merge(ii-1)
       enddo

       !$omp parallel do default(none) &
       !$omp shared(should_merge,idx_merge,scan_array,xyzh_merge,vxyzu_merge,kk) &
       !$omp shared(xyzh,vxyzu,npart) &
       !$omp private(ii,mm) &
       !$omp reduction(+:npart_regions)
       do ii = 1,npart
         if (should_merge(ii) == 1) then
            mm = scan_array(ii) + 1
            idx_merge(mm) = ii
            xyzh_merge(1:4,mm) = xyzh(1:4,ii)
            vxyzu_merge(1:3,mm) = vxyzu(1:3,ii)
            npart_regions(kk) = npart_regions(kk) + 1
         endif
       enddo
       !$omp end parallel do

       if (apr_verbose) print*,nmerge,'particles selected for merge'
       ! Now send them to be merged
       if (nmerge > 1) call merge_with_special_tree(nmerge,idx_merge,xyzh_merge(:,1:nmerge),&
                                            vxyzu_merge(:,1:nmerge),kk,xyzh,vxyzu,apr_level,nkilled,&
                                            nrelax,relaxlist,npartnew)
       nmerge_total = nmerge_total + nkilled ! actually merged
       if (apr_verbose) then
          print*,'merged: ',nkilled,kk
          print*,'npart: ',npartnew - nkilled
       endif
       npart_regions(kk) = npart_regions(kk) - nkilled
    enddo
 enddo
 ! update npart as required
 npart = npartnew
 npart_regions(1) = npartnew - sum(npart_regions(2:apr_max))
 if (apr_verbose) print*,'particles at each level:',npart_regions(:)

 ! If we need to relax, do it here
 if (nrelax > 0 .and. do_relax) call relax_particles(npart,n_ref,xyzh_ref,force_ref,nrelax,relaxlist)
 ! Turn it off now because we only want to do this on first splits
 if (.not. gr) do_relax = .false.

 ! As we may have killed particles, time to do an array shuffle
 call shuffle_part(npart)

 ! Tidy up
 if (do_relax) then
    deallocate(xyzh_ref,force_ref,pmass_ref)
 endif
 deallocate(relaxlist,should_merge,idx_merge,scan_array,rneighs,idx_split)

 if (apr_verbose) print*,'total particles at end of apr: ',npart

 ! summary variables
 print_apr = .true.
 iosum_apr(1) = ntrack
 iosum_apr(2) = apr_max
 iosum_apr(3) = iosum_apr(3) + nsplit_total
 iosum_apr(4) = iosum_apr(4) + nmerge_total
 do ii = 1,apr_max
    iosum_apr(ii+4) = count(apr_level(1:npart) == ii)
 enddo

end subroutine update_apr

!-----------------------------------------------------------------------
!+
!  routine to split one particle into two
!+
!-----------------------------------------------------------------------
subroutine splitpart(i,i_new,rneigh)
 use part,         only:xyzh
 use physcon,      only:pi
 use vectorutils, only:cross_product3D,rotatevec
 use get_apr_level, only:split_dir_func
 use dim, only:ind_timesteps
 integer, intent(in) :: i,i_new
 real, optional :: rneigh
 real :: sep

 if (adjusted_split) then
    sep = min(sep_factor*xyzh(4,i),0.35*rneigh)
    sep = sep/xyzh(4,i)  ! for consistency later on
 else
    sep = sep_factor
 endif

 call split_dir_func(i,i_new,sep)

end subroutine splitpart

!-----------------------------------------------------------------------
!+
!  Take in all particles that *might* be merged at this apr_level
!  and use our special tree to merge what has left the region
!+
!-----------------------------------------------------------------------
subroutine merge_with_special_tree(nmerge,mergelist,xyzh_merge,vxyzu_merge,current_apr,&
                                     xyzh,vxyzu,apr_level,nkilled,nrelax,relaxlist,npartnew)
 use neighkdtree,   only:build_tree,ncells,leaf_is_active,get_cell_location
 use mpiforce,      only:cellforce
 use kdtree,        only:inodeparts,inoderange
 use part,          only:kill_particle,npartoftype,igas
 use part,          only:combine_two_particles
 use dim,           only:ind_timesteps,maxvxyzu
 use get_apr_level, only:get_apr,put_in_smallest_bin
 integer,         intent(inout) :: nmerge,nkilled,nrelax,relaxlist(:),npartnew
 integer(kind=1), intent(inout) :: apr_level(:)
 integer,         intent(in)    :: current_apr,mergelist(:)
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,            intent(inout) :: xyzh_merge(:,:),vxyzu_merge(:,:)
 integer :: remainder,icell,n_cell,apri,m
 integer :: eldest,tuther
 real    :: com(3)
 type(cellforce)        :: cell

 ! First ensure that we're only sending in a multiple of 2 to the tree
 remainder = modulo(nmerge,2)
 nmerge = nmerge - remainder

 call build_tree(nmerge,nmerge,xyzh_merge(:,1:nmerge),vxyzu_merge(:,1:nmerge),&
                      for_apr=.true.)
 ! Now use the centre of mass of each cell to check whether it should
 ! be merged or not
 com = 0.
 over_cells: do icell=1,int(ncells)
    if (leaf_is_active(icell) == 0) cycle over_cells !--skip empty cells
    n_cell = inoderange(2,icell)-inoderange(1,icell)+1

    call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
    com(1) = cell%xpos(1)
    com(2) = cell%xpos(2)
    com(3) = cell%xpos(3)

    call get_apr(com(1:3),icentre,apri)

    ! If the apr level based on the com is lower than the current level,
    ! we merge!
    if (apri < current_apr) then

       eldest = mergelist(inodeparts(inoderange(1,icell)))
       tuther = mergelist(inodeparts(inoderange(1,icell) + 1)) !as in kdtree

       ! merge by averaging everything
       call combine_two_particles(eldest,tuther)

       xyzh(4,eldest) = (0.5*(xyzh(4,eldest) + xyzh(4,tuther)))*(2.0**(1./3.))
       apr_level(eldest) = apr_level(eldest) - int(1,kind=1)
       if (ind_timesteps) call put_in_smallest_bin(eldest)

       ! add it to the shuffling list if needed
       if (do_relax) then
          nrelax = nrelax + 1
          relaxlist(nrelax) = eldest
       endif

       ! discard tuther (t'other)
       call kill_particle(tuther,npartoftype)
       nkilled = nkilled + 2 ! this refers to the number of children killed
       ! If this particle was on the shuffle list previously, take it off
       do m = 1,nrelax
          if (relaxlist(m) == tuther) relaxlist(m) = 0
       enddo
    endif

 enddo over_cells

end subroutine merge_with_special_tree

end module apr