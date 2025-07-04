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
! :Dependencies: apr_region, dim, get_apr_level, io, io_summary, kdtree,
!   linklist, mpiforce, part, physcon, quitdump, random, relaxem,
!   timestep_ind, utils_apr, vectorutils
!
 use dim, only:use_apr
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
 logical :: directional = .true.

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
 integer,         intent(inout) :: ierr
 integer(kind=1), intent(inout) :: apr_level(:)
 logical :: previously_set
 integer :: i

 ! the resolution levels are in addition to the base resolution
 apr_max = apr_max_in + 1

 ! if we're reading in a file that already has the levels set,
 ! don't override these
 previously_set = .false.
 if (sum(int(apr_level(1:npart))) > npart) then
    previously_set = .true.
    do_relax = .false.
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
    ntrack_max = 1000
    ntrack = 0 ! to start with
 elseif (apr_type == -1) then
    ntrack_max = 2
 else
    ntrack_max = 1
 endif

 if (ntrack_max > 1) directional = .false. ! no directional splitting for creating/multiple regions

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
 use dim,        only:maxp,ind_timesteps,maxvxyzu
 use part,       only:ntot,isdead_or_accreted,igas,aprmassoftype,&
                    shuffle_part,iphase,iactive,maxp
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
 integer :: ii,jj,kk,npartnew,nsplit_total,apri,npartold,ll
 integer :: n_ref,nrelax,nmerge,nkilled,apr_current,nmerge_total
 real, allocatable :: xyzh_ref(:,:),force_ref(:,:),pmass_ref(:)
 real, allocatable :: xyzh_merge(:,:),vxyzu_merge(:,:)
 integer, allocatable :: relaxlist(:),mergelist(:),iclosest
 real :: get_apr_in(3)

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

 if (apr_verbose) print*,'started splitting'

 do jj = 1,apr_max-1
    do ll = 1,ntrack ! for multiple regions
       icentre = ll
       npartold = npartnew ! to account for new particles as they are being made
       split_over_active: do ii = 1,npartold

          ! only do this on active particles
          if (ind_timesteps) then
             if (.not.iactive(iphase(ii))) cycle split_over_active
          endif

          apr_current = apr_level(ii)
          get_apr_in(1) = xyzh(1,ii)
          get_apr_in(2) = xyzh(2,ii)
          get_apr_in(3) = xyzh(3,ii)
          ! this is the refinement level it *should* have based
          ! on it's current position
          call get_apr(get_apr_in,icentre,apri)
          ! if the level it should have is greater than the
          ! level it does have, increment it up one
          if (apri > apr_current) then
             call splitpart(ii,npartnew)
             if (do_relax .and. (apri == top_level)) then
                nrelax = nrelax + 2
                relaxlist(nrelax-1) = ii
                relaxlist(nrelax)   = npartnew
             endif
             nsplit_total = nsplit_total + 1
          endif
       enddo split_over_active
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
 allocate(mergelist(npart),xyzh_merge(4,npart),vxyzu_merge(maxvxyzu,npart))
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

       merge_over_active: do ii = 1,npart
          ! note that here we only do this process for particles that are not already counted in the blending region
          if ((apr_level(ii) == kk) .and. (.not.isdead_or_accreted(xyzh(4,ii)))) then ! avoid already dead particles
             if (ind_timesteps) then
                if (.not.iactive(iphase(ii))) cycle merge_over_active
             endif
             if (ntrack > 1) call find_closest_region(xyzh(1:3,ii),iclosest)
             if (iclosest == ll) then
                nmerge = nmerge + 1
                mergelist(nmerge) = ii
                xyzh_merge(1:4,nmerge) = xyzh(1:4,ii)
                vxyzu_merge(1:3,nmerge) = vxyzu(1:3,ii)
                npart_regions(kk) = npart_regions(kk) + 1
             endif
          endif
       enddo merge_over_active
       if (apr_verbose) print*,nmerge,'particles selected for merge'
       ! Now send them to be merged
       if (nmerge > 1) call merge_with_special_tree(nmerge,mergelist(1:nmerge),xyzh_merge(:,1:nmerge),&
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
 do_relax = .false.

 ! As we may have killed particles, time to do an array shuffle
 call shuffle_part(npart)

 ! Tidy up
 if (do_relax) then
    deallocate(xyzh_ref,force_ref,pmass_ref)
 endif
 deallocate(mergelist,relaxlist)

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
subroutine splitpart(i,npartnew)
 use part,         only:copy_particle_all,apr_level,xyzh,vxyzu,npartoftype,igas
 use part,         only:set_particle_type
 use physcon,      only:pi
 use dim,          only:ind_timesteps
 use random,       only:ran2
 use vectorutils, only:cross_product3D,rotatevec
 use utils_apr,  only:apr_region_is_circle,icentre
 integer, intent(in) :: i
 integer, intent(inout) :: npartnew
 integer :: j,npartold,next_door
 real :: theta,dx,dy,dz,x_add,y_add,z_add,sep,rneigh
 real :: v(3),u(3),w(3),a,b,c,mag_v
 integer, save :: iseed = 4
 integer(kind=1) :: aprnew

 if (adjusted_split) then
    call closest_neigh(i,next_door,rneigh)
    sep = min(sep_factor*xyzh(4,i),0.35*rneigh)
    sep = sep/xyzh(4,i)  ! for consistency later on
 else
    sep = sep_factor
 endif

 ! Calculate the plane that the particle must be split along
 ! to be tangential to the splitting region. Particles are split
 ! on this plane but rotated randomly on it.
 dx = xyzh(1,i) - apr_centre(1,icentre)
 dy = xyzh(2,i) - apr_centre(2,icentre)
 if (.not.apr_region_is_circle) then
    dz = xyzh(3,i) - apr_centre(3,icentre)       ! for now, let's split about the CoM

    if (directional) then
       ! Calculate a vector, v, that lies on the plane
       u = (/1.0,0.5,1.0/)
       w = (/dx,dy,dz/)
       call cross_product3D(u,w,v)

       ! rotate it around the normal to the plane by a random amount
       theta = ran2(iseed)*2.*pi
       call rotatevec(v,w,theta)
    else
       ! No directional splitting, so just create a unit vector in a random direction
       a = ran2(iseed) - 0.5
       b = ran2(iseed) - 0.5
       c = ran2(iseed) - 0.5
       v = (/a, b, c/)
    endif

    mag_v = sqrt(dot_product(v,v))
    if (mag_v > tiny(mag_v)) then
       v = v/mag_v
    else
       v = 0.
    endif
 else
    dz = 0.
    u = 0.
    w = 0.
    v = 0.
    theta = atan2(dy,dx) + 0.5*pi
    v(1) = cos(theta)
    v(2) = sin(theta)
 endif

 ! Now apply it
 x_add = sep*v(1)*xyzh(4,i)
 y_add = sep*v(2)*xyzh(4,i)
 z_add = sep*v(3)*xyzh(4,i)

 npartold = npartnew
 npartnew = npartold + 1
 npartoftype(igas) = npartoftype(igas) + 1
 aprnew = apr_level(i) + int(1,kind=1) ! to prevent compiler warnings


 !--create the new particle
 do j=npartold+1,npartnew
    call copy_particle_all(i,j,new_part=.true.)
    xyzh(1,j) = xyzh(1,i) + x_add
    xyzh(2,j) = xyzh(2,i) + y_add
    xyzh(3,j) = xyzh(3,i) + z_add
    vxyzu(:,j) = vxyzu(:,i)
    xyzh(4,j) = xyzh(4,i)*(0.5**(1./3.))
    apr_level(j) = aprnew
    if (ind_timesteps) call put_in_smallest_bin(j)
 enddo


 ! Edit the old particle that was sent in and kept
 xyzh(1,i) = xyzh(1,i) - x_add
 xyzh(2,i) = xyzh(2,i) - y_add
 xyzh(3,i) = xyzh(3,i) - z_add
 apr_level(i) = aprnew
 xyzh(4,i) = xyzh(4,i)*(0.5**(1./3.))
 if (ind_timesteps) call put_in_smallest_bin(i)


end subroutine splitpart

!-----------------------------------------------------------------------
!+
!  Take in all particles that *might* be merged at this apr_level
!  and use our special tree to merge what has left the region
!+
!-----------------------------------------------------------------------
subroutine merge_with_special_tree(nmerge,mergelist,xyzh_merge,vxyzu_merge,current_apr,&
                                     xyzh,vxyzu,apr_level,nkilled,nrelax,relaxlist,npartnew)
 use linklist, only:set_linklist,ncells,ifirstincell,get_cell_location
 use mpiforce, only:cellforce
 use kdtree,   only:inodeparts,inoderange
 use part,     only:kill_particle,npartoftype,igas
 use part,     only:combine_two_particles
 use dim,      only:ind_timesteps,maxvxyzu
 use get_apr_level, only:get_apr
 integer,         intent(inout) :: nmerge,nkilled,nrelax,relaxlist(:),npartnew
 integer(kind=1), intent(inout) :: apr_level(:)
 integer,         intent(in)    :: current_apr,mergelist(:)
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,            intent(inout) :: xyzh_merge(:,:),vxyzu_merge(:,:)
 integer :: remainder,icell,i,n_cell,apri,m
 integer :: eldest,tuther
 real    :: com(3)
 type(cellforce)        :: cell

 ! First ensure that we're only sending in a multiple of 2 to the tree
 remainder = modulo(nmerge,2)
 nmerge = nmerge - remainder

 call set_linklist(nmerge,nmerge,xyzh_merge(:,1:nmerge),vxyzu_merge(:,1:nmerge),&
                      for_apr=.true.)
 ! Now use the centre of mass of each cell to check whether it should
 ! be merged or not
 com = 0.
 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)
    if (i == 0) cycle over_cells !--skip empty cells
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

!-----------------------------------------------------------------------
!+
!  Find the closest neighbour to a particle (needs replacing)
!+
!-----------------------------------------------------------------------
subroutine closest_neigh(i,next_door,rmin)
 use part, only:xyzh,npart
 integer, intent(in)  :: i
 integer, intent(out) :: next_door
 real,    intent(out) :: rmin
 real :: dx,dy,dz,rtest
 integer :: j

 ! DP note: this is not MPI safe...
 rmin = huge(rmin)
 next_door = 0
 do j = 1,npart
    if (j == i) cycle
    dx = xyzh(1,i) - xyzh(1,j)
    dy = xyzh(2,i) - xyzh(2,j)
    dz = xyzh(3,i) - xyzh(3,j)
    rtest = dx**2 + dy**2 + dz**2
    if (rtest < rmin) then
       next_door = j
       rmin = rtest
    endif
 enddo

 rmin = sqrt(rmin)

end subroutine closest_neigh

!-----------------------------------------------------------------------
!+
!  routine to put a particle on the shortest timestep
!+
!-----------------------------------------------------------------------
subroutine put_in_smallest_bin(i)
 use timestep_ind, only:nbinmax
 use part,         only:ibin
 integer, intent(in) :: i

 ibin(i) = nbinmax

end subroutine put_in_smallest_bin

end module apr
