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
! :Runtime parameters:
!   - apr_drad   : *size of step to next region*
!   - apr_max    : *number of additional refinement levels (3 -> 2x resolution)*
!   - apr_rad    : *radius of innermost region*
!   - apr_type   : *1: static, 2: sink, 3: clumps, 4: sequential sinks, 5: com*
!   - ref_dir    : *increase (1) or decrease (-1) resolution*
!   - track_part : *number of sink to track*
!
! :Dependencies: apr_region, dim, infile_utils, io, kdtree, linklist,
!   mpiforce, part, physcon, ptmass, quitdump, random, relaxem,
!   timestep_ind, vectorutils
!
 use dim, only:use_apr
 implicit none

 public :: init_apr,update_apr,read_options_apr,write_options_apr
 public :: create_or_update_apr_clump
 public :: use_apr

 ! default values for runtime parameters
 integer, public :: apr_max_in = 3
 integer, public :: ref_dir = 1
 integer, public :: apr_type = 1
 integer, public :: apr_max = 4
 real,    public :: apr_rad = 1.0
 real,    public :: apr_drad = 0.1
 real,    public :: apr_centre(3) = 0.

 private
 integer :: top_level = 1, ntrack = 0, track_part = 0
 real, allocatable    :: apr_regions(:)
 integer, allocatable :: npart_regions(:)
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
 use part,       only:npart,massoftype,aprmassoftype
 use apr_region, only:set_apr_centre,set_apr_regions
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

 ! initiliase the regions
 call set_apr_centre(apr_type,apr_centre,ntrack,track_part)
 if (.not.allocated(apr_regions)) allocate(apr_regions(apr_max),npart_regions(apr_max))
 call set_apr_regions(ref_dir,apr_max,apr_regions,apr_rad,apr_drad)
 npart_regions = 0

 ! now set the aprmassoftype array, this stores all the masses for the different resolution levels
 do i = 1,apr_max
    aprmassoftype(:,i) = massoftype(:)/(2.**(i-1))
 enddo

 ierr = 0

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
                    shuffle_part,iphase,iactive,poten,&
                    maxp,xyzmh_ptmass
 use quitdump,   only:quit
 use relaxem,    only:relax_particles
 use apr_region, only:dynamic_apr,set_apr_centre
 use io,         only:fatal
 real,    intent(inout)         :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:)
 integer, intent(inout)         :: npart
 integer(kind=1), intent(inout) :: apr_level(:)
 integer :: ii,jj,kk,npartnew,nsplit_total,apri,npartold
 integer :: n_ref,nrelax,nmerge,nkilled,apr_current
 real, allocatable :: xyzh_ref(:,:),force_ref(:,:),pmass_ref(:)
 real, allocatable :: xyzh_merge(:,:),vxyzu_merge(:,:)
 integer, allocatable :: relaxlist(:),mergelist(:)
 real :: get_apr_in(3),radi,radi_max

 if (npart >= 0.9*maxp) then
    call fatal('apr','maxp is not large enough; set --maxp on the command line to something larger than ',var='maxp',ival=maxp)
 endif

 ! if the centre of the region can move, update it
 if (dynamic_apr) then
    if (ntrack > 0) then
       call create_or_update_apr_clump(npart,xyzh,vxyzu,poten,apr_level,&
                                        xyzmh_ptmass,aprmassoftype)
    else
       call set_apr_centre(apr_type,apr_centre,ntrack,track_part)
    endif
 endif

 ! If this routine doesn't need to be used, just skip it
 if (apr_max == 1) return

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


 do jj = 1,apr_max-1
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
       call get_apr(get_apr_in,apri)
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
 do jj = 1,apr_max-1
    kk = apr_max - jj + 1             ! to go from apr_max -> 2
    mergelist = -1 ! initialise
    nmerge = 0
    nkilled = 0
    xyzh_merge = 0.
    vxyzu_merge = 0.
    radi_max = 0.

    merge_over_active: do ii = 1,npart
       ! note that here we only do this process for particles that are not already counted in the blending region
       if ((apr_level(ii) == kk) .and. (.not.isdead_or_accreted(xyzh(4,ii)))) then ! avoid already dead particles
          if (ind_timesteps) then
             if (.not.iactive(iphase(ii))) cycle merge_over_active
          endif
          nmerge = nmerge + 1
          mergelist(nmerge) = ii
          xyzh_merge(1:4,nmerge) = xyzh(1:4,ii)
          vxyzu_merge(:,nmerge) = vxyzu(:,ii)
          npart_regions(kk) = npart_regions(kk) + 1
       endif
       radi = sqrt(dot_product(xyzh(1:3,ii),xyzh(1:3,ii)))
       if (radi > radi_max) radi_max = radi
    enddo merge_over_active
    ! Now send them to be merged
    if (nmerge > 1) call merge_with_special_tree(nmerge,mergelist(1:nmerge),xyzh_merge(:,1:nmerge),&
                                         vxyzu_merge(:,1:nmerge),kk,xyzh,vxyzu,apr_level,nkilled,&
                                         nrelax,relaxlist,npartnew)
    if (apr_verbose) then
       print*,'merged: ',nkilled,kk
       print*,'npart: ',npartnew - nkilled
    endif
    npart_regions(kk) = npart_regions(kk) - nkilled
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

end subroutine update_apr

!-----------------------------------------------------------------------
!+
!  routine to return the adaptive particle refinement level based on position
!  and the boundaries set by the apr_* arrays
!+
!-----------------------------------------------------------------------
subroutine get_apr(pos,apri)
 use io, only:fatal
 use apr_region, only:apr_region_is_circle
 real, intent(in)     :: pos(3)
 integer, intent(out) :: apri
 integer :: jj, kk
 real :: dx,dy,dz,r

 apri = -1 ! to prevent compiler warnings

 do jj = 1,apr_max
    if (ref_dir == 1) then
       kk = apr_max - jj + 1       ! going from apr_max -> 1
    else
       kk = jj                    ! going from 1 -> apr_max
    endif
    dx = pos(1) - apr_centre(1)
    dy = pos(2) - apr_centre(2)
    dz = pos(3) - apr_centre(3)

    if (apr_region_is_circle) dz = 0.

    r = sqrt(dx**2 + dy**2 + dz**2)

    if (r < apr_regions(kk)) then
       apri = kk
       return
    endif
 enddo

 if (apri == -1) call fatal('apr_region, get_apr','could not find apr level')

end subroutine get_apr

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
 use apr_region,  only:apr_region_is_circle
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
 dx = xyzh(1,i) - apr_centre(1)
 dy = xyzh(2,i) - apr_centre(2)
 if (.not.apr_region_is_circle) then
    dz = xyzh(3,i) - apr_centre(3)

    ! Calculate a vector, v, that lies on the plane
    u = (/1.0,0.5,1.0/)
    w = (/dx,dy,dz/)
    call cross_product3D(u,w,v)

    ! rotate it around the normal to the plane by a random amount
    theta = ran2(iseed)*2.*pi
    call rotatevec(v,w,theta)

    if (.not.directional) then
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
 use part,     only:kill_particle,npartoftype,aprmassoftype,igas
 use dim,      only:ind_timesteps,maxvxyzu
 integer,         intent(inout) :: nmerge,nkilled,nrelax,relaxlist(:),npartnew
 integer(kind=1), intent(inout) :: apr_level(:)
 integer,         intent(in)    :: current_apr,mergelist(:)
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,            intent(inout) :: xyzh_merge(:,:),vxyzu_merge(:,:)
 integer :: remainder,icell,i,n_cell,apri,m
 integer :: eldest,tuther
 real    :: com(3),pmassi,v2eldest,v2tuther,v_mag,vcom(maxvxyzu),vcom_mag
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
    call get_apr(com(1:3),apri)

    ! If the apr level based on the com is lower than the current level,
    ! we merge!
    if (apri < current_apr) then
       eldest = mergelist(inodeparts(inoderange(1,icell)))
       tuther = mergelist(inodeparts(inoderange(1,icell) + 1)) !as in kdtree

       ! keep eldest, reassign it to have the com properties
       xyzh(1,eldest) = cell%xpos(1)
       xyzh(2,eldest) = cell%xpos(2)
       xyzh(3,eldest) = cell%xpos(3)

       ! calculate the magnitude of the velocity to conserve kinetic energy (for now!)
       ! direction is in com, magnitude set by conservation
       pmassi = aprmassoftype(igas,apr_level(eldest))
       v2eldest = dot_product(vxyzu(1:3,eldest),vxyzu(1:3,eldest))
       v2tuther = dot_product(vxyzu(1:3,tuther),vxyzu(1:3,tuther))
       v_mag = sqrt(0.5*(v2eldest + v2tuther))
       vcom = 0.5*(vxyzu(:,eldest) + vxyzu(:,tuther))
       vcom_mag = sqrt(dot_product(vcom(1:3),vcom(1:3)))
       !vxyzu(1:3,eldest) = vcom/vcom_mag * v_mag

       ! or instead, just set it from the com
       vxyzu(:,eldest) = vcom

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
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_apr(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_apr'
 logical :: igotall1,igotall2

 imatch  = .true.
 igotall1 = .true.
 igotall2 = .true.
 select case(trim(name))
 case('apr_max')
    read(valstring,*,iostat=ierr) apr_max_in
    ngot = ngot + 1
    if (apr_max_in  <  0) call fatal(label,'apr_max < 0 in input options')
 case('ref_dir')
    read(valstring,*,iostat=ierr) ref_dir
    ngot = ngot + 1
 case('apr_type')
    read(valstring,*,iostat=ierr) apr_type
    ngot = ngot + 1
 case('apr_rad')
    read(valstring,*,iostat=ierr) apr_rad
    ngot = ngot + 1
    if (apr_rad  <  tiny(apr_rad)) call fatal(label,'apr_rad too small in input options')
 case('apr_drad')
    read(valstring,*,iostat=ierr) apr_drad
    ngot = ngot + 1
    if (apr_drad  <  tiny(apr_drad)) call fatal(label,'apr_drad too small in input options')
 case default
    imatch = .false.
    select case(apr_type)
    case(1)
       call read_options_apr1(name,valstring,imatch,igotall1,ierr)
    case(2,4)
       call read_options_apr2(name,valstring,imatch,igotall2,ierr)
    end select
 end select
 igotall = (ngot >= 5) .and. igotall1 .and. igotall2

end subroutine read_options_apr

!-----------------------------------------------------------------------
!+
!  extra subroutines for reading in different styles of apr zones
!+
!-----------------------------------------------------------------------
subroutine read_options_apr1(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_apr1'

 imatch  = .true.
 select case(trim(name))
 case('apr_centre(1)')
    read(valstring,*,iostat=ierr) apr_centre(1)
    ngot = ngot + 1
 case('apr_centre(2)')
    read(valstring,*,iostat=ierr) apr_centre(2)
    ngot = ngot + 1
 case('apr_centre(3)')
    read(valstring,*,iostat=ierr) apr_centre(3)
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 3)

end subroutine read_options_apr1

subroutine read_options_apr2(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_apr2'

 imatch  = .true.
 select case(trim(name))
 case('track_part')
    read(valstring,*,iostat=ierr) track_part
    ngot = ngot + 1
    if (track_part  <  1) call fatal(label,'track_part not chosen in input options')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)

end subroutine read_options_apr2

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_apr(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for adaptive particle refinement'
 call write_inopt(apr_max_in,'apr_max','number of additional refinement levels (3 -> 2x resolution)',iunit)
 call write_inopt(ref_dir,'ref_dir','increase (1) or decrease (-1) resolution',iunit)
 call write_inopt(apr_type,'apr_type','1: static, 2: sink, 3: clumps, 4: sequential sinks, 5: com',iunit)

 select case (apr_type)
 case (1)
    call write_inopt(apr_centre(1),'apr_centre(1)','centre of region x position',iunit)
    call write_inopt(apr_centre(2),'apr_centre(2)','centre of region y position',iunit)
    call write_inopt(apr_centre(3),'apr_centre(3)','centre of region z position',iunit)
 case(2,4)
    call write_inopt(track_part,'track_part','number of sink to track',iunit)
 case default
    ! write nothing
 end select

 call write_inopt(apr_rad,'apr_rad','radius of innermost region',iunit)
 call write_inopt(apr_drad,'apr_drad','size of step to next region',iunit)

end subroutine write_options_apr

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

!-----------------------------------------------------------------------
!+
!  Create a new apr region that is centred on a dense clump
!  (This is work in progress)
!+
!-----------------------------------------------------------------------
subroutine create_or_update_apr_clump(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,aprmassoftype)
 use apr_region, only:set_apr_centre
 use part, only:igas,rhoh
 use ptmass, only:rho_crit_cgs
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: apr_level(:)
 real, intent(in) :: xyzh(:,:), vxyzu(:,:), aprmassoftype(:,:),xyzmh_ptmass(:,:)
 real(kind=4), intent(in) :: poten(:)
 integer :: nbins, ii, ibin, nmins, jj, apri
 integer, allocatable :: counter(:), minima(:), min_particle(:)
 real, allocatable :: radius(:), ave_poten(:)
 real :: rin, rout, dbin, dx, dy, dz, rad, gradleft, gradright
 real :: minpoten, pmassi, rhoi

 ! set up arrays
 nbins = 100
 allocate(counter(nbins),radius(nbins),ave_poten(nbins),&
    minima(nbins),min_particle(nbins))

 ! Currently hardwired but this is problematic
 rin = 10.
 rout = 100.
 dbin = (rout-rin)/real(nbins-1)
 do ii = 1,nbins
    radius(ii) = rin + real(ii-1)*dbin
 enddo

 ave_poten = 0.
 counter = 0
 ! Create an azimuthally averaged potential energy vs. radius profile
 do ii = 1,npart
    dx = xyzh(1,ii) - xyzmh_ptmass(1,1)
    dy = xyzh(2,ii) - xyzmh_ptmass(2,1)
    dz = xyzh(3,ii) - xyzmh_ptmass(3,1)
    rad = sqrt(dx**2 + dy**2 + dz**2)
    pmassi = aprmassoftype(igas,apr_level(ii))

    ibin = int((rad - radius(1))/dbin + 1)
    if ((ibin > nbins) .or. (ibin < 1)) cycle

    ave_poten(ibin) = ave_poten(ibin) + poten(ii)/pmassi
    counter(ibin) = counter(ibin) + 1
 enddo

 ! average with the number of particles in the bin
 do ii = 1,nbins
    if (counter(ii) > 0) then
       ave_poten(ii) = ave_poten(ii)/counter(ii)
    else
       ave_poten(ii) = 0.
    endif
 enddo

 ! Identify what radius the local minima are at
 minima = 0
 nmins = 0
 do ii = 2, nbins-1
    gradleft = (ave_poten(ii) - ave_poten(ii-1))/(radius(ii) - radius(ii-1))
    gradright = (ave_poten(ii+1) - ave_poten(ii))/(radius(ii+1) - radius(ii))
    if (gradleft * gradright < 0.) then
       nmins = nmins + 1
       minima(nmins) = ii
    endif
 enddo
 if (nmins == 0) return

 ! Identify the particles in these minima that have the lowest potential energy
 ! this is quite inefficient, in future should save these above into the bins so
 ! you just need to cycle through the subset? Don't know if this is faster
 minpoten = 1.0
 do jj = 1,nmins
    do ii = 1,npart
       dx = xyzh(1,ii) - xyzmh_ptmass(1,1)
       dy = xyzh(2,ii) - xyzmh_ptmass(2,1)
       dz = xyzh(3,ii) - xyzmh_ptmass(3,1)
       rad = sqrt(dx**2 + dy**2 + dz**2)
       pmassi = aprmassoftype(igas,apr_level(ii))

       ibin = int((rad - radius(1))/dbin + 1)
       if ((ibin == (minima(jj))) .or. &
        (ibin - 1 == (minima(jj))) .or. &
        (ibin + 1 == (minima(jj)))) then
          if ((poten(ii)/pmassi) < minpoten) then
             minpoten = poten(ii)/pmassi
             min_particle(jj) = ii
          endif
       endif
    enddo
 enddo

 ! For the moment, force there to only be one minimum
 ! and let it be the lowest
 nmins = 1

 ! Check they are not already within a region of low potential energy
 ! If they are, replace the existing particle as the one to be tracked
 over_mins: do jj = 1,nmins
    ii = min_particle(jj)
    ! check that the particle at the lowest potential energy has also met the
    ! density criteria
    pmassi = aprmassoftype(igas,apr_level(ii))
    rhoi = rhoh(xyzh(4,ii),pmassi)
    if (rhoi < rho_crit_cgs) cycle over_mins

    ! get the refinement level of the particle in the middle of the potential
    call get_apr(xyzh(1:3,ii),apri)
    if ((ref_dir == -1) .and. (apri == apr_max) .and. (ntrack<1)) then
       ! it's a newly identified clump, time to derefine it
       ntrack = ntrack + 1
       track_part = ii
    else
       ! it's an existing clump, update the position of it's centre
       track_part = ii
    endif
 enddo over_mins
 if (ntrack > 0) call set_apr_centre(apr_type,apr_centre,ntrack,track_part)
 print*,'tracking ',track_part,ntrack

 ! tidy up
 deallocate(counter,ave_poten,radius,minima,min_particle)

end subroutine create_or_update_apr_clump

end module apr
