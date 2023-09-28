!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr
  !
  ! Contains everything for live adaptive particle refinement
  !
  ! :References: None
  !
  ! :Owner: Rebecca Nealon
  !
  ! :Runtime parameters:
  !   - apr_max_in        : number of refinement levels (3 -> 2x resolution)
  !   - ref_dir           : increase (1) or decrease (-1) resolution from the base resolution
  !   - apr_type          : choice of region, defined in apr_region.f90
  !
  ! :Dependencies: None
  !
  implicit none

  public :: init_apr,update_apr,read_options_apr,write_options_apr,hacky_write

  private
  integer :: apr_max, apr_max_in = 3
  integer :: ref_dir = 1, top_level = 1, apr_type = 1
  real    :: apr_centre(3)
  real, allocatable    :: apr_regions(:)
  integer, allocatable :: npart_regions(:)
  real    :: sep_factor = 0.25
  logical :: apr_verbose = .false.
  logical :: do_relax = .true.
  logical :: adjusted_split = .true.

contains

  !-----------------------------------------------------------------------
  !+
  !  Initialising all the apr arrays and properties
  !+
  !-----------------------------------------------------------------------
  subroutine init_apr(apr_level,apr_weight,ierr)
    use dim, only:maxp_hard
    use part, only:npart,massoftype
    use apr_region, only:set_apr_region
    integer, intent(inout) :: ierr,apr_level(:)
    real, intent(inout)    :: apr_weight(:)
    logical :: previously_set
    real    :: apr_rad

    ! the resolution levels are in addition to the base resolution
    apr_max = apr_max_in + 1

    ! if we're reading in a file that already has the levels set,
    ! don't override these
    previously_set = .false.
    if (sum(apr_level(1:npart)) > npart) then
      previously_set = .true.
      do_relax = .false.
    endif

    if (.not.previously_set) then
      ! initialise the base resolution level
      if (ref_dir == 1) then
        apr_level = 1
      else
        apr_level = apr_max
      endif
    endif

    ! initiliase the regions
    call set_apr_region(apr_type,apr_centre,apr_rad)
    allocate(apr_regions(apr_max),npart_regions(apr_max))
    if (ref_dir == 1) then
      apr_regions(apr_max) = apr_rad
      apr_regions(1) = 2.0    ! TBD: this should be replaced with a routine that automagically calculates the steps safely
!      apr_regions(2) = 0.4
!      apr_regions(3) = 0.3
      top_level = apr_max
    else
      apr_regions(1) = apr_rad
      apr_regions(apr_max) = 2.0    ! TBD: this should be replaced with a routine that automagically calculates the steps safely
      apr_regions(2) = 0.3
      top_level = 1
    endif
    npart_regions = 0

    ! if we are derefining we make sure that
    ! massoftype(igas) is associated with the
    ! largest particle
    if (ref_dir == -1) then
      massoftype(:) = massoftype(:) * 2.**(apr_max -1)
    endif

    ! set the initial weightings
    apr_weight(:) = 1.0

    ierr = 0

  end subroutine init_apr

  !-----------------------------------------------------------------------
  !+
  !  Subroutine to check if particles need to be split or merged
  !+
  !-----------------------------------------------------------------------
  subroutine update_apr(npart,xyzh,vxyzu,fxyzu,apr_level,apr_weight)
    use dim,      only:maxp_hard
    use part,     only:ntot,isdead_or_accreted,igas,massoftypefunc,&
                       shuffle_part
    use quitdump, only:quit
    use relaxem,  only:relax_particles
    real, intent(inout) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),apr_weight(:)
    integer, intent(inout) :: npart,apr_level(:)
    integer :: ii,jj,kk,npartnew,nsplit_total,apri,npartold
    integer :: n_ref,nrelax,nmerge,nkilled,apr_current
    real, allocatable :: xyzh_ref(:,:),force_ref(:,:),pmass_ref(:)
    real, allocatable :: xyzh_merge(:,:),vxyzu_merge(:,:)
    integer, allocatable :: relaxlist(:),mergelist(:)
    real :: xi,yi,zi

    ! If this routine doesn't need to be used, just skip it
    if (apr_max == 1) return

    print*,apr_weight(1:10)
    read*

    ! Before adjusting the particles, if we're going to
    ! relax them then let's save the reference particles
    if (do_relax) then
      allocate(xyzh_ref(4,maxp_hard),force_ref(3,maxp_hard),pmass_ref(maxp_hard),relaxlist(maxp_hard))
      relaxlist = -1

      n_ref = 0
      xyzh_ref = 0.
      force_ref = 0.
      pmass_ref = 0.

      do ii = 1,npart
        if (.not.isdead_or_accreted(xyzh(4,ii))) then ! ignore dead particles
          n_ref = n_ref + 1
          xyzh_ref(1:4,n_ref) = xyzh(1:4,ii)
          pmass_ref(n_ref) = massoftypefunc(igas,apr_level(ii))
          force_ref(1:3,n_ref) = fxyzu(1:3,ii)*pmass_ref(n_ref)
        endif
      enddo
    endif

    ! Do any particles need to be split?
    npartnew = npart
    npartold = npart
    nsplit_total = 0
    nrelax = 0
    apri = 0 ! to avoid compiler errors


    do jj = 1,apr_max
      npartold = npartnew ! to account for new particles as they are being made

      do ii = 1,npartold
        apr_current = apr_level(ii)
        xi = xyzh(1,ii)
        yi = xyzh(2,ii)
        zi = xyzh(3,ii)
        ! this is the refinement level it *should* have based
        ! on it's current position
        call get_apr((/xi,yi,zi/),apri)
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
    allocate(mergelist(npart),xyzh_merge(4,npart),vxyzu_merge(4,npart))
    npart_regions = 0
    do jj = 1,apr_max-1
      kk = apr_max - jj + 1             ! to go from apr_max -> 2
      mergelist = -1 ! initialise
      nmerge = 0
      nkilled = 0
      xyzh_merge = 0.
      vxyzu_merge = 0.

      do ii = 1,npart
        ! note that here we only do this process for particles that are not already counted in the blending region
        if ((apr_level(ii) == kk) .and. (.not.isdead_or_accreted(xyzh(4,ii)))) then ! avoid already dead particles
          nmerge = nmerge + 1
          mergelist(nmerge) = ii
          xyzh_merge(1:4,nmerge) = xyzh(1:4,ii)
          vxyzu_merge(1:3,nmerge) = vxyzu(1:3,ii)
          npart_regions(kk) = npart_regions(kk) + 1
        endif
      enddo

      ! Now send them to be merged
      if (nmerge > 0) call merge_with_special_tree(nmerge,mergelist,xyzh_merge(:,1:nmerge),&
      vxyzu_merge(:,1:nmerge),kk,xyzh,apr_level,nkilled,&
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
      deallocate(xyzh_ref,force_ref,pmass_ref,relaxlist)
    endif
    deallocate(mergelist)

  end subroutine update_apr

  !-----------------------------------------------------------------------
  !+
  !  routine to return the adaptive particle refinement level based on position
  !  and the boundaries set by the apr_* arrays
  !+
  !-----------------------------------------------------------------------
  subroutine get_apr(pos,apri)
    real, intent(in)     :: pos(3)
    integer, intent(out) :: apri
    integer :: jj, kk
    real :: dx,dy,dz,r

    do jj = 1,apr_max
      if (ref_dir == 1) then
        kk = apr_max - jj + 1       ! going from apr_max -> 1
      else
        kk = jj                    ! going from 1 -> apr_max
      endif
      dx = pos(1) - apr_centre(1)
      dy = pos(2) - apr_centre(2)
      dz = pos(3) - apr_centre(3)
      r = sqrt(dx**2 + dy**2)
      if (r < apr_regions(kk)) then
        apri = kk
        return
      endif
    enddo

    print*,'get_apr did not find apri'

  end subroutine get_apr

  !-----------------------------------------------------------------------
  !+
  !  routine to split one particle into two
  !+
  !-----------------------------------------------------------------------
  subroutine splitpart(i,npartnew)
    use part,    only:copy_particle_all,apr_level,xyzh,vxyzu,npartoftype,igas
    use part,    only:set_particle_type
    use physcon, only:pi
    integer, intent(in) :: i
    integer, intent(inout) :: npartnew
    integer :: j,npartold,aprnew,next_door
    real :: theta,dx,dy,x_add,y_add,sep,rneigh

    if (adjusted_split) then
      call closest_neigh(i,next_door,rneigh)
      sep = min(sep_factor*xyzh(4,i),0.35*rneigh)
      sep = sep/xyzh(4,i)  ! for consistency later on
    else
      sep = sep_factor
    endif

    ! calculate the angle from this particle to the centre of the region
    ! we will split and then rotate the particle positions through this angle
    dx = xyzh(1,i) - apr_centre(1)
    dy = xyzh(2,i) - apr_centre(2)
    theta = atan2(dy,dx) + 0.5*pi
    x_add = sep*cos(theta)*xyzh(4,i)
    y_add = sep*sin(theta)*xyzh(4,i)

    npartold = npartnew
    npartnew = npartold + 1
    npartoftype(igas) = npartoftype(igas) + 1
    aprnew = apr_level(i) + 1

    !--create the new particle
    do j=npartold+1,npartnew
      call copy_particle_all(i,j,new_part=.true.)
      xyzh(1,j) = xyzh(1,i) + x_add
      xyzh(2,j) = xyzh(2,i) + y_add
      vxyzu(:,j) = vxyzu(:,i)
      apr_level(j) = aprnew
      xyzh(4,j) = xyzh(4,i)*(0.5**(1./3.))
    enddo

    ! Edit the old particle that was sent in and kept
    xyzh(1,i) = xyzh(1,i) - x_add
    xyzh(2,i) = xyzh(2,i) - y_add
    apr_level(i) = aprnew
    xyzh(4,i) = xyzh(4,i)*(0.5**(1./3.))

  end subroutine splitpart

  !-----------------------------------------------------------------------
  !+
  !  Take in all particles that *might* be merged at this apr_level
  !  and use our special tree to merge what has left the region
  !+
  !-----------------------------------------------------------------------
  subroutine merge_with_special_tree(nmerge,mergelist,xyzh_merge,vxyzu_merge,current_apr,&
                                     xyzh,apr_level,nkilled,nrelax,relaxlist,npartnew)
    use linklist, only:set_linklist,ncells,ifirstincell,get_cell_location
    use mpiforce, only:cellforce
    use kdtree,   only:inodeparts,inoderange
    use part,     only:kill_particle,npartoftype
    integer, intent(inout) :: nmerge,apr_level(:),nkilled,nrelax,relaxlist(:),npartnew
    integer, intent(in)    :: current_apr,mergelist(:)
    real, intent(inout)    :: xyzh(:,:)
    real, intent(inout)    :: xyzh_merge(:,:),vxyzu_merge(:,:)
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
      if (i <= 0) cycle over_cells !--skip empty cells AND inactive cells
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
        tuther = mergelist(inodeparts(inoderange(2,icell)))

        ! keep eldest, reassign it
        xyzh(1,eldest) = cell%xpos(1)
        xyzh(2,eldest) = cell%xpos(2)
        xyzh(3,eldest) = cell%xpos(3)
        xyzh(4,eldest) = xyzh(4,eldest)*(2.0**(1./3.))
        apr_level(eldest) = apr_level(eldest) - 1

        ! add it to the shuffling list if needed
        if (do_relax) then
          nrelax = nrelax + 1
          relaxlist(nrelax) = eldest
        endif

        ! discard tuther
        call kill_particle(tuther,npartoftype)
        nkilled = nkilled + 2
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

    imatch  = .true.
    select case(trim(name))
    case('apr_max')
      read(valstring,*,iostat=ierr) apr_max_in
      ngot = ngot + 1
      if (apr_max_in  <  1) call fatal(label,'apr_max < 1 in input options')
    case('ref_dir')
      read(valstring,*,iostat=ierr) ref_dir
      ngot = ngot + 1
    case('apr_type')
      read(valstring,*,iostat=ierr) apr_type
      ngot = ngot + 1
    case default
      imatch = .false.
    end select

    igotall = (ngot >= 1)

  end subroutine read_options_apr

  !-----------------------------------------------------------------------
  !+
  !  Writes input options to the input file.
  !+
  !-----------------------------------------------------------------------
  subroutine write_options_apr(iunit)
    use infile_utils, only:write_inopt
    integer, intent(in) :: iunit

    call write_inopt(apr_max_in,'apr_max','number of refinement levels (3 -> 2x resolution)',iunit)
    call write_inopt(ref_dir,'ref_dir','increase (1) or decrease (-1) resolution',iunit)
    call write_inopt(apr_type,'apr_type','1: static',iunit)

  end subroutine write_options_apr

  !-----------------------------------------------------------------------
  !+
  !  Hacky routine that just writes out the raw position, mass, density and h
  !+
  !-----------------------------------------------------------------------
  subroutine hacky_write(ifile)
    use part, only:igas,rhoh,massoftypefunc,Bxyz, &
                   npart,xyzh,apr_level
    use dim,  only:mhd
    character(len=*), intent(in) :: ifile
    integer :: ii,iunit=24
    character(len=120) :: mydumpfile,rootname
    real :: rhoi,pmass
    real :: hfact

    hfact = 1.2 ! straight from the *.in file

    rootname = 'raw'
    write(mydumpfile,"(a,'.raw')") trim(ifile)

    open(iunit,file=mydumpfile,status='replace',form='formatted')
    ! Now we write the header, and here we just hard-wire what we need to
  !  write(iunit,iostat=ierr) time,npart,npart,1.6,1.2,2,3, &
  !   6,1,3,3,-1.0,-1.0,1.0,1.0,12,'cartesian'
  !   if (ierr/=0) print*,'AH SOMETHING WRONG IN WRITE'
    write(iunit,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
    1,'x', &
    2,'y', &
    3,'h', &
    4,'rho', &
    5,'mass',&
    6,'apri',&
    7,'Bx',&
    8,'By'

    do ii=1,npart
      pmass = massoftypefunc(igas,apr_level(ii))
      rhoi = rhoh(xyzh(4,ii),pmass)
      if (.not.mhd) then
        write(iunit,'(8(es18.10,1X))') xyzh(1:2,ii), xyzh(4,ii), rhoi, pmass, real(apr_level(ii)), 0.0, 0.
      else
        write(iunit,'(8(es18.10,1X))') xyzh(1:2,ii), xyzh(4,ii), rhoi, pmass, real(apr_level(ii)), Bxyz(1,ii), Bxyz(2,ii)
      endif
    enddo

    close(iunit)

  end subroutine hacky_write

  subroutine closest_neigh(i,next_door,rmin)
    use part, only:xyzh,npart
    integer, intent(in)  :: i
    integer, intent(out) :: next_door
    real,    intent(out) :: rmin
    real :: dx,dy,dz,rtest
    integer :: j

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


end module apr
