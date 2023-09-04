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
  !   - apr_max           : number of resolution levels
  !   - ref_dir           : increase (0) or decrease (1) resolution from the base resolution
  !   - [x,y,z]_centre    : centre coordinates of the region to be more highly resolved
  !   - apr_rad           : radius of the region to be more highly resolved
  !
  ! :Dependencies: None
  !
  implicit none

  public :: init_apr,update_apr,read_options_apr,write_options_apr

  private
  integer :: apr_max = 2
  integer :: ref_dir = 0
  real    :: x_centre = 0.0, y_centre = 0.0, z_centre = 0.0
  real    :: apr_rad = 0.2
  real, allocatable    :: apr_regions(:)
  real    :: sep_factor = 0.2
  integer, save :: looped_through = 0
  logical :: do_relax = .false.

contains

  !-----------------------------------------------------------------------
  !+
  !  Initialising all the apr arrays and properties
  !+
  !-----------------------------------------------------------------------
  subroutine init_apr(apr_level,ierr)
    use dim, only:maxp_hard
    integer, intent(inout) :: ierr,apr_level(:)

    ! initialise the base resolution level
    if (ref_dir == 0) then
      apr_level = 1
    else
      apr_level = apr_max
    endif

    ! initiliase the regions
    allocate(apr_regions(apr_max))
    apr_regions(apr_max) = apr_rad
    apr_regions(1) = 2.0    ! TBD: this should be replaced with a routine that automagically calculates the steps safely

    ierr = 0

  end subroutine init_apr

  !-----------------------------------------------------------------------
  !+
  !  Subroutine to check if particles need to be split or merged
  !+
  !-----------------------------------------------------------------------
  subroutine update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
    use dim, only:maxp_hard
    use part,             only:ntot,isdead_or_accreted,igas,apr_massoftype,&
                               shuffle_part
    use quitdump,         only:quit
    use relaxem,          only:relax_particles
    real, intent(inout) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:)
    integer, intent(inout) :: npart,apr_level(:)
    integer :: ii,jj,kk,npartnew,nsplit_total,apri,npartold
    integer :: n_ref,nrelax,nmerge,nkilled
    real, allocatable :: xyzh_ref(:,:),force_ref(:,:),pmass_ref(:)
    real, allocatable :: xyzh_merge(:,:),vxyzu_merge(:,:)
    integer, allocatable :: relaxlist(:),mergelist(:)

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
          pmass_ref(n_ref) = apr_massoftype(igas,apr_level(ii))
          force_ref(1:3,n_ref) = fxyzu(1:3,ii)*pmass_ref(n_ref)
        endif
      enddo
    endif

    ! Do any particles need to be split?
    npartnew = npart
    npartold = npart
    nsplit_total = 0
    nrelax = 0

    do jj = 1,apr_max
      do ii = 1,npartold
        ! this is the refinement level it *should* have based
        ! on it's current position
        call get_apr(xyzh(1:2,ii),apri)
        ! if the level it should have is greater than the
        ! level it does have, increment it up one
        if (apri > apr_level(ii)) then
          call splitpart(ii,npartnew)
          ntot = npartnew
          if (do_relax) then
            nrelax = nrelax + 2
            relaxlist(nrelax-1) = ii
            relaxlist(nrelax)   = npartnew
          endif
        endif
      enddo
    enddo

    ! Take into account all the added particles
    print*,'split: ',(npartnew-npartold)
    npart = npartnew

    ! Do any particles need to be merged?
    allocate(mergelist(npart),xyzh_merge(4,npart),vxyzu_merge(4,npart))
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
        endif
      enddo
      ! Now send them to be merged
      if (nmerge > 0) call merge_with_special_tree(nmerge,mergelist,xyzh_merge(:,1:nmerge),&
                                              vxyzu_merge(:,1:nmerge),kk,xyzh,apr_level,nkilled,&
                                              nrelax,relaxlist,npartnew)
      print*,'merged: ',nkilled,kk
    enddo
    ! update npart as required
    npart = npartnew

    ! If we need to relax, do it here
    if (nrelax > 0) call relax_particles(npart,n_ref,xyzh_ref,force_ref,nrelax,relaxlist)
    ! Turn it off now because we only want to do this on first splits
    do_relax = .false.

    ! As we may have killed particles, time to do an array shuffle
    call shuffle_part(npart)

    ! Do my hacky write
    looped_through = looped_through + 1
    call nasty_write(looped_through,npart,xyzh(:,:),fxyzu(:,:),apr_level(:))


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
    real, intent(in)     :: pos(2)
    integer, intent(out) :: apri
    integer :: jj, kk
    real :: dx,dy,r

    do jj = 1,apr_max
      kk = apr_max - jj + 1       ! going from apr_max -> 1
      dx = pos(1) - x_centre
      dy = pos(2) - y_centre
      !  dz = pos(3) - z_centre
      r = sqrt(dx**2 + dy**2)
      if (r < apr_regions(kk)) then
        apri = kk
        return
      endif
    enddo

    print*,'function get_apr did not find a level'

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
    integer :: j,npartold,aprnew
    real :: theta,dx,dy,x_add,y_add

    ! calculate the angle from this particle to the centre of the region
    ! we will split and then rotate the particle positions through this angle
    dx = xyzh(1,i) - x_centre
    dy = xyzh(2,i) - y_centre
    theta = atan2(dy,dx) + 0.5*pi
    x_add = sep_factor*cos(theta)*xyzh(4,i)
    y_add = sep_factor*sin(theta)*xyzh(4,i)

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
    use part,     only:kill_particle,npartoftype,get_partinfo,iphase
    integer, intent(inout) :: nmerge,apr_level(:),nkilled,nrelax,relaxlist(:),npartnew
    integer, intent(in)    :: current_apr,mergelist(:)
    real, intent(inout)    :: xyzh(:,:)
    real, intent(inout)    :: xyzh_merge(:,:),vxyzu_merge(:,:)
    integer :: remainder,icell,i,n_cell,apri,m
    integer :: eldest,tuther,iamtypei
    logical :: iactive,isgas,isdust
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
      call get_apr(com(1:2),apri)
      ! If the apr level based on the com is lower than the current level,
      ! we merge!
      if (apri < current_apr) then
        eldest = mergelist(inodeparts(inoderange(1,icell)))
        tuther = mergelist(inodeparts(inoderange(2,icell)))
        call get_partinfo(iphase(tuther),iactive,isgas,isdust,iamtypei)

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
        nkilled = nkilled + 1
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
      read(valstring,*,iostat=ierr) apr_max
      ngot = ngot + 1
      if (apr_max  <  1) call fatal(label,'apr_max < 1 in input options')
    case('ref_dir')
      read(valstring,*,iostat=ierr) ref_dir
      ngot = ngot + 1
    case('x_centre')
      read(valstring,*,iostat=ierr) x_centre
      ngot = ngot + 1
    case('y_centre')
      read(valstring,*,iostat=ierr) y_centre
      ngot = ngot + 1
    case('z_centre')
      read(valstring,*,iostat=ierr) z_centre
      ngot = ngot + 1
    case('apr_rad')
      read(valstring,*,iostat=ierr) apr_rad
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

    call write_inopt(apr_max,'apr_max','maximum number of resolution levels',iunit)
    call write_inopt(ref_dir,'ref_dir','increase (0) or decrease (1) resolution',iunit)
    call write_inopt(x_centre,'x_centre','x pos of region',iunit)
    call write_inopt(y_centre,'y_centre','y pos of region',iunit)
    call write_inopt(z_centre,'z_centre','z pos of region',iunit)
    call write_inopt(apr_rad,'apr_rad','radius of region',iunit)

  end subroutine write_options_apr

  !-----------------------------------------------------------------------
  !+
  !  Hacky routine that just writes out the raw position, mass, density and h
  !+
  !-----------------------------------------------------------------------
  subroutine nasty_write(ifile,npart,xyzh,fxyzu,apr_level)
    use part, only:igas,rhoh,apr_massoftype
    real, intent(in) :: xyzh(:,:),fxyzu(:,:)
    integer, intent(in) :: apr_level(:),ifile,npart
    integer :: ii,iunit=24,apri
    character(len=120) :: mydumpfile,rootname
    real :: rhoi,hfact,pmass

    hfact = 1.2 ! straight from the *.in file

    rootname = 'rawinfo'
    write(mydumpfile,"(a,'_',i5.5,'.dat')") trim(rootname),ifile

    open(iunit,file=mydumpfile,status='replace')
    write(iunit,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
    1,'x', &
    2,'y', &
    3,'h', &
    4,'rho', &
    5,'mass',&
    6,'fx',&
    7,'fy',&
    8,'fz',&
    9,'apri'

    do ii = 1,npart
      apri = apr_level(ii)
      pmass = apr_massoftype(igas,apri)
      rhoi = rhoh(xyzh(4,ii),pmass)
      write(iunit,'(9(es18.10,1X))') xyzh(1:2,ii),xyzh(4,ii),rhoi,pmass,fxyzu(1:3,ii),real(apri)
    enddo

    close(iunit)

  end subroutine nasty_write

end module apr
