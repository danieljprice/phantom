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
  real    :: sep_factor = 0.35
  integer, save :: looped_through = 0

contains

  !-----------------------------------------------------------------------
  !+
  !  Initialising all the apr arrays and properties
  !+
  !-----------------------------------------------------------------------
  subroutine init_apr(apr_level,ierr)
    use dim, only:maxp_hard
    use memory,          only:allocate_memory
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
  subroutine update_apr(npart,xyzh,vxyzu,apr_level)
    use part, only:ntot,npartoftype
    use quitdump,         only:quit
    real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
    integer, intent(inout) :: npart,apr_level(:)
    integer :: ii,jj,npartnew,nsplit_total,apri,npartold

    ! Do any particles need to be split?
    npartnew = npart
    npartold = npart
    nsplit_total = 0

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
        endif
      enddo
    enddo
    ! Take into account all the added particles
    print*,'split: ',(npartnew-npartold)
    npart = npartnew

    looped_through = looped_through + 1

    ! Do any particles need to be merged?

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
      xyzh(4,j) = xyzh(4,i)*(2)**(1./3.)
    enddo

    ! Edit the old particle that was sent in and kept
    xyzh(1,i) = xyzh(1,i) - x_add
    xyzh(2,i) = xyzh(2,i) - y_add
    apr_level(i) = aprnew
    xyzh(4,i) = xyzh(4,i)*2.0**(1./3.)

  end subroutine splitpart

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

end module apr
