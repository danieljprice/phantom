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
  !   - apr_max      : number of resolution levels
  !   - x_centre,y_centre : centre of the region to be more highly resolved
  !   - apr_rad      : radius of the region to be more highly resolved
  !
  ! :Dependencies: None
  !
  implicit none

  public :: init_apr,update_apr,read_options_apr,write_options_apr

  private
  integer :: apr_max = 1
  real    :: x_centre = 0.0, y_centre = 0.0
  real    :: apr_rad = 0.2

contains

  !--------------------------------------------------------------------
  ! Subroutine to initialise arrays
  !--------------------------------------------------------------------
  subroutine init_apr(apr_level,ierr)
    integer, intent(inout) :: ierr,apr_level(:)

    apr_level = 3

    ierr = 0

  end subroutine init_apr

  !--------------------------------------------------------------------
  ! Subroutine to check if particles need to be split or merged
  !--------------------------------------------------------------------
  subroutine update_apr(npart,xyzh,vxyzu,apr_level)
    real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
    integer, intent(inout) :: npart,apr_level(:)
    integer :: ii

    ! Do any particles need to be split?


    ! Do any particles need to be merged?

  end subroutine update_apr

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
   case('x_centre')
      read(valstring,*,iostat=ierr) x_centre
      ngot = ngot + 1
    case('y_centre')
       read(valstring,*,iostat=ierr) x_centre
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
  call write_inopt(x_centre,'x_centre','x pos of region',iunit)
  call write_inopt(y_centre,'y_centre','y pos of region',iunit)
  call write_inopt(apr_rad,'apr_rad','radius of region',iunit)

end subroutine write_options_apr

end module apr
