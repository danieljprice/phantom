!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_apr
!
! Utility routines for APR
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters:
!   - apr_drad     : *size of step to next region*
!   - apr_max      : *number of additional refinement levels (3 -> 2x resolution)*
!   - apr_rad      : *radius of innermost region*
!   - apr_type     : *1: static, 2: sink, 3: clumps, 4: sequential sinks, 5: com, 6: vertical*
!   - ref_dir      : *increase (1) or decrease (-1) resolution*
!   - rho_crit_cgs : *density above which apr zones are created (g/cm^3)*
!   - track_part   : *number of sink to track*
!
! :Dependencies: infile_utils, io, part, ptmass
!

 use ptmass, only:rho_crit_cgs

 implicit none

 public :: read_options_apr,write_options_apr

 ! default values for runtime parameters are stored here
 integer :: apr_max_in = 3, ref_dir = 1, apr_type = 1, apr_max = 4
 integer :: top_level = 1, ntrack = 1, ntrack_max = 10
 integer :: icentre = 1, track_part_in = 1
 integer, allocatable :: npart_regions(:), track_part(:)
 real :: apr_rad = 1.0, apr_drad = 0.1, apr_centre_in(3) = 0.
 real, allocatable :: apr_regions(:), apr_centre(:,:)
 real, save :: apr_H(2,100)  ! we enforce this to be 100

 logical :: apr_region_is_circle = .false.

contains

!-----------------------------------------------------------------------
!+
!  Find inner and outer radius - assuming that the disc is around the
!  central star for now
!+
!-----------------------------------------------------------------------
subroutine find_inner_and_outer_radius(npart,xyzh,rmin,rmax)
 use part, only: xyzmh_ptmass, isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:)
 real, intent(out)   :: rmin,rmax
 integer :: ii
 real    :: rmin_test, rmax_test, xi, yi, zi, r2_test

 rmin_test = huge(rmin_test)
 rmax_test = tiny(rmax_test) ! just big and small initial guesses


 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(npart,xyzh,xyzmh_ptmass,rmin_test,rmax_test) &
 !$omp private(ii,xi,yi,zi,r2_test)
 do ii = 1,npart
    if (isdead_or_accreted(xyzh(4,ii))) cycle
    xi = xyzh(1,ii) - xyzmh_ptmass(1,1)
    yi = xyzh(2,ii) - xyzmh_ptmass(2,1)
    zi = xyzh(3,ii) - xyzmh_ptmass(3,1)
    r2_test = xi**2 + yi**2 + zi**2

    if (r2_test < rmin_test) rmin_test = r2_test
    if (r2_test > rmax_test) rmax_test = r2_test

 enddo
 !$omp end parallel do

 rmin = sqrt(rmin_test)
 rmax = sqrt(rmax_test)

end subroutine find_inner_and_outer_radius

!-----------------------------------------------------------------------
!+
!  routine to find the closest apr centre to a position
!+
!-----------------------------------------------------------------------

subroutine find_closest_region(pos,iclosest)
 real, intent(in) :: pos(3)
 integer, intent(out) :: iclosest
 real :: r2,rtest,dx,dy,dz
 integer :: ii

 rtest = huge(rtest)

 iclosest = -1

 if (ntrack == 0) return ! nothing to do here

 do ii = 1,ntrack
    dx = pos(1) - apr_centre(1,ii)
    dy = pos(2) - apr_centre(2,ii)
    dz = pos(3) - apr_centre(3,ii)
    r2 = dx**2 + dy**2 + dz**2
    if (r2 < rtest) then
       iclosest = ii
       rtest = r2
    endif
 enddo

end subroutine find_closest_region


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
 logical :: igotall1,igotall2,igotall3

 imatch   = .true.
 igotall1 = .true.
 igotall2 = .true.
 igotall3 = .true.
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
    case(3)
       call read_options_apr3(name,valstring,imatch,igotall3,ierr)
    end select
 end select
 igotall = (ngot >= 5) .and. igotall1 .and. igotall2 .and. igotall3

end subroutine read_options_apr

!-----------------------------------------------------------------------
!+
!  extra subroutines for reading in different styles of apr zones
!  subroutine 1: when you need to send in the x,y,z position
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
    read(valstring,*,iostat=ierr) apr_centre_in(1)
    ngot = ngot + 1
 case('apr_centre(2)')
    read(valstring,*,iostat=ierr) apr_centre_in(2)
    ngot = ngot + 1
 case('apr_centre(3)')
    read(valstring,*,iostat=ierr) apr_centre_in(3)
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 3)

end subroutine read_options_apr1

!-----------------------------------------------------------------------
!+
!  extra subroutines for reading in different styles of apr zones
!  subroutine 2: when you need to track a particular particle
!+
!-----------------------------------------------------------------------

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
    read(valstring,*,iostat=ierr) track_part_in
    ngot = ngot + 1
    if (track_part_in  <  1) call fatal(label,'track_part not chosen in input options')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)

end subroutine read_options_apr2

!-----------------------------------------------------------------------
!+
!  extra subroutines for reading in different styles of apr zones
!  subroutine 3: setting a threshold density
!+
!-----------------------------------------------------------------------

subroutine read_options_apr3(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_apr3'

 imatch  = .true.
 select case(trim(name))
 case('rho_crit_cgs')
    read(valstring,*,iostat=ierr) rho_crit_cgs
    if (rho_crit_cgs < 0.) call fatal(label,'rho_crit < 0')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)

end subroutine read_options_apr3

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
 call write_inopt(apr_type,'apr_type','1: static, 2: sink, 3: clumps, 4: sequential sinks, 5: com, 6: vertical',iunit)

 select case (apr_type)
 case (1)
    call write_inopt(apr_centre_in(1),'apr_centre(1)','centre of region x position',iunit)
    call write_inopt(apr_centre_in(2),'apr_centre(2)','centre of region y position',iunit)
    call write_inopt(apr_centre_in(3),'apr_centre(3)','centre of region z position',iunit)
 case(2,4)
    call write_inopt(track_part_in,'track_part','number of sink to track',iunit)
 case(3)
    call write_inopt(rho_crit_cgs,'rho_crit_cgs','density above which apr zones are created (g/cm^3)',iunit)
 case default
    ! write nothing
 end select

 call write_inopt(apr_rad,'apr_rad','radius of innermost region',iunit)
 call write_inopt(apr_drad,'apr_drad','size of step to next region',iunit)

end subroutine write_options_apr

end module utils_apr
