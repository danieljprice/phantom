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
!   - split_dir    : *1: tangent to boundary, 2: along trajectory, 3: purely randomly*
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
 integer :: split_dir = 1
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
 use part, only:xyzmh_ptmass, isdead_or_accreted
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
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_apr(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for adaptive particle refinement'
 call write_inopt(apr_max_in,'apr_max','number of additional refinement levels (3 -> 2x resolution)',iunit)
 call write_inopt(ref_dir,'ref_dir','increase (1) or decrease (-1) resolution',iunit)
 call write_inopt(split_dir,'split_dir','1: tangent to boundary, 2: along trajectory, 3: purely randomly',iunit)
 call write_inopt(apr_type,'apr_type','1: static, 2: sink, 3: clumps, 4: sequential sinks, 5: com, 6: vertical',iunit)

 select case (apr_type)
 case(1)
    call write_inopt(apr_centre_in(1),'apr_centre(1)','centre of region x position',iunit)
    call write_inopt(apr_centre_in(2),'apr_centre(2)','centre of region y position',iunit)
    call write_inopt(apr_centre_in(3),'apr_centre(3)','centre of region z position',iunit)
 case(2,4)
    call write_inopt(track_part_in,'track_part','number of sink to track',iunit)
 case(3)
    call write_inopt(rho_crit_cgs,'rho_crit_cgs','density above which apr zones are created (g/cm^3)',iunit)
 end select

 call write_inopt(apr_rad,'apr_rad','radius of innermost region',iunit)
 call write_inopt(apr_drad,'apr_drad','size of step to next region',iunit)

end subroutine write_options_apr

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_apr(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(apr_max_in,'apr_max',db,errcount=nerr,min=0)
 call read_inopt(ref_dir,'ref_dir',db,errcount=nerr,min=-1,max=1)
 call read_inopt(split_dir,'split_dir',db,errcount=nerr,min=1,max=3)
 call read_inopt(apr_type,'apr_type',db,errcount=nerr,min=1,max=6)

 select case(apr_type)
 case(1)
    call read_inopt(apr_centre_in(1),'apr_centre(1)',db,errcount=nerr)
    call read_inopt(apr_centre_in(2),'apr_centre(2)',db,errcount=nerr)
    call read_inopt(apr_centre_in(3),'apr_centre(3)',db,errcount=nerr)
 case(2,4)
    call read_inopt(track_part_in,'track_part',db,errcount=nerr,min=1)
 case(3)
    call read_inopt(rho_crit_cgs,'rho_crit_cgs',db,errcount=nerr,min=0.)
 end select

 call read_inopt(apr_rad,'apr_rad',db,errcount=nerr,min=tiny(apr_rad))
 call read_inopt(apr_drad,'apr_drad',db,errcount=nerr,min=tiny(apr_drad))

end subroutine read_options_apr

!-----------------------------------------------------------------------
!+
!  Writes tracking file for APR regions
!+
!-----------------------------------------------------------------------
subroutine write_aprtrack(tdump,dumpfile)
 use io, only:iaprdump,error
 real,             intent(in) :: tdump
 character(len=*), intent(in) :: dumpfile
 integer :: ierr, i, j
 character(len=10) :: filename
 character(len=3)  :: padded_ntrack
 character(len=11) :: label
 character(len=256) :: fmt
 integer :: dumpfile_int, dump_length, start_pos
 logical :: iexist

 if (ntrack == 0) return ! nothing to do here

 ! clever formatting
 dump_length = len_trim(dumpfile)
 start_pos = dump_length - 5 + 1
 read(dumpfile(start_pos:dump_length), *) dumpfile_int
 ! dynamically make the formatting string to accomodate the right number of regions
 fmt = '(ES18.10,2X,I16.5,2X,3(ES18.10,1X)'
 do j = 1,apr_max_in - 1
    fmt = trim(fmt) // ',ES18.10,1X'
 enddo
 fmt = trim(fmt) // ',ES18.10)'

 do i = 1,ntrack
    write(padded_ntrack, '(I3.3)') i
    filename = 'apr_' // padded_ntrack // '.ev'

    ! check if the file exists or not
    inquire(file=filename,exist=iexist)
    if (.not.iexist .or. (tdump < tiny(tdump))) then
       ! create a new file
       open(unit=iaprdump,file=filename,status='replace',form='formatted',iostat=ierr)
       write(iaprdump, '("# APR info for region ",i3)') i
       write(iaprdump,"('#',5(1x,'[',i2.2,1x,a11,']',2x))",advance="no") &
          1,'time', &
          2,'dump', &
          3,'x centre', &
          4,'y centre', &
          5,'z centre'
       do j = 1,apr_max-1
          write(label, '(A7,I0)') 'radius_', j  ! the different radii

          write(iaprdump, "(1x,'[',i2.2,1x,a11,']',2x)", advance="no") 5 + j, label
       enddo
       write(iaprdump,*)
    else
       ! append the existing file
       open(unit=iaprdump,file=filename,status='old',form='formatted',position='append',iostat=ierr)
    endif

    if (ierr /= 0) then
       call error('write_aprtrack','could not open APR tracking file for writing')
    else
       if (ref_dir == 1) then
          write(iaprdump,fmt) tdump,dumpfile_int,apr_centre(1:3,i),apr_regions(2:apr_max)
       else
          write(iaprdump,fmt) tdump,dumpfile_int,apr_centre(1:3,i),apr_regions(1:apr_max-1)
       endif

    endif
    close(unit=iaprdump)
 enddo

end subroutine write_aprtrack

end module utils_apr
