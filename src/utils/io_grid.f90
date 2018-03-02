!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: io_grid
!
!  DESCRIPTION:
!  module for read/write of gridded data to/from file
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: fileutils, hdf5utils, io
!+
!--------------------------------------------------------------------------
module io_grid
 implicit none
 integer, parameter, public :: nformats = 1
 character(len=4), dimension(nformats) :: labelformat = &
      (/'hdf5'/)

 public :: read_grid_header,read_grid_column

contains

!----------------------------------------------------------------
!+
!  routine to print available file formats
!+
!----------------------------------------------------------------
subroutine print_grid_formats()
 integer :: i

 do i=1,nformats
    if (i==1) then
       print "(i2,') ',a)",i,trim(labelformat(i))//'    (default if no format specified)'
    else
       print "(i2,') ',a)",i,trim(labelformat(i))
    endif
 enddo

end subroutine print_grid_formats

!----------------------------------------------------------------
!+
!  routine to determine file format based on filenames
!
!  returns the stripped base filename for the file, without the
!  extension, in outfile. e.g. file_longd.struct would return
!  'file' as the stripped filename.
!+
!----------------------------------------------------------------
integer function get_grid_format(filename,outfile)
 use fileutils, only:lcase
 character(len=*), intent(in)  :: filename
 character(len=*), intent(out) :: outfile

 outfile = trim(filename)
 get_grid_format = 0
 if (index(lcase(filename),'.h5') /= 0) then
    get_grid_format = 1
 endif

end function get_grid_format

!----------------------------------------------------------------
!+
!  interface routine to read gridded data of any format
!+
!----------------------------------------------------------------
subroutine read_grid_header(gridformat,filename,nx,ny,nz,ncols,time,ierr)
 use hdf5utils, only:read_grid_hdf5_header
 use io,        only:cstring
 character(len=*), intent(in)  :: gridformat,filename
 integer,          intent(out) :: nx,ny,nz,ncols,ierr
 real,             intent(out) :: time

 ierr = 0
 time = 0.
 nx   = 0
 ny   = 0
 nz   = 0
 ncols = 0
 select case(gridformat)
 case('hdf5')
    call read_grid_hdf5_header(cstring(filename),nx,ny,nz,ncols,ierr)
 case default
    print*,' input file format for '//trim(filename)//' not recognised, skipping...'
    ierr = 1
    return
 end select

end subroutine read_grid_header

!----------------------------------------------------------------
!+
!  interface routine to read gridded data of any format
!+
!----------------------------------------------------------------
subroutine read_grid_column(gridformat,filename,icolumn,nx,ny,nz,datcol,ierr)
 use hdf5utils, only:read_grid_hdf5_column
 use io,        only:cstring
 character(len=*), intent(in)  :: gridformat,filename
 integer,          intent(in)  :: icolumn,nx,ny,nz
 real,             intent(out) :: datcol(:)
 integer,          intent(out) :: ierr

 select case(gridformat)
 case('hdf5')
    call read_grid_hdf5_column(cstring(filename),icolumn,nx,ny,nz,datcol,ierr)
 case default
    print*,' input file format for '//trim(filename)//' not recognised, skipping...'
    ierr = 1
    return
 end select

end subroutine read_grid_column

end module io_grid
