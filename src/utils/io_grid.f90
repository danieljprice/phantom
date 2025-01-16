!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module io_grid
!
! module for read/write of gridded data to/from file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: fileutils, hdf5utils, io, readwrite_griddata
!
 implicit none
 integer, parameter :: doub_prec = kind(0.d0)
 integer, parameter, public :: nformats = 4
 character(len=6), dimension(nformats) :: labelformat = &
      (/'ascii ','binary','stream','hdf5  '/)

 public :: read_grid_header,read_grid_column
 integer, parameter, private :: io_unit = 85

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
       print "(1x,a)",'"'//trim(labelformat(i))//'"    (default if no format specified)'
    else
       print "(1x,a)",'"'//trim(labelformat(i))//'"'
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
    get_grid_format = 4
 elseif (index(lcase(filename),'.gridstream') /= 0) then
    get_grid_format = 3
 elseif (index(lcase(filename),'.grid') /= 0) then
    get_grid_format = 2
 endif

end function get_grid_format

!----------------------------------------------------------------
!+
!  interface routine to read gridded data of any format
!+
!----------------------------------------------------------------
subroutine read_grid_header(gridformat,filename,nx,ny,nz,ncols,time,ierr)
 use hdf5utils,          only:read_grid_hdf5_header
 use readwrite_griddata, only:isgridformat,open_gridfile_r
 use io,                 only:cstring
 character(len=*), intent(in)  :: gridformat,filename
 integer,          intent(out) :: nx,ny,nz,ncols,ierr
 real(doub_prec),  intent(out) :: time
 integer :: npix(3)

 ierr = 0
 time = 0.d0
 nx   = 0
 ny   = 0
 nz   = 0
 ncols = 0
 select case(gridformat)
 case('hdf5')
    call read_grid_hdf5_header(cstring(filename),nx,ny,nz,ncols,ierr)
 case default
    if (isgridformat('grid'//gridformat)) then
       call open_gridfile_r(io_unit,filename,'grid'//trim(gridformat),3,ncols,npix,time,ierr)
       nx = npix(1); ny = npix(2); nz = npix(3)
    else
       print*,' input file format for '//trim(filename)//' not recognised, skipping...'
       ierr = 1
    endif
 end select

end subroutine read_grid_header

!----------------------------------------------------------------
!+
!  interface routine to read gridded data of any format
!+
!----------------------------------------------------------------
subroutine read_grid_column(gridformat,filename,icolumn,nx,ny,nz,datcol,ierr)
 use hdf5utils,          only:read_grid_hdf5_column
 use readwrite_griddata, only:read_gridcolumn,isgridformat
 use io,                 only:cstring
 character(len=*), intent(in)  :: gridformat,filename
 integer,          intent(in)  :: icolumn,nx,ny,nz
 real(doub_prec),  intent(out) :: datcol(:)
 integer,          intent(out) :: ierr

 select case(gridformat)
 case('hdf5')
    call read_grid_hdf5_column(cstring(filename),icolumn,nx,ny,nz,datcol,ierr)
 case default
    if (isgridformat('grid'//gridformat)) then
       call read_gridcolumn(io_unit,datcol,nx*ny*nz,ierr)
    else
       print*,' input file format for '//trim(filename)//' not recognised, skipping...'
       ierr = 1
    endif
 end select

end subroutine read_grid_column

end module io_grid
