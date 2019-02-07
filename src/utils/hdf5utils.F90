!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: hdf5utils
!
!  DESCRIPTION:
!  This module provides Fortran interfaces to
!  the c routines that handle HDF5 calls
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module hdf5utils
 !interface to the c versions
 interface
  subroutine write_grid_hdf5(filename,datasetname,datarr,nx,ny,nz,ierr) bind(c)
   use, intrinsic :: iso_c_binding
   character(kind=c_char), intent(in)  :: filename
   character(kind=c_char), intent(in)  :: datasetname
   integer(kind=c_int),    intent(in)  :: nx,ny,nz
   real(kind=c_float),     intent(in)  :: datarr(nx*ny*nz)
   integer(kind=c_int),    intent(out) :: ierr
  end subroutine write_grid_hdf5
 end interface

 interface
  subroutine read_grid_hdf5_header(filename,nx,ny,nz,ncol,ierr) bind(c)
   use, intrinsic :: iso_c_binding
   character(kind=c_char), intent(in)  :: filename
   integer(kind=c_int),    intent(out) :: nx,ny,nz,ncol,ierr
  end subroutine read_grid_hdf5_header
 end interface

 !--generic interface for reading both real and double precision arrays
 interface read_grid_hdf5_column
  subroutine read_grid_hdf5_column(filename,ireadcol,nx,ny,nz,datarr,ierr) bind(c)
   use, intrinsic :: iso_c_binding
   character(kind=c_char), intent(in)  :: filename
   integer(kind=c_int),    intent(in)  :: ireadcol,nx,ny,nz
   real(kind=c_double),    intent(out) :: datarr(nx*ny*nz)
   integer(kind=c_int),    intent(out) :: ierr
  end subroutine read_grid_hdf5_column

  subroutine read_grid_hdf5_column_float(filename,ireadcol,nx,ny,nz,datarr,ierr) bind(c)
   use, intrinsic :: iso_c_binding
   character(kind=c_char), intent(in)  :: filename
   integer(kind=c_int),    intent(in)  :: ireadcol,nx,ny,nz
   real(kind=c_float),     intent(out) :: datarr(nx*ny*nz)
   integer(kind=c_int),    intent(out) :: ierr
  end subroutine read_grid_hdf5_column_float
 end interface read_grid_hdf5_column

end module hdf5utils
