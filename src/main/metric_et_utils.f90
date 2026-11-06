!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric_et_utils
!
! Utilities for handling tabulated metrics from the Einstein Toolkit
! (or from a numerical relativity code), and for interpolating the metric
! and its spatial derivatives from the grid to arbitrary positions
!
! The grid is vertex-centred: node (i,j,k) is at gridorigin + (i-1,j-1,k-1)*dxgrid.
! For a periodic grid, node gridsize+1 is identified with node 1.
!
! :References: Magnall et al. (2023), Phys. Rev. D 108, 103534
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: table_utils
!
 implicit none

 real, allocatable, public :: gcovgrid(:,:,:,:,:)
 real, allocatable, public :: gcongrid(:,:,:,:,:)
 real, allocatable, public :: sqrtggrid(:,:,:)
 real, allocatable, public :: metricderivsgrid(:,:,:,:,:,:)
 real,              public :: dxgrid(3), gridorigin(3)
 integer,           public :: gridsize(3)
 logical,           public :: gridinit = .false.

 ! Declaration of grid limits and dimensions
 integer,           public :: nx,ny,nz

 logical,            public :: periodic_grid = .false.
 character(len=120), public :: metric_file   = 'tabulated_metric.dat'

 ! metric evolved by the numerical relativity code (see extern_nr),
 ! on a periodic grid of ngrid_evolve nodes per dimension
 logical,            public :: evolve_metric = .false.
 integer,            public :: ngrid_evolve  = 32
 ! timestep limit from the expansion rate of the evolved metric
 real,               public :: dtmax_evolve  = huge(1.)

 ! routine that advances an evolved metric on the grid to a given time
 abstract interface
  subroutine advance_metric_func(time)
   real, intent(in) :: time
  end subroutine advance_metric_func
 end interface
 procedure(advance_metric_func), pointer, public :: advance_metric => null()

 public :: allocate_grid,initialize_grid,print_metric_grid
 public :: write_tabulated_metric,read_tabulated_metric
 public :: interpolate_gcov,interpolate_gcov_derivs

 private

contains

!---------------------------------------------------------------
!+
!  allocate memory for the metric grid
!+
!---------------------------------------------------------------
subroutine allocate_grid(nxin,nyin,nzin,dx,dy,dz,originx,originy,originz)
 integer, intent(in) :: nxin,nyin,nzin
 real,    intent(in) :: dx,dy,dz,originx,originy,originz

 nx = nxin
 ny = nyin
 nz = nzin
 gridsize = (/nx,ny,nz/)
 dxgrid   = (/dx,dy,dz/)
 gridorigin = (/originx,originy,originz/)

 if (allocated(gcovgrid)) deallocate(gcovgrid)
 if (allocated(gcongrid)) deallocate(gcongrid)
 if (allocated(sqrtggrid)) deallocate(sqrtggrid)
 if (allocated(metricderivsgrid)) deallocate(metricderivsgrid)
 allocate(gcovgrid(0:3,0:3,nx,ny,nz))
 allocate(gcongrid(0:3,0:3,nx,ny,nz))
 allocate(sqrtggrid(nx,ny,nz))

 !metric derivs are stored in the form
 ! mu comp, nu comp, deriv, gridx,gridy,gridz
 ! Note that this is only the spatial derivs of
 ! the metric and we will need an additional array
 ! for time derivs
 allocate(metricderivsgrid(0:3,0:3,3,nx,ny,nz))

end subroutine allocate_grid

!---------------------------------------------------------------
!+
!  initialise a uniform, non-periodic grid with n(1:3) nodes
!  spanning xmin(1:3) to xmax(1:3) inclusive
!+
!---------------------------------------------------------------
subroutine initialize_grid(n,xmin,xmax)
 integer, intent(in) :: n(3)
 real,    intent(in) :: xmin(3),xmax(3)
 real :: dx(3)

 ! Calculate the step size in each direction
 dx = (xmax - xmin)/(n - 1)
 call allocate_grid(n(1),n(2),n(3),dx(1),dx(2),dx(3),xmin(1),xmin(2),xmin(3))
 periodic_grid = .false.
 gridinit = .true.

end subroutine initialize_grid

!---------------------------------------------------------------
!+
!  print information about the metric grid
!+
!---------------------------------------------------------------
subroutine print_metric_grid()

 print*, "Grid spacing (x,y,z) is : ", dxgrid
 print*, "Grid origin (x,y,z) is: ", gridorigin
 print*, "Covariant metric tensor of the grid is: ", gcovgrid(:,:,1,1,1)

end subroutine print_metric_grid

!---------------------------------------------------------------
!+
!  write tabulated metric to file
!+
!---------------------------------------------------------------
subroutine write_tabulated_metric(filename, ierr)
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer :: iunit

 ! Open the file for writing
 open(newunit=iunit,file=filename,status='replace',form='unformatted',action='write',iostat=ierr)
 if (ierr /= 0) then
    ierr = 1
    return
 endif

 ! Write the dimensions of the grid
 write(iunit) gridsize

 ! Write the grid origin and spacing
 write(iunit) gridorigin
 write(iunit) dxgrid

 ! Write the metric values to the file
 write(iunit) gcovgrid
 write(iunit) gcongrid
 write(iunit) sqrtggrid
 write(iunit) metricderivsgrid

 ! Close the file
 close(iunit)
 ierr = 0

end subroutine write_tabulated_metric

!---------------------------------------------------------------
!+
!  read tabulated metric from file
!+
!---------------------------------------------------------------
subroutine read_tabulated_metric(filename, ierr)
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer :: iunit

 ! Open the file for reading
 open(newunit=iunit,file=filename,status='old',form='unformatted',action='read',iostat=ierr)
 if (ierr /= 0) return

 ! Read the dimensions of the grid
 read(iunit) gridsize

 ! Read the grid origin and spacing
 read(iunit) gridorigin
 read(iunit) dxgrid

 call allocate_grid(gridsize(1),gridsize(2),gridsize(3),&
        dxgrid(1),dxgrid(2),dxgrid(3),&
        gridorigin(1),gridorigin(2),gridorigin(3))

 ! Read the metric values from the file
 read(iunit) gcovgrid
 read(iunit) gcongrid
 read(iunit) sqrtggrid
 read(iunit) metricderivsgrid

 gridinit = .true.

 ! Close the file
 close(iunit)
 ierr = 0

end subroutine read_tabulated_metric

!-----------------------------------------------------------------------
!+
!  interpolate the covariant metric to a given position, using tricubic
!  Hermite interpolation of the metric and its tabulated first derivatives
!+
!-----------------------------------------------------------------------
pure subroutine interpolate_gcov(position,gcov)
 use table_utils, only:interp3d_hermite
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)

 call interp3d_hermite(16,gridsize,gcovgrid,metricderivsgrid,gridorigin,dxgrid,&
                       periodic_grid,position,gcov)

end subroutine interpolate_gcov

!-----------------------------------------------------------------------
!+
!  spatial derivatives of the covariant metric, dgcov(:,:,i) = d g_munu / dx^i,
!  at a given position: the exact derivatives of the Hermite interpolant,
!  so that they are consistent with the metric from interpolate_gcov
!+
!-----------------------------------------------------------------------
pure subroutine interpolate_gcov_derivs(position,dgcov)
 use table_utils, only:interp3d_hermite
 real, intent(in)  :: position(3)
 real, intent(out) :: dgcov(0:3,0:3,3)
 real :: gcov(0:3,0:3)

 call interp3d_hermite(16,gridsize,gcovgrid,metricderivsgrid,gridorigin,dxgrid,&
                       periodic_grid,position,gcov,dgcov)

end subroutine interpolate_gcov_derivs

end module metric_et_utils
