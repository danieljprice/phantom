!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric_et_utils
!
! Utilities for handling tabulated metrics from the Einstein Toolkit
!
! :References: Magnall et al. (2023), Phys. Rev. D 108, 103534
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

 real, allocatable :: gcovgrid(:,:,:,:,:)
 real, allocatable :: gcongrid(:,:,:,:,:)
 real, allocatable :: sqrtggrid(:,:,:)
 real, allocatable :: metricderivsgrid(:,:,:,:,:,:)
 real              :: dxgrid(3), gridorigin(3)
 integer           :: gridsize(3)
 logical           :: gridinit = .false.

 ! Declaration of grid limits and dimensions
 integer, public :: nx,ny,nz
 real, parameter :: xmin = -10.0, xmax = 10.0
 real, parameter :: ymin = -10.0, ymax = 10.0
 real, parameter :: zmin = -10.0, zmax = 10.0
 real, parameter :: mass = 1.0  ! Mass of the central object

contains

!---------------------------------------------------------------
!+
!  allocate memory for the metric grid
!+
!---------------------------------------------------------------
subroutine allocate_grid(nxin,nyin,nzin,dx,dy,dz,originx,originy,originz)
 integer, intent(in) :: nxin,nyin,nzin
 real, intent(in) :: dx,dy,dz,originx,originy,originz

 nx = nxin
 ny = nyin
 nz = nzin
 gridsize(1) = nx
 gridsize(2) = ny
 gridsize(3) = nz

 dxgrid(1) = dx
 dxgrid(2) = dy
 dxgrid(3) = dz

 gridorigin(1) = originx
 gridorigin(2) = originy
 gridorigin(3) = originz

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
!  initialise a metric grid with a uniform grid
!  (currently size is hardwired but just for testing...)
!+
!---------------------------------------------------------------
subroutine initialize_grid()
 ! Local variable declarations
 real :: dx, dy, dz, x0(3)

 nx = 100
 ny = 100
 nz = 100

 ! Calculate the step size in each direction
 dx = (xmax - xmin) / (nx - 1)
 dy = (ymax - ymin) / (ny - 1)
 dz = (zmax - zmin) / (nz - 1)

 x0 = [0.,0.,0.]
 call allocate_grid(nx,ny,nz,dx,dy,dz,x0(1),x0(2),x0(3))

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
subroutine write_tabulated_metric(metric_file, ierr)
 character(len=*), intent(in) :: metric_file
 integer, intent(out) :: ierr
 integer :: iunit

 ! Open the file for writing
 open(newunit=iunit,file=metric_file,status='replace',form='unformatted',action='write',iostat=ierr)
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
subroutine read_tabulated_metric(metric_file, ierr)
 character(len=*), intent(in) :: metric_file
 integer, intent(out) :: ierr
 integer :: iunit

 ! Open the file for reading
 open(newunit=iunit,file=metric_file,status='old',form='unformatted',action='read',iostat=ierr)
 if (ierr /= 0) return

 ! Read the dimensions of the grid
 read(iunit) gridsize

 ! Read the grid origin and spacing
 read(iunit) gridorigin
 read(iunit) dxgrid

 nx = gridsize(1)
 ny = gridsize(2)
 nz = gridsize(3)

 call allocate_grid(nx,ny,nz,&
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

end module metric_et_utils
