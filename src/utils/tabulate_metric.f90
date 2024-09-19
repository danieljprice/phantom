!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program tabulate_metric
!
! tabulate_metric
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: tabulate_metric [no arguments]
!
! :Dependencies: metric, metric_et_utils
!
 use metric_et_utils
 !use metric

 implicit none

 integer :: ierr
 character(len=64) :: metric_file = 'tabuled_metric.dat'


 ! Init grid and tabulated metric
 call initialize_grid()

 ! Fill and interpolate metric in the grid
 call fill_grid()

 ! Write Data in file
 call write_tabulated_metric(metric_file, ierr)

 if (ierr /= 0) then
    print *, 'Error writing metric data to file'
 else
    print *, 'Metric data successfully written to file'
 endif

contains

subroutine fill_grid()
 use metric
 integer :: i, j, k
 real :: dx, dy, dz
 real :: position(3)
 real :: gcov(0:3,0:3)
 real :: gcon(0:3,0:3)
 real :: sqrtg
 real :: dgcovdx(0:3,0:3)
 real :: dgcovdy(0:3,0:3)
 real :: dgcovdz(0:3,0:3)
 ! Triple loop to fill the grid
 dx = (xmax - xmin) / (nx - 1)
 dy = (ymax - ymin) / (ny - 1)
 dz = (zmax - zmin) / (nz - 1)

 do i = 1, nx
    do j = 1, ny
       do k = 1, nz
          ! Calculate the current position in the grid
          position(1) = xmin + (i - 1) * dx
          position(2) = ymin + (j - 1) * dy
          position(3) = zmin + (k - 1) * dz
          ! Store the calculated values in the grid arrays
          call get_metric_cartesian(position,gcov,gcon,sqrtg)
          !call get_metric_derivs(position,dgcovdx, dgcovdy, dgcovdz)
          call metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
          gcovgrid(:,:,i,j,k) = gcov
          gcongrid(:,:,i,j,k) = gcon
          sqrtggrid(i,j,k) = sqrtg
          metricderivsgrid(:,:,1,i,j,k) = dgcovdx
          metricderivsgrid(:,:,2,i,j,k) = dgcovdy
          metricderivsgrid(:,:,3,i,j,k) = dgcovdz
       enddo
    enddo
 enddo
end subroutine fill_grid

end program tabulate_metric

