!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program tabulate_metric
!
! Tabulates an analytic metric (compiled with METRIC=...) on a uniform
! grid and writes it to a file that can be read with METRIC=et
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: tabulate_metric outfile nx ny nz xmin xmax ymin ymax zmin zmax [spin]
!
! :Dependencies: metric, metric_et_utils, metric_tools
!
 use metric,          only:a,mass1
 use metric_tools,    only:get_metric,get_metric_derivs
 use metric_et_utils, only:initialize_grid,write_tabulated_metric,gcovgrid,gcongrid,&
                           sqrtggrid,metricderivsgrid,gridorigin,dxgrid
 implicit none
 integer, parameter :: nargs_min = 10
 real,    parameter :: rmin_fac = 0.5  ! nodes with r < rmin_fac*mass1 are set to Minkowski
 character(len=120) :: outfile,string
 integer :: n(3),i,j,k,ierr,iarg
 real    :: xmin(3),xmax(3),lim(6),position(3),gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
 real    :: dgcovdx(0:3,0:3),dgcovdy(0:3,0:3),dgcovdz(0:3,0:3)

 if (command_argument_count() < nargs_min) then
    print "(a)",'Usage: tabulate_metric outfile nx ny nz xmin xmax ymin ymax zmin zmax [spin]'
    stop
 endif
 call get_command_argument(1,outfile)
 do iarg=1,3
    call get_command_argument(iarg+1,string)
    read(string,*) n(iarg)
 enddo
 do iarg=1,6
    call get_command_argument(iarg+4,string)
    read(string,*) lim(iarg)
 enddo
 xmin = lim(1:5:2)
 xmax = lim(2:6:2)
 if (command_argument_count() > nargs_min) then
    call get_command_argument(nargs_min+1,string)
    read(string,*) a
 endif
 print "(a,es10.3,a,3(i5),a,6(1x,f8.2))",' tabulating metric with spin a = ',a,' on grid ',n,' over ',lim

 ! Init grid and tabulated metric
 call initialize_grid(n,xmin,xmax)

 ! Fill the metric in the grid (triple loop over the grid)
 !$omp parallel do default(none) shared(n,gridorigin,dxgrid,mass1) &
 !$omp shared(gcovgrid,gcongrid,sqrtggrid,metricderivsgrid) &
 !$omp private(i,j,k,position,gcov,gcon,sqrtg,dgcovdx,dgcovdy,dgcovdz)
 do k=1,n(3)
    do j=1,n(2)
       do i=1,n(1)
          ! Calculate the current position in the grid
          position = gridorigin + (/i-1,j-1,k-1/)*dxgrid
          if (norm2(position) < rmin_fac*mass1) then
             ! avoid the coordinate singularity at the origin; never used by orbits
             gcov = 0.
             gcov(0,0) = -1.
             gcov(1,1) = 1.
             gcov(2,2) = 1.
             gcov(3,3) = 1.
             gcon = gcov
             sqrtg = 1.
             dgcovdx = 0.
             dgcovdy = 0.
             dgcovdz = 0.
          else
             call get_metric(position,gcov,gcon,sqrtg)
             call get_metric_derivs(position,dgcovdx,dgcovdy,dgcovdz)
          endif
          ! Store the calculated values in the grid arrays
          gcovgrid(:,:,i,j,k) = gcov
          gcongrid(:,:,i,j,k) = gcon
          sqrtggrid(i,j,k)    = sqrtg
          metricderivsgrid(:,:,1,i,j,k) = dgcovdx
          metricderivsgrid(:,:,2,i,j,k) = dgcovdy
          metricderivsgrid(:,:,3,i,j,k) = dgcovdz
       enddo
    enddo
 enddo
 !$omp end parallel do

 ! Write Data in file
 call write_tabulated_metric(outfile,ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR writing metric data to '//trim(outfile)
 else
    print "(a)",' metric data written to '//trim(outfile)
 endif

end program tabulate_metric
