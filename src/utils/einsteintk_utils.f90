module einsteintk_utils
      implicit none 
      real, allocatable :: gcovgrid(:,:,:,:,:)
      real, allocatable :: gcongrid(:,:,:,:,:)
      real, allocatable :: sqrtggrid(:,:,:)
      real, allocatable :: tmunugrid(:,:,:,:,:)
      real, allocatable :: metricderivsgrid(:,:,:,:,:,:)
      real              :: dxgrid(3), gridorigin(3), boundsize(3)
      integer           :: gridsize(3)
      logical           :: gridinit = .false.
      character(len=128)  :: logfilestor,evfilestor,dumpfilestor,infilestor
contains
      subroutine init_etgrid(nx,ny,nz,dx,dy,dz,originx,originy,originz)
      integer, intent(in) :: nx,ny,nz
      real,    intent(in) :: dx,dy,dz,originx,originy,originz
      !integer, intent(in) :: boundsizex, boundsizey, boundsizez
      
      gridsize(1) = nx
      gridsize(2) = ny
      gridsize(3) = nz 

      dxgrid(1) = dx
      dxgrid(2) = dy
      dxgrid(3) = dz
      
      gridorigin(1) = originx
      gridorigin(2) = originy
      gridorigin(3) = originz

      ! How mmany grid points is the boundary?
      ! boundsize(1) = boundsizex
      ! boundsize(2) = boundsizey
      ! boundsize(3) = boundsizez
      
      allocate(gcovgrid(0:3,0:3,nx,ny,nz))
      allocate(gcongrid(0:3,0:3,nx,ny,nz))
      allocate(sqrtggrid(nx,ny,nz))

      ! Will need to delete this at somepoint 
      ! For now it is the simplest way 
      allocate(tmunugrid(0:3,0:3,nx,ny,nz))

      ! metric derivs are stored in the form 
      ! mu comp, nu comp, deriv, gridx,gridy,gridz 
      ! Note that this is only the spatial derivs of 
      ! the metric and we will need an additional array 
      ! for time derivs
      allocate(metricderivsgrid(0:3,0:3,3,nx,ny,nz)) 
      
      gridinit = .true.

    end subroutine init_etgrid
    
    subroutine print_etgrid()
      ! Subroutine for printing quantities of the ET grid 

      print*, "Grid spacing (x,y,z) is : ", dxgrid
      print*, "Grid origin (x,y,z) is: ", gridorigin
      !print*, "Grid size is: ", sizeof(gcovgrid)
      print*, "Covariant metric tensor of the grid is: ", gcovgrid(:,:,1,1,1)
      !print*, "Contravariant metric tensor of the grid is: ", gcongrid
      !print*, "Negative sqrtg of the grid is: ", sqrtggrid

    end subroutine print_etgrid
end module einsteintk_utils
