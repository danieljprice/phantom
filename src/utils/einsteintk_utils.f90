module einsteintk_utils
      implicit none 
      real, allocatable :: gcovgrid(:,:,:,:,:)
      real, allocatable :: gcongrid(:,:,:,:,:)
      real, allocatable :: sqrtggrid(:,:,:)
      real, allocatable :: tmunugrid(:,:,:,:,:)
      real, allocatable :: rhostargrid(:,:,:)
      real, allocatable :: pxgrid(:,:,:,:)
      real, allocatable :: metricderivsgrid(:,:,:,:,:,:)
      real              :: dxgrid(3), gridorigin(3), boundsize(3)
      integer           :: gridsize(3)
      logical           :: gridinit = .false.
      logical           :: exact_rendering
      character(len=128)  :: logfilestor,evfilestor,dumpfilestor,infilestor
contains
      subroutine init_etgrid(nx,ny,nz,dx,dy,dz,originx,originy,originz)
      integer, intent(in) :: nx,ny,nz
      real,    intent(in) :: dx,dy,dz,originx,originy,originz
      
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

      ! Will need to delete this at somepoint 
      ! For now it is the simplest way 
      allocate(tmunugrid(0:3,0:3,nx,ny,nz))

      allocate(pxgrid(3,nx,ny,nz))

      allocate(rhostargrid(nx,ny,nz))

      ! metric derivs are stored in the form 
      ! mu comp, nu comp, deriv, gridx,gridy,gridz 
      ! Note that this is only the spatial derivs of 
      ! the metric and we will need an additional array 
      ! for time derivs
      allocate(metricderivsgrid(0:3,0:3,3,nx,ny,nz)) 
      
      gridinit = .true.
      !exact_rendering = exact

    end subroutine init_etgrid
    
    subroutine print_etgrid()
      ! Subroutine for printing quantities of the ET grid 

      print*, "Grid spacing (x,y,z) is : ", dxgrid
      print*, "Grid origin (x,y,z) is: ", gridorigin
      print*, "Covariant metric tensor of the grid is: ", gcovgrid(:,:,1,1,1)

    end subroutine print_etgrid

    subroutine get_particle_rhs(i,vx,vy,vz,fx,fy,fz,e_rhs)
      use part,   only: vxyzu,fxyzu,fext
      integer, intent(in) :: i
      real, intent(out) :: vx,vy,vz,fx,fy,fz,e_rhs 

      !vxyz
      vx = vxyzu(1,i)
      vy = vxyzu(2,i) 
      vz = vxyzu(3,i)
      
      ! dp/dt 
      !print*, "fext: ", fext(:,i)
      !print*, "fxyzu: ", fxyzu(:,i)
      !fx = fxyzu(1,i) + fext(1,i)
      !print*, "fx: ", fx 
      !fy = fxyzu(2,i) + fext(2,i)
      !fz = fxyzu(3,i) + fext(3,i)
      fx = fext(1,i)
      fy = fext(2,i)
      fz = fext(3,i)


      ! de/dt
      e_rhs = 0. 

    end subroutine get_particle_rhs

    subroutine get_particle_val(i,x,y,z,px,py,pz,e)
      use part,   only: xyzh, pxyzu
      integer, intent(in) :: i
      real, intent(out) :: x,y,z,px,py,pz,e 

      !xyz
      x = xyzh(1,i) 
      y = xyzh(2,i) 
      z = xyzh(3,i)
      
      ! p
      px = pxyzu(1,i)
      py = pxyzu(2,i)
      pz = pxyzu(3,i)

      ! e
      ! ???
      e = pxyzu(4,i) 

    end subroutine get_particle_val

    subroutine set_particle_val(i,x,y,z,px,py,pz,e)
      use part, only: xyzh, pxyzu
      integer, intent(in) :: i
      real, intent(in) :: x,y,z,px,py,pz,e
      ! Subroutine for setting the particle values in phantom 
      ! using the values stored in einstein toolkit before a dump 
      
      !xyz
      xyzh(1,i) = x 
      xyzh(2,i) = y 
      xyzh(3,i) = z 

      ! p 
      pxyzu(1,i) = px 
      pxyzu(2,i) = py 
      pxyzu(3,i) = pz 
      pxyzu(4,i) = e 


    end subroutine set_particle_val

    subroutine get_phantom_dt(dtout)
      use part, only:xyzh
      real, intent(out) :: dtout
      real, parameter :: safety_fac = 0.2
      real :: minh 

      ! Get the smallest smoothing length
      minh = minval(xyzh(4,:))

      ! Courant esque condition from Rosswog 2021+
      ! Since c is allways one in our units 
      dtout = safety_fac*minh
      print*, "dtout phantom: ", dtout


    end subroutine get_phantom_dt

    subroutine set_rendering(flag)
      logical, intent(in) :: flag 

      exact_rendering = flag

    end subroutine set_rendering
    
    ! Do I move this to tmunu2grid?? 
    ! I think yes 
    

    ! Moved to einsteintk_wrapper.f90 to fix dependency issues 
    
    ! subroutine get_metricderivs_all(dtextforce_min)
    !   use part, only:npart, xyzh,vxyzu,metrics,metricderivs,dens,fext
    !   use timestep, only:bignumber,C_force
    !   use extern_gr, only:get_grforce
    !   use metric_tools, only:pack_metricderivs
    !   real, intent(out) :: dtextforce_min
    !   integer :: i
    !   real :: pri,dtf

    !   pri = 0.
    !   dtextforce_min = bignumber

    !   !$omp parallel do default(none) &
    !   !$omp shared(npart, xyzh,metrics,metricderivs,vxyzu,dens,C_force,fext) &
    !   !$omp firstprivate(pri) & 
    !   !$omp private(i,dtf) & 
    !   !$omp reduction(min:dtextforce_min)
    !   do i=1, npart 
    !     call pack_metricderivs(xyzh(1:3,i),metricderivs(:,:,:,i))
    !     call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i), &
    !          vxyzu(1:3,i),dens(i),vxyzu(4,i),pri,fext(1:3,i),dtf)
    !     dtextforce_min = min(dtextforce_min,C_force*dtf)
    !   enddo 
    !   !$omp end parallel do 
    ! end subroutine get_metricderivs_all
end module einsteintk_utils
