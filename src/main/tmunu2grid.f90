module tmunu2grid
    implicit none 

contains
    subroutine get_tmunugrid_all(npart,xyzh,vxyzu,tmunus)
        use einsteintk_utils, only: dxgrid, gridorigin,boundsize,gridsize,gcovgrid,tmunugrid
        use interpolations3D, only: interpolate3D
        use boundary,         only: xmin,ymin,zmin,xmax,ymax,zmax
        use part, only: massoftype,igas,rhoh
        integer, intent(in) :: npart
        real, intent(in)    :: xyzh(:,:), vxyzu(:,:), tmunus(:,:,:)
        real                      :: weight,h,rho,pmass
        real                      :: xmininterp(3)
        integer                   :: ngrid(3)
        real,allocatable          :: datsmooth(:,:,:), dat(:)
        integer                   :: nnodes,i,k,j, ilower, iupper, jlower, jupper, klower, kupper
        logical                   :: normalise

        if (.not. allocated(datsmooth)) allocate (datsmooth(gridsize(1),gridsize(2),gridsize(3)))
        if (.not. allocated(dat)) allocate (dat(npart))
        ! All particles have equal weighting in the interp
        ! Here we calculate the weight for the first particle
        ! Get the smoothing length 
        h = xyzh(4,1)
        ! Get pmass
        pmass = massoftype(igas)
        ! Get density
        rho = rhoh(h,pmass)
        
        call get_weight(pmass,h,rho,weight)
        !print*, "Weighting for particle smoothing is: ", weight
        !weight = 1.
        ! For now we can set this to the origin, but it might need to be
        ! set to the grid origin of the CCTK_grid since we have boundary points
        ! TODO This should also be the proper phantom values and not a magic number
        !xmin(:) = gridorigin(:) - 0.5*dxgrid(:) ! We move the origin back by 0.5*dx to make a pseudo cell-centered grid
        xmininterp(1) =  xmin 
        xmininterp(2) =  ymin 
        xmininterp(3) =  zmin 
        
        !print*, "xmin: ", xmin
        !print*, "xmax: ", xmax
        call get_particle_domain(gridorigin(1),xmin,xmax,dxgrid(1),ilower,iupper)
        call get_particle_domain(gridorigin(2),ymin,ymax,dxgrid(2),jlower,jupper)
        call get_particle_domain(gridorigin(3),zmin,zmax,dxgrid(3),klower,kupper)
        !print*, "ivals: ", ilower, iupper
        ! nnodes is just the size of the mesh 
        ! might not be needed
        ! We note that this is not actually the size of the einstein toolkit grid
        ! As we want our periodic boundary to be on the particle domain not the 
        ! ET grid domain 
        ngrid(1) = (iupper-ilower)
        ngrid(2) = (jupper-jlower)
        ngrid(3) = (kupper-klower) 
        nnodes   = (iupper-ilower)*(jupper-jlower)*(kupper-klower)
        ! Do we want to normalise interpolations? 
        normalise = .true.

        
        !print*, "ngrid: ", ngrid
        
        !print*,"tmunu val:  ", tmunus(:,:,1)
        ! tt component

        tmunugrid = 0.
        do k=1,4
            do j=1,4
                do i=1, npart
                    dat(i) = tmunus(k,j,i)
                    !  if (dat(i) < 1.0 .and. i > 4) then 
                    !      print*, "dat: ", dat(i)
                    !      print*, "i is: ", i
                    !      stop 
                    !  endif
                enddo 
                !print*, "gcov: ", gcovgrid(:,:,1,1,1)
                !print*, "tmunugrid:  ", tmunugrid(:,:,1,1,1)
                ! print*, "k,j :", k, j
                ! print*, "Dat: ", dat(1:30)

                ! Get the position of the first grid cell x,y,z 
                ! print*, "x position of 1, 1, 1", gridorigin(:)
                ! print*, "x position of 1,1,1 calculated (cell centered)", xmin(1) + (1.-0.5)*dxgrid(1)
                ! Call to interpolate 3D 
                call interpolate3D(xyzh,weight,npart, &
                             xmininterp,tmunugrid(k-1,j-1,ilower:iupper,jlower:jupper,klower:kupper), &
                             nnodes,dxgrid,normalise,dat,ngrid)

                !print*, "Interpolated grid values are: ", datsmooth(4:38,4:38,4:38)
            enddo 
        enddo 
        ! do i=4,35
        !     do j=4,35
        !         do k=4,35
        !             if (tmunugrid(0,0,i,j,k) > 1.0008253314232896) then 
        !                 print*, "tmunugrid: ", tmunugrid(0,0,i,j,k) 
        !                 print*, "i,j,k: ", i,j,k
        !                 print*, "grid position i : ", gridorigin(1) + i*dxgrid(1)
        !                 print*, "grid position j : ", gridorigin(2) + j*dxgrid(2) 
        !                 print*, "grid position k : ", gridorigin(3) + k*dxgrid(3) 
                        
        !                 !stop 
        !             endif 
        !         enddo 
        !     enddo 
        ! enddo     
        !print*, "tmunugrid: ", tmunugrid(0,0,5,5,5:35) 
        !stop
    end subroutine get_tmunugrid_all

    subroutine get_weight(pmass,h,rhoi,weight)
        real, intent(in)  :: pmass,h,rhoi
        real, intent(out) :: weight 
        
        weight = (pmass*h**3.)/rhoi 

    end subroutine get_weight

    subroutine get_dat(tmunus,dat)
        real, intent(in)  :: tmunus
        real, intent(out) :: dat

    end subroutine get_dat

    subroutine get_particle_domain(gridorigin,xmin,xmax,dxgrid,ilower,iupper)
        real,    intent(in)  :: gridorigin, xmin,xmax, dxgrid
        integer, intent(out) :: ilower, iupper


        ilower = int((xmin - gridorigin)/dxgrid) + 1 ! +1 since our arrays start at 1 not 0 
        iupper = int((xmax - gridorigin)/dxgrid) + 1

    end subroutine get_particle_domain

end module tmunu2grid