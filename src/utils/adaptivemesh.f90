!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module adaptivemesh
!
! module to construct/hold adaptive mesh structure
! used for interpolating particles to adaptive grid
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim
!
 use dim, only:periodic
 implicit none
 !--this controls the number of cells on each level
 integer, parameter :: nsub = 2
 integer, parameter :: ndim = 3
 integer, parameter :: maxlevels = 20
 integer, parameter :: minlevels = 7
 real, parameter    :: overrefinefac = 1.25
 !
 !--memory allowed for tree:
 !  bear in mind that total number of cells is maxmeshes*nsub**ndim
 !  so use maxp/nsub**ndim to get number of cells = maxp
 !
 integer :: maxmeshes = 0 ! zero until allocated
 !
 !--resolution of the root grid (2^ifirstlevel)^ndim
 !
 integer, parameter :: ifirstlevel = 1

 integer, parameter :: maxchildren = nsub**ndim
 !
 !--grid array stores children in first 8 cells, then parent
 !  (currently we do not store the parent as it is not needed)
 !
 integer, allocatable :: gridnodes(:,:)
! integer, parameter :: ilocparent = maxchildren + 1

 integer :: maxlevel

contains

subroutine build_mesh(xyzh,np,nmesh,xmin,dxmax)
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(in)  :: np
 integer, intent(out) :: nmesh
 real,    intent(in), optional :: xmin(ndim),dxmax(ndim)
 real :: xminp(ndim),dxmaxp(ndim)
 integer :: i,ierr,nattempts

 if (present(xmin)) then
    if (.not.present(dxmax)) stop 'error in call to build_mesh: dxmax not present'
    !
    !--if boundaries are sent in, use these as the root node size
    !
    xminp(:)  = xmin(:)
    dxmaxp(:) = dxmax(:)
 else
    !
    !--otherwise, set the root node size to fit all the particles
    !
    xminp(:)  = huge(0.)
    dxmaxp(:) = -huge(0.)
    do i=1,np
       xminp(:)  = min(xyzh(1:ndim,i),xminp(:))
       dxmaxp(:) = max(xyzh(1:ndim,i),dxmaxp(:))
    enddo
    dxmaxp(:) = dxmaxp(:) - xminp(:) + 1.e-6
 endif
 print*,'build_mesh: xmin = ',xminp(:)
 print*,'            xmax = ',xminp(:) + dxmaxp(:)
 !
 !--build the mesh
 !
 nattempts = 0
 ierr = 1
 maxmeshes = max(int(overrefinefac*(2*np + 1)),2**minlevels) ! allow twice the number of particles

 do while(ierr /= 0 .and. nattempts < 10)
    nattempts = nattempts + 1
    !
    !--initialise the tree structure
    !  (root node only)
    !
    maxlevel = ifirstlevel
    nmesh   = 1
    !
    !--allocate memory
    !
    print*,' allocating mesh with max=',maxmeshes
    allocate(gridnodes(maxchildren,maxmeshes))
    gridnodes(:,1) = -1
    ierr = 0

    !$omp parallel do default(none) schedule(guided,10) &
    !$omp shared(np,xyzh,xminp,dxmaxp,nmesh) &
    !$omp private(i) &
    !$omp reduction(+:ierr)
    do i=1,np
       call refine_mesh(xyzh(:,i),1,ifirstlevel,xminp,dxmaxp,nmesh,ierr)
    enddo
    !$omp end parallel do
    if (ierr /= 0) then
       maxmeshes = maxmeshes*2
       call free_mesh
    endif
 enddo
 if (ierr /= 0) stop 'could not allocate enough memory for mesh (after 10 attempts)'

 print "(a,i10,a,i10,a,i2,a,i6,a)",&
       ' build_mesh: nmeshes = ',nmesh,', ncells = ',nmesh*(nsub**ndim),', max level = ',maxlevel, &
       ' (effective resolution = ',nsub**maxlevel,'^3)'

end subroutine build_mesh

recursive subroutine refine_mesh(xyzhi,imesh,level,xminl,dxmax,nmesh,ierr)
 real,    intent(in)    :: xyzhi(4)
 integer, intent(in)    :: imesh
 integer, intent(in)    :: level
 real,    intent(in)    :: xminl(ndim)
 real,    intent(in)    :: dxmax(ndim)
 integer, intent(inout) :: nmesh
 integer, intent(out)   :: ierr
 real :: dxcell(ndim),xminnew(ndim)
 real                     :: hmincell,radkern,dlevel
 integer                  :: icell,isubmesh,ipix,jpix,kpix,isublevel
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 !character(len=maxlevels), parameter :: string = '1234567890'

 ierr = 0
 if (level > maxlevels) then
    print "(2(a,i2),a)",' INTERNAL ERROR: level ',level,' > maxlevels (',maxlevels,') in refine_mesh'
    stop
 endif
 dlevel    = 1./real(nsub**level)
 dxcell(:) = dxmax(:)*dlevel
 !print*,string(1:level)//'level ',level,' effective grid size = ',nint(dxmax(:)/dxcell(:))

 hmincell = overrefinefac*0.5*minval(dxcell(:))

 !--work out which daughter cell the particle contributes to
 radkern = 2.*xyzhi(4)
 ipixmin = int((xyzhi(1) - radkern - xminl(1))/dxcell(1)) + 1
 jpixmin = int((xyzhi(2) - radkern - xminl(2))/dxcell(2)) + 1
 kpixmin = int((xyzhi(3) - radkern - xminl(3))/dxcell(3)) + 1

 ipixmax = int((xyzhi(1) + radkern - xminl(1))/dxcell(1)) + 1
 jpixmax = int((xyzhi(2) + radkern - xminl(2))/dxcell(2)) + 1
 kpixmax = int((xyzhi(3) + radkern - xminl(3))/dxcell(3)) + 1
 !print*,string(1:level)//'x range = ',xminl(1),xminl(1) + (nsub)*dxcell(1),ipixmin,ipixmax
 !print*,string(1:level)//'y range = ',xminl(2),xminl(2) + (nsub)*dxcell(2),jpixmin,jpixmax
 !print*,string(1:level)//'z range = ',xminl(3),xminl(3) + (nsub)*dxcell(3),kpixmin,kpixmax

 if (ipixmin < 1)    ipixmin = 1  ! make sure they only contribute
 if (jpixmin < 1)    jpixmin = 1  ! to pixels in the image
 if (kpixmin < 1)    kpixmin = 1
 if (ipixmax > nsub) ipixmax = nsub
 if (jpixmax > nsub) jpixmax = nsub
 if (kpixmax > nsub) kpixmax = nsub

 if (periodic) then
    !--these errors should never happen with periodic bc's
    if (ipixmax < 1) stop 'INTERNAL ERROR: ipixmax < 1'
    if (jpixmax < 1) stop 'INTERNAL ERROR: jpixmax < 1'
    if (kpixmax < 1) stop 'INTERNAL ERROR: kpixmax < 1'
    if (ipixmin > nsub) stop 'INTERNAL ERROR: ipixmin > nsub'
    if (jpixmin > nsub) stop 'INTERNAL ERROR: jpixmin > nsub'
    if (kpixmin > nsub) stop 'INTERNAL ERROR: kpixmin > nsub'
 endif
 ipix = int((xyzhi(1) - xminl(1))/dxcell(1)) + 1
 jpix = int((xyzhi(2) - xminl(2))/dxcell(2)) + 1
 kpix = int((xyzhi(3) - xminl(3))/dxcell(3)) + 1

 do kpix=kpixmin,kpixmax
    do jpix=jpixmin,jpixmax
       do ipix=ipixmin,ipixmax

          icell = ((kpix-1)*nsub + (jpix-1))*nsub + ipix

          !print*,string(1:level)//'%node ',imesh,' daughter cell = ',icell,' (',ipix,jpix,kpix,' )'
          !if (icell <= 0 .or. icell > nsub**ndim) stop 'refine_mesh: error in cell id'

          isubmesh = gridnodes(icell,imesh) !grid(imesh)%daughter(icell)
          if (isubmesh > 0) then
             !
             !--if cell is already refined, descend one level
             !  and repeat the process
             !
             xminnew(1)  = xminl(1) + (ipix-1)*dxcell(1)
             xminnew(2)  = xminl(2) + (jpix-1)*dxcell(2)
             xminnew(3)  = xminl(3) + (kpix-1)*dxcell(3)
             !print*,string(1:level)//'%already refined' !' new xmin = ',xminl(:)

             isublevel = level + 1
             if (isublevel > maxlevels) then ! should never happen
                print*,'INTERNAL ERROR: sublevel = ',isublevel, &
                       ' mesh ',imesh,' grid = ',gridnodes(:,imesh) !grid(imesh)%daughter(:)
             endif
             call refine_mesh(xyzhi,isubmesh,isublevel,xminnew,dxmax,nmesh,ierr)

             !
             !--this line specifies the refinement criterion
             !
          elseif (((xyzhi(4) < hmincell .and. xyzhi(4) > tiny(0.)) &
                   .or.(level < minlevels)) &
                  .and. level < maxlevels) then
             !
             !--if node is not already refined and the refinement
             !  criterion is satisfied (here if h < 0.5*dx)
             !

             !$omp critical
             nmesh = nmesh + 1
             if (nmesh > maxmeshes) then
                !print*,'ERROR: nmesh > maxmeshes (',maxmeshes,'): change parameter and recompile'
                ierr = 1
             else
                !print*,string(1:level)//'%refining: adding node ',nmesh,' on level ',level+1
                gridnodes(icell,imesh) = nmesh
                gridnodes(:,nmesh)     = -1
                !
                !--could stop here if the refinement criterion was
                !  based on the number of particles in the cell
                !  However, with hmin we should proceed to the next
                !  level to see if this is refined or not
                !
                xminnew(1)  = xminl(1) + (ipix-1)*dxcell(1)
                xminnew(2)  = xminl(2) + (jpix-1)*dxcell(2)
                xminnew(3)  = xminl(3) + (kpix-1)*dxcell(3)

                isublevel = level + 1
                isubmesh  = nmesh      ! be careful to pass by value here, not by reference...
             endif
             !$omp end critical

             if (ierr /= 0) return

             call refine_mesh(xyzhi,isubmesh,isublevel,xminnew,dxmax,nmesh,ierr)
             maxlevel = max(level+1,maxlevel)
          endif

       enddo
    enddo
 enddo

end subroutine refine_mesh

subroutine free_mesh()

 if (allocated(gridnodes)) deallocate(gridnodes)

end subroutine free_mesh

end module adaptivemesh
