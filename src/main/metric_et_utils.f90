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
! An evolved metric may also have refined blocks (an adaptive mesh, see
! extern_nr): an octree of blocks of nb^3 nodes, each refining an octant of
! its parent by a factor 2. The metric at a point is then interpolated from
! the finest block containing it, which holds the metric on its (nb+2)^3
! nodes from one node below its first node to one beyond its last.
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
 ! mixed spatial derivatives (xy, xz, yz, xyz) for Hermite interpolation
 real, allocatable          :: metriccrossgrid(:,:,:,:,:,:)
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
 ! solve the full Einstein constraints for the initial particle data
 ! (else use the linearised longitudinal-gauge slice)
 logical,            public :: nr_solve_init = .true.
 ! deposit of T_munu from the particles to the grid: 0 = normalised SPH
 ! interpolation at the nodes, 1 = SPH kernel sampled at the nodes,
 ! 2 = zero-second-moment kernel sampled at the nodes (smoothing error
 ! O((kh)^4)), 3 = normalised SPH kernel integrated exactly over the cells
 integer,            public :: nr_deposit = 2
 ! minimum smoothing length for the deposit at the nodes, in grid spacings
 ! (below 1, particles with h below the grid spacing have their kernel
 ! integrated exactly over the cells, which is accurate but slow)
 real,               public :: nr_deposit_hmin = 1.
 ! pressureless matter: deposit the coordinate densities rho* and rho* u_i
 ! twice per advance, from which the NR code computes T_munu at each stage,
 ! instead of processing the particles at every stage
 logical,            public :: nr_dust_source = .true.
 ! timestep limit from the expansion rate of the evolved metric
 real,               public :: dtmax_evolve  = huge(1.)
 ! cosmology with an evolved metric: redshift at which the scale factor is 1
 ! (< 0 if not a cosmological run), and the proper time (age) and mean scale
 ! factor of the grid, written to the dump header (proper time < 0 if unset)
 real,               public :: redshift_init = -1.
 real,               public :: proper_time   = -1.
 real,               public :: scale_factor  = 1.

 ! refined blocks: level, first node x0, spacing, children in the octants
 ! (x,y,z) -> 1 + ix + 2 iy + 4 iz (-1 = none), and the metric on the nodes
 type, public :: metric_block
  integer :: level = 0
  integer :: child(8) = -1
  real    :: x0(3) = 0., dx(3) = 0.
  real, allocatable :: gcov(:,:,:,:,:),derivs(:,:,:,:,:,:),cross(:,:,:,:,:,:)
 end type metric_block
 type(metric_block), allocatable, public :: mblocks(:)
 integer, public :: nmblocks = 0, mblocksize = 0
 ! level-1 block at each level-1 block position (-1 = none)
 integer, allocatable :: map1(:,:,:)
 ! refinement of the evolved metric: particles need a grid spacing below
 ! their smoothing length divided by nr_refine_h
 real,               public :: nr_refine_h = 1.

 ! routine that advances an evolved metric on the grid to a given time
 abstract interface
  subroutine advance_metric_func(time)
   real, intent(in) :: time
  end subroutine advance_metric_func
 end interface
 procedure(advance_metric_func), pointer, public :: advance_metric => null()

 public :: allocate_grid,initialize_grid,print_metric_grid
 public :: write_tabulated_metric,read_tabulated_metric,set_grid_ready
 public :: interpolate_gcov,interpolate_gcov_derivs
 public :: allocate_blocks,set_block_ready,find_block

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
 if (allocated(metriccrossgrid)) deallocate(metriccrossgrid)
 allocate(gcovgrid(0:3,0:3,nx,ny,nz))
 allocate(gcongrid(0:3,0:3,nx,ny,nz))
 allocate(sqrtggrid(nx,ny,nz))

 !metric derivs are stored in the form
 ! mu comp, nu comp, deriv, gridx,gridy,gridz
 ! Note that this is only the spatial derivs of
 ! the metric and we will need an additional array
 ! for time derivs
 allocate(metricderivsgrid(0:3,0:3,3,nx,ny,nz))
 allocate(metriccrossgrid(0:3,0:3,4,nx,ny,nz))

end subroutine allocate_grid

!---------------------------------------------------------------
!+
!  mark the grid as ready for interpolation once gcovgrid and
!  metricderivsgrid hold the metric, computing the mixed
!  derivatives used by the Hermite interpolation
!+
!---------------------------------------------------------------
subroutine set_grid_ready()
 use table_utils, only:hermite_cross_derivs_3d

 call hermite_cross_derivs_3d(16,gridsize,metricderivsgrid,dxgrid,periodic_grid,metriccrossgrid)
 gridinit = .true.

end subroutine set_grid_ready

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

 call set_grid_ready()

 ! Close the file
 close(iunit)
 ierr = 0

end subroutine read_tabulated_metric

!---------------------------------------------------------------
!+
!  allocate nblocks refined blocks of nb^3 nodes (none if nblocks = 0)
!+
!---------------------------------------------------------------
subroutine allocate_blocks(nblocks,nb)
 integer, intent(in) :: nblocks,nb
 integer :: i,m

 if (allocated(mblocks)) deallocate(mblocks)
 if (allocated(map1)) deallocate(map1)
 nmblocks = nblocks
 mblocksize = nb
 if (nblocks <= 0) return
 m = nb + 2
 allocate(mblocks(nblocks))
 do i=1,nblocks
    allocate(mblocks(i)%gcov(0:3,0:3,m,m,m),mblocks(i)%derivs(0:3,0:3,3,m,m,m))
    allocate(mblocks(i)%cross(0:3,0:3,4,m,m,m))
 enddo
 allocate(map1(0:2*gridsize(1)/nb-1,0:2*gridsize(2)/nb-1,0:2*gridsize(3)/nb-1))
 map1 = -1

end subroutine allocate_blocks

!---------------------------------------------------------------
!+
!  mark block id ready for interpolation once its level, position,
!  children, metric and metric derivatives are set: compute the mixed
!  derivatives and enter a level-1 block into the lookup table
!+
!---------------------------------------------------------------
subroutine set_block_ready(id)
 use table_utils, only:hermite_cross_derivs_3d
 integer, intent(in) :: id
 integer :: c(3),m(3)

 m = mblocksize + 2
 call hermite_cross_derivs_3d(16,m,mblocks(id)%derivs,mblocks(id)%dx,.false.,mblocks(id)%cross)
 if (mblocks(id)%level == 1) then
    c = nint((mblocks(id)%x0 - gridorigin)/(mblocksize*mblocks(id)%dx))
    map1(c(1),c(2),c(3)) = id
 endif

end subroutine set_block_ready

!-----------------------------------------------------------------------
!+
!  finest refined block containing a position (0 if none, i.e. the
!  position is only covered by the root grid)
!+
!-----------------------------------------------------------------------
pure integer function find_block(position) result(id)
 real, intent(in) :: position(3)
 real    :: rel(3),size1(3)
 integer :: c(3),cid,oct(3)

 id = 0
 if (nmblocks <= 0) return
 size1 = 0.5*mblocksize*dxgrid
 c = floor((position - gridorigin)/size1)
 if (periodic_grid) then
    c = modulo(c,2*gridsize/mblocksize)
 elseif (any(c < 0 .or. c >= 2*gridsize/mblocksize)) then
    return
 endif
 cid = map1(c(1),c(2),c(3))
 do while (cid > 0)
    id = cid
    rel = block_offset(position,id)
    oct = 0
    where (rel >= 0.5*mblocksize*mblocks(id)%dx) oct = 1
    cid = mblocks(id)%child(1 + oct(1) + 2*oct(2) + 4*oct(3))
 enddo

end function find_block

!-----------------------------------------------------------------------
!+
!  position relative to the first node of block id (wrapped into the
!  periodic box)
!+
!-----------------------------------------------------------------------
pure function block_offset(position,id) result(rel)
 real,    intent(in) :: position(3)
 integer, intent(in) :: id
 real :: rel(3)

 rel = position - mblocks(id)%x0
 if (periodic_grid) rel = modulo(rel,gridsize*dxgrid)

end function block_offset

!-----------------------------------------------------------------------
!+
!  tricubic Hermite interpolation of the metric (and optionally of its
!  derivatives) at a position, from the finest block containing it or
!  else from the root grid
!+
!-----------------------------------------------------------------------
pure subroutine interpolate_metric_at(position,gcov,dgcov)
 use table_utils, only:interp3d_hermite
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: dgcov(0:3,0:3,3)
 integer :: id,m(3)

 id = find_block(position)
 if (id > 0) then
    m = mblocksize + 2
    associate(b => mblocks(id))
       call interp3d_hermite(16,m,b%gcov,b%derivs,b%cross,b%x0 - b%dx,b%dx,.false.,&
                             b%x0 + block_offset(position,id),gcov,dgcov)
    end associate
 else
    call interp3d_hermite(16,gridsize,gcovgrid,metricderivsgrid,metriccrossgrid,gridorigin,dxgrid,&
                          periodic_grid,position,gcov,dgcov)
 endif

end subroutine interpolate_metric_at

!-----------------------------------------------------------------------
!+
!  interpolate the covariant metric to a given position, using tricubic
!  Hermite interpolation of the metric and its tabulated first derivatives
!+
!-----------------------------------------------------------------------
pure subroutine interpolate_gcov(position,gcov)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)

 call interpolate_metric_at(position,gcov)

end subroutine interpolate_gcov

!-----------------------------------------------------------------------
!+
!  spatial derivatives of the covariant metric, dgcov(:,:,i) = d g_munu / dx^i,
!  at a given position: the exact derivatives of the Hermite interpolant,
!  so that they are consistent with the metric from interpolate_gcov
!+
!-----------------------------------------------------------------------
pure subroutine interpolate_gcov_derivs(position,dgcov)
 real, intent(in)  :: position(3)
 real, intent(out) :: dgcov(0:3,0:3,3)
 real :: gcov(0:3,0:3)

 call interpolate_metric_at(position,gcov,dgcov)

end subroutine interpolate_gcov_derivs

end module metric_et_utils
