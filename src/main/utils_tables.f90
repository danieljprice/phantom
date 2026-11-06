!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module table_utils
!
! This module contains low-level utilities related to tabulated data
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

 public :: yinterp, linspace, logspace, diff, flip_array, interpolator
 public :: find_nearest_index, interp_1d, linear_interpolator_one_d, interpolate_1d
 public :: differentiate, interp3d_hermite

 private

contains

!----------------------------------------------------
!+
!  Function to linearly interpolate values from
!  from an array. Given x it finds the position
!  in the table of x values (xtab) and interpolates
!  the value of y from the table of y values (ytab)
!  to return the y value corresponding to the
!  position x
!+
!----------------------------------------------------
pure real function yinterp(ytab,xtab,x)
 real, intent(in) :: ytab(:),xtab(:)
 real, intent(in) :: x
 integer :: i,n
 real    :: xprev,xi,dx,dydx

 n = size(xtab)

 yinterp = 0.
 if (n <= 0) return

 if (x <= xtab(1)) then
    yinterp = ytab(1)
    return
 elseif (x >= xtab(n)) then
    yinterp = ytab(n)
    return
 endif

 xprev = xtab(1)
 over_tab: do i=2,n
    xi = xtab(i)
    if (x >= xprev .and. x < xi) then
       dx      = x - xprev
       dydx    = (ytab(i) - ytab(i-1))/(xi - xprev)
       yinterp = ytab(i-1) + dx*dydx
       exit over_tab
    endif
    xprev = xi
 enddo over_tab

end function yinterp

!--------------------------------------------------------
!+
!  Function to fill an array with equally spaced points
!  e.g. call linspace(rgrid,rmin,rmax)
!
!  the argument dx is optional giving the grid spacing
!+
!--------------------------------------------------------
pure subroutine linspace(x,xmin,xmax,dx)
 real, intent(out) :: x(:)
 real, intent(in)  :: xmin,xmax
 real, intent(out), optional :: dx
 integer :: i, n
 real    :: dxi

 n  = size(x)
 dxi = (xmax - xmin)/real(n-1)
 do i=1,n
    x(i) = xmin + (i-1)*dxi
 enddo
 if (present(dx)) dx = dxi

end subroutine linspace

!--------------------------------------------------------
!+
!  Function to fill an array with equally log-spaced points
!+
!--------------------------------------------------------
pure subroutine logspace(x,xmin,xmax)
 real, intent(out) :: x(:)
 real, intent(in)  :: xmin,xmax
 integer :: i, n
 real    :: dx

 n  = size(x)

 dx = log10(xmax/xmin)/real(n-1)
 do i=1,n
    x(i) = log10(xmin) + (i-1)*dx
 enddo

 x = 10.**x

end subroutine logspace

!----------------------------------------------------------------
!+
!  Finds index of the array value closest to a given value for an
!  ordered array
!+
!----------------------------------------------------------------
subroutine interpolator(array,value,valueidx)
 real,    intent(in)  :: array(:)
 real,    intent(in)  :: value
 integer, intent(out) :: valueidx

 valueidx = minloc(abs(array - value), dim = 1)

end subroutine interpolator

!----------------------------------------------------------------
!+
!  Reverses the elements of a 1-d array
!+
!----------------------------------------------------------------
subroutine flip_array(array)
 real, intent(inout) :: array(:)
 real, allocatable   :: flipped_array(:)
 integer             :: i

 allocate(flipped_array(size(array)))
 do i = 1, size(array)
    flipped_array(i) = array(size(array) - i + 1)
 enddo
 array = flipped_array

end subroutine flip_array

!----------------------------------------------------------------
!+
!  Takes a n-dim array and produces a (n-1)-dim array with the
!  ith element being the (i+1)th element minus the ith element of
!  the original array
!+
!----------------------------------------------------------------
subroutine diff(array, darray)
 real,              intent(in)  :: array(:)
 real, allocatable, intent(out) :: darray(:)
 integer                        :: i

 allocate(darray(size(array)-1))
 do i = 1, size(array)-1
    darray(i) = array(i+1) - array(i)
 enddo

end subroutine diff

!-----------------------------------------------------------------------
!+
!  Find index of nearest lower value in array
!+
!-----------------------------------------------------------------------
subroutine find_nearest_index(arr,val,indx)
 real,    intent(in)  :: arr(:), val
 integer, intent(out) :: indx
 integer              :: istart,istop,i

 istart = 1
 istop  = size(arr)
 indx = istart
 if (val >= arr(istop)) then
    indx = istop-1   ! -1 to avoid array index overflow
 elseif (val <= arr(istart)) then
    indx = istart
 else
    i = istart
    do while (i <= istop)
       if (arr(i)>val) then
          indx = i-1
          exit
       endif
       i = i+1
    enddo
 endif

end subroutine find_nearest_index

!-----------------------------------------------------------------------
!+
!  1D linear interpolation routine
!+
!-----------------------------------------------------------------------
pure real function interp_1d(x,x1,x2,y1,y2)
 real, intent(in) :: x, x1, x2, y1, y2

 interp_1d = y1 + (x-x1)*(y2-y1)/(x2-x1)

end function interp_1d

!-----------------------------------------------------------------------
!+
!  similar but just interpolates between two values
!  val0 and val1 where u is the fraction of the way
!+
!-----------------------------------------------------------------------
pure subroutine linear_interpolator_one_d(val0,val1,u,val)
 real, intent(out) :: val
 real, intent(in)  :: val0,val1,u

 val=(1.-u)*val0+u*val1

end subroutine linear_interpolator_one_d

function interpolate_1d(x,datax,datay,dydx) result(y)
 real, intent(in) :: datax(:),datay(:)
 real, intent(in) :: dydx(:) !--input so that it does not need to be recalculated at all calls
 real, intent(in) :: x
 real, allocatable :: x0(:),y0(:)
 real             :: y,dx,dxtest
 logical          :: uniform
 integer          :: i,Nsize,Nref

 Nsize=size(datax(:))
 Nref=0
 allocate(x0(Nsize))
 allocate(y0(Nsize))

 x0 = datax(:)
 y0 = datay(:)
 dx = x0(2)-x0(1)
 dxtest = x0(Nsize)-x0(Nsize-1)

 if (abs(dx-dxtest)<1.E-4) then
    uniform=.true. !--uniform dx
 else
    uniform=.false.
 endif

 if (uniform) then
    Nref=floor((x-x0(1))/dx+1)
    if (Nref<1 .or. Nref>Nsize) then
       print*,'cannot interpolate a value out of the grid: x = ',x,y0(:),Nref
       stop
    endif
 else
    do i=1,Nsize
       if (x>=x0(i)) then
          cycle
       else
          Nref=i-1
          exit
       endif
    enddo
 endif
 y = dydx(Nref)*(x-x0(Nref))+y0(Nref)

!    print*, x,x0(Nref),y,y0(Nref),Nref,dx,dxtest,abs(dx-dxtest)<1.E-4,uniform
 deallocate(x0)
 deallocate(y0)

end function interpolate_1d

subroutine differentiate(y,x,dydx)
 !-- based on numpy gradient
 real,              intent(in)    :: y(:),x(:)
 real, allocatable, intent(inout) :: dydx(:) !will be deallocated in grids_for_setup.f90:deallocate_sigma()
 real, allocatable :: dx(:),dx1(:),dx2(:)
 real, allocatable :: a(:),b(:),c(:)
 integer :: Nsize, Nsizedx

 Nsize = size(x)
 allocate(dydx(Nsize))
 Nsizedx = Nsize-1
 allocate(dx(Nsizedx),dx1(Nsizedx),dx2(Nsizedx))
 allocate(a(Nsize),b(Nsize),c(Nsize))

 dx = x(2:)-x(1:Nsize-1)
 dx1 = dx(1:Nsizedx-1)
 dx2 = dx(2:)

 !--calculate non-edge values
 a = -(dx2(:))/(dx1(:) * (dx1(:) + dx2(:)))
 b = (dx2(:) - dx1(:)) / (dx1(:) * dx2(:))
 c = dx1(:) / (dx2(:) * (dx1(:) + dx2(:)))

 dydx(2:Nsize-1) = a(:) * y(1:Nsize-2) + b(:) * y(2:Nsize-1) + c(:) * y(3:)

 !--calculate edge value 1

 dx1(1) = dx(1)
 dx2(1) = dx(2)
 a(1) = -(2. * dx1(1) + dx2(1))/(dx1(1) * (dx1(1) + dx2(1)))
 b(1) = (dx1(1) + dx2(1)) / (dx1(1) * dx2(1))
 c(1) = - dx1(1) / (dx2(1) * (dx1(1) + dx2(1)))

 dydx(1) = a(1) * y(1) + b(1) * y(2) + c(1) * y(3)

 !-- calculate edge value Nsize

 dx1(1) = dx(Nsizedx-1)
 dx2(1) = dx(Nsizedx)
 a(1) = (dx2(1)) / (dx1(1) * (dx1(1) + dx2(1)))
 b(1) = - (dx2(1) + dx1(1)) / (dx1(1) * dx2(1))
 c(1) = (2. * dx2(1) + dx1(1)) / (dx2(1) * (dx1(1) + dx2(1)))

 dydx(Nsize) = a(1) * y(Nsize-2) + b(1) * y(Nsize-1) + c(1) * y(Nsize)

 deallocate(dx,dx1,dx2)
 deallocate(a,b,c)

!   do i=1,size(x)
 !     print*,i,x(i)
 ! enddo

end subroutine differentiate

!-----------------------------------------------------------------------
!+
!  Tricubic Hermite interpolation of nc fields tabulated on a uniform,
!  vertex-centred 3D grid (node (i,j,k) at x0 + (i-1,j-1,k-1)*dx), given
!  the values f and their first derivatives df at the nodes. Mixed
!  derivatives are taken by centred differences of df.
!
!  The interpolant is C1 and exact for polynomials up to cubic order;
!  values converge at 4th order and gradients at 3rd order. The returned
!  gradient is the exact derivative of the interpolant, so a force
!  computed from it is consistent with the interpolated field (this is
!  what gives energy conservation for geodesics in a tabulated metric).
!
!  For a periodic grid, node n+1 is identified with node 1; otherwise
!  positions outside the grid are extrapolated from the edge cell and
!  one-sided differences are used at the edges.
!+
!-----------------------------------------------------------------------
pure subroutine interp3d_hermite(nc,n,f,df,x0,dx,periodic,x,val,grad)
 integer, intent(in)  :: nc,n(3)
 real,    intent(in)  :: f(nc,n(1),n(2),n(3)),df(nc,3,n(1),n(2),n(3))
 real,    intent(in)  :: x0(3),dx(3),x(3)
 logical, intent(in)  :: periodic
 real,    intent(out) :: val(nc)
 real,    intent(out), optional :: grad(nc,3)
 real    :: t(3),h(0:1,0:1,3),dh(0:1,0:1,3),d(nc)
 integer :: i0(3),ijk(3),cx,cy,cz,a,b,c,idim

 call cell_index_3d(x,x0,dx,n,periodic,i0,t)
 do idim=1,3
    call hermite_basis(t(idim),h(:,:,idim),dh(:,:,idim))
 enddo

 val = 0.
 if (present(grad)) grad = 0.
 do cz=0,1
    do cy=0,1
       do cx=0,1
          ijk(1) = grid_index(i0(1)+cx,n(1),periodic)
          ijk(2) = grid_index(i0(2)+cy,n(2),periodic)
          ijk(3) = grid_index(i0(3)+cz,n(3),periodic)
          do c=0,1
             do b=0,1
                do a=0,1
                   call hermite_nodal_data_3d(nc,n,f,df,dx,periodic,a,b,c,ijk,d)
                   ! convert derivatives to units of the grid spacing
                   d = dx(1)**a * dx(2)**b * dx(3)**c * d
                   val = val + h(a,cx,1)*h(b,cy,2)*h(c,cz,3)*d
                   if (present(grad)) then
                      grad(:,1) = grad(:,1) + dh(a,cx,1)*h(b,cy,2)*h(c,cz,3)*d
                      grad(:,2) = grad(:,2) + h(a,cx,1)*dh(b,cy,2)*h(c,cz,3)*d
                      grad(:,3) = grad(:,3) + h(a,cx,1)*h(b,cy,2)*dh(c,cz,3)*d
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
 enddo
 if (present(grad)) then
    do idim=1,3
       grad(:,idim) = grad(:,idim)/dx(idim)
    enddo
 endif

end subroutine interp3d_hermite

!-----------------------------------------------------------------------
!+
!  lower node index i0 of the grid cell containing x, and the fractional
!  position t within the cell (for a non-periodic grid the edge cell is
!  used for positions outside the grid)
!+
!-----------------------------------------------------------------------
pure subroutine cell_index_3d(x,x0,dx,n,periodic,i0,t)
 real,    intent(in)  :: x(3),x0(3),dx(3)
 integer, intent(in)  :: n(3)
 logical, intent(in)  :: periodic
 integer, intent(out) :: i0(3)
 real,    intent(out) :: t(3)
 real    :: s
 integer :: idim

 do idim=1,3
    s = (x(idim) - x0(idim))/dx(idim)
    i0(idim) = floor(s) + 1
    if (.not.periodic) i0(idim) = min(max(i0(idim),1),max(n(idim)-1,1))
    t(idim) = s - (i0(idim) - 1)
 enddo

end subroutine cell_index_3d

!-----------------------------------------------------------------------
!+
!  map a node index onto the grid: wrap if periodic, clamp otherwise
!+
!-----------------------------------------------------------------------
pure integer function grid_index(i,n,periodic)
 integer, intent(in) :: i,n
 logical, intent(in) :: periodic

 if (periodic) then
    grid_index = modulo(i-1,n) + 1
 else
    grid_index = min(max(i,1),n)
 endif

end function grid_index

!-----------------------------------------------------------------------
!+
!  neighbouring nodes im, ip used for a centred difference at node i, and
!  their separation in units of the grid spacing (one-sided differences
!  at the edges of a non-periodic grid)
!+
!-----------------------------------------------------------------------
pure subroutine fd_neighbours(i,n,periodic,im,ip,nsep)
 integer, intent(in)  :: i,n
 logical, intent(in)  :: periodic
 integer, intent(out) :: im,ip,nsep

 if (periodic) then
    im = grid_index(i-1,n,periodic)
    ip = grid_index(i+1,n,periodic)
    nsep = 2
 else
    im = max(i-1,1)
    ip = min(i+1,n)
    nsep = ip - im
 endif

end subroutine fd_neighbours

!-----------------------------------------------------------------------
!+
!  1D cubic Hermite basis functions h(a,c) and dh/dt on t in [0,1], where
!  c=0,1 is the end of the interval and a=0 multiplies the value, a=1 the
!  derivative (in units of the interval length)
!+
!-----------------------------------------------------------------------
pure subroutine hermite_basis(t,h,dh)
 real, intent(in)  :: t
 real, intent(out) :: h(0:1,0:1),dh(0:1,0:1)
 real :: t2,t3

 t2 = t*t
 t3 = t2*t
 h(0,0) = 2.*t3 - 3.*t2 + 1.
 h(1,0) = t3 - 2.*t2 + t
 h(0,1) = -2.*t3 + 3.*t2
 h(1,1) = t3 - t2
 dh(0,0) = 6.*t2 - 6.*t
 dh(1,0) = 3.*t2 - 4.*t + 1.
 dh(0,1) = -6.*t2 + 6.*t
 dh(1,1) = 3.*t2 - 2.*t

end subroutine hermite_basis

!-----------------------------------------------------------------------
!+
!  nodal data for the tricubic Hermite interpolant at node ijk: the
!  derivative d^(a+b+c) f / dx^a dy^b dz^c for a,b,c in {0,1}. First
!  derivatives come from df; mixed derivatives are centred differences
!  of df in the remaining directions
!+
!-----------------------------------------------------------------------
pure subroutine hermite_nodal_data_3d(nc,n,f,df,dx,periodic,a,b,c,ijk,d)
 integer, intent(in)  :: nc,n(3),a,b,c,ijk(3)
 real,    intent(in)  :: f(nc,n(1),n(2),n(3)),df(nc,3,n(1),n(2),n(3)),dx(3)
 logical, intent(in)  :: periodic
 real,    intent(out) :: d(nc)
 integer :: idir,i,j,k,jm,jp,km,kp,nsepj,nsepk

 i = ijk(1)
 j = ijk(2)
 k = ijk(3)
 if (a+b+c == 0) then
    d = f(:,i,j,k)
    return
 endif
 ! direction of the tabulated first derivative to start from
 if (a==1) then
    idir = 1
 elseif (b==1) then
    idir = 2
 else
    idir = 3
 endif

 if (a==1 .and. b==1 .and. c==1) then
    ! d^3 f/dxdydz = d^2 (df/dx)/dydz
    call fd_neighbours(j,n(2),periodic,jm,jp,nsepj)
    call fd_neighbours(k,n(3),periodic,km,kp,nsepk)
    d = (df(:,1,i,jp,kp) - df(:,1,i,jp,km) - df(:,1,i,jm,kp) + df(:,1,i,jm,km)) &
        /(nsepj*dx(2)*nsepk*dx(3))
 elseif (idir==1 .and. b==1) then
    ! d^2 f/dxdy = d(df/dx)/dy
    call fd_neighbours(j,n(2),periodic,jm,jp,nsepj)
    d = (df(:,1,i,jp,k) - df(:,1,i,jm,k))/(nsepj*dx(2))
 elseif (idir < 3 .and. c==1) then
    ! d^2 f/dxdz or d^2 f/dydz = d(df/dx or df/dy)/dz
    call fd_neighbours(k,n(3),periodic,km,kp,nsepk)
    d = (df(:,idir,i,j,kp) - df(:,idir,i,j,km))/(nsepk*dx(3))
 else
    d = df(:,idir,i,j,k)
 endif

end subroutine hermite_nodal_data_3d

end module table_utils
