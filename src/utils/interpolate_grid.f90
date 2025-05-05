module interpolate_grid

implicit none

public

contains

function interpolate_1d(x,datax,datay,dydx) result(y)
    real, intent(in) :: datax(:),datay(:)
    real, intent(in) :: dydx(:) !--input so that it does not need to be recalculated at all calls             
    real, intent(in) :: x
    real, dimension(:), allocatable :: x0,y0 
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
    
    if(abs(dx-dxtest)<1.E-4) then  
       uniform=.true. !--uniform dx 
    else
       uniform=.false.
    endif  

    if(uniform) then
       Nref=floor((x-x0(1))/dx+1)
       if (Nref<1 .or. Nref>Nsize) then
          print*,'cannot interpolate a value out of the grid: x = ',x,y0(:),Nref
          stop
       endif
    else
       do i=1,Nsize
          if(x>=x0(i)) then
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
   real, intent(in) :: y(:),x(:)
   real, dimension(:), allocatable :: dx,dx1,dx2
   real, intent(inout), dimension(:), allocatable :: dydx !will be deallocated in grids_for_setup.f90:deallocate_sigma()
   real, dimension(:), allocatable :: a,b,c
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

end module interpolate_grid
