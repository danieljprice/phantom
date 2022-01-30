module interpolate_grid

implicit none

public

contains

function interpolate_1d(x,dataset,dydx) result(y)
    real, intent(in) :: dataset(:,:)
    real, intent(in) :: dydx(:) !--input so that it does not need to be recalculated at all calls             
    real, intent(in) :: x
    real, dimension(:), allocatable :: x0,y0 
    real             :: y,dx,dxtest
    logical          :: uniform=.false.
    integer          :: i,Nsize,Nref

    Nsize=size(dataset(:,1))
    allocate(x0(Nsize))
    allocate(y0(Nsize))

    x0 = dataset(:,1)
    y0 = dataset(:,2) 
    dx = x0(2)-x0(1)
    dxtest = x0(Nsize)-x0(Nsize-1)

    if(abs(dx-dxtest)<1.E-4) uniform=.true. !--uniform dx   

    if(uniform) then
       Nref=floor((x-x0(1))/dx+1)
       if (Nref<1.) then
          print*,'cannot interpolate a value out of the grid'
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

!    print*, x,x0(Nref),y,y0(Nref),Nref
    deallocate(x0)
    deallocate(y0)
  
end function interpolate_1d

subroutine differentiate(y,x,dydx) 
   !-- based on numpy gradient
   real, intent(in) :: y(:),x(:)
   real,  dimension(size(y)-1) :: dx,dx1,dx2
   real, dimension(:), allocatable :: dydx
   real, dimension(size(y))    :: a,b,c
   integer :: Nsize, Nsizedx,i

   Nsize = size(x)
   allocate(dydx(Nsize))
   Nsizedx = Nsize-1

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
   
!   do i=1,size(x)
 !     print*,i,x(i)
  ! enddo

end subroutine differentiate 

end module interpolate_grid
