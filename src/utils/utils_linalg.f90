module linearalgebra



implicit none

  INTEGER, PARAMETER :: int64 = SELECTED_INT_KIND(16)
  INTEGER, PARAMETER :: int32 = SELECTED_INT_KIND(8)

  INTEGER, PARAMETER :: real64 = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: real32 = SELECTED_REAL_KIND(6)


contains

function inverse(a_, n) result(x)

    implicit none

    save

    integer(kind=int64), intent(in) :: &
         n
    real(kind=real64), dimension(n, n), intent(in) :: &
         a_

    real(kind=real64), dimension(n, n) :: &
         x

    real(kind=real64), dimension(n, n) :: &
         a

    ! this subroutine finds the inverse of matrix a

    integer(kind=int64) :: &
         k,i,j,l,n1
    real(kind=real64) :: &
         ajji,r

    if (n == 1) then
       x(1,1) = 1.d0 / a_(1,1)
       return
    endif

    a(:,:) = a_(:,:)
    x(:,:) = 0.d0
    do j = 1, n
       x(j,j) = 1.d0
    end do

    n1 = n-1

    ! no conditioning
    ! reduce matrix to upper triangular form

    do j=1,n-1
       ajji = 1.D0 / a(j,j)
       do i=j+1,n
          r = -a(i,j) * ajji
          do k=j+1,n
             a(i,k) = a(i,k) + r * a(j,k)
          enddo
          do k=1,n
             x(i,k) = x(i,k) + r * x(j,k)
          enddo
       end do
    end do

    ! use back substitution to find the vector solution

    do l=1, n
       x(n,l) = x(n,l) / a(n,n)
    end do

    do l=1, n
       do i=n-1, 1, -1
          r = 0.d0
          do j=i+1,n
             r = r + a(i,j) * x(j,l)
          end do
          ! r = dot_product(a(i,i+1:n), x(i+1:n))
          x(i,l) = (x(i,l) - r) / a(i,i)
       end do
    end do

  end function inverse



end module linearalgebra
