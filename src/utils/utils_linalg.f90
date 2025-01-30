!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module linalg
!
! linalg
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: None
!

 implicit none
contains
function inverse(matrix,n)
 integer, intent(in)::n
 real, dimension(n,n),intent(in) :: matrix
 real, dimension(n,2*n) :: a,temp
 integer::i,j,k
 real :: ratio,divisor
 real, dimension(n,n) :: inverse

 a(:,:) = 0.
 inverse(:,:) = 0.

 !Augmenting Identity Matrix of Order n
 do i=1,n
    do j=1,n
       a(i,j) = matrix(i,j)
    enddo
 enddo

 !we should do a row swap if the elements on diagonal are 0

 do i=1,n
    do j=1,n
       if (i==j) then
          a(i,j+n) = 1.
       endif
    enddo
 enddo

 !we swap the rows if we enounter 0 as diagonal element
 temp = a(:,:)
 do i=1,n
    if (a(i,i)==0.0) then
       do j=1,n
          if (j  /=  i) then
             if (a(j,i)  /=  0.) then
                a(j,:) = temp(i,:)
                a(i,:) = temp(j,:)
                temp = a(:,:)
                exit
             endif
          endif
       enddo
    endif
 enddo

 !Applying Guass Jordan Elimination
 do i=1,n
    do j=1,n
       if (i  /=  j) then
          ratio = a(j,i)/a(i,i)
          do k=1,2*n
             a(j,k) = a(j,k) - ratio*a(i,k)
          enddo
       endif
    enddo
 enddo

 !dividing by the diagonal elements to get identity matrix
 do i=1,n
    divisor = a(i,i)
    do j=1,2*n
       a(i,j) = a(i,j)/divisor
    enddo
 enddo
 do i=1,n
    do j=1,n
       inverse(i,j)=a(i,j+n)
    enddo
 enddo

end function inverse

end module linalg
