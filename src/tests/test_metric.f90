module testmetric
 implicit none
contains

!----------------------------------------------------------------
!+
!  Tests of the metric
!+
!----------------------------------------------------------------
subroutine test_metric_i(gcov,gcon,ntests,npass)
 use testutils, only:checkvalbuf
 use utils_gr,  only:dot_product_gr
 integer, intent(inout)   :: ntests,npass
 real,    intent(in)      :: gcov(0:3,0:3),gcon(0:3,0:3)
 real, dimension(0:3,0:3) :: gg
 real, parameter          :: tol = 6.e-11
 real                     :: sum,errmax
 integer                  :: i,j,error,ncheck

 ntests = ntests+1
 error  = 0

 ! Product of metric and its inverse
 gg = 0.
 gg = matmul(gcov,gcon)
 sum = 0
 do j=0,3
    do i=0,3
       sum = sum + gg(i,j)
    enddo
 enddo

! Check to see that the product is 4 (trace of identity)
 call checkvalbuf(sum,4.,tol,'[F]: gddgUU ',error,ncheck,errmax)

 if (error/=0) then
    print*, 'gdown*gup /= Identity'
    do i=0,3
       write(*,*) gg(i,:)
    enddo
 else
    npass = npass+1
 endif

end subroutine test_metric_i

subroutine test_u0(gcov,gcon,v,ntests,npass)
 use testutils, only:checkvalbuf
 use utils_gr,  only:dot_product_gr
 integer, intent(inout) :: ntests,npass
 real,    intent(in)    :: gcov(0:3,0:3),gcon(0:3,0:3),v(3)
 integer ::error,ncheck
 real    :: gvv,v4(4)

 ntests = ntests+1
 error  = 0

 v4  = (/1.,v/)
 gvv = dot_product_gr(v4,v4,gcov)

 ! Check to see if U0 is imaginary (i.e. bad velocity)
 call checkvalbuf(int(sign(1.,gvv)),-1,0,'[F]: sign of dot_product_gr(v4,v4,gcov))',error,ncheck)

 if (error/=0) then
    print*,'Warning, U0 is imaginary: dot_product_gr(v4,v4,gcov)=',gvv,' > 0'
 else
    npass = npass + 1
 endif

end subroutine test_u0

end module testmetric
