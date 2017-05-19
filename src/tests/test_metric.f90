module testmetric
implicit none
contains
!----------------------------------------------------------------
!+
!  Subroutine to test the metric for different values of x and v
!  It should be used in the test suite.
!+
!----------------------------------------------------------------
subroutine test_metric(ntests,npass)
   use metric_tools, only: get_metric
   use metric, only: metric_type
   use utils_gr, only: get_metric3plus1, dot_product_gr
   integer, intent(inout) :: ntests,npass
   real :: x(1:3), v(1:3), v4(0:3)
   real :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
   integer :: ierr, ix,iy,iz,vx,vy,vz
   ierr = 0

   write(*,'(/,a)') '--> testing metric'
   write(*,'(a,/)') '    metric type = '//trim(metric_type)

   do ix=0,50
      do iy=0,50
         do iz=0,50
            do vx=0,10
               do vy=0,10
                  do vz=0,10
                     x = (/0.1*ix,0.1*iy,0.1*iz/)
                     call get_metric(x,gcov,gcon,sqrtg)
                     v = (/0.1*vx,0.1*vy,0.1*vz/)
                     v4(0) = 1.
                     v4(1:3) = v(:)
                     ! Only allow valid combinations of position and velocity to be tested.
                     ! i.e. Not faster than the speed of light locally.
                     if (dot_product_gr(v4,v4,gcov) < 0.) then
                        ! print*,dot_product_gr(v4,v4,gcov)
                        call test_metric_i(x,v,ntests,npass)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   write(*,'(/,a,/)') '<-- metric test complete'

end subroutine test_metric


!----------------------------------------------------------------
!+
!  Subroutine containing the actual tests of the metric
!+
!----------------------------------------------------------------
subroutine test_metric_i(x,v,ntests,npass)
   use testutils, only: checkvalbuf
   use metric_tools, only: get_metric
   use utils_gr, only: dot_product_gr
   integer, intent(inout) :: ntests,npass
   real, intent(in) :: x(1:3), v(1:3)
   real, dimension(0:3,0:3) :: gcov,gcon,gg
   real :: sqrtg, v4(0:3), sum
   real, parameter :: tol=6.e-11
   ! real, parameter :: tol=4.e-13
   integer :: i,j, nerrors,ncheck,n_error
   real :: errmax

   ntests=ntests+1
   nerrors = 0

   call get_metric(x,gcov,gcon,sqrtg)
   gg = 0.
   gg = matmul(gcov,gcon)
   sum = 0
   do j=0,3
      do i=0,3
         sum = sum + gg(i,j)
      enddo
   enddo

   n_error = 0
   call checkvalbuf(sum,4.,tol,'[F]: gddgUU ',n_error,ncheck,errmax)
   nerrors = nerrors+n_error
   if (n_error>0) then
      print*, 'gdown*gup /= Identity'
      do i=0,3
         write(*,*) gg(i,:)
      enddo
   endif

   v4=(/1.,v/)
   call checkvalbuf(int(sign(1.,dot_product_gr(v4,v4,gcov))),-1,0,'[F]: &
   &sign of dot_product_gr(v4,v4,gcov))',nerrors,ncheck)
   if (dot_product_gr(v4,v4,gcov) > 0.) then
      nerrors = nerrors + 1
      print*,'Warning: Bad combination of position and velocity... &
      &dot_product_gr(v4,v4,gcov)=',dot_product_gr(v4,v4,gcov),' > 0'
   endif

   if (nerrors/=0) then
      print*,'-- Metric test failed with:'
      print*,'    x =',x
      print*,'    v =',v
      print*,''
   else
      npass  = npass + 1
   endif

end subroutine test_metric_i

end module testmetric
