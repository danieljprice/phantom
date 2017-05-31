module testcons2prim
implicit none
contains
subroutine test_cons2prim(ntests,npass)
   use eos,          only: gamma
   use cons2prim_gr, only: get_u
   use utils_gr,     only: dot_product_gr
   use metric_tools, only: get_metric
   integer, intent(inout) :: ntests,npass
   integer :: j,ix,iy,iz,vx,vy,vz, ncombinations, ncombpassed
   real :: x(1:3),v(1:3),v4(0:3)
   real :: dens,u,p,gcov(0:3,0:3),gcon(0:3,0:3),sqrtg

   write(*,'(/,a,/)') '--> testing conservative2primitive solver'

   ntests = ntests + 1
   ncombinations = 0
   ncombpassed = 0

   gamma = 5./3.
   dens = 10.
!$omp parallel do default (none) &
!$omp shared(dens) &
!$omp reduction(+:ncombinations,ncombpassed) &
!$omp private(ix,iy,iz,vx,vy,vz,x,gcov,gcon,sqrtg,v,v4,j,u,p)
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
                        do j=0,10
                           p = j*.1
                           call get_u(u,p,dens)
                           call test_cons2prim_i(x,v,dens,u,p,ncombinations,ncombpassed)
                        enddo
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
!$omp end parallel do
   print*,ncombinations,' combinations tried, out of which ',ncombpassed,'passed.'
   if (ncombpassed==ncombinations) npass = npass + 1
   write(*,'(/,a,/)') '<-- conservative2primitive test complete'

end subroutine test_cons2prim

subroutine test_cons2prim_i(x,v,dens,u,p,ntests,npass)
   use cons2prim_gr, only: conservative2primitive,primitive2conservative
   use testutils, only: checkval,checkvalbuf
   ! use checks, only: check
   use testmetric, only: test_metric_i
   real, intent(in) :: x(1:3),v(1:3),dens,u,p
   integer, intent(inout) :: ntests,npass
   real :: rho,pmom(1:3),en
   real :: v_out(1:3),dens_out,u_out,p_out
   real, parameter :: tol = 4.e-12
   ! real, parameter :: tol = 4.e-10
   integer :: nerrors, ierr,metricpass, j
   integer :: ncheck,dummy
   real :: errmax

   ntests = ntests + 1
   nerrors = 0

   ! Used for initial guess in conservative2primitive
   v_out    = v
   dens_out = dens
   u_out    = u
   p_out    = p

   metricpass = 0
   call test_metric_i(x,v,dummy,metricpass)
   if (metricpass==0) then
      print*,'Warning: Metric test failed so cons2prim may also fail...'
   endif

   call primitive2conservative(x,v,dens,u,P,rho,pmom,en,'entropy')
   call conservative2primitive(x,v_out,dens_out,u_out,p_out,rho,pmom,en,ierr,'entropy')

   ! call checkval(ierr,0,0,n_error,'ierr = 0 for convergence')
   call checkvalbuf(ierr,0,0,'[F]: ierr (convergence)',nerrors,ncheck)
   ! nerrors = nerrors + n_error

   ! call checkval(3,v_out,v,tol,n_error,'v_out = v')
   do j=1,3
      call checkvalbuf(v_out(j),v(j),tol,'[F]: v_out',nerrors,ncheck,errmax)
      ! nerrors = nerrors + n_error
   enddo

   ! call checkval(dens_out,dens,tol,n_error,'dens_out = dens')
   call checkvalbuf(dens_out,dens,tol,'[F]: dens_out',nerrors,ncheck,errmax)
   ! nerrors = nerrors + n_error

   ! call checkval(u_out,u,tol,n_error,'u_out = u')
   call checkvalbuf(u_out,u,tol,'[F]: u_out',nerrors,ncheck,errmax)
   ! nerrors = nerrors + n_error

   ! call checkval(p_out,p,tol,n_error,'p_out = p')
   call checkvalbuf(p_out,p,tol,'[F]: p_out',nerrors,ncheck,errmax)
   ! nerrors = nerrors + n_error

   if (nerrors/=0) then
      print*,'-- cons2prim test failed with'
      print*,'  - IN:'
      print*,'     x    =',x
      print*,'     v    =',v
      print*,'     dens =',dens
      print*,'     u    =',u
      print*,'     p    =',p
      print*,'  - OUT:'
      print*,'     v    =',v_out
      print*,'     dens =',dens_out
      print*,'     u    =',u_out
      print*,'     p    =',p_out
      print*,''
   else
      npass = npass + 1
   endif
end subroutine test_cons2prim_i

end module testcons2prim
