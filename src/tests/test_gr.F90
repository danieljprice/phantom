module testgr
 implicit none

 public :: test_gr

 private

contains

subroutine test_gr(ntests,npass)
 use io,              only:id,master
 use testutils,       only:checkvalbuf,checkvalbuf_end,checkval
 use eos,             only:gamma,equationofstate,ieos
 use utils_gr,        only:dot_product_gr
 use metric_tools,    only:get_metric
 use metric,          only:metric_type
 integer, intent(inout) :: ntests,npass
 real    :: radii(5),theta(5),phi(5),vx(5),vy(5),vz(5)
 real    :: utherm(4),density(4)
 real    :: position(3),v(3),v4(0:3),sqrtg,gcov(0:3,0:3),gcon(0:3,0:3)
 real    :: pi,ri,thetai,phii,vxi,vyi,vzi,x,y,z,p,dens,u,pondens,spsound
 integer :: i,j,k,l,m,n,ii,jj,count
 integer :: ncomb_metric,npass_metric,ncomb_cons2prim,npass_cons2prim

 if (id==master) write(*,"(/,a,/)") '--> TESTING GENERAL RELATIVITY'

 write(*,'(/,a)') '--> testing metric and cons2prim'
 write(*,'(a,/)') '    metric type = '//trim(metric_type)

 ntests = ntests + 1
 ncomb_metric = 0
 npass_metric = 0
 ncomb_cons2prim = 0
 npass_cons2prim = 0

 gamma = 5./3.

 pi     = 4.*atan(1.)

 radii  = (/2.1,2.5,3.0,5.0,10.0/)
 theta  = (/0.,pi/4.,pi/2.,3.*pi/4.,pi/)
 phi    = (/0.,pi/4.,pi/2.,pi,3*pi/2./)

 vx = (/0.,0.25,0.5,0.75,1./)
 vy = vx
 vz = vx

 utherm   = (/0.,2.,10.,100./)
 density  = (/1.,2.,10.,100./)

 do i=1,size(radii)
    ri = radii(i)
    do j=1,size(theta)
       thetai = theta(j)
       do k=1,size(phi)
          phii = phi(k)
          x = ri*sin(thetai)*cos(phii)
          y = ri*sin(thetai)*sin(phii)
          z = ri*cos(thetai)
          position = (/x,y,z/)
          do l=1,size(vx)
             vxi=vx(l)
             do m=1,size(vy)
                vyi=vy(m)
                do n=1,size(vz)
                   vzi=vz(n)

                   call get_metric(position,gcov,gcon,sqrtg)
                   v = (/vxi,vyi,vzi/)
                   v4(0) = 1.
                   v4(1:3) = v(:)

                   ! Only allow valid combinations of position and velocity to be tested.
                   ! i.e. Not faster than the speed of light locally (U0 real, not imaginary).
                   if (dot_product_gr(v4,v4,gcov) < 0.) then
                      count = npass_metric
                      call test_metric_i(gcov,gcon,ncomb_metric,npass_metric)
                      if (npass_metric/=count+1) print*,'Warning: Metric test failed so cons2prim may also fail...'
                      do ii=1,size(utherm)
                         u = utherm(ii)
                         do jj=1,size(density)
                            dens = density(jj)
                            call equationofstate(ieos,pondens,spsound,dens,x,y,z,u)
                            p = pondens*dens
                            call test_cons2prim_i(position,v,dens,u,p,ncomb_cons2prim,npass_cons2prim)
                         enddo
                      enddo
                   endif

                enddo
             enddo
          enddo
       enddo
    enddo
 enddo


 print*,ncomb_metric,' combinations tried, out of which ',npass_metric,'passed the metric test'
 print*,ncomb_cons2prim,' combinations tried, out of which ',npass_cons2prim,'passed the cons2prim test'
 if (npass_metric==ncomb_metric .and. npass_cons2prim==ncomb_cons2prim) npass = npass + 1

 if (id==master) write(*,"(/,a)") '<-- GR TESTS COMPLETE'

end subroutine test_gr

! Indivdual test subroutines start here

!----------------------------------------------------------------
!+
!  Test of the metric
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

!----------------------------------------------------------------
!+
!  Test of U0 (velocity)
!+
!----------------------------------------------------------------
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

subroutine test_cons2prim_i(x,v,dens,u,p,ntests,npass)
 use cons2prim,    only:conservative_to_primitive,primitive_to_conservative
 use testutils,    only:checkval,checkvalbuf
 use metric_tools, only:get_grpacki
 real, intent(in) :: x(1:3),v(1:3),dens,u,p
 integer, intent(inout) :: ntests,npass
 real :: grpacki(0:3,0:3,2)
 real :: rho,pmom(1:3),en
 real :: v_out(1:3),dens_out,u_out,p_out
 real, parameter :: tol = 4.e-12
 integer :: nerrors, ierr, j
 integer :: ncheck
 real :: errmax

 ntests = ntests + 1
 nerrors = 0

 ! Used for initial guess in conservative2primitive
 v_out    = v
 dens_out = dens
 u_out    = u
 p_out    = p

 call get_grpacki(x,grpacki)
 call primitive_to_conservative(x,grpacki,v,dens,u,P,rho,pmom,en)
 call conservative_to_primitive(x,grpacki,v_out,dens_out,u_out,p_out,rho,pmom,en,ierr)

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

end module testgr
