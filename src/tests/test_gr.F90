module testgr
 implicit none
 public :: test_gr

 private

contains

subroutine test_gr(ntests,npass)
 use io,              only:id,master
 use testutils,       only:checkvalbuf,checkvalbuf_end,checkval
 use testmetric,      only:test_metric_i
 use testcons2prim,   only:test_cons2prim_i
 use eos,             only:gamma,equationofstate,ieos
 use utils_gr,        only:dot_product_gr
 use metric_tools,    only:get_metric
 use metric,          only:metric_type
 integer, intent(inout) :: ntests,npass
 real :: radii(5),theta(5),phi(5),vx(5),vy(5),vz(5)
 real :: utherm(4),density(4)
 real :: position(3),v(3),v4(0:3),sqrtg,gcov(0:3,0:3),gcon(0:3,0:3)
 real :: pi,ri,thetai,phii,vxi,vyi,vzi,x,y,z,p,dens,u,pondens,spsound
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
                   ! i.e. Not faster than the speed of light locally.
                   if (dot_product_gr(v4,v4,gcov) < 0.) then
                      count = npass_metric
                      call test_metric_i(position,gcov,gcon,v,ncomb_metric,npass_metric,checkxv=.false.)
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

end module testgr
