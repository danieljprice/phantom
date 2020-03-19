module testgr
 implicit none

 public :: test_gr

 private

contains

subroutine test_gr(ntests,npass)
 use io,  only:id,master
 integer, intent(inout)   :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING GENERAL RELATIVITY'
 call test_combinations(ntests,npass)
 call test_precession(ntests,npass)
 call test_inccirc(ntests,npass)
 if (id==master) write(*,"(/,a)") '<-- GR TESTS COMPLETE'

end subroutine test_gr

! Indivdual test subroutines start here
subroutine test_precession(ntests,npass)
 use testutils,    only:checkval,update_test_scores
 use metric_tools, only:imetric,imet_kerr,imet_schwarzschild
#ifdef KERR
 use metric,       only:a
#else
 real :: a
#endif
 integer, intent(inout) :: ntests,npass
 integer :: nerr(6),norbits,nstepsperorbit
 real    :: dt,period,x0,vy0,tmax,angtol,postol
 real    :: angmom(3),angmom0(3),xyz(3),vxyz(3)

 write(*,'(/,a)') '--> testing step_extern_gr (precession)'
 if (imetric /= imet_kerr .and. imetric /= imet_schwarzschild) then
    write(*,'(/,a)') '   Skipping test! Metric is not Kerr (or Schwarzschild).'
    return
 endif

 a              = 0.
 x0             = 90.
 vy0            = 0.0521157
 xyz            = (/x0,0. ,0./)
 vxyz           = (/0.,vy0,0./)
 period         = 2390. ! approximate
 norbits        = 4
 tmax           = norbits*period
 nstepsperorbit = 1000
 dt             = 0.239 !period/nstepsperorbit

 call integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom)

 angtol = 1.08e-15
 postol = 1.4e-5
 call checkval(angmom(1),angmom0(1),angtol,nerr(1),'error in angmomx')
 call checkval(angmom(2),angmom0(2),angtol,nerr(2),'error in angmomy')
 call checkval(angmom(3),angmom0(3),angtol,nerr(3),'error in angmomz')
 call checkval(xyz(1), 77.606726748045929,postol,nerr(4),'error in final x position')
 call checkval(xyz(2),-45.576259888019351,postol,nerr(5),'error in final y position')
 call checkval(xyz(3),0.0                ,postol,nerr(6),'error in final z position')

 call update_test_scores(ntests,nerr,npass)

end subroutine test_precession

subroutine test_inccirc(ntests,npass)
 use physcon,   only:pi
 use testutils, only:checkval,update_test_scores
 use metric_tools, only:imetric,imet_kerr
#ifdef KERR
 use metric,    only:a
#else
 real :: a
#endif
 integer, intent(inout) :: ntests,npass
 integer :: nerr(6),norbits,nstepsperorbit
 real    :: dt,period,tmax
 real    :: angmom(3),angmom0(3),xyz(3),vxyz(3)
 real :: m,omega,phi,q,r,rdot,rho2,theta,thetadot,vx,vy,vz,x1,y1,z1
 real :: R2,rfinal

 write(*,'(/,a)') '--> testing step_extern_gr (inclined circular orbit)'

 if (imetric /= imet_kerr) then
    write(*,'(/,a)') '   Skipping test! Metric is not Kerr.'
    return
 endif

 a        = 1.
 r        = 10.
 theta    = 45.*pi/180. ! convert to radians
 phi      = 0.
 m        = 1.
 q        = sqrt(r**2 - a**2*cos(theta)**2)
 rho2     = r**2 + a**2*cos(theta)**2
 omega    = q*sqrt(m)/(sin(theta)*(rho2*sqrt(r)+a*q*sqrt(m)*sin(theta))) !shakura 1987
 rdot     = 0.
 thetadot = 0.

 ! Cartesian coordinates
 x1 = sqrt(r**2+a**2)*sin(theta)*cos(phi)
 y1 = sqrt(r**2+a**2)*sin(theta)*sin(phi)
 z1 = r*cos(theta)
 vx = r/sqrt(r**2+a**2)*sin(theta)*cos(phi)*rdot + sqrt(r**2+a**2)*(cos(theta)*cos(phi)*thetadot-sin(theta)*sin(phi)*omega)
 vy = r/sqrt(r**2+a**2)*sin(theta)*sin(phi)*rdot + sqrt(r**2+a**2)*(cos(theta)*sin(phi)*thetadot+sin(theta)*cos(phi)*omega)
 vz = cos(theta)*rdot-r*sin(theta)*thetadot

 xyz  = (/x1,y1,z1/)
 vxyz = (/vx,vy,vz/)

 period         = 2390. ! approximate
 norbits        = 4
 tmax           = norbits*period
 nstepsperorbit = 1000
 dt             = 0.239 !period/nstepsperorbit

 call integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom)

 R2     = dot_product(xyz,xyz)
 rfinal = sqrt(0.5*(R2-a**2) + 0.5*sqrt((R2-a**2)**2 + 4.*a**2*xyz(3)**2))

 nerr = 0
 call checkval(angmom(1),angmom0(1),6.e-10,nerr(1),'error in angmomx')
 call checkval(angmom(2),angmom0(2),6.e-10,nerr(2),'error in angmomy')
 call checkval(angmom(3),angmom0(3),6.e-10,nerr(3),'error in angmomz')
 call checkval(rfinal   ,r         ,5.08e-6,nerr(4),'error in final r position')

 call update_test_scores(ntests,nerr,npass)

end subroutine test_inccirc

subroutine integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom)
 use io,             only:iverbose
 use part,           only:igas,npartoftype,massoftype,set_particle_type,get_ntypes
 use step_lf_global, only:step_extern_gr
 use eos,            only:ieos
 use cons2prim,      only:prim2consall
 use metric_tools,   only:init_metric,unpack_metric
 use extern_gr,      only:get_grforce_all
 real, intent(in) :: tmax,dt
 real, intent(inout) :: xyz(3), vxyz(3)
 real, intent(out)   :: angmom0(3),angmom(3)
 integer :: nsteps,ntypes,npart
 real    :: time,dtextforce,massi,blah
 real    :: xyzh(4,1),vxyzu(4,1),fext(3,1),pxyzu(4,1),dens(1),metrics(0:3,0:3,2,1),metricderivs(0:3,0:3,3,1)

 npart        = 1

 xyzh         = 0.
 vxyzu        = 0.
 pxyzu        = 0.
 fext         = 0.
 metrics      = 0.
 metricderivs = 0.

 xyzh(1:3,1)  = xyz(:)
 vxyzu(1:3,1) = vxyz(:)
 xyzh(4,:)    = 1.
 vxyzu(4,:)   = 0.
 massi        = 1.e-10
 call set_particle_type(1,igas)

 npartoftype(igas) = npart
 massoftype(igas)  = massi
 ntypes            = get_ntypes(npartoftype)

 !
 ! initialise runtime parameters
 !
 ieos           = 11
 iverbose       = 1
 time           = 0
 blah           = dt

 call init_metric(npart,xyzh,metrics,metricderivs)
 call prim2consall(npart,xyzh,metrics,vxyzu,dens,pxyzu,use_dens=.false.)
 call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,dens,fext,dtextforce)
 call calculate_angmom(xyzh(1:3,1),metrics(:,:,:,1),massi,vxyzu(1:3,1),angmom0)

 nsteps = 0
 do while (time <= tmax)
    nsteps = nsteps + 1
    time   = time   + dt
    dtextforce = blah
    call step_extern_gr(npart,ntypes,dt,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,time)
 enddo

 call calculate_angmom(xyzh(1:3,1),metrics(:,:,:,1),massi,vxyzu(1:3,1),angmom)

 xyz(:)  = xyzh(1:3,1)
 vxyz(:) = vxyzu(1:3,1)

end subroutine integrate_geodesic

subroutine calculate_angmom(xyzi,metrici,massi,vxyzi,angmomi)
 use metric_tools, only:unpack_metric
 use vectorutils,  only:cross_product3D
 use utils_gr,     only:dot_product_gr
 real, intent(in)  :: xyzi(3),metrici(:,:,:),massi,vxyzi(3)
 real, intent(out) :: angmomi(3)
 real              :: alpha_gr,beta_gr_UP(3),bigvi(3),fourvel_space(3),lorentzi,v2i,gammaijdown(3,3)

 call unpack_metric(metrici,betaUP=beta_gr_UP,alpha=alpha_gr,gammaijdown=gammaijdown)
 bigvi    = (vxyzi+beta_gr_UP)/alpha_gr
 v2i      = dot_product_gr(bigvi,bigvi,gammaijdown)
 lorentzi = 1./sqrt(1.-v2i)
 fourvel_space = (lorentzi/alpha_gr)*vxyzi
 call cross_product3D(xyzi,fourvel_space,angmomi) ! position cross with four-velocity
 angmomi=angmomi*massi

end subroutine calculate_angmom

subroutine test_combinations(ntests,npass)
 use physcon,         only:pi
 use testutils,       only:checkvalbuf,checkvalbuf_end,checkval
 use eos,             only:gamma,equationofstate,ieos
 use utils_gr,        only:dot_product_gr
 use metric_tools,    only:get_metric
 use metric,          only:metric_type
 integer, intent(inout) :: ntests,npass
 real    :: radii(5),theta(5),phi(5),vx(5),vy(5),vz(5)
 real    :: utherm(4),density(4)
 real    :: position(3),v(3),v4(0:3),sqrtg,gcov(0:3,0:3),gcon(0:3,0:3)
 real    :: ri,thetai,phii,vxi,vyi,vzi,x,y,z,p,dens,u,pondens,spsound
 integer :: i,j,k,l,m,n,ii,jj,count
 integer :: ncomb_metric,npass_metric,ncomb_cons2prim,npass_cons2prim
 write(*,'(/,a)') '--> testing metric and cons2prim with combinations of variables'
 write(*,'(a,/)') '    metric type = '//trim(metric_type)

 ntests = ntests + 1
 ncomb_metric = 0
 npass_metric = 0
 ncomb_cons2prim = 0
 npass_cons2prim = 0

 gamma = 5./3.

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

end subroutine test_combinations

!----------------------------------------------------------------
!+
!  Test of the metric
!+
!----------------------------------------------------------------
subroutine test_metric_i(gcov,gcon,ntests,npass)
 use testutils, only:checkvalbuf
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

subroutine test_cons2prim_i(x,v,dens,u,p,ntests,npass)
 use cons2primsolver, only:conservative2primitive,primitive2conservative,ien_entropy
 use testutils,       only:checkval,checkvalbuf,update_test_scores
 use metric_tools,    only:pack_metric
 use eos,             only:gamma
 real, intent(in) :: x(1:3),v(1:3),dens,u,p
 integer, intent(inout) :: ntests,npass
 real :: metrici(0:3,0:3,2)
 real :: rho,pmom(1:3),en
 real :: v_out(1:3),dens_out,u_out,p_out
 real, parameter :: tol = 4.e-12
 integer :: nerrors(1), ierr, j
 integer :: ncheck
 real :: errmax

 nerrors = 0

 ! Used for initial guess in conservative2primitive
 v_out    = v
 dens_out = dens
 u_out    = u
 p_out    = p

 call pack_metric(x,metrici)
 call primitive2conservative(x,metrici,v,dens,u,P,rho,pmom,en,ien_entropy,gamma)
 call conservative2primitive(x,metrici,v_out,dens_out,u_out,p_out,rho,pmom,en,ierr,ien_entropy,gamma)

 ! call checkval(ierr,0,0,n_error,'ierr = 0 for convergence')
 call checkvalbuf(ierr,0,0,'[F]: ierr (convergence)',nerrors(1),ncheck)
 ! nerrors = nerrors + n_error

 ! call checkval(3,v_out,v,tol,n_error,'v_out = v')
 do j=1,3
    call checkvalbuf(v_out(j),v(j),tol,'[F]: v_out',nerrors(1),ncheck,errmax)
    ! nerrors = nerrors + n_error
 enddo

 ! call checkval(dens_out,dens,tol,n_error,'dens_out = dens')
 call checkvalbuf(dens_out,dens,tol,'[F]: dens_out',nerrors(1),ncheck,errmax)
 ! nerrors = nerrors + n_error

 ! call checkval(u_out,u,tol,n_error,'u_out = u')
 call checkvalbuf(u_out,u,tol,'[F]: u_out',nerrors(1),ncheck,errmax)
 ! nerrors = nerrors + n_error

 ! call checkval(p_out,p,tol,n_error,'p_out = p')
 call checkvalbuf(p_out,p,tol,'[F]: p_out',nerrors(1),ncheck,errmax)
 ! nerrors = nerrors + n_error

 call update_test_scores(ntests,nerrors,npass)
 if (nerrors(1)/=0) then
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
 endif
end subroutine test_cons2prim_i

end module testgr
