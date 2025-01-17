!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testgr
!
! Unit tests of General Relativity
!
! :References: Liptai & Price (2019), MNRAS
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: cons2prim, cons2primsolver, eos, extern_gr, inverse4x4,
!   io, metric, metric_tools, part, physcon, substepping, testutils, units,
!   utils_gr, vectorutils
!
 use testutils, only:checkval,checkvalbuf,checkvalbuf_end,update_test_scores
 implicit none

 public :: test_gr

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests for General Relativity
!+
!-----------------------------------------------------------------------
subroutine test_gr(ntests,npass)
 use io,      only:id,master
 use units,   only:set_units
 use physcon, only:solarm
 integer, intent(inout)   :: ntests,npass

 call set_units(mass=1.d6*solarm,G=1.d0,c=1.d0)
 if (id==master) write(*,"(/,a,/)") '--> TESTING GENERAL RELATIVITY'
 call test_combinations_all(ntests,npass)
 call test_precession(ntests,npass)
 call test_inccirc(ntests,npass)
 if (id==master) write(*,"(/,a)") '<-- GR TESTS COMPLETE'

end subroutine test_gr

!-----------------------------------------------------------------------
!+
!   Test of orbital precession in the Kerr metric
!+
!-----------------------------------------------------------------------
subroutine test_precession(ntests,npass)
 use metric_tools, only:imetric,imet_kerr,imet_schwarzschild
 use metric,       only:a
 integer, intent(inout) :: ntests,npass
 integer :: nerr(6),norbits,nstepsperorbit
 real    :: dt,period,x0,vy0,tmax,angtol,postol
 real    :: angmom(3),angmom0(3),xyz(3),vxyz(3)

 write(*,'(/,a)') '--> testing substep_gr (precession)'
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

!-----------------------------------------------------------------------
!+
!   Test of inclined circular orbit in the Kerr metric
!+
!-----------------------------------------------------------------------
subroutine test_inccirc(ntests,npass)
 use physcon,      only:pi
 use metric_tools, only:imetric,imet_kerr
 use metric,       only:a
 integer, intent(inout) :: ntests,npass
 integer :: nerr(6),norbits,nstepsperorbit
 real    :: dt,period,tmax
 real    :: angmom(3),angmom0(3),xyz(3),vxyz(3)
 real :: m,omega,phi,q,r,rdot,rho2,theta,thetadot,vx,vy,vz,x1,y1,z1
 real :: R2,rfinal

 write(*,'(/,a)') '--> testing substep_gr (inclined circular orbit)'

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

!-----------------------------------------------------------------------
!+
!   test the geodesic integrator using test particle integration
!   and the substep_gr routine
!+
!-----------------------------------------------------------------------
subroutine integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom)
 use io,             only:iverbose
 use part,           only:igas,npartoftype,massoftype,set_particle_type,get_ntypes,ien_type,&
                          xyzmh_ptmass,vxyz_ptmass,pxyzu_ptmass,metrics_ptmass,&
                          metricderivs_ptmass,fxyz_ptmass,nptmass
 use substepping,    only:substep_gr
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
 ien_type       = 1

 call init_metric(npart,xyzh,metrics,metricderivs)
 call prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,use_dens=.false.,dens=dens)
 call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,fext,dtextforce,dens=dens)
 call calculate_angmom(xyzh(1:3,1),metrics(:,:,:,1),massi,vxyzu(1:3,1),angmom0)

 nsteps = 0
 do while (time <= tmax)
    nsteps = nsteps + 1
    time   = time   + dt
    dtextforce = blah
    !  call substep_gr(npart,ntypes,dt,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,time)
    call substep_gr(npart,nptmass,ntypes,dt,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,time,&
                       xyzmh_ptmass,vxyz_ptmass,pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass,fxyz_ptmass)
 enddo

 call calculate_angmom(xyzh(1:3,1),metrics(:,:,:,1),massi,vxyzu(1:3,1),angmom)

 xyz(:)  = xyzh(1:3,1)
 vxyz(:) = vxyzu(1:3,1)

end subroutine integrate_geodesic

!-----------------------------------------------------------------------
!+
!   compute the angular momentum for the orbit
!+
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
!+
!  Test various combinations of position, velocity and fluid quantities
!+
!-----------------------------------------------------------------------
subroutine test_combinations_all(ntests,npass)
 use eos, only:ieos
 integer, intent(inout) :: ntests,npass
 integer, parameter     :: eos_to_test(2) = (/2,12/)
 integer                :: i

 do i = 1,size(eos_to_test)
    ieos = eos_to_test(i)
    call test_combinations(ntests,npass)
 enddo

end subroutine test_combinations_all

!-----------------------------------------------------------------------
!+
!  Test various combinations of position, velocity and fluid quantities
!+
!-----------------------------------------------------------------------
subroutine test_combinations(ntests,npass)
 use physcon,         only:pi
 use eos,             only:gamma,equationofstate,ieos
 use utils_gr,        only:dot_product_gr
 use metric_tools,    only:get_metric,get_metric_derivs,imetric,imet_kerr
 use metric,          only:metric_type
 integer, intent(inout) :: ntests,npass
 real    :: radii(5),theta(5),phi(5),vx(5),vy(5),vz(5)
 real    :: utherm(7),density(7),errmax,errmaxg,errmaxc,errmaxd
 real    :: position(3),v(3),v4(0:3),sqrtg,gcov(0:3,0:3),gcon(0:3,0:3)
 real    :: ri,thetai,phii,vxi,vyi,vzi,x,y,z,p,t,dens,u,pondens,spsound
 real    :: dgdx1(0:3,0:3),dgdx2(0:3,0:3),dgdx3(0:3,0:3)
 integer :: i,j,k,l,m,n,ii,jj
 integer :: ncheck_metric,nfail_metric,ncheck_cons2prim,nfail_cons2prim
 integer :: ncheckg,nfailg,ncheckd,nfaild
 real, parameter :: tol = 2.e-15
 real, parameter :: tolc = 1.e-12
 real, parameter :: told = 4.e-7

 write(*,'(/,a)')    '--> testing metric and cons2prim with combinations of variables'
 write(*,'(a)')      '    metric type = '//trim(metric_type)
 write(*,'(a,I4,/)') '    eos         = ', ieos

 ntests = ntests + 4
 ncheck_metric = 0
 nfail_metric = 0
 ncheckg = 0
 nfailg  = 0
 ncheck_cons2prim = 0
 nfail_cons2prim = 0
 ncheckd = 0
 nfaild = 0
 errmax = 0.
 errmaxg = 0.
 errmaxc = 0.
 errmaxd = 0.

 ! ieos=12
 gamma = 5./3.

 radii  = (/2.1,2.5,3.0,5.0,10.0/)
 theta  = (/0.,pi/4.,pi/2.,3.*pi/4.,pi/)
 phi    = (/0.,pi/4.,pi/2.,pi,3.*pi/2./)

 vx = (/0.,0.25,0.5,0.75,1./)
 vy = vx
 vz = vx

 utherm   = (/1.e-3,1.,10.,100.,1000.,1.e5,1.e7/)
 density  = (/1.e-10,1.e-5,1.e-3,1.,10.,100.,1000./)

 t = -1. ! initial temperature guess to avoid complier warning

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

          call get_metric(position,gcov,gcon,sqrtg)
          call test_metric_i(gcov,gcon,sqrtg,ncheck_metric,nfail_metric,errmax,ncheckg,nfailg,errmaxg,tol)

          ! Check below is because Kerr metric derivatives currently badly behaved at the poles
          ! Would be nice to remove this...
          if ((imetric /= imet_kerr) .or. (x**2 + y**2 > 1.e-12)) then
             call get_metric_derivs(position,dgdx1,dgdx2,dgdx3)
             call test_metric_derivs_i(position,dgdx1,dgdx2,dgdx3,ncheckd,nfaild,errmaxd,told)
          endif

          do l=1,size(vx)
             vxi=vx(l)
             do m=1,size(vy)
                vyi=vy(m)
                do n=1,size(vz)
                   vzi=vz(n)

                   v = (/vxi,vyi,vzi/)
                   v4(0) = 1.
                   v4(1:3) = v(:)

                   ! Only allow valid combinations of position and velocity to be tested.
                   ! i.e. Not faster than the speed of light locally (U0 real, not imaginary).
                   if (dot_product_gr(v4,v4,gcov) < 0.) then
                      do ii=1,size(utherm)
                         u = utherm(ii)
                         do jj=1,size(density)
                            dens = density(jj)
                            call equationofstate(ieos,pondens,spsound,dens,x,y,z,t,u)
                            p = pondens*dens
                            call test_cons2prim_i(position,v,dens,u,p,ncheck_cons2prim,nfail_cons2prim,errmaxc,tolc)
                         enddo
                      enddo
                   endif

                enddo
             enddo
          enddo
       enddo
    enddo
 enddo

 call checkvalbuf_end('inv * metric = identity',ncheck_metric,nfail_metric,errmax,tol)
 call checkvalbuf_end('sqrt g = -det(g)',ncheckg,nfailg,errmaxg,tol)
 call checkvalbuf_end('d/dx^i g_munu',ncheckd,nfaild,errmaxd,told)
 call checkvalbuf_end('conservative to primitive',ncheck_cons2prim,nfail_cons2prim,errmaxc,tolc)
 if (nfail_metric==0)    npass = npass + 1
 if (nfailg==0)          npass = npass + 1
 if (nfaild==0)          npass = npass + 1
 if (nfail_cons2prim==0) npass = npass + 1

end subroutine test_combinations

!----------------------------------------------------------------
!+
!  Test of the metric
!+
!----------------------------------------------------------------
subroutine test_metric_i(gcov,gcon,sqrtg,ncheck,nfail,errmax,ncheckg,nfailg,errmaxg,tol)
 use inverse4x4, only:inv4x4
 integer, intent(inout)   :: ncheck,nfail,ncheckg,nfailg
 real,    intent(in)      :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg,tol
 real,    intent(inout)   :: errmax,errmaxg
 real, dimension(0:3,0:3) :: gg
 real                     :: sum,det
 integer                  :: i,j

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
 call checkvalbuf(sum,4.,tol,'[F]: gddgUU ',nfail,ncheck,errmax)

 !if (nfail /= 0) then
 !   print*,' metric '
 !   print "(4(es10.3,1x))",gcov
 !   print*,' inverse '
 !   print "(4(es10.3,1x))",gcon
 !   print*,' gg '
 !   print "(4(es10.3,1x))",gg
 !   print*, 'gdown*gup /= Identity'
 !endif

 ! Check that the determinant of the metric matches the one returned
 call inv4x4(gcov,gg,det)
 call checkvalbuf(-det,sqrtg,tol,'sqrt(g) ',nfailg,ncheckg,errmaxg)

end subroutine test_metric_i

!----------------------------------------------------------------
!+
!  Check that analytic metric derivs give similar answer to
!  numerical differences of the metric
!+
!----------------------------------------------------------------
subroutine test_metric_derivs_i(x,dgdx1,dgdx2,dgdx3,ncheck,nfail,errmax,tol)
 use metric_tools, only:numerical_metric_derivs
 real, intent(in) :: x(1:3),dgdx1(0:3,0:3),dgdx2(0:3,0:3),dgdx3(0:3,0:3),tol
 integer, intent(inout) :: ncheck,nfail
 real,    intent(inout) :: errmax
 real :: dgdx_1(0:3,0:3),dgdx_2(0:3,0:3),dgdx_3(0:3,0:3)
 integer :: j,i

 call numerical_metric_derivs(x,dgdx_1,dgdx_2,dgdx_3)
 do j=0,3
    do i=0,3
       call checkvalbuf(dgdx1(i,j),dgdx_1(i,j),tol,'dgcov/dx',nfail,ncheck,errmax)
       call checkvalbuf(dgdx2(i,j),dgdx_2(i,j),tol,'dgcov/dy',nfail,ncheck,errmax)
       call checkvalbuf(dgdx3(i,j),dgdx_3(i,j),tol,'dgcov/dz',nfail,ncheck,errmax)
    enddo
 enddo

end subroutine test_metric_derivs_i

!----------------------------------------------------------------
!+
!  Test of the conservative to primitive variable solver
!+
!----------------------------------------------------------------
subroutine test_cons2prim_i(x,v,dens,u,p,ncheck,nfail,errmax,tol)
 use cons2primsolver, only:conservative2primitive,primitive2conservative
 use part,            only:ien_entropy,ien_etotal,ien_entropy_s
 use metric_tools,    only:pack_metric,unpack_metric
 use eos,             only:ieos,equationofstate,calc_temp_and_ene
 use physcon,         only:radconst

 real, intent(in) :: x(1:3),v(1:3),dens,p,tol
 real,    intent(inout) :: u
 integer, intent(inout) :: ncheck,nfail
 real,    intent(inout) :: errmax
 real :: metrici(0:3,0:3,2)
 real :: rho2,pmom2(1:3),en2,temp
 real :: p2,u2,t2,dens2,gamma2,v2(1:3)
 real :: pondens2,spsound2
 real :: v_out(1:3),dens_out,u_out,p_out,gamma_out
 real :: toli
 integer :: ierr,i,j,nfailprev,ien_type
 real, parameter :: tolg = 1.e-7, tolp = 1.5e-6

 ! perturb the state
 dens2 = 2.*dens
 u2 = 2.*u
 t2 = -1.

 call equationofstate(ieos,pondens2,spsound2,dens2,x(1),x(2),x(3),t2,u2)
 P2 = pondens2 * dens2
 v2 = v

 over_energy_variables: do i = 1,3
    ! Used for initial guess in conservative2primitive
    v_out    = v
    dens_out = dens
    u_out    = u
    p_out    = p
    gamma_out = 1. + p/(dens*u)
    errmax   = 0.
    nfailprev = nfail
    temp = 1.e7 ! arbitrary initial guess
    gamma2 = 1. + P2/(dens2*u2)

    call pack_metric(x,metrici)
    if (ieos == 12 .and. i /= 3) then
       ! entropy_K and etotal cannot use with gasplusrad eos
       cycle
    elseif (i == 1) then
       ien_type = ien_entropy
       toli = 1.5e-11
    elseif (i == 2) then
       ien_type = ien_etotal
       toli = 1.5e-9
    else
       ien_type = ien_entropy_s
       toli = 1.5e-11
    endif

    call primitive2conservative(x,metrici,v,dens2,u2,P2,rho2,pmom2,en2,ien_type)
    call conservative2primitive(x,metrici,v_out,dens_out,u_out,p_out,temp,gamma_out,rho2,pmom2,en2,ierr,ien_type)

    call checkvalbuf(ierr,0,0,'[F]: ierr (convergence)',nfail,ncheck)
    do j=1,3
       call checkvalbuf(v_out(j),v2(j),toli,'[F]: v_out',nfail,ncheck,errmax)
    enddo
    call checkvalbuf(dens_out,dens2,toli,'[F]: dens_out',nfail,ncheck,errmax)
    call checkvalbuf(u_out,u2,toli,'[F]: u_out',nfail,ncheck,errmax)
    call checkvalbuf(p_out,p2,tolp,'[F]: p_out',nfail,ncheck,errmax)
    call checkvalbuf(gamma_out,gamma2,tolg,'[F]: gamma_out',nfail,ncheck,errmax)

    if (nfail > nfailprev .and. nfail < 10) then
       print*,'-- cons2prim test failed with'
       print*,'   ien_type =',ien_type
       print*,'   ieos     =',ieos
       print*,'  - IN:'
       print*,'     x    =',x
       print*,'     v    =',v2
       print*,'     dens =',dens2
       print*,'     u    =',u2
       print*,'     p    =',p2
       print*,'     gamma=',gamma2
       print*,'  - OUT:'
       print*,'     v    =',v_out
       print*,'     dens =',dens_out
       print*,'     u    =',u_out
       print*,'     p    =',p_out
       print*,'     gamma=',gamma_out
       print*,''
    endif
 enddo over_energy_variables

end subroutine test_cons2prim_i

end module testgr
