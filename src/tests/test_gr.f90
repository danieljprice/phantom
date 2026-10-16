!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testgr
!
! Unit tests of General Relativity
!
! :References: Liptai & Price (2019), MNRAS
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: cons2prim, cons2primsolver, eos, eos_idealplusrad,
!   extern_gr, inverse4x4, io, mesa_microphysics, metric, metric_tools,
!   part, physcon, ptmass, substepping, table_utils, testutils,
!   timestep_ind, timing, units, utils_gr, vectorutils
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
 integer, intent(inout) :: ntests,npass

 call set_units(mass=1.d6*solarm,G=1.d0,c=1.d0)
 if (id==master) write(*,"(/,a,/)") '--> TESTING GENERAL RELATIVITY'
 call test_combinations_all(ntests,npass)
 call test_precession(ntests,npass)
 call test_inccirc(ntests,npass)
 call test_rn_charged(ntests,npass)
 call test_u_los(ntests,npass)
 call test_entropy_step(ntests,npass)
 call test_tde_entropy(ntests,npass)
 call test_mesa_forward_backward(ntests,npass)
 if (id==master) write(*,"(/,a)") '<-- GR TESTS COMPLETE'

end subroutine test_gr

!-----------------------------------------------------------------------
!+
!   Test of orbital precession in the Kerr/Schwarzschild/RN metric
!+
!-----------------------------------------------------------------------
subroutine test_precession(ntests,npass)
 use metric_tools, only:imetric,imet_kerr,imet_schwarzschild,imet_rn
 use metric,       only:a,charge
 integer, intent(inout) :: ntests,npass
 integer :: nerr(6),norbits,nstepsperorbit,itest
 real    :: dt,period,x0,vy0,tmax,angtol,postol
 real    :: angmom(3),angmom0(3),xyz(3),vxyz(3)

 write(*,'(/,a)') '--> testing substep_gr (precession)'
 if (imetric /= imet_kerr .and. imetric /= imet_schwarzschild .and. imetric /= imet_rn) then
    write(*,'(/,a)') '   Skipping test! Metric is not Kerr, Schwarzschild or RN.'
    return
 endif

 a              = 0.
 charge         = 0.
 x0             = 90.
 vy0            = 0.0521157
 period         = 2390. ! approximate
 norbits        = 4
 tmax           = norbits*period
 nstepsperorbit = 1000
 dt             = 0.239 !period/nstepsperorbit

 do itest=1,2
    xyz  = (/x0,0. ,0./)
    vxyz = (/0.,vy0,0./)
    call integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom,use_sink=(itest==2))

    angtol = 1.08e-15
    postol = 1.4e-5
    call checkval(angmom(1),angmom0(1),angtol,nerr(1),'error in angmomx')
    call checkval(angmom(2),angmom0(2),angtol,nerr(2),'error in angmomy')
    call checkval(angmom(3),angmom0(3),angtol,nerr(3),'error in angmomz')
    call checkval(xyz(1), 77.606726748045929,postol,nerr(4),'error in final x position')
    call checkval(xyz(2),-45.576259888019351,postol,nerr(5),'error in final y position')
    call checkval(xyz(3),0.0                ,postol,nerr(6),'error in final z position')
 enddo

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
 integer :: nerr(6),norbits,nstepsperorbit,itest
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

 period         = 2390. ! approximate
 norbits        = 4
 tmax           = norbits*period
 nstepsperorbit = 1000
 dt             = 0.239 !period/nstepsperorbit

 do itest=1,2
    xyz  = (/x1,y1,z1/)
    vxyz = (/vx,vy,vz/)
    call integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom,use_sink=(itest==2))

    R2     = dot_product(xyz,xyz)
    rfinal = sqrt(0.5*(R2-a**2) + 0.5*sqrt((R2-a**2)**2 + 4.*a**2*xyz(3)**2))

    nerr = 0
    call checkval(angmom(1),angmom0(1),6.e-10,nerr(1),'error in angmomx')
    call checkval(angmom(2),angmom0(2),6.e-10,nerr(2),'error in angmomy')
    call checkval(angmom(3),angmom0(3),6.e-10,nerr(3),'error in angmomz')
    call checkval(rfinal   ,r         ,5.08e-6,nerr(4),'error in final r position')
 enddo

 call update_test_scores(ntests,nerr,npass)

end subroutine test_inccirc

!-----------------------------------------------------------------------
!+
!   test the geodesic integrator using test particle integration
!   and the substep_gr routine
!+
!-----------------------------------------------------------------------
subroutine integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom,use_sink)
 use io,             only:iverbose,id,master
 use part,           only:igas,npartoftype,massoftype,set_particle_type,get_ntypes,ien_type,&
                          xyzmh_ptmass,vxyz_ptmass,pxyzu_ptmass,metrics_ptmass,&
                          metricderivs_ptmass,fxyz_ptmass,nptmass,&
                          fxyz_ptmass_tree,dsdt_ptmass,dptmass,fsink_old,ibin_wake,gtgrad,group_info, &
                          bin_info,nmatrix,n_group,n_ingroup,n_sing,rho,rhoh
 use substepping,    only:substep_gr
 use eos,            only:ieos,gamma
 use cons2prim,      only:prim2consall
 use metric_tools,   only:init_metric,unpack_metric
 use extern_gr,      only:get_grforce_all
 use timestep_ind,   only:nbinmax
 use ptmass,         only:use_fourthorder,set_integration_precision
 real,    intent(in)    :: tmax,dt
 real,    intent(inout) :: xyz(3), vxyz(3)
 real,    intent(out)   :: angmom0(3),angmom(3)
 logical, intent(in)    :: use_sink
 integer :: nsteps,ntypes,npart
 real    :: time,dtextforce,massi,blah
 real    :: xyzh(4,1),vxyzu(4,1),fext(3,1),pxyzu(4,1),dens(1),metrics(0:3,0:3,2,1),metricderivs(0:3,0:3,3,1)

 npartoftype = 0
 npart = 0
 xyzh         = 0.
 vxyzu        = 0.
 pxyzu        = 0.
 fext         = 0.
 metrics      = 0.
 metricderivs = 0.

 nptmass = 0
 xyzmh_ptmass = 0.
 vxyz_ptmass = 0.
 pxyzu_ptmass = 0.
 metrics_ptmass = 0.
 metricderivs_ptmass = 0.
 fxyz_ptmass = 0.

 if (use_sink) then
    ! use a single sink particle
    nptmass = 1
    xyzmh_ptmass(1:3,1) = xyz(:)
    vxyz_ptmass(1:3,1) = vxyz(:)
    xyzmh_ptmass(4,1)  = 1.e-10
    xyzmh_ptmass(5,1)  = 1.
 else
    ! use a single gas particle
    npart        = 1
    xyzh(1:3,1)  = xyz(:)
    vxyzu(1:3,1) = vxyz(:)
    xyzh(4,1)    = 1.
    vxyzu(4,1)   = 0.
    massi        = 1.e-10
    call set_particle_type(1,igas)
    npartoftype(igas) = npart
    massoftype(igas)  = massi
    rho(1) = rhoh(xyzh(4,1),massi)
 endif

 ntypes = get_ntypes(npartoftype)

 !
 ! initialise runtime parameters
 !
 ieos           = 11  ! zero pressure equation of state
 iverbose       = 1
 time           = 0
 blah           = dt
 ien_type       = 1
 gamma          = 5./3.

 if (use_sink) then
    if (id==master) print*,'   using sink'
    call init_metric(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass)
    call prim2consall(nptmass,xyzmh_ptmass,metrics_ptmass,vxyz_ptmass,pxyzu_ptmass,use_sink=.true.)
    call get_grforce_all(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass,vxyz_ptmass,fxyz_ptmass,dtextforce,use_sink=.true.)
    call calculate_angmom(xyzmh_ptmass(1:3,1),metrics_ptmass(:,:,:,1),xyzmh_ptmass(4,1),vxyz_ptmass(1:3,1),angmom0)
 else
    if (id==master) print*,'   using gas'
    call init_metric(npart,xyzh,metrics,metricderivs)
    call prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,dens=dens)
    call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,fext,dtextforce,dens=dens)
    call calculate_angmom(xyzh(1:3,1),metrics(:,:,:,1),massi,vxyzu(1:3,1),angmom0)
 endif

 nsteps = 0
 use_fourthorder = .false.
 call set_integration_precision()
 do while (time <= tmax)
    nsteps = nsteps + 1
    time   = time   + dt
    dtextforce = blah
    call substep_gr(npart,ntypes,nptmass,dt,dtextforce,time,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext, &
                    xyzmh_ptmass,vxyz_ptmass,pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass,fxyz_ptmass,&
                    fxyz_ptmass_tree,dsdt_ptmass,dptmass,fsink_old,nbinmax,ibin_wake,gtgrad,group_info, &
                    bin_info,nmatrix,n_group,n_ingroup,n_sing)
 enddo

 if (use_sink) then
    call calculate_angmom(xyzmh_ptmass(1:3,1),metrics_ptmass(:,:,:,1),xyzmh_ptmass(4,1),vxyz_ptmass(1:3,1),angmom)
    xyz(:)  = xyzmh_ptmass(1:3,1)
    vxyz(:) = vxyz_ptmass(1:3,1)
 else
    call calculate_angmom(xyzh(1:3,1),metrics(:,:,:,1),massi,vxyzu(1:3,1),angmom)
    xyz(:)  = xyzh(1:3,1)
    vxyz(:) = vxyzu(1:3,1)
 endif

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

 if (ieos == 12) call benchmark_cons2prim_entropy()

end subroutine test_combinations

!----------------------------------------------------------------
!+
!  Test of the metric
!+
!----------------------------------------------------------------
subroutine test_metric_i(gcov,gcon,sqrtg,ncheck,nfail,errmax,ncheckg,nfailg,errmaxg,tol)
 use inverse4x4, only:inv4x4
 integer, intent(inout) :: ncheck,nfail,ncheckg,nfailg
 real,    intent(in)    :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg,tol
 real,    intent(inout) :: errmax,errmaxg
 real :: gg(0:3,0:3)
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
 real,    intent(in)    :: x(1:3),dgdx1(0:3,0:3),dgdx2(0:3,0:3),dgdx3(0:3,0:3),tol
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

 real,    intent(in)    :: x(1:3),v(1:3),dens,p,tol
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

!----------------------------------------------------------------
!+
!  Benchmark cons2prim for ideal gas + radiation with entropy
!+
!----------------------------------------------------------------
subroutine benchmark_cons2prim_entropy()
 use cons2primsolver, only:conservative2primitive,primitive2conservative
 use part,            only:ien_entropy_s
 use metric_tools,    only:pack_metric
 use eos,             only:ieos,gmw
 use eos_idealplusrad, only:get_idealplusrad_enfromtemp,get_idealplusrad_pres
 use timing,          only:print_time
 use units,           only:unit_density,unit_pressure,unit_ergg
 use io,              only:id,master
 real :: metrici(0:3,0:3,2)
 real :: x(3),v(3),dens,u,p,temp,gamma,rho,pmom(3),en
 real :: v_out(3),dens_out,u_out,p_out,gamma_out
 real :: rhocgs,ucgs,prescgs,mu,t1,t2
 integer :: ierr,i,nfail,nrepeat
 integer, parameter :: nrepeat_default = 1000

 if (ieos /= 12) return
 if (id /= master) return

 mu = gmw
 nrepeat = nrepeat_default
 nfail = 0

 x    = (/10.,0.1,0.2/)
 v    = (/0.1,0.2,0.05/)
 dens = 1.
 temp = 1.e7
 rhocgs  = dens*unit_density
 call get_idealplusrad_enfromtemp(rhocgs,temp,mu,ucgs)
 call get_idealplusrad_pres(rhocgs,temp,mu,prescgs)
 u = ucgs/unit_ergg
 p = prescgs/unit_pressure
 gamma = 1. + p/(dens*u)

 call pack_metric(x,metrici)
 call primitive2conservative(x,metrici,v,dens,u,p,rho,pmom,en,ien_entropy_s)

 v_out    = v
 dens_out = 0.9*dens
 u_out    = u
 p_out    = 0.9*p
 gamma_out = gamma
 temp     = 0.9*temp

 call cpu_time(t1)
 do i=1,nrepeat
    ierr = 0
    call conservative2primitive(x,metrici,v_out,dens_out,u_out,p_out,temp,gamma_out,&
                               rho,pmom,en,ierr,ien_entropy_s)
    if (ierr /= 0) nfail = nfail + 1
 enddo
 call cpu_time(t2)

 write(*,'(/,a)') '--> benchmark: cons2prim with ien_entropy_s (ieos=12)'
 if (nfail > 0) write(*,'(a,i0,a,i0)') '   convergence failures: ',nfail,' / ',nrepeat
 call print_time(real((t2-t1)/nrepeat,kind=4),'conservative2primitive')

end subroutine benchmark_cons2prim_entropy

!-----------------------------------------------------------------------
!+
!   Test of circular orbits in the Reissner-Nordstrom metric
!   with various charges (Q/M = 0, 0.5, 0.9, 1.1)
!+
!-----------------------------------------------------------------------
subroutine test_rn_charged(ntests,npass)
 use io,           only:id,master
 use metric_tools, only:imetric,imet_rn
 use metric,       only:mass1,charge
 integer, intent(inout) :: ntests,npass
 integer :: nerr(6),itest,iq
 real    :: dt,period,tmax,r0,v0,mass_val,charge_val
 real    :: angmom(3),angmom0(3),xyz(3),vxyz(3)
 real    :: charges_to_test(4)
 real, parameter :: postol = 1.e-5, angtol = 1.e-10

 if (id==master) write(*,'(/,a)') '--> testing RN charged orbits'

 if (imetric /= imet_rn) then
    if (id==master) write(*,'(/,a)') '   Skipping test! Metric is not RN.'
    return
 endif

 ! Q/M ratios to test: Schwarzschild-like, Sub-extremal, Near-extremal, Naked Singularity
 charges_to_test = (/0.0, 0.5, 0.9, 1.1/)
 mass_val = 1.0
 r0 = 10.0

 do iq=1,size(charges_to_test)
    charge_val = charges_to_test(iq)
    charge = charge_val
    mass1  = mass_val

    ! Circular orbit coordinate velocity for RN metric:
    ! v_circ = sqrt( M/r - Q^2/r^2 )
    v0 = sqrt( mass_val/r0 - charge_val**2/r0**2 )

    period = 2.0 * 3.14159265358979 * r0 / v0
    tmax = 1.0 * period
    dt = period / 1000.0

    if (id==master) write(*,'(a,F4.1,a,F4.2)') '   testing Q/M = ', charge_val/mass_val, ' v0 = ', v0

    do itest=1,2
       xyz  = (/r0, 0., 0./)
       vxyz = (/0., v0, 0./)
       call integrate_geodesic(tmax,dt,xyz,vxyz,angmom0,angmom,use_sink=(itest==2))

       nerr = 0
       call checkval(angmom(1),angmom0(1),angtol,nerr(1),'error in angmomx')
       call checkval(angmom(2),angmom0(2),angtol,nerr(2),'error in angmomy')
       call checkval(angmom(3),angmom0(3),angtol,nerr(3),'error in angmomz')
       call checkval(sqrt(dot_product(xyz,xyz)), r0, postol, nerr(4),'error in final r position')

       call update_test_scores(ntests,nerr,npass)
    enddo
 enddo

end subroutine test_rn_charged

!-----------------------------------------------------------------------
!+
!   Test of the line of sight covariant 4-velocity used in the shock
!   viscosity: it should reduce to Gamma*v in flat space and, for a
!   fixed metric, be a monotonic function of the coordinate velocity
!   so that the shock heating is positive
!+
!-----------------------------------------------------------------------
subroutine test_u_los(ntests,npass)
 use physcon,      only:pi
 use utils_gr,     only:get_u_los,get_u0
 use metric_tools, only:get_metric
 use metric,       only:a
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nv = 9
 integer :: nfail(2),i,j,k,l,ia,ib,ierr,nwrong
 real    :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg,rhat(3),pos(3),vel(nv),a_prev
 real    :: ua,ub,u0,errmax
 real, parameter :: radii(3) = (/3.,5.,10./), theta(3) = (/0.3,pi/2.,2.5/)
 real, parameter :: spin(2) = (/0.,0.9/)

 write(*,'(/,a)') '--> testing line of sight 4-velocity for shock viscosity'
 vel = (/(-0.8 + 0.2*(i-1),i=1,nv)/)

 ! flat space: U_i rhat^i = Gamma v
 gcov = 0.
 gcov(0,0) = -1.
 do i=1,3
    gcov(i,i) = 1.
 enddo
 rhat = (/1.,2.,2./)/3.
 errmax = 0.
 do i=1,nv
    errmax = max(errmax,abs(get_u_los(gcov,vel(i),rhat) - vel(i)/sqrt(1. - vel(i)**2)))
 enddo
 call checkval(errmax,0.,1.e-15,nfail(1),'U_los = Gamma*v in flat space')

 ! curved space: (U_los(v_a) - U_los(v_b))*(v_a - v_b) >= 0 for all pairs of velocities
 a_prev = a
 nwrong = 0
 do l=1,size(spin)
    a = spin(l)
    do i=1,size(radii)
       do j=1,size(theta)
          do k=0,3
             pos = radii(i)*(/sin(theta(j))*cos(0.5*pi*k),sin(theta(j))*sin(0.5*pi*k),cos(theta(j))/)
             rhat = (/cos(0.7*k),sin(0.7*k)*cos(theta(j)),sin(0.7*k)*sin(theta(j))/)
             call get_metric(pos,gcov,gcon,sqrtg)
             do ia=1,nv
                call get_u0(gcov,vel(ia)*rhat,u0,ierr)
                if (ierr /= 0) cycle  ! skip velocities faster than light
                ua = get_u_los(gcov,vel(ia),rhat)
                do ib=1,nv
                   call get_u0(gcov,vel(ib)*rhat,u0,ierr)
                   if (ierr /= 0) cycle
                   ub = get_u_los(gcov,vel(ib),rhat)
                   if ((ua - ub)*(vel(ia) - vel(ib)) < 0.) nwrong = nwrong + 1
                enddo
             enddo
          enddo
       enddo
    enddo
 enddo
 a = a_prev
 call checkval(nwrong,0,0,nfail(2),'U_los monotonic in v for fixed metric')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_u_los

!-----------------------------------------------------------------------
!+
!   Test that one explicit step in the GR entropy variable gives the
!   correct change in u when gas is strongly heated or cooled, for cold
!   ideal gas and for radiation-dominated gas+radiation. With a plain
!   Euler step ds = du/(kT) the heating of cold gas is overestimated,
!   overflowing T(s) for large du/u
!+
!-----------------------------------------------------------------------
subroutine test_entropy_step(ntests,npass)
 use eos,     only:get_entropy,get_p_from_rho_s,correct_entropy_derivative,&
                   get_u_from_rhoT,get_pressure,get_temperature,gmw,gamma,init_eos
 use part,    only:ien_entropy,ien_entropy_s
 use physcon, only:kboltz
 use units,   only:unit_density
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nx = 6, neos = 3
 real,    parameter :: xfrac(nx) = (/0.5,10.,1.e6,-0.5,-1.,-2./)
 integer, parameter :: ieos_test(neos) = (/2,12,10/)
 real,    parameter :: rho_cgs(neos) = (/1.,1.e-10,1./)  ! cold gas, radiation dominated, MESA
 real,    parameter :: temp0(neos) = (/10.,1.e5,1.e6/)
 real,    parameter :: tol(neos) = (/1.e-14,1.e-14,0.07/)  ! MESA: consistency of the two tables
 integer :: nerr(nx*(neos+1)),j,k,ieos,ierr
 real    :: rho,u,pres,temp,s,dsdt,dudt,dt,x,xlim,uexact,u1,pres1,temp1,gmw_prev
 real    :: xyz0(3),vxyzu0(4)

 write(*,'(/,a)') '--> testing entropy step for strongly heated or cooled gas'
 gmw_prev = gmw
 gmw    = 0.6
 gamma  = 5./3.
 dt     = 1.
 nerr   = 0
 xyz0   = 0.  ! position and velocity are not used by these equations of state
 vxyzu0 = 0.
 do k=1,neos
    ieos = ieos_test(k)
    call init_eos(ieos,ierr)
    if (ierr /= 0) then
       write(*,'(a,i2,a)') '   skipping ieos=',ieos,' (eos initialisation failed, e.g. tables not found)'
       cycle
    endif
    rho = rho_cgs(k)/unit_density
    u   = get_u_from_rhoT(rho,temp0(k),ieos)
    vxyzu0(4) = u
    pres = get_pressure(ieos,xyz0,rho,vxyzu0)
    temp = get_temperature(ieos,xyz0,rho,vxyzu0)
    s    = get_entropy(rho,pres,gmw,ieos)
    do j=1,nx
       x    = xfrac(j)
       if (ieos==10 .and. x > 100.) cycle  ! would go off the MESA tables
       dudt = x*u/dt
       ! cooling that would make u negative is limited as in the non-GR code
       xlim = x
       if (x <= -1.) xlim = x/(1. - x)
       uexact = u*(1. + xlim)

       ! s entropy: Euler step computed as in force, corrected, then T(s) and u(T) from the eos
       dsdt = dudt/(kboltz*temp)
       dsdt = correct_entropy_derivative(ien_entropy_s,x,dsdt,dt,ieos,xyz0,rho,u,pres)
       call get_p_from_rho_s(ieos,s + dt*dsdt,rho,gmw,pres1,temp1)
       u1   = get_u_from_rhoT(rho,temp1,ieos)
       call checkval(u1/uexact,1.,tol(k),nerr((k-1)*nx+j),'u/u_exact after s entropy step')
    enddo
 enddo

 ! K entropy: u is proportional to K at fixed density
 do j=1,nx
    x    = xfrac(j)
    xlim = x
    if (x <= -1.) xlim = x/(1. - x)
    u1 = u + dt*correct_entropy_derivative(ien_entropy,x,x*u/dt,dt,2,xyz0,rho,u,pres)
    call checkval(u1/(u*(1. + xlim)),1.,1.e-15,nerr(neos*nx+j),'u/u_exact after K entropy step')
 enddo
 gmw = gmw_prev

 call update_test_scores(ntests,nerr,npass)

end subroutine test_entropy_step

!-----------------------------------------------------------------------
!+
!   Test that shock heating never decreases the entropy of a particle,
!   using a low resolution polytrope passing close to the black hole.
!   Tidal distortion of the star triggers the shock viscosity. With
!   conductivity switched off, K can only increase
!+
!-----------------------------------------------------------------------
subroutine test_tde_entropy(ntests,npass)
 use io,             only:id,master,iverbose,real4
 use dim,            only:maxp,maxalpha
 use part,           only:init_part,npart,npartoftype,xyzh,vxyzu,pxyzu,fext,dens,metrics,&
                          metricderivs,eos_vars,rad,massoftype,hfact,alphaind,xyzmh_ptmass,&
                          vxyz_ptmass,nptmass,nsinkproperties,ien_type,ien_entropy,&
                          isdead_or_accreted,igas
 use mpidomain,      only:i_belong
 use eos,            only:ieos,gamma,X_in,Z_in,ishock_heating
 use options,        only:alpha,alphau,beta,ireconav,iexternalforce
 use setstar,        only:star_t,set_star,set_defaults_star,ipoly
 use setorbit,       only:orbit_t,set_defaults_orbit,set_orbit
 use externalforces, only:initialise_externalforces
 use extern_gr,      only:get_grforce_all
 use cons2prim,      only:prim2consall
 use metric_tools,   only:init_metric
 use deriv,          only:get_derivs_global,get_density_global
 use step_lf_global, only:init_step,step
 use timestep,       only:dtmax,dtcourant,dtforce
 integer, intent(inout) :: ntests,npass
 type(star_t)  :: star
 type(orbit_t) :: orbit
 integer :: nfail(3),ierr,i,nsteps,ndecrease
 integer(kind=8) :: ntot
 real    :: xyzmh_in(nsinkproperties,2),vxyz_in(3,2),xstar(3),vstar(3)
 real    :: t,tmax,dt,dtnew,dtextforce,rhozero,kprev(maxp)

 if (id==master) write(*,'(/,a)') '--> testing entropy does not decrease in a tidal disruption'
 iverbose = 0

 ! polytropic star, with no relaxation since the test only needs tidal distortion
 call init_part()
 ieos  = 2
 gamma = 5./3.
 call set_defaults_star(star)
 star%m = '1*msun'
 star%r = '1*rsun'
 star%iprofile = ipoly
 star%np = 1000
 nptmass = 0
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,relax=.false.,&
               use_var_comp=.false.,write_rho_to_file=.false.,rhozero=rhozero,&
               npart_total=ntot,mask=i_belong,ierr=ierr,write_files=.false.)
 call checkval(ierr,0,0,nfail(1),'set_star runs with ierr = 0')

 ! eccentric orbit with pericentre at 10 M, starting just before pericentre
 call set_defaults_orbit(orbit)
 orbit%elems%a = '50.'
 orbit%e = 0.8
 orbit%f = -60.
 xyzmh_in = 0.
 vxyz_in  = 0.
 i = 0
 call set_orbit(orbit,1.,npart*massoftype(igas),1.e-8,1.e-8,xyzmh_in,vxyz_in,i,.false.,ierr)
 call checkval(ierr,0,0,nfail(2),'set_orbit runs with ierr = 0')
 call update_test_scores(ntests,nfail(1:2),npass)
 xstar = xyzmh_in(1:3,2) - xyzmh_in(1:3,1)
 vstar = vxyz_in(1:3,2) - vxyz_in(1:3,1)
 do i=1,npart
    xyzh(1:3,i)  = xyzh(1:3,i) + xstar
    vxyzu(1:3,i) = vstar
 enddo

 ! evolve K with shock viscosity but no conductivity
 ien_type       = ien_entropy
 ishock_heating = 1
 alpha          = 1.
 alphau         = 0.
 beta           = 2.
 ireconav       = -1
 if (maxalpha==maxp) alphaind(1,:) = real4(alpha)
 iexternalforce = 1

 t = 0.
 tmax  = 200.
 dtmax = tmax
 call get_density_global(2,zero_fxyzu=.true.)
 call init_metric(npart,xyzh,metrics,metricderivs)
 call prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,dens=dens)
 call initialise_externalforces(iexternalforce,ierr)
 call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,fext,dtextforce,dens=dens)
 call get_derivs_global()
 call init_step(npart,t,dtmax)

 ! count particles whose entropy decreases during any step
 ndecrease = 0
 nsteps = 0
 dt = min(dtcourant,dtforce,dtextforce)
 do while (t < tmax)
    kprev(1:npart) = pxyzu(4,1:npart)
    t = t + dt
    call step(npart,npart,t,dt,dtextforce,dtnew)
    nsteps = nsteps + 1
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i)) .and. pxyzu(4,i) < kprev(i)) ndecrease = ndecrease + 1
    enddo
    dt = min(dtnew,tmax-t)
 enddo
 if (id==master) print*,' nsteps = ',nsteps
 call checkval(ndecrease,0,0,nfail(3),'number of times K decreased')
 call update_test_scores(ntests,nfail(3:3),npass)

end subroutine test_tde_entropy

!-----------------------------------------------------------------------
!+
!  forward/backward consistency test between the two sets of MESA EoS
!  tables for the same composition (X_in,Z_in):
!
!   - classical tables, tabulated in (rho, u): give P, T and s
!   - GR tables, tabulated in (rho, s): give P, T and u
!
!  Forward:  (rho,u) -> s from classical tables -> u,P,T from GR tables
!  Backward: (rho,s) -> u from GR tables -> s,P,T from classical tables
!
!  Also checks the round trip through the code-unit routines in the eos
!  module used by the GR cons2prim solver (get_entropy, get_p_from_rho_s,
!  get_u_from_rho_s). Points where either table would extrapolate are
!  skipped.
!+
!-----------------------------------------------------------------------
subroutine test_mesa_forward_backward(ntests,npass)
 use io,                only:id,master
 use eos,               only:ieos,init_eos,X_in,Z_in,gmw,get_entropy,get_p_from_rho_s,get_u_from_rho_s
 use mesa_microphysics, only:getvalue_mesa,getvalue_mesa_gr
 use units,             only:unit_density,unit_pressure,unit_ergg
 use table_utils,       only:linspace
 integer, intent(inout) :: ntests,npass
 integer, parameter     :: ncheck_min = 50
 integer                :: nfail(8),ncheck(8),i,j,ierr,ierr1,ierr2,ieos_prev
 real                   :: errmax(8),tol_u,tol_s,tol_p,tol_t
 real                   :: rhoi,ui,si,p1,t1,u2,p2,t2,s2
 real                   :: rho_code,p_code,s_code,u_code,temp
 real                   :: logrho(16),logu(16),logs(16)

 if (id==master) write(*,"(/,a)") '--> testing forward/backward consistency of MESA EoS tables'

 ieos_prev = ieos
 ieos = 10
 call init_eos(ieos,ierr)
 if (ierr /= 0) then
    if (id==master) write(*,"(a)") ' skipped: MESA EoS tables not available'
    ieos = ieos_prev
    return
 endif
 if (id==master) write(*,"(1x,a,f5.3,a,f5.3)") 'using X = ',X_in,', Z = ',Z_in

 !
 ! tolerances on relative errors, set by the interpolation error between
 ! the two tables (different grids: GR table spacing is 0.27 dex in rho and
 ! 0.012 dex in s; classical table is 0.059 dex in u and 0.2 dex in V).
 ! Measured maxima on this grid are 1.5% in s, 6% in u (an error ds/s is
 ! amplified to du/u ~ s/c_v ds/s ~ 15 ds/s) and 39% in P and T, the latter
 ! in H/He partial ionisation zones where T(u) varies sharply between nodes
 ! (median errors are 0.1% in s, 1.5% in u, 0.05% in P and T)
 !
 tol_s = 0.025
 tol_u = 0.1
 tol_p = 0.5
 tol_t = 0.5

 ! log-spaced grids in cgs units: rho in g/cm^3, u in erg/g, s in erg/g/K
 ! (getvalue_mesa returns s in erg/g/K, but getvalue_mesa_gr expects s*kB*N_A,
 !  the units returned by eos:entropy with ientropy=3)
 ! (rho <= 1 g/cm^3 avoids cold degenerate states below the zero-temperature
 !  energy, where the classical table is not monotonic in u)
 call linspace(logrho,-8.,0.)
 call linspace(logu,12.,17.)
 call linspace(logs,9.,10.5)

 nfail = 0; ncheck = 0; errmax = 0.
 do i=1,size(logrho)
    rhoi = 10.**logrho(i)
    !
    ! forward: (rho,u) -> s from classical tables -> u from GR tables
    !
    do j=1,size(logu)
       ui = 10.**logu(j)
       call getvalue_mesa(rhoi,ui,9,si,ierr1)
       call getvalue_mesa_gr(rhoi,si,2,u2,ierr2)
       if (ierr1 /= 0 .or. ierr2 /= 0) cycle  ! skip points off either table
       call getvalue_mesa(rhoi,ui,2,p1)
       call getvalue_mesa(rhoi,ui,4,t1)
       call getvalue_mesa_gr(rhoi,si,1,p2)
       call getvalue_mesa_gr(rhoi,si,3,t2)
       call checkvalbuf(u2,ui,tol_u,'u(rho,s(rho,u))',nfail(1),ncheck(1),errmax(1),use_rel_tol=.true.)
       call checkvalbuf(p2,p1,tol_p,'P(rho,s(rho,u))',nfail(2),ncheck(2),errmax(2),use_rel_tol=.true.)
       call checkvalbuf(t2,t1,tol_t,'T(rho,s(rho,u))',nfail(3),ncheck(3),errmax(3),use_rel_tol=.true.)
       !
       ! same round trip via the code-unit routines in the eos module
       !
       rho_code = rhoi/unit_density
       p_code   = p1/unit_pressure
       s_code   = get_entropy(rho_code,p_code,gmw,ieos)
       temp     = t1
       call get_p_from_rho_s(ieos,s_code,rho_code,gmw,p2,temp)
       call get_u_from_rho_s(ieos,s_code,rho_code,u_code)
       call checkvalbuf(u_code*unit_ergg,ui,tol_u,'u(rho,S) via eos module',nfail(7),ncheck(7),errmax(7),use_rel_tol=.true.)
       call checkvalbuf(p2,p_code,tol_p,'P(rho,S) via eos module',nfail(8),ncheck(8),errmax(8),use_rel_tol=.true.)
    enddo
    !
    ! backward: (rho,s) -> u from GR tables -> s from classical tables
    !
    do j=1,size(logs)
       si = 10.**logs(j)
       call getvalue_mesa_gr(rhoi,si,2,ui,ierr2)
       call getvalue_mesa(rhoi,ui,9,s2,ierr1)
       if (ierr1 /= 0 .or. ierr2 /= 0) cycle  ! skip points off either table
       call getvalue_mesa_gr(rhoi,si,1,p2)
       call getvalue_mesa_gr(rhoi,si,3,t2)
       call getvalue_mesa(rhoi,ui,2,p1)
       call getvalue_mesa(rhoi,ui,4,t1)
       call checkvalbuf(s2,si,tol_s,'s(rho,u(rho,s))',nfail(4),ncheck(4),errmax(4),use_rel_tol=.true.)
       call checkvalbuf(p1,p2,tol_p,'P(rho,u(rho,s))',nfail(5),ncheck(5),errmax(5),use_rel_tol=.true.)
       call checkvalbuf(t1,t2,tol_t,'T(rho,u(rho,s))',nfail(6),ncheck(6),errmax(6),use_rel_tol=.true.)
    enddo
 enddo
 call checkvalbuf_end('forward:  u(rho,s(rho,u)) = u',ncheck(1),nfail(1),errmax(1),tol_u)
 call checkvalbuf_end('forward:  P(rho,s(rho,u)) = P',ncheck(2),nfail(2),errmax(2),tol_p)
 call checkvalbuf_end('forward:  T(rho,s(rho,u)) = T',ncheck(3),nfail(3),errmax(3),tol_t)
 call checkvalbuf_end('backward: s(rho,u(rho,s)) = s',ncheck(4),nfail(4),errmax(4),tol_s)
 call checkvalbuf_end('backward: P(rho,u(rho,s)) = P',ncheck(5),nfail(5),errmax(5),tol_p)
 call checkvalbuf_end('backward: T(rho,u(rho,s)) = T',ncheck(6),nfail(6),errmax(6),tol_t)
 call checkvalbuf_end('eos module: u(rho,S(rho,P)) = u',ncheck(7),nfail(7),errmax(7),tol_u)
 call checkvalbuf_end('eos module: P(rho,S(rho,P)) = P',ncheck(8),nfail(8),errmax(8),tol_p)
 call update_test_scores(ntests,nfail,npass)

 ! make sure that a reasonable number of points were on both tables
 nfail = 0
 call checkval(min(ncheck(1),ncheck(4)) >= ncheck_min,.true.,nfail(1),'enough points on both tables')
 call update_test_scores(ntests,nfail(1:1),npass)
 ieos = ieos_prev

end subroutine test_mesa_forward_backward

end module testgr
