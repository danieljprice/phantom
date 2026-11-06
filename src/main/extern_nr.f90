!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_nr
!
! Coupling of the GRSPH particles to a metric evolved by the numerical
! relativity code (libnr, through nr_api), with METRIC=et and NR=yes.
!
! The metric lives on the periodic, vertex-centred grid of metric_et_utils,
! with ngrid_evolve nodes per dimension spanning the particle domain. The NR
! code advances it whenever Phantom updates the metric to a later time
! (metric_et:update_metric). At every Runge-Kutta stage the NR code asks for
! the stress-energy tensor for the metric of that stage: the metric is
! interpolated from the stage grid to the particles, the primitive variables
! are recovered from the conserved ones, and T_munu is interpolated from the
! particles to the grid nodes (as in Magnall et al. 2023), normalised so the
! rest mass on the grid equals that of the particles. This uses private copies
! of the particle metric and primitives, so the particle arrays are untouched.
!
! The initial slice is flat with unit lapse and zero shift, and the extrinsic
! curvature K_ij = -sqrt(8 pi rho_E/3) delta_ij solves the Hamiltonian
! constraint at each node: initial data for an expanding, (nearly)
! homogeneous universe, as in the FLRW tests.
!
! :References: Magnall et al. (2023), Phys. Rev. D 108, 103534
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, cons2prim, extern_gr, interpolations3D, io,
!   metric_et_utils, metric_tools, nr_api, part, physcon
!
 use nr_api, only:nr_dp
 implicit none

 public :: init_extern_nr

 private

 ! private copies of the particle metric, primitives and stress-energy tensor
 real, allocatable :: metrics_nr(:,:,:,:),metricderivs_nr(:,:,:,:),vxyzu_nr(:,:)
 real, allocatable :: dens_nr(:),eos_vars_nr(:,:),tmunus_nr(:,:,:)
 ! metric on the grid, in the real kind of the NR interface
 real(nr_dp), allocatable :: gcov_nr(:,:,:,:,:),gcon_nr(:,:,:,:,:),sqrtg_nr(:,:,:)
 real(nr_dp), allocatable :: dgcov_nr(:,:,:,:,:,:)
 ! T_munu on the grid from the last call to get_tmunu_grid
 real(nr_dp), allocatable :: tmunu_last(:,:,:,:,:)
 ! file for the time evolution of the metric
 integer :: iunit_nr = 0

contains

!-----------------------------------------------------------------------
!+
!  set up the grid and the initial data, and start evolving the metric
!+
!-----------------------------------------------------------------------
subroutine init_extern_nr(time)
 use io,              only:fatal,iprint,id,master,fileprefix
 use nr_api,          only:nr_linked,nr_init,nr_formulation
 use boundary,        only:xmin,ymin,zmin,dxbound,dybound,dzbound
 use metric_et_utils, only:allocate_grid,ngrid_evolve,periodic_grid,advance_metric
 use physcon,         only:pi
 real, intent(in) :: time
 integer :: n(3),i,ierr
 real(nr_dp) :: dx(3),origin(3)
 real(nr_dp), allocatable :: alpha(:,:,:),beta(:,:,:,:),g3(:,:,:,:,:),kij(:,:,:,:,:),tmunu(:,:,:,:,:)

 if (.not.nr_linked()) call fatal('extern_nr','evolve_metric = T requires Phantom compiled with NR=yes')
 n = ngrid_evolve
 dx = [dxbound,dybound,dzbound]/n
 origin = [xmin,ymin,zmin]
 call allocate_grid(n(1),n(2),n(3),dx(1),dx(2),dx(3),origin(1),origin(2),origin(3))
 periodic_grid = .true.
 allocate(gcov_nr(0:3,0:3,n(1),n(2),n(3)),gcon_nr(0:3,0:3,n(1),n(2),n(3)),sqrtg_nr(n(1),n(2),n(3)))
 allocate(dgcov_nr(0:3,0:3,3,n(1),n(2),n(3)),tmunu(0:3,0:3,n(1),n(2),n(3)))
 allocate(alpha(n(1),n(2),n(3)),beta(3,n(1),n(2),n(3)))
 allocate(g3(3,3,n(1),n(2),n(3)),kij(3,3,n(1),n(2),n(3)))

 ! flat initial slice with unit lapse and zero shift
 gcov_nr = 0.
 gcov_nr(0,0,:,:,:) = -1.
 alpha = 1.
 beta  = 0.
 g3    = 0.
 kij   = 0.
 do i=1,3
    gcov_nr(i,i,:,:,:) = 1.
    g3(i,i,:,:,:) = 1.
 enddo
 dgcov_nr = 0.
 ! extrinsic curvature from the Hamiltonian constraint, with energy density T_00
 call get_tmunu_grid(real(time,kind=nr_dp),gcov_nr,dgcov_nr,tmunu)
 do i=1,3
    kij(i,i,:,:,:) = -sqrt(8.*pi*tmunu(0,0,:,:,:)/3.)
 enddo

 call nr_init(n,origin,dx,real(time,kind=nr_dp),alpha,beta,g3,kij,get_tmunu_grid,ierr)
 if (ierr /= 0) call fatal('extern_nr','could not initialise the NR evolution')
 call load_nr_metric()
 advance_metric => advance_nr_metric

 if (id==master) then
    write(iprint,"(a,a,a,i4,a)") ' evolving the metric with NR (',trim(nr_formulation()),') on a ',&
                                 ngrid_evolve,'^3 periodic grid'
    open(newunit=iunit_nr,file=trim(fileprefix)//'_nr.ev',status='replace',form='formatted')
    write(iunit_nr,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") 1,'time',2,'nsteps',3,'a',4,'alpha',&
          5,'rho',6,'rho_E',7,'constraints'
 endif
 call write_nr_ev(real(time,kind=nr_dp),0)

end subroutine init_extern_nr

!-----------------------------------------------------------------------
!+
!  advance the metric to a later time (called by metric_et:update_metric)
!+
!-----------------------------------------------------------------------
subroutine advance_nr_metric(time)
 use nr_api, only:nr_time,nr_evolve_to
 real, intent(in) :: time
 integer :: nsteps

 ! ignore requests that differ from the current time only by round-off
 if (time - nr_time() > 1.e-12*abs(time)) then
    call nr_evolve_to(real(time,kind=nr_dp),nsteps)
    call load_nr_metric()
    call write_nr_ev(nr_time(),nsteps)
 endif

end subroutine advance_nr_metric

!-----------------------------------------------------------------------
!+
!  copy the current metric of the NR evolution to the grid
!+
!-----------------------------------------------------------------------
subroutine load_nr_metric()
 use nr_api,          only:nr_get_metric
 use metric_et_utils, only:gcovgrid,gcongrid,sqrtggrid,metricderivsgrid,gridinit

 call nr_get_metric(gcov_nr,gcon_nr,sqrtg_nr,dgcov_nr)
 gcovgrid  = real(gcov_nr)
 gcongrid  = real(gcon_nr)
 sqrtggrid = real(sqrtg_nr)
 metricderivsgrid = real(dgcov_nr)
 gridinit = .true.

end subroutine load_nr_metric

!-----------------------------------------------------------------------
!+
!  stress-energy tensor T_munu of the particles at the grid nodes, for
!  the metric gcov (with spatial derivatives metricderivs) on the grid
!+
!-----------------------------------------------------------------------
subroutine get_tmunu_grid(t,gcov,metricderivs,tmunu)
 use part,            only:npart,xyzh,pxyzu,vxyzu,dens
 use metric_et_utils, only:gcovgrid,metricderivsgrid,gridinit
 use metric_tools,    only:init_metric
 use cons2prim,       only:cons2primall
 use extern_gr,       only:get_tmunu_all
 real(nr_dp), intent(in)  :: t,gcov(0:,0:,:,:,:),metricderivs(0:,0:,:,:,:,:)
 real(nr_dp), intent(out) :: tmunu(0:,0:,:,:,:)

 ! metric at the particles, interpolated from the grid
 gcovgrid = real(gcov)
 metricderivsgrid = real(metricderivs)
 gridinit = .true.
 call allocate_work(npart)
 vxyzu_nr = vxyzu(:,1:npart)
 dens_nr  = dens(1:npart)
 call init_metric(npart,xyzh,metrics_nr,metricderivs_nr)

 ! primitive variables and T_munu on the particles
 call cons2primall(npart,xyzh,metrics_nr,pxyzu,vxyzu_nr,dens_nr,eos_vars_nr)
 call get_tmunu_all(npart,xyzh,metrics_nr,vxyzu_nr,metricderivs_nr,dens_nr,tmunus_nr)

 call deposit_tmunu(npart,tmunus_nr,tmunu)
 if (.not.allocated(tmunu_last)) allocate(tmunu_last,mold=tmunu)
 tmunu_last = tmunu

end subroutine get_tmunu_grid

!-----------------------------------------------------------------------
!+
!  interpolate T_munu from the particles to the grid nodes, using the
!  normalised SPH interpolant with weights m/(rho* h^3), and scale it so
!  that the rest mass on the grid, integral of rho*, equals that of the
!  particles
!+
!-----------------------------------------------------------------------
subroutine deposit_tmunu(npart,tmunus,tmunu)
 use part,             only:xyzh,rho,massoftype,igas,isdead_or_accreted
 use metric_et_utils,  only:gridsize,dxgrid,gridorigin
 use interpolations3D, only:interpolate3D_vecexact
 integer,     intent(in)  :: npart
 real,        intent(in)  :: tmunus(:,:,:)
 real(nr_dp), intent(out) :: tmunu(0:,0:,:,:,:)
 integer, parameter :: ndat = 17
 real, allocatable :: dat(:,:),weights(:),datsmooth(:,:,:,:)
 integer, allocatable :: itype(:)
 integer :: i,mu,nu,nalive
 real :: cfac

 allocate(dat(npart,ndat),weights(npart),itype(npart))
 allocate(datsmooth(ndat,gridsize(1),gridsize(2),gridsize(3)))
 nalive = 0
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       itype(i) = -1
       weights(i) = 0.
       dat(i,:) = 0.
    else
       nalive = nalive + 1
       itype(i) = igas
       weights(i) = massoftype(igas)/(rho(i)*xyzh(4,i)**3)
       do nu=0,3
          do mu=0,3
             dat(i,1 + mu + 4*nu) = tmunus(mu+1,nu+1,i)
          enddo
       enddo
       dat(i,ndat) = rho(i)
    endif
 enddo

 ! pixel i of the interpolation is at xmin + i*dx, so shifting the origin
 ! by one cell puts pixel i on grid node i
 call interpolate3D_vecexact(xyzh,weights,dat,ndat,itype,npart,gridorigin(1)-dxgrid(1),&
                             gridorigin(2)-dxgrid(2),gridorigin(3)-dxgrid(3),datsmooth,&
                             gridsize(1),gridsize(2),gridsize(3),dxgrid(1),dxgrid(2),dxgrid(3),&
                             .true.,.true.,.true.,.true.)

 cfac = nalive*massoftype(igas)/(sum(datsmooth(ndat,:,:,:))*product(dxgrid))
 do nu=0,3
    do mu=0,3
       tmunu(mu,nu,:,:,:) = cfac*datsmooth(1 + mu + 4*nu,:,:,:)
    enddo
 enddo

end subroutine deposit_tmunu

!-----------------------------------------------------------------------
!+
!  allocate the private copies of the particle arrays
!+
!-----------------------------------------------------------------------
subroutine allocate_work(npart)
 use part, only:vxyzu,eos_vars
 integer, intent(in) :: npart

 if (allocated(dens_nr)) then
    if (size(dens_nr) == npart) return
    deallocate(metrics_nr,metricderivs_nr,vxyzu_nr,dens_nr,eos_vars_nr,tmunus_nr)
 endif
 allocate(metrics_nr(0:3,0:3,2,npart),metricderivs_nr(0:3,0:3,3,npart))
 allocate(vxyzu_nr(size(vxyzu,1),npart),dens_nr(npart),eos_vars_nr(size(eos_vars,1),npart))
 allocate(tmunus_nr(0:3,0:3,npart))

end subroutine allocate_work

!-----------------------------------------------------------------------
!+
!  write the grid-averaged scale factor a = det(gamma)^(1/6), lapse,
!  rest-mass density of the particles, energy density
!  rho_E = n^mu n^nu T_munu and L2 norm of the constraints to the _nr.ev file
!+
!-----------------------------------------------------------------------
subroutine write_nr_ev(t,nsteps)
 use io,     only:id,master
 use nr_api, only:nr_constraints
 real(nr_dp), intent(in) :: t
 integer,     intent(in) :: nsteps
 real(nr_dp) :: cnorm,abar,alphabar,rhoE,alpha
 integer :: i,j,k,mu,nu,nnodes

 ! constraints with the source for the current metric (sets tmunu_last and dens_nr)
 call nr_constraints(cnorm)
 abar = 0.
 alphabar = 0.
 rhoE = 0.
 nnodes = size(sqrtg_nr)
 do k=1,size(sqrtg_nr,3)
    do j=1,size(sqrtg_nr,2)
       do i=1,size(sqrtg_nr,1)
          ! g^00 = -1/alpha^2, sqrt(-g) = alpha sqrt(gamma) and n^mu = -alpha g^(mu 0)
          alpha = 1./sqrt(-gcon_nr(0,0,i,j,k))
          alphabar = alphabar + alpha
          abar = abar + (sqrtg_nr(i,j,k)/alpha)**(1./3.)
          do nu=0,3
             do mu=0,3
                rhoE = rhoE + alpha**2*gcon_nr(mu,0,i,j,k)*gcon_nr(nu,0,i,j,k)*tmunu_last(mu,nu,i,j,k)
             enddo
          enddo
       enddo
    enddo
 enddo
 if (id==master) then
    write(iunit_nr,"(es18.10,1x,i10,5(1x,es18.10))") t,nsteps,abar/nnodes,alphabar/nnodes,&
                    sum(dens_nr)/size(dens_nr),rhoE/nnodes,cnorm
    flush(iunit_nr)
 endif

end subroutine write_nr_ev

end module extern_nr
