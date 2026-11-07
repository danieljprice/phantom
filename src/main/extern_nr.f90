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
! Initial data are those of a flat FLRW universe at a = 1 with a linear
! scalar perturbation in the longitudinal gauge (e.g. Macpherson et al. 2017),
!
!   ds^2 = -(1 + 2 phi) dt^2 + (1 - 2 phi) dx^i dx_i,  phi' = 0,
!
! so that alpha = sqrt(1 + 2 phi), beta = 0, gamma_ij = (1 - 2 phi) delta_ij and
! K_ij = -H (1 - 2 phi)/sqrt(1 + 2 phi) delta_ij, with the conformal Hubble rate
! H = sqrt(8 pi rho_E/3) from the mean energy density. For dust, phi follows
! from the conserved density of the particles, rho* = rho sqrt(gamma), through
! the linearised Hamiltonian constraint
!
!   (del^2 - 3 H^2 - 12 pi rho) phi = 4 pi rho delta*,   delta* = rho*/<rho*> - 1,
!
! solved spectrally on the periodic grid (phi = 0 for a uniform universe). The
! particle velocities must satisfy the momentum constraint, v = -H grad phi/(4 pi rho),
! as in setup_flrw.
!
! :References: Magnall et al. (2023), Phys. Rev. D 108, 103534
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, cons2prim, einsteintk_utils, extern_gr,
!   interpolations3D, io, metric_et_utils, metric_tools, nr_api, part,
!   physcon, timestep, utils_dft
!
 use nr_api, only:nr_dp
 implicit none

 public :: init_extern_nr

 private

 ! private copies of the particle metric, primitives and stress-energy tensor
 real, allocatable :: metrics_nr(:,:,:,:),metricderivs_nr(:,:,:,:),vxyzu_nr(:,:)
 real, allocatable :: dens_nr(:),eos_vars_nr(:,:),tmunus_nr(:,:,:),rhosave_nr(:)
 real, allocatable :: xyzh_nr(:,:),pxyzu_nr(:,:)
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
 use metric_et_utils, only:allocate_grid,ngrid_evolve,periodic_grid,advance_metric,nr_solve_init
 real, intent(in) :: time
 integer :: n(3),i,ierr
 real(nr_dp) :: dx(3),origin(3),hubble
 real(nr_dp), allocatable :: alpha(:,:,:),beta(:,:,:,:),g3(:,:,:,:,:),kij(:,:,:,:,:),tmunu(:,:,:,:,:)
 real(nr_dp), allocatable :: rhostar(:,:,:)
 real, allocatable :: phi(:,:,:)

 if (.not.nr_linked()) call fatal('extern_nr','evolve_metric = T requires Phantom compiled with NR=yes')
 n = ngrid_evolve
 dx = [dxbound,dybound,dzbound]/n
 origin = [xmin,ymin,zmin]
 call allocate_grid(n(1),n(2),n(3),real(dx(1)),real(dx(2)),real(dx(3)),&
                    real(origin(1)),real(origin(2)),real(origin(3)))
 periodic_grid = .true.
 allocate(gcov_nr(0:3,0:3,n(1),n(2),n(3)),gcon_nr(0:3,0:3,n(1),n(2),n(3)),sqrtg_nr(n(1),n(2),n(3)))
 allocate(dgcov_nr(0:3,0:3,3,n(1),n(2),n(3)),tmunu(0:3,0:3,n(1),n(2),n(3)))
 allocate(alpha(n(1),n(2),n(3)),beta(3,n(1),n(2),n(3)))
 allocate(g3(3,3,n(1),n(2),n(3)),kij(3,3,n(1),n(2),n(3)))
 allocate(rhostar(n(1),n(2),n(3)),phi(n(1),n(2),n(3)))

 ! T_munu and the conserved density deposited on the flat slice
 gcov_nr = 0.
 gcov_nr(0,0,:,:,:) = -1.
 do i=1,3
    gcov_nr(i,i,:,:,:) = 1.
 enddo
 dgcov_nr = 0.
 call get_tmunu_rhostar_grid(real(time,kind=nr_dp),gcov_nr,dgcov_nr,tmunu,rhostar)

 if (nr_solve_init) then
    call solved_initial_slice(time,n,origin,dx,tmunu,rhostar,g3,kij,hubble)
 else
    call linear_initial_slice(dx,tmunu,rhostar,g3,kij,hubble)
 endif
 ! lapse of the longitudinal gauge, alpha^2 = 1 + 2 phi with gamma_xx = 1 - 2 phi
 phi = real(0.5*(1. - g3(1,1,:,:,:)))
 alpha = sqrt(1. + 2.*phi)
 beta  = 0.

 call nr_init(n,origin,dx,real(time,kind=nr_dp),alpha,beta,g3,kij,get_tmunu_grid,ierr)
 if (ierr /= 0) call fatal('extern_nr','could not initialise the NR evolution')
 call load_nr_metric()
 advance_metric => advance_nr_metric

 if (id==master) then
    write(iprint,"(a,a,a,i4,a)") ' evolving the metric with NR (',trim(nr_formulation()),') on a ',&
                                 ngrid_evolve,'^3 periodic grid'
    write(iprint,"(a,es12.5,a,es10.3)") ' initial conformal Hubble rate = ',hubble,&
                                        ', max |phi| = ',maxval(abs(phi))
    open(newunit=iunit_nr,file=trim(fileprefix)//'_nr.ev',status='replace',form='formatted')
    write(iunit_nr,"('#',15(1x,'[',i2.2,1x,a11,']',2x))") 1,'time',2,'nsteps',3,'a',4,'alpha',&
          5,'rho',6,'rho_E',7,'constraints',8,'phi amp',9,'delta amp',10,'ham L1',11,'rho max',&
          12,'mom L1',13,'8 pi |S|',14,'rms delta*',15,'min rho*'
 endif
 call write_nr_ev(real(time,kind=nr_dp),0)

end subroutine init_extern_nr

!-----------------------------------------------------------------------
!+
!  initial slice satisfying the full (discrete) Einstein constraints for
!  the deposited particles: the conserved density D = rho* and momentum
!  density D u_i = -T_0i on the flat slice do not depend on the metric,
!  and the NR code solves for a conformally flat metric on the slice of
!  the longitudinal gauge (uniform alpha K, alpha^2 gamma_xx = 1). The
!  matter source in each iteration is T_munu deposited from the particles
!  in the trial metric, as the evolution will compute it
!+
!-----------------------------------------------------------------------
subroutine solved_initial_slice(time,n,origin,dx,tmunu,rhostar,g3,kij,hubble)
 use io,     only:fatal,iprint,id,master
 use nr_api, only:nr_solve_constraints
 real,        intent(in)  :: time
 integer,     intent(in)  :: n(3)
 real(nr_dp), intent(in)  :: origin(3),dx(3),tmunu(0:,0:,:,:,:),rhostar(:,:,:)
 real(nr_dp), intent(out) :: g3(:,:,:,:,:),kij(:,:,:,:,:),hubble
 real(nr_dp), allocatable :: mom(:,:,:,:),alpha(:,:,:)
 real(nr_dp) :: kmean,hrel,mrel
 integer :: i,niter,ierr

 allocate(mom(3,n(1),n(2),n(3)),alpha(n(1),n(2),n(3)))
 do i=1,3
    mom(i,:,:,:) = -tmunu(0,i,:,:,:)
 enddo
 call nr_solve_constraints(n,origin,dx,rhostar,mom,alpha,g3,kij,kmean,niter,hrel,mrel,ierr,&
                           real(time,kind=nr_dp),get_tmunu_grid)
 if (id==master) write(iprint,"(a,i0,a,es10.3,a,es10.3)") ' initial data: constraints solved in ',niter,&
    ' iterations, H/16 pi E = ',hrel,', |M|/8 pi |P| = ',mrel
 if (ierr /= 0) call fatal('extern_nr','constraint solve for the initial data did not converge')
 hubble = -kmean/3.

end subroutine solved_initial_slice

!-----------------------------------------------------------------------
!+
!  initial slice of the linearised longitudinal gauge: the potential phi
!  from the linearised Hamiltonian constraint for the deposited rho*,
!  (del^2 - 3H^2 - 12 pi rho*) phi = 4 pi (rho* - <rho*>)
!+
!-----------------------------------------------------------------------
subroutine linear_initial_slice(dx,tmunu,rhostar,g3,kij,hubble)
 use physcon,   only:pi
 use utils_dft, only:solve_helmholtz_periodic
 real(nr_dp), intent(in)  :: dx(3),tmunu(0:,0:,:,:,:),rhostar(:,:,:)
 real(nr_dp), intent(out) :: g3(:,:,:,:,:),kij(:,:,:,:,:),hubble
 real, allocatable :: phi(:,:,:)
 real(nr_dp) :: rhobar
 integer :: i

 hubble = sqrt(8.*pi*sum(tmunu(0,0,:,:,:))/size(rhostar)/3.)
 rhobar = sum(rhostar)/size(rhostar)
 allocate(phi,mold=real(rhostar))
 call solve_helmholtz_periodic(real(4.*pi*(rhostar - rhobar)),real(dx),real(3.*hubble**2 + 12.*pi*rhobar),phi)
 g3  = 0.
 kij = 0.
 do i=1,3
    g3(i,i,:,:,:)  = 1. - 2.*phi
    kij(i,i,:,:,:) = -hubble*(1. - 2.*phi)/sqrt(1. + 2.*phi)
 enddo

end subroutine linear_initial_slice

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
!  copy the current metric of the NR evolution to the grid, and set the
!  timestep limit for the particles from the expansion rate
!+
!-----------------------------------------------------------------------
subroutine load_nr_metric()
 use nr_api,          only:nr_get_metric,nr_dt_expansion
 use metric_et_utils, only:gcovgrid,gcongrid,sqrtggrid,metricderivsgrid,set_grid_ready,dtmax_evolve
 use part,            only:npart,divcurlv,xyzh,isdead_or_accreted
 real :: divvmax
 integer :: i

 call nr_get_metric(gcov_nr,gcon_nr,sqrtg_nr,dgcov_nr)
 !
 ! the particles also need small steps while the metric or the matter that
 ! sources it changes rapidly: limit both |alpha K| dt and the fractional
 ! change of rho* per step, |div v| dt, to courant_expansion
 !
 divvmax = 0.
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) divvmax = max(divvmax,abs(real(divcurlv(1,i))))
 enddo
 dtmax_evolve = real(min(nr_dt_expansion(real(divvmax,kind=nr_dp)),real(huge(dtmax_evolve),kind=nr_dp)))
 gcovgrid  = real(gcov_nr)
 gcongrid  = real(gcon_nr)
 sqrtggrid = real(sqrtg_nr)
 metricderivsgrid = real(dgcov_nr)
 call set_grid_ready()

end subroutine load_nr_metric

!-----------------------------------------------------------------------
!+
!  stress-energy tensor T_munu of the particles at the grid nodes, for
!  the metric gcov (with spatial derivatives metricderivs) on the grid
!  (the routine called by the NR code at each Runge-Kutta stage)
!+
!-----------------------------------------------------------------------
subroutine get_tmunu_grid(t,gcov,metricderivs,tmunu)
 real(nr_dp), intent(in)  :: t,gcov(0:,0:,:,:,:),metricderivs(0:,0:,:,:,:,:)
 real(nr_dp), intent(out) :: tmunu(0:,0:,:,:,:)

 call get_tmunu_rhostar_grid(t,gcov,metricderivs,tmunu)

end subroutine get_tmunu_grid

!-----------------------------------------------------------------------
!+
!  T_munu of the particles at the grid nodes at time t for the metric on
!  the grid, and optionally the conserved density rho*.
!
!  The particles are only advanced at the start of each hydro step (time),
!  while the NR code integrates across the step, so their state is
!  extrapolated to t: positions with the velocity, the kernel-summed rho*
!  with the continuity equation, rho*(t) = rho* exp(-div v (t - time)), with
!  h kept consistent with it (h ~ rho*^(-1/3)), and the conserved momentum
!  with the external (gravitational) force. This makes the source second
!  order in the step. cons2prim reads rho* from the particle array, so this
!  is replaced by the extrapolated values and restored afterwards
!+
!-----------------------------------------------------------------------
subroutine get_tmunu_rhostar_grid(t,gcov,metricderivs,tmunu,rhostar)
 use part,            only:npart,xyzh,pxyzu,vxyzu,fext,dens,rho,divcurlv
 use metric_et_utils, only:gcovgrid,metricderivsgrid,set_grid_ready
 use metric_tools,    only:init_metric
 use cons2prim,       only:cons2primall
 use extern_gr,       only:get_tmunu_all
 use timestep,        only:time
 real(nr_dp), intent(in)  :: t,gcov(0:,0:,:,:,:),metricderivs(0:,0:,:,:,:,:)
 real(nr_dp), intent(out) :: tmunu(0:,0:,:,:,:)
 real(nr_dp), intent(out), optional :: rhostar(:,:,:)
 real :: dtlag

 call allocate_work(npart)
 dtlag = real(t) - time

 ! particle state extrapolated to t
 rhosave_nr = rho(1:npart)
 call extrapolate_particles(npart,dtlag,xyzh,vxyzu,pxyzu,fext,divcurlv,rhosave_nr,&
                            xyzh_nr,pxyzu_nr,rho)

 ! metric at the particles, interpolated from the grid
 gcovgrid = real(gcov)
 metricderivsgrid = real(metricderivs)
 call set_grid_ready()
 vxyzu_nr = vxyzu(:,1:npart)
 dens_nr  = dens(1:npart)
 call init_metric(npart,xyzh_nr,metrics_nr,metricderivs_nr)

 ! primitive variables and T_munu on the particles
 call cons2primall(npart,xyzh_nr,metrics_nr,pxyzu_nr,vxyzu_nr,dens_nr,eos_vars_nr)
 call get_tmunu_all(npart,xyzh_nr,metrics_nr,vxyzu_nr,metricderivs_nr,dens_nr,tmunus_nr)

 call deposit_tmunu(npart,xyzh_nr,tmunus_nr,tmunu,rhostar)
 rho(1:npart) = rhosave_nr
 if (.not.allocated(tmunu_last)) allocate(tmunu_last,mold=tmunu)
 tmunu_last = tmunu

end subroutine get_tmunu_rhostar_grid

!-----------------------------------------------------------------------
!+
!  extrapolate the particle positions, smoothing lengths, conserved
!  density and conserved momentum by a time dt, using the velocity, the
!  velocity divergence and the external force at the start of the step.
!  Positions are wrapped back into the periodic domain
!+
!-----------------------------------------------------------------------
subroutine extrapolate_particles(npart,dt,xyzh,vxyzu,pxyzu,fext,divcurlv,rho0,xyzh_t,pxyzu_t,rho_t)
 use boundary, only:xmin,ymin,zmin,dxbound,dybound,dzbound
 integer, intent(in)  :: npart
 real,    intent(in)  :: dt,xyzh(:,:),vxyzu(:,:),pxyzu(:,:),fext(:,:),rho0(:)
 real(kind=4), intent(in) :: divcurlv(:,:)
 real,    intent(out) :: xyzh_t(:,:),pxyzu_t(:,:)
 real,    intent(inout) :: rho_t(:)
 real :: fac
 integer :: i

 !$omp parallel do default(none) &
 !$omp shared(npart,dt,xyzh,vxyzu,pxyzu,fext,divcurlv,rho0,xyzh_t,pxyzu_t,rho_t) &
 !$omp shared(xmin,ymin,zmin,dxbound,dybound,dzbound) &
 !$omp private(i,fac)
 do i=1,npart
    fac = exp(-real(divcurlv(1,i))*dt)
    rho_t(i) = rho0(i)*fac
    xyzh_t(1,i) = xmin + modulo(xyzh(1,i) + vxyzu(1,i)*dt - xmin,dxbound)
    xyzh_t(2,i) = ymin + modulo(xyzh(2,i) + vxyzu(2,i)*dt - ymin,dybound)
    xyzh_t(3,i) = zmin + modulo(xyzh(3,i) + vxyzu(3,i)*dt - zmin,dzbound)
    xyzh_t(4,i) = xyzh(4,i)/fac**(1./3.)
    pxyzu_t(:,i) = pxyzu(:,i)
    pxyzu_t(1:3,i) = pxyzu(1:3,i) + fext(1:3,i)*dt
 enddo
 !$omp end parallel do

end subroutine extrapolate_particles

!-----------------------------------------------------------------------
!+
!  interpolate T_munu from the particles to the grid nodes, using the
!  normalised SPH interpolant with weights m/(rho* h^3), and scale it so
!  that the rest mass on the grid, integral of rho*, equals that of the
!  particles. Optionally returns the (normalised) rho* on the grid
!+
!-----------------------------------------------------------------------
subroutine deposit_tmunu(npart,xyzh,tmunus,tmunu,rhostar)
 use part,             only:rho,massoftype,igas,isdead_or_accreted,iphase,iamtype,maxphase,maxp
 use metric_et_utils,  only:gridsize,dxgrid,gridorigin,nr_deposit,nr_deposit_hmin
 use interpolations3D, only:interpolate3D_vecexact,interpolate3D_nodes
 use einsteintk_utils, only:set_rendering
 integer,     intent(in)  :: npart
 real,        intent(in)  :: xyzh(:,:)
 real,        intent(in)  :: tmunus(:,:,:)
 real(nr_dp), intent(out) :: tmunu(0:,0:,:,:,:)
 real(nr_dp), intent(out), optional :: rhostar(:,:,:)
 integer, parameter :: ndat = 17
 real, allocatable :: dat(:,:),weights(:),datsmooth(:,:,:,:)
 integer, allocatable :: itype(:)
 integer :: i,mu,nu,itypei
 real :: cfac,mtot,pmassi

 allocate(dat(npart,ndat),weights(npart),itype(npart))
 allocate(datsmooth(ndat,gridsize(1),gridsize(2),gridsize(3)))
 mtot = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,tmunus,rho,massoftype,iphase,dat,weights,itype,maxphase,maxp) &
 !$omp private(i,mu,nu,itypei,pmassi) &
 !$omp reduction(+:mtot)
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       itype(i) = -1
       weights(i) = 0.
       dat(i,:) = 0.
    else
       itypei = igas
       if (maxphase==maxp) itypei = iamtype(iphase(i))
       pmassi = massoftype(itypei)
       mtot = mtot + pmassi
       itype(i) = itypei
       weights(i) = pmassi/(rho(i)*xyzh(4,i)**3)
       do nu=0,3
          do mu=0,3
             dat(i,1 + mu + 4*nu) = tmunus(mu+1,nu+1,i)
          enddo
       enddo
       dat(i,ndat) = rho(i)
    endif
 enddo
 !$omp end parallel do

 if (nr_deposit == 1 .or. nr_deposit == 2) then
    ! kernel sampled at the nodes, optionally with zero second moment, with
    ! h >= nr_deposit_hmin dx, or integrated exactly over the cells if the
    ! smoothing length is still below a grid spacing
    call interpolate3D_nodes(xyzh,weights,dat,ndat,itype,npart,gridorigin,dxgrid,gridsize,&
                             nr_deposit==2,nr_deposit_hmin*maxval(dxgrid),.true.,datsmooth)
 else
    ! normalised SPH interpolation, with the kernel integrated exactly over
    ! the cells if nr_deposit = 3. Pixel i of the interpolation is at
    ! xmin + i*dx, so shifting the origin by one cell puts pixel i on node i
    call set_rendering(nr_deposit == 3)
    call interpolate3D_vecexact(xyzh,weights,dat,ndat,itype,npart,gridorigin(1)-dxgrid(1),&
                                gridorigin(2)-dxgrid(2),gridorigin(3)-dxgrid(3),datsmooth,&
                                gridsize(1),gridsize(2),gridsize(3),dxgrid(1),dxgrid(2),dxgrid(3),&
                                .true.,.true.,.true.,.true.)
 endif

 cfac = mtot/(sum(datsmooth(ndat,:,:,:))*product(dxgrid))
 do nu=0,3
    do mu=0,3
       tmunu(mu,nu,:,:,:) = cfac*datsmooth(1 + mu + 4*nu,:,:,:)
    enddo
 enddo
 if (present(rhostar)) rhostar = cfac*datsmooth(ndat,:,:,:)

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
    deallocate(metrics_nr,metricderivs_nr,vxyzu_nr,dens_nr,eos_vars_nr,tmunus_nr,rhosave_nr)
    deallocate(xyzh_nr,pxyzu_nr)
 endif
 allocate(metrics_nr(0:3,0:3,2,npart),metricderivs_nr(0:3,0:3,3,npart))
 allocate(vxyzu_nr(size(vxyzu,1),npart),dens_nr(npart),eos_vars_nr(size(eos_vars,1),npart))
 allocate(tmunus_nr(0:3,0:3,npart),rhosave_nr(npart))
 allocate(xyzh_nr(4,npart),pxyzu_nr(size(vxyzu,1),npart))

end subroutine allocate_work

!-----------------------------------------------------------------------
!+
!  write the grid-averaged scale factor a = det(gamma)^(1/6), lapse,
!  rest-mass density of the particles, energy density
!  rho_E = n^mu n^nu T_munu, L2 norm of the constraints, amplitude of the
!  metric perturbation phi (from det(gamma)^(1/3) = a^2 (1 - 2 phi)),
!  amplitude of the density contrast of the particles, L1 norm of the
!  Hamiltonian constraint and maximum rest-mass density to the _nr.ev file
!+
!-----------------------------------------------------------------------
subroutine write_nr_ev(t,nsteps)
 use physcon, only:pi
 use io,      only:id,master
 use nr_api,  only:nr_constraints
 real(nr_dp), intent(in) :: t
 integer,     intent(in) :: nsteps
 real(nr_dp) :: cnorm,hnorm,mnorm,abar,alphabar,rhoE,alpha,rhobar,phimin,phimax,phi,smag
 integer :: i,j,k,mu,nu,nnodes

 ! constraints with the source for the current metric (sets tmunu_last and dens_nr)
 call nr_constraints(cnorm,hnorm,mnorm)
 abar = 0.
 alphabar = 0.
 rhoE = 0.
 smag = 0.
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
          ! |S_i| with S_i = -n^mu T_mu i = alpha g^(mu 0) T_mu i
          smag = smag + alpha*sqrt(sum(matmul(gcon_nr(:,0,i,j,k),tmunu_last(:,1:3,i,j,k))**2))
       enddo
    enddo
 enddo
 abar = abar/nnodes
 phimin = huge(phimin)
 phimax = -huge(phimax)
 do k=1,size(sqrtg_nr,3)
    do j=1,size(sqrtg_nr,2)
       do i=1,size(sqrtg_nr,1)
          alpha = 1./sqrt(-gcon_nr(0,0,i,j,k))
          phi = 0.5*(1. - (sqrtg_nr(i,j,k)/alpha)**(2./3.)/abar**2)
          phimin = min(phimin,phi)
          phimax = max(phimax,phi)
       enddo
    enddo
 enddo
 rhobar = sum(dens_nr)/size(dens_nr)
 if (id==master) then
    write(iunit_nr,"(es18.10,1x,i10,13(1x,es18.10))") t,nsteps,abar,alphabar/nnodes,&
                    rhobar,rhoE/nnodes,cnorm,0.5*(phimax - phimin),&
                    0.5*(maxval(dens_nr) - minval(dens_nr))/rhobar,hnorm,maxval(dens_nr),&
                    mnorm,8.*pi*smag/nnodes,sqrt(sum((dens_nr/rhobar - 1.)**2)/size(dens_nr)),&
                    minval(dens_nr)/rhobar
    flush(iunit_nr)
 endif

end subroutine write_nr_ev

end module extern_nr
