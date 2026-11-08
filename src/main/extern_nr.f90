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
! With amr_maxlevel > 0 the NR grid is refined where the particles need
! it: a particle with smoothing length h needs a grid spacing below
! h/nr_refine_h, so that its kernel is resolved by the nodes it is
! deposited to (and the NR code may also refine on the curvature of the
! metric). The refined blocks are rebuilt from the particles before each
! advance, the dust source is deposited on every block with the smoothing
! length floored at that block's grid spacing, and the particles take the
! metric from the finest block containing them.
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
 use nr_api,          only:nr_linked,nr_init,nr_formulation,nr_maxlevel
 use boundary,        only:xmin,ymin,zmin,dxbound,dybound,dzbound
 use metric_et_utils, only:allocate_grid,ngrid_evolve,periodic_grid,advance_metric,nr_solve_init,&
                            proper_time,scale_factor,redshift_init,nr_dust_source,nr_deposit
 use part,            only:npart,eos_vars,igasP
 real, intent(in) :: time
 integer :: n(3),i,ierr
 logical :: changed
 real(nr_dp) :: dx(3),origin(3),hubble,alphamean,amean
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
 ! the dust source needs pressureless matter
 if (nr_dust_source .and. npart > 0) then
    if (maxval(abs(eos_vars(igasP,1:npart))) > 0.) then
       nr_dust_source = .false.
       if (id==master) write(iprint,"(a)") ' nr_dust_source = F: the matter has pressure'
    endif
 endif
 ! refined blocks from the particles (one more level per pass)
 if (nr_maxlevel() > 0) then
    if (.not.nr_dust_source) call fatal('extern_nr','mesh refinement needs the dust source (nr_dust_source)')
    if (nr_deposit /= 1 .and. nr_deposit /= 2) call fatal('extern_nr','mesh refinement needs nr_deposit = 1 or 2')
    do i=1,nr_maxlevel()
       call regrid_nr(changed)
       if (.not.changed) exit
    enddo
 endif
 call load_nr_metric()
 advance_metric => advance_nr_metric
 ! proper time (age) at the start, unless restarting from a dump that has it:
 ! for a cosmological run that of the Einstein-de Sitter universe, 2/(3 H)
 if (proper_time < 0.) then
    proper_time = 0.
    if (redshift_init >= 0.) proper_time = real(2./(3.*hubble))
 endif
 call grid_means(alphamean,amean)
 scale_factor = real(amean)

 if (id==master) then
    write(iprint,"(a,a,a,i4,a)") ' evolving the metric with NR (',trim(nr_formulation()),') on a ',&
                                 ngrid_evolve,'^3 periodic grid'
    write(iprint,"(a,es12.5,a,es10.3)") ' initial conformal Hubble rate = ',hubble,&
                                        ', max |phi| = ',maxval(abs(phi))
    open(newunit=iunit_nr,file=trim(fileprefix)//'_nr.ev',status='replace',form='formatted')
    write(iunit_nr,"('#',18(1x,'[',i2.2,1x,a11,']',2x))") 1,'time',2,'nsteps',3,'a',4,'alpha',&
          5,'rho',6,'rho_E',7,'constraints',8,'phi amp',9,'delta amp',10,'ham L1',11,'rho max',&
          12,'mom L1',13,'8 pi |S|',14,'rms delta*',15,'min rho*',16,'proper time',17,'redshift',&
          18,'nblocks'
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
 use nr_api,          only:nr_time,nr_evolve_to,nr_maxlevel
 use metric_et_utils, only:proper_time,scale_factor,nr_dust_source
 real, intent(in) :: time
 real(nr_dp) :: t0,alpha0,alpha1,abar
 integer :: nsteps
 logical :: changed

 ! ignore requests that differ from the current time only by round-off
 if (time - nr_time() > 1.e-12*abs(time)) then
    t0 = nr_time()
    call grid_means(alpha0,abar)
    if (nr_maxlevel() > 0) call regrid_nr(changed)
    if (nr_dust_source) call set_dust_source(time)
    call nr_evolve_to(real(time,kind=nr_dp),nsteps)
    call load_nr_metric()
    ! proper time of the mean lapse, d tau = alpha dt (trapezoidal rule)
    call grid_means(alpha1,abar)
    proper_time  = real(proper_time + 0.5*(alpha0 + alpha1)*(nr_time() - t0))
    scale_factor = real(abar)
    call write_nr_ev(nr_time(),nsteps)
 endif

end subroutine advance_nr_metric

!-----------------------------------------------------------------------
!+
!  matter source for the NR evolution up to time t1, for dust: the
!  coordinate densities rho* and rho* u_i deposited from the particles at
!  their current time and extrapolated to t1. The NR code interpolates
!  them linearly in time and computes T_munu in the metric of each stage,
!  so the particles are deposited twice per advance instead of being
!  processed at every Runge-Kutta stage of every NR step
!+
!-----------------------------------------------------------------------
subroutine set_dust_source(t1)
 use part,     only:npart,xyzh,pxyzu,vxyzu,fext,rho,divcurlv
 use timestep, only:time
 use nr_api,   only:nr_set_dust_source,nr_nblocks,nr_blocksize
 real, intent(in) :: t1
 real(nr_dp), allocatable :: dens0(:,:,:),mom0(:,:,:,:),dens1(:,:,:),mom1(:,:,:,:)
 real, allocatable :: bdat0(:,:,:,:,:),bdat1(:,:,:,:,:)
 integer :: n(3),nb,id

 call allocate_work(npart)
 n = shape(sqrtg_nr)
 nb = nr_blocksize()
 allocate(dens0(n(1),n(2),n(3)),dens1(n(1),n(2),n(3)),mom0(3,n(1),n(2),n(3)),mom1(3,n(1),n(2),n(3)))
 allocate(bdat0(5,nb,nb,nb,nr_nblocks()),bdat1(5,nb,nb,nb,nr_nblocks()))
 call deposit_dust(npart,xyzh,pxyzu,dens0,mom0,bdat0)
 ! particles extrapolated to t1 (rho* temporarily replaced, as in the callback)
 rhosave_nr = rho(1:npart)
 call extrapolate_particles(npart,t1 - time,xyzh,vxyzu,pxyzu,fext,divcurlv,rhosave_nr,&
                            xyzh_nr,pxyzu_nr,rho)
 call deposit_dust(npart,xyzh_nr,pxyzu_nr,dens1,mom1,bdat1)
 rho(1:npart) = rhosave_nr
 call nr_set_dust_source(real(time,kind=nr_dp),dens0,mom0,real(t1,kind=nr_dp),dens1,mom1)
 do id=1,nr_nblocks()
    call nr_set_dust_source(real(time,kind=nr_dp),real(bdat0(5,:,:,:,id),kind=nr_dp),&
                            real(bdat0(1:3,:,:,:,id),kind=nr_dp),real(t1,kind=nr_dp),&
                            real(bdat1(5,:,:,:,id),kind=nr_dp),real(bdat1(1:3,:,:,:,id),kind=nr_dp),id)
 enddo

end subroutine set_dust_source

!-----------------------------------------------------------------------
!+
!  rebuild the refined blocks of the NR grid: each particle whose
!  smoothing length h is below nr_refine_h times the root grid spacing
!  needs refinement down to the level whose spacing is below h/nr_refine_h,
!  over the radius of its kernel
!+
!-----------------------------------------------------------------------
subroutine regrid_nr(changed)
 use part,            only:npart,xyzh,isdead_or_accreted
 use nr_api,          only:nr_regrid,nr_maxlevel
 use metric_et_utils, only:dxgrid,nr_refine_h
 use kernel,          only:radkern
 logical, intent(out) :: changed
 real(nr_dp), allocatable :: xpt(:,:),rpt(:)
 integer, allocatable :: lpt(:)
 real :: hres
 integer :: i,np

 ! particles needing refinement have h < hres
 hres = nr_refine_h*minval(dxgrid)
 np = 0
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i)) .and. xyzh(4,i) < hres) np = np + 1
 enddo
 allocate(xpt(3,np),rpt(np),lpt(np))
 np = 0
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i)) .and. xyzh(4,i) < hres) then
       np = np + 1
       xpt(:,np) = xyzh(1:3,i)
       rpt(np) = radkern*xyzh(4,i)
       lpt(np) = min(nr_maxlevel(),ceiling(log(hres/xyzh(4,i))/log(2.)))
    endif
 enddo
 call nr_regrid(changed,xpt,rpt,lpt)

end subroutine regrid_nr

!-----------------------------------------------------------------------
!+
!  mean lapse and mean scale factor gamma^(1/6) over the grid nodes
!+
!-----------------------------------------------------------------------
subroutine grid_means(alphabar,abar)
 real(nr_dp), intent(out) :: alphabar,abar
 real(nr_dp) :: alpha
 integer :: i,j,k

 alphabar = 0.
 abar = 0.
 do k=1,size(sqrtg_nr,3)
    do j=1,size(sqrtg_nr,2)
       do i=1,size(sqrtg_nr,1)
          ! g^00 = -1/alpha^2 and sqrt(-g) = alpha sqrt(gamma)
          alpha = 1./sqrt(-gcon_nr(0,0,i,j,k))
          alphabar = alphabar + alpha
          abar = abar + (sqrtg_nr(i,j,k)/alpha)**(1./3.)
       enddo
    enddo
 enddo
 alphabar = alphabar/size(sqrtg_nr)
 abar = abar/size(sqrtg_nr)

end subroutine grid_means

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
 call load_nr_blocks()

end subroutine load_nr_metric

!-----------------------------------------------------------------------
!+
!  copy the metric of the refined blocks of the NR grid (if any)
!+
!-----------------------------------------------------------------------
subroutine load_nr_blocks()
 use nr_api,          only:nr_nblocks,nr_blocksize,nr_get_block,nr_get_block_metric
 use metric_et_utils, only:allocate_blocks,set_block_ready,mblocks
 real(nr_dp), allocatable :: gcov(:,:,:,:,:),gcon(:,:,:,:,:),sqrtg(:,:,:),dgcov(:,:,:,:,:,:)
 real(nr_dp) :: origin(3),dx(3)
 integer :: id,m,level,parent,child(8)

 call allocate_blocks(nr_nblocks(),nr_blocksize())
 m = nr_blocksize() + 2
 allocate(gcov(0:3,0:3,m,m,m),gcon(0:3,0:3,m,m,m),sqrtg(m,m,m),dgcov(0:3,0:3,3,m,m,m))
 do id=1,nr_nblocks()
    call nr_get_block(id,level,origin,dx,parent,child)
    call nr_get_block_metric(id,gcov,gcon,sqrtg,dgcov)
    mblocks(id)%level  = level
    mblocks(id)%child  = child
    mblocks(id)%x0     = real(origin)
    mblocks(id)%dx     = real(dx)
    mblocks(id)%gcov   = real(gcov)
    mblocks(id)%derivs = real(dgcov)
    call set_block_ready(id)
 enddo

end subroutine load_nr_blocks

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
 use part,            only:npart,xyzh,pxyzu,vxyzu,fext,dens,rho,divcurlv,eos_vars
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
 ! cons2prim takes gamma (and other inputs) from eos_vars
 eos_vars_nr = eos_vars(:,1:npart)
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
!  interpolate T_munu from the particles to the grid nodes, normalised so
!  that the rest mass on the grid, integral of rho*, equals that of the
!  particles. Optionally returns the (normalised) rho* on the grid
!+
!-----------------------------------------------------------------------
subroutine deposit_tmunu(npart,xyzh,tmunus,tmunu,rhostar)
 use part, only:rho
 integer,     intent(in)  :: npart
 real,        intent(in)  :: xyzh(:,:)
 real,        intent(in)  :: tmunus(:,:,:)
 real(nr_dp), intent(out) :: tmunu(0:,0:,:,:,:)
 real(nr_dp), intent(out), optional :: rhostar(:,:,:)
 integer, parameter :: ndat = 17
 real, allocatable :: dat(:,:),datgrid(:,:,:,:)
 integer :: i,mu,nu

 allocate(dat(npart,ndat))
 !$omp parallel do default(none) shared(npart,tmunus,rho,dat) private(i,mu,nu)
 do i=1,npart
    do nu=0,3
       do mu=0,3
          dat(i,1 + mu + 4*nu) = tmunus(mu+1,nu+1,i)
       enddo
    enddo
    dat(i,ndat) = rho(i)
 enddo
 !$omp end parallel do
 call deposit_grid(npart,xyzh,dat,ndat,datgrid)
 do nu=0,3
    do mu=0,3
       tmunu(mu,nu,:,:,:) = datgrid(1 + mu + 4*nu,:,:,:)
    enddo
 enddo
 if (present(rhostar)) rhostar = datgrid(ndat,:,:,:)

end subroutine deposit_tmunu

!-----------------------------------------------------------------------
!+
!  coordinate densities of dust on the grid: the conserved density
!  D = rho* and the momentum density D u_i (with u_i = p_i for dust),
!  normalised so that the rest mass on the grid equals that of the
!  particles, and optionally the same on the refined blocks,
!  bdat(1:3,:,:,:,id) = D u_i and bdat(5,:,:,:,id) = D (bdat(4,...) is
!  the partition of unity, see deposit_blocks)
!+
!-----------------------------------------------------------------------
subroutine deposit_dust(npart,xyzh,pxyzu,dens,mom,bdat)
 use part, only:rho
 integer,     intent(in)  :: npart
 real,        intent(in)  :: xyzh(:,:),pxyzu(:,:)
 real(nr_dp), intent(out) :: dens(:,:,:),mom(:,:,:,:)
 real,        intent(out), optional :: bdat(:,:,:,:,:)
 real, allocatable :: dat(:,:),datgrid(:,:,:,:)
 integer :: i

 ! D u_i, a field of ones (whose deposit is the partition of unity) and D
 allocate(dat(npart,5))
 !$omp parallel do default(none) shared(npart,pxyzu,rho,dat) private(i)
 do i=1,npart
    dat(i,1:3) = rho(i)*pxyzu(1:3,i)
    dat(i,4) = 1.
    dat(i,5) = rho(i)
 enddo
 !$omp end parallel do
 call deposit_grid(npart,xyzh,dat,5,datgrid)
 dens = datgrid(5,:,:,:)
 mom  = datgrid(1:3,:,:,:)
 if (present(bdat)) then
    if (size(bdat,5) > 0) call deposit_blocks(npart,xyzh,dat,5,4,sum(datgrid(4,:,:,:))/size(dens),bdat)
 endif

end subroutine deposit_dust

!-----------------------------------------------------------------------
!+
!  interpolate the particle data dat(npart,ndat), whose last column is
!  rho*, to the grid nodes with weights m/(rho* h^3) (the deposit chosen
!  by nr_deposit), and scale all fields so that the rest mass on the grid,
!  integral of rho*, equals that of the particles
!+
!-----------------------------------------------------------------------
subroutine deposit_grid(npart,xyzh,dat,ndat,datgrid)
 use metric_et_utils,  only:gridsize,dxgrid,gridorigin,nr_deposit,nr_deposit_hmin
 use interpolations3D, only:interpolate3D_vecexact,interpolate3D_nodes
 use einsteintk_utils, only:set_rendering
 integer, intent(in)  :: npart,ndat
 real,    intent(in)  :: xyzh(:,:),dat(:,:)
 real, allocatable, intent(out) :: datgrid(:,:,:,:)
 real, allocatable :: weights(:)
 integer, allocatable :: itype(:)
 real :: cfac,mtot

 allocate(datgrid(ndat,gridsize(1),gridsize(2),gridsize(3)))
 call deposit_weights(npart,xyzh,weights,itype,mtot)

 if (nr_deposit == 1 .or. nr_deposit == 2) then
    ! kernel sampled at the nodes, optionally with zero second moment, with
    ! h >= nr_deposit_hmin dx, or integrated exactly over the cells if the
    ! smoothing length is still below a grid spacing
    call interpolate3D_nodes(xyzh,weights,dat,ndat,itype,npart,gridorigin,dxgrid,gridsize,&
                             nr_deposit==2,nr_deposit_hmin*maxval(dxgrid),.true.,datgrid)
 else
    ! normalised SPH interpolation, with the kernel integrated exactly over
    ! the cells if nr_deposit = 3. Pixel i of the interpolation is at
    ! xmin + i*dx, so shifting the origin by one cell puts pixel i on node i
    call set_rendering(nr_deposit == 3)
    call interpolate3D_vecexact(xyzh,weights,dat,ndat,itype,npart,gridorigin(1)-dxgrid(1),&
                                gridorigin(2)-dxgrid(2),gridorigin(3)-dxgrid(3),datgrid,&
                                gridsize(1),gridsize(2),gridsize(3),dxgrid(1),dxgrid(2),dxgrid(3),&
                                .true.,.true.,.true.,.true.)
 endif

 cfac = mtot/(sum(datgrid(ndat,:,:,:))*product(dxgrid))
 datgrid = cfac*datgrid

end subroutine deposit_grid

!-----------------------------------------------------------------------
!+
!  interpolate the particle data dat(npart,ndat) to the nodes of every
!  refined block of the NR grid as deposit_grid does to the root grid
!  (nr_deposit = 1 or 2), with the smoothing length floored at
!  nr_deposit_hmin times the block's grid spacing. The particles are first
!  sorted into lists of the blocks that their kernels overlap.
!
!  The kernel sampled at the nodes is a partition of unity only to a
!  fraction of a per cent, by an amount that depends on h/dx, while the
!  root grid is normalised by its total mass. Column ione of dat holds
!  ones, whose deposit sum_b (m_b/rho*_b) W_b is the partition of unity,
!  and pnorm is its mean on the (normalised) root grid: each level is
!  scaled so that the mean of its partition of unity is also pnorm
!+
!-----------------------------------------------------------------------
subroutine deposit_blocks(npart,xyzh,dat,ndat,ione,pnorm,bdat)
 use nr_api,           only:nr_nblocks,nr_blocksize,nr_get_block,nr_maxlevel
 use metric_et_utils,  only:gridsize,nr_deposit,nr_deposit_hmin
 use interpolations3D, only:interpolate3D_nodes
 integer, intent(in)  :: npart,ndat,ione
 real,    intent(in)  :: xyzh(:,:),dat(:,:),pnorm
 real,    intent(out) :: bdat(:,:,:,:,:)
 real, allocatable :: weights(:),xyzh_b(:,:),dat_b(:,:),w_b(:)
 integer, allocatable :: itype(:),nin(:),first(:),plist(:),itype_b(:),lev(:)
 real(nr_dp) :: origin(3),dx(3)
 real :: mtot,x0(3),dxl(3),psum(0:20),nsum(0:20)
 integer :: id,nblk,nb,m,j,level,parent,child(8),nbv(3),nper(3)

 nblk = nr_nblocks()
 nb = nr_blocksize()
 call deposit_weights(npart,xyzh,weights,itype,mtot)
 allocate(nin(nblk),first(nblk))
 ! count the particles overlapping each block, then list them
 call list_block_particles(npart,xyzh,itype,nr_maxlevel(),nb,nin,first)
 allocate(plist(sum(nin)))
 call list_block_particles(npart,xyzh,itype,nr_maxlevel(),nb,nin,first,plist)
 nbv = nb
 allocate(lev(nblk))
 psum = 0.
 nsum = 0.
 !$omp parallel do default(none) schedule(dynamic) &
 !$omp shared(nblk,nb,nbv,nin,first,plist,xyzh,weights,dat,ndat,ione,itype,bdat,gridsize,lev) &
 !$omp shared(nr_deposit,nr_deposit_hmin) &
 !$omp private(id,m,j,level,parent,child,origin,dx,x0,dxl,nper,xyzh_b,w_b,dat_b,itype_b) &
 !$omp reduction(+:psum,nsum)
 do id=1,nblk
    call nr_get_block(id,level,origin,dx,parent,child)
    lev(id) = level
    x0 = real(origin)
    dxl = real(dx)
    nper = gridsize*2**level
    m = nin(id)
    allocate(xyzh_b(4,m),w_b(m),dat_b(m,ndat),itype_b(m))
    do j=1,m
       xyzh_b(:,j) = xyzh(1:4,plist(first(id)+j-1))
       w_b(j)      = weights(plist(first(id)+j-1))
       dat_b(j,:)  = dat(plist(first(id)+j-1),1:ndat)
       itype_b(j)  = itype(plist(first(id)+j-1))
    enddo
    call interpolate3D_nodes(xyzh_b,w_b,dat_b,ndat,itype_b,m,x0,dxl,nbv,nr_deposit==2,&
                             nr_deposit_hmin*maxval(dxl),.true.,bdat(:,:,:,:,id),nper)
    psum(level) = psum(level) + sum(bdat(ione,:,:,:,id))
    nsum(level) = nsum(level) + nb**3
    deallocate(xyzh_b,w_b,dat_b,itype_b)
 enddo
 !$omp end parallel do
 ! normalise each level to the mean partition of unity of the root grid
 do id=1,nblk
    bdat(:,:,:,:,id) = bdat(:,:,:,:,id)*pnorm*nsum(lev(id))/psum(lev(id))
 enddo

end subroutine deposit_blocks

!-----------------------------------------------------------------------
!+
!  particles whose kernels (with the smoothing length floored as in the
!  deposit) overlap each refined block: without plist, count them (nin)
!  and set the offsets first; with plist, list them from first(id)
!+
!-----------------------------------------------------------------------
subroutine list_block_particles(npart,xyzh,itype,maxlevel,nb,nin,first,plist)
 use nr_api,          only:nr_find_block
 use metric_et_utils, only:gridsize,dxgrid,gridorigin,nr_deposit_hmin
 use kernel,          only:radkern
 integer, intent(in)    :: npart,itype(:),maxlevel,nb
 real,    intent(in)    :: xyzh(:,:)
 integer, intent(inout) :: nin(:),first(:)
 integer, intent(out), optional :: plist(:)
 integer :: i,l,id,ci,cj,ck,c(3),clo(3),chi(3),nl(3)
 real :: dxl(3),size_l(3),r
 logical :: counting

 counting = .not.present(plist)
 nin = 0
 do i=1,npart
    if (itype(i) < 0) cycle
    do l=1,maxlevel
       dxl = dxgrid/2.**l
       size_l = nb*dxl
       nl = gridsize*2**l/nb
       r = radkern*max(xyzh(4,i),nr_deposit_hmin*maxval(dxl))
       clo = floor((xyzh(1:3,i) - r - gridorigin)/size_l)
       chi = floor((xyzh(1:3,i) + r - gridorigin)/size_l)
       ! no block twice when the kernel spans the whole (periodic) level
       chi = min(chi,clo + nl - 1)
       do ck=clo(3),chi(3)
          do cj=clo(2),chi(2)
             do ci=clo(1),chi(1)
                c = [ci,cj,ck]
                id = nr_find_block(l,c)
                if (id <= 0) cycle
                nin(id) = nin(id) + 1
                if (.not.counting) plist(first(id) + nin(id) - 1) = i
             enddo
          enddo
       enddo
    enddo
 enddo
 if (counting) then
    first(1) = 1
    do id=2,size(nin)
       first(id) = first(id-1) + nin(id-1)
    enddo
 endif

end subroutine list_block_particles

!-----------------------------------------------------------------------
!+
!  deposit weights m/(rho* h^3) and types of the particles (-1 for dead
!  or accreted particles), and their total mass
!+
!-----------------------------------------------------------------------
subroutine deposit_weights(npart,xyzh,weights,itype,mtot)
 use part, only:rho,massoftype,igas,isdead_or_accreted,iphase,iamtype,maxphase,maxp
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    allocatable, intent(out) :: weights(:)
 integer, allocatable, intent(out) :: itype(:)
 real,    intent(out) :: mtot
 integer :: i,itypei
 real :: pmassi

 allocate(weights(npart),itype(npart))
 mtot = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,rho,massoftype,iphase,weights,itype,maxphase,maxp) &
 !$omp private(i,itypei,pmassi) &
 !$omp reduction(+:mtot)
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       itype(i) = -1
       weights(i) = 0.
    else
       itypei = igas
       if (maxphase==maxp) itypei = iamtype(iphase(i))
       pmassi = massoftype(itypei)
       mtot = mtot + pmassi
       itype(i) = itypei
       weights(i) = pmassi/(rho(i)*xyzh(4,i)**3)
    endif
 enddo
 !$omp end parallel do

end subroutine deposit_weights

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
 use physcon,         only:pi
 use io,              only:id,master
 use nr_api,          only:nr_constraints,nr_nblocks
 use metric_et_utils, only:proper_time,redshift_init
 real(nr_dp), intent(in) :: t
 integer,     intent(in) :: nsteps
 real(nr_dp) :: cnorm,hnorm,mnorm,abar,alphabar,rhoE,alpha,rhobar,phimin,phimax,phi,smag,z
 real(nr_dp), allocatable :: tmunu(:,:,:,:,:)
 integer :: i,j,k,mu,nu,nnodes

 ! constraints with the source for the current metric (sets tmunu_last and
 ! dens_nr, which with refined blocks are set here, as the NR code then uses
 ! the dust source on every level)
 if (nr_nblocks() > 0) then
    allocate(tmunu,mold=tmunu_last)
    call get_tmunu_grid(t,gcov_nr,dgcov_nr,tmunu)
 endif
 call nr_constraints(cnorm,hnorm,mnorm)
 call grid_means(alphabar,abar)
 rhoE = 0.
 smag = 0.
 nnodes = size(sqrtg_nr)
 do k=1,size(sqrtg_nr,3)
    do j=1,size(sqrtg_nr,2)
       do i=1,size(sqrtg_nr,1)
          ! n^mu = -alpha g^(mu 0)
          alpha = 1./sqrt(-gcon_nr(0,0,i,j,k))
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
 z = -1.
 if (redshift_init >= 0.) z = (1. + redshift_init)/abar - 1.
 if (id==master) then
    write(iunit_nr,"(es18.10,1x,i10,15(1x,es18.10),1x,i10)") t,nsteps,abar,alphabar,&
                    rhobar,rhoE/nnodes,cnorm,0.5*(phimax - phimin),&
                    0.5*(maxval(dens_nr) - minval(dens_nr))/rhobar,hnorm,maxval(dens_nr),&
                    mnorm,8.*pi*smag/nnodes,sqrt(sum((dens_nr/rhobar - 1.)**2)/size(dens_nr)),&
                    minval(dens_nr)/rhobar,proper_time,z,nr_nblocks()
    flush(iunit_nr)
 endif

end subroutine write_nr_ev

end module extern_nr
