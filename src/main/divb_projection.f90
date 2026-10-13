!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module divb_projection
!
! Adjoint projection of the magnetic field onto a discretely
! divergence-free subspace. The divergence D is the difference operator
! used for cleaning; its adjoint G is the symmetric SPH gradient used
! for grad psi. Both use the kernel gradient corrected for smoothing
! lengths set from number density, so the pair matches the force loop
! also for unequal particle masses (adaptive particle refinement).
! The scalar system D M_V^-1 D^T pi = D B is solved by Jacobi-
! preconditioned conjugate gradient, once per call. The solve stops at
! the practical bound on h|div B|/|B| of Tsukamoto (2026), or earlier if
! the relative residual tolerance is met.
!
! Single MPI rank only: the neighbour list has no remote particles and
! the reductions are local.
!
! :References:
!   Tsukamoto (2026), arXiv:2606.28197
!   Tricco & Price (2012), J. Comput. Phys. 231, 7214-7236
!   Tricco, Price & Bate (2016), J. Comput. Phys. 322, 326-344
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, derivutils, dim, io, kernel, neighcompact, part,
!   pcg, sortutils
!
 implicit none
 public :: project_divb,get_divb,divb_adjoint_mismatch

 private

 ! compacted neighbour list and edge weights, rebuilt on every call:
 ! d_ij = cedge(1,ic)*rhat(:,ic) and d_ji = -cedge(2,ic)*rhat(:,ic),
 ! stored in single precision to save memory; D and G read the same
 ! values, so they remain exact adjoints
 integer, allocatable :: ivar(:,:),ijvar(:)
 real(kind=4), allocatable :: rhat(:,:),cedge(:,:)
 real,    allocatable :: pmass(:),aval(:),qsum(:,:),honb(:)
 logical, allocatable :: active(:)
 integer              :: ncompact = 0

 ! conjugate-gradient vectors and scratch for the stopping test
 real,    allocatable :: pimult(:),svec(:),zvec(:),pvec(:),wvec(:),gvec(:,:)
 real,    allocatable :: chi(:),dchi(:)

 ! chi just after the previous projection, and the step it was taken
 real,    allocatable :: chi_prev(:)
 integer              :: nsteps_prev = -1

 ! state of the current solve, read by the stopping test
 real :: rnorm0 = 0.,rnorm = 0.,tol_rel = 0.,chi_top_ref = 0.,chi_fred = 0.1

 ! practical stop on chi = h|div B|/|B| (Tsukamoto 2026): RMS below
 ! 1e-5, and the worst 1% of the change since the last projection
 ! reduced by 0.1/N_interval
 real, parameter :: chi_abs   = 1.e-5
 real, parameter :: chi_ftop  = 0.01
 real, parameter :: chi_fred0 = 0.1

 ! floor on |B| in chi, as a fraction of the rms |B|, so that particles
 ! with little or no field do not dominate the practical stop
 real, parameter :: bfloor_frac = 1.e-6

contains

!-----------------------------------------------------------------------
!+
!  Project B so that the discrete divergence is driven to the solver
!  tolerance, then set the cleaning scalar psi and its time derivative
!  to zero. Bxyz returns the projected primitive field.
!+
!-----------------------------------------------------------------------
subroutine project_divb(npart,xyzh,rho,gradh,Bevol,Bxyz,dBevol,divcurlB,tol,maxit,nsteps,ierr,&
                        divb_before,divb_after,emag_ratio,nits)
 use io,   only:id,master,iprint,warning
 use pcg,  only:pcg_solve
 integer,      intent(in)    :: npart,maxit,nsteps
 real,         intent(in)    :: xyzh(:,:),rho(:),tol
 real(kind=4), intent(in)    :: gradh(:,:)
 real,         intent(inout) :: Bevol(:,:),Bxyz(:,:),dBevol(:,:)
 real(kind=4), intent(inout) :: divcurlB(:,:)
 integer,      intent(out)   :: ierr
 real,         intent(out), optional :: divb_before,divb_after,emag_ratio
 integer,      intent(out), optional :: nits
 integer :: i,niter,ninterval,ierr_pcg
 real    :: eb0,eb1,ratio

 ierr   = 0
 niter  = 0
 rnorm0 = 0.
 rnorm  = 0.
 ratio  = 1.
 if (npart > 0) call prepare_stencil(npart,xyzh,rho,gradh,ierr)
 if (npart <= 0 .or. ierr /= 0) then
    call set_optional_outputs(divb_before,divb_after,emag_ratio,rnorm0,rnorm,ratio)
    if (present(nits)) nits = 0
    call release_storage
    return
 endif

 ! primitive field from the evolved variable B/rho
 !$omp parallel do default(none) private(i) shared(npart,Bxyz,Bevol,rho)
 do i=1,npart
    if (rho(i) > 0.) Bxyz(1:3,i) = Bevol(1:3,i)*rho(i)
 enddo
 !$omp end parallel do

 call apply_div(npart,Bxyz,svec)
 eb0 = magnetic_energy(npart,Bxyz,rho)
 call set_h_on_B(npart,xyzh,Bxyz)

 ! steps since the previous projection set how far the worst
 ! particles must be cleaned before the practical stop
 ninterval = nsteps
 if (nsteps_prev >= 0) ninterval = nsteps - nsteps_prev
 call reset_chi_prev(npart)
 call set_chi_reference(npart,max(ninterval,1))

 tol_rel = tol
 rnorm0  = norm_v(npart,svec,rho)
 rnorm   = rnorm0
 if (rnorm0 > tiny(rnorm0)) then
    call pcg_solve(npart,apply_laplacian,apply_jacobi,projection_converged,maxit, &
                   pimult,svec,zvec,pvec,wvec,niter,ierr_pcg)
    if (ierr_pcg /= 0) call warning('divb_proj','conjugate gradient direction lost, stopping')
 else
    pimult(1:npart) = 0.
 endif

 ! B <- B - G pi, the symmetric-gradient correction
 call apply_grad(npart,pimult,gvec,rho)
 !$omp parallel do default(none) private(i) shared(npart,active,Bxyz,gvec,Bevol,rho)
 do i=1,npart
    if (.not. active(i)) cycle
    Bxyz(1:3,i)  = Bxyz(1:3,i) - gvec(1:3,i)
    Bevol(1:3,i) = Bxyz(1:3,i)/rho(i)
 enddo
 !$omp end parallel do

 ! psi was the cleaning field of the pre-projection divergence; zero
 ! its rate too, otherwise the next half-kick re-sources it
 Bevol(4,1:npart)  = 0.
 dBevol(4,1:npart) = 0.

 ! store the difference divergence of the projected field
 call apply_div(npart,Bxyz,svec)
 !$omp parallel do default(none) private(i) shared(npart,active,divcurlB,svec)
 do i=1,npart
    if (active(i)) divcurlB(1,i) = real(svec(i),kind=4)
 enddo
 !$omp end parallel do
 rnorm = norm_v(npart,svec,rho)
 call set_h_on_B(npart,xyzh,Bxyz)
 call store_chi_prev(npart)
 nsteps_prev = nsteps

 eb1 = magnetic_energy(npart,Bxyz,rho)
 if (eb0 > 0.) ratio = eb1/eb0
 if (id == master) then
    write(iprint,'(a,i6,a,es10.3,a,es10.3,a,es10.3)') &
       ' div B projection:',niter,' iterations, |div B|_V ',rnorm0,' -> ',rnorm,', dE_B/E_B ',ratio - 1.
 endif
 if (niter >= maxit .and. rnorm > tol*(rnorm0 + tiny(rnorm0))) then
    call warning('divb_proj','iteration cap reached before div B tolerance')
 endif
 call set_optional_outputs(divb_before,divb_after,emag_ratio,rnorm0,rnorm,ratio)
 if (present(nits)) nits = niter
 call release_storage

contains

!-----------------------------------------------------------------------
!+
!  w = L p = D G p
!+
!-----------------------------------------------------------------------
subroutine apply_laplacian(p,w)
 real, intent(in)  :: p(:)
 real, intent(out) :: w(:)

 call apply_grad(npart,p,gvec,rho)
 call apply_div(npart,gvec,w)

end subroutine apply_laplacian

!-----------------------------------------------------------------------
!+
!  Stopping test: the relative residual tolerance is met, or both
!  practical tests are met: the rms of chi is below the absolute floor,
!  and either the worst 1% of |chi - chi_prev| has been reduced by
!  chi_fred or that top tail is itself below the floor.
!+
!-----------------------------------------------------------------------
logical function projection_converged(r)
 real, intent(in) :: r(:)
 real    :: rms_all,rms_top,rms_dtop
 logical :: top_ok

 rnorm = norm_v(npart,r,rho)
 projection_converged = (rnorm <= tol_rel*(rnorm0 + tiny(rnorm0)))
 if (projection_converged) return

 call measure_chi(r,rms_all,rms_top,rms_dtop)
 top_ok = (rms_dtop <= chi_fred*chi_top_ref) .or. (rms_top < chi_abs)
 projection_converged = top_ok .and. (rms_all < chi_abs)

end function projection_converged

end subroutine project_divb

!-----------------------------------------------------------------------
!+
!  Difference divergence (D B)_i of the primitive field on the current
!  particle distribution (zero on particles outside the active set)
!+
!-----------------------------------------------------------------------
subroutine get_divb(npart,xyzh,rho,gradh,Bxyz,divb,ierr)
 integer,      intent(in)  :: npart
 real,         intent(in)  :: xyzh(:,:),rho(:),Bxyz(:,:)
 real(kind=4), intent(in)  :: gradh(:,:)
 real,         intent(out) :: divb(:)
 integer,      intent(out) :: ierr

 divb(1:npart) = 0.
 call prepare_stencil(npart,xyzh,rho,gradh,ierr)
 if (ierr == 0) call apply_div(npart,Bxyz,divb)
 call release_storage

end subroutine get_divb

!-----------------------------------------------------------------------
!+
!  Relative mismatch in the adjoint identity
!  sum_i pi_i (D X)_i = sum_i V_i (G pi)_i . X_i
!  for given test fields X and pi, which must vanish outside the
!  active set (e.g. on boundary particles)
!+
!-----------------------------------------------------------------------
subroutine divb_adjoint_mismatch(npart,xyzh,rho,gradh,xvec,piscal,mismatch,ierr)
 integer,      intent(in)  :: npart
 real,         intent(in)  :: xyzh(:,:),rho(:),xvec(:,:),piscal(:)
 real(kind=4), intent(in)  :: gradh(:,:)
 real,         intent(out) :: mismatch
 integer,      intent(out) :: ierr
 integer :: i
 real    :: lhs,rhs,scale

 mismatch = 0.
 call prepare_stencil(npart,xyzh,rho,gradh,ierr)
 if (ierr /= 0) then
    call release_storage
    return
 endif
 call apply_div(npart,xvec,svec)
 call apply_grad(npart,piscal,gvec,rho)

 lhs = 0.
 rhs = 0.
 do i=1,npart
    if (.not. active(i)) cycle
    lhs = lhs + piscal(i)*svec(i)
    rhs = rhs + pmass(i)/rho(i)*dot_product(gvec(1:3,i),xvec(1:3,i))
 enddo
 scale = abs(lhs) + abs(rhs)
 if (scale > 0.) mismatch = abs(lhs - rhs)/scale
 call release_storage

end subroutine divb_adjoint_mismatch

!-----------------------------------------------------------------------
!+
!  Allocate the workspace, build the symmetric neighbour list, the
!  edge weights and the Jacobi preconditioner
!+
!-----------------------------------------------------------------------
subroutine prepare_stencil(npart,xyzh,rho,gradh,ierr)
 use io, only:warning,fatal
 integer,      intent(in)  :: npart
 real,         intent(in)  :: xyzh(:,:),rho(:)
 real(kind=4), intent(in)  :: gradh(:,:)
 integer,      intent(out) :: ierr

 call release_storage
 allocate(active(npart),pmass(npart),aval(npart),qsum(3,npart),honb(npart), &
          pimult(npart),svec(npart),zvec(npart),pvec(npart),wvec(npart),gvec(3,npart), &
          chi(npart),dchi(npart),ivar(3,npart),stat=ierr)
 if (ierr /= 0) call fatal('divb_proj','out of memory: set use_divb_proj = F')

 call build_neighbour_list(npart,xyzh)
 if (ncompact <= 0) then
    call warning('divb_proj','empty neighbour list; build the tree before projecting')
    ierr = 1
    return
 endif
 call build_weights(npart,xyzh,rho,gradh)
 call build_preconditioner(rho)

end subroutine prepare_stencil

!-----------------------------------------------------------------------
!+
!  Fill the compacted neighbour list from the current tree, sizing the
!  edge storage from the actual neighbour count
!+
!-----------------------------------------------------------------------
subroutine build_neighbour_list(npart,xyzh)
 use io,           only:fatal
 use kernel,       only:radkern
 use part,         only:hfact
 use neighcompact, only:compact_neigh_counts,get_compacted_neighbour_list
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 integer :: nneigh_average,icompactmax,ncompactlocal,nedges,ierr
 integer(kind=8) :: nedges8

 ! first guess of 1.5 times the mean neighbour number; retry if short
 call compact_neigh_counts(npart,radkern,hfact,nneigh_average,icompactmax)
 nedges8 = int(1.5*nneigh_average,kind=8)*int(npart,kind=8)
 if (nedges8 > huge(nedges)) call fatal('divb_proj','too many neighbour pairs: set use_divb_proj = F')
 nedges = max(int(nedges8),1)
 do
    allocate(ijvar(nedges),stat=ierr)
    if (ierr /= 0) call fatal('divb_proj','out of memory: set use_divb_proj = F')
    call get_compacted_neighbour_list(xyzh,ivar,ijvar,ncompact,ncompactlocal,nedges)
    if (nedges <= size(ijvar)) exit
    deallocate(ijvar)
 enddo
 ncompact = ncompactlocal
 allocate(rhat(3,max(nedges,1)),cedge(2,max(nedges,1)),stat=ierr)
 if (ierr /= 0) call fatal('divb_proj','out of memory: set use_divb_proj = F')

end subroutine build_neighbour_list

!-----------------------------------------------------------------------
!+
!  Particle volumes, the active set and the directed edge coefficients.
!  cedge(1) is m_j/rho_i times the corrected kernel gradient at h_i,
!  cedge(2) the same coefficient for the reversed pair, so the gradient
!  gather is exactly the adjoint of the divergence.
!+
!-----------------------------------------------------------------------
subroutine build_weights(npart,xyzh,rho,gradh)
 use dim,      only:periodic,maxphase,maxp,igradomega,igradzeta
 use boundary, only:mod_boundary
 use derivutils, only:get_grkern_corrected
 use part,     only:isdead_or_accreted,iphase,iamtype,iamboundary,igas,&
                    get_pmass
 integer,      intent(in) :: npart
 real,         intent(in) :: xyzh(:,:),rho(:)
 real(kind=4), intent(in) :: gradh(:,:)
 integer :: n,i,j,k,ic,itype
 real    :: dx,dy,dz,rij,hi1,hj1,qi,qj,pmassi,pmassj,zetaomegai,zetaomegaj
 real    :: grki,grkj,grkt
 logical :: use_gas

 ! masses of gas and boundary particles; other types carry no field
 use_gas = .not.(maxphase==maxp)
 !$omp parallel do default(none) private(i,itype) &
 !$omp shared(npart,xyzh,rho,pmass,active,iphase,use_gas,maxphase,maxp)
 do i=1,npart
    pmass(i) = 0.
    active(i) = .false.
    if (isdead_or_accreted(xyzh(4,i)) .or. rho(i) <= 0.) cycle
    itype = igas
    if (maxphase==maxp) itype = iamtype(iphase(i))
    if (itype /= igas .and. .not.iamboundary(itype)) cycle
    pmass(i) = get_pmass(i,use_gas)
 enddo
 !$omp end parallel do

 ! rows of the list are gas or boundary; boundary particles are fixed
 do n=1,ncompact
    i = ivar(3,n)
    active(i) = pmass(i) > 0.
    if (maxphase==maxp) active(i) = active(i) .and. .not.iamboundary(iamtype(iphase(i)))
 enddo

 !$omp parallel do default(none) &
 !$omp private(n,i,j,k,ic,dx,dy,dz,rij,hi1,hj1,qi,qj,pmassi,pmassj,zetaomegai,zetaomegaj,grki,grkj,grkt) &
 !$omp shared(ncompact,ivar,ijvar,active,pmass,rho,xyzh,gradh,rhat,cedge)
 do n=1,ncompact
    i = ivar(3,n)
    hi1 = 1./xyzh(4,i)
    pmassi = pmass(i)
    zetaomegai = real(gradh(igradzeta,i))*real(gradh(igradomega,i))
    do k=1,ivar(1,n)
       ic = ivar(2,n) + k
       j  = ijvar(ic)
       rhat(:,ic)  = 0.
       cedge(:,ic) = 0.
       if (.not.active(i) .or. pmass(j) <= 0.) cycle
       dx = xyzh(1,i) - xyzh(1,j)
       dy = xyzh(2,i) - xyzh(2,j)
       dz = xyzh(3,i) - xyzh(3,j)
       if (periodic) call mod_boundary(dx,dy,dz)
       rij = sqrt(dx*dx + dy*dy + dz*dz)
       if (rij <= 0.) cycle
       rhat(:,ic) = real((/dx,dy,dz/)/rij,kind=4)
       pmassj = pmass(j)

       ! d_ij = m_j/rho_i * grad_i W_ij(h_i), corrected for h(n)
       qi = rij*hi1
       call get_grkern_corrected(qi*qi,qi,hi1**4,zetaomegai,pmassj,grki,grkt)
       cedge(1,ic) = real(pmassj*grki/rho(i),kind=4)

       ! same coefficient on the reversed pair; grad_j W points along -rhat
       hj1 = 1./xyzh(4,j)
       qj = rij*hj1
       zetaomegaj = real(gradh(igradzeta,j))*real(gradh(igradomega,j))
       call get_grkern_corrected(qj*qj,qj,hj1**4,zetaomegaj,pmassi,grkj,grkt)
       cedge(2,ic) = real(pmassi*grkj/rho(j),kind=4)
    enddo
 enddo
 !$omp end parallel do

end subroutine build_weights

!-----------------------------------------------------------------------
!+
!  Jacobi preconditioner: diagonal of L = D M_V^{-1} D^T, plus the row
!  sums of d_ij that the gradient reapplies on every iteration
!+
!-----------------------------------------------------------------------
subroutine build_preconditioner(rho)
 real, intent(in) :: rho(:)
 integer :: n,i,j,k,ic
 real    :: sumd2

 aval = 0.
 qsum = 0.
 !$omp parallel do default(none) private(n,i,j,k,ic,sumd2) &
 !$omp shared(ncompact,ivar,ijvar,active,pmass,rho,cedge,rhat,qsum,aval)
 do n=1,ncompact
    i = ivar(3,n)
    if (.not. active(i)) cycle
    sumd2 = 0.
    do k=1,ivar(1,n)
       ic = ivar(2,n) + k
       j  = ijvar(ic)
       qsum(:,i) = qsum(:,i) + real(cedge(1,ic))*rhat(:,ic)
       if (active(j)) sumd2 = sumd2 + real(cedge(1,ic))**2*rho(j)/pmass(j)
    enddo
    aval(i) = dot_product(qsum(:,i),qsum(:,i))*rho(i)/pmass(i) + sumd2
 enddo
 !$omp end parallel do

end subroutine build_preconditioner

!-----------------------------------------------------------------------
!+
!  (D X)_i = sum_j d_ij . (X_j - X_i) on active particles
!+
!-----------------------------------------------------------------------
subroutine apply_div(npart,vec,divs)
 integer, intent(in)  :: npart
 real,    intent(in)  :: vec(:,:)
 real,    intent(out) :: divs(:)
 integer :: n,i,j,k,ic
 real    :: acc

 divs(1:npart) = 0.
 !$omp parallel do default(none) private(n,i,j,k,ic,acc) &
 !$omp shared(ncompact,ivar,ijvar,cedge,rhat,active,vec,divs)
 do n=1,ncompact
    i = ivar(3,n)
    if (.not. active(i)) cycle
    acc = 0.
    do k=1,ivar(1,n)
       ic = ivar(2,n) + k
       j  = ijvar(ic)
       acc = acc + cedge(1,ic)*dot_product(rhat(:,ic),vec(1:3,j) - vec(1:3,i))
    enddo
    divs(i) = acc
 enddo
 !$omp end parallel do

end subroutine apply_div

!-----------------------------------------------------------------------
!+
!  (G pi)_i = (1/V_i) * [sum_j pi_j d_ji - pi_i sum_j d_ij]
!  which is -rho times the symmetric gradient of psi = pi/V
!+
!-----------------------------------------------------------------------
subroutine apply_grad(npart,piscal,g,rho)
 integer, intent(in)  :: npart
 real,    intent(in)  :: piscal(:),rho(:)
 real,    intent(out) :: g(:,:)
 integer :: n,i,j,k,ic
 real    :: sumj(3)

 g(:,1:npart) = 0.
 !$omp parallel do default(none) private(n,i,j,k,ic,sumj) &
 !$omp shared(ncompact,ivar,ijvar,cedge,rhat,qsum,active,pmass,rho,piscal,g)
 do n=1,ncompact
    i = ivar(3,n)
    if (.not. active(i)) cycle
    sumj = 0.
    do k=1,ivar(1,n)
       ic = ivar(2,n) + k
       j  = ijvar(ic)
       if (active(j)) sumj = sumj - piscal(j)*cedge(2,ic)*rhat(:,ic)
    enddo
    g(:,i) = (sumj - piscal(i)*qsum(:,i))*rho(i)/pmass(i)
 enddo
 !$omp end parallel do

end subroutine apply_grad

!-----------------------------------------------------------------------
!+
!  z = P^{-1} r with P = diag(a_i). Particles with a vanishing
!  diagonal (no neighbours) are left out of the solve.
!+
!-----------------------------------------------------------------------
subroutine apply_jacobi(r,z)
 real, intent(in)  :: r(:)
 real, intent(out) :: z(:)
 integer :: i

 !$omp parallel do default(none) private(i) shared(active,aval,r,z)
 do i=1,size(r)
    if (active(i) .and. aval(i) > 0.) then
       z(i) = r(i)/aval(i)
    else
       z(i) = 0.
    endif
 enddo
 !$omp end parallel do

end subroutine apply_jacobi

!-----------------------------------------------------------------------
!+
!  Volume-weighted norm ||q||_V on the active set
!+
!-----------------------------------------------------------------------
real function norm_v(npart,q,rho)
 integer, intent(in) :: npart
 real,    intent(in) :: q(:),rho(:)
 integer :: i
 real    :: acc

 acc = 0.
 !$omp parallel do default(none) private(i) shared(npart,active,pmass,rho,q) reduction(+:acc)
 do i=1,npart
    if (active(i)) acc = acc + pmass(i)/rho(i)*q(i)*q(i)
 enddo
 !$omp end parallel do
 norm_v = sqrt(acc)

end function norm_v

!-----------------------------------------------------------------------
!+
!  Discrete magnetic energy sum_i V_i |B_i|^2 on the active set
!+
!-----------------------------------------------------------------------
real function magnetic_energy(npart,Bxyz,rho)
 integer, intent(in) :: npart
 real,    intent(in) :: Bxyz(:,:),rho(:)
 integer :: i
 real    :: acc

 acc = 0.
 !$omp parallel do default(none) private(i) shared(npart,active,pmass,rho,Bxyz) reduction(+:acc)
 do i=1,npart
    if (active(i)) acc = acc + pmass(i)/rho(i)*dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
 enddo
 !$omp end parallel do
 magnetic_energy = acc

end function magnetic_energy

!-----------------------------------------------------------------------
!+
!  h_i/|B_i|, so that chi_i = h_i |div B|_i/|B_i| = |div B|_i * honb_i
!+
!-----------------------------------------------------------------------
subroutine set_h_on_B(npart,xyzh,Bxyz)
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:),Bxyz(:,:)
 integer :: i,nact
 real    :: b2sum,bfloor

 ! floor |B| at a fraction of its rms over the active set
 b2sum = 0.
 nact  = 0
 !$omp parallel do default(none) private(i) shared(npart,active,Bxyz) reduction(+:b2sum,nact)
 do i=1,npart
    if (active(i)) then
       b2sum = b2sum + dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
       nact  = nact + 1
    endif
 enddo
 !$omp end parallel do
 bfloor = bfloor_frac*sqrt(b2sum/max(nact,1)) + tiny(0.)

 !$omp parallel do default(none) private(i) shared(npart,active,honb,xyzh,Bxyz,bfloor)
 do i=1,npart
    honb(i) = 0.
    if (active(i)) honb(i) = xyzh(4,i)/max(norm2(Bxyz(1:3,i)),bfloor)
 enddo
 !$omp end parallel do

end subroutine set_h_on_B

!-----------------------------------------------------------------------
!+
!  Start chi_prev from zero if the particle set has changed size
!  since the last projection (e.g. refinement, accretion, injection).
!  Reordering at fixed npart is not detected.
!+
!-----------------------------------------------------------------------
subroutine reset_chi_prev(npart)
 integer, intent(in) :: npart

 if (allocated(chi_prev)) then
    if (size(chi_prev) == npart) return
    deallocate(chi_prev)
 endif
 allocate(chi_prev(npart))
 chi_prev = 0.

end subroutine reset_chi_prev

!-----------------------------------------------------------------------
!+
!  Record the opening top-tail divergence and the reduction factor
!  0.1/N_interval used by the practical stop
!+
!-----------------------------------------------------------------------
subroutine set_chi_reference(npart,ninterval)
 integer, intent(in) :: npart,ninterval
 real :: rms_all,rms_top

 call measure_chi(svec(1:npart),rms_all,rms_top,chi_top_ref)
 chi_fred = chi_fred0/real(ninterval)

end subroutine set_chi_reference

!-----------------------------------------------------------------------
!+
!  RMS of chi_i = h_i |r_i|/|B_i| over the active set, and the RMS of
!  the top 1% of chi and of |chi - chi_prev|, using the uncorrected
!  field for |B|
!+
!-----------------------------------------------------------------------
subroutine measure_chi(r,rms_all,rms_top,rms_dtop)
 real,    intent(in)  :: r(:)
 real,    intent(out) :: rms_all,rms_top,rms_dtop
 integer :: i,nact,ntop
 real    :: chii,acc

 nact = 0
 acc  = 0.
 do i=1,size(r)
    if (.not. active(i)) cycle
    chii = abs(r(i))*honb(i)
    nact = nact + 1
    chi(nact)  = chii
    dchi(nact) = abs(chii - chi_prev(i))
    acc = acc + chii*chii
 enddo
 rms_all  = 0.
 rms_top  = 0.
 rms_dtop = 0.
 if (nact <= 0) return
 rms_all  = sqrt(acc/real(nact))
 ntop     = min(max(1,ceiling(chi_ftop*real(nact))),nact)
 rms_top  = rms_of_largest(nact,ntop,chi)
 rms_dtop = rms_of_largest(nact,ntop,dchi)

end subroutine measure_chi

!-----------------------------------------------------------------------
!+
!  Root-mean-square of the ntop largest entries of a(1:n); a is reordered
!+
!-----------------------------------------------------------------------
real function rms_of_largest(n,ntop,a)
 use sortutils, only:select_inplace
 integer, intent(in)    :: n,ntop
 real,    intent(inout) :: a(:)

 ! partition so that the ntop largest values are in a(n-ntop+1:n)
 call select_inplace(n,n-ntop+1,a)
 rms_of_largest = sqrt(sum(a(n-ntop+1:n)**2)/real(ntop))

end function rms_of_largest

!-----------------------------------------------------------------------
!+
!  chi just after this projection, kept for the next interval
!+
!-----------------------------------------------------------------------
subroutine store_chi_prev(npart)
 integer, intent(in) :: npart

 chi_prev(1:npart) = abs(svec(1:npart))*honb(1:npart)

end subroutine store_chi_prev

!-----------------------------------------------------------------------
!+
!  Drop the workspace; only chi_prev is kept between calls
!+
!-----------------------------------------------------------------------
subroutine release_storage()

 if (allocated(ivar))   deallocate(ivar)
 if (allocated(ijvar))  deallocate(ijvar)
 if (allocated(rhat))   deallocate(rhat)
 if (allocated(cedge))  deallocate(cedge)
 if (allocated(active)) deallocate(active)
 if (allocated(pmass))  deallocate(pmass)
 if (allocated(aval))   deallocate(aval)
 if (allocated(qsum))   deallocate(qsum)
 if (allocated(honb))   deallocate(honb)
 if (allocated(pimult)) deallocate(pimult)
 if (allocated(svec))   deallocate(svec)
 if (allocated(zvec))   deallocate(zvec)
 if (allocated(pvec))   deallocate(pvec)
 if (allocated(wvec))   deallocate(wvec)
 if (allocated(gvec))   deallocate(gvec)
 if (allocated(chi))    deallocate(chi)
 if (allocated(dchi))   deallocate(dchi)
 ncompact = 0

end subroutine release_storage

!-----------------------------------------------------------------------
!+
!  Optional diagnostics returned to the unit test
!+
!-----------------------------------------------------------------------
subroutine set_optional_outputs(divb_before,divb_after,emag_ratio,rbefore,rafter,ratio)
 real, intent(out), optional :: divb_before,divb_after,emag_ratio
 real, intent(in) :: rbefore,rafter,ratio

 if (present(divb_before)) divb_before = rbefore
 if (present(divb_after))  divb_after  = rafter
 if (present(emag_ratio))  emag_ratio  = ratio

end subroutine set_optional_outputs

end module divb_projection
