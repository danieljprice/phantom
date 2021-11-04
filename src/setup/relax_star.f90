!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module relaxstar
!
! Automated relaxation of stellar density profile,
!   iterating towards hydrostatic equilibrium
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - maxits   : *maximum number of relaxation iterations*
!   - tol_dens : *% error in density to stop relaxation*
!   - tol_ekin : *tolerance on ekin/epot to stop relaxation*
!
! :Dependencies: checksetup, damping, deriv, energies, eos, fileutils,
!   infile_utils, initial, io, io_summary, memory, options, part, physcon,
!   ptmass, readwrite_dumps, sortutils, step_lf_global, table_utils, units
!
 implicit none
 public :: relax_star,write_options_relax,read_options_relax

 real,    private :: tol_ekin = 1.e-7 ! criteria for being converged
 real,    private :: tol_dens = 1.   ! allow 1% RMS error in density
 integer, private :: maxits = 1000

 real,    private :: gammaprev,hfactprev
 integer, private :: ieos_prev

 private

contains

!----------------------------------------------------------------
!+
!  relax a star to hydrostatic equilibrium. We run the main
!  code but with a fake equation of state, low neighbour number
!  and fixing the entropy as a function of r
!
!  IN:
!    rho(nt)   - tabulated density as function of r (in code units)
!    pr(nt)    - tabulated pressure as function of r (in code units)
!    r(nt)     - radius for each point in the table
!
!  IN/OUT:
!    xyzh(:,:) - positions and smoothing lengths of all particles
!+
!----------------------------------------------------------------
subroutine relax_star(nt,rho,pr,r,npart,xyzh)
 use table_utils, only:yinterp
 use deriv,       only:get_derivs_global
 use part,        only:vxyzu
 use step_lf_global, only:init_step,step
 use initial,       only:initialise
 use memory,      only:allocate_memory
 use energies,    only:compute_energies,ekin,epot,etherm
 use checksetup,  only:check_setup
 use io,          only:error,warning
 use fileutils,   only:getnextfilename
 use readwrite_dumps, only:write_fulldump,init_readwrite_dumps
 use eos, only:gamma
 use physcon,     only:pi
 use options,     only:iexternalforce
 use io_summary,  only:summary_initialise
 integer, intent(in)    :: nt
 integer, intent(inout) :: npart
 real,    intent(in)    :: rho(nt),pr(nt),r(nt)
 real,    intent(inout) :: xyzh(:,:)
 integer :: nits,nerr,nwarn,iunit
 real    :: t,dt,dtmax,rmserr,rstar,mstar,tdyn
 real    :: entrop(nt),utherm(nt),mr(nt),rmax,dtext,dtnew
 logical :: converged,use_step
 logical, parameter :: fix_entrop = .true. ! fix entropy instead of thermal energy
 logical, parameter :: write_files = .true.
 character(len=20) :: filename
 !
 ! save settings and set a bunch of options
 !
 rstar = maxval(r)
 mr = get_mr(rho,r)
 mstar = mr(nt)
 tdyn  = 2.*pi*sqrt(rstar**3/(32.*mstar))
 print*,'rstar  = ',rstar,' mstar = ',mstar, ' tdyn = ',tdyn
 call set_options_for_relaxation(tdyn)
 call summary_initialise()
 !
 ! check particle setup is sensible
 !
 call check_setup(nwarn,nerr)
 if (nerr > 0) then
    call error('relax_star','cannot relax star because particle setup contains errors')
    call restore_original_options()
    return
 endif
 use_step = .false.
 if (iexternalforce > 0) then
    call warning('relax_star','asynchronous shifting not implemented with external forces: evolving in time instead')
    use_step = .true.
 endif
 !
 ! define utherm(r) based on P(r) and rho(r)
 ! and use this to set the thermal energy of all particles
 !
 entrop = pr/rho**gamma
 utherm = pr/(rho*(gamma-1.))
 if (any(utherm <= 0.)) then
    call error('relax_star','relax-o-matic needs non-zero pressure array set in order to work')
    call restore_original_options()
    return
 endif
 call reset_u_and_get_errors(npart,xyzh,vxyzu,nt,mr,rho,utherm,entrop,fix_entrop,rmax,rmserr)
 !
 ! compute derivatives the first time around (needed if using actual step routine)
 !
 t = 0.
 call allocate_memory(2*npart)
 call get_derivs_global()
 call reset_u_and_get_errors(npart,xyzh,vxyzu,nt,mr,rho,utherm,entrop,fix_entrop,rmax,rmserr)
 call compute_energies(t)
 !
 ! perform sanity checks
 !
 if (etherm > abs(epot)) then
    call error('relax_star','cannot relax star because it is unbound (etherm > epot)')
    print*,' Etherm = ',etherm,' Epot = ',Epot
    print*
    call restore_original_options()
    return
 endif
 print "(/,3(a,1pg11.3),/,a,0pf6.2,a,es11.3,a)",&
   ' RELAX-A-STAR-O-MATIC: Etherm:',etherm,' Epot:',Epot, ' R*:',maxval(r), &
   '       WILL stop WHEN: dens error < ',tol_dens,'% AND Ekin/Epot < ',tol_ekin,' OR Iter=0'

 filename = 'relax_00000'
 if (write_files) then
    call init_readwrite_dumps()
    call write_fulldump(t,filename)
    open(newunit=iunit,file='relax.ev',status='replace')
    write(iunit,"(a)") '# nits,rmax,etherm,epot,ekin/epot,L2_{err}'
 endif
 converged = .false.
 dt = epsilon(0.) ! To avoid error in sink-gas substepping
 dtext = huge(dtext)
 if (use_step) then
    dtmax = tdyn
    call init_step(npart,t,dtmax)
 endif
 nits = 0
 do while (.not. converged)
    nits = nits + 1
    !
    ! shift particles by one "timestep"
    !
    t = t + dt
    if (use_step) then
       call step(npart,npart,t,dt,dtext,dtnew)
       dt = dtnew
    else
       call shift_particles(npart,xyzh,vxyzu,dt)
    endif
    !
    ! reset thermal energy and calculate information
    !
    call reset_u_and_get_errors(npart,xyzh,vxyzu,nt,mr,rho,utherm,entrop,fix_entrop,rmax,rmserr)
    !
    ! compute energies and check for convergence
    !
    call compute_energies(t)
    converged = ((ekin > 0. .and. ekin/abs(epot) < tol_ekin .and. &
                 rmserr < 0.01*tol_dens) .or. nits >= maxits)
    !
    ! print information to screen
    !
    if (use_step) then
       print "(a,es10.3,a,2pf6.2,2(a,1pg11.3))",' Relaxing star: t/dyn:',t/tdyn,', dens error:',rmserr,'%, R*:',rmax, &
        ' Ekin/Epot:',ekin/abs(epot)
    else
       print "(a,i4,a,i4,a,2pf6.2,2(a,1pg11.3))",' Relaxing star: Iter',nits,'/',maxits, &
             ', dens error:',rmserr,'%, R*:',rmax,' Ekin/Epot:',ekin/abs(epot)
    endif
    !
    ! additional diagnostic output, mainly for debugging/checking
    !
    if (write_files) then
       !
       ! write information to the relax.ev file
       !
       write(iunit,*) nits,rmax,etherm,epot,ekin/abs(epot),rmserr
       !
       ! write dump files
       !
       if (mod(nits,5)==0) then
          filename = getnextfilename(filename)
          call write_fulldump(t,filename)
          call flush(iunit)
       endif
    endif
 enddo
 if (write_files) close(iunit)
 !
 ! unfake some things
 !
 call restore_original_options()

end subroutine relax_star

!----------------------------------------------------------------
!+
!  shift particles: this is like timestepping but done
!  asynchronously. Each particle shifts by dx = 0.5*dt^2*a
!  where dt is the local courant timestep, i.e. h/c_s
!+
!----------------------------------------------------------------
subroutine shift_particles(npart,xyzh,vxyzu,dtmin)
 use deriv, only:get_derivs_global
 use part,  only:fxyzu,fext,xyzmh_ptmass,nptmass,rhoh,massoftype,igas
 use ptmass,only:get_accel_sink_gas
 use eos,   only:get_spsound
 use options, only:ieos
 integer, intent(in) :: npart
 real, intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real, intent(out)   :: dtmin
 real :: dx(3),dti,phi,rhoi,cs,hi
 integer :: i,nlargeshift
!
! shift particles asynchronously
!
 dtmin = huge(dtmin)
 nlargeshift = 0
 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(npart,xyzh,vxyzu,fxyzu,fext,xyzmh_ptmass,nptmass,massoftype,ieos) &
 !$omp private(i,dx,dti,phi,cs,rhoi,hi) &
 !$omp reduction(min:dtmin) &
 !$omp reduction(+:nlargeshift)
 do i=1,npart
    fext(1:3,i) = 0.
    if (nptmass > 0) then
       call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                               xyzmh_ptmass,fext(1,i),fext(2,i),fext(3,i),phi)
    endif
    hi = xyzh(4,i)
    rhoi = rhoh(hi,massoftype(igas))
    cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
    dti = 0.3*hi/cs   ! h/cs
    dx  = 0.5*dti**2*(fxyzu(1:3,i) + fext(1:3,i))
    if (dot_product(dx,dx) > hi**2) then
       dx = dx / sqrt(dot_product(dx,dx)) * hi  ! Avoid large shift in particle position
       nlargeshift = nlargeshift + 1
    endif
    xyzh(1:3,i) = xyzh(1:3,i) + dx(:)
    vxyzu(1:3,i) = dx(:)/dti ! fake velocities, so we can measure kinetic energy
    dtmin = min(dtmin,dti)   ! used to print a "time" in the output (but it is fake)
 enddo
 !$omp end parallel do
 if (nlargeshift > 0) print*,'Warning: Restricted dx for ', nlargeshift, 'particles'
 !
 ! get forces on particles
 !
 call get_derivs_global()

end subroutine shift_particles

!----------------------------------------------------------------
!+
!  reset the thermal energy to be exactly p(r)/((gam-1)*rho(r))
!  according to the desired p(r) and rho(r)
!  also compute error between true rho(r) and desired rho(r)
!+
!----------------------------------------------------------------
subroutine reset_u_and_get_errors(npart,xyzh,vxyzu,nt,mr,rho,utherm,entrop,fix_entrop,rmax,rmserr)
 use table_utils, only:yinterp
 use sortutils,   only:find_rank,r2func
 use part,        only:rhoh,massoftype,igas,maxvxyzu,iorder=>ll
 use eos,         only:gamma
 integer, intent(in) :: npart,nt
 real, intent(in)    :: xyzh(:,:),mr(nt),rho(nt),utherm(nt),entrop(nt)
 real, intent(inout) :: vxyzu(:,:)
 real, intent(out)   :: rmax,rmserr
 logical, intent(in) :: fix_entrop
 real :: ri,rhor,rhoi,rho1,mstar,massri
 integer :: i

 rho1 = yinterp(rho,mr,0.)
 rmax = 0.
 rmserr = 0.
 call find_rank(npart,r2func,xyzh(1:3,:),iorder)
 mstar = mr(nt)
 do i = 1,npart
    ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    massri = mstar * real(iorder(i)-1) / real(npart)
    rhor = yinterp(rho,mr,massri) ! analytic rho(r)
    rhoi = rhoh(xyzh(4,i),massoftype(igas)) ! actual rho
    if (maxvxyzu >= 4) then
       if (fix_entrop) then
          vxyzu(4,i) = (yinterp(entrop,mr,massri)*rhoi**(gamma-1.))/(gamma-1.)
       else
          vxyzu(4,i) = yinterp(utherm,mr,massri)
       endif
    endif
    rmserr = rmserr + (rhor - rhoi)**2
    rmax   = max(rmax,ri)
 enddo
 rmserr = sqrt(rmserr/npart)/rho1

end subroutine reset_u_and_get_errors

!----------------------------------------------------------------
!+
!  set code options specific to relaxation calculations
!+
!----------------------------------------------------------------
subroutine set_options_for_relaxation(tdyn)
 use eos,  only:ieos,gamma
 use part, only:hfact,maxvxyzu
 use damping, only:damp,tdyn_s
 use options, only:idamp
 use units,   only:utime
 real, intent(in) :: tdyn

 gammaprev = gamma
 hfactprev = hfact
 ieos_prev = ieos
 !
 ! turn on settings appropriate to relaxation
 !
 !gamma = 2.
 !hfact = 0.8 !0.7
 if (maxvxyzu >= 4) ieos = 2
 if (tdyn > 0.) then
    idamp = 2
    tdyn_s = tdyn*utime
 else
    idamp = 1
    damp = 0.05
 endif

end subroutine set_options_for_relaxation

!--------------------------------------------------
!+
!  get mass coordinate m(r) = \int 4.*pi*rho*r^2 dr
!+
!--------------------------------------------------
function get_mr(rho,r) result(mr)
 use physcon, only:pi
 real, intent(in)  :: rho(:),r(:)
 real :: mr(size(r))
 integer :: i

 mr(1) = 0.
 do i=2,size(rho)
    mr(i) = mr(i-1) + 4./3.*pi*rho(i) * (r(i) - r(i-1)) * (r(i)**2 + r(i)*r(i-1) + r(i-1)**2)
 enddo

end function get_mr

!----------------------------------------------------------------
!+
!  restore previous settings
!+
!----------------------------------------------------------------
subroutine restore_original_options
 use eos,     only:ieos,gamma
 use damping, only:damp
 use options, only:idamp
 use part,    only:hfact,vxyzu,npart

 gamma = gammaprev
 hfact = hfactprev
 ieos  = ieos_prev
 idamp = 0
 damp = 0.
 vxyzu(1:3,1:npart) = 0.

end subroutine restore_original_options

!----------------------------------------------------------------
!+
!  write relaxation options to .setup file
!+
!----------------------------------------------------------------
subroutine write_options_relax(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(tol_ekin,'tol_ekin','tolerance on ekin/epot to stop relaxation',iunit)
 call write_inopt(tol_dens,'tol_dens','% error in density to stop relaxation',iunit)
 call write_inopt(maxits,'maxits','maximum number of relaxation iterations',iunit)

end subroutine write_options_relax

!----------------------------------------------------------------
!+
!  read relaxation options from .setup file
!+
!----------------------------------------------------------------
subroutine read_options_relax(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(tol_ekin,'tol_ekin',db,errcount=nerr)
 call read_inopt(tol_dens,'tol_dens',db,errcount=nerr)
 call read_inopt(maxits,'maxits',db,errcount=nerr)

end subroutine read_options_relax

end module relaxstar
