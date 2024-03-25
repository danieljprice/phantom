!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module ptmass
!
! This module contains everything to do with
!  sink / point mass particles
!
!  These are treated quite differently to SPH particles,
!  are not included in the neighbour lists and in principle
!  should be stored (as identical copies) on every MPI processor
!  NOTE: only certain types of particles are allowed to be accreted onto
!        sink particles (during creation or normal accretion).  The list
!        of 'accretable' particles is given in (and can be modified in)
!        in function 'is_accretable' in the 'part' module.
!
! :References: Bate, Bonnell & Price (1995), MNRAS 277, 362-376 [BBP95]
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - f_acc           : *particles < f_acc*h_acc accreted without checks*
!   - f_crit_override : *unconditional sink formation if rho > f_crit_override*rho_crit*
!   - h_acc           : *accretion radius for new sink particles*
!   - h_soft_sinkgas  : *softening length for new sink particles*
!   - h_soft_sinksink : *softening length between sink particles*
!   - icreate_sinks   : *allow automatic sink particle creation*
!   - r_crit          : *critical radius for point mass creation (no new sinks < r_crit from existing sink)*
!   - r_merge_cond    : *sinks will merge if bound within this radius*
!   - r_merge_uncond  : *sinks will unconditionally merge within this separation*
!   - rho_crit_cgs    : *density above which sink particles are created (g/cm^3)*
!
! :Dependencies: boundary, dim, eos, eos_barotropic, eos_piecewise,
!   extern_geopot, externalforces, fastmath, infile_utils, io, io_summary,
!   kdtree, kernel, linklist, mpidomain, mpiutils, options, part,
!   ptmass_heating, units, vectorutils
!
 use part, only:nsinkproperties,gravity,is_accretable,&
                ihsoft,ihacc,ispinx,ispiny,ispinz,imacc,iJ2,iReff
 use io,   only:iscfile,iskfile,id,master
 implicit none

 public :: init_ptmass, finish_ptmass
 public :: pt_write_sinkev, pt_close_sinkev
 public :: get_accel_sink_gas, get_accel_sink_sink
 public :: merge_sinks
 public :: ptmass_predictor, ptmass_corrector
 public :: ptmass_not_obscured
 public :: ptmass_accrete, ptmass_create
 public :: write_options_ptmass, read_options_ptmass
 public :: update_ptmass
 public :: calculate_mdot
 public :: ptmass_calc_enclosed_mass
 public :: ptmass_boundary_crossing

 ! settings affecting routines in module (read from/written to input file)
 integer, public :: icreate_sinks = 0
 real,    public :: rho_crit_cgs  = 1.e-10
 real,    public :: r_crit = 5.e-3
 real,    public :: h_acc  = 1.e-3
 real,    public :: f_acc  = 0.8
 real,    public :: h_soft_sinkgas  = 0.0
 real,    public :: h_soft_sinksink = 0.0
 real,    public :: r_merge_uncond  = 0.0    ! sinks will unconditionally merge if they touch
 real,    public :: r_merge_cond    = 0.0    ! sinks will merge if bound within this radius
 real,    public :: f_crit_override = 0.0    ! 1000.
 ! Note for above: if f_crit_override > 0, then will unconditionally make a sink when rho > f_crit_override*rho_crit_cgs
 ! This is a dangerous parameter since failure to form a sink might be indicative of another problem.
 ! This is a hard-coded parameter due to this danger, but will appear in the .in file if set > 0.

 ! additional public variables
 integer, public :: ipart_rhomax
 real,    public :: r_crit2,rho_crit
 real,    public :: r_merge2        = 0.0 ! initialise to prevent test failure
 real,    public :: r_merge_uncond2 = 0.0 ! initialise to prevent test failure
 real,    public :: r_merge_cond2   = 0.0 ! initialise to prevent test failure

 ! calibration of timestep control on sink-sink and sink-gas orbital integration
 ! this is hardwired because can be adjusted by changing C_force
 ! just means that with the default setting of C_force the orbits are accurate
 real, parameter :: dtfacphi  = 0.05
 real, parameter :: dtfacphi2 = dtfacphi*dtfacphi

 ! parameters to control output regarding sink particles
 logical, private, parameter :: record_created   = .false. ! verbose tracking of why sinks are not created
 logical, private            :: write_one_ptfile = .true.  ! default logical to determine if we are writing one or nptmass data files
 logical, private            :: l_crit_override  = .false. ! logical to determine the printing of f_crit_override to the .in file
 character(len=50), private  :: pt_prefix = 'Sink'
 character(len=50), private  :: pt_suffix = '00.sink'      ! will be overwritten to .ev for write_one_ptfile = .false.

 integer, public, parameter :: ndptmass = 13
 integer, public, parameter :: &
       idxmsi           =  1, &
       idymsi           =  2, &
       idzmsi           =  3, &
       idmsi            =  4, &
       idspinxsi        =  5, &
       idspinysi        =  6, &
       idspinzsi        =  7, &
       idvxmsi          =  8, &
       idvymsi          =  9, &
       idvzmsi          = 10, &
       idfxmsi          = 11, &
       idfymsi          = 12, &
       idfzmsi          = 13

 private

contains
!----------------------------------------------------------------
!+
!  if (tofrom==.true.)  Acceleration from/to gas particles due to sink particles;
!                       required in initial.F90 & step_leapfrog.F90 to update all accelerations
!  if (tofrom==.false.) Acceleration on gas due to sink particles (but not vice-versa);
!                       this is typically used to calculate phi (in compute_energies in
!                       energies.F90); in this case, fxi,fyi,fzi should be dummy input
!                       variables that do not affect the sink's motion.
!+
!----------------------------------------------------------------
subroutine get_accel_sink_gas(nptmass,xi,yi,zi,hi,xyzmh_ptmass,fxi,fyi,fzi,phi, &
                                       pmassi,fxyz_ptmass,dsdt_ptmass,fonrmax,dtphi2)
#ifdef FINVSQRT
 use fastmath,      only:finvsqrt
#endif
 use kernel,        only:kernel_softening,radkern
 use vectorutils,   only:unitvec
 use extern_geopot, only:get_geopot_force
 integer,           intent(in)    :: nptmass
 real,              intent(in)    :: xi,yi,zi,hi
 real,              intent(inout) :: fxi,fyi,fzi,phi
 real,              intent(in)    :: xyzmh_ptmass(nsinkproperties,nptmass)
 real,    optional, intent(in)    :: pmassi
 real,    optional, intent(inout) :: fxyz_ptmass(4,nptmass),dsdt_ptmass(3,nptmass)
 real,    optional, intent(out)   :: fonrmax,dtphi2
 real                             :: ftmpxi,ftmpyi,ftmpzi
 real                             :: dx,dy,dz,rr2,ddr,dr3,f1,f2,pmassj,J2,shat(3),Rsink
 real                             :: hsoft,hsoft1,hsoft21,q2i,qi,psoft,fsoft
 real                             :: fxj,fyj,fzj,dsx,dsy,dsz
 integer                          :: j
 logical                          :: tofrom
 !
 ! Determine if acceleration is from/to gas, or to gas
 !
 if (present(pmassi) .and. present(fxyz_ptmass) .and. present(fonrmax)) then
    tofrom  = .true.
    fonrmax = 0.
 else
    tofrom  = .false.
 endif

 ftmpxi = 0.  ! use temporary summation variable
 ftmpyi = 0.  ! (better for round-off, plus we need this bit of
 ftmpzi = 0.  ! the force to calculate the dtphi timestep)
 phi    = 0.
 f2     = 0.

 do j=1,nptmass
    dx     = xi - xyzmh_ptmass(1,j)
    dy     = yi - xyzmh_ptmass(2,j)
    dz     = zi - xyzmh_ptmass(3,j)
    pmassj = xyzmh_ptmass(4,j)
    hsoft  = xyzmh_ptmass(ihsoft,j)
    J2     = xyzmh_ptmass(iJ2,j)
    if (hsoft > 0.0) hsoft = max(hsoft,hi)
    if (pmassj < 0.0) cycle

    rr2    = dx*dx + dy*dy + dz*dz + epsilon(rr2)
#ifdef FINVSQRT
    ddr    = finvsqrt(rr2)
#else
    ddr    = 1./sqrt(rr2)
#endif
    dsx = 0.
    dsy = 0.
    dsz = 0.
    fxj = 0.
    fyj = 0.
    fzj = 0.
    if (rr2 < (radkern*hsoft)**2) then
       !
       ! if the sink particle is given a softening length, soften the
       ! force and potential if r < radkern*hsoft
       !
       hsoft1 = 1.0/hsoft
       hsoft21= hsoft1**2
       q2i    = rr2*hsoft21
       qi     = sqrt(q2i)
       call kernel_softening(q2i,qi,psoft,fsoft)  ! Note: psoft < 0

       ! acceleration of gas due to point mass particle
       f1     = pmassj*fsoft*hsoft21*ddr
       ftmpxi = ftmpxi - dx*f1
       ftmpyi = ftmpyi - dy*f1
       ftmpzi = ftmpzi - dz*f1
       phi    = phi + pmassj*psoft*hsoft1  ! potential (spline-softened)

       ! acceleration of sink from gas
       if (tofrom) f2 = pmassi*fsoft*hsoft21*ddr
    else
       ! no softening on the sink-gas interaction
       dr3  = ddr*ddr*ddr

       ! acceleration of gas due to point mass particle
       f1     = pmassj*dr3
       ftmpxi = ftmpxi - dx*f1
       ftmpyi = ftmpyi - dy*f1
       ftmpzi = ftmpzi - dz*f1
       phi    = phi    - pmassj*ddr      ! potential (GM/r)

       ! acceleration of sink from gas
       if (tofrom) f2 = pmassi*dr3

       ! additional accelerations due to oblateness
       if (abs(J2) > 0.) then
          shat = unitvec(xyzmh_ptmass(ispinx:ispinz,j))
          Rsink = xyzmh_ptmass(iReff,j)
          call get_geopot_force(dx,dy,dz,ddr,f1,Rsink,J2,shat,ftmpxi,ftmpyi,ftmpzi,phi,dsx,dsy,dsz,fxj,fyj,fzj)
       endif
    endif

    if (tofrom) then
       ! backreaction of gas onto sink
       fxyz_ptmass(1,j) = fxyz_ptmass(1,j) + dx*f2 + fxj*pmassi/pmassj
       fxyz_ptmass(2,j) = fxyz_ptmass(2,j) + dy*f2 + fyj*pmassi/pmassj
       fxyz_ptmass(3,j) = fxyz_ptmass(3,j) + dz*f2 + fzj*pmassi/pmassj

       ! backreaction torque of gas onto oblate sink
       dsdt_ptmass(1,j) = dsdt_ptmass(1,j) + pmassi*dsx
       dsdt_ptmass(2,j) = dsdt_ptmass(2,j) + pmassi*dsy
       dsdt_ptmass(3,j) = dsdt_ptmass(3,j) + pmassi*dsz

       ! timestep is sqrt(separation/force)
       fonrmax = max(f1,f2,fonrmax)
    endif
 enddo
 !
 ! external force timestep based on sqrt(phi)/accel
 !
 if (present(dtphi2)) then
    if (abs(phi) > epsilon(phi)) then
       f2     = ftmpxi*ftmpxi + ftmpyi*ftmpyi + ftmpzi*ftmpzi
       !dtphi is sqrt of this, but for optimisation we take the sqrt outside of the loop
       dtphi2 = dtfacphi2*abs(phi)/f2
    else
       dtphi2 = huge(dtphi2)
    endif
 endif
 !
 ! add temporary sums to existing force on gas particle
 !
 fxi = fxi + ftmpxi
 fyi = fyi + ftmpyi
 fzi = fzi + ftmpzi

end subroutine get_accel_sink_gas

!----------------------------------------------------------------
!+
!  Compute force on sink particles due to other sinks and
!  from external potentials
!+
!----------------------------------------------------------------
subroutine get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,phitot,dtsinksink,&
            iexternalforce,ti,merge_ij,merge_n,dsdt_ptmass)
#ifdef FINVSQRT
 use fastmath,       only:finvsqrt
#endif
 use externalforces, only:externalforce
 use extern_geopot,  only:get_geopot_force
 use kernel,         only:kernel_softening,radkern
 use vectorutils,    only:unitvec
 integer, intent(in)  :: nptmass
 real,    intent(in)  :: xyzmh_ptmass(nsinkproperties,nptmass)
 real,    intent(out) :: fxyz_ptmass(4,nptmass)
 real,    intent(out) :: phitot,dtsinksink
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: ti
 integer, intent(out) :: merge_ij(:),merge_n
 real,    intent(out) :: dsdt_ptmass(3,nptmass)
 real    :: xi,yi,zi,pmassi,pmassj,fxi,fyi,fzi,phii
 real    :: ddr,dx,dy,dz,rr2,rr2j,dr3,f1,f2
 real    :: hsoft1,hsoft21,q2i,qi,psoft,fsoft
 real    :: fextx,fexty,fextz,phiext !,hsofti
 real    :: fterm,pterm,potensoft0,dsx,dsy,dsz
 real    :: J2i,rsinki,shati(3)
 real    :: J2j,rsinkj,shatj(3)
 integer :: i,j

 dtsinksink = huge(dtsinksink)
 fxyz_ptmass(:,:) = 0.
 dsdt_ptmass(:,:) = 0.
 phitot   = 0.
 merge_n  = 0
 merge_ij = 0
 if (nptmass <= 1) return
 !
 !--get self-contribution to the potential if sink-sink softening is used
 !
 if (h_soft_sinksink > 0.) then
    hsoft1 = 1.0/h_soft_sinksink
    hsoft21= hsoft1**2
    call kernel_softening(0.,0.,potensoft0,fterm)
 else
    hsoft1 = 0.  ! to avoid compiler warnings
    hsoft21 = 0.
    potensoft0 = 0.
 endif
 !
 !--compute N^2 forces on point mass particles due to each other
 !
 !$omp parallel do default(none) &
 !$omp shared(nptmass,xyzmh_ptmass,fxyz_ptmass,merge_ij,r_merge2,dsdt_ptmass) &
 !$omp shared(iexternalforce,ti,h_soft_sinksink,potensoft0,hsoft1,hsoft21) &
 !$omp private(i,xi,yi,zi,pmassi,pmassj) &
 !$omp private(dx,dy,dz,rr2,rr2j,ddr,dr3,f1,f2) &
 !$omp private(fxi,fyi,fzi,phii,dsx,dsy,dsz) &
 !$omp private(fextx,fexty,fextz,phiext) &
 !$omp private(q2i,qi,psoft,fsoft) &
 !$omp private(fterm,pterm,J2i,J2j,shati,shatj,rsinki,rsinkj) &
 !$omp reduction(min:dtsinksink) &
 !$omp reduction(+:phitot,merge_n)
 do i=1,nptmass
    xi     = xyzmh_ptmass(1,i)
    yi     = xyzmh_ptmass(2,i)
    zi     = xyzmh_ptmass(3,i)
    pmassi = xyzmh_ptmass(4,i)
    !hsofti = xyzmh_ptmass(5,i)
    if (pmassi < 0.) cycle
    J2i    = xyzmh_ptmass(iJ2,i)

    fxi    = 0.
    fyi    = 0.
    fzi    = 0.
    phii   = 0.
    dsx    = 0.
    dsy    = 0.
    dsz    = 0.
    do j=1,nptmass
       if (i==j) cycle
       dx     = xi - xyzmh_ptmass(1,j)
       dy     = yi - xyzmh_ptmass(2,j)
       dz     = zi - xyzmh_ptmass(3,j)
       pmassj = xyzmh_ptmass(4,j)
       !hsoftj = xyzmh_ptmass(5,j)
       if (pmassj < 0.) cycle
       J2j = xyzmh_ptmass(iJ2,j)

       rr2  = dx*dx + dy*dy + dz*dz + epsilon(rr2)

#ifdef FINVSQRT
       ddr  = finvsqrt(rr2)
#else
       ddr  = 1./sqrt(rr2)
#endif

       if (rr2 < (radkern*h_soft_sinksink)**2) then
          !
          ! if the sink particle is given a softening length, soften the
          ! force and potential if r < radkern*h_soft_sinksink
          !
          q2i    = rr2*hsoft21
          qi     = sqrt(q2i)
          call kernel_softening(q2i,qi,psoft,fsoft)  ! Note: psoft < 0

          ! acceleration of sink1 from sink2
          fterm = fsoft*hsoft21*ddr
          f1    = pmassj*fterm
          fxi   = fxi - dx*f1
          fyi   = fyi - dy*f1
          fzi   = fzi - dz*f1
          pterm = psoft*hsoft1
          phii  = phii + pmassj*pterm ! potential (spline-softened)
       else
          ! no softening on the sink-sink interaction
          dr3   = ddr*ddr*ddr

          ! acceleration of sink1 from sink2
          f1    = pmassj*dr3
          fxi   = fxi - dx*f1
          fyi   = fyi - dy*f1
          fzi   = fzi - dz*f1
          pterm = -ddr
          phii  = phii + pmassj*pterm    ! potential (GM/r)

          ! additional acceleration due to oblateness of sink particles j and i
          if (abs(J2j) > 0.) then
             shatj = unitvec(xyzmh_ptmass(ispinx:ispinz,j))
             rsinkj = xyzmh_ptmass(iReff,j)
             call get_geopot_force(dx,dy,dz,ddr,f1,rsinkj,J2j,shatj,fxi,fyi,fzi,phii)
          endif
          if (abs(J2i) > 0.) then
             shati = unitvec(xyzmh_ptmass(ispinx:ispinz,i))
             rsinki = xyzmh_ptmass(iReff,i)
             call get_geopot_force(dx,dy,dz,ddr,f1,rsinki,J2i,shati,fxi,fyi,fzi,phii,dsx,dsy,dsz)
          endif
       endif
       if (rr2 < r_merge2) then
          if (merge_ij(i)==0) then
             merge_n = merge_n + 1
             merge_ij(i) = j
          else
             ! if we have already identified a nearby sink, replace the tag with the nearest sink
             dx   = xi - xyzmh_ptmass(1,merge_ij(i))
             dy   = yi - xyzmh_ptmass(2,merge_ij(i))
             dz   = zi - xyzmh_ptmass(3,merge_ij(i))
             rr2j = dx*dx + dy*dy + dz*dz + epsilon(rr2j)
             if (rr2 < rr2j) merge_ij(i) = j
          endif
       endif
    enddo
    phitot = phitot + 0.5*pmassi*phii  ! total potential (G M_1 M_2/r)

    !
    !--apply external forces
    !
    if (iexternalforce > 0) then
       call externalforce(iexternalforce,xi,yi,zi,0.,ti,fextx,fexty,fextz,phiext,ii=-i)
       fxi = fxi + fextx
       fyi = fyi + fexty
       fzi = fzi + fextz
       phii   = phii + phiext
       phitot = phitot + phiext
    endif
    !
    !--self-contribution to the potential if sink-sink softening is used
    !  Note: we do NOT add this for sink-sink interactions because the
    !  positions are assumed to be UNCORRELATED, hence the self-contribution
    !  is not important. Other particles (e.g. gas) are assumed to have
    !  correlated positions, so the self-contribution is important
    !
    !pterm = 0.5*pmassi*pmassi*potensoft0*hsoft1
    !phii = phii + pterm
    !phitot = phitot + pterm
    !
    !--store sink-sink forces (only)
    !
    fxyz_ptmass(1,i) = fxyz_ptmass(1,i) + fxi
    fxyz_ptmass(2,i) = fxyz_ptmass(2,i) + fyi
    fxyz_ptmass(3,i) = fxyz_ptmass(3,i) + fzi
    fxyz_ptmass(4,i) = fxyz_ptmass(4,i) + phii
    dsdt_ptmass(1,i) = dsdt_ptmass(1,i) + pmassi*dsx
    dsdt_ptmass(2,i) = dsdt_ptmass(2,i) + pmassi*dsy
    dsdt_ptmass(3,i) = dsdt_ptmass(3,i) + pmassi*dsz
 enddo
 !$omp end parallel do

 !
 !--sink-sink timestep based on sqrt(phi)/accel
 !  minimum is taken over all sink particles
 !
 do i=1,nptmass
    fxi  = fxyz_ptmass(1,i)
    fyi  = fxyz_ptmass(2,i)
    fzi  = fxyz_ptmass(3,i)
    phii = fxyz_ptmass(4,i)
    f2   = fxi*fxi + fyi*fyi + fzi*fzi
    !print*,'phi = ',phii,' accel = ',sqrt(f2)
    !
    !--we use an additional tolerance here on the sink-sink timestep
    !  so that with the default C_force of ~0.25 we get a few
    !  hundred steps per orbit
    !
    if (f2 > 0. .and. nptmass > 1) then
       dtsinksink = min(dtsinksink,dtfacphi*sqrt(abs(phii)/f2))
    endif
 enddo

end subroutine get_accel_sink_sink
!----------------------------------------------------------------
!+
!  Update position of sink particles if they cross the periodic boundary
!+
!----------------------------------------------------------------
subroutine ptmass_boundary_crossing(nptmass,xyzmh_ptmass)
 use boundary,  only:cross_boundary
 use mpidomain, only:isperiodic
 integer, intent(in)    :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 integer                :: i,ncross

 ncross = 0
 do i = 1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) call cross_boundary(isperiodic,xyzmh_ptmass(:,i),ncross)
 enddo

end subroutine ptmass_boundary_crossing

!----------------------------------------------------------------
!+
!  predictor step for the point masses
!  (called from inside a parallel section)
!+
!----------------------------------------------------------------
subroutine ptmass_predictor(nptmass,dt,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass)
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: dt
 real,    intent(inout) :: xyzmh_ptmass(nsinkproperties,nptmass)
 real,    intent(inout) :: vxyz_ptmass(3,nptmass)
 real,    intent(in)    :: fxyz_ptmass(4,nptmass),dsdt_ptmass(3,nptmass)
 real    :: vxhalfi,vyhalfi,vzhalfi
 integer :: i

 !$omp parallel do schedule(static) default(none) &
 !$omp shared(nptmass,dt,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass) &
 !$omp private(i,vxhalfi,vyhalfi,vzhalfi)
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       vxhalfi = vxyz_ptmass(1,i) + 0.5*dt*fxyz_ptmass(1,i)
       vyhalfi = vxyz_ptmass(2,i) + 0.5*dt*fxyz_ptmass(2,i)
       vzhalfi = vxyz_ptmass(3,i) + 0.5*dt*fxyz_ptmass(3,i)
       xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dt*vxhalfi
       xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dt*vyhalfi
       xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dt*vzhalfi
       vxyz_ptmass(1,i) = vxhalfi
       vxyz_ptmass(2,i) = vyhalfi
       vxyz_ptmass(3,i) = vzhalfi
       xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) + 0.5*dt*dsdt_ptmass(1,i)
       xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) + 0.5*dt*dsdt_ptmass(2,i)
       xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) + 0.5*dt*dsdt_ptmass(3,i)
    endif
 enddo
 !$omp end parallel do

end subroutine ptmass_predictor

!----------------------------------------------------------------
!+
!  corrector step for the point masses
!  (called from inside a parallel section)
!+
!----------------------------------------------------------------
subroutine ptmass_corrector(nptmass,dt,vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,dsdt_ptmass,iexternalforce)
 use externalforces, only:update_vdependent_extforce_leapfrog,is_velocity_dependent
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: dt
 real,    intent(inout) :: vxyz_ptmass(3,nptmass), xyzmh_ptmass(nsinkproperties,nptmass)
 real,    intent(in)    :: fxyz_ptmass(4,nptmass)
 real,    intent(in)    :: dsdt_ptmass(3,nptmass)
 integer, intent(in)    :: iexternalforce
 real :: vxhalfi,vyhalfi,vzhalfi
 real :: fxi,fyi,fzi,fextv(3)
 integer :: i

 !
 ! handle special case of velocity-dependent external forces
 ! in the leapfrog integrator
 !
 if (is_velocity_dependent(iexternalforce)) then
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,dsdt_ptmass,dt,nptmass,iexternalforce) &
    !$omp private(vxhalfi,vyhalfi,vzhalfi,fxi,fyi,fzi,fextv) &
    !$omp private(i)
    do i=1,nptmass
       if (xyzmh_ptmass(4,i) > 0.) then
          vxhalfi = vxyz_ptmass(1,i)
          vyhalfi = vxyz_ptmass(2,i)
          vzhalfi = vxyz_ptmass(3,i)
          fxi = fxyz_ptmass(1,i)
          fyi = fxyz_ptmass(2,i)
          fzi = fxyz_ptmass(3,i)
          call update_vdependent_extforce_leapfrog(iexternalforce,&
               vxhalfi,vyhalfi,vzhalfi,fxi,fyi,fzi,fextv,dt,&
               xyzmh_ptmass(1,i),xyzmh_ptmass(2,i),xyzmh_ptmass(3,i))
          fxi = fxi + fextv(1)
          fyi = fyi + fextv(2)
          fzi = fzi + fextv(3)
          vxyz_ptmass(1,i) = vxhalfi + 0.5*dt*fxi
          vxyz_ptmass(2,i) = vyhalfi + 0.5*dt*fyi
          vxyz_ptmass(3,i) = vzhalfi + 0.5*dt*fzi
          xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) + 0.5*dt*dsdt_ptmass(1,i)
          xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) + 0.5*dt*dsdt_ptmass(2,i)
          xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) + 0.5*dt*dsdt_ptmass(3,i)
       endif
    enddo
    !$omp end parallel do
 else
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,nptmass) &
    !$omp private(i)
    do i=1,nptmass
       if (xyzmh_ptmass(4,i) > 0.) then
          vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + 0.5*dt*fxyz_ptmass(1,i)
          vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + 0.5*dt*fxyz_ptmass(2,i)
          vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + 0.5*dt*fxyz_ptmass(3,i)
          xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) + 0.5*dt*dsdt_ptmass(1,i)
          xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) + 0.5*dt*dsdt_ptmass(2,i)
          xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) + 0.5*dt*dsdt_ptmass(3,i)
       endif
    enddo
    !$omp end parallel do
 endif

end subroutine ptmass_corrector

!----------------------------------------------------------------
!+
!  Determines if particle pairs are obscured by a sink particle; if
!  so, then they should not contribute to one another's properties
!  (except for gravity)
!  The input to this function assumes particle i is at 0,0,0.
!
!  The first version is implemented in sphNG, and the second method
!  takes the size of the sink particle in to account to determine if
!  it is obscuring the particle.  This is done by finding point c and
!  determining if it is within the accretion radius of the sink.
!  Point c is the point along ij that is closest to the sink particle
!  (therefore line{cs} \perpto line{ij}).
!  Tests shows that sphNG version is more stable than the geometrical
!  version, which is more stable than nothing.
!+
!----------------------------------------------------------------
logical function ptmass_not_obscured(xj,yj,zj,xsink,ysink,zsink,r_sink)
 real, intent(in) :: xj,yj,zj,xsink,ysink,zsink,r_sink
 real             :: xc,yc,zc
 real             :: sep1_ij,sep_ic,sep2_cs
 !
 !  sphNG method
 ptmass_not_obscured = .false.                            ! == add_contribution from force.F90
 return
 !
 !  Full geometrical method
 sep1_ij = 1.0/sqrt(xj*xj + yj*yj + zj*zj)                ! 1.0/separation between i & j
 sep_ic  = (xj*xsink + yj*ysink + zj*zsink )*sep1_ij      ! separation between i & c = cos(a)*sep_is, where cos(a) = (ij*isink)/(|ij||isink|)
 xc      = xj*sep1_ij * sep_ic                            ! x-coordinate of point c
 yc      = yj*sep1_ij * sep_ic                            ! y-coordinate of point c
 zc      = zj*sep1_ij * sep_ic                            ! z-coordinate of point c
 sep2_cs = (xc-xsink)**2 + (yc-ysink)**2 + (zc-zsink)**2  ! separation^2 between c & sink
 if (sep2_cs < r_sink*r_sink ) then                       ! determine if the closest point along ij is within twice the sink radius
    ptmass_not_obscured = .false.
 else
    ptmass_not_obscured = .true.
 endif
 !
end function ptmass_not_obscured
!----------------------------------------------------------------
!+
!  accrete particles onto point masses
!+
!----------------------------------------------------------------
!----------------------------------------------------------------
! Routine updated by CJN 12/06/11
! and again by CJN on 30/03/14
! and again by JHW on 09/12/14
!
! Also should include a thermal energy component of sinks for
! calculating conserved quantities - otherwise accreted particle
! energy is thrown away.
!
! Now includes check to ensure that the particle is actually bound
! to the point mass and not just passing through its neighbourhood:
!      (a) specific angular momentum of particle must be less than that
!          required for it to form a circular orbit at hacc
!      (b) particle must be bound
!      (c) particle must be more bound to current point mass than any other
! Since the order of accretion should not matter, the sink's original
! characteristics will be used in the checks.  However, the most updated
! values will be used to update the sink's characteristics since the order
! in which particles is added is irrelevant.
!----------------------------------------------------------------
subroutine ptmass_accrete(is,nptmass,xi,yi,zi,hi,vxi,vyi,vzi,fxi,fyi,fzi, &
                          itypei,pmassi,xyzmh_ptmass,vxyz_ptmass,accreted, &
                          dptmass,time,facc,nbinmax,ibin_wakei,nfaili)

!$ use omputils, only:ipart_omp_lock
 use part,       only: ihacc
 use kernel,     only: radkern2
 use io,         only: iprint,iverbose,fatal
 use io_summary, only: iosum_ptmass,maxisink,print_acc
 integer,           intent(in)    :: is,nptmass,itypei
 real,              intent(in)    :: xi,yi,zi,pmassi,vxi,vyi,vzi,fxi,fyi,fzi,time,facc
 real,              intent(inout) :: hi
 real,              intent(in)    :: xyzmh_ptmass(nsinkproperties,nptmass)
 real,              intent(in)    :: vxyz_ptmass(3,nptmass)
 logical,           intent(out)   :: accreted
 real,              intent(inout) :: dptmass(:,:)
 integer(kind=1),   intent(in)    :: nbinmax
 integer(kind=1),   intent(inout) :: ibin_wakei
 integer, optional, intent(out)   :: nfaili
 integer            :: i,ifail
 real               :: dx,dy,dz,r2,dvx,dvy,dvz,v2,hacc
 logical, parameter :: iofailreason=.false.
 integer            :: j
 real               :: mpt,drdv,angmom2,angmomh2,epart,dxj,dyj,dzj,dvxj,dvyj,dvzj,rj2,vj2,epartj
 logical            :: mostbound

 accreted = .false.
 ifail    = 0
 !
 ! Verify particle is 'accretable'
 if (.not. is_accretable(itypei) ) then
    if (present(nfaili)) nfaili = 5
    if (iverbose >= 1 .and. iofailreason) &
       write(iprint,"(/,a)") 'ptmass_accrete: FAILED: particle is not an accretable type'
    return
 endif
 !
 sinkloop : do i=is,nptmass
    hacc = xyzmh_ptmass(ihacc,i)
    mpt  = xyzmh_ptmass(4,i)
    if (mpt < 0.) cycle
    dx = xi - xyzmh_ptmass(1,i)
    dy = yi - xyzmh_ptmass(2,i)
    dz = zi - xyzmh_ptmass(3,i)
    r2 = dx*dx + dy*dy + dz*dz
    dvx = vxi - vxyz_ptmass(1,i)
    dvy = vyi - vxyz_ptmass(2,i)
    dvz = vzi - vxyz_ptmass(3,i)
    v2 = dvx*dvx + dvy*dvy + dvz*dvz
!
!  See if particle passes conditions to be accreted
!
    if (r2 < (facc*hacc)**2) then
       ! accrete indiscriminately
       accreted = .true.
       ifail    = -1
    elseif (r2 < hacc**2) then
       ibin_wakei = nbinmax
       drdv = dx*dvx + dy*dvy + dz*dvz
       ! compare specific angular momentum
       angmom2  = r2*v2 - drdv*drdv
       angmomh2 = mpt*hacc
       if (angmom2 < angmomh2) then
          ! check if bound
          epart = 0.5*v2 - mpt/sqrt(r2)
          if (epart < 0.) then
             ! check to ensure it is most bound to this particle
             mostbound = .true.
             j = 1
             do while (mostbound .and. j <= nptmass)
                if (j /= i) then
                   dxj    = xi - xyzmh_ptmass(1,j)
                   dyj    = yi - xyzmh_ptmass(2,j)
                   dzj    = zi - xyzmh_ptmass(3,j)
                   rj2    = dxj*dxj + dyj*dyj + dzj*dzj
                   dvxj   = vxi - vxyz_ptmass(1,j)
                   dvyj   = vyi - vxyz_ptmass(2,j)
                   dvzj   = vzi - vxyz_ptmass(3,j)
                   vj2    = dvxj*dvxj + dvyj*dvyj + dvzj*dvzj
                   epartj = 0.5*vj2 - xyzmh_ptmass(4,j)/sqrt(rj2)
                   if (epartj < epart) mostbound = .false.
                endif
                j = j + 1
             enddo
             if ( mostbound ) then
                accreted = .true.
                ifail = -2
             else
                ifail = 4
             endif
          else
             ifail = 3
          endif
       else
          ifail = 2
       endif
    else
       ifail = 1
       if (r2 < radkern2*hi*hi) ibin_wakei = nbinmax
    endif
    if (iverbose >= 1 .and. iofailreason) then
       !--Forced off since output will be unreasonably large
       select case(ifail)
       case(4)
          write(iprint,"(/,a)") 'ptmass_accrete: FAILED: particle is not most bound to this sink'
       case(3)
          write(iprint,"(/,a,Es9.2)") 'ptmass_accrete: FAILED: particle is not bound: e = ',epart
       case(2)
          write(iprint,"(/,a,Es9.2,a,Es9.2)") 'ptmass_accrete: FAILED: angular momentum is too large: ' &
                                              ,angmom2,' > ',angmomh2
       case(1)
          write(iprint,"(/,a)") 'ptmass_accrete: FAILED: r2 > hacc**2'
       case(-1)
          write(iprint,"(/,a)") 'ptmass_accrete: PASSED indiscriminately: particle will be accreted'
       case(-2)
          write(iprint,"(/,a)") 'ptmass_accrete: PASSED: particle will be accreted'
       case default
          write(iprint,"(/,a)") 'ptmass_accrete: FAILED: unknown reason'
       end select
    endif
    if (present(nfaili)) nfaili = ifail
!
! if accreted==true, then checks all passed => accrete particle
!
    if ( accreted ) then
!$     call omp_set_lock(ipart_omp_lock(i))

! Set new position for the sink particles
       dptmass(idxmsi,i) = dptmass(idxmsi,i) + xi*pmassi
       dptmass(idymsi,i) = dptmass(idymsi,i) + yi*pmassi
       dptmass(idzmsi,i) = dptmass(idzmsi,i) + zi*pmassi

! Set new mass and increment accreted mass
       dptmass(idmsi,i) = dptmass(idmsi,i) + pmassi

! Set new spin angular momentum; this component is the angular momentum
! of the accreted particles about the origin
       dptmass(idspinxsi,i) = dptmass(idspinxsi,i) + pmassi*(yi*vzi - zi*vyi)
       dptmass(idspinysi,i) = dptmass(idspinysi,i) + pmassi*(zi*vxi - xi*vzi)
       dptmass(idspinzsi,i) = dptmass(idspinzsi,i) + pmassi*(xi*vyi - yi*vxi)

! Set new velocities for the sink particles
       dptmass(idvxmsi,i) = dptmass(idvxmsi,i) + vxi*pmassi
       dptmass(idvymsi,i) = dptmass(idvymsi,i) + vyi*pmassi
       dptmass(idvzmsi,i) = dptmass(idvzmsi,i) + vzi*pmassi

! Set new accelerations for the sink particles
       dptmass(idfxmsi,i) = dptmass(idfxmsi,i) + fxi*pmassi
       dptmass(idfymsi,i) = dptmass(idfymsi,i) + fyi*pmassi
       dptmass(idfzmsi,i) = dptmass(idfzmsi,i) + fzi*pmassi

! Track values for summary
       print_acc = .true.
       if (nptmass > maxisink) then
          iosum_ptmass(1,1) = iosum_ptmass(1,1) + 1
          if (ifail == -1) iosum_ptmass(2,1) = iosum_ptmass(2,1) + 1
       else
          iosum_ptmass(1,i) = iosum_ptmass(1,i) + 1
          if (ifail == -1) iosum_ptmass(2,i) = iosum_ptmass(2,i) + 1
       endif

!$     call omp_unset_lock(ipart_omp_lock(i))
       hi = -abs(hi)

! avoid possibility that two sink particles try to accrete the same gas particle by exiting the loop
       exit sinkloop
    endif
 enddo sinkloop

end subroutine ptmass_accrete

!-----------------------------------------------------------------------
!+
!  Update ptmass position, spin, velocity, acceleration, and mass
!  of sink particles once all particles are accreted
!  Regarding Spin Angular Momentum, S:
!  If calculated serially, then for particle i,
!  S = S + (m_i M_sink)/(M_sink+m_i) [(r_i-r_sink) x (v_i-v_sink)]
!  This assumes that the sink properties will be updated before the next
!  accretion event
!  To be compatible with parallel construction, this is equivalent to
!  S = S + sum_i (L_i) + L_{sink,before all accretion} - L_{sink, after all accretion}
!  where the angular momenta are calculated about the origin, x=y=z=0.
!  The latter is used; sum_i (L_i) is calculated in ptmass_accrete, and
!  two angular momentum terms are calculated here.
!+
!-----------------------------------------------------------------------
subroutine update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)
 real,    intent(in)    :: dptmass(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 integer, intent(in)    :: nptmass

 real                   :: newptmass(nptmass),newptmass1(nptmass)

 ! Add angular momentum of sink particle using old properties (taken about the origin)
 xyzmh_ptmass(ispinx,1:nptmass) =xyzmh_ptmass(ispinx,1:nptmass)+xyzmh_ptmass(4,1:nptmass) &
                                *(xyzmh_ptmass(2,1:nptmass)*vxyz_ptmass(3,1:nptmass)      &
                                - xyzmh_ptmass(3,1:nptmass)*vxyz_ptmass(2,1:nptmass))
 xyzmh_ptmass(ispiny,1:nptmass) =xyzmh_ptmass(ispiny,1:nptmass)+xyzmh_ptmass(4,1:nptmass) &
                                *(xyzmh_ptmass(3,1:nptmass)*vxyz_ptmass(1,1:nptmass)      &
                                - xyzmh_ptmass(1,1:nptmass)*vxyz_ptmass(3,1:nptmass))
 xyzmh_ptmass(ispinz,1:nptmass) =xyzmh_ptmass(ispinz,1:nptmass)+xyzmh_ptmass(4,1:nptmass) &
                                *(xyzmh_ptmass(1,1:nptmass)*vxyz_ptmass(2,1:nptmass)      &
                                - xyzmh_ptmass(2,1:nptmass)*vxyz_ptmass(1,1:nptmass))
 ! Calculate new masses
 newptmass(1:nptmass)           =xyzmh_ptmass(4,1:nptmass)+dptmass(idmsi,1:nptmass)
 newptmass1(1:nptmass)          =1./newptmass(1:nptmass)
 ! Update position and accreted mass
 xyzmh_ptmass(1,1:nptmass)      =(dptmass(idxmsi,1:nptmass)+xyzmh_ptmass(1,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 xyzmh_ptmass(2,1:nptmass)      =(dptmass(idymsi,1:nptmass)+xyzmh_ptmass(2,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 xyzmh_ptmass(3,1:nptmass)      =(dptmass(idzmsi,1:nptmass)+xyzmh_ptmass(3,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 xyzmh_ptmass(imacc, 1:nptmass) = xyzmh_ptmass(imacc,1:nptmass)+dptmass(idmsi,    1:nptmass)
 ! Add angular momentum contribution from the gas particles
 xyzmh_ptmass(ispinx,1:nptmass) =xyzmh_ptmass(ispinx,1:nptmass)+dptmass(idspinxsi,1:nptmass)
 xyzmh_ptmass(ispiny,1:nptmass) =xyzmh_ptmass(ispiny,1:nptmass)+dptmass(idspinysi,1:nptmass)
 xyzmh_ptmass(ispinz,1:nptmass) =xyzmh_ptmass(ispinz,1:nptmass)+dptmass(idspinzsi,1:nptmass)
 ! Update velocity, force, and final mass
 vxyz_ptmass(1,1:nptmass)       =(dptmass(idvxmsi,1:nptmass)+vxyz_ptmass(1,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 vxyz_ptmass(2,1:nptmass)       =(dptmass(idvymsi,1:nptmass)+vxyz_ptmass(2,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 vxyz_ptmass(3,1:nptmass)       =(dptmass(idvzmsi,1:nptmass)+vxyz_ptmass(3,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 fxyz_ptmass(1,1:nptmass)       =(dptmass(idfxmsi,1:nptmass)+fxyz_ptmass(1,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 fxyz_ptmass(2,1:nptmass)       =(dptmass(idfymsi,1:nptmass)+fxyz_ptmass(2,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 fxyz_ptmass(3,1:nptmass)       =(dptmass(idfzmsi,1:nptmass)+fxyz_ptmass(3,1:nptmass)*xyzmh_ptmass(4,1:nptmass))*newptmass1
 xyzmh_ptmass(4,1:nptmass)      =newptmass(1:nptmass)
 ! Subtract angular momentum of sink particle using new properties (taken about the origin)
 xyzmh_ptmass(ispinx,1:nptmass) =xyzmh_ptmass(ispinx,1:nptmass)-xyzmh_ptmass(4,1:nptmass) &
                                *(xyzmh_ptmass(2,1:nptmass)*vxyz_ptmass(3,1:nptmass)      &
                                - xyzmh_ptmass(3,1:nptmass)*vxyz_ptmass(2,1:nptmass))
 xyzmh_ptmass(ispiny,1:nptmass) =xyzmh_ptmass(ispiny,1:nptmass)-xyzmh_ptmass(4,1:nptmass) &
                                *(xyzmh_ptmass(3,1:nptmass)*vxyz_ptmass(1,1:nptmass)      &
                                - xyzmh_ptmass(1,1:nptmass)*vxyz_ptmass(3,1:nptmass))
 xyzmh_ptmass(ispinz,1:nptmass) =xyzmh_ptmass(ispinz,1:nptmass)-xyzmh_ptmass(4,1:nptmass) &
                                *(xyzmh_ptmass(1,1:nptmass)*vxyz_ptmass(2,1:nptmass)      &
                                - xyzmh_ptmass(2,1:nptmass)*vxyz_ptmass(1,1:nptmass))

end subroutine update_ptmass

!-------------------------------------------------------------------------
!+
! Subroutine to automatically create and insert a sink particle
! once certain conditions are met
!
! Conditions are given in section 2.2.2 of BBP95 and in the Phantom paper
!+
!-------------------------------------------------------------------------
subroutine ptmass_create(nptmass,npart,itest,xyzh,vxyzu,fxyzu,fext,divcurlv,poten,&
                         massoftype,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,time)
 use part,   only:ihacc,ihsoft,igas,iamtype,get_partinfo,iphase,iactive,maxphase,rhoh, &
                  ispinx,ispiny,ispinz,fxyz_ptmass_sinksink,eos_vars,igasP,igamma
 use dim,    only:maxp,maxneigh,maxvxyzu,maxptmass,ind_timesteps
 use kdtree, only:getneigh
 use kernel, only:kernel_softening,radkern
 use io,     only:id,iprint,fatal,iverbose,nprocs
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 use part,     only:ibin,ibin_wake
 use linklist, only:getneigh_pos,ifirstincell,listneigh=>listneigh_global
 use eos,           only:gamma
 use eos_barotropic,only:gamma_barotropic
 use eos_piecewise, only:gamma_pwp
 use options,  only:ieos
 use units,    only:unit_density
 use io_summary, only:summary_variable_rhomax,summary_ptmass_fail, &
                      inosink_notgas,inosink_divv,inosink_h,inosink_active, &
                      inosink_therm,inosink_grav,inosink_Etot,inosink_poten,inosink_max
 use mpiutils, only:reduceall_mpi,bcast_mpi,reduceloc_mpi
 integer,         intent(inout) :: nptmass
 integer,         intent(in)    :: npart,itest
 real,            intent(inout) :: xyzh(:,:)
 real,            intent(in)    :: vxyzu(:,:),fxyzu(:,:),fext(:,:),massoftype(:)
 real(4),         intent(in)    :: divcurlv(:,:),poten(:)
 real,            intent(inout) :: xyzmh_ptmass(:,:)
 real,            intent(inout) :: vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 real,            intent(in)    :: time
 integer(kind=1)    :: iphasei,ibin_wakei,ibin_itest
 integer            :: nneigh
 integer, parameter :: maxcache      = 12000
 integer, parameter :: nneigh_thresh = 1024 ! approximate epot if neigh>neigh_thresh; (-ve for off)
 real, save :: xyzcache(maxcache,3)
 real    :: dptmass(ndptmass,nptmass+1)
 real    :: xi,yi,zi,hi,hi1,hi21,xj,yj,zj,hj1,hj21,xk,yk,zk,hk1
 real    :: rij2,rik2,rjk2,dx,dy,dz
 real    :: vxi,vyi,vzi,dv2,dvx,dvy,dvz,rhomax
 real    :: alpha_grav,alphabeta_grav,radxy2,radxz2,radyz2
 real    :: etot,epot,ekin,etherm,erot,erotx,eroty,erotz
 real    :: rcrossvx,rcrossvy,rcrossvz,fxj,fyj,fzj
 real    :: pmassi,pmassj,pmassk,rhoj
 real    :: q2i,qi,psofti,psoftj,psoftk,fsoft,epot_mass,epot_rad,pmassgas1
 real    :: hcheck,hcheck2,f_acc_local
 real(4) :: divvi,potenj_min,poteni
 integer :: ifail,nacc,j,k,n,nk,itype,itypej,itypek,ifail_array(inosink_max),id_rhomax,nneigh_act
 logical :: accreted,iactivej,isgasj,isdustj,calc_exact_epot,ForceCreation

 ifail       = 0
 ifail_array = 0
 poteni      = 0._4
 potenj_min  = huge(poteni)
!
! find the location of the maximum density across
! all MPI threads
!
 rhomax = 0.
 if (itest > 0 .and. itest <= npart) then
    iphasei = iphase(itest)
    itype   = iamtype(iphasei)
    rhomax  = rhoh(xyzh(4,itest),massoftype(itype))
 endif
 call reduceloc_mpi('max',rhomax,id_rhomax)
 ForceCreation = (f_crit_override > 0. .and. rhomax > f_crit_override*rho_crit)
!
! get properties of particle on the thread
! where it belongs
!
 if (id == id_rhomax) then
    if (itest < 0 .or. itest > npart) call fatal('ptmass','index out of range testing for sink creation')
    if (ForceCreation) then
       write(iprint,"(/,1x,a,2(Es18.6,a))") 'ptmass_create: WARNING! rhomax = ',rhomax*unit_density,' > ', &
                                             f_crit_override*rho_crit_cgs,' = f_crit_override*rho_crit  (cgs units)'
       write(iprint,"(/,1x,a)")             'ptmass_create: WARNING! Forcing sink formation despite tests not passing!'
    endif
    xi = xyzh(1,itest)
    yi = xyzh(2,itest)
    zi = xyzh(3,itest)
    hi = xyzh(4,itest)
    vxi = vxyzu(1,itest)
    vyi = vxyzu(2,itest)
    vzi = vxyzu(3,itest)
    iphasei = iphase(itest)
    divvi = divcurlv(1,itest)
    if (ind_timesteps) ibin_itest = ibin(itest)
    if (gravity) poteni = poten(itest)
 endif
!
! broadcast properties of the particle being tested to all threads
!
 call bcast_mpi(xi,id_rhomax)
 call bcast_mpi(yi,id_rhomax)
 call bcast_mpi(zi,id_rhomax)
 call bcast_mpi(hi,id_rhomax)
 call bcast_mpi(vxi,id_rhomax)
 call bcast_mpi(vyi,id_rhomax)
 call bcast_mpi(vzi,id_rhomax)
 call bcast_mpi(iphasei,id_rhomax)
 call bcast_mpi(divvi,id_rhomax)
 if (ind_timesteps) call bcast_mpi(ibin_itest,id_rhomax)
 if (gravity) call bcast_mpi(poteni,id_rhomax)
 !
 ! determine radius in which to check the criteria
 !
 hcheck      = radkern*hi               ! = h_acc in previous versions of Phantom; current method is faster
 f_acc_local = max(f_acc,hcheck/h_acc)  ! = 1.0   in previous versions of Phantom; current method is faster
 hcheck2     = hcheck*hcheck
 !
 ! initialise variables
 !
 hi1  = 1.0/hi
 hi21 = hi1**2
 if (maxphase==maxp) then
    itype = iamtype(iphasei)
 else
    itype = igas
 endif
 pmassi = massoftype(itype)
 pmassj = massoftype(igas)
 pmassk = pmassj
 itypej = igas
 itypek = igas
 iactivej = .true.
 pmassgas1 = 1.0/pmassj

 if (id==id_rhomax) call summary_variable_rhomax(itest,rhoh(hi,pmassi)*real(unit_density),iprint,nptmass)

 if (iverbose >= 1 .and. id==id_rhomax) &
    write(iprint,"(a,i10,a,i2,a)",advance='no') &
     ' ptmass_create: Testing particle i=',itest,' on thread ',id,' for ptmass creation...'

 ! CHECK 0: make sure particle is a gas particle (sanity check, should be unnecessary)
 if (.not. is_accretable(itype)) then
    if (iverbose >= 1) write(iprint,"(/,1x,a)") 'ptmass_create: FAILED because not a gas particle'
    call summary_ptmass_fail(inosink_notgas)
    if (.not. record_created) return
    ifail_array(inosink_notgas) = 1
 endif

 ! CHECK 1: divv < 0
 if (divvi > 0._4) then
    if (iverbose >= 1) write(iprint,"(/,1x,a)") 'ptmass_create: FAILED because div v > 0'
    call summary_ptmass_fail(inosink_divv)
    if (.not. record_created .and. .not.ForceCreation) return
    ifail_array(inosink_divv) = 1
 endif

 ! CHECK 2: 2h < h_acc
 if (hi > 0.5*h_acc) then
    if (iverbose >= 1) write(iprint,"(/,1x,2(a,es10.3),a)") 'ptmass_create: FAILED because 2h > h_acc (',2*hi,' > ',h_acc,')'
    call summary_ptmass_fail(inosink_h)
    if (.not. record_created) return
    ifail_array(inosink_h) = 1
 endif

 ekin   = 0.
 epot   = -epsilon(epot)
 etherm = 0.
 erot   = 0.
 erotx  = 0.
 eroty  = 0.
 erotz  = 0.
 epot_mass  = 0.
 epot_rad   = 0.
 nneigh_act = 0

 ! CHECK 3: all neighbours are all active ( & perform math for checks 4-6)
 ! find neighbours within the checking radius of hcheck
 call getneigh_pos((/xi,yi,zi/),0.,hcheck,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
 ! determine if we should approximate epot
 calc_exact_epot = .true.
 if ((nneigh_thresh > 0 .and. nneigh > nneigh_thresh) .or. (nprocs > 1)) calc_exact_epot = .false.
!$omp parallel default(none) &
!$omp shared(nprocs) &
!$omp shared(maxp,maxphase) &
!$omp shared(nneigh,listneigh,xyzh,xyzcache,vxyzu,massoftype,iphase,pmassgas1,calc_exact_epot,hcheck2,eos_vars) &
!$omp shared(itest,id,id_rhomax,ifail,xi,yi,zi,hi,vxi,vyi,vzi,hi1,hi21,itype,pmassi,ieos,gamma,poten) &
#ifdef PERIODIC
!$omp shared(dxbound,dybound,dzbound) &
#endif
!$omp shared(ibin_wake,ibin_itest) &
!$omp private(n,j,xj,yj,zj,hj1,hj21,psoftj,rij2,nk,k,xk,yk,zk,hk1,psoftk,rjk2,psofti,rik2) &
!$omp private(dx,dy,dz,dvx,dvy,dvz,dv2,isgasj,isdustj) &
!$omp private(rhoj,q2i,qi,fsoft,rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2) &
!$omp firstprivate(pmassj,pmassk,itypej,iactivej,itypek) &
!$omp reduction(+:nneigh_act,ekin,erotx,eroty,erotz,etherm,epot,epot_mass,epot_rad) &
!$omp reduction(min:potenj_min)
!$omp do
 over_neigh: do n=1,nneigh
    j = listneigh(n)
    !
    ! get mass and particle type to immediately determine if active and accretable
    if (maxphase==maxp) then
       call get_partinfo(iphase(j),iactivej,isgasj,isdustj,itypej)
       pmassj = massoftype(itypej)
       if (.not. is_accretable(itypej) ) cycle over_neigh ! Verify particle is 'accretable'
    endif

    if (n <= maxcache) then
       xj = xyzcache(n,1)
       yj = xyzcache(n,2)
       zj = xyzcache(n,3)
    else
       xj = xyzh(1,j)
       yj = xyzh(2,j)
       zj = xyzh(3,j)
    endif
    dx = xi - xj
    dy = yi - yj
    dz = zi - zj
#ifdef PERIODIC
    if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
    if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
    if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
    rij2 = dx*dx + dy*dy + dz*dz
    if (rij2 < hcheck2) then

       if (ind_timesteps) then
          ibin_wake(j) = max(ibin_wake(j),ibin_itest)
          if (.not.iactivej .or. ifail==inosink_active) then
             ifail = inosink_active
             cycle over_neigh
          endif
       endif

       nneigh_act = nneigh_act + 1

       dvx = vxi - vxyzu(1,j)
       dvy = vyi - vxyzu(2,j)
       dvz = vzi - vxyzu(3,j)

       hj1  = 1.0/xyzh(4,j)
       hj21 = hj1**2

       ! kinetic energy
       dv2  = dvx*dvx + dvy*dvy + dvz*dvz
       ekin = ekin + pmassj*dv2

       ! rotational energies around each axis
       rcrossvx = (dy*dvz - dz*dvy)
       rcrossvy = (dz*dvx - dx*dvz)
       rcrossvz = (dx*dvy - dy*dvx)

       radxy2 = dx*dx + dy*dy
       radyz2 = dy*dy + dz*dz
       radxz2 = dx*dx + dz*dz

       if (radyz2 > 0.) erotx = erotx + pmassj*rcrossvx*rcrossvx/radyz2
       if (radxz2 > 0.) eroty = eroty + pmassj*rcrossvy*rcrossvy/radxz2
       if (radxy2 > 0.) erotz = erotz + pmassj*rcrossvz*rcrossvz/radxy2

       ! thermal energy (for gas only)
       if (itypej==igas) then
          rhoj = rhoh(xyzh(4,j),pmassj)
          if (maxvxyzu >= 4) then
             etherm = etherm + pmassj*vxyzu(4,j)
          else
             if (ieos==2 .and. gamma > 1.001) then
                etherm = etherm + pmassj*(eos_vars(igasP,j)/rhoj)/(gamma - 1.)
             elseif (ieos==5 .and. gamma > 1.001) then
                etherm = etherm + pmassj*(eos_vars(igasP,j)/rhoj)/(eos_vars(igamma,j) - 1.)
             elseif (ieos==8) then
                etherm = etherm + pmassj*(eos_vars(igasP,j)/rhoj)/(gamma_barotropic(rhoj) - 1.)
             elseif (ieos==9) then
                etherm = etherm + pmassj*(eos_vars(igasP,j)/rhoj)/(gamma_pwp(rhoj) - 1.)
             else
                etherm = etherm + pmassj*1.5*(eos_vars(igasP,j)/rhoj)
             endif
          endif
       endif

       ! gravitational potential energy of clump
       if (gravity) then
          potenj_min = min(potenj_min,poten(j))
          if (calc_exact_epot) then
             if (nprocs > 1) call fatal('ptmass_create', 'cannot use calc_exact_epot with MPI')
             ! Calculate potential energy exactly
             !
             ! add contribution of i-j (since, e.g., rij2 is already calculated)
             !
             q2i    = rij2*hi21
             qi     = sqrt(q2i)
             call kernel_softening(q2i,qi,psofti,fsoft)
             q2i    = rij2*hj21
             qi     = sqrt(q2i)
             call kernel_softening(q2i,qi,psoftj,fsoft)
             epot   = epot + 0.5*pmassi*pmassj*(psofti*hi1 + psoftj*hj1)
             !
             ! add contribution of k-j for all k >= j (to avoid double counting, but include self-contribution)
             !
             over_neigh_k: do nk=n,nneigh
                k = listneigh(nk)
                if (k==itest .and. id==id_rhomax) cycle over_neigh_k ! contribution already added
                if (maxphase==maxp) then
                   itypek = iamtype(iphase(k))
                   pmassk = massoftype(itypek)
                   if (.not. is_accretable(itypek) ) cycle over_neigh_k
                endif

                if (nk <= maxcache) then
                   xk = xyzcache(nk,1)
                   yk = xyzcache(nk,2)
                   zk = xyzcache(nk,3)
                else
                   xk = xyzh(1,k)
                   yk = xyzh(2,k)
                   zk = xyzh(3,k)
                endif
                dx = xi - xk
                dy = yi - yk
                dz = zi - zk
#ifdef PERIODIC
                if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
                if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
                if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
                rik2 = dx*dx + dy*dy + dz*dz
                if (rik2 < hcheck2) then
                   dx = xj - xk
                   dy = yj - yk
                   dz = zj - zk
                   hk1 = 1.0/xyzh(4,k)
#ifdef PERIODIC
                   if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
                   if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
                   if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
                   rjk2   = dx*dx + dy*dy + dz*dz
                   ! Since neither j or k is dominant, sum half the contribution at each particle
                   ! Due to the construction of the loop, both j & k contributions are added here
                   q2i    = rjk2*hj21
                   qi     = sqrt(q2i)
                   call kernel_softening(q2i,qi,psoftj,fsoft)
                   q2i    = rjk2*hk1**2
                   qi     = sqrt(q2i)
                   call kernel_softening(q2i,qi,psoftk,fsoft)
                   epot   = epot + 0.5*pmassj*pmassk*(psoftj*hj1 + psoftk*hk1)
                endif
             enddo over_neigh_k
          else
             ! Calculate mass to approximate potential energy
             epot_mass = epot_mass + pmassj*pmassgas1             ! to avoid rounding errors
             ! Calculate the mass-weighted average distance to the sink particle candidate
             epot_rad  = epot_rad  + pmassj*pmassgas1*sqrt(rij2)
          endif
       endif
    endif
 enddo over_neigh
!$omp enddo
!$omp end parallel

 if (.not. calc_exact_epot) then
    epot_mass = reduceall_mpi('+', epot_mass)
    epot_rad  = reduceall_mpi('+', epot_rad)
    epot_mass = epot_mass + pmassi*pmassgas1  !self-contribution of the candidate particle
 endif
 !
 !--Update tracking array & reset ifail if required
 !  Note that if ifail_array(inosink_notgas,inosink_divv,inosink_h)==1 and record_created==.false.,
 !  this subroutine will already have been exited, and this loop will never be reached
 if ( record_created .or. ForceCreation) then
    if ( ifail==inosink_active ) then
       ifail_array(inosink_active) = 1
    elseif (ifail_array(inosink_notgas)==1) then
       ifail = inosink_notgas
    elseif (ifail_array(inosink_divv)==1) then
       ifail = inosink_divv
    elseif (ifail_array(inosink_h)==1) then
       ifail = inosink_h
    endif
 endif
 !
 ! communicate failure on any MPI thread to all threads
 !
 ifail = int(reduceall_mpi('max',ifail))
 ifail_array = int(reduceall_mpi('max',ifail_array))
 !
 ! Continue checks (non-sensical for ifail==1 since energies not completely calculated)
 !
 if (ifail==0 .or. ((ForceCreation .or. record_created) .and. ifail_array(inosink_active) == 0) ) then
    ! finish computing energies
    ekin  = 0.5*ekin
    erotx = 0.5*erotx
    eroty = 0.5*eroty
    erotz = 0.5*erotz
    erot  = sqrt(erotx*erotx + eroty*eroty + erotz*erotz)

    ekin = reduceall_mpi('+', ekin)
    erot = reduceall_mpi('+', erot)

    if (gravity) then
       if (.not. calc_exact_epot) then
          ! Approximate the potential enegy by approximating a uniform density sphere.
          ! If half the mass is in a sphere of radius epot_rad, then the total mass
          ! (assuming isotropy) should be in a sphere with twice this volume,
          ! (i.e. epot_rad -> epot_rad*2^(1/3)).
          epot_rad = epot_rad/epot_mass*1.25992
          epot     = -0.6*(epot_mass/pmassgas1)**2/epot_rad
       endif

       ! CHECK 4: ratio of thermal to gravitational energy alpha <= 1/2 (Eq. 2.9 of BBP95)
       alpha_grav = abs(etherm/epot)
       if (alpha_grav > 0.5) then
          ifail     = inosink_therm
          ifail_array(inosink_therm) = 1
       endif

       ! CHECK 5: ratio of thermal to grav plus ratio of rotational to grav energy <= 1 (Eq. 2.10 of BBP95)
       alphabeta_grav = alpha_grav + abs(erot/epot)
       if (alphabeta_grav > 1.0) then
          ifail     = inosink_grav
          ifail_array(inosink_grav) = 1
       endif

       ! CHECK 7: particle i is at minimum in potential
       if (poteni > potenj_min) then
          ifail = inosink_poten
          ifail_array(inosink_poten) = 1
       endif
    else
       alpha_grav     = 0.0
       alphabeta_grav = 0.0
    endif

    ! CHECK 6: total energy of clump is < 0
    etot = ekin + etherm + epot
    if (etot > 0.) then
       ifail     = inosink_Etot
       ifail_array(inosink_Etot) = 1
    endif
 else
    alpha_grav     = 0.0
    alphabeta_grav = 0.0
    etot           = 0.0
 endif

 ! communicate failure to all MPI threads
 ifail = int(reduceall_mpi('max',ifail))
 ifail_array = int(reduceall_mpi('max',ifail_array))

 ! override failure if the candidate particle is too dense! (some critera still apply)
 if (ForceCreation) then
    if (ifail > 0 .and. is_accretable(itype) .and. hi < 0.5*h_acc) then
       if (id==id_rhomax) then
          write(iprint,"(/,1x,a)")'ptmass_create: OVERRIDING sink failure creation given high density'
          ! list all failure modes that are overridden
          if (ifail_array(inosink_therm)==1) then
             write(iprint,"(/,1x,a,es10.3)") &
             'ptmass_create: FAILURE OVERRIDED when thermal energy/grav energy > 0.5: alpha_grav = ',alpha_grav
          endif
          if (ifail_array(inosink_grav)==1) then
             write(iprint,"(/,1x,a,2es10.3)") &
             'ptmass_create: FAILURE OVERRIDED when alpha_grav + beta_grav > 1, alpha, beta = ',alpha_grav, abs(erot/epot)
          endif
          if (ifail_array(inosink_Etot)==1) then
             write(iprint,"(/,1x,a,es10.3)") &
            'ptmass_create: FAILURE OVERRIDED when total energy > 0, etot = ',etot
          endif
          if (ifail_array(inosink_poten)==1) then
             write(iprint,"(/,1x,a,'phi = ',es10.3,' min =',es10.3)") &
             'ptmass_create: FAILURE OVERRIDED when not at potential minimum ',poteni,potenj_min
          endif
          if (ifail_array(inosink_divv)==1) then
             write(iprint,"(/,1x,a,es10.3)") 'ptmass_create: FAILURE OVERRIDED when  div v > 0', divvi
          endif
       endif
       ifail       = 0
       ifail_array = 0
    endif
 endif

 if (iverbose >= 1 .and. id==id_rhomax) then
    select case(ifail)
    case(0)
       write(iprint,"(1x,a)") 'ptmass_create: OK'
    case(inosink_active)
       write(iprint,"(/,1x,a)") &
       'ptmass_create: FAILED because not all particles within h_acc are active'
    case(inosink_therm)
       write(iprint,"(/,1x,a,es10.3)") &
       'ptmass_create: FAILED because thermal energy/grav energy > 0.5: alpha_grav = ',alpha_grav
    case(inosink_grav)
       write(iprint,"(/,1x,a,2es10.3)") &
       'ptmass_create: FAILED because alpha_grav + beta_grav > 1, alpha, beta = ',alpha_grav, abs(erot/epot)
    case(inosink_Etot)
       write(iprint,"(/,1x,a,es11.3)") &
       'ptmass_create: FAILED because total energy > 0, etot = ',etot
    case(inosink_poten)
       write(iprint,"(/,1x,a,'phi = ',es10.3,' min =',es10.3)") &
       'ptmass_create: FAILED because not at potential minimum ',poteni,potenj_min
    case default
       write(iprint,"(/,1x,a)") 'ptmass_create: FAILED (unknown reason)'
    end select
 endif
 !
 ! create new point mass, at position of original particle but with zero mass. Then accrete particles within hacc to form sink
 !
 if (ifail==0) then
    nptmass = nptmass + 1
    if (nptmass > maxptmass) call fatal('ptmass_create','nptmass > maxptmass')
    n = nptmass
    xyzmh_ptmass(:,n)      = 0.              ! zero all quantities by default
    xyzmh_ptmass(1:3,n)    = (/xi,yi,zi/)
    xyzmh_ptmass(4,n)      = 0.              ! zero mass
    xyzmh_ptmass(ihacc,n)  = h_acc
    xyzmh_ptmass(ihsoft,n) = h_soft_sinkgas
    vxyz_ptmass(:,n)       = 0.              ! zero velocity, get this by accreting
    itypej = igas                            ! default particle type to be accreted
    pmassj = massoftype(igas)                ! default particle mass to be accreted
    !
    ! accrete neighbours (including self)
    !
    nacc       = 0
    dptmass    = 0.
    ibin_wakei = 0 ! dummy argument that has no meaning in this situation
    do n=1,nneigh
       j = listneigh(n)
       if (maxphase==maxp) then
          itypej = iamtype(iphase(j))
          pmassj = massoftype(itypej)
       endif
       fxj = fxyzu(1,j) + fext(1,j)
       fyj = fxyzu(2,j) + fext(2,j)
       fzj = fxyzu(3,j) + fext(3,j)
       call ptmass_accrete(nptmass,nptmass,xyzh(1,j),xyzh(2,j),xyzh(3,j),xyzh(4,j),&
                           vxyzu(1,j),vxyzu(2,j),vxyzu(3,j),fxj,fyj,fzj, &
                           itypej,pmassj,xyzmh_ptmass,vxyz_ptmass,accreted, &
                           dptmass,time,f_acc_local,ibin_wakei,ibin_wakei)

       if (accreted) nacc = nacc + 1
    enddo

    ! perform reduction just for this sink
    dptmass(:,nptmass) = reduceall_mpi('+',dptmass(:,nptmass))
    nacc = int(reduceall_mpi('+', nacc))

    ! update ptmass position, spin, velocity, acceleration, and mass
    fxyz_ptmass(:,nptmass) = 0.0
    fxyz_ptmass_sinksink(:,nptmass) = 0.0
    call update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)

    if (id==id_rhomax) then
       write(iprint,"(a,i3,a,4(es10.3,1x),a,i6,a,es10.3)") ' created ptmass #',nptmass,&
       ' at (x,y,z,t)=(',xyzmh_ptmass(1:3,nptmass),time,') by accreting ',nacc,' particles: M=',xyzmh_ptmass(4,nptmass)
    endif
    if (nacc <= 0) call fatal('ptmass_create',' created ptmass but failed to accrete anything')
    !
    ! open new file to track new sink particle details & and update all sink-tracking files;
    ! fxyz_ptmass, fxyz_ptmass_sinksink are total force on sinks and sink-sink forces.
    !
    if (write_one_ptfile) then
       if (nptmass==1) call pt_open_sinkev(0)  ! otherwise file is already open
    else
       call pt_open_sinkev(nptmass)
    endif
    call pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink)
 else
    !
    ! record failure reason for summary
    !
    call summary_ptmass_fail(ifail)
 endif
 ! print details to file, if requested
 if (record_created) then
    write(iscfile,'(es18.10,1x,3(i18,1x),8(es18.9,1x),8(i18,1x))') &
       time,nptmass+1,itest,nneigh_act,rhoh(hi,pmassi),divvi,alpha_grav,alphabeta_grav,etot,epot,ekin,etherm,ifail_array
    call flush(iscfile)
 endif

end subroutine ptmass_create

!-----------------------------------------------------------------------
!+
!  Merge sinks
!  If sinks are within r_merge_uncond, they will be automatically merged
!  If sinks are within r_merge_cond, they will merge if they are bound
!  A system is bound if
!     Ekin + Epot < 0
!     0.5*mu*dv^2 - G*m1*m2/r < 0
!  where
!     mu = m1*m2/(m1+m2)
!  is the reduced mass.  Therefore, a system is bound if
!     0.5*m1*m2/(m1+m2) dv^2 - G*m1*m2/dr < 0
!  which can be rearranged to
!     0.5*dv^2 - G*(m1+m2)/dr < 0
!  to remove a division.  Therefore, in code units, we use
!     Ekin = 0.5*dv^2
!     Epot = -(m1+m2)/dr
!
!  The merging is similar to that in update_ptmass.
!  We do not remove merged sinks from the list, but tag them with a
!  negative mass.
!+
!-----------------------------------------------------------------------
subroutine merge_sinks(time,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,merge_ij)
 use io,    only:iprint,warning,iverbose,id,master
 real,    intent(in)    :: time
 integer, intent(in)    :: nptmass,merge_ij(nptmass)
 real,    intent(inout) :: xyzmh_ptmass(nsinkproperties,nptmass)
 real,    intent(inout) :: vxyz_ptmass(3,nptmass),fxyz_ptmass(4,nptmass)
 integer :: i,j
 real    :: rr2,xi,yi,zi,mi,vxi,vyi,vzi,xj,yj,zj,mj,vxj,vyj,vzj,Epot,Ekin
 real    :: mij,mij1
 logical :: lmerge
 character(len=15) :: typ

 do i=1,nptmass
    if (merge_ij(i) > 0 .and. xyzmh_ptmass(4,i) > 0.) then
       j = merge_ij(i)
       if (merge_ij(j) == i .and. xyzmh_ptmass(4,j) > 0.) then
          lmerge = .false.
          xi  = xyzmh_ptmass(1,i)
          yi  = xyzmh_ptmass(2,i)
          zi  = xyzmh_ptmass(3,i)
          mi  = xyzmh_ptmass(4,i)
          xj  = xyzmh_ptmass(1,j)
          yj  = xyzmh_ptmass(2,j)
          zj  = xyzmh_ptmass(3,j)
          mj  = xyzmh_ptmass(4,j)
          vxi = vxyz_ptmass(1,i)
          vyi = vxyz_ptmass(2,i)
          vzi = vxyz_ptmass(3,i)
          vxj = vxyz_ptmass(1,j)
          vyj = vxyz_ptmass(2,j)
          vzj = vxyz_ptmass(3,j)
          rr2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
          if (rr2 < r_merge_uncond2) then
             lmerge = .true.
             typ    = 'unconditionally'
          elseif (rr2 < r_merge_cond2) then
             Ekin = 0.5*( (vxi-vxj)**2 + (vyi-vyj)**2 + (vzi-vzj)**2 )
             Epot = -(mi+mj)/rr2
             if (Ekin + Epot < 0.) lmerge = .true.
             typ    = 'conditionally'
          endif
          if (lmerge) then
             ! Add angular momentum of sink particle i using old properties (taken about the origin)
             xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) + mi*(yi*vzi - zi*vyi)
             xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) + mi*(zi*vxi - xi*vzi)
             xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) + mi*(xi*vyi - yi*vxi)
             ! Calculate new masses
             mij  = mi + mj
             mij1 = 1.0/mij
             ! Update quantities
             xyzmh_ptmass(1:3,i)    = (xyzmh_ptmass(1:3,i)*mi + xyzmh_ptmass(1:3,j)*mj)*mij1
             xyzmh_ptmass(4,i)      = mij
             xyzmh_ptmass(imacc,i)  = xyzmh_ptmass(imacc,i)  + xyzmh_ptmass(imacc,j)
             xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) + xyzmh_ptmass(ispinx,j) + mj*(yj*vzj - zj*vyj)
             xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) + xyzmh_ptmass(ispiny,j) + mj*(zj*vxj - xj*vzj)
             xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) + xyzmh_ptmass(ispinz,j) + mj*(xj*vyj - yj*vxj)
             vxyz_ptmass(1:3,i)     = (vxyz_ptmass(1:3,i)*mi + vxyz_ptmass(1:3,j)*mj)*mij1
             fxyz_ptmass(1:3,i)     = (fxyz_ptmass(1:3,i)*mi + fxyz_ptmass(1:3,j)*mj)*mij1
             ! Subtract angular momentum of sink particle using new properties (taken about the origin)
             xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) &
                                    - mij*(xyzmh_ptmass(2,i)*vxyz_ptmass(3,i) - xyzmh_ptmass(3,i)*vxyz_ptmass(2,i))
             xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) &
                                    - mij*(xyzmh_ptmass(3,i)*vxyz_ptmass(1,i) - xyzmh_ptmass(1,i)*vxyz_ptmass(3,i))
             xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) &
                                    - mij*(xyzmh_ptmass(1,i)*vxyz_ptmass(2,i) - xyzmh_ptmass(2,i)*vxyz_ptmass(1,i))
             ! Kill sink j by setting negative mass
             xyzmh_ptmass(4,j)      = -abs(mj)
             ! print success
             write(iprint,"(/,1x,3a,I8,a,I8,a,F10.4)") 'merge_sinks: ',typ,' merged sinks ',i,' & ',j,' at time = ',time
          elseif (id==master .and. iverbose>=1) then
             write(iprint,"(/,1x,a,I8,a,I8,a,F10.4)") &
             'merge_sinks: failed to conditionally merge sinks ',i,' & ',j,' at time = ',time
          endif
       elseif (xyzmh_ptmass(4,j) > 0. .and. id==master .and. iverbose>=1) then
          write(iprint,"(/,1x,a,I8,a,I8,a,F10.4)") &
          'merge_sinks: There is a mismatch in sink indicies and relative proximity for ',i,' & ',j,' at time = ',time
       endif
    endif
 enddo

end subroutine merge_sinks

!-----------------------------------------------------------------------
!+
!  Open files to track sink particle data
!+
!-----------------------------------------------------------------------
subroutine init_ptmass(nptmass,logfile)
 integer,          intent(in) :: nptmass
 character(len=*), intent(in) :: logfile
 integer                      :: i,idot
 character(len=150)           :: filename

 if (id /= master) return ! only do this on master thread
 !
 !--Extract prefix & suffix
 !
 idot = index(logfile,'.')
 if (idot==0) idot = len_trim(logfile) + 1
 pt_prefix = logfile(1:idot-3)

 !
 !--Define file name components and finalise suffix & open files
 !
 if (icreate_sinks > 0) then
    write_one_ptfile = .true.
    write(pt_suffix,'(2a)') logfile(len(trim(pt_prefix))+1:idot-1),".sink"
    if (nptmass > 0) call pt_open_sinkev(0)
 else
    write_one_ptfile = .false.
    write(pt_suffix,'(2a)') logfile(len(trim(pt_prefix))+1:idot-1),".ev"
    do i = 1,nptmass
       call pt_open_sinkev(i)
    enddo
 endif
 !
 !--Open file for tracking sink creation (if required)
 !
 if (record_created) then
    filename = trim(pt_prefix)//"SinkCreated"//trim(pt_suffix)
    open(unit=iscfile,file=trim(filename),form='formatted',status='replace')
    write(iscfile,'("# Data of particles attempting to be converted into sinks.  Columns 13-20: 0 = T, 1 = F")')
    write(iscfile,"('#',20(1x,'[',i2.2,1x,a11,']',2x))") &
           1,'time', &
           2,'nptmass+1', &
           3,'itest',     &
           4,'neigh',     &
           5,'rho',       &
           6,'div v',     &
           7,'alpha',     &
           8,'alphabeta', &
           9,'etot',      &
          10,'epot',      &
          11,'ekin',      &
          12,'etherm',    &
          13,'is gas',    &
          14,'div v < 0', &
          15,'2h < h_acc',&
          16,'all active',&
          17,'alpha < 0', &
          18,'a+b <= 1',  &
          19,'etot < 0',  &
          20,'pot_min'
 else
    iscfile = -abs(iscfile)
 endif

end subroutine init_ptmass
!-----------------------------------------------------------------------
!+
!  finalise ptmass stuff, free memory, close files
!+
!-----------------------------------------------------------------------
subroutine finish_ptmass(nptmass)
 integer, intent(in) :: nptmass

 call pt_close_sinkev(nptmass)

end subroutine finish_ptmass
!-----------------------------------------------------------------------
!+
!  write open sink data files
!+
!-----------------------------------------------------------------------
subroutine pt_open_sinkev(num)
 integer, intent(in) :: num
 integer             :: iunit
 character(len=200)  :: filename

 if (id /= master) return ! only do this on master thread

 if (write_one_ptfile) then
    write(filename,'(2a)') trim(pt_prefix),trim(pt_suffix)
 else
    write(filename,'(2a,I4.4,2a)') trim(pt_prefix),"Sink",num,"N",trim(pt_suffix)
 endif
 iunit = iskfile+num
 open(unit=iunit,file=trim(filename),form='formatted',status='replace')
 if (write_one_ptfile) then
    write(iunit,'(a)') 'To extract one file per sink: make sinks; ./phantomsinks '
 endif
 write(iunit,"('#',20(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',    &
          2,'x',       &
          3,'y',       &
          4,'z',       &
          5,'mass',    &
          6,'vx',      &
          7,'vy',      &
          8,'vz',      &
          9,'spinx',   &
         10,'spiny',   &
         11,'spinz',   &
         12,'macc',    &  ! total mass accreted
         13,'fx',      &
         14,'fy',      &
         15,'fz',      &
         16,'fssx',    &
         17,'fssy',    &
         18,'fssz',    &
         19,'sink ID', &
         20,'nptmass'

end subroutine pt_open_sinkev
!-----------------------------------------------------------------------
!+
!  close sink data files
!+
!-----------------------------------------------------------------------
subroutine pt_close_sinkev(nptmass)
 integer, intent(in) :: nptmass
 integer             :: i,iunit

 if (id == master) then ! only on master thread
    if (write_one_ptfile) then
       close(iskfile)
    else
       do i = 1,nptmass
          iunit = iskfile+i
          close(iunit)
       enddo
    endif
 endif

end subroutine pt_close_sinkev
!-----------------------------------------------------------------------
!+
!  write sink data to files
!+
!-----------------------------------------------------------------------
subroutine pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink)
 integer, intent(in) :: nptmass
 real,    intent(in) :: time, xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),fxyz_ptmass_sinksink(:,:)
 integer             :: i,iunit

 if (id /= master) return ! only do this on master thread

 iunit = iskfile
 do i = 1,nptmass
    if (.not. write_one_ptfile) iunit = iskfile+i
    if (xyzmh_ptmass(4,i) > 0.) then
       write(iunit,"(18(1pe18.9,1x),2(I18,1x))") &
       time, xyzmh_ptmass(1:4,i),vxyz_ptmass(1:3,i), &
       xyzmh_ptmass(ispinx,i),xyzmh_ptmass(ispiny,i),xyzmh_ptmass(ispinz,i), &
       xyzmh_ptmass(imacc,i),fxyz_ptmass(1:3,i),fxyz_ptmass_sinksink(1:3,i),i,nptmass
       if (i==nptmass .or. (.not. write_one_ptfile)) call flush(iunit)
    endif
 enddo

end subroutine pt_write_sinkev
!-----------------------------------------------------------------------
!+
!  compute mass accretion rate
!+
!-----------------------------------------------------------------------
subroutine calculate_mdot(nptmass,time,xyzmh_ptmass)
 use part,        only: imdotav,imacc,i_tlast,i_mlast
 integer, intent(in) :: nptmass
 real,    intent(in) :: time
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 integer             :: i
 real                :: dt

 do i=1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       dt = time - xyzmh_ptmass(i_tlast,i)
       xyzmh_ptmass(imdotav,i) = (xyzmh_ptmass(imacc,i) - xyzmh_ptmass(i_mlast,i))/dt
       xyzmh_ptmass(i_mlast,i) = xyzmh_ptmass(imacc,i)
       xyzmh_ptmass(i_tlast,i) = time
    endif
 enddo
end subroutine calculate_mdot

!-----------------------------------------------------------------------
!+
!  calculate (weighted) sum of particle mass enclosed in sink softening radius
!+
!-----------------------------------------------------------------------
subroutine ptmass_calc_enclosed_mass(nptmass,npart,xyzh)
 use part, only:sink_has_heating,imassenc,ihsoft,massoftype,igas,xyzmh_ptmass,isdead_or_accreted
 use ptmass_heating, only:isink_heating,heating_kernel
 use kernel,         only:radkern2
 integer, intent(in) :: nptmass,npart
 real,    intent(in) :: xyzh(:,:)
 integer             :: i,j
 real                :: wi,q2,x0,y0,z0,hsoft21

 do i = 1,nptmass
    if (.not. sink_has_heating(xyzmh_ptmass(:,i))) cycle
    wi = 0.
    x0 = xyzmh_ptmass(1,i)
    y0 = xyzmh_ptmass(2,i)
    z0 = xyzmh_ptmass(3,i)
    hsoft21 = 1./xyzmh_ptmass(ihsoft,i)**2

    !$omp parallel do default (none) &
    !$omp reduction(+:wi) &
    !$omp shared(npart,xyzh,x0,y0,z0,i,hsoft21,isink_heating) &
    !$omp private(j,q2)
    do j = 1,npart
       if (.not. isdead_or_accreted(xyzh(4,j))) then
          q2 = ((xyzh(1,j)-x0)**2 + (xyzh(2,j)-y0)**2 + (xyzh(3,j)-z0)**2)*hsoft21
          if (q2 < radkern2) wi = wi + heating_kernel(q2,isink_heating)  ! wj = 1 for uniform heating
       endif
    enddo
    !$omp end parallel do
    xyzmh_ptmass(imassenc,i) = wi * massoftype(igas)
 enddo

end subroutine ptmass_calc_enclosed_mass

!-----------------------------------------------------------------------
!+
!  writes sink particle options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_ptmass(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling sink particles'
 if (gravity) then
    call write_inopt(icreate_sinks,'icreate_sinks','allow automatic sink particle creation',iunit)
    if (icreate_sinks > 0) then
       call write_inopt(rho_crit_cgs,'rho_crit_cgs','density above which sink particles are created (g/cm^3)',iunit)
       call write_inopt(r_crit,'r_crit','critical radius for point mass creation (no new sinks < r_crit from existing sink)', &
                        iunit)
       call write_inopt(h_acc, 'h_acc' ,'accretion radius for new sink particles',iunit)
       if (f_crit_override > 0. .or. l_crit_override) then
          call write_inopt(f_crit_override,'f_crit_override' ,'unconditional sink formation if rho > f_crit_override*rho_crit',&
                           iunit)
       endif
       call write_inopt(h_soft_sinkgas,'h_soft_sinkgas','softening length for new sink particles', iunit)
    endif
 endif
 call write_inopt(h_soft_sinksink,'h_soft_sinksink','softening length between sink particles',iunit)
 call write_inopt(f_acc,'f_acc','particles < f_acc*h_acc accreted without checks',iunit)
 call write_inopt(r_merge_uncond,'r_merge_uncond','sinks will unconditionally merge within this separation',iunit)
 call write_inopt(r_merge_cond,'r_merge_cond','sinks will merge if bound within this radius',iunit)

end subroutine write_options_ptmass

!-----------------------------------------------------------------------
!+
!  reads sink particle options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ptmass(name,valstring,imatch,igotall,ierr)
 use io,         only:warning,fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 real                          :: h_soft   ! to ensure backwards compatibility
 character(len=30), parameter  :: label = 'read_options_ptmass'

 ! none of the options apply if no gravity
 if (.not.gravity) then
    igotall = .true.
 endif

 imatch  = .true.
 select case(trim(name))
 case('icreate_sinks')
    read(valstring,*,iostat=ierr) icreate_sinks
    ngot = ngot + 1
    if (icreate_sinks < 0) call fatal(label,'sink creation option out of range')
 case('rho_crit_cgs')
    read(valstring,*,iostat=ierr) rho_crit_cgs
    if (rho_crit_cgs < 0.) call fatal(label,'rho_crit < 0')
    ngot = ngot + 1
 case('r_crit')
    read(valstring,*,iostat=ierr) r_crit
    if (r_crit < 0.) call fatal(label,'r_crit < 0')
    if (icreate_sinks==1 .and. r_crit < 2.0*h_acc) then
       call warning(label,'Strongly suggest r_crit >= 2.0*h_acc')
    endif
    ngot = ngot + 1
 case('h_acc')
    read(valstring,*,iostat=ierr) h_acc
    if (h_acc <= 0.) call fatal(label,'h_acc < 0')
    ngot = ngot + 1
 case('f_crit_override')
    read(valstring,*,iostat=ierr) f_crit_override
    if (f_crit_override < 0.) f_crit_override = 0.  ! reset to zero since a negative value does not make sense
    if (f_crit_override > 0. .and. f_crit_override < 100. ) call fatal(label,'Give star formation a chance! Reset to > 100')
    l_crit_override = .true.
 case('h_soft')  ! to ensure backwards compatibility
    read(valstring,*,iostat=ierr) h_soft
    if (h_soft > 0.) call fatal(label,'h_soft has been renamed to h_soft_sinkgas.  Please modify in-file before retrying')
 case('h_soft_sinkgas')
    read(valstring,*,iostat=ierr) h_soft_sinkgas
    if (h_soft_sinkgas < 0.) call fatal(label,'h_soft_sinkgas < 0')
    ngot = ngot + 1
 case('h_soft_sinksink')
    read(valstring,*,iostat=ierr) h_soft_sinksink
    if (h_soft_sinksink < 0.) call fatal(label,'h_soft_sinksink < 0')
    ngot = ngot + 1
 case('f_acc')
    read(valstring,*,iostat=ierr) f_acc
    if (f_acc < 0.0) call fatal(label,'f_acc < 0')
    if (f_acc > 1.0) call fatal(label,'f_acc > 1')
    ngot = ngot + 1
 case('r_merge_uncond')
    read(valstring,*,iostat=ierr) r_merge_uncond
    if (icreate_sinks==1 .and. r_merge_uncond < 2.0*h_acc) then
       call warning(label,'Strongly suggest r_merge_uncond >= 2.0*h_acc')
    endif
    ngot = ngot + 1
 case('r_merge_cond')
    read(valstring,*,iostat=ierr) r_merge_cond
    if (r_merge_cond > 0. .and. r_merge_cond < r_merge_uncond) call fatal(label,'0 < r_merge_cond < r_merge_uncond')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
 if (icreate_sinks > 0) then
    igotall = (ngot >= 8)
 else
    igotall = (ngot >= 4)
 endif

end subroutine read_options_ptmass
!-----------------------------------------------------------------------
end module ptmass
