!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: ptmass
!
!  DESCRIPTION:
!  This module contains everything to do with
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
!  REFERENCES: Bate, Bonnell & Price (1995), MNRAS 277, 362-376 [BBP95]
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    f_acc           -- particles < f_acc*h_acc accreted without checks
!    h_acc           -- accretion radius for new sink particles
!    h_soft_sinkgas  -- softening length for new sink particles
!    h_soft_sinksink -- softening length between sink particles
!    icreate_sinks   -- allow automatic sink particle creation
!    r_crit          -- critical radius for point mass creation (no new sinks < r_crit from existing sink)
!    rho_crit_cgs    -- density above which sink particles are created (g/cm^3)
!
!  DEPENDENCIES: boundary, dim, eos, externalforces, fastmath,
!    infile_utils, io, io_summary, kdtree, kernel, linklist, mpiutils,
!    options, part, units
!+
!--------------------------------------------------------------------------
module ptmass
 use dim,  only:maxptmass
 use part, only:nsinkproperties,gravity,is_accretable
 use io,   only:iscfile,ipafile,iskfile
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"
 public :: init_ptmass, finish_ptmass
 public :: pt_write_sinkev, pt_close_sinkev
 public :: get_accel_sink_gas, get_accel_sink_sink
 public :: ptmass_predictor, ptmass_corrector
 public :: ptmass_not_obscured
 public :: ptmass_accrete, ptmass_create
 public :: write_options_ptmass, read_options_ptmass
 public :: update_ptmass

 ! settings affecting routines in module (read from/written to input file)
 integer, public :: icreate_sinks = 0
 real,    public :: rho_crit_cgs  = 1.e-10
 real,    public :: r_crit = 5.e-3
 real,    public :: h_acc  = 1.e-3
 real,    public :: f_acc  = 0.8
 real,    public :: h_soft_sinkgas  = 0.0
 real,    public :: h_soft_sinksink = 0.0

 ! additional public variables
 integer, public :: ipart_rhomax
 real,    public :: r_crit2,rho_crit

 integer,         public :: rhomax_ipart
 integer(kind=1), public :: rhomax_iphase,rhomax_ibin
 real,            public :: rhomax_xyzh(4),rhomax_vxyz(3)
 real(kind=4),    public :: rhomax_divv

 ! calibration of timestep control on sink-sink and sink-gas orbital integration
 ! this is hardwired because can be adjusted by changing C_force
 ! just means that with the default setting of C_force the orbits are accurate
 real, parameter :: dtfacphi  = 0.05
 real, parameter :: dtfacphi2 = dtfacphi*dtfacphi

 ! parameters to control output regarding sink particles
 logical, private, parameter :: record_created  = .false.  ! verbose tracking of why sinks are not created
 logical, private, parameter :: record_accreted = .false.  ! verbose tracking of particle accretion
 character(len=50), private  :: pt_prefix = 'Sink'
 character(len=50), private  :: pt_suffix = '00.ev'

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
!  if (tofrom==.true.)  acceleration from/to gas particles due to sink particles
!  if (tofrom==.false.) acceleration on gas due to sink particles (but not vice-versa)
!+
!----------------------------------------------------------------
subroutine get_accel_sink_gas(nptmass,xi,yi,zi,hi,xyzmh_ptmass,fxi,fyi,fzi,phi, &
                                       pmassi,fxyz_ptmass,fonrmax,dtphi2)
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use kernel,   only:kernel_softening,radkern
 use part,     only:ihacc,ihsoft
 integer,           intent(in)    :: nptmass
 real,              intent(in)    :: xi,yi,zi,hi
 real,              intent(inout) :: fxi,fyi,fzi,phi
 real,              intent(in)    :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real,    optional, intent(in)    :: pmassi
 real,    optional, intent(inout) :: fxyz_ptmass(4,maxptmass)
 real,    optional, intent(out)   :: fonrmax,dtphi2
 real                             :: ftmpxi,ftmpyi,ftmpzi
 real                             :: dx,dy,dz,rr2,ddr,dr3,f1,f2,pmassj
 real                             :: hsoft,hsoft1,hsoft21,q2i,qi,psoft,fsoft
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
    if (hsoft > 0.0)  hsoft = max(hsoft,hi)

    rr2    = dx*dx + dy*dy + dz*dz + epsilon(rr2)
#ifdef FINVSQRT
    ddr    = finvsqrt(rr2)
#else
    ddr    = 1./sqrt(rr2)
#endif
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
    endif

    if (tofrom) then
       ! backreaction of gas onto sink
       fxyz_ptmass(1,j) = fxyz_ptmass(1,j) + dx*f2
       fxyz_ptmass(2,j) = fxyz_ptmass(2,j) + dy*f2
       fxyz_ptmass(3,j) = fxyz_ptmass(3,j) + dz*f2
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

 return
end subroutine get_accel_sink_gas

!----------------------------------------------------------------
!+
!  Compute force on sink particles due to other sinks and
!  from external potentials
!+
!----------------------------------------------------------------
subroutine get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,phitot,dtsinksink,&
            iexternalforce,ti)
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use externalforces, only:externalforce
 use kernel,   only:kernel_softening,radkern
 !$ use omputils, only:ipart_omp_lock
 integer, intent(in)  :: nptmass
 real,    intent(in)  :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real,    intent(out) :: fxyz_ptmass(4,maxptmass)
 real,    intent(out) :: phitot,dtsinksink
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: ti
 real    :: xi,yi,zi,pmassi,pmassj,fxi,fyi,fzi,phii
 real    :: ddr,dx,dy,dz,rr2,dr3,f1,f2,term
 real    :: hsoft,hsoft1,hsoft21,q2i,qi,psoft,fsoft
 real    :: fextx,fexty,fextz,phiext !,hsofti
 real    :: fterm, pterm
 integer :: i,j

 dtsinksink = huge(dtsinksink)
 fxyz_ptmass(:,:) = 0.
 phitot = 0.
 !
 !--compute N^2 forces on point mass particles due to each other
 !
 !$omp parallel do default(none) &
 !$omp shared(nptmass,xyzmh_ptmass,fxyz_ptmass,ipart_omp_lock) &
 !$omp shared(iexternalforce,ti,h_soft_sinksink) &
 !$omp private(i,xi,yi,zi,pmassi,pmassj) &
 !$omp private(dx,dy,dz,rr2,ddr,dr3,f1,f2) &
 !$omp private(fxi,fyi,fzi,phii,term) &
 !$omp private(fextx,fexty,fextz,phiext) &
 !$omp private(hsoft,hsoft1,hsoft21,q2i,qi,psoft,fsoft) &
 !$omp private(fterm,pterm) &
 !$omp reduction(min:dtsinksink) &
 !$omp reduction(+:phitot)
 do i=1,nptmass
    xi     = xyzmh_ptmass(1,i)
    yi     = xyzmh_ptmass(2,i)
    zi     = xyzmh_ptmass(3,i)
    pmassi = xyzmh_ptmass(4,i)
    !hsofti = xyzmh_ptmass(5,i)
    fxi    = 0.
    fyi    = 0.
    fzi    = 0.
    phii   = 0.
    do j=i+1,nptmass
       dx     = xi - xyzmh_ptmass(1,j)
       dy     = yi - xyzmh_ptmass(2,j)
       dz     = zi - xyzmh_ptmass(3,j)
       pmassj = xyzmh_ptmass(4,j)
       !hsoftj = xyzmh_ptmass(5,j)

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
          hsoft1 = 1.0/h_soft_sinksink
          hsoft21= hsoft1**2
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

          ! acceleration of sink2 from sink1
          f2    = pmassi*fterm
          term  = pmassi*pmassj*psoft*hsoft1
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

          ! acceleration of sink2 from sink1
          f2   = pmassi*dr3
          term = -pmassi*pmassj*ddr
       endif

       phitot = phitot + term  ! potential (G M_1 M_2/r)

       !$ call omp_set_lock(ipart_omp_lock(j))
       fxyz_ptmass(1,j) = fxyz_ptmass(1,j) + dx*f2
       fxyz_ptmass(2,j) = fxyz_ptmass(2,j) + dy*f2
       fxyz_ptmass(3,j) = fxyz_ptmass(3,j) + dz*f2
       fxyz_ptmass(4,j) = fxyz_ptmass(4,j) + pmassi*pterm
       !$ call omp_unset_lock(ipart_omp_lock(j))

    enddo

    !
    !--apply external forces
    !
    if (iexternalforce > 0) then
       call externalforce(iexternalforce,xi,yi,zi,0.,ti,fextx,fexty,fextz,phiext)
       fxi = fxi + fextx
       fyi = fyi + fexty
       fzi = fzi + fextz
       phii   = phii + phiext
       phitot = phitot + phiext
    endif

    !
    !--store sink-sink forces (only)
    !
    !$ call omp_set_lock(ipart_omp_lock(i))
    fxyz_ptmass(1,i) = fxyz_ptmass(1,i) + fxi
    fxyz_ptmass(2,i) = fxyz_ptmass(2,i) + fyi
    fxyz_ptmass(3,i) = fxyz_ptmass(3,i) + fzi
    fxyz_ptmass(4,i) = fxyz_ptmass(4,i) + phii ! Note: No self contribution to the potential for sink-sink softening.
    !$ call omp_unset_lock(ipart_omp_lock(i))

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
    if (f2  >  0) then
       dtsinksink = min(dtsinksink,dtfacphi*sqrt(abs(phii)/f2))
    endif
 enddo

end subroutine get_accel_sink_sink

!----------------------------------------------------------------
!+
!  predictor step for the point masses
!  (called from inside a parallel section)
!+
!----------------------------------------------------------------
subroutine ptmass_predictor(nptmass,dt,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: dt
 real,    intent(inout) :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real,    intent(inout) :: vxyz_ptmass(3,maxptmass)
 real,    intent(in)    :: fxyz_ptmass(4,maxptmass)
 real    :: vxhalfi,vyhalfi,vzhalfi
 integer :: i

 !$omp parallel do schedule(static) default(none) &
 !$omp shared(nptmass,dt,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass) &
 !$omp private(i,vxhalfi,vyhalfi,vzhalfi)
 do i=1,nptmass
    vxhalfi = vxyz_ptmass(1,i) + 0.5*dt*fxyz_ptmass(1,i)
    vyhalfi = vxyz_ptmass(2,i) + 0.5*dt*fxyz_ptmass(2,i)
    vzhalfi = vxyz_ptmass(3,i) + 0.5*dt*fxyz_ptmass(3,i)
    xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dt*vxhalfi
    xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dt*vyhalfi
    xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dt*vzhalfi
    vxyz_ptmass(1,i) = vxhalfi
    vxyz_ptmass(2,i) = vyhalfi
    vxyz_ptmass(3,i) = vzhalfi
 enddo
 !$omp end parallel do

end subroutine ptmass_predictor

!----------------------------------------------------------------
!+
!  corrector step for the point masses
!  (called from inside a parallel section)
!+
!----------------------------------------------------------------
subroutine ptmass_corrector(nptmass,dt,vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,iexternalforce)
 use externalforces, only:update_vdependent_extforce_leapfrog,is_velocity_dependent
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: dt
 real,    intent(inout) :: vxyz_ptmass(3,maxptmass)
 real,    intent(in)    :: fxyz_ptmass(4,maxptmass), xyzmh_ptmass(nsinkproperties,maxptmass)
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
    !$omp shared(vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,dt,nptmass,iexternalforce) &
    !$omp private(vxhalfi,vyhalfi,vzhalfi,fxi,fyi,fzi,fextv) &
    !$omp private(i)
    do i=1,nptmass
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
    enddo
    !$omp end parallel do
 else
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(vxyz_ptmass,fxyz_ptmass,dt,nptmass) &
    !$omp private(i)
    do i=1,nptmass
       vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + 0.5*dt*fxyz_ptmass(1,i)
       vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + 0.5*dt*fxyz_ptmass(2,i)
       vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + 0.5*dt*fxyz_ptmass(3,i)
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
 real,              intent(in)    :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real,              intent(in)    :: vxyz_ptmass(3,maxptmass)
 logical,           intent(out)   :: accreted
 real,              intent(inout) :: dptmass(ndptmass,maxptmass)
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
       !$ call omp_set_lock(ipart_omp_lock(i))

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

       !$ call omp_unset_lock(ipart_omp_lock(i))
       hi = -abs(hi)

       if (record_accreted) then
          !$omp critical(trackacc)
          call fatal('ptmass', 'track_accreted has been deprecated because it relied on OpenMP-unsafe code')
          call track_accreted(time,i,dx,dy,dz,dvx,dvy,dvz,xyzmh_ptmass(1,i),xyzmh_ptmass(2,i), &
            xyzmh_ptmass(3,i),xyzmh_ptmass(4,i),pmassi,vxyz_ptmass(1:3,i),xyzmh_ptmass(1,i), &
            xyzmh_ptmass(2,i),xyzmh_ptmass(3,i),vxyz_ptmass(1,i),vxyz_ptmass(2,i),vxyz_ptmass(3,i))
          !$omp end critical(trackacc)
       endif

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
 use part,  only:ispinx,ispiny,ispinz,imacc
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
subroutine ptmass_create(nptmass,npart,itest,xyzh,vxyzu,fxyzu,fext,divcurlv,massoftype,&
                         xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,time,&
                         xyzhi,vxyzi,iphasei_test,divvi_test,ibini)
 use part,   only:ihacc,ihsoft,igas,iamtype,get_partinfo,iphase,iactive,maxphase,rhoh, &
                  ispinx,ispiny,ispinz,fxyz_ptmass_sinksink
 use dim,    only:maxp,maxneigh,maxvxyzu
 use kdtree, only:getneigh
 use kernel, only:kernel_softening
 use io,     only:iprint,fatal,iverbose,nprocs
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
#ifdef IND_TIMESTEPS
 use part,     only:ibin,ibin_wake
#endif
 use linklist, only:getneigh_pos,ifirstincell
 use eos,      only:equationofstate,gamma,gamma_pwp,utherm
 use options,  only:ieos
 use units,    only:unit_density
 use io_summary, only:summary_variable_rhomax,summary_ptmass_fail
 use mpiutils, only:reduceall_mpi
 integer,         intent(inout) :: nptmass
 integer,         intent(in)    :: npart,itest
 real,            intent(inout) :: xyzh(:,:)
 real,            intent(in)    :: vxyzu(:,:), fxyzu(:,:), fext(:,:), massoftype(:)
 real(4),         intent(in)    :: divcurlv(:,:)
 real,            intent(inout) :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real,            intent(inout) :: vxyz_ptmass(3,maxptmass),fxyz_ptmass(4,maxptmass)
 real,            intent(in)    :: time
 real, optional,  intent(in)    :: xyzhi(4),vxyzi(3)
 real(kind=4),    optional, intent(in) :: divvi_test
 integer(kind=1), optional, intent(in) :: iphasei_test,ibini
 integer(kind=1)    :: iphasei,ibin_wakei
 integer            :: nneigh
 integer            :: listneigh(maxneigh)
 integer, parameter :: maxcache      = 12000
 integer, parameter :: nneigh_thresh = 1024 ! approximate epot if neigh>neigh_thresh; (-ve for off)
#ifdef IND_TIMESTEPS
 integer(kind=1)    :: ibin_itest
#endif
 real    :: xyzcache(maxcache,3)
 real    :: dptmass(ndptmass,nptmass+1)
 real    :: xi,yi,zi,hi,hi1,hi21,xj,yj,zj,hj1,hj21,xk,yk,zk,hk1
 real    :: rij2,rik2,rjk2,dx,dy,dz,h_acc2
 real    :: vxi,vyi,vzi,dv2,dvx,dvy,dvz
 real    :: alpha_grav,alphabeta_grav,radxy2,radxz2,radyz2
 real    :: etot,epot,ekin,etherm,erot,erotx,eroty,erotz
 real    :: rcrossvx,rcrossvy,rcrossvz,fxj,fyj,fzj
 real    :: pmassi,pmassj,pmassk,ponrhoj,rhoj,spsoundj
 real    :: q2i,qi,psofti,psoftj,psoftk,fsoft,epot_mass,epot_rad,pmassgas1
 real(4) :: divvi
 integer :: ifail,nacc,j,k,n,nk,itype,itypej,itypek,ifail7(7)
 logical :: accreted,iactivej,isdustj,iactivek,isdustk,calc_exact_epot

 ifail  = 0
 ifail7 = 0

 ! itest's characteristics
 if (present(xyzhi)) then
    xi = xyzhi(1)
    yi = xyzhi(2)
    zi = xyzhi(3)
    hi = xyzhi(4)
    vxi = vxyzi(1)
    vyi = vxyzi(2)
    vzi = vxyzi(3)
    iphasei= iphasei_test
    divvi = divvi_test
#ifdef IND_TIMESTEPS
    ibin_itest = ibini
#endif
 else
    if (itest <= 0 .or. itest > npart) return
    xi = xyzh(1,itest)
    yi = xyzh(2,itest)
    zi = xyzh(3,itest)
    hi = xyzh(4,itest)
    vxi = vxyzu(1,itest)
    vyi = vxyzu(2,itest)
    vzi = vxyzu(3,itest)
    iphasei = iphase(itest)
    divvi = divcurlv(1,itest)
#ifdef IND_TIMESTEPS
    ibin_itest = ibin(itest)
#endif
 endif
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

 call summary_variable_rhomax(itest,rhoh(hi,pmassi)*real(unit_density),iprint,nptmass)

 if (iverbose >= 1) write(iprint,"(a,i10,a)",advance='no') ' ptmass_create: Testing ',itest,' for ptmass creation...'

 ! CHECK 0: make sure particle is a gas particle (sanity check, should be unnecessary)
 if (itype /= igas) then
    if (iverbose >= 1) write(iprint,"(/,1x,a)") 'ptmass_create: FAILED because not a gas particle'
    call summary_ptmass_fail(5)
    if (.not. record_created) return
    ifail7(5) = 1
 endif

 ! CHECK 1: divv < 0
 if (divvi > 0._4) then
    if (iverbose >= 1) write(iprint,"(/,1x,a)") 'ptmass_create: FAILED because div v > 0'
    call summary_ptmass_fail(6)
    if (.not. record_created) return
    ifail7(6) = 1
 endif

 ! CHECK 2: 2h < h_acc
 if (hi > 0.5*h_acc) then
    if (iverbose >= 1) write(iprint,"(/,1x,a,es10.3,a,es10.3,a)") 'ptmass_create: FAILED because 2h > h_acc (',h_acc,')'
    call summary_ptmass_fail(7)
    if (.not. record_created) return
    ifail7(7) = 1
 endif

 h_acc2 = h_acc*h_acc
 ekin   = 0.
 epot   = -epsilon(epot)
 etherm = 0.
 erot   = 0.
 erotx  = 0.
 eroty  = 0.
 erotz  = 0.
 epot_mass = 0.
 epot_rad  = 0.

 ! CHECK 3: all neighbours within h_acc are all active ( & perform math for checks 4-6)
 ! find neighbours within h_acc
 call getneigh_pos((/xi,yi,zi/),0.,h_acc,3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
 ! determine if we should approximate epot
 calc_exact_epot = .true.
 if ((nneigh_thresh > 0 .and. nneigh > nneigh_thresh) .or. (nprocs > 1)) calc_exact_epot = .false.
!$omp parallel default(none) &
!$omp shared(nprocs) &
!$omp shared(maxp,maxphase) &
!$omp shared(nneigh,listneigh,h_acc2,xyzh,xyzcache,vxyzu,massoftype,iphase,pmassgas1,calc_exact_epot) &
!$omp shared(itest,ifail,xi,yi,zi,hi,vxi,vyi,vzi,hi1,hi21,itype,pmassi,ieos,gamma) &
#ifdef PERIODIC
!$omp shared(dxbound,dybound,dzbound) &
#endif
#ifdef IND_TIMESTEPS
!$omp shared(ibin_wake,ibin_itest) &
#endif
!$omp private(n,j,xj,yj,zj,hj1,hj21,psoftj,rij2,nk,k,xk,yk,zk,hk1,psoftk,rjk2,psofti,rik2) &
!$omp private(dx,dy,dz,dvx,dvy,dvz,dv2,isdustj,iactivek,isdustk) &
!$omp private(rhoj,ponrhoj,spsoundj,q2i,qi,fsoft,rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2) &
!$omp firstprivate(pmassj,pmassk,itypej,iactivej,itypek) &
!$omp reduction(+:ekin,erotx,eroty,erotz,etherm,epot,epot_mass,epot_rad)
!$omp do
 over_neigh: do n=1,nneigh
    j = listneigh(n)
    !
    ! get mass and particle type to immediately determine if active and accretable
    if (maxphase==maxp) then
       call get_partinfo(iphase(j),iactivej,isdustj,itypej)
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
    if (rij2 < h_acc2) then

#ifdef IND_TIMESTEPS
       ibin_wake(j) = max(ibin_wake(j),ibin_itest)
       if (.not.iactivej .or. ifail==1) then
          ifail = 1
          cycle over_neigh
       endif
#endif

       dvx = vxi - vxyzu(1,j)
       dvy = vyi - vxyzu(2,j)
       dvz = vzi - vxyzu(3,j)

       hj1  = 1.0/xyzh(4,j)
       hj21 = hj1**2

       ! kinetic energy
       dv2 = dvx*dvx + dvy*dvy + dvz*dvz
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
             etherm = etherm + pmassj*utherm(vxyzu(4,j),rhoj)
          else
             call equationofstate(ieos,ponrhoj,spsoundj,rhoj,xj,yj,zj)
             if (ieos==2 .and. gamma > 1.001) then
                etherm = etherm + pmassj*ponrhoj/(gamma - 1.)
             elseif (ieos==9) then
                etherm = etherm + pmassj*ponrhoj/(gamma_pwp(rhoj) - 1.)
             else
                etherm = etherm + pmassj*1.5*ponrhoj
             endif
          endif
       endif

       ! gravitational potential energy of clump
       if (gravity) then
          if (calc_exact_epot) then
             if (nprocs > 1) call fatal('ptmass_create', 'cannot use calc_exact_epot with MPI')
             ! Calculate potential energy exactly
             !
             ! add contribution of i-j (since, e.g., rij2 is already calculated)
             q2i    = rij2*hi21
             qi     = sqrt(q2i)
             call kernel_softening(q2i,qi,psofti,fsoft)
             q2i    = rij2*hj21
             qi     = sqrt(q2i)
             call kernel_softening(q2i,qi,psoftj,fsoft)
             epot   = epot + 0.5*pmassi*pmassj*(psofti*hi1 + psoftj*hj1)
             !
             ! add contribution of k-j for all k >= j (to avoid double counting, but include self-contribution)
             over_neigh_k: do nk=n,nneigh
                k = listneigh(nk)
                if (k==itest) cycle over_neigh_k ! contribution already added
                if (maxphase==maxp) then
                   call get_partinfo(iphase(k),iactivek,isdustk,itypek)
                   pmassk = massoftype(itypek)
                   if (.not. is_accretable(itypek) ) cycle over_neigh_k
                endif
                !
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
                if (rik2 < h_acc2) then
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
 !  Note that if ifail7(5,6,7)==1 and record_created==.false., this subroutine will
 !  already have been exited, and this loop will never be reached
 if ( record_created ) then
    if ( ifail==1 ) then
       ifail7(1) = 1
    elseif (ifail7(5)==1) then
       ifail = 5
    elseif (ifail7(6)==1) then
       ifail = 6
    elseif (ifail7(7)==1) then
       ifail = 7
    endif
 endif

 !
 !--Continue checks (non-sensical for ifail==1 since energies not completely calculated)
 if (ifail==0 .or. (record_created .and. ifail > 1)) then
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
          ifail     = 2
          ifail7(2) = 1
       endif

       ! CHECK 5: ratio of thermal to grav plus ratio of rotational to grav energy <= 1 (Eq. 2.10 of BBP95)
       alphabeta_grav = alpha_grav + abs(erot/epot)
       if (alphabeta_grav > 1.0) then
          ifail     = 3
          ifail7(3) = 1
       endif
    else
       alpha_grav     = 0.0
       alphabeta_grav = 0.0
    endif

    ! CHECK 6: total energy of clump is < 0
    etot = ekin + etherm + epot
    if (etot > 0.) then
       ifail     = 4
       ifail7(4) = 1
    endif
 else
    alpha_grav     = 0.0
    alphabeta_grav = 0.0
    etot           = 0.0
 endif

 if (iverbose >= 1) then
    select case(ifail)
    case(4)
       write(iprint,"(/,a,es10.3)") &
       'ptmass_create: FAILED because total energy > 0, etot = ',etot
    case(3)
       write(iprint,"(/,a,2es10.3)") &
       'ptmass_create: FAILED because alpha_grav + beta_grav > 1, alpha, beta = ',alpha_grav, abs(erot/epot)
    case(2)
       write(iprint,"(/,a,es10.3)") &
       'ptmass_create: FAILED because thermal energy/grav energy > 0.5: alpha_grav = ',alpha_grav
    case(1)
       write(iprint,"(/,a)") &
       'ptmass_create: FAILED because not all particles within h_acc are active'
    case(0)
       write(iprint,"(a)") 'ptmass_create: OK'
    case default
       write(iprint,"(/,a)") 'ptmass_create: FAILED (unknown reason)'
    end select
 endif

 !
 ! create new point mass, at position of original
 ! particle but with zero mass. Then accrete particles
 ! within hacc to form sink
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
                           dptmass,time,1.0,ibin_wakei,ibin_wakei)

       if (accreted) nacc = nacc + 1
    enddo

    ! perform reduction just for this sink
    dptmass(:,nptmass) = reduceall_mpi('+',dptmass(:,nptmass))
    nacc = int(reduceall_mpi('+', nacc))

    ! update ptmass position, spin, velocity, acceleration, and mass
    call update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)

    write(iprint,"(a,i3,a,4(es10.3,1x),a,i6,a,es10.3)") ' created ptmass #',nptmass,&
     ' at (x,y,z,t)=(',xyzmh_ptmass(1:3,nptmass),time,') by accreting ',nacc,' particles: M=',xyzmh_ptmass(4,nptmass)
    if (nacc <= 0) call fatal('ptmass_create',' created ptmass but failed to accrete anything')
    !
    ! open new file to track new sink particle details & and update all sink-tracking files;
    ! fxyz_ptmass, fxyz_ptmass_sinksink are total force on sinks and sink-sink forces.
    !
    fxyz_ptmass = 0.0
    fxyz_ptmass_sinksink = 0.0
    call pt_open_sinkev(nptmass)
    call pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink)
 else
    !
    ! record failure reason for summary
    !
    call summary_ptmass_fail(ifail)
 endif
 ! print details to file, if requested
 if (record_created) then
    write(iscfile,'(es18.10,1x,3(i18,1x),5(es18.9,1x),7(i18,1x))') &
       time,nptmass+1,itest,nneigh,rhoh(hi,pmassi),divvi,alpha_grav,alphabeta_grav,etot,ifail7
    call flush(iscfile)
 endif

end subroutine ptmass_create

!-----------------------------------------------------------------------
!+
!  Open files to track sink particle data
!+
!-----------------------------------------------------------------------
subroutine init_ptmass(nptmass,logfile,dumpfile)
 integer,          intent(in) :: nptmass
 character(len=*), intent(in) :: logfile,dumpfile
 integer                      :: i,idot,idash
 character(len=150)           :: filename
 !
 !--Extract prefix
 !
 idash = index(dumpfile,'_')
 write(pt_prefix,"(a)") dumpfile(1:idash-1)
 !
 !--Extract suffix
 !
 idot = index(logfile,'.')
 if (idot==0) idot = len_trim(logfile) + 1
 write(pt_suffix,'(2a)') logfile(len(trim(pt_prefix))+1:idot-1),".ev"
 !
 !--Open file for each sink particle
 !
 do i = 1,nptmass
    call pt_open_sinkev(i)
 enddo
 !
 !--Open file for tracking sink creation (if required)
 !
 if (record_created) then
    filename = trim(pt_prefix)//"SinkCreated"//trim(pt_suffix)
    open(unit=iscfile,file=trim(filename),form='formatted',status='replace')
    write(iscfile,'("# Data of particles attempting to be converted into sinks.  Columns 10-17: 0 = T, 1 = F")')
    write(iscfile,"('#',16(1x,'[',i2.2,1x,a11,']',2x))") &
           1,'time', &
           2,'nptmass+1', &
           3,'itest',     &
           4,'neigh',     &
           5,'rho',       &
           6,'div v',     &
           7,'alpha',     &
           8,'alphabeta', &
           9,'etot',      &
          10,'all active',&
          11,'alpha < 0', &
          12,'a+b <= 1',  &
          13,'etot < 0',  &
          14,'is gas',    &
          15,'div v < 0', &
          16,'2h < h_acc'
 else
    iscfile = -abs(iscfile)
 endif
 !
 !--Open file for tracking particle accretion (if required)
 !
 if (record_accreted) then
    filename = trim(pt_prefix)//"ParticleAccretion"//trim(pt_suffix)
    open(unit=ipafile,file=trim(filename),form='formatted',status='replace')
    write(ipafile,'("# Data of accreted particles")')
    write(ipafile,"('#',26(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',  &
          2,'sinki', &
          3,'angx',  &
          4,'angy',  &
          5,'angz',  &
          6,'mpart', &
          7,'msink', &
          8,'mu',    &
          9,'dx',    &
          10,'dy',   &
          11,'dz',   &
          12,'dvx',  &
          13,'dvy',  &
          14,'dvz',  &
          15,'xsink_new', &
          16,'ysink_new', &
          17,'zsink_new', &
          18,'xsink_old', &
          19,'ysink_old', &
          20,'zsink_old', &
          21,'vxsink_new',&
          22,'vysink_new',&
          23,'vzsink_new',&
          24,'vxsink_old',&
          25,'vysink_old',&
          26,'vzsink_old'
 else
    ipafile = -abs(ipafile)
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

 write(filename,'(2a,I4.4,2a)') trim(pt_prefix),"Sink",num,"N",trim(pt_suffix)
 iunit = iskfile+num
 open(unit=iunit,file=trim(filename),form='formatted',status='replace')
 write(iunit,"('#',18(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',  &
          2,'x',     &
          3,'y',     &
          4,'z',     &
          5,'mass',  &
          6,'vx',    &
          7,'vy',    &
          8,'vz',    &
          9,'spinx', &
         10,'spiny', &
         11,'spinz', &
         12,'macc',  &
         13,'fx',    &
         14,'fy',    &
         15,'fz',    &
         16,'fssx',  &
         17,'fssy',  &
         18,'fssz'

end subroutine pt_open_sinkev

!-----------------------------------------------------------------------
!+
!  close sink data files
!+
!-----------------------------------------------------------------------
subroutine pt_close_sinkev(nptmass)
 integer, intent(in) :: nptmass
 integer             :: i,iunit

 do i = 1,nptmass
    iunit = iskfile+i
    close(iunit)
 enddo

end subroutine pt_close_sinkev

!-----------------------------------------------------------------------
!+
!  write sink data to files
!+
!-----------------------------------------------------------------------
subroutine pt_write_sinkev(nptmass,time,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink)
 use part,        only: ispinx,ispiny,ispinz,imacc
 integer, intent(in) :: nptmass
 real,    intent(in) :: time, xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),fxyz_ptmass_sinksink(:,:)
 integer             :: i,iunit

 do i = 1,nptmass
    iunit = iskfile+i
    write(iunit,"(18(1pe18.9,1x))") &
    time, xyzmh_ptmass(1:4,i),vxyz_ptmass(1:3,i), &
    xyzmh_ptmass(ispinx,i),xyzmh_ptmass(ispiny,i),xyzmh_ptmass(ispinz,i), &
    xyzmh_ptmass(imacc,i),fxyz_ptmass(1:3,i),fxyz_ptmass_sinksink(1:3,i)
    call flush(iunit)
 enddo

end subroutine pt_write_sinkev

!-----------------------------------------------------------------------
!+
! writes accreted particle properties to file
! Author: CJN.  Modified: JHW (July 2015)
!+
!-----------------------------------------------------------------------
subroutine track_accreted(time,sinki,dx,dy,dz,dvx,dvy,dvz,xsink,ysink,zsink,msink,mpart,v_new,xs,ys,zs,vxs,vys,vzs)
 integer, intent(in) :: sinki
 real,    intent(in) :: time,dx,dy,dz,dvx,dvy,dvz,xsink,ysink,zsink,msink,mpart
 real,    intent(in) :: v_new(3),xs,ys,zs, vxs,vys,vzs
 real                :: spinm
 real                :: pos(3),vel(3),ang(3)

 pos(1) = dx
 pos(2) = dy
 pos(3) = dz

 vel(1) = dvx
 vel(2) = dvy
 vel(3) = dvz

 spinm  = msink*mpart/(msink+mpart)
 ang(1) = (pos(2)*vel(3) - pos(3)*vel(2))
 ang(2) = (pos(3)*vel(1) - pos(1)*vel(3))
 ang(3) = (pos(1)*vel(2) - pos(2)*vel(1))
 ang(:) = spinm*ang(:)

 write(ipafile,'(es18.10,1x,i18,1x,24(es18.10,1x))') time,sinki,ang(1),ang(2),ang(3),mpart,msink,spinm, &
      dx,dy,dz,dvx,dvy,dvz,xsink,ysink,zsink,xs,ys,zs, v_new,vxs,vys,vzs
 call flush(ipafile)

end subroutine track_accreted

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
       call write_inopt(r_crit,'r_crit','critical radius for point mass creation (no new sinks < r_crit from existing sink)',iunit)
       call write_inopt(h_acc, 'h_acc' ,'accretion radius for new sink particles',iunit)
       call write_inopt(h_soft_sinkgas,'h_soft_sinkgas','softening length for new sink particles', iunit)
    endif
 endif
 call write_inopt(h_soft_sinksink,'h_soft_sinksink','softening length between sink particles',iunit)
 call write_inopt(f_acc,'f_acc','particles < f_acc*h_acc accreted without checks',iunit)

end subroutine write_options_ptmass

!-----------------------------------------------------------------------
!+
!  reads sink particle options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ptmass(name,valstring,imatch,igotall,ierr)
 use io,         only:fatal
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
    ngot = ngot + 1
 case('h_acc')
    read(valstring,*,iostat=ierr) h_acc
    if (h_acc <= 0.) call fatal(label,'h_acc < 0')
    ngot = ngot + 1
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
 case default
    imatch = .false.
 end select

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
 if (icreate_sinks > 0) then
    igotall = (ngot >= 6)
 else
    igotall = (ngot >= 2)
 endif

end subroutine read_options_ptmass
!-----------------------------------------------------------------------
end module ptmass
