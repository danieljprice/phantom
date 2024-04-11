!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module step_extern
!
! Computes sub-steps in the RESPA algorithm
!
!   Multiple option of sub stepping can be choosed depending on
!   the physics and the precision needed
!
!   Only Hydro : step_extern_sph
!   Hydro + GR : step_extern_sph_gr step_extern_gr
!   2nd order with all fast physics implemented : step extern
!   4th order (Work in progress, only gravitionnal interaction
!   sink-sink and sink-gas) : step_extern_FSI step_extern_PEFRL
!
! :References:
!     Verlet (1967), Phys. Rev. 159, 98-103
!     Tuckerman, Berne & Martyna (1992), J. Chem. Phys. 97, 1990-2001
!     Rantala + (2020) (2023),Chin (2007a)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary_dyn, chem, cons2prim, cons2primsolver, cooling,
!   cooling_ism, damping, deriv, dim, dust_formation, eos, extern_gr,
!   externalforces, growth, io, io_summary, krome_interface, metric_tools,
!   mpiutils, options, part, ptmass, ptmass_radiation, timestep,
!   timestep_ind, timestep_sts, timing, units
!
 implicit none

 public :: step_extern_lf
 public :: step_extern_gr
 public :: step_extern_sph
 public :: step_extern_sph_gr
 public :: step_extern_FSI

 real,parameter :: dk(3) = (/1./6.,2./3.,1./6./)
 real,parameter :: ck(2) = (/0.5,0.5/)

 private

contains

subroutine step_extern_sph_gr(dt,npart,xyzh,vxyzu,dens,pxyzu,metrics)
 use part,            only:isdead_or_accreted,igas,massoftype,rhoh,eos_vars,igasP,&
                              ien_type,eos_vars,igamma,itemp
 use cons2primsolver, only:conservative2primitive
 use eos,             only:ieos
 use io,              only:warning
 use metric_tools,    only:pack_metric
 use timestep,        only:xtol
 real,    intent(in)    :: dt
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:),dens(:),metrics(:,:,:,:)
 real,    intent(in)    :: pxyzu(:,:)
 real,    intent(out)   :: vxyzu(:,:)
 integer, parameter :: nitermax = 50
 integer :: i,niter,ierr
 real    :: xpred(1:3),vold(1:3),diff
 logical :: converged
 real    :: rhoi,pri,tempi,gammai

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,dens,dt,xtol) &
 !$omp shared(pxyzu,metrics,ieos,massoftype,ien_type,eos_vars) &
 !$omp private(i,niter,diff,xpred,vold,converged,ierr) &
 !$omp private(pri,rhoi,tempi,gammai)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then

       !-- unpack and compute values for initial guess in cons2prim
       pri    = eos_vars(igasP,i)
       tempi  = eos_vars(itemp,i)
       gammai = eos_vars(igamma,i)
       rhoi   = rhoh(xyzh(4,i),massoftype(igas))

       call conservative2primitive(xyzh(1:3,i),metrics(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i),&
                                      pri,tempi,gammai,rhoi,pxyzu(1:3,i),pxyzu(4,i),ierr,ien_type)
       if (ierr > 0) call warning('cons2primsolver [in step_extern_sph_gr (a)]','enthalpy did not converge',i=i)
       !
       ! main position update
       !
       xpred = xyzh(1:3,i) + dt*vxyzu(1:3,i)
       vold  = vxyzu(1:3,i)
       converged = .false.
       niter = 0
       do while (.not. converged .and. niter<=nitermax)
          niter = niter + 1
          call conservative2primitive(xyzh(1:3,i),metrics(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i),&
                                         pri,tempi,gammai,rhoi,pxyzu(1:3,i),pxyzu(4,i),ierr,ien_type)
          if (ierr > 0) call warning('cons2primsolver [in step_extern_sph_gr (b)]','enthalpy did not converge',i=i)
          xyzh(1:3,i) = xpred + 0.5*dt*(vxyzu(1:3,i)-vold)
          diff = maxval(abs(xyzh(1:3,i)-xpred)/xpred)
          if (diff < xtol) converged = .true.
          ! UPDATE METRIC HERE
          call pack_metric(xyzh(1:3,i),metrics(:,:,:,i))
       enddo
       if (niter > nitermax) call warning('step_extern_sph_gr','Reached max number of x iterations. x_err ',val=diff)

       ! repack values
       eos_vars(igasP,i)  = pri
       eos_vars(itemp,i)  = tempi
       eos_vars(igamma,i) = gammai
    endif
 enddo
 !$omp end parallel do

end subroutine step_extern_sph_gr

subroutine step_extern_gr(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,time)
 use dim,            only:maxptmass,maxp,maxvxyzu
 use io,             only:iverbose,id,master,iprint,warning,fatal
 use externalforces, only:externalforce,accrete_particles,update_externalforce
 use options,        only:iexternalforce
 use part,           only:maxphase,isdead_or_accreted,iamboundary,igas,iphase,iamtype,&
                             massoftype,rhoh,ien_type,eos_vars,igamma,itemp,igasP
 use io_summary,     only:summary_variable,iosumextr,iosumextt,summary_accrete
 use timestep,       only:bignumber,C_force,xtol,ptol
 use eos,            only:equationofstate,ieos
 use cons2primsolver,only:conservative2primitive
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 use damping,        only:calc_damp,apply_damp,idamp
 integer, intent(in)    :: npart,ntypes
 real,    intent(in)    :: dtsph,time
 real,    intent(inout) :: dtextforce
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:),pxyzu(:,:),dens(:),metrics(:,:,:,:),metricderivs(:,:,:,:)
 integer :: i,itype,nsubsteps,naccreted,its,ierr,nlive
 real    :: timei,t_end_step,hdt,pmassi
 real    :: dt,dtf,dtextforcenew,dtextforce_min
 real    :: pri,spsoundi,pondensi,tempi,gammai
 real, save :: pprev(3),xyz_prev(3),fstar(3),vxyz_star(3),xyz(3),pxyz(3),vxyz(3),fexti(3)
 !$omp threadprivate(pprev,xyz_prev,fstar,vxyz_star,xyz,pxyz,vxyz,fexti)
 real    :: x_err,pmom_err,accretedmass,damp_fac
 ! real, save :: dmdt = 0.
 logical :: last_step,done,converged,accreted
 integer, parameter :: itsmax = 50
 integer :: pitsmax,xitsmax
 real    :: perrmax,xerrmax
 real :: rhoi,hi,eni,uui,densi

 pitsmax = 0
 xitsmax = 0
 perrmax = 0.
 xerrmax = 0.

 !
 ! determine whether or not to use substepping
 !
 if (dtextforce < dtsph) then
    dt = dtextforce
    last_step = .false.
 else
    dt = dtsph
    last_step = .true.
 endif

 timei = time
 itype          = igas
 pmassi         = massoftype(igas)
 t_end_step     = timei + dtsph
 nsubsteps      = 0
 dtextforce_min = huge(dt)
 done           = .false.
 substeps: do while (timei <= t_end_step .and. .not.done)
    hdt           = 0.5*dt
    timei         = timei + dt
    nsubsteps     = nsubsteps + 1
    dtextforcenew = bignumber

    call calc_damp(time, damp_fac)

    if (.not.last_step .and. iverbose > 1 .and. id==master) then
       write(iprint,"(a,f14.6)") '> external forces only : t=',timei
    endif
    !---------------------------
    ! predictor during substeps
    !---------------------------
    !
    ! predictor step for external forces, also recompute external forces
    !
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,vxyzu,fext,iphase,ntypes,massoftype) &
    !$omp shared(maxphase,maxp,eos_vars) &
    !$omp shared(dt,hdt,xtol,ptol) &
    !$omp shared(ieos,pxyzu,dens,metrics,metricderivs,ien_type) &
    !$omp private(i,its,spsoundi,tempi,rhoi,hi,eni,uui,densi) &
    !$omp private(converged,pmom_err,x_err,pri,ierr,gammai) &
    !$omp firstprivate(pmassi,itype) &
    !$omp reduction(max:xitsmax,pitsmax,perrmax,xerrmax) &
    !$omp reduction(min:dtextforcenew)
    predictor: do i=1,npart
       xyz(1) = xyzh(1,i)
       xyz(2) = xyzh(2,i)
       xyz(3) = xyzh(3,i)
       hi = xyzh(4,i)
       if (.not.isdead_or_accreted(hi)) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             pmassi = massoftype(itype)
          endif

          its       = 0
          converged = .false.
          !
          ! make local copies of array quantities
          !
          pxyz(1:3) = pxyzu(1:3,i)
          eni       = pxyzu(4,i)
          vxyz(1:3) = vxyzu(1:3,i)
          uui       = vxyzu(4,i)
          fexti     = fext(:,i)

          pxyz      = pxyz + hdt*fexti

          !-- unpack thermo variables for the first guess in cons2prim
          densi     = dens(i)
          pri       = eos_vars(igasP,i)
          gammai    = eos_vars(igamma,i)
          tempi     = eos_vars(itemp,i)
          rhoi      = rhoh(hi,massoftype(igas))

          ! Note: grforce needs derivatives of the metric,
          ! which do not change between pmom iterations
          pmom_iterations: do while (its <= itsmax .and. .not. converged)
             its   = its + 1
             pprev = pxyz
             call conservative2primitive(xyz,metrics(:,:,:,i),vxyz,densi,uui,pri,&
                                            tempi,gammai,rhoi,pxyz,eni,ierr,ien_type)
             if (ierr > 0) call warning('cons2primsolver [in step_extern_gr (a)]','enthalpy did not converge',i=i)
             call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i),vxyz,densi,uui,pri,fstar)
             pxyz = pprev + hdt*(fstar - fexti)
             pmom_err = maxval(abs(pxyz - pprev))
             if (pmom_err < ptol) converged = .true.
             fexti = fstar
          enddo pmom_iterations
          if (its > itsmax ) call warning('step_extern_gr',&
                                    'max # of pmom iterations',var='pmom_err',val=pmom_err)
          pitsmax = max(its,pitsmax)
          perrmax = max(pmom_err,perrmax)

          call conservative2primitive(xyz,metrics(:,:,:,i),vxyz,densi,uui,pri,tempi,&
                                         gammai,rhoi,pxyz,eni,ierr,ien_type)
          if (ierr > 0) call warning('cons2primsolver [in step_extern_gr (b)]','enthalpy did not converge',i=i)
          xyz = xyz + dt*vxyz
          call pack_metric(xyz,metrics(:,:,:,i))

          its        = 0
          converged  = .false.
          vxyz_star = vxyz
          ! Note: since particle positions change between iterations
          !  the metric and its derivatives need to be updated.
          !  cons2prim does not require derivatives of the metric,
          !  so those can updated once the iterations are complete
          !  in order to reduce the number of computations.
          xyz_iterations: do while (its <= itsmax .and. .not. converged)
             its         = its+1
             xyz_prev    = xyz
             call conservative2primitive(xyz,metrics(:,:,:,i),vxyz_star,densi,uui,&
                                            pri,tempi,gammai,rhoi,pxyz,eni,ierr,ien_type)
             if (ierr > 0) call warning('cons2primsolver [in step_extern_gr (c)]','enthalpy did not converge',i=i)
             xyz  = xyz_prev + hdt*(vxyz_star - vxyz)
             x_err = maxval(abs(xyz-xyz_prev))
             if (x_err < xtol) converged = .true.
             vxyz = vxyz_star
             ! UPDATE METRIC HERE
             call pack_metric(xyz,metrics(:,:,:,i))
          enddo xyz_iterations
          call pack_metricderivs(xyz,metricderivs(:,:,:,i))
          if (its > itsmax ) call warning('step_extern_gr','Reached max number of x iterations. x_err ',val=x_err)
          xitsmax = max(its,xitsmax)
          xerrmax = max(x_err,xerrmax)

          ! re-pack arrays back where they belong
          xyzh(1:3,i) = xyz(1:3)
          pxyzu(1:3,i) = pxyz(1:3)
          vxyzu(1:3,i) = vxyz(1:3)
          vxyzu(4,i) = uui
          fext(:,i)  = fexti
          dens(i) = densi
          eos_vars(igasP,i)  = pri
          eos_vars(itemp,i)  = tempi
          eos_vars(igamma,i) = gammai

          ! Skip remainder of update if boundary particle; note that fext==0 for these particles
          if (iamboundary(itype)) cycle predictor
       endif
    enddo predictor
    !$omp end parallel do

    if (iverbose >= 2 .and. id==master) then
       write(iprint,*)                '------ Iterations summary: -------------------------------'
       write(iprint,"(a,i2,a,f14.6)") 'Most pmom iterations = ',pitsmax,' | max error = ',perrmax
       write(iprint,"(a,i2,a,f14.6)") 'Most xyz  iterations = ',xitsmax,' | max error = ',xerrmax
       write(iprint,*)
    endif

    !
    ! corrector step on gas particles (also accrete particles at end of step)
    !
    accretedmass = 0.
    naccreted    = 0
    nlive = 0
    dtextforce_min = bignumber
    !$omp parallel default(none) &
    !$omp shared(npart,xyzh,metrics,metricderivs,vxyzu,fext,iphase,ntypes,massoftype,hdt,timei) &
    !$omp shared(maxphase,maxp) &
    !$omp private(i,accreted) &
    !$omp shared(ieos,dens,pxyzu,iexternalforce,C_force) &
    !$omp private(pri,pondensi,spsoundi,tempi,dtf) &
    !$omp firstprivate(itype,pmassi) &
    !$omp reduction(min:dtextforce_min) &
    !$omp reduction(+:accretedmass,naccreted,nlive) &
    !$omp shared(idamp,damp_fac)
    !$omp do
    accreteloop: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             pmassi = massoftype(itype)
             !  if (itype==iboundary) cycle accreteloop
          endif

          call equationofstate(ieos,pondensi,spsoundi,dens(i),xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
          pri = pondensi*dens(i)
          call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i),pri,fext(1:3,i),dtf)
          dtextforce_min = min(dtextforce_min,C_force*dtf)

          if (idamp > 0) then
             call apply_damp(fext(1,i), fext(2,i), fext(3,i), vxyzu(1:3,i), xyzh(1:3,i), damp_fac)
          endif

          !
          ! correct v to the full step using only the external force
          !
          pxyzu(1:3,i) = pxyzu(1:3,i) + hdt*fext(1:3,i)
          ! Do we need call cons2prim here ??

          if (iexternalforce > 0) then
             call accrete_particles(iexternalforce,xyzh(1,i),xyzh(2,i), &
                                       xyzh(3,i),xyzh(4,i),pmassi,timei,accreted,i)
             if (accreted) then
                accretedmass = accretedmass + pmassi
                naccreted = naccreted + 1
             endif
          endif
          nlive = nlive + 1
       endif
    enddo accreteloop
    !$omp enddo
    !$omp end parallel

    if (npart > 2 .and. nlive < 2) then
       call fatal('step','all particles accreted',var='nlive',ival=nlive)
    endif

    if (iverbose >= 2 .and. id==master .and. naccreted /= 0) write(iprint,"(a,es10.3,a,i4,a)") &
          'Step: at time ',timei,', ',naccreted,' particles were accreted. Mass accreted = ',accretedmass

    dtextforcenew = min(dtextforce_min,dtextforcenew)
    dtextforce    = dtextforcenew

    if (last_step) then
       done = .true.
    else
       dt = dtextforce
       if (timei + dt > t_end_step) then
          dt = t_end_step - timei
          last_step = .true.
       endif
    endif

 enddo substeps

 if (nsubsteps > 1) then
    if (iverbose>=1 .and. id==master) then
       write(iprint,"(a,i6,a,f8.2,a,es10.3,a,es10.3)") &
              ' using ',nsubsteps,' substeps (dthydro/dtextf = ',dtsph/dtextforce_min,'), dt = ',dtextforce_min,' dtsph = ',dtsph
    endif
    call summary_variable('ext',iosumextr ,nsubsteps,dtsph/dtextforce_min)
    call summary_variable('ext',iosumextt ,nsubsteps,dtextforce_min,1.0/dtextforce_min)
 endif

end subroutine step_extern_gr

 !----------------------------------------------------------------
 !+
 !  This is the equivalent of the routine below when no external
 !  forces, sink particles or cooling are used
 !+
 !----------------------------------------------------------------
subroutine step_extern_sph(dt,npart,xyzh,vxyzu)
 use part, only:isdead_or_accreted
 real,    intent(in)    :: dt
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: vxyzu(:,:)
 integer :: i

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,dt) &
 !$omp private(i)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       !
       ! main position update
       !
       xyzh(1,i) = xyzh(1,i) + dt*vxyzu(1,i)
       xyzh(2,i) = xyzh(2,i) + dt*vxyzu(2,i)
       xyzh(3,i) = xyzh(3,i) + dt*vxyzu(3,i)
    endif
 enddo
 !$omp end parallel do

end subroutine step_extern_sph

 !----------------------------------------------------------------
 !+
 !  This is the equivalent of the routine below with no cooling
 !  and external forces except ptmass. (4th order scheme)
 !+
 !----------------------------------------------------------------
subroutine step_extern_FSI(dtextforce,dtsph,time,npart,nptmass,xyzh,vxyzu,fext, &
                           xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fsink_old,dsdt_ptmass)
 use part,           only: isdead_or_accreted,igas,massoftype
 use io,             only:iverbose,id,master,iprint,warning,fatal
 use io_summary,     only:summary_variable,iosumextr,iosumextt
 real,    intent(in)    :: dtsph,time
 integer, intent(in)    :: npart,nptmass
 real,    intent(inout) :: dtextforce
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(4,nptmass),dsdt_ptmass(3,nptmass),fsink_old(4,nptmass)
 real    :: dt,t_end_step,dtextforce_min
 real    :: pmassi,timei
 logical :: done,last_step
 integer :: nsubsteps

 !
 ! determine whether or not to use substepping
 !
 if (dtextforce < dtsph) then
    dt = dtextforce
    last_step = .false.
 else
    dt = dtsph
    last_step = .true.
 endif

 timei = time
 pmassi         = massoftype(igas)
 t_end_step     = timei + dtsph
 nsubsteps      = 0
 dtextforce_min = huge(dt)
 done           = .false.

 substeps: do while (timei <= t_end_step .and. .not.done)
    timei = timei + dt
    if (abs(dt) < tiny(0.)) call fatal('step_extern','dt <= 0 in sink-gas substepping',var='dt',val=dt)
    nsubsteps = nsubsteps + 1
    call get_force_4th(nptmass,npart,nsubsteps,pmassi,timei,dtextforce,&
                      xyzh,fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass)
    call kick_4th (dk(1),dt,npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,fext,fxyz_ptmass,dsdt_ptmass)
    call drift_4th(ck(1),dt,npart,nptmass,xyzh,xyzmh_ptmass,vxyzu,vxyz_ptmass,dsdt_ptmass)
    call get_force_4th(nptmass,npart,nsubsteps,pmassi,timei,dtextforce,&
                      xyzh,fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass) ! Direct calculation of the force and force gradient
    fsink_old=fxyz_ptmass
    call get_gradf_extrap_4th(nptmass,npart,nsubsteps,pmassi,timei,dtextforce,&
                             xyzh,fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,fsink_old) ! extrapolation method Omelyan
    !call get_gradf_4th(nptmass,npart,pmassi,dt,xyzh,fext,xyzmh_ptmass,fxyz_ptmass,fsink_old)
    call kick_4th (dk(2),dt,npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,fext,fxyz_ptmass,dsdt_ptmass)
    call drift_4th(ck(2),dt,npart,nptmass,xyzh,xyzmh_ptmass,vxyzu,vxyz_ptmass,dsdt_ptmass)
    call get_force_4th(nptmass,npart,nsubsteps,pmassi,timei,dtextforce,&
                      xyzh,fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass)
    call kick_4th (dk(3),dt,npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,fext,fxyz_ptmass,dsdt_ptmass)
    if (iverbose >= 2 ) write(iprint,*) "nsubsteps : ",nsubsteps,"time,dt : ",timei,dt

    dtextforce_min = min(dtextforce_min,dtextforce)

    if (last_step) then
       done = .true.
    else
       dt = dtextforce
       if (timei + dt > t_end_step) then
          dt = t_end_step - timei
          last_step = .true.
       endif
    endif
 enddo substeps

 if (nsubsteps > 1) then
    if (iverbose>=1 .and. id==master) then
       write(iprint,"(a,i6,a,f8.2,a,es10.3,a,es10.3)") &
              ' using ',nsubsteps,' substeps (dthydro/dtextf = ',dtsph/dtextforce_min,'), dt = ',dtextforce_min,' dtsph = ',dtsph
    endif
    call summary_variable('ext',iosumextr ,nsubsteps,dtsph/dtextforce_min)
    call summary_variable('ext',iosumextt ,nsubsteps,dtextforce_min,1.0/dtextforce_min)
 endif


end subroutine step_extern_FSI


 !----------------------------------------------------------------
 !+
 !  drift routine for the 4th order scheme
 !+
 !----------------------------------------------------------------

subroutine drift_4th(ck,dt,npart,nptmass,xyzh,xyzmh_ptmass,vxyzu,vxyz_ptmass,dsdt_ptmass)
 use part,     only:isdead_or_accreted,ispinx,ispiny,ispinz
 use io  ,     only:id,master
 use mpiutils, only:bcast_mpi
 real,    intent(in)    :: dt,ck
 integer, intent(in)    :: npart,nptmass
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),dsdt_ptmass(:,:)
 real    :: ckdt
 integer :: i

 ckdt = ck*dt

 ! Drift gas particles

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,ckdt) &
 !$omp private(i)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       xyzh(1,i) = xyzh(1,i) + ckdt*vxyzu(1,i)
       xyzh(2,i) = xyzh(2,i) + ckdt*vxyzu(2,i)
       xyzh(3,i) = xyzh(3,i) + ckdt*vxyzu(3,i)
    endif
 enddo
 !$omp end parallel do

 ! Drift sink particles
 if(nptmass>0) then
    if(id==master) then
       !$omp parallel do default(none) &
       !$omp shared(nptmass,xyzmh_ptmass,vxyz_ptmass,dsdt_ptmass,ckdt) &
       !$omp private(i)
       do i=1,nptmass
          if (xyzmh_ptmass(4,i) > 0.) then
             xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + ckdt*vxyz_ptmass(1,i)
             xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + ckdt*vxyz_ptmass(2,i)
             xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + ckdt*vxyz_ptmass(3,i)
             xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) + ckdt*dsdt_ptmass(1,i)
             xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) + ckdt*dsdt_ptmass(2,i)
             xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) + ckdt*dsdt_ptmass(3,i)
          endif
       enddo
       !$omp end parallel do
    endif
    call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
 endif
end subroutine drift_4th


 !----------------------------------------------------------------
 !+
 !  kick routine for the 4th order scheme
 !+
 !----------------------------------------------------------------

subroutine kick_4th(dk,dt,npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,fext,fxyz_ptmass,dsdt_ptmass)
 use part, only: isdead_or_accreted,ispinx,ispiny,ispinz
 use io  ,     only:id,master
 use mpiutils, only:bcast_mpi
 real,    intent(in)    :: dt,dk
 integer, intent(in)    :: npart,nptmass
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:),fext(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),dsdt_ptmass(:,:)
 integer :: i
 real    :: dkdt

 dkdt = dk*dt

 ! Kick gas particles

 !$omp parallel do default(none) &
 !$omp shared(npart,fext,xyzh,vxyzu,dkdt) &
 !$omp private(i)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       vxyzu(1,i) = vxyzu(1,i) + dkdt*fext(1,i)
       vxyzu(2,i) = vxyzu(2,i) + dkdt*fext(2,i)
       vxyzu(3,i) = vxyzu(3,i) + dkdt*fext(3,i)
    endif
 enddo
 !$omp end parallel do

 ! Kick sink particles
 if (nptmass>0) then
    if(id==master) then
       !$omp parallel do default(none) &
       !$omp shared(nptmass,xyzmh_ptmass,fxyz_ptmass,vxyz_ptmass,dsdt_ptmass,dkdt) &
       !$omp private(i)
       do i=1,nptmass
          if (xyzmh_ptmass(4,i) > 0.) then
             vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + dkdt*fxyz_ptmass(1,i)
             vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + dkdt*fxyz_ptmass(2,i)
             vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + dkdt*fxyz_ptmass(3,i)
             xyzmh_ptmass(ispinx,i) = xyzmh_ptmass(ispinx,i) + dkdt*dsdt_ptmass(1,i)
             xyzmh_ptmass(ispiny,i) = xyzmh_ptmass(ispiny,i) + dkdt*dsdt_ptmass(2,i)
             xyzmh_ptmass(ispinz,i) = xyzmh_ptmass(ispinz,i) + dkdt*dsdt_ptmass(3,i)
          endif
       enddo
       !$omp end parallel do

    endif
    call bcast_mpi(vxyz_ptmass(:,1:nptmass))
 endif

end subroutine kick_4th

 !----------------------------------------------------------------
 !+
 !  force routine for the 4th order scheme
 !+
 !----------------------------------------------------------------

subroutine get_force_4th(nptmass,npart,nsubsteps,pmassi,timei,dtextforce,xyzh,fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass)
 use options,        only:iexternalforce
 use dim,            only:maxptmass
 use io,             only:iverbose,master,id,iprint,warning,fatal
 use part,           only:epot_sinksink,fxyz_ptmass_sinksink,dsdt_ptmass_sinksink
 use ptmass,         only:get_accel_sink_gas,get_accel_sink_sink,merge_sinks
 use timestep,       only:bignumber,C_force
 use mpiutils,       only:bcast_mpi,reduce_in_place_mpi,reduceall_mpi
 integer,           intent(in)    :: nptmass,npart,nsubsteps
 real,              intent(inout) :: xyzh(:,:),fext(:,:)
 real,              intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(4,nptmass),dsdt_ptmass(3,nptmass)
 real,              intent(inout) :: dtextforce
 real,              intent(in)    :: timei,pmassi
 integer         :: merge_ij(nptmass)
 integer         :: merge_n
 integer         :: i
 real            :: dtf,dtextforcenew,dtsinkgas,dtphi2,fonrmax
 real            :: fextx,fexty,fextz
 real            :: fonrmaxi,phii,dtphi2i

 dtextforcenew = bignumber
 dtsinkgas     = bignumber
 dtphi2        = bignumber
 fonrmax = 0

 if (nptmass>0) then
    if (id==master) then
       call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                               dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass)
       if (merge_n > 0) then
          call merge_sinks(timei,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,merge_ij)
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                  dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass)
          fxyz_ptmass_sinksink=fxyz_ptmass
          dsdt_ptmass_sinksink=dsdt_ptmass
          if (iverbose >= 2) write(iprint,*) 'dt(sink-sink) = ',C_force*dtf
       endif
    else
       fxyz_ptmass(:,:) = 0.
       dsdt_ptmass(:,:) = 0.
    endif
    call bcast_mpi(epot_sinksink)
    call bcast_mpi(dtf)
    dtextforcenew = min(dtextforcenew,C_force*dtf)
 endif
 if (iverbose >= 3 ) write(iprint,*) "dt_sink_sink",dtextforcenew
 !$omp parallel default(none) &
 !$omp shared(npart,nptmass,xyzh,xyzmh_ptmass,fext) &
 !$omp private(fextx,fexty,fextz) &
 !$omp private(fonrmaxi,dtphi2i,phii,pmassi,dtf) &
 !$omp reduction(min:dtextforcenew,dtphi2) &
 !$omp reduction(max:fonrmax) &
 !$omp reduction(+:fxyz_ptmass,dsdt_ptmass)
 !$omp do
 do i=1,npart
    fextx = 0.
    fexty = 0.
    fextz = 0.
    call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                  fextx,fexty,fextz,phii,pmassi,fxyz_ptmass,dsdt_ptmass,fonrmaxi,dtphi2i)
    fonrmax = max(fonrmax,fonrmaxi)
    dtphi2  = min(dtphi2,dtphi2i)
    fext(1,i) = fextx
    fext(2,i) = fexty
    fext(3,i) = fextz
 enddo
 !$omp enddo
 !$omp end parallel

 if (nptmass > 0) then
    call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))
    call reduce_in_place_mpi('+',dsdt_ptmass(:,1:nptmass))
 endif

 if(nptmass>0) then
    if (fonrmax > 0.) then
       dtsinkgas = min(dtsinkgas,C_force*1./sqrt(fonrmax),C_force*sqrt(dtphi2))
    endif
    if (iverbose >= 3 ) write(iprint,*) nsubsteps,'dt,(ext/sink-sink) = ',dtextforcenew,', dt(sink-gas) = ',dtsinkgas
    dtextforcenew = min(dtextforcenew,dtsinkgas)
    dtextforce = dtextforcenew
 endif

 dtextforcenew = reduceall_mpi('min',dtextforcenew)

end subroutine get_force_4th


 !----------------------------------------------------------------
 !+
 !  grad routine for the 4th order scheme (FSI)
 !+
 !----------------------------------------------------------------


subroutine get_gradf_4th(nptmass,npart,pmassi,dt,xyzh,fext,xyzmh_ptmass,fxyz_ptmass,fsink_old)
 use dim,            only:maxptmass
 use ptmass,         only:get_gradf_sink_gas,get_gradf_sink_sink
 use mpiutils,       only:reduce_in_place_mpi
 use io,             only:id,master
 integer, intent(in) :: nptmass,npart
 real, intent(inout) :: xyzh(:,:),fext(:,:)
 real, intent(inout) :: xyzmh_ptmass(:,:),fxyz_ptmass(4,nptmass)
 real, intent(in)    :: fsink_old(4,nptmass)
 real, intent(inout) :: dt
 real, intent(in)    :: pmassi
 real            :: fextx,fexty,fextz
 integer         :: i


 if (nptmass>0) then
    if(id==master) then
       call get_gradf_sink_sink(nptmass,dt,xyzmh_ptmass,fxyz_ptmass,fsink_old)
    else
       fxyz_ptmass(:,:) = 0.
    endif
 endif

 !$omp parallel default(none) &
 !$omp shared(npart,nptmass,xyzh,xyzmh_ptmass,fext,dt,pmassi,fsink_old) &
 !$omp private(fextx,fexty,fextz) &
 !$omp reduction(+:fxyz_ptmass)
 !$omp do
 do i=1,npart
    fextx = fext(1,i)
    fexty = fext(2,i)
    fextz = fext(3,i)
    call get_gradf_sink_gas(nptmass,dt,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                xyzmh_ptmass,fextx,fexty,fextz,pmassi,fxyz_ptmass,fsink_old)
    fext(1,i) = fext(1,i)+ fextx
    fext(2,i) = fext(2,i)+ fexty
    fext(3,i) = fext(3,i)+ fextz
 enddo
 !$omp enddo
 !$omp end parallel

 if (nptmass > 0) then
    call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))
    !call reduce_in_place_mpi('+',dsdt_ptmass(:,1:nptmass))
 endif

end subroutine get_gradf_4th

 !----------------------------------------------------------------
 !+
 !  grad routine for the 4th order scheme (FSI), extrapolation method
 !+
 !----------------------------------------------------------------


subroutine get_gradf_extrap_4th(nptmass,npart,nsubsteps,pmassi,timei,dtextforce,xyzh, &
                               fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,fsink_old)
 use options,        only:iexternalforce
 use dim,            only:maxptmass
 use part,           only:epot_sinksink
 use io,             only:master,id
 use ptmass,         only:get_accel_sink_gas,get_accel_sink_sink,merge_sinks
 use timestep,       only:bignumber
 use mpiutils,       only:reduce_in_place_mpi
 integer,           intent(in)    :: nptmass,npart,nsubsteps
 real,              intent(inout) :: xyzh(:,:),fext(:,:)
 real,              intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(4,nptmass),dsdt_ptmass(3,nptmass)
 real,              intent(in)    :: fsink_old(4,nptmass)
 real,              intent(inout) :: dtextforce
 real,              intent(in)    :: timei,pmassi,dt
 integer         :: merge_ij(nptmass)
 integer         :: merge_n
 integer         :: i
 real            :: dtf
 real            :: fextx,fexty,fextz,xi,yi,zi
 real            :: fonrmaxi,phii,dtphi2i,extrapfac

 extrapfac = (1/24.)*dt**2

 if (nptmass>0) then
    if (id==master) then
       call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                 dtf,iexternalforce,timei,merge_ij,merge_n, &
                                 dsdt_ptmass,extrapfac,fsink_old)
       if (merge_n > 0) then
          call merge_sinks(timei,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,merge_ij)
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                    dtf,iexternalforce,timei,merge_ij,merge_n, &
                                    dsdt_ptmass,extrapfac,fsink_old)
       endif
    else
       fxyz_ptmass(:,:) = 0.
       dsdt_ptmass(:,:) = 0.
    endif
 endif


 !$omp parallel default(none) &
 !$omp shared(npart,nptmass,xyzh,xyzmh_ptmass,fext,extrapfac,fsink_old) &
 !$omp private(fextx,fexty,fextz,xi,yi,zi) &
 !$omp private(fonrmaxi,dtphi2i,phii,pmassi,dtf) &
 !$omp reduction(+:fxyz_ptmass,dsdt_ptmass)
 !$omp do
 do i=1,npart
    fextx = 0.
    fexty = 0.
    fextz = 0.
    xi = xyzh(1,i) + extrapfac*fext(1,i)
    yi = xyzh(2,i) + extrapfac*fext(2,i)
    zi = xyzh(3,i) + extrapfac*fext(3,i)
    call get_accel_sink_gas(nptmass,xi,yi,zi,xyzh(4,i),xyzmh_ptmass,&
                    fextx,fexty,fextz,phii,pmassi,fxyz_ptmass, &
                    dsdt_ptmass,fonrmaxi,dtphi2i,extrapfac,fsink_old)
    fext(1,i) = fextx
    fext(2,i) = fexty
    fext(3,i) = fextz
 enddo
 !$omp enddo
 !$omp end parallel

 if (nptmass > 0) then
    call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))
    call reduce_in_place_mpi('+',dsdt_ptmass(:,1:nptmass))
 endif
end subroutine get_gradf_extrap_4th



! NOTE: The chemistry and cooling here is implicitly calculated.  That is,
!       dt is *passed in* to the chemistry & cooling routines so that the
!       output will be at the correct time of time + dt.  Since this is
!       implicit, there is no cooling timestep.  Explicit cooling is
!       calculated in force and requires a cooling timestep.

subroutine cooling_abundances_update(i,pmassi,xyzh,vxyzu,eos_vars,abundance,nucleation,dust_temp, &
                                     divcurlv,abundc,abunde,abundo,abundsi,dt,dphot0,idK2,idmu,idkappa, &
                                     idgamma,imu,igamma,nabn,dphotflag,nabundances)
 use dim,             only:h2chemistry,do_nucleation,use_krome,update_muGamma,store_dust_temperature
 use options,         only:icooling
 use chem,            only:update_abundances,get_dphot
 use dust_formation,  only:evolve_dust
 use cooling,         only:energ_cooling,cooling_in_step
 use part,            only:rhoh
#ifdef KROME
 use part,            only: T_gas_cool
 use krome_interface, only: update_krome
#endif
 real,         intent(inout) :: vxyzu(:,:),xyzh(:,:)
 real,         intent(inout) :: eos_vars(:,:),abundance(:,:)
 real,         intent(inout) :: nucleation(:,:),dust_temp(:)
 real(kind=4), intent(in)    :: divcurlv(:,:)
 real,         intent(inout) :: abundc,abunde,abundo,abundsi
 real,         intent(in)    :: dt,dphot0,pmassi
 integer,      intent(in)    :: idK2,idmu,idkappa,idgamma,imu,igamma
 integer,      intent(in)    :: i,nabn,dphotflag,nabundances

 real :: dudtcool,rhoi,ui,dphot
 real :: abundi(nabn)

 dudtcool = 0.
 rhoi = rhoh(xyzh(4,i),pmassi)
 !
 ! CHEMISTRY
 !
 if (h2chemistry) then
    !
    ! Get updated abundances of all species, updates 'chemarrays',
    !
    dphot = get_dphot(dphotflag,dphot0,xyzh(1,i),xyzh(2,i),xyzh(3,i))
    call update_abundances(vxyzu(4,i),rhoi,abundance(:,i),nabundances,&
               dphot,dt,abundi,nabn,eos_vars(imu,i),abundc,abunde,abundo,abundsi)
 endif
#ifdef KROME
 ! evolve chemical composition and determine new internal energy
 ! Krome also computes cooling function but only associated with chemical processes
 ui = vxyzu(4,i)
 call update_krome(dt,xyzh(:,i),ui,rhoi,abundance(:,i),eos_vars(igamma,i),eos_vars(imu,i),T_gas_cool(i))
 dudtcool = (ui-vxyzu(4,i))/dt
#else
 !evolve dust chemistry and compute dust cooling
 if (do_nucleation) then
    call evolve_dust(dt, xyzh(:,i), vxyzu(4,i), nucleation(:,i), dust_temp(i), rhoi)
    eos_vars(imu,i)    = nucleation(idmu,i)
    eos_vars(igamma,i) = nucleation(idgamma,i)
 endif
 !
 ! COOLING
 !
 if (icooling > 0 .and. cooling_in_step) then
    if (h2chemistry) then
       !
       ! Call cooling routine, requiring total density, some distance measure and
       ! abundances in the 'abund' format
       !
       call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,&
                 dust_temp(i),eos_vars(imu,i), eos_vars(igamma,i),abund_in=abundi)
    elseif (store_dust_temperature) then
       ! cooling with stored dust temperature
       if (do_nucleation) then
          call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,&
                    dust_temp(i),nucleation(idmu,i),nucleation(idgamma,i),nucleation(idK2,i),nucleation(idkappa,i))
       elseif (update_muGamma) then
          call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,&
                    dust_temp(i),eos_vars(imu,i), eos_vars(igamma,i))
       else
          call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,dust_temp(i))
       endif
    else
       ! cooling without stored dust temperature
       call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool)
    endif
 endif
#endif
 ! update internal energy
 if (cooling_in_step .or. use_krome) vxyzu(4,i) = vxyzu(4,i) + dt * dudtcool


end subroutine cooling_abundances_update



subroutine external_force_update(xi,yi,zi,hi,vxi,vyi,vzi,timei,i,dtextforcenew,dtf,dt, &
                                 fextx,fexty,fextz,extf_is_velocity_dependent,iexternalforce)
 use timestep,       only:C_force
 use externalforces, only: externalforce,update_vdependent_extforce_leapfrog
 real,    intent(in) :: xi,yi,zi,hi,vxi,vyi,vzi,timei,dt
 real, intent(inout) :: dtextforcenew,dtf,fextx,fexty,fextz
 integer, intent(in) :: iexternalforce,i
 logical, intent(in) :: extf_is_velocity_dependent
 real :: fextxi,fextyi,fextzi,poti
 real :: fextv(3)

 call externalforce(iexternalforce,xi,yi,zi,hi, &
   timei,fextxi,fextyi,fextzi,poti,dtf,i)
 dtextforcenew = min(dtextforcenew,C_force*dtf)

 fextx = fextx + fextxi
 fexty = fexty + fextyi
 fextz = fextz + fextzi
!
!  Velocity-dependent external forces require special handling
!  in leapfrog (corrector is implicit)
!
 if (extf_is_velocity_dependent) then
    fextxi = fextx
    fextyi = fexty
    fextzi = fextz
    call update_vdependent_extforce_leapfrog(iexternalforce,vxi,vyi,vzi, &
                                             fextxi,fextyi,fextzi,fextv,dt,xi,yi,zi)
    fextx = fextx + fextv(1)
    fexty = fexty + fextv(2)
    fextz = fextz + fextv(3)
 endif


end subroutine external_force_update



 !----------------------------------------------------------------
 !+
 !  Substepping of external and sink particle forces.
 !  Also updates position of all particles even if no external
 !  forces applied. This is the internal loop of the RESPA
 !  algorithm over the "fast" forces.
 !+
 !----------------------------------------------------------------
subroutine step_extern_lf(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,fext,fxyzu,time,nptmass, &
                          xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,nbinmax,ibin_wake)
 use dim,            only:maxptmass,maxp,maxvxyzu,store_dust_temperature,use_krome,itau_alloc,&
                            do_nucleation,update_muGamma,h2chemistry,ind_timesteps
 use io,             only:iverbose,id,master,iprint,warning,fatal
 use externalforces, only:externalforce,accrete_particles,update_externalforce, &
                             update_vdependent_extforce_leapfrog,is_velocity_dependent
 use ptmass,         only:ptmass_predictor,ptmass_corrector,ptmass_accrete, &
                             get_accel_sink_gas,get_accel_sink_sink,merge_sinks,f_acc,pt_write_sinkev, &
                             idxmsi,idymsi,idzmsi,idmsi,idspinxsi,idspinysi,idspinzsi, &
                             idvxmsi,idvymsi,idvzmsi,idfxmsi,idfymsi,idfzmsi, &
                             ndptmass,update_ptmass
 use options,        only:iexternalforce,icooling
 use part,           only:maxphase,abundance,nabundances,epot_sinksink,eos_vars,&
                             isdead_or_accreted,iamboundary,igas,iphase,iamtype,massoftype,rhoh,divcurlv, &
                             fxyz_ptmass_sinksink,dsdt_ptmass_sinksink,dust_temp,tau,&
                             nucleation,idK2,idmu,idkappa,idgamma,imu,igamma
 use chem,           only:update_abundances,get_dphot
 use cooling_ism,    only:dphot0,energ_cooling_ism,dphotflag,abundsi,abundo,abunde,abundc,nabn
 use io_summary,     only:summary_variable,iosumextr,iosumextt,summary_accrete,summary_accrete_fail
 use timestep,       only:bignumber,C_force
 use timestep_sts,   only:sts_it_n
 use mpiutils,       only:bcast_mpi,reduce_in_place_mpi,reduceall_mpi
 use damping,        only:calc_damp,apply_damp,idamp
 use ptmass_radiation,only:get_rad_accel_from_ptmass,isink_radiation
 use cooling,        only:energ_cooling,cooling_in_step
 use dust_formation, only:evolve_dust,calc_muGamma
 use units,          only:unit_density
#ifdef KROME
 use part,            only: T_gas_cool
 use krome_interface, only: update_krome
#endif
 integer,         intent(in)    :: npart,ntypes,nptmass
 real,            intent(in)    :: dtsph,time
 real,            intent(inout) :: dtextforce
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:),fxyzu(:,:)
 real,            intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),dsdt_ptmass(:,:)
 integer(kind=1), intent(in)    :: nbinmax
 integer(kind=1), intent(inout) :: ibin_wake(:)
 integer         :: i,itype,nsubsteps,naccreted,nfail,nfaili,merge_n,nlive
 integer         :: merge_ij(nptmass)
 integer(kind=1) :: ibin_wakei
 real            :: timei,hdt,fextx,fexty,fextz,fextxi,fextyi,fextzi,phii,pmassi
 real            :: dtphi2,dtphi2i,vxhalfi,vyhalfi,vzhalfi,fxi,fyi,fzi
 real            :: dudtcool,fextv(3),poti,ui,rhoi,mui,gammai,ph,ph_tot
 real            :: dt,dtextforcenew,dtsinkgas,fonrmax,fonrmaxi
 real            :: dtf,accretedmass,t_end_step,dtextforce_min
 real, allocatable :: dptmass(:,:) ! dptmass(ndptmass,nptmass)
 real            :: damp_fac,dphot
 real, save      :: dmdt = 0.
 real            :: abundi(nabn),gmwvar
 logical         :: accreted,extf_is_velocity_dependent
 logical         :: last_step,done


 !
 ! determine whether or not to use substepping
 !
 if (dtextforce < dtsph) then
    dt = dtextforce
    last_step = .false.
 else
    dt = dtsph
    last_step = .true.
 endif

 timei = time
 extf_is_velocity_dependent = is_velocity_dependent(iexternalforce)
 accretedmass   = 0.
 itype          = igas
 pmassi         = massoftype(igas)
 t_end_step     = timei + dtsph
 nsubsteps      = 0
 dtextforce_min = huge(dt)
 done           = .false.
 ! allocate memory for dptmass array (avoids ifort bug)
 allocate(dptmass(ndptmass,nptmass))

 substeps: do while (timei <= t_end_step .and. .not.done)
    hdt           = 0.5*dt
    timei         = timei + dt
    if (abs(dt) < tiny(0.)) call fatal('step_extern','dt <= 0 in sink-gas substepping',var='dt',val=dt)
    nsubsteps     = nsubsteps + 1
    dtextforcenew = bignumber
    dtsinkgas     = bignumber
    dtphi2        = bignumber

    call calc_damp(time, damp_fac)

    if (.not.last_step .and. iverbose > 1 .and. id==master) then
       write(iprint,"(a,f14.6)") '> external/ptmass forces only : t=',timei
    endif
    !
    ! update time-dependent external forces
    !
    call update_externalforce(iexternalforce,timei,dmdt)

    !---------------------------
    ! predictor during substeps
    !---------------------------
    !
    ! point mass predictor step
    !
    if (nptmass > 0) then
       if (id==master) then
          call ptmass_predictor(nptmass,dt,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass)
          !
          ! get sink-sink forces (and a new sink-sink timestep.  Note: fxyz_ptmass is zeroed in this subroutine)
          ! pass sink-sink forces to variable fxyz_ptmass_sinksink for later writing.
          !
          if (iexternalforce==14) call update_externalforce(iexternalforce,timei,dmdt)
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                      dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass)
          if (merge_n > 0) then
             call merge_sinks(timei,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,merge_ij)
             call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                         dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass)
          endif
          fxyz_ptmass_sinksink=fxyz_ptmass
          dsdt_ptmass_sinksink=dsdt_ptmass
          if (iverbose >= 2) write(iprint,*) 'dt(sink-sink) = ',C_force*dtf
       else
          fxyz_ptmass(:,:) = 0.
          dsdt_ptmass(:,:) = 0.
       endif
       call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
       call bcast_mpi(vxyz_ptmass(:,1:nptmass))
       call bcast_mpi(epot_sinksink)
       call bcast_mpi(dtf)
       dtextforcenew = min(dtextforcenew,C_force*dtf)
    endif

    !
    ! predictor step for sink-gas and external forces, also recompute sink-gas and external forces
    !
    fonrmax = 0.
    !$omp parallel default(none) &
    !$omp shared(maxp,maxphase) &
    !$omp shared(npart,xyzh,vxyzu,fext,abundance,iphase,ntypes,massoftype) &
    !$omp shared(eos_vars,dust_temp,store_dust_temperature) &
    !$omp shared(dt,hdt,timei,iexternalforce,extf_is_velocity_dependent,cooling_in_step,icooling) &
    !$omp shared(xyzmh_ptmass,vxyz_ptmass,idamp,damp_fac) &
    !$omp shared(nptmass,nsubsteps,C_force,divcurlv,dphotflag,dphot0) &
    !$omp shared(abundc,abundo,abundsi,abunde) &
    !$omp shared(nucleation,do_nucleation,update_muGamma,h2chemistry,unit_density) &
    !$omp private(dphot,abundi,gmwvar,ph,ph_tot) &
    !$omp private(ui,rhoi, mui, gammai) &
    !$omp private(i,dudtcool,fxi,fyi,fzi,phii) &
    !$omp private(fextx,fexty,fextz,fextxi,fextyi,fextzi,poti,fextv,accreted) &
    !$omp private(fonrmaxi,dtphi2i,dtf) &
    !$omp private(vxhalfi,vyhalfi,vzhalfi) &
    !$omp firstprivate(pmassi,itype) &
#ifdef KROME
    !$omp shared(T_gas_cool) &
#endif
    !$omp reduction(min:dtextforcenew,dtphi2) &
    !$omp reduction(max:fonrmax) &
    !$omp reduction(+:fxyz_ptmass,dsdt_ptmass)
    !$omp do
    predictor: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype  = iamtype(iphase(i))
             pmassi = massoftype(itype)
          endif
          !
          ! predict v to the half step
          !
          vxyzu(1:3,i) = vxyzu(1:3,i) + hdt*fext(1:3,i)
          !
          ! main position update
          !
          xyzh(1,i) = xyzh(1,i) + dt*vxyzu(1,i)
          xyzh(2,i) = xyzh(2,i) + dt*vxyzu(2,i)
          xyzh(3,i) = xyzh(3,i) + dt*vxyzu(3,i)
          !
          ! Skip remainder of update if boundary particle; note that fext==0 for these particles
          if (iamboundary(itype)) cycle predictor
          !
          ! compute and add sink-gas force
          !
          fextx = 0.
          fexty = 0.
          fextz = 0.
          if (nptmass > 0) then
             call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                         fextx,fexty,fextz,phii,pmassi,fxyz_ptmass,dsdt_ptmass,fonrmaxi,dtphi2i)
             fonrmax = max(fonrmax,fonrmaxi)
             dtphi2  = min(dtphi2,dtphi2i)
          endif
          !
          ! compute and add external forces
          !
          if (iexternalforce > 0) then
             call external_force_update(xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i), &
                                    vxyzu(1,i),vxyzu(1,i),vxyzu(1,i),timei,i, &
                                    dtextforcenew,dtf,dt,fextx,fexty,fextz, &
                                    extf_is_velocity_dependent,iexternalforce)
          endif

          if (idamp > 0) then
             call apply_damp(fextx, fexty, fextz, vxyzu(1:3,i), xyzh(1:3,i), damp_fac)
          endif

          if (maxvxyzu >= 4 .and. itype==igas) then
             call cooling_abundances_update(i,pmassi,xyzh,vxyzu,eos_vars,abundance,nucleation,dust_temp, &
                                            divcurlv,abundc,abunde,abundo,abundsi,dt,dphot0,idK2,idmu,idkappa, &
                                            idgamma,imu,igamma,nabn,dphotflag,nabundances)
          endif
          fext(1,i) = fextx
          fext(2,i) = fexty
          fext(3,i) = fextz
       endif
    enddo predictor
    !$omp enddo
    !$omp end parallel

    if (nptmass > 0 .and. isink_radiation > 0) then
       if (itau_alloc == 1) then
          call get_rad_accel_from_ptmass(nptmass,npart,xyzh,xyzmh_ptmass,fext,tau)
       else
          call get_rad_accel_from_ptmass(nptmass,npart,xyzh,xyzmh_ptmass,fext)
       endif
    endif

    !
    ! reduction of sink-gas forces from each MPI thread
    !
    if (nptmass > 0) then
       call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))
       call reduce_in_place_mpi('+',dsdt_ptmass(:,1:nptmass))
    endif
    !---------------------------
    ! corrector during substeps
    !---------------------------
    !
    ! corrector step on sinks (changes velocities only, does not change position)
    !
    if (nptmass > 0) then
       if (id==master) then
          call ptmass_corrector(nptmass,dt,vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,dsdt_ptmass,iexternalforce)
       endif
       call bcast_mpi(vxyz_ptmass(:,1:nptmass))
    endif

    !
    ! corrector step on gas particles (also accrete particles at end of step)
    !
    accretedmass = 0.
    nfail        = 0
    naccreted    = 0
    nlive        = 0
    ibin_wakei   = 0
    dptmass(:,:) = 0.

    !$omp parallel default(none) &
    !$omp shared(maxp,maxphase) &
    !$omp shared(npart,xyzh,vxyzu,fext,iphase,ntypes,massoftype,hdt,timei,nptmass,sts_it_n) &
    !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,f_acc) &
    !$omp shared(iexternalforce) &
    !$omp shared(nbinmax,ibin_wake) &
    !$omp reduction(+:dptmass) &
    !$omp private(i,accreted,nfaili,fxi,fyi,fzi) &
    !$omp firstprivate(itype,pmassi,ibin_wakei) &
    !$omp reduction(+:accretedmass,nfail,naccreted,nlive)
    !$omp do
    accreteloop: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             pmassi = massoftype(itype)
             if (iamboundary(itype)) cycle accreteloop
          endif
          !
          ! correct v to the full step using only the external force
          !
          vxyzu(1:3,i) = vxyzu(1:3,i) + hdt*fext(1:3,i)

          if (iexternalforce > 0) then
             call accrete_particles(iexternalforce,xyzh(1,i),xyzh(2,i), &
                                       xyzh(3,i),xyzh(4,i),pmassi,timei,accreted)
             if (accreted) accretedmass = accretedmass + pmassi
          endif
          !
          ! accretion onto sink particles
          ! need position, velocities and accelerations of both gas and sinks to be synchronised,
          ! otherwise will not conserve momentum
          ! Note: requiring sts_it_n since this is supertimestep with the most active particles
          !
          if (nptmass > 0 .and. sts_it_n) then
             fxi = fext(1,i)
             fyi = fext(2,i)
             fzi = fext(3,i)
             if (ind_timesteps) ibin_wakei = ibin_wake(i)

             call ptmass_accrete(1,nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                                    vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),fxi,fyi,fzi,&
                                    itype,pmassi,xyzmh_ptmass,vxyz_ptmass,&
                                    accreted,dptmass,timei,f_acc,nbinmax,ibin_wakei,nfaili)
             if (accreted) then
                naccreted = naccreted + 1
                cycle accreteloop
             else
                if (ind_timesteps) ibin_wake(i) = ibin_wakei
             endif
             if (nfaili > 1) nfail = nfail + 1
          endif
          nlive = nlive + 1
       endif
    enddo accreteloop
    !$omp enddo
    !$omp end parallel

    if (npart > 2 .and. nlive < 2) then
       call fatal('step','all particles accreted',var='nlive',ival=nlive)
    endif

    !
    ! reduction of sink particle changes across MPI
    !
    if (nptmass > 0) then
       call reduce_in_place_mpi('+',dptmass(:,1:nptmass))

       naccreted = int(reduceall_mpi('+',naccreted))
       nfail = int(reduceall_mpi('+',nfail))

       if (id==master) call update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)

       call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
       call bcast_mpi(vxyz_ptmass(:,1:nptmass))
       call bcast_mpi(fxyz_ptmass(:,1:nptmass))
    endif

    if (iverbose >= 2 .and. id==master .and. naccreted /= 0) write(iprint,"(a,es10.3,a,i4,a,i4,a)") &
          'Step: at time ',timei,', ',naccreted,' particles were accreted amongst ',nptmass,' sink(s).'

    if (nptmass > 0) then
       call summary_accrete_fail(nfail)
       call summary_accrete(nptmass)
       ! only write to .ev during substeps if no gas particles present
       if (npart==0) call pt_write_sinkev(nptmass,timei,xyzmh_ptmass,vxyz_ptmass, &
                                             fxyz_ptmass,fxyz_ptmass_sinksink)
    endif
    !
    ! check if timestep criterion was violated during substeps
    !
    if (nptmass > 0) then
       if (fonrmax > 0.) then
          dtsinkgas = min(dtsinkgas,C_force*1./sqrt(fonrmax),C_force*sqrt(dtphi2))
       endif
       if (iverbose >= 2) write(iprint,*) nsubsteps,'dt(ext/sink-sink) = ',dtextforcenew,', dt(sink-gas) = ',dtsinkgas
       dtextforcenew = min(dtextforcenew,dtsinkgas)
    endif

    dtextforcenew = reduceall_mpi('min',dtextforcenew)

    dtextforce_min = min(dtextforce_min,dtextforcenew)
    dtextforce = dtextforcenew

    if (last_step) then
       done = .true.
    else
       dt = dtextforce
       if (timei + dt > t_end_step) then
          dt = t_end_step - timei
          last_step = .true.
       endif
    endif
 enddo substeps

 deallocate(dptmass)

 if (nsubsteps > 1) then
    if (iverbose>=1 .and. id==master) then
       write(iprint,"(a,i6,a,f8.2,a,es10.3,a,es10.3)") &
              ' using ',nsubsteps,' substeps (dthydro/dtextf = ',dtsph/dtextforce_min,'), dt = ',dtextforce_min,' dtsph = ',dtsph
    endif
    call summary_variable('ext',iosumextr ,nsubsteps,dtsph/dtextforce_min)
    call summary_variable('ext',iosumextt ,nsubsteps,dtextforce_min,1.0/dtextforce_min)
 endif

end subroutine step_extern_lf


end module step_extern
