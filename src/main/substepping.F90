!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module substepping
!
! Computes sub-steps in the RESPA algorithm
!
!   Multiple option of sub stepping can be choosed depending on
!   the physics and the precision needed
!
!   Only Hydro : substep_sph
!   Hydro + GR : substep_sph_gr substep_gr
!   2nd order with all fast physics implemented  : substep (use_fourthorder = false)
!   4th order without vdep forces and oblateness : substep (not yet implemented)
!
! :References:
!     Verlet (1967), Phys. Rev. 159, 98-103
!     Tuckerman, Berne & Martyna (1992), J. Chem. Phys. 97, 1990-2001
!     Rantala + (2020) (2023),Chin (2007a)
!
! :Owner: Alison Young
!
! :Runtime parameters: None
!
! :Dependencies: chem, cons2primsolver, cooling, cooling_ism, damping, dim,
!   dust_formation, eos, extern_gr, externalforces, io, io_summary,
!   krome_interface, metric_tools, mpiutils, options, part, ptmass,
!   ptmass_radiation, subgroup, timestep, timestep_sts, timing
!
 implicit none


 public :: substep_gr
 public :: substep_sph
 public :: substep_sph_gr
 public :: substep
 public :: get_force
 public :: combine_forces_gr,combine_forces_gr_one

 private

contains

subroutine substep_sph_gr(dt,npart,xyzh,vxyzu,dens,pxyzu,metrics)
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
       if (ierr > 0) call warning('cons2primsolver [in substep_sph_gr (a)]','enthalpy did not converge',i=i)
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
          if (ierr > 0) call warning('cons2primsolver [in substep_sph_gr (b)]','enthalpy did not converge',i=i)
          xyzh(1:3,i) = xpred + 0.5*dt*(vxyzu(1:3,i)-vold)
          diff = maxval(abs(xyzh(1:3,i)-xpred)/xpred)
          if (diff < xtol) converged = .true.
          ! UPDATE METRIC HERE
          call pack_metric(xyzh(1:3,i),metrics(:,:,:,i))
       enddo
       if (niter > nitermax) call warning('substep_sph_gr','Reached max number of x iterations. x_err ',val=diff)

       ! repack values
       eos_vars(igasP,i)  = pri
       eos_vars(itemp,i)  = tempi
       eos_vars(igamma,i) = gammai
    endif
 enddo
 !$omp end parallel do

end subroutine substep_sph_gr

subroutine substep_gr(npart,nptmass,ntypes,dtsph,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,time,&
                       xyzmh_ptmass,vxyz_ptmass,pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass,fxyz_ptmass)
 use dim,            only:maxptmass,maxvxyzu,use_apr
 use io,             only:iverbose,id,master,iprint,warning,fatal
 use part,           only:isdead_or_accreted,iamboundary,igas,iamtype,&
                             massoftype,rhoh,igamma,itemp,igasP
 use io_summary,     only:summary_variable,iosumextr,iosumextt,summary_accrete
 use timestep,       only:bignumber
 use eos,            only:equationofstate
 use cons2primsolver,only:conservative2primitive
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 use damping,        only:calc_damp,apply_damp
 integer, intent(in)    :: npart,ntypes,nptmass
 real,    intent(in)    :: dtsph,time
 real,    intent(inout) :: dtextforce
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:),pxyzu(:,:),dens(:),metrics(:,:,:,:),metricderivs(:,:,:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 real,    intent(inout) :: pxyzu_ptmass(:,:),metrics_ptmass(:,:,:,:),metricderivs_ptmass(:,:,:,:)

 integer :: itype,nsubsteps,naccreted,nlive,nlive_sinks,naccreted_sinks
 real    :: timei,t_end_step,hdt,pmassi
 real    :: dt,dtextforcenew,dtextforce_min
 real    :: accretedmass,damp_fac
 ! real, save :: dmdt = 0.
 logical :: last_step,done
 integer, parameter :: itsmax = 50
 integer :: pitsmax,xitsmax
 real    :: perrmax,xerrmax

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


    call predict_gr(xyzh,vxyzu,ntypes,pxyzu,fext,npart,nptmass,dt,timei,hdt, &
                    dens,metrics,metricderivs,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,&
                    metrics_ptmass,metricderivs_ptmass,pxyzu_ptmass,pitsmax,perrmax, &
                    xitsmax,xerrmax,dtextforcenew)


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
    call accrete_gr(xyzh,vxyzu,dens,fext,metrics,metricderivs,nlive,naccreted,&
                    pxyzu,accretedmass,hdt,npart,nptmass,&
                    ntypes,dtextforce_min,timei,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,&
                    metrics_ptmass,metricderivs_ptmass,&
                    pxyzu_ptmass,nlive_sinks,naccreted_sinks)

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

end subroutine substep_gr

 !----------------------------------------------------------------
 !+
 !  This is the equivalent of the routine below when no external
 !  forces, sink particles or cooling are used
 !+
 !----------------------------------------------------------------
subroutine substep_sph(dt,npart,xyzh,vxyzu)
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

end subroutine substep_sph

!----------------------------------------------------------------
 !+
 !  Substepping of external and sink particle forces.
 !  Also updates position of all particles even if no external
 !  forces applied. This is the internal loop of the RESPA
 !  algorithm over the "fast" forces.
 !  (Here it can be FSI or Leapfrog)
 !+
 !----------------------------------------------------------------
subroutine substep(npart,ntypes,nptmass,dtsph,dtextforce,time,xyzh,vxyzu,fext, &
                   xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dptmass, &
                   sf_ptmass,fsink_old,nbinmax,ibin_wake,gtgrad,group_info, &
                   bin_info,nmatrix,n_group,n_ingroup,n_sing,isionised)
 use io,             only:iverbose,id,master,iprint,fatal
 use options,        only:iexternalforce
 use part,           only:fxyz_ptmass_sinksink,ndptmass
 use io_summary,     only:summary_variable,iosumextr,iosumextt
 use externalforces, only:is_velocity_dependent
 use ptmass,         only:use_fourthorder,use_regnbody,ck,dk,ptmass_check_stars,icreate_sinks
 use subgroup,     only:group_identify
 integer,         intent(in)    :: npart,ntypes
 integer,         intent(inout) :: n_group,n_ingroup,n_sing,nptmass
 integer,         intent(inout) :: group_info(:,:)
 real,            intent(in)    :: dtsph,time
 real,            intent(inout) :: dtextforce
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:)
 real,            intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),dsdt_ptmass(:,:)
 real,            intent(inout) :: dptmass(ndptmass,nptmass),fsink_old(:,:),gtgrad(:,:),bin_info(:,:)
 integer(kind=1), intent(in)    :: nbinmax
 integer        , intent(inout) :: sf_ptmass(:,:)
 integer(kind=1), intent(inout) :: ibin_wake(:),nmatrix(nptmass,nptmass)
 logical,         intent(in)    :: isionised(:)
 logical :: extf_vdep_flag,done,last_step,accreted
 integer :: force_count,nsubsteps
 real    :: timei,time_par,dt,dtgroup,t_end_step
 real    :: dtextforce_min
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
 time_par = time
 extf_vdep_flag = is_velocity_dependent(iexternalforce)
 t_end_step     = timei + dtsph
 nsubsteps      = 0
 dtextforce_min = huge(dt)
 done           = .false.
 accreted       = .false.

 substeps: do while (timei <= t_end_step .and. .not.done)
    force_count = 0
    timei = timei + dt
    if (abs(dt) < tiny(0.)) call fatal('substepping','dt <= 0 in sink-gas substepping',var='dt',val=dt)
    nsubsteps     = nsubsteps + 1

    if (.not.last_step .and. iverbose > 1 .and. id==master) then
       write(iprint,"(a,f14.6)") '> external/ptmass forces only : t=',timei
    endif
!
! Main integration scheme
!
    call kick(dk(1),dt,npart,nptmass,ntypes,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass, &
              fext,fxyz_ptmass,dsdt_ptmass,dptmass)

    call drift(ck(1),dt,time_par,npart,nptmass,ntypes,xyzh,xyzmh_ptmass,&
               vxyzu,vxyz_ptmass,fxyz_ptmass,gtgrad,n_group,n_ingroup,&
               group_info,bin_info)

    call get_force(nptmass,npart,nsubsteps,ntypes,time_par,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass, &
                   vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,dk(2),force_count,extf_vdep_flag,sf_ptmass,&
                   bin_info,group_info,nmatrix,isionised=isionised)

    if (use_fourthorder) then !! FSI 4th order scheme

       ! FSI extrapolation method (Omelyan 2006)
       call get_force(nptmass,npart,nsubsteps,ntypes,time_par,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass, &
                      vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,dk(2),force_count,extf_vdep_flag,sf_ptmass, &
                      bin_info,group_info,nmatrix,fsink_old)

       call kick(dk(2),dt,npart,nptmass,ntypes,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                 fext,fxyz_ptmass,dsdt_ptmass,dptmass)

       call drift(ck(2),dt,time_par,npart,nptmass,ntypes,xyzh,xyzmh_ptmass,&
                  vxyzu,vxyz_ptmass,fxyz_ptmass,gtgrad,n_group,n_ingroup,&
                  group_info,bin_info)

       call get_force(nptmass,npart,nsubsteps,ntypes,time_par,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass, &
                      vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,dk(3),force_count,extf_vdep_flag,sf_ptmass, &
                      bin_info,group_info,nmatrix,isionised=isionised)

       ! the last kick phase of the scheme will perform the accretion loop after velocity update

       call kick(dk(3),dt,npart,nptmass,ntypes,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,fext, &
                 fxyz_ptmass,dsdt_ptmass,dptmass,ibin_wake,nbinmax,timei, &
                 fxyz_ptmass_sinksink,accreted)

       if (use_regnbody) then
          call group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,nmatrix,&
                              dtext=dtgroup)
          call get_force(nptmass,npart,nsubsteps,ntypes,time_par,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass, &
                         vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,dk(3),force_count,extf_vdep_flag,sf_ptmass, &
                         bin_info,group_info,nmatrix)
       elseif (accreted) then
          call get_force(nptmass,npart,nsubsteps,ntypes,time_par,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass, &
                         vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,dk(3),force_count,extf_vdep_flag,sf_ptmass,&
                         bin_info,group_info,nmatrix)
       endif
    else  !! standard leapfrog scheme
       ! the last kick phase of the scheme will perform the accretion loop after velocity update
       call kick(dk(2),dt,npart,nptmass,ntypes,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,fext, &
                fxyz_ptmass,dsdt_ptmass,dptmass,ibin_wake,nbinmax,timei, &
                fxyz_ptmass_sinksink,accreted)
       if (accreted) then
          call get_force(nptmass,npart,nsubsteps,ntypes,time_par,dtextforce,xyzh,vxyzu,fext,xyzmh_ptmass, &
                         vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,dk(2),force_count,extf_vdep_flag,sf_ptmass,&
                         bin_info,group_info,nmatrix)
       endif
    endif

    dtextforce_min = min(dtextforce_min,dtextforce)

    if (last_step) then
       done = .true.
    else
       dt = dtextforce
       dtgroup = dtextforce
       if (timei + dt > t_end_step) then
          dt = t_end_step - timei
          last_step = .true.
       endif
    endif
 enddo substeps

 if (icreate_sinks == 2) call ptmass_check_stars(xyzmh_ptmass,sf_ptmass,nptmass,timei)

 if (nsubsteps > 1) then
    if (iverbose >=1 .and. id==master) then
       write(iprint,"(a,i6,3(a,es10.3))") ' using ',nsubsteps,' substeps '//&
             '(dthydro/dtextf =',dtsph/dtextforce_min,'), dt =',dtextforce_min,' dtsph =',dtsph
    endif
    call summary_variable('ext',iosumextr ,nsubsteps,dtsph/dtextforce_min)
    call summary_variable('ext',iosumextt ,nsubsteps,dtextforce_min,1.0/dtextforce_min)
 endif

end subroutine substep

 !----------------------------------------------------------------
 !+
 !  drift routine for the whole system (part and ptmass)
 !+
 !----------------------------------------------------------------

subroutine drift(cki,dt,time_par,npart,nptmass,ntypes,xyzh,xyzmh_ptmass,vxyzu, &
                 vxyz_ptmass,fxyz_ptmass,gtgrad,n_group,n_ingroup,group_info, &
                 bin_info)
 use part, only: isdead_or_accreted,ispinx,ispiny,ispinz,igarg
 use ptmass,   only:ptmass_drift,use_regnbody
 use subgroup, only:evolve_groups
 use io  ,     only:id,master
 use mpiutils, only:bcast_mpi
 real,    intent(in)    :: dt,cki
 integer, intent(in)    :: npart,nptmass,ntypes
 real,    intent(inout) :: time_par
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:),gtgrad(:,:),bin_info(:,:)
 integer, intent(in)    :: n_ingroup,n_group
 integer, intent(inout) :: group_info(:,:)
 integer :: i
 real    :: ckdt

 ckdt = cki*dt

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
 if (nptmass>0) then
    if (id==master) then
       if (use_regnbody) then
          call ptmass_drift(nptmass,ckdt,xyzmh_ptmass,vxyz_ptmass,group_info,n_ingroup)
       else
          call ptmass_drift(nptmass,ckdt,xyzmh_ptmass,vxyz_ptmass)
       endif
    endif
    call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
 endif

 if (use_regnbody) then
    call evolve_groups(n_group,nptmass,time_par,time_par+cki*dt,group_info,bin_info, &
                       xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad)
 endif

 time_par = time_par + ckdt !! update time for external potential in force routine

end subroutine drift

 !----------------------------------------------------------------
 !+
 !  kick routine for the whole system (part and ptmass)
 !+
 !----------------------------------------------------------------

subroutine kick(dki,dt,npart,nptmass,ntypes,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass, &
                fext,fxyz_ptmass,dsdt_ptmass,dptmass,ibin_wake, &
                nbinmax,timei,fxyz_ptmass_sinksink,accreted)
 use part,           only:isdead_or_accreted,massoftype,iamtype,iamboundary,iphase,ispinx,ispiny,ispinz,igas,ndptmass
 use part,           only:apr_level,aprmassoftype
 use ptmass,         only:f_acc,ptmass_accrete,pt_write_sinkev,update_ptmass,ptmass_kick
 use externalforces, only:accrete_particles
 use options,        only:iexternalforce
 use io  ,           only:id,master,fatal,iprint,iverbose
 use io_summary,     only:summary_accrete,summary_accrete_fail
 use mpiutils,       only:bcast_mpi,reduce_in_place_mpi,reduceall_mpi
 use dim,            only:ind_timesteps,maxp,maxphase,use_apr
 use timestep_sts,   only:sts_it_n
 use timing,         only:get_timings,increment_timer,itimer_acc
 real,                      intent(in)    :: dt,dki
 integer,                   intent(in)    :: npart,nptmass,ntypes
 real,                      intent(inout) :: xyzh(:,:)
 real,                      intent(inout) :: vxyzu(:,:),fext(:,:)
 real,                      intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),dsdt_ptmass(:,:)
 real,                      intent(inout) :: dptmass(ndptmass,nptmass)
 real,            optional, intent(inout) :: fxyz_ptmass_sinksink(:,:)
 real,            optional, intent(in)    :: timei
 integer(kind=1), optional, intent(inout) :: ibin_wake(:)
 integer(kind=1), optional, intent(in)    :: nbinmax
 logical        , optional, intent(inout)   :: accreted
 real(kind=4)    :: t1,t2,tcpu1,tcpu2
 integer(kind=1) :: ibin_wakei
 logical         :: is_accretion
 integer         :: i,itype,nfaili
 integer         :: naccreted,nfail,nlive
 real            :: dkdt,pmassi,fxi,fyi,fzi,accretedmass

 if (present(timei) .and. present(ibin_wake) .and. present(nbinmax)) then
    is_accretion = .true.
 else
    is_accretion = .false.
 endif

 itype = igas
 pmassi = massoftype(igas)

 dkdt = dki*dt

 ! Kick sink particles
 if (nptmass>0) then
    if (id==master) then
       call ptmass_kick(nptmass,dkdt,vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,dsdt_ptmass)
    endif
    call bcast_mpi(vxyz_ptmass(:,1:nptmass))
    call bcast_mpi(xyzmh_ptmass(ispinx,1:nptmass))
    call bcast_mpi(xyzmh_ptmass(ispiny,1:nptmass))
    call bcast_mpi(xyzmh_ptmass(ispinz,1:nptmass))
 endif


 ! Kick gas particles

 if (.not.is_accretion) then
    !$omp parallel do default(none) &
    !$omp shared(maxp,maxphase) &
    !$omp shared(iphase,ntypes) &
    !$omp shared(npart,fext,xyzh,vxyzu,dkdt) &
    !$omp firstprivate(itype) &
    !$omp private(i)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             if (iamboundary(itype)) cycle
          endif
          vxyzu(1,i) = vxyzu(1,i) + dkdt*fext(1,i)
          vxyzu(2,i) = vxyzu(2,i) + dkdt*fext(2,i)
          vxyzu(3,i) = vxyzu(3,i) + dkdt*fext(3,i)
       endif
    enddo
    !$omp end parallel do

 else
    call get_timings(t1,tcpu1)
    accretedmass = 0.
    nfail        = 0
    naccreted    = 0
    nlive        = 0
    ibin_wakei   = 0
    dptmass(:,1:nptmass) = 0.
    !$omp parallel do default(none) &
    !$omp shared(maxp,maxphase) &
    !$omp shared(npart,xyzh,vxyzu,fext,dkdt,iphase,ntypes,massoftype,timei,nptmass,sts_it_n) &
    !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,f_acc,apr_level,aprmassoftype) &
    !$omp shared(iexternalforce) &
    !$omp shared(nbinmax,ibin_wake) &
    !$omp private(i,accreted,nfaili,fxi,fyi,fzi) &
    !$omp firstprivate(itype,pmassi,ibin_wakei) &
    !$omp reduction(+:accretedmass) &
    !$omp reduction(+:nfail) &
    !$omp reduction(+:naccreted) &
    !$omp reduction(+:nlive) &
    !$omp reduction(+:dptmass)
    accreteloop: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             if (iamboundary(itype)) cycle accreteloop
             if (use_apr) then
                pmassi = aprmassoftype(itype,apr_level(i))
             else
                pmassi = massoftype(itype)
             endif
          elseif (use_apr) then
             pmassi = aprmassoftype(igas,apr_level(i))
          endif
          !
          ! correct v to the full step using only the external force
          !
          vxyzu(1,i) = vxyzu(1,i) + dkdt*fext(1,i)
          vxyzu(2,i) = vxyzu(2,i) + dkdt*fext(2,i)
          vxyzu(3,i) = vxyzu(3,i) + dkdt*fext(3,i)

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
                              itype,pmassi,xyzmh_ptmass,vxyz_ptmass,accreted, &
                              dptmass,timei,f_acc,nbinmax,ibin_wakei,nfaili)
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
    !$omp end parallel do

    call get_timings(t2,tcpu2)
    call increment_timer(itimer_acc,t2-t1,tcpu2-tcpu1)

    if (npart > 2 .and. nlive < 2) then
       call fatal('step','all particles accreted',var='nlive',ival=nlive)
    endif

!
! reduction of sink particle changes across MPI
!
    accreted = .false.
    if (nptmass > 0) then
       naccreted = int(reduceall_mpi('+',naccreted))
       nfail = int(reduceall_mpi('+',nfail))
       if (naccreted > 0) then
          accreted = .true.
          call reduce_in_place_mpi('+',dptmass(:,1:nptmass))
          if (id==master) call update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)
       endif
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
 endif


end subroutine kick

!----------------------------------------------------------------
!+
!  force routine for the whole system. First is computed the
!  sink/sink interaction and extf on sink, then comes forces
!  on gas. sink/gas, extf and dampening. Finally there is an
!  update of abundances and temp depending on cooling method
!  during the last force calculation of the substep.
!+
!----------------------------------------------------------------
subroutine get_force(nptmass,npart,nsubsteps,ntypes,timei,dtextforce,xyzh,vxyzu, &
                     fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,dt,dki, &
                     force_count,extf_vdep_flag,sf_ptmass,bin_info,group_info,&
                     nmatrix,fsink_old,isionised)
 use io,              only:iverbose,master,id,iprint,warning,fatal
 use dim,             only:maxp,maxvxyzu,itau_alloc,gr,use_apr,maxptmass
 use ptmass,          only:get_accel_sink_gas,get_accel_sink_sink,merge_sinks, &
                           ptmass_vdependent_correction,n_force_order,use_regnbody,&
                           icreate_sinks
 use options,         only:iexternalforce
 use part,            only:maxphase,abundance,nabundances,epot_sinksink,eos_vars,&
                           isdead_or_accreted,iamboundary,igas,iphase,iamtype,massoftype,divcurlv, &
                           fxyz_ptmass_sinksink,dsdt_ptmass_sinksink,dust_temp,tau,&
                           nucleation,idK2,idmu,idkappa,idgamma,imu,igamma,n_group,n_ingroup,n_sing,&
                           apr_level,aprmassoftype
 use cooling_ism,     only:dphot0,dphotflag,abundsi,abundo,abunde,abundc,nabn
 use timestep,        only:bignumber,C_force
 use mpiutils,        only:bcast_mpi,reduce_in_place_mpi,reduceall_mpi
 use damping,         only:apply_damp,idamp,calc_damp
 use externalforces,  only:update_externalforce
 use ptmass_radiation,only:get_rad_accel_from_ptmass,isink_radiation
 use subgroup,        only:group_identify
 use timing,          only:get_timings,increment_timer,itimer_gasf,itimer_sinksink
 integer,                  intent(in)    :: npart,nsubsteps,ntypes
 integer,                  intent(inout) :: force_count,nptmass
 integer,                  intent(inout) :: sf_ptmass(:,:)
 real,                     intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:)
 real,                     intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,                     intent(inout) :: fxyz_ptmass(4,maxptmass),dsdt_ptmass(3,maxptmass)
 real,                     intent(inout) :: dtextforce
 real,                     intent(in)    :: timei,dki,dt
 logical,                  intent(in)    :: extf_vdep_flag
 real,                     intent(inout) :: bin_info(:,:)
 integer,                  intent(inout) :: group_info(:,:)
 integer(kind=1),          intent(inout) :: nmatrix(:,:)
 real,           optional, intent(inout) :: fsink_old(4,maxptmass)
 logical,        optional, intent(in)    :: isionised(:)
 real(kind=4)    :: t1,t2,tcpu1,tcpu2
 integer         :: merge_ij(nptmass)
 integer         :: merge_n
 integer         :: i,itype
 real, save      :: dmdt = 0.
 real            :: dtf,dtextforcenew,dtsinkgas,dtphi2,fonrmax
 real            :: fextx,fexty,fextz,xi,yi,zi,pmassi,damp_fac
 real            :: fonrmaxi,phii,dtphi2i
 real            :: dkdt,extrapfac
 logical         :: extrap,last

 if (present(fsink_old)) then
    fsink_old = fxyz_ptmass
    extrap  = .true.
 else
    extrap  = .false.
 endif

 force_count   = force_count + 1
 extrapfac     = (1./24.)*dt**2
 dkdt          = dki*dt
 itype         = igas
 pmassi        = massoftype(igas)
 dtextforcenew = bignumber
 dtsinkgas     = bignumber
 dtphi2        = bignumber
 fonrmax       = 0
 last          = (force_count == n_force_order)

 !
 ! update time-dependent external forces
 !
 call calc_damp(timei, damp_fac)
 call update_externalforce(iexternalforce,timei,dmdt)
 !
 ! Sink-sink interactions (loop over ptmass in get_accel_sink_sink)
 !
 if (nptmass > 0) then
    call get_timings(t1,tcpu1)
    if (id==master) then
       if (extrap) then
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                   dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass, &
                                   group_info,bin_info,extrapfac,fsink_old)
          if (merge_n > 0) then
             call merge_sinks(timei,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,sf_ptmass,merge_ij)
             if (use_regnbody) call group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,&
                                                   vxyz_ptmass,group_info,bin_info,nmatrix,dtext=dt)
             call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                      dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass, &
                                      group_info,bin_info,extrapfac,fsink_old)
          endif
       else
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                   dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass, &
                                   group_info,bin_info)
          if (merge_n > 0) then
             call merge_sinks(timei,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,sf_ptmass,merge_ij)
             if (use_regnbody) call group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,&
                                                   vxyz_ptmass,group_info,bin_info,nmatrix,dtext=dt)
             call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                                      dtf,iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass, &
                                      group_info,bin_info)
          endif
       endif
       if (iverbose >= 2) write(iprint,*) 'dt(sink-sink) = ',C_force*dtf
       if (last) then
          fxyz_ptmass_sinksink(:,1:nptmass) = fxyz_ptmass (:,1:nptmass)
          dsdt_ptmass_sinksink(:,1:nptmass) = dsdt_ptmass (:,1:nptmass)
       endif
    else
       fxyz_ptmass(:,1:nptmass) = 0.
       dsdt_ptmass(:,1:nptmass) = 0.
    endif
    call bcast_mpi(epot_sinksink)
    call bcast_mpi(dtf)
    if (icreate_sinks==2) then
       call bcast_mpi(nptmass)
       call bcast_mpi(sf_ptmass)
    endif
    dtextforcenew = min(dtextforcenew,C_force*dtf)
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_sinksink,t2-t1,tcpu2-tcpu1)
 endif

 !
 !-- Forces on gas particles (Sink/gas,extf,damp,cooling,rad pressure)
 !
 call get_timings(t1,tcpu1)

 !$omp parallel default(none) &
 !$omp shared(maxp,maxphase) &
 !$omp shared(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,fext) &
 !$omp shared(eos_vars,dust_temp,idamp,damp_fac,abundance,iphase,ntypes,massoftype) &
 !$omp shared(dkdt,dt,timei,iexternalforce,extf_vdep_flag,last,aprmassoftype,apr_level) &
 !$omp shared(divcurlv,dphotflag,dphot0,nucleation,extrap) &
 !$omp shared(abundc,abundo,abundsi,abunde,extrapfac,fsink_old) &
 !$omp shared(isink_radiation,itau_alloc,tau,isionised) &
 !$omp private(fextx,fexty,fextz,xi,yi,zi) &
 !$omp private(i,fonrmaxi,dtphi2i,phii,dtf) &
 !$omp firstprivate(pmassi,itype) &
 !$omp reduction(min:dtextforcenew,dtphi2) &
 !$omp reduction(max:fonrmax) &
 !$omp reduction(+:fxyz_ptmass,dsdt_ptmass,bin_info)
 !$omp do
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (ntypes > 1 .and. maxphase==maxp) then
          itype  = iamtype(iphase(i))
          if (use_apr) then
             pmassi = aprmassoftype(itype,apr_level(i))
          else
             pmassi = massoftype(itype)
          endif
       endif
       fextx = 0.
       fexty = 0.
       fextz = 0.
       if (extrap) then
          xi = xyzh(1,i) + extrapfac*fext(1,i)
          yi = xyzh(2,i) + extrapfac*fext(2,i)
          zi = xyzh(3,i) + extrapfac*fext(3,i)
       else
          xi = xyzh(1,i)
          yi = xyzh(2,i)
          zi = xyzh(3,i)
       endif
       if (nptmass > 0) then
          if (extrap) then
             call get_accel_sink_gas(nptmass,xi,yi,zi,xyzh(4,i),xyzmh_ptmass,&
                                     fextx,fexty,fextz,phii,pmassi,fxyz_ptmass, &
                                     dsdt_ptmass,fonrmaxi,dtphi2i,bin_info,&
                                     extrapfac,fsink_old)
          else
             call get_accel_sink_gas(nptmass,xi,yi,zi,xyzh(4,i),xyzmh_ptmass,&
                                     fextx,fexty,fextz,phii,pmassi,fxyz_ptmass,&
                                     dsdt_ptmass,fonrmaxi,dtphi2i,bin_info)
             fonrmax = max(fonrmax,fonrmaxi)
             dtphi2  = min(dtphi2,dtphi2i)
          endif
       endif

       !
       ! compute and add external forces
       !
       if (iexternalforce > 0) then
          call get_external_force_gas(xi,yi,zi,xyzh(4,i),vxyzu(1,i), &
                                      vxyzu(2,i),vxyzu(3,i),timei,i, &
                                      dtextforcenew,dtf,dkdt,fextx,fexty, &
                                      fextz,extf_vdep_flag,iexternalforce)
       endif
       !
       ! damping
       !
       if (idamp > 0) then
          call apply_damp(fextx, fexty, fextz, vxyzu(1:3,i), (/xi,yi,zi/), damp_fac)
       endif
       !
       ! Radiation pressure force with isink_radiation
       !
       if (nptmass > 0 .and. isink_radiation > 0) then
          if (extrap) then
             if (itau_alloc == 1) then
                call get_rad_accel_from_ptmass(nptmass,npart,i,xi,yi,zi,xyzmh_ptmass,fextx,fexty,fextz, &
                                               tau=tau,fsink_old=fsink_old,extrapfac=extrapfac)
             else
                call get_rad_accel_from_ptmass(nptmass,npart,i,xi,yi,zi,xyzmh_ptmass,fextx,fexty,fextz, &
                                               fsink_old=fsink_old,extrapfac=extrapfac)
             endif
          else
             if (itau_alloc == 1) then
                call get_rad_accel_from_ptmass(nptmass,npart,i,xi,yi,zi,xyzmh_ptmass,fextx,fexty,fextz,tau)
             else
                call get_rad_accel_from_ptmass(nptmass,npart,i,xi,yi,zi,xyzmh_ptmass,fextx,fexty,fextz)
             endif
          endif
       endif

       fext(1,i) = fextx
       fext(2,i) = fexty
       fext(3,i) = fextz
       !
       ! temperature and abundances update (only done during the last force calculation of the substep)
       !
       if (maxvxyzu >= 4 .and. itype==igas .and. last) then
          call cooling_abundances_update(i,pmassi,xyzh,vxyzu,eos_vars,abundance,nucleation,dust_temp, &
                                         divcurlv,abundc,abunde,abundo,abundsi,dt,dphot0,isionised(i))
       endif
    endif
 enddo
 !$omp enddo
 !$omp end parallel
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_gasf,t2-t1,tcpu2-tcpu1)


 if (nptmass > 0) then
    call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))
    call reduce_in_place_mpi('+',dsdt_ptmass(:,1:nptmass))
    if (id==master .and. extf_vdep_flag) then
       call ptmass_vdependent_correction(nptmass,dkdt,vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,iexternalforce)
    endif
 endif

 if (last) then
    if (nptmass > 0) then
       if (fonrmax > 0.) then
          dtsinkgas = min(dtsinkgas,C_force*1./sqrt(fonrmax),C_force*sqrt(dtphi2))
       endif
       if (iverbose >= 2) write(iprint,*) nsubsteps,'dt(ext/sink-sink) = ',dtextforcenew,', dt(sink-gas) = ',dtsinkgas
       dtextforcenew = min(dtextforcenew,dtsinkgas)
    endif

    dtextforcenew = reduceall_mpi('min',dtextforcenew)
    dtextforce = dtextforcenew
 endif

end subroutine get_force

!-----------------------------------------------------------------------------------
!+
! Update of abundances and internal energy using cooling method (see cooling module)
! NOTE: The chemistry and cooling here is implicitly calculated.  That is,
!       dt is *passed in* to the chemistry & cooling routines so that the
!       output will be at the correct time of time + dt.  Since this is
!       implicit, there is no cooling timestep.  Explicit cooling is
!       calculated in force and requires a cooling timestep.
!+
!------------------------------------------------------------------------------------
subroutine cooling_abundances_update(i,pmassi,xyzh,vxyzu,eos_vars,abundance,nucleation,dust_temp, &
                                     divcurlv,abundc,abunde,abundo,abundsi,dt,dphot0,isionisedi)
 use dim,             only:h2chemistry,do_nucleation,use_krome,update_muGamma,store_dust_temperature
 use part,            only:idK2,idmu,idkappa,idgamma,imu,igamma,nabundances
 use cooling_ism,     only:nabn,dphotflag
 use options,         only:icooling
 use chem,            only:update_abundances,get_dphot
 use dust_formation,  only:evolve_dust,calc_muGamma
 use cooling,         only:energ_cooling,cooling_in_step
 use part,            only:rhoh
#ifdef KROME
 use part,            only: T_gas_cool
 use krome_interface, only: update_krome
 real                       :: ui
#endif
 real,         intent(inout) :: vxyzu(:,:),xyzh(:,:)
 real,         intent(inout) :: eos_vars(:,:),abundance(:,:)
 real,         intent(inout) :: nucleation(:,:),dust_temp(:)
 real(kind=4), intent(in)    :: divcurlv(:,:)
 real,         intent(inout) :: abundc,abunde,abundo,abundsi
 real(kind=8), intent(in)    :: dphot0
 real,         intent(in)    :: dt,pmassi
 logical,      intent(in)    :: isionisedi
 integer,      intent(in)    :: i

 real :: dudtcool,rhoi,dphot,pH,pH_tot
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
 elseif (update_muGamma) then
    call calc_muGamma(rhoi, dust_temp(i),eos_vars(imu,i),eos_vars(igamma,i), pH, pH_tot)
 endif
 !
 ! COOLING
 !
 if (icooling > 0 .and. cooling_in_step .and. icooling/=9) then
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
 if (isionisedi .or. icooling == 9) dudtcool = 0.
 if (cooling_in_step .or. use_krome) vxyzu(4,i) = vxyzu(4,i) + dt * dudtcool


end subroutine cooling_abundances_update

 !----------------------------------------------------------------
 !+
 !  routine for external force applied on gas particle
 !+
 !----------------------------------------------------------------

subroutine get_external_force_gas(xi,yi,zi,hi,vxi,vyi,vzi,timei,i,dtextforcenew,dtf,dkdt, &
                                 fextx,fexty,fextz,extf_is_velocity_dependent,iexternalforce)
 use timestep,       only:C_force
 use externalforces, only: externalforce,update_vdependent_extforce
 real,    intent(in) :: xi,yi,zi,hi,vxi,vyi,vzi,timei,dkdt
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
    call update_vdependent_extforce(iexternalforce,vxi,vyi,vzi, &
                                             fextxi,fextyi,fextzi,fextv,dkdt,xi,yi,zi)
    fextx = fextx + fextv(1)
    fexty = fexty + fextv(2)
    fextz = fextz + fextv(3)
 endif

end subroutine get_external_force_gas

 !----------------------------------------------------------------
 !+
 !  routine for prediction substep in GR case
 !+
 !----------------------------------------------------------------
subroutine predict_gr(xyzh,vxyzu,ntypes,pxyzu,fext,npart,nptmass,dt,timei,hdt, &
                     dens,metrics,metricderivs,&
                     xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,&
                     metrics_ptmass,metricderivs_ptmass,pxyzu_ptmass,&
                     pitsmax,perrmax,&
                     xitsmax,xerrmax,dtextforcenew)

 use dim,            only:maxptmass,maxp,maxvxyzu,use_apr
 use io,             only:master,warning,fatal
 use externalforces, only:externalforce,accrete_particles,update_externalforce
 use part,           only:maxphase,isdead_or_accreted,iamboundary,igas,iphase,iamtype,&
                          massoftype,rhoh,ien_type,eos_vars,igamma,itemp,igasP,&
                          aprmassoftype,apr_level,epot_sinksink,fxyz_ptmass_sinksink,dsdt_ptmass,&
                          fsink_old
 use io_summary,     only:summary_variable,iosumextr,iosumextt,summary_accrete
 use timestep,       only:bignumber,xtol,ptol
 use eos,            only:equationofstate,ieos
 use cons2primsolver,only:conservative2primitive
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 use damping,        only:calc_damp,apply_damp
 use options,        only:iexternalforce
 use ptmass,         only:get_accel_sink_sink,get_accel_sink_gas


 real,    intent(inout)   :: xyzh(:,:),vxyzu(:,:),fext(:,:),pxyzu(:,:),dens(:),metrics(:,:,:,:),metricderivs(:,:,:,:)
 real,    intent(inout)   :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),pxyzu_ptmass(:,:)
 real,    intent(inout)   :: metrics_ptmass(:,:,:,:),metricderivs_ptmass(:,:,:,:)
 real,    intent(in)      :: dt,hdt,timei
 integer, intent(in)      :: npart,ntypes,nptmass
 integer, intent(inout)   :: pitsmax,xitsmax
 real,    intent(inout)   :: perrmax,xerrmax,dtextforcenew

 integer :: i,its,ierr,itype
 integer, parameter :: itsmax = 50
 real    :: pmassi
 real    :: pri,spsoundi,tempi,gammai
 real, save :: pprev(3),xyz_prev(3),fstar(3),vxyz_star(3),xyz(3),pxyz(3),vxyz(3),fexti(3)
 !$omp threadprivate(pprev,xyz_prev,fstar,vxyz_star,xyz,pxyz,vxyz,fexti)
 real    :: x_err,pmom_err
 ! real, save :: dmdt = 0.
 logical :: converged
 real :: rhoi,hi,eni,uui,densi,poti
 real :: bin_info(6,nptmass)
 real :: dtphi2,dtsinksink,fonrmax
 integer :: merge_ij(nptmass),merge_n

 pmassi  = massoftype(igas)
 itype   = igas
 fsink_old = fxyz_ptmass
 fxyz_ptmass = 0.
 !----------------------------------------------
 ! calculate acceleration sink-sink
 !----------------------------------------------
 call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass_sinksink,epot_sinksink,dtsinksink,&
                         iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass)
 !----------------------------------------------
 ! predictor during substeps for gas particles
 !----------------------------------------------
 !
 ! predictor step for external forces, also recompute external forces
 !
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,fext,iphase,ntypes,massoftype) &
 !$omp shared(nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass) &
 !$omp shared(maxphase,maxp,eos_vars) &
 !$omp shared(bin_info,dtphi2,poti,fonrmax) &
 !$omp shared(fxyz_ptmass_sinksink) &
 !$omp shared(dt,hdt,xtol,ptol,aprmassoftype,apr_level) &
 !$omp shared(ieos,pxyzu,dens,metrics,metricderivs,ien_type) &
 !$omp shared(pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass) &
 !$omp shared(dtsinksink,epot_sinksink,merge_ij,merge_n,dsdt_ptmass,iexternalforce) &
 !$omp private(i,its,spsoundi,tempi,rhoi,hi,eni,uui,densi) &
 !$omp private(converged,pmom_err,x_err,pri,ierr,gammai) &
 !$omp firstprivate(pmassi,itype) &
 !$omp reduction(max:xitsmax,pitsmax,perrmax,xerrmax) &
 !$omp reduction(min:dtextforcenew)
 predictor: do i=1,npart
    xyz(1) = xyzh(1,i)
    xyz(2) = xyzh(2,i)
    xyz(3) = xyzh(3,i)
    hi     = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (ntypes > 1 .and. maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (use_apr) then
             pmassi = aprmassoftype(itype,apr_level(i))
          else
             pmassi = massoftype(itype)
          endif
       elseif (use_apr) then
          pmassi = aprmassoftype(igas,apr_level(i))
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

          if (ierr > 0) call warning('cons2primsolver [in substep_gr (a)]','enthalpy did not converge',i=i)
          call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i),vxyz,densi,uui,pri,fstar)

          ! calculate force between sink-gas particles
          call get_accel_sink_gas(nptmass,xyz(1),xyz(2),xyz(3),hi,xyzmh_ptmass, &
                                  fstar(1),fstar(2),fstar(3),poti,pmassi,fxyz_ptmass,&
                                  dsdt_ptmass,fonrmax,dtphi2,bin_info)

          pxyz = pprev + hdt*(fstar - fexti)
          pmom_err = maxval(abs(pxyz - pprev))
          if (pmom_err < ptol) converged = .true.
          fexti = fstar
       enddo pmom_iterations
       if (its > itsmax ) call warning('substep_gr',&
                              'max # of pmom iterations',var='pmom_err',val=pmom_err)
       pitsmax = max(its,pitsmax)
       perrmax = max(pmom_err,perrmax)

       call conservative2primitive(xyz,metrics(:,:,:,i),vxyz,densi,uui,pri,tempi,&
                                    gammai,rhoi,pxyz,eni,ierr,ien_type)
       if (ierr > 0) call warning('cons2primsolver [in substep_gr (b)]','enthalpy did not converge',i=i)
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
          if (ierr > 0) call warning('cons2primsolver [in substep_gr (c)]','enthalpy did not converge',i=i)
          xyz  = xyz_prev + hdt*(vxyz_star - vxyz)
          x_err = maxval(abs(xyz-xyz_prev))
          if (x_err < xtol) converged = .true.
          vxyz = vxyz_star
          ! UPDATE METRIC HERE
          call pack_metric(xyz,metrics(:,:,:,i))
       enddo xyz_iterations
       call pack_metricderivs(xyz,metricderivs(:,:,:,i))
       if (its > itsmax ) call warning('substep_gr','Reached max number of x iterations. x_err ',val=x_err)
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

 call predict_gr_sink(xyzmh_ptmass,vxyz_ptmass,ntypes,pxyzu_ptmass,fsink_old,fxyz_ptmass,fxyz_ptmass_sinksink,nptmass,&
                      dt,timei,hdt,metrics_ptmass,metricderivs_ptmass,dtextforcenew,pitsmax,perrmax,&
                      xitsmax,xerrmax)

end subroutine predict_gr

 !----------------------------------------------------------------
 !+
 !  routine for prediction substep in GR case
 !+
 !----------------------------------------------------------------
subroutine predict_gr_sink(xyzmh_ptmass,vxyz_ptmass,ntypes,pxyzu_ptmass,fsink_old,fxyz_ptmass,&
                     fxyz_ptmass_sinksink,nptmass,dt,timei,hdt, &
                     metrics_ptmass,metricderivs_ptmass,dtextforcenew,pitsmax,perrmax, &
                     xitsmax,xerrmax)
 use dim,            only:maxptmass
 use io,             only:master,warning,fatal
 use part,           only:epot_sinksink,dsdt_ptmass
 use timestep,       only:bignumber,xtol,ptol
 use cons2primsolver,only:conservative2primitive
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 use ptmass,         only:get_accel_sink_sink

 real,    intent(inout)   :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),fxyz_ptmass_sinksink(:,:),pxyzu_ptmass(:,:)
 real,    intent(inout)   :: metrics_ptmass(:,:,:,:),metricderivs_ptmass(:,:,:,:),dtextforcenew,fsink_old(:,:)
 real,    intent(in)      :: dt,hdt,timei
 integer, intent(in)      :: nptmass,ntypes
 integer, intent(inout)   :: pitsmax,xitsmax
 real,    intent(inout)   :: perrmax,xerrmax

 integer :: i,its,ierr
 integer, parameter :: itsmax = 50
 real    :: pmassi,xyzhi(4)
 real    :: pri,tempi,gammai
 real, save :: pprev(3),xyz_prev(3),fstar(3),vxyz_star(3),xyz(3),pxyz(3),vxyz(3),fexti(3)
 !$omp threadprivate(pprev,xyz_prev,fstar,vxyz_star,xyz,pxyz,vxyz,fexti)
 real    :: x_err,pmom_err
 ! real, save :: dmdt = 0.
 logical :: converged
 real :: rhoi,hi,eni,uui,densi
 integer :: merge_ij(2),merge_n
 real    :: dtsinksink

 !---------------------------
 ! predictor during substeps
 !---------------------------
 !
 ! predictor step for external forces, also recompute external forces
 !
 !$omp parallel do default(none) &
 !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink) &
 !$omp shared(dt,hdt,xtol,ptol,nptmass) &
 !$omp shared(pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass,fsink_old) &
 !$omp shared(dtsinksink,epot_sinksink,merge_ij,merge_n,dsdt_ptmass) &
 !$omp private(i,its,tempi,rhoi,hi,eni,uui,densi,xyzhi) &
 !$omp private(converged,pmom_err,x_err,pri,ierr,gammai,pmassi) &
 !$omp reduction(max:xitsmax,pitsmax,perrmax,xerrmax) &
 !$omp reduction(min:dtextforcenew)

 predictor: do i=1,nptmass
    xyzhi(1) = xyzmh_ptmass(1,i)
    xyzhi(2) = xyzmh_ptmass(2,i)
    xyzhi(3) = xyzmh_ptmass(3,i)
    pmassi   = xyzmh_ptmass(4,i)
    hi       = xyzmh_ptmass(5,i)

    xyz(1)   = xyzhi(1)
    xyz(2)   = xyzhi(2)
    xyz(3)   = xyzhi(3)
    xyzhi(4) = hi
    if (pmassi >= 0) then
       its       = 0
       converged = .false.
       !
       ! make local copies of array quantities
       !
       pxyz(1:3) = pxyzu_ptmass(1:3,i)
       eni       = 0.
       vxyz(1:3) = vxyz_ptmass(1:3,i)
       uui       = 0.
       fexti     = fsink_old(1:3,i)
       pxyz      = pxyz + hdt*fexti

       !-- unpack thermo variables for the first guess in cons2prim
       densi     = 1.
       pri       = 0.
       gammai    = 0.
       tempi     = 0.
       rhoi      = 1.
       ! Note: grforce needs derivatives of the metric,
       ! which do not change between pmom iterations
       pmom_iterations: do while (its <= itsmax .and. .not. converged)
          its   = its + 1
          pprev = pxyz
          call conservative2primitive(xyz,metrics_ptmass(:,:,:,i),vxyz,densi,uui,pri,&
                                     tempi,gammai,rhoi,pxyz,eni,ierr,1)

          if (ierr > 0) call warning('cons2primsolver [in substep_gr (a)]','enthalpy did not converge',i=i)

          call get_grforce(xyzhi,metrics_ptmass(:,:,:,i),metricderivs_ptmass(:,:,:,i),vxyz,densi,uui,pri,fstar)

          ! add forces from curvature on sink, from sink sink interaction and sink-gas interaction
          fstar = fstar + fxyz_ptmass_sinksink(1:3,i) + fxyz_ptmass(1:3,i)

          pxyz = pprev + hdt*(fstar - fexti)
          pmom_err = maxval(abs(pxyz - pprev))
          if (pmom_err < ptol) converged = .true.
          fexti = fstar
       enddo pmom_iterations
       if (its > itsmax ) call warning('substep_gr',&
                              'max # of pmom iterations',var='pmom_err',val=pmom_err)
       pitsmax = max(its,pitsmax)
       perrmax = max(pmom_err,perrmax)

       call conservative2primitive(xyz,metrics_ptmass(:,:,:,i),vxyz,densi,uui,pri,tempi,&
                                    gammai,rhoi,pxyz,eni,ierr,1)

       if (ierr > 0) call warning('cons2primsolver [in substep_gr (b)]','enthalpy did not converge',i=i)
       xyz = xyz + dt*vxyz
       call pack_metric(xyz,metrics_ptmass(:,:,:,i))

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
          call conservative2primitive(xyz,metrics_ptmass(:,:,:,i),vxyz_star,densi,uui,&
                                       pri,tempi,gammai,rhoi,pxyz,eni,ierr,1)
          if (ierr > 0) call warning('cons2primsolver [in substep_gr (c)]','enthalpy did not converge',i=i)
          xyz  = xyz_prev + hdt*(vxyz_star - vxyz)
          x_err = maxval(abs(xyz-xyz_prev))
          if (x_err < xtol) converged = .true.
          vxyz = vxyz_star
          ! UPDATE METRIC HERE
          call pack_metric(xyz,metrics_ptmass(:,:,:,i))
       enddo xyz_iterations
       call pack_metricderivs(xyz,metricderivs_ptmass(:,:,:,i))
       if (its > itsmax ) call warning('substep_gr','Reached max number of x iterations. x_err ',val=x_err)
       xitsmax = max(its,xitsmax)
       xerrmax = max(x_err,xerrmax)

       ! re-pack arrays back where they belong
       xyzmh_ptmass(1:3,i) = xyz(1:3)
       pxyzu_ptmass(1:3,i) = pxyz(1:3)
       vxyz_ptmass(1:3,i)  = vxyz(1:3)
       fxyz_ptmass(1:3,i)  = fexti

    endif
 enddo predictor
!$omp end parallel do

end subroutine predict_gr_sink

 !----------------------------------------------------------------
 !+
 !  routine for accretion step in GR case
 !+
 !----------------------------------------------------------------
subroutine accrete_gr(xyzh,vxyzu,dens,fext,metrics,metricderivs,nlive,naccreted,&
                      pxyzu,accretedmass,hdt,npart,nptmass,ntypes,dtextforce_min,timei,&
                      xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,&
                      metrics_ptmass,metricderivs_ptmass,&
                      pxyzu_ptmass,nlive_sinks,naccreted_sinks)

 use dim,            only:maxptmass,maxp,maxvxyzu,use_apr
 use io,             only:master,warning,fatal
 use externalforces, only:externalforce,accrete_particles,update_externalforce
 use options,        only:iexternalforce
 use part,           only:maxphase,isdead_or_accreted,iamboundary,igas,iphase,iamtype,&
                          massoftype,rhoh,igamma,itemp,igasP,aprmassoftype,apr_level,&
                          fxyz_ptmass_sinksink,fsink_old
 use io_summary,     only:summary_variable,iosumextr,iosumextt,summary_accrete
 use timestep,       only:bignumber,C_force
 use eos,            only:equationofstate,ieos
 use cons2primsolver,only:conservative2primitive
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 use damping,        only:calc_damp,apply_damp,idamp
 use part,           only:epot_sinksink
 use ptmass,         only:get_accel_sink_sink,get_accel_sink_gas

 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:),pxyzu(:,:),dens(:),metrics(:,:,:,:),metricderivs(:,:,:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),pxyzu_ptmass(:,:)
 real,    intent(inout) :: metrics_ptmass(:,:,:,:),metricderivs_ptmass(:,:,:,:)
 integer, intent(in)    :: npart,ntypes,nptmass
 integer, intent(inout) :: nlive,naccreted
 integer, intent(inout) :: nlive_sinks,naccreted_sinks
 real,    intent(inout) :: accretedmass
 real,    intent(in)    :: hdt,timei
 real,    intent(inout) :: dtextforce_min

 logical :: accreted
 integer :: i,itype
 real    :: pmassi
 real    :: dtf
 real    :: pri,spsoundi,pondensi,tempi
 real    :: damp_fac
 ! real, save :: dmdt = 0.
 integer, parameter :: itsmax = 50
 real :: bin_info(6,nptmass),dsdt_ptmass(3,nptmass)
 real :: dtphi2,dtsinksink,fonrmax,poti
 integer :: merge_ij(nptmass),merge_n

 pmassi     = massoftype(igas)
 itype      = igas
 fsink_old = fxyz_ptmass
 fxyz_ptmass = 0.
 !----------------------------------------------
 ! calculate acceleration sink-sink
 !----------------------------------------------
 call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass_sinksink,epot_sinksink,dtsinksink,&
                         iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass)

 dtextforce_min = min(dtextforce_min,C_force*dtsinksink)

 !$omp parallel default(none) &
 !$omp shared(npart,xyzh,metrics,metricderivs,vxyzu,fext,iphase,ntypes,massoftype,hdt,timei) &
 !$omp shared(nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass) &
 !$omp shared(pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass) &
 !$omp shared(maxphase,maxp,aprmassoftype,apr_level) &
 !$omp shared(dtsinksink,epot_sinksink,merge_ij,merge_n,dsdt_ptmass) &
 !$omp shared(bin_info,dtphi2,poti,fonrmax) &
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
          if (use_apr) then
          else
             pmassi = massoftype(itype)
          endif
          !  if (itype==iboundary) cycle accreteloop
       elseif (use_apr) then
          pmassi = aprmassoftype(igas,apr_level(i))
       endif

       call equationofstate(ieos,pondensi,spsoundi,dens(i),xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
       pri = pondensi*dens(i)

       call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i),pri,fext(1:3,i),dtf)
       call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass, &
                                  fext(1,i),fext(2,i),fext(3,i),poti,pmassi,fxyz_ptmass,&
                                  dsdt_ptmass,fonrmax,dtphi2,bin_info)

       dtextforce_min = min(dtextforce_min,C_force*dtf,C_force*sqrt(dtphi2))

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
 call accrete_gr_sink(xyzmh_ptmass,vxyz_ptmass,fsink_old,fxyz_ptmass,fxyz_ptmass_sinksink,&
                      metrics_ptmass,metricderivs_ptmass,nlive_sinks,naccreted_sinks,pxyzu_ptmass,&
                      accretedmass,hdt,nptmass,dtextforce_min,timei,dtsinksink)
end subroutine accrete_gr

 !----------------------------------------------------------------
 !+
 !  routine for accretion step in GR case
 !+
 !----------------------------------------------------------------
subroutine accrete_gr_sink(xyzmh_ptmass,vxyz_ptmass,fsink_old,fxyz_ptmass,fxyz_ptmass_sinksink,metrics_ptmass,metricderivs_ptmass,&
                            nlive_sinks,naccreted_sinks,&
                            pxyzu_ptmass,accretedmass,hdt,nptmass,dtextforce_min,timei,dtsinksink)
 use part,           only:ihsoft
 use externalforces, only:accrete_particles
 use options,        only:iexternalforce
 use timestep,       only:bignumber,C_force
 use cons2primsolver,only:conservative2primitive
 use extern_gr,      only:get_grforce
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),pxyzu_ptmass(:,:),fsink_old(:,:)
 real,    intent(inout) :: metrics_ptmass(:,:,:,:),metricderivs_ptmass(:,:,:,:),fxyz_ptmass_sinksink(:,:)
 integer, intent(in)    :: nptmass
 integer, intent(inout) :: nlive_sinks,naccreted_sinks
 real,    intent(inout) :: accretedmass
 real,    intent(in)    :: hdt,timei,dtsinksink
 real,    intent(inout) :: dtextforce_min

 logical :: accreted
 integer :: i
 real    :: xyzhi(4),pmassi,densi,pri
 real    :: dtf,hsofti
 ! real, save :: dmdt = 0.
 integer, parameter :: itsmax = 50

 !$omp parallel default(none) &
 !$omp shared(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass,vxyz_ptmass,fxyz_ptmass,hdt,timei) &
 !$omp shared(dtsinksink,fxyz_ptmass_sinksink) &
 !$omp private(i,accreted) &
 !$omp shared(pxyzu_ptmass,iexternalforce,C_force) &
 !$omp private(dtf,xyzhi,hsofti,pmassi,pri,densi) &
 !$omp reduction(min:dtextforce_min) &
 !$omp reduction(+:accretedmass,naccreted_sinks,nlive_sinks)
 !$omp do
 accreteloop: do i=1,nptmass
    pri = 0.
    densi = 1.

    ! add this force due to the curvature of the metric.
    xyzhi(1:3) = xyzmh_ptmass(1:3,i)

    ! if a sink particle is already eaten by the black hole, skip it...
    pmassi = xyzmh_ptmass(4,i)
    if (pmassi < 0.) cycle accreteloop
    !
    ! the smoothing length is used inside get_grforce to set the timestep based
    ! on h/abs(dp/dt), but for sink particles this is meaningless unless
    ! a softening length is set
    !
    hsofti = xyzmh_ptmass(ihsoft,i)
    xyzhi(4) = xyzmh_ptmass(5,i)
    if (hsofti > 0.) xyzhi(4) = hsofti
    ! add force from sink-gas interaction to the sink-sink interaction array
    fxyz_ptmass_sinksink(:,i) = fxyz_ptmass_sinksink(:,i) + fxyz_ptmass(:,i)
    ! calculate force due to curvature
    call get_grforce(xyzhi,metrics_ptmass(:,:,:,i),metricderivs_ptmass(:,:,:,i),vxyz_ptmass(1:3,i),&
                     densi,0.,pri,fxyz_ptmass(1:3,i),dtf)
    ! get total force on sinks
    fxyz_ptmass(:,i) = fxyz_ptmass(:,i) + fxyz_ptmass_sinksink(:,i)

    dtextforce_min = min(dtextforce_min,C_force*dtf)
    !
    ! correct v to the full step using only the external force
    !
    pxyzu_ptmass(1:3,i) = pxyzu_ptmass(1:3,i) + hdt*fxyz_ptmass(1:3,i)
    if (pmassi < 0.) then
       !
       ! sending the mass twice here is deliberate, as an accreted sink particle is indicated by
       ! a negative mass, unlike gas particles which are flagged with a negative smoothing length
       !
       call accrete_particles(iexternalforce,xyzmh_ptmass(1,i),xyzmh_ptmass(2,i), &
                                 xyzmh_ptmass(3,i),xyzmh_ptmass(4,i),xyzmh_ptmass(4,i),timei,accreted)
       if (accreted) then
          accretedmass = accretedmass + abs(xyzmh_ptmass(4,i))
          naccreted_sinks = naccreted_sinks + 1
       endif
    endif
    nlive_sinks = nlive_sinks + 1

 enddo accreteloop
 !$omp enddo
 !$omp end parallel

end subroutine accrete_gr_sink

subroutine combine_forces_gr(nptmass,fsinks,fgr)
 real, intent(in)    :: fsinks(:,:)
 integer, intent(in) :: nptmass

 real, intent(inout) :: fgr(:,:)

 integer :: i

 do i=1,nptmass
    fgr(:,i) = fsinks(:,i) + fgr(:,i)
 enddo
end subroutine combine_forces_gr


subroutine combine_forces_gr_one(fsink,fgr)
 real, intent(in)    :: fsink(:)
 real, intent(inout) :: fgr(:)

 fgr(:) = fgr(:) + fsink(:)

end subroutine combine_forces_gr_one

end module substepping
