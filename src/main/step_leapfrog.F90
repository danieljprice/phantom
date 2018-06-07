!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: step_lf_global
!
!  DESCRIPTION:
!   Computes one (hydro) timestep
!
!   Change this subroutine to change the timestepping algorithm
!
!   This version uses a Velocity Verlet (leapfrog) integrator with
!   substepping (operator splitting) of external/sink particle forces,
!   following the Reversible RESPA algorithm of Tuckerman et al. 1992
!
!  REFERENCES:
!     Verlet (1967), Phys. Rev. 159, 98-103
!     Tuckerman, Berne & Martyna (1992), J. Chem. Phys. 97, 1990-2001
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: chem, coolfunc, deriv, dim, eos, externalforces, growth,
!    io, io_summary, mpiutils, options, part, ptmass, timestep,
!    timestep_ind, timestep_sts
!+
!--------------------------------------------------------------------------
module step_lf_global
 use dim,  only:maxp,maxvxyzu,maxBevol,ndusttypes
 use part, only:vpred,Bpred,dustpred,ppred
 use timestep_ind, only:maxbins,itdt,ithdt,itdt1,ittwas
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"
 real               :: ibin_dts(4,0:maxbins)

contains

!------------------------------------------------------------
!+
!  initialisation routine necessary for individual timesteps
!+
!------------------------------------------------------------
subroutine init_step(npart,time,dtmax)
#ifdef IND_TIMESTEPS
 use timestep_ind, only:get_dt
 use part,         only:ibin,twas
#endif
 integer, intent(in) :: npart
 real,    intent(in) :: time,dtmax
#ifdef IND_TIMESTEPS
 integer             :: i
!
! twas is set so that at start of step we predict
! forwards to half of current timestep
!
 !$omp parallel do schedule(static) private(i)
 do i=1,npart
    twas(i) = time + 0.5*get_dt(dtmax,ibin(i))
 enddo
 !
 ! For each ibin option, calculate dt, dt/2, 1/dt and twas
 do i=0,maxbins
    ibin_dts(itdt,  i) = get_dt(dtmax,int(i,kind=1))
    ibin_dts(ithdt, i) = 0.5*ibin_dts(itdt,i)
    ibin_dts(itdt1, i) = 1.0/ibin_dts(itdt,i)
    ibin_dts(ittwas,i) = time + 0.5*get_dt(dtmax,int(i,kind=1))
 enddo
#endif

end subroutine init_step

!------------------------------------------------------------
!+
!  main timestepping routine
!+
!------------------------------------------------------------
subroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)
 use dim,            only:maxp,ndivcurlv,maxvxyzu,maxptmass,maxalpha,nalpha,h2chemistry,use_dustgrowth,gr
 use io,             only:iprint,fatal,iverbose,id,master,warning
 use options,        only:damp,tolv,iexternalforce,icooling,use_dustfrac
 use part,           only:xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol, &
                          isdead_or_accreted,rhoh,dhdrho,&
                          iphase,iamtype,massoftype,maxphase,igas,idust,mhd,maxBevol,&
                          switches_done_in_derivs,iboundary,get_ntypes,npartoftype,&
                          dustfrac,dustevol,ddustfrac,temperature,alphaind,nptmass,store_temperature,&
                          dustprop,ddustprop,dustproppred,pxyzu,dens,metrics,metricderivs
 use eos,            only:get_spsound
 use options,        only:avdecayconst,alpha,ieos,alphamax
 use deriv,          only:derivs
 use timestep,       only:dterr,bignumber
 use mpiutils,       only:reduceall_mpi
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,ibin_wake
 use io_summary,     only:summary_printout,summary_variable,iosumtvi,iowake
 use coolfunc,       only:energ_coolfunc
#ifdef IND_TIMESTEPS
 use timestep,       only:dtmax,dtmax_rat,dtdiff,mod_dtmax_in_step
 use timestep_ind,   only:get_dt,nbinmax,decrease_dtmax,ibinnow
 use timestep_sts,   only:sts_get_dtau_next,use_sts,ibin_sts,sts_it_n
 use part,           only:ibin,ibin_old,twas,iactive
#endif
#ifdef GR
 use metric,         only:imetric
 use metric_tools,   only:imet_minkowski
 use cons2prim,      only:conservative_to_primitive
#endif
#ifdef DUSTGROWTH
 use growth,         only:update_dustprop
#endif
 integer, intent(inout) :: npart
 integer, intent(in)    :: nactive
 real,    intent(in)    :: t,dtsph
 real,    intent(inout) :: dtextforce
 real,    intent(out)   :: dtnew
 integer            :: i,its,np,ntypes,itype,nwake
 real               :: timei,erri,errmax,v2i,errmaxmean
 real               :: vxi,vyi,vzi,eni,vxoldi,vyoldi,vzoldi,hdtsph,pmassi
 real               :: alphaloci,divvdti,source,tdecay1,hi,rhoi,ddenom,spsoundi
 real               :: v2mean,hdti
 real               :: pxi,pyi,pzi,p2i,p2mean
#ifdef IND_TIMESTEPS
 real               :: dtsph_next,dti,time_now
#ifdef MPI
 logical, parameter :: allow_waking = .false.
#else
 logical, parameter :: allow_waking = .true.
#endif
#else
 integer(kind=1), parameter :: nbinmax = 0
#endif
 integer, parameter :: maxits = 30
 logical            :: converged,store_itype
!
! set initial quantities
!
 timei  = t
 hdtsph = 0.5*dtsph
 dterr  = bignumber

! determine twas for each ibin
#ifdef IND_TIMESTEPS
 if (sts_it_n) then
    time_now = timei + dtsph
    do i=0,maxbins
       ibin_dts(ittwas,i) = (int(time_now*ibin_dts(itdt1,i),kind=8) + 0.5)*ibin_dts(itdt,i)
    enddo
 endif
#endif

!--------------------------------------
! velocity predictor step, using dtsph
!--------------------------------------
 itype   = igas
 ntypes  = get_ntypes(npartoftype)
 pmassi  = massoftype(itype)
 store_itype = (maxphase==maxp .and. ntypes > 1)

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,fxyzu,iphase,hdtsph,store_itype) &
#ifdef GR
 !$omp shared(pxyzu) &
#endif
 !$omp shared(Bevol,dBevol,dustevol,ddustfrac,use_dustfrac) &
 !$omp shared(dustprop,ddustprop,dustproppred) &
#ifdef IND_TIMESTEPS
 !$omp shared(ibin,ibin_old,twas,timei) &
#endif
 !$omp firstprivate(itype) &
 !$omp private(i,hdti)
 predictor: do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
#ifdef IND_TIMESTEPS
       if (iactive(iphase(i))) ibin_old(i) = ibin(i) ! only required for ibin_neigh in force.F90
       !
       !--synchronise all particles to their half timesteps
       !
       hdti = twas(i) - timei
#else
       hdti = hdtsph
#endif
       if (store_itype) itype = iamtype(iphase(i))
       if (itype==iboundary) cycle predictor
       !
       ! predict v and u to the half step with "slow" forces
       !
#ifdef GR
       pxyzu(:,i) = pxyzu(:,i) + hdti*fxyzu(:,i)
#else
       vxyzu(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
#endif
       if (itype==idust .and. use_dustgrowth) then
          dustproppred(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)
       endif
       if (itype==igas) then
          if (mhd)          Bevol(:,i)    = Bevol(:,i)        + hdti*dBevol(:,i)
          if (use_dustfrac) dustevol(:,i) = abs(dustevol(:,i) + hdti*ddustfrac(:,i))
       endif
    endif
 enddo predictor
 !omp end parallel do

!----------------------------------------------------------------------
! substepping with external and sink particle forces, using dtextforce
! accretion onto sinks/potentials also happens during substepping
!----------------------------------------------------------------------
#ifdef GR
 if (iexternalforce > 0 .and. imetric /= imet_minkowski) then
    call step_extern_gr(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,t,damp)
 else
    call step_extern_sph_gr(dtsph,npart,xyzh,vxyzu,dens,pxyzu,metrics)
 endif

#else
 if (nptmass > 0 .or. iexternalforce > 0 .or. (h2chemistry .and. icooling > 0) .or. damp > 0.) then
    call step_extern(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,fext,t,damp, &
                     nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nbinmax,ibin_wake)
 else
    call step_extern_sph(dtsph,npart,xyzh,vxyzu)
 endif
#endif

 timei = timei + dtsph
!----------------------------------------------------
! interpolation of SPH quantities needed in the SPH
! force evaluations, using dtsph
!----------------------------------------------------
!$omp parallel do default(none) schedule(guided,1) &
!$omp shared(xyzh,vxyzu,vpred,fxyzu,divcurlv,npart,store_itype) &
!$omp shared(pxyzu,ppred) &
!$omp shared(Bevol,dBevol,Bpred,dtsph,massoftype,iphase) &
!$omp shared(dustevol,ddustprop,dustprop,dustproppred,dustfrac,ddustfrac,dustpred,use_dustfrac) &
!$omp shared(alphaind,ieos,alphamax) &
!$omp shared(temperature) &
#ifdef IND_TIMESTEPS
!$omp shared(twas,timei) &
#endif
!$omp private(hi,rhoi,tdecay1,source,ddenom,hdti) &
!$omp private(i,spsoundi,alphaloci,divvdti) &
!$omp firstprivate(pmassi,itype,avdecayconst,alpha)
 predict_sph: do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (store_itype) then
          itype = iamtype(iphase(i))
          pmassi = massoftype(itype)
          if (itype==iboundary) then
#ifdef GR
             ppred(:,i) = pxyzu(:,i)
#else
             vpred(:,i) = vxyzu(:,i)
#endif
             if (mhd)          Bpred(:,i)  = Bevol (:,i)
             if (use_dustgrowth) dustproppred(:,:) = dustprop(:,:)
             if (use_dustfrac) dustpred(:,i) = dustevol(:,i)
             cycle predict_sph
          endif
       endif
       !
       ! make prediction for h
       !
       if (ndivcurlv >= 1) then
          xyzh(4,i) = xyzh(4,i) - dtsph*dhdrho(xyzh(4,i),pmassi)*rhoh(xyzh(4,i),pmassi)*divcurlv(1,i)
       endif
       !
       ! make a prediction for v and u to the full step for use in the
       ! force evaluation. These have already been updated to the
       ! half step, so only need a half step (0.5*dtsph) here
       !
#ifdef IND_TIMESTEPS
       hdti = timei - twas(i)   ! interpolate to end time
#else
       hdti = 0.5*dtsph
#endif
#ifdef GR
       ppred(:,i) = pxyzu(:,i) + hdti*fxyzu(:,i)
#else
       vpred(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
#endif
       if (use_dustgrowth .and. itype==idust) then
          dustproppred(:,i) = dustproppred(:,i) + hdti*ddustprop(:,i)
       endif
       if (itype==igas) then
          if (mhd) Bpred(:,i) = Bevol (:,i) + hdti*dBevol(:,i)
          if (use_dustfrac) then
             rhoi          = rhoh(xyzh(4,i),pmassi)
             dustpred(:,i) = dustevol(:,i) + hdti*ddustfrac(:,i)
!------------------------------------------------
!--sqrt(rho*epsilon) method
!             dustfrac(:,i) = min(dustpred(:,i)**2/rhoi,1.) ! dustevol = sqrt(rho*eps)
!------------------------------------------------
!--asin(sqrt(epsilon)) method
             dustfrac(:,i) = sin(dustpred(:,i))**2
!------------------------------------------------
          endif
       endif
       !
       ! viscosity switch ONLY (conductivity and resistivity do not use MM97-style switches)
       !
       if (maxalpha==maxp) then
          hi   = xyzh(4,i)
          rhoi = rhoh(hi,pmassi)
          if (store_temperature) then
             spsoundi = get_spsound(ieos,xyzh(:,i),rhoi,vpred(:,i),temperature(i))
          else
             spsoundi = get_spsound(ieos,xyzh(:,i),rhoi,vpred(:,i))
          endif
          tdecay1  = avdecayconst*spsoundi/hi
          ddenom   = 1./(1. + dtsph*tdecay1) ! implicit integration for decay term
          if (nalpha >= 2) then
             ! Cullen and Dehnen (2010) switch
             alphaloci = alphaind(2,i)
             if (alphaind(1,i) < alphaloci) then
                alphaind(1,i) = real(alphaloci,kind=kind(alphaind))
             else
                alphaind(1,i) = real((alphaind(1,i) + dtsph*alphaloci*tdecay1)*ddenom,kind=kind(alphaind))
             endif
          else
             if (ndivcurlv < 1) call fatal('step','alphaind used but divv not stored')
             ! MM97
             source = max(0.0_4,-divcurlv(1,i))
             alphaind(1,i) = real(min((alphaind(1,i) + dtsph*(source + alpha*tdecay1))*ddenom,alphamax),kind=kind(alphaind))
          endif
       endif
    endif
 enddo predict_sph
 !$omp end parallel do
!
! recalculate all SPH forces, and new timestep
!
 if ((iexternalforce /= 0 .or. nptmass > 0) .and. id==master .and. iverbose >= 2) &
   write(iprint,"(a,f14.6,/)") '> full step            : t=',timei

 if (npart > 0) then
    if (gr) vpred = vxyzu ! Need primitive utherm as a guess in cons2prim
    call derivs(1,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,&
                divcurlB,Bpred,dBevol,dustproppred,ddustprop,dustfrac,ddustfrac,temperature,timei,dtsph,dtnew,&
                ppred,dens,metrics)
    if (gr) vxyzu = vpred ! May need primitive variables elsewhere?
 endif
!
! if using super-timestepping, determine what dt will be used on the next loop
!
#ifdef IND_TIMESTEPS
 if ( use_sts ) call sts_get_dtau_next(dtsph_next,dtsph,dtmax,dtdiff,nbinmax)
 if (mod_dtmax_in_step .and. sts_it_n) then
    call decrease_dtmax(npart,maxbins,timei-dtsph,dtmax_rat,dtmax,ibin,ibin_wake,ibin_sts,ibin_dts,mod_dtmax_in_step)
 endif
#endif
!
!-------------------------------------------------------------------------
!  leapfrog corrector step: most of the time we should not need to take
!  any extra iterations, but to be reversible for velocity-dependent
!  forces we must iterate until velocities agree.
!-------------------------------------------------------------------------
 its        = 0
 converged  = .false.
 errmaxmean = 0.0
 iterations: do while (its < maxits .and. .not.converged)
    its     = its + 1
    errmax  = 0.
    v2mean  = 0.
    p2mean  = 0.
    np      = 0
    itype   = igas
    pmassi  = massoftype(igas)
    ntypes  = get_ntypes(npartoftype)
    store_itype = (maxphase==maxp .and. ntypes > 1)
    nwake   = 0
!$omp parallel default(none) &
!$omp shared(xyzh,vxyzu,vpred,fxyzu,npart,hdtsph,store_itype) &
!$omp shared(pxyzu,ppred) &
!$omp shared(Bevol,dBevol,iphase,its) &
!$omp shared(dustevol,ddustfrac,use_dustfrac) &
!$omp shared(dustprop,ddustprop,dustproppred) &
!$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass,massoftype) &
!$omp shared(dtsph,icooling) &
#ifdef IND_TIMESTEPS
!$omp shared(dtmax,ibin,ibin_old,ibin_sts,twas,timei,use_sts,dtsph_next,ibin_wake,sts_it_n) &
!$omp shared(ibin_dts,nbinmax,ibinnow) &
!$omp private(dti,hdti) &
#endif
!$omp private(i,vxi,vyi,vzi,vxoldi,vyoldi,vzoldi) &
!$omp private(pxi,pyi,pzi,p2i) &
!$omp private(erri,v2i,eni) &
!$omp reduction(max:errmax) &
!$omp reduction(+:np,v2mean,p2mean,nwake) &
!$omp firstprivate(pmassi,itype)
!$omp do
    corrector: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (store_itype) itype = iamtype(iphase(i))
          if (itype==iboundary) cycle corrector
#ifdef IND_TIMESTEPS
          !
          !--update active particles
          !
          if (iactive(iphase(i))) then
             ibin_wake(i) = 0       ! cannot wake active particles
             hdti = timei - twas(i) ! = 0.5*get_dt(dtmax,ibin_old(i)) if .not.use_sts, dtmax has not changed, particle was not just woken up
             if (use_sts) then
                if (ibin(i) < ibin_sts(i)) ibin(i) = min(ibin_sts(i), nbinmax ) ! increase ibin if needed for super timestepping
                if (.not.sts_it_n .or. (sts_it_n .and. ibin_sts(i) > ibin(i))) then
                   dti = hdti + 0.5*dtsph_next
                else
                   dti = hdti + ibin_dts(ithdt,ibin(i))
                endif
             else
                dti = hdti + ibin_dts(ithdt,ibin(i))
             endif
#ifdef GR
             pxyzu(:,i) = pxyzu(:,i) + dti*fxyzu(:,i)
#else
             vxyzu(:,i) = vxyzu(:,i) + dti*fxyzu(:,i)
#endif
             if (use_dustgrowth .and. itype==idust) dustproppred(:,i) = dustproppred(:,i) + dti*ddustprop(:,i)
             if (itype==igas) then
                if (mhd)          Bevol(:,i)    = Bevol(:,i)    + dti*dBevol(:,i)
                if (use_dustfrac) dustevol(:,i) = dustevol(:,i) + dti*ddustfrac(:,i)
             endif
             twas(i) = twas(i) + dti
          endif
          !
          !--synchronise all particles
          !
          hdti = timei - twas(i)

#ifdef GR
          pxyzu(:,i) = pxyzu(:,i) + hdti*fxyzu(:,i)
#else
          vxyzu(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
#endif

          if (itype==igas) then
             if (mhd)          Bevol(:,i)  = Bevol(:,i)  + hdti*dBevol(:,i)
             if (use_dustfrac) dustevol(:,i) = dustevol(:,i) + hdti*ddustfrac(:,i)
          endif
          !
          !--Wake inactive particles for next step, if required
          !
          if (sts_it_n .and. ibin_wake(i) > ibin(i) .and. allow_waking) then
             ibin_wake(i) = min(int(nbinmax,kind=1),ibin_wake(i))
             nwake        = nwake + 1
             twas(i)      = ibin_dts(ittwas,ibin_wake(i))
             ibin(i)      = ibin_wake(i)
             ibin_wake(i) = 0 ! reset flag
          endif
#else
          !
          ! For velocity-dependent forces compare the new v
          ! with the predicted v used in the force evaluation.
          ! Determine whether or not we need to iterate.
          !

#ifdef GR
          pxi = pxyzu(1,i) + hdtsph*fxyzu(1,i)
          pyi = pxyzu(2,i) + hdtsph*fxyzu(2,i)
          pzi = pxyzu(3,i) + hdtsph*fxyzu(3,i)
          eni = pxyzu(4,i) + hdtsph*fxyzu(4,i)

          erri = (pxi - ppred(1,i))**2 + (pyi - ppred(2,i))**2 + (pzi - ppred(3,i))**2
          !if (erri > errmax) print*,id,' errmax = ',erri,' part ',i,vxi,vxoldi,vyi,vyoldi,vzi,vzoldi
          errmax = max(errmax,erri)

          p2i = pxi*pxi + pyi*pyi + pzi*pzi
          p2mean = p2mean + p2i
          np = np + 1

          pxyzu(1,i) = pxi
          pxyzu(2,i) = pyi
          pxyzu(3,i) = pzi
          !--this is the energy equation if non-isothermal
          pxyzu(4,i) = eni
#else
          vxi = vxyzu(1,i) + hdtsph*fxyzu(1,i)
          vyi = vxyzu(2,i) + hdtsph*fxyzu(2,i)
          vzi = vxyzu(3,i) + hdtsph*fxyzu(3,i)
          if (maxvxyzu >= 4) eni = vxyzu(4,i) + hdtsph*fxyzu(4,i)

          erri = (vxi - vpred(1,i))**2 + (vyi - vpred(2,i))**2 + (vzi - vpred(3,i))**2
          !if (erri > errmax) print*,id,' errmax = ',erri,' part ',i,vxi,vxoldi,vyi,vyoldi,vzi,vzoldi
          errmax = max(errmax,erri)

          v2i    = vxi*vxi + vyi*vyi + vzi*vzi
          v2mean = v2mean + v2i
          np     = np + 1

          vxyzu(1,i) = vxi
          vxyzu(2,i) = vyi
          vxyzu(3,i) = vzi
          !--this is the energy equation if non-isothermal
          if (maxvxyzu >= 4) then
             vxyzu(4,i) = eni
             if (icooling==3) call energ_coolfunc(vxyzu(4,i),rhoh(xyzh(4,i),massoftype(itype)),dtsph,v2i)
          endif
#endif

          if (itype==igas) then
             !
             ! corrector step for magnetic field and dust
             !
             if (mhd)          Bevol(:,i)  = Bevol(:,i)  + hdtsph*dBevol(:,i)
             if (use_dustfrac) dustevol(:,i) = dustevol(:,i) + hdtsph*ddustfrac(:,i)
          endif
#endif
       endif
    enddo corrector
!$omp enddo
!$omp end parallel

#ifdef GR
    call check_velocity_error(errmax,p2mean,np,its,tolv,dtsph,timei,damp,dterr,errmaxmean,converged)
#else
    call check_velocity_error(errmax,v2mean,np,its,tolv,dtsph,timei,damp,dterr,errmaxmean,converged)
#endif

    if (.not.converged .and. npart > 0) then
       !$omp parallel do private(i) schedule(static)
       do i=1,npart
          if (store_itype) itype = iamtype(iphase(i))
#ifdef IND_TIMESTEPS
          if (iactive(iphase(i))) then
#ifdef GR
             ppred(:,i) = pxyzu(:,i)
#else
             vpred(:,i) = vxyzu(:,i)
#endif
             if (mhd)          Bpred(:,i)  = Bevol(:,i)
             if (use_dustfrac) dustpred(:,i) = dustevol(:,i)
          endif
#else
#ifdef GR
          ppred(:,i) = pxyzu(:,i)
#else
          vpred(:,i) = vxyzu(:,i)
#endif
          if (mhd)          Bpred(:,i)  = Bevol(:,i)
          if (use_dustfrac) dustpred(:,i) = dustevol(:,i)
!
! shift v back to the half step
!
#ifdef GR
          pxyzu(:,i) = pxyzu(:,i) - hdtsph*fxyzu(:,i)
#else
          vxyzu(:,i) = vxyzu(:,i) - hdtsph*fxyzu(:,i)
#endif
          if (itype==igas) then
             if (mhd)          Bevol(:,i)  = Bevol(:,i)  - hdtsph*dBevol(:,i)
             if (use_dustfrac) dustevol(:,i) = dustevol(:,i) - hdtsph*ddustfrac(:,i)
          endif
#endif
       enddo
       !$omp end parallel do
!
!   get new force using updated velocity: no need to recalculate density etc.
!

       if (gr) vpred = vxyzu ! Need primitive utherm as a guess in cons2prim
       call derivs(2,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,divcurlB, &
                     Bpred,dBevol,dustproppred,ddustprop,dustfrac,ddustfrac,&
                     temperature,timei,dtsph,dtnew,ppred,dens,metrics)
       if (gr) vxyzu = vpred ! May need primitive variables elsewhere?
    endif
#ifdef DUSTGROWTH
    call update_dustprop(npart,dustproppred) !--update dustprop values
#endif
 enddo iterations
 ! Summary statements & crash if velocity is not converged
 if (nwake > 1) call summary_variable('wake',iowake,  0,real(nwake))
 if (its   > 1) call summary_variable('tolv',iosumtvi,0,real(its)  )
 if (maxits > 1 .and. its >= maxits) then
    call summary_printout(iprint,nptmass)
    call fatal('step','VELOCITY ITERATIONS NOT CONVERGED!!')
 endif

#ifdef GR
 call conservative_to_primitive(npart,xyzh,metrics,pxyzu,vxyzu,dens)
#endif

 return
end subroutine step

#ifdef GR
subroutine step_extern_sph_gr(dt,npart,xyzh,vxyzu,dens,pxyzu,metrics)
 use part,         only:isdead_or_accreted
 use cons2prim,    only:conservative_to_primitive
 use io,           only:warning
 use metric_tools, only:pack_metric
 real,    intent(in)    :: dt
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:),dens(:),metrics(:,:,:,:)
 real,    intent(in)    :: pxyzu(:,:)
 real,    intent(out)   :: vxyzu(:,:)
 integer, parameter :: nitermax = 50
 real,    parameter ::     xtol = 1.e-15
 integer :: i,niter
 real    :: xpred(1:3),vold(1:3),diff
 logical :: converged

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,dens,dt) &
 !$omp shared(pxyzu,metrics) &
 !$omp private(i,niter,diff,xpred,vold,converged)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call conservative_to_primitive(xyzh(:,i),metrics(:,:,:,i),pxyzu(:,i),vxyzu(:,i),dens(i))
       !
       ! main position update
       !
       xpred = xyzh(1:3,i) + dt*vxyzu(1:3,i)
       vold  = vxyzu(1:3,i)
       converged = .false.
       niter = 0
       do while (.not. converged .and. niter<=nitermax)
          niter = niter + 1
          call conservative_to_primitive(xyzh(:,i),metrics(:,:,:,i),pxyzu(:,i),vxyzu(:,i),dens(i))
          xyzh(1:3,i) = xpred + 0.5*dt*(vxyzu(1:3,i)-vold)
          diff = maxval(abs(xyzh(1:3,i)-xpred)/xpred)
          if (diff < xtol) converged = .true.
          ! UPDATE METRIC HERE
          call pack_metric(xyzh(1:3,i),metrics(:,:,:,i))
       enddo
       if (niter > nitermax) call warning('step_extern_sph_gr','Reached max number of x iterations. x_err ',val=diff)

    endif
 enddo
 !$omp end parallel do

end subroutine step_extern_sph_gr

subroutine step_extern_gr(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,time,damp)
 use dim,            only:maxptmass,maxp,maxvxyzu
 use io,             only:iverbose,id,master,iprint,warning
 use externalforces, only:externalforce,accrete_particles,update_externalforce
 use options,        only:iexternalforce
 use part,           only:maxphase,isdead_or_accreted,iboundary,igas,iphase,iamtype,massoftype,rhoh
 use io_summary,     only:summary_variable,iosumextsr,iosumextst,iosumexter,iosumextet,iosumextr,iosumextt, &
                          summary_accrete,summary_accrete_fail
 use timestep,       only:bignumber
 use eos,            only:equationofstate,ieos
 use cons2prim,      only:conservative_to_primitive
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 integer, intent(in)    :: npart,ntypes
 real,    intent(in)    :: dtsph,time,damp
 real,    intent(inout) :: dtextforce
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:),pxyzu(:,:),dens(:),metrics(:,:,:,:),metricderivs(:,:,:,:)
 integer :: i,itype,nsubsteps,naccreted,its
 real    :: timei,t_end_step,hdt,pmassi
 real    :: dt,dtf,dtextforcenew,dtextforce_min
 real    :: pi,pprev(3),xyz_prev(3),spsoundi,pondensi
 real    :: x_err,pmom_err,fstar(3),vxyzu_star(4),accretedmass
 ! real, save :: dmdt = 0.
 logical :: last_step,done,converged,accreted
 integer, parameter :: itsmax = 50
 real, parameter :: ptol = 1.e-15, xtol = 1.e-15
 integer, save :: pitsmax = 0, xitsmax = 0
!
! determine whether or not to use substepping
!
 if ((iexternalforce > 0) .and. dtextforce < dtsph) then
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

    if (.not.last_step .and. iverbose > 1 .and. id==master) then
       write(iprint,"(a,f14.6)") '> external forces only : t=',timei
    endif

    !---------------------------
    ! predictor during substeps
    !---------------------------
    !
    ! predictor step for external forces, also recompute external forces
    !
    !$omp parallel default(none) &
    !$omp shared(npart,xyzh,vxyzu,fext,iphase,ntypes,massoftype) &
    !$omp shared(dt,hdt) &
    !$omp shared(its,pxyzu,dens,metrics,metricderivs) &
    !$omp private(i,dtf,vxyzu_star,fstar) &
    !$omp private(converged,pprev,pmom_err,xyz_prev,x_err,pi) &
    !$omp firstprivate(pmassi,itype) &
    !$omp reduction(max:xitsmax,pitsmax) &
    !$omp reduction(min:dtextforcenew)
    !$omp do
    predictor: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             pmassi = massoftype(itype)
          endif

          its       = 0
          converged = .false.

          pxyzu(1:3,i) = pxyzu(1:3,i) + hdt*fext(1:3,i)

! Note: grforce needs derivatives of the metric, which do not change between pmom iterations
          pmom_iterations: do while (its <= itsmax .and. .not. converged)
             its   = its + 1
             pprev = pxyzu(1:3,i)
             call conservative_to_primitive(xyzh(:,i),metrics(:,:,:,i),pxyzu(:,i),vxyzu(:,i),dens(i),pi)
             call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i),pi,fstar,dtf)
             pxyzu(1:3,i) = pprev + hdt*(fstar - fext(1:3,i))
             pmom_err = maxval( abs( (pxyzu(1:3,i) - pprev)/pprev ) )
             if (pmom_err < ptol) converged = .true.
             fext(1:3,i) = fstar
          enddo pmom_iterations
          if (its > itsmax ) call warning('step_extern_gr','Reached max number of pmom iterations. pmom_err ',val=pmom_err)

          pitsmax = max(its,pitsmax)

          call conservative_to_primitive(xyzh(:,i),metrics(:,:,:,i),pxyzu(:,i),vxyzu(:,i),dens(i),pi)
          xyzh(1:3,i) = xyzh(1:3,i) + dt*vxyzu(1:3,i)


          its        = 0
          converged  = .false.
          vxyzu_star = vxyzu(:,i)
! Note: since particle positions change between iterations the metric and its derivatives need to be updated.
!       conservative_to_primitive does not require derivatives of the metric, so those can updated once the iterations
!       are complete, in order to reduce the number of computations.
          xyz_iterations: do while (its <= itsmax .and. .not. converged)
             its         = its+1
             xyz_prev    = xyzh(1:3,i)
             call conservative_to_primitive(xyzh(:,i),metrics(:,:,:,i),pxyzu(:,i),vxyzu_star,dens(i))
             xyzh(1:3,i)  = xyz_prev + hdt*(vxyzu_star(1:3) - vxyzu(1:3,i))
             x_err = maxval( abs( (xyzh(1:3,i)-xyz_prev)/xyz_prev ) )
             if (x_err < xtol) converged = .true.
             vxyzu(:,i)   = vxyzu_star
             ! UPDATE METRIC HERE
             call pack_metric(xyzh(1:3,i),metrics(:,:,:,i))
          enddo xyz_iterations
          call pack_metricderivs(xyzh(1:3,i),metricderivs(:,:,:,i))
          if (its > itsmax ) call warning('step_extern_gr','Reached max number of x iterations. x_err ',val=x_err)
          xitsmax = max(its,xitsmax)


          ! Skip remainder of update if boundary particle; note that fext==0 for these particles
          if (itype==iboundary) cycle predictor
       endif
    enddo predictor
    !$omp enddo
    !$omp end parallel

    !
    ! corrector step on gas particles (also accrete particles at end of step)
    !
    accretedmass = 0.
    naccreted    = 0
    dtextforce_min = bignumber

    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,metrics,metricderivs,vxyzu,fext,iphase,ntypes,massoftype,hdt,timei) &
    !$omp private(i,accreted) &
    !$omp shared(ieos,dens,pxyzu,iexternalforce) &
    !$omp private(pi,pondensi,spsoundi,dtf) &
    !$omp firstprivate(itype,pmassi) &
    !$omp reduction(min:dtextforce_min) &
    !$omp reduction(+:accretedmass,naccreted)
    accreteloop: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             pmassi = massoftype(itype)
             !  if (itype==iboundary) cycle accreteloop
          endif

          call equationofstate(ieos,pondensi,spsoundi,dens(i),xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))
          pi = pondensi*dens(i)
          call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i),pi,fext(1:3,i),dtf)
          dtextforce_min = min(dtf,dtextforce_min)
          !
          ! correct v to the full step using only the external force
          !
          pxyzu(1:3,i) = pxyzu(1:3,i) + hdt*fext(1:3,i)

          if (iexternalforce > 0) then
             call accrete_particles(iexternalforce,xyzh(1,i),xyzh(2,i), &
                                    xyzh(3,i),xyzh(4,i),pmassi,timei,accreted,i)
             if (accreted) then
                accretedmass = accretedmass + pmassi
                naccreted = naccreted + 1
             endif
          endif
       endif
    enddo accreteloop
    !$omp end parallel do

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

#endif
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
!  Substepping of external and sink particle forces.
!  Also updates position of all particles even if no external
!  forces applied. This is the internal loop of the RESPA
!  algorithm over the "fast" forces.
!+
!----------------------------------------------------------------
subroutine step_extern(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,fext,time,damp,nptmass, &
                       xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nbinmax,ibin_wake)
 use dim,            only:maxptmass,maxp,maxvxyzu
 use io,             only:iverbose,id,master,iprint,warning
 use externalforces, only:externalforce,accrete_particles,update_externalforce, &
                          update_vdependent_extforce_leapfrog,is_velocity_dependent
 use ptmass,         only:ptmass_predictor,ptmass_corrector,ptmass_accrete, &
                          get_accel_sink_gas,get_accel_sink_sink,f_acc,pt_write_sinkev, &
                          idxmsi,idymsi,idzmsi,idmsi,idspinxsi,idspinysi,idspinzsi, &
                          idvxmsi,idvymsi,idvzmsi,idfxmsi,idfymsi,idfzmsi, &
                          ndptmass,update_ptmass
 use options,        only:iexternalforce
 use part,           only:maxphase,abundance,nabundances,h2chemistry,temperature,store_temperature,epot_sinksink,&
                          isdead_or_accreted,iboundary,igas,iphase,iamtype,massoftype,rhoh,divcurlv, &
                          fxyz_ptmass_sinksink
 use options,        only:icooling
 use chem,           only:energ_h2cooling
 use io_summary,     only:summary_variable,iosumextsr,iosumextst,iosumexter,iosumextet,iosumextr,iosumextt, &
                          summary_accrete,summary_accrete_fail
 use timestep,       only:bignumber,C_force
 use timestep_sts,   only:sts_it_n
 use mpiutils,       only:bcast_mpi,reduce_in_place_mpi,reduceall_mpi
 integer,         intent(in)    :: npart,ntypes,nptmass
 real,            intent(in)    :: dtsph,time,damp
 real,            intent(inout) :: dtextforce
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:)
 real,            intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 integer(kind=1), intent(in)    :: nbinmax
 integer(kind=1), intent(inout) :: ibin_wake(:)
 integer         :: i,itype,nsubsteps,idudtcool,ichem,naccreted,nfail,nfaili
 integer(kind=1) :: ibin_wakei
 real            :: timei,hdt,fextx,fexty,fextz,fextxi,fextyi,fextzi,phii,pmassi
 real            :: dtphi2,dtphi2i,vxhalfi,vyhalfi,vzhalfi,fxi,fyi,fzi,deni
 real            :: dudtcool,fextv(3),fac,poti
 real            :: dt,dtextforcenew,dtsinkgas,fonrmax,fonrmaxi
 real            :: dtf,accretedmass,t_end_step,dtextforce_min
 real            :: dptmass(ndptmass,nptmass)
 real, save      :: dmdt = 0.
 logical         :: accreted,extf_is_velocity_dependent
 logical         :: last_step,done
!
! determine whether or not to use substepping
!
 if ((iexternalforce > 0 .or. nptmass > 0) .and. dtextforce < dtsph) then
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
 fac            = 0.
 t_end_step     = timei + dtsph
 nsubsteps      = 0
 dtextforce_min = huge(dt)
 done           = .false.

 substeps: do while (timei <= t_end_step .and. .not.done)
    hdt           = 0.5*dt
    timei         = timei + dt
    nsubsteps     = nsubsteps + 1
    dtextforcenew = bignumber
    dtsinkgas     = bignumber
    dtphi2        = bignumber

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
          call ptmass_predictor(nptmass,dt,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
          !
          ! get sink-sink forces (and a new sink-sink timestep.  Note: fxyz_ptmass is zeroed in this subroutine)
          ! pass sink-sink forces to variable fxyz_ptmass_sinksink for later writing.
          !
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,dtf,iexternalforce,timei)
          fxyz_ptmass_sinksink=fxyz_ptmass
          if (iverbose >= 2) write(iprint,*) 'dt(sink-sink) = ',C_force*dtf
       else
          fxyz_ptmass(:,:) = 0.
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
    !$omp shared(npart,xyzh,vxyzu,fext,abundance,iphase,ntypes,massoftype) &
    !$omp shared(temperature) &
    !$omp shared(dt,hdt,timei,iexternalforce,extf_is_velocity_dependent,icooling) &
    !$omp shared(xyzmh_ptmass,vxyz_ptmass,damp) &
    !$omp shared(nptmass,f_acc,nsubsteps,C_force,divcurlv) &
    !$omp private(i,ichem,idudtcool,dudtcool,fxi,fyi,fzi,phii) &
    !$omp private(fextx,fexty,fextz,fextxi,fextyi,fextzi,poti,deni,fextv,accreted) &
    !$omp private(fonrmaxi,dtphi2i,dtf) &
    !$omp private(vxhalfi,vyhalfi,vzhalfi) &
    !$omp firstprivate(pmassi,itype) &
    !$omp reduction(+:accretedmass) &
    !$omp reduction(min:dtextforcenew,dtsinkgas,dtphi2) &
    !$omp reduction(max:fonrmax) &
    !$omp reduction(+:fxyz_ptmass)
    !$omp do
    predictor: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
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
          if (itype==iboundary) cycle predictor
          !
          ! compute and add sink-gas force
          !
          fextx = 0.
          fexty = 0.
          fextz = 0.
          if (nptmass > 0) then
             call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                      fextx,fexty,fextz,phii,pmassi,fxyz_ptmass,fonrmaxi,dtphi2i)
             fonrmax = max(fonrmax,fonrmaxi)
             dtphi2  = min(dtphi2,dtphi2i)
          endif
          !
          ! compute and add external forces
          !
          if (iexternalforce > 0) then
             call externalforce(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i), &
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
                vxhalfi = vxyzu(1,i)
                vyhalfi = vxyzu(2,i)
                vzhalfi = vxyzu(3,i)
                fxi = fextx
                fyi = fexty
                fzi = fextz
                call update_vdependent_extforce_leapfrog(iexternalforce,&
                     vxhalfi,vyhalfi,vzhalfi, &
                     fxi,fyi,fzi,fextv,dt,xyzh(1,i),xyzh(2,i),xyzh(3,i))
                fextx = fextx + fextv(1)
                fexty = fexty + fextv(2)
                fextz = fextz + fextv(3)
             endif
          endif
          if (damp > 0.) then
             fextx = fextx - damp*vxyzu(1,i)
             fexty = fexty - damp*vxyzu(2,i)
             fextz = fextz - damp*vxyzu(3,i)
          endif
          fext(1,i) = fextx
          fext(2,i) = fexty
          fext(3,i) = fextz
          !
          ! ARP added:
          ! Chemistry must be updated on each substep, else will severely underestimate.
          !
          if ((maxvxyzu >= 4).and.(icooling > 0).and.h2chemistry) then
             !--Flag determines no cooling, just update abundances.
             ichem     = 1
             idudtcool = 0
             !--Dummy variable to fill for cooling (will remain empty)
             dudtcool  = 0.
             !--Provide a blank dudt_cool element, not needed here
             call energ_h2cooling(vxyzu(4,i),dudtcool,rhoh(xyzh(4,i),pmassi),abundance(:,i), &
                                  nabundances,dt,xyzh(1,i),xyzh(2,i),xyzh(3,i), &
                                  divcurlv(1,i),idudtcool,ichem)
          endif
       endif
    enddo predictor
    !$omp enddo
    !$omp end parallel

    !
    ! reduction of sink-gas forces from each MPI thread
    !
    call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))

    !---------------------------
    ! corrector during substeps
    !---------------------------
    !
    ! corrector step on sinks (changes velocities only, does not change position)
    !
    if (nptmass > 0) then
       if (id==master) then
          call ptmass_corrector(nptmass,dt,vxyz_ptmass,fxyz_ptmass,xyzmh_ptmass,iexternalforce)
       endif
       call bcast_mpi(vxyz_ptmass(:,1:nptmass))
    endif

    !
    ! corrector step on gas particles (also accrete particles at end of step)
    !
    accretedmass = 0.
    nfail        = 0
    naccreted    = 0
    ibin_wakei   = 0
    dptmass(:,1:nptmass) = 0.

    !$omp parallel default(none) &
    !$omp shared(npart,xyzh,vxyzu,fext,iphase,ntypes,massoftype,hdt,timei,nptmass,sts_it_n) &
    !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,f_acc) &
    !$omp shared(iexternalforce) &
    !$omp shared(nbinmax,ibin_wake) &
    !$omp reduction(+:dptmass) &
    !$omp private(i,accreted,nfaili,fxi,fyi,fzi) &
    !$omp firstprivate(itype,pmassi,ibin_wakei) &
    !$omp reduction(+:accretedmass,nfail,naccreted)
    !$omp do
    accreteloop: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             itype = iamtype(iphase(i))
             pmassi = massoftype(itype)
             !if (itype==iboundary) cycle accreteloop
          endif
          !
          ! correct v to the full step using only the external force
          !
          vxyzu(1:3,i) = vxyzu(1:3,i) + hdt*fext(1:3,i)

          if (iexternalforce > 0) then
             call accrete_particles(iexternalforce,xyzh(1,i),xyzh(2,i), &
                                    xyzh(3,i),xyzh(4,i),pmassi,timei,accreted,i)
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
#ifdef IND_TIMESTEPS
             ibin_wakei = ibin_wake(i)
#endif
             call ptmass_accrete(1,nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                                 vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),fxi,fyi,fzi,&
                                 itype,pmassi,xyzmh_ptmass,vxyz_ptmass,&
                                 accreted,dptmass,timei,f_acc,nbinmax,ibin_wakei,nfaili)
             if (accreted) then
                naccreted = naccreted + 1
                cycle accreteloop
#ifdef IND_TIMESTEPS
             else
                ibin_wake(i) = ibin_wakei
#endif
             endif
             if (nfaili > 1) nfail = nfail + 1
          endif
       endif

    enddo accreteloop
    !$omp enddo
    !$omp end parallel

    !
    ! reduction of sink particle changes across MPI
    !
    call reduce_in_place_mpi('+',dptmass(:,1:nptmass))

    naccreted = int(reduceall_mpi('+',naccreted))
    nfail = int(reduceall_mpi('+',nfail))

    if (id==master) call update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)

    call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
    call bcast_mpi(vxyz_ptmass(:,1:nptmass))
    call bcast_mpi(fxyz_ptmass(:,1:nptmass))

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

 if (nsubsteps > 1) then
    if (iverbose>=1 .and. id==master) then
       write(iprint,"(a,i6,a,f8.2,a,es10.3,a,es10.3)") &
           ' using ',nsubsteps,' substeps (dthydro/dtextf = ',dtsph/dtextforce_min,'), dt = ',dtextforce_min,' dtsph = ',dtsph
    endif
    if (nptmass >  0 .and. iexternalforce <= 0) then
       call summary_variable('ext',iosumextsr,nsubsteps,dtsph/dtextforce_min)
       call summary_variable('ext',iosumextst,nsubsteps,dtextforce_min,1.0/dtextforce_min)
    elseif (nptmass <= 0 .and. iexternalforce >  0) then
       call summary_variable('ext',iosumexter,nsubsteps,dtsph/dtextforce_min)
       call summary_variable('ext',iosumextet,nsubsteps,dtextforce_min,1.0/dtextforce_min)
    else
       call summary_variable('ext',iosumextr ,nsubsteps,dtsph/dtextforce_min)
       call summary_variable('ext',iosumextt ,nsubsteps,dtextforce_min,1.0/dtextforce_min)
    endif
 endif

end subroutine step_extern

!-----------------------------------------------------
!+
!  Check error in v^1 compared to the predicted v^*
!  in the leapfrog corrector step. If this is not
!  within some tolerance we iterate the corrector step
!+
!-----------------------------------------------------
subroutine check_velocity_error(errmax,v2mean,np,its,tolv,dt,timei,damp,dterr,errmaxmean,converged)
 use io,         only:id,master,iprint,iverbose,warning
#ifndef IND_TIMESTEPS
 use timestep,   only:dtcourant,dtforce,bignumber
#endif
 use mpiutils,   only:reduceall_mpi
 use io_summary, only:summary_variable,iosumtve,iosumtvv
 real,    intent(inout) :: errmax,v2mean,errmaxmean
 integer, intent(in)    :: np,its
 real,    intent(in)    :: tolv,dt,timei,damp
 real,    intent(out)   :: dterr
 logical, intent(out)   :: converged
 real            :: errtol,vmean
#ifndef IND_TIMESTEPS
 real            :: dtf
#endif
 integer(kind=8) :: nptot
!
!  Check v^1 against the predicted velocity
!  this will only be different because of velocity-dependent forces
!  (i.e., viscosity). If these are not a small fraction of the total force
!  we need to take iterations in order to get this right in leapfrog.
!  Also, if this occurs, we take the timestep down to try to avoid the need
!  for iterations.
!
 nptot = reduceall_mpi('+',np)
 v2mean = reduceall_mpi('+',v2mean)
 errmax = reduceall_mpi('max',errmax)
 if (damp > 0.) call warning('step','damping is ON')

 if (nptot > 0) then
    v2mean = v2mean/real(nptot)
 else
    v2mean = 0.
 endif
 if (v2mean > tiny(v2mean)) then
    errmax = errmax/sqrt(v2mean)
 else
    errmax = 0.
 endif
 errmaxmean = errmaxmean + errmax
 errtol = tolv
 dterr = huge(dterr)
 if (tolv < 1.e2 .and. damp < tiny(damp)) then
#ifndef IND_TIMESTEPS
    dtf = min(dtcourant,dtforce)
    !--if errors are controlling the timestep
    if (dtf > dt .and. dtf < bignumber) then
       errtol = errtol*(dt/dtf)**2
    endif
    if (its==1) then
       if (errtol > tiny(errtol) .and. errmax > epsilon(errmax)) then ! avoid divide-by-zero
          dterr = dt*sqrt(errtol/errmax)
       endif
    endif
#endif
 endif
!
! if the error in the predicted velocity exceeds the tolerance, take iterations
!
! if (maxits > 1 .and. tolv < 1.e2) then
 if (tolv < 1.e2 .and. damp < tiny(damp)) then
    converged = (errmax < tolv)
    if (id==master .and. .not.converged) then
       vmean = sqrt(v2mean)
       call summary_variable('tolv',iosumtve,0,errmax)
       call summary_variable('tolv',iosumtvv,0,vmean)
    endif
    if (id==master .and.((iverbose >= 1.and.(.not.converged.or.its > 1)) .or. iverbose >= 2)) &
       write(iprint,"(a,i2,a,es10.3,a,es10.3,a)") &
            ' velocity-dependent force (tolv) iteration ',its,' errmax = ',errmax,' (vmean = ',sqrt(v2mean),')'
 else
    converged = .true.
 endif
 if (id==master .and. tolv > 1.e2 .and. errmax > 1.0d-2) &
   write(iprint,"(a,es10.3,a,es10.3,a,es10.3)") &
        ' velocity-dependent force (tolv) at time = ',timei,' with errmax = ',errmax,' and vmean = ',sqrt(v2mean)

end subroutine check_velocity_error

end module step_lf_global
