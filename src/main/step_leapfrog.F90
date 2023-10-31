!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module step_lf_global
!
! Computes one (hydro) timestep
!
!   Change this subroutine to change the timestepping algorithm
!
!   This version uses a Velocity Verlet (leapfrog) integrator with
!   substepping (operator splitting) of external/sink particle forces,
!   following the Reversible RESPA algorithm of Tuckerman et al. 1992
!
! :References:
!     Verlet (1967), Phys. Rev. 159, 98-103
!     Tuckerman, Berne & Martyna (1992), J. Chem. Phys. 97, 1990-2001
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary_dyn, chem, cons2prim, cons2primsolver, cooling,
!   cooling_ism, damping, deriv, dim, dust_formation, eos, extern_gr,
!   externalforces, growth, io, io_summary, krome_interface, metric_tools,
!   mpiutils, options, part, ptmass, ptmass_radiation, timestep,
!   timestep_ind, timestep_sts, timing
!
 use dim,  only:maxp,maxvxyzu,do_radiation,ind_timesteps
 use part, only:vpred,Bpred,dustpred,ppred
 use part, only:radpred
 use timestep_ind, only:maxbins,itdt,ithdt,itdt1,ittwas
 implicit none
 real :: ibin_dts(4,0:maxbins)

contains

!------------------------------------------------------------
!+
!  initialisation routine necessary for individual timesteps
!+
!------------------------------------------------------------
subroutine init_step(npart,time,dtmax)
 use timestep_ind, only:get_dt,nbinmax
 use part,         only:ibin,twas,iphase,iamboundary,iamtype
 integer, intent(in) :: npart
 real,    intent(in) :: time,dtmax
 integer             :: i
 !
 ! first time through, move all particles on shortest timestep
 ! then allow them to gradually adjust levels.
 ! Keep boundary particles on level 0 since forces are never calculated
 ! and to prevent boundaries from limiting the timestep
 !
 if (ind_timesteps) then
    if (time < tiny(time)) then
       !$omp parallel do schedule(static) private(i)
       do i=1,npart
          ibin(i) = nbinmax
          if (iamboundary(iamtype(iphase(i)))) ibin(i) = 0
       enddo
    endif
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
    !
    do i=0,maxbins
       ibin_dts(itdt,  i) = get_dt(dtmax,int(i,kind=1))
       ibin_dts(ithdt, i) = 0.5*ibin_dts(itdt,i)
       ibin_dts(itdt1, i) = 1.0/ibin_dts(itdt,i)
       ibin_dts(ittwas,i) = time + 0.5*get_dt(dtmax,int(i,kind=1))
    enddo
 endif

end subroutine init_step

!------------------------------------------------------------
!+
!  main timestepping routine
!+
!------------------------------------------------------------
subroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)
 use dim,            only:maxp,ndivcurlv,maxvxyzu,maxptmass,maxalpha,nalpha,h2chemistry,&
                          use_dustgrowth,use_krome,gr,do_radiation
 use io,             only:iprint,fatal,iverbose,id,master,warning
 use options,        only:iexternalforce,use_dustfrac,implicit_radiation
 use part,           only:xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol, &
                          rad,drad,radprop,isdead_or_accreted,rhoh,dhdrho,&
                          iphase,iamtype,massoftype,maxphase,igas,idust,mhd,&
                          iamboundary,get_ntypes,npartoftypetot,&
                          dustfrac,dustevol,ddustevol,eos_vars,alphaind,nptmass,&
                          dustprop,ddustprop,dustproppred,pxyzu,dens,metrics,ics
 use options,        only:avdecayconst,alpha,ieos,alphamax
 use deriv,          only:derivs
 use timestep,       only:dterr,bignumber,tolv
 use mpiutils,       only:reduceall_mpi
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass
 use io_summary,     only:summary_printout,summary_variable,iosumtvi,iowake, &
                          iosumflrp,iosumflrps,iosumflrc
 use cooling,        only:ufloor
 use boundary_dyn,   only:dynamic_bdy,update_xyzminmax
#ifdef KROME
 use part,           only:gamma_chem
#endif
 use timestep,       only:dtmax,dtmax_ifactor,dtdiff
 use timestep_ind,   only:get_dt,nbinmax,decrease_dtmax,dt_too_small
 use timestep_sts,   only:sts_get_dtau_next,use_sts,ibin_sts,sts_it_n
 use part,           only:ibin,ibin_old,twas,iactive,ibin_wake
#ifdef GR
 use part,           only:metricderivs
 use metric_tools,   only:imet_minkowski,imetric
 use cons2prim,      only:cons2primall
 use extern_gr,      only:get_grforce_all
#else
 use cooling,        only:cooling_in_step
#endif
 use timing,         only:increment_timer,get_timings,itimer_extf
 use growth,         only:check_dustprop
 use damping,        only:idamp
 use cons2primsolver, only:conservative2primitive,primitive2conservative
 use eos,             only:equationofstate

 integer, intent(inout) :: npart
 integer, intent(in)    :: nactive
 real,    intent(in)    :: t,dtsph
 real,    intent(inout) :: dtextforce
 real,    intent(out)   :: dtnew
 integer            :: i,its,np,ntypes,itype,nwake,nvfloorp,nvfloorps,nvfloorc,ialphaloc
 real               :: timei,erri,errmax,v2i,errmaxmean
 real               :: vxi,vyi,vzi,eni,hdtsph,pmassi
 real               :: alphaloci,source,tdecay1,hi,rhoi,ddenom,spsoundi
 real               :: v2mean,hdti
 real(kind=4)       :: t1,t2,tcpu1,tcpu2
 real               :: pxi,pyi,pzi,p2i,p2mean
 real               :: dtsph_next,dti,time_now
 logical, parameter :: allow_waking = .true.
 integer, parameter :: maxits = 30
 logical            :: converged,store_itype
!
! set initial quantities
!
 timei  = t
 hdtsph = 0.5*dtsph
 dterr  = bignumber
! determine twas for each ibin
 if (ind_timesteps .and. sts_it_n) then
    time_now = timei + dtsph
    do i=0,maxbins
       ibin_dts(ittwas,i) = (int(time_now*ibin_dts(itdt1,i),kind=8) + 0.5)*ibin_dts(itdt,i)
    enddo
 endif

!--------------------------------------
! velocity predictor step, using dtsph
!--------------------------------------
 itype   = igas
 ntypes  = get_ntypes(npartoftypetot)
 pmassi  = massoftype(itype)
 store_itype = (maxphase==maxp .and. ntypes > 1)
 ialphaloc = 2
 nvfloorp  = 0

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,fxyzu,iphase,hdtsph,store_itype) &
 !$omp shared(rad,drad,pxyzu) &
 !$omp shared(Bevol,dBevol,dustevol,ddustevol,use_dustfrac) &
 !$omp shared(dustprop,ddustprop,dustproppred,ufloor) &
 !$omp shared(ibin,ibin_old,twas,timei) &
 !$omp firstprivate(itype) &
 !$omp private(i,hdti) &
 !$omp reduction(+:nvfloorp)
 predictor: do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (ind_timesteps) then
          if (iactive(iphase(i))) ibin_old(i) = ibin(i) ! only required for ibin_neigh in force.F90
          !
          !--synchronise all particles to their half timesteps
          !
          hdti = twas(i) - timei
       else
          hdti = hdtsph
       endif
       if (store_itype) itype = iamtype(iphase(i))
       if (iamboundary(itype)) cycle predictor
       !
       ! predict v and u to the half step with "slow" forces
       !
       if (gr) then
          pxyzu(:,i) = pxyzu(:,i) + hdti*fxyzu(:,i)
       else
          vxyzu(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
       endif

       !--floor the thermal energy if requested and required
       if (ufloor > 0.) then
          if (vxyzu(4,i) < ufloor) then
             vxyzu(4,i) = ufloor
             nvfloorp   = nvfloorp + 1
          endif
       endif

       if (itype==idust .and. use_dustgrowth) then
          dustprop(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)
       endif
       if (itype==igas) then
          if (mhd)          Bevol(:,i) = Bevol(:,i) + hdti*dBevol(:,i)
          if (do_radiation) rad(:,i)   = rad(:,i) + hdti*drad(:,i)
          if (use_dustfrac) then
             dustevol(:,i) = abs(dustevol(:,i) + hdti*ddustevol(:,i))
             if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)
          endif
       endif
    endif
 enddo predictor
 !omp end parallel do
 if (use_dustgrowth) call check_dustprop(npart,dustprop(1,:))


!----------------------------------------------------------------------
! substepping with external and sink particle forces, using dtextforce
! accretion onto sinks/potentials also happens during substepping
!----------------------------------------------------------------------
 call get_timings(t1,tcpu1)
#ifdef GR
 if ((iexternalforce > 0 .and. imetric /= imet_minkowski) .or. idamp > 0) then
    call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
    call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,dens,fext,dtextforce)
    call step_extern_gr(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,t)
 else
    call step_extern_sph_gr(dtsph,npart,xyzh,vxyzu,dens,pxyzu,metrics)
 endif

#else
 if (nptmass > 0 .or. iexternalforce > 0 .or. h2chemistry .or. cooling_in_step .or. idamp > 0) then
    call step_extern(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,fext,fxyzu,t, &
                     nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nbinmax,ibin_wake)
 else
    call step_extern_sph(dtsph,npart,xyzh,vxyzu)
 endif
#endif
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_extf,t2-t1,tcpu2-tcpu1)

 timei = timei + dtsph
 nvfloorps  = 0
!----------------------------------------------------
! interpolation of SPH quantities needed in the SPH
! force evaluations, using dtsph
!----------------------------------------------------
!$omp parallel do default(none) schedule(guided,1) &
!$omp shared(maxp,maxphase,maxalpha) &
!$omp shared(xyzh,vxyzu,vpred,fxyzu,divcurlv,npart,store_itype) &
!$omp shared(pxyzu,ppred) &
!$omp shared(Bevol,dBevol,Bpred,dtsph,massoftype,iphase) &
!$omp shared(dustevol,ddustprop,dustprop,dustproppred,dustfrac,ddustevol,dustpred,use_dustfrac) &
!$omp shared(alphaind,ieos,alphamax,ialphaloc) &
!$omp shared(eos_vars,ufloor) &
!$omp shared(twas,timei) &
!$omp shared(rad,drad,radpred)&
!$omp private(hi,rhoi,tdecay1,source,ddenom,hdti) &
!$omp private(i,spsoundi,alphaloci) &
!$omp firstprivate(pmassi,itype,avdecayconst,alpha) &
!$omp reduction(+:nvfloorps)
 predict_sph: do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (store_itype) then
          itype = iamtype(iphase(i))
          pmassi = massoftype(itype)
          if (iamboundary(itype)) then
             if (gr) then
                ppred(:,i) = pxyzu(:,i)
             else
                vpred(:,i) = vxyzu(:,i)
             endif
             if (mhd)          Bpred(:,i)  = Bevol (:,i)
             if (use_dustgrowth) dustproppred(:,i) = dustprop(:,i)
             if (use_dustfrac)   dustpred(:,i) = dustevol(:,i)
             if (do_radiation)   radpred(:,i) = rad(:,i)
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
       if (ind_timesteps) then
          hdti = timei - twas(i)   ! interpolate to end time
       else
          hdti = 0.5*dtsph
       endif

       if (gr) then
          ppred(:,i) = pxyzu(:,i) + hdti*fxyzu(:,i)
       else
          vpred(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
       endif

       !--floor the thermal energy if requested and required
       if (ufloor > 0.) then
          if (vpred(4,i) < ufloor) then
             vpred(4,i) = ufloor
             nvfloorps  = nvfloorps + 1
          endif
       endif

       if (use_dustgrowth .and. itype==idust) then
          dustproppred(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)
       endif
       if (itype==igas) then
          if (mhd) Bpred(:,i) = Bevol (:,i) + hdti*dBevol(:,i)
          if (use_dustfrac) then
             rhoi          = rhoh(xyzh(4,i),pmassi)
             dustpred(:,i) = dustevol(:,i) + hdti*ddustevol(:,i)
             if (use_dustgrowth) dustproppred(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)
          endif
          if (do_radiation) radpred(:,i) = rad(:,i) + hdti*drad(:,i)
       endif
       !
       ! viscosity switch ONLY (conductivity and resistivity do not use MM97-style switches)
       !
       if (maxalpha==maxp) then
          hi   = xyzh(4,i)
          rhoi = rhoh(hi,pmassi)
          spsoundi = eos_vars(ics,i)
          tdecay1  = avdecayconst*spsoundi/hi
          ddenom   = 1./(1. + dtsph*tdecay1) ! implicit integration for decay term
          if (nalpha >= 2) then
             ! Cullen and Dehnen (2010) switch
             alphaloci = alphaind(ialphaloc,i)
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

 if (use_dustgrowth) call check_dustprop(npart,dustproppred(1,:))

!
! recalculate all SPH forces, and new timestep
!
 if ((iexternalforce /= 0 .or. nptmass > 0) .and. id==master .and. iverbose >= 2) &
   write(iprint,"(a,f14.6,/)") '> full step            : t=',timei

 if (npart > 0) then
    if (gr) vpred = vxyzu ! Need primitive utherm as a guess in cons2prim
    dt_too_small = .false.

    call derivs(1,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,&
                divcurlB,Bpred,dBevol,radpred,drad,radprop,dustproppred,ddustprop,&
                dustpred,ddustevol,dustfrac,eos_vars,timei,dtsph,dtnew,&
                ppred,dens,metrics)

    if (do_radiation .and. implicit_radiation) then
       rad = radpred
       vxyzu(4,1:npart) = vpred(4,1:npart)
    endif

    if (gr) vxyzu = vpred ! May need primitive variables elsewhere?
    if (dt_too_small) then
       ! dt < dtmax/2**nbinmax and exit
       ! Perform this here rather than in get_newbin so that we can get some diagnostic info
       ! This is only used for individual timesteps
       call summary_printout(iprint,nptmass)
       call fatal('step','step too small: bin would exceed maximum')
    endif
 endif
!
! if using super-timestepping, determine what dt will be used on the next loop
!
 if (ind_timesteps) then
    if ( use_sts ) call sts_get_dtau_next(dtsph_next,dtsph,dtmax,dtdiff,nbinmax)
    if (dtmax_ifactor /=0 .and. sts_it_n) then
       call decrease_dtmax(npart,maxbins,timei-dtsph,dtmax_ifactor,dtmax,ibin,ibin_wake,ibin_sts,ibin_dts)
    endif
 endif
!
!-------------------------------------------------------------------------
!  leapfrog corrector step: most of the time we should not need to take
!  any extra iterations, but to be reversible for velocity-dependent
!  forces we must iterate until velocities agree.
!-------------------------------------------------------------------------
 its        = 0
 converged  = .false.
 errmaxmean = 0.0
 nwake      = 0
 nvfloorc   = 0
 iterations: do while (its < maxits .and. .not.converged .and. npart > 0)
    its     = its + 1
    errmax  = 0.
    v2mean  = 0.
    p2mean  = 0.
    np      = 0
    itype   = igas
    pmassi  = massoftype(igas)
    ntypes  = get_ntypes(npartoftypetot)
    store_itype = (maxphase==maxp .and. ntypes > 1)
!$omp parallel default(none) &
!$omp shared(xyzh,vxyzu,vpred,fxyzu,npart,hdtsph,store_itype) &
!$omp shared(pxyzu,ppred) &
!$omp shared(Bevol,dBevol,iphase,its) &
!$omp shared(dustevol,ddustevol,use_dustfrac) &
!$omp shared(dustprop,ddustprop,dustproppred) &
!$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass,massoftype) &
!$omp shared(dtsph,ieos,ufloor) &
!$omp shared(ibin,ibin_old,ibin_sts,twas,timei,use_sts,dtsph_next,ibin_wake,sts_it_n) &
!$omp shared(ibin_dts,nbinmax) &
!$omp private(dti,hdti) &
!$omp shared(rad,radpred,drad)&
!$omp private(i,vxi,vyi,vzi) &
!$omp private(pxi,pyi,pzi,p2i) &
!$omp private(erri,v2i,eni) &
!$omp reduction(max:errmax) &
!$omp reduction(+:np,v2mean,p2mean,nwake,nvfloorc) &
!$omp firstprivate(pmassi,itype)
!$omp do
    corrector: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (store_itype) itype = iamtype(iphase(i))
          if (iamboundary(itype)) cycle corrector
          if (ind_timesteps) then
             !
             !--update active particles
             !
             if (iactive(iphase(i))) then
                ibin_wake(i) = 0       ! cannot wake active particles
                hdti = timei - twas(i) ! = 0.5*get_dt(dtmax,ibin_old(i)) if .not.use_sts & dtmax has not changed & particle was not just woken up
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

                if (gr) then
                   pxyzu(:,i) = pxyzu(:,i) + dti*fxyzu(:,i)
                else
                   vxyzu(:,i) = vxyzu(:,i) + dti*fxyzu(:,i)
                endif

                if (use_dustgrowth .and. itype==idust) dustprop(:,i) = dustprop(:,i) + dti*ddustprop(:,i)
                if (itype==igas) then
                   if (mhd)          Bevol(:,i) = Bevol(:,i) + dti*dBevol(:,i)
                   if (do_radiation) rad(:,i)   = rad(:,i)   + dti*drad(:,i)
                   if (use_dustfrac) then
                      dustevol(:,i) = dustevol(:,i) + dti*ddustevol(:,i)
                      if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + dti*ddustprop(:,i)
                   endif
                endif
                twas(i) = twas(i) + dti
             endif
             !
             !--synchronise all particles
             !
             hdti = timei - twas(i)

             if (gr) then
                pxyzu(:,i) = pxyzu(:,i) + hdti*fxyzu(:,i)
             else
                vxyzu(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
             endif

             !--floor the thermal energy if requested and required
             if (ufloor > 0.) then
                if (vxyzu(4,i) < ufloor) then
                   vxyzu(4,i) = ufloor
                   nvfloorc   = nvfloorc + 1
                endif
             endif

             if (itype==idust .and. use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)
             if (itype==igas) then
                if (mhd)          Bevol(:,i) = Bevol(:,i) + hdti*dBevol(:,i)
                if (do_radiation) rad(:,i)   = rad(:,i)   + hdti*drad(:,i)
                if (use_dustfrac) then
                   dustevol(:,i) = dustevol(:,i) + hdti*ddustevol(:,i)
                   if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)
                endif
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
          else  ! not individual timesteps == global timestepping
             !
             ! For velocity-dependent forces compare the new v
             ! with the predicted v used in the force evaluation.
             ! Determine whether or not we need to iterate.
             !

             if (gr) then
                pxi = pxyzu(1,i) + hdtsph*fxyzu(1,i)
                pyi = pxyzu(2,i) + hdtsph*fxyzu(2,i)
                pzi = pxyzu(3,i) + hdtsph*fxyzu(3,i)
                eni = pxyzu(4,i) + hdtsph*fxyzu(4,i)

                erri = (pxi - ppred(1,i))**2 + (pyi - ppred(2,i))**2 + (pzi - ppred(3,i))**2
                errmax = max(errmax,erri)

                p2i = pxi*pxi + pyi*pyi + pzi*pzi
                p2mean = p2mean + p2i
                np = np + 1

                pxyzu(1,i) = pxi
                pxyzu(2,i) = pyi
                pxyzu(3,i) = pzi
                pxyzu(4,i) = eni
             else
                vxi = vxyzu(1,i) + hdtsph*fxyzu(1,i)
                vyi = vxyzu(2,i) + hdtsph*fxyzu(2,i)
                vzi = vxyzu(3,i) + hdtsph*fxyzu(3,i)
                if (maxvxyzu >= 4) eni = vxyzu(4,i) + hdtsph*fxyzu(4,i)

                erri = (vxi - vpred(1,i))**2 + (vyi - vpred(2,i))**2 + (vzi - vpred(3,i))**2
                errmax = max(errmax,erri)

                v2i    = vxi*vxi + vyi*vyi + vzi*vzi
                v2mean = v2mean + v2i
                np     = np + 1

                vxyzu(1,i) = vxi
                vxyzu(2,i) = vyi
                vxyzu(3,i) = vzi
                !--this is the energy equation if non-isothermal
                if (maxvxyzu >= 4) vxyzu(4,i) = eni
             endif

             if (itype==idust .and. use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdtsph*ddustprop(:,i)
             if (itype==igas) then
                !
                ! corrector step for magnetic field and dust
                !
                if (mhd)          Bevol(:,i) = Bevol(:,i)  + hdtsph*dBevol(:,i)
                if (do_radiation) rad(:,i)   = rad(:,i) + hdtsph*drad(:,i)
                if (use_dustfrac) then
                   dustevol(:,i) = dustevol(:,i) + hdtsph*ddustevol(:,i)
                   if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdtsph*ddustprop(:,i)
                endif
             endif
          endif
       endif
    enddo corrector
!$omp enddo
!$omp end parallel
    if (use_dustgrowth) call check_dustprop(npart,dustprop(1,:))

    if (gr) then
       call check_velocity_error(errmax,p2mean,np,its,tolv,dtsph,timei,idamp,dterr,errmaxmean,converged)
    else
       call check_velocity_error(errmax,v2mean,np,its,tolv,dtsph,timei,idamp,dterr,errmaxmean,converged)
    endif

    if (.not.converged .and. npart > 0) then
!$omp parallel do default(none)&
!$omp private(i) &
!$omp shared(npart,hdtsph)&
!$omp shared(store_itype,vxyzu,fxyzu,vpred,iphase) &
!$omp shared(Bevol,dBevol,Bpred,pxyzu,ppred) &
!$omp shared(dustprop,ddustprop,dustproppred,use_dustfrac,dustevol,dustpred,ddustevol) &
!$omp shared(rad,drad,radpred) &
!$omp firstprivate(itype) &
!$omp schedule(static)
       until_converged: do i=1,npart
          if (store_itype) itype = iamtype(iphase(i))
          if (iamboundary(itype)) cycle until_converged

          if (ind_timesteps) then
             if (iactive(iphase(i))) then

                if (gr) then
                   ppred(:,i) = pxyzu(:,i)
                else
                   vpred(:,i) = vxyzu(:,i)
                endif
                if (use_dustgrowth) dustproppred(:,i) = dustprop(:,i)
                if (mhd)          Bpred(:,i)  = Bevol(:,i)
                if (use_dustfrac) dustpred(:,i) = dustevol(:,i)
                if (do_radiation) radpred(:,i) = rad(:,i)
             endif
          else
             if (gr) then
                ppred(:,i) = pxyzu(:,i)
             else
                vpred(:,i) = vxyzu(:,i)
             endif
             if (use_dustgrowth) dustproppred(:,i) = dustprop(:,i)
             if (mhd)          Bpred(:,i)  = Bevol(:,i)
             if (use_dustfrac) dustpred(:,i) = dustevol(:,i)
             if (do_radiation) radpred(:,i) = rad(:,i)
             !
             ! shift v back to the half step
             !
             if (gr) then
                pxyzu(:,i) = pxyzu(:,i) - hdtsph*fxyzu(:,i)
             else
                vxyzu(:,i) = vxyzu(:,i) - hdtsph*fxyzu(:,i)
             endif
             if (itype==idust .and. use_dustgrowth) dustprop(:,i) = dustprop(:,i) - hdtsph*ddustprop(:,i)
             if (itype==igas) then
                if (mhd)          Bevol(:,i)  = Bevol(:,i)  - hdtsph*dBevol(:,i)
                if (use_dustfrac) then
                   dustevol(:,i) = dustevol(:,i) - hdtsph*ddustevol(:,i)
                   if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) - hdtsph*ddustprop(:,i)
                endif
                if (do_radiation) rad(:,i) = rad(:,i) - hdtsph*drad(:,i)
             endif
          endif
       enddo until_converged
!$omp end parallel do

       if (use_dustgrowth) call check_dustprop(npart,dustprop(1,:))

!
!   get new force using updated velocity: no need to recalculate density etc.
!
       if (gr) vpred = vxyzu ! Need primitive utherm as a guess in cons2prim
       call derivs(2,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,divcurlB, &
                     Bpred,dBevol,radpred,drad,radprop,dustproppred,ddustprop,dustpred,ddustevol,dustfrac,&
                     eos_vars,timei,dtsph,dtnew,ppred,dens,metrics)
       if (gr) vxyzu = vpred ! May need primitive variables elsewhere?
       if (do_radiation .and. implicit_radiation) then
          rad = radpred
          vxyzu(4,1:npart) = vpred(4,1:npart)
       endif
    endif
 enddo iterations
 
 ! MPI reduce summary variables
 nwake     = int(reduceall_mpi('+', nwake))
 nvfloorp  = int(reduceall_mpi('+', nvfloorp))
 nvfloorps = int(reduceall_mpi('+', nvfloorps))
 nvfloorc  = int(reduceall_mpi('+', nvfloorc))

 if (dynamic_bdy) call update_xyzminmax(dtsph)

 ! Summary statements & crash if velocity is not converged
 if (nwake    > 0) call summary_variable('wake', iowake,    0,real(nwake)    )
 if (nvfloorp > 0) call summary_variable('floor',iosumflrp, 0,real(nvfloorp) )
 if (nvfloorps> 0) call summary_variable('floor',iosumflrps,0,real(nvfloorps))
 if (nvfloorc > 0) call summary_variable('floor',iosumflrc, 0,real(nvfloorc) )
 if (its      > 1) call summary_variable('tolv', iosumtvi,  0,real(its)      )
 if (maxits   > 1 .and. its >= maxits) then
    call summary_printout(iprint,nptmass)
    call fatal('step','VELOCITY ITERATIONS NOT CONVERGED!!')
 endif

#ifdef GR
 call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
#endif

end subroutine step

#ifdef GR
subroutine step_extern_sph_gr(dt,npart,xyzh,vxyzu,dens,pxyzu,metrics)
 use part,            only:isdead_or_accreted,igas,massoftype,rhoh,eos_vars,igasP,&
                           ien_type,eos_vars,igamma,itemp
 use cons2primsolver, only:conservative2primitive
 use eos,             only:ieos,get_pressure
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
    !$omp parallel do default(none) schedule(runtime)  &
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
    !$omp do schedule(runtime)
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
subroutine step_extern(npart,ntypes,dtsph,dtextforce,xyzh,vxyzu,fext,fxyzu,time,nptmass, &
                       xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nbinmax,ibin_wake)
 use dim,            only:maxptmass,maxp,maxvxyzu,store_dust_temperature,use_krome,itau_alloc,do_nucleation
 use io,             only:iverbose,id,master,iprint,warning,fatal
 use externalforces, only:externalforce,accrete_particles,update_externalforce, &
                          update_vdependent_extforce_leapfrog,is_velocity_dependent
 use ptmass,         only:ptmass_predictor,ptmass_corrector,ptmass_accrete, &
                          get_accel_sink_gas,get_accel_sink_sink,merge_sinks,f_acc,pt_write_sinkev, &
                          idxmsi,idymsi,idzmsi,idmsi,idspinxsi,idspinysi,idspinzsi, &
                          idvxmsi,idvymsi,idvzmsi,idfxmsi,idfymsi,idfzmsi, &
                          ndptmass,update_ptmass
 use options,        only:iexternalforce,icooling
 use part,           only:maxphase,abundance,nabundances,h2chemistry,eos_vars,epot_sinksink,&
                          isdead_or_accreted,iamboundary,igas,iphase,iamtype,massoftype,rhoh,divcurlv, &
                          fxyz_ptmass_sinksink,dust_temp,tau,nucleation,idK2,idmu,idkappa,idgamma
 use chem,           only:update_abundances,get_dphot
 use cooling_ism,    only:dphot0,energ_cooling_ism,dphotflag,abundsi,abundo,abunde,abundc,nabn
 use io_summary,     only:summary_variable,iosumextr,iosumextt,summary_accrete,summary_accrete_fail
 use timestep,       only:bignumber,C_force
 use timestep_sts,   only:sts_it_n
 use mpiutils,       only:bcast_mpi,reduce_in_place_mpi,reduceall_mpi
 use damping,        only:calc_damp,apply_damp,idamp
 use ptmass_radiation,only:get_rad_accel_from_ptmass,isink_radiation
 use cooling,        only:energ_cooling,cooling_in_step
 use dust_formation, only:evolve_dust
#ifdef KROME
 use part,            only: gamma_chem,mu_chem,dudt_chem,T_gas_cool
 use krome_interface, only: update_krome
#endif
 integer,         intent(in)    :: npart,ntypes,nptmass
 real,            intent(in)    :: dtsph,time
 real,            intent(inout) :: dtextforce
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fext(:,:),fxyzu(:,:)
 real,            intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 integer(kind=1), intent(in)    :: nbinmax
 integer(kind=1), intent(inout) :: ibin_wake(:)
 integer         :: i,itype,nsubsteps,naccreted,nfail,nfaili,merge_n,nlive
 integer         :: merge_ij(nptmass)
 integer(kind=1) :: ibin_wakei
 real            :: timei,hdt,fextx,fexty,fextz,fextxi,fextyi,fextzi,phii,pmassi
 real            :: dtphi2,dtphi2i,vxhalfi,vyhalfi,vzhalfi,fxi,fyi,fzi
 real            :: dudtcool,fextv(3),poti,ui,rhoi
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
          call ptmass_predictor(nptmass,dt,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
          !
          ! get sink-sink forces (and a new sink-sink timestep.  Note: fxyz_ptmass is zeroed in this subroutine)
          ! pass sink-sink forces to variable fxyz_ptmass_sinksink for later writing.
          !
          if (iexternalforce==14) call update_externalforce(iexternalforce,timei,dmdt)
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,dtf,iexternalforce,timei,merge_ij,merge_n)
          if (merge_n > 0) then
             call merge_sinks(timei,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,merge_ij)
             call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,dtf,iexternalforce,timei,merge_ij,merge_n)
          endif
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
    !$omp shared(maxp,maxphase) &
    !$omp shared(npart,xyzh,vxyzu,fext,abundance,iphase,ntypes,massoftype) &
    !$omp shared(eos_vars,dust_temp,store_dust_temperature) &
    !$omp shared(dt,hdt,timei,iexternalforce,extf_is_velocity_dependent,cooling_in_step,icooling) &
    !$omp shared(xyzmh_ptmass,vxyz_ptmass,idamp,damp_fac) &
    !$omp shared(nptmass,nsubsteps,C_force,divcurlv,dphotflag,dphot0) &
    !$omp shared(abundc,abundo,abundsi,abunde) &
    !$omp shared(nucleation,do_nucleation) &
#ifdef KROME
    !$omp shared(gamma_chem,mu_chem,dudt_chem) &
#endif
    !$omp private(dphot,abundi,gmwvar) &
    !$omp private(ui,rhoi) &
    !$omp private(i,dudtcool,fxi,fyi,fzi,phii) &
    !$omp private(fextx,fexty,fextz,fextxi,fextyi,fextzi,poti,fextv,accreted) &
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
          if (idamp > 0) then
             call apply_damp(fextx, fexty, fextz, vxyzu(1:3,i), xyzh(1:3,i), damp_fac)
          endif
          fext(1,i) = fextx
          fext(2,i) = fexty
          fext(3,i) = fextz

          if (maxvxyzu >= 4 .and. itype==igas) then
             ! NOTE: The chemistry and cooling here is implicitly calculated.  That is,
             !       dt is *passed in* to the chemistry & cooling routines so that the
             !       output will be at the correct time of time + dt.  Since this is
             !       implicit, there is no cooling timestep.  Explicit cooling is
             !       calculated in force and requires a cooling timestep.

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
                call update_abundances(vxyzu(4,i),rhoi,abundance(:,i),&
                      nabundances,dphot,dt,abundi,nabn,gmwvar,abundc,abunde,abundo,abundsi)
             endif
#ifdef KROME
             ! evolve chemical composition and determine new internal energy
             ! Krome also computes cooling function but only associated with chemical processes
             ui = vxyzu(4,i)
             call update_krome(dt,xyzh(:,i),ui,rhoi,abundance(:,i),gamma_chem(i),mu_chem(i),T_gas_cool(i))
             dudt_chem(i) = (ui-vxyzu(4,i))/dt
             dudtcool     = dudt_chem(i)
#else
             !evolve dust chemistry and compute dust cooling
             if (do_nucleation) call evolve_dust(dt, xyzh(:,i), vxyzu(4,i), nucleation(:,i), dust_temp(i), rhoi)
             !
             ! COOLING
             !
             if (icooling > 0 .and. cooling_in_step) then
                if (h2chemistry) then
                   !
                   ! Call cooling routine, requiring total density, some distance measure and
                   ! abundances in the 'abund' format
                   !
                   call energ_cooling_ism(vxyzu(4,i),rhoi,divcurlv(1,i),gmwvar,abundi,dudtcool)
                elseif (store_dust_temperature) then
                   ! cooling with stored dust temperature
                   if (do_nucleation) then
                      call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),dudtcool,rhoi,dt,&
                           dust_temp(i),nucleation(idmu,i),nucleation(idgamma,i),nucleation(idK2,i),nucleation(idkappa,i))
                   else
                      call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),dudtcool,rhoi,dt,dust_temp(i))
                   endif
                else
                   ! cooling without stored dust temperature
                   call energ_cooling(xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),dudtcool,rhoi,dt)
                endif
             endif
#endif
             ! update internal energy
             if (cooling_in_step .or. use_krome) vxyzu(4,i) = vxyzu(4,i) + dt * dudtcool
          endif
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
    if (nptmass > 0) call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))

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

end subroutine step_extern

!-----------------------------------------------------
!+
!  Check error in v^1 compared to the predicted v^*
!  in the leapfrog corrector step. If this is not
!  within some tolerance we iterate the corrector step
!+
!-----------------------------------------------------
subroutine check_velocity_error(errmax,v2mean,np,its,tolv,dt,timei,idamp,dterr,errmaxmean,converged)
 use io,         only:id,master,iprint,iverbose,warning
 use timestep,   only:dtcourant,dtforce,bignumber
 use mpiutils,   only:reduceall_mpi
 use io_summary, only:summary_variable,iosumtve,iosumtvv
 real,    intent(inout) :: errmax,v2mean,errmaxmean
 integer, intent(in)    :: np,its,idamp
 real,    intent(in)    :: tolv,dt,timei
 real,    intent(out)   :: dterr
 logical, intent(out)   :: converged
 real            :: errtol,vmean
 real            :: dtf
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
 if (idamp > 0) call warning('step','damping is ON')

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
 if (tolv < 1.e2 .and. idamp == 0) then
    if (.not.ind_timesteps) then
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
    endif
 endif
!
! if the error in the predicted velocity exceeds the tolerance, take iterations
!
! if (maxits > 1 .and. tolv < 1.e2) then
 if (tolv < 1.e2 .and. idamp == 0) then
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
