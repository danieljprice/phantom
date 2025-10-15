!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Dependencies: boundary_dyn, cons2prim, cons2primsolver, cooling,
!   cooling_radapprox, damping, deriv, dim, eos, extern_gr, growth, io,
!   io_summary, metric_tools, mpiutils, options, part, porosity, ptmass,
!   substepping, timestep, timestep_ind, timestep_sts, timing
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
                          use_dustgrowth,use_krome,gr,do_radiation,use_apr,use_sinktree
 use io,             only:iprint,fatal,iverbose,id,master,warning
 use options,        only:iexternalforce,use_dustfrac,implicit_radiation
 use part,           only:xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol, &
                          rad,drad,radprop,isdead_or_accreted,rhoh,dhdrho,&
                          iphase,iamtype,massoftype,maxphase,igas,idust,mhd,&
                          iamboundary,get_ntypes,npartoftypetot,apr_level,&
                          dustfrac,dustevol,ddustevol,eos_vars,alphaind,nptmass,&
                          dustprop,ddustprop,dustproppred,pxyzu,dens,metrics,ics,&
                          filfac,filfacpred,mprev,filfacprev,aprmassoftype,isionised,&
                          epot_sinksink,fxyz_ptmass_tree
 use options,        only:avdecayconst,alpha,ieos,alphamax
 use deriv,          only:derivs
 use timestep,       only:dterr,bignumber,tolv,C_force
 use mpiutils,       only:reduceall_mpi,bcast_mpi
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass, &
                          dsdt_ptmass,fsink_old,ibin_wake,dptmass,sf_ptmass, &
                          pxyzu_ptmass,metrics_ptmass
 use part,           only:n_group,n_ingroup,n_sing,gtgrad,group_info,bin_info,nmatrix
 use io_summary,     only:summary_printout,summary_variable,iosumtvi,iowake, &
                          iosumflrp,iosumflrps,iosumflrc
 use boundary_dyn,   only:dynamic_bdy,update_xyzminmax
 use timestep,       only:dtmax,dtmax_ifactor,dtdiff
 use timestep_ind,   only:get_dt,nbinmax,decrease_dtmax,dt_too_small
 use timestep_sts,   only:sts_get_dtau_next,use_sts,ibin_sts,sts_it_n
 use part,           only:ibin,ibin_old,twas,iactive,ibin_wake
 use part,           only:metricderivs,metricderivs_ptmass,fxyz_ptmass_sinksink
 use metric_tools,   only:imet_minkowski,imetric
 use cons2prim,      only:cons2primall,cons2primall_sink
 use extern_gr,      only:get_grforce_all
 use cooling,        only:ufloor,cooling_in_step,Tfloor
 use cooling_radapprox,only:radcool_evolve_ui
 use timing,         only:increment_timer,get_timings,itimer_substep
 use growth,         only:check_dustprop
 use options,        only:use_porosity,icooling
 use porosity,       only:get_filfac
 use damping,        only:idamp
 use cons2primsolver, only:conservative2primitive,primitive2conservative
 use eos,             only:equationofstate
 use substepping,     only:substep,substep_gr, &
                           substep_sph_gr,substep_sph,combine_forces_gr
 use ptmass,         only:get_accel_sink_sink,get_accel_sink_gas,ptmass_kick

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
 real               :: dtsinksink
 real               :: fonrmax,poti,dtphi2
 real               :: fext_gas(4,npart)
 integer            :: merge_ij(nptmass)
 integer            :: merge_n
 real(kind=4)       :: t1,t2,tcpu1,tcpu2
 real               :: pxi,pyi,pzi,p2i,p2mean
 real               :: dtsph_next,dti,time_now
 logical, parameter :: allow_waking = .true.
 integer, parameter :: maxits = 30
 logical            :: converged,store_itype

!
! set initial quantities
!
 fext_gas = 0.
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
 !$omp shared(dustprop,ddustprop,dustproppred,ufloor,icooling,Tfloor) &
 !$omp shared(mprev,filfacprev,filfac,use_porosity) &
 !$omp shared(ibin,ibin_old,twas,timei) &
 !$omp firstprivate(itype) &
 !$omp private(i,hdti) &
 !$omp reduction(+:nvfloorp)
 predictor: do i=1,npart
    ! print *, "predictor, i=", i
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
          if (icooling == 9) then
             vxyzu(1:3,i) = vxyzu(1:3,i) + hdti*fxyzu(1:3,i)
             call radcool_evolve_ui(vxyzu(4,i),hdti,i,Tfloor,xyzh(4,i))
          else
             vxyzu(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
          endif
       endif

       !--floor the thermal energy if requested and required
       if (ufloor > 0. .and. icooling /= 9) then
          if (vxyzu(4,i) < ufloor) then
             vxyzu(4,i) = ufloor
             nvfloorp   = nvfloorp + 1
          endif
       endif
       if (use_porosity) then
          mprev(i) = dustprop(1,i)
          filfacprev(i) = filfac(i)
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
 !$omp end parallel do

 !
 ! 1st ptmass kick (sink-gas)
 !
 if (use_sinktree .and. nptmass>0) then
    if (id==master) then
       call ptmass_kick(nptmass,hdtsph,vxyz_ptmass,fxyz_ptmass_tree,xyzmh_ptmass,dsdt_ptmass,.true.)
    endif
    call bcast_mpi(vxyz_ptmass(:,1:nptmass))
 endif

 if (use_dustgrowth) then
    if (use_porosity) then
       call get_filfac(npart,xyzh,mprev,filfac,dustprop,hdti)
    endif
    call check_dustprop(npart,dustprop,filfac,mprev,filfacprev)
 endif

!----------------------------------------------------------------------
! substepping with external and sink particle forces, using dtextforce
! accretion onto sinks/potentials also happens during substepping
!----------------------------------------------------------------------
 call get_timings(t1,tcpu1)
 if (gr) then
    call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
    call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,fext,dtextforce,dens=dens)
    ! first calculate all the force arrays on sink particles
    if (nptmass > 0) then

       call cons2primall_sink(nptmass,xyzmh_ptmass,metrics_ptmass,pxyzu_ptmass,vxyz_ptmass)
       call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass_sinksink,epot_sinksink,dtsinksink,&
                            iexternalforce,timei,merge_ij,merge_n,dsdt_ptmass)
       call get_grforce_all(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass,&
                            vxyz_ptmass,fxyz_ptmass,dtextforce,use_sink=.true.)
       do i=1,nptmass
          fxyz_ptmass(1:3,i) = fxyz_ptmass(1:3,i) + fxyz_ptmass_sinksink(1:3,i)
       enddo
       do i=1,npart
          call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass, &
                                  fext(1,i),fext(2,i),fext(3,i),poti,pmassi,fxyz_ptmass,&
                                  dsdt_ptmass,fonrmax,dtphi2,bin_info)
       enddo
    endif

    if ((iexternalforce > 0 .and. imetric /= imet_minkowski) .or. idamp > 0 .or. nptmass > 0 .or. &
        (nptmass > 0 .and. imetric == imet_minkowski)) then

       ! for now use the minimum of the two timesteps as dtextforce
       dtextforce = min(dtextforce, C_force*dtsinksink, C_force*sqrt(dtphi2))
       call substep_gr(npart,nptmass,ntypes,dtsph,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,t,&
                       xyzmh_ptmass,vxyz_ptmass,pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass,fxyz_ptmass)
    else
       call substep_sph_gr(dtsph,npart,xyzh,vxyzu,dens,pxyzu,metrics)
    endif
 else
    if (nptmass > 0 .or. iexternalforce > 0 .or. h2chemistry .or. cooling_in_step .or. idamp > 0) then
       call substep(npart,ntypes,nptmass,dtsph,dtextforce,t,xyzh,vxyzu,&
                    fext,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_tree,dsdt_ptmass,&
                    dptmass,sf_ptmass,fsink_old,nbinmax,ibin_wake,gtgrad, &
                    group_info,bin_info,nmatrix,n_group,n_ingroup,n_sing,isionised)
    else
       call substep_sph(dtsph,npart,xyzh,vxyzu)
    endif
 endif
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_substep,t2-t1,tcpu2-tcpu1)

 timei = timei + dtsph
 nvfloorps  = 0



!----------------------------------------------------
! interpolation of SPH quantities needed in the SPH
! force evaluations, using dtsph
!----------------------------------------------------
!$omp parallel do default(none) schedule(guided,1) &
!$omp shared(maxp,maxphase,maxalpha) &
!$omp shared(xyzh,vxyzu,vpred,fxyzu,divcurlv,npart,store_itype) &
!$omp shared(pxyzu,ppred,apr_level,aprmassoftype) &
!$omp shared(Bevol,dBevol,Bpred,dtsph,massoftype,iphase) &
!$omp shared(dustevol,ddustprop,dustprop,dustproppred,dustfrac,ddustevol,dustpred,use_dustfrac) &
!$omp shared(filfac,filfacpred,use_porosity) &
!$omp shared(alphaind,ieos,alphamax,ialphaloc) &
!$omp shared(eos_vars,ufloor,icooling,Tfloor) &
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
          if (use_apr) then
             pmassi = aprmassoftype(itype,apr_level(i))
          else
             pmassi = massoftype(itype)
          endif
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
          if (icooling == 9) then
             vpred(1:3,i) = vxyzu(1:3,i) + hdti*fxyzu(1:3,i)
             call radcool_evolve_ui(vxyzu(4,i),hdti,i,Tfloor,xyzh(4,i),vpred(4,i))
          else
             vpred(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
          endif
       endif

       !--floor the thermal energy if requested and required
       if (ufloor > 0.) then
          if (vpred(4,i) < ufloor) then
             vpred(4,i) = ufloor
             nvfloorps  = nvfloorps + 1
          endif
       endif
       if (use_porosity) filfacpred(i) = filfac(i)
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
 if (use_dustgrowth) then
    if (use_porosity) then
       call get_filfac(npart,xyzh,dustprop(1,:),filfacpred,dustproppred,hdti)
    endif
    call check_dustprop(npart,dustproppred(:,:),filfacpred,dustprop(1,:),filfac)
 endif
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
                dustpred,ddustevol,filfacpred,dustfrac,eos_vars,timei,dtsph,dtnew,&
                ppred,dens,metrics,apr_level)
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
    pmassi  = massoftype(igas) ! this does not appear to be used below
    ntypes  = get_ntypes(npartoftypetot)
    store_itype = (maxphase==maxp .and. ntypes > 1)
!$omp parallel default(none) &
!$omp shared(xyzh,vxyzu,vpred,fxyzu,npart,hdtsph,store_itype) &
!$omp shared(pxyzu,ppred) &
!$omp shared(Bevol,dBevol,iphase,its) &
!$omp shared(dustevol,ddustevol,use_dustfrac) &
!$omp shared(dustprop,ddustprop,dustproppred) &
!$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass,massoftype) &
!$omp shared(dtsph,ieos,ufloor,icooling,Tfloor) &
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
                   if (icooling == 9) then
                      vxyzu(1:3,i) = vxyzu(1:3,i) + dti*fxyzu(1:3,i)
                      call radcool_evolve_ui(vxyzu(4,i),dti,i,Tfloor,xyzh(4,i))
                   else
                      vxyzu(:,i) = vxyzu(:,i) + dti*fxyzu(:,i)
                   endif
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
                if (icooling == 9) then
                   vxyzu(1:3,i) = vxyzu(1:3,i) + hdti*fxyzu(1:3,i)
                   call radcool_evolve_ui(vxyzu(4,i),hdti,i,Tfloor,xyzh(4,i))
                else
                   vxyzu(:,i) = vxyzu(:,i) + hdti*fxyzu(:,i)
                endif
             endif

             !--floor the thermal energy if requested and required
             if (ufloor > 0.) then
                if (vxyzu(4,i) < ufloor .and. icooling /= 9) then
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
                if (maxvxyzu >= 4) then
                   if (icooling == 9) then
                      call radcool_evolve_ui(vxyzu(4,i),hdtsph,i,Tfloor,xyzh(4,i),eni)
                   else
                      eni = vxyzu(4,i) + hdtsph*fxyzu(4,i)
                   endif
                endif
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
    if (use_dustgrowth) then
       if (use_porosity) then
          call get_filfac(npart,xyzh,mprev,filfac,dustprop,dtsph)
       endif
       call check_dustprop(npart,dustprop,filfac,mprev,filfacprev)
    endif

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
!$omp shared(filfac,filfacpred,use_porosity) &
!$omp shared(rad,drad,radpred,icooling,Tfloor,xyzh) &
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
                if (use_porosity) filfacpred(i) = filfac(i)
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
             if (use_porosity) filfacpred(i) = filfac(i)
             if (mhd)          Bpred(:,i)  = Bevol(:,i)
             if (use_dustfrac) dustpred(:,i) = dustevol(:,i)
             if (do_radiation) radpred(:,i) = rad(:,i)
             !
             ! shift v back to the half step
             !
             if (gr) then
                pxyzu(:,i) = pxyzu(:,i) - hdtsph*fxyzu(:,i)
             else
                if (icooling == 9) then
                   call radcool_evolve_ui(vxyzu(4,i),-hdtsph,i,Tfloor,xyzh(4,i))
                   vxyzu(1:3,i) = vxyzu(1:3,i) - hdtsph*fxyzu(1:3,i)
                else
                   vxyzu(:,i) = vxyzu(:,i) - hdtsph*fxyzu(:,i)
                endif
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

       if (use_dustgrowth) then
          if (use_porosity) then
             call get_filfac(npart,xyzh,mprev,filfac,dustprop,dtsph)
          endif
          call check_dustprop(npart,dustprop,filfac,mprev,filfacprev)
       endif
!
!   get new force using updated velocity: no need to recalculate density etc.
!
       if (gr) vpred = vxyzu ! Need primitive utherm as a guess in cons2prim
       call derivs(2,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,divcurlB, &
                     Bpred,dBevol,radpred,drad,radprop,dustproppred,ddustprop,dustpred,ddustevol,filfacpred,&
                     dustfrac,eos_vars,timei,dtsph,dtnew,ppred,dens,metrics,apr_level)
       if (gr) vxyzu = vpred ! May need primitive variables elsewhere?
       if (do_radiation .and. implicit_radiation) then
          rad = radpred
          vxyzu(4,1:npart) = vpred(4,1:npart)
       endif
    endif
    if (icooling == 9 .and. iverbose >=2) then
       print *, "end of iteration", maxval(vpred(4,:)), minval(vpred(4,:))
       print *, "end of iteration", maxval(vxyzu(4,:)), minval(vxyzu(4,:))
    endif
 enddo iterations

 !
 ! 2nd ptmass kick (no need to predict vel ptmass as they are not coupled to any vel dep force)
 !
 if (use_sinktree .and. nptmass>0) then
    if (id==master) then
       call ptmass_kick(nptmass,hdtsph,vxyz_ptmass,fxyz_ptmass_tree,xyzmh_ptmass,dsdt_ptmass,.true.)
    endif
    call bcast_mpi(vxyz_ptmass(:,1:nptmass))
 endif


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

 if (gr) call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)

end subroutine step

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
