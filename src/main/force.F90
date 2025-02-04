!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module forces
!
! This module is the "guts" of the code
!   Calculates force and rates of change for all particles
!
! :References:
!
!  Code paper:
!      Price et al. (2018), PASA 35, e031
!  Hydro:
!      Price (2012), J. Comp. Phys. 231, 759-794
!      Lodato & Price (2010), MNRAS 405, 1212-1226
!      Price & Federrath (2010), MNRAS 406, 1659-1674
!  MHD:
!      Tricco & Price (2012), J. Comp. Phys. 231, 7214-7236
!      Tricco, Price & Bate (2016), MNRAS 322, 326-344
!      Wurster, Price & Ayliffe (2014), MNRAS 444, 1104-1112
!      Wurster, Price & Bate (2016), MNRAS 457, 1037-1061
!  Dust:
!      Laibe & Price (2012a), MNRAS 420, 2345-2364
!      Laibe & Price (2012b), MNRAS 420, 2365-2376
!      Price & Laibe (2015), MNRAS 451, 5332-5345
!      Hutchison, Price & Laibe (2018), MNRAS 476, 2186-2198
!      Ballabio et al. (2018), MNRAS 477, 2766-2771
!      Mentiplay, Price, Pinte & Laibe (2020), MNRAS 499, 3806-3818
!      Price & Laibe (2020), MNRAS 495, 3929-3934
!  Radiation:
!      Whitehouse & Bate (2004), MNRAS 353, 1078-1094
!  GR:
!      Liptai & Price (2019), MNRAS 485, 819-842
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: boundary, cooling, dim, dust, eos, eos_shen,
!   eos_stamatellos, fastmath, growth, io, io_summary, kdtree, kernel,
!   linklist, metric_tools, mpiderivs, mpiforce, mpimemory, mpiutils,
!   nicil, omputils, options, part, physcon, ptmass, ptmass_heating,
!   radiation_utils, timestep, timestep_ind, timestep_sts, timing, units,
!   utils_gr, viscosity
!
 use dim, only:maxfsum,maxxpartveciforce,maxp,ndivcurlB,ndivcurlv,&
               maxdusttypes,maxdustsmall,do_radiation
 use mpiforce, only:cellforce,stackforce
 use linklist, only:ifirstincell
 use kdtree,   only:inodeparts,inoderange
 use part,     only:iradxi,ifluxx,ifluxy,ifluxz,ikappa,ien_type,ien_entropy,ien_etotal,ien_entropy_s

 implicit none

 integer, parameter :: maxcellcache = 1000

 public :: force, reconstruct_dv, get_drag_terms ! latter to avoid compiler warning

 !--indexing for xpartveci array
 integer, parameter ::       &
       ixi             = 1,  &
       iyi             = 2,  &
       izi             = 3,  &
       ihi             = 4,  &
       ivxi            = 5,  &
       ivyi            = 6,  &
       ivzi            = 7,  &
       ieni            = 8,  &
       iBevolxi        = 9,  &
       iBevolyi        = 10, &
       iBevolzi        = 11, &
       ipsi            = 12, &
       igradhi1        = 13, &
       igradhi2        = 14, &
       ialphai         = 15, &
       ialphaBi        = 16, &
       ivwavei         = 17, &
       irhoi           = 18, &
       irhogasi        = 19, &
       ispsoundi       = 20, &
       itempi          = 21, &
       isxxi           = 22, &
       isxyi           = 23, &
       isxzi           = 24, &
       isyyi           = 25, &
       isyzi           = 26, &
       iszzi           = 27, &
       ivisctermisoi   = 28, &
       ivisctermanisoi = 29, &
       ipri            = 30, &
       ipro2i          = 31, &
       ietaohmi        = 32, &
       ietahalli       = 33, &
       ietaambii       = 34, &
       ijcbcbxi        = 35, &
       ijcbcbyi        = 36, &
       ijcbcbzi        = 37, &
       ijcbxi          = 38, &
       ijcbyi          = 39, &
       ijcbzi          = 40, &
       idivBi          = 41, &
       icurlBxi        = 42, &
       icurlByi        = 43, &
       icurlBzi        = 44, &
       igrainmassi     = 45, &
       igraindensi     = 46, &
       ifilfaci        = 47, &
       ifxi_drag       = 48, &
       ifyi_drag       = 49, &
       ifzi_drag       = 50, &
       idti            = 51, &
       idvxdxi         = 52, &
       idvzdzi         = 60, &
 !--dust arrays initial index
       idustfraci      = 61, &
 !--dust arrays final index
       idustfraciend   = 61 + (maxdusttypes - 1), &
       itstop          = 62 + (maxdusttypes - 1), &
       itstopend       = 62 + 2*(maxdusttypes - 1), &
 !--final dust index
       lastxpvdust     = 62 + 2*(maxdusttypes - 1), &
       iradxii         = lastxpvdust + 1, &
       iradfxi         = lastxpvdust + 2, &
       iradfyi         = lastxpvdust + 3, &
       iradfzi         = lastxpvdust + 4, &
       iradkappai      = lastxpvdust + 5, &
       iradlambdai     = lastxpvdust + 6, &
       iradrbigi       = lastxpvdust + 7, &
 !--final radiation index
       lastxpvrad      = lastxpvdust + 7, &
 !--gr primitive density
       idensGRi        = lastxpvrad + 1, &
 !--gr metrics
       imetricstart    = idensGRi + 1, &
       imetricend      = imetricstart + 31

 !--indexing for fsum array
 integer, parameter ::   &
       ifxi           = 1,  &
       ifyi           = 2,  &
       ifzi           = 3,  &
       ipot           = 4,  &
       idrhodti       = 5,  &
       idudtdissi     = 6,  &
       idendtdissi    = 7,  &
       idivBsymi      = 8,  &
       idBevolxi      = 9,  &
       idBevolyi      = 10, &
       idBevolzi      = 11, &
       idivBdiffi     = 12, &
 !--dust array indexing
       ifdragxi       = 13, &
       ifdragyi       = 14, &
       ifdragzi       = 15, &
       iddustevoli    = 16, &
       iddustevoliend = 16 +   (maxdustsmall-1), &
       idudtdusti     = 17 +   (maxdustsmall-1), &
       idudtdustiend  = 17 + 2*(maxdustsmall-1), &
       ideltavxi      = 18 + 2*(maxdustsmall-1), &
       ideltavxiend   = 18 + 3*(maxdustsmall-1), &
       ideltavyi      = 19 + 3*(maxdustsmall-1), &
       ideltavyiend   = 19 + 4*(maxdustsmall-1), &
       ideltavzi      = 20 + 4*(maxdustsmall-1), &
       ideltavziend   = 20 + 5*(maxdustsmall-1), &
       idvix          = 21 + 5*(maxdustsmall-1), &
       idviy          = 22 + 5*(maxdustsmall-1), &
       idviz          = 23 + 5*(maxdustsmall-1), &
       idensgasi      = 24 + 5*(maxdustsmall-1), &
       icsi           = 25 + 5*(maxdustsmall-1), &
       idradi         = 25 + 5*(maxdustsmall-1) + 1

 private

contains

!----------------------------------------------------------------
!+
!  compute all forces and rates of change on the particles
!+
!----------------------------------------------------------------
subroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&
                 rad,drad,radprop,dustprop,dustgasprop,dustfrac,ddustevol,fext,fxyz_drag,&
                 ipart_rhomax,dt,stressmax,eos_vars,dens,metrics,apr_level)

 use dim,          only:maxvxyzu,maxneigh,mhd,mhd_nonideal,lightcurve,mpi,use_dust,use_apr
 use io,           only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs
 use linklist,     only:ncells,get_neighbour_list,get_hmaxcell,get_cell_location,listneigh
 use options,      only:iresistive_heating
 use part,         only:rhoh,dhdrho,rhoanddhdrho,alphaind,iactive,gradh,&
                        hrho,iphase,igas,maxgradh,dvdx,eta_nimhd,deltav,poten,iamtype,&
                        dragreg,filfac,fxyz_dragold,aprmassoftype
 use timestep,     only:dtcourant,dtforce,dtrad,bignumber,dtdiff
 use io_summary,   only:summary_variable, &
                        iosumdtf,iosumdtd,iosumdtv,iosumdtc,iosumdto,iosumdth,iosumdta, &
                        iosumdgs,iosumdge,iosumdgr,iosumdtfng,iosumdtdd,iosumdte,iosumdtB,iosumdense
#ifdef FINVSQRT
 use fastmath,     only:finvsqrt
#endif
 use physcon,      only:pi
 use viscosity,    only:irealvisc,shearfunc,dt_viscosity
#ifdef IND_TIMESTEPS
 use timestep_ind, only:nbinmax,ibinnow,get_newbin
 use timestep_sts, only:nbinmaxsts
 use timestep,     only:nsteps,time
#else
 use timestep,     only:C_cour,C_force
#endif
 use part,         only:divBsymm,isdead_or_accreted,ngradh,gravity,ibin_wake
 use mpiutils,     only:reduce_mpi,reduceall_mpi,reduceloc_mpi,bcast_mpi
#ifdef GRAVITY
 use kernel,       only:kernel_softening
 use kdtree,       only:expand_fgrav_in_taylor_series
 use linklist,     only:get_distance_from_centre_of_mass
 use part,         only:xyzmh_ptmass,nptmass,massoftype,maxphase,is_accretable,ihacc
 use ptmass,       only:icreate_sinks,rho_crit,r_crit2,h_acc
 use units,        only:unit_density
#endif
#ifdef DUST
 use kernel,       only:wkern_drag,cnormk_drag
#endif
 use dust,         only:drag_implicit
 use nicil,        only:nimhd_get_jcbcb
 use mpiderivs,    only:send_cell,recv_cells,check_send_finished,init_cell_exchange,&
                        finish_cell_exchange,recv_while_wait,reset_cell_counters,cell_counters
 use mpimemory,    only:reserve_stack,reset_stacks,get_cell,write_cell
 use mpimemory,    only:stack_remote  => force_stack_1
 use mpimemory,    only:stack_waiting => force_stack_2
 use io_summary,   only:iosumdtr
 use timing,       only:increment_timer,get_timings,itimer_force_local,itimer_force_remote
 use omputils,     only:omp_thread_num,omp_num_threads

 integer,      intent(in)    :: icall,npart
 integer(kind=1), intent(in) :: apr_level(:)
 real,         intent(in)    :: xyzh(:,:)
 real,         intent(inout) :: vxyzu(:,:)
 real,         intent(in)    :: dustfrac(:,:)
 real,         intent(in)    :: dustprop(:,:)
 real,         intent(inout) :: dustgasprop(:,:)
 real,         intent(in)    :: fext(:,:)
 real,         intent(inout) :: fxyz_drag(:,:)
 real,         intent(in)    :: eos_vars(:,:)
 real,         intent(out)   :: fxyzu(:,:),ddustevol(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real,         intent(out)   :: dBevol(:,:)
 real(kind=4), intent(inout) :: divcurlv(:,:)
 real(kind=4), intent(in)    :: divcurlB(:,:)
 real,         intent(in)    :: dt,stressmax
 integer,      intent(out)   :: ipart_rhomax ! test this particle for point mass creation
 real,         intent(in)    :: rad(:,:)
 real,         intent(out)   :: drad(:,:)
 real,         intent(inout) :: radprop(:,:)
 real,         intent(in)    :: dens(:), metrics(:,:,:,:)

 real, save :: xyzcache(maxcellcache,4)
!$omp threadprivate(xyzcache)
 integer :: i,icell,nneigh
 integer :: nstokes,nsuper,ndrag,ndustres,ndense
 real    :: dtmini,dtohm,dthall,dtambi,dtvisc
 real    :: dustresfacmean,dustresfacmax
#ifdef GRAVITY
 real    :: potensoft0,dum,dx,dy,dz,fxi,fyi,fzi,poti,epoti
 real    :: rhomax,rhomax_thread
 logical :: use_part
 integer :: ipart_rhomax_thread,j,id_rhomax
 real    :: hi,pmassi,rhoi
 logical :: iactivei,iamdusti
 integer :: iamtypei
#endif
#ifdef DUST
 real                   :: frac_stokes,frac_super
#endif
 logical :: realviscosity,useresistiveheat
#ifndef IND_TIMESTEPS
 real    :: dtmaxi,minglobdt
#else
 integer :: nbinmaxnew,nbinmaxstsnew,ncheckbin
 integer :: ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd
 integer :: ndtvisc,ndtohm,ndthall,ndtambi,ndtdust,ndtrad,ndtclean
 real    :: dtitmp,dtrat,dtmaxi
 real    :: dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean,dtcleanfacmean
 real    :: dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax ,dtcleanfacmax
 real    :: dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean,dtdustfacmean ,dtradfacmean
 real    :: dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax, dtdustfacmax  ,dtradfacmax
#endif
 integer(kind=1)           :: ibinnow_m1

 type(cellforce)           :: cell,xsendbuf,xrecvbuf(nprocs)
 integer                   :: mpitype
 logical                   :: remote_export(nprocs),do_export,idone(nprocs),thread_complete(omp_num_threads)
 integer                   :: irequestsend(nprocs),irequestrecv(nprocs)
 integer                   :: ncomplete_mpi

 real(kind=4)              :: t1,t2,tcpu1,tcpu2

#ifdef IND_TIMESTEPS
 nbinmaxnew      = 0
 nbinmaxstsnew   = 0
 ndtforce        = 0
 ndtforceng      = 0
 ndtcool         = 0
 ndtdrag         = 0
 ndtdragd        = 0
 ndtdust         = 0
 ncheckbin       = 0
 ndtvisc         = 0
 ndtohm          = 0
 ndthall         = 0
 ndtambi         = 0
 ndtrad          = 0
 ndtclean        = 0
 dtfrcfacmean    = 0.0
 dtfrcngfacmean  = 0.0
 dtdragfacmean   = 0.0
 dtdragdfacmean  = 0.0
 dtcoolfacmean   = 0.0
 dtviscfacmean   = 0.0
 dtohmfacmean    = 0.0
 dthallfacmean   = 0.0
 dtambifacmean   = 0.0
 dtdustfacmean   = 0.0
 dtfrcfacmax     = 0.0
 dtfrcngfacmax   = 0.0
 dtdragfacmax    = 0.0
 dtdragdfacmax   = 0.0
 dtcoolfacmax    = 0.0
 dtviscfacmax    = 0.0
 dtohmfacmax     = 0.0
 dthallfacmax    = 0.0
 dtambifacmax    = 0.0
 dtdustfacmax    = 0.0
 dtradfacmean    = 0.0
 dtradfacmax     = 0.0
 dtcleanfacmean  = 0.0
 dtcleanfacmax   = 0.0
 ibinnow_m1      = ibinnow - 1_1
#else
 ibinnow_m1      = 0
#endif
 dustresfacmean  = 0.0
 dustresfacmax   = 0.0
 dtmaxi          = 0.
 dtcourant       = bignumber
 dtforce         = bignumber
 dtvisc          = bignumber
 dtmini          = bignumber
 dtohm           = bignumber
 dthall          = bignumber
 dtambi          = bignumber
 dtrad           = bignumber

 if (iverbose >= 3 .and. id==master) write(iprint,*) 'forces: cell cache =',maxcellcache

 realviscosity    = (irealvisc > 0)
 useresistiveheat = (iresistive_heating > 0)
 if (ndivcurlv < 1) call fatal('force','divv not stored but it needs to be')

 !--dust/gas stuff
 ndrag         = 0
 nstokes       = 0
 nsuper        = 0
 ndustres      = 0
 !store the force vector needed for implicit drag
 if (use_dust .and. drag_implicit) then
    !$omp parallel do default(none) shared(fext,npart,fxyz_drag,fxyz_dragold) private(i)
    do i=1,npart
       fxyz_dragold(1,i) = fxyz_drag(1,i)
       fxyz_dragold(2,i) = fxyz_drag(2,i)
       fxyz_dragold(3,i) = fxyz_drag(3,i)
    enddo
    !$omp end parallel do
 endif

 ! sink particle creation
 ndense        = 0
 ipart_rhomax  = 0
#ifdef GRAVITY
 rhomax        = 0.
#endif

 if (mpi) then
    call reset_stacks
    call reset_cell_counters(cell_counters)
 endif

!
!-- verification for non-ideal MHD
!
 if (mhd_nonideal .and. ndivcurlB < 4) call fatal('force','non-ideal MHD needs curl B stored, but ndivcurlB < 4')
!
!--check that compiled options are compatible with this routine
!
 if (maxgradh /= maxp) call fatal('force','need storage of gradh (maxgradh=maxp)')

!$omp parallel default(none) &
!$omp shared(maxp) &
!$omp shared(ncells,ifirstincell) &
!$omp shared(xyzh) &
!$omp shared(dustprop) &
!$omp shared(dragreg) &
!$omp shared(filfac) &
!$omp shared(dustgasprop) &
!$omp shared(fxyz_drag) &
!$omp shared(fxyz_dragold) &
!$omp shared(fext) &
!$omp shared(vxyzu) &
!$omp shared(fxyzu) &
!$omp shared(divcurlv) &
!$omp shared(iphase) &
!$omp shared(dvdx) &
!$omp shared(gradh) &
!$omp shared(divcurlb) &
!$omp shared(bevol) &
!$omp shared(rad,radprop,drad) &
!$omp shared(eta_nimhd) &
!$omp shared(alphaind) &
!$omp shared(stressmax) &
!$omp shared(divBsymm) &
!$omp shared(dBevol) &
!$omp shared(eos_vars) &
!$omp shared(dt) &
!$omp shared(nprocs,icall) &
!$omp shared(poten) &
!$omp private(icell,i) &
!$omp private(cell) &
!$omp private(remote_export) &
!$omp private(idone) &
!$omp private(nneigh) &
!$omp private(mpitype) &
!$omp shared(dens) &
!$omp shared(metrics) &
!$omp shared(apr_level,aprmassoftype) &
#ifdef GRAVITY
!$omp shared(massoftype,npart,maxphase) &
!$omp private(hi,pmassi,rhoi) &
!$omp private(iactivei,iamdusti,iamtypei) &
!$omp private(dx,dy,dz,poti,fxi,fyi,fzi,potensoft0,dum,epoti) &
!$omp shared(xyzmh_ptmass,nptmass) &
!$omp shared(rhomax,ipart_rhomax,icreate_sinks,rho_crit,r_crit2,h_acc) &
!$omp private(rhomax_thread,ipart_rhomax_thread,use_part,j) &
#endif
!$omp shared(id) &
!$omp private(do_export) &
!$omp private(irequestrecv) &
!$omp private(irequestsend) &
!$omp private(xrecvbuf) &
!$omp private(xsendbuf) &
!$omp shared(cell_counters) &
!$omp shared(thread_complete) &
!$omp shared(ncomplete_mpi) &
!$omp shared(stack_remote) &
!$omp shared(stack_waiting) &
#ifdef IND_TIMESTEPS
!$omp shared(nbinmax,nbinmaxsts) &
!$omp private(dtitmp,dtrat) &
!$omp reduction(+:ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd,ncheckbin,ndtvisc,ndtrad,ndtclean) &
!$omp reduction(+:ndtohm,ndthall,ndtambi,ndtdust,dtohmfacmean,dthallfacmean,dtambifacmean,dtdustfacmean) &
!$omp reduction(+:dtfrcfacmean,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean,dtviscfacmean) &
!$omp reduction(+:dtradfacmean,dtcleanfacmean) &
!$omp reduction(max:dtohmfacmax,dthallfacmax,dtambifacmax,dtdustfacmax,dtradfacmax,dtcleanfacmax) &
!$omp reduction(max:dtfrcfacmax,dtfrcngfacmax,dtdragfacmax,dtdragdfacmax,dtcoolfacmax,dtviscfacmax) &
!$omp reduction(max:nbinmaxnew,nbinmaxstsnew) &
#endif
!$omp reduction(+:ndustres,dustresfacmean,ndense) &
!$omp reduction(min:dtrad) &
!$omp reduction(min:dtohm,dthall,dtambi,dtdiff) &
!$omp reduction(min:dtcourant,dtforce,dtvisc) &
!$omp reduction(max:dtmaxi,dustresfacmax) &
!$omp reduction(min:dtmini) &
!$omp shared(dustfrac) &
!$omp shared(ddustevol) &
!$omp shared(deltav) &
!$omp shared(ibin_wake,ibinnow_m1) &
!$omp shared(t1) &
!$omp shared(t2) &
!$omp shared(tcpu1) &
!$omp shared(tcpu2)

 call init_cell_exchange(xrecvbuf,irequestrecv,thread_complete,ncomplete_mpi,mpitype)

 !$omp master
 call get_timings(t1,tcpu1)
 !$omp end master

 !--initialise send requests to 0
 irequestsend = 0

 !$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)

    !--skip empty cells AND inactive cells
    if (i <= 0) cycle over_cells

    cell%icell = icell

    call start_cell(cell,iphase,xyzh,vxyzu,gradh,divcurlv,divcurlB,dvdx,Bevol, &
                    dustfrac,dustprop,fxyz_dragold,eta_nimhd,eos_vars,alphaind,stressmax,&
                    rad,radprop,dens,metrics,apr_level,dt)
    if (cell%npcell == 0) cycle over_cells

    call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)

    !--get the neighbour list and fill the cell cache
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache, &
                           getj=.true.,f=cell%fgrav,remote_export=remote_export)

    cell%owner = id
    do_export = any(remote_export)

    if (mpi) then
       call recv_cells(stack_remote,xrecvbuf,irequestrecv,cell_counters)
       if (do_export) then
          if (stack_waiting%n > 0) then
             !--wait for broadcast to complete, continue to receive whilst doing so
             idone(:) = .false.
             do while(.not.all(idone))
                call check_send_finished(irequestsend,idone)
                call recv_cells(stack_remote,xrecvbuf,irequestrecv,cell_counters)
             enddo
          endif
          call reserve_stack(stack_waiting,cell%waiting_index)
          call send_cell(cell,remote_export,irequestsend,xsendbuf,cell_counters,mpitype)  ! send to remote
       endif
    endif

    call compute_cell(cell,listneigh,nneigh,Bevol,xyzh,vxyzu,fxyzu, &
                      iphase,divcurlv,divcurlB,alphaind,eta_nimhd,eos_vars, &
                      dustfrac,dustprop,fxyz_dragold,gradh,ibinnow_m1,ibin_wake,stressmax,xyzcache,&
                      rad,radprop,dens,metrics,apr_level,dt)

    if (do_export) then
       call write_cell(stack_waiting,cell)
    else
       call finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,dvdx,&
                             divBsymm,divcurlv,dBevol,ddustevol,deltav,dustgasprop,fxyz_drag,fext,dragreg,&
                             filfac,dtcourant,dtforce,dtvisc,dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi, &
#ifdef IND_TIMESTEPS
                             nbinmaxnew,nbinmaxstsnew,ncheckbin, &
                             ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd, &
                             ndtvisc,ndtohm,ndthall,ndtambi,ndtdust,ndtrad,ndtclean, &
                             dtitmp,dtrat, &
                             dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean, &
                             dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax, &
                             dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean ,dtdustfacmean, &
                             dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax  ,dtdustfacmax, &
                             dtradfacmean ,dtcleanfacmean, &
                             dtradfacmax  ,dtcleanfacmax, &
#endif
                             ndustres,dustresfacmax,dustresfacmean, &
                             rad,drad,radprop,dtrad)
    endif

 enddo over_cells
 !$omp enddo

 if (stack_waiting%n > 0) then
    idone(:) = .false.
    do while(.not.all(idone))
       call check_send_finished(irequestsend,idone)
       call recv_cells(stack_remote,xrecvbuf,irequestrecv,cell_counters)
    enddo
 endif

 if (mpi) then
    call recv_while_wait(stack_remote,xrecvbuf,irequestrecv,&
         irequestsend,thread_complete,cell_counters,ncomplete_mpi)
    call reset_cell_counters(cell_counters)
 endif

 !$omp master
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_force_local,t2-t1,tcpu2-tcpu1)
 call get_timings(t1,tcpu1)
 !$omp end master
 !$omp barrier

 igot_remote: if (mpi .and. stack_remote%n > 0) then
    !$omp do schedule(runtime)
    over_remote: do i = 1,stack_remote%n
       cell = get_cell(stack_remote,i)

       call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,maxcellcache, &
                               getj=.true.,f=cell%fgrav,&
                               cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti)

       call compute_cell(cell,listneigh,nneigh,Bevol,xyzh,vxyzu,fxyzu, &
                         iphase,divcurlv,divcurlB,alphaind,eta_nimhd,eos_vars, &
                         dustfrac,dustprop,fxyz_dragold,gradh,ibinnow_m1,ibin_wake,stressmax,xyzcache,&
                         rad,radprop,dens,metrics,apr_level,dt)

       remote_export = .false.
       remote_export(cell%owner+1) = .true. ! use remote_export array to send back to the owner

       idone(:) = .false.
       do while(.not.all(idone))
          call check_send_finished(irequestsend,idone)
          call recv_cells(stack_waiting,xrecvbuf,irequestrecv,cell_counters)
       enddo

       call send_cell(cell,remote_export,irequestsend,xsendbuf,cell_counters,mpitype) ! send the cell back to owner

    enddo over_remote
    !$omp enddo

    !$omp master
    stack_remote%n = 0
    !$omp end master

    idone(:) = .false.
    do while(.not.all(idone))
       call check_send_finished(irequestsend,idone)
       call recv_cells(stack_waiting,xrecvbuf,irequestrecv,cell_counters)
    enddo

 endif igot_remote

 if (mpi) call recv_while_wait(stack_waiting,xrecvbuf,irequestrecv,&
          irequestsend,thread_complete,cell_counters,ncomplete_mpi)

 iam_waiting: if (mpi .and. stack_waiting%n > 0) then
    !$omp do schedule(runtime)
    over_waiting: do i = 1, stack_waiting%n
       cell = get_cell(stack_waiting,i)

       call finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,dvdx, &
                                          divBsymm,divcurlv,dBevol,ddustevol,deltav,dustgasprop,fxyz_drag,fext,dragreg, &
                                          filfac,dtcourant,dtforce,dtvisc,dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi, &
#ifdef IND_TIMESTEPS
                                          nbinmaxnew,nbinmaxstsnew,ncheckbin, &
                                          ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd, &
                                          ndtvisc,ndtohm,ndthall,ndtambi,ndtdust,ndtrad,ndtclean, &
                                          dtitmp,dtrat, &
                                          dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean, &
                                          dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax, &
                                          dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean ,dtdustfacmean, &
                                          dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax  ,dtdustfacmax, &
                                          dtradfacmean ,dtcleanfacmean, &
                                          dtradfacmax  ,dtcleanfacmax, &
#endif
                                          ndustres,dustresfacmax,dustresfacmean, &
                                          rad,drad,radprop,dtrad)

    enddo over_waiting
    !$omp enddo

    stack_waiting%n = 0

 endif iam_waiting

 call finish_cell_exchange(irequestrecv,xsendbuf,mpitype)

!$omp master
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_force_remote,t2-t1,tcpu2-tcpu1)
!$omp end master

#ifdef GRAVITY
 if (icreate_sinks > 0) then
    rhomax_thread = 0.
    ipart_rhomax_thread = 0
    !$omp do schedule(runtime)
    over_parts: do i=1,npart
       hi = xyzh(4,i)
#ifdef IND_TIMESTEPS
       if (iactive(iphase(i)) .and..not.isdead_or_accreted(hi)) then
#else
       if (.not.isdead_or_accreted(hi)) then
#endif
          if (maxphase==maxp) then
             iamtypei = iamtype(iphase(i))
             if (.not.is_accretable(iamtypei)) cycle over_parts
          else
             iamtypei = igas
          endif
          if (use_apr) then
             pmassi = aprmassoftype(iamtypei,apr_level(i))
          else
             pmassi = massoftype(iamtypei)
          endif
          rhoi = rhoh(hi,pmassi)
          if (rhoi > rho_crit) then
             if (rhoi > rhomax_thread) then
                !
                !--find the maximum density on particles outside the
                !  allowed minimum distance from other sink particles
                !
                use_part = .true.
                over_ptmass: do j=1,nptmass
                   if (icreate_sinks==2 .and. xyzmh_ptmass(ihacc,j)<h_acc ) cycle
                   if (xyzmh_ptmass(4,j) > 0. .and.       &
                       (xyzh(1,i) - xyzmh_ptmass(1,j))**2 &
                     + (xyzh(2,i) - xyzmh_ptmass(2,j))**2 &
                     + (xyzh(3,i) - xyzmh_ptmass(3,j))**2 < r_crit2) then
                      ndense   = ndense + 1
                      use_part = .false.
                      exit over_ptmass
                   endif
                enddo over_ptmass
                if (use_part) then
                   rhomax_thread = rhoi
                   ipart_rhomax_thread = i
                endif
             endif
          endif
       endif
    enddo over_parts
    !$omp enddo
    if (rhomax_thread > rho_crit) then
       !$omp critical(rhomaxadd)
       if (rhomax_thread > rhomax) then
          rhomax = rhomax_thread
          ipart_rhomax = ipart_rhomax_thread
       endif
       !$omp end critical(rhomaxadd)
    endif
 endif
#endif
!$omp end parallel

#ifdef IND_TIMESTEPS
 ! check for nbinmaxnew = 0, can happen if all particles
 ! are dead/inactive, e.g. after sink creation or if all
 ! have moved to a higher ibin; the following step on the
 ! higher ibin will yeild non-zero and modify nbinmax
 ! appropriately
 if (ncheckbin==0) then
    nbinmaxnew    = nbinmax
    nbinmaxstsnew = nbinmaxsts
 endif
#endif

#ifdef GRAVITY
 if (reduceall_mpi('max',ipart_rhomax) > 0) then
    call reduceloc_mpi('max',rhomax,id_rhomax)
    if (id /= id_rhomax) ipart_rhomax = -1
 endif
 if (icreate_sinks > 0 .and. ipart_rhomax > 0 .and. iverbose>=1) then
    print*,' got rhomax = ',rhomax*unit_density,' on particle ',ipart_rhomax !,rhoh(xyzh(4,ipart_rhomax))
 endif
 ndense = int(reduce_mpi('+',ndense))
 if (ndense > 0) call summary_variable('dense',iosumdense,ndense,0.)
#endif

#ifdef DUST
 ndrag = int(reduceall_mpi('+',ndrag))
 if (ndrag > 0) then
    nstokes = int(reduce_mpi('+',nstokes))
    nsuper =  int(reduce_mpi('+',nsuper))
    frac_stokes = nstokes/real(ndrag)
    frac_super  = nsuper/real(ndrag)
    if (iverbose >= 1 .and. id==master) then
       if (nstokes > 0) call warning('force','using Stokes drag regime',var='%Stokes',val=100.*frac_stokes)
       if (nsuper > 0)  call warning('force','supersonic Epstein regime',val=100.*frac_super,var='%super')
    endif
    if (nstokes > 0) call summary_variable('dust',iosumdgs,nstokes,100.*frac_stokes)
    if (nsuper  > 0) call summary_variable('dust',iosumdge,nsuper ,100.*frac_super )
 else
    frac_stokes = 0.
    frac_super  = 0.
 endif
 if (ndustres > 0) call summary_variable('dust',iosumdgr,ndustres  ,dustresfacmean/real(ndustres), dustresfacmax )
#endif

#ifdef IND_TIMESTEPS
 nbinmax    = int(reduceall_mpi('max',nbinmaxnew),kind=1)
 nbinmaxsts = int(reduceall_mpi('max',nbinmaxstsnew),kind=1)
 ndtforce   = int(reduce_mpi('+',ndtforce))
 ndtforceng = int(reduce_mpi('+',ndtforceng))
 ndtcool    = int(reduce_mpi('+',ndtcool))
 ndtdrag    = int(reduce_mpi('+',ndtdrag))
 ndtdragd   = int(reduce_mpi('+',ndtdragd))
 ndtdust    = int(reduce_mpi('+',ndtdust))
 ndtrad     = int(reduce_mpi('+',ndtrad))
 ndtclean   = int(reduce_mpi('+',ndtclean))

 !  Print warning statements, if required
 if (iverbose >= 1 .and. id==master) then
    if (ndtforce   > 0) write(iprint,*) 'force controlling timestep on ',ndtforce,' gas particles'
    if (ndtforceng > 0) write(iprint,*) 'force controlling timestep on ',ndtforce,' non-gas particles'
    if (ndtcool    > 0) write(iprint,*) 'cooling controlling timestep on ',ndtcool,' particles'
    if (ndtdrag    > 0) write(iprint,*) 'drag controlling timestep on ',ndtdrag,' gas particles'
    if (ndtdragd   > 0) write(iprint,*) 'drag controlling timestep on ',ndtdrag,' dust particles'
    if (ndtdust    > 0) write(iprint,*) 'dust diffusion controlling timestep on ',ndtdust,' particles'
    if (ndtrad     > 0) write(iprint,*) 'radiation diffusion controlling timestep on ',ndtrad,' particles'
    if (ndtclean   > 0) write(iprint,*) 'B-cleaning controlling timestep on ',ndtclean,' particles'
    if (ndtvisc    > 0) then
       write(iprint,*)   'thread ',id,' WARNING: viscosity           constraining timestep on ',ndtvisc,' particles by factor ', &
                       dtviscfacmean/real(ndtvisc)
    endif
    if (mhd_nonideal) then
       if (ndtohm  > 0) &
        write(iprint,'(a,Es16.9,I8,a,I8,a,2F9.2)') 'WARNING: at (time, step) = ',time,nsteps, &
                                                   ', ohmic resistivity   constraining timestep on ',ndtohm, &
                                                   ' particles by (ave, max) factor of',dtohmfacmean/real(ndtohm),dtohmfacmax
       if (ndthall > 0) &
        write(iprint,'(a,Es16.9,I8,a,I8,a,2F9.2)') 'WARNING: at (time, step) = ',time,nsteps, &
                                                   ', Hall Effect         constraining timestep on ',ndthall, &
                                                   ' particles by (ave, max) factor of',dthallfacmean/real(ndthall),dthallfacmax
       if (ndtambi > 0) &
        write(iprint,'(a,Es16.9,I8,a,I8,a,2F9.2)') 'WARNING: at (time, step) = ',time,nsteps, &
                                                   ', ambipolar diffusion constraining timestep on ',ndtambi, &
                                                   ' particles by (ave, max) factor of',dtambifacmean/real(ndtambi),dtambifacmax
    endif
 endif
 !  Save values for summary
 if (ndtforce   > 0)  call summary_variable('dt',iosumdtf  ,ndtforce  ,dtfrcfacmean  /real(ndtforce)  ,dtfrcfacmax  )
 if (ndtforceng > 0)  call summary_variable('dt',iosumdtfng,ndtforceng,dtfrcngfacmean/real(ndtforceng),dtfrcngfacmax)
 if (ndtcool    > 0)  call summary_variable('dt',iosumdtc  ,ndtcool   ,dtcoolfacmean /real(ndtcool)   ,dtcoolfacmax )
 if (ndtdrag    > 0)  call summary_variable('dt',iosumdtd  ,ndtdrag   ,dtdragfacmean /real(ndtdrag)   ,dtdragfacmax )
 if (ndtdragd   > 0)  call summary_variable('dt',iosumdtdd ,ndtdragd  ,dtdragdfacmean/real(ndtdragd)  ,dtdragdfacmax)
 if (ndtvisc    > 0)  call summary_variable('dt',iosumdtv  ,ndtvisc   ,dtviscfacmean /real(ndtvisc)   ,dtviscfacmax)
 if (ndtdust    > 0)  call summary_variable('dt',iosumdte  ,ndtdust   ,dtdustfacmean /real(ndtdust)   ,dtdustfacmax)
 if (mhd_nonideal) then
    if (ndtohm  > 0)  call summary_variable('dt',iosumdto  ,ndtohm    ,dtohmfacmean  /real(ndtohm)    ,dtohmfacmax  )
    if (ndthall > 0)  call summary_variable('dt',iosumdth  ,ndthall   ,dthallfacmean /real(ndthall)   ,dthallfacmax )
    if (ndtambi > 0)  call summary_variable('dt',iosumdta  ,ndtambi   ,dtambifacmean /real(ndtambi)   ,dtambifacmax )
 endif
 if (ndtrad     > 0)  call summary_variable('dt',iosumdtr  ,ndtrad    ,dtradfacmean  /real(ndtrad)    ,dtradfacmax  )
 if (ndtclean   > 0)  call summary_variable('dt',iosumdtB  ,ndtclean  ,dtcleanfacmean/real(ndtclean)  ,dtcleanfacmax)
#else

 dtcourant = reduceall_mpi('min',dtcourant)
 dtforce   = reduceall_mpi('min',dtforce)
 dtvisc    = reduceall_mpi('min',dtvisc)
 dtmini    = reduce_mpi('min',dtmini)
 dtmaxi    = reduce_mpi('max',dtmaxi)

 if (iverbose >= 2 .and. id==master) write(iprint,*) 'dtmin = ',C_Cour*dtmini, ' dtmax = ',C_cour*dtmaxi, &
    ' dtmax/dtmin = ',dtmaxi/(dtmini + epsilon(0.)),'dtcour/dtf = ',(C_cour*dtcourant)/(C_force*dtforce + epsilon(0.))
 if ( dtforce < dtcourant ) call summary_variable('dt',iosumdtf,0,0.0)
 if ( dtvisc  < dtcourant ) call summary_variable('dt',iosumdtv,0,0.0)
 if ( mhd_nonideal ) then
    ! Note: We are not distinguishing between use_STS and .not.use_STS here since if
    !       use_STS==.true., then dtohm=dtambi=bignumber.
    dtohm    = reduceall_mpi('min',dtohm )
    dthall   = reduceall_mpi('min',dthall)
    dtambi   = reduceall_mpi('min',dtambi)
    if ( dthall < dtcourant ) call summary_variable('dt',iosumdth,0,0.0)
    if ( dtohm  < dtcourant ) call summary_variable('dt',iosumdto,0,0.0)
    if ( dtambi < dtcourant ) call summary_variable('dt',iosumdta,0,0.0)

    minglobdt = min(dtvisc,dtohm,dthall,dtambi)

    if (minglobdt < dtcourant) then
       dtcourant = minglobdt
       if      (abs(dtcourant-dtvisc) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','viscosity constraining timestep')
          call summary_variable('dt',iosumdtv,0,0.0,0.0, .true. )
       elseif (abs(dtcourant-dthall) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','Hall Effect constraining timestep')
          call summary_variable('dt',iosumdth,0,0.0,0.0, .true. )
       elseif (abs(dtcourant-dtohm ) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','ohmic resistivity constraining timestep')
          call summary_variable('dt',iosumdto,0,0.0,0.0, .true. )
       elseif (abs(dtcourant-dtambi) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','ambipolar diffusion constraining timestep')
          call summary_variable('dt',iosumdta,0,0.0,0.0, .true. )
       endif
    endif
 else
    if (dtvisc < dtcourant) then
       dtcourant = dtvisc
       if (iverbose >= 1 .and. id==master) call warning('force','viscosity constraining timestep')
       call summary_variable('dt',iosumdtv,0,0.0,0.0, .true. )
    endif
 endif

 if (do_radiation) then
    dtrad = reduceall_mpi('min',dtrad)
    if (dtrad < dtcourant) then
       call summary_variable('dt',iosumdtr,0,0.0)
       if (iverbose >= 1 .and. id==master) &
          call warning('force','radiation is constraining timestep')
       call summary_variable('dt',iosumdtr,0,0.0,0.0,.true.)
    endif
 endif

 if ( dtforce < dtcourant ) call summary_variable('dt',iosumdtf,0,0.0,0.0, .true. )
#endif

end subroutine force

!----------------------------------------------------------------
!+
!  Internal subroutine that computes the force summations
!
!  MAKE SURE THIS ROUTINE IS INLINED BY THE COMPILER
!+
!----------------------------------------------------------------
subroutine compute_forces(i,iamgasi,iamdusti,xpartveci,hi,hi1,hi21,hi41,gradhi,gradsofti, &
                          beta, &
                          pmassi,listneigh,nneigh,xyzcache,fsum,vsigmax, &
                          ifilledcellcache,realviscosity,useresistiveheat, &
                          xyzh,vxyzu,Bevol,iphasei,iphase,massoftype, &
                          divcurlB,eta_nimhd, eos_vars, &
                          dustfrac,dustprop,fxyz_drag,gradh,divcurlv,alphaind, &
                          alphau,alphaB,bulkvisc,stressmax,&
                          ndrag,nstokes,nsuper,ts_min,ibinnow_m1,ibin_wake,ibin_neighi,&
                          ignoreself,rad,radprop,dens,metrics,apr_level,dt)
#ifdef FINVSQRT
 use fastmath,    only:finvsqrt
#endif
 use kernel,      only:grkern,cnormk,radkern2
 use part,        only:igas,idust,iohm,ihall,iambi,maxphase,iactive,&
                       iamtype,iamdust,get_partinfo,mhd,maxvxyzu,maxdvdx,igasP,ics,iradP,itemp
 use dim,         only:maxalpha,maxp,mhd_nonideal,gravity,gr,use_apr
 use part,        only:rhoh,dvdx,aprmassoftype
 use nicil,       only:nimhd_get_jcbcb,nimhd_get_dBdt
 use eos,         only:ieos,eos_is_non_ideal
 use eos_stamatellos, only:gradP_cool,getopac_opdep
#ifdef GRAVITY
 use kernel,      only:kernel_softening
 use ptmass,      only:ptmass_not_obscured
#endif
#ifdef PERIODIC
 use boundary,    only:dxbound,dybound,dzbound
#endif
 use dim,         only:use_dust,use_dustgrowth,ind_timesteps
 use dust,        only:get_ts,idrag,icut_backreaction,ilimitdustflux,irecon,drag_implicit
 use kernel,      only:wkern_drag,cnormk_drag,wkern,cnormk
 use part,        only:ndustsmall,grainsize,graindens,ndustsmall,grainsize,graindens,filfac
 use options,     only:use_porosity,icooling
 use growth,      only:get_size
 use kernel,      only:wkern,cnormk
#ifdef IND_TIMESTEPS
 use part,        only:ibin_old,iamboundary
 use timestep_ind,only:get_dt
#endif
 use timestep,    only:bignumber
 use options,     only:overcleanfac,use_dustfrac,ireconav,limit_radiation_flux
 use units,       only:get_c_code
 use metric_tools,only:imet_minkowski,imetric
 use utils_gr,    only:get_bigv
 use radiation_utils, only:get_rad_R
 use io,          only:fatal
 integer,         intent(in)    :: i
 logical,         intent(in)    :: iamgasi,iamdusti
 real,            intent(in)    :: xpartveci(:)
 real(kind=8),    intent(in)    :: hi1,hi21,hi41,gradhi,gradsofti
 real,            intent(in)    :: hi,beta
 real,            intent(in)    :: pmassi
 integer,         intent(in)    :: listneigh(:)
 integer,         intent(in)    :: nneigh
 real,            intent(in)    :: xyzcache(:,:)
 real,            intent(out)   :: fsum(maxfsum)
 real,            intent(out)   :: vsigmax
 logical,         intent(in)    :: ifilledcellcache
 logical,         intent(in)    :: realviscosity,useresistiveheat
 real,            intent(in)    :: xyzh(:,:)
 real,            intent(inout) :: vxyzu(:,:)
 real,            intent(in)    :: Bevol(:,:)
 real(kind=4),    intent(in)    :: divcurlB(:,:)
 real,            intent(in)    :: dustfrac(:,:)
 real,            intent(in)    :: dustprop(:,:)
 real,            intent(in)    :: fxyz_drag(:,:)
 integer(kind=1), intent(in)    :: iphasei
 integer(kind=1), intent(in)    :: iphase(:)
 real,            intent(in)    :: massoftype(:)
 real,            intent(in)    :: eta_nimhd(:,:)
 real,            intent(in)    :: eos_vars(:,:)
 real(kind=4),    intent(in)    :: alphaind(:,:)
 real(kind=4),    intent(in)    :: gradh(:,:),divcurlv(:,:)
 real,            intent(in)    :: alphau,alphaB,bulkvisc,stressmax
 integer,         intent(inout) :: ndrag,nstokes,nsuper
 real,            intent(out)   :: ts_min
 integer(kind=1), intent(out)   :: ibin_wake(:),ibin_neighi
 integer(kind=1), intent(in)    :: ibinnow_m1
 logical,         intent(in)    :: ignoreself
 real,            intent(in)    :: rad(:,:),dens(:),metrics(:,:,:,:)
 real,            intent(inout) :: radprop(:,:)
 integer(kind=1), intent(in)    :: apr_level(:)
 real,            intent(in)    :: dt
 integer :: j,n,iamtypej
 logical :: iactivej,iamgasj,iamdustj
 real    :: rij2,q2i,qi,xj,yj,zj,dx,dy,dz,runix,runiy,runiz,rij1,hfacgrkern
 real    :: grkerni,grgrkerni,dvx,dvy,dvz,projv,denij,vsigi,vsigu,dudtdissi
 real    :: projBi,projBj,dBx,dBy,dBz,dB2,projdB
 real    :: dendissterm,dBdissterm,dudtresist,dpsiterm,pmassonrhoi
 real    :: gradpi,projsxi,projsyi,projszi
 real    :: gradp,projsx,projsy,projsz,Bxj,Byj,Bzj,Bj,Bj1,psij
 real    :: grkernj,grgrkernj,autermj,avBtermj,vsigj,spsoundj,tempj
 real    :: gradpj,pro2j,projsxj,projsyj,projszj,sxxj,sxyj,sxzj,syyj,syzj,szzj,dBrhoterm
 real    :: visctermisoj,visctermanisoj,enj,hj,mrhoj5,alphaj,pmassj,rho1j
 real    :: rhoj,prj,rhoav1
 real    :: hj1,hj21,q2j,qj,vwavej,divvj
 real    :: dvdxi(9),dvdxj(9)
#ifdef GRAVITY
 real    :: fmi,fmj,dsofti,dsoftj
 logical :: add_contribution
#else
 logical, parameter :: add_contribution = .true.
#endif
 real    :: phi,phii,phij,fgrav,fgravi,fgravj,termi
#ifdef KROME
 real    :: gammaj
#endif
 integer :: iregime,idusttype,l
 real    :: dragterm,dragheating,wdrag,dv2,tsijtmp
 real    :: grkernav,tsj(maxdusttypes),dustfracterms(maxdusttypes),term
 real    :: projvstar,projf_drag,epstsj,sdrag1,sdrag2!,rhogas1i
 real    :: winter
 real    :: dBevolx,dBevoly,dBevolz,divBsymmterm,divBdiffterm
 real    :: rho21i,rho21j,Bxi,Byi,Bzi,psii,pmjrho21grkerni,pmjrho21grkernj
 real    :: auterm,avBterm,mrhoi5,vsigB
 real    :: jcbcbj(3),jcbj(3),dBnonideal(3),dBnonidealj(3),divBi,curlBi(3),curlBj(3)
 real    :: vsigavi,vsigavj
 real    :: dustfraci(maxdusttypes),dustfracj(maxdusttypes),tsi(maxdusttypes)
 real    :: sqrtrhodustfraci(maxdusttypes),sqrtrhodustfracj(maxdusttypes)
 real    :: dustfracisum,dustfracjsum,rhogasj,epstsi!,rhogas1j
 real    :: vwavei,rhoi,rho1i,spsoundi,tempi
 real    :: sxxi,sxyi,sxzi,syyi,syzi,szzi
 real    :: visctermiso,visctermaniso
 real    :: pri,pro2i
 real    :: etaohmi,etahalli,etaambii
 real    :: jcbcbi(3),jcbi(3)
 real    :: alphai,grainmassi,graindensi,filfaci
 logical :: usej
 integer :: iamtypei
 real    :: radFi(3),radFj(3),radRj,radDFWi,radDFWj,c_code,radkappai,radkappaj,&
            radDi,radDj,radeni,radenj,radlambdai,radlambdaj
 real    :: xi,yi,zi,densi,eni,metrici(0:3,0:3,2)
 real    :: vxi,vyi,vzi,vxj,vyj,vzj,projvi,projvj
 real    :: fxi_drag,fyi_drag,fzi_drag,dti
 real    :: qrho2i,qrho2j
 integer :: ii,ia,ib,ic
 real    :: densj
 real    :: bigvi(1:3),bigv2i,alphagri,lorentzi
 real    :: veli(3),vij,vijstar
 real    :: bigv2j,alphagrj,enthi,enthj
 real    :: dlorentzv,lorentzj,lorentzi_star,lorentzj_star,projbigvi,projbigvj
 real    :: bigvj(1:3),velj(3),metricj(0:3,0:3,2),projbigvstari,projbigvstarj
 real    :: radPj,fgravxi,fgravyi,fgravzi
 real    :: gradpx,gradpy,gradpz,gradP_cooli,gradP_coolj

 ! unpack
 xi            = xpartveci(ixi)
 yi            = xpartveci(iyi)
 zi            = xpartveci(izi)
 vxi           = xpartveci(ivxi)
 vyi           = xpartveci(ivyi)
 vzi           = xpartveci(ivzi)
 eni           = xpartveci(ieni)
 vwavei        = xpartveci(ivwavei)
 rhoi          = xpartveci(irhoi)
 rho1i         = 1./rhoi
 spsoundi      = xpartveci(ispsoundi)
 tempi         = xpartveci(itempi)
 sxxi          = xpartveci(isxxi)
 sxyi          = xpartveci(isxyi)
 sxzi          = xpartveci(isxzi)
 syyi          = xpartveci(isyyi)
 syzi          = xpartveci(isyzi)
 szzi          = xpartveci(iszzi)
 visctermiso   = xpartveci(ivisctermisoi)
 visctermaniso = xpartveci(ivisctermanisoi)
 pri           = xpartveci(ipri)
 pro2i         = xpartveci(ipro2i)
 etaohmi       = xpartveci(ietaohmi)
 etahalli      = xpartveci(ietahalli)
 etaambii      = xpartveci(ietaambii)
 jcbcbi(1)     = xpartveci(ijcbcbxi)
 jcbcbi(2)     = xpartveci(ijcbcbyi)
 jcbcbi(3)     = xpartveci(ijcbcbzi)
 jcbi(1)       = xpartveci(ijcbxi)
 jcbi(2)       = xpartveci(ijcbyi)
 jcbi(3)       = xpartveci(ijcbzi)
 alphai        = xpartveci(ialphai)
 divBi         = xpartveci(idivBi)
 curlBi(1)     = xpartveci(icurlBxi)
 curlBi(2)     = xpartveci(icurlByi)
 curlBi(3)     = xpartveci(icurlBzi)

 if (gr) then
    densi      = xpartveci(idensGRi)
    ii = imetricstart
    do ic = 1,2
       do ib = 0,3
          do ia = 0,3
             metrici(ia,ib,ic) = xpartveci(ii)
             ii = ii + 1
          enddo
       enddo
    enddo
 endif

 if (use_dustgrowth) then
    grainmassi = xpartveci(igrainmassi)
    graindensi = xpartveci(igraindensi)
    filfaci    = xpartveci(ifilfaci)
 endif
 fxi_drag = xpartveci(ifxi_drag)
 fyi_drag = xpartveci(ifyi_drag)
 fzi_drag = xpartveci(ifzi_drag)
 if (ind_timesteps) then
    dti      = xpartveci(idti)
 endif
 dvdxi(1:9)    = xpartveci(idvxdxi:idvzdzi)

 if (gr) then
    veli = [vxi,vyi,vzi]
    call get_bigv(metrici,veli,bigvi,bigv2i,alphagri,lorentzi)
 endif

 fsum(:) = 0.
 vsigmax = 0.
 pmassonrhoi = pmassi*rho1i
 hfacgrkern  = hi41*cnormk*gradhi

 iamtypei = iamtype(iphasei)

 ! default settings for active/phase if iphase not used
 iactivej = .true.
 iamtypej = igas
 iamgasj  = .true.
 iamdustj = .false.

 ! to find max ibin of all of i's neighbours
 ibin_neighi = 0_1

 ! dust
 ts_min = bignumber

 ! radiation
 c_code = get_c_code()

 ! various pre-calculated quantities
 if (mhd) then
    Bxi  = xpartveci(iBevolxi)*rhoi
    Byi  = xpartveci(iBevolyi)*rhoi
    Bzi  = xpartveci(iBevolzi)*rhoi
    psii = xpartveci(ipsi)
 else
    Bxi  = 0.0
    Byi  = 0.0
    Bzi  = 0.0
    psii = 0.0
 endif
 if (use_dust .and. use_dustfrac) then
    dustfraci(:) = xpartveci(idustfraci:idustfraciend)
    dustfracisum = sum(dustfraci(:))
    if (ilimitdustflux) then
       tsi(:) = min(xpartveci(itstop:itstopend),hi/spsoundi) ! flux limiter from Ballabio et al. (2018)
    else
       tsi(:) = xpartveci(itstop:itstopend)
    endif
    epstsi = sum(dustfraci(:)*tsi(:))
!--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
    sqrtrhodustfraci(:) = sqrt(dustfraci(:)/(1.-dustfraci(:)))
 else
    dustfraci(:) = 0.
    dustfracisum = 0.
    tsi(:)       = 0.
    epstsi       = 0.
    sqrtrhodustfraci(:) = 0.
 endif
 rho21i = rho1i*rho1i
 mrhoi5  = 0.5*pmassi*rho1i
 !avterm  = mrhoi5*alphai       !  artificial viscosity parameter
 auterm  = mrhoi5*alphau       !  artificial thermal conductivity parameter
 avBterm = mrhoi5*alphaB*rho1i
!
!--initialise the following to zero for the case
!
 usej      = .false.
 grkernj   = 0.
 alphaj    = alphai
 divvj     = 0.
 dvdxj(:)  = 0.
 rhoj      = 0.
 rho1j     = 0.
 mrhoj5    = 0.
 gradpj    = 0.
 projsxj   = 0.
 projsyj   = 0.
 projszj   = 0.
 dudtresist = 0.
 dpsiterm   = 0.
 fgravi = 0.
 fgravj = 0.
 phii   = 0.
 phij   = 0.
 phi    = 0.
 dBnonideal(:) = 0.0
 Bxj = 0.
 Byj = 0.
 Bzj = 0.
 psij = 0.
 visctermisoj = 0.
 visctermanisoj = 0.
 fgravxi = 0.
 fgravyi = 0.
 fgravzi = 0.
 if (icooling == 9) then
    gradP_cool(i) = 0.
    gradpx = 0.
    gradpy = 0.
    gradpz = 0.
    gradP_cooli=0.
    gradP_coolj=0.
 endif

 loop_over_neighbours2: do n = 1,nneigh

    j = abs(listneigh(n))
    if ((ignoreself) .and. (i==j)) cycle loop_over_neighbours2

    if (ifilledcellcache .and. n <= maxcellcache) then
       ! positions from cache are already mod boundary
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
    q2i = rij2*hi21
    !--hj is in the cell cache but not in the neighbour cache
    !  as not accessed during the density summation
    if (ifilledcellcache .and. n <= maxcellcache) then
       hj1 = xyzcache(n,4)
    else
       hj1 = 1./xyzh(4,j)
    endif
    hj21 = hj1*hj1
    q2j  = rij2*hj21
    is_sph_neighbour: if (q2i < radkern2 .or. q2j < radkern2) then
#ifdef GRAVITY
       !  Determine if neighbouring particle is hidden by a sink particle;
       !  if so, do not add contribution.
       add_contribution = .true.
       !k = 1
       !do while (k <= nptmass .and. add_contribution)
       !   xkpt = xyzmh_ptmass(1,k)
       !   ykpt = xyzmh_ptmass(2,k)
       !   zkpt = xyzmh_ptmass(3,k)
       !   vpos = (xkpt-xpartveci(ixi))*(xkpt-xj) &
       !        + (ykpt-xpartveci(iyi))*(ykpt-yj) &
       !        + (zkpt-xpartveci(izi))*(zkpt-zj)
       !   if (vpos < 0.0) then
       !      add_contribution = ptmass_not_obscured(-dx,-dy,-dz,  &
       !                          xkpt-xpartveci(ixi),ykpt-xpartveci(iyi),zkpt-xpartveci(izi), &
       !                          xyzmh_ptmass(ihacc,k))
       !   endif
       !   k = k + 1
       !enddo
#endif

       if (rij2 > epsilon(rij2)) then
#ifdef FINVSQRT
          rij1 = finvsqrt(rij2)
#else
          rij1 = 1./sqrt(rij2)
#endif
          qi   = (rij2*rij1)*hi1  ! this is qi = rij*hi1
       else
          rij1 = 0.
          qi   = 0.
       endif

       if (q2i < radkern2) then
          grkerni = grkern(q2i,qi)*hfacgrkern
#ifdef GRAVITY
          call kernel_softening(q2i,qi,phii,fmi)
          phii   = phii*hi1
          fmi    = fmi*hi21
          dsofti = gradsofti*grkerni
          fgravi = fmi + dsofti
#endif
       else
          grkerni = 0.
#ifdef GRAVITY
          phii   = -rij1
          fmi    = rij1*rij1
          fgravi = fmi
#endif
       endif

       runix = dx*rij1
       runiy = dy*rij1
       runiz = dz*rij1
       !
       !--compute the contribution neighbours with h_j and grad W (h_j)
       !
       if (q2j < radkern2) then
          qj = (rij2*rij1)*hj1
          grkernj = grkern(q2j,qj)*hj21*hj21*cnormk*gradh(1,j) ! ndim + 1
#ifdef GRAVITY
          call kernel_softening(q2j,qj,phij,fmj)
          fmj    = fmj*hj21
          dsoftj = gradh(2,j)*grkernj
          fgravj = fmj + dsoftj
#endif
          usej = .true.
       else
          grkernj = 0.
#ifdef GRAVITY
          fmj    = rij1*rij1
          fgravj = fmj
#endif
          usej = .false.
       endif
       if (mhd) usej = .true.
       if (use_dust) usej = .true.
       if (maxvxyzu >= 4 .and. .not.gravity) usej = .true.

       !--get individual timestep/ multiphase information (querying iphase)
       if (maxphase==maxp) then
          call get_partinfo(iphase(j),iactivej,iamgasj,iamdustj,iamtypej)
#ifdef IND_TIMESTEPS
          ! Particle j is a neighbour of an active particle;
          ! flag it to see if it needs to be woken up next step.
          if (.not.iamboundary(iamtypej)) then
             ibin_wake(j)  = max(ibinnow_m1,ibin_wake(j))
             ibin_neighi = max(ibin_neighi,ibin_old(j))
          endif
          if (use_dust .and. drag_implicit) then
             dti = min(dti,get_dt(dt,ibin_old(j)))
          endif
#endif
       endif
       if (use_apr) then
          pmassj = aprmassoftype(iamtypej,apr_level(j))
       else
          pmassj = massoftype(iamtypej)
       endif

       fgrav = 0.5*(pmassj*fgravi + pmassi*fgravj)

       !  If particle is hidden by the sink, treat the neighbour as
       !  not gas; gravitational contribution will be added after the
       !  isgas if-statement
       if (.not. add_contribution) then
          iamgasj = .false.
          usej    = .false.
       endif

       !--get dv : needed for timestep and av term
       vxj = vxyzu(1,j)
       vyj = vxyzu(2,j)
       vzj = vxyzu(3,j)
       projvi = vxi*runix + vyi*runiy + vzi*runiz
       projvj = vxj*runix + vyj*runiy + vzj*runiz
       dvx = vxi - vxj
       dvy = vyi - vyj
       dvz = vzi - vzj

       ! do the projection here
       projv = dvx*runix + dvy*runiy + dvz*runiz

       if (iamgasj .and. maxvxyzu >= 4) then
          enj = vxyzu(4,j)
          if (eos_is_non_ideal(ieos)) then  ! only do this if eos requires temperature in physical units
             tempj = eos_vars(itemp,j)
             denij = 0.5*(eni/tempi + enj/tempj)*(tempi - tempj)  ! dU = c_V * dT
          else
             denij = eni - enj
          endif
       else
          denij = 0.
          enj = 0.
       endif

       if (gr) then
          !-- Get velocites required in GR shock capturing
          velj  = [vxj,vyj,vzj]
          metricj(0:3,0:3,1:2) = metrics(:,:,:,j)
          call get_bigv(metricj,velj,bigvj,bigv2j,alphagrj,lorentzj)

          ! Reduce problem to 1D along the line of sight
          projbigvi = bigvi(1)*runix + bigvi(2)*runiy + bigvi(3)*runiz !dot_product(bigvi,rij)
          projbigvj = bigvj(1)*runix + bigvj(2)*runiy + bigvj(3)*runiz

          ! Reconstruction of velocity
          if (maxdvdx==maxp) dvdxj(:) = dvdx(:,j)
          projbigvstari = projbigvi
          projbigvstarj = projbigvj
          if (ireconav >= 0) call reconstruct_dv_gr(projbigvi,projbigvj,runix,runiy,runiz,rij2*rij1,&
                                 dvdxi,dvdxj,projbigvstari,projbigvstarj,ireconav)

          ! Relativistic version of vi-vj
          vij = abs((projbigvi-projbigvj)/(1.-projbigvi*projbigvj))
          vijstar = abs((projbigvstari-projbigvstarj)/(1.-projbigvstari*projbigvstarj))
       endif

       if (iamgasi .and. iamgasj) then
          !--work out vsig for timestepping and av
          if (gr) then
             ! Relativistic version vij + csi    (could put a beta here somewhere?)
             vsigi     = (beta*vij+spsoundi)/(1.+beta*vij*spsoundi)
             vsigavi   = alphai*vsigi
          else
             vsigi     = max(vwavei - beta*projv,0.)
             vsigavi   = max(alphai*vwavei - beta*projv,0.)
          endif
          if (vsigi > vsigmax) vsigmax = vsigi

          if (mhd) then
             hj   = xyzh(4,j)
             rhoj = rhoh(hj,pmassj)
             Bxj  = Bevol(1,j)*rhoj
             Byj  = Bevol(2,j)*rhoj
             Bzj  = Bevol(3,j)*rhoj
             psij = Bevol(4,j)

             dBx = Bxi - Bxj
             dBy = Byi - Byj
             dBz = Bzi - Bzj
             projBi = Bxi*runix + Byi*runiy + Bzi*runiz
             if (usej) projBj = Bxj*runix + Byj*runiy + Bzj*runiz
             projdB = dBx*runix + dBy*runiy + dBz*runiz
             dB2 = dBx*dBx + dBy*dBy + dBz*dBz
             divBdiffterm = -pmassj*projdB*grkerni
          endif
       else
          !-- v_sig for pairs of particles that are not gas-gas
          vsigi = max(-projv,0.0)
          if (vsigi > vsigmax) vsigmax = vsigi
          vsigavi = 0.
          dBx = 0.; dBy = 0.; dBz = 0.; dB2 = 0.
          projBi = 0.; projBj = 0.; divBdiffterm = 0.
       endif

       !--get terms required for particle j
       if (usej) then
          hj       = 1./hj1
          rhoj     = rhoh(hj,pmassj)
          rho1j    = 1./rhoj
          rho21j   = rho1j*rho1j

          if (maxdvdx==maxp) dvdxj(:) = dvdx(:,j)

          if (iamgasj) then
             if (ndivcurlv >= 1) divvj = divcurlv(1,j)
             if (use_dustfrac) then
                dustfracj(:) = dustfrac(:,j)
                dustfracjsum = sum(dustfracj(:))
                rhogasj      = rhoj*(1. - dustfracjsum)
                !--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
                sqrtrhodustfracj(:) = sqrt(dustfracj(:)/(1.-dustfracj(:)))
             else
                dustfracj(:) = 0.
                dustfracjsum = 0.
                rhogasj      = rhoj
                sqrtrhodustfracj(:) = 0.
             endif

             if (maxalpha==maxp) alphaj  = alphaind(1,j)

             prj = eos_vars(igasP,j)
             spsoundj = eos_vars(ics,j)
             radPj = 0.
             if (do_radiation) radPj = radprop(iradP,j)
             !
             !--calculate j terms (which were precalculated outside loop for i)
             !
             call get_stress(prj,spsoundj,rhoj,rho1j,xj,yj,zj,pmassj,Bxj,Byj,Bzj, &
                        pro2j,vwavej, &
                        sxxj,sxyj,sxzj,syyj,syzj,szzj,visctermisoj,visctermanisoj, &
                        realviscosity,divvj,bulkvisc,dvdxj,stressmax,radPj)

             mrhoj5   = 0.5*pmassj*rho1j
             autermj  = mrhoj5*alphau
             avBtermj = mrhoj5*alphaB*rho1j

             if (gr) then
                ! Relativistic version vij + csi
                vsigj   = (beta*vij+spsoundj)/(1.+beta*vij*spsoundj)
                vsigavj = alphaj*vsigj
             else
                vsigj = max(vwavej - beta*projv,0.)
                vsigavj = max(alphaj*vwavej - beta*projv,0.)
             endif
             if (vsigj > vsigmax) vsigmax = vsigj
          else
             vsigj = max(-projv,0.)
             if (vsigj > vsigmax) vsigmax = vsigj
             vsigavj = 0.; vwavej = 0.; avBtermj = 0.; autermj = 0. ! avoid compiler warnings
             sxxj = 0.; sxyj = 0.; sxzj = 0.; syyj = 0.; syzj = 0.; szzj = 0.; pro2j = 0.; prj = 0.
             dustfracj = 0.; dustfracjsum = 0.; sqrtrhodustfracj = 0.
          endif
       else ! set to zero terms which are used below without an if (usej)
          !rhoj      = 0.
          rho1j     = 0.
          rho21j    = 0.

          mrhoj5    = 0.
          autermj   = 0.
          avBtermj  = 0.
          psij = 0.

          gradpj    = 0.
          projsxj   = 0.
          projsyj   = 0.
          projszj   = 0.
          projBj = 0.
          prj   = 0.
          pro2j = 0.
          vwavej = 0.
          vsigavj = 0.
          spsoundj = 0.
          dustfracj = 0.
          dustfracjsum = 0.
          sqrtrhodustfracj = 0.
          dvdxj(:) = 0.
          sxxj = 0.; sxyj = 0.; sxzj = 0.; syyj = 0.; syzj = 0.; szzj = 0.
       endif

       ifgas: if (iamgasi .and. iamgasj) then

          !
          !--artificial viscosity term
          !
          qrho2i = 0.
          qrho2j = 0.
          if (gr) then
             densj = dens(j)
             enthi  = 1.+eni+pri/densi
             enthj  = 1.+enj+prj/densj
          endif

!------------------
#ifdef DISC_VISCOSITY
!------------------
          !
          !--This is for "physical" disc viscosity
          !  (We multiply by h/rij, use cs for the signal speed, apply to both approaching/receding,
          !   with beta viscosity only applied to approaching pairs)
          !
          if (gr) then
             lorentzi_star = 1./sqrt(1.-projbigvstari**2)
             lorentzj_star = 1./sqrt(1.-projbigvstarj**2)
             dlorentzv = lorentzi_star*projbigvstari - lorentzj_star*projbigvstarj
             if (projv < 0.) then
                qrho2i = -0.5*rho1i*vsigavi*enthi*dlorentzv*hi*rij1
                if (usej) qrho2j = -0.5*rho1j*vsigavj*enthj*dlorentzv*hj*rij1
             else
                qrho2i = -0.5*rho1i*alphai*spsoundi*enthi*dlorentzv*hi*rij1
                if (usej) qrho2j = -0.5*rho1j*alphaj*spsoundj*enthj*dlorentzv*hj*rij1
             endif
             if (ien_type == ien_etotal) then ! total energy
                dudtdissi = - pmassj*((pro2i + qrho2i)*projvj*grkerni + &
                                      (pro2j + qrho2j)*projvi*grkernj)
             else
                dudtdissi = -0.5*pmassj*rho1i*alphai*spsoundi*enthi*dlorentzv*hi*rij1*projv*grkerni
             endif
          else
             if (projv < 0.) then
                qrho2i = - 0.5*rho1i*(alphai*spsoundi - beta*projv)*hi*rij1*projv
                if (usej) qrho2j = - 0.5*rho1j*(alphaj*spsoundj - beta*projv)*hj*rij1*projv
             else
                qrho2i = - 0.5*rho1i*alphai*spsoundi*hi*rij1*projv
                if (usej) qrho2j = - 0.5*rho1j*alphaj*spsoundj*hj*rij1*projv
             endif
             dudtdissi = -0.5*pmassj*rho1i*alphai*spsoundi*hi*rij1*projv**2*grkerni
          endif

!--DISC_VISCOSITY--
#else
!------------------
          if (projv < 0.) then
             if (gr) then
                lorentzi_star = 1./sqrt(1.-projbigvstari**2)
                lorentzj_star = 1./sqrt(1.-projbigvstarj**2)
                dlorentzv = lorentzi_star*projbigvstari - lorentzj_star*projbigvstarj
                qrho2i = -0.5*rho1i*vsigavi*enthi*dlorentzv
                if (usej) qrho2j = -0.5*rho1j*vsigavj*enthj*dlorentzv
             else
                qrho2i = - 0.5*rho1i*vsigavi*projv
                if (usej) qrho2j = - 0.5*rho1j*vsigavj*projv
             endif
          endif
          !--energy conservation from artificial viscosity
          if (ien_type == ien_etotal) then ! total energy
             dudtdissi = - pmassj*((pro2i + qrho2i)*projvj*grkerni + &
                                   (pro2j + qrho2j)*projvi*grkernj)
          else
             dudtdissi = pmassj*qrho2i*projv*grkerni
          endif
!--DISC_VISCOSITY--
#endif
!------------------

          !--add av term to pressure
          gradpi = pmassj*(pro2i + qrho2i)*grkerni
          if (usej) gradpj = pmassj*(pro2j + qrho2j)*grkernj
          !-- calculate grad P from gas pressure alone for cooling
          ! .not. mhd required to past jetnimhd test
          if (icooling == 9 .and. .not. mhd) then
             gradP_cooli =  pmassj*pri*rho1i*rho1i*grkerni
             gradP_coolj = 0.
             if (usej) then
                gradp_coolj =  pmassj*prj*rho1j*rho1j*grkernj
             endif
          endif

          !--artificial thermal conductivity (need j term)
          if (maxvxyzu >= 4) then
             if (gr) then
                denij = alphagri*eni/lorentzi - alphagrj*enj/lorentzj
                if (imetric==imet_minkowski) then  ! Eq 60 in LP19
                   rhoav1 = 2./(enthi*densi + enthj*densj)
                   vsigu = min(1.,sqrt(abs(pri-prj)*rhoav1))
                else
                   vsigu = abs(vij)  ! Eq 61 in LP19
                endif
             elseif (gravity) then
                vsigu = abs(projv)
             else
                rhoav1 = 2./(rhoi + rhoj)
                vsigu = sqrt(abs(pri - prj)*rhoav1)
             endif
             dendissterm = vsigu*denij*(auterm*grkerni + autermj*grkernj)
          endif

          if (mhd) then
             !
             ! artificial resistivity
             !
             vsigB = sqrt((dvx - projv*runix)**2 + (dvy - projv*runiy)**2 + (dvz - projv*runiz)**2)
             dBdissterm = (avBterm*grkerni + avBtermj*grkernj)*vsigB

             !--energy dissipation due to artificial resistivity
             if (useresistiveheat) dudtresist = -0.5*dB2*dBdissterm

             pmjrho21grkerni = pmassj*rho21i*grkerni
             pmjrho21grkernj = pmassj*rho21j*grkernj

             termi        = pmjrho21grkerni*projBi
             divBsymmterm =  termi + pmjrho21grkernj*projBj
             dBrhoterm    = -termi
             !
             ! grad psi term for divergence cleaning
             ! cleaning evolving d/dt (psi/c_h) as in Tricco, Price & Bate (2016)
             !
             dpsiterm = overcleanfac*(pmjrho21grkerni*psii*vwavei + pmjrho21grkernj*psij*vwavej)
             ! coefficient for associated cleaning timestep
             Bj = sqrt(Bxj**2 + Byj**2 + Bzj**2)
             if (Bj > 0.0) then
                Bj1 = 1.0/Bj
             else
                Bj1 = 0.0
             endif
             !
             ! non-ideal MHD terms
             !
             if (mhd_nonideal) then
                call nimhd_get_dBdt(dBnonideal,etaohmi,etahalli,etaambii,curlBi,jcbi,jcbcbi,runix,runiy,runiz)
                dBnonideal = dBnonideal*pmjrho21grkerni
                curlBj = divcurlB(2:4,j)
                call nimhd_get_jcbcb(jcbcbj,jcbj,curlBj,Bxj,Byj,Bzj,Bj1)
                call nimhd_get_dBdt(dBnonidealj,eta_nimhd(iohm,j),eta_nimhd(ihall,j),eta_nimhd(iambi,j) &
                                   ,curlBj,jcbj,jcbcbj,runix,runiy,runiz)
                dBnonideal = dBnonideal + dBnonidealj*pmjrho21grkernj
             endif
             !
             ! dB/dt evolution equation
             !
             dBevolx = dBrhoterm*dvx + dBdissterm*dBx - dpsiterm*runix - dBnonideal(1)
             dBevoly = dBrhoterm*dvy + dBdissterm*dBy - dpsiterm*runiy - dBnonideal(2)
             dBevolz = dBrhoterm*dvz + dBdissterm*dBz - dpsiterm*runiz - dBnonideal(3)
          endif
          !
          !--get projection of anisotropic part of stress tensor
          !  in direction of particle pair
          !
          projsxi = (sxxi*runix + sxyi*runiy + sxzi*runiz)*grkerni
          projsyi = (sxyi*runix + syyi*runiy + syzi*runiz)*grkerni
          projszi = (sxzi*runix + syzi*runiy + szzi*runiz)*grkerni
          if (usej) then
             projsxj = (sxxj*runix + sxyj*runiy + sxzj*runiz)*grkernj
             projsyj = (sxyj*runix + syyj*runiy + syzj*runiz)*grkernj
             projszj = (sxzj*runix + syzj*runiy + szzj*runiz)*grkernj
          endif
          !
          !--physical viscosity term (direct second derivatives)
          !
          if (realviscosity .and. maxdvdx /= maxp) then
             grgrkerni = -2.*grkerni*rij1
             gradpi = gradpi + visctermiso*projv*grgrkerni
             projsxi = projsxi + visctermaniso*dvx*grgrkerni
             projsyi = projsyi + visctermaniso*dvy*grgrkerni
             projszi = projszi + visctermaniso*dvz*grgrkerni
             dudtdissi = dudtdissi + grgrkerni*(visctermiso*projv**2 &
                                 + visctermaniso*(dvx*dvx + dvy*dvy + dvz*dvz))
             if (usej) then
                grgrkernj = -2.*grkernj*rij1
                gradpj = gradpj + visctermisoj*projv*grgrkernj
                projsxj = projsxj + visctermanisoj*dvx*grgrkernj
                projsyj = projsyj + visctermanisoj*dvy*grgrkernj
                projszj = projszj + visctermanisoj*dvz*grgrkernj
             endif
          endif

          !--terms used in force
          gradp = gradpi + gradpj
          projsx = projsxi + projsxj
          projsy = projsyi + projsyj
          projsz = projszi + projszj

          fsum(ifxi) = fsum(ifxi) - runix*(gradp + fgrav) - projsx
          fsum(ifyi) = fsum(ifyi) - runiy*(gradp + fgrav) - projsy
          fsum(ifzi) = fsum(ifzi) - runiz*(gradp + fgrav) - projsz
          fsum(ipot) = fsum(ipot) + pmassj*phii ! no need to symmetrise (see PM07)
          if (icooling == 9) then
             gradpx = gradpx + runix*(gradP_cooli + gradP_coolj)
             gradpy = gradpy + runiy*(gradP_cooli + gradP_coolj)
             gradpz = gradpz + runiz*(gradP_cooli + gradP_coolj)
          endif

          !--calculate divv for use in du, h prediction, av switch etc.
          fsum(idrhodti) = fsum(idrhodti) + projv*grkerni

          if (maxvxyzu >= 4) then
             !--viscous heating
             fsum(idudtdissi) = fsum(idudtdissi) + dudtdissi + dudtresist
             !--energy dissipation due to conductivity
             fsum(idendtdissi) = fsum(idendtdissi) + dendissterm
          endif

          !--add contribution to particle i's force
          if (mhd) then
             !--div B in symmetric form (for source term subtraction)
             fsum(idivBsymi) = fsum(idivBsymi) + divBsymmterm
             fsum(idBevolxi) = fsum(idBevolxi) + dBevolx
             fsum(idBevolyi) = fsum(idBevolyi) + dBevoly
             fsum(idBevolzi) = fsum(idBevolzi) + dBevolz
             !--div B in difference form for dpsi/dt evolution
             fsum(idivBdiffi) = fsum(idivBdiffi) + divBdiffterm
          endif

          if (do_radiation) then
             radDFWj = 0.
             if (usej) then
                radFj(1:3) = radprop(ifluxx:ifluxz,j)
                radkappaj = radprop(ikappa,j)
                radenj = rad(iradxi,j)
                radRj = get_rad_R(rhoj,radenj,radFj,radkappaj)

                if (limit_radiation_flux) then
                   radlambdaj = (2. + radRj)/(6. + 3*radRj + radRj*radRj)
                else
                   radlambdaj = 1./3.
                endif

                radDj = c_code*radlambdaj/radkappaj/rhoj

                radDFWj = pmassj*radDj*grkernj*rho21j &
                          *(radFj(1)*runix + radFj(2)*runiy + radFj(3)*runiz)
             endif
             radFi(1:3) = xpartveci(iradfxi:iradfzi)
             radkappai = xpartveci(iradkappai)
             radeni = xpartveci(iradxii)
             radlambdai = xpartveci(iradlambdai)

             radDi = c_code*radlambdai/radkappai/rhoi

             ! TWO FIRST DERIVATIVES !
             radDFWi = pmassj*radDi*grkerni*rho21i*&
                (radFi(1)*runix + radFi(2)*runiy + radFi(3)*runiz)

             fsum(idradi) = fsum(idradi) + (radDFWi + radDFWj)
             ! BROOKSHAW !
             ! fsum(idradi) = fsum(idradi) + pmassj/rhoi/rhoj*&
             !                                ((4*radDi*radDj)/(radDi+radDj))*&
             !                                (rhoi*radeni-rhoj*radenj)*grkerni*rij1
          endif

          if (use_dust .and. use_dustfrac) then
             tsj = 0.
             do l=1,ndustsmall
                ! get stopping time - for one fluid dust we do not know deltav, but it is small by definition
                if (use_dustgrowth) then !- only work for ndustsmall=1 though
                   if (use_porosity) then
                      call get_ts(idrag,l,get_size(grainmassi,graindensi,filfaci),&
                      graindensi*filfaci,rhogasj,rhoj*dustfracjsum,spsoundj,0.,tsj(l),iregime)
                   else
                      call get_ts(idrag,l,get_size(grainmassi,graindensi),&
                      graindensi,rhogasj,rhoj*dustfracjsum,spsoundj,0.,tsj(l),iregime)
                   endif
                else
                   call get_ts(idrag,l,grainsize(l),graindens(l),rhogasj,rhoj*dustfracjsum,spsoundj,0.,tsj(l),iregime)
                endif
             enddo
             if (ilimitdustflux) tsj(:)   = min(tsj(:),hj/spsoundj) ! flux limiter from Ballabio et al. (2018)
             epstsj   = sum(dustfracj(:)*tsj(:))

             !! Check that weighted sums of Tsj and tilde(Tsj) are equal (see Hutchison et al. 2017)
             !if (ndustsmall>1) then
             !   if (abs(sum(dustfracj*tsj) - sum(dustfracj*(tsj-epstsj))/(1. - dustfracjsum)) > 1e-14) &
             !      print*,'Stop! tsj or epstsj in force is incorrect!'
             !else
             !   if (abs(tsj(1) - (tsj(1)-epstsj)/(1.-dustfracj(1)))>1e-14) &
             !      print*,'Warning! error in tsj is = ',tsj(1)-(tsj(1)-epstsj)/(1.-dustfracjsum)
             !endif

             do l=1,ndustsmall
                if (dustfraci(l) > 0. .or. dustfracj(l) > 0.) then
                   ! define average of kernels
                   grkernav = 0.5*(grkerni + grkernj)

                   !--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
                   dustfracterms(l) = pmassj*sqrtrhodustfracj(l)*rho1j     &
                                      *((tsi(l)-epstsi)*(1.-dustfraci(l))/(1.-dustfracisum)   &
                                        +(tsj(l)-epstsj)*(1.-dustfracj(l))/(1.-dustfracjsum)) &
                                      *(pri - prj)*grkernav*rij1

                   fsum(iddustevoli+(l-1)) = fsum(iddustevoli+(l-1)) - dustfracterms(l)

                   !--sqrt(rho*epsilon) method and sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
                   if (maxvxyzu >= 4) fsum(idudtdusti+(l-1)) = fsum(idudtdusti+(l-1)) - sqrtrhodustfraci(l)*dustfracterms(l)*denij
                endif
                ! Equation 270 in Phantom paper
                if (dustfraci(l) < 1.) then
                   term = tsi(l)/(1. - dustfracisum)*pmassj*(pro2i*grkerni + pro2j*grkernj)
                   fsum(ideltavxi+(l-1)) = fsum(ideltavxi+(l-1)) + term*runix
                   fsum(ideltavyi+(l-1)) = fsum(ideltavyi+(l-1)) + term*runiy
                   fsum(ideltavzi+(l-1)) = fsum(ideltavzi+(l-1)) + term*runiz
                endif
             enddo
          endif

       else !ifgas
          !
          !  gravity between particles of different types, or between gas pairs that are hidden by a sink
          !
          fsum(ifxi) = fsum(ifxi) - fgrav*runix
          fsum(ifyi) = fsum(ifyi) - fgrav*runiy
          fsum(ifzi) = fsum(ifzi) - fgrav*runiz
          fsum(ipot) = fsum(ipot) + pmassj*phii ! no need to symmetrise (see PM07)

          !
          ! gas-dust: compute drag terms
          !
          if (use_dust .and. idrag>0) then
             if (iamgasi .and. iamdustj .and. icut_backreaction==0) then
                projvstar = projv
                if (irecon >= 0) call reconstruct_dv(projv,dx,dy,dz,runix,runiy,runiz,dvdxi,dvdxj,projvstar,irecon)
                dv2 = projvstar**2 ! dvx*dvx + dvy*dvy + dvz*dvz
                if (q2i < q2j) then
                   wdrag = wkern_drag(q2i,qi)*hi21*hi1*cnormk_drag
                else
                   wdrag = wkern_drag(q2j,qj)*hj21*hj1*cnormk_drag
                endif
                if (use_dustgrowth) then
                   if (use_porosity) then
                      call get_ts(idrag,1,get_size(dustprop(1,j),dustprop(2,j),filfac(j)),&
                      dustprop(2,j)*filfac(j),rhoi,rhoj,spsoundi,dv2,tsijtmp,iregime)
                   else
                      call get_ts(idrag,1,get_size(dustprop(1,j),dustprop(2,j)),&
                      dustprop(2,j),rhoi,rhoj,spsoundi,dv2,tsijtmp,iregime)
                   endif
                else
                   !--the following works for large grains only (not hybrid large and small grains)
                   idusttype = iamtypej - idust + 1
                   call get_ts(idrag,idusttype,grainsize(idusttype),graindens(idusttype),rhoi,rhoj,spsoundi,dv2,tsijtmp,iregime)
                endif
                ndrag = ndrag + 1
                if (iregime > 2)  nstokes = nstokes + 1
                if (iregime == 2) nsuper = nsuper + 1
                if (drag_implicit) then
                   if (ind_timesteps) then
                      call get_drag_terms(tsijtmp,dti,sdrag1,sdrag2,fxi_drag,fyi_drag,fzi_drag,&
                                       runix,runiy,runiz,fxyz_drag(:,j),projf_drag)
                   else
                      call get_drag_terms(tsijtmp,dt,sdrag1,sdrag2,fxi_drag,fyi_drag,fzi_drag,&
                                       runix,runiy,runiz,fxyz_drag(:,j),projf_drag)
                   endif
                   ! implicit drag, slightly modified version of Loren-Aguilar & Bate 2014, 2015
                   call get_drag_terms(tsijtmp,dt,sdrag1,sdrag2,fxi_drag,fyi_drag,fzi_drag,&
                                       runix,runiy,runiz,fxyz_drag(:,j),projf_drag)
                   dragterm = 3.*pmassj*(sdrag1*projvstar - sdrag2*projf_drag)/(rhoi + rhoj)*wdrag
                else
                   ! explicit drag, with timestep condition
                   dragterm = 3.*pmassj/((rhoi + rhoj)*tsijtmp)*projvstar*wdrag
                endif
                ts_min = min(ts_min,tsijtmp)
                !store acceleration without drag here fr the next iteration
                fsum(ifdragxi) = fsum(ifdragxi) - dragterm*runix
                fsum(ifdragyi) = fsum(ifdragyi) - dragterm*runiy
                fsum(ifdragzi) = fsum(ifdragzi) - dragterm*runiz
                if (maxvxyzu >= 4) then
                   !--energy dissipation due to drag
                   dragheating = dragterm*projv
                   fsum(idudtdissi) = fsum(idudtdissi) + dragheating
                endif
             elseif (iamdusti .and. iamgasj) then
                projvstar = projv
                if (irecon >= 0) call reconstruct_dv(projv,dx,dy,dz,runix,runiy,runiz,dvdxi,dvdxj,projvstar,irecon)
                dv2 = projvstar**2 !dvx*dvx + dvy*dvy + dvz*dvz
                if (q2i < q2j) then
                   wdrag = wkern_drag(q2i,qi)*hi21*hi1*cnormk_drag
                else
                   wdrag = wkern_drag(q2j,qj)*hj21*hj1*cnormk_drag
                endif
                if (use_dustgrowth) then
                   if (use_porosity) then
                      call get_ts(idrag,1,get_size(grainmassi,graindensi,filfaci),&
                      graindensi*filfaci,rhoj,rhoi,spsoundj,dv2,tsijtmp,iregime)
                   else
                      call get_ts(idrag,1,get_size(grainmassi,graindensi),graindensi,rhoj,rhoi,spsoundj,dv2,tsijtmp,iregime)
                   endif
                   if (q2i < q2j) then
                      winter = wkern(q2i,qi)*hi21*hi1*cnormk
                   else
                      winter = wkern(q2j,qj)*hj21*hj1*cnormk
                   endif
                   !--following quantities are weighted by mass rather than mass/density
                   fsum(idensgasi) = fsum(idensgasi) + pmassj*winter
                   fsum(idvix)     = fsum(idvix)     + pmassj*dvx*winter
                   fsum(idviy)     = fsum(idviy)     + pmassj*dvy*winter
                   fsum(idviz)     = fsum(idviz)     + pmassj*dvz*winter
                   fsum(icsi)      = fsum(icsi)      + pmassj*spsoundj*winter
                else
                   !--the following works for large grains only (not hybrid large and small grains)
                   idusttype = iamtypei - idust + 1
                   call get_ts(idrag,idusttype,grainsize(idusttype),graindens(idusttype),rhoj,rhoi,spsoundj,dv2,tsijtmp,iregime)
                endif
                if (drag_implicit) then
                   ! implicit drag, slightly modified version of Loren-Aguilar & Bate 2014, 2015
                   if (ind_timesteps) then
                      call get_drag_terms(tsijtmp,dti,sdrag1,sdrag2,fxi_drag,fyi_drag,fzi_drag,&
                                       runix,runiy,runiz,fxyz_drag(:,j),projf_drag)
                   else
                      call get_drag_terms(tsijtmp,dt,sdrag1,sdrag2,fxi_drag,fyi_drag,fzi_drag,&
                                       runix,runiy,runiz,fxyz_drag(:,j),projf_drag)
                   endif
                   dragterm = 3.*pmassj*(sdrag1*projvstar - sdrag2*projf_drag)/(rhoi + rhoj)*wdrag
                else
                   ! explicit drag, with timestep condition
                   dragterm = 3.*pmassj/((rhoi + rhoj)*tsijtmp)*projvstar*wdrag
                endif
                ts_min = min(ts_min,tsijtmp)
                ndrag = ndrag + 1
                if (iregime > 2)  nstokes = nstokes + 1
                if (iregime == 2) nsuper = nsuper + 1
                fsum(ifdragxi) = fsum(ifdragxi) - dragterm*runix ! + because projv is opposite
                fsum(ifdragyi) = fsum(ifdragyi) - dragterm*runiy
                fsum(ifdragzi) = fsum(ifdragzi) - dragterm*runiz
             endif
          endif
       endif ifgas

       !--self gravity contribution to total energy equation
       if (gr .and. gravity .and. ien_type == ien_etotal) then
          fgravxi = fgravxi - runix*fgrav
          fgravyi = fgravyi - runiy*fgrav
          fgravzi = fgravzi - runiz*fgrav
       endif

#ifdef GRAVITY
    else !is_sph_neighbour
       !
       !--if particle is a trial neighbour, but not an SPH neighbour
       !  then compute the 1/r^2 force contribution
       !  (no softening here, as by definition we
       !   are outside the kernel radius)
       !
#ifdef FINVSQRT
       rij1 = finvsqrt(rij2)
#else
       rij1 = 1./sqrt(rij2)
#endif
       fgrav  = rij1*rij1*rij1
       if (maxphase==maxp) then
          iamtypej = iamtype(iphase(j))
       endif
       if (use_apr) then
          pmassj = aprmassoftype(iamtypej,apr_level(j))
       else
          pmassj = massoftype(iamtypej)
       endif
       phii   = -rij1
       fgravj = fgrav*pmassj
       fsum(ifxi) = fsum(ifxi) - dx*fgravj
       fsum(ifyi) = fsum(ifyi) - dy*fgravj
       fsum(ifzi) = fsum(ifzi) - dz*fgravj
       fsum(ipot) = fsum(ipot) + pmassj*phii

       !--self gravity contribution to total energy equation
       if (gr .and. gravity .and. ien_type == ien_etotal) then
          fgravxi = fgravxi - dx*fgravj
          fgravyi = fgravyi - dy*fgravj
          fgravzi = fgravzi - dz*fgravj
       endif
#endif
    endif is_sph_neighbour

 enddo loop_over_neighbours2

 if (icooling == 9) gradP_cool(i) = sqrt(gradpx*gradpx + gradpy*gradpy + gradpz*gradpz)

 if (gr .and. gravity .and. ien_type == ien_etotal) then
    fsum(idudtdissi) = fsum(idudtdissi) + vxi*fgravxi + vyi*fgravyi + vzi*fgravzi
 endif

end subroutine compute_forces

!----------------------------------------------------------------
!+
!  Internal subroutine that computes pressure and other derived
!  quantities necessary to get a force, given that we have rho.
!+
!----------------------------------------------------------------
subroutine get_stress(pri,spsoundi,rhoi,rho1i,xi,yi,zi, &
                 pmassi,Bxi,Byi,Bzi, &
                 pro2i,vwavei, &
                 sxxi,sxyi,sxzi,syyi,syzi,szzi,visctermiso,visctermaniso, &
                 realviscosity,divvi,bulkvisc,dvdx,stressmax, &
                 radPi)

 use dim,             only:maxdvdx,maxp
 use part,            only:mhd,strain_from_dvdx
 use viscosity,       only:shearfunc

 real,    intent(in)    :: pri,spsoundi,rhoi,rho1i,xi,yi,zi,pmassi
 real,    intent(in)    :: Bxi,Byi,Bzi
 real,    intent(out)   :: pro2i,vwavei
 real,    intent(out)   :: sxxi,sxyi,sxzi,syyi,syzi,szzi
 real,    intent(out)   :: visctermiso,visctermaniso
 logical, intent(in)    :: realviscosity
 real,    intent(in)    :: divvi,bulkvisc,stressmax
 real,    intent(in)    :: dvdx(9)
 real,    intent(in)    :: radPi

 real :: Bro2i,Brhoxi,Brhoyi,Brhozi
 real :: stressiso,term,graddivvcoeff,del2vcoeff,strain(6)
 real :: shearvisc,etavisc,valfven2i

 sxxi = 0.
 sxyi = 0.
 sxzi = 0.
 syyi = 0.
 syzi = 0.
 szzi = 0.
 visctermiso   = 0.
 visctermaniso = 0.
 stressiso     = 0.

 if (realviscosity) then
    !--get shear viscosity coefficient from function
    shearvisc = shearfunc(xi,yi,zi,spsoundi)
    etavisc   = rhoi*shearvisc
!
!--add physical viscosity terms to stress tensor
!  (construct S^ij/ rho^2 for use in the force equation)
!
    if (maxdvdx==maxp) then
       strain = strain_from_dvdx(dvdx)
       !--get stress (multiply by coefficient for use in second derivative)
       term = -shearvisc*pmassi*rho1i  ! shearvisc = eta/rho, so this is eta/rho**2
       sxxi = term*strain(1)
       sxyi = term*strain(2)
       sxzi = term*strain(3)
       syyi = term*strain(4)
       syzi = term*strain(5)
       szzi = term*strain(6)
       stressiso = (2./3.*shearvisc - bulkvisc)*divvi*rho1i + stressmax
    else
       graddivvcoeff = 0.5*(bulkvisc*rhoi + etavisc/3.)   ! 0.5 here is because we
       del2vcoeff    = 0.5*etavisc                   ! average between particle pairs

       !--construct isotropic and anisotropic terms from above
       visctermiso   = 2.5*graddivvcoeff*pmassi*rho1i*rho1i
       visctermaniso = (del2vcoeff - 0.5*graddivvcoeff)*pmassi*rho1i*rho1i
    endif
 endif

 if (mhd) then
!
!--construct useful terms based on the B-field
!
    Brhoxi = Bxi*rho1i
    Brhoyi = Byi*rho1i
    Brhozi = Bzi*rho1i

    Bro2i     = Brhoxi*Brhoxi + Brhoyi*Brhoyi + Brhozi*Brhozi
    valfven2i = Bro2i*rhoi
    vwavei    = sqrt(spsoundi*spsoundi + valfven2i)

    !--MHD terms in stress tensor
    sxxi  = sxxi - pmassi*Brhoxi*Brhoxi
    sxyi  = sxyi - pmassi*Brhoxi*Brhoyi
    sxzi  = sxzi - pmassi*Brhoxi*Brhozi
    syyi  = syyi - pmassi*Brhoyi*Brhoyi
    syzi  = syzi - pmassi*Brhoyi*Brhozi
    szzi  = szzi - pmassi*Brhozi*Brhozi
!
!--construct total isotropic pressure term (gas + magnetic + stress)
!
    pro2i = (pri + radPi)*rho1i*rho1i + stressiso + 0.5*Bro2i

 else
!
!--construct m*p/(rho^2 \Omega) in force equation using pressure
!
    pro2i  = (pri + radPi)*rho1i*rho1i + stressiso
    vwavei = spsoundi

    if (do_radiation) then
       vwavei = sqrt(vwavei*vwavei + 4.*radPi/(3.*rhoi))  ! Commercon et al. (2011)
    endif
 endif

end subroutine get_stress

!----------------------------------------------------------------

subroutine start_cell(cell,iphase,xyzh,vxyzu,gradh,divcurlv,divcurlB,dvdx,Bevol, &
                     dustfrac,dustprop,fxyz_drag,eta_nimhd,eos_vars,alphaind,stressmax,&
                     rad,radprop,dens,metrics,apr_level,dt)

 use io,        only:fatal
 use options,   only:alpha,use_dustfrac,limit_radiation_flux
 use dim,       only:maxp,ndivcurlv,ndivcurlB,maxdvdx,maxalpha,maxvxyzu,mhd,mhd_nonideal,&
                use_dustgrowth,gr,use_dust,ind_timesteps,use_apr
 use part,      only:iamgas,maxphase,rhoanddhdrho,igas,massoftype,get_partinfo,&
                     iohm,ihall,iambi,ndustsmall,iradP,igasP,ics,itemp,aprmassoftype
 use viscosity, only:irealvisc,bulkvisc
 use dust,      only:get_ts,idrag
 use options,   only:use_porosity
 use part,      only:grainsize,graindens,filfac
 use growth,    only:get_size
 use part,        only:ibin_old
 use timestep_ind,    only:get_dt
 use nicil,           only:nimhd_get_jcbcb
 use radiation_utils, only:get_rad_R
 type(cellforce),    intent(inout) :: cell
 integer(kind=1),    intent(in)    :: iphase(:)
 real,               intent(in)    :: xyzh(:,:)
 real,               intent(inout) :: vxyzu(:,:)
 real(kind=4),       intent(in)    :: gradh(:,:)
 real(kind=4),       intent(in)    :: divcurlv(:,:)
 real(kind=4),       intent(in)    :: divcurlB(:,:)
 real(kind=4),       intent(in)    :: dvdx(:,:)
 real,               intent(in)    :: Bevol(:,:)
 real,               intent(in)    :: dustfrac(:,:)
 real,               intent(in)    :: dustprop(:,:)
 real,               intent(in)    :: fxyz_drag(:,:)
 real,               intent(in)    :: eta_nimhd(:,:)
 real(kind=4),       intent(in)    :: alphaind(:,:)
 real,               intent(in)    :: stressmax
 real,               intent(in)    :: rad(:,:)
 real,               intent(inout) :: radprop(:,:)
 real,               intent(in)    :: dens(:)
 real,               intent(in)    :: metrics(:,:,:,:)
 real,               intent(in)    :: eos_vars(:,:)
 real,               intent(in)    :: dt
 integer(kind=1),    intent(in)    :: apr_level(:)
 real         :: radRi
 real         :: radPi

 real         :: divcurlvi(ndivcurlv)
 real         :: dvdxi(9),curlBi(3),jcbcbi(3),jcbi(3)
 real         :: hi,rhoi,rho1i,dhdrhoi,pmassi,eni
 real(kind=8) :: hi1
 real         :: dustfraci(maxdusttypes),dustfracisum,rhogasi,pro2i,pri,spsoundi,tempi
 real         :: sxxi,sxyi,sxzi,syyi,syzi,szzi,visctermiso,visctermaniso
 real         :: tstopi(maxdusttypes)
 real         :: Bxi,Byi,Bzi,B2i,Bi1
 real         :: vwavei,alphai
 integer      :: i,j,iamtypei,ip,ii,ia,ib,ic
 real         :: densi
 integer      :: iregime
 real         :: dti

 logical :: iactivei,iamgasi,iamdusti,realviscosity

 realviscosity = (irealvisc > 0)

 cell%npcell = 0
 over_parts: do ip = inoderange(1,cell%icell),inoderange(2,cell%icell)
    i = inodeparts(ip)

    if (i < 0) then
       cycle over_parts
    endif

    if (maxphase==maxp) then
       call get_partinfo(iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif
    if (.not.iactivei) then ! handles boundaries + case where first particle in cell is inactive
       cycle over_parts
    endif

    if (use_apr) then
       pmassi = aprmassoftype(iamtypei,apr_level(i))
    else
       pmassi = massoftype(iamtypei)
    endif

    hi = xyzh(4,i)
    if (hi < 0.) call fatal('force','negative smoothing length',i,var='h',val=hi)

    !
    !--compute density and related quantities from the smoothing length
    !
    call rhoanddhdrho(hi,hi1,rhoi,rho1i,dhdrhoi,pmassi)
    !
    !--velocity gradients, used for reconstruction and physical viscosity
    !
    if (maxdvdx==maxp) dvdxi(:) = dvdx(:,i)

    if (iamgasi) then
       if (ndivcurlv >= 1) divcurlvi(:) = real(divcurlv(:,i),kind=kind(divcurlvi))
       if (maxvxyzu >= 4) then
          eni = vxyzu(4,i)
       else
          eni = 0.0
       endif

       !
       ! one-fluid dust properties
       !
       if (use_dustfrac) then
          dustfraci(:) = dustfrac(:,i)
          dustfracisum = sum(dustfraci)
          rhogasi      = rhoi*(1. - dustfracisum)
          do j=1,ndustsmall
             if (dustfraci(j) > 1. .or. dustfraci(j) < 0.) call fatal('force','invalid eps',var='dustfrac',val=dustfraci(j))
          enddo
          if (dustfracisum > 1.) call fatal('force','invalid eps',var='sum of dustfrac',val=dustfracisum)
       else
          dustfraci(:) = 0.
          dustfracisum = 0.
          rhogasi      = rhoi
       endif

       if (mhd) then
          Bxi = Bevol(1,i) * rhoi ! B/rho -> B (conservative to primitive)
          Byi = Bevol(2,i) * rhoi
          Bzi = Bevol(3,i) * rhoi
          if (mhd_nonideal) then
             B2i = Bxi**2 + Byi**2 + Bzi**2
             if (B2i > 0.0) then
                Bi1 = 1.0/sqrt(B2i)
             else
                Bi1 = 0.0
             endif
             curlBi = divcurlB(2:4,i)
             call nimhd_get_jcbcb(jcbcbi,jcbi,curlBi,Bxi,Byi,Bzi,Bi1)
          endif
       endif

       if (gr) densi = dens(i)
       pri = eos_vars(igasP,i)
       spsoundi = eos_vars(ics,i)
       tempi = eos_vars(itemp,i)
       radPi = 0.
       if (do_radiation) radPi = radprop(iradP,i)
       !
       ! calculate terms required in the force evaluation
       !
       call get_stress(pri,spsoundi,rhoi,rho1i, &
                  xyzh(1,i),xyzh(2,i),xyzh(3,i), &
                  pmassi, &
                  Bxi,Byi,Bzi, &
                  pro2i, &
                  vwavei,sxxi,sxyi,sxzi,syyi,syzi,szzi, &
                  visctermiso,visctermaniso,realviscosity,divcurlvi(1),bulkvisc,dvdxi,stressmax, &
                  radPi)

       !
       ! get stopping time - for one fluid dust we don't know deltav, but as small by definition we assume=0
       !
       tstopi = 0.
       if (use_dust .and. use_dustfrac .and. iamgasi) then
          do j=1,ndustsmall
             if (use_dustgrowth) then
                if (use_porosity) then
                   call get_ts(idrag,j,get_size(dustprop(1,i),dustprop(2,i),filfac(j)),&
                   dustprop(2,i)*filfac(j),rhogasi,rhoi*dustfracisum,spsoundi,0.,tstopi(j),iregime)
                else
                   call get_ts(idrag,j,get_size(dustprop(1,i),dustprop(2,i)),&
                   dustprop(2,i),rhogasi,rhoi*dustfracisum,spsoundi,0.,tstopi(j),iregime)
                endif
             else
                call get_ts(idrag,j,grainsize(j),graindens(j),rhogasi,rhoi*dustfracisum,spsoundi,0.,tstopi(j),iregime)
             endif
          enddo
       endif

    else ! not a gas particle
       vwavei = 0.
       rhogasi = 0.
       pri = 0.
       pro2i = 0.
       sxxi = 0.
       sxyi = 0.
       sxzi = 0.
       syyi = 0.
       syzi = 0.
       szzi = 0.
       visctermiso = 0.
       visctermaniso = 0.
       dustfraci = 0.
       tempi = 0.
       densi = 0.
    endif

    !
    !-- cell packing
    !
    cell%npcell                                    = cell%npcell + 1
    cell%arr_index(cell%npcell)                    = ip
    cell%iphase(cell%npcell)                       = iphase(i)
    cell%xpartvec(ixi,cell%npcell)                 = xyzh(1,i)
    cell%xpartvec(iyi,cell%npcell)                 = xyzh(2,i)
    cell%xpartvec(izi,cell%npcell)                 = xyzh(3,i)
    cell%xpartvec(ihi,cell%npcell)                 = xyzh(4,i)
    cell%xpartvec(igradhi1,cell%npcell)            = gradh(1,i)
#ifdef GRAVITY
    cell%xpartvec(igradhi2,cell%npcell)            = gradh(2,i)
#endif
    cell%xpartvec(ivxi,cell%npcell)                = vxyzu(1,i)
    cell%xpartvec(ivyi,cell%npcell)                = vxyzu(2,i)
    cell%xpartvec(ivzi,cell%npcell)                = vxyzu(3,i)
    if (maxvxyzu >= 4) then
       cell%xpartvec(ieni,cell%npcell)             = vxyzu(4,i)
    endif
    if (mhd) then
       if (iamgasi) then
          cell%xpartvec(iBevolxi,cell%npcell)      = Bevol(1,i)
          cell%xpartvec(iBevolyi,cell%npcell)      = Bevol(2,i)
          cell%xpartvec(iBevolzi,cell%npcell)      = Bevol(3,i)
          cell%xpartvec(ipsi,cell%npcell)          = Bevol(4,i)

          cell%xpartvec(idivBi,  cell%npcell)      = divcurlB(1,i)
          cell%xpartvec(icurlBxi,cell%npcell)      = divcurlB(2,i)
          cell%xpartvec(icurlByi,cell%npcell)      = divcurlB(3,i)
          cell%xpartvec(icurlBzi,cell%npcell)      = divcurlB(4,i)
       else
          cell%xpartvec(iBevolxi:ipsi,cell%npcell) = 0. ! to avoid compiler warning
       endif
    endif
    if (use_apr) then
       cell%apr(cell%npcell)                        = apr_level(i)
    else
       cell%apr(cell%npcell)                        = 1
    endif

    alphai = alpha
    if (maxalpha==maxp) then
       alphai = alphaind(1,i)
    endif
    cell%xpartvec(ialphai,cell%npcell)            = alphai
    cell%xpartvec(irhoi,cell%npcell)              = rhoi
    cell%xpartvec(idustfraci:idustfraciend,cell%npcell) = dustfraci
    if (iamgasi) then
       cell%xpartvec(ivwavei,cell%npcell)         = vwavei
       cell%xpartvec(irhogasi,cell%npcell)        = rhogasi
       cell%xpartvec(ipri,cell%npcell)            = pri
       cell%xpartvec(ispsoundi,cell%npcell)       = spsoundi
       cell%xpartvec(itempi,cell%npcell)          = tempi
       cell%xpartvec(isxxi,cell%npcell)           = sxxi
       cell%xpartvec(isxyi,cell%npcell)           = sxyi
       cell%xpartvec(isxzi,cell%npcell)           = sxzi
       cell%xpartvec(isyyi,cell%npcell)           = syyi
       cell%xpartvec(isyzi,cell%npcell)           = syzi
       cell%xpartvec(iszzi,cell%npcell)           = szzi
       cell%xpartvec(ivisctermisoi,cell%npcell)   = visctermiso
       cell%xpartvec(ivisctermanisoi,cell%npcell) = visctermaniso
       cell%xpartvec(ipri,cell%npcell)            = pri
       cell%xpartvec(ipro2i,cell%npcell)          = pro2i
       if (use_dustfrac) then
          cell%xpartvec(itstop:itstopend,cell%npcell) = tstopi
       endif
       if (use_dust .and. ind_timesteps) then
          dti = get_dt(dt,ibin_old(i))
          cell%xpartvec(idti,cell%npcell) = dti
       endif

       if (do_radiation) then
          radRi = get_rad_R(rhoi,rad(iradxi,i),radprop(ifluxx:ifluxz,i),radprop(ikappa,i))
          cell%xpartvec(iradxii,cell%npcell)         = rad(iradxi,i)
          cell%xpartvec(iradfxi:iradfzi,cell%npcell) = radprop(ifluxx:ifluxz,i)
          cell%xpartvec(iradkappai,cell%npcell)      = radprop(ikappa,i)
          if (limit_radiation_flux) then
             cell%xpartvec(iradlambdai,cell%npcell) = (2. + radRi)/(6. + 3.*radRi + radRi*radRi)
          else
             cell%xpartvec(iradlambdai,cell%npcell) = 1./3.
          endif
          cell%xpartvec(iradrbigi,cell%npcell)       = radRi
       endif

    endif
    if (mhd_nonideal) then
       cell%xpartvec(ietaohmi,cell%npcell)        = eta_nimhd(iohm,i)
       cell%xpartvec(ietahalli,cell%npcell)       = eta_nimhd(ihall,i)
       cell%xpartvec(ietaambii,cell%npcell)       = eta_nimhd(iambi,i)
       cell%xpartvec(ijcbcbxi,cell%npcell)        = jcbcbi(1)
       cell%xpartvec(ijcbcbyi,cell%npcell)        = jcbcbi(2)
       cell%xpartvec(ijcbcbzi,cell%npcell)        = jcbcbi(3)
       cell%xpartvec(ijcbxi,cell%npcell)          = jcbi(1)
       cell%xpartvec(ijcbyi,cell%npcell)          = jcbi(2)
       cell%xpartvec(ijcbzi,cell%npcell)          = jcbi(3)
    endif

    if (gr) then
       cell%xpartvec(idensGRi,cell%npcell)        = densi
       ii = imetricstart
       do ic = 1,2
          do ib = 1,4
             do ia = 1,4
                cell%xpartvec(ii,cell%npcell)     = metrics(ia,ib,ic,i)
                ii = ii + 1
             enddo
          enddo
       enddo
    endif

    if (use_dustgrowth) then
       cell%xpartvec(igrainmassi,cell%npcell)     = dustprop(1,i)
       cell%xpartvec(igraindensi,cell%npcell)     = dustprop(2,i)
       cell%xpartvec(ifilfaci,cell%npcell)        = filfac(i)
    endif
    if (use_dust) then
       cell%xpartvec(ifxi_drag,cell%npcell)       = fxyz_drag(1,i)
       cell%xpartvec(ifyi_drag,cell%npcell)       = fxyz_drag(2,i)
       cell%xpartvec(ifzi_drag,cell%npcell)       = fxyz_drag(3,i)
    endif
    cell%xpartvec(idvxdxi:idvzdzi,cell%npcell)    = dvdx(1:9,i)
 enddo over_parts

end subroutine start_cell

subroutine compute_cell(cell,listneigh,nneigh,Bevol,xyzh,vxyzu,fxyzu, &
                        iphase,divcurlv,divcurlB,alphaind,eta_nimhd, eos_vars, &
                        dustfrac,dustprop,fxyz_drag,gradh,ibinnow_m1,ibin_wake,stressmax,xyzcache,&
                        rad,radprop,dens,metrics,apr_level,dt)
 use io,          only:error,id
 use dim,         only:maxvxyzu,use_apr
 use options,     only:beta,alphau,alphaB,iresistive_heating
 use part,        only:get_partinfo,iamgas,mhd,igas,maxphase,massoftype,aprmassoftype
 use viscosity,   only:irealvisc,bulkvisc

 type(cellforce), intent(inout)  :: cell

 integer,         intent(in)     :: listneigh(:)
 integer,         intent(in)     :: nneigh
 real,            intent(in)     :: Bevol(:,:)
 real,            intent(in)     :: xyzh(:,:)
 real,            intent(inout)  :: vxyzu(:,:)
 real,            intent(in)     :: fxyzu(:,:)
 integer(kind=1), intent(in)     :: iphase(:)
 real(kind=4),    intent(in)     :: divcurlv(:,:)
 real(kind=4),    intent(in)     :: divcurlB(:,:)
 real(kind=4),    intent(in)     :: alphaind(:,:)
 real,            intent(in)     :: eta_nimhd(:,:)
 real,            intent(in)     :: dustfrac(:,:)
 real,            intent(in)     :: dustprop(:,:)
 real,            intent(in)     :: fxyz_drag(:,:)
 real,            intent(in)     :: eos_vars(:,:)
 real(kind=4),    intent(in)     :: gradh(:,:)
 integer(kind=1), intent(inout)  :: ibin_wake(:)
 integer(kind=1), intent(in)     :: ibinnow_m1
 real,            intent(in)     :: stressmax
 real,            intent(in)     :: xyzcache(:,:)
 real,            intent(in)     :: rad(:,:)
 real,            intent(inout)  :: radprop(:,:)
 real,            intent(in)     :: dens(:),metrics(:,:,:,:)
 real,            intent(in)     :: dt
 integer(kind=1), intent(in)     :: apr_level(:)

 real                            :: hi
 real(kind=8)                    :: hi1,hi21,hi31,hi41
 real(kind=8)                    :: gradhi,gradsofti
 real                            :: pmassi

 integer                         :: iamtypei

 logical                         :: iactivei
 logical                         :: iamgasi
 logical                         :: iamdusti

 logical                         :: realviscosity
 logical                         :: useresistiveheat
 logical                         :: ignoreself

 integer                         :: i,ip

 realviscosity    = (irealvisc > 0)
 useresistiveheat = (iresistive_heating > 0)

 over_parts: do ip = 1,cell%npcell

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(ip),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif

    if (.not.iactivei) then ! handles case where first particle in cell is inactive
       cycle over_parts     ! also boundary particles are inactive
    endif

    i = inodeparts(cell%arr_index(ip))

    if (use_apr) then
       pmassi = aprmassoftype(iamtypei,cell%apr(ip))
    else
       pmassi = massoftype(iamtypei)
    endif

    hi    = cell%xpartvec(ihi,ip)
    hi1   = 1./hi
    hi21  = hi1*hi1
    hi31  = hi1*hi21
    hi41  = hi21*hi21

    if (cell%xpartvec(igradhi1,ip) > 0.) then
       gradhi = cell%xpartvec(igradhi1,ip)
    else
       call error('force','stored gradh is zero, resetting to 1')
       gradhi = 1.
    endif
#ifdef GRAVITY
    gradsofti = cell%xpartvec(igradhi2,ip)
#endif

    !
    !--loop over current particle's neighbours (includes self)
    !
    ignoreself = (cell%owner == id)
    call compute_forces(i,iamgasi,iamdusti,cell%xpartvec(:,ip),hi,hi1,hi21,hi41,gradhi,gradsofti, &
                         beta, &
                         pmassi,listneigh,nneigh,xyzcache,cell%fsums(:,ip),cell%vsigmax(ip), &
                         .true.,realviscosity,useresistiveheat, &
                         xyzh,vxyzu,Bevol,cell%iphase(ip),iphase,massoftype, &
                         divcurlB,eta_nimhd,eos_vars, &
                         dustfrac,dustprop,fxyz_drag,gradh,divcurlv,alphaind, &
                         alphau,alphaB,bulkvisc,stressmax, &
                         cell%ndrag,cell%nstokes,cell%nsuper,cell%tsmin(ip),ibinnow_m1,ibin_wake,cell%ibinneigh(ip), &
                         ignoreself,rad,radprop,dens,metrics,apr_level,dt)

 enddo over_parts

end subroutine compute_cell

subroutine finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,dvdx,&
                                         divBsymm,divcurlv,dBevol,ddustevol,deltav,dustgasprop,fxyz_drag,fext,dragreg, &
                                         filfac,dtcourant,dtforce,dtvisc,dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi, &
#ifdef IND_TIMESTEPS
                                         nbinmaxnew,nbinmaxstsnew,ncheckbin, &
                                         ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd, &
                                         ndtvisc,ndtohm,ndthall,ndtambi,ndtdust,ndtrad,ndtclean, &
                                         dtitmp,dtrat, &
                                         dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean, &
                                         dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax,  &
                                         dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean ,dtdustfacmean, &
                                         dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax  ,dtdustfacmax,  &
                                         dtradfacmean ,dtcleanfacmean, &
                                         dtradfacmax  ,dtcleanfacmax, &
#endif
                                         ndustres,dustresfacmax,dustresfacmean,&
                                         rad,drad,radprop,dtrad)

 use io,             only:fatal,warning
 use dim,            only:mhd,mhd_nonideal,lightcurve,use_dust,maxdvdx,use_dustgrowth,gr,use_krome,&
                          store_dust_temperature,do_nucleation,update_muGamma,h2chemistry,use_apr
 use eos,            only:ieos,iopacity_type
 use options,        only:alpha,ipdv_heating,ishock_heating,psidecayfac,overcleanfac, &
                          use_dustfrac,damp,icooling,implicit_radiation
 use part,           only:rhoanddhdrho,iboundary,igas,maxphase,maxvxyzu,nptmass,xyzmh_ptmass,eos_vars, &
                          massoftype,get_partinfo,tstop,strain_from_dvdx,ithick,iradP,sinks_have_heating,&
                          luminosity,nucleation,idK2,idkappa,dust_temp,pxyzu,ndustsmall,imu,&
                          igamma,aprmassoftype
 use cooling,        only:energ_cooling,cooling_in_step
 use ptmass_heating, only:energ_sinkheat
 use dust,           only:drag_implicit
#ifdef IND_TIMESTEPS
 use part,           only:ibin
 use timestep_ind,   only:get_newbin,check_dtmin
 use timestep_sts,   only:sts_it_n,ibin_sts
#endif
 use viscosity,      only:bulkvisc,dt_viscosity,irealvisc,shearfunc
 use kernel,         only:kernel_softening
 use linklist,       only:get_distance_from_centre_of_mass
 use kdtree,         only:expand_fgrav_in_taylor_series
 use nicil,          only:nicil_get_dudt_nimhd,nicil_get_dt_nimhd
 use timestep,       only:C_cour,C_cool,C_force,C_rad,C_ent,bignumber,dtmax
 use timestep_sts,   only:use_sts
 use units,          only:unit_ergg,unit_density,get_c_code
 use eos_shen,       only:eos_shen_get_dTdu
 use metric_tools,   only:unpack_metric
 use utils_gr,       only:get_u0
 use io,             only:error
 use growth,         only:get_size
 use dust,           only:idrag,get_ts
 use physcon,        only:fourpi
 use options,        only:use_porosity
 use part,           only:Omega_k
 use io,             only:warning
 use physcon,        only:c,kboltz
 use eos_stamatellos, only:duSPH
 integer,            intent(in)    :: icall
 type(cellforce),    intent(inout) :: cell
 real,               intent(inout) :: fxyzu(:,:)
 real,               intent(in)    :: xyzh(:,:)
 real,               intent(inout) :: vxyzu(:,:)
 real,               intent(in)    :: dt
 real(kind=4),       intent(in)    :: dvdx(:,:)
 real(kind=4),       intent(out)   :: poten(:)
 real(kind=4),       intent(out)   :: divBsymm(:)
 real(kind=4),       intent(out)   :: divcurlv(:,:)
 real,               intent(out)   :: dBevol(:,:)
 real,               intent(out)   :: ddustevol(:,:)
 real,               intent(out)   :: deltav(:,:,:)
 real,               intent(out)   :: dustgasprop(:,:)
 real,               intent(in)    :: filfac(:)
 integer,            intent(out)   :: dragreg(:)
 real,               intent(inout) :: fxyz_drag(:,:)
 real,               intent(in)    :: fext(:,:)
 real,               intent(inout) :: dtcourant,dtforce,dtvisc
 real,               intent(inout) :: dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi
#ifdef IND_TIMESTEPS
 integer,            intent(inout) :: nbinmaxnew,nbinmaxstsnew,ncheckbin
 integer,            intent(inout) :: ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd
 integer,            intent(inout) :: ndtvisc,ndtohm,ndthall,ndtambi,ndtdust,ndtrad,ndtclean
 real,               intent(inout) :: dtitmp,dtrat
 real,               intent(inout) :: dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean ,dtdragdfacmean,dtcoolfacmean
 real,               intent(inout) :: dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax  ,dtdragdfacmax ,dtcoolfacmax
 real,               intent(inout) :: dtviscfacmean,dtohmfacmean  ,dthallfacmean ,dtambifacmean ,dtdustfacmean
 real,               intent(inout) :: dtviscfacmax ,dtohmfacmax   ,dthallfacmax  ,dtambifacmax  ,dtdustfacmax
 real,               intent(inout) :: dtradfacmean ,dtcleanfacmean
 real,               intent(inout) :: dtradfacmax  ,dtcleanfacmax
#endif
 integer,            intent(inout) :: ndustres
 real,               intent(inout) :: dustresfacmean,dustresfacmax
 real,               intent(in)    :: rad(:,:),radprop(:,:)
 real,               intent(out)   :: drad(:,:)
 real,               intent(inout) :: dtrad
 real    :: c_code,dtradi,radlambdai,radkappai
 real    :: xpartveci(maxxpartveciforce),fsum(maxfsum)
 real    :: rhoi,rho1i,rhogasi,hi,hi1,pmassi,tempi,gammai
 real    :: Bxyzi(3),curlBi(3),dvdxi(9),straini(6)
 real    :: xi,yi,zi,B2i,f2i,divBsymmi,betai,frac_divB,divBi,vcleani
 real    :: pri,spsoundi,drhodti,divvi,shearvisc,fac,pdv_work
 real    :: psii,dtau
 real    :: eni,dudtnonideal
 real    :: dustfraci(maxdusttypes),dustfracisum
 real    :: tstopi(maxdusttypes),tseff,dtdustdenom
 real    :: etaambii,etahalli,etaohmi
 real    :: vsigmax,vwavei,fxyz4
 real    :: dTdui,dTdui_cgs,rho_cgs
 real    :: dudt_radi
#ifdef GRAVITY
 real    :: potensoft0,dum,dx,dy,dz,fxi,fyi,fzi,poti,epoti
#endif
 real    :: vsigdtc,dtc,dtf,dti,dtcool,dtdiffi,ts_min,dtent
 real    :: dtohmi,dtambii,dthalli,dtvisci,dtdrag,dtdusti,dtclean
 integer :: iamtypei
 logical :: iactivei,iamgasi,iamdusti,realviscosity
#ifdef IND_TIMESTEPS
 integer(kind=1)       :: ibin_neighi
 logical               :: allow_decrease,dtcheck
 character(len=16)     :: dtchar
#endif
 real    :: tstopint,gmassi,gdensi
 integer :: ireg
 integer               :: ip,i
 real                  :: densi, vxi,vyi,vzi,u0i,dudtcool,dudtheat
 real                  :: posi(3),veli(3),gcov(0:3,0:3),metrici(0:3,0:3,2)
 integer               :: ii,ia,ib,ic,ierror
 eni = 0.
 realviscosity = (irealvisc > 0)

 over_parts: do ip = 1,cell%npcell

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(ip),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif

    if (.not.iactivei) then ! handles case where first particle in cell is inactive
       cycle over_parts
    endif

    if (use_apr) then
       pmassi = aprmassoftype(iamtypei,cell%apr(ip))
    else
       pmassi = massoftype(iamtypei)
    endif

    i = inodeparts(cell%arr_index(ip))

    fsum(:)       = cell%fsums(:,ip)
    xpartveci(:)  = cell%xpartvec(:,ip)
#ifdef IND_TIMESTEPS
    ibin_neighi   = cell%ibinneigh(ip)
#endif

    dtc     = dtmax
    dtf     = bignumber
    dtcool  = bignumber
    dtvisci = bignumber
    dtohmi  = bignumber
    dthalli = bignumber
    dtambii = bignumber
    dtdiffi = bignumber
    dtclean = bignumber
    dtdusti = bignumber
    dtdrag  = bignumber
    dtradi  = bignumber
    dtent   = bignumber

    xi         = xpartveci(ixi)
    yi         = xpartveci(iyi)
    zi         = xpartveci(izi)
    hi         = xpartveci(ihi)
    hi1        = 1./hi
    spsoundi   = xpartveci(ispsoundi)
    tempi      = xpartveci(itempi)
    vsigmax    = cell%vsigmax(ip)
    ts_min     = cell%tsmin(ip)
    tstopi     = 0.
    dustfraci  = 0.
    dustfracisum = 0.
    gammai = eos_vars(igamma,i)

    vxi = xpartveci(ivxi)
    vyi = xpartveci(ivyi)
    vzi = xpartveci(ivzi)

    u0i = 1.
    if (gr) then
       veli = (/vxi,vyi,vzi/)
       posi = (/xi,yi,zi/)

       densi = xpartveci(idensGRi)
       ii = imetricstart
       do ic = 1,2
          do ib = 0,3
             do ia = 0,3
                metrici(ia,ib,ic) = xpartveci(ii)
                ii = ii + 1
             enddo
          enddo
       enddo

       call unpack_metric(metrici,gcov=gcov)
       call get_u0(gcov,veli,u0i,ierror)
       if (ierror > 0) call error('get_u0 in force','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
    endif

    if (iamgasi) then
       rhoi    = xpartveci(irhoi)
       rho1i   = 1./rhoi
       rhogasi = xpartveci(irhogasi)
       pri     = xpartveci(ipri)
       vwavei  = xpartveci(ivwavei)

       if (mhd) then
          Bxyzi(1) = xpartveci(iBevolxi) * rhoi
          Bxyzi(2) = xpartveci(iBevolyi) * rhoi
          Bxyzi(3) = xpartveci(iBevolzi) * rhoi
          B2i      = Bxyzi(1)**2 + Bxyzi(2)**2 + Bxyzi(3)**2
          divBi    = xpartveci(idivBi)
       endif

       if (mhd_nonideal) then
          curlBi(1) = xpartveci(icurlBxi)
          curlBi(2) = xpartveci(icurlByi)
          curlBi(3) = xpartveci(icurlBzi)
          etaohmi   = xpartveci(ietaohmi)
          etaambii  = xpartveci(ietaambii)
          etahalli  = xpartveci(ietahalli)
       endif

       if (maxvxyzu >= 4) then
          eni = xpartveci(ieni)
       endif
       if (maxdvdx == maxp) then
          dvdxi(:) = xpartveci(idvxdxi:idvzdzi)
       else
          dvdxi(:) = 0.
       endif

       if (use_dustfrac) then
          dustfraci(:) = xpartveci(idustfraci:idustfraciend)
          dustfracisum = sum(dustfraci(:))
          tstopi(:)    = xpartveci(itstop:itstopend)
       endif
    else
       rho1i = 0.
       vwavei = 0.
    endif

#ifdef GRAVITY
    !--add self-contribution
    call kernel_softening(0.,0.,potensoft0,dum)
    epoti = 0.5*pmassi*(fsum(ipot) + pmassi*potensoft0*hi1)
    !
    !--add contribution from distant nodes, expand these in Taylor series about node centre
    !  use xcen directly, -1 is placeholder
    !
    call get_distance_from_centre_of_mass(cell%icell,xi,yi,zi,dx,dy,dz)
    call expand_fgrav_in_taylor_series(cell%fgrav,dx,dy,dz,fxi,fyi,fzi,poti)
    fsum(ifxi) = fsum(ifxi) + fxi
    fsum(ifyi) = fsum(ifyi) + fyi
    fsum(ifzi) = fsum(ifzi) + fzi
    if (gr .and. ien_type == ien_etotal) then
       fsum(idudtdissi) = fsum(idudtdissi) + vxi*fxi + vyi*fyi + vzi*fzi
    endif
    epoti = epoti + 0.5*pmassi*poti
    poten(i) = real(epoti,kind=kind(poten))
#endif

    if (mhd .and. iamgasi) then
       !
       !--for MHD, need to make the force stable when beta < 1.  In this regime,
       !  subtract off the B(div B)/rho term (Borve, Omang & Trulsen 2001, 2005);
       !  outside of this regime, do nothing, but (smoothly) transition between
       !  regimes.  Tests in Feb 2018 suggested that beginning the transition
       !  too close to beta = 1 lead to some inaccuracies, and that the optimal
       !  transition range was 2-10 or 2-5, depending on the test.
       !
       divBsymmi  = fsum(idivBsymi)
       if (B2i > 0.0) then
          betai = 2.0*pri/B2i
          if (betai < 2.0) then
             frac_divB = 1.0
          elseif (betai < 10.0) then
             frac_divB = (10.0 - betai)*0.125
          else
             frac_divB = 0.0
          endif
       else
          frac_divB = 0.0
       endif
       fsum(ifxi) = fsum(ifxi) - Bxyzi(1)*divBsymmi*frac_divB
       fsum(ifyi) = fsum(ifyi) - Bxyzi(2)*divBsymmi*frac_divB
       fsum(ifzi) = fsum(ifzi) - Bxyzi(3)*divBsymmi*frac_divB
       divBsymm(i) = real(rhoi*divBsymmi,kind=kind(divBsymm)) ! for output store div B as rho*div B
    endif

    f2i = fsum(ifxi)**2 + fsum(ifyi)**2 + fsum(ifzi)**2

#ifdef DRIVING
    ! force is first initialised in driving routine
    fxyzu(1,i) = fxyzu(1,i) + fsum(ifxi)
    fxyzu(2,i) = fxyzu(2,i) + fsum(ifyi)
    fxyzu(3,i) = fxyzu(3,i) + fsum(ifzi)
#else
    fxyzu(1,i) = fsum(ifxi)
    fxyzu(2,i) = fsum(ifyi)
    fxyzu(3,i) = fsum(ifzi)
#endif
    if (use_dust) then
       if (drag_implicit) then
          fxyz_drag(1,i) = fsum(ifdragxi)
          fxyz_drag(2,i) = fsum(ifdragyi)
          fxyz_drag(3,i) = fsum(ifdragzi)
       endif
       fxyzu(1,i) = fxyzu(1,i) + fsum(ifdragxi)
       fxyzu(2,i) = fxyzu(2,i) + fsum(ifdragyi)
       fxyzu(3,i) = fxyzu(3,i) + fsum(ifdragzi)
    endif

    drhodti = pmassi*fsum(idrhodti)

    isgas: if (iamgasi) then
       divvi = -drhodti*rho1i
       if (ndivcurlv >= 1) divcurlv(1,i) = real(divvi,kind=kind(divcurlv)) ! store divv from forces

       if (maxvxyzu >= 4 .or. lightcurve) then
          if (maxdvdx == maxp .and. realviscosity) then
             shearvisc = shearfunc(xi,yi,zi,spsoundi)
             straini   = strain_from_dvdx(dvdxi(:))
             fsum(idudtdissi) = fsum(idudtdissi) + (bulkvisc - 2./3.*shearvisc)*divvi**2 &
                           + 0.5*shearvisc*(straini(1)**2 + 2.*(straini(2)**2 + straini(3)**2 + straini(5)**2) &
                           + straini(4)**2 + straini(6)**2)
          endif
          fxyz4 = 0.
          if (ien_type == ien_etotal) then
             fxyz4 = fxyz4 + fsum(idudtdissi) + fsum(idendtdissi)
          elseif (ien_type == ien_entropy_s) then
             fxyz4 = fxyz4 + real(u0i/tempi*(fsum(idudtdissi) + fsum(idendtdissi))/kboltz)
          elseif (ien_type == ien_entropy) then ! here eni is the entropy
             if (gr .and. ishock_heating > 0) then
                fxyz4 = fxyz4 + (gammai - 1.)*densi**(1.-gammai)*u0i*fsum(idudtdissi)
             elseif (ishock_heating > 0) then
                fxyz4 = fxyz4 + (gammai - 1.)*rhoi**(1.-gammai)*fsum(idudtdissi)
             endif
             ! add conductivity for GR
             if (gr) then
                fxyz4 = fxyz4 + (gammai - 1.)*densi**(1.-gammai)*u0i*fsum(idendtdissi)
             endif
#ifdef GR
#ifdef ISENTROPIC
             fxyz4 = 0.
#endif
             if (lightcurve) then
                luminosity(i) = real(pmassi*u0i*(fsum(idendtdissi)+fsum(idudtdissi)),kind=kind(luminosity))
             endif
#endif
          elseif (ieos==16) then ! here eni is the temperature
             if (abs(damp) < tiny(damp)) then
                rho_cgs = rhoi * unit_density
                call eos_shen_get_dTdu(rho_cgs,eni,0.05,dTdui_cgs)
                dTdui = real(dTdui_cgs / unit_ergg)
                !use cgs
                fxyz4 = fxyz4 + dTdui*(pri*rho1i*rho1i*drhodti + fsum(idudtdissi))
             else
                fxyz4 = 0.
             endif
          else ! eni is the internal energy
             fac = rhoi/rhogasi
             pdv_work = pri*rho1i*rho1i*drhodti
             if (ipdv_heating > 0) then
                fxyz4 = fxyz4 + fac*pdv_work
             endif
             if (ishock_heating > 0) then
                if (fsum(idudtdissi) < -epsilon(0.)) &
                   call warning('force','-ve entropy derivative',i,var='dudt_diss',val=fsum(idudtdissi))
                fxyz4 = fxyz4 + fac*fsum(idudtdissi)
             endif
             !
             !--store pdV work and shock heating in separate array needed for some applications
             !  this is a kind of luminosity if it were all radiated
             !
             if (lightcurve) then
                pdv_work = pri*rho1i*rho1i*drhodti
                if (pdv_work > tiny(pdv_work)) then ! pdv_work < 0 is possible, and we want to ignore this case
                   dudt_radi = fac*pdv_work + fac*fsum(idudtdissi)
                else
                   dudt_radi = fac*fsum(idudtdissi)
                endif
                luminosity(i) = real(pmassi*dudt_radi,kind=kind(luminosity))
             endif
             if (mhd_nonideal) then
                call nicil_get_dudt_nimhd(dudtnonideal,etaohmi,etaambii,rhoi,curlBi,Bxyzi)
                fxyz4 = fxyz4 + fac*dudtnonideal
             endif
             !--add conductivity and resistive heating
             fxyz4 = fxyz4 + fac*fsum(idendtdissi)
             if (icooling > 0 .and. dt > 0. .and. .not. cooling_in_step) then
                if (h2chemistry) then
                   !
                   ! Call cooling routine, requiring total density, some distance measure and
                   ! abundances in the 'abund' format
                   !
                   call energ_cooling(xi,yi,zi,vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,&
                        dust_temp(i),eos_vars(imu,i), eos_vars(igamma,i))
                elseif (store_dust_temperature) then
                   ! cooling with stored dust temperature
                   if (do_nucleation) then
                      call energ_cooling(xi,yi,zi,vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,&
                           dust_temp(i),eos_vars(imu,i),eos_vars(igamma,i),nucleation(idK2,i),nucleation(idkappa,i))
                   elseif (update_muGamma) then
                      call energ_cooling(xi,yi,zi,vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,&
                           dust_temp(i),eos_vars(imu,i),eos_vars(igamma,i))
                   else
                      call energ_cooling(xi,yi,zi,vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,dust_temp(i))
                   endif
                else
                   ! cooling without stored dust temperature
                   call energ_cooling(xi,yi,zi,vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool)
                endif
                fxyz4 = fxyz4 + fac*dudtcool
             endif
             !  if (nuclear_burning) then
             !     call energ_nuclear(xi,yi,zi,vxyzu(4,i),dudtnuc,rhoi,0.,Tgas=tempi)
             !     fxyz4 = fxyz4 + fac*dudtnuc
             !  endif
             if (sinks_have_heating(nptmass,xyzmh_ptmass)) then
                call energ_sinkheat(nptmass,xyzmh_ptmass,xi,yi,zi,dudtheat)
                fxyz4 = fxyz4 + fac*dudtheat
             endif
             ! extra terms in du/dt from one fluid dust
             if (use_dustfrac) then
                !fxyz4 = fxyz4 + 0.5*fac*rho1i*fsum(idudtdusti)
                fxyz4 = fxyz4 + 0.5*fac*rho1i*sum(fsum(idudtdusti:idudtdustiend))
             endif
          endif
          if (do_radiation .and. implicit_radiation) then
             luminosity(i) = real(pmassi*fxyz4,kind=kind(luminosity))
             !fxyzu(4,i) = 0.
          else
             if (maxvxyzu >= 4) fxyzu(4,i) = fxyz4
             if (icooling == 9) then
                call energ_cooling(xi,yi,zi,vxyzu(4,i),rhoi,dt,divcurlv(1,i),dudtcool,duhydro=fxyz4,ipart=i)
                dusph(i) = fxyz4
             endif
          endif
       endif

       if (mhd) then
          !
          ! sum returns d(B/rho)/dt, just what we want!
          !
          dBevol(1,i) = fsum(idBevolxi)
          dBevol(2,i) = fsum(idBevolyi)
          dBevol(3,i) = fsum(idBevolzi)
          !
          ! hyperbolic/parabolic cleaning terms (dpsi/dt) from Tricco & Price (2012)
          !
          if (psidecayfac > 0.) then
             vcleani = overcleanfac*vwavei
             dtau = psidecayfac*vcleani*hi1
             !
             ! we clean using the difference operator for div B
             !
             psii = xpartveci(ipsi)

             ! new cleaning evolving d/dt (psi/c_h)
             dBevol(4,i) = -vcleani*fsum(idivBdiffi)*rho1i - psii*dtau - 0.5*psii*divvi
             dtclean   = C_cour*hi/(vcleani + tiny(0.))
          endif
       endif

       if (use_dustfrac) then
          !--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
          ddustevol(:,i) = 0.5*(fsum(iddustevoli:iddustevoliend)*rho1i/((1.-dustfraci(1:maxdustsmall))**2.))
          deltav(1,:,i)  = fsum(ideltavxi:ideltavxiend)
          deltav(2,:,i)  = fsum(ideltavyi:ideltavyiend)
          deltav(3,:,i)  = fsum(ideltavzi:ideltavziend)
       endif
       ! timestep based on Courant condition
       vsigdtc = max(vsigmax,vwavei)
       if (vsigdtc > tiny(vsigdtc)) then
          dtc = C_cour*hi/(vsigdtc*max(alpha,1.0))
       endif

       ! cooling timestep dt < fac*u/(du/dt)
       if (maxvxyzu >= 4 .and. .not. gr) then ! not with gr which uses entropy
          if (eni + dtc*fxyzu(4,i) < epsilon(0.) .and. eni > epsilon(0.)) dtcool = C_cool*abs(eni/fxyzu(4,i))
       endif

       ! s entropy timestep to avoid too large s entropy leads to infinite temperature
       if (gr .and. ien_type == ien_entropy_s) then
          dtent = C_ent*abs(pxyzu(4,i)/fxyzu(4,i))
       endif

       ! timestep based on non-ideal MHD
       if (mhd_nonideal) then
          call nicil_get_dt_nimhd(dtohmi,dthalli,dtambii,hi,etaohmi,etahalli,etaambii)
          if ( use_STS ) then
             dtdiffi = min(dtohmi,dtambii)
             dtdiff  = min(dtdiff,dtdiffi)
             dtohmi  = bignumber
             dtambii = bignumber
          endif
       endif

       ! timestep from physical viscosity
       dtvisci = dt_viscosity(xi,yi,zi,hi,spsoundi)

       ! Check to ensure we have enough resolution for gas-dust pairs, where
       ! ts_min is already minimised over all dust neighbours for gas particle i
       if (ts_min < bignumber) then
          if (hi > ts_min*spsoundi) then
             ndustres       = ndustres + 1
             dustresfacmean = dustresfacmean + hi/(ts_min*spsoundi)
             dustresfacmax  = max(dustresfacmax, hi/(ts_min*spsoundi))
          endif
       endif

    else ! not gas
       if (use_dustgrowth .and. iamdusti) then
          !- return interpolations to their respective arrays
          dustgasprop(2,i) = fsum(idensgasi) !- rhogas
          !- interpolations are mass weigthed, divide result by rhog,i
          dustgasprop(4,i) = sqrt(fsum(idvix)**2 + fsum(idviy)**2 + fsum(idviz)**2)/dustgasprop(2,i) !- |dv|
          dustgasprop(1,i) = fsum(icsi)/dustgasprop(2,i) !- sound speed

          !- get the Stokes number with get_ts using the interpolated quantities
          rhoi             = xpartveci(irhoi)
          gdensi           = xpartveci(igraindensi)
          gmassi           = xpartveci(igrainmassi)
          if (use_porosity) then
             call get_ts(idrag,1,get_size(gmassi,gdensi,filfac(i)),gdensi*filfac(i),&
             dustgasprop(2,i),rhoi,dustgasprop(1,i),dustgasprop(4,i)**2,tstopint,ireg)
             dragreg(i) = ireg
          else
             call get_ts(idrag,1,get_size(gmassi,gdensi),gdensi,&
             dustgasprop(2,i),rhoi,dustgasprop(1,i),dustgasprop(4,i)**2,tstopint,ireg)
          endif
          dustgasprop(3,i) = tstopint * Omega_k(i) !- Stokes number
       endif

       if (maxvxyzu > 4) fxyzu(4,i) = 0.
       ! timestep based on Courant condition for non-gas particles
       vsigdtc = vsigmax
       if (vsigdtc > tiny(vsigdtc)) then
          dtc = C_cour*hi/vsigdtc
       endif
    endif isgas

    ! initialise timestep to Courant timestep & perform sanity check
    dti = dtc
    if (dtc < tiny(dtc) .or. dtc > huge(dtc)) call fatal('force','invalid dtc',var='dtc',val=dtc)

    ! timestep based on force condition
    if (abs(f2i) > epsilon(f2i)) then
       dtf = C_force*sqrt(hi/sqrt(f2i))
    endif

    ! one fluid dust timestep
    if (use_dustfrac .and. iamgasi) then
       if (minval(dustfraci(1:ndustsmall)) > 0. .and. spsoundi > 0. .and. dustfracisum > epsilon(0.)) then
          tseff = (1.-dustfracisum)/dustfracisum*sum(dustfraci(1:ndustsmall)*tstopi(1:ndustsmall))
          dtdustdenom = dustfracisum*tseff*spsoundi**2
          if (dtdustdenom > tiny(dtdustdenom)) then
             dtdusti = C_force*hi*hi/dtdustdenom
          endif
       endif
    endif

    ! stopping time and timestep based on it (when using dust-as-particles)
    if (use_dust .and. use_dustfrac) then
       tstop(:,i) = tstopi(:)
    elseif (use_dust .and. .not.use_dustfrac) then
       tstop(:,i) = ts_min
       if (drag_implicit) then
          dtdrag = 90.*ts_min
       else
          dtdrag = 0.9*ts_min
       endif
    endif

    if (do_radiation .and. iamgasi .and. .not.implicit_radiation) then
       if (radprop(ithick,i) < 0.5) then
          drad(iradxi,i) = 0.
       else
          if (iopacity_type == 0) then ! infinite opacity equals no radiation diffusion
             drad(iradxi,i) = radprop(iradP,i)*drhodti*rho1i*rho1i
             dtradi = bignumber
          else
             drad(iradxi,i) = fsum(idradi) + radprop(iradP,i)*drhodti*rho1i*rho1i
             c_code     = get_c_code()
             radkappai  = xpartveci(iradkappai)
             radlambdai = xpartveci(iradlambdai)
             ! eq30 Whitehouse & Bate 2004
             dtradi = C_rad*hi*hi*rhoi*radkappai/c_code/radlambdai
             ! additional timestep constraint to ensure that
             ! radiation energy is positive after the integration
             if ((rad(iradxi,i) + dtradi*drad(iradxi,i)) < 0) then
                if (rad(iradxi,i) > 0.) dtradi = -rad(iradxi,i)/drad(iradxi,i)/1e1
                call warning('force','radiation may become negative, limiting timestep')
             endif
          endif
       endif
    endif

#ifdef IND_TIMESTEPS
    !-- The new timestep for particle i
    dtitmp = min(dtf,dtcool,dtclean,dtvisci,dtdrag,dtohmi,dthalli,dtambii,dtdusti,dtradi)

    if (dtitmp < dti + tiny(dtitmp) .and. dtitmp < dtmax) then
       dti     = dtitmp
       if (dti < tiny(dti) .or. dti > huge(dti)) call fatal('force','invalid dti',var='dti',val=dti) ! sanity check
       dtcheck = .true.
       dtrat   = dtc/dti
       if ( iamgasi ) then
          call check_dtmin(dtcheck,dti,dtf    ,dtrat,ndtforce  ,dtfrcfacmean  ,dtfrcfacmax  ,dtchar,'dt_gasforce' )
          call check_dtmin(dtcheck,dti,dtcool ,dtrat,ndtcool   ,dtcoolfacmean ,dtcoolfacmax ,dtchar,'dt_cool'     )
          call check_dtmin(dtcheck,dti,dtvisci,dtrat,ndtvisc   ,dtviscfacmean ,dtviscfacmax ,dtchar,'dt_visc'     )
          call check_dtmin(dtcheck,dti,dtdrag ,dtrat,ndtdrag   ,dtdragfacmean ,dtdragfacmax ,dtchar,'dt_gasdrag'  )
          call check_dtmin(dtcheck,dti,dtohmi ,dtrat,ndtohm    ,dtohmfacmean  ,dtohmfacmax  ,dtchar,'dt_ohm'      )
          call check_dtmin(dtcheck,dti,dthalli,dtrat,ndthall   ,dthallfacmean ,dthallfacmax ,dtchar,'dt_hall'     )
          call check_dtmin(dtcheck,dti,dtambii,dtrat,ndtambi   ,dtambifacmean ,dtambifacmax ,dtchar,'dt_ambi'     )
          call check_dtmin(dtcheck,dti,dtdusti,dtrat,ndtdust   ,dtdustfacmean ,dtdustfacmax ,dtchar,'dt_dust'     )
          call check_dtmin(dtcheck,dti,dtradi ,dtrat,ndtrad    ,dtradfacmean  ,dtradfacmax  ,dtchar,'dt_radiation')
          call check_dtmin(dtcheck,dti,dtclean,dtrat,ndtclean  ,dtcleanfacmean,dtcleanfacmax,dtchar,'dt_clean'    )
       else
          call check_dtmin(dtcheck,dti,dtf    ,dtrat,ndtforceng,dtfrcngfacmean,dtfrcngfacmax,dtchar,'dt_force'    )
          call check_dtmin(dtcheck,dti,dtdrag ,dtrat,ndtdragd  ,dtdragdfacmean,dtdragdfacmax,dtchar,'dt_drag'     )
       endif
       if (dtcheck) call fatal('force','unknown dti',var='dti',val=dti)
    else
       dtchar = 'dt_courant'
    endif
    !
    allow_decrease = ((icall < 2) .and. sts_it_n)
    call get_newbin(dti,dtmax,ibin(i),allow_decrease,dtchar=dtchar) ! get new timestep bin based on dti
    !
    ! Saitoh-Makino limiter, do not allow timestep to be more than 1 bin away from neighbours
    !
    ibin(i) = max(ibin(i),ibin_neighi-1_1)
    !
    ! find the new maximum value of ibin
    nbinmaxnew = max(nbinmaxnew,int(ibin(i)))
    ncheckbin  = ncheckbin + 1

    ! ibin_sts: based entirely upon the diffusive timescale
    if ( use_sts ) then
       ibin_sts(i) = 0 ! we actually want dtdiff, and this is just a tracer; should reduce the number of sts active particles for speed
       call get_newbin(dtdiffi,dtmax,ibin_sts(i),allow_decrease,.false.)
       nbinmaxstsnew = max(nbinmaxstsnew,int(ibin_sts(i)))
    endif

#else
    ! global timestep needs to be minimum over all particles

    dtcourant = min(dtcourant,dtc)
    dtforce   = min(dtforce,dtf,dtcool,dtdrag,dtdusti,dtclean,dtent)
    dtvisc    = min(dtvisc,dtvisci)
    if (mhd_nonideal .and. iamgasi) then
       dtohm  = min(dtohm,  dtohmi  )
       dthall = min(dthall, dthalli )
       dtambi = min(dtambi, dtambii )
    endif
    dtmini  = min(dtmini,dti)
    dtmaxi  = max(dtmaxi,dti)
    dtrad   = min(dtrad,dtradi)
#endif
 enddo over_parts
end subroutine finish_cell_and_store_results

!-----------------------------------------------------------------------------
!+
!  Apply reconstruction to velocity gradients
!  As described in Price & Laibe (2020), MNRAS 495, 3929-3934
!+
!-----------------------------------------------------------------------------
subroutine reconstruct_dv(projv,dx,dy,dz,rx,ry,rz,dvdxi,dvdxj,projvstar,ilimiter)
 real, intent(in)  :: projv,dx,dy,dz,rx,ry,rz,dvdxi(9),dvdxj(9)
 real, intent(out) :: projvstar
 integer, intent(in) :: ilimiter
 real :: slopei,slopej,slope,sep

 sep = 0.5
 ! CAUTION: here we use dx, not the unit vector to
 ! define the projected slope. This is fine as
 ! long as the slope limiter is linear, otherwise
 ! one should use the unit vector
 slopei = dx*(rx*dvdxi(1) + ry*dvdxi(4) + rz*dvdxi(7)) &
        + dy*(rx*dvdxi(2) + ry*dvdxi(5) + rz*dvdxi(8)) &
        + dz*(rx*dvdxi(3) + ry*dvdxi(6) + rz*dvdxi(9))

 slopej = dx*(rx*dvdxj(1) + ry*dvdxj(4) + rz*dvdxj(7)) &
        + dy*(rx*dvdxj(2) + ry*dvdxj(5) + rz*dvdxj(8)) &
        + dz*(rx*dvdxj(3) + ry*dvdxj(6) + rz*dvdxj(9))

 if (ilimiter > 0) then
    slope = slope_limiter(slopei,slopej)
    projvstar = projv - 2.*sep*slope
 else
    !
    !--reconstruction with no slope limiter
    !  (mainly useful for testing purposes)
    !
    projvstar = projv - sep*(slopei + slopej)
 endif
 ! apply entropy condition
 !if (projvstar*projv < 0.) projvstar = sign(1.0,projv)*min(abs(projv),abs(projvstar))
 !projvstar = sign(1.0,projv)*min(abs(projv),abs(projvstar))

end subroutine reconstruct_dv

!-----------------------------------------------------------------------------
!+
!  Apply reconstruction to GR velocity gradients
!+
!-----------------------------------------------------------------------------
subroutine reconstruct_dv_gr(projvi,projvj,rx,ry,rz,dr,dvdxi,dvdxj,projvstari,projvstarj,ilimiter)
 real, intent(in)  :: projvi,projvj,rx,ry,rz,dr,dvdxi(9),dvdxj(9)
 real, intent(out) :: projvstari,projvstarj
 integer, intent(in) :: ilimiter
 real :: slopei,slopej,slope

 ! CAUTION: here we use dx, not the unit vector to
 ! define the projected slope. This is fine as
 ! long as the slope limiter is linear, otherwise
 ! one should use the unit
 slopei = rx*(rx*dvdxi(1) + ry*dvdxi(4) + rz*dvdxi(7)) &
        + ry*(rx*dvdxi(2) + ry*dvdxi(5) + rz*dvdxi(8)) &
        + rz*(rx*dvdxi(3) + ry*dvdxi(6) + rz*dvdxi(9))

 slopej = rx*(rx*dvdxj(1) + ry*dvdxj(4) + rz*dvdxj(7)) &
        + ry*(rx*dvdxj(2) + ry*dvdxj(5) + rz*dvdxj(8)) &
        + rz*(rx*dvdxj(3) + ry*dvdxj(6) + rz*dvdxj(9))

 if (ilimiter > 0) then
    slope = dr*slope_limiter_gr(slopei,slopej)
 else
    !
    !--reconstruction with no slope limiter
    !  (mainly useful for testing purposes)
    !
    slope = 0.5*dr*(slopei + slopej)
 endif
 projvstari = (projvi - slope) / (1. - projvi*slope)
 projvstarj = (projvj + slope) / (1. + projvj*slope)

end subroutine reconstruct_dv_gr

!-----------------------------------------------------------------------------
!+
!  Slope limiter used for velocity gradient reconstruction,
!  as described in Price & Laibe (2020), MNRAS 495, 3929-3934
!+
!-----------------------------------------------------------------------------
real function slope_limiter(sl,sr) result(s)
 real, intent(in) :: sl,sr
! integer, intent(in) :: ilimiter

 s = 0.

 ! Van Leer monotonised central (MC)
 if (sl*sr > 0.) s = sign(1.0,sl)*min(abs(0.5*(sl + sr)),2.*abs(sl),2.*abs(sr))

 ! Van Leer
 !if (sl*sr > 0.) s = 2.*sl*sr/(sl + sr)

 ! minmod
 !if (sl > 0. .and. sr > 0.) then
 !   s = min(abs(sl),abs(sr))
 !elseif (sl < 0. .and. sr < 0.) then
 !   s = -min(abs(sl),abs(sr))
 !endif

end function slope_limiter

!-----------------------------------------------------------------------------
!+
!  Slope limiter used for velocity gradient reconstruction,
!  as described in Price & Laibe (2020), MNRAS 495, 3929-3934
!+
!-----------------------------------------------------------------------------
real function slope_limiter_gr(sl,sr) result(s)
 real, intent(in) :: sl,sr

 s = 0.
 ! Van Leer is the only slope limiter we found that works for relativistic shocks
 if (sl*sr > 0.) s = 2.*sl*sr/(sl + sr)

end function slope_limiter_gr

!-----------------------------------------------------------------------------
!+
!  Compute additional drag terms needed for the implicit scheme
!  As described in Loren-Anguilar & Bate (2015), MNRAS 454, 4114-4119
!+
!-----------------------------------------------------------------------------
subroutine get_drag_terms(tsijtmp,dt,sdrag1,sdrag2,fxi_drag,fyi_drag,fzi_drag,&
                          runix,runiy,runiz,fxyz_drag,projf_drag)
 real,      intent(in)    :: tsijtmp,dt,fxi_drag,fyi_drag,fzi_drag,runix,runiy,runiz
 real,      intent(in)    :: fxyz_drag(:)
 real,      intent(out)   :: sdrag1,sdrag2,projf_drag

 projf_drag = (fxi_drag - fxyz_drag(1))*runix + (fyi_drag - fxyz_drag(2))*runiy + (fzi_drag - fxyz_drag(3))*runiz

 if (dt > epsilon(0.)) then
    sdrag1 = (1. - exp(-dt/tsijtmp))/dt
    sdrag2 = ((dt+tsijtmp)*(1. - exp(-dt/tsijtmp)) - dt)/dt
 else
    sdrag1 = 1./tsijtmp
    sdrag2 = 0.
 endif
end subroutine get_drag_terms

end module forces
