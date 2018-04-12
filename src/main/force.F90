!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: forces
!
!  DESCRIPTION:
!   This module is the "guts" of the code
!   Calculates force and rates of change for all particles
!
!  REFERENCES:
!    Price (2012), J. Comp. Phys.
!    Lodato & Price (2010), MNRAS
!    Price & Federrath (2010), MNRAS
!    Tricco & Price (2012), J. Comp. Phys.
!
!  OWNER: Conrad Chan
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, chem, cooling, dim, dust, eos, fastmath, growth,
!    io, io_summary, kdtree, kernel, linklist, mpiderivs, mpiforce,
!    mpiutils, nicil, options, part, physcon, ptmass, stack, timestep,
!    timestep_ind, timestep_sts, units, viscosity
!+
!--------------------------------------------------------------------------
module forces
 use dim,      only:maxfsum,maxxpartveciforce,maxBevol,maxp,ndivcurlB,ndivcurlv,ndusttypes
 use mpiforce, only:cellforce,stackforce

 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 integer, parameter :: maxcellcache = 50000

 public :: force

 !--indexing for xpartveci array
 integer, parameter :: &
       ixi  = 1, &
       iyi  = 2, &
       izi  = 3, &
       ihi  = 4, &
       ivxi = 5, &
       ivyi = 6, &
       ivzi = 7, &
       ieni = 8, &
       iBevolxi = 9, &
       iBevolyi = 10, &
       iBevolzi = 11, &
       ipsi = 12, &
       igradhi1    = 13, &
       igradhi2    = 14, &
       ialphai     = 15, &
       ialphaBi    = 16, &
       ivwavei     = 17, &
       irhoi       = 18, &
       irhogasi    = 19, &
       ispsoundi   = 20, &
       isxxi       = 21, &
       isxyi       = 22, &
       isxzi       = 23, &
       isyyi       = 24, &
       isyzi       = 25, &
       iszzi       = 26, &
       ivisctermisoi   = 27, &
       ivisctermanisoi = 28, &
       ipri        = 29, &
       ipro2i      = 30, &
       ietaohmi    = 31, &
       ietahalli   = 32, &
       ietaambii   = 33, &
       ijcbcbxi    = 34, &
       ijcbcbyi    = 35, &
       ijcbcbzi    = 36, &
       ijcbxi      = 37, &
       ijcbyi      = 38, &
       ijcbzi      = 39, &
       iponrhoi    = 40, &
       icurlBxi    = 41, &
       icurlByi    = 42, &
       icurlBzi    = 43, &
       !--dust arrays initial index
       igrainsizei = 44, &
       igraindensi = 45, &
       idustfraci  = 46, &
       itstop      = 47 + (ndusttypes-1), &
       !--dust arrays final index
       idustfraciend = itstop-1, &
       itstopend     = maxxpartveciforce

 !--indexing for fsum array
 integer, parameter :: &
       ifxi        = 1, &
       ifyi        = 2, &
       ifzi        = 3, &
       ipot        = 4, &
       idrhodti    = 5, &
       idudtdissi  = 6, &
       idendtdissi = 7, &
       idivBsymi   = 8, &
       idBevolxi   = 9, &
       idBevolyi   = 10, &
       idBevolzi   = 11, &
       idivBdiffi  = 12, &
       !--dust arrays initial index
       iddustfraci = 13, &
       idudtdusti  = 14 +   (ndusttypes-1), &
       ideltavxi   = 15 + 2*(ndusttypes-1), &
       ideltavyi   = 16 + 3*(ndusttypes-1), &
       ideltavzi   = 17 + 4*(ndusttypes-1), &
       !--dust arrays final index
       iddustfraciend = idudtdusti-1, &
       idudtdustiend  = ideltavxi -1, &
       ideltavxiend   = ideltavyi -1, &
       ideltavyiend   = ideltavzi -1, &
       ideltavziend   = maxfsum

 private

contains

!----------------------------------------------------------------
!+
!  compute all forces and rates of change on the particles
!+
!----------------------------------------------------------------
subroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,&
                 ipart_rhomax,dt,stressmax,temperature)
 use dim,          only:maxvxyzu,maxalpha,maxneigh,maxstrain,&
                        switches_done_in_derivs,mhd,mhd_nonideal,lightcurve
 use io,           only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs
 use linklist,     only:ncells,ifirstincell,get_neighbour_list,get_hmaxcell,get_cell_location
 use options,      only:iresistive_heating
 use part,         only:rhoh,dhdrho,rhoanddhdrho,alphaind,nabundances,ll,get_partinfo,iactive,gradh,&
                        hrho,iphase,maxphase,igas,iboundary,maxgradh,straintensor, &
                        eta_nimhd,deltav,poten
 use timestep,     only:dtcourant,dtforce,bignumber,dtdiff
 use io_summary,   only:summary_variable, &
                        iosumdtf,iosumdtd,iosumdtv,iosumdtc,iosumdto,iosumdth,iosumdta, &
                        iosumdgs,iosumdge,iosumdgr,iosumdtfng,iosumdtdd,iosumdte
#ifdef FINVSQRT
 use fastmath,     only:finvsqrt
#endif
 use physcon,      only:pi
 use viscosity,    only:irealvisc,shearfunc,dt_viscosity
#ifdef IND_TIMESTEPS
 use timestep_ind, only:nbinmax,ibinnow,get_newbin
 use timestep_sts, only:nbinmaxsts,ibin_sts
 use part,         only:ibin
 use timestep,     only:nsteps,time
#else
 use timestep,     only:C_cour,C_force
#endif
 use part,         only:divBsymm,isdead_or_accreted,h2chemistry,ngradh,gravity,ibin_wake
 use mpiutils,     only:reduce_mpi,reduceall_mpi,reduceloc_mpi,bcast_mpi
 use cooling,      only:energ_cooling
 use chem,         only:energ_h2cooling
#ifdef GRAVITY
 use kernel,       only:kernel_softening
 use kdtree,       only:expand_fgrav_in_taylor_series
 use linklist,     only:get_distance_from_centre_of_mass
 use part,         only:xyzmh_ptmass,nptmass,massoftype
 use ptmass,       only:icreate_sinks,rho_crit,r_crit2,&
                        rhomax_xyzh,rhomax_vxyz,rhomax_iphase,rhomax_divv,rhomax_ipart,rhomax_ibin
 use units,        only:unit_density
#endif
#ifdef DUST
 !use dust,         only:get_ts
 use kernel,       only:wkern_drag,cnormk_drag
#endif
 use nicil,        only:nimhd_get_jcbcb,nimhd_get_dt,nimhd_get_dBdt,nimhd_get_dudt
#ifdef LIGHTCURVE
 use part,         only:luminosity
#endif
#ifdef MPI
 use mpiderivs,    only:send_cell,recv_cells,check_send_finished,init_cell_exchange,finish_cell_exchange, &
                       recv_while_wait,reset_cell_counters
 use stack,        only:reserve_stack
 use stack,        only:stack_remote => force_stack_1
 use stack,        only:stack_waiting => force_stack_2
#endif

 integer,      intent(in)    :: icall,npart
 real,         intent(in)    :: xyzh(:,:)
 real,         intent(inout) :: vxyzu(:,:)
 real,         intent(in)    :: dustfrac(:,:),dustprop(:,:)
 real,         intent(inout) :: temperature(:)
 real,         intent(out)   :: fxyzu(:,:),ddustfrac(:,:),ddustprop(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real,         intent(out)   :: dBevol(:,:)
 real(kind=4), intent(inout) :: divcurlv(:,:)
 real(kind=4), intent(in)    :: divcurlB(:,:)
 real,         intent(in)    :: dt,stressmax
 integer,      intent(out)   :: ipart_rhomax ! test this particle for point mass creation

 real, save :: xyzcache(maxcellcache,4)
 integer, save :: listneigh(maxneigh)
!$omp threadprivate(xyzcache,listneigh)
 integer :: i,icell,nneigh
 integer :: nstokes,nsuper,ndrag,ndustres
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
 real    :: frac_stokes, frac_super
#endif
 logical :: realviscosity,useresistiveheat
#ifndef IND_TIMESTEPS
 real    :: dtmaxi
#else
 integer :: nbinmaxnew,nbinmaxstsnew,ncheckbin
 integer :: ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd
 integer :: ndtvisc,ndtohm,ndthall,ndtambi,ndtdust
 real    :: dtitmp,dtrat,dtmaxi
 real    :: dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean
 real    :: dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax
 real    :: dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean,dtdustfacmean
 real    :: dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax, dtdustfacmax
 logical :: allow_decrease,dtcheck
#endif
 integer(kind=1)           :: ibinnow_m1

 logical                   :: remote_export(nprocs)
 type(cellforce)           :: cell

#ifdef MPI
 logical                   :: do_export

 integer                   :: irequestsend(nprocs),irequestrecv(nprocs)
 type(cellforce)           :: xrecvbuf(nprocs),xsendbuf
#endif

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
 ibinnow_m1      = ibinnow - 1_1
#else
 ibinnow_m1      = 0
#endif
 dtmaxi          = 0.

 dustresfacmean  = 0.0
 dustresfacmax   = 0.0
 dtcourant       = bignumber
 dtforce         = bignumber
 dtvisc          = bignumber
 dtmini          = bignumber
 dtohm           = bignumber
 dthall          = bignumber
 dtambi          = bignumber
 if (iverbose >= 3 .and. id==master) write(iprint,*) 'forces: cell cache =',maxcellcache

 realviscosity    = (irealvisc > 0)
 useresistiveheat = (iresistive_heating > 0)
 if (ndivcurlv < 1) call fatal('force','divv not stored but it needs to be')
 if (switches_done_in_derivs) call fatal('force','need switches_done_in_derivs=.false.')

 !--dust/gas stuff
 ndrag         = 0
 nstokes       = 0
 nsuper        = 0
 ndustres      = 0

 ! sink particle creation
 ipart_rhomax  = 0
#ifdef GRAVITY
 rhomax        = 0.
#endif

#ifdef MPI
 call init_cell_exchange(xrecvbuf,irequestrecv)
 stack_waiting%n = 0
 stack_remote%n = 0
 call reset_cell_counters
#endif

!
!-- verification for non-ideal MHD
 if (mhd_nonideal .and. ndivcurlB < 4) call fatal('force','non-ideal MHD needs curl B stored, but ndivcurlB < 4')
!
!--check that compiled options are compatible with this routine
!
 if (maxgradh /= maxp) call fatal('force','need storage of gradh (maxgradh=maxp)')

!$omp parallel default(none) &
!$omp shared(ncells,ll,ifirstincell) &
!$omp shared(xyzh) &
!$omp shared(dustprop) &
!$omp shared(ddustprop) &
!$omp shared(vxyzu) &
!$omp shared(fxyzu) &
!$omp shared(divcurlv) &
!$omp shared(iphase) &
!$omp shared(straintensor) &
!$omp shared(gradh) &
!$omp shared(divcurlb) &
!$omp shared(bevol) &
!$omp shared(eta_nimhd) &
!$omp shared(alphaind) &
!$omp shared(stressmax) &
!$omp shared(divBsymm) &
!$omp shared(dBevol) &
!$omp shared(temperature) &
!$omp shared(dt) &
!$omp shared(nprocs,icall) &
!$omp shared(poten) &
!$omp private(icell,i) &
!$omp private(cell) &
!$omp private(remote_export) &
!$omp private(nneigh) &
#ifdef GRAVITY
!$omp shared(massoftype,npart) &
!$omp private(hi,pmassi,rhoi) &
!$omp private(iactivei,iamdusti,iamtypei) &
!$omp private(dx,dy,dz,poti,fxi,fyi,fzi,potensoft0,dum,epoti) &
!$omp shared(xyzmh_ptmass,nptmass) &
!$omp shared(rhomax,ipart_rhomax,icreate_sinks,rho_crit,r_crit2) &
!$omp private(rhomax_thread,ipart_rhomax_thread,use_part,j) &
#endif
#ifdef MPI
!$omp shared(id) &
!$omp private(do_export) &
!$omp shared(irequestrecv,irequestsend) &
!$omp shared(stack_remote,stack_waiting) &
!$omp shared(xsendbuf,xrecvbuf) &
#endif
#ifdef IND_TIMESTEPS
!$omp shared(ibin,ibin_sts,nbinmax,nbinmaxsts) &
!$omp private(allow_decrease,dtitmp,dtcheck,dtrat) &
!$omp reduction(+:ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd,ncheckbin,ndtvisc) &
!$omp reduction(+:ndtohm,ndthall,ndtambi,ndtdust,dtohmfacmean,dthallfacmean,dtambifacmean,dtdustfacmean) &
!$omp reduction(+:dtfrcfacmean,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean,dtviscfacmean) &
!$omp reduction(max:dtohmfacmax,dthallfacmax,dtambifacmax,dtdustfacmax) &
!$omp reduction(max:dtfrcfacmax,dtfrcngfacmax,dtdragfacmax,dtdragdfacmax,dtcoolfacmax,dtviscfacmax) &
!$omp reduction(max:nbinmaxnew,nbinmaxstsnew) &
#endif
!$omp reduction(min:dtohm,dthall,dtambi,dtdiff) &
!$omp reduction(min:dtcourant,dtforce,dtvisc) &
!$omp reduction(max:dtmaxi) &
!$omp reduction(min:dtmini) &
!$omp reduction(+:ndustres,dustresfacmean) &
!$omp reduction(max:dustresfacmax) &
!$omp shared(dustfrac) &
!$omp shared(ddustfrac) &
!$omp shared(deltav) &
!$omp shared(ibin_wake,ibinnow_m1)

!$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)

    !--skip empty cells AND inactive cells
    if (i <= 0) cycle over_cells

    cell%icell = icell

    call start_cell(cell,iphase,xyzh,vxyzu,gradh,divcurlv,divcurlB,straintensor,Bevol, &
                         dustfrac,dustprop,eta_nimhd,temperature,alphaind,stressmax)
    if (cell%npcell == 0) cycle over_cells

    call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
    !
    !--get the neighbour list and fill the cell cache
    !

    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true., &
#ifdef GRAVITY
                           f=cell%fgrav, &
#endif
                           remote_export=remote_export)
#ifdef MPI
    cell%owner                   = id
    cell%remote_export(1:nprocs) = remote_export
    do_export = any(remote_export)

!$omp critical
    call recv_cells(stack_remote,xrecvbuf,irequestrecv)
!$omp end critical

    if (do_export) then
!$omp critical
       if (stack_waiting%n > 0) call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
       call reserve_stack(stack_waiting,cell%waiting_index)
       ! export the cell: direction 0 for exporting
       call send_cell(cell,0,irequestsend,xsendbuf)
!$omp end critical
    endif
#endif

    call compute_cell(cell,listneigh,nneigh,Bevol,xyzh,vxyzu,fxyzu, &
                      iphase,divcurlv,divcurlB,alphaind,eta_nimhd, temperature, &
                      dustfrac,gradh,ibinnow_m1,ibin_wake,stressmax,xyzcache)

#ifdef MPI
    if (do_export) then
       stack_waiting%cells(cell%waiting_index) = cell
    else
#endif
       call finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,straintensor,&
                             divBsymm,divcurlv,dBevol,ddustfrac,deltav, &
                             dtcourant,dtforce,dtvisc,dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi, &
#ifdef IND_TIMESTEPS
                             nbinmaxnew,nbinmaxstsnew,ncheckbin, &
                             ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd, &
                             ndtvisc,ndtohm,ndthall,ndtambi,ndtdust, &
                             dtitmp,dtrat, &
                             dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean, &
                             dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax, &
                             dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean ,dtdustfacmean, &
                             dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax  ,dtdustfacmax, &
#endif
                             ndustres,dustresfacmax,dustresfacmean)

#ifdef MPI
    endif
#endif
 enddo over_cells
!$omp enddo

#ifdef MPI
!$omp barrier

!$omp single
 if (stack_waiting%n > 0) call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
 call recv_while_wait(stack_remote,xrecvbuf,irequestrecv,irequestsend)
 call reset_cell_counters
!$omp end single

 igot_remote: if (stack_remote%n > 0) then
!$omp do schedule(runtime)
    over_remote: do i = 1,stack_remote%n
       cell = stack_remote%cells(i)

       call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true., &
#ifdef GRAVITY
                         f=cell%fgrav, local_gravity=.true., &
#endif
                         cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti)

       call compute_cell(cell,listneigh,nneigh,Bevol,xyzh,vxyzu,fxyzu, &
                         iphase,divcurlv,divcurlB,alphaind,eta_nimhd, temperature, &
                         dustfrac,gradh,ibinnow_m1,ibin_wake,stressmax,xyzcache)

       cell%remote_export(id+1) = .false.

!$omp critical
       call recv_cells(stack_waiting,xrecvbuf,irequestrecv)
       call check_send_finished(stack_waiting,irequestsend,irequestrecv,xrecvbuf)
       call send_cell(cell,1,irequestsend,xsendbuf)
!$omp end critical
    enddo over_remote
!$omp enddo
!$omp barrier
!$omp single
    stack_remote%n = 0
    call check_send_finished(stack_waiting,irequestsend,irequestrecv,xrecvbuf)
!$omp end single
 endif igot_remote
!$omp barrier
!$omp single
 call recv_while_wait(stack_waiting,xrecvbuf,irequestrecv,irequestsend)
!$omp end single

 iam_waiting: if (stack_waiting%n > 0) then
!$omp do schedule(runtime)
    over_waiting: do i = 1, stack_waiting%n
       cell = stack_waiting%cells(i)

       if (any(cell%remote_export(1:nprocs))) then
          print*,id,cell%remote_export(1:nprocs)
          call fatal('force', 'not all results returned from remote processor')
       endif

       call finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,straintensor, &
                                          divBsymm,divcurlv,dBevol,ddustfrac,deltav, &
                                          dtcourant,dtforce,dtvisc,dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi, &
#ifdef IND_TIMESTEPS
                                          nbinmaxnew,nbinmaxstsnew,ncheckbin, &
                                          ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd, &
                                          ndtvisc,ndtohm,ndthall,ndtambi,ndtdust, &
                                          dtitmp,dtrat, &
                                          dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean, &
                                          dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax, &
                                          dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean ,dtdustfacmean, &
                                          dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax  ,dtdustfacmax, &
#endif
                                          ndustres,dustresfacmax,dustresfacmean)

    enddo over_waiting
!$omp enddo
!$omp barrier
!$omp single
    stack_waiting%n = 0
!$omp end single
 endif iam_waiting

!$omp single
 call finish_cell_exchange(irequestrecv,xsendbuf)
!$omp end single
#endif

#ifdef GRAVITY
 if (icreate_sinks > 0) then
    rhomax_thread = 0.
    ipart_rhomax_thread = 0
!$omp do schedule(runtime)
    do i=1,npart
       hi = xyzh(4,i)
#ifdef IND_TIMESTEPS
       if (iactive(iphase(i)) .and..not.isdead_or_accreted(hi)) then
#else
       if (.not.isdead_or_accreted(hi)) then
#endif
          if (maxphase==maxp) then
             call get_partinfo(iphase(i),iactivei,iamdusti,iamtypei)
          else
             iamtypei = igas
          endif
          pmassi = massoftype(iamtypei)
          rhoi = rhoh(hi,pmassi)
          if (rhoi > rho_crit) then
             if (rhoi > rhomax_thread) then
                !
                !--find the maximum density on particles outside the
                !  allowed minimum distance from other sink particles
                !
                use_part = .true.
                over_ptmass: do j=1,nptmass
                   if ((xyzh(1,i) - xyzmh_ptmass(1,j))**2 &
                     + (xyzh(2,i) - xyzmh_ptmass(2,j))**2 &
                     + (xyzh(3,i) - xyzmh_ptmass(3,j))**2 < r_crit2) then
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
    enddo
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

#ifdef IND_TIMESTEPS
 ! check for nbinmaxnew = 0, can happen if all particles
 ! are dead/inactive, e.g. after sink creation
 if (ncheckbin==0) then
    nbinmaxnew    = nbinmax
    nbinmaxstsnew = nbinmaxsts
 endif
#endif
!$omp end parallel

#ifdef GRAVITY
 if (reduceall_mpi('max',ipart_rhomax) > 0) then
    call reduceloc_mpi('max',rhomax,id_rhomax)
    if (id == id_rhomax) then
       rhomax_ipart  = ipart_rhomax
       rhomax_xyzh   = xyzh(1:4,ipart_rhomax)
       rhomax_vxyz   = vxyzu(1:3,ipart_rhomax)
       rhomax_iphase = iphase(ipart_rhomax)
       rhomax_divv   = divcurlv(1,ipart_rhomax)
#ifdef IND_TIMESTEPS
       rhomax_ibin = ibin(ipart_rhomax)
#endif
    else
       ipart_rhomax = -1
    endif
    call bcast_mpi(rhomax_ipart,id_rhomax)
    call bcast_mpi(rhomax_xyzh,id_rhomax)
    call bcast_mpi(rhomax_vxyz,id_rhomax)
    call bcast_mpi(rhomax_iphase,id_rhomax)
    call bcast_mpi(rhomax_divv,id_rhomax)
#ifdef IND_TIMESTEPS
    call bcast_mpi(rhomax_ibin,id_rhomax)
#endif
 endif
 if (icreate_sinks > 0 .and. ipart_rhomax > 0 .and. iverbose>=1) then
    print*,' got rhomax = ',rhomax*unit_density,' on particle ',ipart_rhomax !,rhoh(xyzh(4,ipart_rhomax))
 endif
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
 if (ndustres > 0) call summary_variable('dust',iosumdgr,ndustres,dustresfacmean /real(ndustres),dustresfacmax )
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

 !  Print warning statements, if required
 if (iverbose >= 1 .and. id==master) then
    if (ndtforce   > 0) write(iprint,*) 'force controlling timestep on ',ndtforce,' gas particles'
    if (ndtforceng > 0) write(iprint,*) 'force controlling timestep on ',ndtforce,' non-gas particles'
    if (ndtcool    > 0) write(iprint,*) 'cooling controlling timestep on ',ndtcool,' particles'
    if (ndtdrag    > 0) write(iprint,*) 'drag controlling timestep on ',ndtdrag,' gas particles'
    if (ndtdragd   > 0) write(iprint,*) 'drag controlling timestep on ',ndtdrag,' dust particles'
    if (ndtdust    > 0) write(iprint,*) 'dust diffusion controlling timestep on ',ndtdust,' particles'
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
    if (min(dtvisc,dtohm,dthall,dtambi) < dtcourant) then
       dtcourant = min(dtvisc,dtohm,dthall,dtambi)
       if      (abs(dtcourant-dtvisc) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','viscosity constraining Courant timestep')
          call summary_variable('dt',iosumdtv,0,0.0,0.0, .true. )
       else if (abs(dtcourant-dthall) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','Hall Effect constraining Courant timestep')
          call summary_variable('dt',iosumdth,0,0.0,0.0, .true. )
       else if (abs(dtcourant-dtohm ) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','ohmic resistivity constraining Courant timestep')
          call summary_variable('dt',iosumdto,0,0.0,0.0, .true. )
       else if (abs(dtcourant-dtambi) < tiny(dtcourant) ) then
          if (iverbose >= 1 .and. id==master) call warning('force','ambipolar diffusion constraining Courant timestep')
          call summary_variable('dt',iosumdta,0,0.0,0.0, .true. )
       endif
    endif
 else
    if (dtvisc < dtcourant) then
       dtcourant = dtvisc
       if (iverbose >= 1 .and. id==master) call warning('force','viscosity constraining Courant timestep')
       call summary_variable('dt',iosumdtv,0,0.0,0.0, .true. )
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
                          xyzh,vxyzu,Bevol,iphase,massoftype, &
                          divcurlB,eta_nimhd, temperature, &
                          dustfrac,gradh,divcurlv,alphaind, &
                          alphau,alphaB,bulkvisc,stressmax,&
                          ndrag,nstokes,nsuper,ts_min,ibinnow_m1,ibin_wake,ibin_neighi,&
                          ignoreself)
#ifdef FINVSQRT
 use fastmath,    only:finvsqrt
#endif
 use kernel,      only:grkern,cnormk,radkern2
 use part,        only:igas,idust,iboundary,iohm,ihall,iambi
 use part,        only:maxphase,iactive,iamtype,iamdust,get_partinfo
 use part,        only:mhd,maxvxyzu,maxBevol,maxstrain
 use dim,         only:maxalpha,maxp,mhd_nonideal,gravity,store_temperature
 use part,        only:rhoh,maxgradh,straintensor
 use nicil,       only:nimhd_get_jcbcb,nimhd_get_dBdt
#ifdef GRAVITY
 use kernel,      only:kernel_softening
 use ptmass,      only:ptmass_not_obscured
#endif
#ifdef PERIODIC
 use boundary,    only:dxbound,dybound,dzbound
#endif
 use dim,         only:use_dust,use_dustgrowth
 use dust,        only:grainsize,graindens
#ifdef DUST
 use dust,        only:get_ts,grainsize,graindens,idrag,icut_backreaction,ilimitdustflux
 use kernel,      only:wkern_drag,cnormk_drag
 use part,        only:dustprop
#endif
#ifdef DUSTGROWTH
 use growth,  only:get_vrelonvfrag
 use eos,         only:get_temperature,ieos
 use part,        only:xyzmh_ptmass
#endif
#ifdef IND_TIMESTEPS
 use part,        only:ibin_old
#endif
 use timestep,    only:bignumber
 use options,     only:overcleanfac,use_dustfrac
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
 integer(kind=1), intent(in)    :: iphase(:)
 real,            intent(in)    :: massoftype(:)
 real,            intent(in)    :: eta_nimhd(:,:)
 real,            intent(inout) :: temperature(:)
 real(kind=4),    intent(in)    :: alphaind(:,:)
 real(kind=4),    intent(in)    :: gradh(:,:),divcurlv(:,:)
 real,            intent(in)    :: alphau,alphaB,bulkvisc,stressmax
 integer,         intent(inout) :: ndrag,nstokes,nsuper
 real,            intent(out)   :: ts_min
 integer(kind=1), intent(out)   :: ibin_wake(:),ibin_neighi
 integer(kind=1), intent(in)    :: ibinnow_m1
 logical,         intent(in)    :: ignoreself
 integer :: l,j,n,iamtypej
 logical :: iactivej,iamgasj,iamdustj
 real    :: rij2,q2i,qi,xj,yj,zj,dx,dy,dz,runix,runiy,runiz,rij1,hfacgrkern
 real    :: grkerni,grgrkerni,dvx,dvy,dvz,projv,denij,vsigi,vsigu,dudtdissi
 real    :: projBi,projBj,dBx,dBy,dBz,dB2,projdB
 real    :: dendissterm,dBdissterm,dudtresist,dpsiterm,pmassonrhoi
 real    :: gradpi,projsxi,projsyi,projszi
 real    :: gradp,projsx,projsy,projsz,Bxj,Byj,Bzj,Bj,Bj1,psij
 real    :: dpsitermj,grkernj,grgrkernj,autermj,avBtermj,vsigj,spsoundj
 real    :: gradpj,pro2j,projsxj,projsyj,projszj,sxxj,sxyj,sxzj,syyj,syzj,szzj,psitermj,dBrhoterm
 real    :: visctermisoj,visctermanisoj,enj,tempj,hj,mrhoj5,alphaj,pmassj,rho1j,vsigBj
 real    :: rhoj,ponrhoj,prj,rhoav1
 real    :: hj1,hj21,q2j,qj,vwavej,divvj
 real    :: strainj(6)
#ifdef GRAVITY
 real    :: fmi,fmj,dsofti,dsoftj
 logical :: add_contribution
#else
 logical, parameter :: add_contribution = .true.
#endif
 real    :: phi,phii,phij,fgrav,fgravi,fgravj,termi
#ifdef DUST
 integer :: iregime
 real    :: dragterm,dragheating,wdrag,tsij(ndusttypes),dv2
 real    :: grkernav,tsj(ndusttypes),dustfracterms(ndusttypes),term
 !real    :: Dav(ndusttypes),vsigeps,depsdissterm(ndusttypes)
#ifdef DUSTGROWTH
 real    :: Sti, Ti, ri
#endif
#endif
 real    :: dBevolx,dBevoly,dBevolz,divBsymmterm,divBdiffterm
 real    :: rho21i,rho21j,Bxi,Byi,Bzi,psii,pmjrho21grkerni,pmjrho21grkernj
 real    :: auterm,avBterm,mrhoi5,vsigB
 real    :: jcbcbj(3),jcbj(3),dBnonideal(3),dBnonidealj(3),curlBi(3),curlBj(3)
 real    :: vsigavi,vsigavj
 real    :: dustfraci(ndusttypes),dustfracj(ndusttypes),tsi(ndusttypes)
 real    :: sqrtrhodustfraci(ndusttypes),sqrtrhodustfracj(ndusttypes)
 real    :: dustfracisum,dustfracjsum,epstsi,epstsj,rhogas1i,rhogasj,rhogas1j
 real    :: vwavei,rhoi,rho1i,spsoundi
 real    :: sxxi,sxyi,sxzi,syyi,syzi,szzi
 real    :: visctermiso,visctermaniso
 real    :: pri,pro2i
 real    :: etaohmi,etahalli,etaambii
 real    :: jcbcbi(3),jcbi(3)
 real    :: alphai,grainsizei,graindensi
 logical :: usej

 ! unpack
 vwavei        = xpartveci(ivwavei)
 rhoi          = xpartveci(irhoi)
 rho1i         = 1./rhoi
 spsoundi      = xpartveci(ispsoundi)
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
 curlBi(1)     = xpartveci(icurlBxi)
 curlBi(2)     = xpartveci(icurlByi)
 curlBi(3)     = xpartveci(icurlBzi)
 if (use_dustgrowth) then
    grainsizei = xpartveci(igrainsizei)
    graindensi = xpartveci(igraindensi)
 endif

 fsum(:) = 0.
 vsigmax = 0.
 pmassonrhoi = pmassi*rho1i
 hfacgrkern  = hi41*cnormk*gradhi

 ! default settings for active/phase if iphase not used
 iactivej = .true.
 iamtypej = igas
 iamgasj  = .true.
 iamdustj = .false.

 ! to find max ibin of all of i's neighbours
 ibin_neighi = 0_1

 ! dust
 ts_min = bignumber

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
 if (use_dustfrac) then
    dustfraci(:) = xpartveci(idustfraci:idustfraciend)
    dustfracisum = sum(dustfraci(:))
    if (ilimitdustflux) then
       tsi(:) = min(xpartveci(itstop:itstopend),hi/spsoundi) ! flux limiter from Ballabio et al. (2018)
    else
       tsi(:) = xpartveci(itstop:itstopend)
    endif
    epstsi = sum(dustfraci(:)*tsi(:))
!--sqrt(rho*epsilon) method
    sqrtrhodustfraci(:) = sqrt(rhoi*dustfraci(:))
!--asin(sqrt(epsilon)) method
!    sqrtrhodustfraci(:) = asin(sqrt(dustfraci(:)))
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
 alphaj     = alphai
 strainj(:) = 0.
 rhoj      = 0.
 rho1j     = 0.
 mrhoj5    = 0.
 gradpj    = 0.
 projsxj   = 0.
 projsyj   = 0.
 projszj   = 0.
 dpsitermj = 0.
 psitermj  = 0.
 vsigBj    = 0.
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
 loop_over_neighbours2: do n = 1,nneigh

    j = abs(listneigh(n))
    if ((ignoreself) .and. (i==j)) cycle loop_over_neighbours2

    if (ifilledcellcache .and. n <= maxcellcache) then
       ! positions from cache are already mod boundary
       xj = xyzcache(n,1)
       yj = xyzcache(n,2)
       zj = xyzcache(n,3)
       dx = xpartveci(ixi) - xj
       dy = xpartveci(iyi) - yj
       dz = xpartveci(izi) - zj
    else
       xj = xyzh(1,j)
       yj = xyzh(2,j)
       zj = xyzh(3,j)
       dx = xpartveci(ixi) - xj
       dy = xpartveci(iyi) - yj
       dz = xpartveci(izi) - zj
    endif
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
    q2j = rij2*hj21
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
       !         + (ykpt-xpartveci(iyi))*(ykpt-yj) &
       !         + (zkpt-xpartveci(izi))*(zkpt-zj)
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
          qi = (rij2*rij1)*hi1  ! this is qi = rij*hi1
       else
          rij1 = 0.
          qi = 0.
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
       if (use_dustfrac) usej = .true.
       if (maxvxyzu >= 4 .and. .not.gravity) usej = .true.

       !--get individual timestep/ multiphase information (querying iphase)
       if (maxphase==maxp) then
          call get_partinfo(iphase(j),iactivej,iamdustj,iamtypej)
          iamgasj = (iamtypej==igas .or. iamtypej==iboundary)
#ifdef IND_TIMESTEPS
          ! Particle j is a neighbour of an active particle;
          ! flag it to see if it needs to be woken up next step.
          if (iamtypej /= iboundary) then
! #ifndef MPI
             ibin_wake(j)  = max(ibinnow_m1,ibin_wake(j))
! #endif
             ibin_neighi = max(ibin_neighi,ibin_old(j))
          endif
#endif
       endif
       pmassj = massoftype(iamtypej)

       fgrav = 0.5*pmassj*(fgravi + fgravj)

       !  If particle is hidden by the sink, treat the neighbour as
       !  not gas; gravitational contribution will be added after the
       !  isgas if-statement
       if (.not. add_contribution) then
          iamgasj = .false.
          usej    = .false.
       endif

       !--get dv : needed for timestep and av term
       dvx = xpartveci(ivxi) - vxyzu(1,j)
       dvy = xpartveci(ivyi) - vxyzu(2,j)
       dvz = xpartveci(ivzi) - vxyzu(3,j)

       projv = dvx*runix + dvy*runiy + dvz*runiz

       if (iamgasj .and. maxvxyzu >= 4) then
          enj   = vxyzu(4,j)
          tempj = 0.0
          if (store_temperature) then
             tempj = temperature(j)
          endif
          denij = xpartveci(ieni) - enj
       else
          denij = 0.
          tempj = 0.0
       endif
       if (iamgasi .and. iamgasj) then
          !--work out vsig for timestepping and av
          vsigi   = max(vwavei - beta*projv,0.)
          vsigavi = max(alphai*vwavei - beta*projv,0.)
          if (vsigi > vsigmax) vsigmax = vsigi

          if (mhd) then
             hj   = xyzh(4,j)
             rhoj = rhoh(hj,pmassj)
             Bxj  = Bevol(1,j)*rhoj
             Byj  = Bevol(2,j)*rhoj
             Bzj  = Bevol(3,j)*rhoj

             if (maxBevol >= 4) psij = Bevol(4,j)
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

          if (iamgasj) then
             if (realviscosity .and. maxstrain==maxp) then
                divvj = divcurlv(1,j)
                strainj(:) = straintensor(:,j)
             else
                divvj = 0.
                strainj(:) = 0.
             endif
             if (use_dustfrac) then
                dustfracj(:) = dustfrac(:,j)
                dustfracjsum = sum(dustfracj(:))
                rhogasj      = rhoj*(1. - dustfracjsum)
                rhogas1j     = 1./rhogasj
!--sqrt(rho*epsilon) method
                sqrtrhodustfracj(:) = sqrt(rhoj*dustfracj(:))
!--asin(sqrt(epsilon)) method
!                sqrtrhodustfracj(:) = asin(sqrt(dustfracj(:)))
             else
                dustfracj(:) = 0.
                dustfracjsum = 0.
                rhogasj      = rhoj
                sqrtrhodustfracj(:) = 0.
             endif

             if (maxalpha==maxp)  alphaj  = alphaind(1,j)
             !
             !--calculate j terms (which were precalculated outside loop for i)
             !
             call get_P(rhoj,rho1j,xj,yj,zj,pmassj,enj,tempj,Bxj,Byj,Bzj,dustfracj, &
                        ponrhoj,pro2j,prj,spsoundj,vwavej, &
                        sxxj,sxyj,sxzj,syyj,syzj,szzj,visctermisoj,visctermanisoj, &
                        realviscosity,divvj,bulkvisc,strainj,stressmax)

             if (store_temperature) then
                vxyzu(4,j)     = enj
                temperature(j) = tempj
             endif

             mrhoj5   = 0.5*pmassj*rho1j
             autermj  = mrhoj5*alphau
             avBtermj = mrhoj5*alphaB*rho1j

             vsigj = max(vwavej - beta*projv,0.)
             vsigavj = max(alphaj*vwavej - beta*projv,0.)
             if (vsigj > vsigmax) vsigmax = vsigj
          else
             vsigj = max(-projv,0.)
             if (vsigj > vsigmax) vsigmax = vsigj
             vsigavj = 0.; vwavej = 0.; avBtermj = 0.; autermj = 0. ! avoid compiler warnings
             sxxj = 0.; sxyj = 0.; sxzj = 0.; syyj = 0.; syzj = 0.; szzj = 0.; pro2j = 0.; prj = 0.
             dustfracj = 0.; sqrtrhodustfracj = 0.
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
          vwavej = 0.
          vsigavj = 0.
          spsoundj = 0.
          dustfracj = 0.
          dustfracjsum = 0.
          sqrtrhodustfracj = 0.
       endif

       ifgas: if (iamgasi .and. iamgasj) then

          !
          !--artificial viscosity term
          !
#ifdef DISC_VISCOSITY
          !
          !--This is for "physical" disc viscosity
          !  (We multiply by h/rij, use cs for the signal speed, apply to both approaching/receding,
          !   with beta viscosity only applied to approaching pairs)
          !
          if (projv < 0.) then
             gradpi = pmassj*(pro2i - 0.5*rho1i*(alphai*spsoundi - beta*projv)*hi*rij1*projv)*grkerni
             if (usej) gradpj = pmassj*(pro2j - 0.5*rho1j*(alphaj*spsoundj - beta*projv)*hj*rij1*projv)*grkernj
          else
             gradpi = pmassj*(pro2i - 0.5*rho1i*alphai*spsoundi*hi*rij1*projv)*grkerni
             if (usej) gradpj = pmassj*(pro2j - 0.5*rho1j*alphaj*spsoundj*hj*rij1*projv)*grkernj
          endif
          dudtdissi = -0.5*pmassj*rho1i*alphai*spsoundi*hi*rij1*projv**2*grkerni
#else
          if (projv < 0.) then
             !--add av term to pressure
             gradpi = pmassj*(pro2i - 0.5*rho1i*vsigavi*projv)*grkerni
             if (usej) gradpj = pmassj*(pro2j - 0.5*rho1j*vsigavj*projv)*grkernj

             !--energy conservation from artificial viscosity (don't need j term)
             dudtdissi = -0.5*pmassj*rho1i*vsigavi*projv**2*grkerni
          else
             gradpi = pmassj*pro2i*grkerni
             if (usej) gradpj = pmassj*pro2j*grkernj
             dudtdissi = 0.
          endif
#endif
          !--artificial thermal conductivity (need j term)
          if (maxvxyzu >= 4) then
             if (gravity) then
                vsigu = abs(projv)
             else
                rhoav1 = 2./(rhoi + rhoj)
                vsigu = sqrt(abs(pri - prj)*rhoav1)  !abs(projv) !sqrt(abs(denij))
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
             ! new cleaning evolving d/dt (psi/c_h) as in Tricco, Price & Bate (2016)
             !
             if (maxBevol >= 4) dpsiterm = overcleanfac*(pmjrho21grkerni*psii*vwavei + pmjrho21grkernj*psij*vwavej)
             !
             ! non-ideal MHD terms
             !
             if (mhd_nonideal) then
                call nimhd_get_dBdt(dBnonideal,etaohmi,etahalli,etaambii,curlBi, &
                                  jcbi,jcbcbi,runix,runiy,runiz)
                dBnonideal = dBnonideal*pmjrho21grkerni
                if (usej) then
                   Bj  = sqrt(Bxj**2 + Byj**2 + Bzj**2)
                   if (Bj > 0.0) then
                      Bj1 = 1.0/Bj
                   else
                      Bj1 = 0.0
                   endif
                   curlBj = divcurlB(2:4,j)
                   call nimhd_get_jcbcb(jcbcbj,jcbj,curlBj,Bxj,Byj,Bzj,Bj1)
                   call nimhd_get_dBdt(dBnonidealj,eta_nimhd(iohm,j),eta_nimhd(ihall,j),eta_nimhd(iambi,j) &
                                      ,curlBj,jcbj,jcbcbj,runix,runiy,runiz)
                   dBnonideal = dBnonideal + dBnonidealj*pmjrho21grkernj
                endif
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
          if (realviscosity .and. maxstrain /= maxp) then
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
#ifdef DUST
          if (use_dustfrac) then
             do l = 1,ndusttypes
                   ! get stopping time - for one fluid dust we do not know deltav, but it is small by definition
                   call get_ts(idrag,grainsize(l),graindens,rhogasj,rhoj*dustfracjsum,spsoundj,0.,tsj(l),iregime)
             enddo
             if (ilimitdustflux) tsj(:)   = min(tsj(:),hj/spsoundj) ! flux limiter from Ballabio et al. (2018)
             epstsj   = sum(dustfracj(:)*tsj(:))
             rhogas1i = rho1i/(1.-dustfracisum)
             rhogas1j = 1./rhogasj

             !! Check that weighted sums of Tsj and tilde(Tsj) are equal (see Hutchison et al. 2017)
             !if (ndusttypes>1) then
             !   if (abs(sum(dustfracj*tsj) - sum(dustfracj*(tsj-epstsj))/(1. - dustfracjsum)) > 1e-14) &
             !      print*,'Stop! tsj or epstsj in force is incorrect!'
             !else
             !   if (abs(tsj(1) - (tsj(1)-epstsj)/(1.-dustfracj(1)))>1e-14) &
             !      print*,'Warning! error in tsj is = ',tsj(1)-(tsj(1)-epstsj)/(1.-dustfracjsum)
             !endif

             do l = 1,ndusttypes
                if (dustfraci(l) > 0. .or. dustfracj(l) > 0.) then
                   ! define averages of diffusion coefficient and kernels
                   !Dav(l)   = dustfraci(l)*tsi(l) + dustfracj(l)*tsj(l)
                   grkernav = 0.5*(grkerni + grkernj)

                   ! these are equations (43) and (45) from Price & Laibe (2015)
                   ! but note there is a sign error in the term in eqn (45) in the paper
                   !dustfracterm(l)  = pmassj*rho1j*Dav(:)*(pri - prj)*grkernav*rij1
!--sqrt(rho*epsilon) method
                   dustfracterms(l) = pmassj*sqrtrhodustfracj(l)*rho1j                     &
                                      *((tsi(l)-epstsi)*rhogas1i+(tsj(l)-epstsj)*rhogas1j) &
                                      *(pri - prj)*grkernav*rij1
!--asin(sqrt(epsilon)) method
!                   dustfracterms(l) = pmassj*sin(sqrtrhodustfracj(l))     &
!                                      *( (tsi(l)-epstsi)*rhogas1i*rho1j   &
!                                        +(tsj(l)-epstsj)*rhogas1j*rho1i ) &
!                                      *(pri - prj)*grkernav*rij1
!                   if (sqrtrhodustfraci(l) == 0.) then
!                      dustfracterms(l) = dustfracterms(l)/(2.*cos(sqrtrhodustfraci(l)))
!                   else
!                      dustfracterms(l) = dustfracterms(l)*sin(sqrtrhodustfraci(l)) &
!                                         /sin(2.*sqrtrhodustfraci(l))
!                      if (sin(2.*sqrtrhodustfraci(l)) == 0. ) stop 'dividing by zero'
!                   endif

                   !vsigeps = 0.5*(spsoundi + spsoundj) !abs(projv)
                   !depsdissterm(l) = pmassj*sqrtrhodustfracj(l)*rho1j*grkernav*vsigeps !(auterm*grkerni + autermj*grkernj)*vsigeps
                   !dustfracterms(l) = dustfracterms(l) - depsdissterm(l)*(dustfraci(l) - dustfracj(l))! &
                   !!*1.e-1/(dustfraci(l) + dustfracj(l))

                   fsum(iddustfraci+(l-1)) = fsum(iddustfraci+(l-1)) - dustfracterms(l)
                   !fsum(iddustfraci+(l-1)) = fsum(iddustfraci+(l-1)) - dustfracterm(l)
!--sqrt(rho*epsilon) method
                   if (maxvxyzu >= 4) fsum(idudtdusti+(l-1)) = fsum(idudtdusti+(l-1)) - sqrtrhodustfraci(l)*dustfracterms(l)*denij
!--asin(sqrt(epsilon)) method
!                   if (maxvxyzu >= 4) fsum(idudtdusti+(l-1)) = fsum(idudtdusti+(l-1)) &
!                                      - dustfracterms(l)*sin(2.*sqrtrhodustfraci(l))/rho1i*denij
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
#endif

        else !ifgas
        !
        !  gravity between particles of different types, or between gas pairs that are hidden by a sink
        !
           fsum(ifxi) = fsum(ifxi) - fgrav*runix
           fsum(ifyi) = fsum(ifyi) - fgrav*runiy
           fsum(ifzi) = fsum(ifzi) - fgrav*runiz
           fsum(ipot) = fsum(ipot) + pmassj*phii ! no need to symmetrise (see PM07)
#ifdef DUST
          !
          ! gas-dust: compute drag terms
          !
          if (idrag>0) then
             if (iamgasi .and. iamdustj .and. icut_backreaction==0) then
                dv2 = dvx*dvx + dvy*dvy + dvz*dvz
                if (q2i < q2j) then
                   wdrag = wkern_drag(q2i,qi)*hi21*hi1*cnormk_drag
                else
                   wdrag = wkern_drag(q2j,qj)*hj21*hj1*cnormk_drag
                endif
                if (use_dustgrowth) then
                   call get_ts(idrag,dustprop(1,j),dustprop(2,j),rhoi,rhoj,spsoundi,dv2,tsij(1),iregime)
                else
                   do l = 1,ndusttypes
                      call get_ts(idrag,grainsize(l),graindens,rhoi,rhoj,spsoundi,dv2,tsij(l),iregime)
                   enddo
                endif
                ndrag = ndrag + 1
                if (iregime > 2)  nstokes = nstokes + 1
                if (iregime == 2) nsuper = nsuper + 1
                dragterm = sum(3.*pmassj/((rhoi + rhoj)*tsij(:))*projv*wdrag)
                ts_min = min(ts_min,minval(tsij(:)))
                fsum(ifxi) = fsum(ifxi) - dragterm*runix
                fsum(ifyi) = fsum(ifyi) - dragterm*runiy
                fsum(ifzi) = fsum(ifzi) - dragterm*runiz
                if (maxvxyzu >= 4) then
                   !--energy dissipation due to drag
                   dragheating = dragterm*projv
                   fsum(idudtdissi) = fsum(idudtdissi) + dragheating
                endif
             elseif (iamdusti .and. iamgasj) then
                dv2 = dvx*dvx + dvy*dvy + dvz*dvz
                if (q2i < q2j) then
                   wdrag = wkern_drag(q2i,qi)*hi21*hi1*cnormk_drag
                else
                   wdrag = wkern_drag(q2j,qj)*hj21*hj1*cnormk_drag
                endif
                if (use_dustgrowth) then
                   call get_ts(idrag,grainsizei,graindensi,rhoj,rhoi,spsoundj,dv2,tsij(1),iregime)
                else
                   do l = 1,ndusttypes
                      call get_ts(idrag,grainsize(l),graindens,rhoj,rhoi,spsoundj,dv2,tsij(l),iregime)
                   enddo
                endif
#ifdef DUSTGROWTH
                if (usej) then
                   dustprop(5,i) = dustprop(5,i) + 3*pmassj/rhoj*projv*wdrag !--interpolate vd-vg for the dust particle i
                   ri = sqrt(xyzh(1,i)**2+xyzh(2,i)**2)
                   Sti = tsij*sqrt(xyzmh_ptmass(4,1)/ri**3) !--G=1 in code units
                   Ti = get_temperature(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
                   call get_vrelonvfrag(xyzh(:,i),dustprop(:,i),spsoundi,Sti,Ti) !--store vrel and vrel/vfrag
                endif
#endif
                dragterm = sum(3.*pmassj/((rhoi + rhoj)*tsij(:))*projv*wdrag)
                ts_min = min(ts_min,minval(tsij(:)))
                ndrag = ndrag + 1
                if (iregime > 2)  nstokes = nstokes + 1
                if (iregime == 2) nsuper = nsuper + 1
                fsum(ifxi) = fsum(ifxi) - dragterm*runix ! + because projv is opposite
                fsum(ifyi) = fsum(ifyi) - dragterm*runiy
                fsum(ifzi) = fsum(ifzi) - dragterm*runiz
             endif
          endif
#endif
       endif ifgas
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
       pmassj = massoftype(iamtypej)
       phii   = -rij1
       fgravj = fgrav*pmassj
       fsum(ifxi) = fsum(ifxi) - dx*fgravj
       fsum(ifyi) = fsum(ifyi) - dy*fgravj
       fsum(ifzi) = fsum(ifzi) - dz*fgravj
       fsum(ipot) = fsum(ipot) + pmassj*phii
#endif
    endif is_sph_neighbour
 enddo loop_over_neighbours2

 return
end subroutine compute_forces

!----------------------------------------------------------------
!+
!  Internal subroutine that computes pressure and other derived
!  quantities necessary to get a force, given that we have rho.
!+
!----------------------------------------------------------------
subroutine get_P(rhoi,rho1i,xi,yi,zi,pmassi,eni,tempi,Bxi,Byi,Bzi,dustfraci, &
                 ponrhoi,pro2i,pri,spsoundi,vwavei, &
                 sxxi,sxyi,sxzi,syyi,syzi,szzi,visctermiso,visctermaniso, &
                 realviscosity,divvi,bulkvisc,strain,stressmax)

 use dim,       only:maxvxyzu,maxstrain,maxp,store_temperature
 use part,      only:mhd
 use eos,       only:equationofstate
 use options,   only:ieos
 use viscosity, only:shearfunc
 real,    intent(in)    :: rhoi,rho1i,xi,yi,zi,pmassi
 real,    intent(inout) :: eni,tempi
 real,    intent(in)    :: Bxi,Byi,Bzi,dustfraci(:)
 real,    intent(out)   :: ponrhoi,pro2i,pri,spsoundi,vwavei
 real,    intent(out)   :: sxxi,sxyi,sxzi,syyi,syzi,szzi
 real,    intent(out)   :: visctermiso,visctermaniso
 logical, intent(in)    :: realviscosity
 real,    intent(in)    :: divvi,bulkvisc,stressmax
 real,    intent(in)    :: strain(6)

 real :: Bro2i,Brhoxi,Brhoyi,Brhozi,rhogasi,gasfrac
 real :: stressiso,term,graddivvcoeff,del2vcoeff
 real :: shearvisc,etavisc,valfven2i,p_on_rhogas
!
!--get pressure (actually pr/dens) and sound speed from equation of state
!
 gasfrac = (1. - sum(dustfraci))  ! rhogas/rho
 rhogasi = rhoi*gasfrac       ! rhogas = (1-eps)*rho
 if (maxvxyzu >= 4) then
    if (store_temperature) then
       call equationofstate(ieos,p_on_rhogas,spsoundi,rhogasi,xi,yi,zi,eni,tempi)
    else
       call equationofstate(ieos,p_on_rhogas,spsoundi,rhogasi,xi,yi,zi,eni)
    endif
 else
    call equationofstate(ieos,p_on_rhogas,spsoundi,rhogasi,xi,yi,zi)
 endif

 pri     = p_on_rhogas*rhogasi
 ponrhoi = p_on_rhogas*gasfrac

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
    if (maxstrain==maxp) then
       !--get stress (multiply by coefficient for use in second derivative)
       term = -shearvisc*pmassi*rho1i ! shearvisc = eta/rho, so this is eta/rho**2
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
    pro2i = ponrhoi*rho1i + stressiso + 0.5*Bro2i

 else
!
!--construct m*p/(rho^2 \Omega) in force equation using pressure
!
    pro2i  = ponrhoi*rho1i + stressiso
    vwavei = spsoundi
 endif

 return
end subroutine get_P

#ifdef IND_TIMESTEPS
!----------------------------------------------------------------
!+
!  Checks which timestep is the limiting dt.  Book keeping is done here
!+
!----------------------------------------------------------------
subroutine check_dtmin(dtcheck,dti,dtopt,dtrat,ndtopt,dtoptfacmean,dtoptfacmax,dtchar_out,dtchar_in)
 integer, intent(inout) :: ndtopt
 real,    intent(in)    :: dti,dtopt,dtrat
 real,    intent(inout) :: dtoptfacmean,dtoptfacmax
 logical, intent(inout) :: dtcheck
 character(len=*), intent(out)   :: dtchar_out
 character(len=*), intent(in)    :: dtchar_in
 !
 if (.not. dtcheck) return
 !
 if ( abs(dti-dtopt) < tiny(dti)) then
    dtcheck      = .false.
    ndtopt       = ndtopt + 1
    dtoptfacmean = dtoptfacmean + dtrat
    dtoptfacmax  = max(dtoptfacmax, dtrat)
    dtchar_out   = dtchar_in
 endif
 !
 return
end subroutine check_dtmin
#endif

!----------------------------------------------------------------

subroutine start_cell(cell,iphase,xyzh,vxyzu,gradh,divcurlv,divcurlB,straintensor,Bevol, &
                     dustfrac,dustprop,eta_nimhd,temperature,alphaind,stressmax)

 use io,        only:fatal
 use options,   only:alpha,use_dustfrac
 use dim,       only:maxp,ndivcurlv,ndivcurlB,maxstrain,maxalpha,maxvxyzu,mhd,mhd_nonideal,&
                use_dustgrowth,store_temperature
 use part,      only:iamgas,maxphase,iboundary,rhoanddhdrho,igas,massoftype,get_partinfo,&
                     iohm,ihall,iambi
 use viscosity, only:irealvisc,bulkvisc
#ifdef DUST
 use dust,      only:get_ts,grainsize,graindens,idrag
#endif
 use nicil,     only:nimhd_get_dt,nimhd_get_jcbcb
 use kdtree,    only:inodeparts,inoderange
 type(cellforce),    intent(inout) :: cell
 integer(kind=1),    intent(in)    :: iphase(:)
 real,               intent(in)    :: xyzh(:,:)
 real,               intent(inout) :: vxyzu(:,:)
 real(kind=4),       intent(in)    :: gradh(:,:)
 real(kind=4),       intent(in)    :: divcurlv(:,:)
 real(kind=4),       intent(in)    :: divcurlB(:,:)
 real(kind=4),       intent(in)    :: straintensor(:,:)
 real,               intent(in)    :: Bevol(:,:)
 real,               intent(in)    :: dustfrac(:,:),dustprop(:,:)
 real,               intent(in)    :: eta_nimhd(:,:)
 real,               intent(inout) :: temperature(:)
 real(kind=4),       intent(in)    :: alphaind(:,:)
 real,               intent(in)    :: stressmax

 real         :: divcurlvi(ndivcurlv)
 real         :: straini(6),curlBi(3),jcbcbi(3),jcbi(3)
 real         :: hi,rhoi,rho1i,dhdrhoi,pmassi,eni,tempi
 real(kind=8) :: hi1
 real         :: dustfraci(ndusttypes),dustfracisum,rhogasi,ponrhoi,pro2i,pri,spsoundi
 real         :: sxxi,sxyi,sxzi,syyi,syzi,szzi,visctermiso,visctermaniso
#ifdef DUST
 real         :: tstopi(ndusttypes)
#endif
 real         :: Bxi,Byi,Bzi,Bi,B2i,Bi1
 real         :: vwavei,alphai

 integer      :: i,j,iamtypei,ip

#ifdef DUST
 integer :: iregime
#endif

 logical :: iactivei,iamgasi,iamdusti,realviscosity

 realviscosity = (irealvisc > 0)

 cell%npcell = 0
 over_parts: do ip = inoderange(1,cell%icell),inoderange(2,cell%icell)
    i = inodeparts(ip)

    if (i < 0) then
       cycle over_parts
    endif

    if (maxphase==maxp) then
       call get_partinfo(iphase(i),iactivei,iamdusti,iamtypei)
       iamgasi = (iamtypei==igas)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif
    if (.not.iactivei) then ! handles case where first particle in cell is inactive
       cycle over_parts
    endif
    if (iamtypei==iboundary) then ! do not compute forces on boundary parts
       cycle over_parts
    endif

    pmassi = massoftype(iamtypei)
    hi = xyzh(4,i)
    if (hi < 0.) call fatal('force','negative smoothing length',i,var='h',val=hi)

    !
    !--compute density and related quantities from the smoothing length
    !
    call rhoanddhdrho(hi,hi1,rhoi,rho1i,dhdrhoi,pmassi)

    if (iamgasi) then
       if (ndivcurlv >= 1) divcurlvi(:) = real(divcurlv(:,i),kind=kind(divcurlvi))
       if (realviscosity .and. maxstrain==maxp) straini(:) = straintensor(:,i)
       if (maxvxyzu >= 4) then
          eni   = vxyzu(4,i)
          tempi = 0.0
          if (store_temperature) then
             tempi = temperature(i)
          endif
       else
          eni   = 0.0
          tempi = 0.0
       endif

       !
       ! one-fluid dust properties
       !
       if (use_dustfrac) then
          dustfraci(:) = dustfrac(:,i)
          dustfracisum = sum(dustfraci)
          rhogasi      = rhoi*(1. - dustfracisum)
          do j = 1,ndusttypes
             if (dustfraci(j) > 1. .or. dustfraci(j) < 0.) call fatal('force','invalid eps',var='dustfrac',val=dustfraci(j))
          enddo
       else
          dustfraci(:) = 0.
          rhogasi      = rhoi
       endif

       if (mhd) then
          Bxi = Bevol(1,i) * rhoi ! B/rho -> B (conservative to primitive)
          Byi = Bevol(2,i) * rhoi
          Bzi = Bevol(3,i) * rhoi
       endif

       !
       ! calculate terms required in the force evaluation
       !
       call get_P(rhoi,rho1i, &
                  xyzh(1,i),xyzh(2,i),xyzh(3,i), &
                  pmassi, &
                  eni, tempi, &
                  Bxi,Byi,Bzi, &
                  dustfraci(:), &
                  ponrhoi,pro2i,pri,spsoundi, &
                  vwavei,sxxi,sxyi,sxzi,syyi,syzi,szzi, &
                  visctermiso,visctermaniso,realviscosity,divcurlvi(1),bulkvisc,straini,stressmax)
#ifdef DUST
       !
       ! get stopping time - for one fluid dust we don't know deltav, but as small by definition we assume=0
       !
       if (use_dustfrac .and. iamgasi) then
          do j = 1,ndusttypes
             call get_ts(idrag,grainsize(j),graindens,rhogasi,rhoi*dustfracisum,spsoundi,0.,tstopi(j),iregime)
          enddo
       endif
#endif

       if (store_temperature) then
          vxyzu(4,i)     = eni
          temperature(i) = tempi
       endif

       if (mhd_nonideal) then
          B2i = Bxi**2 + Byi**2 + Bzi**2
          Bi  = sqrt(B2i)
          if (Bi > 0.0) then
             Bi1 = 1.0/Bi
          else
             Bi1 = 0.0
          endif
          curlBi = divcurlB(2:4,i)
          call nimhd_get_jcbcb(jcbcbi,jcbi,curlBi,Bxi,Byi,Bzi,Bi1)
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

          if (maxBevol >= 4) then
             cell%xpartvec(ipsi,cell%npcell)       = Bevol(4,i)
          endif
          if (maxBevol < 3 .or. maxBevol > 4) call fatal('densityiterate','error in maxBevol setting')

          cell%xpartvec(icurlBxi,cell%npcell)      = divcurlB(2,i)
          cell%xpartvec(icurlByi,cell%npcell)      = divcurlB(3,i)
          cell%xpartvec(icurlBzi,cell%npcell)      = divcurlB(4,i)
       else
          cell%xpartvec(iBevolxi:ipsi,cell%npcell) = 0. ! to avoid compiler warning
       endif
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
       cell%xpartvec(iponrhoi,cell%npcell)        = ponrhoi
       cell%xpartvec(ispsoundi,cell%npcell)       = spsoundi
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
#ifdef DUST
       if (use_dustfrac) then
          cell%xpartvec(itstop:itstopend,cell%npcell)       = tstopi
       endif
#endif
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
#ifdef DUSTGROWTH
    cell%xpartvec(igrainsizei,cell%npcell)        = dustprop(1,i)
    cell%xpartvec(igraindensi,cell%npcell)        = dustprop(2,i) !--Add dustprop(3,i) for vrelonvfrag
#endif
 enddo over_parts

end subroutine start_cell

subroutine compute_cell(cell,listneigh,nneigh,Bevol,xyzh,vxyzu,fxyzu, &
                        iphase,divcurlv,divcurlB,alphaind,eta_nimhd, temperature, &
                        dustfrac,gradh,ibinnow_m1,ibin_wake,stressmax,xyzcache)
 use io,          only:error
#ifdef MPI
 use io,          only:id
#endif
 use dim,         only:maxvxyzu
 use options,     only:beta,alphau,alphaB,iresistive_heating
 use part,        only:get_partinfo,iamgas,iboundary,mhd,igas,maxphase,massoftype
 use viscosity,   only:irealvisc,bulkvisc
 use kdtree,      only:inodeparts

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
 real,            intent(inout)  :: temperature(:)
 real(kind=4),    intent(in)     :: gradh(:,:)
 integer(kind=1), intent(inout)  :: ibin_wake(:)
 integer(kind=1), intent(in)     :: ibinnow_m1
 real,            intent(in)     :: stressmax
 real,            intent(in)     :: xyzcache(:,:)

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
       call get_partinfo(cell%iphase(ip),iactivei,iamdusti,iamtypei)
       iamgasi = (iamtypei==igas)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif

    if (.not.iactivei) then ! handles case where first particle in cell is inactive
       cycle over_parts
    endif
    if (iamtypei==iboundary) then ! do not compute forces on boundary parts
       cycle over_parts
    endif

    i = inodeparts(cell%arr_index(ip))

    pmassi = massoftype(iamtypei)

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
#ifdef MPI
    if (cell%owner == id) then
       ignoreself = .true.
    else
       ignoreself = .false.
    endif
#else
    ignoreself = .true.
#endif

    call compute_forces(i,iamgasi,iamdusti,cell%xpartvec(:,ip),hi,hi1,hi21,hi41,gradhi,gradsofti, &
                         beta, &
                         pmassi,listneigh,nneigh,xyzcache,cell%fsums(:,ip),cell%vsigmax(ip), &
                         .true.,realviscosity,useresistiveheat, &
                         xyzh,vxyzu,Bevol,iphase,massoftype, &
                         divcurlB,eta_nimhd, temperature, &
                         dustfrac,gradh,divcurlv,alphaind, &
                         alphau,alphaB,bulkvisc,stressmax, &
                         cell%ndrag,cell%nstokes,cell%nsuper,cell%dtdrag(ip),ibinnow_m1,ibin_wake,cell%ibinneigh(ip), &
                         ignoreself)

 enddo over_parts

end subroutine compute_cell

subroutine finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,straintensor,&
                                         divBsymm,divcurlv,dBevol,ddustfrac,deltav, &
                                         dtcourant,dtforce,dtvisc,dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi, &
#ifdef IND_TIMESTEPS
                                         nbinmaxnew,nbinmaxstsnew,ncheckbin, &
                                         ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd, &
                                         ndtvisc,ndtohm,ndthall,ndtambi,ndtdust, &
                                         dtitmp,dtrat, &
                                         dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean, &
                                         dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax, &
                                         dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean ,dtdustfacmean, &
                                         dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax  ,dtdustfacmax, &
#endif
                                         ndustres,dustresfacmax,dustresfacmean)
 use io,             only:fatal
#ifdef FINVSQRT
 use fastmath,       only:finvsqrt
#endif
 use dim,            only:mhd,mhd_nonideal,lightcurve,use_dust,maxstrain
 use eos,            only:use_entropy,gamma
 use options, only:ishock_heating,icooling,psidecayfac,overcleanfac,alpha,ipdv_heating,use_dustfrac
 use part,           only:h2chemistry,rhoanddhdrho,abundance,iboundary,igas,maxphase,maxvxyzu,nabundances, &
                          massoftype,get_partinfo,tstop
#ifdef IND_TIMESTEPS
 use part,           only:ibin
 use timestep_ind,   only:get_newbin
 use timestep_sts,   only:sts_it_n,ibin_sts
#endif
 use viscosity,      only:bulkvisc,dt_viscosity,irealvisc,shearfunc
 use kernel,         only:kernel_softening
 use linklist,       only:get_distance_from_centre_of_mass
 use kdtree,         only:expand_fgrav_in_taylor_series,inodeparts
 use nicil,          only:nimhd_get_dudt,nimhd_get_dt
 use cooling,        only:energ_cooling
 use chem,           only:energ_h2cooling
 use timestep,       only:C_cour,C_cool,C_force,bignumber,dtmax
 use timestep_sts,   only:use_sts
#ifdef LIGHTCURVE
 use part,           only:luminosity
#endif

 integer,            intent(in)    :: icall
 type(cellforce),    intent(inout) :: cell
 real,               intent(inout) :: fxyzu(:,:)
 real,               intent(in)    :: xyzh(:,:)
 real,               intent(inout) :: vxyzu(:,:)
 real,               intent(in)    :: dt
 real(kind=4),       intent(in)    :: straintensor(:,:)
 real(kind=4),       intent(out)   :: poten(:)
 real(kind=4),       intent(out)   :: divBsymm(:)
 real(kind=4),       intent(out)   :: divcurlv(:,:)
 real,               intent(out)   :: dBevol(:,:)
 real,               intent(out)   :: ddustfrac(:,:)
 real,               intent(out)   :: deltav(:,:,:)

 real,               intent(inout) :: dtcourant,dtforce,dtvisc
 real,               intent(inout) :: dtohm,dthall,dtambi,dtdiff,dtmini,dtmaxi
#ifdef IND_TIMESTEPS
 integer,            intent(inout) :: nbinmaxnew,nbinmaxstsnew,ncheckbin
 integer,            intent(inout) :: ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd
 integer,            intent(inout) :: ndtvisc,ndtohm,ndthall,ndtambi,ndtdust
 real,               intent(inout) :: dtitmp,dtrat
 real,               intent(inout) :: dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean
 real,               intent(inout) :: dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax
 real,               intent(inout) :: dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean ,dtdustfacmean
 real,               intent(inout) :: dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax  ,dtdustfacmax
#endif
 integer,            intent(inout) :: ndustres
 real,               intent(inout) :: dustresfacmean,dustresfacmax

 real    :: xpartveci(maxxpartveciforce),fsum(maxfsum)
 real    :: rhoi,rho1i,rhogasi,hi,hi1,pmassi
 real    :: Bxyzi(maxBevol),curlBi(3),straini(6)
 real    :: xi,yi,zi,B2i,f2i,divBsymmi,betai,frac_divB,vcleani
 real    :: ponrhoi,spsoundi,drhodti,divvi,shearvisc,fac,pdv_work
 real    :: psii,dtau
 real    :: eni,dudtnonideal
 real    :: dustfraci(ndusttypes),dustfracisum
 real    :: tstopi(ndusttypes),tseff,dtdustdenom
 real    :: etaambii,etahalli,etaohmi
 real    :: vsigmax,vwavei,fxyz4
#ifdef LIGHTCURVE
 real    :: dudt_radi
#endif
#ifdef GRAVITY
 real    :: potensoft0,dum,dx,dy,dz,fxi,fyi,fzi,poti,epoti
#endif
 real    :: vsigdtc,dtc,dtf,dti,dtcool,dtdiffi
 real    :: dtohmi,dtambii,dthalli,dtvisci,dtdrag,dtdusti,dtclean
 integer :: idudtcool,ichem,iamtypei
 logical :: iactivei,iamgasi,iamdusti,realviscosity
#ifdef IND_TIMESTEPS
 integer(kind=1)       :: ibin_neighi
 logical               :: allow_decrease,dtcheck
 character(len=16)     :: dtchar
#endif
 integer               :: ip,i

 eni = 0.
 realviscosity = (irealvisc > 0)

 over_parts: do ip = 1,cell%npcell

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(ip),iactivei,iamdusti,iamtypei)
       iamgasi = (iamtypei==igas)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif

    if (.not.iactivei) then ! handles case where first particle in cell is inactive
       cycle over_parts
    endif
    if (iamtypei==iboundary) then ! do not compute forces on boundary parts
       cycle over_parts
    endif

    pmassi = massoftype(iamtypei)

    i = inodeparts(cell%arr_index(ip))

    fsum(:)       = cell%fsums(:,ip)
    xpartveci(:)  = cell%xpartvec(:,ip)
#ifdef IND_TIMESTEPS
    ibin_neighi   = cell%ibinneigh(ip)
#endif

    xi         = xpartveci(ixi)
    yi         = xpartveci(iyi)
    zi         = xpartveci(izi)
    hi         = xpartveci(ihi)
    hi1        = 1./hi
    spsoundi   = xpartveci(ispsoundi)
    vsigmax    = cell%vsigmax(ip)
    dtdrag     = cell%dtdrag(ip)
    tstopi     = 0.

    if (iamgasi) then
       rhoi    = xpartveci(irhoi)
       rho1i   = 1./rhoi
       rhogasi = xpartveci(irhogasi)
       ponrhoi = xpartveci(iponrhoi)
       vwavei  = xpartveci(ivwavei)

       if (mhd) then
          Bxyzi(1) = xpartveci(iBevolxi) * rhoi
          Bxyzi(2) = xpartveci(iBevolyi) * rhoi
          Bxyzi(3) = xpartveci(iBevolzi) * rhoi
          B2i      = Bxyzi(1)**2 + Bxyzi(2)**2 + Bxyzi(3)**2
       endif

       if (mhd_nonideal) then
          curlBi(1) = xpartveci(icurlBxi)
          curlBi(2) = xpartveci(icurlByi)
          curlBi(3) = xpartveci(icurlBzi)
          etaohmi   = xpartveci(ietaohmi)
          etaambii  = xpartveci(ietaambii)
          etahalli  = xpartveci(ietahalli)
       endif

       straini(:) = 0.
       if (maxvxyzu >= 4) then
          if (maxstrain == maxp .and. realviscosity) then
             straini(:) = straintensor(:,i)
          endif
          eni = xpartveci(ieni)
       endif

       if (use_dustfrac) then
          dustfraci(:) = xpartveci(idustfraci:idustfraciend)
          tstopi(:)    = xpartveci(itstop:itstopend)
       else
          dustfraci = 0.
          tstopi    = 0.
       endif

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
          betai = 2.0*ponrhoi*rhoi/B2i
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
    drhodti = pmassi*fsum(idrhodti)

    isgas: if (iamgasi) then
       divvi = -drhodti*rho1i
       if (ndivcurlv >= 1) divcurlv(1,i) = real(divvi,kind=kind(divcurlv)) ! store divv from forces

       if (maxvxyzu >= 4 .or. lightcurve) then
          if (maxstrain == maxp .and. realviscosity) then
             shearvisc = shearfunc(xi,yi,zi,spsoundi)
             fsum(idudtdissi) = fsum(idudtdissi) + (bulkvisc - 2./3.*shearvisc)*divvi**2 &
                           + 0.5*shearvisc*(straini(1)**2 + 2.*(straini(2)**2 + straini(3)**2 + straini(5)**2) &
                           + straini(4)**2 + straini(6)**2)
          endif
          fxyz4 = 0.
          if (use_entropy) then
             if (ishock_heating > 0) then
                fxyz4 = fxyz4 + (gamma - 1.)*rhoi**(1.-gamma)*fsum(idudtdissi)
             endif
          else
             fac = rhoi/rhogasi
             pdv_work = ponrhoi*rho1i*drhodti
             if (ipdv_heating > 0) then
                fxyz4 = fxyz4 + fac*pdv_work
             endif
             if (ishock_heating > 0) then
                fxyz4 = fxyz4 + fac*fsum(idudtdissi)
             endif
#ifdef LIGHTCURVE
             if (lightcurve) then
                pdv_work = ponrhoi*rho1i*drhodti
                if (pdv_work > tiny(pdv_work)) then ! pdv_work < 0 is possible, and we want to ignore this case
                   dudt_radi = fac*pdv_work + fac*fsum(idudtdissi)
                else
                   dudt_radi = fac*fsum(idudtdissi)
                endif
                luminosity(i) = pmassi*dudt_radi
             endif
#endif
             if (mhd_nonideal) then
                call nimhd_get_dudt(dudtnonideal,etaohmi,etaambii,rhoi,curlBi,Bxyzi(1:3))
                fxyz4 = fxyz4 + fac*dudtnonideal
             endif
             !--add conductivity and resistive heating
             fxyz4 = fxyz4 + fac*fsum(idendtdissi)
             if (icooling > 0) then
                if (h2chemistry) then
                   idudtcool = 1
                   ichem = 0
                   call energ_h2cooling(vxyzu(4,i),fxyz4,rhoi,&
                        abundance(:,i),nabundances,dt,xyzh(1,i),xyzh(2,i),xyzh(3,i),&
                        divcurlv(1,i),idudtcool,ichem)
                else
                   !call energ_cooling(icooling,vxyzu(4,i),fxyz4,xyzh(1,i),xyzh(2,i),xyzh(3,i))
                   call energ_cooling(icooling,vxyzu(4,i),fxyz4,xyzh(1,i),xyzh(2,i),xyzh(3,i),rhoi,vxyzu(:,i),dt)
                endif
             endif
             ! extra terms in du/dt from one fluid dust
             if (use_dustfrac) then
                !fxyz4 = fxyz4 + 0.5*fac*rho1i*fsum(idudtdusti)
                fxyz4 = fxyz4 + 0.5*fac*rho1i*sum(fsum(idudtdusti:idudtdustiend))
             endif
          endif
          if (maxvxyzu >= 4) fxyzu(4,i) = fxyz4
       endif

       dtclean = bignumber
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
          if (maxBevol >= 4 .and. psidecayfac > 0.) then
             vcleani = overcleanfac*vwavei
             dtau = psidecayfac*vcleani*hi1
             !
             ! we clean using the difference operator for div B
             !
             psii = xpartveci(ipsi)

             ! new cleaning evolving d/dt (psi/c_h)
             dBevol(4,i) = -vcleani*fsum(idivBdiffi)*rho1i - psii*dtau - 0.5*psii*divvi

             ! timestep from cleaning (should only matter if overcleaning applied)
             ! the factor of 0.5 is empirical, from checking when overcleaning with ind. timesteps is stable
             dtclean = 0.5*C_cour*hi/(vcleani + epsilon(0.))
          endif
       endif

       if (use_dustfrac) then
!--sqrt(rho*epsilon) method
          ddustfrac(:,i) = 0.5*(fsum(iddustfraci:iddustfraciend)-sqrt(rhoi*dustfraci(:))*divvi)
!--asin(sqrt(epsilon)) method
!          ddustfrac(:,i) = fsum(iddustfraci:iddustfraciend)
          deltav(1,:,i)  = fsum(ideltavxi:ideltavxiend)
          deltav(2,:,i)  = fsum(ideltavyi:ideltavyiend)
          deltav(3,:,i)  = fsum(ideltavzi:ideltavziend)
       endif

       ! timestep based on Courant condition
       vsigdtc = max(vsigmax,vwavei)
       if (vsigdtc > tiny(vsigdtc)) then
          dtc = min(C_cour*hi/(vsigdtc*max(alpha,1.0)),dtclean)
       else
          dtc = min(dtmax,dtclean)
       endif

       ! cooling timestep dt < fac*u/(du/dt)
       if (maxvxyzu >= 4 .and. icooling > 0) then
          dtcool = C_cool*abs(eni/fxyzu(4,i))
       else
          dtcool = bignumber
       endif

       ! timestep based on non-ideal MHD
       if (mhd_nonideal) then
          call nimhd_get_dt(dtohmi,dthalli,dtambii,hi,etaohmi,etahalli,etaambii)
          if ( use_STS ) then
             dtdiffi = min(dtohmi,dtambii)
             dtdiff  = min(dtdiff,dtdiffi)
             dtohmi  = bignumber
             dtambii = bignumber
          endif
       else
          dtohmi  = bignumber
          dthalli = bignumber
          dtambii = bignumber
          dtdiffi = bignumber
       endif

       ! timestep from physical viscosity
       dtvisci = dt_viscosity(xi,yi,zi,hi,spsoundi)

       ! Check to ensure we have enough resolution for gas-dust pairs, where
       ! dtdrag is already minimised over all dust neighbours for gas particle i
       if (dtdrag < bignumber) then
          if (hi > dtdrag*spsoundi) then
             ndustres       = ndustres + 1
             dustresfacmean = dustresfacmean + hi/(dtdrag*spsoundi)
             dustresfacmax  = max(dustresfacmax, hi/(dtdrag*spsoundi))
          endif
       endif

    else ! not gas

       if (maxvxyzu > 4) fxyzu(4,i) = 0.
       ! timestep based on Courant condition for non-gas particles
       vsigdtc = vsigmax
       if (vsigdtc > tiny(vsigdtc)) then
          dtc = C_cour*hi/vsigdtc
       else
          dtc = dtmax
       endif
       dtcool  = bignumber
       dtvisci = bignumber
       dtohmi  = bignumber
       dthalli = bignumber
       dtambii = bignumber
       dtdiffi = bignumber

    endif isgas

    ! initialise timestep to Courant timestep & perform sanity check
    dti = dtc
    if (dtc < tiny(dtc) .or. dtc > huge(dtc)) call fatal('force','invalid dtc',var='dtc',val=dtc)

    ! timestep based on force condition
    if (abs(f2i) > epsilon(f2i)) then
#ifdef FINVSQRT
       dtf = C_force*sqrt(hi*finvsqrt(f2i))
#else
       dtf = C_force*sqrt(hi/sqrt(f2i))
#endif
    else
       dtf = bignumber
    endif

    ! one fluid dust timestep
    if (use_dustfrac .and. iamgasi .and. minval(dustfraci) > 0. .and. spsoundi > 0.) then
       tseff = (1.-dustfracisum)/dustfracisum*sum(dustfraci(:)*tstopi(:))
       dtdustdenom = dustfracisum*tseff*spsoundi**2
       if (dtdustdenom > tiny(dtdustdenom)) then
          dtdusti = C_force*hi*hi/dtdustdenom
       else
          dtdusti = bignumber
       endif
    else
       dtdusti = bignumber
    endif

    ! stopping time
    if (use_dust .and. use_dustfrac) then
       tstop(:,i) = tstopi(:)
    elseif (use_dust .and. .not.use_dustfrac) then
       tstop(:,i) = dtdrag
    endif

#ifdef IND_TIMESTEPS
    !-- The new timestep for particle i
    dtitmp = min(dtf,dtcool,dtvisci,dtdrag,dtohmi,dthalli,dtambii,dtdusti)
    if (dtitmp < dti .and. dtitmp < dtmax) then
       dti     = dtitmp
       if (dti < tiny(dti) .or. dti > huge(dti)) call fatal('force','invalid dti',var='dti',val=dti) ! sanity check
       dtcheck = .true.
       dtrat   = dtc/dti
       if ( iamgasi ) then
          call check_dtmin(dtcheck,dti,dtf    ,dtrat,ndtforce  ,dtfrcfacmean  ,dtfrcfacmax ,dtchar,'dt_gasforce')
          call check_dtmin(dtcheck,dti,dtcool ,dtrat,ndtcool   ,dtcoolfacmean ,dtcoolfacmax,dtchar,'dt_cool'    )
          call check_dtmin(dtcheck,dti,dtvisci,dtrat,ndtvisc   ,dtviscfacmean ,dtviscfacmax,dtchar,'dt_visc'    )
          call check_dtmin(dtcheck,dti,dtdrag ,dtrat,ndtdrag   ,dtdragfacmean ,dtdragfacmax,dtchar,'dt_gasdrag' )
          call check_dtmin(dtcheck,dti,dtohmi ,dtrat,ndtohm    ,dtohmfacmean  ,dtohmfacmax ,dtchar,'dt_ohm'     )
          call check_dtmin(dtcheck,dti,dthalli,dtrat,ndthall   ,dthallfacmean ,dthallfacmax,dtchar,'dt_hall'    )
          call check_dtmin(dtcheck,dti,dtambii,dtrat,ndtambi   ,dtambifacmean ,dtambifacmax,dtchar,'dt_ambi'    )
          call check_dtmin(dtcheck,dti,dtdusti,dtrat,ndtdust   ,dtdustfacmean ,dtdustfacmax,dtchar,'dt_dust'    )
       else
          call check_dtmin(dtcheck,dti,dtf    ,dtrat,ndtforceng,dtfrcngfacmean,dtfrcngfacmax,dtchar,'dt_force'  )
          call check_dtmin(dtcheck,dti,dtdrag ,dtrat,ndtdragd  ,dtdragdfacmean,dtdragdfacmax,dtchar,'dt_drag'   )
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
    ! find the new maximum number of bins
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
    dtforce   = min(dtforce,dtf,dtcool,dtdrag,dtdusti)
    dtvisc    = min(dtvisc,dtvisci)
    if (mhd_nonideal .and. iamgasi) then
       dtohm  = min(dtohm,  dtohmi  )
       dthall = min(dthall, dthalli )
       dtambi = min(dtambi, dtambii )
    endif
    dtmini = min(dtmini,dti)
    dtmaxi = max(dtmaxi,dti)
#endif

 enddo over_parts
end subroutine finish_cell_and_store_results

#ifdef MPI
pure subroutine combine_cells(cella, cellb)
 type(cellforce),   intent(inout)        :: cella
 type(cellforce),   intent(in)           :: cellb

 integer                                 :: i

 do i = 1,cella%npcell
    cella%fsums(:,i) = cella%fsums(:,i) + cellb%fsums(:,i)
 enddo

 cella%ndrag   = cella%ndrag     + cellb%ndrag
 cella%nstokes = cella%nstokes   + cellb%nstokes
 cella%nsuper  = cella%nsuper    + cellb%nsuper

 cella%remote_export = (cella%remote_export .and. cellb%remote_export)

end subroutine combine_cells
#endif

end module forces
