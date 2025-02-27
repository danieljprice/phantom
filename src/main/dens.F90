!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module densityforce
!
! This module is the "guts" of the code
!  Calculates density by iteration with smoothing length
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, io, io_summary, kdtree, kernel, linklist,
!   mpidens, mpiderivs, mpimemory, mpiutils, omputils, options, part,
!   timestep, timing, viscosity
!
 use dim,     only:maxdvdx,maxp,maxrhosum,maxdustlarge
 use dim,     only:calculate_density,calculate_divcurlB
 use kdtree,  only:inodeparts,inoderange
 use kernel,  only:cnormk,wab0,gradh0,dphidh0,radkern2
 use mpidens, only:celldens,stackdens
 use timing,  only:getused,printused,print_time

 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 public :: densityiterate,get_neighbour_stats

 !--indexing for xpartveci array
 integer, parameter :: &
       ixi  = 1, &
       iyi  = 2, &
       izi  = 3, &
       ivxi = 4, &
       ivyi = 5, &
       ivzi = 6, &
       ieni = 7, &
       iBevolxi = 8, &
       iBevolyi = 9, &
       iBevolzi = 10, &
       ipsi = 11, &
       ifxi = 12, &
       ifyi = 13, &
       ifzi = 14, &
       iradxii = 15

 !--indexing for rhosum array
 integer, parameter :: &
       irhoi            = 1, &
       igradhi          = 2, &
       igradsofti       = 3, &
       idivvi           = 4, &
       idvxdxi          = 5, &
       idvxdyi          = 6, &
       idvxdzi          = 7, &
       idvydxi          = 8, &
       idvydyi          = 9, &
       idvydzi          = 10, &
       idvzdxi          = 11, &
       idvzdyi          = 12, &
       idvzdzi          = 13, &
       idaxdxi          = 14, &
       idaxdyi          = 15, &
       idaxdzi          = 16, &
       idaydxi          = 17, &
       idaydyi          = 18, &
       idaydzi          = 19, &
       idazdxi          = 20, &
       idazdyi          = 21, &
       idazdzi          = 22, &
       irxxi            = 23, &
       irxyi            = 24, &
       irxzi            = 25, &
       iryyi            = 26, &
       iryzi            = 27, &
       irzzi            = 28, &
       idivBi           = 29, &
       idBxdxi          = 30, &
       idBxdyi          = 31, &
       idBxdzi          = 32, &
       idBydxi          = 33, &
       idBydyi          = 34, &
       idBydzi          = 35, &
       idBzdxi          = 36, &
       idBzdyi          = 37, &
       idBzdzi          = 38, &
       irhodusti        = 39, &
       irhodustiend     = 39 + (maxdustlarge - 1), &
       iradfxi          = irhodustiend + 1, &
       iradfyi          = irhodustiend + 2, &
       iradfzi          = irhodustiend + 3


 !--kernel related parameters
 !real, parameter    :: cnormk = 1./pi, wab0 = 1., gradh0 = -3.*wab0, radkern2 = 4F.0
 integer, parameter :: isizecellcache = 1000
 integer, parameter :: isizeneighcache = 0
 integer, parameter :: maxdensits = 100

 !--statistics which can be queried later
 integer, private         :: maxneighact,nrelink
 integer(kind=8), private :: nneightry,maxneightry,nneighact,ncalc
 integer(kind=8), private :: nptot = -1

 private

contains

!----------------------------------------------------------------
!+
!  this is the main routine for the whole code
!+
!----------------------------------------------------------------
subroutine densityiterate(icall,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                          fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level)
 use dim,       only:maxp,maxneigh,ndivcurlv,ndivcurlB,maxalpha,mhd_nonideal,nalpha,&
                     use_dust,fast_divcurlB,mpi,gr,use_apr
 use io,        only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs
 use linklist,  only:ifirstincell,ncells,get_neighbour_list,get_hmaxcell,&
                     listneigh,get_cell_location,set_hmaxcell,sync_hmax_mpi
 use part,      only:mhd,rhoh,dhdrho,rhoanddhdrho,ll,get_partinfo,iactive,&
                     hrho,iphase,igas,idust,iamgas,periodic,all_active,dustfrac
 use mpiutils,  only:reduceall_mpi,barrier_mpi,reduce_mpi,reduceall_mpi
 use mpimemory, only:reserve_stack,swap_stacks,reset_stacks,write_cell
 use mpimemory, only:stack_remote  => dens_stack_1
 use mpimemory, only:stack_waiting => dens_stack_2
 use mpimemory, only:stack_redo    => dens_stack_3
 use mpiderivs, only:send_cell,recv_cells,check_send_finished,init_cell_exchange,&
                     finish_cell_exchange,recv_while_wait,reset_cell_counters,cell_counters
 use timestep,  only:rhomaxnow
 use part,      only:ngradh
 use viscosity, only:irealvisc
 use io_summary,only:summary_variable,iosumhup,iosumhdn
 use timing,    only:increment_timer,get_timings,itimer_dens_local,itimer_dens_remote
 use omputils,  only:omp_thread_num,omp_num_threads
 integer,      intent(in)    :: icall,npart,nactive
 integer(kind=1), intent(in) :: apr_level(:)
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(in)    :: vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real(kind=4), intent(out)   :: divcurlv(:,:)
 real(kind=4), intent(out)   :: divcurlB(:,:)
 real(kind=4), intent(out)   :: alphaind(:,:)
 real(kind=4), intent(inout) :: gradh(:,:)  ! requires in for icall = 3
 real,         intent(out)   :: stressmax
 real,         intent(in)    :: rad(:,:)
 real,         intent(inout) :: radprop(:,:)
 real(kind=4), intent(out)   :: dvdx(:,:)

 real,   save :: xyzcache(isizecellcache,3)
!$omp threadprivate(xyzcache)

 integer :: i,icell
 integer :: nneigh,np,npcell
 integer :: nwarnup,nwarndown,nwarnroundoff

 logical :: getdv,realviscosity,getdB,converged
 logical :: iactivei,iamgasi,iamdusti
 integer :: iamtypei

 real    :: rhomax

 logical                   :: redo_neighbours

 integer                   :: j,k,l
 integer                   :: irequestsend(nprocs),irequestrecv(nprocs)

 type(celldens)            :: cell,xsendbuf,xrecvbuf(nprocs)
 integer                   :: mpitype

 integer                   :: n_remote_its,nlocal
 integer                   :: ncomplete_mpi
 real                      :: ntotal
 logical                   :: remote_export(nprocs),do_export,idone(nprocs),thread_complete(omp_num_threads)
 logical                   :: iterations_finished

 real(kind=4)              :: t1,t2,tcpu1,tcpu2

 if (mpi) then
    call reset_stacks
    call reset_cell_counters(cell_counters)
 endif

 if (iverbose >= 3 .and. id==master) &
    write(iprint,*) ' cell cache =',isizecellcache,' neigh cache = ',isizeneighcache,' icall = ',icall

 if (icall==0 .or. icall==1) then
    call reset_neighbour_stats(nneightry,nneighact,maxneightry,maxneighact,ncalc,nrelink)
    nwarnup       = 0
    nwarndown     = 0
    nwarnroundoff = 0
    np = 0
 endif

 !
 ! flag for whether or not we need to calculate velocity derivatives
 ! whilst doing the density iterations (needed for viscosity switches
 ! and for physical viscosity)
 !
 realviscosity = (irealvisc > 0)
 getdv = ((maxalpha==maxp .or. ndivcurlv >= 4) .and. (icall <= 1 .or. icall==3)) .or. &
         (maxdvdx==maxp .and. (use_dust .or. realviscosity .or. gr))
 if (getdv .and. ndivcurlv < 1) call fatal('densityiterate','divv not stored but it needs to be')
 getdB = (mhd .and. (ndivcurlB >= 4 .or. mhd_nonideal))

 if ( all_active ) stressmax  = 0.   ! condition is required for independent timestepping

 ! Flag for what to calculate.  To accurately calculate divcurlB, density calculation must be
 ! done first.  If calculating them sequentially, then only calculate the relevant parts for
 ! each calculation.
 ! Both logicals are .true. if fast_divcurlB = .true. set in config.F90
 if (.not. fast_divcurlB) then
    if (icall==3) then
       calculate_density  = .false.
       calculate_divcurlB = .true.
    else
       calculate_density  = .true.
       calculate_divcurlB = .false.
    endif
 endif

 ! number of cells that only have neighbours on this MPI task
 nlocal = 0

 rhomax = 0.0
!$omp parallel default(none) &
!$omp shared(icall) &
!$omp shared(ncells) &
!$omp shared(ll) &
!$omp shared(ifirstincell) &
!$omp shared(xyzh) &
!$omp shared(vxyzu) &
!$omp shared(fxyzu) &
!$omp shared(fext) &
!$omp shared(gradh) &
!$omp shared(iphase) &
!$omp shared(apr_level) &
!$omp shared(Bevol) &
!$omp shared(divcurlv) &
!$omp shared(divcurlB) &
!$omp shared(alphaind) &
!$omp shared(dustfrac) &
!$omp shared(dvdx) &
!$omp shared(id) &
!$omp shared(nprocs) &
!$omp shared(getdB) &
!$omp shared(getdv) &
!$omp shared(realviscosity) &
!$omp shared(iverbose) &
!$omp shared(iprint) &
!$omp shared(rad,radprop) &
!$omp shared(calculate_density) &
!$omp shared(stack_remote) &
!$omp shared(stack_waiting) &
!$omp shared(stack_redo) &
!$omp shared(iterations_finished) &
!$omp shared(n_remote_its) &
!$omp shared(t1) &
!$omp shared(t2) &
!$omp shared(tcpu1) &
!$omp shared(tcpu2) &
!$omp shared(cell_counters) &
!$omp shared(thread_complete) &
!$omp shared(ncomplete_mpi) &
!$omp reduction(+:nlocal) &
!$omp private(do_export) &
!$omp private(j) &
!$omp private(k) &
!$omp private(l) &
!$omp private(ntotal) &
!$omp private(remote_export) &
!$omp private(nneigh) &
!$omp private(npcell) &
!$omp private(cell) &
!$omp private(iamgasi) &
!$omp private(iamtypei) &
!$omp private(iactivei) &
!$omp private(iamdusti) &
!$omp private(converged) &
!$omp private(redo_neighbours) &
!$omp private(irequestsend) &
!$omp private(xsendbuf) &
!$omp private(xrecvbuf) &
!$omp private(irequestrecv) &
!$omp private(idone) &
!$omp private(mpitype) &
!$omp reduction(+:ncalc) &
!$omp reduction(+:np) &
!$omp reduction(max:maxneighact) &
!$omp reduction(max:maxneightry) &
!$omp reduction(+:nneighact) &
!$omp reduction(+:nneightry) &
!$omp reduction(+:nrelink) &
!$omp reduction(+:stressmax) &
!$omp reduction(max:rhomax) &
!$omp private(i)


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

    !--get the neighbour list and fill the cell cache
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
                           remote_export=remote_export)
    do_export = any(remote_export)
    cell%icell  = icell
    cell%owner  = id
    cell%nits   = 0
    cell%nneigh = 0

    call start_cell(cell,iphase,xyzh,vxyzu,fxyzu,fext,Bevol,rad,apr_level)
    call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
    call get_hmaxcell(icell,cell%hmax)

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
          call reserve_stack(stack_waiting,cell%waiting_index)  ! make a reservation on the stack
          call send_cell(cell,remote_export,irequestsend,xsendbuf,cell_counters,mpitype)  ! send the cell to remote
       endif
    endif

    call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,rad,apr_level)

    if (do_export) then
       call write_cell(stack_waiting,cell)
    else
       converged = (.not. calculate_density)
       local_its: do while (.not. converged)
          call finish_cell(cell,converged)
          call compute_hmax(cell,redo_neighbours)
          if (icall == 0) converged = .true.
          if (.not. converged) then
             if (redo_neighbours) then
                call set_hmaxcell(cell%icell,cell%hmax)
                call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
                                      cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti, &
                                      remote_export=remote_export)

                if (any(remote_export)) then
                   do_export = .true.
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
                nrelink = nrelink + 1
             endif

             call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,rad,apr_level)

             if (do_export) then
                call write_cell(stack_waiting,cell)
                exit local_its
             endif

          endif
       enddo local_its
       if (.not. do_export) then
          call store_results(icall,cell,getdv,getdB,realviscosity,stressmax,xyzh,gradh,divcurlv, &
               divcurlB,alphaind,dvdx,vxyzu,&
               dustfrac,rhomax,nneightry,nneighact,maxneightry,maxneighact,np,ncalc,radprop)
          nlocal = nlocal + 1
       endif
    endif
 enddo over_cells
 !$omp enddo

 ! if any cells were sent
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
 endif

 !$omp master
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_dens_local,t2-t1,tcpu2-tcpu1)
 call get_timings(t1,tcpu1)

 if (iverbose>=6) then
    ntotal = real(nlocal) + real(stack_waiting%n)
    if (ntotal > 0) then
       write(iprint,*) id,'domain decomposition efficiency: local cells / ncells = ',real(nlocal)/ntotal
    else
       write(iprint,*) id,'domain decomposition efficiency: local cells / ncells = 0'
    endif
 endif

 n_remote_its = 0
 iterations_finished = .false.
 if (.not.mpi) iterations_finished = .true.
 !$omp end master
 !$omp barrier

 remote_its: do while(.not. iterations_finished)

    !$omp master
    n_remote_its = n_remote_its + 1
    !$omp end master
    call reset_cell_counters(cell_counters)
    !$omp barrier

    igot_remote: if (stack_remote%n > 0) then
       !$omp do schedule(runtime)
       over_remote: do i = 1,stack_remote%n
          cell = stack_remote%cells(i)

          ! icell is unused (-1 here)
          call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
                                  cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti)

          call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,rad,apr_level)

          remote_export = .false.
          remote_export(cell%owner+1) = .true. ! use remote_export array to send back to the owner

          ! communication happened while computing contributions to remote cells
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
    call reset_cell_counters(cell_counters)
    !$omp barrier

    iam_waiting: if (mpi .and. stack_waiting%n > 0) then
       !$omp do schedule(runtime)
       over_waiting: do i = 1, stack_waiting%n
          cell = stack_waiting%cells(i)

          if (calculate_density) then
             call finish_cell(cell,converged)
             call compute_hmax(cell,redo_neighbours)
          else
             converged = .true.
          endif

          ! check for incoming cells (if converged, this may not be checked until enxt cell)
          call recv_cells(stack_remote,xrecvbuf,irequestrecv,cell_counters)

          if (.not. converged) then
             call set_hmaxcell(cell%icell,cell%hmax)
             call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
                                    cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti, &
                                    remote_export=remote_export)

             idone(:) = .false.
             do while(.not.all(idone))
                call check_send_finished(irequestsend,idone)
                call recv_cells(stack_remote,xrecvbuf,irequestrecv,cell_counters)
             enddo
             call reserve_stack(stack_redo,cell%waiting_index)
             call send_cell(cell,remote_export,irequestsend,xsendbuf,cell_counters,mpitype) ! send the cell to remote

             call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,rad,apr_level)

             call write_cell(stack_redo,cell)
          else
             call store_results(icall,cell,getdv,getdB,realviscosity,stressmax,xyzh,gradh,divcurlv, &
                  divcurlB,alphaind,dvdx,vxyzu, &
                  dustfrac,rhomax,nneightry,nneighact,maxneightry,maxneighact,np,ncalc,radprop)
          endif

       enddo over_waiting
       !$omp enddo

       !$omp master
       stack_waiting%n = 0
       !$omp end master

       idone(:) = .false.
       do while(.not.all(idone))
          call check_send_finished(irequestsend,idone)
          call recv_cells(stack_remote,xrecvbuf,irequestrecv,cell_counters)
       enddo
    endif iam_waiting

    if (mpi) call recv_while_wait(stack_remote,xrecvbuf,irequestrecv,&
             irequestsend,thread_complete,cell_counters,ncomplete_mpi)

    !$omp master
    if (reduceall_mpi('max',stack_redo%n) > 0) then
       call swap_stacks(stack_waiting, stack_redo)
    else
       iterations_finished = .true.
    endif
    stack_redo%n = 0
    !$omp end master
    !$omp barrier

 enddo remote_its

 !$omp master
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_dens_remote,t2-t1,tcpu2-tcpu1)
 !$omp end master

 if (mpi) call finish_cell_exchange(irequestrecv,xsendbuf,mpitype)

 !$omp end parallel

 if (mpi) call sync_hmax_mpi

 if (calculate_density) then
    !--reduce values
    if (realviscosity .and. maxdvdx==maxp) then
       stressmax = reduceall_mpi('max',stressmax)
    endif
    rhomax    = reduceall_mpi('max',rhomax)
    rhomaxnow = rhomax

    if (realviscosity .and. maxdvdx==maxp .and. stressmax > 0. .and. iverbose > 0 .and. id==master) then
       call warning('force','applying negative stress correction',var='max',val=-stressmax)
    endif
    !
    !--warnings
    !
    if (icall==1) then
       if (nwarnup   > 0) call summary_variable('hupdn',iosumhup,0,real(nwarnup  ))
       if (nwarndown > 0) call summary_variable('hupdn',iosumhdn,0,real(nwarndown))
       if (iverbose  >=1) call reduce_and_print_warnings(nwarnup,nwarndown,nwarnroundoff)
    endif
    !
    !--diagnostics
    !
    if (icall==0 .or. icall==1) call reduce_and_print_neighbour_stats(np)
 endif

end subroutine densityiterate

!----------------------------------------------------------------
!+
!  Internal subroutine that computes the contribution to
!  the density sums from a list of neighbours
!
!  MAKE SURE THIS ROUTINE IS INLINED BY THE COMPILER
!+
!----------------------------------------------------------------
pure subroutine get_density_sums(i,xpartveci,hi,hi1,hi21,iamtypei,iamgasi,iamdusti,apri,&
                                 listneigh,nneigh,nneighi,dxcache,xyzcache,rhosum,&
                                 ifilledcellcache,ifilledneighcache,getdv,getdB,&
                                 realviscosity,xyzh,vxyzu,Bevol,fxyzu,fext,ignoreself,rad,apr_level)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 use kernel,   only:get_kernel,get_kernel_grav1
 use part,     only:iphase,iamgas,iamdust,iamtype,maxphase,ibasetype,igas,idust,rhoh
 use part,     only:massoftype,iradxi,aprmassoftype
 use dim,      only:ndivcurlv,gravity,maxp,nalpha,use_dust,do_radiation,use_apr,maxpsph
 use options,  only:implicit_radiation
 integer,      intent(in)    :: i
 real,         intent(in)    :: xpartveci(:)
 real(kind=8), intent(in)    :: hi,hi1,hi21
 integer,      intent(in)    :: iamtypei,apri
 logical,      intent(in)    :: iamgasi,iamdusti
 integer,      intent(in)    :: listneigh(:)
 integer(kind=1), intent(in) :: apr_level(:)
 integer,      intent(in)    :: nneigh
 integer,      intent(out)   :: nneighi
 real,         intent(inout) :: dxcache(:,:)
 real,         intent(in)    :: xyzcache(:,:)
 real,         intent(out)   :: rhosum(:)
 logical,      intent(in)    :: ifilledcellcache,ifilledneighcache
 logical,      intent(in)    :: getdv,realviscosity
 logical,      intent(in)    :: getdB
 real,         intent(in)    :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real,         intent(in)    :: Bevol(:,:)
 logical,      intent(in)    :: ignoreself
 real,         intent(in)    :: rad(:,:)
 integer(kind=1)             :: iphasej
 integer                     :: iamtypej
 integer                     :: j,n,iloc
 real                        :: dx,dy,dz,runix,runiy,runiz
 real                        :: rij2,rij,rij1,q2i,qi,q2prev,rij1grkern
 real                        :: wabi,grkerni,dwdhi,dphidhi
 real                        :: projv,dvx,dvy,dvz,dax,day,daz
 real                        :: projdB,dBx,dBy,dBz,fxi,fyi,fzi,fxj,fyj,fzj
 real                        :: rhoi, rhoj,pmassi,pmassj
 logical                     :: same_type,gas_gas,iamdustj
 real                        :: dradenij

 rhosum(:) = 0.
 if (ignoreself) then
    nneighi = 1   ! self
 else
    nneighi = 0
 endif

 ! defaults for type determination
 ! these are determined from iphase if multiple phases are used
 same_type = .true.
 gas_gas   = .true.

 dphidhi   = 0.
 dx = 0. ! to avoid compiler warnings
 dy = 0.
 dz = 0.
 dvx = 0.
 dvy = 0.
 dvz = 0.
 if (nalpha > 1) then
    fxi = xpartveci(ifxi)
    fyi = xpartveci(ifyi)
    fzi = xpartveci(ifzi)
 endif

 loop_over_neigh: do n = 1,nneigh

    j = listneigh(n)
    !--do self contribution separately to avoid problems with 1/sqrt(0.)
    if ((ignoreself) .and. (j==i)) cycle loop_over_neigh
    if (j > maxpsph) cycle loop_over_neigh

    if (ifilledneighcache .and. n <= isizeneighcache) then
       rij2 = dxcache(1,n)
    else
       if (ifilledcellcache .and. n <= isizecellcache) then
          ! positions from cache are already mod boundary
          dx = xpartveci(ixi) - xyzcache(n,1)
          dy = xpartveci(iyi) - xyzcache(n,2)
          dz = xpartveci(izi) - xyzcache(n,3)
       else
          dx = xpartveci(ixi) - xyzh(1,j)
          dy = xpartveci(iyi) - xyzh(2,j)
          dz = xpartveci(izi) - xyzh(3,j)
       endif
#ifdef PERIODIC
       if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
       if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
       if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
       rij2 = dx*dx + dy*dy + dz*dz
       if (n <= isizeneighcache) dxcache(1,n) = rij2
    endif

    q2i = rij2*hi21
!
!--do interaction if r/h < compact support size
!
    if (q2i < radkern2) then
       if (ifilledneighcache .and. n <= isizeneighcache) then
          q2prev = dxcache(2,n)
          if (q2prev < radkern2) then
             rij = dxcache(3,n)
          else
             rij = sqrt(rij2)
          endif
       else
          rij = sqrt(rij2)
       endif

       qi = rij*hi1
       !--kernel and gradient
       if (gravity) then
          call get_kernel_grav1(q2i,qi,wabi,grkerni,dphidhi)
       else
          call get_kernel(q2i,qi,wabi,grkerni)
       endif

       if (n <= isizeneighcache) then
          !   could possibly ONLY store q2i if q2i>q2prev so that
          !   the maximum number of sqrts are stored
          dxcache(2,n) = q2i ! if h decreasing we don
          dxcache(3,n) = rij
          dxcache(4,n) = grkerni
          !--can ONLY fill this on first pass
          if (.not.ifilledneighcache) then
             dxcache(5,n) = dx
             dxcache(6,n) = dy
             dxcache(7,n) = dz
          endif
       endif

       !
       ! Density, gradh and div v are only computed using
       ! neighbours of the same type
       !
       if (maxphase==maxp) then
          iphasej   = iphase(j)
          iamtypej  = iamtype(iphasej)
          iamdustj  = iamdust(iphasej)
          same_type = ((iamtypei == iamtypej) .or. (ibasetype(iamtypej)==iamtypei))
          gas_gas   = (iamgasi .and. same_type)  ! this ensure that boundary particles are included in gas_gas calculations
       endif

       ! adjust masses for apr
       ! this defaults to massoftype if apr_level=1
       if (use_apr) then
          pmassi = aprmassoftype(iamtypei,apri)
          pmassj = aprmassoftype(iamtypej,apr_level(j))
       else
          pmassi = massoftype(iamtypei)
          pmassj = massoftype(iamtypej)
       endif

       sametype: if (same_type) then
          dwdhi = (-qi*grkerni - 3.*wabi)
          rhosum(irhoi)      = rhosum(irhoi) + wabi*pmassj
          rhosum(igradhi)    = rhosum(igradhi) + dwdhi*pmassj
          rhosum(igradsofti) = rhosum(igradsofti) + dphidhi*pmassj
          nneighi            = nneighi + 1
          !
          ! calculate things needed for viscosity switches
          ! and real viscosity
          !
          if ((getdv .or. getdB .or. do_radiation) .and. calculate_divcurlB) then
             rij1 = 1./(rij + epsilon(rij))
             if (ifilledneighcache .and. n <= isizeneighcache) then
                !--dx,dy,dz are either in neighbour cache or have been calculated
                dx = dxcache(5,n)
                dy = dxcache(6,n)
                dz = dxcache(7,n)
             endif
             rij1grkern = rij1*grkerni
             runix = dx*rij1grkern*pmassi
             runiy = dy*rij1grkern*pmassi
             runiz = dz*rij1grkern*pmassi

             if (getdv) then
                !--get dv and den
                dvx = xpartveci(ivxi) - vxyzu(1,j)
                dvy = xpartveci(ivyi) - vxyzu(2,j)
                dvz = xpartveci(ivzi) - vxyzu(3,j)
                projv = dvx*runix + dvy*runiy + dvz*runiz
                rhosum(idivvi) = rhosum(idivvi) + projv

                if (maxdvdx > 0 .or. ndivcurlv > 1 .or. nalpha > 1) then
                   rhosum(idvxdxi) = rhosum(idvxdxi) + dvx*runix
                   rhosum(idvxdyi) = rhosum(idvxdyi) + dvx*runiy
                   rhosum(idvxdzi) = rhosum(idvxdzi) + dvx*runiz
                   rhosum(idvydxi) = rhosum(idvydxi) + dvy*runix
                   rhosum(idvydyi) = rhosum(idvydyi) + dvy*runiy
                   rhosum(idvydzi) = rhosum(idvydzi) + dvy*runiz
                   rhosum(idvzdxi) = rhosum(idvzdxi) + dvz*runix
                   rhosum(idvzdyi) = rhosum(idvzdyi) + dvz*runiy
                   rhosum(idvzdzi) = rhosum(idvzdzi) + dvz*runiz

                   if (nalpha > 1 .and. gas_gas) then
                      !--divergence of acceleration for Cullen & Dehnen switch
                      fxj = fxyzu(1,j) + fext(1,j)
                      fyj = fxyzu(2,j) + fext(2,j)
                      fzj = fxyzu(3,j) + fext(3,j)
                      dax = fxi - fxj
                      day = fyi - fyj
                      daz = fzi - fzj

                      rhosum(idaxdxi) = rhosum(idaxdxi) + dax*runix
                      rhosum(idaxdyi) = rhosum(idaxdyi) + dax*runiy
                      rhosum(idaxdzi) = rhosum(idaxdzi) + dax*runiz
                      rhosum(idaydxi) = rhosum(idaydxi) + day*runix
                      rhosum(idaydyi) = rhosum(idaydyi) + day*runiy
                      rhosum(idaydzi) = rhosum(idaydzi) + day*runiz
                      rhosum(idazdxi) = rhosum(idazdxi) + daz*runix
                      rhosum(idazdyi) = rhosum(idazdyi) + daz*runiy
                      rhosum(idazdzi) = rhosum(idazdzi) + daz*runiz
                   endif
                   rhosum(irxxi) = rhosum(irxxi) - dx*runix
                   rhosum(irxyi) = rhosum(irxyi) - dx*runiy
                   rhosum(irxzi) = rhosum(irxzi) - dx*runiz
                   rhosum(iryyi) = rhosum(iryyi) - dy*runiy
                   rhosum(iryzi) = rhosum(iryzi) - dy*runiz
                   rhosum(irzzi) = rhosum(irzzi) - dz*runiz

                endif
             endif

             if (getdB .and. gas_gas) then
                ! we need B instead of B/rho, so used our estimated h here
                ! either it is close enough to be converged,
                ! or worst case it runs another iteration and re-calculates
                rhoi = rhoh(real(hi),  pmassi)
                rhoj = rhoh(xyzh(4,j), pmassj)
                dBx = xpartveci(iBevolxi)*rhoi - Bevol(1,j)*rhoj
                dBy = xpartveci(iBevolyi)*rhoi - Bevol(2,j)*rhoj
                dBz = xpartveci(iBevolzi)*rhoi - Bevol(3,j)*rhoj
                projdB = dBx*runix + dBy*runiy + dBz*runiz

                ! difference operator of divB
                rhosum(idivBi) = rhosum(idivBi) + projdB

                rhosum(idBxdxi) = rhosum(idBxdxi) + dBx*runix
                rhosum(idBxdyi) = rhosum(idBxdyi) + dBx*runiy
                rhosum(idBxdzi) = rhosum(idBxdzi) + dBx*runiz
                rhosum(idBydxi) = rhosum(idBydxi) + dBy*runix
                rhosum(idBydyi) = rhosum(idBydyi) + dBy*runiy
                rhosum(idBydzi) = rhosum(idBydzi) + dBy*runiz
                rhosum(idBzdxi) = rhosum(idBzdxi) + dBz*runix
                rhosum(idBzdyi) = rhosum(idBzdyi) + dBz*runiy
                rhosum(idBzdzi) = rhosum(idBzdzi) + dBz*runiz
             endif

             if (do_radiation .and. gas_gas .and. .not. implicit_radiation) then
                rhoi = rhoh(real(hi), pmassi)
                rhoj = rhoh(xyzh(4,j), pmassj)
                dradenij = rad(iradxi,j)*rhoj - xpartveci(iradxii)*rhoi
                rhosum(iradfxi) = rhosum(iradfxi) + dradenij*runix
                rhosum(iradfyi) = rhosum(iradfyi) + dradenij*runiy
                rhosum(iradfzi) = rhosum(iradfzi) + dradenij*runiz
             endif

          endif
       elseif (use_dust .and. (iamgasi  .and. iamdustj)) then
          iloc = irhodusti + iamtypej - idust
          rhosum(iloc) = rhosum(iloc) + wabi
       endif sametype

    elseif (n <= isizeneighcache) then
       ! q2prev > radkern2 from cache indicates rij has NOT been calculated for this pair
       dxcache(2,n) = q2i
       if (.not.ifilledneighcache) then
          dxcache(5,n) = dx
          dxcache(6,n) = dy
          dxcache(7,n) = dz
       endif
    endif
 enddo loop_over_neigh

end subroutine get_density_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract the matrix used in exact linear
!  interpolations from the summations calculated during
!  the density loop
!+
!----------------------------------------------------------------
pure subroutine calculate_rmatrix_from_sums(rhosum,denom,rmatrix,idone)
 real,    intent(in)  :: rhosum(:)
 real,    intent(out) :: denom
 real,    intent(out) :: rmatrix(6)
 logical, intent(out) :: idone
 real :: rxxi,rxyi,rxzi,ryyi,ryzi,rzzi

 rxxi = rhosum(irxxi)
 rxyi = rhosum(irxyi)
 rxzi = rhosum(irxzi)
 ryyi = rhosum(iryyi)
 ryzi = rhosum(iryzi)
 rzzi = rhosum(irzzi)

 denom = rxxi*ryyi*rzzi + 2.*rxyi*rxzi*ryzi &
        - rxxi*ryzi*ryzi - ryyi*rxzi*rxzi - rzzi*rxyi*rxyi

 rmatrix(1) = ryyi*rzzi - ryzi*ryzi    ! xx
 rmatrix(2) = rxzi*ryzi - rzzi*rxyi    ! xy
 rmatrix(3) = rxyi*ryzi - rxzi*ryyi    ! xz
 rmatrix(4) = rzzi*rxxi - rxzi*rxzi    ! yy
 rmatrix(5) = rxyi*rxzi - rxxi*ryzi    ! yz
 rmatrix(6) = rxxi*ryyi - rxyi*rxyi    ! zz
 idone = .true.

 return
end subroutine calculate_rmatrix_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract div and curl v from the sums
!  calculated during the density loop
!+
!----------------------------------------------------------------
pure subroutine calculate_divcurlv_from_sums(rhosum,termnorm,divcurlvi,ndivcurlv,denom,rmatrix)
 use part, only:nalpha
 integer, intent(in)  :: ndivcurlv
 real,    intent(in)  :: rhosum(:),denom,rmatrix(6)
 real,    intent(in)  :: termnorm
 real,    intent(out) :: divcurlvi(5)
 real :: div_a
 real :: gradaxdx,gradaxdy,gradaxdz,gradaydx,gradaydy,gradaydz,gradazdx,gradazdy,gradazdz
 real :: ddenom,gradvxdxi,gradvxdyi,gradvxdzi
 real :: gradvydxi,gradvydyi,gradvydzi,gradvzdxi,gradvzdyi,gradvzdzi
 real :: dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi
 logical, parameter :: use_exact_linear = .true.

 !--divergence of the velocity field
 if (ndivcurlv >= 1) divcurlvi(1) = -rhosum(idivvi)*termnorm

 !--curl of the velocity field
 if (ndivcurlv >= 4) then
    divcurlvi(2) = -(rhosum(idvzdyi) - rhosum(idvydzi))*termnorm
    divcurlvi(3) = -(rhosum(idvxdzi) - rhosum(idvzdxi))*termnorm
    divcurlvi(4) = -(rhosum(idvydxi) - rhosum(idvxdyi))*termnorm
 endif

 !--time derivative of div v, needed for Cullen-Dehnen switch
 if (nalpha >= 2) then
    !--Divvdt For switch
    if (use_exact_linear .and. abs(denom) > tiny(denom)) then
       ddenom = 1./denom
       call exactlinear(gradaxdx,gradaxdy,gradaxdz,rhosum(idaxdxi),rhosum(idaxdyi),rhosum(idaxdzi),rmatrix,ddenom)
       call exactlinear(gradaydx,gradaydy,gradaydz,rhosum(idaydxi),rhosum(idaydyi),rhosum(idaydzi),rmatrix,ddenom)
       call exactlinear(gradazdx,gradazdy,gradazdz,rhosum(idazdxi),rhosum(idazdyi),rhosum(idazdzi),rmatrix,ddenom)
       div_a = -(gradaxdx + gradaydy + gradazdz)

       call exactlinear(gradvxdxi,gradvxdyi,gradvxdzi, &
                         rhosum(idvxdxi),rhosum(idvxdyi),rhosum(idvxdzi),rmatrix,ddenom)
       call exactlinear(gradvydxi,gradvydyi,gradvydzi, &
                         rhosum(idvydxi),rhosum(idvydyi),rhosum(idvydzi),rmatrix,ddenom)
       call exactlinear(gradvzdxi,gradvzdyi,gradvzdzi, &
                         rhosum(idvzdxi),rhosum(idvzdyi),rhosum(idvzdzi),rmatrix,ddenom)

       dvxdxi = -gradvxdxi
       dvxdyi = -gradvxdyi
       dvxdzi = -gradvxdzi
       dvydxi = -gradvydxi
       dvydyi = -gradvydyi
       dvydzi = -gradvydzi
       dvzdxi = -gradvzdxi
       dvzdyi = -gradvzdyi
       dvzdzi = -gradvzdzi
    else
       div_a = -termnorm*(rhosum(idaxdxi) + rhosum(idaydyi) + rhosum(idazdzi))
       dvxdxi = -termnorm*rhosum(idvxdxi)
       dvxdyi = -termnorm*rhosum(idvxdyi)
       dvxdzi = -termnorm*rhosum(idvxdzi)
       dvydxi = -termnorm*rhosum(idvydxi)
       dvydyi = -termnorm*rhosum(idvydyi)
       dvydzi = -termnorm*rhosum(idvydzi)
       dvzdxi = -termnorm*rhosum(idvzdxi)
       dvzdyi = -termnorm*rhosum(idvzdyi)
       dvzdzi = -termnorm*rhosum(idvzdzi)
    endif
    divcurlvi(5) = div_a - (dvxdxi**2 + dvydyi**2 + dvzdzi**2 + &
                             2.*(dvxdyi*dvydxi + dvxdzi*dvzdxi + dvydzi*dvzdyi))
 endif

end subroutine calculate_divcurlv_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract div, curl and grad B from the sums
!  calculated during the density loop
!+
!----------------------------------------------------------------
pure subroutine calculate_divcurlB_from_sums(rhosum,termnorm,divcurlBi,ndivcurlB)
 integer, intent(in)  :: ndivcurlB
 real,    intent(in)  :: rhosum(:)
 real,    intent(in)  :: termnorm
 real,    intent(out) :: divcurlBi(ndivcurlB)

 ! we need these for adaptive resistivity switch
 if (ndivcurlB >= 1) divcurlBi(1) = -rhosum(idivBi)*termnorm
 if (ndivcurlB >= 4) then
    divcurlBi(2) = -(rhosum(idBzdyi) - rhosum(idBydzi))*termnorm
    divcurlBi(3) = -(rhosum(idBxdzi) - rhosum(idBzdxi))*termnorm
    divcurlBi(4) = -(rhosum(idBydxi) - rhosum(idBxdyi))*termnorm
 endif

end subroutine calculate_divcurlB_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract velocity gradients from summations
!  calculated during the density loop.
!+
!----------------------------------------------------------------
subroutine calculate_strain_from_sums(rhosum,termnorm,denom,rmatrix,dvdx)
 real, intent(in)  :: rhosum(:)
 real, intent(in)  :: termnorm,denom
 real, intent(in)  :: rmatrix(6)
 real, intent(out) :: dvdx(9)

 real :: ddenom,gradvxdxi,gradvxdyi,gradvxdzi
 real :: gradvydxi,gradvydyi,gradvydzi,gradvzdxi,gradvzdyi,gradvzdzi
 real :: dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi

! if (abs(denom) > tiny(denom)) then ! do exact linear first derivatives
 if (.false.) then ! do exact linear first derivatives
    ddenom = 1./denom
    call exactlinear(gradvxdxi,gradvxdyi,gradvxdzi, &
                     rhosum(idvxdxi),rhosum(idvxdyi),rhosum(idvxdzi),rmatrix,ddenom)
    call exactlinear(gradvydxi,gradvydyi,gradvydzi, &
                     rhosum(idvydxi),rhosum(idvydyi),rhosum(idvydzi),rmatrix,ddenom)
    call exactlinear(gradvzdxi,gradvzdyi,gradvzdzi, &
                     rhosum(idvzdxi),rhosum(idvzdyi),rhosum(idvzdzi),rmatrix,ddenom)

    !print*,'dvxdxi = ',-rhosum(idvxdxi)*termnorm,gradvxdxi
    dvxdxi = -gradvxdxi
    dvxdyi = -gradvxdyi
    dvxdzi = -gradvxdzi
    dvydxi = -gradvydxi
    dvydyi = -gradvydyi
    dvydzi = -gradvydzi
    dvzdxi = -gradvzdxi
    dvzdyi = -gradvzdyi
    dvzdzi = -gradvzdzi
 else

    !--these make rho*dv/dx_i
    dvxdxi = -rhosum(idvxdxi)*termnorm
    dvxdyi = -rhosum(idvxdyi)*termnorm
    dvxdzi = -rhosum(idvxdzi)*termnorm
    dvydxi = -rhosum(idvydxi)*termnorm
    dvydyi = -rhosum(idvydyi)*termnorm
    dvydzi = -rhosum(idvydzi)*termnorm
    dvzdxi = -rhosum(idvzdxi)*termnorm
    dvzdyi = -rhosum(idvzdyi)*termnorm
    dvzdzi = -rhosum(idvzdzi)*termnorm
 endif

 dvdx(:) = (/dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi/)

end subroutine calculate_strain_from_sums

!----------------------------------------------------------------
!+
!  Internal subroutine to get maximum of stress tensor
!  (to avoid tensile instability)
!+
!----------------------------------------------------------------
pure subroutine get_max_stress(dvdx,divvi,rho1i,stressmax,shearvisc,bulkvisc)
 use part, only:strain_from_dvdx
 real, intent(in)    :: dvdx(9), divvi, rho1i, shearvisc, bulkvisc
 real, intent(inout) :: stressmax
 real :: strainmax,stressiso,strain(6)

 strain = strain_from_dvdx(dvdx)

 ! shearvisc = eta/rho, so this is eta/rho**2
 strainmax = -shearvisc*rho1i*maxval(strain) ! 1/rho*L^2/T*1/L*L/T = 1/rho*L^2/T^2
 stressiso = (2./3.*shearvisc - bulkvisc)*divvi*rho1i ! NB: divv -ve at high Mach no.

 ! we initialise stressmax to zero, as for the purpose of preventing the
 ! tensile instability we only care if the total stress is negative
 ! if stress tensor is positive, don't need correction (stressmax=0)
 stressmax = max(stressmax,-(stressiso + strainmax))
 stressmax = 0.

end subroutine get_max_stress

!----------------------------------------------------------------
!+
!  Internal subroutine that inverts the matrix to get an
!  exact linear derivative
!+
!----------------------------------------------------------------
pure subroutine exactlinear(gradAx,gradAy,gradAz,dAx,dAy,dAz,rmatrix,ddenom)
 real, intent(out) :: gradAx,gradAy,gradAz
 real, intent(in)  :: dAx,dAy,dAz
 real, intent(in)  :: rmatrix(6)
 real, intent(in)  :: ddenom
 !
 !--we return the gradient as the following matrix inversion:
 !  gradAx =(dAx*termxx + dAy*termxy + dAz*termxz)*ddenom
 !  gradAy =(dAx*termxy + dAy*termyy + dAz*termyz)*ddenom
 !  gradAz =(dAx*termxz + dAy*termyz + dAz*termzz)*ddenom
 !
 gradAx =(dAx*rmatrix(1) + dAy*rmatrix(2) + dAz*rmatrix(3))*ddenom
 gradAy =(dAx*rmatrix(2) + dAy*rmatrix(4) + dAz*rmatrix(5))*ddenom
 gradAz =(dAx*rmatrix(3) + dAy*rmatrix(5) + dAz*rmatrix(6))*ddenom

end subroutine exactlinear

!----------------------------------------------------------------
!+
!  subroutine to reduce and print warnings across processors
!  related to h-rho iterations
!+fxyzu
!----------------------------------------------------------------
subroutine reduce_and_print_warnings(nwarnup,nwarndown,nwarnroundoff)
 use mpiutils, only:reduce_mpi
 use io,       only:id,master,iprint
 integer, intent(inout) :: nwarnup,nwarndown,nwarnroundoff

 nwarnup       = int(reduce_mpi('+',nwarnup))
 nwarndown     = int(reduce_mpi('+',nwarndown))
 nwarnroundoff = int(reduce_mpi('+',nwarnroundoff))

#ifndef NOWARNRESTRICTEDHJUMP
 if (id==master .and. nwarnup > 0) then
    write(iprint,*) ' WARNING: restricted h jump (up) ',nwarnup,' times'
 endif
 if (id==master .and. nwarndown > 0) then
    write(iprint,*) ' WARNING: restricted h jump (down) ',nwarndown,' times'
 endif
#endif
 if (id==master .and. nwarnroundoff > 0) then
    write(iprint,*) ' WARNING: denom in exact linear gradients zero on ',nwarnroundoff,' particles'
 endif

end subroutine reduce_and_print_warnings

!----------------------------------------------------------------
!+
!  query function to return neighbour statistics
!  (must be called *after* density evaluation
!   and will be correct ONLY on the master thread)
!+
!----------------------------------------------------------------
subroutine get_neighbour_stats(trialmean,actualmean,maxtrial,maxactual,nrhocalc,nactualtot)
 real,            intent(out) :: trialmean,actualmean
 integer,         intent(out) :: maxtrial,maxactual
 integer(kind=8), intent(out) :: nrhocalc,nactualtot

 if (nptot > 0) then
    trialmean    = nneightry/real(nptot)
    actualmean   = nneighact/real(nptot)
    maxtrial     = int(maxneightry)
    maxactual    = maxneighact
    nrhocalc     = ncalc
    nactualtot   = nneighact
 else ! densityforce has not been called
    trialmean = -1; actualmean   = -1
    maxtrial  = -1; maxactual    = -1
    nrhocalc  = -1; nactualtot   = -1
 endif

end subroutine get_neighbour_stats

subroutine reset_neighbour_stats(nneightry,nneighact,maxneightry,maxneighact,ncalc,nrelink)
 integer,         intent(out) :: maxneighact,nrelink
 integer(kind=8), intent(out) :: ncalc,nneightry,nneighact,maxneightry

 nneightry = 0
 nneighact = 0
 maxneightry = 0
 maxneighact = 0
 ncalc = 0_8
 nneighact = 0
 nrelink = 0_8

end subroutine reset_neighbour_stats

!----------------------------------------------------------------
!+
!  function to collate neighbour-finding statistics across
!  processors and print the results
!+
!----------------------------------------------------------------
subroutine reduce_and_print_neighbour_stats(np)
 use mpiutils, only:reduce_mpi
 use io,       only:iprint,id,master,iverbose
 integer, intent(in) :: np

 nptot       = reduce_mpi('+',np)
 nneightry   = reduce_mpi('+',nneightry)
 nneighact   = reduce_mpi('+',nneighact)
 maxneightry = reduce_mpi('max',maxneightry)
 maxneighact = int(reduce_mpi('max',maxneighact))
 nrelink     = int(reduce_mpi('+',nrelink))
 ncalc       = reduce_mpi('+',ncalc)

 if (id==master .and. iverbose >= 2 .and. nptot > 0 .and. nneighact > 0) then
    write(iprint,"(1x,a,f11.2,2(a,f7.2))") 'trial neigh mean  :',nneightry/real(nptot), &
                 ', real neigh mean = ',nneighact/real(nptot), &
                 ' ratio try/act= ',nneightry/real(nneighact)
    write(iprint,"(1x,a,i11,a,i8)")   'trial neigh max   :',maxneightry,', max real neigh = ',maxneighact
    write(iprint,"(1x,a,i11,a,f7.3)") 'n neighbour calls :',nrelink, ', mean per part   = ',nrelink/real(nptot) + 1
    write(iprint,"(1x,a,i11,a,f7.3)") 'n density calcs   :',ncalc,', mean per part   = ',ncalc/real(nptot)
 endif

end subroutine reduce_and_print_neighbour_stats
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
pure subroutine compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext, &
                             xyzcache,rad,apr_level)
 use part,        only:get_partinfo,iamgas,igas,maxphase
 use viscosity,   only:irealvisc
 use io,          only:id
 use dim,         only:mpi,use_apr

 type(celldens),  intent(inout)  :: cell

 integer,         intent(in)     :: listneigh(:)
 integer,         intent(in)     :: nneigh
 logical,         intent(in)     :: getdv
 logical,         intent(in)     :: getdB
 real,            intent(in)     :: Bevol(:,:)
 real,            intent(in)     :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real,            intent(in)     :: xyzcache(isizecellcache,3)
 real,            intent(in)     :: rad(:,:)
 integer(kind=1), intent(in)     :: apr_level(:)

 real                            :: dxcache(7,isizeneighcache)

 real(kind=8)                    :: hi
 real(kind=8)                    :: hi1,hi21,hi31,hi41

 integer                         :: iamtypei
 logical                         :: iactivei,iamgasi,iamdusti

 logical                         :: realviscosity
 logical                         :: ignoreself
 integer                         :: nneighi,apri
 integer                         :: i,lli

 realviscosity = (irealvisc > 0)

 over_parts: do i = 1,cell%npcell
    lli = inodeparts(cell%arr_index(i))
    ! note: only active particles have been sent here
    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamgasi  = .true.
    endif

    hi    = cell%h(i)
    hi1   = 1./hi
    hi21  = hi1*hi1
    hi31  = hi1*hi21
    hi41  = hi21*hi21

    if (use_apr) then
       apri = cell%apr(i)
    else
       apri = 1
    endif


    ignoreself = (cell%owner == id)

    call get_density_sums(lli,cell%xpartvec(:,i),hi,hi1,hi21,iamtypei,iamgasi,iamdusti,&
                          apri,listneigh,nneigh,nneighi,dxcache,xyzcache,&
                          cell%rhosums(:,i),.true.,.false.,getdv,getdB,realviscosity,&
                          xyzh,vxyzu,Bevol,fxyzu,fext,ignoreself,rad,apr_level)

    cell%nneightry = nneigh
    cell%nneigh(i) = nneighi
 enddo over_parts

end subroutine compute_cell
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
pure subroutine compute_hmax(cell,redo_neighbours)
 use kernel, only:radkern
 type(celldens), intent(inout) :: cell
 logical,         intent(out)  :: redo_neighbours
 real                          :: hmax_old,hmax

 redo_neighbours = .false.
 if (cell%npcell > 0) then
    hmax_old = cell%hmax
    hmax     = 1.01*maxval(cell%h(1:cell%npcell))
    if (hmax > hmax_old) redo_neighbours = .true.
    cell%hmax  = hmax
    cell%rcuti = radkern*hmax
 endif
end subroutine compute_hmax
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine start_cell(cell,iphase,xyzh,vxyzu,fxyzu,fext,Bevol,rad,apr_level)
 use io,          only:fatal
 use dim,         only:maxp,maxvxyzu,do_radiation,use_apr,maxpsph
 use part,        only:maxphase,get_partinfo,mhd,igas,iamgas,&
                       iamboundary,ibasetype,iradxi

 type(celldens),     intent(inout) :: cell
 integer(kind=1),    intent(in)    :: iphase(:)
 real,               intent(in)    :: xyzh(:,:)
 real,               intent(in)    :: vxyzu(:,:)
 real,               intent(in)    :: fxyzu(:,:)
 real,               intent(in)    :: fext(:,:)
 real,               intent(in)    :: Bevol(:,:)
 real,               intent(in)    :: rad(:,:)
 integer(kind=1),            intent(in)    :: apr_level(:)

 integer :: i,ip
 integer :: iamtypei
 logical :: iactivei,iamgasi,iamdusti

 cell%npcell = 0
 over_parts: do ip = inoderange(1,cell%icell),inoderange(2,cell%icell)
    i = inodeparts(ip)

    if (i < 0 .or. i > maxpsph) then
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
    if (.not.iactivei) then ! skip boundary particles + inactive particles
       cycle over_parts
    endif

    cell%npcell = cell%npcell + 1

    cell%arr_index(cell%npcell)               = ip
    cell%iphase(cell%npcell)                  = iphase(i)

    cell%xpartvec(ixi,cell%npcell)            = xyzh(1,i)
    cell%xpartvec(iyi,cell%npcell)            = xyzh(2,i)
    cell%xpartvec(izi,cell%npcell)            = xyzh(3,i)

    cell%h(cell%npcell)                       = xyzh(4,i)
    cell%h_old(cell%npcell)                   = xyzh(4,i)

    cell%xpartvec(ivxi,cell%npcell)           = vxyzu(1,i)
    cell%xpartvec(ivyi,cell%npcell)           = vxyzu(2,i)
    cell%xpartvec(ivzi,cell%npcell)           = vxyzu(3,i)

    if (maxvxyzu >= 4) then
       cell%xpartvec(ieni,cell%npcell)        = vxyzu(4,i)
    endif

    cell%xpartvec(ifxi,cell%npcell)           = fxyzu(1,i) + fext(1,i)
    cell%xpartvec(ifyi,cell%npcell)           = fxyzu(2,i) + fext(2,i)
    cell%xpartvec(ifzi,cell%npcell)           = fxyzu(3,i) + fext(3,i)

    if (mhd) then
       if (iamgasi) then
          cell%xpartvec(iBevolxi,cell%npcell) = Bevol(1,i)
          cell%xpartvec(iBevolyi,cell%npcell) = Bevol(2,i)
          cell%xpartvec(iBevolzi,cell%npcell) = Bevol(3,i)
          cell%xpartvec(ipsi,cell%npcell)     = Bevol(4,i)
       else
          cell%xpartvec(iBevolxi:ipsi,cell%npcell)   = 0. ! to avoid compiler warning
       endif
    endif

    if (do_radiation) cell%xpartvec(iradxii,cell%npcell) = rad(iradxi,i)

    if (use_apr) then
       cell%apr(cell%npcell)                     = apr_level(i)
    else
       cell%apr(cell%npcell)                     = 1
    endif

 enddo over_parts

end subroutine start_cell
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine finish_cell(cell,cell_converged)
 use dim,      only:use_apr
 use io,       only:iprint,fatal
 use part,     only:get_partinfo,iamgas,maxphase,massoftype,igas,hrho,aprmassoftype
 use options,  only:tolh

 type(celldens),  intent(inout) :: cell
 logical,         intent(out)   :: cell_converged
 real                           :: rhosum(maxrhosum)
 real                           :: dhdrhoi,rhohi,omegai
 real                           :: rhoi
 real(kind=8)                   :: gradhi
 real                           :: func,dfdh1,hi,hi_old,hnew
 real                           :: pmassi, xyzh(4)
 integer                        :: i,iamtypei,apri !,nwarnup,nwarndown
 logical                        :: iactivei,iamgasi,iamdusti,converged

 cell%nits = cell%nits + 1
 cell_converged = .true.
 over_parts: do i = 1,cell%npcell
    hi = cell%h(i)
    hi_old = cell%h_old(i)
    rhosum = cell%rhosums(:,i)

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamgasi  = .true.
    endif
    !if (.not.iactivei) print*,' ERROR: should be no inactive particles here',iamtypei,iactivei

    apri = cell%apr(i)
    if (use_apr) then
       pmassi = aprmassoftype(iamtypei,apri)
    else
       pmassi = massoftype(iamtypei)
    endif

    call finish_rhosum(rhosum,pmassi,hi,.true.,rhoi=rhoi,rhohi=rhohi,&
                       gradhi=gradhi,dhdrhoi_out=dhdrhoi,omegai_out=omegai)

    func = rhohi - rhoi
    if (omegai > tiny(omegai)) then
       dfdh1 = dhdrhoi/omegai
    else
       dfdh1 = dhdrhoi/abs(omegai + epsilon(omegai))
    endif
    hnew = hi - func*dfdh1
    if (hnew > 1.2*hi) then
       ! nwarnup   = nwarnup + 1
       hnew      = 1.2*hi
    elseif (hnew < 0.8*hi) then
       ! nwarndown = nwarndown + 1
       hnew      = 0.8*hi
    endif

    converged = ((abs(hnew-hi)/hi_old) < tolh .and. omegai > 0. .and. hi > 0.)
    if (cell_converged) cell_converged = converged

    if ((.not. converged) .and. (cell%nits >= maxdensits)) then
       xyzh(1) = cell%xpartvec(ixi,i)
       xyzh(2) = cell%xpartvec(iyi,i)
       xyzh(3) = cell%xpartvec(izi,i)
       write(iprint,*) 'ERROR: density iteration failed after ',cell%nits,' iterations'
       write(iprint,*) 'hnew = ',hnew,' hi_old = ',hi_old,' nneighi = ',cell%nneigh(i)
       write(iprint,*) 'rhoi = ',rhoi,' gradhi = ',gradhi
       write(iprint,*) 'error = ',abs(hnew-hi)/hi_old,' tolh = ',tolh
       write(iprint,*) 'itype = ',iamtypei
       write(iprint,*) 'x,y,z = ',xyzh(1:3)
       write(iprint,*) 'vx,vy,vz = ',cell%xpartvec(ivxi:ivzi,i)
       call fatal('densityiterate','could not converge in density',inodeparts(cell%arr_index(i)),'error',abs(hnew-hi)/hi_old)
    endif

    if (converged) then
       cell%h(i) = hi
    else
       cell%h(i) = hnew
    endif

 enddo over_parts

end subroutine finish_cell
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
pure subroutine finish_rhosum(rhosum,pmassi,hi,iterating,rhoi,rhohi,gradhi,gradsofti,dhdrhoi_out,omegai_out)
 use part,  only:rhoh,dhdrho
 real,          intent(in)              :: rhosum(maxrhosum)
 real,          intent(in)              :: pmassi
 real,          intent(in)              :: hi
 logical,       intent(in)              :: iterating !false for the last bit where we are computing the final result
 real,          intent(out)             :: rhoi
 real(kind=8),  intent(out)             :: gradhi
 real,          intent(out),  optional  :: rhohi
 real(kind=8),  intent(out),  optional  :: gradsofti
 real,          intent(out),  optional  :: dhdrhoi_out
 real,          intent(out),  optional  :: omegai_out

 real           :: omegai,dhdrhoi
 real(kind=8)   :: hi1,hi21,hi31,hi41

 hi1   = 1./hi
 hi21  = hi1*hi1
 hi31  = hi1*hi21
 hi41  = hi21*hi21

 rhoi   = cnormk*(rhosum(irhoi) + wab0*pmassi)*hi31
 gradhi = cnormk*(rhosum(igradhi) + gradh0*pmassi)*hi41

 dhdrhoi = dhdrho(hi,pmassi)
 omegai = 1. - dhdrhoi*gradhi
 gradhi = 1./omegai

 if (iterating) then
    rhohi = rhoh(hi,pmassi)
    dhdrhoi_out = dhdrhoi
    omegai_out = omegai
 else
    gradsofti = (rhosum(igradsofti) + dphidh0*pmassi)*hi21 ! NB: no cnormk in gradsoft
    gradsofti = gradsofti*dhdrhoi
 endif

end subroutine finish_rhosum
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine store_results(icall,cell,getdv,getdb,realviscosity,stressmax,xyzh,&
                         gradh,divcurlv,divcurlB,alphaind,dvdx,vxyzu,&
                         dustfrac,rhomax,nneightry,nneighact,maxneightry,&
                         maxneighact,np,ncalc,radprop)
 use part,        only:hrho,rhoh,get_partinfo,iamgas,&
                       mhd,maxphase,massoftype,igas,ndustlarge,ndustsmall,xyzh_soa,&
                       maxgradh,idust,ifluxx,ifluxz,ithick,aprmassoftype
 use io,          only:fatal,real4
 use dim,         only:maxp,ndivcurlv,ndivcurlB,nalpha,use_dust,do_radiation,use_apr
 use options,     only:use_dustfrac,implicit_radiation
 use viscosity,   only:bulkvisc,shearparam
 use linklist,    only:set_hmaxcell
 use kernel,      only:radkern
 use kdtree,      only:inodeparts

 integer,         intent(in)    :: icall
 type(celldens),  intent(in)    :: cell
 logical,         intent(in)    :: getdv
 logical,         intent(in)    :: getdB
 logical,         intent(in)    :: realviscosity
 real,            intent(inout) :: stressmax
 real,            intent(inout) :: xyzh(:,:)
 real(kind=4),    intent(inout) :: gradh(:,:)
 real(kind=4),    intent(inout) :: divcurlv(:,:)
 real(kind=4),    intent(inout) :: divcurlB(:,:)
 real(kind=4),    intent(inout) :: alphaind(:,:)
 real(kind=4),    intent(inout) :: dvdx(:,:)
 real,            intent(in)    :: vxyzu(:,:)
 real,            intent(out)   :: dustfrac(:,:)
 real,            intent(inout) :: rhomax
 integer(kind=8), intent(inout) :: nneightry
 integer(kind=8), intent(inout) :: nneighact
 integer(kind=8), intent(inout) :: maxneightry
 integer,         intent(inout) :: maxneighact
 integer,         intent(inout) :: np
 integer(kind=8), intent(inout) :: ncalc
 real,            intent(inout) :: radprop(:,:)

 real         :: rhosum(maxrhosum)

 integer      :: iamtypei,i,lli,l,apri
 logical      :: iactivei,iamgasi,iamdusti
 logical      :: igotrmatrix
 real         :: hi,hi1,hi21,hi31,hi41
 real         :: pmassi,rhoi
 real(kind=8) :: gradhi,gradsofti
 real         :: divcurlvi(5),rmatrix(6),dvdxi(9)
 real         :: divcurlBi(ndivcurlB)
 real         :: rho1i,term,denom,rhodusti(maxdustlarge)

 do i = 1,cell%npcell
    lli = inodeparts(cell%arr_index(i))
    hi = cell%h(i)
    rhosum = cell%rhosums(:,i)

    if (hi < 0.) call fatal('densityiterate','hi < 0 after iterations',lli,var='h',val=hi)

    hi1   = 1./hi
    hi21  = hi1*hi1
    hi31  = hi1*hi21
    hi41  = hi21*hi21

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamgasi  = .true.
       iamdusti = .false.
    endif

    apri = cell%apr(i)
    if (use_apr) then
       pmassi = aprmassoftype(iamtypei,apri)
    else
       pmassi = massoftype(iamtypei)
    endif

    if (calculate_density) then
       call finish_rhosum(rhosum,pmassi,hi,.false.,rhoi=rhoi,gradhi=gradhi,gradsofti=gradsofti)

       !
       !--store final results of density iteration
       !
       xyzh(4,lli) = hrho(rhoi,pmassi)
       xyzh_soa(cell%arr_index(i),4) = xyzh(4,lli)

       if (xyzh(4,lli) < 0.) call fatal('densityiterate','setting negative h from hrho',i,var='rhoi',val=real(rhoi))

       if (maxgradh==maxp) then
          gradh(1,lli) = real(gradhi,kind=kind(gradh))
#ifdef GRAVITY
          gradh(2,lli) = real(gradsofti,kind=kind(gradh))
#endif
       endif
       rhomax = max(rhomax,real(rhoi))
    else
       rhoi = rhoh(hi,pmassi)
    endif

    if (calculate_divcurlB) then
       gradhi = gradh(1,lli)
       rho1i  = 1./rhoi
       if (use_dust .and. .not. use_dustfrac) then
          !
          ! for 2-fluid dust compute dust density on gas particles
          ! and store it in dustfrac as dust-to-gas ratio
          ! so that rho times dustfrac gives dust density
          !
          dustfrac(:,lli) = 0.
          if (iamgasi) then
             do l=1,ndustlarge
                rhodusti(l) = cnormk*massoftype(idust+l-1)*(rhosum(irhodusti+l-1))*hi31 !TDB fix apr here
                dustfrac(ndustsmall+l,lli) = rhodusti(l)*rho1i ! dust-to-gas ratio
             enddo
          endif
       endif
       !
       ! store divv and curl v and related quantities
       !
       igotrmatrix = .false.

       term = cnormk*gradhi*rho1i*hi41
       if (getdv) then
          call calculate_rmatrix_from_sums(rhosum,denom,rmatrix,igotrmatrix)
          call calculate_divcurlv_from_sums(rhosum,term,divcurlvi,ndivcurlv,denom,rmatrix)
          divcurlv(1:ndivcurlv,lli) = real(divcurlvi(1:ndivcurlv),kind=kind(divcurlv)) ! save to global memory
          if (nalpha >= 3) alphaind(3,lli) = real4(divcurlvi(5))
       else ! we always need div v for h prediction
          if (ndivcurlv >= 1) divcurlv(1,lli) = -real4(rhosum(idivvi)*term)
          if (nalpha >= 2) alphaind(2,lli) = 0.
       endif
       !
       ! store div B, curl B and related quantities
       !
       if (mhd .and. iamgasi) then
          if (getdB) then
             call calculate_divcurlB_from_sums(rhosum,term,divcurlBi,ndivcurlB)
             divcurlB(:,lli) = real(divcurlBi(:),kind=kind(divcurlB))
          else
             divcurlBi(:) = 0.
          endif
       endif
       !
       !--get strain tensor from summations
       !
       if (maxdvdx==maxp .and. getdv) then
          if (.not.igotrmatrix) call calculate_rmatrix_from_sums(cell%rhosums(:,i),denom,rmatrix,igotrmatrix)
          call calculate_strain_from_sums(cell%rhosums(:,i),term,denom,rmatrix,dvdxi)
          ! check for negative stresses to prevent tensile instability
          if (realviscosity) call get_max_stress(dvdxi,divcurlvi(1),rho1i,stressmax,shearparam,bulkvisc)
          ! store strain tensor
          dvdx(:,lli) = real(dvdxi(:),kind=kind(dvdx))
       endif

       if (do_radiation .and. iamgasi .and. .not. implicit_radiation) then
          radprop(ifluxx:ifluxz,lli) = cell%rhosums(iradfxi:iradfzi,i)*term
       endif
    endif

    if (calculate_density) then
       ! stats
       nneightry = nneightry + cell%nneightry
       nneighact = nneighact + cell%nneigh(i)
       maxneightry = max(int(maxneightry),cell%nneightry)
       maxneighact = max(maxneighact,cell%nneigh(i))
    endif
 enddo
 np = np + cell%npcell
 ncalc = ncalc + cell%npcell * cell%nits

end subroutine store_results

end module densityforce
