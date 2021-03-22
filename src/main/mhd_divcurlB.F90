!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module mhd_divcurlB
!
! This module calculates divB & curlB.  Since the evolved variable is B/rho,
!  this needs to be performed after the density is updated to prevent race
!  conditions.  For simplicity, we have copied dens.F90, and stripped out
!  everything that is not related to the calculation of divcurlB.
!  To calculate divcurlB in dens (permitting the race condition), set
!  mhd_racecondition = .true. in config.F90.
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, fastmath, io, io_summary, kdtree, kernel,
!   linklist, mpidens, mpiderivs, mpiutils, options, part, stack, timestep,
!   timing, viscosity
!
 use dim,     only:maxp,maxrhosum
 use kdtree,  only:inodeparts,inoderange
 use kernel,  only:cnormk,radkern2
 use mpidens, only:celldens,stackdens
 use timing,  only:getused,printused,print_time

 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 public :: calculate_divcurlB

 !--indexing for xpartveci array
 !  using the same array as dens.F90, but re-indexing to use only the first 6 entries
 integer, parameter :: &
       ixi      = 1, &
       iyi      = 2, &
       izi      = 3, &
       iBevolxi = 4, &
       iBevolyi = 5, &
       iBevolzi = 6

 !--indexing for rhosum array
 !  using the same array as dens.F90, but re-indexing to use only the first 10 entries
 integer, parameter :: &
       idivBi           =  1, &
       idBxdxi          =  2, &
       idBxdyi          =  3, &
       idBxdzi          =  4, &
       idBydxi          =  5, &
       idBydyi          =  6, &
       idBydzi          =  7, &
       idBzdxi          =  8, &
       idBzdyi          =  9, &
       idBzdzi          = 10

 integer, parameter :: isizecellcache = 50000

 private

contains

!----------------------------------------------------------------
!+
!  this is the main routine to calculate divcurlB
!+
!----------------------------------------------------------------
subroutine calculate_divcurlB(icall,npart,nactive,xyzh,divcurlB,Bevol,gradh)
 use dim,       only:maxneigh,ndivcurlB
 use io,        only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs
 use linklist,  only:ifirstincell,ncells,get_neighbour_list,get_hmaxcell,&
                     listneigh,get_cell_location,set_hmaxcell,sync_hmax_mpi
 use part,      only:rhoh,dhdrho,rhoanddhdrho,ll,get_partinfo,iactive,&
                     hrho,iphase,igas,iamgas,periodic
 use mpiutils,  only:reduceall_mpi,barrier_mpi,reduce_mpi,reduceall_mpi
#ifdef MPI
 use stack,     only:reserve_stack,swap_stacks
 use stack,     only:stack_remote => dens_stack_1
 use stack,     only:stack_waiting => dens_stack_2
 use stack,     only:stack_redo => dens_stack_3
 use mpiderivs, only:send_cell,recv_cells,check_send_finished,init_cell_exchange, &
                     finish_cell_exchange,recv_while_wait,reset_cell_counters
#endif
 integer,      intent(in)    :: icall,npart,nactive
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real(kind=4), intent(out)   :: divcurlB(:,:)
 real(kind=4), intent(in)    :: gradh(:,:)

 real,   save :: xyzcache(isizecellcache,3)
!$omp threadprivate(xyzcache)

 integer :: i,icell
 integer :: nneigh
 logical :: iactivei,iamgasi
 integer :: iamtypei
 logical :: remote_export(nprocs)
 integer                   :: npcell
 type(celldens)            :: cell
 logical                   :: redo_neighbours
#ifdef MPI
 integer                   :: j,k,l
 integer                   :: irequestsend(nprocs),irequestrecv(nprocs)
 type(celldens)            :: xrecvbuf(nprocs),xsendbuf
 integer                   :: mpiits,nlocal
 real                      :: ntotal
 logical                   :: iterations_finished
 logical                   :: do_export

 call init_cell_exchange(xrecvbuf,irequestrecv)
 stack_waiting%n = 0
 stack_remote%n  = 0
 stack_redo%n    = 0
#endif

#ifdef MPI
 ! number of local only cells
 nlocal = 0
 call reset_cell_counters
#endif

!$omp parallel default(none) &
!$omp shared(icall,ncells,ll,ifirstincell,xyzh,gradh,iphase,Bevol,divcurlB,id,nprocs,iverbose,iprint) &
#ifdef MPI
!$omp shared(xrecvbuf,xsendbuf,irequestrecv,irequestsend,stack_remote,stack_waiting) &
!$omp shared(stack_redo,iterations_finished,mpiits) &
!$omp reduction(+:nlocal) &
!$omp private(j,k,l,ntotal,do_export) &
#endif
!$omp private(i,remote_export,nneigh,npcell,cell,iamgasi,iamtypei,iactivei,redo_neighbours)
!$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)

    !--skip empty cells AND inactive cells
    if (i <= 0) cycle over_cells
    !
    !--get the neighbour list and fill the cell cache
    !
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
                           remote_export=remote_export)
#ifdef MPI
    if (any(remote_export)) then
       do_export = .true.
    else
       do_export = .false.
    endif
#endif
    cell%icell                   = icell
#ifdef MPI
    cell%owner                   = id
#endif
    cell%nits                    = 0
    cell%nneigh                  = 0
    cell%remote_export(1:nprocs) = remote_export

    call start_cell(cell,iphase,xyzh,Bevol)

    call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
    call get_hmaxcell(icell,cell%hmax)

#ifdef MPI
!$omp critical
    call recv_cells(stack_remote,xrecvbuf,irequestrecv)
!$omp end critical

    if (do_export) then
!$omp critical
       if (stack_waiting%n > 0) call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
       ! make a reservation on the stack
       call reserve_stack(stack_waiting,cell%waiting_index)
       ! export the cell: direction 0 for exporting
       call send_cell(cell,0,irequestsend,xsendbuf)
!$omp end critical
    endif
#endif

    call compute_cell(cell,listneigh,nneigh,Bevol,xyzh,xyzcache)

#ifdef MPI
    if (do_export) then
       ! write directly to stack
       stack_waiting%cells(cell%waiting_index) = cell
    else
       if (.not. do_export) then
#endif
          call store_results(icall,cell,xyzh,gradh,divcurlB)
#ifdef MPI
          nlocal = nlocal + 1
       endif
    endif
#endif
 enddo over_cells
!$omp enddo

#ifdef MPI
!$omp barrier

!$omp single
 if (stack_waiting%n > 0) call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
 call recv_while_wait(stack_remote,xrecvbuf,irequestrecv,irequestsend)
!$omp end single

!$omp single
 if (iverbose>=6) then
    ntotal = real(nlocal) + real(stack_waiting%n)
    if (ntotal > 0) then
       write(iprint,*) id,'local ratio = ',real(nlocal)/ntotal
    else
       write(iprint,*) id,'local ratio = 0'
    endif
 endif

 mpiits = 0
 iterations_finished = .false.
!$omp end single
 remote_its: do while(.not. iterations_finished)
!$omp single
    mpiits = mpiits + 1
    call reset_cell_counters
!$omp end single

    igot_remote: if (stack_remote%n > 0) then
!$omp do schedule(runtime)
       over_remote: do i = 1,stack_remote%n
          cell = stack_remote%cells(i)

          ! icell is unused (-1 here)
          call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
                                  cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti)

          call compute_cell(cell,listneigh,nneigh,Bevol,xyzh,xyzcache)

          cell%remote_export(id+1) = .false.

          ! communication happened while computing contributions to remote cells
!$omp critical
          call recv_cells(stack_waiting,xrecvbuf,irequestrecv)
          call check_send_finished(stack_waiting,irequestsend,irequestrecv,xrecvbuf)
          ! direction return (1)
          call send_cell(cell,1,irequestsend,xsendbuf)
!$omp end critical
       enddo over_remote
!$omp enddo
!$omp barrier
!$omp single
       ! reset remote stack
       stack_remote%n = 0
       ! ensure send has finished
       call check_send_finished(stack_waiting,irequestsend,irequestrecv,xrecvbuf)
!$omp end single
    endif igot_remote
!$omp barrier
!$omp single
    call recv_while_wait(stack_waiting,xrecvbuf,irequestrecv,irequestsend)
    call reset_cell_counters
!$omp end single
    iam_waiting: if (stack_waiting%n > 0) then
!$omp do schedule(runtime)
       over_waiting: do i = 1, stack_waiting%n
          cell = stack_waiting%cells(i)

          if (any(cell%remote_export(1:nprocs))) then
             print*,id,cell%remote_export(1:nprocs)
             print*,id,'mpiits',mpiits
             print*,id,'owner',cell%owner
             print*,id,'icell',cell%icell
             print*,id,'npcell',cell%npcell
             print*,id,'xpos',cell%xpos
             print*,id,'stackpos',i
             print*,id,'waitindex',cell%waiting_index
             call fatal('mhd_divcurlB', 'not all results returned from remote processor')
          endif

          ! communication happened while finishing cell
!$omp critical
          call recv_cells(stack_remote,xrecvbuf,irequestrecv)
!$omp end critical
          call store_results(icall,cell,xyzh,gradh,divcurlB)

       enddo over_waiting
!$omp enddo
!$omp barrier
!$omp single
       ! reset stacks
       stack_waiting%n = 0
       call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
!$omp end single
    endif iam_waiting
!$omp barrier
!$omp single
    call recv_while_wait(stack_remote,xrecvbuf,irequestrecv,irequestsend)
!$omp end single

!$omp single
    if (reduceall_mpi('max',stack_redo%n) > 0) then
       call swap_stacks(stack_waiting, stack_redo)
    else
       iterations_finished = .true.
    endif
    stack_redo%n = 0
!$omp end single

 enddo remote_its

#endif
!$omp end parallel
#ifdef MPI
 call finish_cell_exchange(irequestrecv,xsendbuf)
 call sync_hmax_mpi
#endif

end subroutine calculate_divcurlB

!----------------------------------------------------------------
!+
!  Internal subroutine that computes the contribution to divcurlB
!
!  MAKE SURE THIS ROUTINE IS INLINED BY THE COMPILER
!+
!----------------------------------------------------------------
pure subroutine get_density_sums(i,xpartveci,hi,hi1,hi21,iamtypei,iamgasi,&
                                 listneigh,nneigh,xyzcache,rhosum,xyzh,Bevol,ignoreself)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 use kernel,   only:get_kernel
 use part,     only:iphase,iamgas,iamtype,maxphase,ibasetype,igas,rhoh,massoftype
 use dim,      only:maxp
 integer,      intent(in)    :: i
 real,         intent(in)    :: xpartveci(:)
 real(kind=8), intent(in)    :: hi,hi1,hi21
 integer,      intent(in)    :: iamtypei
 logical,      intent(in)    :: iamgasi
 integer,      intent(in)    :: listneigh(:)
 integer,      intent(in)    :: nneigh
 real,         intent(in)    :: xyzcache(:,:)
 real,         intent(out)   :: rhosum(:)
 real,         intent(in)    :: xyzh(:,:)
 real,         intent(in)    :: Bevol(:,:)
 logical,      intent(in)    :: ignoreself
 integer(kind=1)             :: iphasej
 integer                     :: iamtypej
 integer                     :: j,n
 real                        :: dx,dy,dz,runix,runiy,runiz
 real                        :: rij2,rij,rij1,q2i,qi,rij1grkern
 real                        :: wabi,grkerni,projdB,dBx,dBy,dBz,rhoi,rhoj
 logical                     :: same_type,gas_gas

 rhosum(:) = 0.

 ! defaults for type determination
 ! these are determined from iphase if multiple phases are used
 same_type = .true.
 gas_gas   = .true.
 dx = 0. 
 dy = 0.
 dz = 0.
 loop_over_neigh: do n = 1,nneigh

    j = listneigh(n)
    !--do self contribution separately to avoid problems with 1/sqrt(0.)
    if ((ignoreself) .and. (j==i)) cycle loop_over_neigh

    !--calculate separation
    dx = xpartveci(ixi) - xyzh(1,j)
    dy = xpartveci(iyi) - xyzh(2,j)
    dz = xpartveci(izi) - xyzh(3,j)
#ifdef PERIODIC
    if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
    if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
    if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
    rij2 = dx*dx + dy*dy + dz*dz
    q2i  = rij2*hi21

    if (q2i < radkern2) then
       rij = sqrt(rij2)
       qi  = rij*hi1
       !--kernel
       call get_kernel(q2i,qi,wabi,grkerni)
       !
       ! Density, gradh and div v are only computed using
       ! neighbours of the same type
       !
       if (maxphase==maxp) then
          iphasej   = iphase(j)
          iamtypej  = iamtype(iphasej)
          same_type = ((iamtypei == iamtypej) .or. (ibasetype(iamtypej)==iamtypei))
          gas_gas   = (iamgasi .and. same_type)
       endif

       if (gas_gas) then
          rij1 = 1./(rij + epsilon(rij))
          rij1grkern = rij1*grkerni
          runix = dx*rij1grkern
          runiy = dy*rij1grkern
          runiz = dz*rij1grkern

          rhoi = rhoh(real(hi),  massoftype(igas))
          rhoj = rhoh(xyzh(4,j), massoftype(igas))
          dBx = xpartveci(iBevolxi)*rhoi - Bevol(1,j)*rhoj
          dBy = xpartveci(iBevolyi)*rhoi - Bevol(2,j)*rhoj
          dBz = xpartveci(iBevolzi)*rhoi - Bevol(3,j)*rhoj
          projdB = dBx*runix + dBy*runiy + dBz*runiz

          ! difference operator of divB
          rhosum(idivBi)  = rhosum(idivBi)  + projdB
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

    endif
 enddo loop_over_neigh

end subroutine get_density_sums

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
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
pure subroutine compute_cell(cell,listneigh,nneigh,Bevol,xyzh,xyzcache)
 use part,        only:get_partinfo,iamgas,igas,maxphase
#ifdef MPI
 use io,          only:id
#endif

 type(celldens),  intent(inout)  :: cell
 integer,         intent(in)     :: listneigh(:)
 integer,         intent(in)     :: nneigh
 real,            intent(in)     :: Bevol(:,:)
 real,            intent(in)     :: xyzh(:,:)
 real,            intent(in)     :: xyzcache(isizecellcache,3)
 real(kind=8)                    :: hi
 real(kind=8)                    :: hi1,hi21,hi41
 integer                         :: i,lli,iamtypei
 logical                         :: iactivei,iamgasi,iamdusti,ignoreself

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
    hi41  = hi21*hi21

#ifdef MPI
    if (cell%owner == id) then
       ignoreself = .true.
    else
       ignoreself = .false.
    endif
#else
    ignoreself = .true.
#endif
    call get_density_sums(lli,cell%xpartvec(:,i),hi,hi1,hi21,iamtypei,iamgasi,&
                          listneigh,nneigh,xyzcache,cell%rhosums(:,i),xyzh,Bevol,ignoreself)

 enddo over_parts

end subroutine compute_cell
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine start_cell(cell,iphase,xyzh,Bevol)
 use io,          only:fatal
 use dim,         only:maxp
 use part,        only:maxphase,get_partinfo,igas,iamgas,iamboundary,ibasetype

 type(celldens),     intent(inout) :: cell
 integer(kind=1),    intent(in)    :: iphase(:)
 real,               intent(in)    :: xyzh(:,:)
 real,               intent(in)    :: Bevol(:,:)
 integer :: i,ip
 integer :: iamtypei
 logical :: iactivei,iamgasi,iamdusti

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
    if (.not.iactivei) then ! skip boundary particles + inactive particles
       cycle over_parts
    endif

    cell%npcell = cell%npcell + 1

    cell%arr_index(cell%npcell)            = ip
    cell%iphase(cell%npcell)               = iphase(i)

    cell%xpartvec(:  ,cell%npcell)         = 0. ! to avoid compiler warning
    cell%xpartvec(ixi,cell%npcell)         = xyzh(1,i)
    cell%xpartvec(iyi,cell%npcell)         = xyzh(2,i)
    cell%xpartvec(izi,cell%npcell)         = xyzh(3,i)

    cell%h(cell%npcell)                    = xyzh(4,i)

    if (iamgasi) then
       cell%xpartvec(iBevolxi,cell%npcell) = Bevol(1,i)
       cell%xpartvec(iBevolyi,cell%npcell) = Bevol(2,i)
       cell%xpartvec(iBevolzi,cell%npcell) = Bevol(3,i)
    endif

 enddo over_parts

end subroutine start_cell
!--------------------------------------------------------------------------
!+
! Store div B & curl B
!+
!--------------------------------------------------------------------------
subroutine store_results(icall,cell,xyzh,gradh,divcurlB)
 use part,        only:rhoh,get_partinfo,iamgas,maxphase,massoftype,igas
 use dim,         only:maxp,ndivcurlB

 integer,         intent(in)    :: icall
 type(celldens),  intent(in)    :: cell
 real,            intent(in)    :: xyzh(:,:)
 real(kind=4),    intent(in)    :: gradh(:,:)
 real(kind=4),    intent(inout) :: divcurlB(:,:)
 integer      :: iamtypei,i,lli
 logical      :: iactivei,iamgasi,iamdusti
 real         :: hi,pmassi,rhoi,term
 real         :: divcurlBi(ndivcurlB),rhosum(maxrhosum)

 do i = 1,cell%npcell
    lli  = inodeparts(cell%arr_index(i))

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif

    if (iamgasi) then
       hi     = cell%h(i)
       pmassi = massoftype(iamtypei)
       rhoi   = rhoh(hi,pmassi)
       rhosum = cell%rhosums(:,i)
       term   = cnormk*pmassi*gradh(1,lli)/(rhoi*hi**4)
       call calculate_divcurlB_from_sums(rhosum,term,divcurlBi,ndivcurlB)
       divcurlB(:,lli) = real(divcurlBi(:),kind=kind(divcurlB))
    endif
 enddo

end subroutine store_results
!--------------------------------------------------------------------------

end module mhd_divcurlB
