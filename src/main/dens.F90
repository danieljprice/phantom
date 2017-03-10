!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: densityforce
!
!  DESCRIPTION:
!  This module is the "guts" of the code
!  Calculates density by iteration with smoothing length
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, domain, eos, fastmath, io, io_summary,
!    kernel, linklist, mpiderivs, mpiutils, nicil, options, part, timestep,
!    timing, viscosity
!+
!--------------------------------------------------------------------------
module densityforce
 use dim,      only:maxstrain,maxvxyzu,maxp,ndusttypes
 use part,     only:maxBevol,mhd
 use part,     only:straintensor
 use kernel,   only:cnormk,wab0,gradh0,dphidh0,radkern2
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 public :: densityiterate,get_neighbour_stats,vwave,get_alphaloc

 !--indexing for xpartveci array
 integer, parameter :: maxxpartveci = 11
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
       ipsi = 11

 !--indexing for rhosum array
 integer, parameter :: maxrhosum = 39
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
       irhodusti        = 39

 !--kernel related parameters
 !real, parameter    :: cnormk = 1./pi, wab0 = 1., gradh0 = -3.*wab0, radkern2 = 4.0
 integer, parameter :: isizecellcache = 50000
 integer, parameter :: isizeneighcache = 12000

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
                          fxyzu,fext,alphaind,gradh)
 use dim,       only:maxp,maxneigh,ndivcurlv,ndivcurlB,maxvxyzu,maxalpha, &
                     mhd_nonideal,nalpha,use_dust,use_dustfrac
 use boundary,  only:dxbound,dybound,dzbound
 use eos,       only:get_spsound,get_temperature
 use io,        only:iprint,fatal,iverbose,id,master,real4,warning,error
 use linklist,  only:ncells,ifirstincell,get_neighbour_list,get_hmaxcell,set_hmaxcell
 use options,   only:ieos,tolh,alphaB,alpha,alphamax
 use part,      only:mhd,maxBevol,rhoh,dhdrho,rhoanddhdrho,massoftype,&
                     ll,get_partinfo,iactive,maxgradh,&
                     hrho,iphase,maxphase,igas,idust,iboundary,iamgas,periodic,&
                     set_boundaries_to_active,all_active,n_R,n_electronT
#ifdef FINVSQRT
 use fastmath,  only:finvsqrt
#endif
#ifdef MPI
 use domain,    only:cellbelong
 use io,        only:nprocs
 use mpiutils,  only:fill_buffer,reduceall_mpi
 use mpiderivs, only:init_results_exchange,finish_results_exchange,&
                     send_results,recv_density_results,check_send_finished
 use timing,    only:getused
 use linklist,  only:update_hmax_remote
#endif
 use timestep,  only:rho_dtthresh,mod_dtmax,mod_dtmax_now
 use nicil,     only:nicil_get_ion_n,nicil_translate_error
 use part,      only:ngradh,dustfrac
 use viscosity, only:irealvisc,bulkvisc,shearparam
 use io_summary,only:summary_variable,iosumhup,iosumhdn
 integer,      intent(in)    :: icall,npart,nactive
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(in)    :: vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real(kind=4), intent(in)    :: Bevol(:,:)
 real(kind=4), intent(out)   :: divcurlv(:,:)
 real(kind=4), intent(out)   :: divcurlB(:,:)
 real(kind=4), intent(out)   :: alphaind(:,:)
 real(kind=4), intent(out)   :: gradh(:,:)
 real,         intent(out)   :: stressmax

 integer, parameter :: maxdensits = 50

 real :: xpartveci(maxxpartveci)
 real :: rhosum(maxrhosum)
 integer, save :: listneigh(maxneigh)
 real,   save :: xyzcache(3,isizecellcache)
 real,   save :: dxcache(7,isizeneighcache)
!$omp threadprivate(dxcache,xyzcache,listneigh)

 integer :: i,itsdensity,icell,nneightemp
 integer :: nneigh,nneighi,np
 integer :: nwarnup,nwarndown,nwarnroundoff
 integer :: ierr
 real :: hi,hi_old,hminbisect,hmaxbisect
 real(kind=8) :: hi1,hi21,hi31,hi41
 real :: hnew
 real :: omegai,dhdrhoi,rhohi,func,dfdh1
 real :: rho1i,pmassi,rhodusti(ndusttypes)
 real :: hmaxcelli,hcut
 real :: spsoundi
 real :: divcurlvi(5)
 real :: divcurlBi(ndivcurlB),gradBi,Bi,alphaBi
 real :: term
 real :: denom
 real         :: straini(6),rmatrix(6)
 real(kind=8) :: rhoi
 real(kind=8) :: gradhi,gradsofti
 real         :: Bxi,Byi,Bzi,psii
 real         :: temperaturei
 real         :: rhomax,xi_limiter

 integer :: iamtypei,isize
 logical :: ifilledcellcache,ifilledneighcache,moreincell,activecell,activeneighbsonly
 logical :: iactivei,iamgasi,iamdusti,getdv,realviscosity,igotrmatrix,bisection
 logical :: getdB,converged
#ifdef MPI
 integer, parameter :: maxrecv = 3 + ndivcurlv + ndivcurlB + 6
 integer      :: nactive_thisproc,n,nactivetot
 real         :: xrecvbuf(maxrecv,nprocs),xbuf(maxrecv)
 integer      :: ireqrecv(nprocs),ireqsend(nprocs),nrecv
 real(kind=4) :: t1,t2,t3
 logical      :: have_sent
#endif

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
 getdv = ((maxalpha==maxp .or. ndivcurlv >= 4) .and. icall <= 1) .or. &
         (realviscosity .and. maxstrain==maxp)
 if (getdv .and. ndivcurlv < 1) call fatal('densityiterate','divv not stored but it needs to be')
 getdB = (mhd .and. (ndivcurlB >= 4 .or. mhd_nonideal))

 ! set stress terms to zero initially (make sure they are firstprivate)
 straini(:) = 0.
 if ( all_active ) stressmax  = 0.   ! condition is required for independent timestepping

#ifdef MPI
 nactive_thisproc = 0
 nrecv = 0
 call getused(t1)
 have_sent = .false.
 call init_results_exchange(nprocs,xrecvbuf,ireqrecv)
#endif
 rhomax = 0.0

!$omp parallel default(none) &
!$omp shared(ncells,ll,ifirstincell,npart,set_boundaries_to_active) &
!$omp shared(tolh,iprint) &
!$omp shared(xyzh,vxyzu,divcurlv,massoftype,iphase,fxyzu,fext) &
!$omp shared(gradh,straintensor,dxbound,dybound,dzbound) &
!$omp shared(irealvisc,realviscosity,bulkvisc,shearparam,getdv,getdB,icall) &
!$omp shared(ieos,n_R,n_electronT,dustfrac) &
!$omp private(isize) &
!$omp firstprivate(straini) &
!$omp private(icell,i,hmaxcelli,iamtypei) &
!$omp private(ifilledcellcache,moreincell,activeneighbsonly) &
!$omp private(hi,hi_old,ifilledneighcache,converged) &
!$omp private(itsdensity,hi1,hi21,hi31,hi41,rhoi,rho1i,gradhi,gradsofti,pmassi) &
!$omp private(nneighi,nneigh,nneightemp) &
!$omp private(hcut,rhohi,dhdrhoi,omegai,func,dfdh1,bisection,hminbisect,hmaxbisect) &
!$omp private(hnew,spsoundi,temperaturei) &
!$omp private(divcurlvi) &
!$omp private(divcurlBi) &
!$omp private(activecell,iactivei,iamgasi,iamdusti) &
!$omp private(igotrmatrix) &
!$omp private(rhosum,xpartveci,rhodusti) &
!$omp private(term,xi_limiter) &
!$omp private(ierr) &
!$omp reduction(+:ncalc,nrelink,nneighact,nneightry,nwarnup,nwarndown) &
!$omp reduction(+:np) &
!$omp reduction(max:maxneightry,maxneighact) &
!$omp shared(Bevol) &
!$omp shared(divcurlB,alphaind,alphaB,alpha,alphamax) &
!$omp private(psii,Bxi,Byi,Bzi,gradBi,Bi,alphaBi) &
!$omp private(denom,rmatrix) &
!$omp reduction(max:stressmax,rhomax) &
#ifdef MPI
!$omp reduction(+:nwarnroundoff) &
!$omp shared(id,ireqrecv,ireqsend,xrecvbuf,xbuf,nprocs,have_sent) &
!$omp private(n) &
!$omp reduction(+:nactive_thisproc,nrecv)
#else
!$omp reduction(+:nwarnroundoff)
#endif
!$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)

    !--skip empty cells AND inactive cells
    if (i <= 0) cycle over_cells

#ifdef MPI
    !
    ! skip leaf cells that belong on other procs
    !
    if (.not.cellbelong(icell,id,nprocs,ncells)) cycle over_cells
#endif

    activecell = .true.

    !--if only one particle in cell, do not waste time filling the cell cache
    !  (actually on for this call we fill the cell cache anyway mainly
    !   because for periodic bcs this is a more efficient way of getting
    !   the modulus across the boundary)
    moreincell = (ll(i) /= 0)

    call get_hmaxcell(icell,hmaxcelli)
    !
    !--get the neighbour list and fill the cell cache
    !
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,isizecellcache,.false.,getj=.false.)
    ifilledcellcache = .true.

    !--go through cell list
    over_parts: do while (i /= 0)

    if (i < 0) then    ! i < 0 indicates inactive
       i = ll(abs(i))   ! this is a quicker way of skipping inactives
       cycle over_parts ! but does not apply to first particle in cell
    endif

    if (i > npart) then ! skip mpi ghosts
       i = ll(i)
       cycle over_parts
    endif

    !--get individual timestep/ multiphase information (querying iphase)
    if (maxphase==maxp) then
       call get_partinfo(iphase(i),iactivei,iamdusti,iamtypei)
       iamgasi = iamgas(iphase(i))
    else
       iactivei = .true.
       iamtypei = igas
       iamgasi  = .true.
    endif
    if (.not.iactivei) then ! handles case where first particle in active cell is inactive
       i = ll(i)
       cycle over_parts
    endif
    ! Individual timesteps:  boundary particles are active only on the first step,
    !                        thus treat like gas to initialise the arrays; else they will
    !                        be skipped with the above command
    ! Global timesteps:      boundary particles are always active, thus use the logical
    !                        deriv_first_call to control initialisation
    if (iamtypei==iboundary) then
       if (set_boundaries_to_active) then
          iactivei = .true.
          iamtypei = igas
          iamgasi  = .true.
       else
          i = ll(i)
          cycle over_parts
       endif
    endif

    pmassi = massoftype(iamtypei)
    xpartveci(ixi)  = xyzh(1,i)
    xpartveci(iyi)  = xyzh(2,i)
    xpartveci(izi)  = xyzh(3,i)
    hi              = xyzh(4,i)
    hi_old          = hi
    if (hi < 0.) call fatal('densityiterate','negative smoothing length',i,var='h',val=hi)
    xpartveci(ivxi) = vxyzu(1,i)
    xpartveci(ivyi) = vxyzu(2,i)
    xpartveci(ivzi) = vxyzu(3,i)
    if (maxvxyzu >= 4) xpartveci(ieni) = vxyzu(4,i)

    if (mhd) then
       if (iamgasi) then
          xpartveci(iBevolxi) = Bevol(1,i)
          xpartveci(iBevolyi) = Bevol(2,i)
          xpartveci(iBevolzi) = Bevol(3,i)
          if (maxBevol >= 4) xpartveci(ipsi) = Bevol(4,i)
          if (maxBevol < 3 .or. maxBevol > 4) call fatal('densityiterate','error in maxBevol setting')
       else
          xpartveci(iBevolxi:ipsi) = 0. ! to avoid compiler warning
       endif
    endif

    ifilledneighcache = .false.

    itsdensity = 0

    bisection = .false. ! start with Newton-Raphson iterations
    iterate: do itsdensity = 1, maxdensits
       ncalc = ncalc + 1
       hi1   = 1./hi
       hi21  = hi1*hi1
       hi31  = hi1*hi21
       hi41  = hi21*hi21
!
!--if h exceeds current cell max, re-calculate the list of neighbours for this cell
!  (common to all particles in the cell)
!
       if (hi > hmaxcelli) then
          !write(iprint,*) 'relinking cell ',icell,' part ',i,' hi = ',hi,' hmaxcell =',hmaxcelli
          hcut = 1.01*hi
          call set_hmaxcell(icell,hcut)
          isize = isizecellcache
          if (.not.moreincell) isize = 0 ! do not fill the cell cache for last particle in cell

          call get_neighbour_list(icell,listneigh,nneightemp,xyzh,xyzcache,&
               isize,activeneighbsonly,hmaxcelli,getj=.false.)

          hmaxcelli = hcut
          if (nneightemp >= 0) then  ! less than zero means neighbour list unchanged
             nneigh = nneightemp
             nrelink = nrelink + 1
             if (moreincell) then
                ifilledcellcache = .true.
             else
                ifilledcellcache = .false.
             endif
             ifilledneighcache = .false.
          endif
          !write(iprint,*) 'new nneigh = ',nneigh
       endif
!
!--get contribution to density sums from
!  current particle's neighbours (NB: neighbour list includes self)
!
       call get_density_sums(i,xpartveci,hi1,hi21,iamtypei,iamgasi,listneigh,nneigh, &
                             nneighi,dxcache,xyzcache,rhosum,&
                             ifilledcellcache,ifilledneighcache,getdv,getdB,realviscosity, &
                             xyzh,vxyzu,Bevol,fxyzu,fext)

       if (moreincell) ifilledcellcache = .true.

       !--compute rho and gradhi adding self-contribution
       rhoi      = cnormk*pmassi*(rhosum(irhoi) + wab0)*hi31
       gradhi    = cnormk*pmassi*(rhosum(igradhi) + gradh0)*hi41
       gradsofti = pmassi*(rhosum(igradsofti) + dphidh0)*hi21 ! NB: no cnormk in gradsoft

       rhohi   = rhoh(hi,pmassi)
       dhdrhoi = dhdrho(hi,pmassi)
       omegai  = 1. - dhdrhoi*gradhi

       func = rhohi - rhoi

       if (bisection) then
          if (func < 0.) then
             hmaxbisect = hi
          else
             hminbisect = hi
          endif
          hnew = 0.5*(hminbisect + hmaxbisect)
!          write(iprint,*) 'bisection (',itsdensity,'): hnew = ',&
!                          hnew,i,hi,func,hminbisect,hmaxbisect,nneighi
       else
          if (omegai > tiny(omegai)) then
             dfdh1 = dhdrhoi/omegai
          else
             dfdh1 = dhdrhoi/abs(omegai + epsilon(omegai))
          endif

          hnew = hi - func*dfdh1
          !
          !--stop Newton-Raphson from taking huge (probably wrong) jumps
          !
          if (hnew > 1.2*hi) then
             nwarnup   = nwarnup + 1
             hnew      = 1.2*hi
          elseif (hnew < 0.8*hi) then
             nwarndown = nwarndown + 1
             hnew      = 0.8*hi
          endif
       endif

       gradhi = 1./omegai
       gradsofti = gradsofti*dhdrhoi

       if (icall==0 .and. itsdensity==1) then
          converged = .true.   ! icall = 0 does not do any iterations
       else
          converged = ((abs(hnew-hi)/hi_old) < tolh .and. omegai > 0. .and. hi > 0.)
       endif

       if (converged) then
          np          = np + 1
          nneighact   = nneighact + nneighi
          nneightry   = nneightry + nneigh
          maxneightry = max(nneigh,int(maxneightry))
          maxneighact = max(nneighi,maxneighact)
          exit iterate
       elseif (itsdensity < maxdensits) then
          if ((omegai <= 0. .or. hi <= 0.).and. .not.bisection) then
             !--start bisection iterations
             bisection = .true.
             hminbisect = 0.
             hmaxbisect = 1.e6
             if (periodic) hmaxbisect = min(hmaxbisect,0.25*max(dxbound,dybound,dzbound))
             !--reset h to original value
             hnew = hi_old
          endif
          if (itsdensity >= maxdensits-5 .and. itsdensity > 10) then
             write(iprint,*) 'error in (hnew - hi)/hi_old = ',abs(hnew-hi)/hi_old
             call warning('densityiterate','density iteration failing',i,'hnew',hnew)
          endif
          hi = hnew
       endif
    enddo iterate

    if (itsdensity >= maxdensits) then
       write(iprint,*) 'ERROR: density iteration failed after ',itsdensity,' iterations'
       write(iprint,*) 'hnew = ',hnew,' hi_old = ',hi_old,' nneighi = ',nneighi
       write(iprint,*) 'rhoi = ',rhoi,' gradhi = ',gradhi
       write(iprint,*) 'error = ',abs(hnew-hi)/hi_old,' tolh = ',tolh
       write(iprint,*) 'itype = ',iamtypei
       write(iprint,*) 'x,y,z = ',xyzh(1:3,i)
       write(iprint,*) 'v_x,v_y,v_z = ',vxyzu(1:3,i)
       if (maxvxyzu >= 4) write(iprint,*) 'u = ',vxyzu(4,i)
       write(iprint,*) 'c_s         = ',get_spsound(ieos,xyzh(:,i),real(rhoi),vxyzu(:,i))
       write(iprint,*) 'temperature = ',get_temperature(ieos,xyzh(:,i),real(rhoi),vxyzu(:,i))
       call fatal('densityiterate','could not converge in density',i,'error',abs(hnew-hi)/hi_old)
    endif
!
!--store final results of density iteration
!
    if (hi < 0.) call fatal('densityiterate','hi < 0 after iterations',i,var='h',val=hi)
    xyzh(4,i) = hrho(rhoi,pmassi)   ! we use h to get rho for inactives, so save "rho" rather than h.
    if (xyzh(4,i) < 0.) call fatal('densityiterate','setting negative h from hrho',i,var='rhoi',val=real(rhoi))
    if (maxgradh==maxp) then
       gradh(1,i) = real(gradhi,kind=kind(gradh))
#ifdef GRAVITY
       gradh(2,i) = real(gradsofti,kind=kind(gradh))
#endif
    endif
    rho1i  = 1./rhoi
    rhomax = max(rhomax,real(rhoi))
    if (use_dust .and. .not. use_dustfrac .and. iamgasi) then
       !
       ! for 2-fluid dust compute dust density on gas particles
       ! and store it in dustfrac as dust-to-gas ratio
       ! so that rho times dustfrac gives dust density
       !
       rhodusti(:) = cnormk*massoftype(idust)*(rhosum(irhodusti))*hi31
       dustfrac(:,i) = rhodusti(:)*rho1i ! dust-to-gas ratio
       if ( ndusttypes>1 ) call fatal('dens','2-fluid not compatible with ndusttypes > 1')
    endif
!
! store divv and curl v and related quantities
!
    igotrmatrix = .false.
    term = cnormk*pmassi*gradhi*rho1i*hi41
    if (getdv) then
       call calculate_rmatrix_from_sums(rhosum,denom,rmatrix,igotrmatrix)
       call calculate_divcurlv_from_sums(rhosum,term,divcurlvi,xi_limiter,ndivcurlv,denom,rmatrix)
       divcurlv(1:ndivcurlv,i) = real(divcurlvi(1:ndivcurlv),kind=kind(divcurlv)) ! save to global memory
       !
       ! Cullen & Dehnen (2010) viscosity switch, set alphaloc
       !
       if (nalpha >= 2 .and. iamgasi) then
          spsoundi = get_spsound(ieos,xyzh(:,i),real(rhoi),vxyzu(:,i))
          alphaind(2,i) = real4(get_alphaloc(divcurlvi(5),spsoundi,hi,xi_limiter,alpha,alphamax))
       endif
    else ! we always need div v for h prediction
       if (ndivcurlv >= 1) divcurlv(1,i) = -real4(rhosum(idivvi)*term)
       if (nalpha >= 2) alphaind(2,i) = 0.
    endif
!
! store div B, curl B and related quantities
!
    if (mhd .and. iamgasi) then
       Bxi = xpartveci(iBevolxi)
       Byi = xpartveci(iBevolyi)
       Bzi = xpartveci(iBevolzi)
       if (maxBevol >= 4) psii = xpartveci(ipsi)

       if (getdB) then
          term = cnormk*pmassi*gradhi*rho1i*hi41
          call calculate_divcurlB_from_sums(rhosum,term,divcurlBi,gradBi,ndivcurlB)
          divcurlB(:,i) = real(divcurlBi(:),kind=kind(divcurlB))
       else
          divcurlBi(:) = 0.
       endif
    endif
!
!--get strain tensor from summations
!
    realvisc: if (realviscosity .and. maxstrain==maxp .and. iamgasi) then
       if (.not.igotrmatrix) &
          call calculate_rmatrix_from_sums(rhosum,denom,rmatrix,igotrmatrix)

       term = cnormk*pmassi*gradhi*rho1i*hi41
       call calculate_strain_from_sums(rhosum,term,denom,rmatrix,straini)

       ! check for negative stresses to prevent tensile instability
       call get_max_stress(straini,divcurlvi(1),rho1i,stressmax,shearparam,bulkvisc)

       ! store strain tensor
       straintensor(:,i) = real(straini(:),kind=kind(straintensor))
    endif realvisc
!
!--calculate number densities
!
    if (mhd_nonideal .and. iamgasi) then
       temperaturei = get_temperature(ieos,xpartveci(ixi:izi),real(rhoi),vxyzu(:,i))
       call nicil_get_ion_n(real(rhoi),temperaturei,n_R(:,i),n_electronT(i),ierr)
       if (ierr/=0) then
          call nicil_translate_error(ierr)
          call fatal('densityiterate','error in calcuating number densities for non-ideal MHD')
       endif
    endif
!
#ifdef MPI
    nactive_thisproc = nactive_thisproc + 1
!$omp critical
    if (have_sent) then
       ! have to make sure last send completed before sending another
       call check_send_finished(nprocs,ireqsend,ireqrecv,xrecvbuf,nrecv,getdv,getdB,realviscosity)
    endif

    ! broadcast particle result to all procs
    n = 0
    call fill_buffer(xbuf,xyzh(4,i),n)
    call fill_buffer(xbuf,gradhi,n)
    call fill_buffer(xbuf,gradsofti,n)
    call fill_buffer(xbuf,divcurlvi,n)
    if (mhd .and. getdB) then
       call fill_buffer(xbuf,divcurlBi,n)
    endif
    if (realviscosity .and. maxstrain==maxp) call fill_buffer(xbuf,straini,n)
    if (n > size(xbuf)) call fatal('densityiterate','mpi buffer size exceeded: increase ipartbufsize and recompile')
    call send_results(nprocs,xbuf,n,i,ireqsend)
    have_sent = .true.

    ! check for receives
    call recv_density_results(nprocs,xrecvbuf,ireqrecv,nrecv,getdv,getdB,realviscosity)
!$omp end critical
#endif

    !--go to next particle in cell link list
    i = ll(i)

    enddo over_parts
 enddo over_cells
!$omp enddo
!$omp end parallel

#ifdef MPI
 nactivetot = int(reduceall_mpi('+',nactive_thisproc))
 call getused(t2)

 ! sanity check
 if (iverbose >= 2) then
    print*,' thread ',id,' nactive = ',nactive_thisproc,' of ',nactivetot,' time = ',t2-t1,'s'
 endif
 !if (nactivetot /= nactive) call fatal('densityiterate','MPI error: nactive /= nsent',var='nactive',ival=nactivetot)
!
!--check for remaining receives onto this proc.
!
 call finish_results_exchange(nprocs,xrecvbuf,ireqrecv,nrecv,nactivetot-nactive_thisproc,getdv,getdB,realviscosity)
 call getused(t3)

 if (nactivetot - nrecv < 0) call fatal('densityiterate','received more than nactive',var='nrecv',ival=nrecv)
 if (iverbose >= 2) then
    print "(1x,i2,a,f6.2,a,f6.2,a)", &
      id,' time spent waiting for receives = ',t3-t2,'s = ',100.*(t3-t2)/((t3-t1) + epsilon(t1)),'%'
 endif

 ! reduce max stress across MPI procs
 if (realviscosity .and. maxstrain==maxp) then
    stressmax = reduceall_mpi('max',stressmax)
 endif
!
!--update tree with new hmax
!
 call update_hmax_remote(ncells)
#endif

 if (realviscosity .and. maxstrain==maxp .and. stressmax > 0. .and. iverbose > 0 .and. id==master) then
    call warning('force','applying negative stress correction',var='max',val=-stressmax)
 endif
!
!--determine if we need to decrease dtmax at the next opportunity
!
 if (mod_dtmax .and. rhomax > rho_dtthresh) mod_dtmax_now = .true.
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

end subroutine densityiterate

!----------------------------------------------------------------
!+
!  Internal subroutine that computes the contribution to
!  the density sums from a list of neighbours
!
!  MAKE SURE THIS ROUTINE IS INLINED BY THE COMPILER
!+
!----------------------------------------------------------------
pure subroutine get_density_sums(i,xpartveci,hi1,hi21,iamtypei,iamgasi,listneigh,nneigh,nneighi, &
                             dxcache,xyzcache,rhosum, &
                             ifilledcellcache,ifilledneighcache,getdv,getdB, &
                             realviscosity,xyzh,vxyzu,Bevol,fxyzu,fext)
#ifdef PERIODIC
  use boundary, only:dxbound,dybound,dzbound
#endif
#ifdef FINVSQRT
  use fastmath, only:finvsqrt
#endif
  use kernel,   only:get_kernel,get_kernel_grav1
  use part,     only:iphase,iamgas,iamtype,maxphase,iboundary,idust
  use dim,      only:ndivcurlv,gravity,maxp,nalpha,use_dustfrac,use_dust
  integer,      intent(in)    :: i
  real,         intent(in)    :: xpartveci(:)
  real(kind=8), intent(in)    :: hi1,hi21
  integer,      intent(in)    :: iamtypei
  logical,      intent(in)    :: iamgasi
  integer,      intent(in)    :: listneigh(:)
  integer,      intent(in)    :: nneigh
  integer,      intent(out)   :: nneighi
  real,         intent(inout) :: dxcache(:,:)
  real,         intent(in)    :: xyzcache(:,:)
  real,         intent(out)   :: rhosum(:)
  logical,      intent(in)    :: ifilledcellcache
  logical,      intent(inout) :: ifilledneighcache
  logical,      intent(in)    :: getdv,realviscosity
  logical,      intent(in)    :: getdB
  real,         intent(in)    :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:)
  real(kind=4), intent(in)    :: Bevol(:,:)
  integer(kind=1)             :: iphasej
  integer                     :: iamtypej
  integer                     :: j,n
  real                        :: dx,dy,dz,runix,runiy,runiz
  real                        :: rij2,rij,rij1,q2i,qi,q2prev,rij1grkern
  real                        :: wabi,grkerni,dwdhi,dphidhi
  real                        :: projv,dvx,dvy,dvz,dax,day,daz
  real                        :: projdB,dBx,dBy,dBz,fxi,fyi,fzi,fxj,fyj,fzj
  logical                     :: same_type,gas_gas

  rhosum(:) = 0.
  nneighi   = 1   ! self

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
     fxi = fxyzu(1,i) + fext(1,i)
     fyi = fxyzu(2,i) + fext(2,i)
     fzi = fxyzu(3,i) + fext(3,i)
  endif

  loop_over_neigh: do n = 1,nneigh

     j = listneigh(n)
     !--do self contribution separately to avoid problems with 1/sqrt(0.)
     if (j==i) cycle loop_over_neigh

     if (ifilledneighcache .and. n <= isizeneighcache) then
        rij2 = dxcache(1,n)
     else
        if (ifilledcellcache .and. n <= isizecellcache) then
           ! positions from cache are already mod boundary
           dx = xpartveci(ixi) - xyzcache(1,n)
           dy = xpartveci(iyi) - xyzcache(2,n)
           dz = xpartveci(izi) - xyzcache(3,n)
        else
           dx = xpartveci(ixi) - xyzh(1,j)
           dy = xpartveci(iyi) - xyzh(2,j)
           dz = xpartveci(izi) - xyzh(3,j)
#ifdef PERIODIC
           if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
           if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
           if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
        endif

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
           same_type = ((iamtypei == iamtypej) .or. (iamtypej==iboundary))
           gas_gas   = (iamgasi .and. same_type)  ! this ensure that boundary particles are included in gas_gas calculations
        endif

        sametype: if (same_type) then
           dwdhi = (-qi*grkerni - 3.*wabi)
           rhosum(irhoi)      = rhosum(irhoi) + wabi
           rhosum(igradhi)    = rhosum(igradhi) + dwdhi
           rhosum(igradsofti) = rhosum(igradsofti) + dphidhi
           nneighi            = nneighi + 1
           !
           ! calculate things needed for viscosity switches
           ! and real viscosity
           !
           if (getdv .or. getdB) then

              rij1 = 1./(rij + epsilon(rij))
              if (ifilledneighcache .and. n <= isizeneighcache) then
              !--dx,dy,dz are either in neighbour cache or have been calculated
                 dx = dxcache(5,n)
                 dy = dxcache(6,n)
                 dz = dxcache(7,n)
              endif
              rij1grkern = rij1*grkerni
              runix = dx*rij1grkern
              runiy = dy*rij1grkern
              runiz = dz*rij1grkern

              if (getdv) then
                 !--get dv and den
                 dvx = xpartveci(ivxi) - vxyzu(1,j)
                 dvy = xpartveci(ivyi) - vxyzu(2,j)
                 dvz = xpartveci(ivzi) - vxyzu(3,j)
                 projv = dvx*runix + dvy*runiy + dvz*runiz
                 rhosum(idivvi) = rhosum(idivvi) + projv

                 if (realviscosity .or. ndivcurlv > 1 .or. nalpha > 1) then
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
                       fyj = fxyzu(1,j) + fext(1,j)
                       fzj = fxyzu(1,j) + fext(1,j)
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
                 dBx = xpartveci(iBevolxi) - real(Bevol(1,j))
                 dBy = xpartveci(iBevolyi) - real(Bevol(2,j))
                 dBz = xpartveci(iBevolzi) - real(Bevol(3,j))
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

           endif
        elseif (use_dust .and. iamgasi .and. iamtypej==idust) then
           rhosum(irhodusti) = rhosum(irhodusti) + wabi
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

  ifilledneighcache = .true.

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
pure subroutine calculate_divcurlv_from_sums(rhosum,termnorm,divcurlvi,xi_limiter,ndivcurlv,denom,rmatrix)
  use part, only:nalpha
  integer, intent(in)  :: ndivcurlv
  real,    intent(in)  :: rhosum(:),denom,rmatrix(6)
  real,    intent(in)  :: termnorm
  real,    intent(out) :: divcurlvi(5),xi_limiter
  real :: div_a
  real :: gradaxdx,gradaxdy,gradaxdz,gradaydx,gradaydy,gradaydz,gradazdx,gradazdy,gradazdz
  real :: ddenom,gradvxdxi,gradvxdyi,gradvxdzi
  real :: gradvydxi,gradvydyi,gradvydzi,gradvzdxi,gradvzdyi,gradvzdzi
  real :: dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi
  real :: fac,traceS,txy,txz,tyz,txx,tyy,tzz!,Ri
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
     if (use_exact_linear) then
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
     !if (divcurlvi(1) < 0.) then
     !   Ri = -1.
     !else
     !   Ri = 1.
     !endif
     txx = dvxdxi - divcurlvi(1)/3.
     tyy = dvydyi - divcurlvi(1)/3.
     tzz = dvzdzi - divcurlvi(1)/3.
     txy = 0.5*(dvxdyi + dvydxi)
     txz = 0.5*(dvxdzi + dvzdxi)
     tyz = 0.5*(dvydzi + dvzdyi)
     fac    = max(-divcurlvi(1),0.)**2 !(2.*(1. - Ri)**4*divcurlvi(1))**2
     !traceS = txx**2 + tyy**2 + tzz**2 + 2.*(txy**2 + txz**2 + tyz**2)
     !traceS = txy**2 + txz**2 + tyz**2
     traceS = (dvzdyi - dvydzi)**2 + (dvxdzi - dvzdxi)**2 + (dvydxi - dvxdyi)**2
     if (fac + traceS > 0.) then
        xi_limiter = fac/(fac + traceS)
     else
        xi_limiter = 1.
     endif
  else
     xi_limiter = 1.
  endif

end subroutine calculate_divcurlv_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract div, curl and grad B from the sums
!  calculated during the density loop
!+
!----------------------------------------------------------------
pure subroutine calculate_divcurlB_from_sums(rhosum,termnorm,divcurlBi,gradBi,ndivcurlB)
  integer, intent(in)  :: ndivcurlB
  real,    intent(in)  :: rhosum(:)
  real,    intent(in)  :: termnorm
  real,    intent(out) :: divcurlBi(ndivcurlB),gradBi

  ! we need these for adaptive resistivity switch
  if (ndivcurlB >= 1) divcurlBi(1) = -rhosum(idivBi)*termnorm
  if (ndivcurlB >= 4) then
     divcurlBi(2) = -(rhosum(idBzdyi) - rhosum(idBydzi))*termnorm
     divcurlBi(3) = -(rhosum(idBxdzi) - rhosum(idBzdxi))*termnorm
     divcurlBi(4) = -(rhosum(idBydxi) - rhosum(idBxdyi))*termnorm
  endif
  gradBi = rhosum(idBxdxi) * rhosum(idBxdxi) &
         + rhosum(idBxdyi) * rhosum(idBxdyi) &
         + rhosum(idBxdzi) * rhosum(idBxdzi) &
         + rhosum(idBydxi) * rhosum(idBydxi) &
         + rhosum(idBydyi) * rhosum(idBydyi) &
         + rhosum(idBydzi) * rhosum(idBydzi) &
         + rhosum(idBzdxi) * rhosum(idBzdxi) &
         + rhosum(idBzdyi) * rhosum(idBzdyi) &
         + rhosum(idBzdzi) * rhosum(idBzdzi)
  gradBi = sqrt(termnorm * termnorm * gradBi)

end subroutine calculate_divcurlB_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract the strain tensor from summations
!  calculated during the density loop.
!+
!----------------------------------------------------------------
subroutine calculate_strain_from_sums(rhosum,termnorm,denom,rmatrix,strain)
  real, intent(in)  :: rhosum(:)
  real, intent(in)  :: termnorm,denom
  real, intent(in)  :: rmatrix(6)
  real, intent(out) :: strain(6)

  real :: ddenom,gradvxdxi,gradvxdyi,gradvxdzi
  real :: gradvydxi,gradvydyi,gradvydzi,gradvzdxi,gradvzdyi,gradvzdzi
  real :: dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi

!  if (abs(denom) > tiny(denom)) then ! do exact linear first derivatives
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

  strain(1) = (dvxdxi + dvxdxi)
  strain(2) = (dvxdyi + dvydxi)
  strain(3) = (dvxdzi + dvzdxi)
  strain(4) = (dvydyi + dvydyi)
  strain(5) = (dvydzi + dvzdyi)
  strain(6) = (dvzdzi + dvzdzi)

  return
end subroutine calculate_strain_from_sums

!----------------------------------------------------------------
!+
!  Internal subroutine to get maximum of stress tensor
!  (to avoid tensile instability)
!+
!----------------------------------------------------------------
pure subroutine get_max_stress(strain,divvi,rho1i,stressmax,shearvisc,bulkvisc)
 real, intent(in)    :: strain(6), divvi, rho1i, shearvisc, bulkvisc
 real, intent(inout) :: stressmax
 real :: strainmax,stressiso

 ! shearvisc = eta/rho, so this is eta/rho**2
 strainmax = -shearvisc*rho1i*maxval(strain) ! 1/rho*L^2/T*1/L*L/T = 1/rho*L^2/T^2
 stressiso = (2./3.*shearvisc - bulkvisc)*divvi*rho1i ! NB: divv -ve at high Mach no.

 ! we initialise stressmax to zero, as for the purpose of preventing the
 ! tensile instability we only care if the total stress is negative
 ! if stress tensor is positive, don't need correction (stressmax=0)
 stressmax = max(stressmax,-(stressiso + strainmax))

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

 return
end subroutine exactlinear

!----------------------------------------------------------------
!+
!  query function to return the wave speed
!  (called from step for decay timescale in alpha switches)
!+
!----------------------------------------------------------------
real function vwave(xyzhi,pmassi,ieos,vxyzui,Bxyzi)
 use eos,  only:equationofstate
 use part, only:maxp,mhd,maxvxyzu,rhoh
 real,         intent(in) :: xyzhi(4),pmassi
 real,         intent(in) :: vxyzui(maxvxyzu)
 integer,      intent(in) :: ieos
 real(kind=4), intent(in), optional :: Bxyzi(3)
 real :: spsoundi,hi,rhoi,ponrhoi,valfven2i

 hi = xyzhi(4)
 rhoi = rhoh(hi,pmassi)
 if (maxvxyzu==4) then
    call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xyzhi(1),xyzhi(2),xyzhi(3),vxyzui(4))
 else
    call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xyzhi(1),xyzhi(2),xyzhi(3))
 endif
 if (present(Bxyzi)) then
    valfven2i = (Bxyzi(1)**2 + Bxyzi(2)**2 + Bxyzi(3)**2)/rhoi
    vwave = sqrt(valfven2i**2 + spsoundi**2)
 else
    vwave = spsoundi
 endif

end function vwave

!-------------------------------------------------------------------------------
!+
!  function to return alphaloc from known values of d(divv)/dt and sound speed
!  for use in Cullen & Dehnen (2010) switch
!+
!-------------------------------------------------------------------------------
pure real function get_alphaloc(divvdti,spsoundi,hi,xi_limiter,alphamin,alphamax)
 !use kernel, only:radkern
 real, intent(in) :: divvdti,spsoundi,hi,xi_limiter,alphamin,alphamax
 real :: source
 real :: temp

 source = 10.*hi**2*xi_limiter*max(-divvdti,0.)
 temp = spsoundi**2 !+ source
 if (temp > epsilon(temp)) then
    get_alphaloc = max(min(source/temp,alphamax),alphamin)
 else
    get_alphaloc = alphamin
 endif

end function get_alphaloc

!----------------------------------------------------------------
!+
!  subroutine to reduce and print warnings across processors
!  related to h-rho iterations
!+
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

 return
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

 nptot = reduce_mpi('+',np)
 nneightry   = reduce_mpi('+',nneightry)
 nneighact   = reduce_mpi('+',nneighact)
 maxneightry = reduce_mpi('max',maxneightry)
 maxneighact = int(reduce_mpi('max',maxneighact))
 nrelink     = int(reduce_mpi('+',nrelink))
 ncalc       = reduce_mpi('+',ncalc)
 if (id==master .and. iverbose >= 2 .and. nptot > 0 .and. nneighact > 0) then
    write(iprint,"(1x,a,f11.2,2(a,f7.2))") 'trial neigh mean  :',nneightry/real(nptot),&
                 ', real neigh mean = ',nneighact/real(nptot), &
                 ' ratio try/act= ',nneightry/real(nneighact)
    write(iprint,"(1x,a,i11,a,i8)")   'trial neigh max   :',maxneightry,', max real neigh = ',maxneighact
    write(iprint,"(1x,a,i11,a,f7.3)") 'n neighbour calls :',nrelink, ', mean per part   = ',nrelink/real(nptot) + 1
    write(iprint,"(1x,a,i11,a,f7.3)") 'n density calcs   :',ncalc,', mean per part   = ',ncalc/real(nptot)
 endif

end subroutine reduce_and_print_neighbour_stats

end module densityforce
