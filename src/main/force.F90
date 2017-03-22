!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, chem, cooling, dim, domain, dust, eos, fastmath,
!    io, io_summary, kdtree, kernel, linklist, mpiderivs, mpiutils, nicil,
!    options, part, physcon, ptmass, timestep, timestep_ind, timestep_sts,
!    timing, units, viscosity
!+
!--------------------------------------------------------------------------
module forces
 use dim, only:ndusttypes
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 integer, parameter :: maxcellcache = 50000

 public :: force

 !--indexing for xpartveci array
 integer, parameter :: xpartvecvars = 13 ! Number of xpartvec variables
 integer, parameter :: xpartvecarrs = 2  ! Number of variables that are arrays
 integer, parameter :: maxxpartveci = xpartvecvars + xpartvecarrs*(ndusttypes-1) ! Total number of values
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
       ipsi     = 11, &
       !--dust arrays initial index
       idustfraci = 12, &
       itstop     = 13 + (ndusttypes-1), &
       !--dust arrays final index
       idustfraciend = itstop-1, &
       itstopend     = maxxpartveci

 !--indexing for fsum array
 integer, parameter :: fsumvars = 17 ! Number of fsum variables
 integer, parameter :: fsumarrs = 5  ! Number of variables that are arrays
 integer, parameter :: maxfsum = fsumvars + fsumarrs*(ndusttypes-1) ! Total number of values
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
subroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,dustfrac,ddustfrac,&
                 ipart_rhomax,dt,stressmax)
 use dim,      only:maxp,ndivcurlv,ndivcurlB,maxvxyzu,maxalpha,maxneigh,maxstrain,&
                    switches_done_in_derivs,mhd,maxBevol,mhd_nonideal,use_dustfrac,lightcurve
 use eos,          only:use_entropy,gamma,equationofstate,get_temperature_from_ponrho
 use io,           only:iprint,fatal,iverbose,id,master,real4,warning,error
 use linklist,     only:ncells,ifirstincell,get_neighbour_list
 use options,      only:alpha,alphau,alphaB,beta, &
                        ipdv_heating,ishock_heating,iresistive_heating,icooling, &
                        psidecayfac,overcleanfac
 use part,         only:rhoh,dhdrho,rhoanddhdrho,massoftype,&
                        alphaind,abundance,nabundances, &
                        ll,get_partinfo,iactive,gradh,&
                        hrho,iphase,maxphase,igas,iboundary,maxgradh,straintensor, &
                        n_R,n_electronT,deltav
 use timestep,     only:dtcourant,dtforce,C_cour,C_force,C_cool,dtmax,bignumber,dtdiff
 use io_summary,   only:summary_variable, &
                        iosumdtf,iosumdtd,iosumdtv,iosumdtc,iosumdto,iosumdth,iosumdta, &
                        iosumdgs,iosumdge,iosumdgr,iosumdtfng,iosumdtdd,iosumdte
#ifdef FINVSQRT
 use fastmath,     only:finvsqrt
#endif
 use physcon,      only:pi
 use viscosity,    only:irealvisc,bulkvisc,shearfunc,dt_viscosity
#ifdef IND_TIMESTEPS
 use timestep_ind, only:nbinmax,ibinnow,get_newbin
 use timestep_sts, only:nbinmaxsts,sts_modify_ibin,sts_it_n
 use part,         only:ibin,ibinsink
 use timestep_sts, only:ibinsts
 use timestep,     only:nsteps,time
#endif
 use part,         only:divBsymm,Bextx,Bexty,Bextz,isdead_or_accreted,h2chemistry,ngradh,gravity,ibin_wake
 use timestep_sts, only:use_sts
#ifdef MPI
 use domain,       only:cellbelong
 use io,           only:nprocs
 use mpiutils,     only:fill_buffer
 use mpiderivs,    only:init_results_exchange,finish_results_exchange,&
                        send_results,recv_force_results,check_send_finished
 use timing,       only:getused
#endif
 use mpiutils,     only:reduce_mpi,reduceall_mpi
 use cooling,      only:energ_cooling
 use chem,         only:energ_h2cooling
#ifdef GRAVITY
 use kernel,       only:kernel_softening
 use kdtree,       only:get_distance_from_centre_of_mass,expand_fgrav_in_taylor_series
 use part,         only:xyzmh_ptmass,nptmass,poten
 use ptmass,       only:icreate_sinks,rho_crit,r_crit2
 use units,        only:unit_density
#endif
#ifdef DUST
 use dust,         only:get_ts,grainsize,graindens,idrag
#endif
#ifdef LIGHTCURVE
 use part,         only:luminosity
#endif
 use nicil,        only:nicil_get_eta,nimhd_get_jcbcb,nimhd_get_dt,nimhd_get_dBdt, &
                        nimhd_get_dudt,nicil_translate_error
 integer,      intent(in)    :: icall,npart
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(inout) :: vxyzu(:,:)
 real,         intent(in)    :: dustfrac(:,:)
 real,         intent(out)   :: fxyzu(:,:), ddustfrac(:,:)
 real(kind=4), intent(in)    :: Bevol(:,:)
 real(kind=4), intent(out)   :: dBevol(:,:)
 real(kind=4), intent(inout) :: divcurlv(:,:)
 real(kind=4), intent(in)    :: divcurlB(:,:)
 real,         intent(in)    :: dt,stressmax
 integer,      intent(out)   :: ipart_rhomax ! test this particle for point mass creation
 real :: xpartveci(maxxpartveci),tseff
 real :: fsum(maxfsum)
 real, save :: xyzcache(4,maxcellcache)
 integer, save :: listneigh(maxneigh)
!$omp threadprivate(xyzcache,listneigh)
 integer :: i,j,icell,nneigh
 integer :: nstokes,nsuper,ndrag,ndustres
 integer :: ierr
 real :: hi
 real(kind=8) :: hi1,hi21,hi31,hi41
 real :: dhdrhoi
 real :: rho1i,pro2i,pri,pmassi
 real :: vsigmax,dti,dtc,dtf,dtcool,dtdrag,dtdusti
 real :: dtmini,vsigdtc
 real :: drhodti,f2i,ponrhoi,spsoundi,vwavei,temperaturei
 real :: divvi,vcleani
 real :: divcurlvi(ndivcurlv)
 real :: divcurlBi(ndivcurlB)
 real :: dudtnonideal
 real :: dtvisci,dtohmi,dthalli,dtambii,dtdiffi
 real :: dtohm,dthall,dtambi,dtclean
 real :: jcbcbi(3),jcbi(3),etaohmi,etahalli,etaambii
 real :: source
 real :: sxxi,sxyi,sxzi,syyi,syzi,szzi,dtvisc
 real :: dustresfacmean
 real :: dustresfacmax
 real :: visctermiso,visctermaniso
 real :: alphai
 real :: straini(6)
 real :: rhoi,rhogasi,dustfraci(ndusttypes),dustfracisum,fac,csi
 real(kind=8)       :: gradhi,gradsofti
 real               :: Bxi,Byi,Bzi,B2i,Bi,Bi1,divBsymmi,psii,dtau,frac_divB,betai
 real :: pdv_work,dudt_radi,shearvisc,fxyz4
#ifdef GRAVITY
 real    :: potensoft0,dum,dx,dy,dz,fxi,fyi,fzi,poti,epoti
 real    :: fgrav(20)
 real    :: rhomax,rhomax_thread
 logical :: use_part
 integer :: ipart_rhomax_thread
#endif
#ifdef DUST
 real :: frac_stokes, frac_super
 integer :: iregime
#endif
#ifdef MPI
 integer, parameter :: maxrecv = 2*maxvxyzu + 1 + maxBevol + 1
 integer      :: nactive_thisproc,n,nactivetot
 real         :: xrecvbuf(maxrecv,nprocs),xsendbuf(maxrecv)
 integer      :: ireqrecv(nprocs),ireqsend(nprocs),nrecv
 real(kind=4) :: t1,t2,t3
 logical      :: have_sent
#endif
 integer :: iamtypei,idudtcool,ichem
 integer(kind=1) :: ibin_neigh
 logical :: ifilledcellcache,moreincell
 logical :: iactivei,iamdusti,iamgasi,realviscosity,useresistiveheat
#ifndef IND_TIMESTEPS
 !integer(kind=1), save :: ibin_wake(1)
 real    :: dtmaxi

 dtmaxi = 0.
#else
 !!integer(kind=1), save :: ibin_wake(maxp)
 integer(kind=1)       :: ibinnow_m1
 integer :: nbinmaxnew,nbinmaxstsnew,ncheckbin
 integer :: ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd
 integer :: ndtvisc,ndtohm,ndthall,ndtambi,ndtdust
 real    :: dtitmp,dtrat
 real    :: dtfrcfacmean ,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean
 real    :: dtfrcfacmax  ,dtfrcngfacmax ,dtdragfacmax ,dtdragdfacmax ,dtcoolfacmax
 real    :: dtviscfacmean,dtohmfacmean  ,dthallfacmean,dtambifacmean,dtdustfacmean
 real    :: dtviscfacmax ,dtohmfacmax   ,dthallfacmax ,dtambifacmax, dtdustfacmax
 logical :: allow_decrease,dtcheck,check_ibinsink
 character(len=16) :: dtchar

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
 if (all(ibinsink(1:npart)==0)) then
    check_ibinsink = .false.
 else
    check_ibinsink = .true.
endif
#endif

 dustresfacmean  = 0.0
 dustresfacmax   = 0.0
 dtcourant   = bignumber
 dtforce     = bignumber
 dtvisc      = bignumber
 dtmini      = bignumber
 dtohm       = bignumber
 dthall      = bignumber
 dtambi      = bignumber
 ibin_wake   = 0
 if (iverbose >= 3 .and. id==master) write(iprint,*) 'forces: cell cache =',maxcellcache

 realviscosity    = (irealvisc > 0)
 useresistiveheat = (iresistive_heating > 0)
 if (ndivcurlv < 1) call fatal('force','divv not stored but it needs to be')
 if (switches_done_in_derivs) call fatal('force','need switches_done_in_derivs=.false.')
 !--set stress terms to zero initially (make sure they are firstprivate)
 straini(:)    = 0.
 !--set alphai if not evolved (should be firstprivate)
 alphai        = alpha
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
 nactive_thisproc = 0
 nrecv = 0
 call getused(t1)
 have_sent = .false.
 call init_results_exchange(nprocs,xrecvbuf,ireqrecv)
#endif
!
!-- verification for non-ideal MHD
 if (mhd_nonideal .and. ndivcurlB < 4) call fatal('force','non-ideal MHD needs curl B stored, but ndivcurlB < 4')
!
!--check that compiled options are compatible with this routine
!
 if (maxgradh /= maxp) call fatal('force','need storage of gradh (maxgradh=maxp)')
!
!$omp parallel default(none) &
!$omp shared(ncells,ll,ifirstincell,npart,icall) &
!$omp shared(beta,gamma,icooling) &
!$omp shared(xyzh,vxyzu,fxyzu,divcurlv,massoftype,iphase,abundance,straintensor,deltav) &
!$omp shared(dt,dtmax,gradh,ishock_heating,ipdv_heating) &
!$omp shared(alphaind,alpha,alphau,alphaB) &
!$omp shared(irealvisc,realviscosity,useresistiveheat,bulkvisc,C_cour,C_force,C_cool,stressmax) &
!$omp shared(n_R,n_electronT,use_sts) &
!$omp firstprivate(straini,alphai) &
!$omp private(icell,i,j,iamtypei,iamgasi,ierr) &
!$omp private(ifilledcellcache,moreincell) &
!$omp private(hi,hi1,hi21,hi31,hi41,rhoi,rho1i,gradhi,gradsofti,pmassi) &
!$omp private(nneigh,idudtcool,ichem) &
!$omp private(dhdrhoi) &
!$omp private(ponrhoi,spsoundi,vwavei,pro2i,pri,temperaturei) &
!$omp private(vsigmax,dti,vsigdtc,dtcool,dtdrag) &
!$omp private(divcurlvi,divvi,source,drhodti) &
!$omp private(divcurlBi) &
!$omp private(jcbi,jcbcbi,etaohmi,etahalli,etaambii) &
!$omp private(dudtnonideal,dtohmi,dthalli,dtambii,dtdiffi) &
!$omp private(dtc,iactivei,iamdusti) &
!$omp private(sxxi,sxyi,sxzi,syyi,syzi,szzi) &
!$omp private(xpartveci,tseff,fsum) &
!$omp private(dtvisci,visctermiso,visctermaniso) &
!$omp private(pdv_work,dudt_radi,shearvisc,fxyz4) &
#ifdef LIGHTCURVE
!$omp shared(luminosity) &
#endif
#ifdef GRAVITY
!$omp private(fgrav,dx,dy,dz,poti,fxi,fyi,fzi,potensoft0,dum,epoti) &
!$omp shared(xyzmh_ptmass,nptmass,poten) &
!$omp shared(rhomax,ipart_rhomax,icreate_sinks,rho_crit,r_crit2) &
!$omp private(rhomax_thread,ipart_rhomax_thread,use_part) &
#endif
#ifdef MPI
!$omp shared(id,ireqrecv,ireqsend,xrecvbuf,xsendbuf,nprocs,have_sent) &
!$omp private(n) &
!$omp reduction(+:nactive_thisproc,nrecv) &
#endif
!$omp shared(Bevol,divBsymm,dBevol,Bextx,Bexty,Bextz,psidecayfac,overcleanfac) &
!$omp shared(divcurlB) &
!$omp private(dtau,vcleani,dtclean) &
!$omp private(divBsymmi,psii,Bxi,Byi,Bzi,B2i,Bi,Bi1,frac_divB,betai) &
#ifdef IND_TIMESTEPS
!$omp shared(ibin,ibinsts,ibinsink,nbinmax,nbinmaxsts,sts_it_n,check_ibinsink) &
!$omp private(allow_decrease,dtitmp,dtcheck,dtrat,dtchar) &
!$omp reduction(+:ndtforce,ndtforceng,ndtcool,ndtdrag,ndtdragd,ncheckbin,ndtvisc) &
!$omp reduction(+:ndtohm,ndthall,ndtambi,ndtdust,dtohmfacmean,dthallfacmean,dtambifacmean,dtdustfacmean) &
!$omp reduction(+:dtfrcfacmean,dtfrcngfacmean,dtdragfacmean,dtdragdfacmean,dtcoolfacmean,dtviscfacmean) &
!$omp reduction(max:dtohmfacmax,dthallfacmax,dtambifacmax,dtdustfacmax) &
!$omp reduction(max:dtfrcfacmax,dtfrcngfacmax,dtdragfacmax,dtdragdfacmax,dtcoolfacmax,dtviscfacmax) &
!$omp reduction(max:nbinmaxnew,nbinmaxstsnew) &
#else
!$omp reduction(min:dtohm,dthall,dtambi) &
!$omp reduction(min:dtcourant,dtforce,dtmini,dtvisc) &
!$omp reduction(max:dtmaxi) &
#endif
!$omp shared(ibin_wake) &
#ifdef DUST
!$omp shared(idrag,graindens,grainsize) &
!$omp private(iregime) &
#endif
!$omp shared(dustfrac,ddustfrac) &
!$omp private(rhogasi,dustfraci,dustfracisum,fac,csi,dtdusti,ibin_neigh) &
!$omp reduction(+:nstokes,nsuper,ndrag,ndustres,dustresfacmean) &
!$omp reduction(min:dtdiff) &
!$omp reduction(max:dustresfacmax) &
!$omp private(f2i,dtf)
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

    moreincell = (ll(i) /= 0)
    !
    !--get the neighbour list and fill the cell cache
    !
#ifdef GRAVITY
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,.false.,getj=.true.,f=fgrav)
#else
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,.false.,getj=.true.)
#endif
    ifilledcellcache = .true.

    over_parts: do while(i /= 0)

       if (i < 0) then  ! i < 0 indicates inactive
          i = ll(abs(i)) ! this is a quicker way of skipping inactives
          cycle over_parts ! but cannot be used for first in cell
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
          i = ll(i)
          cycle over_parts
       endif
       if (iamtypei==iboundary) then ! do not compute forces on boundary parts
          i = ll(i)
          cycle over_parts
       endif
       !
       !--fill the temporary array for this particle
       !
       pmassi = massoftype(iamtypei)
       xpartveci(:) = 0.
       xpartveci(ixi)  = xyzh(1,i)
       xpartveci(iyi)  = xyzh(2,i)
       xpartveci(izi)  = xyzh(3,i)
       hi              = xyzh(4,i)
       if (hi < 0.) call fatal('force','negative smoothing length',i,var='h',val=hi)
       xpartveci(ivxi) = vxyzu(1,i)
       xpartveci(ivyi) = vxyzu(2,i)
       xpartveci(ivzi) = vxyzu(3,i)
       if (maxvxyzu >= 4) xpartveci(ieni) = vxyzu(4,i)

       if (mhd .and. iamgasi) then
          Bxi = Bevol(1,i)
          Byi = Bevol(2,i)
          Bzi = Bevol(3,i)
          B2i = Bxi**2 + Byi**2 + Bzi**2
          xpartveci(iBevolxi) = Bxi
          xpartveci(iBevolyi) = Byi
          xpartveci(iBevolzi) = Bzi
          if (maxBevol >= 4) xpartveci(ipsi) = Bevol(4,i)
          if (maxBevol < 3 .or. maxBevol > 4) call fatal('force','error in maxBevol setting')
       endif

       !
       !--compute density and related quantities from the smoothing length
       !
       call rhoanddhdrho(hi,hi1,rhoi,rho1i,dhdrhoi,pmassi)
       hi21 = hi1*hi1
       hi31 = hi21*hi1
       hi41 = hi21*hi21
       !
       ! retrieve other quantities from their stored values
       !
       if (gradh(1,i) > 0.) then   ! gradh stored from previous timestep
          gradhi = gradh(1,i)       ! but check that it is non-zero...
       else
          call error('force','stored gradh is zero, resetting to 1')
          gradhi = 1.
       endif
#ifdef GRAVITY
       gradsofti = gradh(2,i)
#endif
       if (iamgasi) then
          if (realviscosity .and. maxstrain==maxp) straini(:) = straintensor(:,i)
          if (ndivcurlv >= 1) divcurlvi(:) = divcurlv(:,i)
          if (ndivcurlB >= 1) divcurlBi(:) = divcurlB(:,i)

          !
          ! one-fluid dust properties
          !
          if (use_dustfrac .and. iamgasi) then
             dustfraci(:) = dustfrac(:,i)
             dustfracisum = sum(dustfraci)
             rhogasi      = rhoi*(1.-dustfracisum)
             do j = 1,ndusttypes
                if (dustfraci(j) > 1. .or. dustfraci(j) < 0.) call fatal('force','invalid eps',var='dustfrac',val=dustfraci(j))
             enddo
          else
             dustfraci(:) = 0.
             rhogasi   = rhoi
          endif
          xpartveci(idustfraci:idustfraciend) = dustfraci(:)

          !
          ! calculate terms required in the force evaluation
          !
          call get_P(rhoi,rho1i,xpartveci(ixi),xpartveci(iyi),xpartveci(izi), &
                     pmassi,xpartveci(ieni),Bxi,Byi,Bzi,dustfraci(:),ponrhoi,pro2i,pri,spsoundi, &
                     vwavei,sxxi,sxyi,sxzi,syyi,syzi,szzi, &
                     visctermiso,visctermaniso,realviscosity,divcurlvi(1),bulkvisc,straini,stressmax)

#ifdef DUST
          !
          ! get stopping time - for one fluid dust we don't know deltav, but as small by definition we assume=0
          !
          if (use_dustfrac .and. iamgasi) then
             do j = 1,ndusttypes
                call get_ts(idrag,grainsize(j),graindens,rhogasi,rhoi*dustfracisum,spsoundi,0.,xpartveci(itstop+(j-1)),iregime)
             enddo
          endif
#endif
          !
          ! calculate terms required for non-ideal MHD
          !
          if (mhd_nonideal) then
             Bi  = sqrt( B2i )
             if (Bi > 0.0) then
               Bi1 = 1.0/Bi
             else
               Bi1 = 0.0
             endif
             temperaturei = get_temperature_from_ponrho(ponrhoi)
             call nicil_get_eta(etaohmi,etahalli,etaambii,Bi,rhoi,temperaturei &
                               ,n_R(:,i),n_electronT(i),ierr)
             if (ierr/=0) then
               call nicil_translate_error(ierr)
               call fatal('force','error in calcuating eta(i) for non-ideal MHD')
             endif
             call nimhd_get_jcbcb(jcbcbi,jcbi,divcurlBi(2:4),Bxi,Byi,Bzi,Bi1)
          endif
          !
          ! viscosity and resistivity switches
          !
          if (maxalpha==maxp)  alphai  = alphaind(1,i)
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
       endif
!
!--loop over current particle's neighbours (includes self)
!
       call compute_forces(i,iamgasi,iamdusti,xpartveci,hi,hi1,hi21,hi41,gradhi,gradsofti, &
                           pro2i,pri,spsoundi,vwavei,beta, &
                           visctermiso,visctermaniso,sxxi,sxyi,sxzi,syyi,syzi,szzi, &
                           pmassi,rhoi,rho1i,listneigh,nneigh,xyzcache,fsum,vsigmax, &
                           ifilledcellcache,realviscosity,useresistiveheat, &
                           xyzh,vxyzu,Bevol,iphase,massoftype, &
                           etaohmi,etahalli,etaambii,jcbcbi,jcbi,divcurlB,n_R,n_electronT, &
                           dustfrac,gradh,divcurlv,alphaind, &
                           alphai,alphau,alphaB,bulkvisc,stressmax,npart,&
                           ndrag,nstokes,nsuper,dtdrag,ibin_wake,ibin_neigh)

#ifdef GRAVITY
       !--add self-contribution
       call kernel_softening(0.,0.,potensoft0,dum)
       epoti = 0.5*pmassi*(fsum(ipot) + pmassi*potensoft0*hi1)
       !
       !--add contribution from distant nodes, expand these in Taylor series about node centre
       !
       call get_distance_from_centre_of_mass(icell,xpartveci(ixi),xpartveci(iyi),xpartveci(izi),dx,dy,dz)
       call expand_fgrav_in_taylor_series(fgrav,dx,dy,dz,fxi,fyi,fzi,poti)
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
          !  regimes
          !
          divBsymmi  = fsum(idivBsymi)
          if (B2i > 0.0) then
             betai = 2.0*ponrhoi*rhoi/B2i
             if (betai < 1.0) then
                frac_divB = 1.0
             elseif (betai < 2.0) then
                frac_divB = 2.0 - betai
             else
                frac_divB = 0.0
             endif
          else
             frac_divB = 0.0
          endif
          fsum(ifxi) = fsum(ifxi) - Bxi*divBsymmi*frac_divB
          fsum(ifyi) = fsum(ifyi) - Byi*divBsymmi*frac_divB
          fsum(ifzi) = fsum(ifzi) - Bzi*divBsymmi*frac_divB
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
       !--contributions to thermal energy/entropy equation

       divvi = -drhodti*rho1i
       if (ndivcurlv >= 1) divcurlv(1,i) = real(divvi,kind=kind(divcurlv)) ! store divv from forces

       if (maxvxyzu >= 4 .or. lightcurve) then
          if (maxstrain == maxp .and. realviscosity) then
             shearvisc = shearfunc(xpartveci(ixi),xpartveci(iyi),xpartveci(izi),spsoundi)
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
                call nimhd_get_dudt(dudtnonideal,etaohmi,etaambii,rhoi,real(divcurlBi(2:4)),real(Bevol(1:3,i)))
                fxyz4 = fxyz4 + fac*dudtnonideal
             endif
             !--add conductivity and resistive heating
             fxyz4 = fxyz4 + fac*fsum(idendtdissi)
             if (icooling > 0) then
                if (h2chemistry) then
                   idudtcool = 1
                   ichem = 0
                   call energ_h2cooling(vxyzu(4,i),fxyz4,rhoi,&
                        abundance(:,i),nabundances,dt,xyzh(1,i),xyzh(2,i),xyzh(3,i),idudtcool,ichem)
                else
                   call energ_cooling(icooling,vxyzu(4,i),fxyz4,xyzh(1,i),xyzh(2,i),xyzh(3,i))
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
          ! sum returns d(B/rho)/dt, convert this to dB/dt
          !
          dBevol(1,i) = real(rhoi*fsum(idBevolxi) - Bxi*divvi,kind=kind(dBevol))
          dBevol(2,i) = real(rhoi*fsum(idBevolyi) - Byi*divvi,kind=kind(dBevol))
          dBevol(3,i) = real(rhoi*fsum(idBevolzi) - Bzi*divvi,kind=kind(dBevol))

          !
          ! hyperbolic/parabolic cleaning terms (dpsi/dt) from Tricco & Price (2012)
          !
          if (maxBevol >= 4 .and. psidecayfac > 0.) then
             vcleani = overcleanfac*vwavei
             fsum(idivBdiffi) = fsum(idivBdiffi)*rho1i
             dtau = psidecayfac*vcleani*hi1
             !
             ! we clean using the difference operator for div B
             !
             psii = xpartveci(ipsi)

             ! original cleaning as in Tricco & Price (2012)
             !dBevol(4,i) = real(-vwavei*vwavei*fsum(idivBdiffi) - psii*dtau - 0.5*psii*divvi,kind=kind(dBevol))

             ! new cleaning evolving d/dt (psi/c_h)
             dBevol(4,i) = real(-vcleani*fsum(idivBdiffi) - psii*dtau - 0.5*psii*divvi,kind=kind(dBevol))

             ! timestep from cleaning (should only matter if overcleaning applied)
             ! the factor of 0.5 is empirical, from checking when overcleaning with ind. timesteps is stable
             dtclean = 0.5*C_cour*hi/(vcleani + epsilon(0.))
          endif
       endif

       if (use_dustfrac) then
          ddustfrac(:,i) = 0.5*(fsum(iddustfraci:iddustfraciend)-sqrt(rhoi*dustfraci(:))*divvi)
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
       if (maxvxyzu >= 4 .and. icooling==1) then
          dtcool = C_cool*abs(vxyzu(4,i)/fxyzu(4,i))
       else
          dtcool = bignumber
       endif

       ! timestep based on non-ideal MHD
       if (mhd_nonideal) then
          call nimhd_get_dt(dtohmi,dthalli,dtambii,hi,etaohmi,etahalli,etaambii)
          if ( use_sts ) then
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
       dtvisci = dt_viscosity(xpartveci(ixi),xpartveci(iyi),xpartveci(izi),hi,spsoundi)

       ! Check to ensure we have enough resolution for gas-dust pairs, where
       ! dtdrag is already minimised over all dust neighbours for gas particle i
       if (dtdrag < bignumber) then
          if (hi > dtdrag*spsoundi) then
             ndustres       = ndustres + 1
             dustresfacmean = dustresfacmean + hi/(dtdrag*spsoundi)
             dustresfacmax  = max(dustresfacmax, hi/(dtdrag*spsoundi))
          endif
       endif

else

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
          tseff = (1.-dustfracisum)/dustfracisum*sum(dustfraci(:)*xpartveci(itstop:itstopend))
          dtdusti = C_force*hi*hi/(dustfracisum*tseff*spsoundi**2)
       else
          dtdusti = bignumber
       endif

       !
#ifdef IND_TIMESTEPS
       !-- The new timestep for particle i
       dtitmp = min(dtf,dtcool,dtvisci,dtdrag,dtohmi,dthalli,dtambii,dtdusti)
       if (dtitmp < dti .and. dtitmp < dtmax) then
         dti     = dtitmp
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
       ibin(i) = max(ibin(i),ibin_neigh-1_1)
       !
       ! This will to keep particles near sink candidates awake (add 1 to smoothly move the particle)
       !
       if (check_ibinsink) then
          if (ibinsink(i) > ibin(i)) ibin(i) = ibin(i) + 1_1
          ! This is a more useful way to manage this array in the cases where
          ! multiple particles are being tested for sink creation
          ibinsink(i) = max(ibinsink(i) - 1_1, 0_1)
       endif
       !
       ! find the new maximum number of bins
       nbinmaxnew = max(nbinmaxnew,int(ibin(i)))
       ncheckbin  = ncheckbin + 1

       ! ibinsts: based entirely upon the diffusive timescale
       if ( use_sts ) then
          ibinsts(i) = 0 ! we actually want dtdiff, and this is just a tracer; should reduce the number of sts active particles for speed
          call get_newbin(dtdiffi,dtmax,ibinsts(i),allow_decrease,.false.)
          nbinmaxstsnew = max(nbinmaxstsnew,int(ibinsts(i)))
          ibinsts(i) = -ibinsts(i) ! set as negative to tag as active; will be reset to positive shortly
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

#ifdef MPI
       nactive_thisproc = nactive_thisproc + 1
!$omp critical(send)
       if (have_sent) then
          ! have to make sure last send completed before sending another
          call check_send_finished(nprocs,ireqsend,ireqrecv,xrecvbuf,nrecv)
       endif

       ! broadcast particle result to all procs
       n = 0
       call fill_buffer(xsendbuf,fxyzu(:,i),n)
#ifdef IND_TIMESTEPS
       call fill_buffer(xsendbuf,ibin(i),n)
#endif
       if (ndivcurlv >= 1) call fill_buffer(xsendbuf,divcurlv(1,i),n)
       if (mhd) then
          call fill_buffer(xsendbuf,divBsymm(i),n)
          call fill_buffer(xsendbuf,dBevol(:,i),n)
       endif
       !if (gravity) call fill_buffer(xsendbuf,poten(i),n)
       if (n > size(xsendbuf)) call fatal('force','mpi buffer size exceeded: increase ipartbufsize and recompile')
       call send_results(nprocs,xsendbuf,n,i,ireqsend)
       have_sent = .true.

       ! check for receives
       call recv_force_results(nprocs,xrecvbuf,ireqrecv,nrecv)
!$omp end critical(send)
#endif

       !--go to next particle in cell link list
       i = ll(i)

    enddo over_parts

 enddo over_cells

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
!$omp end parallel

#ifdef IND_TIMESTEPS
 ! check for nbinmaxnew = 0, can happen if all particles
 ! are dead/inactive, e.g. after sink creation
 if (ncheckbin==0) then
    nbinmaxnew    = nbinmax
    nbinmaxstsnew = nbinmaxsts
 endif
#endif

#ifdef GRAVITY
 if (icreate_sinks > 0 .and. ipart_rhomax > 0 .and. iverbose>=1) then
      print*,' got rhomax = ',rhomax*unit_density,' on particle ',ipart_rhomax !,rhoh(xyzh(4,ipart_rhomax))
 endif
#endif

#ifdef DUST
 ndrag = reduce_mpi('+',ndrag)
 if (ndrag > 0) then
    nstokes = reduce_mpi('+',nstokes)
    nsuper =  reduce_mpi('+',nsuper)
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

 nbinmaxsts = int(reduceall_mpi('max',nbinmaxstsnew),kind=1)
 nbinmax    = int(reduceall_mpi('max',nbinmaxnew),kind=1)
 ndtforce   = int(reduce_mpi('+',ndtforce))
 ndtforceng = int(reduce_mpi('+',ndtforceng))
 ndtcool    = int(reduce_mpi('+',ndtcool))
 ndtdrag    = int(reduce_mpi('+',ndtdrag))
 ndtdragd   = int(reduce_mpi('+',ndtdragd))
 ndtdust    = int(reduce_mpi('+',ndtdust))

 !  If use_sts, increase ibin as required for ibinsts > ibin
 if ( use_sts .and. sts_it_n ) call sts_modify_ibin(npart,ibin,nbinmax)

 !--Saitoh-Makino limiter
 !  currently switched OFF as not working with
 !  new individual timestepping routine (28/06/16 - DJP)
 if (.false.) then
 ibinnow_m1 = min(ibinnow - 1_1,nbinmax) ! min to prevent failures if nbinmax decreases
!$omp parallel default(none) &
!$omp shared(npart,ibin,ibin_wake,ibinnow_m1,iphase) &
!$omp private(i)
!$omp do schedule(runtime)
   do i=1,npart
      if (ibin_wake(i)==1) then
         if (ibin(i) < ibinnow_m1) then
            ibin(i) = ibinnow_m1
         endif
      endif
   enddo
!$omp enddo
!$omp end parallel
 endif

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
 dtmaxi    = reduce_mpi('max',dtmini)
 if (iverbose >= 2 .and. id==master) write(iprint,*) 'dtmin = ',C_Cour*dtmini, ' dtmax = ',C_cour*dtmaxi, &
    ' dtmax/dtmin = ',dtmaxi/(dtmini + epsilon(0.)),'dtcour/dtf = ',(C_cour*dtcourant)/(C_force*dtforce + epsilon(0.))
 if ( dtforce < dtcourant ) call summary_variable('dt',iosumdtf,0,0.0)
 if ( dtvisc  < dtcourant ) call summary_variable('dt',iosumdtv,0,0.0)
 if ( mhd_nonideal ) then
    ! Note: We are not distinguishing between use_sts and .not.use_sts here since if
    !       use_sts==.true., then dtohm=dtambi=bignumber.
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

#ifdef MPI
 nactivetot = int(reduceall_mpi('+',nactive_thisproc))
! if (nactivetot /= nactive) call fatal('force','MPI error: nactive/=nsent in force',var='nactive',ival=nactive)
!
!--check for remaining receives onto this proc.
!
 call getused(t2)
 if (iverbose >= 2) then
    write(iprint,*) ' thread ',id,' nactive = ',nactive_thisproc,' of ',nactivetot,npart,' time = ',t2-t1,'s'
 endif
 call finish_results_exchange(nprocs,xrecvbuf,ireqrecv,nrecv,nactivetot-nactive_thisproc)
 call getused(t3)

 if (nactivetot - nrecv < 0) call fatal('force','received more than nactive',var='nrecv',ival=nrecv)
 if (iverbose >= 2) then
    print "(1x,i2,a,f6.2,a,f6.2,a)", &
      id,' time spent waiting for receives = ',t3-t2,'s = ',100.*(t3-t2)/((t3-t1) + epsilon(t1)),'%'
 endif
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
                           pro2i,pri,spsoundi,vwavei,beta, &
                           visctermiso,visctermaniso,sxxi,sxyi,sxzi,syyi,syzi,szzi, &
                           pmassi,rhoi,rho1i,listneigh,nneigh,xyzcache,fsum,vsigmax, &
                           ifilledcellcache,realviscosity,useresistiveheat, &
                           xyzh,vxyzu,Bevol,iphase,massoftype, &
                           etaohmi,etahalli,etaambii,jcbcbi,jcbi,divcurlB,n_R,n_electronT, &
                           dustfrac,gradh,divcurlv,alphaind, &
                           alphai,alphau,alphaB,bulkvisc,stressmax,npart,&
                           ndrag,nstokes,nsuper,ts_min,ibin_wake,ibin_neigh)
#ifdef FINVSQRT
  use fastmath,    only:finvsqrt
#endif
  use kernel,      only:grkern,cnormk,radkern2
  use part,        only:igas,maxphase,iactive,iamtype,iamdust,idust,get_partinfo,iboundary
  use part,        only:mhd,maxvxyzu,maxBevol,maxstrain
  use dim,         only:maxalpha,maxp,mhd_nonideal,gravity,use_dust,use_dustfrac,lightcurve
  use part,        only:rhoh,maxgradh,straintensor
  use eos,         only:get_temperature_from_ponrho
  use nicil,       only:nicil_get_eta,nimhd_get_jcbcb,nimhd_get_dBdt
#ifdef GRAVITY
  use kernel,      only:kernel_softening
  use part,        only:nptmass,xyzmh_ptmass,ihacc
  use ptmass,      only:ptmass_not_obscured
#endif
#ifdef PERIODIC
  use boundary,    only:dxbound,dybound,dzbound
#endif
#ifdef DUST
  use dust,        only:get_ts,grainsize,graindens,idrag,icut_backreaction
  use kernel,      only:wkern_drag,cnormk_drag
#endif
#ifdef IND_TIMESTEPS
  use part,        only:ibinold
#endif
  use timestep,    only:bignumber
  use options,     only:overcleanfac
  integer,         intent(in)  :: i
  logical,         intent(in)  :: iamgasi,iamdusti
  real,            intent(in)  :: xpartveci(:)
  real(kind=8),    intent(in)  :: hi1,hi21,hi41,gradhi,gradsofti
  real,            intent(in)  :: hi,pro2i,pri,spsoundi,vwavei,beta
  real,            intent(in)  :: visctermiso,visctermaniso
  real,            intent(in)  :: sxxi,sxyi,sxzi,syyi,syzi,szzi
  real,            intent(in)  :: pmassi,rhoi,rho1i
  integer,         intent(in)  :: listneigh(:)
  integer,         intent(in)  :: nneigh
  real,            intent(in)  :: xyzcache(:,:)
  real,            intent(out) :: fsum(maxfsum)
  real,            intent(out) :: vsigmax
  logical,         intent(in)  :: ifilledcellcache
  logical,         intent(in)  :: realviscosity,useresistiveheat
  real,            intent(in)  :: xyzh(:,:),vxyzu(:,:)
  real(kind=4),    intent(in)  :: Bevol(:,:)
  real(kind=4),    intent(in)  :: divcurlB(:,:)
  real,            intent(in)  :: dustfrac(:,:)
  integer(kind=1), intent(in)  :: iphase(:)
  real,            intent(in)  :: massoftype(:)
  real,            intent(in)  :: jcbcbi(:),jcbi(:),n_R(:,:),n_electronT(:)
  real,            intent(in)  :: etaohmi,etahalli,etaambii
  real(kind=4),    intent(in)  :: alphaind(:,:)
  real(kind=4),    intent(in)  :: gradh(:,:),divcurlv(:,:)
  real,            intent(in)  :: alphai,alphau,alphaB,bulkvisc,stressmax
  integer,         intent(in)  :: npart
  integer, intent(inout) :: ndrag,nstokes,nsuper
  real,            intent(out) :: ts_min
  integer(kind=1), intent(out) :: ibin_wake(:),ibin_neigh
  integer         :: j,l,n,iamtypej,ierr
  logical         :: iactivej,iamgasj,iamdustj
  real :: rij2,q2i,qi,xj,yj,zj,dx,dy,dz,runix,runiy,runiz,rij1,hfacgrkern
  real :: grkerni,grgrkerni,dvx,dvy,dvz,projv,denij,vsigi,vsigu,dudtdissi
  real :: projBi,projBj,dBx,dBy,dBz,dB2,projdB
  real :: dendissterm,dBdissterm,dudtresist,dpsiterm,pmassonrhoi
  real :: gradpi,projsxi,projsyi,projszi
  real :: gradp,projsx,projsy,projsz,Bxj,Byj,Bzj,Bj,Bj1,psij
  real :: dpsitermj,grkernj,grgrkernj,autermj,avBtermj,vsigj,spsoundj
  real :: gradpj,pro2j,projsxj,projsyj,projszj,sxxj,sxyj,sxzj,syyj,syzj,szzj,psitermj,dBrhoterm
  real :: visctermisoj,visctermanisoj,enj,hj,mrhoj5,alphaj,pmassj,rho1j
  real :: rhoj,ponrhoj,prj,rhoav1
  real :: hj1,hj21,q2j,qj,vwavej,divvj
  real :: strainj(6)
  real :: dustfracisum,rhogasj
#ifdef GRAVITY
  integer :: k
  real    :: fmi,fmj,dsofti,dsoftj
  real    :: xkpt,ykpt,zkpt,vpos
  logical :: add_contribution
#else
  logical, parameter :: add_contribution = .true.
#endif
  real :: phi,phii,phij,fgrav,fgravi,fgravj,termi
#ifdef DUST
  integer :: iregime
  real    :: dragterm,dragheating,wdrag,tsij(ndusttypes),dv2
  real    :: Dav(ndusttypes),grkernav,tsj(ndusttypes),dustfracterms(ndusttypes)
  real    :: rhogas1i,rhogas1j
#endif
  real :: dBevolx,dBevoly,dBevolz,divBsymmterm,divBdiffterm
  real :: rho21i,rho21j,Bxi,Byi,Bzi,psii,pmjrho21grkerni,pmjrho21grkernj
  real :: auterm,avBterm,mrhoi5,vsigB
  real :: etaohmj,etahallj,etaambij,temperaturej
  real :: jcbcbj(3),jcbj(3),dBnonideal(3),dBnonidealj(3)
  real :: vsigavi,vsigavj,term
  real :: dustfraci(ndusttypes),dustfracj(ndusttypes),tsi(ndusttypes),sqrtrhodustfraci(ndusttypes),sqrtrhodustfracj(ndusttypes) !,vsigeps,depsdissterm
  real :: dustfracjsum,epstsi,epstsj
  logical :: usej

  fsum(:) = 0.
  vsigmax = 0.
  pmassonrhoi = pmassi*rho1i
  hfacgrkern  = hi41*cnormk*gradhi

  ! default settings for active/phase if iphase not used
  iactivej = .true.
  iamtypej = igas
  iamgasj  = .true.
  iamdustj = .false.
  ibin_neigh = 0_1

  ! dust
  ts_min = bignumber

  ! various pre-calculated quantities
  Bxi  = xpartveci(iBevolxi)
  Byi  = xpartveci(iBevolyi)
  Bzi  = xpartveci(iBevolzi)
  psii = xpartveci(ipsi)
  if (use_dustfrac) then
     dustfraci(:) = xpartveci(idustfraci:idustfraciend)
     dustfracisum = sum(dustfraci(:))
     tsi(:)       = xpartveci(itstop:itstopend)
     epstsi       = sum(dustfraci*tsi)
     sqrtrhodustfraci(:) = sqrt(rhoi*dustfraci(:))
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
  visctermisoj = 0.
  visctermanisoj = 0.

  loop_over_neighbours2: do n = 1,nneigh

     j = abs(listneigh(n))
     if (i==j) cycle loop_over_neighbours2

     if (ifilledcellcache .and. n <= maxcellcache) then
        ! positions from cache are already mod boundary
        xj = xyzcache(1,n)
        yj = xyzcache(2,n)
        zj = xyzcache(3,n)
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
#ifdef PERIODIC
        if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
        if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
        if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
     endif
     rij2 = dx*dx + dy*dy + dz*dz
     q2i = rij2*hi21

     !--hj is in the cell cache but not in the neighbour cache
     !  as not accessed during the density summation
     if (ifilledcellcache .and. n <= maxcellcache) then
        hj1 = xyzcache(4,n)
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
        k = 1
        do while (k <= nptmass .and. add_contribution)
           xkpt = xyzmh_ptmass(1,k)
           ykpt = xyzmh_ptmass(2,k)
           zkpt = xyzmh_ptmass(3,k)
           vpos = (xkpt-xpartveci(ixi))*(xkpt-xj) &
                + (ykpt-xpartveci(iyi))*(ykpt-yj) &
                + (zkpt-xpartveci(izi))*(zkpt-zj)
           if (vpos < 0.0) then
              add_contribution = ptmass_not_obscured(-dx,-dy,-dz,  &
                                 xkpt-xpartveci(ixi),ykpt-xpartveci(iyi),zkpt-xpartveci(izi), &
                                 xyzmh_ptmass(ihacc,k))
           endif
           k = k + 1
        enddo
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
        if (use_dust .or. use_dustfrac) usej = .true.
        if (maxvxyzu >= 4 .and. .not.gravity) usej = .true.

        !--get individual timestep/ multiphase information (querying iphase)
        if (maxphase==maxp) then
           call get_partinfo(iphase(j),iactivej,iamdustj,iamtypej)
           iamgasj = (iamtypej==igas .or. iamtypej==iboundary)
#ifdef IND_TIMESTEPS
           ! Particle j is a neighbour of an active particle;
           ! flag it to see if it needs to be woken up next step.
           if (iamtypej/=iboundary) then
              ibin_wake(j) = 1
              ibin_neigh = max(ibin_neigh,ibinold(j))
           endif
#endif
        endif
        pmassj = massoftype(iamtypej)

        fgrav = 0.5*pmassj*(fgravi + fgravj)

        !--If particle is hidden by the sink, treat the neighbour as
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
           denij = xpartveci(ieni) - enj
        else
           denij = 0.
        endif
        if (iamgasi .and. iamgasj) then
           !--work out vsig for timestepping and av
           vsigi   = max(vwavei - beta*projv,0.)
           vsigavi = max(alphai*vwavei - beta*projv,0.)
           if (vsigi > vsigmax) vsigmax = vsigi

           if (mhd) then
              Bxj = Bevol(1,j)
              Byj = Bevol(2,j)
              Bzj = Bevol(3,j)

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
                 sqrtrhodustfracj(:) = sqrt(rhoj*dustfracj(:))
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
              call get_P(rhoj,rho1j,xj,yj,zj,pmassj,enj,Bxj,Byj,Bzj,dustfracj, &
                         ponrhoj,pro2j,prj,spsoundj,vwavej, &
                         sxxj,sxyj,sxzj,syyj,syzj,szzj,visctermisoj,visctermanisoj, &
                         realviscosity,divvj,bulkvisc,strainj,stressmax)

              mrhoj5   = 0.5*pmassj*rho1j
              autermj  = mrhoj5*alphau
              avBtermj = mrhoj5*alphaB*rho1j

              vsigj = max(vwavej - beta*projv,0.)
              vsigavj = max(alphaj*vwavej - beta*projv,0.)
              if (vsigj > vsigmax) vsigmax = vsigj
           else
              vsigj = max(-projv,0.)
              if (vsigj > vsigmax) vsigmax = vsigj
              vsigavj = 0.
           endif
        else ! set to zero terms which are used below without an if (usej)
           rhoj      = 0.
           rho1j     = 0.
           rho21j    = 0.

           mrhoj5    = 0.
           autermj   = 0.
           avBtermj  = 0.

           gradpj    = 0.
           projsxj   = 0.
           projsyj   = 0.
           projszj   = 0.
           projBj = 0.
           vwavej = 0.
           vsigavj = 0.
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
           !--artificial resistivity
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
           !--grad psi term for divergence cleaning
           ! original cleaning as in Tricco & Price (2012)
           !if (maxBevol >= 4) dpsiterm = pmjrho21grkerni*psii + pmjrho21grkernj*psij

           ! new cleaning evolving d/dt (psi/c_h) as in Tricco, Price & Bate (2016)
           if (maxBevol >= 4) dpsiterm = overcleanfac*(pmjrho21grkerni*psii*vwavei + pmjrho21grkernj*psij*vwavej)
           !
           ! non-ideal MHD terms
           if (mhd_nonideal) then
              call nimhd_get_dBdt(dBnonideal,etaohmi,etahalli,etaambii,real(divcurlB(2:4,i)) &
                                 ,jcbi,jcbcbi,runix,runiy,runiz)
              dBnonideal = dBnonideal*pmjrho21grkerni
              if (usej) then
                 Bj  = sqrt(Bxj**2 + Byj**2 + Bzj**2)
                 if (Bj > 0.0) then
                   Bj1 = 1.0/Bj
                 else
                   Bj1 = 0.0
                 endif
                temperaturej = get_temperature_from_ponrho(ponrhoj)
                call nicil_get_eta(etaohmj,etahallj,etaambij,Bj,rhoj,temperaturej &
                                  ,n_R(:,j),n_electronT(j),ierr)
                call nimhd_get_jcbcb(jcbcbj,jcbj,real(divcurlB(2:4,j)),Bxj,Byj,Bzj,Bj1)
                call nimhd_get_dBdt(dBnonidealj,etaohmj,etahallj,etaambij,real(divcurlB(2:4,j)),jcbj,jcbcbj,runix,runiy,runiz)
                dBnonideal = dBnonideal + dBnonidealj*pmjrho21grkernj
              endif
           endif
           !
           ! dB/dt evolution equation
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

        if (maxvxyzu >= 4 .or. lightcurve) then
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
           epstsj   = sum(dustfracj*tsj)
           rhogas1i = rho1i/(1.-dustfracisum)
           rhogas1j = 1./rhogasj
           
           ! Check that weighted sums of Tsj and tilde(Tsj) are equal (see Hutchison et al. 2017)
           if (ndusttypes>1) then
              if (abs(sum(dustfracj*tsj) - sum(dustfracj*(tsj-epstsj))/(1. - dustfracjsum)) > 1e-14) print*,'Stop! tsj or epstsj in force is incorrect!'
           endif

           do l = 1,ndusttypes
              if (dustfraci(l) > 0. .or. dustfracj(l) > 0.) then
                 ! define averages of diffusion coefficient and kernels
                 !Dav(l)   = dustfraci(l)*tsi(l) + dustfracj(l)*tsj(l)
                 grkernav = 0.5*(grkerni + grkernj)

                 ! these are equations (43) and (45) from Price & Laibe (2015)
                 ! but note there is a sign error in the term in eqn (45) in the paper
                 !dustfracterm(l)  = pmassj*rho1j*Dav(:)*(pri - prj)*grkernav*rij1
                 dustfracterms(l) = pmassj*sqrtrhodustfracj(l)*rho1j*((tsi(l)-epstsi)*rhogas1i+(tsj(l)-epstsj)*rhogas1j)*(pri - prj)*grkernav*rij1
                 
                 !vsigeps = 0.5*(spsoundi + spsoundj) !abs(projv)
                 !depsdissterm = pmassj*sqrtrhodustfracj*rho1j*grkernav*vsigeps !(auterm*grkerni + autermj*grkernj)*vsigeps
                 !dustfracterms = dustfracterms - depsdissterm*(dustfraci - dustfracj)
                 fsum(iddustfraci+(l-1)) = fsum(iddustfraci+(l-1)) - dustfracterms(l)
                 !fsum(iddustfraci+(l-1)) = fsum(iddustfraci+(l-1)) - dustfracterm(l)
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
        if (idrag>0 .and. add_contribution) then
           if (iamgasi .and. iamdustj .and. icut_backreaction==0) then
              dv2 = dvx*dvx + dvy*dvy + dvz*dvz
              if (q2i < q2j) then
                 wdrag = wkern_drag(q2i,qi)*hi21*hi1*cnormk_drag
              else
                 wdrag = wkern_drag(q2j,qj)*hj21*hj1*cnormk_drag
              endif
              do l = 1,ndusttypes
                 call get_ts(idrag,grainsize(l),graindens,rhoi,rhoj,spsoundi,dv2,tsij(l),iregime)
              enddo
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
              do l = 1,ndusttypes
                 call get_ts(idrag,grainsize(l),graindens,rhoj,rhoi,spsoundj,dv2,tsij(l),iregime)
              enddo
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
subroutine get_P(rhoi,rho1i,xi,yi,zi,pmassi,eni,Bxi,Byi,Bzi,dustfraci, &
                 ponrhoi,pro2i,pri,spsoundi,vwavei, &
                 sxxi,sxyi,sxzi,syyi,syzi,szzi,visctermiso,visctermaniso, &
                 realviscosity,divvi,bulkvisc,strain,stressmax)

  use dim,       only:maxvxyzu,maxstrain,maxp
  use part,      only:mhd
  use eos,       only:equationofstate
  use options,   only:ieos
  use viscosity, only:shearfunc
  real,    intent(in)  :: rhoi,rho1i,xi,yi,zi,pmassi,eni
  real,    intent(in)  :: Bxi,Byi,Bzi,dustfraci(:)
  real,    intent(out) :: ponrhoi,pro2i,pri,spsoundi,vwavei
  real,    intent(out) :: sxxi,sxyi,sxzi,syyi,syzi,szzi
  real,    intent(out) :: visctermiso,visctermaniso
  logical, intent(in)  :: realviscosity
  real,    intent(in)  :: divvi,bulkvisc,stressmax
  real,    intent(in)  :: strain(6)

  real :: Bro2i,Brhoxi,Brhoyi,Brhozi,rhogasi,gasfrac
  real :: stressiso,term,graddivvcoeff,del2vcoeff
  real :: shearvisc,etavisc,valfven2i,p_on_rhogas
!
!--get pressure (actually pr/dens) and sound speed from equation of state
!
  gasfrac = (1. - sum(dustfraci))  ! rhogas/rho
  rhogasi = rhoi*gasfrac           ! rhogas = (1-eps)*rho
  if (maxvxyzu >= 4) then
     call equationofstate(ieos,p_on_rhogas,spsoundi,rhogasi,xi,yi,zi,eni)
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
  integer,          intent(inout) :: ndtopt
  real,             intent(in)    :: dti,dtopt,dtrat
  real,             intent(inout) :: dtoptfacmean,dtoptfacmax
  logical,          intent(inout) :: dtcheck
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
end module forces
