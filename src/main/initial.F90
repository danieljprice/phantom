!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module initial
!
! This module initialises (and ends) the run
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: HIIRegion, analysis, apr, boundary, boundary_dyn,
!   centreofmass, checkconserved, checkoptions, checksetup, cons2prim,
!   cooling, cpuinfo, damping, densityforce, deriv, dim, dust,
!   dust_formation, einsteintk_utils, energies, eos, evwrite, extern_gr,
!   externalforces, fastmath, fileutils, forcing, growth, inject, io,
!   io_summary, krome_interface, linklist, metric, metric_et_utils,
!   metric_tools, mf_write, mpibalance, mpidomain, mpimemory, mpitree,
!   mpiutils, nicil, nicil_sup, omputils, options, part, partinject,
!   porosity, ptmass, radiation_utils, readwrite_dumps, readwrite_infile,
!   subgroup, substepping, timestep, timestep_ind, timestep_sts, timing,
!   tmunu2grid, units, writeheader
!

 implicit none
 public :: initialise,finalise,startrun,endrun
 real(kind=4), private :: twall_start, tcpu_start

 private

contains

!----------------------------------------------------------------
!+
!  short initialisation routine that should be called
!  by any utility which will subsequently call derivs
!+
!----------------------------------------------------------------
subroutine initialise()
 use dim,              only:mpi,gr
 use io,               only:fatal,die,id,master,nprocs,ievfile
#ifdef FINVSQRT
 use fastmath,         only:testsqrt
#endif
 use omputils,         only:init_omp,info_omp
 use options,          only:set_default_options
 use io_summary,       only:summary_initialise
 use boundary,         only:set_boundary
 use writeheader,      only:write_codeinfo
 use evwrite,          only:init_evfile
 use mpidomain,        only:init_domains
 use cpuinfo,          only:print_cpuinfo
 use checkoptions,     only:check_compile_time_settings
 use readwrite_dumps,  only:init_readwrite_dumps
 use metric,           only:metric_type
 use metric_et_utils,  only:read_tabulated_metric,gridinit
 integer :: ierr

!
!--write 'PHANTOM' and code version
!
 if (id==master) call write_codeinfo(6)
!
!--check that it is OK to use fast sqrt functions
!  on this architecture
!
#ifdef FINVSQRT
 if (id==master) write(*,"(1x,a)") 'checking fast inverse sqrt...'
 call testsqrt(ierr,(id==master))
 if (ierr /= 0) call die
 if (id==master) write(*,"(1x,a,/)") 'done'
#else
 if (id==master) write(*,"(1x,a)") 'Using NATIVE inverse sqrt'
#endif

!
!--set default options (incl. units)
!
 call set_default_options
 call set_boundary
 call init_evfile(ievfile,'testlog',.false.)
!
!--initialise values for summary array
!
 call summary_initialise
!
!--check compile-time settings are OK
!
 call check_compile_time_settings(ierr)
 if (ierr /= 0) call fatal('initialise','incompatible compile-time settings')
!
!--initialise openMP things if required
!
 if (id==master) call print_cpuinfo()
 if (id==master) call info_omp
 call init_omp
!
!--initialise MPI domains
!
 call init_domains(nprocs)
!
!--initialise metric if tabulated
!
 if (gr  .and. metric_type=='et') then
    call read_tabulated_metric('tabuled_metric.dat',ierr)
    if (ierr == 0) gridinit = .true.
 endif

 call init_readwrite_dumps()

end subroutine initialise

!----------------------------------------------------------------
!+
!  routine which starts a Phantom run
!+
!----------------------------------------------------------------
subroutine startrun(infile,logfile,evfile,dumpfile,noread)
 use mpiutils,         only:reduceall_mpi,barrier_mpi,reduce_in_place_mpi
 use dim,              only:maxp,maxalpha,maxvxyzu,maxptmass,maxdusttypes,itau_alloc,itauL_alloc,&
                            nalpha,mhd,mhd_nonideal,do_radiation,gravity,use_dust,mpi,do_nucleation,&
                            use_dustgrowth,ind_timesteps,idumpfile,update_muGamma,use_apr,use_sinktree,gr,&
                            maxpsph
 use deriv,            only:derivs
 use evwrite,          only:init_evfile,write_evfile,write_evlog
 use energies,         only:compute_energies
 use io,               only:idisk1,iprint,ievfile,error,iwritein,flush_warnings,&
                            die,fatal,id,master,nprocs,real4,warning,iverbose
 use externalforces,   only:externalforce,initialise_externalforces,update_externalforce,&
                            externalforce_vdependent
 use options,          only:iexternalforce,icooling,use_dustfrac,rhofinal1,rhofinal_cgs
 use readwrite_infile, only:read_infile,write_infile
 use readwrite_dumps,  only:read_dump,write_fulldump
 use part,             only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,tau, tau_lucy, &
                            npartoftype,maxtypes,ndusttypes,alphaind,ntot,ndim,update_npartoftypetot,&
                            maxphase,iphase,isetphase,iamtype,igas,idust,imu,igamma,massoftype, &
                            nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,fxyz_ptmass_sinksink,&
                            epot_sinksink,get_ntypes,isdead_or_accreted,dustfrac,ddustevol,&
                            nden_nimhd,dustevol,rhoh,gradh,apr_level,aprmassoftype,&
                            Bevol,Bxyz,dustprop,filfac,ddustprop,ndustsmall,iboundary,eos_vars,dvdx, &
                            n_group,n_ingroup,n_sing,nmatrix,group_info,bin_info,isionised,shortsinktree,&
                            fxyz_ptmass_tree
 use part,             only:pxyzu,dens,metrics,rad,radprop,drad,ithick
 use densityforce,     only:densityiterate
 use linklist,         only:set_linklist
 use boundary_dyn,     only:dynamic_bdy,init_dynamic_bdy
 use substepping,      only:combine_forces_gr
#ifdef GR
 use part,             only:metricderivs,metricderivs_ptmass,metrics_ptmass,pxyzu_ptmass
 use cons2prim,        only:prim2consall
 use eos,              only:ieos
 use extern_gr,        only:get_grforce_all,get_tmunu_all,get_tmunu_all_exact
 use metric_tools,     only:init_metric,imet_minkowski,imetric
 use einsteintk_utils
 use tmunu2grid
#endif
 use units,            only:utime,umass,unit_Bfield
 use eos,              only:gmw,gamma
 use nicil,            only:nicil_initialise
 use nicil_sup,        only:use_consistent_gmw
 use ptmass,           only:init_ptmass,get_accel_sink_gas,get_accel_sink_sink, &
                            h_acc,r_crit,r_crit2,rho_crit,rho_crit_cgs,icreate_sinks, &
                            r_merge_uncond,r_merge_cond,r_merge_uncond2,r_merge_cond2,r_merge2, &
                            use_regnbody
 use timestep,         only:time,dt,dtextforce,C_force,dtmax,dtmax_user,idtmax_n
 use timing,           only:get_timings
 use timestep_ind,     only:ibinnow,maxbins,init_ibin,istepfrac
 use timing,           only:get_timings
 use part,             only:ibin,ibin_old,ibin_wake,alphaind
 use readwrite_dumps,  only:dt_read_in
 use timestep,         only:dtcourant,dtforce
#ifdef STS_TIMESTEPS
 use timestep,         only:dtdiff
#endif
 use timestep_sts,     only:sts_initialise
#ifdef DRIVING
 use forcing,          only:init_forcing
#endif
 use dust,             only:init_drag
 use growth,           only:init_growth
 use porosity,         only:init_porosity,init_filfac
 use options,          only:use_porosity
#ifdef MFLOW
 use mf_write,         only:mflow_write,mflow_init
 use io,               only:imflow
#endif
#ifdef VMFLOW
 use mf_write,         only:vmflow_write,vmflow_init
 use io,               only:ivmflow
#endif
#ifdef BINPOS
 use mf_write,         only:binpos_write,binpos_init
 use io,               only:ibinpos,igpos
#endif
 use dust_formation,   only:init_nucleation,set_abundances
#ifdef INJECT_PARTICLES
 use inject,           only:init_inject,inject_particles
 use partinject,       only:update_injected_particles
 use timestep_ind,     only:nbinmax
#endif
 use apr,             only:init_apr
#ifdef KROME
 use krome_interface,  only:initialise_krome
#endif
#ifdef LIVE_ANALYSIS
 use analysis,         only:do_analysis
 use part,             only:igas
 use fileutils,        only:numfromfile
 use io,               only:ianalysis
 use radiation_utils,  only:set_radiation_and_gas_temperature_equal
#endif
 use mpibalance,       only:balancedomains
 use part,             only:ibelong
 use writeheader,      only:write_codeinfo,write_header
 use eos,              only:ieos,init_eos
 use checksetup,       only:check_setup
 use cooling,          only:init_cooling
 use cpuinfo,          only:print_cpuinfo
 use units,            only:udist,unit_density
 use centreofmass,     only:get_centreofmass
 use energies,         only:etot,angtot,totmom,mdust,xyzcom,mtot
 use checkconserved,   only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in,mtot_in
 use fileutils,        only:make_tags_unique
 use damping,          only:idamp
 use subgroup,         only:group_identify,init_subgroup,update_kappa
 use HIIRegion,        only:iH2R,initialize_H2R,update_ionrates
 character(len=*), intent(in)  :: infile
 character(len=*), intent(out) :: logfile,evfile,dumpfile
 logical,          intent(in), optional :: noread
 integer         :: ierr,i,j,nerr,nwarn,ialphaloc,irestart,merge_n,merge_ij(maxptmass),boundi,boundf
 real            :: poti,hfactfile
 real            :: hi,pmassi,rhoi1
 real            :: dtsinkgas,dtsinksink,fonrmax,dtphi2,dtnew_first,dtinject
 real            :: stressmax,xmin,ymin,zmin,xmax,ymax,zmax,dx,dy,dz,tolu,toll
 real            :: dummy(3)
 real            :: gmw_nicil
#ifndef GR
 real            :: dtf,fextv(3)
#endif
 integer         :: itype,iposinit,ipostmp,ntypes,nderivinit
 logical         :: iexist,read_input_files
 character(len=len(dumpfile)) :: dumpfileold
 character(len=7) :: dust_label(maxdusttypes)
#ifdef INJECT_PARTICLES
 character(len=len(dumpfile)) :: file1D
 integer :: npart_old
#endif


 read_input_files = .true.
 if (present(noread)) read_input_files = .not.noread

 if (read_input_files) then
!
!--do preliminary initialisation
!
    call initialise
!
!--read parameters from the infile
!
    call read_infile(infile,logfile,evfile,dumpfile)
!
!--initialise log output
!
    if (iprint /= 6 .and. id==master) then
       open(unit=iprint,file=logfile,form='formatted',status='replace')
!
!--write opening "splash screen" to logfile
!
       call write_codeinfo(iprint)
       call print_cpuinfo(iprint)
    endif
    if (id==master) write(iprint,"(a)") ' starting run '//trim(infile)

    if (id==master) call write_header(1,infile,evfile,logfile,dumpfile)
!
!--read particle setup from dumpfile
!
    call read_dump(trim(dumpfile),time,hfactfile,idisk1,iprint,id,nprocs,ierr)
    if (ierr /= 0) call fatal('initial','error reading dumpfile')
    call check_setup(nerr,nwarn,restart=.true.) ! sanity check what has been read from file
    if (nwarn > 0) then
       print "(a)"
       call warning('initial','WARNINGS from particle data in file',var='# of warnings',ival=nwarn)
    endif
    if (nerr > 0)  call fatal('initial','errors in particle data from file',var='errors',ival=nerr)
!
!--if starting from a restart dump, rename the dumpefile to that of the previous non-restart dump
!
    irestart = index(dumpfile,'.restart')
    if (irestart > 0) write(dumpfile,'(2a,I5.5)') dumpfile(:irestart-1),'_',idumpfile
 endif
!
!--reset dtmax (required only to permit restart dumps)
!
 dtmax_user = dtmax           ! the user defined dtmax
 if (idtmax_n < 1) idtmax_n = 1
 dtmax      = dtmax/idtmax_n  ! dtmax required to satisfy the walltime constraints
!
!--Initialise dynamic boundaries in the first instance
!
 if (dynamic_bdy) call init_dynamic_bdy(1,npart,nptmass,dtmax)
!
!--initialise values for non-ideal MHD
!
 if (mhd_nonideal) then
    call nicil_initialise(real(utime),real(udist),real(umass),real(unit_Bfield),ierr,iprint,iprint)
    if (ierr/=0) call fatal('initial','error initialising nicil (the non-ideal MHD library)')

    call use_consistent_gmw(ierr,gmw,gmw_nicil)
    if (ierr/=0) write(iprint,'(2(a,Es18.7))') &
       ' initial: Modifying mean molecular weight from ',gmw,' to ',gmw_nicil
 endif
 nden_nimhd = 0.0
!
!--Initialise and verify parameters for super-timestepping
!
#ifdef STS_TIMESTEPS
 call sts_initialise(ierr,dtdiff)
 if (ierr > 0) call fatal('initial','supertimestep: nu > 1 or < 0 or NaN.')
#endif
!
!--initialise the equation of state
!  (must be done AFTER the units are known & AFTER mu is calculated in non-ideal MHD)
!
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('initial','error initialising equation of state')
!
!--get total number of particles (on all processors)
!
 ntot = reduceall_mpi('+',npart)
 call update_npartoftypetot
 if (id==master) write(iprint,"(a,i12)") ' npart total   = ',ntot
 if (npart > 0) then
    if (id==master .and. maxalpha==maxp)  write(iprint,*) 'mean alpha  initial: ',sum(alphaind(1,1:npart))/real(npart)
 endif

#ifdef DRIVING
!
!--initialise turbulence driving
!
 if (id==master) write(iprint,*) 'waiting on input for turbulent driving...'
 call init_forcing(dumpfile,infile,time)
#endif

 if (use_dust) then
    call init_drag(ierr)
    if (ierr /= 0) call fatal('initial','error initialising drag coefficients')
    if (use_dustgrowth) then
       call init_growth(ierr)
       if (ierr /= 0) call fatal('initial','error initialising growth variables')
       if (use_porosity) then
          call init_porosity(ierr)
          if (ierr /= 0) call fatal('initial','error initialising porosity variables')
          call init_filfac(npart,xyzh,vxyzu)
       endif
    endif
 endif
!
!--initialise cooling function
!  this will initialise all cooling variables, including if h2chemistry = true
 if (icooling > 0) call init_cooling(id,master,iprint,ierr)

 if (idamp > 0 .and. idamp < 3 .and. any(abs(vxyzu(1:3,:)) > tiny(0.)) .and. abs(time) < tiny(time)) then
    call error('setup','damping on: setting non-zero velocities to zero')
    vxyzu(1:3,:) = 0.
 endif

 ! initialise apr if it is being used
 if (use_apr) then
    call init_apr(apr_level,ierr)
 else
    apr_level(:) = 1
 endif

!
!--The code works in B/rho as its conservative variable, but writes B to dumpfile
!  So we now convert our primitive variable read, B, to the conservative B/rho
!  This necessitates computing the density sum.
!
 if (mhd .or. use_dustfrac) then
    if (npart > 0) then
       call set_linklist(npart,npart,xyzh,vxyzu)
       fxyzu = 0.
       call densityiterate(2,npart,npart,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                              fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level)
    endif

    ! now convert to B/rho
    do i=1,npart
       itype      = iamtype(iphase(i))
       hi         = xyzh(4,i)
       pmassi     = massoftype(itype)
       rhoi1      = 1.0/rhoh(hi,pmassi)
       if (mhd) then
          Bevol(1,i) = Bxyz(1,i) * rhoi1
          Bevol(2,i) = Bxyz(2,i) * rhoi1
          Bevol(3,i) = Bxyz(3,i) * rhoi1
       endif
       if (use_dustfrac) then
          !--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
          dustevol(:,i) = 0.
          dustevol(1:ndustsmall,i) = sqrt(dustfrac(1:ndustsmall,i)/(1.-dustfrac(1:ndustsmall,i)))
       endif
    enddo
 endif

 if (ind_timesteps) then
    ibin(:)       = 0
    ibin_old(:)   = 0
    ibin_wake(:)  = 0
    if (dt_read_in) call init_ibin(npart,dtmax)
    istepfrac     = 0
    ibinnow       = 0
 else
    dtcourant = huge(dtcourant)
    dtforce   = huge(dtforce)
 endif
 dtinject  = huge(dtinject)

!
!--balance domains prior to starting calculation
!  (make sure this is called AFTER iphase has been set)
!
 if (mpi) then
    do i=1,npart
       ibelong(i) = id
    enddo
    call balancedomains(npart)
    if (use_sinktree) then
       ibelong((maxpsph)+1:maxp) = -1
       boundi = (maxpsph)+(nptmass / nprocs)*id
       boundf = (maxpsph)+(nptmass / nprocs)*(id+1)
       if (id == nprocs-1) boundf = boundf + mod(nptmass,nprocs)
       ibelong(boundi+1:boundf) = id
    endif
 endif

 if (use_sinktree) then
    shortsinktree = 1 ! init shortsinktree to 1 to avoid any problem if nptmass change during the calculation
    fxyz_ptmass_tree = 0.
 endif

!
!--get timestep for external forces
!
 dtextforce = huge(dtextforce)
 fext(:,:)  = 0.

#ifdef GR
#ifdef PRIM2CONS_FIRST
 ! COMPUTE METRIC HERE
 call init_metric(npart,xyzh,metrics,metricderivs)
 ! -- The conserved quantites (momentum and entropy) are being computed
 ! -- directly from the primitive values in the starting dumpfile.
 call prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,use_dens=.false.,dens=dens)
 write(iprint,*) ''
 call warning('initial','using preprocessor flag -DPRIM2CONS_FIRST')
 write(iprint,'(a,/)') ' This means doing prim2cons BEFORE the initial density calculation for this simulation.'
#endif
 ! --- Need rho computed by sum to do primitive to conservative, since dens is not read from file
 if (npart>0) then
    call set_linklist(npart,npart,xyzh,vxyzu)
    fxyzu = 0.
    call densityiterate(2,npart,npart,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                              fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level)
 endif
#ifndef PRIM2CONS_FIRST
 call init_metric(npart,xyzh,metrics,metricderivs)
 call prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,use_dens=.false.,dens=dens)
#endif
 if (iexternalforce > 0 .and. imetric /= imet_minkowski) then
    call initialise_externalforces(iexternalforce,ierr)
    if (ierr /= 0) call fatal('initial','error in external force settings/initialisation')
    call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,fext,dtextforce,dens=dens)
 endif
#else
 if (iexternalforce > 0) then
    call initialise_externalforces(iexternalforce,ierr)
    call update_externalforce(iexternalforce,time,0.)
    if (ierr /= 0) call fatal('initial','error in external force settings/initialisation')
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,vxyzu,fext,time,iexternalforce,C_force) &
    !$omp private(i,poti,dtf,fextv) &
    !$omp reduction(min:dtextforce)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          call externalforce(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i), &
                              xyzh(4,i),time,fext(1,i),fext(2,i),fext(3,i),poti,dtf,i)
          dtextforce = min(dtextforce,C_force*dtf)
          ! add velocity-dependent part
          call externalforce_vdependent(iexternalforce,xyzh(1:3,i),vxyzu(1:3,i),fextv,poti)
          fext(1:3,i) = fext(1:3,i) + fextv
       endif
    enddo
    !$omp end parallel do
 endif
#endif

 if (iexternalforce > 0) then
    dtextforce = reduceall_mpi('min',dtextforce)
    if (id==master) write(iprint,*) 'dt(extforce)  = ',dtextforce
 endif

!
!-- Set external force to zero on boundary particles
!
 if (maxphase==maxp) then
!$omp parallel do default(none) &
!$omp shared(npart,fext,iphase) private(i)
    do i=1,npart
       if (iamtype(iphase(i))==iboundary) fext(:,i)=0.
    enddo
!$omp end parallel do
 endif
!
!--get timestep and forces for sink particles
!
 dtsinkgas = huge(dtsinkgas)
 r_crit2   = r_crit*r_crit
 rho_crit  = real(rho_crit_cgs/unit_density)
 r_merge_uncond2 = r_merge_uncond**2
 r_merge_cond2   = r_merge_cond**2
 r_merge2        = max(r_merge_uncond2,r_merge_cond2)
 if (rhofinal_cgs > 0.) then
    rhofinal1 = real(unit_density/rhofinal_cgs)
 else
    rhofinal1 = 0.0
 endif
 if (iH2R > 0 .and. id==master) then
    call initialize_H2R
 else
    isionised = .false.
 endif
 if (nptmass > 0) then
    if (id==master) write(iprint,"(a,i12)") ' nptmass       = ',nptmass
    if (iH2R > 0) call update_ionrates(nptmass,xyzmh_ptmass,h_acc)
    if (.not. gr) then
       ! compute initial sink-sink forces and get timestep
       if (use_regnbody) then
          call init_subgroup
          call group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,nmatrix)
       endif
       call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,dtsinksink,&
                               iexternalforce,time,merge_ij,merge_n,dsdt_ptmass,&
                               group_info,bin_info)
    endif
#ifdef GR
    ! calculate metric derivatives and the external force caused by the metric on the sink particles
    ! this will also return the timestep for sink-sink
    call init_metric(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass)
    call prim2consall(nptmass,xyzmh_ptmass,metrics_ptmass,&
                     vxyz_ptmass,pxyzu_ptmass,use_dens=.false.,use_sink=.true.)
    call get_grforce_all(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass,&
                     vxyz_ptmass,fxyz_ptmass,dtextforce,use_sink=.true.)
    ! sinks in GR, provide external force due to metric to determine the sink total force
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,dtsinksink,&
                             iexternalforce,time,merge_ij,merge_n,dsdt_ptmass)
#endif
    dtsinksink = C_force*dtsinksink
    if (id==master) write(iprint,*) 'dt(sink-sink) = ',dtsinksink
    dtextforce = min(dtextforce,dtsinksink)

    ! compute initial sink-gas forces and get timestep
    pmassi = massoftype(igas)
    ntypes = get_ntypes(npartoftype)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (ntypes > 1 .and. maxphase==maxp) then
             if (use_apr) then
                pmassi = aprmassoftype(iamtype(iphase(i)),apr_level(i))
             else
                pmassi = massoftype(iamtype(iphase(i)))
             endif
          elseif (use_apr) then
             pmassi = aprmassoftype(igas,apr_level(i))
          endif
          if (.not.use_sinktree) then
             call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass, &
                                     fext(1,i),fext(2,i),fext(3,i),poti,pmassi,fxyz_ptmass,&
                                     dsdt_ptmass,fonrmax,dtphi2,bin_info)
             dtsinkgas = min(dtsinkgas,C_force*1./sqrt(fonrmax),C_force*sqrt(dtphi2))
          endif
       endif
    enddo
    !
    ! reduction of sink-gas forces from each MPI thread
    !
    call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))

    if (id==master) write(iprint,*) 'dt(sink-gas)  = ',dtsinkgas

    dtextforce = min(dtextforce,dtsinkgas)
    !  Reduce dt over MPI tasks
    dtsinkgas = reduceall_mpi('min',dtsinkgas)
    dtextforce = reduceall_mpi('min',dtextforce)
    if (use_regnbody) call update_kappa(xyzmh_ptmass,vxyz_ptmass,bin_info,group_info,n_group)
 endif
 call init_ptmass(nptmass,logfile)
 if (gravity .and. icreate_sinks > 0 .and. id==master) then
    write(iprint,*) 'Sink radius and critical densities:'
    write(iprint,*) ' h_acc                    == ',h_acc*udist,'cm'
    write(iprint,*) ' h_fact*(m/rho_crit)^(1/3) = ',hfactfile*(massoftype(igas)/rho_crit)**(1./3.)*udist,'cm'
    write(iprint,*) ' rho_crit         == ',rho_crit_cgs,'g cm^{-3}'
    write(iprint,*) ' m(h_fact/h_acc)^3 = ', massoftype(igas)*(hfactfile/h_acc)**3*unit_density,'g cm^{-3}'
    if (r_merge_uncond < 2.0*h_acc) then
       write(iprint,*) ' WARNING! Sink creation is on, but but merging is off!  Suggest setting r_merge_uncond >= 2.0*h_acc'
    endif
    dsdt_ptmass = 0. ! could introduce NaN in ptmass spins if not initialised (no get_accel done before creating sink)
    fxyz_ptmass = 0.
    fxyz_ptmass_sinksink = 0.
 endif
 if (abs(time) <= tiny(0.)) then
    !initialize nucleation array at the start of the run only
    if (do_nucleation) call init_nucleation
    !initialize optical depth array tau
    if (itau_alloc == 1) tau = 0.
    !initialize Lucy optical depth array tau_lucy
    if (itauL_alloc == 1) tau_lucy = 2./3.
 endif
 if (update_muGamma) then
    eos_vars(igamma,:) = gamma
    eos_vars(imu,:) = gmw
    call set_abundances !to get mass_per_H
 endif
!
!--inject particles at t=0, and get timestep constraint on this
!
#ifdef INJECT_PARTICLES
 call init_inject(ierr)
 if (ierr /= 0) call fatal('initial','error initialising particle injection')
 !rename wind profile filename
 inquire(file='wind_profile1D.dat',exist=iexist)
 if (iexist) then
    i = len(trim(dumpfile))
    if (dumpfile(i-2:i) == 'tmp') then
       file1D = dumpfile(1:i-9) // '1D.dat'
    else
       file1D = dumpfile(1:i-5) // '1D.dat'
    endif
    call rename('wind_profile1D.dat',trim(file1D))
 endif
 npart_old = npart
 call inject_particles(time,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                       npart,npart_old,npartoftype,dtinject)
 call update_injected_particles(npart_old,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject)
#endif

!
!--set initial chemical abundance values
!
#ifdef KROME
 call initialise_krome()
 dtextforce = min(dtextforce,dtmax/2.0**10)  ! Required since a cooling timestep is not initialised for implicit cooling
#endif
!
!--calculate (all) derivatives the first time around
!
 dtnew_first   = dtmax  ! necessary in case ntot = 0
 nderivinit    = 1
 ! call derivs twice with Cullen-Dehnen switch to update accelerations
 if (maxalpha==maxp .and. nalpha >= 0) nderivinit = 2
 if (do_radiation) nderivinit = 1

 !$omp parallel do default(none) &
 !$omp shared(npart,eos_vars,fxyzu) &
 !$omp private(i)
 do i=1,npart
    eos_vars(3,i) = -1.0 ! initial guess for temperature overridden in eos
    fxyzu(:,i) = 0.      ! so that div_a is 0 in first call to viscosity switch
 enddo
 !$omp end parallel do

 do j=1,nderivinit
    if (ntot > 0) call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
                              rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,filfac,&
                              dustfrac,eos_vars,time,0.,dtnew_first,pxyzu,dens,metrics,apr_level)
#ifdef LIVE_ANALYSIS
    call do_analysis(dumpfile,numfromfile(dumpfile),xyzh,vxyzu, &
                     massoftype(igas),npart,time,ianalysis)
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,dustevol,&
                ddustevol,filfac,dustfrac,eos_vars,time,0.,dtnew_first,pxyzu,dens,metrics,apr_level)

    if (do_radiation) call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)
#endif
 enddo

 if (nalpha >= 2) then
    ialphaloc = 2
    !$omp parallel do private(i)
    do i=1,npart
       alphaind(1,i) = max(alphaind(1,i),alphaind(ialphaloc,i)) ! set alpha = max(alphaloc,alpha)
    enddo
 endif
!
!--set initial timestep
!
 if (.not.ind_timesteps) then
    dt = min(dtnew_first,dtinject)
    if (id==master) then
       write(iprint,*) 'dt(forces)    = ',dtforce
       write(iprint,*) 'dt(courant)   = ',dtcourant
       write(iprint,*) 'dt initial    = ',dt
    endif
 endif
!
!--initialise dynamic boundaries in the second instance
!
 if (dynamic_bdy) call init_dynamic_bdy(2,npart,nptmass,dtmax)
!
!--Calculate current centre of mass
!
 call get_centreofmass(xyzcom,dummy,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
!
!--write second header to logfile/screen
!
 if (id==master .and. read_input_files) call write_header(2,infile,evfile,logfile,dumpfile,ntot)

 call init_evfile(ievfile,evfile,.true.)
 call write_evfile(time,dt)
 if (id==master) call write_evlog(iprint)
#ifdef MFLOW
 call mflow_init(imflow,evfile,infile) !take evfile in input to create string.mf
 call mflow_write(time, dt)
#endif

#ifdef VMFLOW
 call vmflow_init(ivmflow,evfile,infile) !take evfile in input to create string_v.mflowv
 call vmflow_write(time, dt)
#endif

#ifdef BINPOS
 call binpos_init(ibinpos,evfile) !take evfile in input to create string.binpos
 call binpos_write(time, dt)
#endif
!
!--Determine the maximum separation of particles
 xmax = -0.5*huge(xmax)
 ymax = -0.5*huge(ymax)
 zmax = -0.5*huge(zmax)
 xmin = -xmax
 ymin = -ymax
 zmin = -zmax
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh) &
 !$omp private(i) &
 !$omp reduction(min:xmin,ymin,zmin) &
 !$omp reduction(max:xmax,ymax,zmax)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       xmin = min(xmin,xyzh(1,i))
       ymin = min(ymin,xyzh(2,i))
       zmin = min(zmin,xyzh(3,i))
       xmax = max(xmax,xyzh(1,i))
       ymax = max(ymax,xyzh(2,i))
       zmax = max(zmax,xyzh(3,i))
    endif
 enddo
 !$omp end parallel do

 xmin = reduceall_mpi('min',xmin)
 ymin = reduceall_mpi('min',ymin)
 zmin = reduceall_mpi('min',zmin)
 xmax = reduceall_mpi('max',xmax)
 ymax = reduceall_mpi('max',ymax)
 zmax = reduceall_mpi('max',zmax)

 dx = abs(xmax - xmin)
 dy = abs(ymax - ymin)
 dz = abs(zmax - zmin)
!
!--Print box sizes and masses
!
 if (id==master .and. iverbose >= 1) then
    if (get_conserv > 0.0) then
       write(iprint,'(1x,a)') 'Initial mass and extent of particle distribution (in code units):'
    else
       write(iprint,'(1x,a)') 'Mass and extent of the particle distribution:'
    endif
    write(iprint,'(2x,a,es18.6)') '     Total mass : ', mtot
    write(iprint,'(2x,a,es18.6)') 'x(max) - x(min) : ', dx
    write(iprint,'(2x,a,es18.6)') 'y(max) - y(min) : ', dy
    write(iprint,'(2x,a,es18.6)') 'z(max) - z(min) : ', dz
    write(iprint,'(a)') ' '
 endif
!
!--Set initial values for continual verification of conservation laws
!  get_conserve=0.5: update centre of mass only; get_conserve=1: update all; get_conserve=-1: update none
!
 if (get_conserv > 0.0) then
    etot_in   = etot
    angtot_in = angtot
    totmom_in = totmom
    mdust_in  = mdust
    mtot_in   = mtot
    if (id==master .and. iverbose >= 1) then
       write(iprint,'(1x,a)') 'Setting initial values to verify conservation laws:'
    endif
 else
    if (id==master .and. iverbose >= 1) then
       write(iprint,'(1x,a)') 'Reading initial values to verify conservation laws from previous run:'
    endif
 endif
 if (id==master) then
    if (iverbose >= 1) then
       write(iprint,'(1x,a,es18.6)') 'Initial total energy:     ', etot_in
       write(iprint,'(1x,a,es18.6)') 'Initial angular momentum: ', angtot_in
       write(iprint,'(1x,a,es18.6)') 'Initial linear momentum:  ', totmom_in
    endif
    if (use_dust) then
       dust_label = 'dust'
       call make_tags_unique(ndusttypes,dust_label)
       do i=1,ndusttypes
          if (mdust_in(i) > 0.) write(iprint,'(1x,a,es18.6)') 'Initial '//trim(dust_label(i))//' mass:     ',mdust_in(i)
       enddo
       write(iprint,'(1x,a,es18.6)') 'Initial total dust mass:', sum(mdust_in(:))
    endif
    if (use_apr) then
       write(iprint,'(1x,a,es18.6)') 'Initial total mass:  ', mtot_in
    endif
 endif
!
!--Print warnings of units if values are not reasonable
!
 tolu = 1.0e2
 toll = 1.0e-2
 if (get_conserv > 0.0) then
    get_conserv = -1.
    if (id==master) then
       if (abs(etot_in) > tolu ) call warning('initial',&
          'consider changing units to reduce abs(total energy)',var='etot',val=etot_in)
       if (mtot > tolu .or. mtot < toll)  call warning('initial',&
          'consider changing units so total mass is closer to unity',var='mtot',val=mtot)
       ! if (dx > tolu .or. dx < toll .or. dy > tolu .or. dy < toll .or. dz > tolu .or. dz < toll) &
       !call warning('initial','consider changing code-units to have box length closer to unity')
    endif
 endif
!
!--write initial conditions to output file
!  if the input file ends in .tmp or .init
!
 iposinit = index(dumpfile,'.init')
 ipostmp  = index(dumpfile,'.tmp')
 if (iposinit > 0 .or. ipostmp > 0) then
#ifdef HDF5
    dumpfileold = trim(dumpfile)//'.h5'
#else
    dumpfileold = dumpfile
#endif
    if (iposinit > 0) then
       dumpfile = trim(dumpfile(1:iposinit-1))
    else
       dumpfile = trim(dumpfile(1:ipostmp-1))
    endif
    call write_fulldump(time,trim(dumpfile))
    if (id==master) call write_infile(infile,logfile,evfile,trim(dumpfile),iwritein,iprint)
    !
    !  delete temporary dump file
    !
    call barrier_mpi() ! Ensure all procs have read temp file before deleting
    inquire(file=trim(dumpfileold),exist=iexist)
    if (id==master .and. iexist) then
       write(iprint,"(/,a,/)") ' ---> DELETING temporary dump file '//trim(dumpfileold)//' <---'
       open(unit=idisk1,file=trim(dumpfileold),status='old')
       close(unit=idisk1,status='delete')
    endif
 endif

 if (id==master) then
    call flush_warnings()
    call flush(iprint)
!
!--get starting cpu time
!
    call get_timings(twall_start,tcpu_start)
 endif

end subroutine startrun

!----------------------------------------------------------------
!+
!  Reset or deallocate things that were allocated in initialise
!+
!----------------------------------------------------------------
subroutine finalise()
 use dim, only: mpi
 use mpitree, only:finish_tree_comms
 use mpimemory, only:deallocate_mpi_memory

 if (mpi) then
    call finish_tree_comms()
    call deallocate_mpi_memory()
 endif

end subroutine finalise

!----------------------------------------------------------------
!+
!  This module ends the run (prints footer and closes log).
!  Only called by master thread.
!+
!----------------------------------------------------------------

subroutine endrun
 use io,       only:iprint,ievfile,iscfile,imflow,ivmflow,ibinpos,igpos
 use timing,   only:printused
 use part,     only:nptmass
 use eos,      only:ieos,finish_eos
 use ptmass,   only:finish_ptmass
 integer           :: ierr
 character(len=10) :: finishdate, finishtime

 call finalise()
 call finish_eos(ieos,ierr)

 write (iprint,"(/,'>',74('_'),'<')")
!
!--print time and date of finishing
!
 call date_and_time(finishdate,finishtime)
 finishdate = finishdate(7:8)//'/'//finishdate(5:6)//'/'//finishdate(1:4)
 finishtime = finishtime(1:2)//':'//finishtime(3:4)//':'//finishtime(5:)
 write(iprint,"(/,' Run finished on ',a,' at ',a,/)") finishdate,finishtime
!
!--print out total code timings:
!
 call printused(twall_start,'Total wall time:',iprint)

 write(iprint,40)
40 format(/, &
   6x,' |   |           |               | |   _|       | |         ',/, &
   6x,' __| __ \   _` | __|  __|   _` | | |  |    _ \  | |  /  __| ',/, &
   6x,' |   | | | (   | |  \__ \  (   | | |  __| (   | |   < \__ \ ',/, &
   6x,'\__|_| |_|\__,_|\__|____/ \__,_|_|_| _|  \___/ _|_|\_\____/ ',/)

 write (iprint,"('>',74('_'),'<')")
!
!--close ev, log& ptmass-related files
!
 close(unit=ievfile)
 close(unit=imflow)  ! does not matter if not open
 close(unit=ivmflow)
 close(unit=ibinpos)
 close(unit=igpos)
 if (iprint /= 6) close(unit=iprint)

 if (iscfile > 0) close(unit=iscfile)

 call finish_ptmass(nptmass)

end subroutine endrun

end module initial
