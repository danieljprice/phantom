!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_infile
!
! This module contains all routines required for
!  reading and writing of input file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - C_cour             : *Courant number*
!   - C_force            : *dt_force number*
!   - X                  : *hydrogen mass fraction for MESA opacity table*
!   - Z                  : *metallicity for MESA opacity table*
!   - alpha              : *shock viscosity parameter*
!   - alphaB             : *shock resistivity parameter*
!   - alphamax           : *MAXIMUM shock viscosity parameter*
!   - alphau             : *shock conductivity parameter*
!   - avdecayconst       : *decay time constant for viscosity switches*
!   - beta               : *beta viscosity*
!   - bulkvisc           : *magnitude of bulk viscosity*
!   - calc_erot          : *include E_rot in the ev_file*
!   - cv_type            : *how to get cv and mean mol weight (0=constant,1=mesa)*
!   - dtmax              : *time between dumps*
!   - dtmax_dratio       : *dynamic dtmax: density ratio controlling decrease (<=0 to ignore)*
!   - dtmax_max          : *dynamic dtmax: maximum allowed dtmax (=dtmax if <= 0)*
!   - dtmax_min          : *dynamic dtmax: minimum allowed dtmax*
!   - dtwallmax          : *maximum wall time between dumps (hhh:mm, 000:00=ignore)*
!   - dumpfile           : *dump file to start from*
!   - flux_limiter       : *limit radiation flux*
!   - hfact              : *h in units of particle spacing [h = hfact(m/rho)^(1/3)]*
!   - ien_type           : *energy variable (0=auto, 1=entropy, 2=energy, 3=entropy_s)*
!   - implicit_radiation : *use implicit integration (Whitehouse, Bate & Monaghan 2005)*
!   - iopacity_type      : *opacity method (0=inf,1=mesa,2=constant,-1=preserve)*
!   - ipdv_heating       : *heating from PdV work (0=off, 1=on)*
!   - irealvisc          : *physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)*
!   - ireconav           : *use reconstruction in shock viscosity (-1=off,0=no limiter,1=Van Leer)*
!   - iresistive_heating : *resistive heating (0=off, 1=on)*
!   - ishock_heating     : *shock heating (0=off, 1=on)*
!   - itsmax_rad         : *max number of iterations allowed in implicit solver*
!   - iverbose           : *verboseness of log (-1=quiet 0=default 1=allsteps 2=debug 5=max)*
!   - kappa_cgs          : *constant opacity value in cm2/g*
!   - logfile            : *file to which output is directed*
!   - nfulldump          : *full dump every n dumps*
!   - nmax               : *maximum number of timesteps (0=just get derivs and stop)*
!   - nmaxdumps          : *stop after n full dumps (-ve=ignore)*
!   - nout               : *write dumpfile every n dtmax (-ve=ignore)*
!   - overcleanfac       : *factor to increase cleaning speed (decreases time step)*
!   - psidecayfac        : *div B diffusion parameter*
!   - ptol               : *tolerance on pmom iterations*
!   - rhofinal_cgs       : *maximum allowed density (cgs) (<=0 to ignore)*
!   - rkill              : *deactivate particles outside this radius (<0 is off)*
!   - shearparam         : *magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)*
!   - tmax               : *end time*
!   - tol_rad            : *tolerance on backwards Euler implicit solve of dxi/dt*
!   - tolh               : *tolerance on h-rho iterations*
!   - tolv               : *tolerance on v iterations in timestepping*
!   - twallmax           : *maximum wall time (hhh:mm, 000:00=ignore)*
!   - use_mcfost         : *use the mcfost library*
!   - xtol               : *tolerance on xyz iterations*
!
! :Dependencies: HIIRegion, apr, boundary_dyn, cooling, damping, dim, dust,
!   dust_formation, eos, externalforces, forcing, gravwaveutils, growth,
!   infile_utils, inject, io, linklist, metric, nicil_sup, options, part,
!   porosity, ptmass, ptmass_radiation, radiation_implicit,
!   radiation_utils, timestep, viscosity
!
 use timestep,  only:dtmax_dratio,dtmax_max,dtmax_min
 use options,   only:nfulldump,nmaxdumps,twallmax,iexternalforce,tolh, &
                     alpha,alphau,alphaB,beta,avdecayconst,damp,rkill, &
                     ipdv_heating,ishock_heating,iresistive_heating,ireconav, &
                     icooling,psidecayfac,overcleanfac,alphamax,calc_erot,rhofinal_cgs, &
                     use_mcfost,use_Voronoi_limits_file,Voronoi_limits_file,use_mcfost_stellar_parameters,&
                     exchange_radiation_energy,limit_radiation_flux,iopacity_type,mcfost_computes_Lacc,&
                     mcfost_uses_PdV,implicit_radiation,mcfost_keep_part,ISM, mcfost_dust_subl
 use timestep,  only:dtwallmax,tolv,xtol,ptol
 use viscosity, only:irealvisc,shearparam,bulkvisc
 use part,      only:hfact,ien_type
 use io,        only:iverbose
 use dim,       only:do_radiation,nucleation,use_dust,use_dustgrowth,mhd_nonideal,compiled_with_mcfost
 implicit none
 logical :: incl_runtime2 = .false.

contains

!-----------------------------------------------------------------
!+
!  writes an input file
!+
!-----------------------------------------------------------------
subroutine write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
 use timestep,        only:tmax,dtmax,dtmax_user,nmax,nout,C_cour,C_force,C_ent
 use io,              only:fatal
 use infile_utils,    only:write_inopt
#ifdef DRIVING
 use forcing,         only:write_options_forcing
#endif
 use externalforces,  only:write_options_externalforces
 use damping,         only:write_options_damping
 use linklist,        only:write_inopts_link
 use dust,            only:write_options_dust
 use growth,          only:write_options_growth
 use porosity,        only:write_options_porosity
#ifdef INJECT_PARTICLES
 use inject,          only:write_options_inject,inject_type,update_injected_par
#endif
 use apr,             only:write_options_apr
 use dust_formation,  only:write_options_dust_formation
 use nicil_sup,       only:write_options_nicil
 use metric,          only:write_options_metric
 use eos,             only:write_options_eos,ieos,X_in,Z_in
 use ptmass,          only:write_options_ptmass
 use ptmass_radiation,only:write_options_ptmass_radiation
 use cooling,         only:write_options_cooling
 use gravwaveutils,   only:write_options_gravitationalwaves
 use radiation_utils,    only:kappa_cgs
 use radiation_implicit, only:tol_rad,itsmax_rad,cv_type
 use dim,                only:maxvxyzu,maxptmass,gravity,sink_radiation,gr,&
                              nalpha,use_apr
 use part,               only:maxp,mhd,maxalpha,nptmass
 use boundary_dyn,       only:write_options_boundary
 use HIIRegion,          only:write_options_H2R
 character(len=*), intent(in) :: infile,logfile,evfile,dumpfile
 integer,          intent(in) :: iwritein,iprint
 integer                      :: ierr
 character(len=10)            :: startdate,starttime

 if (iwritein /= iprint) then
    open(unit=iwritein,iostat=ierr,file=infile,status='replace',form='formatted')
    if (ierr /= 0) call fatal('write_infile','error creating input file, exiting...')
 endif

 if (iwritein /= iprint) then
!
! get date and time to timestamp the header
!
    call date_and_time(startdate,starttime)
    startdate = startdate(7:8)//'/'//startdate(5:6)//'/'//startdate(1:4)
    starttime = starttime(1:2)//':'//starttime(3:4)//':'//starttime(5:)
!
! write header to the file
!
    write(iwritein,"(a)") '# Runtime options file for Phantom, written '//startdate//' '//starttime
    write(iwritein,"(a)") '# Options not present assume their default values'
    write(iwritein,"(a)") '# This file is updated automatically after a full dump'
 endif

 write(iwritein,"(/,a)") '# job name'
 call write_inopt(trim(logfile),'logfile','file to which output is directed',iwritein)
 call write_inopt(trim(dumpfile),'dumpfile','dump file to start from',iwritein)

 write(iwritein,"(/,a)") '# options controlling run time and input/output'
 if (dtmax_user < 0.) dtmax_user = dtmax ! this should only ever be true for phantomsetup
 call write_inopt(tmax,'tmax','end time',iwritein)
 call write_inopt(dtmax_user,'dtmax','time between dumps',iwritein)
 call write_inopt(nmax,'nmax','maximum number of timesteps (0=just get derivs and stop)',iwritein)
 call write_inopt(nout,'nout','write dumpfile every n dtmax (-ve=ignore)',iwritein)
 call write_inopt(nmaxdumps,'nmaxdumps','stop after n full dumps (-ve=ignore)',iwritein)
 call write_inopt(real(twallmax),'twallmax','maximum wall time (hhh:mm, 000:00=ignore)',iwritein,time=.true.)
 call write_inopt(real(dtwallmax),'dtwallmax','maximum wall time between dumps (hhh:mm, 000:00=ignore)',iwritein,time=.true.)
 call write_inopt(nfulldump,'nfulldump','full dump every n dumps',iwritein)
 call write_inopt(iverbose,'iverbose','verboseness of log (-1=quiet 0=default 1=allsteps 2=debug 5=max)',iwritein)

 if (incl_runtime2 .or. rhofinal_cgs > 0.0 .or. dtmax_dratio > 1.0 .or. calc_erot) then
    write(iwritein,"(/,a)") '# options controlling run time and input/output: supplementary features'
    call write_inopt(rhofinal_cgs,'rhofinal_cgs','maximum allowed density (cgs) (<=0 to ignore)',iwritein)
    call write_inopt(dtmax_dratio,'dtmax_dratio','dynamic dtmax: density ratio controlling decrease (<=0 to ignore)',iwritein)
    call write_inopt(dtmax_max,'dtmax_max','dynamic dtmax: maximum allowed dtmax (=dtmax if <= 0)',iwritein)
    call write_inopt(dtmax_min,'dtmax_min','dynamic dtmax: minimum allowed dtmax',iwritein)
    call write_inopt(calc_erot,'calc_erot','include E_rot in the ev_file',iwritein)
 endif

 write(iwritein,"(/,a)") '# options controlling accuracy'
 call write_inopt(C_cour,'C_cour','Courant number',iwritein)
 call write_inopt(C_force,'C_force','dt_force number',iwritein)
 call write_inopt(tolv,'tolv','tolerance on v iterations in timestepping',iwritein,exp=.true.)
 call write_inopt(hfact,'hfact','h in units of particle spacing [h = hfact(m/rho)^(1/3)]',iwritein)
 call write_inopt(tolh,'tolh','tolerance on h-rho iterations',iwritein,exp=.true.)
 if (gr) then
    call write_inopt(C_ent,'C_ent','restrict timestep when ds/dt is too large (not used if ien_type != 3)',iwritein)
    call write_inopt(xtol,'xtol','tolerance on xyz iterations',iwritein)
    call write_inopt(ptol,'ptol','tolerance on pmom iterations',iwritein)
 endif

 call write_inopts_link(iwritein)

 write(iwritein,"(/,a)") '# options controlling hydrodynamics, shock capturing'
 if (maxalpha==maxp .and. nalpha > 0) then
    call write_inopt(alpha,'alpha','MINIMUM shock viscosity parameter',iwritein)
    call write_inopt(alphamax,'alphamax','MAXIMUM shock viscosity parameter',iwritein)
 else
    call write_inopt(alpha,'alpha','shock viscosity parameter',iwritein)
 endif
 if (maxvxyzu >= 4) then
    call write_inopt(alphau,'alphau','shock conductivity parameter',iwritein)
 endif
 if (mhd) then
    call write_inopt(alphaB,'alphaB','shock resistivity parameter',iwritein)
    call write_inopt(psidecayfac,'psidecayfac','div B diffusion parameter',iwritein)
    call write_inopt(overcleanfac,'overcleanfac','factor to increase cleaning speed (decreases time step)',iwritein)
 endif
 call write_inopt(beta,'beta','beta viscosity',iwritein)
 if (maxalpha==maxp .and. maxp > 0) then
    call write_inopt(avdecayconst,'avdecayconst','decay time constant for viscosity switches',iwritein)
 endif
 if (gr) then
    call write_inopt(ireconav,'ireconav','use reconstruction in shock viscosity (-1=off,0=no limiter,1=Van Leer)',iwritein)
 endif
 call write_options_damping(iwritein)

 !
 ! thermodynamics
 !
 call write_options_eos(iwritein)
 if (maxvxyzu >= 4 .and. (ieos==2 .or. ieos==5 .or. ieos==10 .or. ieos==15 .or. ieos==12 .or. ieos==16 &
      .or. ieos==17 .or. ieos==21 .or. ieos==22 .or. ieos==24) ) then
    call write_inopt(ipdv_heating,'ipdv_heating','heating from PdV work (0=off, 1=on)',iwritein)
    call write_inopt(ishock_heating,'ishock_heating','shock heating (0=off, 1=on)',iwritein)
    if (mhd) then
       call write_inopt(iresistive_heating,'iresistive_heating','resistive heating (0=off, 1=on)',iwritein)
    endif
    if (gr) then
       call write_inopt(ien_type,'ien_type','energy variable (0=auto, 1=entropy, 2=energy, 3=entropy_s)',iwritein)
    endif
 endif

 if (maxvxyzu >= 4) call write_options_cooling(iwritein)

 if (compiled_with_mcfost) then
    call write_inopt(use_mcfost,'use_mcfost','use the mcfost library',iwritein)
    if (use_Voronoi_limits_file) call write_inopt(Voronoi_limits_file,'Voronoi_limits_file',&
         'Limit file for the Voronoi tesselation',iwritein)
    call write_inopt(use_mcfost_stellar_parameters,'use_mcfost_stars',&
         'Fix the stellar parameters to mcfost values or update using sink mass',iwritein)
    call write_inopt(mcfost_computes_Lacc,'mcfost_computes_Lacc',&
         'Should mcfost compute the accretion luminosity',iwritein)
    call write_inopt(mcfost_uses_PdV,'mcfost_uses_PdV',&
         'Should mcfost use the PdV work and shock heating?',iwritein)
    call write_inopt(mcfost_keep_part,'mcfost_keep_part',&
         'Fraction of particles to keep for MCFOST',iwritein)
    call write_inopt(ISM,'ISM',&
         'ISM heating : 0 -> no ISM radiation field, 1 -> ProDiMo, 2 -> Bate & Keto',iwritein)
    call write_inopt(mcfost_dust_subl,'mcfost_dust_subl',&
         'Should mcfost do dust sublimation (experimental!)',iwritein)
 endif

 ! only write sink options if they are used, or if self-gravity is on
 if (nptmass > 0 .or. gravity) call write_options_ptmass(iwritein)

 call write_options_externalforces(iwritein,iexternalforce)

 write(iwritein,"(/,a)") '# options controlling physical viscosity'
 call write_inopt(irealvisc,'irealvisc','physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)',iwritein)
 call write_inopt(shearparam,'shearparam','magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)',iwritein)
 call write_inopt(bulkvisc,'bulkvisc','magnitude of bulk viscosity',iwritein)

#ifdef DRIVING
 call write_options_forcing(iwritein)
#endif

 if (use_dust) call write_options_dust(iwritein)
 if (use_dustgrowth) then
    call write_options_growth(iwritein)
    call write_options_porosity(iwritein)
 endif

 write(iwritein,"(/,a)") '# options for injecting/removing particles'
#ifdef INJECT_PARTICLES
 call write_options_inject(iwritein)
 if (inject_type=='sim') call update_injected_par()
#endif
 call write_inopt(rkill,'rkill','deactivate particles outside this radius (<0 is off)',iwritein)

 if (nucleation) call write_options_dust_formation(iwritein)

 if (sink_radiation) then
    write(iwritein,"(/,a)") '# options controling radiation pressure from sink particles'
    call write_options_ptmass_radiation(iwritein)
 endif

 if (mhd_nonideal) call write_options_nicil(iwritein)

 if (do_radiation) then
    write(iwritein,"(/,a)") '# options for radiation'
    call write_inopt(implicit_radiation,'implicit_radiation','use implicit integration (Whitehouse, Bate & Monaghan 2005)',iwritein)
    call write_inopt(exchange_radiation_energy,'gas-rad_exchange','exchange energy between gas and radiation',iwritein)
    call write_inopt(limit_radiation_flux,'flux_limiter','limit radiation flux',iwritein)
    call write_inopt(iopacity_type,'iopacity_type','opacity method (0=inf,1=mesa,2=constant,-1=preserve)',iwritein)
    if (iopacity_type == 1) then
       call write_inopt(X_in,'X','hydrogen mass fraction for MESA opacity table',iwritein)
       call write_inopt(Z_in,'Z','metallicity for MESA opacity table',iwritein)
    elseif (iopacity_type == 2) then
       call write_inopt(kappa_cgs,'kappa_cgs','constant opacity value in cm2/g',iwritein)
    endif
    if (implicit_radiation) then
       call write_inopt(tol_rad,'tol_rad','tolerance on backwards Euler implicit solve of dxi/dt',iwritein)
       call write_inopt(itsmax_rad,'itsmax_rad','max number of iterations allowed in implicit solver',iwritein)
       call write_inopt(cv_type,'cv_type','how to get cv and mean mol weight (0=constant,1=mesa)',iwritein)
    endif
 endif
 if (gr) call write_options_metric(iwritein)
 call write_options_gravitationalwaves(iwritein)
 call write_options_boundary(iwritein)

 if (use_apr) call write_options_apr(iwritein)

 call write_options_H2R(iwritein)

 if (iwritein /= iprint) close(unit=iwritein)
 if (iwritein /= iprint) write(iprint,"(/,a)") ' input file '//trim(infile)//' written successfully.'

end subroutine write_infile

!-----------------------------------------------------------------
!+
!  reads parameters for the run from the input file
!+
!-----------------------------------------------------------------
subroutine read_infile(infile,logfile,evfile,dumpfile)
 use dim,             only:maxvxyzu,maxptmass,gravity,sink_radiation,nucleation,&
                           itau_alloc,store_dust_temperature,gr,do_nucleation,use_apr
 use timestep,        only:tmax,dtmax,nmax,nout,C_cour,C_force,C_ent
 use eos,             only:read_options_eos,ieos,eos_requires_isothermal
 use io,              only:ireadin,iwritein,iprint,warn,die,error,fatal,id,master,fileprefix
 use infile_utils,    only:read_next_inopt,contains_loop,write_infile_series
#ifdef DRIVING
 use forcing,         only:read_options_forcing,write_options_forcing
#endif
 use externalforces,  only:read_options_externalforces
 use linklist,        only:read_inopts_link
 use dust,            only:read_options_dust
 use growth,          only:read_options_growth
 use options,         only:use_porosity
 use porosity,        only:read_options_porosity
 use metric,          only:read_options_metric
#ifdef INJECT_PARTICLES
 use inject,          only:read_options_inject
#endif
 use apr,             only:read_options_apr
 use dust_formation,  only:read_options_dust_formation,idust_opacity
 use nicil_sup,       only:read_options_nicil
 use part,            only:mhd,nptmass
 use cooling,         only:read_options_cooling
 use ptmass,          only:read_options_ptmass
 use ptmass_radiation,   only:read_options_ptmass_radiation,isink_radiation,&
                              alpha_rad,iget_tdust,iray_resolution
 use radiation_utils,    only:kappa_cgs
 use radiation_implicit, only:tol_rad,itsmax_rad,cv_type
 use damping,         only:read_options_damping
 use gravwaveutils,   only:read_options_gravitationalwaves
 use boundary_dyn,    only:read_options_boundary
 use HIIRegion,       only:read_options_H2R
 character(len=*), parameter   :: label = 'read_infile'
 character(len=*), intent(in)  :: infile
 character(len=*), intent(out) :: logfile,evfile,dumpfile
 character(len=len(infile)+4)  :: infilenew
 character(len=10) :: cline
 character(len=20) :: name
 character(len=120) :: valstring
 integer :: ierr,ireaderr,line,idot,ngot,nlinesread
 real    :: ratio
 logical :: imatch,igotallrequired,igotallturb,igotalllink,igotloops
 logical :: igotallbowen,igotallcooling,igotalldust,igotallextern,igotallinject,igotallgrowth,igotallporosity
 logical :: igotallionise,igotallnonideal,igotalleos,igotallptmass,igotalldamping,igotallapr
 logical :: igotallprad,igotalldustform,igotallgw,igotallgr,igotallbdy,igotallH2R
 integer, parameter :: nrequired = 1

 ireaderr = 0
 ierr     = 0
 line     = 0
 idot     = index(infile,'.in') - 1
 if (idot <= 1) idot = len_trim(infile)
 logfile    = infile(1:idot)//'01.log'
 dumpfile   = infile(1:idot)//'_00000.tmp'
 fileprefix = infile(1:idot)
 ngot            = 0
 igotallturb     = .true.
 igotalldust     = .true.
 igotallgrowth   = .true.
 igotallporosity = .true.
 igotalllink     = .true.
 igotallextern   = .true.
 igotallinject   = .true.
 igotallapr      = .true.
 igotalleos      = .true.
 igotallcooling  = .true.
 igotalldamping  = .true.
 igotloops       = .false.
 igotallionise   = .true.
 igotallnonideal = .true.
 igotallbowen    = .true.
 igotallptmass   = .true.
 igotallprad     = .true.
 igotalldustform = .true.
 igotallgw       = .true.
 igotallgr       = .true.
 igotallbdy      = .true.
 igotallH2R      = .true.
 use_Voronoi_limits_file = .false.

 open(unit=ireadin,err=999,file=infile,status='old',form='formatted')
 do while (ireaderr == 0)
    call read_next_inopt(name,valstring,ireadin,ireaderr,nlinesread)
    if (contains_loop(valstring)) igotloops = .true.
    line = line + nlinesread

    !print*,'name: '//trim(name),' value: '//trim(valstring)
    select case(trim(name))
    case('logfile')
       logfile = trim(valstring(1:min(len(logfile),len(valstring))))
       idot = index(logfile,'.log') - 1
       if (idot <= 1) idot = len_trim(logfile)
       evfile  = logfile(1:idot)//'.ev'
    case('dumpfile')
       dumpfile = trim(valstring(1:min(len(dumpfile),len(valstring))))
       ngot = ngot + 1
    case('tmax')
       read(valstring,*,iostat=ierr) tmax
    case('dtmax')
       read(valstring,*,iostat=ierr) dtmax
    case('nmax')
       read(valstring,*,iostat=ierr) nmax
    case('nout')
       read(valstring,*,iostat=ierr) nout
    case('nmaxdumps')
       read(valstring,*,iostat=ierr) nmaxdumps
    case('twallmax')
       read(valstring,*,iostat=ierr) twallmax
    case('dtwallmax')
       read(valstring,*,iostat=ierr) dtwallmax
    case('iverbose')
       read(valstring,*,iostat=ierr) iverbose
    case('rhofinal_cgs')
       read(valstring,*,iostat=ierr) rhofinal_cgs
       incl_runtime2 = .true.
    case('calc_erot')
       read(valstring,*,iostat=ierr) calc_erot
       incl_runtime2 = .true.
    case('dtmax_dratio')
       read(valstring,*,iostat=ierr) dtmax_dratio
       incl_runtime2 = .true.
    case('dtmax_max')
       read(valstring,*,iostat=ierr) dtmax_max
       if (dtmax_max <= 0.0) dtmax_max = dtmax
       ! to prevent comparison errors from round-off
       ratio = dtmax_max/dtmax
       ratio = int(ratio+0.5)+0.0001
       dtmax_max = dtmax*ratio
    case('dtmax_min')
       read(valstring,*,iostat=ierr) dtmax_min
       ! to prevent comparison errors from round-off
       if (dtmax_min > epsilon(dtmax_min)) then
          ratio = dtmax/dtmax_min
          ratio = int(ratio+0.5)+0.0001
          dtmax_min = dtmax/ratio
       endif
    case('C_cour')
       read(valstring,*,iostat=ierr) C_cour
    case('C_force')
       read(valstring,*,iostat=ierr) C_force
    case('tolv')
       read(valstring,*,iostat=ierr) tolv
    case('C_ent')
       read(valstring,*,iostat=ierr) C_ent
    case('xtol')
       read(valstring,*,iostat=ierr) xtol
    case('ptol')
       read(valstring,*,iostat=ierr) ptol
    case('hfact')
       read(valstring,*,iostat=ierr) hfact
    case('tolh')
       read(valstring,*,iostat=ierr) tolh
    case('rkill')
       read(valstring,*,iostat=ierr) rkill
    case('nfulldump')
       read(valstring,*,iostat=ierr) nfulldump
    case('alpha')
       read(valstring,*,iostat=ierr) alpha
    case('alphamax')
       read(valstring,*,iostat=ierr) alphamax
    case('alphau')
       read(valstring,*,iostat=ierr) alphau
    case('alphaB')
       read(valstring,*,iostat=ierr) alphaB
    case('psidecayfac')
       read(valstring,*,iostat=ierr) psidecayfac
    case('overcleanfac')
       read(valstring,*,iostat=ierr) overcleanfac
    case('beta')
       read(valstring,*,iostat=ierr) beta
    case('ireconav')
       read(valstring,*,iostat=ierr) ireconav
    case('avdecayconst')
       read(valstring,*,iostat=ierr) avdecayconst
    case('ipdv_heating')
       read(valstring,*,iostat=ierr) ipdv_heating
    case('ishock_heating')
       read(valstring,*,iostat=ierr) ishock_heating
    case('iresistive_heating')
       read(valstring,*,iostat=ierr) iresistive_heating
    case('ien_type')
       read(valstring,*,iostat=ierr) ien_type
    case('irealvisc')
       read(valstring,*,iostat=ierr) irealvisc
    case('shearparam')
       read(valstring,*,iostat=ierr) shearparam
    case('bulkvisc')
       read(valstring,*,iostat=ierr) bulkvisc
    case('use_mcfost')
       read(valstring,*,iostat=ierr) use_mcfost
    case('Voronoi_limits_file')
       read(valstring,*,iostat=ierr) Voronoi_limits_file
       use_Voronoi_limits_file = .true.
    case('use_mcfost_stars')
       read(valstring,*,iostat=ierr) use_mcfost_stellar_parameters
    case('mcfost_computes_Lacc')
       read(valstring,*,iostat=ierr) mcfost_computes_Lacc
    case('mcfost_uses_PdV')
       read(valstring,*,iostat=ierr) mcfost_uses_PdV
    case('mcfost_keep_part')
       read(valstring,*,iostat=ierr) mcfost_keep_part
    case('ISM')
       read(valstring,*,iostat=ierr) ISM
    case('mcfost_dust_subl')
       read(valstring,*,iostat=ierr) mcfost_dust_subl
    case('implicit_radiation')
       read(valstring,*,iostat=ierr) implicit_radiation
       if (implicit_radiation) store_dust_temperature = .true.
    case('gas-rad_exchange')
       read(valstring,*,iostat=ierr) exchange_radiation_energy
    case('flux_limiter')
       read(valstring,*,iostat=ierr) limit_radiation_flux
    case('iopacity_type')
       read(valstring,*,iostat=ierr) iopacity_type
    case('cv_type')
       read(valstring,*,iostat=ierr) cv_type
    case('kappa_cgs')
       read(valstring,*,iostat=ierr) kappa_cgs
    case('tol_rad')
       read(valstring,*,iostat=ierr) tol_rad
    case('itsmax_rad')
       read(valstring,*,iostat=ierr) itsmax_rad
    case default
       imatch = .false.
       if (.not.imatch) call read_options_externalforces(name,valstring,imatch,igotallextern,ierr,iexternalforce)
#ifdef DRIVING
       if (.not.imatch) call read_options_forcing(name,valstring,imatch,igotallturb,ierr)
#endif
       if (.not.imatch) call read_inopts_link(name,valstring,imatch,igotalllink,ierr)
       !--Extract if one-fluid dust is used from the fileid
       if (.not.imatch .and. use_dust) call read_options_dust(name,valstring,imatch,igotalldust,ierr)
       if (.not.imatch .and. use_dustgrowth) call read_options_growth(name,valstring,imatch,igotallgrowth,ierr)
       if (.not.imatch .and. use_porosity) call read_options_porosity(name,valstring,imatch,igotallporosity,ierr)
       if (.not.imatch .and. gr) call read_options_metric(name,valstring,imatch,igotallgr,ierr)
#ifdef INJECT_PARTICLES
       if (.not.imatch) call read_options_inject(name,valstring,imatch,igotallinject,ierr)
#endif
       if (.not.imatch .and. use_apr) call read_options_apr(name,valstring,imatch,igotallapr,ierr)
       if (.not.imatch .and. nucleation) call read_options_dust_formation(name,valstring,imatch,igotalldustform,ierr)
       if (.not.imatch .and. sink_radiation) then
          call read_options_ptmass_radiation(name,valstring,imatch,igotallprad,ierr)
       endif
       if (.not.imatch .and. mhd_nonideal) call read_options_nicil(name,valstring,imatch,igotallnonideal,ierr)
       if (.not.imatch) call read_options_eos(name,valstring,imatch,igotalleos,ierr)
       if (.not.imatch .and. maxvxyzu >= 4) call read_options_cooling(name,valstring,imatch,igotallcooling,ierr)
       if (.not.imatch) call read_options_damping(name,valstring,imatch,igotalldamping,ierr)
       if (maxptmass > 0) then
          if (.not.imatch) call read_options_ptmass(name,valstring,imatch,igotallptmass,ierr)
          !
          ! read whatever sink options are present, but make them not compulsory
          ! if there are no sinks used and there is no self-gravity
          !
          if (nptmass==0 .and. .not.gravity) igotallptmass = .true.
       endif
       if (.not.imatch) call read_options_gravitationalwaves(name,valstring,imatch,igotallgw,ierr)
       if (.not.imatch) call read_options_boundary(name,valstring,imatch,igotallbdy,ierr)
       if (.not.imatch) call read_options_H2R(name,valstring,imatch,igotallH2R,ierr)
       if (len_trim(name) /= 0 .and. .not.imatch) then
          call warn('read_infile','unknown variable '//trim(adjustl(name))// &
                     ' in input file, value = '//trim(adjustl(valstring)))
       endif
    end select
    if (ierr /= 0 .and. len_trim(name) /= 0) &
       call fatal('read_infile','error extracting '//trim(adjustl(name))//' from input file')
 enddo
 close(unit=ireadin)

 igotallrequired = (ngot  >=  nrequired) .and. igotalllink   .and. igotallbowen   .and. igotalldust &
                    .and. igotalleos    .and. igotallcooling .and. igotallextern  .and. igotallturb &
                    .and. igotallptmass .and. igotallinject  .and. igotallionise  .and. igotallnonideal &
                    .and. igotallgrowth  .and. igotallporosity .and. igotalldamping .and. igotallprad &
                    .and. igotalldustform .and. igotallgw .and. igotallgr .and. igotallbdy .and. igotallapr

 if (ierr /= 0 .or. ireaderr > 0 .or. .not.igotallrequired) then
    ierr = 1
    if (id==master) then
       if (igotallrequired) then
          write(cline,"(i10)") line
          call error('read_infile','error reading '//trim(infile)//' at line '//trim(adjustl(cline)))
          infilenew = trim(infile)//'.new'
       else
          call error('read_infile','input file '//trim(infile)//' is incomplete for current compilation')
          if (.not.igotalleos) write(*,*) 'missing equation of state options'
          if (.not.igotallcooling) write(*,*) 'missing cooling options'
          if (.not.igotalldamping) write(*,*) 'missing damping options'
          if (.not.igotalllink) write(*,*) 'missing link options'
          if (.not.igotallbowen) write(*,*) 'missing Bowen dust options'
          if (.not.igotalldust) write(*,*) 'missing dust options'
          if (.not.igotallgr) write(*,*) 'missing metric parameters (eg, spin, mass)'
          if (.not.igotallgrowth) write(*,*) 'missing growth options'
          if (.not.igotallporosity) write(*,*) 'missing porosity options'

          if (.not.igotallextern) then
             if (gr) then
                write(*,*) 'missing GR quantities (eg: accretion radius)'
             else
                write(*,*) 'missing external force options'
             endif
          endif
          if (.not.igotallinject) write(*,*) 'missing inject-particle options'
          if (.not.igotallapr) write(*,*) 'missing apr options'
          if (.not.igotallionise) write(*,*) 'missing ionisation options'
          if (.not.igotallnonideal) write(*,*) 'missing non-ideal MHD options'
          if (.not.igotallturb) write(*,*) 'missing turbulence-driving options'
          if (.not.igotallprad) write(*,*) 'missing sink particle radiation options'
          if (.not.igotallptmass) write(*,*) 'missing sink particle options'
          if (.not.igotalldustform) write(*,*) 'missing dusty wind options'
          if (.not.igotallgw) write(*,*) 'missing gravitational wave options'
          infilenew = trim(infile)
       endif
       write(*,"(a)") ' REWRITING '//trim(infilenew)//' with all current and available options...'
       call write_infile(infilenew,logfile,evfile,dumpfile,iwritein,iprint)

       if (trim(infilenew) /= trim(infile)) then
          write(*,"(/,a)") ' useful commands to cut and paste:'
          write(*,*) ' diff '//trim(infilenew)//' '//trim(infile)
          write(*,*) ' mv '//trim(infilenew)//' '//trim(infile)
       else
          write(*,"(/,a)") ' (try again using revised input file)'
       endif
    endif
    call die
 elseif (igotloops) then
    if (id==master) then
       write(iprint,"(1x,a)") 'loops detected in input file:'
       call write_infile_series(ireadin,iwritein,infile,line,ierr)
       if (ierr /= 0) call fatal(label,'error writing input file series')
    endif
    call die
 endif
!
!--check options for possible errors
!
 if (id==master) then
    if (dtmax > tmax) call warn(label,'no output dtmax > tmax',1)
    if (nout > nmax)  call warn(label,'no output nout > nmax',1)
    if (nout==0)     call fatal(label,'nout = 0')
    if (C_cour <= 0.)  call fatal(label,'Courant number < 0')
    if (C_cour > 1.)  call fatal(label,'ridiculously big courant number!!')
    if (C_force <= 0.) call fatal(label,'bad choice for force timestep control')
    if (tolv <= 0.)    call fatal(label,'silly choice for tolv (< 0)')
    if (tolv > 1.e-1) call warn(label,'dangerously large tolerance on v iterations')
    if (xtol <= 0.)    call fatal(label,'silly choice for xtol (< 0)')
    if (xtol > 1.e-1) call warn(label,'dangerously large tolerance on xyz iterations')
    if (ptol <= 0.)    call fatal(label,'silly choice for ptol (< 0)')
    if (ptol > 1.e-1) call warn(label,'dangerously large tolerance on pmom iterations')
    if (nfulldump==0 .or. nfulldump > 10000) call fatal(label,'nfulldump = 0')
    if (nfulldump >= 50) call warn(label,'no full dumps for a long time...',1)
    if (twallmax < 0.)  call fatal(label,'invalid twallmax (use 000:00 to ignore)')
    if (dtwallmax < 0.) call fatal(label,'invalid dtwallmax (use 000:00 to ignore)')
    if (hfact < 1. .or. hfact > 5.) &
                         call warn(label,'ridiculous choice of hfact',4)
    if (tolh > 1.e-3)   call warn(label,'tolh is quite large!',2)
    if (tolh < epsilon(tolh)) call fatal(label,'tolh too small to ever converge')
    !if (damp < 0.)     call fatal(label,'damping < 0')
    !if (damp > 1.)     call warn(label,'damping ridiculously big')
    if (alpha < 0.)    call fatal(label,'stupid choice of alpha')
    if (alphau < 0.)   call fatal(label,'stupid choice of alphau')
    if (alpha > 10.)   call warn(label,'very large alpha, need to change timestep',2)
    if (alphau > 10.)  call warn(label,'very large alphau, check timestep',3)
    if (alphamax < tiny(alphamax)) call warn(label,'alphamax = 0 means no shock viscosity',2)
    if (alphamax < 1.) call warn(label,'alphamax < 1 is dangerous if there are shocks: don''t publish crap',2)
    if (alphamax < 0. .or. alphamax > 100.) call fatal(label,'stupid value for alphamax (generally 0.0-1.0)')
    if (mhd) then
       if (alphaB < 0.)   call warn(label,'stupid choice of alphaB',4)
       if (alphaB > 10.)  call warn(label,'very large alphaB, check timestep',3)
       if (alphaB < 1.)   call warn(label,'alphaB < 1 is not recommended, please don''t publish rubbish',2)
       if (psidecayfac < 0.) call fatal(label,'stupid value for psidecayfac')
       if (psidecayfac > 2.) call warn(label,'psidecayfac set outside recommended range (0.1-2.0)')
       if (overcleanfac < 1.0) call warn(label,'overcleanfac less than 1')
    endif
    if (beta < 0.)     call fatal(label,'beta < 0')
    if (beta > 4.)     call warn(label,'very high beta viscosity set')
    if (.not.compiled_with_mcfost) then
       if (maxvxyzu >= 4 .and. eos_requires_isothermal(ieos)) &
          call fatal(label,'storing thermal energy but eos choice requires ISOTHERMAL=yes')
    endif
    if (irealvisc < 0 .or. irealvisc > 12)  call fatal(label,'invalid setting for physical viscosity')
    if (shearparam < 0.)                     call fatal(label,'stupid value for shear parameter (< 0)')
    if (irealvisc==2 .and. shearparam > 1) call error(label,'alpha > 1 for shakura-sunyaev viscosity')
    if (iverbose > 99 .or. iverbose < -9)   call fatal(label,'invalid verboseness setting (two digits only)')

    if (icooling > 0 .and. .not.(ieos == 2 .or. ieos == 5 .or. ieos == 17 .or. ieos == 22 .or. ieos == 24)) &
         call fatal(label,'cooling requires adiabatic eos (ieos=2)')
    if (icooling > 0 .and. (ipdv_heating <= 0 .or. ishock_heating <= 0)) &
         call fatal(label,'cooling requires shock and work contributions')
    if (((isink_radiation == 1 .or. isink_radiation == 3 ) .and. idust_opacity == 0 ) &
       .and. alpha_rad < 1.d-10 .and. itau_alloc == 0) &
         call fatal(label,'no radiation pressure force! adapt isink_radiation/idust_opacity/alpha_rad')
    if (isink_radiation > 1 .and. idust_opacity == 0 ) &
         call fatal(label,'dust opacity not used! change isink_radiation or idust_opacity')
    if (iget_tdust > 2 .and. iray_resolution < 0 ) &
         call fatal(label,'To get dust temperature with Attenuation or Lucy, set iray_resolution >= 0')
    if (do_nucleation .and. ieos == 5) call error(label,'with nucleation you must use ieos = 2')
 endif
 return

999 continue
 if (id == master) then
    call error('Phantom','input file '//trim(infile)//' not found')
    if (adjustl(infile(1:1)) /= ' ') then
       write(*,*) ' creating one with default options...'
       infilenew = trim(infile)
       call write_infile(infilenew,logfile,evfile,dumpfile,iwritein,iprint)
    endif
 endif
 call die

end subroutine read_infile

end module readwrite_infile
