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
!   - calc_erot : *include E_rot in the ev_file*
!   - curlv     : *output curl v in dump files*
!   - dumpfile  : *dump file to start from*
!   - hfact     : *h in units of particle spacing [h = hfact(m/rho)^(1/3)]*
!   - logfile   : *file to which output is directed*
!   - tolh      : *tolerance on h-rho iterations*
!   - track_lum : *write du/dt to dump files (for a "lightcurve")*
!
! :Dependencies: HIIRegion, boundary_dyn, cooling, damping, dim, dust,
!   dust_formation, eos, externalforces, forcing, gravwaveutils, growth,
!   infile_utils, injection, io, io_control, mcfost_utils, metric,
!   neighkdtree, nicil_sup, options, part, porosity, ptmass,
!   ptmass_radiation, radiation_utils, shock_capturing, timestep,
!   utils_apr, viscosity
!
 use options,   only:iexternalforce
 use part,      only:hfact,tolh
 use dim,       only:do_radiation,nucleation,use_dust,use_dustgrowth,mhd_nonideal,compiled_with_mcfost,&
                     inject_parts,curlv,driving,track_lum,disc_viscosity,isothermal
 implicit none

contains

!-----------------------------------------------------------------
!+
!  writes an input file
!+
!-----------------------------------------------------------------
subroutine write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
 use timestep,         only:write_options_timestep
 use io,               only:fatal
 use infile_utils,     only:write_inopt
 use forcing,          only:write_options_forcing
 use externalforces,   only:write_options_externalforces
 use damping,          only:write_options_damping
 use neighkdtree,      only:write_options_tree
 use dust,             only:write_options_dust
 use growth,           only:write_options_growth
 use porosity,         only:write_options_porosity
 use injection,        only:write_options_injection
 use utils_apr,        only:write_options_apr
 use dust_formation,   only:write_options_dust_formation
 use nicil_sup,        only:write_options_nicil
 use metric,           only:write_options_metric
 use eos,              only:write_options_eos
 use ptmass,           only:write_options_ptmass
 use ptmass_radiation, only:write_options_ptmass_radiation
 use cooling,          only:write_options_cooling
 use gravwaveutils,    only:write_options_gravitationalwaves
 use radiation_utils,  only:write_options_radiation
 use dim,              only:maxvxyzu,maxptmass,gravity,sink_radiation,gr,use_apr
 use part,             only:mhd,nptmass
 use boundary_dyn,     only:write_options_boundary
 use HIIRegion,        only:write_options_H2R
 use viscosity,        only:write_options_viscosity
 use mcfost_utils,     only:write_options_mcfost
 use shock_capturing,  only:write_options_shock_capturing
 use io_control,       only:write_options_iocontrol
 use options,          only:write_options_output
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

 call write_options_iocontrol(iwritein)

 write(iwritein,"(/,a)") '# options controlling accuracy'
 call write_options_timestep(iwritein)
 call write_options_tree(iwritein)
 call write_inopt(hfact,'hfact','h in units of particle spacing [h = hfact(m/rho)^(1/3)]',iwritein)
 call write_inopt(tolh,'tolh','tolerance on h-rho iterations',iwritein,exp=.true.)

 call write_options_shock_capturing(iwritein)
 call write_options_damping(iwritein)
 !
 ! thermodynamics
 !
 call write_options_eos(iwritein)

 if (maxvxyzu >= 4) call write_options_cooling(iwritein)

 if (compiled_with_mcfost) call write_options_mcfost(iwritein)

 ! only write sink options if they are used, or if self-gravity is on
 if (nptmass > 0 .or. gravity) call write_options_ptmass(iwritein)

 call write_options_externalforces(iwritein,iexternalforce)

 if (.not. disc_viscosity) call write_options_viscosity(iwritein)

 if (driving) call write_options_forcing(iwritein)
 if (use_dust) call write_options_dust(iwritein)
 if (use_dustgrowth) then
    call write_options_growth(iwritein)
    call write_options_porosity(iwritein)
 endif

 call write_options_injection(iwritein)

 if (nucleation) call write_options_dust_formation(iwritein)

 if (sink_radiation) call write_options_ptmass_radiation(iwritein)

 if (mhd_nonideal) call write_options_nicil(iwritein)

 if (do_radiation) call write_options_radiation(iwritein)

 if (gr) call write_options_metric(iwritein)
 call write_options_boundary(iwritein)

 if (use_apr) call write_options_apr(iwritein)

 call write_options_H2R(iwritein)

 write(iwritein,"(/,a)") '# optional outputs'
 call write_options_output(iwritein)
 call write_options_gravitationalwaves(iwritein)

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
                           itau_alloc,gr,do_nucleation,use_apr
 use timestep,        only:read_options_timestep
 use eos,             only:read_options_eos,ieos,eos_requires_isothermal
 use io,              only:ireadin,iwritein,iprint,warn,die,error,fatal,id,master,fileprefix
 use infile_utils,    only:read_next_inopt,contains_loop,write_infile_series
 use forcing,         only:read_options_forcing,write_options_forcing
 use externalforces,  only:read_options_externalforces
 use neighkdtree,     only:read_options_tree
 use dust,            only:read_options_dust
 use growth,          only:read_options_growth
 use options,         only:use_porosity
 use porosity,        only:read_options_porosity
 use metric,          only:read_options_metric
 use injection,       only:read_options_injection
 use utils_apr,       only:read_options_apr
 use dust_formation,  only:read_options_dust_formation,idust_opacity
 use nicil_sup,       only:read_options_nicil
 use part,            only:mhd,nptmass
 use cooling,         only:read_options_cooling
 use ptmass,          only:read_options_ptmass
 use ptmass_radiation,only:read_options_ptmass_radiation,isink_radiation,&
                           alpha_rad,iget_tdust,iray_resolution
 use radiation_utils, only:read_options_radiation
 use damping,         only:read_options_damping
 use gravwaveutils,   only:read_options_gravitationalwaves
 use boundary_dyn,    only:read_options_boundary
 use HIIRegion,       only:read_options_H2R
 use viscosity,       only:read_options_viscosity
 use mcfost_utils,    only:read_options_mcfost
 use shock_capturing, only:read_options_shock_capturing
 use io_control,      only:read_options_iocontrol
 use options,         only:read_options_output
 character(len=*), parameter   :: label = 'read_infile'
 character(len=*), intent(in)  :: infile
 character(len=*), intent(out) :: logfile,evfile,dumpfile
 character(len=len(infile)+4)  :: infilenew
 character(len=10) :: cline
 character(len=20) :: name
 character(len=120) :: valstring
 integer :: ierr,ireaderr,line,idot,ngot,nlinesread
 logical :: imatch,igotallrequired,igotallturb,igotalltree,igotloops
 logical :: igotallbowen,igotallcooling,igotalldust,igotallextern,igotallinject,igotallgrowth,igotallporosity
 logical :: igotallionise,igotallnonideal,igotalleos,igotallptmass,igotalldamping,igotallapr
 logical :: igotallprad,igotalldustform,igotallgw,igotallgr,igotallbdy,igotallH2R,igotallviscosity
 logical :: igotalltimestep,igotallmcfost,igotallradiation,igotallshocks,igotalliocontrol,igotalloutput
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

 igotloops       = .false.
 igotallturb     = .true.
 igotalldust     = .true.
 igotallgrowth   = .true.
 igotallporosity = .true.
 igotalltree     = .true.
 igotallextern   = .true.
 igotallinject   = .true.
 igotallapr      = .true.
 igotalleos      = .true.
 igotallcooling  = .true.
 igotalldamping  = .true.
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
 igotallviscosity = .true.
 igotalltimestep  = .true.
 igotallmcfost    = .true.
 igotallradiation = .true.
 igotallshocks    = .true.
 igotalliocontrol = .true.
 igotalloutput    = .true.
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
    case('hfact')
       read(valstring,*,iostat=ierr) hfact
    case('tolh')
       read(valstring,*,iostat=ierr) tolh
    case default
       imatch = .false.
       if (.not.imatch) call read_options_output(name,valstring,imatch,igotalloutput,ierr)
       if (.not.imatch) call read_options_iocontrol(name,valstring,imatch,igotalliocontrol,ierr)
       if (.not.imatch) call read_options_shock_capturing(name,valstring,imatch,igotallshocks,ierr)
       if (.not.imatch .and. do_radiation) call read_options_radiation(name,valstring,imatch,igotallradiation,ierr)
       if (.not.imatch .and. compiled_with_mcfost) call read_options_mcfost(name,valstring,imatch,igotallmcfost,ierr)
       if (.not.imatch) call read_options_timestep(name,valstring,imatch,igotalltimestep,ierr)
       if (.not.imatch) call read_options_viscosity(name,valstring,imatch,igotallviscosity,ierr)
       if (.not.imatch) call read_options_externalforces(name,valstring,imatch,igotallextern,ierr,iexternalforce)
       if (.not.imatch .and. driving) call read_options_forcing(name,valstring,imatch,igotallturb,ierr)
       if (.not.imatch) call read_options_tree(name,valstring,imatch,igotalltree,ierr)
       !--Extract if one-fluid dust is used from the fileid
       if (.not.imatch .and. use_dust) call read_options_dust(name,valstring,imatch,igotalldust,ierr)
       if (.not.imatch .and. use_dustgrowth) call read_options_growth(name,valstring,imatch,igotallgrowth,ierr)
       if (.not.imatch .and. use_porosity) call read_options_porosity(name,valstring,imatch,igotallporosity,ierr)
       if (.not.imatch .and. gr) call read_options_metric(name,valstring,imatch,igotallgr,ierr)
       if (.not.imatch) call read_options_injection(name,valstring,imatch,igotallinject,ierr)
       if (.not.imatch .and. use_apr) call read_options_apr(name,valstring,imatch,igotallapr,ierr)
       if (.not.imatch .and. nucleation) call read_options_dust_formation(name,valstring,imatch,igotalldustform,ierr)
       if (.not.imatch .and. sink_radiation) call read_options_ptmass_radiation(name,valstring,imatch,igotallprad,ierr)
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

 igotallrequired = (ngot  >=  nrequired) .and. igotalltree   .and. igotallbowen   .and. igotalldust &
                    .and. igotalleos    .and. igotallcooling .and. igotallextern  .and. igotallturb &
                    .and. igotallptmass .and. igotallinject  .and. igotallionise  .and. igotallnonideal &
                    .and. igotallgrowth  .and. igotallporosity .and. igotalldamping .and. igotallprad &
                    .and. igotalldustform .and. igotallgw .and. igotallgr .and. igotallbdy .and. igotallapr &
                    .and. igotallviscosity .and. igotalltimestep .and. igotallmcfost .and. igotallradiation &
                    .and. igotallshocks .and. igotalliocontrol .and. igotalloutput

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
          if (.not.igotalltree) write(*,*) 'missing tree options'
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
          if (.not.igotallviscosity) write(*,*) 'missing viscosity options'
          if (.not.igotalltimestep) write(*,*) 'missing timestep options'
          if (.not.igotallmcfost) write(*,*) 'missing mcfost options'
          if (.not.igotallradiation) write(*,*) 'missing radiation options'
          if (.not.igotallshocks) write(*,*) 'missing shock capturing options'
          if (.not.igotalliocontrol) write(*,*) 'missing io control options'
          if (.not.igotalloutput) write(*,*) 'missing output options'
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
    if (hfact < 1. .or. hfact > 5.) call warn(label,'ridiculous choice of hfact',4)
    if (tolh > 1.e-3) call warn(label,'tolh is quite large!',2)
    if (tolh < epsilon(tolh)) call fatal(label,'tolh too small to ever converge')

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
