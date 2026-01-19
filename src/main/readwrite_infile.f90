!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
!   - dumpfile : *dump file to start from*
!   - hfact    : *h in units of particle spacing [h = hfact(m/rho)^(1/3)]*
!   - logfile  : *file to which output is directed*
!   - tolh     : *tolerance on h-rho iterations*
!
! :Dependencies: HIIRegion, boundary_dyn, cooling, damping, dim, dust,
!   dust_formation, eos, externalforces, fileutils, forcing, gravwaveutils,
!   growth, infile_utils, injection, io, io_control, mcfost_utils, metric,
!   mpiutils, neighkdtree, nicil_sup, options, part, porosity, ptmass,
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

 ! injection and related options
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
 if (iwritein /= iprint) write(iprint,"(/,a)") ' input file '//trim(infile)//' written successfully'

end subroutine write_infile

!-----------------------------------------------------------------
!+
!  reads parameters for the run from the input file
!+
!-----------------------------------------------------------------
subroutine read_infile(infile,logfile,evfile,dumpfile)
 use io,           only:ireadin,iwritein,iprint,warn,die,error,fatal,id,master,fileprefix
 use infile_utils, only:open_db_from_file,close_db,inopts,check_and_unroll_infile
 use mpiutils,     only:bcast_mpi
 character(len=*), intent(inout) :: infile
 character(len=*), intent(out) :: logfile,evfile,dumpfile
 character(len=*), parameter   :: label = 'read_infile'
 character(len=len(infile)+4)  :: infilenew
 character(len=len(dumpfile))  :: dumpfilenew
 integer :: ierr,idot,nerr,i
 logical :: igotloops,iexist
 type(inopts), allocatable :: db(:)

 ierr     = 0
 idot     = index(infile,'.in') - 1
 if (idot <= 1) idot = len_trim(infile)
 logfile    = infile(1:idot)//'01.log'
 dumpfile   = infile(1:idot)//'_00000.tmp'
 fileprefix = infile(1:idot)
 infilenew  = infile

 ! if the input file is specified as a dump file (name matches string_NNNNN and file exists)
 ! then set the infile name using the prefix of the dump file, e.g. blah_00010 -> blah.in
 dumpfilenew = ''
 if (is_dumpfile_name(fileprefix)) then
    inquire(file=trim(fileprefix),exist=iexist)
    if (iexist) then
       dumpfilenew = trim(fileprefix)
       ! prefix is the string before the last underscore
       idot = index(dumpfilenew,'_',back=.true.) - 1
       if (idot <= 1) idot = len_trim(dumpfilenew)
       fileprefix = dumpfilenew(1:idot)
       ! set the infile name to prefix.in
       infile = trim(fileprefix)//'.in'
       logfile = trim(fileprefix)//'01.log'
       evfile = trim(fileprefix)//'01.ev'
       write(*,"(a)") ' --> input file '//trim(dumpfilenew)// &
                      ' is a dump file, setting infile name to '//trim(infile)
    endif
 endif

 ! check for loops in the input file and unroll them into a series of input files
 igotloops = .false.
 if (id == master) then
    call check_and_unroll_infile(infile,igotloops,ierr)
    if (ierr /= 0) infilenew = trim(infile)//'.new'
 endif
 call bcast_mpi(igotloops)
 if (igotloops) call die

 ! open the input file as a database
 call open_db_from_file(db,infile,ireadin,ierr)

 ! if the input file is not found, create a default one
 if (ierr /= 0) then
    if (id == master) then
       call error('Phantom','input file '//trim(infile)//' not found')
       if (adjustl(infile(1:1)) /= ' ') then
          write(*,*) ' creating one with default options...'
          infilenew = trim(infile)
          call write_infile(infilenew,logfile,evfile,dumpfile,iwritein,iprint)
       endif
    endif
    call die
 endif

 call read_options_from_db(db,nerr,logfile,dumpfile,evfile)

 ! overwrite the dumpfile name if it was specified on the command line
 if (len_trim(dumpfilenew) > 0) dumpfile = dumpfilenew

 ! warn about unknown variables in the input file
 do i=1,size(db)
    if (len_trim(db(i)%tag) > 0 .and. .not.db(i)%retrieved) then
       call warn(label,'unknown variable '//trim(adjustl(db(i)%tag))// &
                 ' in input file, value = '//trim(adjustl(db(i)%val)))
    endif
 enddo

 call close_db(db)

 ! rewrite the input file if options are missing
 if (nerr > 0) then
    if (id == master) then
       call error(label,'input file '//trim(infile)//' is incomplete for current compilation')
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
 endif

end subroutine read_infile

!-----------------------------------------------------------------
!+
!  read options from the database, lower level routine that
!  can be called externally to parse a manually created database
!  into phantom options, e.g. from libphantom-amuse
!+
!-----------------------------------------------------------------
subroutine read_options_from_db(db,nerr,logfile,dumpfile,evfile)
 use dim,              only:gr,do_radiation,compiled_with_mcfost,mhd_nonideal,&
                            use_apr,sink_radiation,maxptmass,driving,use_dust,&
                            use_dustgrowth,nucleation,mhd_nonideal,maxvxyzu
 use io,               only:warn
 use infile_utils,     only:inopts,read_inopt
 use options,          only:use_porosity
 use part,             only:hfact,tolh
 use eos,              only:read_options_eos
 use io_control,       only:read_options_iocontrol
 use options,          only:read_options_output
 use shock_capturing,  only:read_options_shock_capturing
 use forcing,          only:read_options_forcing,write_options_forcing
 use externalforces,   only:read_options_externalforces
 use neighkdtree,      only:read_options_tree
 use dust,             only:read_options_dust
 use growth,           only:read_options_growth
 use porosity,         only:read_options_porosity
 use metric,           only:read_options_metric
 use injection,        only:read_options_injection
 use utils_apr,        only:read_options_apr
 use dust_formation,   only:read_options_dust_formation
 use nicil_sup,        only:read_options_nicil
 use cooling,          only:read_options_cooling
 use ptmass,           only:read_options_ptmass
 use ptmass_radiation, only:read_options_ptmass_radiation
 use radiation_utils,  only:read_options_radiation
 use damping,          only:read_options_damping
 use gravwaveutils,    only:read_options_gravitationalwaves
 use boundary_dyn,     only:read_options_boundary
 use HIIRegion,        only:read_options_H2R
 use viscosity,        only:read_options_viscosity
 use mcfost_utils,     only:read_options_mcfost
 use shock_capturing,  only:read_options_shock_capturing
 use io_control,       only:read_options_iocontrol
 use options,          only:read_options_output
 use timestep,         only:read_options_timestep
 type(inopts),     intent(inout) :: db(:)
 integer,          intent(inout) :: nerr
 character(len=*), intent(out)   :: logfile,dumpfile,evfile
 character(len=*), parameter :: label = 'read_infile'

 ! parse global options
 nerr = 0  ! errors during read: increment this for compulsory options
 call read_inopt(logfile,'logfile',db,errcount=nerr)
 call read_inopt(dumpfile,'dumpfile',db,errcount=nerr)
 evfile  = logfile2evfile(logfile)

 call read_inopt(hfact,'hfact',db,errcount=nerr,min=1.,max=5.)
 call read_inopt(tolh,'tolh',db,errcount=nerr,min=epsilon(tolh))
 if (tolh > 1.e-3) call warn(label,'tolh is quite large!')

 ! parse options internal to other code modules
 call read_options_iocontrol(db,nerr)
 call read_options_timestep(db,nerr)

 call read_options_shock_capturing(db,nerr)
 call read_options_damping(db,nerr)
 if (.not. disc_viscosity) call read_options_viscosity(db,nerr)
 !
 ! thermodynamics
 !
 call read_options_eos(db,nerr)
 if (maxvxyzu >= 4) call read_options_cooling(db,nerr)
 if (compiled_with_mcfost) call read_options_mcfost(db,nerr)

 if (maxptmass > 0) call read_options_ptmass(db,nerr)
 call read_options_externalforces(db,nerr,iexternalforce)
 call read_options_tree(db,nerr)

 if (driving) call read_options_forcing(db,nerr)
 if (use_dust) call read_options_dust(db,nerr)
 if (use_dustgrowth) call read_options_growth(db,nerr)
 if (use_porosity) call read_options_porosity(db,nerr)

 ! injection and related options
 call read_options_injection(db,nerr)
 if (nucleation) call read_options_dust_formation(db,nerr)
 if (sink_radiation) call read_options_ptmass_radiation(db,nerr)

 if (mhd_nonideal) call read_options_nicil(db,nerr)
 if (do_radiation) call read_options_radiation(db,nerr)
 if (gr) call read_options_metric(db,nerr)
 call read_options_boundary(db,nerr)

 if (use_apr) call read_options_apr(db,nerr)
 call read_options_H2R(db,nerr)
 call read_options_output(db,nerr)
 call read_options_gravitationalwaves(db,nerr)

end subroutine read_options_from_db

!-----------------------------------------------------------------
!+
!  get the name of the.ev file from the .log file
!+
!-----------------------------------------------------------------
function logfile2evfile(logfile) result(evfile)
 character(len=*), intent(in)  :: logfile
 character(len=len(logfile)+4) :: evfile
 integer :: idot

 idot = index(logfile,'.log') - 1
 if (idot <= 1) idot = len_trim(logfile)
 evfile  = logfile(1:idot)//'.ev'

end function logfile2evfile

!-----------------------------------------------------------------
!+
!  query if a file is a dump file name (snap_NNNNN)
!+
!-----------------------------------------------------------------
function is_dumpfile_name(filename) result(is_dumpfile)
 use fileutils, only:is_digit
 character(len=*), intent(in)  :: filename
 logical :: is_dumpfile
 integer :: iunderscore,ndigits,i

 is_dumpfile = .false.
 iunderscore = index(filename,'_',back=.true.)
 if (iunderscore <= 0) return
 if (index(filename,'.in') > 0) return
 if (index(filename,'.setup') > 0) return

 ndigits = 0
 do i=iunderscore+1,len_trim(filename)
    if (is_digit(filename(i:i))) then
       ndigits = ndigits + 1
    else
       exit
    endif
 enddo
 is_dumpfile = ndigits >= 5

end function is_dumpfile_name

end module readwrite_infile
