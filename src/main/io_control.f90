!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module io_control
!
! This module contains options and utility routines controlling
!  run time and the decisions about input/output
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io
!
 use io,       only:iverbose
 use timestep, only:tmax,dtmax
 implicit none

 integer, public :: nfulldump,nmaxdumps
 integer, public :: nmax,nout

 real(kind=4), public :: twallmax

 ! final maximum density
 real, public :: rhofinal_cgs,rhofinal1

 public :: write_options_iocontrol,read_options_iocontrol,set_defaults_iocontrol
 public :: at_simulation_end,check_for_full_dump

 private

contains

!-----------------------------------------------------------------------
!+
!  Set default io control options
!+
!-----------------------------------------------------------------------
subroutine set_defaults_iocontrol()
 use dynamic_dtmax, only:set_defaults_dynamic_dtmax

 tmax    = 10.0
 dtmax   =  1.0

 nmaxdumps = -1
 twallmax  = 0.0             ! maximum wall time for run, in seconds
 rhofinal_cgs = 0.           ! Final maximum density (0 == ignored)


 nmax = -1
 nout = -1
 nfulldump = 10              ! frequency of writing full dumps
 iverbose = 0

 ! Default dynamic dtmax options
 call set_defaults_dynamic_dtmax

end subroutine set_defaults_iocontrol

!-----------------------------------------------------------------------
!+
!  Write io control options to .in file
!+
!-----------------------------------------------------------------------
subroutine write_options_iocontrol(iunit)
 use infile_utils,  only:write_inopt
 use dynamic_dtmax, only:write_options_dynamic_dtmax,dtmax_user
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling run time and input/output'
 if (dtmax_user < 0.) dtmax_user = dtmax ! this should only ever be true for phantomsetup
 call write_inopt(tmax,'tmax','end time',iunit)
 call write_inopt(dtmax_user,'dtmax','time between dumps',iunit)
 call write_inopt(nmax,'nmax','maximum number of timesteps (0=just get derivs and stop)',iunit)
 if (rhofinal_cgs > 0.0) call write_inopt(rhofinal_cgs,'rhofinal_cgs','maximum allowed density (cgs) (<=0 to ignore)',iunit)
 if (nmaxdumps > 0) call write_inopt(nmaxdumps,'nmaxdumps','stop after n full dumps (-ve=ignore)',iunit)

 call write_inopt(nout,'nout','write dumpfile every n dtmax (-ve=ignore)',iunit)
 call write_inopt(real(twallmax),'twallmax','maximum wall time (hhh:mm, 000:00=ignore)',iunit,time=.true.)
 call write_inopt(nfulldump,'nfulldump','full dump every n dumps',iunit)
 call write_inopt(iverbose,'iverbose','verboseness of log (-1=quiet 0=default 1=allsteps 2=debug 5=max)',iunit)

 ! options for dynamic dtmax
 call write_options_dynamic_dtmax(iunit)

end subroutine write_options_iocontrol

!-----------------------------------------------------------------------
!+
!  Read io control options from .in file
!+
!-----------------------------------------------------------------------
subroutine read_options_iocontrol(name,valstring,imatch,igotall,ierr)
 use io,            only:fatal,warn
 use dynamic_dtmax, only:read_options_dynamic_dtmax
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer, intent(out) :: ierr
 logical :: igotalldtmax
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_infile'

 imatch = .true.
 igotall = .true. ! default to true for optional parameters
 igotalldtmax = .true.

 select case(trim(name))
 case('tmax')
    read(valstring,*,iostat=ierr) tmax  
    ngot = ngot + 1
 case('dtmax')
    read(valstring,*,iostat=ierr) dtmax
    ngot = ngot + 1
    if (dtmax > tmax) call warn(label,'no output dtmax > tmax')
 case('nmax')
    read(valstring,*,iostat=ierr) nmax
 case('nout')
    read(valstring,*,iostat=ierr) nout
    if (nout > nmax)  call warn(label,'no output nout > nmax')
    if (nout==0)      call fatal(label,'nout = 0')   
 case('nmaxdumps')
    read(valstring,*,iostat=ierr) nmaxdumps
 case('nfulldump')
    read(valstring,*,iostat=ierr) nfulldump
    if (nfulldump==0 .or. nfulldump > 10000) call fatal(label,'nfulldump = 0')
    if (nfulldump >= 50) call warn(label,'no full dumps for a long time...')   
 case('twallmax')
    read(valstring,*,iostat=ierr) twallmax
    if (twallmax < 0.) call fatal(label,'invalid twallmax (use 000:00 to ignore)')
 case('rhofinal_cgs')
    read(valstring,*,iostat=ierr) rhofinal_cgs
 case('iverbose')
    read(valstring,*,iostat=ierr) iverbose
    if (iverbose > 99 .or. iverbose < -9) call fatal(label,'invalid verboseness setting (two digits only)')
 case default
    imatch = .false.
    if (.not.imatch) call read_options_dynamic_dtmax(name,valstring,imatch,igotalldtmax,dtmax,ierr)
 end select

 igotall = (ngot >= 2) .and. igotalldtmax

end subroutine read_options_iocontrol

!----------------------------------------------------------------
!+
!  check if the simulation has ended
!+
!----------------------------------------------------------------
logical function at_simulation_end(time,nsteps,rhomaxnow)
 integer, intent(in) :: nsteps
 real,    intent(in) :: time,rhomaxnow

 at_simulation_end = (time >= tmax) .or. ((nsteps >= nmax).and.(nmax >= 0)) &
                                    .or. (rhomaxnow*rhofinal1 >= 1.0)

end function at_simulation_end

!----------------------------------------------------------------
!+
!  check if the simulation has ended
!+
!----------------------------------------------------------------
function check_for_full_dump(time,rhomaxnow,nsteps,noutput,&
                             twallused,twallperdump,abortrun) result(fulldump)
 use io, only:id,master,iprint
 integer,      intent(in)  :: nsteps,noutput
 real,         intent(in)  :: time,rhomaxnow
 real(kind=4), intent(in)  :: twallused,twallperdump
 logical,      intent(out) :: abortrun
 logical :: fulldump

 fulldump = (nout <= 0 .and. mod(noutput,nfulldump)==0) .or. (mod(noutput,nout*nfulldump)==0)
!
!--if max wall time is set (> 1 sec) stop the run at the last full dump
!  that will fit into the walltime constraint, based on the wall time between
!  the last two dumps added to the current total walltime used.  The factor of three for
!  changing to full dumps is to account for the possibility that the next step will take longer.
!  If we are about to write a small dump but it looks like we won't make the next dump,
!  write a full dump instead and stop the run
!
 abortrun = .false.
 if (twallmax > 1.) then
    if (fulldump) then
       if ((twallused + abs(nfulldump)*twallperdump) > twallmax) then
          abortrun = .true.
       endif
    else
       if ((twallused + 3.0*twallperdump) > twallmax) then
          fulldump = .true.
          if (id==master) write(iprint,"(1x,a)") '>> PROMOTING DUMP TO FULL DUMP BASED ON WALL TIME CONSTRAINTS... '
          nfulldump = 1  !  also set all future dumps to be full dumps (otherwise gets confusing)
          if ((twallused + twallperdump) > twallmax) abortrun = .true.
       endif
    endif
 endif
!
!--Promote to full dump if this is the final dump
!
 if (at_simulation_end(time,nsteps,rhomaxnow)) fulldump = .true.

end function check_for_full_dump

end module io_control