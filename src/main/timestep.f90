!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module timestep
!
! Options and utility routines related to timestepping
! tolerances and simulation accuracy
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - C_cour       : *Courant factor*
!   - C_force      : *dt_force factor*
!   - overcleanfac : *factor to increase div B cleaning speed (decreases timestep)*
!   - psidecayfac  : *div B diffusion parameter*
!   - ptol         : *tolerance on pmom iterations*
!   - tolv         : *tolerance on v iterations in timestepping*
!   - xtol         : *tolerance on xyz iterations*
!
! :Dependencies: dim, infile_utils, io
!
 implicit none
 real    :: tmax,dtmax
 real    :: C_cour,C_force,C_cool,C_rad,C_ent,tolv,xtol,ptol
 real    :: rhomaxnow
 integer :: nsteps
 ! div B cleaning
 real    :: psidecayfac, overcleanfac

 ! internal global variables
 real    :: dt,dtcourant,dtforce,dtrad,dtextforce,dterr,dtdiff,dtinject,time

 real, parameter :: bignumber = 1.e29

 public :: write_options_timestep, read_options_timestep

contains
!-----------------------------------------------------------------
!+
!  routine to set defaults for timestepping parameters
!+
!-----------------------------------------------------------------
subroutine set_defaults_timestep

 C_cour  = 0.3
 C_force = 0.25
 C_cool  = 0.05
 C_rad   = 0.8  ! see Biriukov & Price (2019)
 C_ent   = 3.
 tolv    = 1.e-2
 xtol    = 1.e-7
 ptol    = 1.e-7

 ! div B cleaning (MHD only)
 psidecayfac       = 1.0     ! ratio of parabolic to hyperbolic cleaning
 overcleanfac      = 1.0     ! factor by which to increase cleaning speed for div B cleaning

end subroutine set_defaults_timestep

!-----------------------------------------------------------------
!+
!  routine to print out the timestep information to the log file
!+
!-----------------------------------------------------------------
subroutine print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,&
                       dtrad,dtprint,dtinj,np)
 integer,         intent(in) :: iprint
 real,            intent(in) :: time,dt,dtforce,dtcourant,dterr,dtmax,dtrad
 real,            intent(in), optional :: dtprint,dtinj
 integer(kind=8), intent(in) :: np
 character(len=20) :: str
 integer(kind=8), save :: nplast = 0

 str = ''
 if (np /= nplast) then
    nplast = np
    write(str,"(i12)") np
    str = ', np = '//trim(adjustl(str))
 endif

 if (abs(dt-dtforce) < tiny(dt)) then
    write(iprint,10) time,dt,'(force)'//trim(str)
 elseif (abs(dt-dtcourant) < tiny(dt)) then
    write(iprint,10) time,dt,'(courant)'//trim(str)
 elseif (abs(dt-dterr) < tiny(dt)) then
    write(iprint,10) time,dt,'(tolv)'//trim(str)
 elseif (abs(dt-dtmax) <= epsilon(dt)) then
    write(iprint,10) time,dt,'(dtmax)'//trim(str)
 elseif (present(dtprint) .and. abs(dt-dtprint) < tiny(dt)) then
    write(iprint,10) time,dt,'(dtprint)'//trim(str)
 elseif (present(dtinj) .and. abs(dt-dtinj) < tiny(dt)) then
    write(iprint,10) time,dt,'(dtinject)'//trim(str)
 elseif (abs(dt-dtrad) < tiny(dt)) then
    write(iprint,10) time,dt,'(radiation)'//trim(str)
 else
    !print*,dt,dtforce,dtcourant,dterr,dtmax
    write(iprint,10) time,dt,'(unknown)'//trim(str)
 endif
10 format(' t = ',g12.5,' dt = ',es10.3,1x,a)

end subroutine print_dtlog

!-----------------------------------------------------------------
!+
!  routine to write timestep accuracy options to input file
!+
!-----------------------------------------------------------------
subroutine write_options_timestep(iunit)
 use infile_utils, only:write_inopt
 use dim,          only:gr,mhd
 integer, intent(in) :: iunit

 call write_inopt(C_cour,'C_cour','Courant factor',iunit)
 call write_inopt(C_force,'C_force','dt_force factor',iunit)
 call write_inopt(tolv,'tolv','tolerance on v iterations in timestepping',iunit,exp=.true.)
 if (gr) then
    call write_inopt(C_ent,'C_ent','restrict timestep when ds/dt is too large (not used if ien_type != 3)',iunit)
    call write_inopt(xtol,'xtol','tolerance on xyz iterations',iunit)
    call write_inopt(ptol,'ptol','tolerance on pmom iterations',iunit)
 endif
 if (mhd) then
    call write_inopt(psidecayfac,'psidecayfac','div B diffusion parameter',iunit)
    call write_inopt(overcleanfac,'overcleanfac','factor to increase div B cleaning speed (decreases timestep)',iunit)
 endif

end subroutine write_options_timestep

!-----------------------------------------------------------------
!+
!  routine to read timestep accuracy options from input file
!+
!-----------------------------------------------------------------
subroutine read_options_timestep(name,valstring,imatch,igotall,ierr)
 use dim, only:mhd
 use io,  only:fatal,warn
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 character(len=30), parameter :: label = 'read_options_timestep'

 imatch  = .true.
 igotall = .true. ! default to true for optional parameters

 select case(trim(name))
 case('C_cour')
    read(valstring,*,iostat=ierr) C_cour
    if (C_cour <= 0.) call fatal(label,'Courant number < 0')
    if (C_cour > 1.)  call fatal(label,'ridiculously big courant number!!')
 case('C_force')
    read(valstring,*,iostat=ierr) C_force
    if (C_force <= 0.) call fatal(label,'bad choice for force timestep control')
 case('tolv')
    read(valstring,*,iostat=ierr) tolv
    if (tolv <= 0.)   call fatal(label,'silly choice for tolv (< 0)')
    if (tolv > 1.e-1) call warn(label,'dangerously large tolerance on v iterations')
 case('C_ent')
    read(valstring,*,iostat=ierr) C_ent
 case('xtol')
    read(valstring,*,iostat=ierr) xtol
    if (xtol <= 0.)   call fatal(label,'silly choice for xtol (< 0)')
    if (xtol > 1.e-1) call warn(label,'dangerously large tolerance on xyz iterations')
 case('ptol')
    read(valstring,*,iostat=ierr) ptol
    if (ptol <= 0.)   call fatal(label,'silly choice for ptol (< 0)')
    if (ptol > 1.e-1) call warn(label,'dangerously large tolerance on pmom iterations')
 case('psidecayfac')
    read(valstring,*,iostat=ierr) psidecayfac
    if (mhd .and. psidecayfac < 0.) call fatal(label,'stupid value for psidecayfac')
    if (mhd .and. psidecayfac > 2.) call warn(label,'psidecayfac set outside recommended range (0.1-2.0)')
case('overcleanfac')
    read(valstring,*,iostat=ierr) overcleanfac
    if (mhd .and. overcleanfac < 1.0) call warn(label,'overcleanfac less than 1')  
 case default
    imatch = .false.
 end select

end subroutine read_options_timestep

end module timestep
