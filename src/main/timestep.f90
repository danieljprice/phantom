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
! :Dependencies: dim, infile_utils
!
 implicit none
 real    :: tmax,dtmax
 real    :: C_cour,C_force,C_cool,C_rad,tolv,xtol,ptol
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
subroutine read_options_timestep(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use dim,          only:gr,mhd
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(C_cour,'C_cour',db,min=0.,max=1.,errcount=nerr)
 call read_inopt(C_force,'C_force',db,errcount=nerr,default=C_force)
 call read_inopt(tolv,'tolv',db,min=0.,max=0.1,errcount=nerr)
 if (gr) then
    call read_inopt(xtol,'xtol',db,min=0.,max=0.1,errcount=nerr)
    call read_inopt(ptol,'ptol',db,min=0.,max=0.1,errcount=nerr)
 endif
 if (mhd) then
    call read_inopt(psidecayfac,'psidecayfac',db,min=0.,max=2.,errcount=nerr,default=psidecayfac)
    call read_inopt(overcleanfac,'overcleanfac',db,min=1.,errcount=nerr,default=overcleanfac)
 endif

end subroutine read_options_timestep

end module timestep
