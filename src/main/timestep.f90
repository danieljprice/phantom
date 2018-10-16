!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: timestep
!
!  DESCRIPTION:
!   Options and utility routines related to timestepping
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES:
!+
!--------------------------------------------------------------------------
module timestep
 implicit none
 real    :: tmax,dtmax
 real    :: C_cour,C_force,C_cool,tolv
 integer :: nmax,nout
 integer :: nsteps
 real, parameter :: bignumber = 1.e29

 real    :: dt, dtcourant, dtforce, dtextforce, dterr, dtdiff, time
 real    :: dtmax_dratio, dtmax_max, dtmax_min, rhomaxnow
 real(kind=4) :: dtwallmax
 integer :: dtmax_ifactor
 logical :: restartonshortest

 public

contains
!-----------------------------------------------------------------
!+
!  routine to set defaults for timestepping parameters
!+
!-----------------------------------------------------------------
subroutine set_defaults_timestep

 C_cour  =  0.3
 C_force =  0.25
 C_cool  =  0.05
 tmax    = 10.0
 dtmax   =  1.0
 tolv    = 1.e-2
 nmax    = -1
 nout    = -1

 dtwallmax = 43200.0         ! maximum wall time between dumps (seconds); default = 12h
 restartonshortest = .false. ! whether or not to restart with all parts on shortest step

 ! Values to control dtmax changing with increasing densities
 dtmax_dratio =  0.          ! dtmax will change if this ratio is exceeded in a timestep (recommend 1.258)
 dtmax_max    = -1.0         ! maximum dtmax allowed (to be reset to dtmax if = -1)
 dtmax_min    =  0.          ! minimum dtmax allowed


end subroutine set_defaults_timestep

!-----------------------------------------------------------------
!+
!  routine to print out the timestep information to the log file
!+
!-----------------------------------------------------------------
subroutine print_dtlog(iprint,time,dt,dtforce,dtcourant,dterr,dtmax,dtprint,np)
 integer, intent(in) :: iprint
 real,    intent(in) :: time,dt,dtforce,dtcourant,dterr,dtmax
 real,    intent(in), optional :: dtprint
 integer, intent(in) :: np
 character(len=20) :: str
 integer, save :: nplast = 0

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
 else
    !print*,dt,dtforce,dtcourant,dterr,dtmax
    write(iprint,10) time,dt,'(unknown)'//trim(str)
 endif
10 format(' t = ',g12.5,' dt = ',es10.3,1x,a)

end subroutine print_dtlog
!----------------------------------------------------------------
!+
!  This will determine if dtmax needs to be decreased
!  This subroutine is called at different times depending on
!  whether individual or global timesteps are being used.
!+
!----------------------------------------------------------------
subroutine check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_ifactor,dtmax_log_dratio,&
                                    rhomaxold,rhomaxnew,nfulldump,change_dtmax_now)
 integer,      intent(in)    :: iprint
 integer,      intent(out)   :: dtmax_ifactor
 integer,      intent(inout) :: nfulldump
 real,         intent(inout) :: dtmax,rhomaxold
 real,         intent(in)    :: rhomaxnew,dtmax_log_dratio
 real(kind=4), intent(in)    :: twallperdump
 logical,      intent(in)    :: change_dtmax_now
 real                        :: ratio
 integer                     :: ipower,ifactor,dtmax_ifactor_time
 integer, parameter          :: ifactor_max_dn = 2**2 ! hardcode to allow at most a decrease of 2 bins per step
 integer, parameter          :: ifactor_max_up = 2**1 ! hardcode to allow at most an increase of 1 bin per step

 ! initialise variables
 dtmax_ifactor      = 0
 dtmax_ifactor_time = 0

 ! modify dtmax based upon wall time constraint, if requested
 if ( dtwallmax > 0.0 ) then
    if (twallperdump > dtwallmax) then
       dtmax_ifactor_time = int(2**(int(log(real(twallperdump/dtwallmax))/log(2.0))+1))
       write(iprint,'(1x,a,2(es10.3,a))') &
          "modifying dtmax: ",dtmax," -> ",dtmax/dtmax_ifactor_time," due to wall time constraint"
       ! set nfulldump = 1 to ensure a full dump within a reasonable wall-time
       if (nfulldump > 1) then
          nfulldump = 1
          write(iprint,'(1x,a)')  &
             "modifying dtmax: nfulldump -> 1 to ensure data is not lost due to decreasing dtmax"
       endif
    endif
 endif

 ! modify dtmax based upon density change (algorithm copied from sphNG)
 if (dtmax_log_dratio > 0.0) then
    ratio   = log10(rhomaxnew/rhomaxold)
    ipower  = -(int(ratio/dtmax_log_dratio))
    if (abs(ratio/dtmax_log_dratio) < 0.5) ipower = 1
    ifactor = 2**abs(ipower)
    if (ipower < 0) then
       ! decrease dtmax
       ifactor = min(ifactor,ifactor_max_dn)
       if (dtmax/ifactor >= dtmax_min ) then
          dtmax_ifactor = ifactor
          write(iprint,'(1x,a,2(es10.3,a))') &
          "modifying dtmax: ",dtmax," -> ",dtmax/ifactor," due to density increase"
       endif
    elseif (ipower > 0) then
       ! increase dtmax
       ifactor = ifactor_max_up
       if (dtmax*ifactor <= dtmax_max .and. ifactor*twallperdump < dtwallmax) then
          dtmax_ifactor = -ifactor
          write(iprint,'(1x,a,2(es10.3,a))') &
          "modifying dtmax: ",dtmax," -> ",dtmax*ifactor," due to density decrease/stabilisation"
       endif
    endif
    rhomaxold = rhomaxnew
 endif

 ! Decreasing due to time constraint trumps change due to density
 if (dtmax_ifactor_time > 0) then
    dtmax_ifactor = max(dtmax_ifactor,dtmax_ifactor_time)
 endif

 if (change_dtmax_now) then
    ! dtmax will be modified now for global timestepping
    if (dtmax_ifactor > 0) then
       dtmax =  dtmax/dtmax_ifactor
    elseif (dtmax_ifactor < 0) then
       dtmax = -dtmax*dtmax_ifactor
    endif
 endif

end subroutine check_dtmax_for_decrease

end module timestep
