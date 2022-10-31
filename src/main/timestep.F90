!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module timestep
!
! Options and utility routines related to timestepping
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io
!
 implicit none
 real    :: tmax,dtmax
 real    :: C_cour,C_force,C_cool,C_rad,C_ent,tolv,xtol,ptol
 integer :: nmax,nout
 integer :: nsteps
 real, parameter :: bignumber = 1.e29
 real    :: dtmax0 = -1.  ! save the original dtmax; require when modifying dtmax due to density or time constraints

 real    :: dt,dtcourant,dtforce,dtrad,dtextforce,dterr,dtdiff,time
 real    :: dtmax_dratio, dtmax_max, dtmax_min, rhomaxnow
 real(kind=4) :: dtwallmax
 integer :: dtmax_ifactor
 logical :: allow_increase = .false.  ! to permit dtmax to increase based upon walltime

 public

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
 tmax    = 10.0
 dtmax   =  1.0
 tolv    = 1.e-2
 xtol    = 1.e-7
 ptol    = 1.e-7
 nmax    = -1
 nout    = -1

 dtwallmax = 43200.0         ! maximum wall time between dumps (seconds); default = 12h

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
!----------------------------------------------------------------
!+
!  This will determine if dtmax needs to be decreased
!  This subroutine is called at different times depending on
!  whether individual or global timesteps are being used.
!+
!----------------------------------------------------------------
subroutine check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_ifactor,dtmax_log_dratio,&
                                    rhomaxold,rhomaxnew,nfulldump,change_dtmax_now)
 use io, only: iverbose
 integer,      intent(in)    :: iprint
 integer,      intent(out)   :: dtmax_ifactor
 integer,      intent(inout) :: nfulldump
 real,         intent(inout) :: dtmax,rhomaxold
 real,         intent(in)    :: rhomaxnew,dtmax_log_dratio
 real(kind=4), intent(in)    :: twallperdump
 logical,      intent(in)    :: change_dtmax_now
 real                        :: ratio,dtmax_global,tempvar,diff
 integer                     :: ipower,ifactor,dtmax_ifactor_time
 integer, parameter          :: ifactor_max_dn = 2**2 ! hardcode to allow at most a decrease of 2 bins per step
 integer, parameter          :: ifactor_max_up = 2**1 ! hardcode to allow at most an increase of 1 bin per step

 ! initialise variables
 dtmax_ifactor      = 0
 dtmax_ifactor_time = 0
 if (dtmax_max > 0.) then
    dtmax_global = min(dtmax_max,dtmax0)
 else
    dtmax_global = dtmax0 ! dtmax never be the default negative value
 endif
 dtmax_global = dtmax_global - epsilon(dtmax_global) ! just to be sure that we are not accidentally increasing dtmax

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
    elseif (twallperdump < 0.5*dtwallmax .and. dtmax < dtmax_global .and. allow_increase) then
       dtmax_ifactor_time = 2
       do while(dtmax_ifactor_time*twallperdump < 0.5*dtwallmax .and. dtmax*dtmax_ifactor_time < dtmax_global)
          dtmax_ifactor_time = dtmax_ifactor_time * 2
       enddo
       write(iprint,'(1x,a,2(es10.3,a))') &
          "modifying dtmax: ",dtmax," -> ",dtmax*dtmax_ifactor_time," due to wall time constraint"
       dtmax_ifactor_time = -dtmax_ifactor_time
    endif
 endif

 ! modify dtmax based upon density change (algorithm copied from sphNG, with slight modifications)
 if (dtmax_log_dratio > 0.0) then
    ratio   = log10(rhomaxnew/rhomaxold)
    ipower  = -(int(ratio/dtmax_log_dratio))
    if (abs(ratio/dtmax_log_dratio) < 0.5) ipower = 1
    if (iverbose > 0) then
       write(iprint,'(1x,a,4es10.3,I6)') &
       "modifying dtmax: inspecting ratio rho_new/rho_old, rho_old, rho_new, ipower: ", &
       10**dtmax_log_dratio,rhomaxnew/rhomaxold,rhomaxold,rhomaxnew,ipower
    endif

    if (ipower > 5) ipower = 5 ! limit the largest increase in step size to 2**5 = 32
    if (ipower == 1) then
       tempvar = time/(2.0*dtmax)
       diff    = tempvar - int(tempvar)
       if (0.25 < diff .and. diff < 0.75) then
          write(iprint,'(1x,a,4es10.3)') 'modifying dtmax: Synct autochange attempt, but sync ',diff
          ipower = 0
       endif
    endif
    ifactor = 2**abs(ipower)
    if (ipower < 0) then
       ! decrease dtmax
       ifactor = min(ifactor,ifactor_max_dn)
       if (dtmax/ifactor >= dtmax_min ) then
          dtmax_ifactor = ifactor
          if (iverbose > 0) then
             write(iprint,'(1x,a,2(es10.3,a),2es10.3,2I6)') &
             "modifying dtmax: ",dtmax," -> ",dtmax/ifactor," due to density increase. rho_old, rho_new, power, ifactor: ", &
             rhomaxold,rhomaxnew,ipower,ifactor
          else
             write(iprint,'(1x,a,2(es10.3,a))') "modifying dtmax: ",dtmax," -> ",dtmax/ifactor," due to density increase."
          endif

       endif
    elseif (ipower > 0) then
       ! increase dtmax
       ifactor = min(ifactor,ifactor_max_up)
       if (dtmax*ifactor <= dtmax_max .and. ifactor*twallperdump < dtwallmax) then
          dtmax_ifactor = -ifactor
          if (iverbose > 0) then
             write(iprint,'(1x,a,2(es10.3,a),a,2es10.3,2I6)') &
             "modifying dtmax: ",dtmax," -> ",dtmax*ifactor," due to density decrease/stabilisation. ", &
             "rho_old, rho_new, power, ifactor: ",rhomaxold,rhomaxnew,ipower,ifactor
          else
             write(iprint,'(1x,a,2(es10.3,a))') &
             "modifying dtmax: ",dtmax," -> ",dtmax*ifactor," due to density decrease/stabilisation."
          endif
       endif
    endif
    rhomaxold = rhomaxnew
 endif

 ! Compare modifications based upon time & density; select largest decrease and/or smallest increase
 if (dtmax_ifactor == 0) then
    dtmax_ifactor = dtmax_ifactor_time
 else
    if (dtmax_ifactor_time /= 0) then
       dtmax_ifactor = max(dtmax_ifactor,dtmax_ifactor_time)
    endif
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
