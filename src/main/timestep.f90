!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
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
 integer :: idtmax_n
 integer :: idtmax_n_next
 integer :: idtmax_frac
 integer :: idtmax_frac_next
 real    :: dtmax_user

 real    :: dt,dtcourant,dtforce,dtrad,dtextforce,dterr,dtdiff,time
 real    :: dtmax_dratio, dtmax_max, dtmax_min, rhomaxnow
 real(kind=4) :: dtwallmax
 integer :: dtmax_ifactor,dtmax_ifactorWT

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

 dtwallmax = 86400.          ! maximum wall time between dumps (seconds); will create 'restart' dumps as required

 ! Values to control dtmax changing with increasing densities
 dtmax_dratio =  0.          ! dtmax will change if this ratio is exceeded in a timestep (recommend 1.258)
 dtmax_max    = -1.0         ! maximum dtmax allowed (to be reset to dtmax if = -1)
 dtmax_min    =  0.          ! minimum dtmax allowed

 idtmax_n = 1
 idtmax_n_next = 1
 idtmax_frac = 0
 idtmax_frac_next = 0
 dtmax_user = -1.

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
!  The following two options are toggled in the .in file and are optional:
!  1. there is a big change in gas density (dtmax itself will change)
!  2. there is too long of walltime between dumps; dtmax will change
!     internally to write restart dumps, but dumps will still be produced
!     at the original frequency
!  Note: when dtmax changes or subdumps are created, all dumps are
!  promoted to full dumps
!+
!----------------------------------------------------------------
subroutine check_dtmax_for_decrease(iprint,dtmax,twallperdump,dtmax_log_dratio,&
                                    rhomaxold,rhomaxnew,nfulldump,change_dtmax_now)
 use io, only: iverbose
 integer,      intent(in)    :: iprint
 integer,      intent(inout) :: nfulldump
 real,         intent(inout) :: dtmax,rhomaxold
 real,         intent(in)    :: rhomaxnew,dtmax_log_dratio
 real(kind=4), intent(in)    :: twallperdump
 logical,      intent(in)    :: change_dtmax_now
 real                        :: ratio,dtmax_global,tempvar,diff
 integer                     :: ipower,ifactor
 integer, parameter          :: ifactor_max_dn = 2**2 ! hardcode to allow at most a decrease of 2 bins per step
 integer, parameter          :: ifactor_max_up = 2**1 ! hardcode to allow at most an increase of 1 bin per step

 ! initialise variables
 dtmax_ifactor   = 0
 dtmax_ifactorWT = 0
 if (dtmax_max > 0.) then
    dtmax_global = min(dtmax_max,dtmax_user)
 else
    dtmax_global = dtmax_user ! dtmax never be the default negative value
 endif
 dtmax_global = dtmax_global - epsilon(dtmax_global) ! just to be sure that we are not accidentally increasing dtmax

 !--Modify dtmax_user based upon density evolution
 !  (algorithm copied from sphNG, with slight modifications)
 !  (this is not permitted if we are between dumps as defined by dtmax_user)
 if (dtmax_log_dratio > 0.0 .and. idtmax_frac==0) then
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
       tempvar = time/(2.0*dtmax_user)
       diff    = tempvar - int(tempvar)
       if (0.25 < diff .and. diff < 0.75) then
          if (iverbose > 0) write(iprint,'(1x,a,4es10.3)') 'modifying dtmax: Synct autochange attempt, but sync ',diff
          ipower = 0
       endif
    endif
    ifactor = 2**abs(ipower)
    if (ipower < 0) then
       ! decrease dtmax
       ifactor = min(ifactor,ifactor_max_dn)
       if (dtmax_user/ifactor >= dtmax_min ) then
          dtmax_ifactor = ifactor
          if (iverbose > 0) then
             write(iprint,'(1x,a,2(es10.3,a),2es10.3,2I6)') &
             "modifying dtmax: ",dtmax_user," -> ",dtmax_user/ifactor, &
             " due to density increase. rho_old, rho_new, power, ifactor: ", rhomaxold,rhomaxnew,ipower,ifactor
          else
             write(iprint,'(1x,a,2(es10.3,a))') "modifying dtmax: ",dtmax_user," -> ",dtmax_user/ifactor, &
             " due to density increase."
          endif

       endif
    elseif (ipower > 0) then
       ! increase dtmax
       ifactor = min(ifactor,ifactor_max_up)
       if (dtmax_user*ifactor <= dtmax_max .and. ifactor*twallperdump < dtwallmax) then
          dtmax_ifactor = -ifactor
          if (iverbose > 0) then
             write(iprint,'(1x,a,2(es10.3,a),a,2es10.3,2I6)') &
             "modifying dtmax: ",dtmax_user," -> ",dtmax_user*ifactor," due to density decrease/stabilisation. ", &
             "rho_old, rho_new, power, ifactor: ",rhomaxold,rhomaxnew,ipower,ifactor
          else
             write(iprint,'(1x,a,2(es10.3,a))') &
             "modifying dtmax: ",dtmax_user," -> ",dtmax_user*ifactor," due to density decrease/stabilisation."
          endif
       endif
    endif
    rhomaxold = rhomaxnew
    ! update dtmax_user; since this is only a diagnostic/in-out variable, we can safely
    !                    update it here for both global & individual timestepping
    if (dtmax_ifactor > 0) then
       dtmax_user =  dtmax_user/dtmax_ifactor
    elseif (dtmax_ifactor < 0) then
       dtmax_user = -dtmax_user*dtmax_ifactor
    endif
 endif

!--Modify dtmax based upon wall time constraint, if requested
!  we will not try this is dtmax has just been modified due to density
 if ( dtwallmax > 0.0 .and. dtmax_ifactor==0) then
    if (twallperdump > dtwallmax) then
       dtmax_ifactor = int(2**(int(log(real(twallperdump/dtwallmax))/log(2.0))+1))
       write(iprint,'(1x,a,I4,a)') &
          "modifying dtmax internally due to wall time constraint.  Increasing to ",idtmax_n*dtmax_ifactor," sub-dumps"
       ! set nfulldump = 1 to ensure a full dump within a reasonable wall-time
       if (nfulldump > 1) then
          nfulldump = 1
          write(iprint,'(1x,a)')  &
             "modifying dtmax: nfulldump -> 1 to ensure data is not lost due to decreasing dtmax"
       endif
    elseif (twallperdump < 0.5*dtwallmax .and. idtmax_n > 1) then
       ! let's increase dtmax only by a factor two, despite the possibility of increasing it by more
       if (idtmax_frac==0 .or. idtmax_frac*2==idtmax_n) then
          dtmax_ifactor = -2
          write(iprint,'(1x,a,I4,a)') &
             "modifying dtmax internally due to wall time constraint.  Decreasing to ",-idtmax_n/dtmax_ifactor," sub-dumps"
       endif
    endif
    dtmax_ifactorWT = dtmax_ifactor
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
