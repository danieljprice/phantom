!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dynamic_dtmax
!
! Module to handle changing the time between dumps
! dynamically based on:
!
!  1. the maximum wall time between dumps
!  2. large density changes, slowing the time between dumps
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - dtmax_dratio : *density ratio controlling decrease (<=0 to ignore)*
!   - dtmax_max    : *maximum allowed dtmax (=dtmax if <= 0)*
!   - dtmax_min    : *minimum allowed dtmax*
!   - dtwallmax    : *maximum wall time between dumps (hhh:mm, 000:00=ignore)*
!
! :Dependencies: infile_utils, io
!
 implicit none
 integer, public :: idtmax_n
 integer, public :: idtmax_frac
 integer, public :: idtmax_n_next
 integer, public :: idtmax_frac_next
 integer, public :: dtmax_ifactor,dtmax_ifactorWT

 real, public :: dtmax_user,dtmax_dratio

 real, public :: dtmax_max,dtmax_min
 real :: dtmax_log_dratio
 real(kind=4) :: dtwallmax

 public :: write_options_dynamic_dtmax, read_options_dynamic_dtmax
 public :: set_defaults_dynamic_dtmax,check_dtmax_for_decrease,check_for_restart_dump
 public :: get_dtmax_initial

 private

contains

!-----------------------------------------------------------------------
!+
!  Set the default values for the dynamic dtmax parameters
!+
!-----------------------------------------------------------------------
subroutine set_defaults_dynamic_dtmax

 dtwallmax = 86400.          ! maximum wall time between dumps (seconds); will create 'restart' dumps as required

 ! Values to control dtmax changing with increasing densities
 dtmax_dratio =  0.          ! dtmax will change if this ratio is exceeded in a timestep (recommend 1.258)
 dtmax_log_dratio = 0.0
 dtmax_max    = -1.0         ! maximum dtmax allowed (to be reset to dtmax if = -1)
 dtmax_min    =  0.          ! minimum dtmax allowed

 idtmax_n = 1
 idtmax_n_next = 1
 idtmax_frac = 0
 idtmax_frac_next = 0
 dtmax_user = -1.

end subroutine set_defaults_dynamic_dtmax

!-----------------------------------------------------------------------
!+
!  Write the options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_dynamic_dtmax(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 if (dtwallmax > 0.0) then
    call write_inopt(real(dtwallmax),'dtwallmax','maximum wall time between dumps (hhh:mm, 000:00=ignore)',iunit,time=.true.)
 endif
 if (dtmax_dratio > 1.0) then
    call write_inopt(dtmax_dratio,'dtmax_dratio','density ratio controlling decrease (<=0 to ignore)',iunit)
    call write_inopt(dtmax_max,'dtmax_max','maximum allowed dtmax (=dtmax if <= 0)',iunit)
    call write_inopt(dtmax_min,'dtmax_min','minimum allowed dtmax',iunit)
 endif

end subroutine write_options_dynamic_dtmax

!-----------------------------------------------------------------------
!+
!  Read the options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_dynamic_dtmax(db,nerr,dtmax)
 use io, only:fatal
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 real,         intent(in)    :: dtmax
 integer :: ierr
 character(len=*), parameter :: label = 'read_options'
 real :: ratio

 ! none of these are compulsory
 call read_inopt(dtwallmax,'dtwallmax',db,min=0._4,errcount=nerr,default=0._4)
 call read_inopt(dtmax_dratio,'dtmax_dratio',db,errcount=nerr,default=dtmax_dratio)
 call read_inopt(dtmax_max,'dtmax_max',db,ierr,errcount=nerr,default=dtmax_max)
 if (ierr == 0) then
    if (dtmax_max <= 0.0) dtmax_max = dtmax
    ! to prevent comparison errors from round-off
    ratio = dtmax_max/dtmax
    ratio = int(ratio+0.5)+0.0001
    dtmax_max = dtmax*ratio
 endif
 call read_inopt(dtmax_min,'dtmax_min',db,ierr,errcount=nerr,default=dtmax_min)
 if (ierr == 0) then
    ! to prevent comparison errors from round-off
    if (dtmax_min > epsilon(dtmax_min)) then
       ratio = dtmax/dtmax_min
       ratio = int(ratio+0.5)+0.0001
       dtmax_min = dtmax/ratio
    endif
 endif

end subroutine read_options_dynamic_dtmax

!-----------------------------------------------------------------------
!+
!  Initialise the dtmax parameters
!+
!-----------------------------------------------------------------------
subroutine get_dtmax_initial(dtmax)
 real, intent(inout) :: dtmax

 dtmax_user = dtmax           ! the user defined dtmax
 if (idtmax_n < 1) idtmax_n = 1
 dtmax = dtmax/idtmax_n  ! dtmax required to satisfy the walltime constraints

end subroutine get_dtmax_initial

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
subroutine check_dtmax_for_decrease(iprint,time,dtmax,twallperdump,rhomaxold,rhomaxnew,nfulldump,change_dtmax_now)
 use io, only:iverbose
 integer,      intent(in)    :: iprint
 integer,      intent(inout) :: nfulldump
 real,         intent(in)    :: time
 real,         intent(inout) :: dtmax,rhomaxold
 real,         intent(in)    :: rhomaxnew
 real(kind=4), intent(in)    :: twallperdump
 logical,      intent(in)    :: change_dtmax_now
 real                        :: ratio,dtmax_global,tempvar,diff,dtmax_log_dratio
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

 if (dtmax_dratio > 0.) then
    dtmax_log_dratio = log10(dtmax_dratio)
 else
    dtmax_log_dratio = 0.0
 endif

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

!-----------------------------------------------------------------
!+
!  routine to check for a restart dump
!+
!-----------------------------------------------------------------
subroutine check_for_restart_dump()

 if (dtmax_ifactorWT == 0) then
    idtmax_n_next    =  idtmax_n
    idtmax_frac_next =  idtmax_frac
 elseif (dtmax_ifactorWT > 0) then
    idtmax_n_next    =  idtmax_n   *dtmax_ifactorWT
    idtmax_frac_next =  idtmax_frac*dtmax_ifactorWT
 elseif (dtmax_ifactorWT < 0) then
    idtmax_n_next    = -idtmax_n   /dtmax_ifactorWT
    idtmax_frac_next = -idtmax_frac/dtmax_ifactorWT
 endif
 idtmax_frac_next = idtmax_frac_next + 1
 idtmax_frac_next = mod(idtmax_frac_next,idtmax_n_next)

end subroutine check_for_restart_dump

end module dynamic_dtmax
