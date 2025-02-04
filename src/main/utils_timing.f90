!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module timing
!
! This module contains utilities for code timings
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, mpiutils
!
 implicit none
 integer, private :: istarttime(6)
 real(kind=4), private :: starttime

 data starttime/-1./

 public :: getused,get_timings,printused
 public :: wallclock,log_timing,print_time

 public :: timer,reset_timer,increment_timer,print_timer
 public :: setup_timers,reduce_timer_mpi,reduce_timers

 integer, parameter, private :: treelabel_len = 30
 type timer
    character(len=10)            :: label
    real(kind=4)                 :: wall
    real(kind=4)                 :: cpu
    real(kind=4)                 :: loadbal
    integer                      :: parent
    integer                      :: treesymbol(5) = -1
    character(len=treelabel_len) :: treelabel
 end type timer

 integer, public, parameter ::   itimer_fromstart     = 1,  &
                                 itimer_lastdump      = 2,  &
                                 itimer_step          = 3,  &
                                 itimer_link          = 4,  &
                                 itimer_balance       = 5,  &
                                 itimer_dens          = 6,  &
                                 itimer_dens_local    = 7,  &
                                 itimer_dens_remote   = 8,  &
                                 itimer_force         = 9,  &
                                 itimer_force_local   = 10, &
                                 itimer_force_remote  = 11, &
                                 itimer_radiation     = 12, &
                                 itimer_rad_save      = 13, &
                                 itimer_rad_neighlist = 14, &
                                 itimer_rad_arrays    = 15, &
                                 itimer_rad_its       = 16, &
                                 itimer_rad_flux      = 17, &
                                 itimer_rad_diff      = 18, &
                                 itimer_rad_update    = 19, &
                                 itimer_rad_store     = 20, &
                                 itimer_cons2prim     = 21, &
                                 itimer_substep       = 22, &
                                 itimer_sinksink      = 23, &
                                 itimer_gasf          = 24, &
                                 itimer_acc           = 25, &
                                 itimer_sg_id         = 26, &
                                 itimer_sg_evol       = 27, &
                                 itimer_HII           = 28, &
                                 itimer_ev            = 29, &
                                 itimer_io            = 30
 integer, public, parameter :: ntimers = 30 ! should be equal to the largest itimer index
 type(timer), public :: timers(ntimers)

 private

contains
!--------------------------------------
!+
!  Routines to handle timer objects
!+
!--------------------------------------
subroutine setup_timers
 !
 ! These timers must be initialised with the correct tree hierarchy,
 ! i.e. children must immediately follow their parents or siblings
 !

 !               timer from array     label          parent
 call init_timer(itimer_fromstart   , 'all',         0            )
 call init_timer(itimer_lastdump    , 'last',        0            )
 call init_timer(itimer_step        , 'step',        0            )
 call init_timer(itimer_HII         , 'HII_regions', 0            )
 call init_timer(itimer_link        , 'tree',        itimer_step  )
 call init_timer(itimer_balance     , 'balance',     itimer_link  )
 call init_timer(itimer_dens        , 'density',     itimer_step  )
 call init_timer(itimer_dens_local  , 'local',       itimer_dens  )
 call init_timer(itimer_dens_remote , 'remote',      itimer_dens  )
 call init_timer(itimer_force       , 'force',       itimer_step  )
 call init_timer(itimer_force_local , 'local',       itimer_force )
 call init_timer(itimer_force_remote, 'remote',      itimer_force )
 call init_timer(itimer_radiation   , 'radiation',   itimer_step  )
 call init_timer(itimer_rad_save    , 'save',        itimer_radiation  )
 call init_timer(itimer_rad_neighlist,'neighlist',   itimer_radiation  )
 call init_timer(itimer_rad_arrays  , 'arrays',      itimer_radiation  )
 call init_timer(itimer_rad_its     , 'its',         itimer_radiation  )
 call init_timer(itimer_rad_flux    , 'flux',        itimer_rad_its    )
 call init_timer(itimer_rad_diff    , 'diff',        itimer_rad_its    )
 call init_timer(itimer_rad_update  , 'update',      itimer_rad_its    )
 call init_timer(itimer_rad_store   , 'store',       itimer_radiation  )
 call init_timer(itimer_cons2prim   , 'cons2prim',   itimer_step  )
 call init_timer(itimer_substep     , 'substep',     itimer_step  )
 call init_timer(itimer_sinksink    , 'sink-sink',   itimer_substep  )
 call init_timer(itimer_gasf        , 'gas_force',   itimer_substep  )
 call init_timer(itimer_acc         , 'accretion',   itimer_substep  )
 call init_timer(itimer_sg_id       , 'subg_id',     itimer_substep  )
 call init_timer(itimer_sg_evol     , 'subg_evol',   itimer_substep  )
 call init_timer(itimer_ev          , 'write_ev',    0            )
 call init_timer(itimer_io          , 'write_dump',  0            )

 ! When the timer is initialised, the tree structure is defined as an arary of integers
 ! because the special ASCII characters take up more than 1 space in a character array,
 ! making alignment difficult. The finish_timer_tree_symbols function is called to
 ! convert that array into a string. This makes text formatting more straightforward,
 ! and the tree diagram generation code easier to understand.
 call finish_timer_tree_symbols

end subroutine setup_timers

subroutine init_timer(itimer,label,parent)
 use io, only:fatal
 integer,          intent(in) :: itimer
 character(len=*), intent(in) :: label
 integer,          intent(in) :: parent

 integer :: i, level, previous_parent

 if (parent >= itimer) call fatal('timing', 'Attempting to initialise timer with non-existent parent')
 !
 !--Innitialise timer variables
 !
 call reset_timer(itimer)
 timers(itimer)%label   = trim(label)
 timers(itimer)%parent  = parent
 timers(itimer)%loadbal = 1._4

 !
 !--Determine ASCII characters for printing the tree
 !
 ! 0: '  '
 ! 1: '└─'
 ! 2: '├─'
 ! 3: '│ '

 ! Get timer level
 call get_timer_level(itimer,level)

 ! Set default character to '└'
 timers(itimer)%treesymbol(level) = 1

 ! Pad with spaces to the correct level
 do i = 1, level-1
    timers(itimer)%treesymbol(i) = 0
 enddo

 ! Get the parent of the previous timer
 if (itimer > 1) then
    previous_parent = timers(itimer-1)%parent
 else
    previous_parent = 0
 endif

 if (previous_parent == parent .and. itimer > 1) then
    ! If sibling is above, replace their '└' with '├'
    timers(itimer-1)%treesymbol(level) = 2
 else
    ! If something else is above, add a connecting line '│'
    ! until the branch is reached
    do i = itimer-1,parent+1,-1
       if (timers(i)%treesymbol(level) == 0) then
          ! Replace empty space with '│'
          timers(i)%treesymbol(level) = 3
       elseif (timers(i)%treesymbol(level) == 1) then
          ! Replace '└' at the branch with ├'
          timers(i)%treesymbol(level) = 2
       endif

    enddo
 endif

end subroutine init_timer

subroutine get_timer_level(itimer,level)
 integer, intent(in)  :: itimer
 integer, intent(out) :: level

 integer :: i

 ! Get the level of this timer, where level 1 is the base level
 level = 0
 i = itimer
 do while (i /= 0)
    i = timers(i)%parent
    level = level + 1
 enddo

end subroutine get_timer_level

subroutine finish_timer_tree_symbols

 character(len=treelabel_len) :: treelabel_new
 integer :: i, j, k

 do i = 1, ntimers
    treelabel_new = ''

    ! Length of tree symbol, to calculate space padding on the right
    k = 0

    ! Convert integer arrays into symbols
    do j = 1, size(timers(i)%treesymbol)
       select case (timers(i)%treesymbol(j))
       case (0)
          ! cannot use spaces because they get trimmed,
          ! so use '.' instead of ' ' as a placeholder
          treelabel_new = trim(treelabel_new) // '..'
          k = k + 2
       case (1)
          treelabel_new = trim(treelabel_new) // '└─'
          k = k + 2
       case (2)
          treelabel_new = trim(treelabel_new) // '├─'
          k = k + 2
       case (3)
          treelabel_new = trim(treelabel_new) // '│.'
          k = k + 2
       case default
          continue
       end select
    enddo

    ! Add label
    treelabel_new = trim(treelabel_new) // timers(i)%label
    k = k + len(trim(timers(i)%label))

    ! Pad with spaces
    ! Each symbol takes up 2 more characters, so to accmodate
    ! 6 symbols, remove 6*2 spaces
    do j = treelabel_len-(6*2), k, -1
       treelabel_new = trim(treelabel_new) // '.'
    enddo

    ! Add ':'
    treelabel_new = trim(treelabel_new) // ':'

    ! Replace '.' with ' '
    do j = 1, treelabel_len
       if (treelabel_new(j:j) == '.') treelabel_new(j:j) = ' '
    enddo

    timers(i)%treelabel = treelabel_new
 enddo

end subroutine finish_timer_tree_symbols

subroutine reset_timer(itimer)
 integer, intent(in)        :: itimer

 timers(itimer)%wall = 0.0_4
 timers(itimer)%cpu  = 0.0_4

end subroutine reset_timer

subroutine increment_timer(itimer,wall,cpu)
 integer,      intent(in) :: itimer
 real(kind=4), intent(in) :: wall, cpu

 timers(itimer)%wall = timers(itimer)%wall + wall
 timers(itimer)%cpu  = timers(itimer)%cpu  + cpu

end subroutine increment_timer

subroutine reduce_timers
 integer :: itimer
 do itimer = 1, ntimers
    call reduce_timer_mpi(itimer)
 enddo
end subroutine reduce_timers

subroutine reduce_timer_mpi(itimer)
 use io,       only:nprocs
 use mpiutils, only:reduceall_mpi
 integer, intent(in) :: itimer
 real(kind=4) :: mean,cpumax,cputot

 cputot = reduceall_mpi('+',timers(itimer)%cpu)

 ! load balance = average time / max time (cpu)
 ! where the average is taken over all tasks except for the max
 ! When every time takes the same time, loadbal = 1
 if (nprocs > 1) then
    cpumax = reduceall_mpi('max',timers(itimer)%cpu)
    mean   = (cputot - cpumax) / (real(nprocs,kind=4) - 1.0_4)
    if (cpumax > 0.0_4) then
       timers(itimer)%loadbal = mean / cpumax
    else
       timers(itimer)%loadbal = 0.0_4 ! to indicate an error
    endif
 else
    timers(itimer)%loadbal = 1.0_4
 endif

 timers(itimer)%cpu  = cputot
 timers(itimer)%wall = reduceall_mpi('max',timers(itimer)%wall)

end subroutine reduce_timer_mpi

!-----------------------------------------------
!+
!  Pretty-print timing entries in a nice table
!+
!-----------------------------------------------
subroutine print_timer(lu,itimer,time_total)
 integer,      intent(in) :: lu
 integer,      intent(in) :: itimer
 real(kind=4), intent(in) :: time_total

 ! Print timings
 if (timers(itimer)%wall > epsilon(0._4)) then
    ! Print label and tree structure
    write(lu, "(1x,a)", advance="no") trim(timers(itimer)%treelabel)
    if (time_total > 7200.0) then
       write(lu,"('  ',f7.2,'h   ',f8.2,'h    ',f6.2,'   ',f6.2,'%','   ',f6.2,'%')")  &
            timers(itimer)%wall/3600.,&
            timers(itimer)%cpu/3600.,&
            timers(itimer)%cpu/timers(itimer)%wall,&
            timers(itimer)%loadbal*100.,&
            timers(itimer)%wall/time_total*100.
    elseif (time_total > 120.0) then
       write(lu,"(f7.2,'min ',f8.2,'min    ',f6.2,'   ',f6.2,'%','   ',f6.2,'%')")  &
            timers(itimer)%wall/60.,&
            timers(itimer)%cpu/60.,&
            timers(itimer)%cpu/timers(itimer)%wall,&
            timers(itimer)%loadbal*100.,&
            timers(itimer)%wall/time_total*100.
    else
       write(lu,"('  ',f7.2,'s   ',f8.2,'s    ',f6.2,'   ',f6.2,'%','   ',f6.2,'%')")  &
            timers(itimer)%wall,&
            timers(itimer)%cpu,&
            timers(itimer)%cpu/timers(itimer)%wall,&
            timers(itimer)%loadbal*100.,&
            timers(itimer)%wall/time_total*100.
    endif
 endif

end subroutine print_timer

!--------------------------------------------------------------------
!+
!  Interface to print timings of deriv routine
!+
!--------------------------------------------------------------------
subroutine log_timing(label,twall,tcpu,start,iunit)
 character(len=*), intent(in) :: label
 character(len=len(label)+16) :: string,stringcpu,stringwall
 real(kind=4),     intent(in) :: twall,tcpu
 logical,          intent(in), optional :: start
 integer,          intent(in), optional :: iunit
 character(len=240), save :: logtiming,logtimingwall,logtimingcpu

 write(stringwall,"(1x,a,' =',f10.2,'s')") trim(label),twall
 write(stringcpu,"(1x,a,' =',f10.2,'s')") trim(label),tcpu
 if (twall > epsilon(twall)) then
    write(string,"(1x,a,' =',f10.2)") trim(label),tcpu/twall
 else
    write(string,"(1x,a,' = ---')") trim(label)
 endif

 if (present(start)) then
    logtiming = ')'
    logtimingwall = ')'
    logtimingcpu = ')'
 endif
 logtimingwall = trim(stringwall)//trim(logtimingwall)
 logtimingcpu = trim(stringcpu)//trim(logtimingcpu)
 logtiming = trim(string)//trim(logtiming)
 if (present(iunit)) then
    write(iunit,"(a)") ' (wall:'//trim(logtimingwall)
    write(iunit,"(a)") ' (cpu :'//trim(logtimingcpu)
    write(iunit,"(a)") ' (cpu/wall:'//trim(logtiming)
 endif

 return
end subroutine log_timing

!--------------------------------------------------------------------
!+
!  sets initial time
!+
!--------------------------------------------------------------------
subroutine initialise_timing
 integer :: iday,imonth,iyear,ihour,imin,isec,imsec,ivalues(8)
 character(len=8)  :: date
 character(len=5)  :: zone
 character(len=10) :: time

 call date_and_time(date,time,zone,ivalues)
 iyear  = ivalues(1)
 imonth = ivalues(2)
 iday   = ivalues(3)
 ihour  = ivalues(5)
 imin   = ivalues(6)
 isec   = ivalues(7)
 imsec  = ivalues(8)
 istarttime(1) = iyear
 istarttime(2) = imonth
 istarttime(3) = iday
 istarttime(4) = ihour
 istarttime(5) = imin
 istarttime(6) = isec
 !istarttime(7) = imsec
 starttime = iday*86400._4 + ihour*3600._4 + imin*60._4 + isec + imsec*0.001_4

 return
end subroutine initialise_timing

!--------------------------------------------------------------------
!+
!  Get time used since begining
!+
!--------------------------------------------------------------------
subroutine getused(tused)
 real(kind=4), intent(out) :: tused
 integer :: i,iday,imonth,ihour,imin,isec,imsec,ivalues(8)
 character(len=8)  :: date
 character(len=5)  :: zone
 character(len=10) :: time

 !tused = wallclockabs()

 !--do self-initialisation the first time it is called
 if (starttime < 0.) call initialise_timing

 call date_and_time(date,time,zone,ivalues)
 iday   = ivalues(3)
 ihour  = ivalues(5)
 imin   = ivalues(6)
 isec   = ivalues(7)
 imsec  = ivalues(8)

 if (ivalues(2) < istarttime(2)) then
    ivalues(2) = ivalues(2) + 12
 endif
 do i = istarttime(2), ivalues(2) - 1
    imonth = mod(i,12)
    if (imonth==4 .or. imonth==6 .or. imonth==9 .or. imonth==11) then
       iday = iday + 30
    elseif (imonth==2) then
       iday = iday + 28
    else
       iday = iday + 31
    endif
 enddo
 tused = iday*86400._4 + ihour*3600._4 + imin*60._4 + isec + imsec*0.001_4 - starttime

 return
end subroutine getused

!--------------------------------------------------------------------
!+
!  print (nicely formatted) time used since input time
!+
!--------------------------------------------------------------------
subroutine printused(t1,string,iunit)
 real(kind=4),     intent(in) :: t1
 character(len=*), intent(in), optional :: string
 integer,          intent(in), optional :: iunit
 character(len=64) :: newstring
 real(kind=4) :: t2
 integer :: lunit

 call getused(t2)

 if (present(string)) then
    newstring = trim(string(1:min(len(newstring),len_trim(string))))
 else
    newstring = 'completed in'
 endif

 if (present(iunit)) then
    lunit = iunit
 else
    lunit = 6
 endif

 call print_time(t2-t1,newstring,lunit)

 return
end subroutine printused

!--------------------------------------------------------------------
!+
!  print a time, nicelly formatted into hours, mins, seconds
!+
!--------------------------------------------------------------------
subroutine print_time(time,string,iunit)
 real(kind=4),     intent(in) :: time
 character(len=*), intent(in), optional :: string
 integer,          intent(in), optional :: iunit
 character(len=64) :: newstring
 integer :: nhr,nmin,lunit
 real(kind=4) :: trem

 trem = time
 nhr = int(trem/3600.)
 if (nhr > 0) trem = trem - nhr*3600._4

 nmin = int(trem/60.)
 if (nmin > 0) trem = trem - nmin*60._4

 if (present(string)) then
    newstring = trim(string(1:min(len(newstring),len_trim(string))))
 else
    newstring = 'completed in'
 endif

 if (present(iunit)) then
    lunit = iunit
 else
    lunit = 6
 endif

 if (nhr > 0) then
    write(lunit,"(1x,a,1x,i3,a,i2,a,f6.2,a,1pe12.4,a)") &
          trim(newstring),nhr,' hr, ',nmin,' min, ',trem,' s (=',time,'s)'
 elseif (nmin > 0) then
    write(lunit,"(1x,a,1x,i2,a,f6.2,a,1pe12.4,a)") &
          trim(newstring),nmin,' min, ',trem,' s (=',time,'s)'
 else
    write(lunit,"(1x,a,1x,f6.2,a)") trim(newstring),trem,' s'
 endif

 return
end subroutine print_time
!--------------------------------------------------------------------
!+
!  Interface to get both wall and cpu time
!+
!--------------------------------------------------------------------
subroutine get_timings(twall,tcpu)
 real(kind=4), intent(out) :: twall,tcpu

 call getused(twall)
 call cpu_time(tcpu)

 return
end subroutine get_timings

!--------------------------------------------------------------------
!+
real function wallclock()
!
! Return the wall clock time, with corrections for "across midnight" situations.
! From a routine originally by Aake Nordlund
!+
!-----------------------------------------------------------------------
 integer, save :: count, count_rate=0, count_max
 real, save :: previous=0., offset=0.

 if (count_rate == 0) then                                                     ! initialized?
    call system_clock(count=count, count_rate=count_rate, count_max=count_max)  ! no -- do it!
    offset = -count/real(count_rate)                                            ! time offset
 else
    call system_clock(count=count)                                              ! yes, get count
 endif
 wallclock = count/real(count_rate) + offset
 if (wallclock < previous) then                                                ! across midnight?
    offset = offset + real(count_max)/real(count_rate)                          ! adjust offset
    wallclock = count/real(count_rate) + offset                                 ! adjust time
 endif
 wallclock = wallclock-previous                                                ! differential
 previous = wallclock+previous                                                 ! remember

end function wallclock

!--------------------------------------------------------------------
!+
!real function wallclockabs()
!
! Return the wall clock time (no corrections for across midnight)
! From a routine originally by Aake Nordlund
! DJP: this version returns non-differential wallclock time
!+
!-----------------------------------------------------------------------
!  integer, save:: count, count_rate=0, count_max
!  real, save:: offset=0.
!
!  if (count_rate == 0) then                                                     ! initialized?
!    call system_clock(count=count, count_rate=count_rate, count_max=count_max)  ! no -- do it!
!    offset = -count/real(count_rate)                                            ! time offset
!  else
!   call system_clock(count=count)                                              ! yes, get count
!  endif
!  wallclockabs  = count/real(count_rate) + offset
!
!end function wallclockabs

end module timing
