!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: timing
!
!  DESCRIPTION:
!   This module contains utilities for code timings
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module timing
 implicit none
 integer, private :: istarttime(6)
 real(kind=4), private :: starttime

 data starttime/-1./

 public :: getused,get_timings,printused
 public :: wallclock,log_timing,print_time

 public :: timer,reset_timer,increment_timer

 type timer
    character(len=20) :: label
    real(kind=4) :: wall
    real(kind=4) :: cpu
 end type

 type(timer), public     :: timer_dens,timer_force,timer_link
 private

contains
!--------------------------------------
!+
!  Routines to handle timer objects
!+
!--------------------------------------
subroutine reset_timer(my_timer,label,wall,cpu)
 type(timer),      intent(out) :: my_timer
 real(kind=4),     intent(in), optional :: wall, cpu
 character(len=*), intent(in), optional :: label

 if (present(wall)) then
    my_timer%wall = wall
 else
    my_timer%wall = 0.
 endif
 if (present(cpu)) then
    my_timer%cpu = cpu
 else
    my_timer%cpu = 0.
 endif
 if (present(label)) then
    my_timer%label = trim(label)
 endif

end subroutine reset_timer

subroutine increment_timer(my_timer,wall,cpu)
 type(timer),  intent(inout) :: my_timer
 real(kind=4), intent(in) :: wall, cpu

 my_timer%wall = my_timer%wall + wall
 my_timer%cpu  = my_timer%cpu + cpu

end subroutine increment_timer

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
