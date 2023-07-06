!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module derivutils
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: io, mpiutils, timing
!
 use timing, only: timers,itimer_dens,itimer_force,itimer_link,itimer_extf,itimer_balance,itimer_cons2prim

 implicit none

 private

 public :: do_timing

contains

!-------------------------------------------------------------
!+
!  interface to timing routines
!+
!-------------------------------------------------------------
subroutine do_timing(label,tlast,tcpulast,start,lunit)
 use io,       only:iverbose,id,master
 use mpiutils, only:reduce_mpi
 use timing,   only:increment_timer,log_timing,get_timings
 character(len=*), intent(in)    :: label
 real(kind=4),     intent(inout) :: tlast,tcpulast
 logical,          intent(in), optional :: start
 integer,          intent(in), optional :: lunit
 real(kind=4) :: t2,tcpu2,tcpu,twall

 call get_timings(t2,tcpu2)
 twall = reduce_mpi('+',t2-tlast)
 tcpu  = reduce_mpi('+',tcpu2-tcpulast)

 if (label=='dens') then
    call increment_timer(itimer_dens,t2-tlast,tcpu2-tcpulast)
 elseif (label=='force') then
    call increment_timer(itimer_force,t2-tlast,tcpu2-tcpulast)
 elseif (label=='link') then
    call increment_timer(itimer_link,t2-tlast,tcpu2-tcpulast)
 elseif (label=='cons2prim') then
    call increment_timer(itimer_cons2prim,t2-tlast,tcpu2-tcpulast)
 endif

 if (iverbose >= 2 .and. id==master) then
    if (present(start)) then
       call log_timing(label,t2-tlast,tcpu,start=.true.)
    elseif (present(lunit)) then
       call log_timing(label,t2-tlast,tcpu,iunit=lunit)
    else
       call log_timing(label,t2-tlast,tcpu)
    endif
 endif
 tlast = t2
 tcpulast = tcpu2

end subroutine do_timing

end module derivutils
