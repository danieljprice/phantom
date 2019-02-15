module derivutils
 use timing, only: timer

 implicit none


 type(timer), public     :: timer_dens,timer_force,timer_link,timer_extf

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
    call increment_timer(timer_dens,t2-tlast,tcpu2-tcpulast)
 else if (label=='force') then
    call increment_timer(timer_force,t2-tlast,tcpu2-tcpulast)
 else if (label=='link') then
    call increment_timer(timer_link,t2-tlast,tcpu2-tcpulast)
 else if (label=='extf') then
    call increment_timer(timer_extf,t2-tlast,tcpu2-tcpulast)
 endif

 if (iverbose >= 2) then
    if (id==master) then
       if (present(start)) then
          call log_timing(label,t2-tlast,tcpu,start=.true.)
       elseif (present(lunit)) then
          call log_timing(label,t2-tlast,tcpu,iunit=lunit)
       else
          call log_timing(label,t2-tlast,tcpu)
       endif
    endif
 endif
 tlast = t2
 tcpulast = tcpu2

end subroutine do_timing
end module derivutils
