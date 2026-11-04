!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module systemutils
!
! Utilities for retrieving command line arguments and flags
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils
!
 implicit none
 public :: get_command_option, get_command_option_real, get_command_option_logical

 private

contains
!-------------------------------------------------------------------
!+
!  find integer-valued option from command line arguments
!  as in --arg=blah
!+
!-------------------------------------------------------------------
function get_command_option(variable,default) result(val)
 character(len=*), intent(in) :: variable
 integer,          intent(in), optional :: default
 character(len=80) :: string
 integer(kind=8)   :: val,ival
 real(kind=8)      :: rval
 integer :: ierr

 val = 0
 if (present(default)) val = default
 if (find_option(variable,string)) then
    read(string,*,iostat=ierr) ival
    if (ierr == 0) then
       val = ival
    else
       ! also allow values such as 1e6
       read(string,*,iostat=ierr) rval
       if (ierr == 0) val = nint(rval,kind=8)
    endif
 endif

end function get_command_option

!-------------------------------------------------------------------
!+
!  find real-valued option from command line arguments
!  as in --arg=blah
!+
!-------------------------------------------------------------------
function get_command_option_real(variable,default) result(val)
 character(len=*), intent(in) :: variable
 real,             intent(in), optional :: default
 character(len=80) :: string
 real(kind=8) :: val,rval
 integer :: ierr

 val = 0.d0
 if (present(default)) val = default
 if (find_option(variable,string)) then
    read(string,*,iostat=ierr) rval
    if (ierr == 0) val = rval
 endif

end function get_command_option_real

!-------------------------------------------------------------------
!+
!  find logical-valued option from command line arguments
!  as in --arg=blah
!+
!-------------------------------------------------------------------
function get_command_option_logical(variable,default) result(val)
 character(len=*), intent(in) :: variable
 logical,          intent(in), optional :: default
 character(len=80) :: string
 logical :: val,lval
 integer :: ierr

 val = .false.
 if (present(default)) val = default
 if (find_option(variable,string)) then
    read(string,*,iostat=ierr) lval
    if (ierr == 0) val = lval
 endif

end function get_command_option_logical

!-------------------------------------------------------------------
!+
!  find the value string of the option with exactly this name
!  on the command line (the last one, if given more than once)
!+
!-------------------------------------------------------------------
logical function find_option(variable,valstring)
 use infile_utils, only:get_option_arg
 character(len=*), intent(in)  :: variable
 character(len=*), intent(out) :: valstring
 character(len=80) :: name,val
 integer :: iarg

 find_option = .false.
 valstring = ''
 do iarg=1,command_argument_count()
    if (get_option_arg(iarg,name,val)) then
       if (trim(name) == trim(variable)) then
          find_option = .true.
          valstring = val
       endif
    endif
 enddo

end function find_option

end module systemutils
