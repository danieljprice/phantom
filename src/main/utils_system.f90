!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
! :Dependencies: None
!
 implicit none
 public :: get_command_option

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
 integer, intent(in), optional :: default
 character(len=80) :: string
 integer(kind=8)   :: val
 integer :: ierr,nargs,ieq,iarg

 val = 0.
 if (present(default)) val = default
 nargs = command_argument_count()
 do iarg=1,nargs
    call get_command_argument(iarg,string)
    ieq = index(string,'=')
    if (string(1:1)=='-' .and. index(string,variable) > 0 .and. ieq > 0) then
       read(string(ieq+1:),*,iostat=ierr) val
    endif
 enddo

end function get_command_option

end module systemutils
