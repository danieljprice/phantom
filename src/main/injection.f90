module injection
!
! Wrapper module as interface to swap-in injection modules,
! contains options common to all injection modules
!
 use dim, only:inject_parts
 implicit none

 real, public :: rkill = -1.

 public :: write_options_injection,read_options_injection

contains

!-----------------------------------------------------------------------
!+
!  Wrapper routine for writing injection options
!+
!-----------------------------------------------------------------------
subroutine write_options_injection(iunit)
 use infile_utils, only:write_inopt
 use inject,       only:write_options_inject,inject_type,update_injected_par
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for injecting/removing particles'
 if (inject_parts) call write_options_inject(iunit)
 if (inject_parts .and. inject_type=='sim') call update_injected_par()
 call write_inopt(rkill,'rkill','deactivate particles outside this radius (<0 is off)',iunit)

end subroutine write_options_injection

!-----------------------------------------------------------------------
!+
!  Wrapper routine for reading injection options
!+
!-----------------------------------------------------------------------
subroutine read_options_injection(name,valstring,imatch,igotall,ierr)
 use inject, only:read_options_inject
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer, intent(out) :: ierr
 integer, save :: ngot = 0
 logical :: igotallinject

 imatch  = .true.
 igotall = .false.
 igotallinject = .false.
 select case(trim(name))
 case('rkill')
    read(valstring,*,iostat=ierr) rkill
    ngot = ngot + 1
 case default
    imatch = .false.
    if (.not.imatch .and. inject_parts) call read_options_inject(name,valstring,imatch,igotallinject,ierr)
 end select

 igotall = (ngot >= 0) .and. igotallinject  ! rkill is not compulsory

end subroutine read_options_injection

end module injection