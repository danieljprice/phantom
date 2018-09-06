!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: strings_utils
!
!  DESCRIPTION:
!   This module contains low-level utilities related to manipulating strings
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module strings_utils
 implicit none

contains

!--------------------------------------------------------------------------
!+
!  create an array of numbered strings
!+
!--------------------------------------------------------------------------
subroutine array_of_numbered_strings(pre_string,post_string,complete_string)
 character(len=*), intent(in)  :: pre_string,post_string
 character(len=*), intent(out) :: complete_string(:)

 integer :: i,N,total_len,int_len
 character(len=20) :: num_string,fmt1,fmt2

 N = size(complete_string)
 int_len = floor(log10(real(N) + tiny(0.))) + 1
 total_len = len(trim(pre_string)) + int_len  + len(trim(post_string))
 if (len(complete_string) < total_len) then
    print*,'output string is not long enough!'
    stop
 endif

 if (N > 1) then
    do i=1,N
       write(fmt2,'(I0)') int_len
       fmt1 = '(I0.'//trim(fmt2)//')'
       write(num_string,fmt1) i
       write(complete_string(i),'(A)') trim(pre_string)//trim(adjustl(num_string))//trim(post_string)
    enddo
 else
    write(complete_string,'(A)') trim(pre_string)//trim(post_string)
 endif

 return
end subroutine array_of_numbered_strings

end module strings_utils
