!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: gitinfo
!
!  DESCRIPTION: writes the git information to the logfile
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io
!+
!--------------------------------------------------------------------------
module gitinfo
 implicit none
 public :: get_and_print_gitinfo
 character(len=7), public :: gitsha = ''

 private

contains
!--------------------------------------------------------------------------
subroutine get_and_print_gitinfo(iunit)
 use io,  only: igit
 integer, intent(in) :: iunit
 integer             :: i,j,k,io_local,lfile
 integer, parameter  :: nfiles = 34
 character(len= 16)  :: fmt
 character(len=128)  :: gitinfo,updatedfiles(nfiles)
 character(len=128)  :: updatedfiles3
 logical             :: iexist

 inquire(file="phantom_version",exist=iexist)
 if (iexist) then
    open(unit=igit,file="phantom_version")
    ! Copy in the nine opening lines of information
    do i = 1,9
       read(igit,'(a)') gitinfo
       write(iunit,'(1x,a)') trim(gitinfo)
       if (index(gitinfo,'ID:')> 0) gitsha = gitinfo(20:26)
    enddo
    i     = 0
    lfile = 0
    io_local = 0
    ! Read the list of modified files
    do while (io_local==0 .and. i < nfiles)
       read(igit,'(a)',iostat=io_local) gitinfo
       if (io_local==0) then
          i = i + 1
          if (len(trim(gitinfo))>38) then
             write(updatedfiles(i),"(2a)") gitinfo(1:35),"..."
          else
             updatedfiles(i) = gitinfo
          endif
          lfile = max(lfile,len(trim(updatedfiles(i))))
       endif
    enddo
    close(igit)
    if (i==0) then
       i = 1
       updatedfiles(i) = "none"
    endif
    if (i==nfiles) updatedfiles(i) = "Plus unlisted files"
    ! Put the files three in a line, and write to iunit
    k = 0
    do j = 1,i
       k = k + 1
       if (k==1) then
          write(updatedfiles3,"(7x,a)") trim(updatedfiles(j))
       else
          write(fmt,"(a,I2,a)") "(a,",3+lfile-len(trim(updatedfiles(j-1))),"x,a)"
          write(updatedfiles3,fmt) trim(updatedfiles3),trim(updatedfiles(j))
       endif
       if (k==3 .or. j==i) then
          write(iunit,"(a)") trim(updatedfiles3)
          k = 0
       endif
    enddo
 else
    write(iunit,'(a)') "phantom_version not found; no git information to display"
 endif
 write(iunit,'(a)') ' '

end subroutine get_and_print_gitinfo
!--------------------------------------------------------------------------
end module gitinfo
