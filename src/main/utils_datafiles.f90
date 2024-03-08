!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module datautils
!
! This module contains utilities to transparently handle
!   finding/reading external data files
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
 public :: find_datafile,download_datafile

 private

contains

!----------------------------------------------------------------
!+
!  Search for a datafile, download it if necessary
!  Returns the full path to the data file
!+
!----------------------------------------------------------------
function find_datafile(filename,dir,env_var,url,verbose) result(filepath)
 character(len=*), intent(in) :: filename
 character(len=*), intent(in), optional :: dir,env_var,url
 logical,          intent(in), optional :: verbose
 character(len=120) :: filepath
 character(len=120) :: mydir,env_dir,my_url
 character(len=40)  :: my_env_var
 logical :: iexist,isverbose
 integer :: ierr

 isverbose = .true.
 if (present(verbose)) isverbose = verbose
!
!  We search for the file by:
!
!  1) checking current directory
!  2) looking in data directory specified by environment variable
!  3) downloading from remote url into data directory (if writeable) or current dir (if not)
!
 inquire(file=trim(filename),exist=iexist)
 if (iexist) then
    filepath = trim(filename)
    if (isverbose) print "(a)",' reading '//trim(filepath)
 else
    ierr = 0
    mydir  = ' '
    if (present(env_var)) then
       my_env_var = env_var
       call get_environment_variable(my_env_var,env_dir)
    elseif (present(dir)) then
       env_dir = dir
    else
       env_dir = './'
    endif
    if (len_trim(env_dir) > 0) then
       mydir = trim(env_dir)
       if (present(dir)) mydir = trim(mydir)//'/'//trim(dir)//'/'
       if (isverbose .and. present(env_var)) then
          print "(a)",' Reading '//trim(filename)//' in '//trim(mydir)//&
                                  ' (from '//trim(my_env_var)//' setting)'
       endif
       filepath = trim(mydir)//trim(filename)
       inquire(file=trim(filepath),exist=iexist)
       if (.not.iexist) then
          if (present(url)) then
             !
             ! try to download the file from a remote url
             !
             my_url = url
             if (present(dir)) then
                if (len_trim(dir) > 0) my_url = trim(url)//'/'//trim(dir)//'/'
             else
                my_url = url
             endif
             call download_datafile(trim(my_url),trim(mydir),trim(filename),filepath,ierr)
             if (ierr == 0) then
                inquire(file=trim(filepath),exist=iexist)
                if (.not.iexist) then
                   if (isverbose) print "(a)",' ERROR: downloaded file '//trim(filename)// &
                                              ' does not exist in '//trim(filepath)
                   ierr = 1
                   filepath = trim(filename)
                else
                   if (isverbose) print "(a)",' DOWNLOADED '//trim(filename)//' TO '//trim(filepath)
                endif
             else
                filepath = trim(filename)
                if (isverbose) print "(a)",' ERROR downloading file'
             endif
          endif
       endif
    else
       if (present(dir)) then
          if (len_trim(dir) > 0) mydir = trim(dir)//'/'
       endif
       if (isverbose) then
          print "(3(/,1x,a),/)",'* DATAFILE NOT FOUND: '//trim(filename)//' (not in current directory)  *', &
                                '* PLEASE TYPE "export '//trim(my_env_var)//'=/where/the/code/is" IN YOUR TERMINAL'
       endif
       filepath = trim(filename)
    endif
 endif

end function find_datafile

!---------------------------------------------------------
!+
!  download file from remote url to specified directory
!  (defaults to current dir if specified dir does not
!   have write permissions)
!+
!---------------------------------------------------------
subroutine download_datafile(url,dir,filename,filepath,ierr)
 character(len=*), intent(in)  :: url, dir
 character(len=*), intent(in)  :: filename
 character(len=*), intent(out) :: filepath
 integer,          intent(out) :: ierr

 if (has_write_permission(dir)) then     ! download to data/ directory
    call retrieve_remote_file(url,trim(adjustl(filename)),trim(dir),filepath,ierr)
 elseif (has_write_permission('')) then  ! download to current directory
    print*,'ERROR: cannot write to '//trim(dir)//', writing to current directory'
    call retrieve_remote_file(url,trim(adjustl(filename)),'',filepath,ierr)
 else  ! return an error
    filepath = trim(filename)
    print*,'ERROR: cannot write to '//trim(dir)//' or current directory'
    ierr = 1
 endif

end subroutine download_datafile

!---------------------------------------------------------
!+
!  use wget to retrieve file and check that it succeeds
!+
!---------------------------------------------------------
subroutine retrieve_remote_file(url,file,dir,localfile,ierr)
 character(len=*), intent(in)  :: url,file,dir
 character(len=*), intent(out) :: localfile
 integer,          intent(out) :: ierr
 integer :: ilen,iunit!,ierr1
 logical :: iexist
 character(len=*), parameter :: cmd = 'curl -k'

 print "(80('-'))"
 print "(a)",'  Downloading '//trim(file)//' from '//trim(url)
 print "(80('-'))"

 ierr = 0
 ! check that wget utility exists
 !call execute_command_line('type -p wget > /dev/null',wait=.true.,exitstat=ierr,cmdstat=ierr1)
 call system('type -p curl > /dev/null')

 if (ierr /= 0) then
    print "(a)",' ERROR: curl utility does not exist'
 else
    if (len_trim(dir) > 0) then
       !call execute_command_line(trim(cmd)//' '//trim(url)//trim(file)//' -O '//trim(dir)//trim(file),wait=.true.,&
       !                          exitstat=ierr,cmdstat=ierr1)
       call system(trim(cmd)//' '//trim(url)//trim(file)//' -o '//trim(dir)//trim(file))
       localfile = trim(dir)//trim(file)
    else
       !call execute_command_line(trim(cmd)//' '//trim(url)//trim(file),wait=.true.,exitstat=ierr,cmdstat=ierr1)
       call system(trim(cmd)//' '//trim(url)//trim(file)//' -o '//trim(file))
       localfile = trim(file)
    endif
 endif
 print "(80('-'))"
 !print*,'EXITSTAT=',ierr

 ! check that the file has actually downloaded correctly
 inquire(file=trim(localfile),exist=iexist,size=ilen)

 if (ierr==0 .and. .not.iexist) then
    print*,'ERROR: downloaded file does not exist'
    ierr = 1
    return
 endif
 if (iexist .and. ilen==0) then
    print*,'ERROR: downloaded file is empty'
    ! delete file if it is empty
    open(newunit=iunit,file=trim(localfile),status='old',iostat=ierr)
    close(iunit,status='delete')
    ierr = 2
 endif

end subroutine retrieve_remote_file

!---------------------------------------------------------
!+
!  function to check if a directory has write permissions
!+
!---------------------------------------------------------
logical function has_write_permission(dir)
 character(len=*), intent(in) :: dir
 integer :: iunit,ierr

 has_write_permission = .true.
 open(newunit=iunit,file=trim(dir)//'data.tmp.abcd',action='write',iostat=ierr)
 if (ierr /= 0) has_write_permission = .false.

 close(iunit,status='delete',iostat=ierr)
 if (ierr /= 0) has_write_permission = .false.

end function has_write_permission

end module datautils
