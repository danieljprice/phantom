!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: datautils
!
!  DESCRIPTION:
!   This module contains utilities to transparently handle
!   finding/reading external data files
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
module datautils
 implicit none
 public :: find_datafile

 private

contains

!----------------------------------------------------------------
!+
!  Search for a datafile, download it if necessary
!  Returns the full path to the data file
!+
!----------------------------------------------------------------
function find_datafile(filename,dir,env_var,url) result(filepath)
 character(len=*), intent(in) :: filename
 character(len=*), intent(in), optional :: dir,env_var,url
 character(len=120) :: filepath
 character(len=120) :: mydir,env_dir,my_url
 character(len=40)  :: my_env_var
 logical :: iexist
 integer :: ierr
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
 else
    ierr = 0
    mydir  = ' '
    if (present(env_var)) then
       my_env_var = env_var
    else
       my_env_var = 'DATA_DIR'
    endif
    call get_environment_variable(my_env_var,env_dir)
    if (len_trim(env_dir) > 0) then
       mydir = trim(env_dir)
       if (present(dir)) mydir = trim(mydir)//'/'//trim(dir)//'/'
       print "(a)",' Reading '//trim(filename)//' in '//trim(mydir)//' (from '//trim(my_env_var)//' setting)'
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
                   print "(a)",' ERROR: downloaded file '//trim(filename)//' does not exist in '//trim(filepath)
                   ierr = 1
                   filepath = trim(filename)
                else
                   print "(a)",' DOWNLOADED '//trim(filename)//' TO '//trim(filepath)
                endif
             else
                filepath = trim(filename)
                print "(a)",' ERROR downloading file'
             endif
          endif
       endif
    else
       if (present(dir)) then
          if (len_trim(dir) > 0) mydir = trim(dir)//'/'
       endif
       print "(3(/,1x,a),/)",'* Looking for '//trim(filename)//' in current directory  *', &
                             '* Set '//trim(my_env_var)//' environment variable to read files *', &
                             '* from the $('//trim(my_env_var)//')/'//trim(mydir)//' directory  *'
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

 print "(80('-'))"
 print "(a)",'  Downloading '//trim(file)//' from '//trim(url)
 print "(80('-'))"

 ierr = 0
 ! check that wget utility exists
 !call execute_command_line('type -p wget > /dev/null',wait=.true.,exitstat=ierr,cmdstat=ierr1)
 call system('type -p wget > /dev/null')

 if (ierr /= 0) then
    print "(a)",' ERROR: wget utility does not exist'
 else
    if (len_trim(dir) > 0) then
       !call execute_command_line('wget '//trim(url)//trim(file)//' -O '//trim(dir)//trim(file),wait=.true.,&
       !                          exitstat=ierr,cmdstat=ierr1)
       call system('wget '//trim(url)//trim(file)//' -O '//trim(dir)//trim(file))
       localfile = trim(dir)//trim(file)
    else
       !call execute_command_line('wget '//trim(url)//trim(file),wait=.true.,exitstat=ierr,cmdstat=ierr1)
       call system('wget '//trim(url)//trim(file))
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
 if (ierr /= 0) then
    has_write_permission = .false.
 endif
 close(iunit,status='delete')

end function has_write_permission

end module datautils
