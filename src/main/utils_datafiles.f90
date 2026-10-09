!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Dependencies: iso_c_binding
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
function find_datafile(filename,dir,env_var,url,url_fallback,verbose) result(filepath)
 character(len=*), intent(in) :: filename
 character(len=*), intent(in), optional :: dir,env_var,url,url_fallback
 logical,          intent(in), optional :: verbose
 character(len=:), allocatable :: filepath,mydir,env_dir,my_url,my_env_var
 logical :: iexist,isverbose,ok
 integer :: ierr,env_len

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
    my_env_var = ' '
    if (present(env_var)) then
       my_env_var = env_var
       call get_environment_variable(my_env_var,length=env_len)
       allocate(character(len=env_len) :: env_dir)
       call get_environment_variable(my_env_var,env_dir,status=ierr)
       if (ierr /= 0) env_dir = ''
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
       if (iexist) then
          ! verify sidecar checksum if present; delete corrupt copies so we can re-download
          call verify_sidecar_md5(filepath,ok)
          if (.not.ok) then
             call delete_if_exists(filepath)
             iexist = .false.
             if (isverbose) print "(a)",' ERROR: checksum mismatch for '//trim(filepath)// &
                                        '; deleted, will try to re-download'
          endif
       endif
       if (.not.iexist) then
          if (present(url)) then
             !
             ! try to download the file from a remote url, then optional fallback
             !
             my_url = url
             call download_datafile(trim(my_url),trim(mydir),trim(filename),filepath,ierr)
             if (ierr /= 0 .and. present(url_fallback)) then
                if (len_trim(url_fallback) > 0) then
                   if (isverbose) print "(a)",' trying fallback URL...'
                   call download_datafile(trim(url_fallback),trim(mydir),trim(filename),filepath,ierr)
                endif
             endif
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
                if (isverbose) print "(a)",' ERROR downloading file: file not found on server'
             endif
          endif
       endif
    else
       if (present(dir)) then
          if (len_trim(dir) > 0) mydir = trim(dir)//'/'
       endif
       if (isverbose) then
          if (present(env_var)) then
             print "(3(/,1x,a),/)",'* DATAFILE NOT FOUND: '//trim(filename)//' (not in current directory)  *', &
                                   '* PLEASE TYPE "export '//trim(my_env_var)//'=/where/the/code/is" IN YOUR TERMINAL'
          else
             print "(a)",' ERROR: datafile not found: '//trim(filename)
          endif
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
 character(len=:), allocatable, intent(out) :: filepath
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
!  use curl to retrieve file and check that it succeeds
!+
!---------------------------------------------------------
subroutine retrieve_remote_file(url,file,dir,localfile,ierr)
 character(len=*), intent(in)  :: url,file,dir
 character(len=:), allocatable, intent(out) :: localfile
 integer,          intent(out) :: ierr
 integer :: ilen,ierr1,cmdstat
 logical :: iexist,ishtml
 character(len=:), allocatable :: cmdline
 character(len=64)  :: expected_md5,actual_md5
 character(len=*), parameter :: curlcmd = 'curl -fL'

 print "(80('-'))"
 print "(a)",'  Downloading '//trim(file)//' from '//trim(url)
 print "(80('-'))"

 ierr = 0
 expected_md5 = ' '
 cmdstat = 0

 if (len_trim(dir) > 0) then
    localfile = trim(dir)//trim(file)
 else
    localfile = trim(file)
 endif

 ! for Zenodo URLs, fetch expected MD5 from the record API before downloading
 if (index(url,'zenodo.org') > 0) then
    call get_zenodo_md5(url,file,expected_md5,ierr1)
    if (ierr1 /= 0) then
       print "(a)",' WARNING: could not retrieve Zenodo checksum; will skip MD5 verification'
       expected_md5 = ' '
    endif
 endif

 cmdline = trim(curlcmd)//' '//trim(url)//trim(file)//' -o '//trim(localfile)
 call execute_command_line(trim(cmdline),wait=.true.,exitstat=ierr,cmdstat=cmdstat)
 if (cmdstat /= 0) ierr = 1
 print "(80('-'))"

 ! check that the file has actually downloaded correctly
 inquire(file=trim(localfile),exist=iexist,size=ilen)

 if (ierr /= 0) then
    print "(a)",' ERROR: file not found on server (or curl failed / is not available)'
    call delete_if_exists(localfile)
    return
 endif

 if (.not.iexist) then
    print "(a)",' ERROR: downloaded file does not exist'
    ierr = 1
    return
 endif

 if (ilen == 0) then
    print "(a)",' ERROR: downloaded file is empty'
    call delete_if_exists(localfile)
    ierr = 2
    return
 endif

 ! reject HTML error pages that some servers return with HTTP 200
 call file_looks_like_html(localfile,ishtml)
 if (ishtml) then
    print "(a)",' ERROR: file not found on server (downloaded content looks like HTML)'
    call delete_if_exists(localfile)
    ierr = 3
    return
 endif

 ! verify MD5 against Zenodo metadata when available
 if (len_trim(expected_md5) > 0) then
    call compute_md5(localfile,actual_md5,ierr1)
    if (ierr1 /= 0) then
       print "(a)",' WARNING: could not compute MD5 of downloaded file'
    elseif (.not.md5_equal(expected_md5,actual_md5)) then
       print "(a)",' ERROR: checksum mismatch for '//trim(file)
       print "(a)",'   expected md5: '//trim(expected_md5)
       print "(a)",'   got md5:      '//trim(actual_md5)
       call delete_if_exists(localfile)
       ierr = 4
       return
    else
       call write_md5_sidecar(localfile,expected_md5)
    endif
 endif

end subroutine retrieve_remote_file

!---------------------------------------------------------
!+
!  delete a file if it exists
!+
!---------------------------------------------------------
subroutine delete_if_exists(path)
 character(len=*), intent(in) :: path
 integer :: iunit,ierr
 logical :: iexist

 inquire(file=trim(path),exist=iexist)
 if (.not.iexist) return
 open(newunit=iunit,file=trim(path),status='old',iostat=ierr)
 if (ierr == 0) close(iunit,status='delete')

end subroutine delete_if_exists

!---------------------------------------------------------
!+
!  true if the start of the file looks like an HTML page
!+
!---------------------------------------------------------
subroutine file_looks_like_html(path,ishtml)
 character(len=*), intent(in)  :: path
 logical,          intent(out) :: ishtml
 integer :: iunit,ierr
 character(len=256) :: line

 ishtml = .false.
 open(newunit=iunit,file=trim(path),status='old',action='read',iostat=ierr)
 if (ierr /= 0) return
 do
    read(iunit,'(a)',iostat=ierr) line
    if (ierr /= 0) exit
    if (len_trim(line) == 0) cycle
    line = adjustl(line)
    if (index(line,'<!DOCTYPE') == 1 .or. index(line,'<!doctype') == 1 .or. &
        index(line,'<html') == 1 .or. index(line,'<HTML') == 1 .or. &
        index(line,'<head') == 1 .or. index(line,'<HEAD') == 1) then
       ishtml = .true.
    endif
    exit
 enddo
 close(iunit)

end subroutine file_looks_like_html

!---------------------------------------------------------
!+
!  extract Zenodo record id from a download URL and fetch
!  the expected MD5 for filename from the record API JSON
!+
!---------------------------------------------------------
subroutine get_zenodo_md5(url,filename,md5hex,ierr)
 character(len=*), intent(in)  :: url,filename
 character(len=*), intent(out) :: md5hex
 integer,          intent(out) :: ierr
 character(len=:), allocatable :: recid,apiurl,tmpfile,cmdline
 integer :: i1,i2,ierr1,cmdstat

 md5hex = ' '
 ierr = 1
 i1 = index(url,'/records/')
 if (i1 <= 0) return
 i1 = i1 + len('/records/')
 i2 = index(url(i1:),'/')
 if (i2 <= 0) then
    recid = trim(url(i1:))
 else
    recid = trim(url(i1:i1+i2-2))
 endif
 if (len_trim(recid) == 0) return

 apiurl = 'https://zenodo.org/api/records/'//trim(recid)
 tmpfile = 'zenodo_api_tmp_'//process_suffix()//'.json'
 cmdline = 'curl -fL '//trim(apiurl)//' -o '//trim(tmpfile)
 call execute_command_line(trim(cmdline),wait=.true.,exitstat=ierr1,cmdstat=cmdstat)
 if (cmdstat /= 0 .or. ierr1 /= 0) then
    call delete_if_exists(tmpfile)
    return
 endif

 call extract_md5_from_zenodo_json(tmpfile,filename,md5hex,ierr)
 call delete_if_exists(tmpfile)

end subroutine get_zenodo_md5

!---------------------------------------------------------
!+
!  scan Zenodo API JSON for "key":"<filename>" and nearby
!  "checksum":"md5:<hex>" (no JSON library required)
!+
!---------------------------------------------------------
subroutine extract_md5_from_zenodo_json(jsonfile,filename,md5hex,ierr)
 character(len=*), intent(in)  :: jsonfile,filename
 character(len=*), intent(out) :: md5hex
 integer,          intent(out) :: ierr
 integer :: iunit,ios,n,i,j,k,keypos,cpos
 character(len=:), allocatable :: buf,keystr
 logical :: iexist

 md5hex = ' '
 ierr = 1
 inquire(file=trim(jsonfile),exist=iexist,size=n)
 if (.not.iexist .or. n <= 0) return

 ! read the whole file as a stream; do NOT assign buf=' ' afterwards —
 ! that would reallocate a deferred-length character to len=1 and corrupt memory
 allocate(character(len=n) :: buf)
 open(newunit=iunit,file=trim(jsonfile),status='old',action='read', &
      access='stream',form='unformatted',iostat=ios)
 if (ios /= 0) then
    deallocate(buf)
    return
 endif
 read(iunit,iostat=ios) buf
 close(iunit)
 if (ios /= 0) then
    deallocate(buf)
    return
 endif

 keystr = '"key":"'//trim(filename)//'"'
 keypos = index(buf,trim(keystr))
 if (keypos <= 0) then
    ! also try with spaces after colon
    keystr = '"key": "'//trim(filename)//'"'
    keypos = index(buf,trim(keystr))
 endif
 if (keypos <= 0) then
    deallocate(buf)
    return
 endif

 ! search near the key for md5: (checksum usually follows key in the same JSON object)
 i = max(1, keypos - 200)
 j = min(len(buf), keypos + len_trim(keystr) + 400)
 cpos = index(buf(i:j),'md5:')
 if (cpos <= 0) then
    deallocate(buf)
    return
 endif
 j = i + cpos + 3  ! position after 'md5:'
 k = 0
 md5hex = ' '
 do while (j <= len(buf) .and. k < 32)
    if ((buf(j:j) >= '0' .and. buf(j:j) <= '9') .or. &
        (buf(j:j) >= 'a' .and. buf(j:j) <= 'f') .or. &
        (buf(j:j) >= 'A' .and. buf(j:j) <= 'F')) then
       k = k + 1
       md5hex(k:k) = buf(j:j)
       j = j + 1
    else
       exit
    endif
 enddo
 deallocate(buf)
 if (k == 32) ierr = 0

end subroutine extract_md5_from_zenodo_json

!---------------------------------------------------------
!+
!  compute MD5 of a file using openssl (portable)
!+
!---------------------------------------------------------
subroutine compute_md5(path,md5hex,ierr)
 character(len=*), intent(in)  :: path
 character(len=*), intent(out) :: md5hex
 integer,          intent(out) :: ierr
 character(len=:), allocatable :: cmdline,tmpfile
 character(len=256) :: line
 integer :: iunit,ios,cmdstat,ierr1,i,j

 md5hex = ' '
 ierr = 1
 tmpfile = 'datafile_md5_tmp_'//process_suffix()//'.txt'
 ! Read via stdin so the output line does not contain a potentially long path.
 cmdline = 'openssl dgst -md5 < '//trim(path)//' > '//trim(tmpfile)//' 2>/dev/null'
 call execute_command_line(trim(cmdline),wait=.true.,exitstat=ierr1,cmdstat=cmdstat)
 if (cmdstat /= 0 .or. ierr1 /= 0) then
    call delete_if_exists(tmpfile)
    return
 endif
 open(newunit=iunit,file=trim(tmpfile),status='old',action='read',iostat=ios)
 if (ios /= 0) then
    call delete_if_exists(tmpfile)
    return
 endif
 read(iunit,'(a)',iostat=ios) line
 close(iunit)
 call delete_if_exists(tmpfile)
 if (ios /= 0) return

 ! openssl output is typically "MD5(filename)= hex" or "MD5(filename)= hex"
 i = index(line,'=')
 if (i <= 0) return
 line = adjustl(line(i+1:))
 j = 0
 md5hex = ' '
 do i = 1,len_trim(line)
    if ((line(i:i) >= '0' .and. line(i:i) <= '9') .or. &
        (line(i:i) >= 'a' .and. line(i:i) <= 'f') .or. &
        (line(i:i) >= 'A' .and. line(i:i) <= 'F')) then
       j = j + 1
       if (j <= len(md5hex)) md5hex(j:j) = line(i:i)
    endif
 enddo
 if (j == 32) ierr = 0

end subroutine compute_md5

!---------------------------------------------------------
!+
!  case-insensitive compare of two MD5 hex strings
!+
!---------------------------------------------------------
logical function md5_equal(a,b)
 character(len=*), intent(in) :: a,b
 character(len=32) :: aa,bb
 integer :: i
 character :: ca,cb

 md5_equal = .false.
 if (len_trim(a) < 32 .or. len_trim(b) < 32) return
 aa = a(1:32)
 bb = b(1:32)
 do i = 1,32
    ca = aa(i:i)
    cb = bb(i:i)
    if (ca >= 'A' .and. ca <= 'Z') ca = achar(iachar(ca) + 32)
    if (cb >= 'A' .and. cb <= 'Z') cb = achar(iachar(cb) + 32)
    if (ca /= cb) return
 enddo
 md5_equal = .true.

end function md5_equal

!---------------------------------------------------------
!+
!  write a one-line .md5 sidecar next to the data file
!+
!---------------------------------------------------------
subroutine write_md5_sidecar(path,md5hex)
 character(len=*), intent(in) :: path,md5hex
 integer :: iunit,ierr

 open(newunit=iunit,file=trim(path)//'.md5',status='replace',action='write',iostat=ierr)
 if (ierr /= 0) return
 write(iunit,'(a)',iostat=ierr) trim(md5hex)
 close(iunit)

end subroutine write_md5_sidecar

!---------------------------------------------------------
!+
!  if path.md5 exists, verify the file matches; ok=.true.
!  if no sidecar or checksum matches
!+
!---------------------------------------------------------
subroutine verify_sidecar_md5(path,ok)
 character(len=*), intent(in)  :: path
 logical,          intent(out) :: ok
 character(len=64) :: expected,actual
 integer :: iunit,ierr
 logical :: iexist

 ok = .true.
 inquire(file=trim(path)//'.md5',exist=iexist)
 if (.not.iexist) return

 open(newunit=iunit,file=trim(path)//'.md5',status='old',action='read',iostat=ierr)
 if (ierr /= 0) return
 read(iunit,'(a)',iostat=ierr) expected
 close(iunit)
 if (ierr /= 0) return

 call compute_md5(path,actual,ierr)
 if (ierr /= 0) return
 ok = md5_equal(expected,actual)

end subroutine verify_sidecar_md5

!---------------------------------------------------------
!+
!  suffix to keep temporary files separate between processes
!+
!---------------------------------------------------------
function process_suffix() result(suffix)
 use iso_c_binding, only:c_int
 character(len=:), allocatable :: suffix
 character(len=32) :: pid_string
 interface
    function process_id() bind(C,name='getpid') result(pid)
     import c_int
     integer(c_int) :: pid
    end function process_id
 end interface

 write(pid_string,'(i0)') process_id()
 suffix = trim(pid_string)

end function process_suffix

!---------------------------------------------------------
!+
!  function to check if a directory has write permissions
!+
!---------------------------------------------------------
logical function has_write_permission(dir)
 character(len=*), intent(in) :: dir
 integer :: iunit,ierr

 has_write_permission = .true.
 open(newunit=iunit,file=trim(dir)//'data.tmp.'//process_suffix(),action='write',iostat=ierr)
 if (ierr /= 0) then
    has_write_permission = .false.
    return
 endif

 close(iunit,status='delete',iostat=ierr)
 if (ierr /= 0) has_write_permission = .false.

end function has_write_permission

end module datautils
