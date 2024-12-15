!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module ephemeris
!
! Download and read data files from the JPL horizons ephemeris service
!
! :References:
! https://ssd-api.jpl.nasa.gov/doc/horizons.html
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: datautils, fileutils
!
 implicit none
 integer, parameter :: nelem = 9

 public :: get_ephemeris,nelem

 private

contains

!-----------------------------------------------------------------------
!+
!  read data files from the JPL Horizons ephemeris server
!+
!-----------------------------------------------------------------------
function get_ephemeris(object,ierr) result(elems)
 use datautils, only:download_datafile
 character(len=*), intent(in)  :: object  ! name of the solar system object
 integer,          intent(out) :: ierr
 real :: elems(nelem)
 character(len=512) :: url
 character(len=30)  :: localfile,filepath
 logical :: iexist
 integer :: ierr2

 ierr = 0
 localfile = trim(object)//'.txt'
 inquire(file=localfile,exist=iexist)
 if (.not.iexist) then
    call construct_horizons_api_url(trim(object),url,ierr)
    if (ierr /= 0) then
       print*,' ERROR: could not query horizons database for object='//trim(object)
       return
    endif
    call download_datafile(url=url,dir=localfile,filename='',filepath=filepath,ierr=ierr)
 else
    filepath = localfile
 endif

 call read_ephemeris_file(filepath,elems,ierr2)
 if (ierr2 /= 0) ierr = ierr2

end function get_ephemeris

!-----------------------------------------------------------------------
!+
!  construct API call to the JPL Horizons ephemeris server
!+
!-----------------------------------------------------------------------
subroutine construct_horizons_api_url(object,url,ierr)
 character(len=*), intent(in)  :: object  ! name of the solar system object
 character(len=*), intent(out) :: url     ! url for query
 integer,          intent(out) :: ierr
 character(len=3)  :: cmd
 character(len=10) :: start_epoch,end_epoch
 integer           :: values(8),year,month,day

 ierr = 0
 select case(trim(adjustl(object)))
 case('pluto')
    cmd = '999' ! pluto barycentre
 case('neptune')
    cmd = '899' ! neptune barycentre
 case('uranus')
    cmd = '799' ! uranus barycentre
 case('saturn')
    cmd = '699' ! saturn barycentre
 case('jupiter')
    cmd = '599' ! jupiter barycentre
 case('mars')
    cmd = '499' ! mars barycentre
 case('earth')
    cmd = '399' ! earth-moon barycentre
 case('venus')
    cmd = '299' ! venus barycentre
 case('mercury')
    cmd = '199' ! mercury barycentre
 case('sun')
    cmd = '10'  ! el sol
 case('moon')
    cmd = '301' ! luna
 case default
    cmd = trim(adjustl(object))  ! can just request numbers directly
    ierr = 1
 end select

 !if (present(epoch)) then
!    start_epoch = epoch
 !else
 call date_and_time(values=values)
 year = values(1); month = values(2); day = values(3)
 write(start_epoch,"(i4.4,'-',i2.2,'-',i2.2)") year,month,day
 write(end_epoch,"(i4.4,'-',i2.2,'-',i2.2)") year,month,day+1
 !endif

 url = "'https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='"//trim(cmd)// &
       "'&OBJ_DATA='YES'&MAKE_EPHEM='YES'&EPHEM_TYPE='ELEMENTS'&CENTER='500@10'&START_TIME='"&
        //trim(start_epoch)//"'&STOP_TIME='"//trim(end_epoch)// &
       "'&STEP_SIZE='1DAYS'&REF_SYSTEM='ICRF'&REF_PLANE='ECLIPTIC''"

end subroutine construct_horizons_api_url

!-----------------------------------------------------------------------
!+
!  read required information from the ephemeris data file
!+
!-----------------------------------------------------------------------
subroutine read_ephemeris_file(file,elems,ierr)
 character(len=*), intent(in)  :: file
 real,             intent(out) :: elems(nelem)
 integer, intent(out) :: ierr
 integer :: iu,j
 character(len=80)  :: line
 character(len=*), parameter :: tag(nelem) = &
    (/'A                    ',&
      'EC                   ',&
      'IN                   ',&
      'OM                   ',&
      'W                    ',&
      'TA                   ',&
      'GM km^3/s^2          ',&
      'Vol. Mean Radius km  ',&
      'Density g/cm^3       '/)
 logical :: got_elem(nelem)
 character(len=len(tag)) :: alt_tag(2)

 ! give default parameters
 got_elem(:) = .false.
 elems(:) = 0.

 open(newunit=iu,file=file,status='old',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(file)
    return
 endif
 print "(/,a)",' > reading ephemeris from '//trim(file)
 do while(ierr == 0)
    read(iu,"(a)",iostat=ierr) line
    do j=1,nelem
       if (.not.got_elem(j)) then
          alt_tag(:) = tag(j)
          ! give alternatives as not all variables have consistent naming
          if (j==9) alt_tag = ['Density g cm^-3  ','Density R=1195 km']
          call read_value(line,tag(j),elems(j),got_elem(j),alt_tag=alt_tag)
       endif
    enddo
 enddo
 ierr = 0
 print "(a,/)"
 do j=1,nelem
    if (.not.got_elem(j)) then
      if (j <= 8) then
         ierr = ierr + 1
      endif
      print*,'ERROR: could not find '//trim(tag(j))//' in '//trim(file)
    endif
 enddo
 close(iu)

end subroutine read_ephemeris_file

!-----------------------------------------------------------------------
!+
!  utility routine to read a var = blah string from the ephemeris file
!+
!-----------------------------------------------------------------------
subroutine read_value(line,tag,val,got_val,alt_tag)
 use fileutils, only:string_delete,string_replace,lcase
 character(len=*), intent(in)    :: line, tag
 real,             intent(out)   :: val
 logical,          intent(inout) :: got_val
 character(len=*), intent(in), optional :: alt_tag(:)
 character(len=len(line)) :: string
 integer :: ieq,j,k,ierr,iplus

 val = 0.
 if (.not.got_val) then
    ! make a copy of the line string
    string = ' '//lcase(line)

    call string_delete(string,'(planet)') ! pluto has GM(planet) km^3/s^2

    ! delete double spaces
    call string_replace(string,'  ',' ')
    call string_delete(string,'(')
    call string_delete(string,')')
    call string_delete(string,',')

    ! add space at start of string if necessary
    if (string(1:1) /= ' ') string = ' '//trim(string)

!    print*,' string = '//trim(string)//' looking for '//trim(tag)

    ! search for ' TAG = value' in the line, including spaces
    ! to avoid confusion of A = val with TA= val
    j = max(index(string,' '//trim(lcase(tag))//' ='),&
        index(string,' '//trim(lcase(tag))//'='))

    ! try to match alternative tags if present
    if (j <= 0 .and. present(alt_tag)) then
       do k=1,size(alt_tag)
          j = max(index(string,' '//trim(lcase(alt_tag(k)))//' ='),&
                  index(string,' '//trim(lcase(alt_tag(k)))//'='))
          if (j/=0) exit
       enddo
    endif

    ! if we get a match read the value from what is after the equals sign
    if (j > 0) then
       string = string(j:)  ! start

       ! find the equals sign
       ieq   = index(string,'=')

       ! strip error bars e.g. 3389.92+-0.04
       iplus = index(string(ieq+1:),'+-')-2
       if (iplus<=0) iplus = len_trim(string(ieq+1:))

       read(string(ieq+1:ieq+1+iplus),*,iostat=ierr) val
       if (ierr == 0) then
          !write(*,"(1x,a,g0)",advance='no') trim(tag)//' = ',val
          ! mark this quantity as having already been read
          got_val = .true.
       else
          print*,' error reading '//trim(tag)//' from '//trim(string(ieq+1:ieq+1+iplus))
       endif
    endif
 endif

end subroutine read_value

end module ephemeris
