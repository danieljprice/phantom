!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: io_structurefn
!
!  DESCRIPTION:
!  module for read/write of structure functions to/from file
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: fileutils
!+
!--------------------------------------------------------------------------
module io_structurefn
 implicit none
 integer, parameter, public :: nformats = 4
 character(len=41), dimension(nformats) :: labelformat = &
      (/'ascii  _longd.struct _trans.struct       ', &
        'split  _longd1.struct _trans1.struct etc.', &
        'kitp   .sfn                              ', &
        'cfeder _SF_*.dat _sf_*                   '/)

 character(len=12), dimension(2) :: sftype = &
      (/'longitudinal','transverse  '/)

 integer, parameter, public :: mfile=60,mlabel=60
 integer, parameter, public :: power_unit = 14

 public :: read_structurefn,write_structurefn

contains

!----------------------------------------------------------------
!+
!  routine to print available file formats
!+
!----------------------------------------------------------------
subroutine print_sfn_formats()
 integer :: i

 do i=1,nformats
    if (i==1) then
       print "(i2,') ',a)",i,trim(labelformat(i))//'    (default if no [outformat] specified)'
    else
       print "(i2,') ',a)",i,trim(labelformat(i))
    endif
 enddo

end subroutine print_sfn_formats

!----------------------------------------------------------------
!+
!  routine to determine file format based on filenames
!
!  returns the stripped base filename for the file, without the
!  extension, in outfile. e.g. file_longd.struct would return
!  'file' as the stripped filename.
!+
!----------------------------------------------------------------
integer function get_sfn_format(filename,outfile)
 use fileutils, only:lcase
 character(len=*), intent(in)  :: filename
 character(len=*), intent(out) :: outfile
 integer :: l

 get_sfn_format = 0
 l = len_trim(filename)
 outfile = filename
 if (l > 12) then
    if (filename(l-12:l)=='_longd.struct') then
       get_sfn_format = 1
       outfile = filename(1:index(filename,'_longd.struct'))
    elseif (filename(l-12:l)=='_trans.struct') then
       get_sfn_format = 1
       outfile = filename(1:index(filename,'_trans.struct'))
    endif
 endif
 if (get_sfn_format==0) then
    l = index(filename,'.struct')
    if (l > 0) then
       get_sfn_format = 2
       outfile = filename(1:l-1)
    else
       l = index(filename,'.sfn')
       if (l > 0) then
          get_sfn_format = 3
          outfile = filename(1:l-1)
       else
          l = index(lcase(filename),'_sf_')
          if (l > 0) then
             get_sfn_format = 4
             outfile = filename(1:l-1)
          endif
          l = index(lcase(filename),'_sf_v')
          if (l > 0) then
             get_sfn_format = 4
             outfile = filename(1:l-1)
          endif
          l = index(lcase(filename),'_sf_rho')
          if (l > 0) then
             get_sfn_format = 4
             outfile = filename(1:l-1)
          endif
          l = index(lcase(filename),'_sf_sqrtrho')
          if (l > 0) then
             get_sfn_format = 4
             outfile = filename(1:l-1)
          endif
          l = index(lcase(filename),'_sf_lnrho')
          if (l > 0) then
             get_sfn_format = 4
             outfile = filename(1:l-1)
          endif
          l = index(lcase(filename),'_sf_unk')
          if (l > 0) then
             get_sfn_format = 4
             outfile = filename(1:l-1)
          endif
       endif
    endif
 endif

end function get_sfn_format

!----------------------------------------------------------------
!+
!  interface routine to read structure function of
!  any format
!+
!----------------------------------------------------------------
subroutine read_structurefn(infile,xlag,func,maxlag,maxorder,n_lag,n_order,rho_power,iunitr,ierr)
 character(len=*), intent(in)  :: infile
 real,             intent(out) :: xlag(:)
 real(kind=8),     intent(out) :: func(:,:,:)
 real,             intent(out) :: rho_power
 integer,          intent(out) :: n_lag,n_order,ierr
 integer,          intent(in)  :: maxlag,maxorder,iunitr

 character(len=len(infile)+10)  :: outfile
 integer :: informat
 real    :: xlagmax

 ierr = 0
 informat = get_sfn_format(infile,outfile)
 select case(informat)
 case(1)
    print*,' read from '//trim(labelformat(informat))//' not implemented'
    ierr = 1
    return
 case(2)
    print*,' read from '//trim(labelformat(informat))//' not implemented'
    ierr = 2
    return
 case(3)
    call read_sf(trim(infile),xlag,func,maxlag,maxorder,n_lag,n_order,rho_power,iunitr,ierr)
 case(4)
    call read_cfstruct(trim(infile),xlag,func,maxlag,maxorder,n_lag,n_order,rho_power,iunitr,ierr)
    !print*,' lag = ',xlag(1:10)
    xlagmax = real(int(maxval(xlag(1:n_lag)-0.5)*2.0/sqrt(3.)))
    if (xlagmax > 2.0) then
       print*,' converting lag values to actual distances (divide by ',xlagmax,')'
       xlag(1:n_lag) = xlag(1:n_lag)/xlagmax
    endif
 case default
    print*,' input file format for '//trim(infile)//' not recognised, skipping...'
    ierr = 3
    return
 end select

end subroutine read_structurefn

!----------------------------------------------------------------
!+
!  interface routine to write structure function of
!  any format
!+
!----------------------------------------------------------------
subroutine write_structurefn(ioutformat,outfile,n_lag,xlag,n_order,&
                             func,rho_power,ncount,orders,time,iunit,ierr)
 integer,          intent(in)  :: ioutformat
 character(len=*), intent(in)  :: outfile
 integer,          intent(in)  :: n_lag, n_order
 real,             intent(in)  :: xlag(n_lag)
 real(kind=8),     intent(in)  :: func(2,n_order,n_lag)
 real,             intent(in)  :: rho_power, time
 integer(kind=8),  intent(in)  :: ncount(n_lag)
 real,             intent(in)  :: orders(n_order)
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr

 character(len=mfile) :: origin

 ierr = 0


 select case(ioutformat)
 case(4)
    !
    !--write to Christoph Federrath format ascii files
    !
    call write_cfstruct(trim(outfile),n_lag,xlag(1:n_lag),n_order,&
                        ncount(1:n_lag),func,rho_power,time,iunit)
 case(3)
    !
    !--write to Aake format .sfn file
    !
    origin = trim(outfile)
    print "(a)",' writing to '//trim(outfile)//'.sfn'
    call openw_sf(trim(outfile)//'.sfn',origin,n_lag,xlag(1:n_lag),n_order,1)
    call write_sf(n_lag,n_order,func,orders(1:n_order),rho_power)
    close(power_unit)
 case(2)
    !
    !--split ascii format (one file per order)
    !
    call write_structfiles_split(trim(outfile),n_lag,xlag(1:n_lag),n_order,func,&
         rho_power,time,iunit)

 case default
    !
    !--default _longd.struct and _trans.struct format
    !
    call write_structfiles(trim(outfile),n_lag,xlag(1:n_lag),n_order,func,&
         rho_power,time,iunit)
 end select

end subroutine write_structurefn

!----------------------------------------------------------------
!+
!  routines to output structure functions in Aake format
!+
!----------------------------------------------------------------
subroutine openw_sf (file,origin,n_lag,lag,n_order,n_rho_power)
 character(len=*),     intent(in) :: file
 character(len=mfile), intent(in) :: origin
 integer,              intent(in) :: n_lag,n_order,n_rho_power
 real,                 intent(in) :: lag(n_lag)
 namelist /structurefn/n_lag,n_order,n_rho_power,origin
!
 open (power_unit,file=trim(file),status='unknown',form='formatted')                 ! open unit
 write (power_unit,structurefn)                                                ! dimensions info
 write (power_unit,'(1x,8g15.7)') lag                                          ! lag vector
end subroutine

subroutine write_sf (n_lag, n_order, f, orders, rho_power)
 integer,      intent(in) :: n_lag, n_order
 real,         intent(in) :: rho_power, orders(n_order)
 real(kind=8), intent(in) :: f(2,n_order,n_lag)
 integer :: i_order, i_sf
 real order
 namelist /sf/ i_order, order, i_sf, rho_power
!
 do i_sf = 1,2                                                                 ! make human-readable
    do i_order = 1,n_order                                                        ! make human-readable
       order = orders(i_order)
       write (power_unit,sf)                                                       ! namelist before each fn
       write (power_unit,'(8g15.7)') f(i_sf,i_order,:)  !f(:,i_order,i_sf)         ! f(lag;order;direction)
    enddo
 enddo
end subroutine
!----------------------------------------------------------------
!+
!  routines to read structure functions in Aake format
!+
!----------------------------------------------------------------
subroutine read_sf(sfnfile,xlag,func,maxlag,maxorder,n_lag,n_order,rho_power,iunitr,ierr)
 character(len=*), intent(in)  :: sfnfile
 real,             intent(out) :: xlag(:)
 real(kind=8),     intent(out) :: func(:,:,:)
 real,             intent(out) :: rho_power
 integer,          intent(out) :: n_lag,n_order,ierr
 integer,          intent(in)  :: maxlag,maxorder,iunitr
 character(len=120)                  :: origin
 integer :: i_order,n_rho_power,i_sf,ilag
 real    :: order,xlagmax

 namelist /structurefn/n_lag,n_order,n_rho_power,origin
 namelist /sf/ i_order, order, i_sf, rho_power

 open(unit=iunitr,file=sfnfile,form='formatted',status='old',iostat=ierr)
 if (ierr /= 0) then
    print*,' ERROR: could not open file '//trim(sfnfile)
    return
 endif
 print*,'reading file '//trim(sfnfile)
 read(iunitr,NML=structurefn,iostat=ierr)
 if (ierr  /=  0) then
    print*,'warning: error reading header'
    rewind(iunitr)
    do ilag=1,6
       read(iunitr,*)
    enddo
 endif
 print*,'nlag = ',n_lag,' norder = ',n_order,' n_rho_power = ',n_rho_power

 if (n_lag > maxlag) then
    print*,' WARNING: n_lag = ',n_lag,' exceeds array limits in read_sf'
    print*,'          reading only first ',maxlag,' lines'
    n_lag = maxlag
 endif
 if (n_order > maxorder) then
    print*,' WARNING: n_order = ',n_order,' exceeds array limits in read_sf'
    print*,'          reading only first ',maxorder,' orders'
    n_order = maxorder
 endif

 read(iunitr,*) xlag(1:n_lag)
 print*,' lag = ',xlag(1:n_lag)
 xlagmax = maxval(xlag(1:n_lag))
 if (xlagmax > 1.0) xlag(1:n_lag) = xlag(1:n_lag)*0.5/xlagmax
 ierr = 0
 do while (ierr == 0)
    read(iunitr,NML=sf,iostat=ierr)
    print*,'reading order = ',i_order,order,' i_sf = ',i_sf,' origin = '//trim(origin)
    print*,'rho_power = ',rho_power
    read(iunitr,*,iostat=ierr) func(i_sf,i_order,1:n_lag)
 enddo
 ierr = 0
 close(iunitr)

end subroutine read_sf

!----------------------------------------------------------------
!+
!  routine to write structure functions to file
!  in simple ascii format (all sfns in one file)
!+
!----------------------------------------------------------------
subroutine write_structfiles(basename,n_lag,lag,n_order,f,rho_power,time,iunit)
 character(len=*), intent(in) :: basename
 integer,          intent(in) :: n_lag, n_order
 real,             intent(in) :: lag(n_lag)
 real(kind=8),     intent(in) :: f(2,n_order,n_lag)
 real,             intent(in) :: rho_power, time
 integer,          intent(in) :: iunit
 integer :: i_sf, ilag, ierr, j
 character(len=120) :: filename

 do i_sf = 1,2                                                                 ! make human-readable
    if (i_sf==1) then
       filename = trim(basename)//'_longd.struct'
    else
       filename = trim(basename)//'_trans.struct'
    endif
    open(iunit,file=trim(filename),status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR writing to '//trim(filename)
    else
       write(*,"(a)",advance='no') ' writing to '//trim(filename)//'...'
!
!--write file header
!
       write (iunit,"(a,f10.3,a)",iostat=ierr) &
             '#  structure function file, rho_power = ',rho_power,': created from file '//trim(basename)
       if (ierr /= 0) print "(a)",'ERROR!'

       write (iunit,"(a,8g12.4)",iostat=ierr) &
             '#  time in file = ',time
       if (ierr /= 0) print "(a)",'ERROR!'
!
!--write lag and structure function to file
!
       do ilag=1,n_lag
          write (iunit,'(20(8g15.7))',iostat=ierr) lag(ilag),(f(i_sf,j,ilag),j=1,n_order)
          if (ierr /= 0) print "(a)",'ERROR!'
       enddo
       print*,'ok'
       close(iunit)
    endif
 enddo

end subroutine write_structfiles

!----------------------------------------------------------------
!+
!  routine to write structure functions to sequence of files
!  in simple ascii format (one file per structure function)
!+
!----------------------------------------------------------------
subroutine write_structfiles_split(basename,n_lag,lag,n_order,f,rho_power,time,iunit)
 character(len=*), intent(in) :: basename
 integer,          intent(in) :: n_lag, n_order
 real,             intent(in) :: lag(n_lag)
 real(kind=8),     intent(in) :: f(2,n_order,n_lag)
 real,             intent(in) :: rho_power, time
 integer,          intent(in) :: iunit
 integer :: i_sf, ilag, ierr, i_order
 character(len=120) :: filename

 do i_sf = 1,2
    do i_order=1,n_order
       if (i_sf==1) then
          write(filename,"(a,i2.2,a)") trim(basename)//'_longd',i_order,'.struct'
       else
          write(filename,"(a,i2.2,a)") trim(basename)//'_trans',i_order,'.struct'
       endif
       open(iunit,file=trim(filename),status='replace',form='formatted',iostat=ierr)
       if (ierr /= 0) then
          print "(a)",' ERROR writing to '//trim(filename)
       else
          write(*,"(a)",advance='no') ' writing to '//trim(filename)//'...'
          !
          !--write file header
          !
          write (iunit,"(a,i2,a,f10.3,a)",iostat=ierr) &
                '#  structure function file, order ',i_order,' rho_power = ',rho_power,&
                ': created from file'//trim(basename)
          if (ierr /= 0) print "(a)",'ERROR!'

          write (iunit,"(a,8g12.4)",iostat=ierr) &
                '#  time in file = ',time
          if (ierr /= 0) print "(a)",'ERROR!'
          !
          !--write lag and structure function to file
          !
          do ilag=1,n_lag
             write (iunit,'(20(8g15.7))',iostat=ierr) lag(ilag),f(i_sf,i_order,ilag)
             !if (i_order==1 .and. i_sf==1) print*,' writing line ',ilag,'order 1 sf = ',f(1,1,ilag)
          enddo
          print*,'ok'
       endif
       close(iunit)
    enddo
 enddo

end subroutine write_structfiles_split

!----------------------------------------------------------------
!+
!  routine to write structure functions in Christoph Federrath's
!  analysis code format (ascii)
!+
!----------------------------------------------------------------
subroutine write_cfstruct(basename,n_lag,lag,n_order,n_count,f,rho_power,time,iunit)
 use fileutils, only:lcase
 character(len=*), intent(in) :: basename
 integer,          intent(in) :: n_lag, n_order
 integer(kind=8),  intent(in) :: n_count(n_lag)
 real,             intent(in) :: lag(n_lag)
 real(kind=8),     intent(in) :: f(2,n_order,n_lag)
 real,             intent(in) :: rho_power, time
 integer,          intent(in) :: iunit
 integer            :: i_sf, ilag, ierr, j
 character(len=120) :: filename
 real, parameter    :: tol = 1.e-3
!
!--give file the correct extension if one is not already there
!
 if (index(lcase(basename),'_sf_') > 0) then
    filename = trim(basename)
 else
    if (abs(rho_power-1./3.) <= tol) then      ! rho_power = 0.33
       filename = trim(basename)//'_SF_rho3.dat'
    elseif (abs(rho_power-0.5) <= tol) then   ! rho_power = 0.5
       filename = trim(basename)//'_SF_sqrtrho.dat'
    elseif (abs(rho_power-1.) <= tol) then    ! rho_power = 1
       filename = trim(basename)//'_SF_rho.dat'
    elseif (abs(rho_power+1.) <= tol) then    ! rho_power = -1
       filename = trim(basename)//'_SF_lnrho.dat'
    elseif (abs(rho_power) <= tol) then       ! rho_power = 0
       filename = trim(basename)//'_SF_vels.dat'
    else                                        ! rho_power = unknown
       filename = trim(basename)//'_SF_unknown.dat'
    endif
 endif

 open(iunit,file=trim(filename),status='replace',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR writing to '//trim(filename)
 else
    write(*,"(a)",advance='no') ' writing to '//trim(filename)//'...'
!
!--write lag and structure function to file
!
    do ilag=1,n_lag
       write (iunit,'(i8,1x,50(8g15.7))',iostat=ierr) ilag,lag(ilag),lag(ilag), &
             ((n_count(j),f(i_sf,j,ilag),i_sf=1,2),j=1,n_order)
       if (ierr /= 0) print "(a)",'ERROR!'
    enddo
    print*,'ok'
 endif
 close(iunit)

 return
end subroutine write_cfstruct

!----------------------------------------------------------------
!+
!  routine to read structure functions from Christoph Federrath's
!  analysis code format (ascii)
!+
!----------------------------------------------------------------
subroutine read_cfstruct(filename,xlag,sf,maxlag,maxorder,nlag,norder,rho_power,iunitr,ierr)
 use fileutils, only:get_ncolumns,skip_header,lcase
 character(len=*), intent(in)  :: filename
 real(kind=8),     intent(out) :: sf(:,:,:)
 real,             intent(out) :: xlag(:)
 integer,          intent(out) :: norder,nlag,ierr
 integer,          intent(in)  :: maxlag,maxorder,iunitr
 real,             intent(out) :: rho_power

 integer :: ilag,idum,iorder,ncol,nhead
 real    :: dum
!
!--work out how many columns there are in the file
!
 open(unit=iunitr,file=trim(filename),form='formatted',status='old',iostat=ierr)
 if (ierr /= 0) then
    print*,' ERROR opening '//trim(filename)//' for read'
    return
 endif
 idum = index(lcase(filename),'_sf_')
 select case(trim(lcase(filename(idum+4:))))
 case('v','vels.dat')
    rho_power = 0.
 case('rho3','rho3.dat')
    rho_power = 1./3.
 case('sqrtrho','sqrtrho.dat')
    rho_power = 0.5
 case('rho','rho.dat')
    rho_power = 1.0
 case('lnrho')
    rho_power = -1.0
 case default
    print*,' ERROR! unknown rho_power from string '//trim(lcase(filename(idum+5:)))
    rho_power = 66.
 end select

 print*,'reading file '//trim(filename)//' rho_power = ',rho_power
 call get_ncolumns(iunitr,ncol,nhead)
 norder = (ncol - 3)/4
 print*,'inferring norder = ',norder,' from number of columns'
 if (norder <= 1 .or. norder > 100) then
    print*,' ERROR calculating number of structure function orders (got ',norder,') - aborting'
    ierr = 1
    close(iunitr)
    return
 endif
 call skip_header(iunitr,nhead)

 ilag = 0
 ierr = 0
 do while (ierr==0 .and. ilag < maxlag)
    ilag = ilag + 1
    read(iunitr,*,iostat=ierr) idum,xlag(ilag),dum, &
    (dum,sf(1,iorder,ilag),dum,sf(2,iorder,ilag),iorder=1,norder)
    !print*,' got line ',ilag,'order 1 sf = ',sf(1,1,ilag),' ierr = ',ierr
 enddo
 close(iunitr)
 !print*,'ilag = ',ilag,ierr
 if (ierr > 0) then
    print*,'WARNING! errors reading structure functions from '//trim(filename)
 elseif (ierr < 0) then  ! reached end of file, the expected scenario
    nlag = ilag - 1
    ierr = 0
 else  ! run off array limits, but no errors
    nlag = ilag
 endif
 if (nlag >= maxlag) then
    print*,'WARNING! array size too small, read only ',maxlag,' lines'
    nlag = maxlag
 endif

end subroutine read_cfstruct

end module io_structurefn
