!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: struct2struct
!
!  DESCRIPTION: program to convert between structure function file formats
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: struct2struct [outformat] infiles(*.sfn,*.struct,*_sf_*.dat)
!
!  DEPENDENCIES: io_structurefn
!+
!--------------------------------------------------------------------------
program struct2struct
 use io_structurefn, only:get_sfn_format,nformats,labelformat,&
     write_structurefn,read_structurefn
 implicit none
 integer :: nargs,i,ierr
 integer, parameter :: iunitr = 1, iunitw = 2
 integer, parameter :: maxlag = 4096, maxorder=10
 character(len=120) :: infile,outfile
 integer :: n_lag, n_order,iorder
 integer :: informat,ioutformat,istart
 real    :: rho_power,time
 real :: xlag(maxlag)
 real :: orders(maxorder)
 integer(kind=8) :: ncount(maxlag)
 real(kind=8) :: sf(2,maxorder,maxlag)
 logical :: iexist

 ioutformat = 0
 nargs = command_argument_count()
 if (nargs < 1) then
    call print_usage
    stop
 endif

 !
 !--the first command line argument is either a filename or the
 !  output file format. We assume it is a filename if the file exists
 !
 call get_command_argument(1,infile)

 inquire(file=trim(infile),exist=iexist)
 if (iexist .and. get_sfn_format(infile,outfile) > 0) then
    ioutformat = 1
    print "(a)",' Using default output format = '//trim(labelformat(ioutformat))
    istart = 1
 else
    istart = 2
    read(infile,*,iostat=ierr) ioutformat
    if (ierr /= 0 .or. ioutformat <= 0 .or. ioutformat > nformats) then
       do i=1,nformats
          if (trim(infile(1:6))==trim(labelformat(i)(1:6))) then
             ioutformat = i
          endif
       enddo
    endif
 endif
 if (ioutformat <= 0 .or. ioutformat > nformats) then
    call print_usage
    stop
 endif

 print "(/,a,/)",' struct2struct: we structure your functions'

 over_args: do i=istart,nargs
    call get_command_argument(i,infile)
    informat = get_sfn_format(infile,outfile)

    if (informat == ioutformat) then
       print "(a)",' ERROR: input format = output format for '//trim(infile)
       print "(a)",'        We don''t allow this to prevent accidental overwrites.'
       print "(a)",'        ...doing nothing...'
       cycle over_args
    endif
    time = 0.
    n_lag = 0
    n_order = 0
    ncount(:) = 0
    do iorder=1,maxorder
       orders(iorder) = iorder
    enddo

    call read_structurefn(trim(infile),xlag,sf,maxlag,maxorder,n_lag,n_order,rho_power,iunitr,ierr)
    !
    !--check for errors/inconsistencies in read
    !
    if (ierr /= 0) then
       print*,' error reading structure function from '//trim(infile)//', skipping'
       cycle over_args
    endif

    print*,'nlag = ',n_lag,' norder = ',n_order

    if (ierr==0) then
       call write_structurefn(ioutformat,trim(outfile),n_lag,xlag(1:n_lag),n_order,&
                              sf,rho_power,ncount(1:n_lag),orders(1:n_order),time,iunitw,ierr)
    else
       print*,' ERRORS READING FILE: NO OUTPUT WRITTEN'
    endif
 enddo over_args

 print "(/,a,/)",' struct2struct: may your structure functions scale happily'

contains

subroutine print_usage()
 use io_structurefn, only:print_sfn_formats

 print "(a)",' struct2struct (c) 2009 Daniel Price'
 print "(a)",' Converts between structure function file formats'
 print "(/,a,/)",' Usage: struct2struct [outformat] infiles(*.sfn,*.struct,*_sf_*.dat)'
 print "(a)",' possible values for [outformat]:'
 call print_sfn_formats()
 print "(a)",' input format determined from the file extension'

end subroutine print_usage

end program struct2struct
