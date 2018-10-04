!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: time_average_struct
!
!  DESCRIPTION: program to convert between structure function file formats
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: time_average_struct [outformat] infiles(*.sfn,*.struct,*_sf_*.dat)
!
!  DEPENDENCIES: io_structurefn
!+
!--------------------------------------------------------------------------
program time_average_struct
 use io_structurefn, only:get_sfn_format,nformats,labelformat,&
     read_structurefn,write_structurefn
 implicit none
 integer :: nargs,i,ierr
 integer, parameter :: iunitr = 1, iunitw = 2
 integer, parameter :: maxlag = 4096, maxorder=10
 character(len=120) :: infile,outfile
 integer :: n_lag, n_order,iorder,n_lag_prev,n_order_prev
 integer :: ioutformat,istart,nfiles
 real    :: rho_power,time
 real :: xlag(maxlag)
 real :: xlagprev(maxlag)
 real :: orders(maxorder)
 integer(kind=8) :: ncount(maxlag)
 real(kind=8) :: sf(2,maxorder,maxlag)
 real(kind=8) :: sf_av(2,maxorder,maxlag)
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
    ioutformat = get_sfn_format(infile,outfile)
    print "(a)",' Using default output format (same as input format) = '//trim(labelformat(ioutformat))
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

 print "(/,a,/)",' time_average_struct: we average your structure functions'

 n_order_prev = 0
 n_lag_prev   = 0
 xlagprev(:)  = 0.
 nfiles = 0
 sf_av = 0.d0
 over_args: do i=istart,nargs
    call get_command_argument(i,infile)

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
       n_order = n_order_prev
       n_lag   = n_lag_prev
       xlag    = xlagprev
       cycle over_args
    endif
    if (i > istart) then
       if (n_order /= n_order_prev) then
          print*,' error: n_order different in '//trim(infile)//', skipping'
          cycle over_args
       elseif (n_lag /= n_lag_prev) then
          print*,' error: n_lag different in '//trim(infile)//', skipping'
          cycle over_args
       endif
       if (any(abs(xlag(1:n_lag)-xlagprev(1:n_lag)) > tiny(xlag))) then
          print*,' WARNING: lag values differ in '//trim(infile)//': bins may be different!'
          stop
       endif
    endif

    sf_av(:,1:n_order,1:n_lag) = sf_av(:,1:n_order,1:n_lag) + sf(:,1:n_order,1:n_lag)
    nfiles = nfiles + 1

    n_order_prev = n_order
    n_lag_prev   = n_lag
    xlagprev(:)  = xlag(:)

    print*,'nlag = ',n_lag,' norder = ',n_order

 enddo over_args
!
!--divide sum by nfiles
!
 if (nfiles > 0) then
    sf_av = sf_av/real(nfiles)
 else
    stop ' ERROR: no files read '
 endif

 print*,' writing lag = ',xlag(1:10)
 print*,' writing n_order = ',n_order,' n_lag = ',n_lag
 outfile = 'time_average'
 call write_structurefn(ioutformat,trim(outfile),n_lag,xlag(1:n_lag),n_order,&
                        sf_av,rho_power,ncount(1:n_lag),orders(1:n_order),time,iunitw,ierr)

 print "(/,a,/)",' time_average_struct: may your structure functions scale happily'

contains

subroutine print_usage()
 use io_structurefn, only:print_sfn_formats

 print "(a)",' time_average_struct (c) 2009 Daniel Price'
 print "(a)",' produces time average from a series of structure function files'
 print "(/,a,/)",' Usage: time_average_struct [outformat] infiles(*.sfn,*.struct,*_sf_*.dat)'
 print "(a)",' possible values for [outformat]:'
 call print_sfn_formats()
 print "(a)",' input format determined from the file extension'

end subroutine print_usage

end program time_average_struct
