!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: get_struct_slope
!
!  DESCRIPTION: program to calculate structure function slopes
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: get_struct_slope infile(*.sfn,*.struct,*_sf_*.dat)
!
!  DEPENDENCIES: io_structurefn, leastsquares
!+
!--------------------------------------------------------------------------
program get_struct_slope
 use io_structurefn, only:read_structurefn,sftype
 use leastsquares,   only:fit_slope
 implicit none
 integer            :: nargs,i,ierr
 integer, parameter :: iunitr = 1, iunitw = 23
 integer, parameter :: maxlag = 4096, maxorder=10
 character(len=120) :: infile
 integer            :: n_lag,n_order,iorder,i_sf,itries
 logical            :: iexist
 real               :: rho_power,xminfit,xmaxfit
 real :: xlag(maxlag)
 real(kind=8) :: sf(2,maxorder,maxlag)
 real :: err,errslope,erryint,slope,yint,slope3longd,slope3trans
 real :: errmax,dxminfit,xminfiti,xmaxfiti
 real :: errmaxsf(2)

 nargs = command_argument_count()
 if (nargs < 1) then
    call print_usage
    stop
 endif

 print "(/,a,/)",' get_struct_slope: your slope is our business'

 do i=1,4
    select case(i)
    case(1)
       infile = 'slope_longd.out'
    case(2)
       infile = 'slope_trans.out'
    case(3)
       infile = 'slope_longd.func'
    case(4)
       infile = 'slope_trans.func'
    end select
    inquire(file=trim(infile),exist=iexist)
    if (iexist) then
       print*,'ERROR: '//trim(infile)//' already exists: delete or rename this file and try again'
       ierr = 1
    else
       open(unit=iunitw+i,file=trim(infile),status='new',form='formatted',iostat=ierr)
       if (ierr /= 0) then
          print*,' ERROR: could not open '//trim(infile)//' for output'
       endif
    endif
 enddo
 if (ierr /= 0) then
    print*,'use "rm slope_longd.out slope_trans.out slope_longd.func slope_trans.func"'
    stop
 endif

 over_args: do i=1,nargs
    call get_command_argument(i,infile)

    n_lag = 0
    n_order = 0

    call read_structurefn(trim(infile),xlag,sf,maxlag,maxorder,n_lag,n_order,rho_power,iunitr,ierr)
    !
    !--check for errors/inconsistencies in read
    !
    if (ierr /= 0) then
       print*,' error reading structure function from '//trim(infile)//', skipping'
       cycle over_args
    endif
    print*,'nlag = ',n_lag,' norder = ',n_order

    write(iunitw+1,"('# ',a5,5(a12,1x),a,a)") 'order','slope','yint','errslope','erryint','rfac',&
                                              'f(x) for splash','  [from '//trim(infile)//']'
    write(iunitw+2,"('# ',a5,5(a12,1x),a,a)") 'order','slope','yint','errslope','erryint','rfac',&
                                              'f(x) for splash','  [from '//trim(infile)//']'

    xminfit = 0.10 !5 !01 !0.022
    xmaxfit = 0.15 !0.22
    dxminfit = 0.005

    if (index(infile,'128') /= 0) then
       !xminfit = xminfit*4.
       print*,' adjusting xminfit for 128^3, xminfit = ',xminfit
    elseif (index(infile,'256') /= 0) then
       !xminfit = xminfit*2.
       print*,' adjusting xminfit for 256^3, xminfit = ',xminfit
    endif

    !--get slope of 3rd order structure fn
    call fit_slope(n_lag,xlag(1:n_lag),sf(1,3,1:n_lag),&
                   slope3longd,yint,err,errslope,erryint,xmin=xminfit,xmax=xmaxfit,logplot=.true.)
    call fit_slope(n_lag,xlag(1:n_lag),sf(2,3,1:n_lag),&
                   slope3trans,yint,err,errslope,erryint,xmin=xminfit,xmax=xmaxfit,logplot=.true.)

    !--get slope of 1st order structure fns in the default range
    call fit_slope(n_lag,xlag(1:n_lag),sf(1,1,1:n_lag),&
                   slope,yint,err,errslope,erryint,xmin=xminfit,xmax=xmaxfit,logplot=.true.)
    errmaxsf(1) = 1.-err    ! record the error of the fit
    call fit_slope(n_lag,xlag(1:n_lag),sf(2,1,1:n_lag),&
                   slope,yint,err,errslope,erryint,xmin=xminfit,xmax=xmaxfit,logplot=.true.)
    !--record the error of the fit
    errmaxsf(2) = 1.-err   ! record the error of the fit
!    print*,'ERRMAX = ',errmaxsf(:)

    do i_sf=1,2
       !
       !--start by fitting the full allowable range
       !
       xminfiti = xminfit
       xmaxfiti = xmaxfit
       errmax   = errmaxsf(i_sf)
       do iorder=1,n_order
          !
          !--try making the fitting range shorter until the error is less than or equal
          !  to the error in the fit for the 1st order structure function
          !
          err = 0.
          itries = 1
          !xminfiti = xminfit
          print "(1x,a,i2)",trim(sftype(i_sf))//' order ',iorder
          !do while(xminfiti > 0.005)
          call fit_slope(n_lag,xlag(1:n_lag),sf(i_sf,iorder,1:n_lag),&
                         slope,yint,err,errslope,erryint,xmin=xminfiti,xmax=xmaxfiti,logplot=.true.)
          !   xminfiti = xminfiti - dxminfit
          !   write(i_sf*100+iorder,*) xmaxfiti-xminfiti,slope,1.-err,xminfiti
          !   print*,'iorder = ',iorder,' slope = ',slope,'err = ',1.-err,' compared to ',errmax,' trying xminfit = ',xminfiti
          !enddo
          print*,' best fit slope = ',slope,' +/- ',errslope,' overall goodness of fit = ',100.*err,'%'
          print*,' y intercept    = ',yint,' +/- ',erryint,' unlogged = ',10**(yint),' angle = ',180*atan(slope)/3.1415936d0
          !print "(a,es10.4,a,es10.4)",' function string: 10**(m*log10(x) + yint),m=',slope,',yint=',yint
          !print*,' writing to output file'
          if (i_sf==1) then
             write(iunitw+1,"(i5,2x,6(es12.4,1x),2(a,es10.4))",iostat=ierr) iorder,slope,yint,errslope,erryint,err,&
               slope/slope3longd,'10**(m*log10(x) + yint),m=',slope,',yint=',yint
             write(iunitw+3,"(2(a,es10.4))",iostat=ierr) '10**(m*log10(x) + yint),m=',slope,',yint=',yint
          else
             write(iunitw+2,"(i5,2x,6(es12.4,1x),2(a,es10.4))",iostat=ierr) iorder,slope,yint,errslope,erryint,err,&
               slope/slope3trans,'10**(m*log10(x) + yint),m=',slope,',yint=',yint
             write(iunitw+4,"(2(a,es10.4))",iostat=ierr) '10**(m*log10(x) + yint),m=',slope,',yint=',yint
          endif
       enddo
    enddo
 enddo over_args

 do i=1,4
    close(iunitw+i)
 enddo
 print "(/,a,/)",' get_struct_slope: may your structure functions scale happily'

contains

subroutine print_usage()
 use io_structurefn, only:print_sfn_formats

 print "(a)",' get_struct_slope (c) 2009 Daniel Price'
 print "(a)",' calculates best fit slope for structure functions'
 print "(/,a,/)",' Usage: get_struct_slope infile(*.sfn,*.struct,*_sf_*.dat)'
 print "(a)",' input format is determined from the file extension as follows:'
 call print_sfn_formats

end subroutine print_usage

end program get_struct_slope
