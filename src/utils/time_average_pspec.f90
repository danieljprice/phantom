!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: time_average_pspec
!
!  DESCRIPTION:
!  Utility to get the time averaged power spectra from the
!  .pow files produced the phantom2power utility
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: time_average_pspec compfac *.pow
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
program time_average_pspec
 implicit none
 integer, parameter :: iunit = 10, iout = 7, iprint = 6
 integer, parameter :: maxbins = 4096, maxfiles= 200
 integer :: nargs,iarg,nbins,nbinsprev,nfiles,i,ierr
 integer :: nheaderlines,ncolumns
 character(len=120) :: filename,line
 real :: xval(maxbins),xvalprev(maxbins),ave(maxbins),var(maxbins)
 real :: pspecval(maxbins,maxfiles)
 real :: dum,val,compfac

 nargs = command_argument_count()
 if (nargs < 2) then
    print "(a)",'A utility for computing the time averaged Power Spectrum'
    print "(a)",'from a bunch of ascii .pow files as output by phantom2power'
    print "(/,a)",'Usage: time_average_pspec compfac *.pow'
    stop
 elseif (nargs > maxfiles) then
    print "(a)",' ERROR: number of files exceeds array limits, edit and recompile'
    stop
 endif

 nbinsprev = 0
 xvalprev = 0.
 ave = 0.
 pspecval = 0.
 xval = 0.
 nfiles = 0
 call get_command_argument(1,filename)
 read(filename,*) compfac
 if (compfac < 0. .or. compfac > 5.) then
    stop 'error in specifying the compensation factor'
 else
    print*,' COMPENSATION FACTOR = ',compfac
 endif

 over_files: do iarg=2,nargs
    call get_command_argument(iarg,filename)

    open(unit=iunit,file=filename,iostat=ierr,status='old')
    if (ierr /= 0) then
       print "(a)",' ERROR: cannot open '//trim(filename)//' for read'
    else
       nbins = 0
       !--specify formatting of ascii file - should really determine this automatically
       !  but this works so long as the formatting of the .pow file does not change
       nheaderlines = 3
       ncolumns = 8
       !--skip header lines
       do i=1,nheaderlines
          read(iunit,"(a)",iostat=ierr)
       enddo
       !--read the columns in the file
       do while (ierr==0 .and. nbins < maxbins)
          if (ierr==0) then
             nbins = nbins + 1
             !--read all the columns to deal with crappy ifort ascii output which is spread over multiple lines
             read(iunit,*,iostat=ierr) xval(nbins),val,(dum,i=3,ncolumns)
             pspecval(nbins,nfiles+1) = val*xval(nbins)**compfac
          endif
       enddo
       if (ierr==0) print "(a,i4,a,i4,a)",' ERROR! number of bins ',nbins,' exceeds maximum (',maxbins,')'
       write(iprint,*) 'nbins = ',nbins

       !--error checks
       if (nbins <= 0) then
          print "(a)",' ERROR: no data read from file, skipping'
       endif

       if (nbinsprev > 0) then
          if (nbins /= nbinsprev) then
             print "(a)",' ERROR: number of bins has changed between files, skipping file'
             cycle over_files
          elseif (.not.all(xval(1:nbins)==xvalprev(1:nbins))) then
             print "(a)",' ERROR: location of bins has changed between files, skipping file'
             cycle over_files
          endif
          xvalprev = xval
       endif
       nfiles = nfiles + 1
       ave(1:nbins) = ave(1:nbins) + pspecval(1:nbins,nfiles)
    endif
 enddo over_files

 !--compute average
 write(iprint,*) 'nfiles = ',nfiles
 ave(1:nbins) = ave(1:nbins)/real(nfiles)

 !--compute standard deviation
 var = 0.
 do iarg=1,nfiles
    var(1:nbins) = var(1:nbins) + (pspecval(1:nbins,iarg) -  ave(1:nbins))**2
 enddo
 var(1:nbins) = var(1:nbins)/nfiles

 open(unit=iout,file='averaged_pspec.pow',status='replace',form='formatted')
 write(iout,"(a)") '# [01  xval ] [02  pspec(average)] [03 st. dev] [04 variance]'
 do i=1,nbins
    if (ave(i) > 0.) then ! only spit out non-zero bins
       write(iout,"(4(es14.8,2x))") xval(i),ave(i),sqrt(var(i)),var(i)
    endif
 enddo
 close(iout)

end program time_average_pspec
