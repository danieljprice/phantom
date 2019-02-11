!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: grid2pdf
!
!  DESCRIPTION: This program is a utility for calculating volume-weighted
! PDFs (similar to phantom2pdf, but from already gridded data)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: grid2pdf format filename(s)
!
!  DEPENDENCIES: dim, io, io_grid, pdfs, rhomach
!+
!--------------------------------------------------------------------------
program grid2pdf
 use dim,              only:tagline
 use io,               only:set_io_unit_numbers,iprint,idisk1,real4
 use pdfs,             only:pdf_calc,pdf_write
 use rhomach,          only:get_rhomach_grid,write_rhomach
 use io_grid,          only:print_grid_formats,read_grid_header,read_grid_column
 implicit none
 integer :: npixx,npixy,npixz,nbins,ncols
 integer :: nargs,ierr,ifile,i,j,k,n
 integer, parameter :: iunit = 71
 real, allocatable :: xval(:),pdf(:)
 real, allocatable :: datgrid(:,:,:,:) !rho, vx, vy, vz
 real, allocatable :: rhogrid(:)
 real :: time
 real :: binspacing,pdfmin,pdfmax,rhologmin,rhologmax
 real :: rmsv,rmsvmw,rhomean,xtmp
 real :: rhomeanvw,rhomeanmw,rhovarvw,rhovarmw
 real :: smeanvw,smeanmw,svarvw,svarmw
 character(len=120) :: gridfile,fileout,tagline
 character(len=20)  :: gridformat
 logical :: volweighted,iexist

 call set_io_unit_numbers
 iprint = 6
 npixx = 256
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 2) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: grid2pdf format filename(s)'
    print "(a)",' Available formats:'
    call print_grid_formats()
    stop
 endif
 print "(/,a,/)",' Grid2pdf: we are pleased to do you service'
 tagline = 'Grid2pdf: (c) Daniel Price 2007-2010 (uses SPLASH pdf module)'

 call get_command_argument(1,gridformat)

 over_files: do ifile=2,nargs

    call get_command_argument(ifile,gridfile)
!
!--read header to determine size of array to allocate
!
    call read_grid_header(gridformat,gridfile,npixx,npixy,npixz,ncols,time,ierr)
    if (ierr /= 0) then
       print*,'error reading dumpfile'
       cycle over_files
    endif
    print*,'finished header read, nx,ny,nz = ',npixx,npixy,npixz,' ncols = ',ncols,' ierr = ',ierr

    print "(a,2(i4,' x'),i4,a)",' (allocating memory for ',npixx,npixy,npixz,' grid)'
    if (.not.allocated(datgrid)) allocate(datgrid(1,npixx,npixy,npixz))
    if (.not.allocated(rhogrid)) allocate(rhogrid(npixx*npixy*npixz))
    rhogrid(:) = 0.
    !
    !--read grid data from file
    !
    call read_grid_column(gridformat,gridfile,1,npixx,npixy,npixz,rhogrid,ierr)

    rhomean = sum(rhogrid)/size(rhogrid)
    rhogrid = rhogrid/rhomean
    print*,'rhomean = ',rhomean
    n = 0
    do k=1,npixz
       do j=1,npixy
          do i=1,npixx
             n = n + 1
             datgrid(1,i,j,k) = rhogrid(n)
          enddo
       enddo
    enddo
    !print*,'got rhogrid = ',rhogrid(1:10)

    !
    !--calculate rhomach quantities from the gridded data
    !
    rhologmin = -15.
    rhologmax = 15.
    call get_rhomach_grid(datgrid,npixx,npixy,npixz,rhologmin,&
                       rhomeanvw,rhovarvw,smeanvw,svarvw,rmsv,rhogrid)
    !
    !--if we can extract the file number,
    !  look for a file called "mach.evol" and lookup rmsv in this file
    !
    i = index(gridfile,'_dens')
    if (i > 0) then
       read(gridfile(i-4:i-1),*,iostat=ierr) n
       if (ierr==0 .and. n > 10 .and. n < 500) then
          print*,' got dump number = ',n
          inquire(file='Mach.evol',exist=iexist)
          if (iexist) then
             open(unit=iunit,file='Mach.evol',status='old',form='formatted',iostat=ierr)
             read(iunit,*,iostat=ierr) ! skip header line
             j = 0
             do while (ierr==0 .and. j < n)
                read(iunit,*,iostat=ierr) j,xtmp
                if (j==n .and. ierr==0) then
                   rmsv = xtmp
                   time = real(n)
                   print*,'FOUND entry in Mach.evol: got ',j,'rmsv = ',rmsv
                endif
             enddo
             close(unit=iunit)
          else
             print*,' error: Mach.evol file does not exist'
          endif
       else
          print*,' could not extract dump number from filename '//gridfile(i-4:i-1)//': no Mach.evol lookup'
       endif
    endif

    !
    !--fill in the blanks
    !
    rhomeanmw = 0.
    rhovarmw  = 0.
    svarmw    = 0.
    smeanmw   = 0.
    rmsvmw    = 0.
    !
    !--write line to rhomach file
    !
    write(fileout,"(a,i3.3,a)") 'rhomach_'//trim(gridfile)//'.out'
    call write_rhomach(trim(fileout),time,rhomeanvw,rhomeanmw,rhovarvw,rhovarmw,&
                    rmsv,rmsvmw,smeanvw,smeanmw,svarvw,svarmw)
!
!--allocate memory for PDF calculation
!
    binspacing = 0.1
    nbins = nint((rhologmax - rhologmin)/binspacing)
    print "(a,i3,a)",' (allocating memory for ',nbins,' PDF bins)'
    if (.not.allocated(xval)) allocate(xval(nbins),pdf(nbins))

    !
    !--calculate PDF of lnrho
    !
    call pdf_calc(npixx*npixy*npixz,rhogrid,rhologmin,rhologmax,nbins,xval,pdf,pdfmin,pdfmax,.true.,volweighted,ierr)

    xval = exp(xval)
    call pdf_write(nbins,xval,pdf,'lnrhogrid',volweighted,trim(gridfile),trim(tagline))
    if (allocated(rhogrid)) deallocate(rhogrid)

 enddo over_files

 if (allocated(datgrid)) deallocate(datgrid)
 if (allocated(xval)) deallocate(xval)
 if (allocated(pdf)) deallocate(pdf)
!---------------------------------------------

 print "(/,a,/)",' Grid2pdf: we hope you had a nice griddy time'

end program grid2pdf
