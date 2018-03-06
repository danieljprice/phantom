!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom2pdf
!
!  DESCRIPTION: This program is a utility for calculating volume-weighted
! PDFs (by interpolation to a 3D grid), direct from Phantom dumps
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom2pdf [npixx] dumpfilename(s)
!
!  DEPENDENCIES: boundary, dim, interpolations3D, io, part, pdfs,
!    readwrite_dumps, rhomach
!+
!--------------------------------------------------------------------------
program phantom2pdf
 use dim,              only:tagline
 use part,             only:xyzh,vxyzu,npart,massoftype,rhoh,hfact
 use io,               only:set_io_unit_numbers,iprint,idisk1,real4
 use boundary,         only:xmin,ymin,zmin,dxbound
 use readwrite_dumps,  only:read_dump
 use interpolations3D, only:interpolate3D
 use pdfs,             only:pdf_calc,pdf_write
 use rhomach,          only:get_rhomach_grid,write_rhomach
 implicit none
 integer :: npixx,npixy,npixz,nbins
 integer :: nargs,ierr,i,istart,ifile
 integer, parameter :: iunit = 71
 real, allocatable :: xval(:),pdf(:)
 real, allocatable :: datgrid(:,:,:,:) !rho, vx, vy, vz
 real, allocatable :: rhogrid(:)
 real :: weight,time,pixwidth,totvol
 real :: rhomin,rhomax,binspacing,pdfmin,pdfmax,rhologmin,rhologmax
 real :: dvol,vx2,vy2,vz2,rhoi,rmsv,rmsvmw
 real :: rhomeanvw,rhomeanmw,rhovarvw,rhovarmw
 real :: smeanvw,smeanmw,svarvw,svarmw
 character(len=120) :: dumpfile,fileout
 character(len=20)  :: string
 logical :: volweighted

 call set_io_unit_numbers
 iprint = 6
 npixx = 256
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantom2pdf [npixx] dumpfilename(s)'
    print "(a)",'        (default npixx=256)'
    stop
 endif
 print "(/,a,/)",' Phantom2pdf: we are pleased to do you service'

 call get_command_argument(1,string)
 read(string,*,iostat=ierr) npixx
 istart = 2
 if (ierr /= 0 .or. npixx < 1 .or. npixx > 100000) then
    print*,' npixx not read; assuming 256 '
    npixx = 256
    istart = 1
 endif
 npixy = npixx; npixz = npixx

 print "(a,i4,a)",' (allocating memory for ',npixx,'^3 grid)'
 allocate(datgrid(4,npixx,npixy,npixz))

 over_files: do ifile=istart,nargs

    call get_command_argument(ifile,dumpfile)
!
!--read particle setup from dumpfile
!
    call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) then
       print*,'error reading dumpfile'
       cycle over_files
    endif
!
! calculate quantities on the particles to check against the values on the grid
!
    rhomin = huge(rhomin)
    rhomax = 0.
    rhomeanvw = 0.
    totvol = 0.
    rmsv   = 0.
    print*,'hfact = ',hfact
    do i=1,npart
       rhoi = rhoh(xyzh(4,i),massoftype(1))
       dvol = massoftype(1)/rhoi
       rhomin = min(rhomin,rhoi)
       rhomax = max(rhomax,rhoi)
       vx2 = vxyzu(1,i)**2
       vy2 = vxyzu(2,i)**2
       vz2 = vxyzu(3,i)**2
       rmsv = rmsv + dvol*(vx2 + vy2 + vz2)
       totvol = totvol + dvol
       rhomeanvw = rhomeanvw + dvol*rhoi
    enddo
    rmsv = sqrt(rmsv/totvol)
    rhomeanvw = rhomeanvw/totvol
    print*,'On parts: Total Volume =',totvol,' rms v = ',rmsv
    print*,'On parts: Total Mass = ',massoftype(1)*npart,' mean dens = ',rhomeanvw
    print*,'On parts: max. dens = ',rhomax,' min. dens = ',rhomin

!---------------------------------------------
!
!--interpolate to the 3D grid
!
    pixwidth = dxbound/npixx
    weight = 1.d0/hfact**3
    print*,'using hfact = ',hfact,' weight = ',weight

    call interpolate3D(xyzh,weight,massoftype(1),vxyzu,npart,xmin,ymin,zmin,datgrid, &
                    npixx,npixy,npixz,pixwidth,.true.,0)
    !
    !--calculate rhomach quantities from the gridded data
    !
    rhologmin = -15.
    rhologmax = 15.
    call get_rhomach_grid(datgrid,npixx,npixy,npixz,rhologmin,&
                       rhomeanvw,rhovarvw,smeanvw,svarvw,rmsv,rhogrid)
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
    write(fileout,"(a,i3.3,a)") 'rhomach_'//trim(dumpfile)//'_grid',npixx,'.out'
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
    call pdf_write(nbins,xval,pdf,'lnrhogrid',volweighted,trim(dumpfile),trim(tagline))
    if (allocated(rhogrid)) deallocate(rhogrid)

 enddo over_files

 if (allocated(datgrid)) deallocate(datgrid)
 if (allocated(xval)) deallocate(xval)
 if (allocated(pdf)) deallocate(pdf)
!---------------------------------------------

 print "(/,a,/)",' Phantom2pdf: we thank you for your kind custom'

end program phantom2pdf
