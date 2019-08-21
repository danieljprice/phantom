!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Computes the specific energy distribution
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: infile_utils, io, part, physcon, sortutils
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'tde'
 public :: do_analysis

 private

 integer :: nbins
 real    :: mh
 integer, parameter :: nmaxbins = 5000
 real    :: rmax

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,         only:warning,iprint
 use dump_utils, only:read_array_from_file
 use prompting,  only:prompt
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 character(len=120) :: output
 character(len=20)  :: tdeprefix,tdeparams
 real, dimension(nmaxbins) :: ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,lbins,dndl
 integer :: i,iline,ierr
 logical :: iexist
 real(4) :: luminosity(npart)
 logical, save :: first = .true.

 call read_array_from_file(123,dumpfile,'luminosity',luminosity(1:npart),ierr)
 if (ierr/=0) then
    write(*,*)
    write(*,'("WARNING: could not read luminosity from file. It will be set to zero")')
    write(*,*)
    luminosity = 0.
 else if (first) then
    rmax  = 700.
    call prompt('Enter value for rmax when binning luminosity',rmax,0.)
    first = .false.
 endif

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'hist',numfile
 write(*,'("Output file name is ",A)') output

 ! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G == 1")')

 ! Read black hole mass from .tdeparams file
 iline = index(dumpfile,'_')
 tdeprefix = dumpfile(1:iline-1)
 tdeparams = trim(tdeprefix)//'.tdeparams'
 inquire(file=tdeparams, exist=iexist)
 ierr =1
 if (iexist) then
    call read_tdeparams(tdeparams,mh,iunit,ierr)
    if (ierr /= 0) call warning('analysis','could not open/read '//trim(tdeparams))
 else
    call warning('analysis','file '//trim(tdeparams)//' does not exist')
 endif

 if (.not.iexist .or. ierr/=0) then
    write(iprint,*) 'Assuming default black hole mass'
    mh = 1.
 endif

! Print out the parameters
 write(*,*)
 write(*,'("Parameters are:")')
 write(*,*) 'mh (black hole mass)    = ',mh
 write(*,*)

 nbins = int(sqrt(real(npart)))
 if (nbins > nmaxbins) nbins = nmaxbins

 call tde_analysis(npart,xyzh,vxyzu,real(luminosity),ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,lbins,dndl)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'e',      &
       2,'dn/de',  &
       3,'tr',     &
       4,'dndt',   &
       5,'r',      &
       6,'dlumdr', &
       7,'lumcdf', &
       8,'lum',    &
       9,'dndlum'

 do i = 1,nbins
    write(iunit,'(9(es18.10,1X))') &
       ebins(i),   &
       dnde(i),    &
       tbins(i),   &
       dndt(i),    &
       rbins(i),   &
       dlumdr(i),  &
       lumcdf(i),  &
       lbins(i),   &
       dndl(i)
 enddo

end subroutine do_analysis

!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine tde_analysis(npart,xyzh,vxyzu,luminosity,ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,lbins,dndl)
 use part, only:isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:),vxyzu(:,:),luminosity(:)
 real, intent(out), dimension(nmaxbins) :: ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,lbins,dndl
 integer :: i
 real    :: eps(npart),tr(npart),r(npart),v2,trmin

 !
 !-- Compute the specific energy and return time of each particle, store in an array
 !
 tr  = 0.
 eps = 0.
 do i=1,npart
    r(i)   = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       v2     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       eps(i) = v2/2. - mh/r(i)                                  !-- Specific energy
       if (eps(i)<0.) then
          tr(i) = treturn(mh,eps(i))                             !-- Return time, only set if energy is negative
       else
          tr(i) = 0.
       endif
    endif
 enddo

 ! Create a histogram of the enegies, return times, luminosity as a function of radius, and luminosity
 call hist(npart,eps,       ebins,dnde,  minval(eps),       maxval(eps),       nbins)
 trmin = treturn(mh,minval(eps))
 call hist(npart,tr,        tbins,dndt,  trmin,             trmin*100.,        nbins)
 call hist(npart,r,         rbins,dlumdr,0.,                rmax,              nbins,luminosity)
 call hist(npart,luminosity,lbins,dndl,  minval(luminosity),maxval(luminosity),nbins)

 !-- Calculate the luminosity CDF as a function of radius
 do i=1,nbins
    lumcdf(i) = sum(dlumdr(1:i))
 enddo

 !-- Normalise the CDF
 if (.not.lumcdf(nbins)<tiny(0.)) then
    lumcdf = lumcdf/lumcdf(nbins)
 else
    lumcdf = 0.
 endif

end subroutine tde_analysis

!--------------------------------------------------------------------------------------------------------------------

!
!-- Function to calculate return time from energy
!
real function treturn(mass,en)
 use physcon,        only:pi
 real, intent(in) :: mass,en

 treturn = 2.*pi*mass/(2.*abs(en))**1.5

end function treturn

!
!-- General function to compute a histogram
!
subroutine hist(np,xarray,xhist,yhist,xmin,xmax,n_bins,yarray)
 use sortutils, only:indexx
 use io,        only:warning
 integer, intent(in) :: np,n_bins
 real, intent(in)    :: xarray(np),xmin,xmax
 real, intent(in), optional :: yarray(np)
 real, intent(out)   :: xhist(n_bins),yhist(n_bins)
 integer :: indx(np),i,ibin,j,nbinned
 logical :: binned
 real    :: dx,xleft,xright,xi

! Sort the array (return sorted indices)
 call indexx(np,xarray,indx)

! Spacing between bins
 dx = (xmax-xmin)/n_bins

! Initialise x and y arrays of histogram arrays
 yhist  = 0.
 xhist  = 0.

!
!--- Count how many particles in each bin
!

! Start at first bin
 ibin   = 1
! Lower and Upper values of starting bin
 xleft  = xmin
 xright = xmin + dx

 nbinned = 0

! Loop over particles in sorted order
 do i=1,np
    j  = indx(i)                        ! sorted index
    xi = xarray(j)
    binned = .false.
    if (xi<xmin) cycle    ! skip particles below xmin
    if (xi/=xi)  cycle    ! skip particles with NaN

    do while (.not.binned.and.ibin<=n_bins)

       if (xleft<=xi.and.xi<=xright) then ! add to bin

          if (present(yarray)) then ! sum the quantity
             yhist(ibin) = yhist(ibin) + yarray(j)

          else ! just count the number
             yhist(ibin) = yhist(ibin) + 1.

          endif

          nbinned = nbinned + 1
          binned  = .true.

       else ! go to next bin
          ibin    = ibin + 1
          xleft   = xright                ! increase the lower bin value
          xright  = xright  + dx          ! increase the upper bin value

       endif

    enddo

 enddo

 !--- Warn the user if not all the particles were placed in bins.
 if (nbinned/=np) then
    call warning('analysis '//analysistype,'Not all particles were binned! n_not_binned',ival=np-nbinned)
 endif

! Create the x axis array for the histogram, use the bin centre as the bin value
 xhist(1) = xmin+dx/2.
 do i=2,n_bins
    xhist(i) = xhist(i-1) + dx
 enddo

end subroutine hist

!----------------------------------------------------------------
!+
!  Read tde information from .tdeparams file
!+
!----------------------------------------------------------------
subroutine read_tdeparams(filename,mhole,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: mhole
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)

! Read in parameters from the file .tdeparams
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(mhole,'mh',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_tdeparams

end module
