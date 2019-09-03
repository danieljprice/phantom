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

 real    :: mh   = 1.
 real    :: rmax = 700.

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,         only:warning
 use dump_utils, only:read_array_from_file
 use prompting,  only:prompt
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 character(len=120) :: output
 character(len=20)  :: filename
 real, dimension(:), allocatable :: ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,lbins,dndl,angbins,dndang,vbins,dndv
 integer :: i,ierr
 logical :: iexist
 real(4) :: luminosity(npart)
 logical, save :: first = .true.

 call read_array_from_file(123,dumpfile,'luminosity',luminosity(1:npart),ierr)
 if (ierr/=0) then
    write(*,*)
    write(*,'("WARNING: could not read luminosity from file. It will be set to zero")')
    write(*,*)
    luminosity = 0.
 endif

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'hist',numfile
 write(*,'("Output file name is ",A)') output

 ! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G == 1")')

 ! Read black hole mass from params file
 filename = 'analysis_'//trim(analysistype)//'.params'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_tdeparams(filename,ierr)
 if (.not.iexist.or.ierr/=0) then
    nbins = int(sqrt(real(npart)))
    call write_tdeparams(filename)
    print*,' Edit '//trim(filename)//' and rerun phantomanalysis'
    stop
 endif

 allocate(ebins(nbins),dnde(nbins),tbins(nbins),dndt(nbins),rbins(nbins),dlumdr(nbins),lumcdf(nbins),&
          lbins(nbins),dndl(nbins),angbins(nbins),dndang(nbins),vbins(nbins),dndv(nbins))

 if (first) then
   ! Print out the parameters
    write(*,*)
    write(*,'("Parameters are:")')
    write(*,*) 'mh (black hole mass)    = ',mh
    write(*,*)
    first = .false.
 endif

 call tde_analysis(npart,xyzh,vxyzu,real(luminosity),ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,&
                   lbins,dndl,angbins,dndang,vbins,dndv)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'e',      &
       2,'dn/de',  &
       3,'tr',     &
       4,'dndt',   &
       5,'r',      &
       6,'dlumdr', &
       7,'lumcdf', &
       8,'lum',    &
       9,'dndlum', &
       10,'angm',  &
       11,'dndang',&
       12,'v'     ,&
       13,'dndv'

 do i = 1,nbins
    write(iunit,'(13(es18.10,1X))') &
       ebins(i),   &
       dnde(i),    &
       tbins(i),   &
       dndt(i),    &
       rbins(i),   &
       dlumdr(i),  &
       lumcdf(i),  &
       lbins(i),   &
       dndl(i),    &
       angbins(i), &
       dndang(i),  &
       vbins(i),   &
       dndv(i)
 enddo

 deallocate(ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,lbins,dndl,angbins,dndang,vbins,dndv)

end subroutine do_analysis

!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine tde_analysis(npart,xyzh,vxyzu,luminosity,ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,&
                        lbins,dndl,angbins,dndang,vbins,dndv)
 use part,        only:isdead_or_accreted
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:),vxyzu(:,:),luminosity(:)
 real, intent(out), dimension(:) :: ebins,dnde,tbins,dndt,rbins,dlumdr,lumcdf,lbins,dndl,angbins,dndang,vbins,dndv
 integer :: i
 real, dimension(npart) :: eps,tr,r,Langm,vel
 real :: v2,trmin,Li(3)

 !
 !-- Compute the specific energy and return time of each particle, store in an array
 !
 tr  = 0.
 eps = 0.
 do i=1,npart
    r(i)   = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    v2     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Li)
    Langm(i) = sqrt(dot_product(Li,Li))
    eps(i) = v2/2. - mh/r(i)                                  !-- Specific energy
    if (eps(i)<0.) then
       tr(i) = treturn(mh,eps(i))                             !-- Return time, only set if energy is negative
    else
       tr(i) = 0.
    endif
    vel(i) = sqrt(v2)
 enddo

 ! Create a histogram of the enegies, return times, luminosity as a function of radius, and luminosity
 call hist(npart,eps,       ebins,dnde,  minval(eps),       maxval(eps),       nbins)
 trmin = treturn(mh,minval(eps))
 call hist(npart,tr,        tbins,dndt,  trmin,             trmin*100.,        nbins)
 call hist(npart,r,         rbins,dlumdr,0.,                rmax,              nbins,luminosity)
 call hist(npart,luminosity,lbins,dndl,  minval(luminosity),maxval(luminosity),nbins)
 call hist(npart,Langm,     angbins,dndang,minval(Langm),maxval(Langm),nbins)
 call hist(npart,vel,       vbins,dndv,minval(vel),maxval(vel),nbins)

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
!  Read/write tde information from/to params file
!+
!----------------------------------------------------------------
subroutine write_tdeparams(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing analysis options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# options when performing TDE analysis'
 call write_inopt(mh,'mh','black hole mass in code units',iunit)
 call write_inopt(nbins,'nbins','number of bins',iunit)
  call write_inopt(rmax,'rmax','',iunit)
 close(iunit)

end subroutine write_tdeparams

subroutine read_tdeparams(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading analysis options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mh,'mh',db,min=0.,errcount=nerr)
 call read_inopt(nbins,'nbins',db,min=0,errcount=nerr)
 call read_inopt(rmax,'rmax',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of params file: re-writing...'
    ierr = nerr
 endif

end subroutine read_tdeparams

end module
