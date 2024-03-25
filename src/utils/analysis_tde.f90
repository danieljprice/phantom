!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Computes the the distribution of different quantities in a TDE simulation
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - angmax : *max angular momentum*
!   - angmin : *min angular momentum*
!   - emax   : *max energy*
!   - emin   : *min energy*
!   - lummax : *max luminosity*
!   - lummin : *min luminosity*
!   - mh     : *black hole mass in code units*
!   - nbins  : *number of bins*
!   - rmax   : *max radius*
!   - rmin   : *min radius*
!   - trmax  : *max return time*
!   - trmin  : *min return time*
!   - vmax   : *max velocity*
!   - vmin   : *min velocity*
!
! :Dependencies: dump_utils, infile_utils, io, physcon, prompting,
!   readwrite_dumps, sortutils, vectorutils
!
 implicit none
 character(len=3), parameter, public :: analysistype = 'tde'
 public :: do_analysis

 private

 real, dimension(:), allocatable :: ebins,dmde,tbins,dmdt,rbins,dlumdr,lumcdf,lbins,dmdl,angbins,dmdang,vbins,dmdv

 !---- These can be changed in the params file
 integer :: nbins
 real :: mh     = 1.
 !--- If min=max then lims are set dynamically
 real :: emin   = 0.
 real :: emax   = 0.
 real :: trmin  = 0.
 real :: trmax  = 0.
 real :: rmin   = 0.
 real :: rmax   = 0.
 real :: lummin = 0.
 real :: lummax = 0.
 real :: angmin = 0.
 real :: angmax = 0.
 real :: vmin   = 0.
 real :: vmax   = 0.

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,         only:warning
 use dump_utils, only:read_array_from_file
 use prompting,  only:prompt
 use readwrite_dumps, only: opened_full_dump
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 character(len=120) :: output
 character(len=20)  :: filename
 integer :: i,ierr
 logical :: iexist
 real(4) :: luminosity(npart)


 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

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
 print*, trim(analysistype), 'trim(analysistype)'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_tdeparams(filename,ierr)
 if (.not.iexist.or.ierr/=0) then
    nbins = int(sqrt(real(npart)))
    call write_tdeparams(filename)
    print*,' Edit '//trim(filename)//' and rerun phantomanalysis'
    stop
 endif

 if (.not.allocated(ebins)) then
    allocate(ebins(nbins),dmde(nbins),tbins(nbins),dmdt(nbins),rbins(nbins),dlumdr(nbins),lumcdf(nbins),&
             lbins(nbins),dmdl(nbins),angbins(nbins),dmdang(nbins),vbins(nbins),dmdv(nbins))
 endif

 call tde_analysis(npart,pmass,xyzh,vxyzu,real(luminosity))

 open(iunit,file=output)
 write(iunit,'("# ",es20.12,"   # TIME")') time
 write(iunit,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'e',      &
       2,'dmde',  &
       3,'tr',     &
       4,'dmdt',   &
       5,'r',      &
       6,'dlumdr', &
       7,'lumcdf', &
       8,'lum',    &
       9,'dmdlum', &
       10,'angm',  &
       11,'dmdang',&
       12,'v'     ,&
       13,'dmdv'

 do i = 1,nbins
    write(iunit,'(13(es18.10,1X))') &
       ebins(i),   &
       dmde(i),    &
       tbins(i),   &
       dmdt(i),    &
       rbins(i),   &
       dlumdr(i),  &
       lumcdf(i),  &
       lbins(i),   &
       dmdl(i),    &
       angbins(i), &
       dmdang(i),  &
       vbins(i),   &
       dmdv(i)
 enddo


end subroutine do_analysis

!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine tde_analysis(npart,pmass,xyzh,vxyzu,luminosity)
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
 real, intent(in)    :: pmass,xyzh(:,:),vxyzu(:,:),luminosity(:)
 integer :: i
 real, dimension(npart) :: eps,tr,r,Langm,vel,mass
 real :: v2,Li(3),de,dt,dr,dl,dang,dv

 !
 !-- Compute the specific energy and return time of each particle, store in an array
 !
 tr  = 0.
 eps = 0.
 do i=1,npart
    r(i)   = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    v2     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Li) !Should not multiply by particle mass??

    Langm(i) = sqrt(dot_product(Li,Li))
    eps(i) = v2/2. - mh/r(i)                                  !-- Specific energy
    if (eps(i)<0.) then
       tr(i) = treturn(mh,eps(i))                             !-- Return time, only set if energy is negative
    else
       tr(i) = 0.
    endif
    vel(i) = sqrt(v2)
    mass(i) = pmass
 enddo

 !---- If min=max then do dynamic limits for each dump
 if (abs(emin-emax)<tiny(1.)) then
    emin   = minval(eps)
    emax   = maxval(eps)
 endif
 if (abs(trmin-trmax)<tiny(1.)) then
    trmin  = treturn(mh,minval(eps))
    trmax  = 100.*trmin
 endif
 if (abs(rmin-rmax)<tiny(1.)) then
    rmin   = minval(r)
    rmax   = maxval(r)
 endif
 if (abs(lummin-lummax)<tiny(1.)) then
    lummin = minval(luminosity)
    lummax = maxval(luminosity)
 endif
 if (abs(angmin-angmax)<tiny(1.)) then
    angmin = minval(Langm)
    angmax = maxval(Langm)
 endif
 if (abs(vmin-vmax)<tiny(1.)) then
    vmin   = minval(vel)
    vmax   = maxval(vel)
 endif

 ! Create a histogram of the enegies, return times, luminosity as a function of radius, and luminosity
 ! --- Note: The arrays returned by hist are NOT normalised.
 !        e.g. the dmde returned is actually just mass per bin, and is not mass per energy
 !        Rescaling to have correct dimensions is done further down.
 !        This is useful if comparing distributions that have different bin widths.
 !
 !        Counting the mass in each bin is also more useful than just simply the number of particles
 !        per bin, since it allows for comparison at different resolutions (i.e. different particle mass)
 call hist(npart, eps       , ebins  , dmde  , emin  , emax  , nbins, mass)
 call hist(npart, tr        , tbins  , dmdt  , trmin , trmax , nbins, mass)
 call hist(npart, r         , rbins  , dlumdr, rmin  , rmax  , nbins, luminosity)
 call hist(npart, luminosity, lbins  , dmdl  , lummin, lummax, nbins, mass)
 call hist(npart, Langm     , angbins, dmdang, angmin, angmax, nbins, mass)
 call hist(npart, vel       , vbins  , dmdv  , vmin  , vmax  , nbins, mass)

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

 de   = (emax - emin)/nbins
 dt   = (trmax - trmin)/nbins
 dr   = (rmax - rmin)/nbins
 dl   = (lummax - lummin)/nbins
 dang = (angmax - angmin)/nbins
 dv   = (vmax - vmin)/nbins

 !--- Rescale to have correct dimensions
 dmde   = dmde/de
 dmdt   = dmdt/dt
 dlumdr = dlumdr/dr
 dmdang = dmdang/dang
 dmdv   = dmdv/dv

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
subroutine hist(np,xarray,xhist,yhist,xmin,xmax,n_bins,weights)
 use sortutils, only:indexx
 use io,        only:warning
 integer, intent(in) :: np,n_bins
 real, intent(in)    :: xarray(np),xmin,xmax
 real, intent(in), optional :: weights(np)
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

          if (present(weights)) then ! sum the quantity
             yhist(ibin) = yhist(ibin) + weights(j)

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
 write(iunit,"(a,/)") '# options when performing TDE analysis'
 call write_inopt(mh,'mh','black hole mass in code units',iunit)
 call write_inopt(nbins,'nbins','number of bins',iunit)

 write(iunit,"(/,a)") '# limits for binning (set min=max to have dynamic limits)'
 call write_inopt(emin,'emin','min energy',iunit)
 call write_inopt(emax,'emax','max energy',iunit)

 call write_inopt(trmin,'trmin','min return time',iunit)
 call write_inopt(trmax,'trmax','max return time',iunit)

 call write_inopt(rmin,'rmin','min radius',iunit)
 call write_inopt(rmax,'rmax','max radius',iunit)

 call write_inopt(lummin,'lummin','min luminosity',iunit)
 call write_inopt(lummax,'lummax','max luminosity',iunit)

 call write_inopt(angmin,'angmin','min angular momentum',iunit)
 call write_inopt(angmax,'angmax','max angular momentum',iunit)

 call write_inopt(vmin,'vmin','min velocity',iunit)
 call write_inopt(vmax,'vmax','max velocity',iunit)

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
 call read_inopt(mh,   'mh',   db,min=0.,errcount=nerr)
 call read_inopt(nbins,'nbins',db,min=0, errcount=nerr)

 call read_inopt(emin,'emin',db,         errcount=nerr)
 call read_inopt(emax,'emax',db,min=emin,errcount=nerr)

 call read_inopt(trmin,'trmin',db,min=0.,   errcount=nerr)
 call read_inopt(trmax,'trmax',db,min=trmin,errcount=nerr)

 call read_inopt(rmin,'rmin',db,min=0.,  errcount=nerr)
 call read_inopt(rmax,'rmax',db,min=rmin,errcount=nerr)

 call read_inopt(lummin,'lummin',db,min=0.,    errcount=nerr)
 call read_inopt(lummax,'lummax',db,min=lummin,errcount=nerr)

 call read_inopt(angmin,'angmin',db,min=0.,    errcount=nerr)
 call read_inopt(angmax,'angmax',db,min=angmin,errcount=nerr)

 call read_inopt(vmin,'vmin',db,min=0.,  errcount=nerr)
 call read_inopt(vmax,'vmax',db,min=vmin,errcount=nerr)

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of params file: re-writing...'
    ierr = nerr
 endif

end subroutine read_tdeparams

end module analysis
