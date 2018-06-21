!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
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
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'tde'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 character(len=120) :: output
 integer, parameter :: nmaxbins = 5000
 real,    dimension(nmaxbins) :: ebins,dnde,tbins,dndt
 real,    dimension(npart)    :: eps,tr
 integer :: nbins,i
 real    :: mass1,r,v2,trmin

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'hist',numfile
 write(*,'("Output file name is ",A)') output

 ! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G == 1 and mass1 == 1e6")')

 mass1 = 1.e6
 nbins = int(sqrt(real(npart)))

!
!-- Compute the specific energy and return time of each particles, store in an array
!
 tr  = 0.
 eps = 0.
 do i=1,npart
    r      = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    v2     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    eps(i) = v2/2. - mass1/r                                     !-- Specific energy
    if (eps(i)<0.) then
       tr(i) = treturn(mass1,eps(i))                !-- Return time, only set if energy is negative
    else
       tr(i) = 0.
    endif
 enddo

! Create a histogram of the enegies and return times
 call hist(npart,eps,ebins,dnde,minval(eps),maxval(eps),nbins)
 trmin = treturn(mass1,minval(eps))
 call hist(npart,tr,tbins,dndt,trmin,trmin*100.,nbins)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',4(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'e',    &
       2,'dn/de',&
       3,'tr',   &
       4,'dndt'

 do i = 1,nbins
    write(iunit,'(4(es18.10,1X))') ebins(i),dnde(i),tbins(i),dndt(i)
 enddo

end subroutine do_analysis

real function treturn(mass,en)
 use physcon,        only:pi
 real, intent(in) :: mass,en

  treturn = 2.*pi*mass/(2.*abs(en))**1.5

end function treturn

subroutine hist(np,xarray,xhist,yhist,xmin,xmax,nbins)
 use sortutils, only:indexx
 integer, intent(in) :: np,nbins
 real, intent(in)    :: xarray(np),xmin,xmax
 real, intent(out)   :: xhist(nbins),yhist(nbins)
 integer :: indx(np),i,ibin,j
 real    :: dx,xright

! Sort the array (return sorted indices)
 call indexx(np,xarray,indx)

! Spacing between bins
 dx = (xmax-xmin)/nbins

! Initialise x and y arrays of histogram arrays
 yhist  = 0.
 xhist  = 0.

!
!--- Count how many particles in each bin
!

! Start at first bin
 ibin   = 1
! Upper value of starting bin
 xright = xmin + dx

! Loop over particles in sorted order
 do i=1,np
    j = indx(i)                        ! sorted index
    if (xarray(j)<xright) then
       yhist(ibin) = yhist(ibin) + 1.  ! add to bin
    else
       ibin    = ibin + 1              ! go to next bin
       xright  = xright  + dx          ! increase the upper bin upper value
    endif
 enddo

! Create the x axis array for the histogram, use the bin centre as the bin value
 xhist(1) = xmin+dx/2.
 do i=2,nbins
    xhist(i) = xhist(i-1) + dx
 enddo

end subroutine hist

end module
