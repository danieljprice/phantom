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
!  Analysis routine to calculate the mass weighted density PDF
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, part, pdfs, readwrite_dumps
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'MWpdf'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,              only:maxptmass, maxp
 use part,             only:massoftype,rhoh,hfact,isdead_or_accreted
 use readwrite_dumps,  only:read_dump
 use pdfs,             only:get_pdf,write_pdf
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num, npart, iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:), particlemass, time
 character(len=3)  :: variable='rho'
 integer           :: nbins, i
 real, allocatable :: pdf_rho(:), xbin(:)
 real :: totmass,rhomin,rhomax,binspacing,rhologmin,rhologmax
 real :: vx2,vy2,vz2,rhoi,rmsv
 real :: rhomean, rho(maxp), hi, pmassi

 rhomin = huge(rhomin)
 rhomax = 0.
 rhomean = 0.
 totmass = 0.
 print*,'hfact = ',hfact
 do i=1,npart
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       pmassi = massoftype(1)
       rho(i) = rhoh(hi,pmassi)
       rhoi   = rho(i)
       rhomin = min(rhomin,rhoi)
       rhomax = max(rhomax,rhoi)
       vx2 = vxyzu(1,i)**2
       vy2 = vxyzu(2,i)**2
       vz2 = vxyzu(3,i)**2
       rmsv = rmsv + pmassi*(vx2 + vy2 + vz2)
       totmass = totmass + pmassi
       rhomean = rhomean + pmassi*rhoi
    endif
 enddo

 rmsv = sqrt(rmsv/totmass)
 rhomean = rhomean/totmass
 print*,'On parts:      rms v = ',rmsv
 print*,'On parts: Total Mass = ',totmass,' mean dens = ',rhomean
 print*,'On parts:  max. dens = ',rhomax, ' min. dens = ',rhomin

! rhologmin = log10(rhomin)
! rhologmax = log10(rhomax)

 rhologmin = -10.
 rhologmax = 10.

!
!--check to see if there are densities smaller or larger than the limits of the PDF
!
 if (log10(rhomin)<rhologmin) then
    print*, 'ERROR: you have densities lower than the minimum value of 10^-10'
    print*, 'Bailing out.'
    stop
 endif
 if (log10(rhomax)>rhologmax) then
    print*, 'ERROR: you have densities higher than the maximum value of 10^10'
    print*, 'Bailing out.'
    stop
 endif

!
!--allocate memory for PDF calculation
!
 binspacing = 0.01
 nbins = nint((rhologmax - rhologmin)/binspacing)

 if (.not.allocated(pdf_rho)) allocate(pdf_rho(nbins))
 if (.not.allocated(xbin)) allocate(xbin(nbins))

 !
 !--calculate PDF of rho, and write to text file
 !
! call get_pdf(xbin,pdf,nbins,npart,lnrho,lnrhomin,lnrhomax)
! call get_pdf(xbin,pdf,nbins,npart,rho,rhomin,rhomax)
 call get_pdf(xbin,pdf_rho,nbins,npart,rho,10.**int(rhologmin),10.**int(rhologmax))
 call write_pdf(iunit,dumpfile,variable,nbins,xbin,pdf_rho,time)

 if (allocated(xbin)) deallocate(xbin)
 if (allocated(pdf_rho)) deallocate(pdf_rho)
!---------------------------------------------

end subroutine do_analysis

end module
