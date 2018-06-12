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
!  Analysis routine computing particle tracks
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'tracks'
 public  :: do_analysis

 logical         :: initialised = .false.
 logical         :: list        = .false.

 integer, parameter                     :: maxnr      = 20
 integer, parameter                     :: maxnzoverr = 10
 integer, parameter                     :: maxtrack   = 200
 character(len=17) :: output(maxtrack)
 integer :: itracked(maxtrack)
 real :: radtrack(maxnr)
 real :: zoverrtrack(maxnzoverr)
 integer :: flag(maxnr,maxnzoverr)
 integer                                :: ntracked,dimtracked,nrtracked,nzoverrtracked,nflag
 integer                                :: nbpart
 real, parameter :: rtol      = 1.d-1
 real, parameter :: zoverrtol = 1.d-3
 real            :: rmin,rmax,zmin,zmax

 real :: zstart(maxtrack),vzstart(maxtrack),szstart(maxtrack)
 real, parameter :: grainsize = 1.d-12
 real, parameter :: densgrain = 1.
 real, parameter :: p      = 0.75
 real, parameter :: rin    = 0.01
 real, parameter :: q      = 0.
 real, parameter :: csz    = 0.01
 real, parameter :: sigmaz = 0.01
 real            :: tk,rhog,cs,sz,sigma,lambda,lambdap,lambdam,zpred,zini,vzini
 logical         :: ipredz = .false.
 logical         :: startflag

 save zstart,vzstart,szstart,startflag
 data startflag/.true./

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part, only:iphase
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 real :: costh,sinth,rad,vr,vth
 integer :: i,j,ierr

! --- read the tracked particles in a list or track particular particles---
! --- find a number ntracked of particles and open a related output file---

 if (.not.initialised) then
    if (list) then
       open(unit=223,file='tracklist',status='old',iostat=ierr)
       if (ierr /= 0) stop 'cannot find a file called tracklist'
       ntracked = 0
       do while (ierr==0 .and. ntracked < maxtrack)
          ntracked = ntracked + 1
          read(223,*,iostat=ierr) itracked(ntracked)
          if (ierr /= 0 .or. itracked(ntracked) <= 0 .or. itracked(ntracked) > size(xyzh(1,:))) then
             ntracked = ntracked -1
             if (ierr==0) print *,'Warning: Index of tracked particle exceeds the array dimensions'
          endif
       enddo
       close(223)

    else
       open(unit=224,file='rlist',status='old',iostat=ierr)
       if (ierr /= 0) stop 'cannot find a file called rlist'
       nrtracked = 0
       do while (ierr==0 .and. nrtracked < maxnr)
          nrtracked = nrtracked + 1
          read(224,*,iostat=ierr) radtrack(nrtracked)
          if (ierr /= 0) then
             nrtracked = nrtracked -1
             if (ierr==0) print *,'Warning: Index of tracked particle exceeds the array dimensions'
          endif
       enddo
       close(224)

       open(unit=225,file='zoverrlist',status='old',iostat=ierr)
       if (ierr /= 0) stop 'cannot find a file called zoverrlist'
       nzoverrtracked = 0
       do while (ierr==0 .and. nzoverrtracked < maxnzoverr)
          nzoverrtracked = nzoverrtracked + 1
          read(225,*,iostat=ierr) zoverrtrack(nzoverrtracked)
          if (ierr /= 0) then
             nzoverrtracked = nzoverrtracked -1
             if (ierr==0) print *,'Warning: Index of tracked particle exceeds the array dimensions'
          endif
       enddo
       close(225)

       dimtracked = nrtracked*nzoverrtracked

       print*,'number of radial bins',nrtracked,' number of vertical bins',nzoverrtracked

       do i = 1,nrtracked
          do j = 1,nzoverrtracked
             flag(i,j) = 0
          enddo
       enddo
       ntracked  = 0
       nbpart    = 1

       !--iphase(i)==5 corresponds to a dust particle : only the dust particles are tracked
       do while ( (nbpart < npart) .and. (ntracked /= dimtracked) )
          rad = sqrt(xyzh(1,nbpart)**2+xyzh(2,nbpart)**2)
          do i = 1,nrtracked
             do j = 1,nzoverrtracked
                rmin = radtrack(i) - rtol
                rmax = radtrack(i) + rtol
                zmin = radtrack(i)*(zoverrtrack(j) - zoverrtol)
                zmax = radtrack(i)*(zoverrtrack(j) + zoverrtol)

                if ( (iphase(nbpart)==5).and.(flag(i,j) /= 1).and.(rad >= rmin).and.(rad <= rmax)&
                  .and.(xyzh(3,nbpart) >= zmin).and.(xyzh(3,nbpart) <= zmax)) then

                   flag(i,j) = 1
                   ntracked  = ntracked + 1
                   print*,'nbpart',nbpart,rad,xyzh(3,nbpart)
                   itracked(ntracked) = nbpart

                endif
             enddo
          enddo
          nbpart = nbpart + 1
       enddo

       if (ntracked /= dimtracked) then
          print*,'Warning: some desired particles have not been found: increase the tolerance'
       endif
    endif

    print*,'ntracked',ntracked

    do i = 1,ntracked
       write(output(i),"(a,i8.8,a)") 'track',itracked(i),'.out'
       open(iunit+i-1,file=output(i),status='replace')
       print *,'Writing to file... ',output(i)
    enddo
    initialised = .true.

 else
    do i = 1,ntracked
       open(iunit+i-1,file=output(i),status='old',position='append')
    enddo
 endif
!--------------------------------------------------------------------------

 do i=1,ntracked
    j     = itracked(i)
    rad   = sqrt(xyzh(1,j)**2+xyzh(2,j)**2)
    costh = xyzh(1,j)/rad
    sinth = xyzh(2,j)/rad
    vr    = costh*vxyzu(1,j) + sinth*vxyzu(2,j)
    vth   = sinth*vxyzu(1,j) - costh*vxyzu(2,j)

!--if ipredz = .true., calculate the harmonic oscillator approximation for the solution------
!--DOES NOT WORK YET :)
    if (ipredz) then
       if (startflag) then
          zstart(i)     = xyzh(3,j)
          vzstart(i)    = vxyzu(3,j)
          sigma         = sigmaz*rad**(-p)*(1.0d0 - sqrt(rin/rad))
          szstart(i)    = grainsize*densgrain*2.50662827/(sigma*csz) !2.50662827 = sqrt(2*Pi)
          startflag     = .false.
       endif

       zini    = zstart(i)
       vzini   = vzstart(i)
       sz      = szstart(i)
       tk    = rad**1.5 !Keplerian time in code units
       if (sz > 0.5) then
          lambda  = sqrt(1. - (1./(2*sz))**2)
          zpred   = exp(-time/(2*tk*sz))*(zini*cos(lambda*time/tk) &
                       + 1./lambda*(1./(2.*sz)*zini + vzini*tk)*sin(lambda*time/tk))
       elseif (sz==0.5) then
          zpred   = (zini + vzini*time)*exp(-time/tk)
       else
          lambdap = -1./(2.*sz) + sqrt((1./(2*sz))**2 - 1.)
          lambdam = -1./(2.*sz) - sqrt((1./(2*sz))**2 - 1.)
          zpred   = (zini + (lambdap*zini - vzini*tk)/(lambdam - lambdap))*exp(lambdam*time/tk) + &
                       (-(lambdap*vzini - vzini*tk)/(lambdam - lambdap))*exp(lambdap*time/tk)
       endif
    endif

!--write the output files------------------------------------------------------------------------------
    if (ipredz) then
       write(iunit+i-1,"(20(es16.8,1x))") time,xyzh(1,j),xyzh(2,j),xyzh(3,j), &
                                          zpred,rad, &
                                          vxyzu(1,j),vxyzu(2,j),vxyzu(3,j),vr,vth, &
                                          sigma,sz,tk
    else
       write(iunit+i-1,"(20(es16.8,1x))") time,xyzh(1,j),xyzh(2,j),xyzh(3,j),rad, &
                                          vxyzu(1,j),vxyzu(2,j),vxyzu(3,j),vr,vth
    endif

    close(iunit+i-1)
 enddo

end subroutine do_analysis

end module
