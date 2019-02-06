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
!  Analysis routine for sphNG jet calculations
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, physcon, setup_params
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'jet'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part, only:npartoftype,iamtype,isdead_or_accreted,nptmass,&
                xyzmh_ptmass,igas,rhoh
 use physcon,      only:gg,solarm,au,pi,years
 use setup_params, only:rhozero
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,nout
 real :: macc,mout,mtot,mtotall,mwithinr,maboverho
 real :: xi,yi,zi,hi,rhoi,r,dr,rhatx,rhaty,rhatz,vxi,vyi,vzi,vr
 real :: vr_cut,rcut,rhocut,vmax,tfreefall,tff,tyrs,vmean
 real(kind=8) :: umass,utime,udist,udens
 character(len=30) :: fileout

!
!--global quantities
!
 print*,'npartoftype = ',npartoftype(:)
 print*,'npart = ',npart
 nout = 0
 mout = 0.
 mtot = 0.
 mtotall = 0.
 mwithinr = 0.
 maboverho = 0.
 udist = 1.d16
 umass = solarm
 utime = sqrt(udist**3/(gg*umass))
 udens = umass/udist**3
 print "(a,es10.3)",' udist = ',udist
 print "(a,es10.3)",' umass = ',umass
 print "(a,es10.3)",' utime = ',utime
 print "(a,es10.3)",' udens = ',udens
 tfreefall = sqrt((3.*pi)/(32.*rhozero))
 tff = time/tfreefall
 tyrs = time/(years/utime)
 print*,' t/tff = ',tff,' t(years) = ',tyrs

 !
 !--define mass, radius and velocity cuts
 !  use a small positive but non-zero value of v_cut
 !  to avoid counting material not really part of the outflow
 !  (here 100 m/s in physical units)
 !
 vr_cut = 1.d4/(udist/utime)
 rcut = 100.*(au/udist)
 rhocut = 1.d-14/udens
 vmax = 0.

 print*,'vr cut = ',vr_cut, ' rcut = ',rcut,' rhocut = ',rhocut
 !
 !--compute mass in outflow and accreted mass
 !  as well as maximum velocity
 !
 vmean = 0.
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    rhoi = rhoh(hi,particlemass)
    if (.not. isdead_or_accreted(hi)) then
       r = sqrt(xi*xi + yi*yi + zi*zi)
       dr = 1./r
       rhatx = xi*dr
       rhaty = yi*dr
       rhatz = zi*dr
       vxi = vxyzu(1,i)
       vyi = vxyzu(2,i)
       vzi = vxyzu(3,i)
       vr  = vxi*rhatx + vyi*rhaty + vzi*rhatz
       if (vr > vr_cut) then
          nout = nout + 1
          mout = mout + particlemass
          vmean = vmean + vr
       endif
       mtot = mtot + particlemass
       !--here we count all material that has passed
       !  through a certain radius
       if (r < rcut .or. vr > vr_cut) then
          mwithinr = mwithinr + particlemass
       endif
       !print*,' rhoi = ',rhoi,rhocut
       if (rhoi > rhocut) then
          maboverho = maboverho + particlemass
       endif
       vmax = max(vmax,vr)
    endif
    mtotall = mtotall + particlemass
 enddo
 if (nptmass >= 1) then
    macc = xyzmh_ptmass(4,1)
 else
    macc = 0.
 endif
 mtot = mtot + macc
 mwithinr = mwithinr + macc
 maboverho = maboverho + macc
 vmax = vmax*(udist/utime)/1.d5
 if (nout > 0) vmean = vmean/real(nout)*(udist/utime)/1.d5

 print*,'n(sinks) = ',nptmass,' mtot = ',mtot,mtotall
 print*,'n(vr>',vr_cut,') = ',nout
 print*,'mass (rho>',rhocut,') = ',maboverho
 print*,'mass (r<',rcut,') = ',mwithinr, ' compared to ',macc+mout
 print*,'vmax,mean = ',vmax,vmean,' km/s'

 fileout= 'jetanalysis-macc.ev'
 if (num==1) then
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',10(a14,1x))") 'time','t/tff','t(yrs)','Macc','vmean'
 else
    open(iunit,file=fileout,position='append')
 endif
 print*,'writing to '//trim(fileout)

 fileout= 'jetanalysis-mout.ev'
 if (num==1) then
    open(iunit+1,file=fileout,status='replace')
    write(iunit+1,"('#',10(a14,1x))") 'time','t/tff','t(yrs)','Mout','vdum'
 else
    open(iunit+1,file=fileout,position='append')
 endif
 print*,'writing to '//trim(fileout)

 fileout= 'jetanalysis-mtot.ev'
 if (num==1) then
    open(iunit+2,file=fileout,status='replace')
    write(iunit+2,"('#',10(a14,1x))") 'time','t/tff','t(yrs)','Mtot','vmax'
 else
    open(iunit+2,file=fileout,position='append')
 endif
 print*,'writing to '//trim(fileout)


 write(*,*) 't = ',time,' macc = ',macc,' Mout = ',mout
 write(iunit,"(10(1pe14.6,1x))") time,tff,tyrs,macc,vmean
 write(iunit+1,"(10(1pe14.6,1x))") time,tff,tyrs,mout,-1.
 write(iunit+2,"(10(1pe14.6,1x))") time,tff,tyrs,macc+mout,vmax
 close(iunit)
 close(iunit+1)
 close(iunit+2)

end subroutine do_analysis

end module
