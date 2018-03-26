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
!  Analysis routine for KITP turbulence comparison
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
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
 character(len=20), parameter, public :: analysistype = 'kitp'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,ninbox
 real :: xminbox,xmaxbox,yminbox,ymaxbox,zminbox,zmaxbox
 real :: boxmass,rmsmach,ekin,xmom,ymom,zmom
 real :: xi,yi,zi,vxi,vyi,vzi,v2i
 character(len=30) :: fileout
 integer, parameter :: iunit = 78
!
!--get integrated values of stuff
!
 xminbox = 0.5
 yminbox = 0.5
 zminbox = 0.5
 xmaxbox = xminbox + 1./128.
 ymaxbox = yminbox + 1./128.
 zmaxbox = zminbox + 1./128.

 ninbox = 0
 ekin = 0.
 xmom = 0.
 ymom = 0.
 zmom = 0.
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    if (xi >= xminbox .and. xi <= xmaxbox .and. &
        yi >= yminbox .and. yi <= ymaxbox .and. &
        zi >= zminbox .and. zi <= zmaxbox) then
       ninbox = ninbox + 1
       vxi = vxyzu(1,i)
       vyi = vxyzu(2,i)
       vzi = vxyzu(3,i)
       v2i = vxi*vxi + vyi*vyi + vzi*vzi
       ekin = ekin + v2i
       xmom = xmom + vxi
       ymom = ymom + vyi
       zmom = zmom + vzi
    endif
 enddo

 boxmass = ninbox*particlemass
 if (ninbox > 0) then
    rmsmach = sqrt(ekin/ninbox)
 else
    rmsmach = 0.
 endif
 ekin = 0.5*particlemass*ekin

 fileout= 'kitpanalysis.ev'
 if (num==1) then
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',10(a14,1x))") 'time','mass','ekin','rmsmach','xmom','ymom','zmom'
 else
    open(iunit,file=fileout,position='append')
 endif
 print*,' writing to '//trim(fileout)
 write(*,*) 't = ',time,' mass = ',boxmass,' ekin = ',ekin,' rmsmach = ',rmsmach,' mom = ',xmom,ymom,zmom
 write(iunit,"(10(1pe14.6,1x))") time,boxmass,ekin,rmsmach,xmom,ymom,zmom
 close(iunit)

end subroutine do_analysis

end module
