!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION: Traces the centre of (stellar) mass each galaxy that
!  that is part of a major merger.  Stars initially assigned to each
!  galaxy are assumed to remain in that galaxy for all time.
!
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, physcon, units
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'GalMerger'
 logical, private :: firstcall = .true.

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,iphase,istar,massoftype,igas,&
                         maxp,maxphase,isdead_or_accreted,iamtype
 use units,        only: utime,udist
 use physcon,      only: years,mpc
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i,itype
 real                         :: xi,yi,zi,hi,pmassi,time0
 real                         :: xpos1,ypos1,zpos1,xpos2,ypos2,zpos2,totmass1,totmass2,r
 logical                      :: iexist
 character(len=200)           :: fileout
 !
 !--Open file (appendif exists)
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_stellarCoM.dat'
 inquire(file=fileout,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time', &
          2,'x1',   &
          3,'y1',   &
          4,'z1',   &
          5,'x2',   &
          6,'y2',   &
          7,'z2',   &
          8,'r'
 else
    open(iunit,file=fileout,position='append')
 endif
 !
 xpos1    = 0.d0
 ypos1    = 0.d0
 zpos1    = 0.d0
 totmass1 = 0.d0
 xpos2    = 0.d0
 ypos2    = 0.d0
 zpos2    = 0.d0
 totmass2 = 0.d0
!$omp parallel default(none) &
!$omp shared(npart,xyzh,iphase,massoftype) &
!$omp private(i,itype,xi,yi,zi,hi,pmassi) &
!$omp reduction(+:xpos1,ypos1,zpos1,xpos2,ypos2,zpos2,totmass1,totmass2)
!$omp do
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (itype > 0) then ! avoid problems if called from ICs
             pmassi = massoftype(itype)
          else
             pmassi = massoftype(igas)
          endif
       else
          pmassi = massoftype(igas)
       endif
       if (itype==istar) then
          if (i <= npart/2 ) then
             totmass1 = totmass1 + pmassi
             xpos1    = xpos1    + pmassi*xi
             ypos1    = ypos1    + pmassi*yi
             zpos1    = zpos1    + pmassi*zi
          else
             totmass2 = totmass2 + pmassi
             xpos2    = xpos2    + pmassi*xi
             ypos2    = ypos2    + pmassi*yi
             zpos2    = zpos2    + pmassi*zi
          endif
       endif
    endif
 enddo
!$omp enddo
!$omp end parallel
 !
 xpos1  = xpos1/totmass1 * udist/mpc + 0.5
 ypos1  = ypos1/totmass1 * udist/mpc + 0.5
 zpos1  = zpos1/totmass1 * udist/mpc + 0.5
 xpos2  = xpos2/totmass2 * udist/mpc + 0.5
 ypos2  = ypos2/totmass2 * udist/mpc + 0.5
 zpos2  = zpos2/totmass2 * udist/mpc + 0.5
 r      = sqrt( (xpos1-xpos2)**2 + (ypos1-ypos2)**2 + (zpos1-zpos2)**2 )
 time0  = time*utime/years*1.0d-10
 !
 write(iunit,'(8(es18.10,1x))') time0,xpos1,ypos1,zpos1,xpos2,ypos2,zpos2,r
 close(iunit)
 !
end subroutine do_analysis

end module

