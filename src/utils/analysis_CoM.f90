!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Determines the centre of mass of the system and writes it
!               to file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'centerofmass'
 logical, private :: firstcall = .true.

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass
 use centreofmass, only: get_centreofmass
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real                         :: xpos(3),vpos(3)
 logical                      :: iexist
 character(len=200)           :: fileout
 !
 ! Open file (appendif exists)
 !
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_com.dat'
 inquire(file=fileout,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time', &
          2,'x',    &
          3,'y',    &
          4,'z',    &
          5,'vx',   &
          6,'vy',   &
          7,'vz'
 else
    open(iunit,file=fileout,position='append')
 endif
 !
 ! Call centre of mass subroutine
 !
 call get_centreofmass(xpos,vpos,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 !
 ! Write results to file
 !
 write(iunit,'(7(es18.10,1x))') time,xpos,vpos

 close(iunit)

end subroutine do_analysis

end module analysis
