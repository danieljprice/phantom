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
!  Analysis routine for comparing mass accreted onto sinks
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
 character(len=20), parameter, public :: analysistype = 'multiple_sinks'
 public                               :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,    only:xyzmh_ptmass,nptmass
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,iunit,npart
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: time,particlemass
 character(len=30)            :: dumpprefix,fileout
 integer                      :: i,iline,unitnum,blah

 iline = index(dumpfile,'_')
 dumpprefix = dumpfile(1:iline)
 unitnum = 1001
 do i=1,nptmass
    write(fileout,"(2a,i3.3,a)") 'sink_',trim(dumpprefix),i,'.ev'
    if (num==0) then
       open(unit=unitnum, file=fileout, status='replace')
       write(unitnum,"('#',1x,9('[',i1.1,a17,']',2x))") &
           1,'time',           &
           2,'x co-ordinate',  &
           3,'y co-ordinate',  &
           4,'z co-ordinate',  &
           5,'mass',           &
           6,'mass accreted',  &
           7,'spin ang mom x', &
           8,'spin ang mom y', &
           9,'spin ang mom z'
       write(unitnum,"(9(3x,es18.11e2,1x))") &
           time,                &
           xyzmh_ptmass(1:4,i), &
           xyzmh_ptmass(7:10,i)
       close(unit=unitnum)
       unitnum = unitnum + 1
    else
       open(unit=unitnum, file=fileout, position='append')
       write(unitnum,"(9(3x,es18.11e2,1x))") &
           time,                &
           xyzmh_ptmass(1:4,i), &
           xyzmh_ptmass(7:10,i)
       close(unit=unitnum)
       unitnum = unitnum + 1
    endif
 enddo

end subroutine do_analysis

end module
