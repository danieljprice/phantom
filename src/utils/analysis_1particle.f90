!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to output a single particle as a function of time
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: prompting
!
 implicit none
 character(len=20), parameter, public :: analysistype = '1particle'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use prompting, only:prompt
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, save :: ii = 0, particlei = 1
 character(len=50), save :: filename

 if (ii==0) then
    call prompt('select which particle you want the position and velocity for',particlei)
    write(filename,"(a,i5.5,a)") 'particle_',particlei,'.ev'
    filename=trim(filename)
    open(unit=1234,file=filename,status='replace')
    write(1234,"(a1,a26,6a27)") '#','time','x','y','z','vx','vy','vz'
    ii=1
 else
    open(unit=1234,file=filename,position='append')
 endif

 write(1234,"(7e27.17)") time,xyzh(1:3,particlei),vxyzu(1:3,particlei)

 close(1234)

end subroutine do_analysis

end module analysis
