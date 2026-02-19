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
 use dump_utils, only:read_array_from_file ! here added
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, save :: ii = 0, particlei = 1
 character(len=50), save :: filename
 real(4) :: luminosity(npart) ! here added
 real :: radprop(11,npart) ! here added
 integer :: ierr ! here added

 if (ii==0) then
    call prompt('select which particle you want the position and velocity for',particlei)
    write(filename,"(a,i5.5,a)") 'particle_',particlei,'.ev'
    filename=trim(filename)
    open(unit=1234,file=filename,status='replace')
    write(1234,"(a1,a26,9a27)") '#','time','x','y','z','vx','vy','vz','dotE','radP'
    ii=1
 else
    open(unit=1234,file=filename,position='append')
 endif

 ! read properties from the file
 call read_array_from_file(123,dumpfile,'luminosity',luminosity(1:npart),ierr)
 call read_array_from_file(123,dumpfile,'radiation pressure',radprop(1,1:npart),ierr)
 !call read_array_from_file(123,dumpfile,'radiation pressure',radprop(1:11,1:npart),ierr)

 print*,'radprop',radprop(8,2)
 if (ierr/=0) then
    write(*,*)
    write(*,'("WARNING: could not read radP from file. It will be set to zero")')
    write(*,*)
    radprop = 0.
 endif

 ! write to file (here added)
 write(1234,"(9e27.17)") time,xyzh(1:3,particlei),vxyzu(1:3,particlei),luminosity(particlei),radprop(8,particlei)

 close(1234)

end subroutine do_analysis

end module analysis
