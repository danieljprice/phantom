!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: testbin
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: testbin [no arguments]
!
!  DEPENDENCIES: datafiles, prompting, testbinary
!+
!--------------------------------------------------------------------------
program testbin
 use testbinary, only:test_binary
 use prompting,  only:prompt
 use datafiles,  only:find_phantom_datafile
 implicit none
 integer :: j,ierr,ierr1,itex
 real :: a,e,inc,o,w,f,m1,m2
 character(len=120) :: filename

 m1 = 1.8
 m2 = 0.4
 call prompt('enter primary mass',m1)
 call prompt('enter secondary mass',m2)
 itex = 0

 filename = find_phantom_datafile('orbits.txt','orbits')
 open(unit=2,file=filename,status='old',iostat=ierr)
 if (ierr==0) then
    open(newunit=itex,file='orbits.tex',status='replace',iostat=ierr1)
    print "(a)",' writing to orbits.tex'
 endif
 j = 0
 do while(ierr==0)
    j = j + 1
    read(2,*,iostat=ierr) a,e,inc,o,w,f
    if (ierr==0) call test_binary(m1,m2,a,e,inc,o,w,f,j,itex)
 enddo
 close(2)
 if (itex > 0) close(itex)

end program testbin
