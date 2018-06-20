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
!  Computes the specific energy distribution
!
!  REFERENCES: None
!
!  OWNER: David Liptai
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
 character(len=20), parameter, public :: analysistype = 'tde'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use externalforces, only:mass1
 use sortutils,      only:indexx
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: num,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: particlemass,time
 character(len=120) :: fileout
 integer, parameter :: nmaxbins = 5000
 integer, dimension(nmaxbins) :: counts
 integer, dimension(npart)    :: indx
 real,    dimension(npart)    :: en
 integer :: nbins,i,ibin,j
 real :: de,enup,xyz(3),vxyz(3)

 ! Compute the specific energy of each particles, store in an array
 do i=1,npart
    xyz   = xyzh(1:3,i)
    vxyz  = vxyzu(1:3,i)
    en(i) = dot_product(vxyz,vxyz)/2. - mass1/sqrt(dot_product(xyz,xyz))
 enddo

 ! Sort the array (return sorted indices)
 call indexx(npart, en, indx)

 ! Energy bin width
 nbins = int(sqrt(real(npart)))
 de = (maxval(en)-minval(en))/nbins

 counts = 0
 ibin   = 1
 enup   = minval(en) + de
 do i=1,npart
    counts(ibin) = counts(ibin) + 1
    j = indx(i)
    if (en(j)>enup) then
       enup = enup + de
       ibin = ibin + 1
    endif
 enddo

 write(fileout,"(a)") trim(dumpfile)//'.out'
 print "(a)",' writing to '//trim(fileout)
 open(iunit,file=fileout)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',2(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'spec. en.', &
       2,'n'

 do i = 1,nbins
    write(iunit,'(2(es18.10,1X))') real(i),real(counts(i))
 enddo

end subroutine do_analysis


end module
