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
!  Analysis routine to look at one particle
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, io, physcon
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'particle'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use eos,     only:get_spsound
 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=9) :: output,filename
 integer :: i
 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 integer, parameter :: iecc    = 23
 real :: ecc,G,M,E,vel(3),pos(3),rad,Li(3),Limag,term,mu

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'part',numfile
 write(*,'("Output file name is ",A)') output

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',3x))") &
       1,'x', &
       2,'y', &
       3,'z'

 !CHOOSE YOUR PARTICLE
 i = 99*160000

 if (i > npart) print*,'Particle chosen does not exist.'

 G = 1.0
 M = 1.0
 mu = G*M
 pos = xyzh(1:3,i)
 vel = vxyz(1:3,i)
 rad = sqrt(dot_product(pos,pos))

 E = 0.5*dot_product(vel,vel) - mu/rad

 Li(1) = pmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
 Li(2) = pmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
 Li(3) = pmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))

 Limag = sqrt(dot_product(Li,Li))/pmass
 term = 2.*E*Limag**2/(mu**2)
 ecc = sqrt(1. + term)

 write(iunit,'(3(es18.10,1X))') xyzh(1,i),xyzh(2,i),xyzh(3,i)

 write(filename,"(a,i3.3)")"tracing"!,i
 if (numfile==0) then
    open(unit=iprec,file=filename,status="replace")
    write(iecc,'("# rad and ecc for particle ",es18.10)') i
    write(iecc,"('#',3(1x,'[',i2.2,1x,a11,']',3x))") &
         1,'time', &
         2,'rad', &
         3,'|e|'
 else
    open(unit=iecc,file=filename,status="old",position="append")
 endif
 if (rad > tiny(rad)) then
    write(iecc,'(3(es18.10,1X))') time,rad,ecc
 else
    write(iecc,'(3(es18.10,1X))') time,0.0,0.0
 endif
 close(unit=iecc)

 close(iunit)

end subroutine do_analysis

end module analysis

