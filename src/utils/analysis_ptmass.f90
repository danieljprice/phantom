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
!  Analysis routine for sinks: assumes central sink with label 1,
!                              and all sinks orbit this
!                              assumes G=1.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, io, options, part, physcon, setbinary
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'ptmass'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

 use io,       only:fatal
 use part,     only:xyzmh_ptmass,vxyz_ptmass,ihacc,nptmass
 use physcon,  only:pi
 use setbinary,only:Rochelobe_estimate
 use eos,      only:gamma
 use options,  only:ieos

 character(len=*), intent(in) :: dumpfile
 real,dimension(:,:), intent(in) :: xyzh,vxyzu
 real, intent(in) :: pmass,time
 integer, intent(in) :: npart,iunit,numfile

 integer :: i,nsinks
 real :: G

! Use two variables solely to remove compiler warnings...
 write(*,'("Performing analysis on ",a,"... which is unit ",i5,"...")') trim(dumpfile),iunit

! Check nsinks is >=2
 if (nptmass < 2) call fatal(analysistype,'Not enough sinks...')

! Assuming G=1.
 G=1.0

! Currently assuming all sinks orbit sink1
 do i = 2,nptmass
    call get_binary_params(1,i,xyzmh_ptmass,vxyz_ptmass,time,G)
 enddo

end subroutine do_analysis
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
! Calculate the binary paramters
!-----------------------------------------------------------------------
subroutine get_binary_params(ipri,isec,xyzmh_ptmass,vxyz_ptmass,time,G)
!-----------------------------------------------------------------------
 use io, only:fatal

 implicit none

 integer,intent(in) :: ipri,isec
 real,intent(in) :: time,G
 real,dimension(:,:),intent(in) :: xyzmh_ptmass,vxyz_ptmass

 logical :: exists
 character(len=25) :: output
 integer,parameter :: iunit = 150
 integer :: check
 real :: rbin,mpri,msec,E,Lmag,a,ecc,inc
 real,dimension(3) :: xpri,vpri,xsec,vsec,dr,dv,L

 write(output,'("ptmass_",i0,".dat")') isec

 mpri = xyzmh_ptmass(4,ipri)
 msec = xyzmh_ptmass(4,isec)
! msec = xyzmh_ptmass(4,3)

 xpri(:) = xyzmh_ptmass(1:3,ipri)
 vpri(:) = vxyz_ptmass(1:3,ipri)
 xsec(:) = xyzmh_ptmass(1:3,isec)
! xsec(:) = xyzmh_ptmass(1:3,3)
! xsec(:) = xyzmh_ptmass(1:3,3)
 vsec(:) = vxyz_ptmass(1:3,isec)
! vsec(:) = vxyz_ptmass(1:3,3)
 dr(:) = xpri(:) - xsec(:)
 dv(:) = vpri(:) - vsec(:)
 rbin  = sqrt(dot_product(dr,dr))

! Calculate the binary specific relative ang. mom and energy
 call cross(dr,dv,L)
 Lmag = sqrt(dot_product(L,L))
 E = 0.5*dot_product(dv,dv) - G*(mpri+msec)/rbin

 if (abs(E) < tiny(E)) stop 'binary energy problem'

 if (E < 0) then ! check that system is bound
    call get_ae(Lmag,E,mpri,msec,a,ecc)

    inc = acos(L(3)/sqrt(dot_product(L,L)))

    if (time <= tiny(time)) then
       open(iunit,file=trim(output),status='replace',action='write',iostat=check)
       if (check /= 0) call fatal(analysistype,'unable to open binary.dat file at t=0.0')
       write(iunit,"('#',4(1x,'[',i2.2,1x,a11,']',2x))") &
            1,'time', &
            2,'a', &
            3,'ecc', &
            4,'inc'
    else
       inquire(file=trim(output),exist=exists)
       if (.not. exists) call fatal(analysistype,'t /= 0.0, but the analysis output file does not exist...')
       open(iunit,file=trim(output),status='old',action='write',position='append',iostat=check)
       if (check /= 0) call fatal(analysistype,'unable to open binary.dat file during run')
    endif
    write(iunit,'(4(ES18.10,1X))') time,a,ecc,inc
    close(iunit)
 endif
!-----------------------------------------------------------------------
end subroutine get_binary_params
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
subroutine get_ae(Lmag,E,m1,m2,a,ecc)
!-----------------------------------------------------------------------
! Return the semi-major axis and eccentricity between two objects
!-----------------------------------------------------------------------
 implicit none
 real,intent(out) :: a,ecc
 real,intent(in) :: Lmag,E,m1,m2

 if (Lmag < tiny(Lmag)) stop 'Lmag is zero in get_ae'
 if (abs(E) < tiny(E)) stop 'E is zero in get_ae'

! Hence obtain the binary eccentricity
 ecc = sqrt(1.0 + (2.0*E*Lmag**2)/((m1+m2)**2))

! and semi-major axis
 a = Lmag*Lmag/((m1+m2)*(1.0-ecc*ecc))
!-----------------------------------------------------------------------
end subroutine get_ae
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
subroutine cross(a,b,c)
!-----------------------------------------------------------------------
! Return the vector cross product of two 3d vectors
!-----------------------------------------------------------------------
 implicit none
 real,intent(in),dimension(3)  :: a,b
 real,intent(out),dimension(3) :: c

 c(1) = a(2)*b(3)-b(2)*a(3)
 c(2) = a(3)*b(1)-b(3)*a(1)
 c(3) = a(1)*b(2)-b(1)*a(2)

!-----------------------------------------------------------------------
end subroutine cross
!-----------------------------------------------------------------------
end module analysis
