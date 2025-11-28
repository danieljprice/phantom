!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for sinks: assumes central sink with label 1,
!                              and all sinks orbit this
!                              assumes G=1.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, options, part, physcon, setbinary
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'ptmass'
 public :: do_analysis

 private

contains

!-----------------------------------------------------------------------
!+
!  print the binary parameters for each pair of sink particles
!+
!-----------------------------------------------------------------------
subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,   only:fatal
 use part, only:xyzmh_ptmass,vxyz_ptmass,nptmass
 character(len=*), intent(in) :: dumpfile
 real,     intent(in) :: pmass,time
 real,     intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer,  intent(in) :: numfile,npart,iunit
 integer :: i
 real    :: G

! Use two variables solely to remove compiler warnings...
 write(*,'("Performing analysis on ",a,"... which is unit ",i5,"...")') trim(dumpfile),iunit

! Check nsinks is >=2
 if (nptmass < 2) call fatal(analysistype,'Not enough sinks...')

! Assuming G=1.
 G = 1.0

! Currently assuming all sinks orbit sink1
 do i = 2,nptmass
    call get_binary_params(1,i,xyzmh_ptmass,vxyz_ptmass,time,G)
 enddo

end subroutine do_analysis

!-----------------------------------------------------------------------
!+
!  calculate the orbital parameters for a pair of sink particles
!+
!-----------------------------------------------------------------------
subroutine get_binary_params(ipri,isec,xyzmh_ptmass,vxyz_ptmass,time,G)
 use io,     only:fatal
 use orbits, only:get_orbital_elements
 integer, intent(in) :: ipri,isec
 real,    intent(in) :: time,G
 real,    intent(in) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, parameter :: iunit = 150
 logical :: exists
 integer :: check
 real :: m1,m2,a,ecc,inc,Omega,w,f
 real :: dr(3),dv(3)
 character(len=25) :: output

 write(output,'("ptmass_",i0,".dat")') isec

 m1 = xyzmh_ptmass(4,ipri)
 m2 = xyzmh_ptmass(4,isec)

 dr(:) = xyzmh_ptmass(1:3,ipri) - xyzmh_ptmass(1:3,isec)
 dv(:) = vxyz_ptmass(1:3,ipri) - vxyz_ptmass(1:3,isec)

 call get_orbital_elements(G*(m1+m2),dr,dv,a,ecc,inc,Omega,w,f)

 if (time <= tiny(time)) then
    open(iunit,file=trim(output),status='replace',action='write',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open binary.dat file at t=0.0')
    write(iunit,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time', &
          2,'a', &
          3,'ecc', &
          4,'inc', &
          5,'Omega', &
          6,'w', &
          7,'f'
 else
    inquire(file=trim(output),exist=exists)
    if (.not. exists) call fatal(analysistype,'t /= 0.0, but the analysis output file does not exist...')
    open(iunit,file=trim(output),status='old',action='write',position='append',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open binary.dat file during run')
 endif
 write(iunit,'(7(es18.10,1x))') time,a,ecc,inc,Omega,w,f
 close(iunit)

end subroutine get_binary_params

end module analysis
