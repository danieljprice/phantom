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
! :Dependencies: discanalysisutils, io, part
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
 use discanalysisutils, only:get_binary_params
 character(len=*), intent(in) :: dumpfile
 real,     intent(in) :: pmass,time
 real,     intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer,  intent(in) :: numfile,npart,iunit
 integer :: i
 real    :: G,a,ecc
 character(len=25) :: output

 ! Use two variables solely to remove compiler warnings...
 write(*,'("Performing analysis on ",a,"... which is unit ",i5,"...")') trim(dumpfile),iunit

! Check nsinks is >=2
 if (nptmass < 2) call fatal(analysistype,'Not enough sinks...')

! Assuming G=1.
 G = 1.0

! Currently assuming all sinks orbit sink1
 do i = 2,nptmass
    write(output,'("ptmass_",i0,".dat")') i
    call get_binary_params(1,i,xyzmh_ptmass,vxyz_ptmass,time,a,ecc,G,output)
 enddo

end subroutine do_analysis

end module analysis
