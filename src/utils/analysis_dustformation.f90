!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine computing the time of dust formation for
! each SPH particle
!
! :References: Bermudez-Bustamante et al. (2023), submitted to MNRAS
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, fileutils, part, physcon, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'dustformation'
 public :: do_analysis

 private
 real, allocatable :: t_formation(:)

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,       only:do_nucleation
 use part,      only:nucleation,idK0,idK1
 use units,     only:utime
 use physcon,   only:years
 use fileutils, only:basename
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,n,ntot,iu,j
 real :: K0,K1,rdust
 character(len=30) :: filename

 if (.not.do_nucleation) then
    stop 'ERROR: need DUST_NUCLEATION=yes for this analysis type'
 endif

 if (.not.allocated(t_formation)) then
    allocate(t_formation(npart))
    t_formation = 0
 endif

 if (npart > size(t_formation)) then
    print*,' ERROR npart > npart_in, skipping analysis for '//trim(dumpfile)
    return
 endif

 n = 0
 ntot = 0
 do i=1,npart
    if (t_formation(i) <= 0.) then
       K0 = nucleation(idK0,i)
       K1 = nucleation(idK1,i)
       rdust = 1.28e-04*K1/(K0+tiny(0.))
       if (rdust > 1.e-3) then
          t_formation(i) = time*utime/years
          n = n + 1
       endif
    else
       ntot = ntot + 1
    endif
 enddo

 print*,' time is ',time*utime/years,' years'
 print*,' dust formation just started on ',n,' particles, total particles with dust = ',ntot+n,' / ',npart
 if (n > 1) then
    filename = trim(basename(dumpfile))
    j = index(filename,'_',back=.true.) - 1
    if (j <= 0) j = len_trim(filename)
    filename = filename(1:j)//'.comp'
    print*,' writing to '//trim(filename)
    open(newunit=iu,file=trim(filename),status='replace')
    write(iu,*) '#  t_{formation}'
    do i=1,npart
       write(iu,*) t_formation(i)
    enddo
    close(iu)
 endif

end subroutine do_analysis

end module analysis
