!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: io, metric_tools, part, utils_gr
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'etotgr'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,           only:warning,iprint
 use part,         only:pxyzu,metrics,metricderivs
 use metric_tools, only:init_metric
 use utils_gr,     only:get_u0
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 real :: pdotv,u,etot,e,gcov(0:3,0:3),U0
 integer :: ierr,i

 call init_metric(npart,xyzh,metrics,metricderivs)

 etot = 0.
 do i=1,npart
    pdotv = dot_product(pxyzu(1:3,i),vxyzu(1:3,i))
    gcov  = metrics(:,:,1,i)
    call get_u0(gcov,vxyzu(1:3,i),U0,ierr)
    u = vxyzu(4,i)
    e = pmass*(pdotv + (1.+u)/U0 - 1.)
    etot = etot + e
 enddo

 write(1,*) time,etot

end subroutine do_analysis

end module analysis
