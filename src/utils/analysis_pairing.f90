!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to check for particle pairing, i.e. whether particles
! have a very close neighbour (this can cause problems in post-processing)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: linklist, part
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'pairing'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:set_linklist,getneigh_pos,ifirstincell,listneigh=>listneigh_global
 use part,     only:isdead_or_accreted
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: maxcache = 0
 real, allocatable  :: xyzcache(:,:)
 integer :: nneigh,i,n,j,nwarn,nbin(7),ncheck
 integer(kind=8) :: nneightot
 real :: r2,h2,dx(3)
 real, parameter :: sep_hist(7) = (/0.0001,0.001,0.01,0.1,0.3,1.0,2.0/)

 print*,'> Building tree... '
 call set_linklist(npart,npart,xyzh,vxyzu)

 nwarn = 0
 nbin = 0
 ncheck = 0
 nneightot = 0
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call getneigh_pos(xyzh(1:3,i),0.,xyzh(4,i),3,listneigh,nneigh,xyzh,xyzcache,maxcache,ifirstincell)
       h2 = xyzh(4,i)**2
       ncheck = ncheck + 1
       do n = 1,nneigh
          j = listneigh(n)
          if (j==i) cycle
          dx = xyzh(1:3,i)-xyzh(1:3,j)
          r2 = dot_product(dx,dx)
          where (r2 < sep_hist(:)**2*h2)
             nbin(:) = nbin(:) + 1
          end where
          nneightot = nneightot + 1
          if (r2 < 0.0001*h2 .and. j > i) then
             print*,i,' WARNING: close neighbour ',j,' at r/h=',sqrt(r2/h2),', r = ',sqrt(r2)
             !print*,'xi = ',xyzh(1:4,i)
             !print*,'xj = ',xyzh(1:4,j)
             nwarn = nwarn + 1
          endif
       enddo
    endif
 enddo

 print*,' average number of neighbours = ',nneightot/real(ncheck)
 print "(/,1x,50('-'),/,6x,'r/h',6x,' nneigh ',5x,' percentage',/,1x,50('-'))"
 do i=1,size(sep_hist)
    print "(1x,1pg8.2,3x,i10,5x,1pg8.2,a)",sep_hist(i),nbin(i)/2,nbin(i)/real(nneightot)*100.,' %'
 enddo
 print "(1x,50('-'),/)"

 print*,' WARNING: found pairing (r/h < 0.01) on ',nwarn,' particle pairs'

 !print*,'> Printing kd-tree '
 !do i=1,numnodes
 !   print*,node(i)
 !enddo

end subroutine do_analysis

end module analysis
