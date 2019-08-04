!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testlink
!
!  DESCRIPTION:
!  This module performs unit tests of the link list routines
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, domain, io, kdtree, kernel, linklist,
!    mpiutils, part, random, testutils, timing, unifdis
!+
!--------------------------------------------------------------------------
module testlink
 implicit none
 public :: test_link

 private

contains

subroutine test_link(ntests,npass)
 use dim,      only:maxp,maxneigh,periodic
 use io,       only:id,master,nprocs!,iverbose
 use mpiutils, only:reduceall_mpi
 use part,     only:npart,npartoftype,massoftype,xyzh,vxyzu,hfact,igas,kill_particle
 use kernel,   only:radkern2,radkern
 use unifdis,  only:set_unifdis
 use timing,   only:getused
 use random,   only:ran2
 use domain,   only:i_belong
 use part,            only:maxphase,iphase,isetphase,igas,iactive
 use testutils,       only:checkval,checkvalbuf_start,checkvalbuf,checkvalbuf_end,update_test_scores
 use linklist,        only:set_linklist,get_neighbour_list,ncells,ifirstincell
 use kdtree,          only:inodeparts,inoderange
 use domain,          only:i_belong
#ifdef PERIODIC
 use boundary, only:xmin,xmax,ymin,ymax,zmin,zmax,dybound,dzbound
 use linklist, only:dcellx,dcelly,dcellz
#endif
 use boundary, only:dxbound
 use part,     only:isdead_or_accreted
#ifdef IND_TIMESTEPS
 use part,     only:iphase_soa
#endif
 integer, intent(inout) :: ntests,npass
 real                   :: psep,hzero,totmass,dxboundp,dyboundp,dzboundp
 real                   :: xminp,xmaxp,yminp,ymaxp,zminp,zmaxp
 real                   :: rhozero,hi21,dx,dy,dz,xi,yi,zi,q2,hmin,hmax,hi
 integer                :: i,j,icell,ixyzcachesize,ncellstest,nfailedprev,maxpen
 integer                :: nneigh,nneighexact,nneightry,max1,max2,ncheck1,ncheck2,nwarn
 integer                :: ip
 integer                :: nparttot
 integer(kind=8)        :: nptot
#ifdef IND_TIMESTEPS
 integer                :: npartincell,ierrmax
 logical                :: hasactive
#endif
 integer                :: maxneighi,minneigh,iseed,nlinktest,itest,ndead
 integer(kind=8)        :: meanneigh,i8
 integer :: nfailed(8)
 logical                :: iactivei,iactivej,activecell
 real, allocatable :: xyzcache(:,:)
 integer :: listneigh(maxneigh)
 character(len=1), dimension(3), parameter :: xlabel = (/'x','y','z'/)

 if (id==master) write(*,"(a,/)") '--> TESTING LINKLIST / NEIGHBOUR FINDING'
!
!--set up a random particle distribution
!
 npart = 0
 nptot = 0
#ifdef PERIODIC
 xminp = xmin
 xmaxp = xmax
 yminp = ymin
 ymaxp = ymax
 zminp = zmin
 zmaxp = zmax
#else
 xminp = -1.
 xmaxp = 1.
 yminp = -2.
 ymaxp = 0.
 zminp = 1.
 zmaxp = 3.
#endif
 dxboundp = xmaxp-xminp
 dyboundp = ymaxp-yminp
 dzboundp = zmaxp-zminp
 psep = (xmaxp-xminp)/32.

 call set_unifdis('random',id,master,xminp,xmaxp,yminp,ymaxp,zminp,zmaxp,psep,&
                  hfact,npart,xyzh,periodic,nptot=nptot,mask=i_belong)
 npartoftype(:) = 0
 npartoftype(igas) = npart
 !print*,'thread ',id,' npart = ',npart
 !iverbose = 3

 rhozero = 7.5
 hfact = 1.2
 totmass = rhozero/(dxboundp*dyboundp*dzboundp)
 massoftype(igas) = totmass/reduceall_mpi('+',npart)
 hzero = hfact*(massoftype(1)/rhozero)**(1./3.)

 hmin = 0.01*hzero
 hmax = 0.2/dxboundp !0.25/dxboundp

#ifdef IND_TIMESTEPS
 nlinktest = 3
#else
 nlinktest = 2
#endif

 over_tests: do itest=1,nlinktest

    iseed = -24358
    ip = 0
    do i8=1_8,nptot
       hi = hmin + ran2(iseed)*(hmax - hmin)
       if (i_belong(i8)) then
          ip = ip + 1
          !--give random smoothing lengths
          xyzh(4,ip) = hi
       endif
    enddo
    npart = ip

#ifdef IND_TIMESTEPS
!----------------------------------------------------
! TEST 1: WITH ALL PARTICLES ACTIVE
! TEST 2: SOME PARTICLES DEAD OR ACCRETED
! TEST 3: WITH ONLY A FRACTION OF PARTICLES ACTIVE
!----------------------------------------------------
    do i=1,npart
       if (itest==3) then
          !--partially active
          if (xyzh(4,i) < (hmin + 0.2*(hmax-hmin))) then
             iphase(i) = isetphase(igas,iactive=.true.)
          else
             iphase(i) = isetphase(igas,iactive=.false.)
          endif
       elseif (itest==2) then
          !--mark a number of particles as dead or accreted
          if (xyzh(4,i) > (hmin + 0.2*(hmax-hmin))) xyzh(4,i) = -abs(xyzh(4,i))
          if (mod(i,1000)==0) call kill_particle(i)
          iphase(i) = isetphase(igas,iactive=.true.)
       else
          !--all active
          iphase(i) = isetphase(igas,iactive=.true.)
       endif
    enddo
#else
    if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
#endif
!
!--setup the link list
!
    if (id==master) write(*,"(/,1x,2(a,i1),a,/)") 'Test ',itest,' of ',nlinktest,': building linked list...'
    call set_linklist(npart,npart,xyzh,vxyzu)
!
!--count dead particles
!
    ndead = 0
    do i=1,npart
       if (isdead_or_accreted(xyzh(4,i))) ndead = ndead + 1
    enddo
!
!--check that the number of cells is non-zero
!
    call checkval((ncells>0),.true.,nfailed(1),'ncells > 0')
    call update_test_scores(ntests,nfailed(1:1),npass)
!
!--check the assignment of positive or negative
!  to the head of the cell that specifies whether
!  or not the cell contains active particles
!
#ifdef IND_TIMESTEPS
    if (itest /= 2) then
       ncheck1 = 0
       ncheck2 = 0
       nfailed = 0
       call checkvalbuf_start('active/inactive cells')
       !!$omp parallel do default(none) &
       !!$omp shared(ncells,ifirstincell,iphase,ll) &
       !!$omp private(i,activecell,hasactive,iactivei,npartincell) &
       !!$omp reduction(+:nfail,ncheck1,ncheck2)
       do icell=1,int(ncells)
          if (ifirstincell(icell) < 0) then
             activecell = .false.
          else
             activecell = .true.
          endif
          npartincell = 0
          hasactive   = .false.
          not_empty: if (ifirstincell(icell) /= 0) then
             do i = inoderange(1,icell), inoderange(2,icell)
                npartincell = npartincell + 1
                iactivei = iactive(iphase_soa(i))
                if (iactivei) hasactive = .true.
                if (.not.activecell) then
                   call checkvalbuf(iactivei,.false.,'inactive cell contains active particle',nfailed(1),ncheck1)
                endif
             enddo
             if (activecell .and. npartincell > 0) then
                call checkvalbuf(hasactive,.true.,'active cell has at least one active particle',nfailed(2),ncheck2)
             endif
          endif not_empty
       enddo
       !!$omp end parallel do
       ierrmax = 0
       if (ncheck1 > 0) then
          call checkvalbuf_end('inactive cells have no active particles',ncheck1,nfailed(1),ierrmax,0,npart)
       endif
       call checkvalbuf_end('active cells have at least one active particle',ncheck2,nfailed(2),ierrmax,0,npart)

       call update_test_scores(ntests,nfailed(1:2),npass)
    endif
#endif

    dochecks: if (itest /= 2) then

       ixyzcachesize = 60000*int((radkern/2.0)**3)
       if (.not.allocated(xyzcache)) allocate(xyzcache(ixyzcachesize,3))
!
!--now pick a sample of particles, find their neighbours via "get_neighbour_list" and
!  check it via a direct evaluation
!
       ncellstest = 10
       nfailed(:) = 0
       max1 = 0
       max2 = 0
       ncheck1 = 0
       ncheck2 = 0
       maxneighi = 0
       minneigh  = huge(minneigh)
       meanneigh = 0
       call checkvalbuf_start('neighbour number')
       nwarn = 0
       listneigh(:) = 0

       over_cells: do icell=1,int(ncells)
          i = ifirstincell(icell)
          if (i==0) cycle over_cells

          if (i < 0) then
             activecell = .false.
             i = -i
          else
             activecell = .true.
          endif

          call get_neighbour_list(icell,listneigh,nneightry,xyzh,xyzcache,ixyzcachesize)

          over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
             i = inodeparts(ip)
             iactivei = .true.
#ifdef IND_TIMESTEPS
             i = abs(i)
             iactivei = iactive(iphase(i))
#endif
             xi = xyzh(1,i)
             yi = xyzh(2,i)
             zi = xyzh(3,i)
             hi = xyzh(4,i)
             hi21 = 1./(hi*hi)
             !
             !--first work out the correct answer:
             !  i.e., the actual number of neighbours
             !  by a direct summation over all particles
             !
             nneighexact = 0
             iactivej    = .true.

             do j=1,npart
#ifdef IND_TIMESTEPS
                iactivej = iactive(iphase(j))
#endif
                !
                !--an active cell should return a list of both active
                !  and inactive neighbours. An inactive cell should
                !  get contributions from active neighbours ONLY
                !
                if (activecell .or. iactivej) then
                   dx = xi - xyzh(1,j)
                   dy = yi - xyzh(2,j)
                   dz = zi - xyzh(3,j)
#ifdef PERIODIC
                   if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
                   if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
                   if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
                   q2 = (dx*dx + dy*dy + dz*dz)*hi21
                   if (q2 < radkern2) then
                      nneighexact = nneighexact + 1
                   endif
                endif
             enddo
             maxneighi = max(nneighexact,maxneighi)
             minneigh  = min(nneighexact,minneigh)
             meanneigh = meanneigh + nneighexact

             !
             !--get the number of actual neighbours from
             !  the trial list using the neighbour cache
             !  (with spillover, exactly as in the code)
             !
             nneigh = 0
             do j=1,nneightry
#ifdef IND_TIMESTEPS
                iactivej = iactive(iphase(listneigh(j)))
#endif
                if (activecell .or. iactivej) then
                   if (j <= ixyzcachesize) then
                      dx = xi - xyzcache(j,1)
                      dy = yi - xyzcache(j,2)
                      dz = zi - xyzcache(j,3)
                   else
                      dx = xi - xyzh(1,listneigh(j))
                      dy = yi - xyzh(2,listneigh(j))
                      dz = zi - xyzh(3,listneigh(j))
                   endif
#ifdef PERIODIC
                   if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
                   if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
                   if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
                   q2 = (dx*dx + dy*dy + dz*dz)*hi21
                   if (q2 < radkern2) then
                      nneigh = nneigh + 1
                   endif
                endif
             enddo
             call checkvalbuf(nneigh,nneighexact,0,'nneigh (cached)',nfailed(1),ncheck1,max1)
             !
             !--get the number of actual neighbours using the
             !  un-cached values of xyz
             !
             nneigh = 0
#ifdef PERIODIC
             if (radkern*hi < min(0.5*dxbound-2.*dcellx,0.5*dybound-2.*dcelly,0.5*dzbound-2.*dcellz)) then
#endif
                do j=1,nneightry
#ifdef IND_TIMESTEPS
                   iactivej = iactive(iphase(listneigh(j)))
#endif
                   if (activecell .or. iactivej) then
                      dx = xi - xyzh(1,listneigh(j))
                      dy = yi - xyzh(2,listneigh(j))
                      dz = zi - xyzh(3,listneigh(j))
#ifdef PERIODIC
                      if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
                      if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
                      if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
                      q2 = (dx*dx + dy*dy + dz*dz)*hi21
                      if (q2 < radkern2) then
                         nneigh = nneigh + 1
                      endif
                   endif
                enddo
                nfailedprev = nfailed(2)
                call checkvalbuf(nneigh,nneighexact,0,'nneigh (no cache)',nfailed(2),ncheck2,max2)

                if (nfailed(2) > nfailedprev) then
                   !
                   !--check for double counting in neighbour list
                   !
                   do j=1,nneightry
                      if (nwarn < 20) then
                         if (any(listneigh(1:j-1)==listneigh(j))) then
                            nwarn = nwarn + 1
                            print*,' ERROR: double counting in neighbour list ',listneigh(j) !,listneigh(1:j-1)
                         endif
                      endif
                   enddo
                endif
#ifdef PERIODIC
             endif
#endif
          enddo over_parts
       enddo over_cells

       call checkvalbuf_end('nneigh (cached)',  ncheck1,nfailed(1),max1,0)
       call checkvalbuf_end('nneigh (no cache)',ncheck2,nfailed(2),max2,0,npart)
       write(*,"(1x,2(a,i6),a,f9.2)") 'max nneigh = ',maxneighi,&
    ' min nneigh = ',minneigh,' mean = ',meanneigh/real(ncheck2)

       call update_test_scores(ntests,nfailed(1:2),npass)

    endif dochecks

 enddo over_tests

!
!--check neighbour finding with some pathological configurations
!
 nlinktest = nlinktest + 1
 if (id==master) write(*,"(/,1x,a,i2,a,/)") 'Test ',nlinktest,': building linked list...'
 do maxpen=1,3
    if (id==master) write(*,"(a)") ' particles in a line in '//xlabel(maxpen)//' direction '

    !--particles in a line
    nparttot = 20 * nprocs
    psep  = dxbound/nparttot
    massoftype(1)   = 2.
    npart = 0
    do i=1,nparttot
       if (mod(i,nprocs) == id) then
          npart = npart + 1
          xyzh(:,npart) = 0.
          xyzh(maxpen,npart)  = (i-1)*psep
          xyzh(4,npart)       = hfact*psep
          if (maxphase==maxp) iphase(npart) = isetphase(igas,iactive=.true.)
       endif
    enddo
    npartoftype(:) = 0
    npartoftype(igas) = npart

    call set_linklist(npart,npart,xyzh,vxyzu)
    !
    !--check that the number of cells is non-zero
    !
    call checkval((ncells>0),.true.,nfailed(1),'ncells > 0')
    call update_test_scores(ntests,nfailed(1:1),npass)

 enddo

 if (allocated(xyzcache)) deallocate(xyzcache)

 if (id==master) write(*,"(/,a,/)") '<-- LINKLIST TEST COMPLETE'

end subroutine test_link

end module testlink
