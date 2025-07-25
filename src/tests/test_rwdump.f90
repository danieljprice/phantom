!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testrwdump
!
! Unit test of read/write of particle data to/from dump files
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, dump_utils, eos, io, memory, mpiutils,
!   options, part, physcon, readwrite_dumps, testutils, timing, units
!
 implicit none
 public :: test_rwdump

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of reading/writing dump files
!+
!-----------------------------------------------------------------------
subroutine test_rwdump(ntests,npass)
 use part,            only:npart,npartoftype,massoftype,xyzh,hfact,vxyzu,&
                           Bevol,Bxyz,Bextx,Bexty,Bextz,alphaind,maxalpha,&
                           periodic,maxphase,mhd,maxvxyzu,igas,idust,&
                           maxp,poten,gravity,use_dust,dustfrac,xyzmh_ptmass,&
                           nptmass,nsinkproperties,xyzh_label,xyzmh_ptmass_label,&
                           dustfrac_label,vxyz_ptmass,vxyz_ptmass_label,&
                           vxyzu_label,set_particle_type,iphase,ndustsmall,ndustlarge,ndusttypes,&
                           iorig,copy_particle_all,norig
 use dim,             only:maxp,maxdustsmall,mpi
 use memory,          only:allocate_memory,deallocate_memory
 use testutils,       only:checkval,update_test_scores
 use io,              only:idisk1,id,master,iprint,nprocs
 use readwrite_dumps, only:read_dump,write_fulldump,write_smalldump,read_smalldump,is_small_dump
 use eos,             only:gamma,polyk
 use boundary,        only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax
 use units,           only:set_units,umass,udist
 use physcon,         only:au,solarm
 use mpiutils,        only:barrier_mpi
 use dump_utils,      only:read_array_from_file
 use timing,          only:getused,printused
 use options,         only:use_dustfrac
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(67)
 integer :: i,j,ierr,itest,ngas,ndust,ntot,maxp_old,iu
 real    :: tfile,hfactfile,time,tol,toldp
 real    :: alphawas,Bextxwas,Bextywas,Bextzwas,polykwas
 real    :: xminwas,xmaxwas,yminwas,ymaxwas,zminwas,zmaxwas
 logical :: test_speed
 real(kind=4) :: t1
 real, allocatable :: tmparray(:)

 if (id==master) write(*,"(/,a,/)") '--> TESTING READ/WRITE from dump file'
 test_speed = .false.

 ! This test will reallocate memory, so reset at the end
 maxp_old = maxp

 over_tests: do itest = 1,2

    if (test_speed) then
       ntot = maxp
    else
       ntot = 1001
    endif
    npart = ntot
    ndust = 10
    ngas  = ntot-ndust
    npartoftype(:) = 0
    npartoftype(igas) = ngas
    npartoftype(idust) = ndust
    do i=1,npart
       if (i <= npartoftype(1)) then
          call set_particle_type(i,igas)
       else
          call set_particle_type(i,idust)
       endif
    enddo
    hfact = 1.123
    massoftype(igas) = 30.
    massoftype(idust) = 40.
    gamma = 1.3
    polyk = 2.2
    alphawas = real(0.23_4)
    iu = 4
    do i=1,npart
       iorig(i) = i
       xyzh(1,i) = 1.
       xyzh(2,i) = 2.
       xyzh(3,i) = 3.
       xyzh(4,i) = 4.
       vxyzu(1,i) = 5.
       vxyzu(2,i) = 6.
       vxyzu(3,i) = 7.
       if (maxvxyzu >= 4) vxyzu(iu,i) = 8.
       if (maxalpha==maxp) then
          alphaind(1,i) = real(alphawas,kind=kind(alphaind)) ! 0->1
       endif
       if (mhd) then
          Bxyz(1,i) = 10.
          Bxyz(2,i) = 11.
          Bxyz(3,i) = 12.
          Bevol(4,i) = 13.
       endif
       if (gravity) then
          poten(i) = 15._4
       endif
       if (use_dust) then
          use_dustfrac = .true.
          ndustsmall = maxdustsmall
          ndustlarge = 1
          ndusttypes = ndustsmall + ndustlarge
          dustfrac(:,i) = 0.16_4
       endif
    enddo
    norig   = npart
    nptmass = 10
    do i=1,nptmass
       do j=1,nsinkproperties
          xyzmh_ptmass(j,i) = 100 + j
       enddo
       do j=1,3
          vxyz_ptmass(j,i) = 200 + j
       enddo
    enddo
    time = 20.
    tol = real(tiny(0._4))

    if (periodic) then
       xminwas = -1.
       xmaxwas = 2.
       yminwas = -3.
       ymaxwas = 4.
       zminwas = -5.
       zmaxwas = 6.
       call set_boundary(xminwas,xmaxwas,yminwas,ymaxwas,zminwas,zmaxwas)
    endif
    if (mhd) then
       Bextx = 124.
       Bexty = 125.
       Bextz = 126.
       Bextxwas = Bextx
       Bextywas = Bexty
       Bextzwas = Bextz
    endif
    polykwas = polyk
    call set_units(dist=au,mass=solarm,G=1.d0)
    call copy_particle_all(ngas,2,.false.)   ! Move i=ngas to i=2 (assumes i=2 was killed)
    call copy_particle_all(1,ngas,.true.)    ! Create new i=ngas based upon i=1
!
!--write to file
!
    if (test_speed) call getused(t1)
    if (itest==2) then
       call write_smalldump(time,'test.dump')
    else
       call write_fulldump(time,'test.dump')
    endif
    if (test_speed) call printused(t1)
!
!--scrub variables
!
    massoftype(:) = 0.
    npartoftype(:) = 0
    xyzh  = 0.
    vxyzu = 0.
    polyk = 0.
    if (maxalpha==maxp) alphaind = 0.
    if (mhd) then
       Bevol = 0.
       Bxyz  = 0.
    endif
    if (gravity) then
       poten = 0.
    endif
    if (use_dust) then
       dustfrac(:,:) = 0.
    endif
    gamma = 0.
    Bextx = 0.
    Bexty = 0.
    Bextz = 0.
    call set_units()
    xyzmh_ptmass = 0.
    vxyz_ptmass = 0.
!
!--read from same file
!
    if (test_speed) call getused(t1)
    if (itest==2) then
       if (id==master) write(*,"(/,a)") '--> checking read_dump'
       nfailed = 0
       call deallocate_memory
       call read_dump('test.dump',tfile,hfactfile,idisk1,iprint,id,nprocs,ierr)
       call checkval(ierr,is_small_dump,0,nfailed(1),'read_dump returns is_small_dump error code')
       call update_test_scores(ntests,nfailed,npass)

       if (id==master) write(*,"(/,a)") '--> checking read_smalldump'
       call read_smalldump('test.dump',tfile,hfactfile,idisk1,iprint,id,nprocs,ierr)
       toldp = epsilon(0._4)
    else
       call deallocate_memory
       call read_dump('test.dump',tfile,hfactfile,idisk1,iprint,id,nprocs,ierr)
       toldp = tiny(toldp)
    endif
    if (test_speed) call printused(t1)
!
!--check that values match
!
    nfailed(:) = 0
    if (id==master) write(*,"(/,a)") '--> checking header'
    call checkval(ierr,0,0,nfailed(1),'read without errors')
    call checkval(tfile,time,toldp,nfailed(2),'time')
    call checkval(hfactfile,hfact,toldp,nfailed(3),'hfact')
    call checkval(gamma,1.3,toldp,nfailed(4),'gamma')
    call checkval(polyk,polykwas,toldp,nfailed(5),'polyk')
    call checkval(massoftype(igas),30.,toldp,nfailed(6),'m(gas)')
    call checkval(massoftype(idust),40.,toldp,nfailed(7),'m(dust)')
    call checkval(npartoftype(1),ngas,0,nfailed(8),'npartoftype(1)')
    call checkval(npartoftype(idust),ndust,0,nfailed(9),'npartoftype(2)')
    call checkval(umass,solarm,tiny(umass),nfailed(10),'umass')
    call checkval(udist,au,tiny(udist),nfailed(11),'udist')
    if (periodic) then
       call checkval(xmin,xminwas,tiny(xmin),nfailed(12),'xmin')
       call checkval(xmax,xmaxwas,tiny(xmax),nfailed(13),'xmax')
       call checkval(ymin,yminwas,tiny(ymin),nfailed(14),'ymin')
       call checkval(ymax,ymaxwas,tiny(ymax),nfailed(15),'ymax')
       call checkval(zmin,zminwas,tiny(zmin),nfailed(16),'zmin')
       call checkval(zmax,zmaxwas,tiny(zmax),nfailed(17),'zmax')
       call set_boundary()  ! reset to defaults
    endif
    if (mhd) then
       call checkval(Bextx,Bextxwas,tiny(Bextx),nfailed(18),'Bextx')
       call checkval(Bexty,Bextywas,tiny(Bexty),nfailed(19),'Bexty')
       call checkval(Bextz,Bextzwas,tiny(Bextz),nfailed(20),'Bextz')
    endif
    if (itest==1) then  ! iorig is not dumped to small dumps
       call checkval(iorig(2),    int(ngas   ,kind=8), 0,nfailed(65),'iorig(2)')
       call checkval(iorig(ngas), int(npart+1,kind=8), 0,nfailed(66),'iorig(ngas)')
       call checkval(iorig(npart),int(npart  ,kind=8), 0,nfailed(67),'iorig(N)')
    endif

    call update_test_scores(ntests,nfailed,npass)

    call barrier_mpi()

    if (id==master) write(*,"(/,a)") '--> checking particle arrays'
    if (all(xyzh(4,:) <= 0.)) xyzh(4,:) = 1. ! otherwise indicates dead particles and not checked
    nfailed(:) = 0
    do j=1,4
       call checkval(npart,xyzh(j,:),real(j),tiny(xyzh),nfailed(j),xyzh_label(j))
    enddo
    if (itest /= 2) then
       do j=1,maxvxyzu
          call checkval(npart,vxyzu(j,:),real(4+j),tiny(vxyzu),nfailed(5+j),vxyzu_label(j))
       enddo
       if (maxalpha==maxp) then
          call checkval(npart,alphaind(1,:),alphawas,tol,nfailed(9),'alpha')
       endif
       if (mhd) then
          call checkval(npart,Bxyz(1,:),10.,tol,nfailed(10),'Bx')
          call checkval(npart,Bxyz(2,:),11.,tol,nfailed(11),'By')
          call checkval(npart,Bxyz(3,:),12.,tol,nfailed(12),'Bz')
          call checkval(npart,Bevol(4,:),13.,tol,nfailed(13),'psi')
       endif
       if (gravity) then
          call checkval(npart,poten,15.,tol,nfailed(15),'poten')
       endif
       if (use_dust) then
          do i = 1,ndustsmall
             call checkval(npart,dustfrac(i,:),0.16,real(epsilon(0._4)),nfailed(16),'dustfrac')
          enddo
       endif
    endif
    if (maxphase==maxp) then
       call checkval(ngas,iphase,igas,0,nfailed(17),'particle type 1')
       call checkval(ndust,iphase(npart-ndust+1:npart),idust,0,nfailed(18),'particle type 2')
    endif
    call update_test_scores(ntests,nfailed,npass)

    if (id==master) write(*,"(/,a)") '--> checking sink particle arrays'
    nfailed = 0
    do j=1,nsinkproperties
       call checkval(nptmass,xyzmh_ptmass(j,:),real(100+j),tiny(xyzmh_ptmass),nfailed(j),xyzmh_ptmass_label(j))
    enddo
    if (itest /= 2) then
       do j=1,3
          call checkval(nptmass,vxyz_ptmass(j,:),real(200+j),tiny(vxyz_ptmass),&
                        nfailed(nsinkproperties+j),vxyz_ptmass_label(j))
       enddo
    endif
    call update_test_scores(ntests,nfailed,npass)

    call barrier_mpi()

    ! clean up doggie-doos
    if (use_dust) then
       dustfrac(:,:) = 0.
    endif
    nptmass = 0
    xyzmh_ptmass = 0.
    vxyz_ptmass = 0.

    ! test read of a single array from the file
    if (.not. mpi .and. itest==1) then
       if (id==master) write(*,"(/,a)") '--> checking read of single array from file'
       xyzh(2,:) = 0.
       nfailed = 0
       allocate(tmparray(size(xyzh,dim=2)))
       call read_array_from_file(idisk1,'test.dump','y',tmparray,ierr)
       xyzh(1,:) = tmparray
       deallocate(tmparray)
       call checkval(ierr,0,0,nfailed(1),'error flag')
       call checkval(npart,xyzh(1,:),2.,tiny(xyzh),nfailed(2),'y')
       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1
    endif

    if (id==master) then
       open(unit=idisk1,file='test.dump',status='old')
       close(unit=idisk1,status='delete')
    endif
 enddo over_tests

 call deallocate_memory
 call allocate_memory(int(maxp_old,kind=8))
 if (id==master) write(*,"(/,a)") '<-- READ/WRITE TEST COMPLETE'

end subroutine test_rwdump

end module testrwdump
