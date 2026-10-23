!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testapr
!
! Unit test for adaptive particle refinement
!
! :References:
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: apr, boundary, deriv, dim, energies, io, mpidomain,
!   mpiutils, options, part, physcon, random, testutils, unifdis, utils_apr
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master,fatal
 implicit none
 public :: test_apr,setup_apr_region_for_test

 private

contains

!--------------------------------------------
!+
!  Various tests of the apr module
!+
!--------------------------------------------
subroutine test_apr(ntests,npass)
 use unifdis,      only:set_unifdis
 use boundary,     only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use part,         only:npart,npartoftype,hfact,xyzh,init_part,massoftype,radprop
 use part,         only:isetphase,igas,iphase,vxyzu,fxyzu,apr_level,maxvxyzu,init_rho_from_h
 use mpidomain,    only:i_belong
 use mpiutils,     only:reduceall_mpi
 use dim,          only:periodic,use_apr,do_radiation,use_dust
 use apr,          only:update_apr
 use utils_apr,    only:apr_centre
 use energies,     only:compute_energies,angtot,etot,totmom,ekin,etherm
 use random,       only:ran2
 use deriv,        only:get_density_global
 integer, intent(inout) :: ntests,npass
 real :: psep,rhozero,time,totmass, etotin, totmomin
 real :: angtotin, ekinin, ethermin
 real :: tolang, tolen, tolmom, hcoarse,hfine
 integer :: original_npart,splitted,nfailed(8),i,iseed,ncoarse,nfine

 if (use_apr) then
    if (id==master) write(*,"(/,a)") '--> TESTING APR MODULE'
 else
    if (id==master) write(*,"(/,a)") '--> SKIPPING APR TEST (REQUIRES -DAPR)'
    return
 endif

 ! Tolerances
 tolmom = 2.e-15
 tolang = 8.0e-14
 tolen  = 2.e-15
 nfailed(:) = 0
 iseed = -92757

 ! Set up a uniform box of particles
 call init_part()
 psep = dxbound/20.
 time = 0.
 npartoftype(:) = 0
 npart = 0
 rhozero = 1.0
 totmass = rhozero*dxbound*dybound*dzbound
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact,npart,xyzh,periodic,mask=i_belong)

 original_npart = npart
 npartoftype(igas) = npart
 massoftype(igas) = totmass/reduceall_mpi('+',npart)
 iphase(1:npart) = isetphase(igas,iactive=.true.)
 call init_rho_from_h()

 ! this is to prevent a (reasonable) problem when running this test with DEBUG=yes and radiation
 if (do_radiation) then
    radprop(4,:) = 23.0421 ! just some inconsequential number
 endif

! phases are set explicitly via isetphase/iphase in the APR tests

 ! Set some random velocities
 do i=1,npart
    vxyzu(1:3,i) = (/ran2(iseed),ran2(iseed),ran2(iseed)/)
    if (maxvxyzu > 3) vxyzu(4,i) = ran2(iseed)**2
 enddo

 ! Initialise APR
 call setup_apr_region_for_test()
 apr_centre(:,1:2) = 20. ! just moves the APR region away from the box so you don't have any split or merge
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)

 !
 ! mixed masses on a uniform lattice: neighbour counts (hence h) should match,
 ! unlike h(rho) which would split by m^{1/3}
 !
 write(*,"(/,a)") '--> checking h from number density at mixed APR levels'
 do i=1,npart
    if (mod(i,2)==0) apr_level(i) = 2
 enddo
 call get_density_global(1)
 ncoarse = 0
 nfine   = 0
 hcoarse = 0.
 hfine   = 0.
 do i=1,npart
    if (apr_level(i) == 1) then
       ncoarse = ncoarse + 1
       hcoarse = hcoarse + xyzh(4,i)
    else
       nfine = nfine + 1
       hfine = hfine + xyzh(4,i)
    endif
 enddo
 ncoarse = int(reduceall_mpi('+',ncoarse))
 nfine   = int(reduceall_mpi('+',nfine))
 hcoarse = reduceall_mpi('+',hcoarse)
 hfine   = reduceall_mpi('+',hfine)
 if (ncoarse > 0) hcoarse = hcoarse/ncoarse
 if (nfine > 0) hfine = hfine/nfine
 if (ncoarse > 0 .and. nfine > 0) then
    call checkval(hfine,hcoarse,5.e-2,nfailed(8),'h fine vs coarse (number density)')
 else
    nfailed(8) = 1
    if (id==master) write(*,*) 'ERROR: mixed APR levels not present for h(n) test'
 endif
 call update_test_scores(ntests,nfailed(8:8),npass)
 apr_level(1:npart) = 1
 call init_rho_from_h()

 ! Initialise the energies values
 call compute_energies(0.)
 etotin   = etot
 totmomin = totmom
 angtotin = angtot
 ekinin = ekin
 ethermin = etherm

 ! Now set for a split
 write(*,"(/,a)") '--> conducting a split'
 apr_centre(1:2,1) = 0.25    ! this puts a sphere centred at (0.25,0.25)
 apr_centre(1:2,2) = -0.25   ! and a second sphere at (-0.25,-0.25)
 apr_centre(3,1:2) = 0.      ! and ensures they are in the plane
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)

 ! Check the new conserved values
 call compute_energies(0.)
 call checkval(angtot,angtotin,tolang,nfailed(1),'angular momentum')
 call checkval(totmom,totmomin,tolmom,nfailed(2),'linear momentum')
 call checkval(etot,etotin,tolen,nfailed(3),'total energy')
 call checkval(ekin,ekinin,tolen,nfailed(4),'kinetic energy')
 call checkval(etherm,ethermin,tolen,nfailed(5),'thermal energy')
 call update_test_scores(ntests,nfailed(1:5),npass)

 ! after splitting, the total number of particles should have been updated
 splitted = npart

 ! Move the apr zone out of the box and update again to merge
 write(*,"(/,a)") '--> conducting a merge'
 apr_centre(:,1:2) = 20. ! move the APR zones away again
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)

 ! Check the new conserved values
 call compute_energies(0.)
 !call checkval(angtot,angtotin,tolang,nfailed(1),'angular momentum')
 call checkval(totmom,totmomin,tolmom,nfailed(5),'linear momentum')
 !call checkval(etot,etotin,tolen,nfailed(3),'total energy')
 !call checkval(ekin,ekinin,tolen,nfailed(4),'kinetic energy')
 call checkval(etherm,ethermin,tolen,nfailed(6),'thermal energy')

 ! Check that the original particle number returns
 call checkval(npart,original_npart,0,nfailed(7),'number of particles == original number')
 call update_test_scores(ntests,nfailed(5:7),npass)

 ! nearest neighbour search and refinement by several levels at once
 call test_nearest_neighbour(ntests,npass)
 call test_multilevel_split(ntests,npass)

 ! split and merge with dust as a mixture and as particles
 if (use_dust) then
    call test_apr_dust(ntests,npass,onefluid=.true.)
    call test_apr_dust(ntests,npass,onefluid=.false.)
 endif

 if (id==master) write(*,"(/,a)") '<-- APR TEST COMPLETE'

end subroutine test_apr

!--------------------------------------------------------------------------
!+
!  check the tree-based nearest neighbour distance used to place split
!  particles against a brute force search (with periodic boundaries)
!+
!--------------------------------------------------------------------------
subroutine test_nearest_neighbour(ntests,npass)
 use unifdis,   only:set_unifdis
 use boundary,  only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use part,      only:npart,npartoftype,hfact,xyzh,vxyzu,init_part,massoftype,igas,&
                     set_particle_type
 use dim,       only:periodic
 use mpidomain, only:i_belong
 use deriv,     only:get_density_global
 use apr,       only:get_nearest_neighbour_distance
 integer, intent(inout) :: ntests,npass
 integer :: i,j,nlist,nfailed(1)
 integer, allocatable :: list(:)
 real,    allocatable :: rtree(:),rbrute(:)
 real    :: dx(3),r2

 if (id==master) write(*,"(/,a)") '--> checking nearest neighbour distances for splitting'
 call init_part()
 npart = 0
 call set_unifdis('random',id,master,xmin,xmax,ymin,ymax,zmin,zmax,dxbound/12.,&
                  hfact,npart,xyzh,periodic,verbose=.false.,mask=i_belong)
 npartoftype(igas) = npart
 massoftype(igas)  = dxbound*dybound*dzbound/npart
 do i=1,npart
    call set_particle_type(i,igas)
 enddo
 vxyzu(:,1:npart) = 0.
 call get_density_global(1)

 ! every other particle, so that the list is not just 1..n
 nlist = npart/2
 allocate(list(nlist),rtree(nlist),rbrute(nlist))
 do i=1,nlist
    list(i) = 2*i
 enddo
 call get_nearest_neighbour_distance(nlist,list,npart,xyzh,vxyzu,rtree)
 do i=1,nlist
    rbrute(i) = huge(0.)
    do j=1,npart
       if (j == list(i)) cycle
       dx = xyzh(1:3,list(i)) - xyzh(1:3,j)
       if (periodic) then
          dx(1) = dx(1) - dxbound*nint(dx(1)/dxbound)
          dx(2) = dx(2) - dybound*nint(dx(2)/dybound)
          dx(3) = dx(3) - dzbound*nint(dx(3)/dzbound)
       endif
       r2 = dot_product(dx,dx)
       rbrute(i) = min(rbrute(i),r2)
    enddo
    rbrute(i) = sqrt(rbrute(i))
 enddo
 call checkval(nlist,rtree,rbrute,1.e-14,nfailed(1),'nearest neighbour distance (tree vs brute force)')
 call update_test_scores(ntests,nfailed,npass)
 deallocate(list,rtree,rbrute)

end subroutine test_nearest_neighbour

!--------------------------------------------------------------------------
!+
!  refine by several levels in one call, as happens when restarting with
!  a larger apr_max; mass and momentum are conserved and each particle at
!  level l is one of 2**(l-1) children of an original particle
!+
!--------------------------------------------------------------------------
subroutine test_multilevel_split(ntests,npass)
 use unifdis,   only:set_unifdis
 use boundary,  only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use part,      only:npart,npartoftype,hfact,xyzh,vxyzu,fxyzu,init_part,massoftype,igas,&
                     set_particle_type,apr_level,radprop
 use dim,       only:periodic,maxtypes,maxdusttypes,do_radiation
 use mpidomain, only:i_belong
 use mpiutils,  only:reduceall_mpi
 use random,    only:ran2
 use apr,       only:init_apr,update_apr
 use utils_apr, only:apr_max_in,apr_centre
 integer, intent(inout) :: ntests,npass
 integer :: i,iseed,ierr,original_npart,nfailed(4)
 real    :: mtypein(maxtypes),mdustin(maxdusttypes),momin(3),hmean,nbase
 real    :: mtype(maxtypes),mdust(maxdusttypes),mom(3)

 if (id==master) write(*,"(/,a)") '--> refining by three levels at once'
 call init_part()
 npart = 0
 iseed = -2917
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,dxbound/16.,&
                  hfact,npart,xyzh,periodic,mask=i_belong)
 original_npart = npart
 npartoftype(igas) = npart
 massoftype(igas)  = dxbound*dybound*dzbound/reduceall_mpi('+',npart)
 do i=1,npart
    call set_particle_type(i,igas)
    vxyzu(1:3,i) = (/ran2(iseed),ran2(iseed),ran2(iseed)/)
 enddo
 if (do_radiation) radprop(4,:) = 23.0421

 ! same regions as the other tests, but with three additional levels
 call setup_apr_region_for_test()
 apr_max_in = 3
 call init_apr(apr_level,ierr)
 apr_centre(:,1:2) = 20.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 call get_apr_totals(npart,mtypein,mdustin,momin,hmean)

 apr_centre(1:2,1) = 0.25
 apr_centre(1:2,2) = -0.25
 apr_centre(3,1:2) = 0.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 call get_apr_totals(npart,mtype,mdust,mom,hmean)

 nbase = sum(0.5**(apr_level(1:npart)-1))
 nfailed(:) = 0
 call checkval(maxval(apr_level(1:npart)) == 4,.true.,nfailed(1),'particles refined by three levels')
 call checkval(nint(nbase),original_npart,0,nfailed(2),'number of particles before refinement')
 call checkval(maxtypes,mtype,mtypein,1.e-13,nfailed(3),'mass of each particle type')
 call checkval(3,mom,momin,1.e-13,nfailed(4),'linear momentum')
 call update_test_scores(ntests,nfailed,npass)

 ! restore a single level of refinement for subsequent tests
 apr_max_in = 1
 apr_level(1:npart) = 1

end subroutine test_multilevel_split

!--------------------------------------------------------------------------
!+
!  split and merge with dust, either as a one-fluid mixture with several
!  small grain species or as several types of dust particles; checks that
!  particle type counts, the mass of each type and of each dust species,
!  linear momentum and smoothing lengths survive a split and a merge
!+
!--------------------------------------------------------------------------
subroutine test_apr_dust(ntests,npass,onefluid)
 use unifdis,   only:set_unifdis
 use boundary,  only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use part,      only:npart,npartoftype,hfact,xyzh,init_part,massoftype,radprop,&
                     igas,idust,vxyzu,fxyzu,apr_level,maxvxyzu,set_particle_type,&
                     dustfrac,ndustsmall,ndustlarge,ndusttypes
 use dim,       only:periodic,do_radiation,maxtypes,maxdusttypes,maxdustsmall,maxdustlarge
 use options,   only:use_dustfrac
 use mpidomain, only:i_belong
 use mpiutils,  only:reduceall_mpi
 use apr,       only:update_apr
 use utils_apr, only:apr_centre
 use random,    only:ran2
 use physcon,   only:pi
 integer, intent(inout) :: ntests,npass
 logical, intent(in)    :: onefluid
 integer :: i,j,iseed,itype,npart_previous,original_npart,nfailed(6)
 integer :: npartoftype_in(maxtypes)
 real    :: psep,totmass,hmeanin,hmean,mom(3),momin(3)
 real    :: mtype(maxtypes),mtypein(maxtypes),mdust(maxdusttypes),mdustin(maxdusttypes)
 real, parameter :: tolm = 1.e-13, tolh = 1.e-13
 integer, parameter :: ndust_test = 3

 if (onefluid) then
    if (id==master) write(*,"(/,a)") '--> testing APR with dust as a mixture'
 else
    if (id==master) write(*,"(/,a)") '--> testing APR with dust as particles'
 endif
 iseed = -5135

 ! gas on a uniform lattice with random velocities
 call init_part()
 npartoftype(:) = 0
 npart = 0
 psep = dxbound/20.
 totmass = dxbound*dybound*dzbound
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact,npart,xyzh,periodic,mask=i_belong)
 npartoftype(igas) = npart
 massoftype(igas) = totmass/reduceall_mpi('+',npart)
 do i=1,npart
    call set_particle_type(i,igas)
 enddo

 if (onefluid) then
    ! several small grain species with different, spatially varying, dust fractions
    use_dustfrac = .true.
    ndustsmall = min(ndust_test,maxdustsmall)
    ndustlarge = 0
    do i=1,npart
       do j=1,ndustsmall
          dustfrac(j,i) = 0.01*j*(1. + 0.5*sin(2.*pi*(xyzh(1,i) - xmin)/dxbound))
       enddo
    enddo
 else
    ! several types of dust particles placed randomly
    use_dustfrac = .false.
    ndustsmall = 0
    ndustlarge = min(ndust_test,maxdustlarge)
    do j=1,ndustlarge
       npart_previous = npart
       call set_unifdis('random',id,master,xmin,xmax,ymin,ymax,zmin,zmax,2.*psep,&
                        hfact,npart,xyzh,periodic,inputiseed=-4231-j,verbose=.false.,mask=i_belong)
       itype = idust + j - 1
       do i=npart_previous+1,npart
          call set_particle_type(i,itype)
       enddo
       npartoftype(itype) = npart - npart_previous
       massoftype(itype) = 0.01*j*totmass/reduceall_mpi('+',npartoftype(itype))
    enddo
 endif
 ndusttypes = ndustsmall + ndustlarge
 original_npart = npart
 npartoftype_in = npartoftype

 do i=1,npart
    vxyzu(1:3,i) = (/ran2(iseed),ran2(iseed),ran2(iseed)/)
    if (maxvxyzu > 3) vxyzu(4,i) = ran2(iseed)**2
 enddo
 if (do_radiation) radprop(4,:) = 23.0421

 ! initialise APR with the regions away from the box
 call setup_apr_region_for_test()
 apr_centre(:,1:2) = 20.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 call get_apr_totals(npart,mtypein,mdustin,momin,hmeanin)

 ! split
 if (id==master) write(*,"(/,a)") '--> conducting a split'
 apr_centre(1:2,1) = 0.25
 apr_centre(1:2,2) = -0.25
 apr_centre(3,1:2) = 0.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 nfailed(:) = 0
 call check_apr_totals(npart,mtypein,mdustin,momin,tolm,nfailed(1:4))
 if (.not.onefluid) call checkval(sum(npartoftype(idust:idust+ndustlarge-1)) > &
    sum(npartoftype_in(idust:idust+ndustlarge-1)),.true.,nfailed(5),'dust particles were split')
 call update_test_scores(ntests,nfailed(1:5),npass)

 ! merge
 if (id==master) write(*,"(/,a)") '--> conducting a merge'
 apr_centre(:,1:2) = 20.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 nfailed(:) = 0
 call check_apr_totals(npart,mtypein,mdustin,momin,tolm,nfailed(1:4))
 call checkval(npart,original_npart,0,nfailed(5),'number of particles == original number')
 call get_apr_totals(npart,mtype,mdust,mom,hmean)
 call checkval(hmean,hmeanin,tolh,nfailed(6),'mean smoothing length after merge')
 call update_test_scores(ntests,nfailed(1:6),npass)

 ! restore defaults for subsequent tests
 use_dustfrac = .false.
 ndustsmall = 0
 ndustlarge = 0
 ndusttypes = 0

end subroutine test_apr_dust

!--------------------------------------------------------------------------
!+
!  check particle type counts, mass of each type, mass of each one-fluid
!  dust species and linear momentum against reference values
!+
!--------------------------------------------------------------------------
subroutine check_apr_totals(npart,mtypein,mdustin,momin,tol,nfailed)
 use part, only:npartoftype,iphase,iamtype,xyzh,isdead_or_accreted
 use dim,  only:maxtypes,maxdusttypes
 integer, intent(in)    :: npart
 real,    intent(in)    :: mtypein(maxtypes),mdustin(maxdusttypes),momin(3),tol
 integer, intent(inout) :: nfailed(4)
 integer :: i,ncount(maxtypes),itype
 real    :: mtype(maxtypes),mdust(maxdusttypes),mom(3),hmean

 ncount(:) = 0
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    itype = iamtype(iphase(i))
    ncount(itype) = ncount(itype) + 1
 enddo
 call checkval(all(npartoftype(1:maxtypes) == ncount),.true.,nfailed(1),'npartoftype matches particle types')

 call get_apr_totals(npart,mtype,mdust,mom,hmean)
 call checkval(maxtypes,mtype,mtypein,tol,nfailed(2),'mass of each particle type')
 call checkval(maxdusttypes,mdust,mdustin,tol,nfailed(3),'mass of each dust species')
 call checkval(3,mom,momin,tol,nfailed(4),'linear momentum')

end subroutine check_apr_totals

!--------------------------------------------------------------------------
!+
!  total mass of each particle type, mass of each one-fluid dust species,
!  linear momentum and mean smoothing length of all particles
!+
!--------------------------------------------------------------------------
subroutine get_apr_totals(npart,mtype,mdust,mom,hmean)
 use part,     only:xyzh,vxyzu,iphase,iamtype,dustfrac,ndustsmall,igas,&
                    get_pmass,isdead_or_accreted
 use options,  only:use_dustfrac
 use mpiutils, only:reduceall_mpi
 integer, intent(in)  :: npart
 real,    intent(out) :: mtype(:),mdust(:),mom(3),hmean
 integer :: i,itype,n
 real    :: pmassi

 mtype(:) = 0.
 mdust(:) = 0.
 mom(:)   = 0.
 hmean    = 0.
 n = 0
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    itype  = iamtype(iphase(i))
    pmassi = get_pmass(i,.false.)
    mtype(itype) = mtype(itype) + pmassi
    if (use_dustfrac .and. itype == igas) then
       mdust(1:ndustsmall) = mdust(1:ndustsmall) + pmassi*dustfrac(1:ndustsmall,i)
    endif
    mom(:) = mom(:) + pmassi*vxyzu(1:3,i)
    hmean  = hmean + xyzh(4,i)
    n = n + 1
 enddo
 mtype = reduceall_mpi('+',mtype)
 mdust = reduceall_mpi('+',mdust)
 mom   = reduceall_mpi('+',mom)
 hmean = reduceall_mpi('+',hmean)/reduceall_mpi('+',n)

end subroutine get_apr_totals

!--------------------------------------------
!+
!  Set up an APR region that is used in other tests
!+
!--------------------------------------------
subroutine setup_apr_region_for_test()
 use apr,        only:init_apr,update_apr
 use utils_apr,  only:apr_type,apr_rad,apr_max_in,ref_dir,ntrack
 use part,       only:apr_level
 integer :: ierr

 if (id==master) write(*,"(/,a)") '--> adding an apr region'

 ! set parameters for the region
 apr_max_in  =   1     ! number of additional refinement levels (3 -> 2x resolution)
 ref_dir     =   1     ! increase (1) or decrease (-1) resolution
 apr_type    =  -1     ! choose this so you get the default option which is reserved for the test suite
 apr_rad     =   0.2   ! radius of innermost region
 ntrack      =   2     ! number of regions to track

 ! initialise
 call init_apr(apr_level,ierr)

end subroutine setup_apr_region_for_test

end module testapr
