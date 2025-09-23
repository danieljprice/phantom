!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Dependencies: apr, boundary, dim, energies, io, mpidomain, mpiutils,
!   part, random, testutils, unifdis, utils_apr
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
 use part,         only:isetphase,igas,iphase,vxyzu,fxyzu,apr_level,maxvxyzu,iphase_soa
 use mpidomain,    only:i_belong
 use mpiutils,     only:reduceall_mpi
 use dim,          only:periodic,use_apr,do_radiation
 use apr,          only:update_apr
 use utils_apr,    only:apr_centre
 use energies,     only:compute_energies,angtot,etot,totmom,ekin,etherm
 use random,       only:ran2
 integer, intent(inout) :: ntests,npass
 real :: psep,rhozero,time,totmass, etotin, totmomin
 real :: angtotin, ekinin, ethermin
 real :: tolang, tolen, tolmom
 integer :: original_npart,splitted,nfailed(7),i,iseed

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
 massoftype(igas) = totmass/reduceall_mpi('+',npart)
 iphase(1:npart) = isetphase(igas,iactive=.true.)

 ! this is to prevent a (reasonable) problem when running this test with DEBUG=yes and radiation
 if (do_radiation) then
    radprop(4,:) = 23.0421 ! just some inconsequential number
 endif

 ! and this is because the iphase_soa arrays are not initialised in the test suite
 iphase_soa(:) = 1

 ! Set some random velocities
 do i=1,npart
    vxyzu(1:3,i) = (/ran2(iseed),ran2(iseed),ran2(iseed)/)
    if (maxvxyzu > 3) vxyzu(4,i) = ran2(iseed)**2
 enddo

 ! Initialise APR
 call setup_apr_region_for_test()
 apr_centre(:,1:2) = 20. ! just moves the APR region away from the box so you don't have any split or merge
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)


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

 if (id==master) write(*,"(/,a)") '<-- APR TEST COMPLETE'

end subroutine test_apr

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
