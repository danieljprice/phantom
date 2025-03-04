!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setsolarsystem
!
! Setup asteroid orbits using data from the IAU Minor Planet Center
!
! :References: https://minorplanetcenter.net/data
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: datautils, ephemeris, fileutils, mpc, part, physcon,
!   setbinary, units
!
 implicit none
 public :: set_minor_planets
 public :: add_body,add_sun_and_planets,add_dwarf_planets

 private

contains
!----------------------------------------------------------------
!+
!  set up asteroids as a population of dust particles
!+
!----------------------------------------------------------------
subroutine set_minor_planets(npart,npartoftype,massoftype,xyzh,vxyzu,mtot,itype,sample_orbits)
 use part,      only:set_particle_type,nsinkproperties
 use physcon,   only:au
 use units,     only:udist
 use datautils, only:find_datafile
 use mpc,       only:read_mpc,mpc_entry
 use setbinary, only:set_binary
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 real,    intent(in)    :: mtot
 integer, intent(in)    :: itype
 logical, intent(in), optional :: sample_orbits
 integer :: nbodies,nsample,n,i,j,nptmass,ierr
 real    :: xyzmh_ptmass(nsinkproperties,2),vxyz_ptmass(3,2),hpart
 character(len=120) :: filename
 type(mpc_entry), allocatable :: dat(:)

 nsample = 1
 if (present(sample_orbits)) then
    ! can place many particles evenly sampling the orbit if desired
    if (sample_orbits) nsample = 100
 endif

 nbodies = 0
 filename = find_datafile('Distant.txt',url='https://www.minorplanetcenter.net/iau/MPCORB/')
 call read_mpc(filename,nbodies,dat=dat)
 print "(a,i0,a)",' read orbital data for ',nbodies,' minor planets'

 n = 0
 hpart  = 10.*au/udist/nsample

 do i=1,nbodies
    !
    ! for each solar system object get the xyz positions from the orbital parameters
    !
    !print*,i,'aeiOwM=',dat(i)%a,dat(i)%ecc,dat(i)%inc,dat(i)%O,dat(i)%w,dat(i)%M
    do j=1,nsample
       nptmass = 0
       n = n + 1
       if (nsample==1) then
          call set_binary(mtot,epsilon(0.),dat(i)%a,dat(i)%ecc,0.02,1.e-15,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,incl=dat(i)%inc,&
                    arg_peri=dat(i)%w,posang_ascnode=dat(i)%O,&
                    mean_anomaly=dat(i)%M,verbose=.false.)
       else
          call set_binary(mtot,epsilon(0.),dat(i)%a,dat(i)%ecc,0.02,1.e-15,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,incl=dat(i)%inc,&
                    arg_peri=dat(i)%w,posang_ascnode=dat(i)%O,&
                    f=360.*(n-1)/nsample,verbose=.false.)
       endif
       !
       ! now delete the point masses but set a dust particle as the secondary
       !
       xyzh(1:3,n)  = xyzmh_ptmass(1:3,2)
       xyzh(4,n)    = hpart  ! give a random length scale as the smoothing length
       vxyzu(1:3,n) = vxyz_ptmass(1:3,2)
       call set_particle_type(n,itype)
    enddo

 enddo
 !
 ! set mass of all the minor bodies equal
 !
 print "(a)"
 npart = npart + nbodies*nsample
 npartoftype(itype) = npartoftype(itype) + nbodies*nsample
 massoftype(itype) = 1.e-20

end subroutine set_minor_planets

!----------------------------------------------------------------
!+
!  setup the solar system planets by querying their ephemeris
!  from the JPL server
!+
!----------------------------------------------------------------
subroutine add_sun_and_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)
 integer, intent(inout) :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(in)    :: mtot
 character(len=*), intent(in), optional :: epoch
 integer,          parameter :: nbodies = 10
 character(len=*), parameter :: planet_name(nbodies) = &
     (/'sun    ',&
       'mercury', &
       'venus  ', &
       'earth  ', &
       'moon   ', &
       'mars   ', &
       'jupiter', &
       'saturn ', &
       'uranus ', &
       'neptune'/)
 integer :: i

 do i=1,nbodies
    if (present(epoch)) then
       call add_body(planet_name(i),nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch=epoch)
    else
       call add_body(planet_name(i),nptmass,xyzmh_ptmass,vxyz_ptmass,mtot)
    endif
 enddo

end subroutine add_sun_and_planets

!----------------------------------------------------------------
!+
!  setup the solar system dwarf planets by querying their ephemeris
!  from the JPL server
!+
!----------------------------------------------------------------
subroutine add_dwarf_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)
 integer, intent(inout) :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(in)    :: mtot
 character(len=*), intent(in), optional :: epoch
 integer,          parameter :: nbodies = 5
 character(len=*), parameter :: planet_name(nbodies) = &
       (/'ceres   ',&
         'pluto   ', &
         'eris    ', &
         'makemake', &
         'haumea  '/)
 integer :: i

 do i=1,nbodies
    if (present(epoch)) then
       call add_body(planet_name(i),nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch=epoch)
    else
       call add_body(planet_name(i),nptmass,xyzmh_ptmass,vxyz_ptmass,mtot)
    endif
 enddo

end subroutine add_dwarf_planets

!----------------------------------------------------------------
!+
!  setup a body in the solar system by querying their ephemeris
!  from the JPL server
!+
!----------------------------------------------------------------
subroutine add_body(body_name,nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)
 use ephemeris, only:get_ephemeris,nelem
 use units,     only:umass,udist
 use physcon,   only:gg,km,solarm,solarr,earthm,au
 use setbinary, only:set_binary
 use fileutils, only:lcase
 character(len=*), intent(in)    :: body_name
 integer,          intent(inout) :: nptmass
 real,             intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,             intent(in)    :: mtot
 character(len=*), intent(in), optional :: epoch
 real    :: elems(nelem)
 real    :: xyz_tmp(size(xyzmh_ptmass(:,1)),2),vxyz_tmp(3,2),gm_cgs
 real    :: mbody,rbody,a,e,inc,O,w,f
 integer :: ierr,ntmp
 logical :: got_elem(nelem)

 if (trim(adjustl(lcase(body_name)))=='sun') then
    !
    ! add the Sun, ideally would add here its motion
    ! around the solar system barycentre
    !
    nptmass = nptmass + 1
    xyzmh_ptmass(:,nptmass) = 0.
    xyzmh_ptmass(4,nptmass) = solarm/umass
    xyzmh_ptmass(5,nptmass) = solarr/udist
    vxyz_ptmass(:,nptmass) = 0.
    return
 endif

 if (present(epoch)) then
    elems = get_ephemeris(lcase(body_name),got_elem,ierr,epoch=epoch)
 else
    elems = get_ephemeris(lcase(body_name),got_elem,ierr)
 endif
 if (ierr > 0 .or. .not.any(got_elem)) then
    print "(a)",' ERROR: could not read ephemeris data for '//body_name
    return ! skip if error reading ephemeris file
 endif
 a   = elems(1)*km/udist
 e   = elems(2)
 inc = elems(3)
 O   = elems(4)
 w   = elems(5)
 f   = elems(6)
 gm_cgs = elems(7)*km**3
 mbody  = (gm_cgs/gg)/umass
 rbody = elems(8)*km/udist
 print "(1x,a,1pg10.3)",   '           m/msun = ',mbody*umass/solarm
 print "(1x,a,1pg10.3)",   '         m/mearth = ',mbody*umass/earthm
 print "(1x,a,1pg10.4,a)", '                a = ',a*udist/au,' au'
 print "(1x,a,1pg10.3,a)", '           radius = ',elems(8),' km'
 print "(1x,a,1pg10.3,a)", '          density = ',elems(9),' g/cm^3'
 ntmp = 0
 call set_binary(mtot,0.,a,e,0.01,rbody,&
      xyz_tmp,vxyz_tmp,ntmp,ierr,incl=inc,&
      arg_peri=w,posang_ascnode=O,f=f,verbose=.false.)

 ! add a point mass for each body
 nptmass = nptmass + 1
 xyzmh_ptmass(1:3,nptmass) = xyz_tmp(1:3,2)
 xyzmh_ptmass(4,nptmass) = mbody
 xyzmh_ptmass(5,nptmass) = rbody
 vxyz_ptmass(1:3,nptmass) = vxyz_tmp(1:3,2)
 print "(1x,a,3(1x,1pg10.3))",' x = ',xyz_tmp(1:3,2)
 print "(1x,a,3(1x,1pg10.3))",' v = ',vxyz_tmp(1:3,2)

end subroutine add_body

end module setsolarsystem
