!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    gastemp       -- temperature at injection point in K
!    mdot          -- mass injection rate in grams/second
!    npartperorbit -- particle injection rate in particles/binary orbit
!
!  DEPENDENCIES: infile_utils, io, part, partinject, physcon, random, units
!+
!--------------------------------------------------------------------------
module inject
 implicit none
 character(len=*), parameter, public :: inject_type = 'asteroidwind'

 public :: inject_particles,write_options_inject,read_options_inject

 private

 real :: mdot          = 5.e8       ! mass injection rate in grams/second
 real :: npartperorbit = 100.       ! particle injection rate in particles per orbit
 real :: gastemp       = 3000.      ! gas temperature

contains

subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use io,        only:fatal
 use part,      only:nptmass,massoftype,igas,hfact,ihsoft
 use partinject,only:add_or_update_particle
 use physcon,   only:pi,twopi,gg,kboltz,mass_proton_cgs
 use random,    only:ran2
 use units,     only:udist, umass, utime
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    dimension(3)  :: xyz,vxyz,r1,r2,v2,vhat
 integer :: i,ipart,npinject,seed
 real    :: dmdt,dndt,rasteroid,h,u,speed
 real    :: m1,m2,mu,period,r,q
 real    :: phi,theta

 if (nptmass < 2) call fatal('inject_rochelobe','not enough point masses for asteroid wind injection')
 if (nptmass > 2) call fatal('inject_rochelobe','too many point masses for asteroid wind injection')

 r1        = xyzmh_ptmass(1:3,1)
 r2        = xyzmh_ptmass(1:3,2)
 rasteroid = xyzmh_ptmass(ihsoft,2)
 m1        = xyzmh_ptmass(4,1)
 m2        = xyzmh_ptmass(4,2)
 v2        = vxyz_ptmass(1:3,2)

 speed     = sqrt(dot_product(v2,v2))
 vhat      = v2/speed

 r         = sqrt(dot_product(r1-r2,r1-r2))
 q         = m2/m1
 mu        = 1./(1 + q)
 period    = twopi*sqrt((r*udist)**3/(gg*(m1+m2)*umass))       ! period of orbit in code units

 dmdt      = mdot/(umass/utime)                                ! convert grams/sec to code units
 dndt      = npartperorbit*utime/period                        ! convert particles per orbit into code units

!
!-- Mass of gas particles is set by mass accretion rate and particle injection rate
!
 massoftype(igas) = dmdt/dndt

!
!-- How many particles do we need to inject?
!   (Seems to need at least eight gas particles to not crash) <-- This statement may or may not be true...
!
 if(npartoftype(igas)<8) then
    npinject = 8-npartoftype(igas)
 else
    npinject = max(0, int(0.5 + (time*dmdt/massoftype(igas)) - npartoftype(igas) ))
 endif

!
!-- Randomly inject particles around the asteroids outer 'radius'
!
 do i=1,npinject
    phi       = ran2(seed)*twopi
    theta     = ran2(seed)*pi
    xyz       = r2 + (/rasteroid*cos(phi)*sin(theta),rasteroid*sin(phi)*sin(theta),rasteroid*cos(theta)/)
    vxyz      = 0.999*speed*vhat
    u         = 3.*(kboltz*gastemp/(mu*mass_proton_cgs))/2. * (utime/udist)**2
    h         = hfact*(rasteroid/2.)
    ipart     = npart + 1
    call add_or_update_particle(igas,xyz,vxyz,h,u,ipart,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(mdot         ,'mdot'         ,'mass injection rate in grams/second'              ,iunit)
 call write_inopt(npartperorbit,'npartperorbit','particle injection rate in particles/binary orbit',iunit)
 call write_inopt(gastemp      ,'gastemp'      ,'temperature at injection point in K'              ,iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 select case(trim(name))
 case('mdot')
    read(valstring,*,iostat=ierr) mdot
    ngot = ngot + 1
    if (mdot  <  0.) call fatal(label,'mdot < 0 in input options')
case('npartperorbit')
    read(valstring,*,iostat=ierr) npartperorbit
    ngot = ngot + 1
    if (npartperorbit < 0.) call fatal(label,'npartperorbit < 0 in input options')
 case('gastemp')
    read(valstring,*,iostat=ierr) gastemp
    ngot = ngot + 1
    if (gastemp<= 0.) call fatal(label,'gastemp 0 or negative in input options')
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
