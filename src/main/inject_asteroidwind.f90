!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!    mdot          -- mass injection rate in grams/second
!    npartperorbit -- particle injection rate in particles/binary orbit
!    vlag          -- percentage lag in velocity of wind
!
!  DEPENDENCIES: infile_utils, io, part, partinject, physcon, random, units
!+
!--------------------------------------------------------------------------
module inject
 implicit none
 character(len=*), parameter, public :: inject_type = 'asteroidwind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 private

 real :: mdot          = 5.e8       ! mass injection rate in grams/second
 real :: npartperorbit = 100.       ! particle injection rate in particles per orbit
 real :: vlag          = 0.1        ! percentage lag in velocity of wind

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Inject particles
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,        only:fatal
 use part,      only:nptmass,massoftype,igas,hfact,ihsoft
 use partinject,only:add_or_update_particle
 use physcon,   only:pi,twopi,gg,kboltz,mass_proton_cgs
 use random,    only:ran2
 use units,     only:udist, umass, utime
 use options,   only:iexternalforce
 use externalforces,only:mass1
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real,    dimension(3)  :: xyz,vxyz,r1,r2,v2,vhat
 integer :: i,ipart,npinject,seed,pt
 real    :: dmdt,dndt,rasteroid,h,u,speed
 real    :: m1,m2,mu,period,r,q,semia,spec_energy
 real    :: phi,theta,func,rp,ra,ecc

 if (nptmass < 2 .and. iexternalforce == 0) call fatal('inject_asteroidwind','not enough point masses for asteroid wind injection')
 if (nptmass > 2) call fatal('inject_asteroidwind','too many point masses for asteroid wind injection')

 if (nptmass == 2) then
    pt = 2
    r1 = xyzmh_ptmass(1:3,1)
    m1 = xyzmh_ptmass(4,1)
 else
    pt = 1
    r1 = 0.
    m1 = mass1
 endif

 r2        = xyzmh_ptmass(1:3,pt)
 rasteroid = xyzmh_ptmass(ihsoft,pt)
 m2        = xyzmh_ptmass(4,pt)
 v2        = vxyz_ptmass(1:3,pt)

 speed     = sqrt(dot_product(v2,v2))
 vhat      = v2/speed

 r         = sqrt(dot_product(r1-r2,r1-r2))
 q         = m2/m1
 mu        = 1./(1 + q)
 spec_energy = 0.5*speed**2 - (1.0*m1/r)
 if (iexternalforce == 11) spec_energy = spec_energy - (3.*m1/(r**2))
 semia     = -m1/(2.0*spec_energy)

 rp = 0.42
 ra = 0.98
 ecc = 0.4

 func = (ra*rp/(r**2)) ! - ((1.-ecc)/(1.+ecc)))              ! function to scale dn/dt with r^2
							       ! rp*ra instead of semia**2 is more accurate
							       ! but not by much

 period    = twopi*sqrt((semia*udist)**3/(gg*(m1+m2)*umass))   ! period of orbit
 period    = period/utime                                      ! in code units

 dmdt      = mdot/(umass/utime)                                ! convert grams/sec to code units
 dndt      = npartperorbit/period                        ! convert particles per orbit into code units

!
!-- Mass of gas particles is set by mass accretion rate and particle injection rate
!
 massoftype(igas) = dmdt/dndt

!
!-- How many particles do we need to inject?
!   (Seems to need at least eight gas particles to not crash) <-- This statement may or may not be true...
!
 if (npartoftype(igas)<8) then
    npinject = 8-npartoftype(igas)
 else
!    npinject = max(0, int(0.5 + (time*dndt*func) - npartoftype(igas) ))
   npinject = max(0, int(0.5 + (mod(time,period)/period*dndt*func) ))
   print*,time,npinject
 endif

!
!-- Randomly inject particles around the asteroids outer 'radius'
!
 do i=1,npinject
    phi       = ran2(seed)*twopi
    theta     = ran2(seed)*pi
    xyz       = r2 + (/rasteroid*cos(phi)*sin(theta),rasteroid*sin(phi)*sin(theta),rasteroid*cos(theta)/)
    vxyz      = (1.-vlag/100)*speed*vhat
    u         = 0. ! setup is isothermal so utherm is not stored
    h         = hfact*(rasteroid/2.)
    ipart     = npart + 1
    call add_or_update_particle(igas,xyz,vxyz,h,u,ipart,npart,npartoftype,xyzh,vxyzu)
 enddo

! if (npinject > 0) print*,time,r,func,(mod(time,period)/period)

 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)

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
 call write_inopt(vlag         ,'vlag'         ,'percentage lag in velocity of wind'               ,iunit)

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
 case('vlag')
    read(valstring,*,iostat=ierr) vlag
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
