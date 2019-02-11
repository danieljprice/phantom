!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup for star cluster formation calculations
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, datafiles, eos, io, part, physcon, ptmass,
!    random, setup_params, setvfield, timestep, units, velfield
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!
!  Sets up a star cluster formation calculation
!  following Bate, Bonnell & Bromm (2003). Requires
!  pre-calculated velocity cubes.
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,  only:pi,solarm,pc,years,kboltz,mass_proton_cgs,au
 use velfield, only:set_velfield_from_cubes
 use setup_params, only:rmax,rhozero
 use random,       only:ran2
 use part,         only:igas,set_particle_type
 use io,           only:fatal,master
 use units,        only:umass,udist,utime,set_units
 use setvfield,    only:normalise_vfield
 use timestep,     only:dtmax,tmax
 use centreofmass, only:reset_centreofmass
 use ptmass,       only:h_acc,r_crit,rho_crit_cgs,icreate_sinks
 use datafiles,    only:find_phantom_datafile
 use eos,          only:ieos
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: r2,totmass,xi,yi,zi,epotgrav,t_ff, T,mu
 integer :: ipart,iseed,ierr
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'
 character(len=120) :: filex,filey,filez

 call set_units(dist=0.1d0*pc,mass=solarm,G=1.)

 npartoftype(:) = 0
 npartoftype(1) = size(xyzh(1,:))
 npart = npartoftype(1)
 gamma = 1.0
 T     = 10.0           ! Temperature in Kelvin
 mu    = 2.46           ! Mean molecular weight
 polyk = kboltz*T/(mu*mass_proton_cgs)*(utime/udist)**2

 rmax = 0.1875*(pc/udist)
 r2 = rmax*rmax
 if (id==master) write(*,"(1x,a)") 'Cluster formation setup: '

 totmass = 50.*(solarm/umass)
 massoftype(1) = totmass/real(npart)
 rhozero = totmass/(4./3.*pi*rmax*r2)
 write(*,"(1x,a,f6.3,a)")   '       Rcloud = ',rmax*(udist/pc),' pc'
 write(*,"(1x,a,f6.2,a)")   '       Mcloud = ',totmass*(umass/solarm),' Msun'
 write(*,"(1x,a,es10.3,a)") ' Mean density = ',rhozero*umass/udist**3,' g/cm^3'
 write(*,"(1x,a,es10.3,a,es10.3,a)") 'Particle mass = ',massoftype(1)*(umass/solarm),' Msun'
 hfact = 1.2
 t_ff = sqrt(3.*pi/(32.*rhozero))
 write(*,"(1x,a,es10.3,a,e10.3,a)") 'Freefall time = ',t_ff*(utime/years),' years (',t_ff,' in code units)'
 write(*,"(1x,a,es10.3,a)") '  Sound speed = ',sqrt(polyk)*(udist/utime), ' cm/s'

 ipart = 0
 iseed = -2485
 do while (ipart < npart)
    xi = 2.*rmax*(ran2(iseed)-0.5)
    yi = 2.*rmax*(ran2(iseed)-0.5)
    zi = 2.*rmax*(ran2(iseed)-0.5)
    if ((xi*xi + yi*yi + zi*zi) < r2) then
       ipart = ipart + 1
       xyzh(1,ipart) = xi
       xyzh(2,ipart) = yi
       xyzh(3,ipart) = zi
       xyzh(4,ipart) = hfact*(massoftype(1)/rhozero)**(1./3.)
       call set_particle_type(ipart,igas)
    endif
 enddo
 write(*,"(1x,a)") 'Setting up velocity field on the particles...'
 vxyzu(:,:) = 0.

 filex = find_phantom_datafile(filevx,'velfield')
 filey = find_phantom_datafile(filevy,'velfield')
 filez = find_phantom_datafile(filevz,'velfield')

 call set_velfield_from_cubes(xyzh,vxyzu,npartoftype(igas),filex,filey,filez,&
                              1.,rmax,.false.,ierr)
 if (ierr /= 0) call fatal('setup','error setting up velocity field')

 epotgrav = 3./5.*totmass**2/rmax
 call normalise_vfield(npart,vxyzu,ierr,ke=epotgrav)
 if (ierr /= 0) call fatal('setup','error normalising velocity field')

 ! Setting the centre of mass of the cloud to be zero
 call reset_centreofmass(npart,xyzh,vxyzu)

 ! set options for input file
 tmax  = 2.*t_ff
 dtmax = 0.002*t_ff
 h_acc = 5.*au/udist
 r_crit= 2.*h_acc
 ieos = 8
 icreate_sinks = 1
 rho_crit_cgs  = 1.e-10

end subroutine setpart

end module setup
