!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!  Handles injection of particles in Keplerian shearing flow
!
!  Here is the "philosophy" of this injection routine:
!  Injected particles are queued up at two "injection zones" on either
!  (azimuthal) side of the simulation domain
!
!  All particles are initially boundary types
!
!  These particles will move at a fixed velocity (plus velocity change due to external forces)
!  Once inside the domain, they become "live" gas particles
!
!  When particles leave the simulation domain and enter the opposite injection zone (the exit zone),
!  they become boundary particles again
!  When they leave the injection zone, they are killed
!
!  Some boundary particles are injected along both radial boundaries outside the domain
!  to maintain the correct pressure on particles inside
!
!  The setdisc module is used to create the Keplerian flow for injection
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    HoverR      -- disc aspect ratio at inner sector radius
!    R_in        -- inner total disc radius
!    R_out       -- outer total disc radius
!    Rsect_in    -- inner sector radius (inner injection radius)
!    Rsect_out   -- outer sector radius (outer injection radius)
!    disc_mass   -- total disc mass
!    dr_bound    -- Radial boundary thickness
!    object_mass -- mass of the central object
!    p_index     -- radial surface density profile powerlaw
!    phi_inject  -- azimuthal range of injection zone
!    phimax      -- maximum azimuthal extent (-phimax,phimax)
!    q_index     -- radial sound speed profile powerlaw
!
!  DEPENDENCIES: eos, infile_utils, io, part, partinject, physcon, setdisc
!+
!--------------------------------------------------------------------------
module inject
 implicit none
 character(len=*), parameter, public :: inject_type = 'keplerianshear'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject
 public :: set_injection_parameters

 type injectparams
    real, public :: R_in, R_out, Rsect_in, Rsect_out, width, R_mid, dr_bound, phi_inject
    real, public :: phimax, p_index,q_index, HoverR, object_mass, disc_mass
 end type

 type(injectparams), public :: injp

 integer :: nqueuecrit, ngas_initial, ninjectmax
 logical :: firstrun

 private

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
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use part,      only:igas,iboundary, massoftype
 use physcon,   only:Rg,gg,pi
 use eos,       only:gamma
 use io,        only:master
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 real, parameter :: mu = 1.26 ! Used in Bowen (1988)

 integer :: ninject,nqueue, nkill, nboundary, ndomain, injected, nexit
 logical :: replenish


 ! Subroutine begins:
 massoftype(iboundary) = massoftype(igas)

 !--------------
 ! 1. Determine the state of all particles in the simulation: live, boundary, or dead
 !--------------

 call determine_particle_status(nqueue, nkill, nboundary, ndomain, nexit)

 npartoftype(igas) = ndomain
 npartoftype(iboundary) = nboundary +nkill

 !--------------
 ! 2. Replenish the injection zone if necessary
 !--------------

 firstrun = .not.(time>1.0e-9)
 replenish = .false.

! If this is the first run, injection zone empty and needs to be filled

 if(firstrun) then
    ngas_initial = npartoftype(igas)
    replenish = .true.
    ninjectmax = int(ngas_initial*injp%phi_inject/injp%phimax)

    nqueuecrit = -1
    print*, 'First timestep: initialising queue'
 endif

 if(nqueue< nqueuecrit) then
    print*, 'Queue nearly empty: replenishing'
    replenish=.true. ! If a small number of particles in queue, replenishment required
 endif

 if(npartoftype(igas) < int(0.9*ngas_initial) .and. nqueue < ninjectmax) then
    print*, 'Simulation domain depleted: replenishing'
    replenish = .true.
 endif

 if(injp%phimax>0.9*pi) replenish=.false.

! Number of particles to inject must be proportional to number in simulation



 if(firstrun) then
    ninject = ninjectmax
 else
    ninject = ngas_initial - npartoftype(igas)
 endif
 gamma = 5.0/3.0

 print*, 'Number of particles killed during injection: ', nkill
 print*, 'Number of particles that have left the domain: ',nexit
 print*, 'Number of boundary particles in the simulation: ', nboundary
 print*, 'Number of particles queued up for injection: ', nqueue, nqueuecrit, ninject, ninjectmax
 print*, 'Number of particles in the simulation domain: ', ndomain, npartoftype(igas), ngas_initial



 injected = 0

 if(replenish) then
    call replenish_injection_zone(ninject, time,dtlast, injected)
    print*, injected, ' particles injected '
    if(firstrun) then
       nqueuecrit = int(ninject/2)
       firstrun = .false.
    endif
 else
    print*, 'No injection on this timestep'
 endif

 write(87,*) time, dtlast, ninject, ndomain, nboundary, nkill, nqueue, nqueuecrit, injected
 print*, npart, ' particles in total'
 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)

end subroutine

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(injp%R_in,'R_in','inner total disc radius',iunit)
 call write_inopt(injp%R_out,'R_out','outer total disc radius',iunit)
 call write_inopt(injp%Rsect_in,'Rsect_in','inner sector radius (inner injection radius)',iunit)
 call write_inopt(injp%Rsect_out,'Rsect_out','outer sector radius (outer injection radius)',iunit)
 call write_inopt(injp%dr_bound,'dr_bound','Radial boundary thickness',iunit)
 call write_inopt(injp%phimax, 'phimax', 'maximum azimuthal extent (-phimax,phimax)',iunit)
 call write_inopt(injp%phi_inject, 'phi_inject', 'azimuthal range of injection zone',iunit)
 call write_inopt(injp%p_index, 'p_index', 'radial surface density profile powerlaw',iunit)
 call write_inopt(injp%q_index, 'q_index', 'radial sound speed profile powerlaw',iunit)
 call write_inopt(injp%HoverR, 'HoverR', 'disc aspect ratio at inner sector radius', iunit)
 call write_inopt(injp%disc_mass,'disc_mass', 'total disc mass', iunit)
 call write_inopt(injp%object_mass,'object_mass', 'mass of the central object', iunit)
end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('R_in')
    read(valstring,*,iostat=ierr) injp%R_in
    ngot = ngot + 1
 case('R_out')
    read(valstring,*,iostat=ierr) injp%R_out
    ngot = ngot + 1
 case('Rsect_in')
    read(valstring,*,iostat=ierr) injp%Rsect_in
    ngot = ngot + 1
 case('Rsect_out')
    read(valstring,*,iostat=ierr) injp%Rsect_out
    ngot = ngot + 1
 case('dr_bound')
    read(valstring,*,iostat=ierr) injp%dr_bound
    ngot = ngot+1
 case('phimax')
    read(valstring,*,iostat=ierr) injp%phimax
    injp%phimax = injp%phimax*3.14159264/180.0
    ngot = ngot + 1
 case('phi_inject')
    read(valstring,*,iostat=ierr) injp%phi_inject
    injp%phi_inject = injp%phi_inject*3.141592654/180.0
    ngot = ngot+1
 case('p_index')
    read(valstring,*,iostat=ierr) injp%p_index
    ngot = ngot +1
 case('q_index')
    read(valstring,*,iostat=ierr) injp%q_index
    ngot = ngot +1
 case('HoverR')
    read(valstring,*,iostat=ierr) injp%HoverR
    ngot = ngot +1
 case('disc_mass')
    read(valstring,*,iostat=ierr) injp%disc_mass
    ngot = ngot+1
 case('object_mass')
    read(valstring,*,iostat=ierr) injp%object_mass
    ngot = ngot+1
 end select


 igotall = (ngot >= 11)
 injp%width = injp%Rsect_out - injp%Rsect_in
 injp%R_mid = injp%width/2.0 + injp%Rsect_in


end subroutine

subroutine set_injection_parameters(R_in, R_out, Rsect_in,Rsect_out,dr_bound,&
 phimax,phi_inject,p_index,q_index,HoverR,disc_mass,object_mass)

 real, intent(in) :: R_in, R_out, Rsect_in,Rsect_out, dr_bound
 real, intent(in) :: phimax,p_index,q_index,HoverR, phi_inject
 real, intent(in) :: disc_mass,object_mass

 injp%R_in = R_in
 injp%R_out = R_out
 injp%Rsect_in = Rsect_in
 injp%Rsect_out = Rsect_out
 injp%width = injp%Rsect_out - injp%Rsect_in
 injp%R_mid = injp%width/2.0 + injp%Rsect_in

 injp%phimax = phimax
 injp%p_index = p_index
 injp%q_index = q_index
 injp%HoverR = HoverR
 injp%disc_mass = disc_mass
 injp%object_mass = object_mass

 injp%dr_bound = dr_bound
 injp%phi_inject = phi_inject

 return
end subroutine set_injection_parameters


subroutine determine_particle_status(nqueue, nkill, nboundary, ndomain, nexit)
 use part, only : igas, iboundary, npart, xyzh, kill_particle, set_particle_type

 integer, intent(inout) :: nqueue, nkill,nboundary,ndomain, nexit

 integer :: i
 real :: rpart,phipart, midsep
 real :: rsign, phisign, annulus_halfwidth

 logical :: radial_boundary, buffer_zone,injection_zone
 logical :: exit_zone, dead_zone, simulation_domain

 annulus_halfwidth = 0.5*(injp%Rsect_out-injp%Rsect_in)
 nkill = 0
 nqueue = 0
 nboundary = 0
 ndomain = 0
 nexit = 0

 do i =1,npart

    call calc_polar_coordinates(rpart,phipart,xyzh(1,i), xyzh(2,i))

    rsign = rpart-injp%R_mid
    midsep = abs(rsign)

    rsign = rsign/midsep
    phisign = phipart/abs(phipart)

! Particles outside the simulation and boundary/buffer zones should be killed
    dead_zone = (abs(phipart) > injp%phimax+injp%phi_inject) .or. &
        (midsep > annulus_halfwidth + injp%dr_bound)


! Maintain boundary particles in the radial limits to keep pressure correct
    radial_boundary = (midsep > annulus_halfwidth) .and.  &
        (midsep <= annulus_halfwidth + injp%dr_bound) .and. .not.(dead_zone)

! Two buffer zones, either side of simulation domain in phi
! They contain regions where particles are injected, but injection only occurs in half of both zones
    buffer_zone = (abs(phipart)>injp%phimax) .and. (abs(phipart)<=injp%phimax+injp%phi_inject) &
            .and..not.(dead_zone).and..not.(radial_boundary)

! Can be in a buffer zone without being in an injection zone!
    injection_zone = buffer_zone .and. (rsign*phisign>0.0)
    exit_zone = buffer_zone .and.(rsign*phisign<0.0)

! Or particles can be in the simulation domain!
    simulation_domain = (midsep <=annulus_halfwidth) .and. (abs(phipart) <=injp%phimax)


! Set dead particles
    if(dead_zone) then
       call kill_particle(i)
       nkill = nkill+1
    endif

! Set boundary particles
    if(buffer_zone.or.radial_boundary) then
       call set_particle_type(i,iboundary)
       nboundary = nboundary+1
       if(injection_zone) nqueue = nqueue+1
       if(exit_zone) nexit = nexit+1
    endif

! Set living particles
    if(simulation_domain) then
       call set_particle_type(i,igas)
       ndomain = ndomain+1
    endif

 enddo

 return
end subroutine determine_particle_status


!----
! Subroutine fills the injection zones with boundary particles
!---

subroutine replenish_injection_zone(ninject,time,dtlast,injected)
 use eos,        only:polyk,gamma
 use io,         only:id,master
 use part,       only:igas,iboundary,npart,npartoftype,xyzh,vxyzu,massoftype
 use partinject, only:add_or_update_particle
 use physcon,    only:pi
 use setdisc,    only:set_disc


 integer, intent(inout) :: ninject, injected
 integer :: i, npart_initial, i_part, part_type

 real :: hfact = 1.2
 real :: v0(3),v_subtract(3)
 real :: sig0,sig_in,cs0,vmag,omega0,vmax,time,dtlast
 real :: phi1,phi2,rpart,phipart,phi_rot,omega,tcross
 real :: xyzh_inject(4,ninject)
 real :: vxyzu_inject(4,ninject)


 npart_initial = npart
 vmag   = sqrt(injp%object_mass/injp%R_mid)
 omega0 = sqrt(injp%object_mass/(injp%R_mid*injp%R_mid*injp%R_mid))
 vmax   = sqrt(injp%object_mass/injp%Rsect_in) - vmag

 cs0    = injp%HoverR*sqrt(injp%object_mass)*injp%R_in**(injp%q_index-0.5)
 sig0   = sigma0(injp%disc_mass,injp%R_in,injp%R_out,injp%p_index)
 sig_in = sig0*injp%R_in**injp%p_index


! Set up both injection zones in the upper phi domain initially
 phi1 = injp%phimax
 phi2 = injp%phimax + injp%phi_inject

! call set_disc to establish new particle set

 call set_disc(id,master  = master,                       &
               npart      = ninject,                      &
               rmin       = injp%Rsect_in-injp%dr_bound,  &
               rmax       = injp%Rsect_out+injp%dr_bound, &
               phimin     = phi1,                         &
               phimax     = phi2,                         &
               p_index    = injp%p_index,                 &
               q_index    = injp%q_index,                 &
               HoverR     = injp%HoverR,                  &
               sig_norm   = sig_in,                       &
               star_mass  = injp%object_mass,             &
               polyk      = polyk,                        &
               gamma      = gamma,                        &
               particle_mass = massoftype(1),             &
               hfact      = hfact,                        &
               xyzh       = xyzh_inject,                  &
               vxyzu      = vxyzu_inject,                 &
               writefile  = .false.)

! Rotate components with r < Rmid so that they are in the other injection zone

 do i=1,ninject

    call calc_polar_coordinates(rpart,phipart,xyzh_inject(1,i), xyzh_inject(2,i))

    if(rpart <=injp%R_mid) then
       phi_rot = -2.0*phipart
       call rotate_particle_z(xyzh_inject(:,i), vxyzu_inject(:,i), phi_rot)
    endif

 enddo


! Transform to the corotating frame
 v0 = (/0.0, vmag,0.0/)

 print *, 'Transforming to corotating frame: angular velocity ', omega0

 do i=1,ninject

    call calc_polar_coordinates(rpart,phipart,xyzh_inject(1,i), xyzh_inject(2,i))
    call rotate_vector_z(v0, v_subtract,phipart)
    vxyzu_inject(1:3,i) = vxyzu_inject(1:3,i)-v_subtract(:)

 enddo


! Now add injected particles to the simulation
! If not the first injection, replenishment should occur along lines of constant tcross
! tcross = abs(phimax-phi)/omega(r)

! We fix tcross to its value at Rsect_in

! This implies (phimax-phi) > omega(rpart)*tcross(Rsect_in) to be true for injection

 i_part = npart_initial
 part_type = iboundary

 tcross = time-dtlast

 injected= 0
 do i=1,ninject

    tcross = time-dtlast

    call calc_polar_coordinates(rpart,phipart,xyzh_inject(1,i), xyzh_inject(2,i))

    omega = sqrt(injp%object_mass/(rpart*rpart*rpart)) -omega0

!print*, abs(injp%phimax-phipart), omega, tcross, omega*tcross

    if(abs(injp%phimax-phipart) > omega*tcross) then

       i_part = i_part+1
       call add_or_update_particle(part_type, xyzh_inject(1:3,i), vxyzu_inject(1:3,i), xyzh_inject(4,i), &
    vxyzu_inject(4,i), i_part, npart, npartoftype, xyzh, vxyzu) ! Another Brick in the Wall

       injected = injected+1
    endif

 enddo

 return
end subroutine replenish_injection_zone

!----------------------------------
! Rotates a particle in the z axis
!----------------------------------

subroutine rotate_particle_z(xyz,vxyz,phi)

 real, intent(inout) :: phi, xyz(3), vxyz(3)
 real :: x,y,vx,vy

 x = xyz(1)
 y = xyz(2)

 vx = vxyz(1)
 vy = vxyz(2)

 xyz(1) = x*cos(phi) - y*sin(phi)
 xyz(2) = x*sin(phi) + y*cos(phi)

 vxyz(1) = vx*cos(phi) - vy*sin(phi)
 vxyz(2) = vx*sin(phi) + vy*cos(phi)

 return
end subroutine rotate_particle_z

!----------------------------------
!+
! Rotates a single vector in the z axis
!+
!-----------------------------------

subroutine rotate_vector_z(oldvec,newvec,phi)
 real, intent(inout) :: oldvec(3), newvec(3)
 real, intent(in) :: phi

 newvec(1) = oldvec(1)*cos(phi) - oldvec(2)*sin(phi)
 newvec(2) = oldvec(1)*sin(phi) + oldvec(2)*cos(phi)

 return
end subroutine rotate_vector_z

!
!+
! Helper function to calculate polar co-ordinates from x,y
!+
!

subroutine calc_polar_coordinates(r,phi,x,y)

 real, intent(in) :: x,y
 real,intent(inout) :: r,phi

 r = sqrt(x*x + y*y)
 phi = atan2(y,x)

 return
end subroutine calc_polar_coordinates


!-----------------------------------------------------------------------
!+
!  Simple function to calculate the disc's surface density normalisation
!  for a disc mass, inner and outer radii and the powerlaw index
!+
!-----------------------------------------------------------------------

real function sigma0(Mdisc, Rinner, Router, p_index)
 real, intent(in) :: Mdisc,Rinner, Router, p_index

 real :: exponent

 sigma0 = Mdisc/(2.0*3.141592654)
 exponent = 2.0-p_index
 if(p_index==2.0) then
    sigma0 = sigma0*log(Rinner/Router)
 else
    sigma0 = sigma0*exponent/(Router**exponent - Rinner**exponent)
 endif

end function sigma0

end module inject
