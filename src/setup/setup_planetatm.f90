!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Set up for the de Val Borro et al. planet-disc comparison problem
!
!  REFERENCES: de Val Borro et al. (2006), MNRAS 370, 529
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    HoverRinput  -- H/R at R_in
!    R_in         -- inner radius
!    R_out        -- outer radius
!    accradius1   -- primary accretion radius
!    accradius2   -- secondary accretion radius
!    alphaSS      -- desired alpha_SS
!    mplanet      -- m1/(m1+m2)
!    np           -- number of particles
!    p_indexinput -- surface density profile
!    q_indexinput -- temperature profile
!    sig0         -- disc surface density
!
!  DEPENDENCIES: extern_binary, externalforces, infile_utils, io, options,
!    physcon, prompting, setdisc, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 integer :: np
 integer :: atm_type = 0
 real :: R_in, R_out, HoverRinput, sig0, alphaSS
 real :: p_indexinput, q_indexinput
 real :: Ratm_in,Ratm_out,npart_planet_frac,rho_core_cgs

 private

contains
!----------------------------------------------------------------
!
! spherical density profile as a function of radius
!
!----------------------------------------------------------------
real function atm_dens(r)
 use eos, only:gamma
 real, intent(in) :: r

 select case(atm_type)
  case(0)
     !atm_dens = exp(-(r-r_planet)/scaleheight)
     atm_dens = r**(-3)
  case(1)
     atm_dens = r**(-1./(gamma - 1.))
  case default
     atm_dens = r**(-3)
 end select

end function atm_dens

!----------------------------------------------------------------
!
! This subroutine sets up planet-disc interaction problem
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,         only:set_disc
 use units,           only:set_units,udist,umass,utime
 use physcon,         only:solarm,au,pi,gg
 use io,              only:master
 use options,         only:iexternalforce,alpha
 use timestep,        only:dtmax
 use prompting,       only:prompt
 use extern_binary,   only:accradius1,accradius2 !,binary_posvel
 use extern_binary,   only:binarymassr,eps_soft1,eps_soft2,ramp
 use externalforces,  only:iext_binary,iext_corot_binary,iext_corotate
 use extern_corotate, only:omega_corotate
 use spherical,       only:set_sphere
 use part,            only:set_particle_type

 integer,            intent(in)            :: id
 integer,            intent(out)           :: npart
 integer,            intent(out)           :: npartoftype(:)
 real,               intent(out)           :: xyzh(:,:)
 real,               intent(out)           :: polyk,gamma,hfact
 real,               intent(out)           :: vxyzu(:,:)
 real,               intent(out)           :: massoftype(:)
 real,               intent(inout)         :: time
 character (len=20), intent (in), optional :: fileprefix
 integer :: i
 integer :: nx,npart_disc,npart_planet_atm
 integer(kind=8) :: nptot
 integer :: maxvxyzu,itype
 integer, parameter :: igas=1

 !real :: xbinary(10),vbinary(6)
 real :: a0,Mstar
 real :: vmag,omega0,v_0(3)!,v_subtract(3)
 real :: phipart,r
 real :: udens,rho_core,r_surface,a_orbit
 real :: Mplanet,Mearth,Mjupiter
 real :: xyz_orig(3),psep,vol_sphere
 real :: cs0,cs,G

 logical :: iexist
 character(len=100) :: filename

 !
 !--set code units
 !
 call set_units(dist=0.25*au,mass=solarm,G=1.d0)
 G = gg*umass*utime**2/(udist**3)

 filename=trim(fileprefix)//'.setup'

 !
 !--set default options for the input file
 !
 np = size(xyzh(1,:))
 npart = np
 npartoftype(1) = npart
 maxvxyzu = size(vxyzu(:,1))
 gamma = 1.0
 hfact = 1.2
 time  = 0.
 a0    = 1.
 Mstar = 1.
 !binarymassr = 1.e-3
 !-----------------------
 Mearth   = 5.976e27/umass
 Mjupiter = 1.899e30/umass
 !binarymassr = Mjupiter/(Mstar + Mjupiter) ! Jupiter mass
 !binarymassr = Mearth/(Mstar + Mearth) ! Earth mass
 binarymassr = 15.*Mearth/(Mstar + Mearth) ! Earth mass
 !-----------------------
 Mplanet = Mstar*binarymassr/(1. - binarymassr)
 HoverRinput = 0.05
 accradius1 = 0.0
 accradius2 = 0.075*a0
 R_in  = 0.1*a0
 R_out  = 2.5*a0
 sig0   = 0.002/(pi*a0**2)
 alphaSS = 0.01
 p_indexinput = 1.
 q_indexinput = 0.5
 ramp = .false.!.true.

 !-----------------------
 !iexternalforce = iext_binary
 !iexternalforce = iext_corotate
 iexternalforce = iext_corot_binary
 !-----------------------
 rho_core_cgs = 5.
 Ratm_in  = 1.
 Ratm_out = 3.
 npart_planet_frac = 0.0000001

 if (iexternalforce == iext_corot_binary) then
    udens = umass/udist**3
    rho_core  = rho_core_cgs/udens
    r_surface = (3./(4.*pi)*Mplanet/rho_core)**(1./3.)
    eps_soft1 = r_surface
 else
    eps_soft1 = 0.6*HoverRinput*a0
 endif


 print "(a,/)",'Phantomsetup: routine to setup planet-disc interaction with fixed planet orbit '
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_gwinputfile(filename)
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    !
    !--set default options
    !
    call prompt('Enter total number of gas particles ',np,0,size(xyzh(1,:)))
    call prompt('Enter mplanet/mtot',binarymassr,0.,1.)

    call prompt('Enter accretion radius of the PRIMARY (planet)',accradius1,accradius1,1.)
    call prompt('Enter accretion radius of the SECONDARY (star)',accradius2,accradius2,1.)

    call prompt('Enter softening radius of the PRIMARY (planet)',eps_soft1,0.,1.)
    call prompt('Enter softening radius of the SECONDARY (star)',eps_soft2,0.,1.)

    call prompt('Enter inner disc edge R_in ',R_in,accradius1)
    call prompt('Enter outer disc edge R_out ',R_out,R_in)
    call prompt('Enter H/R at R=R_in ',HoverRinput,0.)
    call prompt('Enter p index of surface density profile Sigma = Sigma0*R^-p',p_indexinput,0.)
    call prompt('Enter q index of temperature profile cs = cs0*R^-q',q_indexinput,0.)
    call prompt('Enter Sigma0 for disc',sig0)
    call prompt('Enter desired value of alpha_SS',alphaSS,0.)

    if (iexternalforce == iext_corot_binary) then
       call prompt('Enter core density in cgs units',rho_core_cgs,0.)
       call prompt('Enter inner atmosphere radius in planet radii',Ratm_in,1.,100.)
       call prompt('Enter outer atmosphere radius in planet radii',Ratm_out,Ratm_in,100.)
       call prompt('Enter atmosphere type (0:isothermal; 1:adiabatic)',atm_type,0,1)
       call prompt('Enter fraction of particles to be used in planet atmosphere',npart_planet_frac,0.,1.)
    endif

    !
    !--write default input file
    !
    call write_gwinputfile(filename)

    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif

 if (iexternalforce == iext_corot_binary) then
    npart_planet_atm = floor(npart_planet_frac*np)
    npart_disc = np - npart_planet_atm
    npart = npart_disc
 else
    npart = np
    npartoftype(1) = npart
 endif

 alpha = alphaSS
 itype = 1

 call set_disc(id,master=master,&
               npart   = npart,&
               rmin    = R_in,   &
               rmax    = R_out,  &
               p_index = p_indexinput,    &
               q_index = q_indexinput,   &
               HoverR  = HoverRinput,  &
               sig_norm = sig0, &
               star_mass = Mstar,  &
               gamma     = gamma,  &
               particle_mass = massoftype(1), &
               hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,polyk=polyk,alpha=alpha, &
               prefix = fileprefix )

 if (iexternalforce == iext_corot_binary .and. npart_planet_atm > 0) then
    !
    ! place particles in sphere
    !
    Ratm_in   = Ratm_in*r_surface
    Ratm_out  = Ratm_out*r_surface

    if (ramp) then
       xyz_orig(:) = (/a0,0.,0./)
    else
       a_orbit = a0 - binarymassr
       xyz_orig(:) = (/a_orbit,0.,0./)
    endif
    vol_sphere  = 4./3.*pi*Ratm_out**3
    nx          = int(npart_planet_atm**(1./3.))
    psep        = vol_sphere**(1./3.)/real(nx)
    nptot       = npart

    call set_sphere('closepacked',id,master,Ratm_in,Ratm_out,psep,hfact,npart,xyzh, &
                    rhofunc=atm_dens,nptot=nptot, &
                    np_requested=npart_planet_atm,xyz_origin=xyz_orig)

    npart_planet_atm = npart-npart_disc
    npartoftype(1) = npart
    do i = npart_disc+1,npart
       !--set the particle type for the atmosphere particles
       call set_particle_type(i,1)
       !-----------------------------------------
       !  Set thermal energy
       !  utherm generally should not be stored
       !  for an isothermal equation of state
       if (maxvxyzu >= 4) then
          if (itype==igas) then
             cs0 = sqrt(1.0/R_in)*(HoverRinput)*sqrt(G*Mstar)*(1.0/R_in)**(-q_indexinput)
             cs = cs_func(cs0,a0,q_indexinput)
             if (gamma > 1.) then
                vxyzu(4,i) = cs**2/(gamma - 1.)/gamma
             else
                vxyzu(4,i) = 1.5*cs**2
             endif
          else
             vxyzu(4,i) = 0.
          endif
       endif
    enddo
 endif

 if (iexternalforce == iext_corotate .or. &
     iexternalforce == iext_corot_binary) then
    !
    !--Change to corotating frame
    !
    ! Calculate velocity at planet
    vmag = sqrt(Mstar/a0)
    omega0 = sqrt(Mstar/a0**3)

    ! v_phi = v_y at y=0
    ! Obtain the true v_phi at any point (r,phi) via rotation in z axis

    v_0 = (/0.0,vmag,0.0/)

    print *, 'Transforming to corotating frame: angular velocity ', omega0

    do i=1,npart_disc
       phipart = atan2(xyzh(2,i),xyzh(1,i))
       r = sqrt(xyzh(1,i)**2+xyzh(2,i)**2)
       !call rotate_z(v_0, v_subtract,phipart)
       !vxyzu(1:3,i) = vxyzu(1:3,i)-r*v_subtract(:)
       vxyzu(1,i) = vxyzu(1,i) - r*(-omega0)*sin(phipart)
       vxyzu(2,i) = vxyzu(2,i) + r*(-omega0)*cos(phipart)
    enddo
    vxyzu(1:3,npart_disc+1:npart) = 0.

    omega_corotate = omega0
 endif

 dtmax = 2.*pi/100.

 !--------------------------------------------------
 ! If you want to translate the disc so it is around the primary uncomment the following lines
 !--------------------------------------------------
! call binary_posvel(time,xbinary,vbinary)
! do i=1,npart
!   xyzh(1,i) = xyzh(1,i) + xbinary(1)
!   xyzh(2,i) = xyzh(2,i) + xbinary(2)
!   xyzh(3,i) = xyzh(3,i) + xbinary(3)
!   vxyzu(1,i) = vxyzu(1,i) + vbinary(1)
!   vxyzu(2,i) = vxyzu(2,i) + vbinary(2)
!   vxyzu(3,i) = vxyzu(3,i) + vbinary(3)
! enddo

 return
end subroutine setpart


subroutine write_gwinputfile(filename)
 use infile_utils,   only:write_inopt
 use extern_binary,  only:accradius1,accradius2,binary_posvel
 use extern_binary,  only:binarymassr
 use options,        only:iexternalforce
 use externalforces, only:iext_corot_binary
 implicit none
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for gwdisc setup routines'

 write(iunit,"(/,a)") '# resolution'

 call write_inopt(np,'np','number of particles',iunit)

 write(iunit,"(/,a)") '# options for binary'

 call write_inopt(binarymassr,'mplanet','m1/(m1+m2)',iunit)
 call write_inopt(accradius1,'accradius1','primary accretion radius',iunit)
 call write_inopt(accradius2,'accradius2','secondary accretion radius',iunit)

 write(iunit,"(/,a)") '# options for accretion disc'

 call write_inopt(R_in,'R_in','inner radius',iunit)
 call write_inopt(R_out,'R_out', 'outer radius',iunit)
 call write_inopt(HoverRinput,'HoverRinput','H/R at R_in',iunit)
 call write_inopt(sig0,'sig0','disc surface density',iunit)
 call write_inopt(p_indexinput,'p_indexinput','surface density profile',iunit)
 call write_inopt(q_indexinput,'q_indexinput','temperature profile',iunit)
 call write_inopt(alphaSS,'alphaSS','desired alpha_SS',iunit)

 if (iexternalforce == iext_corot_binary) then
    write(iunit,"(/,a)") '# options for planet atmosphere'

    call write_inopt(rho_core_cgs,'rho_core','planet core density (cgs units)',iunit)
    call write_inopt(Ratm_in,'Ratm_in','inner atmosphere radius (planet radii)',iunit)
    call write_inopt(Ratm_out,'Ratm_out','outer atmosphere radius (planet radii)',iunit)
    call write_inopt(atm_type,'atm_type','Enter atmosphere type (0:isothermal; 1:adiabatic)',iunit)
    call write_inopt(npart_planet_frac,'Natm/Npart','fraction of particles for planet atmosphere',iunit)
 endif

 close(iunit)

end subroutine write_gwinputfile

subroutine read_gwinputfile(filename)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use extern_binary, only:accradius1,accradius2,binary_posvel
 use extern_binary, only:binarymassr
 use options,        only:iexternalforce
 use externalforces, only:iext_corot_binary
 implicit none
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 21
 integer :: ierr
 type(inopts), dimension(:), allocatable :: db

 print "(a)",'reading setup options from '//trim(filename)

 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(np,'np',db,ierr)
 call read_inopt(binarymassr,'mplanet',db,ierr)
 call read_inopt(accradius1,'accradius1',db,ierr)
 call read_inopt(accradius2,'accradius2',db,ierr)
 call read_inopt(R_in,'R_in',db,ierr)
 call read_inopt(R_out,'R_out',db,ierr)
 call read_inopt(HoverRinput,'HoverRinput',db,ierr)
 call read_inopt(sig0,'sig0',db,ierr)
 call read_inopt(p_indexinput,'p_indexinput',db,ierr)
 call read_inopt(q_indexinput,'q_indexinput',db,ierr)
 call read_inopt(alphaSS,'alphaSS',db,ierr)

 if (iexternalforce == iext_corot_binary) then
    call read_inopt(rho_core_cgs,'rho_core',db,ierr)
    call read_inopt(Ratm_in,'Ratm_in',db,ierr)
    call read_inopt(Ratm_out,'Ratm_out',db,ierr)
    call read_inopt(atm_type,'atm_type',db,ierr)
    call read_inopt(npart_planet_frac,'Natm/Npart',db,ierr)
 endif

 call close_db(db)
 close(iunit)

end subroutine read_gwinputfile

!----------------------------------------------------------------
!
! function to return the sound speed given the radius
!
!----------------------------------------------------------------
pure real function cs_func(cs0,r,q_index)
 real, intent(in) :: cs0,r,q_index

 cs_func = cs0*r**(-q_index)

end function cs_func

!!----------------------------------------------------------------
!!
!! Rotates a vector in the z axis
!!
!!----------------------------------------------------------------
!subroutine rotate_z(oldvec,newvec,phi)
! real, intent(inout) :: oldvec(3), newvec(3)
! real, intent(in) :: phi
!
! newvec(1) = oldvec(1)*cos(phi) - oldvec(2)*sin(phi)
! newvec(2) = oldvec(1)*sin(phi) + oldvec(2)*cos(phi)
!
! return
!end subroutine rotate_z

end module setup
