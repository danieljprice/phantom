!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setBfield
!
!  DESCRIPTION:
!   Interactive setup of magnetic field on the particles
!
!    Can be used to add magnetic field to hydro setups, and
!   is used by utilities like moddump to add magnetic field
!   to purely hydrodynamic dump files before continuing the
!   calculation
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, physcon, prompting, setup_params, units
!+
!--------------------------------------------------------------------------
module setBfield
 implicit none
 public :: set_Bfield

 private

contains

subroutine set_Bfield(npart,npartoftype,xyzh,massoftype,vxyzu,polyk, &
                      Bxyz,Bextx,Bexty,Bextz)
 use units,        only:unit_Bfield
 use setup_params, only:rmax,rhozero
 use physcon,      only:pi
 use io,           only:fatal
 use prompting,    only:prompt
 integer,      intent(in)  :: npart
 integer,      intent(in)  :: npartoftype(:)
 real,         intent(in)  :: xyzh(:,:), vxyzu(:,:)
 real,         intent(in)  :: massoftype(:)
 real,         intent(in)  :: polyk
 real,         intent(out) :: Bxyz(:,:)
 real,         intent(out) :: Bextx,Bexty,Bextz
 integer :: maxp
 integer :: igeom,i
 character(len=1)  :: isetB
 character(len=10) :: string
 real :: totmass,przero,fracx,fracy,fracz,fractot
 real :: Bzero,Bzero2,Bxzero,Byzero,Bzzero
 real :: c1,area,rmasstoflux_crit,rmasstoflux,betazero,valfven
 real :: theta
 logical :: reverse_field_dir,ians

 maxp = size(xyzh(1,:))
!
! print general starting info
!
 print "(/,a,/)",' --- MAGNETIC FIELD SETUP --- '
!
! choose field geometry
!
 print "(2(/,a))",' 1) uniform cartesian field ', &
                  ' 2) uniform toroidal field  '
 igeom = 1
 call prompt('Choose initial magnetic field geometry ',igeom,1,2)

! set initial mean pressure for use in beta calculations
 przero = polyk*rhozero
!
!  set critical mass-to-flux ratio in code units
!  (code units are G=1, mu_0=1 -- see Price & Monaghan 2004a for details)
!  c1 is a normalisation factor taken from Mouschovias & Spitzer 1976
!
 c1 = 0.53
 rmasstoflux_crit = 2./3.*c1*sqrt(5./pi)
 totmass = npart*massoftype(1)
 area = pi*rmax**2
!
! choose field strength
!
 print "(/,' Do you want to enter: ',/,"// &
        "'         magnetic field strength in Gauss (m)',/,"// &
        "' or magnetic field strength in code units (c)',/,"// &
        "'            or Alfven speed in code units (a)',/,"// &
        "'                      or mean plasma beta (b)',/,"// &
        "'                or the mass-to-flux ratio (f)?')"

 isetB = 'm'
 call prompt('Enter your selection ',isetB,list=(/'m','c','a','b','f'/),noblank=.true.)

 Bzero = 0.
 select case(isetB)
 case('m')

    print*,' Code units of mag flux density = ',unit_Bfield
    Bzero = 0.
    call prompt('Enter mag field strength in Gauss ',Bzero)
    Bzero = Bzero/unit_Bfield

 case('c')

    print*,' Enter mag field strength in units of ',unit_Bfield
    write(string,"(es10.3)") unit_Bfield
    call prompt('Enter mag field strength in units of '//trim(string)//' G',Bzero)

 case('a')

    if (rhozero <= 0.) then
       print*,' ERROR: rhozero = ',przero,' from particle setup'
       rhozero = 1.
       call prompt('Enter density in code units used to compute the Alfven speed',rhozero,0.)
    endif
    valfven = 0.
    call prompt('Enter Alfven speed in code units ',valfven,0.)
    Bzero = sqrt(rhozero)*valfven

 case('b')

    if (przero <= 0.) then
       print*,' ERROR: przero = ',przero,' from particle setup'
       przero = 1.
       call prompt('Enter gas pressure in code units used to compute the plasma beta',przero,0.)
    endif
    print*,' mean gas pressure = ',przero
    betazero = 0.
    call prompt('Enter mean plasma beta assuming uniform rho, T (0=hydro)',betazero,0.)

    if (abs(betazero) < tiny(betazero)) then
       Bzero = 0.
    else
       Bzero = sqrt(2.*przero/betazero)
    endif

 case('f')
!
!   set area for mass-to-flux calculation (assumed spherical at the moment)
!
    if (rmax <= 0. .or. rmax > 1e3) then
       rmax = 0.
       call prompt('enter rmax for mass-to-flux ratio calculation ',rmax,0.)
    endif
    area = pi*rmax**2

    print*,' critical mass-to-flux ratio (code units) = ',rmasstoflux_crit
    print "('  < 1 = subcritical,',/,"// &
            "'    1 = critical,',/,"// &
            "'  > 1 = supercritical, ',/,"// &
            "'    0=inf=hydro )')"
    rmasstoflux = 0.
    call prompt('Enter spherical mass-to-flux ratio ',rmasstoflux,0.)

    if (abs(rmasstoflux) < tiny(rmasstoflux)) then
       Bzero = 0.
    else
       Bzero = totmass/(area*rmasstoflux*rmasstoflux_crit)
    endif

 case default
    print*,' error: unknown setting for B, setting B = 0'
    Bzero = 0.
 end select

!
! now actually do setup
!
 select case(igeom)
 case(1)
!
!--uniform cartesian field
!
    print "(a)",' Enter Bx:By:Bz ratio'
    fracx = 1.
    fracy = 0.
    fracz = 0.
    call prompt(' Enter Bx fraction ',fracx,0.)
    call prompt(' Enter By fraction ',fracy,0.)
    call prompt(' Enter Bz fraction ',fracz,0.)
!    read*,fracx,fracy,fracz

    fractot = sqrt(fracx**2 + fracy**2 + fracz**2)
    if (fractot < tiny(0.)) fractot = 1.
    Bxzero = fracx*Bzero/fractot
    Byzero = fracy*Bzero/fractot
    Bzzero = fracz*Bzero/fractot
!
!--spit out actual settings
!
    print "(' Bx_0 = ',es14.5,/,"// &
           "' By_0 = ',es14.5,/,"// &
           "' Bz_0 = ',es14.5)",Bxzero,Byzero,Bzzero
!
!--B
!
    do i=1,npart
       Bxyz(1,i) = Bxzero
       Bxyz(2,i) = Byzero
       Bxyz(3,i) = Bzzero
    enddo

 case(2)
!
!--uniform toroidal field
!
    Bxzero = 0.
    Byzero = 0.
    Bzzero = 0.
    reverse_field_dir = .false.
    if (Bzero >= 0.) then
       ians = .false.
       call prompt('Set field in clockwise direction? (discs in Phantom rotate anti-clockwise)',ians)
       if (.not.ians) reverse_field_dir = .true.
    else
       ians = .false.
       call prompt('Set field in clockwise direction? (discs in Phantom rotate anti-clockwise)',ians)
       if (ians) reverse_field_dir = .true.
    endif
    if (reverse_field_dir) Bzero = -Bzero
    do i=1,npart
       theta=atan2(xyzh(2,i),xyzh(1,i))
       Bxyz(1,i) = Bzero*sin(theta)
       Bxyz(2,i) = -Bzero*cos(theta)
       Bxyz(3,i) = 0.
    enddo

 case default
    call fatal('set_Bfield','unknown field geometry')
 end select

!
!--spit out various information about the magnetic field we have set up
!
 if (rhozero > 0. .and. przero > 0.) then
    Bzero2 = Bzero**2
    valfven = sqrt(Bzero2/rhozero)
    if (Bzero2 > 0.) then
       betazero = przero/(0.5*Bzero2)
    else
       betazero = 0.
    endif
    print "(' Alfven speed = ',es12.4,/,' Plasma beta  = ',es12.4)",valfven,betazero
 endif

!
!--spit out flux to mass ratio (assumes spherical geometry at the moment)
!
 if (area > tiny(area)) then
    if (Bxzero > tiny(0.)) print 10,'x',totmass/(area*Bxzero),totmass/(area*Bxzero)/rmasstoflux_crit
    if (Byzero > tiny(0.)) print 10,'y',totmass/(area*Byzero),totmass/(area*Byzero)/rmasstoflux_crit
    if (Bzzero > tiny(0.)) print 10,'z',totmass/(area*Bzzero),totmass/(area*Bzzero)/rmasstoflux_crit
 endif
10 format (' Mass to flux ratio (',a1,') = ',es12.4,' r/rcrit = ',f9.4)

!
!--spit out field strength in CGS units
!
 print "(' Initial field strength (Bzero) = ',es12.4,' G',"// &
        "' (',es12.4,' in code units)')",Bzero*unit_Bfield,Bzero

!
!--set external field components
! (required for const. stress boundaries and current projection)
!
 Bextx = Bxzero
 Bexty = Byzero
 Bextz = Bzzero
 print "(' Setting external B field = ',3(es12.4,2x),//)",Bextx,Bexty,Bextz

 return
end subroutine set_Bfield

end module setBfield
