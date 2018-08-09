!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: part
!
!  DESCRIPTION:
!  This module contains main particle data for the code and
!  functions used to query particle properties not stored.
!
!  Basically this module defines any quantity that is
!  stored on the particles and defines how that storage
!  is implemented. Thus any routine which requires knowledge
!  of the specifics of this storage should be placed here.
!  (for example, any routine that copies all of the variables
!   stored on a given particle).
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, domain, io, mpiutils
!+
!--------------------------------------------------------------------------
module part
 use dim, only:ndim,maxp,maxsts,ndivcurlv,ndivcurlB,maxvxyzu, &
          maxalpha,maxptmass,maxdvdx, &
          mhd,maxmhd,maxBevol,maxp_h2,maxtemp,periodic, &
          maxgrav,ngradh,maxtypes,h2chemistry,gravity, &
          maxp_dustfrac,use_dust, &
          store_temperature,lightcurve,maxlum,nalpha,maxmhdni, &
          maxne,maxp_growth,maxdustlarge,maxdustsmall,maxdusttypes
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"
!
!--basic storage needed for read/write of particle data
!
 real :: xyzh(4,maxp)
 real :: xyzh_soa(maxp,4)
 real :: vxyzu(maxvxyzu,maxp)
 real(kind=4) :: alphaind(nalpha,maxalpha)
 real(kind=4) :: divcurlv(ndivcurlv,maxp)
 real(kind=4) :: divcurlB(ndivcurlB,maxp)
 real :: Bevol(maxBevol,maxmhd)
 real :: Bxyz(3,maxmhd)
 character(len=*), parameter :: xyzh_label(4) = (/'x','y','z','h'/)
 character(len=*), parameter :: vxyzu_label(4) = (/'vx','vy','vz','u '/)
 character(len=*), parameter :: Bxyz_label(3) = (/'Bx','By','Bz'/)
!
!--storage of dust properties
!
 real :: dustprop(4,maxp_growth)
 real :: St(maxp_growth)
 character(len=*), parameter :: dustprop_label(4) = (/'grainsize ','graindens ','vrel/vfrag','    dv    '/)
!
!--storage in divcurlv
!
 integer, parameter :: idivv = 1
 integer, parameter :: icurlvx = 2
 integer, parameter :: icurlvy = 3
 integer, parameter :: icurlvz = 4
 character(len=*), parameter :: divcurlv_label(4) = &
   (/'divv  ','curlvx','curlvy','curlvz'/)
!
!--storage in divcurlB
!
 integer, parameter :: idivB = 1
 integer, parameter :: icurlBx = 2
 integer, parameter :: icurlBy = 3
 integer, parameter :: icurlBz = 4
 character(len=*), parameter :: divcurlB_label(4) = &
   (/'divB  ','curlBx','curlBy','curlBz'/)
!
!--velocity gradients
!
 real(kind=4) :: dvdx(9,maxdvdx)
!
!--H2 chemistry
!
 integer, parameter :: nabundances = 5
 integer, parameter :: ih2ratio  = 1 ! ratio of H2 to H
 integer, parameter :: iHI       = 2 ! HI abundance
 integer, parameter :: iproton   = 3 ! proton abundance
 integer, parameter :: ielectron = 4 ! electron abundance
 integer, parameter :: iCO       = 5 ! CO abundance
 real :: abundance(nabundances,maxp_h2)
 character(len=*), parameter :: abundance_label(5) = &
   (/'h2ratio','abHIq  ','abhpq  ','abeq   ','abco   '/)
!
!--storage of temperature
!
 real :: temperature(maxtemp)
!
!--one-fluid dust (small grains)
!
 real :: dustfrac(maxdusttypes,maxp_dustfrac)
 character(len=*), parameter :: dustfrac_label(maxdusttypes) = 'dustfrac'
 character(len=*), parameter :: tstop_label(maxdusttypes) = 'tstop'
 real :: dustevol(maxdusttypes,maxp_dustfrac)
 real :: deltav(3,maxdusttypes,maxp_dustfrac)
 character(len=*), parameter :: deltav_label(3) = &
   (/'deltavx','deltavy','deltavz'/)
!
!--sink particles
!
 integer, parameter :: nsinkproperties = 11
 integer, parameter :: ihacc  = 5 ! accretion radius
 integer, parameter :: ihsoft = 6 ! softening radius
 integer, parameter :: imacc  = 7 ! accreted mass
 integer, parameter :: ispinx = 8  ! spin angular momentum x
 integer, parameter :: ispiny = 9  ! spin angular momentum y
 integer, parameter :: ispinz = 10 ! spin angular momentum z
 integer, parameter :: i_tlast = 11 ! time of last injection
 real :: xyzmh_ptmass(nsinkproperties,maxptmass)
 real :: vxyz_ptmass(3,maxptmass)
 real :: fxyz_ptmass(4,maxptmass),fxyz_ptmass_sinksink(4,maxptmass)
 integer :: nptmass = 0   ! zero by default
 real    :: epot_sinksink
 character(len=*), parameter :: xyzmh_ptmass_label(11) = &
  (/'x        ','y        ','z        ','m        ','h        ',&
    'hsoft    ','maccreted','spinx    ','spiny    ','spinz    ','tlast    '/)
 character(len=*), parameter :: vxyz_ptmass_label(3) = (/'vx','vy','vz'/)
!
!--self-gravity
!
 real(kind=4) :: poten(maxgrav)
!
!--Non-ideal MHD
!
 real :: n_R(4,maxmhdni),n_electronT(maxne),eta_nimhd(4,maxmhdni)
 integer, parameter :: iohm  = 1 ! eta_ohm
 integer, parameter :: ihall = 2 ! eta_hall
 integer, parameter :: iambi = 3 ! eta_ambi
 integer, parameter :: iion  = 4 ! ionisation fraction
#ifdef NONIDEALMHD
 character(len=*), parameter :: eta_nimhd_label(4) = (/'eta_{OR}','eta_{HE}','eta_{AD}','ne/n    '/)
#endif
!
!--for analysis routines, do not allocate any more storage
!  than is strictly necessary
!
#ifdef ANALYSIS
 integer, parameter, private :: maxan = 0
 integer, parameter, private :: maxmhdan = 0
 integer, parameter, private :: maxdustan = 0
#else
 integer, parameter, private :: maxan = maxp
 integer, parameter, private :: maxmhdan = maxmhd
 integer, parameter, private :: maxdustan = maxp_dustfrac
#endif
!
!--lightcurves
!
 real(kind=4) :: luminosity(maxlum)
!
!--derivatives (only needed if derivs is called)
!
 real               :: fxyzu(maxvxyzu,maxan)
 real               :: dBevol(maxBevol,maxmhdan)
 real(kind=4)       :: divBsymm(maxmhdan)
 real               :: fext(3,maxan)
 real               :: ddustfrac(maxdusttypes,maxdustan)
 real               :: ddustprop(4,maxp_growth) !--grainsize is the only prop that evolves for now
!
!--storage associated with/dependent on timestepping
!
 real               :: vpred(maxvxyzu,maxan)
 real               :: dustpred(maxdusttypes,maxdustan)
 real               :: Bpred(maxBevol,maxmhdan)
 real               :: dustproppred(4,maxp_growth)
#ifdef IND_TIMESTEPS
 integer(kind=1)    :: ibin(maxan)
 integer(kind=1)    :: ibin_old(maxan)
 integer(kind=1)    :: ibin_wake(maxan)
 real(kind=4)       :: dt_in(maxan)
 real               :: twas(maxan)
#else
 integer(kind=1)    :: ibin_wake(1)
#endif
 integer, parameter :: maxphase = maxan
 integer, parameter :: maxgradh = maxan
 integer(kind=1)    :: iphase(maxphase)
 integer(kind=1)    :: iphase_soa(maxphase)
 logical, public    :: all_active = .true.

 real(kind=4)       :: gradh(ngradh,maxgradh)
 real               :: tstop(maxdusttypes,maxan)
!
!--storage associated with link list
!  (used for dead particle list also)
!
 integer :: ll(maxan)
 real    :: dxi(ndim) ! to track the extent of the particles
!
!--size of the buffer required for transferring particle
!  information between MPI threads
!
 integer, parameter, private :: maxpd =  max(maxp,1) ! avoid divide by zero
 integer, parameter, private :: usedivcurlv = min(ndivcurlv,1)
 integer, parameter :: ipartbufsize = 4 &  ! xyzh
   +maxvxyzu                            &  ! vxyzu
   +maxvxyzu                            &  ! vpred
   +maxvxyzu                            &  ! fxyzu
   +3                                   &  ! fext
   +usedivcurlv                         &  ! divcurlv
   +nalpha*maxalpha/maxpd               &  ! alphaind
   +ngradh*maxgradh/maxpd               &  ! gradh
   +(maxmhd/maxpd)*maxBevol             &  ! Bevol
   +(maxmhd/maxpd)*maxBevol             &  ! Bpred
   +maxphase/maxpd                      &  ! iphase
#ifdef DUST
   +maxdusttypes                        &  ! dustfrac
   +maxdusttypes                        &  ! dustevol
#ifdef DUSTGROWTH
   +1                                   &  ! dustproppred
   +1                                   &  ! ddustprop
#endif
#endif
   +(maxp_h2/maxpd)*nabundances         &  ! abundance
   +(maxgrav/maxpd)                     &  ! poten
   +(maxtemp/maxpd)                     &  ! temperature
#ifdef IND_TIMESTEPS
   +1                                   &  ! ibin
   +1                                   &  ! ibin_old
   +1                                   &  ! ibin_wake
   +1                                   &  ! dt_in
   +1                                   &  ! twas
#endif
   +0

 real            :: hfact,Bextx,Bexty,Bextz
 integer         :: npart
 integer(kind=8) :: ntot
 integer         :: ideadhead = 0

 integer :: npartoftype(maxtypes)
 real    :: massoftype(maxtypes)

 integer :: ndustsmall,ndustlarge,ndusttypes
!
!--labels for each type
!  NOTE: If new particle is added, and it is allowed to be accreted onto
!        a sink particle, add it to the list in 'is_accretable' below
!  NOTE: set_boundaries_to_active = .true., but will be set to .false. at the
!        end of initial.  This will allow boundary particles to always be
!         initialised (even on restarts where not all arrays, e.g. gradh,
!         are not saved)
!
 integer, parameter :: igas        = 1
 integer, parameter :: iboundary   = 3
 integer, parameter :: istar       = 4
 integer, parameter :: idarkmatter = 5
 integer, parameter :: ibulge      = 6
 integer, parameter :: idust       = 7
 integer, parameter :: idustlast   = idust + maxdustlarge - 1
 integer, parameter :: iunknown    = 0
 logical            :: set_boundaries_to_active = .true.
 integer :: i
 character(len=5), dimension(maxtypes), parameter :: &
   labeltype = (/'gas  ','empty','bound','star ','darkm','bulge', &
                 ('dust ', i=idust,idustlast)/)
!
!--generic interfaces for routines
!
 interface hrho
  module procedure hrho4,hrho8,hrho4_pmass,hrho8_pmass,hrhomixed_pmass
 end interface hrho

 private :: hrho4,hrho8,hrho4_pmass,hrho8_pmass,hrhomixed_pmass

contains
!----------------------------------------------------------------
!+
!  this function determines the mass of the particle
!  where use_gas == .not. (maxphase==maxp)
!  currently used only in readwrite_dump
!+
!----------------------------------------------------------------
real function get_pmass(i,use_gas)
 integer, intent(in) :: i
 logical, intent(in) :: use_gas

 if (use_gas) then
    get_pmass = massoftype(igas)
 else
    if (iphase(i) /= 0) then
       get_pmass = massoftype(iamtype(iphase(i)))
    else
       get_pmass = massoftype(igas)
    endif
 endif

end function get_pmass
!
!----------------------------------------------------------------
!+
!  this function gives rho as a function of h
!  (used by output routines to get rho given that we only
!   store h)
!+
!----------------------------------------------------------------
pure real function rhoh(hi,pmassi)
 real, intent(in) :: hi,pmassi

 rhoh = pmassi*(hfact/abs(hi))**3

end function rhoh

!----------------------------------------------------------------
!+
!  this function gives dh/drho as a function of h
!+
!----------------------------------------------------------------
pure real function dhdrho(hi,pmassi)
 real, intent(in) :: hi,pmassi
 real :: rhoi

 rhoi   = rhoh(hi,pmassi)
 dhdrho = -hi/(3.*rhoi)

end function dhdrho
!----------------------------------------------------------------
!+
!  this subroutine does both of the above
!  (routine is an optimisation - also returns divisions)
!+
!----------------------------------------------------------------
subroutine rhoanddhdrho(hi,hi1,rhoi,rho1i,dhdrhoi,pmassi)
 real,         intent(in)  :: hi, pmassi
 real(kind=8), intent(out) :: hi1
 real,         intent(out) :: rhoi
 real,         intent(out) :: rho1i,dhdrhoi
 real, parameter :: third = 1./3.

 hi1 = 1./abs(hi)
 rhoi = pmassi*(hfact*hi1)**3
 rho1i = 1./rhoi
 dhdrhoi = -third*hi*rho1i

end subroutine rhoanddhdrho

!----------------------------------------------------------------
!+
!  this function gives h as a function of rho
!+
!----------------------------------------------------------------
real(kind=4) function hrho4(rhoi)
 real(kind=4), intent(in) :: rhoi

 hrho4 = real(hfact*(massoftype(1)/abs(rhoi))**(1./3.),kind=4)

end function hrho4

real(kind=8) function hrho8(rhoi)
 real(kind=8), intent(in) :: rhoi

 hrho8 = hfact*(massoftype(1)/abs(rhoi))**(1.d0/3.d0)

end function hrho8

real(kind=4) function hrho4_pmass(rhoi,pmassi)
 real(kind=4), intent(in) :: rhoi,pmassi

 hrho4_pmass = real(hfact*(pmassi/abs(rhoi))**(1./3.),kind=4)

end function hrho4_pmass

real(kind=8) function hrho8_pmass(rhoi,pmassi)
 real(kind=8), intent(in) :: rhoi,pmassi

 hrho8_pmass = hfact*(pmassi/abs(rhoi))**(1.d0/3.d0)

end function hrho8_pmass

real(kind=8) function hrhomixed_pmass(rhoi,pmassi)
 real(kind=8), intent(in) :: rhoi
 real(kind=4), intent(in) :: pmassi

 hrhomixed_pmass = hfact*(pmassi/abs(rhoi))**(1.d0/3.d0)

end function hrhomixed_pmass

!----------------------------------------------------------------
!+
!  query function returning whether or not a particle is dead
!  (currently indicated by having a negative smoothing length)
!+
!----------------------------------------------------------------
logical function isdead(i)
 integer, intent(in) :: i

 ! h = 0 indicates a dead particle
 if (abs(xyzh(4,i)) < tiny(xyzh)) then
    isdead = .true.
 else
    isdead = .false.
 endif

end function isdead

pure logical function isdeadh(hi)
 real, intent(in) :: hi

 ! h = 0 indicates a dead particle
 if (abs(hi) < tiny(hi)) then
    isdeadh = .true.
 else
    isdeadh = .false.
 endif

end function isdeadh

pure logical function isdead_or_accreted(hi)
 real, intent(in) :: hi

 ! h <= 0 indicates either dead or accreted
 if (hi < tiny(hi)) then
    isdead_or_accreted = .true.
 else
    isdead_or_accreted = .false.
 endif

end function isdead_or_accreted

!----------------------------------------------------------------
!+
! routine which kills a particle and adds it to the dead list
!+
!----------------------------------------------------------------
subroutine kill_particle(i,npoftype)
 integer, intent(in) :: i
 integer, intent(inout), optional :: npoftype(:)
 integer :: itype

 xyzh(4,i) = 0.
 if (present(npoftype)) then
    ! get the type so we know how to decrement npartoftype
    if (maxphase==maxp) then
       itype = iamtype(iphase(i))
    else
       itype = igas
    endif
    npoftype(itype) = npoftype(itype) - 1
 endif
!$omp critical
 ll(i) = ideadhead
 ideadhead = i
!$omp end critical

end subroutine kill_particle

!----------------------------------------------
!+
!  functions to deconstruct iphase
!  abs(iphase) is the particle type
!  sign(iphase) gives whether it is active/inactive
!+
!----------------------------------------------
pure integer(kind=1) function isetphase(itype,iactive)
 integer, intent(in) :: itype
 logical, intent(in) :: iactive

 if ((set_boundaries_to_active .and. itype==iboundary)   .or. &
     (iactive                  .and. itype/=iboundary) ) then
    isetphase = int(itype,kind=1)
 else
    isetphase = -int(abs(itype),kind=1)
 endif

end function isetphase

pure subroutine get_partinfo(iphasei,isactive,isdust,itype)
 integer(kind=1), intent(in)  :: iphasei
 logical,         intent(out) :: isactive,isdust
 integer,         intent(out) :: itype

! isactive = iactive(iphasei)
! itype = iamtype(iphasei)
! isdust = itype==idust

!--inline versions of above (for speed)
 if (iphasei > 0) then
    isactive = .true.
    itype    = iphasei
 else
    isactive = .false.
    itype    = -iphasei
 endif
#ifdef DUST
 isdust = ((itype>=idust) .and. (itype<=idustlast))
#else
 isdust = .false.
#endif

 return
end subroutine get_partinfo

pure logical function iactive(iphasei)
 integer(kind=1), intent(in) :: iphasei

 iactive = (iphasei > 0)

end function iactive

pure elemental integer function iamtype(iphasei)
 integer(kind=1), intent(in) :: iphasei

 iamtype = abs(iphasei)

end function iamtype

pure elemental function iamtype_int1(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer(kind=1) :: iamtype_int1

 iamtype_int1 = abs(iphasei)

end function iamtype_int1

pure function iamtype_int11(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer(kind=1) :: iamtype_int11

 iamtype_int11 = abs(iphasei)

end function iamtype_int11

pure elemental logical function iamgas(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer :: itype

 itype = iamtype(iphasei)
 iamgas = int(itype)==igas

end function iamgas

pure elemental logical function iamdust(iphasei)
 integer(kind=1), intent(in) :: iphasei
 integer :: itype

 itype = iamtype(iphasei)
 iamdust = ((itype>=idust) .and. (itype<=idustlast))

end function iamdust

pure integer function get_ntypes(noftype)
 integer, intent(in) :: noftype(:)
 integer :: i

 get_ntypes = 0
 do i=1,size(noftype)
    if (noftype(i) > 0) get_ntypes = i
 enddo

end function get_ntypes

!-----------------------------------------------------------------------
!+
!  Determine if particle is of a type that is accretable
!    Modify the if-statement to include all the types of particles
!    that can be accreted onto the sink
!+
!-----------------------------------------------------------------------
pure logical function is_accretable(itype)
 integer, intent(in)  :: itype

 if (itype==igas .or. itype==idust) then
    is_accretable = .true.
 else
    is_accretable = .false.
 endif

end function is_accretable

!----------------------------------------------------------------
!+
!  utility function for setup routines to set initial value
!  of iphase (assumes particle is active)
!+
!----------------------------------------------------------------
subroutine set_particle_type(i,itype)
 use io, only:fatal
 integer, intent(in) :: i,itype

 if (maxphase==maxp) then
    iphase(i) = isetphase(itype,iactive=.true.)
 elseif (itype /= igas) then
    call fatal('set_particle_type','attempt to setup a particle of type > 1, but iphase not allocated')
 endif

end subroutine set_particle_type

!----------------------------------------------------------------
!+
!  utility function to get strain tensor from dvdx array
!+
!----------------------------------------------------------------
pure function strain_from_dvdx(dvdxi) result(strain)
 real, intent(in) :: dvdxi(9)
 real :: strain(6)

 strain(1) = 2.*dvdxi(1)
 strain(2) = dvdxi(2) + dvdxi(4)
 strain(3) = dvdxi(3) + dvdxi(7)
 strain(4) = 2.*dvdxi(5)
 strain(5) = dvdxi(6) + dvdxi(8)
 strain(6) = 2.*dvdxi(9)

end function strain_from_dvdx

!----------------------------------------------------------------
!+
! routine which copies a particle from one location to another
! (prior to a derivs evaluation - so no derivs required)
!+
!----------------------------------------------------------------
subroutine copy_particle(src, dst)
 integer, intent(in) :: src, dst

 xyzh(:,dst)  = xyzh(:,src)
 vxyzu(:,dst) = vxyzu(:,src)
 fext(:,dst)  = fext(:,src)
 if (mhd) then
    Bevol(:,dst) = Bevol(:,src)
    Bxyz(:,dst)  = Bxyz(:,dst)
 endif
 if (ndivcurlv  > 0) divcurlv(:,dst)  = divcurlv(:,src)
 if (maxalpha ==maxp) alphaind(:,dst) = alphaind(:,src)
 if (maxgradh ==maxp) gradh(:,dst)    = gradh(:,src)
 if (maxphase ==maxp) iphase(dst)   = iphase(src)
 if (maxgrav  ==maxp) poten(dst) = poten(src)
#ifdef IND_TIMESTEPS
 ibin(dst)       = ibin(src)
 ibin_old(dst)   = ibin_old(src)
 ibin_wake(dst)  = ibin_wake(src)
 dt_in(dst)      = dt_in(src)
 twas(dst)       = twas(src)
#endif
 if (use_dust) then
    dustfrac(:,dst) = dustfrac(:,src)
    dustevol(:,dst) = dustevol(:,src)
 endif
 if (maxp_h2==maxp) abundance(:,dst) = abundance(:,src)
 if (store_temperature) temperature(dst) = temperature(src)

 return
end subroutine copy_particle

!----------------------------------------------------------------
!+
! routine which copies a particle from one location to another
! (copies everything which is stored on a particle)
!
! Note that link list information CANNOT be copied so link list
! must be rebuilt after a copy operation.
!+
!----------------------------------------------------------------
subroutine copy_particle_all(src,dst)
 integer, intent(in) :: src,dst

 xyzh(:,dst)  = xyzh(:,src)
 vxyzu(:,dst) = vxyzu(:,src)
 vpred(:,dst) = vpred(:,src)
 fxyzu(:,dst) = fxyzu(:,src)
 fext(:,dst)  = fext(:,src)
 if (mhd) then
    Bevol(:,dst)  = Bevol(:,src)
    Bpred(:,dst)  = Bpred(:,src)
    dBevol(:,dst) = dBevol(:,src)
    Bxyz(:,dst)   = Bxyz(:,src)
    divBsymm(dst) = divBsymm(src)
    if (maxmhdni==maxp) then
       n_R(:,dst)       = n_R(:,src)
       n_electronT(dst) = n_electronT(src)
       eta_nimhd(:,dst) = eta_nimhd(:,src)
    endif
 endif
 if (ndivcurlv > 0) divcurlv(:,dst) = divcurlv(:,src)
 if (ndivcurlB > 0) divcurlB(:,dst) = divcurlB(:,src)
 if (maxalpha ==maxp) alphaind(:,dst) = alphaind(:,src)
 if (maxgradh ==maxp) gradh(:,dst) = gradh(:,src)
 if (maxphase ==maxp) iphase(dst) = iphase(src)
 if (maxgrav  ==maxp) poten(dst) = poten(src)
 if (maxlum   ==maxp) luminosity(dst) = luminosity(src)
#ifdef IND_TIMESTEPS
 ibin(dst)       = ibin(src)
 ibin_old(dst)   = ibin_old(src)
 ibin_wake(dst)  = ibin_wake(src)
 dt_in(dst)      = dt_in(src)
 twas(dst)       = twas(src)
#endif
 if (use_dust) then
    dustfrac(:,dst)  = dustfrac(:,src)
    dustevol(:,dst)  = dustevol(:,src)
    dustpred(:,dst)  = dustpred(:,src)
    ddustfrac(:,dst) = ddustfrac(:,src)
    deltav(:,:,dst)  = deltav(:,:,src)
 endif
 if (maxp_h2==maxp) abundance(:,dst) = abundance(:,src)
 if (store_temperature) temperature(dst) = temperature(src)

 return
end subroutine copy_particle_all

!------------------------------------------------------------------
!+
! routine which reorders the particles according to an input list
! (prior to a derivs evaluation - so no derivs required)
! (allocates temporary arrays for each variable, so use with caution)
!+
!------------------------------------------------------------------
subroutine reorder_particles(iorder,np)
 integer, intent(in) :: iorder(:)
 integer, intent(in) :: np

 call copy_array(xyzh(:,1:np), iorder(1:np))
 call copy_array(vxyzu(:,1:np),iorder(1:np))
 call copy_array(fext(:,1:np), iorder(1:np))
 if (mhd) then
    call copy_array(Bevol(:,1:np),iorder(1:np))
    !--also copy the Bfield here, as this routine is used in setup routines
    call copy_array(Bxyz(:,1:np), iorder(1:np))
 endif
 if (ndivcurlv > 0)   call copy_arrayr4(divcurlv(:,1:np),iorder(1:np))
 if (maxalpha ==maxp) call copy_arrayr4(alphaind(:,1:np),iorder(1:np))
 if (maxgradh ==maxp) call copy_arrayr4(gradh(:,1:np),   iorder(1:np))
 if (maxphase ==maxp) call copy_arrayint1(iphase(1:np),  iorder(1:np))
 if (maxgrav  ==maxp) call copy_array1(poten(1:np),      iorder(1:np))
#ifdef IND_TIMESTEPS
 call copy_arrayint1(ibin(1:np),      iorder(1:np))
 call copy_arrayint1(ibin_old(1:np),  iorder(1:np))
 call copy_arrayint1(ibin_wake(1:np), iorder(1:np))
 !call copy_array1(twas(1:np),          iorder(1:np))
#endif

 return
end subroutine reorder_particles

!-----------------------------------------------------------------------
!+
!  routine to compactify the list of particles by removing dead
!  particles from the list
!  (could be openMP parallel if we sent in ndead)
!+
!-----------------------------------------------------------------------
subroutine shuffle_part(np)
 use io, only:fatal
 use domain, only:ibelong
 integer, intent(inout) :: np
 integer :: newpart

 do while (ideadhead /= 0)
    newpart = ideadhead
    if (newpart <= np) then
       if (.not.isdead(np)) then
          ! move particle to new position
          call copy_particle_all(np,newpart)
          ! move ibelong to new position
          ibelong(newpart) = ibelong(np)
          ! update deadhead
          ideadhead = ll(newpart)
       endif
       np = np - 1
    else
       ideadhead = ll(newpart)
    endif
    if (np < 0) call fatal('shuffle','npart < 0')
 enddo

 return
end subroutine shuffle_part

integer function count_dead_particles()
 integer :: i

 i = ideadhead
 count_dead_particles = 0
 do while (i > 0)
    count_dead_particles = count_dead_particles + 1
    i = ll(i)
 enddo

end function count_dead_particles

!-----------------------------------------------------------------------
!+
!  routine to remove dead or accreted particles
!  uses the routines above for efficiency
!+
!-----------------------------------------------------------------------
subroutine delete_dead_or_accreted_particles(npart,npartoftype)
 integer, intent(inout) :: npart,npartoftype(:)
 integer :: i,itype

 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       ! get the type so we know how to decrement npartoftype
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
       else
          itype = igas
       endif
       npartoftype(itype) = npartoftype(itype) - 1
       call kill_particle(i)
    endif
 enddo
 call shuffle_part(npart)

 return
end subroutine delete_dead_or_accreted_particles

!----------------------------------------------------------------
!+
!   change the position and status of a dead particle
!
!+
!----------------------------------------------------------------

subroutine change_status_pos(npart,x,y,z,h,vx,vy,vz)

 integer, intent(in) :: npart
 real, intent (in) :: x,y,z,h
 real, intent (in) :: vx,vy,vz
 integer  :: i,ix

 ix=0

 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       ix=i
       exit
    endif
 enddo

 xyzh(1,ix)=x
 xyzh(2,ix)=y
 xyzh(3,ix)=z
 xyzh(4,ix)=h
 vxyzu(1,ix)=vx
 vxyzu(2,ix)=vy
 vxyzu(3,ix)=vz

 return

end subroutine change_status_pos

!----------------------------------------------------------------
!+
!  pack particle information into a contiguous buffer
!  to send to another processor
!+
!----------------------------------------------------------------
subroutine fill_sendbuf(i,xtemp)
 use io,       only:fatal
 use mpiutils, only:fill_buffer
 integer, intent(in)  :: i
 real,    intent(out) :: xtemp(ipartbufsize)
 integer :: nbuf
!
!--package particle information into one simple wrapper
!
 nbuf = 0
!--NB: could use MPI_PACK here...
 if (i > 0) then
    call fill_buffer(xtemp,xyzh(:,i),nbuf)
    call fill_buffer(xtemp,vxyzu(:,i),nbuf)
    call fill_buffer(xtemp,vpred(:,i),nbuf)
    call fill_buffer(xtemp,fxyzu(:,i),nbuf)
    call fill_buffer(xtemp,fext(:,i),nbuf)
    if (ndivcurlv > 0) then
       call fill_buffer(xtemp,divcurlv(1,i),nbuf)
    endif
    if (maxalpha==maxp) then
       call fill_buffer(xtemp,alphaind(:,i),nbuf)
    endif
    if (maxgradh==maxp) then
       call fill_buffer(xtemp,gradh(:,i),nbuf)
    endif
    if (mhd) then
       call fill_buffer(xtemp,Bevol(:,i),nbuf)
       call fill_buffer(xtemp,Bpred(:,i),nbuf)
    endif
    if (maxphase==maxp) then
       call fill_buffer(xtemp,iphase(i),nbuf)
    endif
    if (use_dust) then
       call fill_buffer(xtemp, dustfrac(:,i),nbuf)
       call fill_buffer(xtemp, dustevol(:,i),nbuf)
    endif
    if (maxp_h2==maxp) then
       call fill_buffer(xtemp, abundance(:,i),nbuf)
    endif
    if (store_temperature) then
       call fill_buffer(xtemp, temperature(i),nbuf)
    endif
    if (maxgrav==maxp) then
       call fill_buffer(xtemp, poten(i),nbuf)
    endif
#ifdef IND_TIMESTEPS
    call fill_buffer(xtemp,ibin(i),nbuf)
    call fill_buffer(xtemp,ibin_old(i),nbuf)
    call fill_buffer(xtemp,ibin_wake(i),nbuf)
    call fill_buffer(xtemp,dt_in(i),nbuf)
    call fill_buffer(xtemp,twas(i),nbuf)
#endif
 endif
 if (nbuf /= ipartbufsize) call fatal('fill_sendbuf','error in send buffer size')

 return
end subroutine fill_sendbuf

!----------------------------------------------------------------
!+
!  unpack particle information from the send buffer
!  after receiving from another processor
!+
!----------------------------------------------------------------
subroutine unfill_buffer(ipart,xbuf)
 use mpiutils, only:unfill_buf
 integer, intent(in) :: ipart
 real,    intent(in) :: xbuf(ipartbufsize)
 integer :: j

 j = 0
 xyzh(:,ipart)          = unfill_buf(xbuf,j,4)
 vxyzu(:,ipart)         = unfill_buf(xbuf,j,maxvxyzu)
 vpred(:,ipart)         = unfill_buf(xbuf,j,maxvxyzu)
 fxyzu(:,ipart)         = unfill_buf(xbuf,j,maxvxyzu)
 fext(:,ipart)          = unfill_buf(xbuf,j,3)
 if (ndivcurlv > 0) then
    divcurlv(1,ipart)  = real(unfill_buf(xbuf,j),kind=kind(divcurlv))
 endif
 if (maxalpha==maxp) then
    alphaind(:,ipart)   = real(unfill_buf(xbuf,j,nalpha),kind(alphaind))
 endif
 if (maxgradh==maxp) then
    gradh(:,ipart)      = real(unfill_buf(xbuf,j,ngradh),kind(gradh))
 endif
 if (mhd) then
    Bevol(:,ipart)      = real(unfill_buf(xbuf,j,maxBevol),kind=kind(Bevol))
    Bpred(:,ipart)      = real(unfill_buf(xbuf,j,maxBevol),kind=kind(Bevol))
 endif
 if (maxphase==maxp) then
    iphase(ipart)       = nint(unfill_buf(xbuf,j),kind=1)
 endif
 if (use_dust) then
    dustfrac(:,ipart)   = unfill_buf(xbuf,j,maxdusttypes)
    dustevol(:,ipart)   = unfill_buf(xbuf,j,maxdusttypes)
 endif
 if (maxp_h2==maxp) then
    abundance(:,ipart)  = unfill_buf(xbuf,j,nabundances)
 endif
 if (store_temperature) then
    temperature(ipart)  = unfill_buf(xbuf,j)
 endif
 if (maxgrav==maxp) then
    poten(ipart)        = real(unfill_buf(xbuf,j),kind=kind(poten))
 endif
#ifdef IND_TIMESTEPS
 ibin(ipart)            = nint(unfill_buf(xbuf,j),kind=1)
 ibin_old(ipart)        = nint(unfill_buf(xbuf,j),kind=1)
 ibin_wake(ipart)       = nint(unfill_buf(xbuf,j),kind=1)
 dt_in(ipart)           = real(unfill_buf(xbuf,j),kind=kind(dt_in))
 twas(ipart)            = unfill_buf(xbuf,j)
#endif

!--just to be on the safe side, set other things to zero
 if (mhd) then
    divBsymm(ipart) = 0.
 endif

 return
end subroutine unfill_buffer

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (rank 2 arrays)
!+
!----------------------------------------------------------------

subroutine copy_array(array,ilist)
 real,    intent(inout) :: array(:,:)
 integer, intent(in)    :: ilist(:)
 real :: arraytemp(size(array(1,:)))
 integer :: i

 do i=1,size(array(:,1))
    arraytemp(:) = array(i,ilist(:))
    array(i,:) = arraytemp
 enddo

 return
end subroutine copy_array

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 2 arrays)
!+
!----------------------------------------------------------------

subroutine copy_arrayr4(array,ilist)
 real(kind=4), intent(inout) :: array(:,:)
 integer,      intent(in)    :: ilist(:)
 real(kind=4) :: arraytemp(size(array(1,:)))
 integer :: i

 do i=1,size(array(:,1))
    arraytemp(:) = array(i,ilist(:))
    array(i,:) = arraytemp
 enddo

 return
end subroutine copy_arrayr4

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 1 arrays)
!+
!----------------------------------------------------------------

subroutine copy_array1(array,ilist)
 real(kind=4), intent(inout) :: array(:)
 integer,      intent(in)    :: ilist(:)
 real(kind=4) :: arraytemp(size(array(:)))

 arraytemp(:) = array(ilist(:))
 array = arraytemp

 return
end subroutine copy_array1

!----------------------------------------------------------------
!+
!  utility to reorder an array
!  (real4, rank 1 arrays)
!+
!----------------------------------------------------------------

subroutine copy_arrayint1(iarray,ilist)
 integer(kind=1), intent(inout) :: iarray(:)
 integer,         intent(in)    :: ilist(:)
 integer(kind=1) :: iarraytemp(size(iarray(:)))

 iarraytemp(:) = iarray(ilist(:))
 iarray = iarraytemp

 return
end subroutine copy_arrayint1

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined box
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_box(xmin, xmax, ymin, ymax, zmin, zmax)
 real, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

 integer :: i
 real :: x, y, z, h

 do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    h = xyzh(4,i)
    if (x  <  xmin .or. x  >  xmax .or. y  <  ymin .or. y  >  ymax .or. z  <  zmin .or. z  >  zmax) then
       xyzh(4,i) = -abs(h)
    endif
 enddo
end subroutine

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined sphere
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_sphere(center, radius)
 real, intent(in) :: center(3), radius

 integer :: i
 real :: r(3), radius_squared

 radius_squared = radius**2
 do i=1,npart
    r = xyzh(1:3,i) - center
    if (dot_product(r,r)  >  radius_squared) then
       xyzh(4,i) = -abs(xyzh(4,i))
    endif
 enddo
end subroutine

!----------------------------------------------------------------
!+
!  Delete particles outside of a defined cylinder
!+
!----------------------------------------------------------------
subroutine delete_particles_outside_cylinder(center, radius, zmax)
 real, intent(in) :: center(3), radius, zmax

 integer :: i
 real :: x, y, z, rcil

 do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    rcil=sqrt((x-center(1))**2+(y-center(2))**2)

    if (rcil>radius .or. abs(z)>zmax) then
       call kill_particle(i)
    endif
 enddo
end subroutine

!----------------------------------------------------------------
!+
!  Delete particles within radius
!+
!----------------------------------------------------------------
subroutine delete_particles_inside_radius(center,radius,npart,npartoftype)
 real, intent(in) :: center(3), radius
 integer, intent(inout) :: npart,npartoftype(:)

 integer :: i,itype
 real :: x,y,z,rcil

 do i=1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    rcil=sqrt((x-center(1))**2+(y-center(2))**2)

    if (rcil<radius) then
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
       else
          itype = igas
       endif
       npartoftype(itype) = npartoftype(itype) - 1
       call kill_particle(i)
    endif
 enddo
 call shuffle_part(npart)

 return
end subroutine

end module part
