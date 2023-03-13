!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cooling_stamatellos
!
! Cooling method of Stamatellos et al. 2007
!
! :References: Stamatellos et al. 2007
!
! :Owner: Alison Young
!
! :Runtime parameters: None
!
! :Dependencies:
!
 
 implicit none
 real, public :: Lstar ! in units of L_sun
 public :: cooling_S07,write_options_cooling_stamatellos,read_options_cooling_stamatellos
 
 contains

!
! Do cooling calculation
!
   subroutine cooling_S07(rhoi,ui,dudti_cool,xi,yi,zi,Tfloor,dudti_sph,dt,i)
     use io,       only:warning
     use physcon,  only:steboltz,pi,solarl
     use units,    only:umass,udist,unit_density,unit_ergg,utime
     use eos_stamatellos, only:getopac_opdep,getintenerg_opdep,gradP_cool,Gpot_cool, &
          iunitst
     use part,       only:eos_vars,igasP
     real,intent(in) :: rhoi,ui,dudti_sph,xi,yi,zi,Tfloor,dt
     integer,intent(in) :: i
     real,intent(out) :: dudti_cool
     real            :: coldensi,kappaBari,kappaParti,Lstar,ri2
     real            :: Tirri,gammai,gmwi,Tmini,Ti,dudt_rad,Teqi
     real            :: tcool,ueqi,umini,tthermi,poti,presi,coldens2i
     
     poti = Gpot_cool(i)
     presi = eos_vars(igasP,i)
!    Tfloor is from input parameters and is background heating
!    Add to stellar heating. Just assuming one star at (0,0,0) for now
     Lstar = 0.0 !in Lsun
     ri2 = xi*xi + yi*yi + zi*zi
     ri2 = ri2 *udist*udist
     Tirri = (Lstar*solarl/(4d0*pi*steboltz)/ri2)**0.25
     Tmini = Tfloor + Tirri

! get opacities & Ti for ui

     call getopac_opdep(ui*unit_ergg,rhoi*unit_density,kappaBari,kappaParti,&
           Ti,gmwi,gammai)
     
     coldensi = sqrt(abs(poti*rhoi)/4.d0/pi)
     coldensi = 0.368d0*coldensi ! n=2 in polytrope formalism Forgan+ 2009
     coldensi = coldensi*umass/udist/udist ! physical units

! testing Lombardi+ method of estimating the mean column density
     coldens2i = 1.014d0 * presi / abs(-gradP_cool(i))! /rhoi ) ! Lombardi+ 2015
     coldens2i = coldens2i *umass/udist/udist ! physical units
!     write(iunitst,'(5E12.5)') coldensi,coldens2i,presi,gradP_cool(i)
! ! USE LOMBARDI METHOD !!
     coldensi = coldens2i
     !!!
     
     tcool = (coldensi**2d0)*kappaBari +(1.d0/kappaParti) ! physical units
     dudt_rad = 4.d0*steboltz*(Tmini**4.d0 - Ti**4.d0)/tcool/unit_ergg*utime! code units
! calculate Teqi
     Teqi = dudti_sph*(coldensi**2.d0*kappaBari + (1.d0/kappaParti))*unit_ergg/utime
     Teqi = Teqi/4.d0/steboltz
     Teqi = Teqi + Tmini**4.d0
     if (Teqi < Tmini**4.d0) then
        Teqi = Tmini
     else
        Teqi = Teqi**0.25d0
     endif
     call getintenerg_opdep(Teqi,rhoi*unit_density,ueqi)
     ueqi = ueqi/unit_ergg
     call getintenerg_opdep(Tmini,rhoi*unit_density,umini)
     umini = umini/unit_ergg
! calculate thermalization timescale
     if ((dudti_sph + dudt_rad) == 0.d0) then
        tthermi = 0d0
     else
        tthermi = abs((ueqi - ui)/(dudti_sph + dudt_rad))
     endif
     
! internal energy update -> put in form where it'll work as dudtcool
     if (tthermi == 0d0) then
        dudti_cool = 0.d0 ! condition if denominator above is zero
     else
        dudti_cool = (ui*exp(-dt/tthermi) + ueqi*(1.d0-exp(-dt/tthermi)) -ui)/dt !code units
     endif
     
     if (isnan(dudti_cool)) then
        print *, "kappaBari=",kappaBari, "kappaParti=",kappaParti
        print *, "poti=",poti, "rhoi=",rhoi, "Ti=", Ti
        print *, "tcool=",tcool,"coldensi=",coldensi,"dudti_sph",dudti_sph
        print *, "Teqi=",Teqi, "dt=",dt,"tthermi=", tthermi,"ueqi=", ueqi
        call warning("In Stamatellos cooling","dudticool=NaN. ui",val=ui)
        stop
     else if (dudti_cool < 0.d0 .and. abs(dudti_cool) > ui/dt) then
     	dudti_cool = (umini - ui)/dt
     endif
     
   end subroutine cooling_S07


 subroutine write_options_cooling_stamatellos(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 !Tfloor handled in cooling.F90
 call write_inopt(Lstar,'Lstar','Luminosity of host star for calculating Tmin (Lsun)',iunit)

end subroutine write_options_cooling_stamatellos

 subroutine read_options_cooling_stamatellos(name,valstring,imatch,igotallstam,ierr)
! use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotallstam
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0


 imatch  = .true.
 igotallstam = .false. ! cooling options are compulsory
 select case(trim(name))
 case('Lstar')
    read(valstring,*,iostat=ierr) Lstar
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if (ngot >= 1) igotallstam = .true.

end subroutine read_options_cooling_stamatellos

  
 end module cooling_stamatellos
