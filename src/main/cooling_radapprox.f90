!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_radapprox
!
! Cooling method of Stamatellos et al. 2007
!
! :References: Stamatellos et al. 2007
!
! :Owner: Alison Young
!
! :Runtime parameters:
!   - EOS_file : *File containing tabulated EOS values*
!   - Lstar    : *Luminosity of host star for calculating Tmin (Lsun)*
!
! :Dependencies: eos_stamatellos, infile_utils, io, part, physcon, units
!

 implicit none
 real  :: Lstar = 0d0 ! in units of L_sun
 real,parameter :: dtcool_crit = 0.0001 ! critical dt_rad/dt_hydro for not applying cooling
 integer :: isink_star ! index of sink to use as illuminating star
 public :: radcool_update_energ,write_options_cooling_radapprox,read_options_cooling_radapprox
 public :: init_star

contains

subroutine init_star()
 use part,    only:nptmass,xyzmh_ptmass
 use io,      only:fatal
 integer :: i,imin=0
 real :: rsink2,rsink2min

 rsink2min = 0d0

 isink_star = 0
 if (nptmass == 0) then
    print *, "NO central star"
 elseif (nptmass == 0) then
    print *, "No stellar heating."
 elseif (nptmass == 1) then
    isink_star = 1
 else
    do i=1,nptmass
       rsink2 = xyzmh_ptmass(1,i)**2 + xyzmh_ptmass(2,i)**2 + xyzmh_ptmass(3,i)**2
       if (i==1 .or. (rsink2 < rsink2min) ) then
          rsink2min = rsink2
          imin = i
       endif
    enddo
    isink_star = imin
 endif
 if (isink_star > 0)  print *, "Using sink no. ", isink_star,&
      "at (xyz)",xyzmh_ptmass(1:3,isink_star)!"as illuminating star."
end subroutine init_star


!
! Do cooling calculation
!
! update energy to return evolved energy array. Called from substep
subroutine radcool_update_energ(i,xi,yi,zi,rhoi,ui,Tfloor,dt,dudti_cool)
 use io,       only:warning
 use physcon,  only:steboltz,pi,solarl,Rg,kb_on_mh,piontwo,rpiontwo
 use units,    only:umass,udist,unit_density,unit_ergg,utime,unit_pressure
 use eos_stamatellos, only:getopac_opdep,getintenerg_opdep,gradP_cool,&
          ttherm_store,teqi_store,opac_store,duSPH
 use part,       only:xyzmh_ptmass,igas
 integer,intent(in) :: i
 real,intent(in) :: xi,yi,zi,rhoi,Tfloor
 real,intent(in) :: ui,dt
 real,intent(out)::dudti_cool
 real            :: coldensi,kappaBari,kappaParti,ri2
 real            :: gmwi,Tmini4,Ti,dudti_rad,Teqi,HLom,du_tot
 real            :: cs2,Om2,Hmod2
 real            :: opaci,ueqi,umini,tthermi,presi,Hcomb

 coldensi = huge(coldensi)
 kappaBari = 0d0
 kappaParti = 0d0
 Teqi = huge(Teqi)
 tthermi = huge(tthermi)
 opaci = epsilon(opaci)
 if (abs(ui) < epsilon(ui)) print *, "ui zero", i

 if (isink_star > 0) then
    ri2 = (xi-xyzmh_ptmass(1,isink_star))**2d0 &
          + (yi-xyzmh_ptmass(2,isink_star))**2d0 &
          + (zi-xyzmh_ptmass(3,isink_star))**2d0
 else
    ri2 = xi**2d0 + yi**2d0 + zi**2d0
 endif

 ! get opacities & Ti for ui
 call getopac_opdep(ui*unit_ergg,rhoi*unit_density,kappaBari,kappaParti,&
           Ti,gmwi)
 presi = kb_on_mh*rhoi*unit_density*Ti/gmwi ! cgs
 presi = presi/unit_pressure !code units

! Modified Lombardi method
 HLom  = presi/abs(gradP_cool(i))/rhoi
 cs2 = presi/rhoi
 if (isink_star > 0 .and. ri2 > 0d0) then
    Om2 = xyzmh_ptmass(4,isink_star)/(ri2**(1.5)) !NB we are using spherical radius here
 else
    Om2 = 0d0
 endif
 Hmod2 = cs2 * piontwo / (Om2 + 8d0*rpiontwo*rhoi)
 Hcomb = 1.d0/sqrt((1.0d0/HLom)**2.0d0 + (1.0d0/Hmod2))
 coldensi = 1.014d0 * Hcomb *rhoi*umass/udist/udist ! physical units

!    Tfloor is from input parameters and is background heating
!    Stellar heating
 if (isink_star > 0 .and. Lstar > 0.d0) then
    Tmini4 = Tfloor**4d0 + exp(-coldensi*kappaBari)*(Lstar*solarl/(16d0*pi*steboltz*ri2*udist*udist))
 else
    Tmini4 = Tfloor**4d0
 endif

 call getintenerg_opdep(Tmini4**(1.0/4.0),rhoi*unit_density,umini)
 umini = umini/unit_ergg
 opaci = (coldensi**2d0)*kappaBari + (1.d0/kappaParti) ! physical units
 opac_store(i) = opaci
 dudti_rad = 4.d0*steboltz*(Tmini4 - Ti**4.d0)/opaci/unit_ergg*utime! code units

 du_tot = duSPH(i)


 ! If radiative cooling is negligible compared to hydrodynamical heating
 ! don't use this method to update energy, just use hydro du/dt. Does it conserve u alright?

 if (abs(du_tot) > epsilon(du_tot) .and.  abs(dudti_rad/du_tot) < dtcool_crit) then
    dudti_cool = du_tot
    if ( (dudti_cool*dt + ui) < umini) then
       dudti_cool = (umini - ui)/dt
    endif
    return
 endif

 Teqi = du_tot * opaci*unit_ergg/utime ! physical units
 Teqi = Teqi/4.d0/steboltz
 Teqi = Teqi + Tmini4
 du_tot = du_tot + dudti_rad
 if (Teqi < Tmini4) then
    Teqi = Tmini4**(1.0/4.0)
 else
    Teqi = Teqi**(1.0/4.0)
 endif
 teqi_store(i) = Teqi

 call getintenerg_opdep(Teqi,rhoi*unit_density,ueqi)
 ueqi = ueqi/unit_ergg

 ! calculate thermalization timescale
 if ((du_tot) == 0.d0) then
    tthermi = 0d0
 else
    tthermi = abs((ueqi - ui)/(du_tot))
 endif

 ttherm_store(i) = tthermi

 ! evolve energy
 if (tthermi == 0d0) then
    dudti_cool = 0d0 ! condition if denominator above is zero
 elseif ( (dt/tthermi) < tiny(ui) ) then
    dudti_cool = 0d0
 else
    dudti_cool = ( ui*exp(-dt/tthermi) + ueqi*(1.d0-exp(-dt/tthermi)) - ui) / dt !code units
 endif

 if (isnan(dudti_cool)) then
    !    print *, "kappaBari=",kappaBari, "kappaParti=",kappaParti
    print *, "rhoi=",rhoi*unit_density, "Ti=", Ti, "Teqi=", Teqi
    print *, "Tmini=",Tmini4**(1.0/4.0), "ri=", ri2**(0.5)
    print *, "opaci=",opaci,"coldensi=",coldensi,"dusph(i)",duSPH(i)
    print *,  "dt=",dt,"tthermi=", tthermi,"umini=", umini
    print *, "dudti_rad=", dudti_rad ,"ueqi=",ueqi,"ui=",ui
    call warning("In Stamatellos cooling","energ=NaN or 0. ui=",val=ui)
    stop
 endif

end subroutine radcool_update_energ


subroutine write_options_cooling_radapprox(iunit)
 use infile_utils, only:write_inopt
 use eos_stamatellos, only: eos_file
 integer, intent(in) :: iunit

 !N.B. Tfloor handled in cooling.F90
 call write_inopt(eos_file,'EOS_file','File containing tabulated EOS values',iunit)
 call write_inopt(Lstar,'Lstar','Luminosity of host star for calculating Tmin (Lsun)',iunit)

end subroutine write_options_cooling_radapprox


subroutine read_options_cooling_radapprox(name,valstring,imatch,igotallra,ierr)
 use io, only:warning,fatal
 use eos_stamatellos, only: eos_file
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotallra
 integer,          intent(out) :: ierr
 integer       :: ieosread
 integer, save :: ngot = 0

 imatch  = .true.
 igotallra = .false. ! cooling options are compulsory
 select case(trim(name))
 case('Lstar')
    read(valstring,*,iostat=ierr) Lstar
    if (Lstar < 0.) call fatal('Lstar','Luminosity cannot be negative')
    ngot = ngot + 1
 case('EOS_file')
    read(valstring,*,iostat=ierr) eos_file
    ngot = ngot + 1
 case('ieos')
    read(valstring,*,iostat=ierr) ieosread
    if (ieosread /= 23) call fatal('ieosread','For icooling=9, you need ieos=23')
 case default
    imatch = .false.
 end select

 if (ngot >= 2) igotallra = .true.

end subroutine read_options_cooling_radapprox

end module cooling_radapprox

