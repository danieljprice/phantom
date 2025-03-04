!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
 real  :: Lstar = 0. ! in units of L_sun
 integer :: isink_star ! index of sink to use as illuminating star
 public :: radcool_update_du,write_options_cooling_radapprox,read_options_cooling_radapprox
 public :: init_star,radcool_evolve_ui

contains

subroutine init_star()
 use part,    only:nptmass,xyzmh_ptmass
 use io,      only:fatal
 integer :: i,imin=0
 real :: rsink2,rsink2min

 rsink2min = 0.

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

subroutine radcool_evolve_ui(ui,dt,i,Tfloor,h,uout)
 use eos_stamatellos, only:ttherm_store,ueqi_store,getintenerg_opdep
 use io,              only:warning
 use units,           only:unit_density,unit_ergg
 use part,            only:rhoh,massoftype,igas
 real, intent(inout) :: ui
 real, intent(in)    :: dt,Tfloor,h
 integer,intent(in)  :: i
 real,optional,intent(out) :: uout
 real :: tthermi,ueqi,utemp,ufloor_cgs,rhoi_cgs
 real :: expdtonttherm

 tthermi = ttherm_store(i)
 ueqi = ueqi_store(i)
 utemp = ui
 rhoi_cgs = rhoh(h,massoftype(igas))*unit_density
 call getintenerg_opdep(Tfloor**(1.0/4.0),rhoi_cgs,ufloor_cgs)

 if (tthermi > epsilon(tthermi) .and. ui /= ueqi) then
    if (dt > 0.) then
       ! evolve energy
       expdtonttherm = exp(-dt/tthermi)
       utemp = ui*expdtonttherm + ueqi*(1.-expdtonttherm)
    elseif (dt < 0.) then
       ! i.e. for the backwards step in the leapfrog integrator
       expdtonttherm = exp(dt/tthermi)
       utemp = (ui - ueqi*(1.-expdtonttherm))/expdtonttherm
    endif

    ! if tthermi ==0 or dt/thermi is neglible then ui doesn't change
    if (isnan(utemp) .or. utemp < epsilon(utemp)) then
       utemp = ui
    endif
 endif
 if (utemp < ufloor_cgs/unit_ergg) utemp = ufloor_cgs/unit_ergg
 if (utemp < 0.) print *, "ERROR! i=",i, ui,ueqi

 if (present(uout)) then
    uout = utemp
 else
    ui = utemp
 endif

end subroutine radcool_evolve_ui


!
! Do cooling calculation
!
! update tthermi and ueqi for ui update
subroutine radcool_update_du(i,xi,yi,zi,rhoi,ui,duhydro,Tfloor)
 use io,       only:warning
 use physcon,  only:steboltz,pi,solarl,kb_on_mh,piontwo,rpiontwo
 use units,    only:umass,udist,unit_density,unit_ergg,utime,unit_pressure
 use eos_stamatellos, only:getopac_opdep,getintenerg_opdep,gradP_cool,&
          ttherm_store,ueqi_store,opac_store
 use part,       only:xyzmh_ptmass,igas,eos_vars,iTemp
 integer,intent(in) :: i
 real,intent(in) :: xi,yi,zi,rhoi
 real,intent(in) :: ui,duhydro,Tfloor
 real            :: coldensi,kappaBari,kappaParti,ri2
 real            :: gmwi,Tmini4,Ti,dudti_rad,Teqi,HLom,du_tot
 real            :: opaci,ueqi,umini,tthermi,presi,Hcomb
 real            :: cs2,Om2,Hmod2,rhoi_cgs,ui_cgs

 coldensi = huge(coldensi)
 kappaBari = 0.
 kappaParti = 0.
 Teqi = huge(Teqi)
 tthermi = huge(tthermi)
 opaci = epsilon(opaci)
 if (abs(ui) < epsilon(ui)) print *, "ui zero", i

 if (isink_star > 0) then
    ri2 = (xi-xyzmh_ptmass(1,isink_star))**2 &
          + (yi-xyzmh_ptmass(2,isink_star))**2 &
          + (zi-xyzmh_ptmass(3,isink_star))**2
 else
    ri2 = xi**2 + yi**2 + zi**2
 endif

 ! get opacities & Ti for ui
 ui_cgs = ui*unit_ergg
 rhoi_cgs = rhoi*unit_density
 call getopac_opdep(ui_cgs,rhoi_cgs,kappaBari,kappaParti,&
           Ti,gmwi)
 eos_vars(iTemp,i) = Ti ! save temperature
 presi = kb_on_mh*rhoi*unit_density*Ti/gmwi ! cgs
 presi = presi/unit_pressure !code units

! Modified Lombardi method
 HLom  = presi/abs(gradP_cool(i))/rhoi
 cs2 = presi/rhoi
 if (isink_star > 0 .and. ri2 > 0.) then
    Om2 = xyzmh_ptmass(4,isink_star)/(ri2**(1.5)) !NB we are using spherical radius here
 else
    Om2 = 0.
 endif
 Hmod2 = cs2 * piontwo / (Om2 + 8.*rpiontwo*rhoi)
 Hcomb = 1./sqrt((1./HLom)**2 + (1./Hmod2))
 coldensi = 1.014 * Hcomb *rhoi*umass/udist/udist ! physical units

!    Tfloor is from input parameters and is background heating
!    Stellar heating
 if (isink_star > 0 .and. Lstar > 0.) then
    Tmini4 = Tfloor**4 + exp(-coldensi*kappaBari)*(Lstar*solarl/(16.*pi*steboltz*ri2*udist*udist))
 else
    Tmini4 = Tfloor**4
 endif

 call getintenerg_opdep(Tmini4**(1.0/4.0),rhoi_cgs,umini)
 umini = umini/unit_ergg
 opaci = (coldensi**2)*kappaBari + (1./kappaParti) ! physical units
 opac_store(i) = opaci
 dudti_rad = 4.*steboltz*(Tmini4 - Ti**4)/opaci/unit_ergg*utime! code units

 du_tot = duhydro

 Teqi = du_tot * opaci*unit_ergg/utime ! physical units
 Teqi = Teqi/4./steboltz
 Teqi = Teqi + Tmini4
 du_tot = du_tot + dudti_rad
 !Check if we need to use the temperature floor
 if (Teqi < Tmini4) then
    Teqi = Tmini4**(1.0/4.0)
 else
    Teqi = Teqi**(1.0/4.0)
 endif

 rhoi_cgs = rhoi*unit_density
 call getintenerg_opdep(Teqi,rhoi_cgs,ueqi)
 ueqi = ueqi/unit_ergg
 ueqi_store(i) = ueqi
 ! calculate thermalization timescale
 if ((du_tot) == 0.) then
    tthermi = 0.
 else
    tthermi = abs((ueqi - ui)/(du_tot))
 endif

 ttherm_store(i) = tthermi

 if (isnan(tthermi) .or. isnan(ueqi)) then
    call warning("In Stamatellos cooling","energ=NaN or 0. ui=",val=ui)
 endif

end subroutine radcool_update_du


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
    if (ieosread /= 24) call fatal('ieosread','For icooling=9, you need ieos=24')
 case default
    imatch = .false.
 end select

 if (ngot >= 2) igotallra = .true.

end subroutine read_options_cooling_radapprox

end module cooling_radapprox

