!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
 integer, public :: od_method = 4 ! default = Modified Lombardi method (Young et al. 2024)
 real  :: Lstar = 0d0 ! in units of L_sun
 integer :: isink_star ! index of sink to use as illuminating star
 integer :: fld_opt = 1 ! by default FLD is switched on
 public :: radcool_update_du,write_options_cooling_radapprox,read_options_cooling_radapprox
 public :: init_star,radcool_evolve_ui

contains

subroutine init_star()
 use part,    only:nptmass,xyzmh_ptmass
 use io,      only:fatal
 integer :: i,imin=0
 real :: rsink2,rsink2min

 rsink2min = 0d0

 isink_star = 0
 if (od_method == 4 .and. nptmass == 0) then
    print *, "NO central star and using od_method = 4"
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
 ! update energy to return evolved energy. Called from substep.
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
    if (dt > 0d0) then
       ! evolve energy
       expdtonttherm = exp(-dt/tthermi)
       utemp = ui*expdtonttherm + ueqi*(1.d0-expdtonttherm)
    elseif (dt < 0d0) then
       ! i.e. for the backwards step in the leapfrog integrator
       expdtonttherm = exp(dt/tthermi)
       utemp = (ui - ueqi*(1-expdtonttherm))/expdtonttherm
    endif

    ! if tthermi ==0 or dt/thermi is neglible then ui doesn't change
    if (isnan(utemp) .or. utemp < epsilon(utemp)) then
       utemp = ui
    endif
 endif
 if (utemp < ufloor_cgs/unit_ergg) utemp = ufloor_cgs/unit_ergg
 if (utemp < 0d0) print *, "ERROR! i=",i, ui,ueqi

 if (present(uout)) then
    uout = utemp
 else
    ui = utemp
 endif

end subroutine radcool_evolve_ui


!
! Calculate equilibrium energy and thermal timescale
subroutine radcool_update_du(i,xi,yi,zi,rhoi,ui,duhydro,Tfloor)
 use io,       only:warning
 use physcon,  only:steboltz,pi,solarl,kb_on_mh,piontwo,rpiontwo
 use units,    only:umass,udist,unit_density,unit_ergg,utime,unit_pressure
 use eos_stamatellos, only:getopac_opdep,getintenerg_opdep,gradP_cool,Gpot_cool,&
          duFLD,doFLD,ttherm_store,ueqi_store,tau_store,du_store
 use part,       only:xyzmh_ptmass,igas,eos_vars,iTemp
 integer,intent(in) :: i
 real,intent(in) :: xi,yi,zi,rhoi
 real,intent(in) :: ui,duhydro,Tfloor
 real            :: coldensi,kappaBari,kappaParti,ri2
 real            :: gmwi,Tmini4,Ti,dudti_rad,Teqi,Hstam,HLom,du_tot
 real            :: cs2,Om2,Hmod2,rhoi_cgs,ui_cgs
 real            :: opaci,ueqi,umini,tthermi,poti,presi,Hcomb,du_FLDi

 coldensi = huge(coldensi)
 poti = Gpot_cool(i)
 du_FLDi = duFLD(i)
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
 ui_cgs = ui*unit_ergg
 rhoi_cgs = rhoi*unit_density
 call getopac_opdep(ui_cgs,rhoi_cgs,kappaBari,kappaParti,&
           Ti,gmwi)
 eos_vars(iTemp,i) = Ti ! save temperature
 presi = kb_on_mh*rhoi*unit_density*Ti/gmwi ! cgs
 presi = presi/unit_pressure !code units

 select case (od_method)
 case (1)
    ! Stamatellos+ 2007 method
    coldensi = sqrt(abs(poti*rhoi)/4.d0/pi) ! G cancels out as G=1 in code
    coldensi = 0.368d0*coldensi ! n=2 in polytrope formalism Forgan+ 2009
    coldensi = coldensi*umass/udist/udist ! physical units
 case (2)
    ! Lombardi+ 2015 method of estimating the mean column density
    coldensi = 1.014d0 * presi / abs(gradP_cool(i))! 1.014d0 * P/(-gradP/rho)
    coldensi = coldensi *umass/udist/udist ! physical units
 case (3)
    ! Combined method
    HStam = sqrt(abs(poti*rhoi)/4.0d0/pi)*0.368d0/rhoi
    HLom  = 1.014d0*presi/abs(gradP_cool(i))/rhoi
    Hcomb = 1.d0/sqrt((1.0d0/HLom)**2.0d0 + (1.0d0/HStam)**2.0d0)
    coldensi = Hcomb*rhoi
    coldensi = coldensi*umass/udist/udist ! physical units
 case (4)
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
 case default
    call warning("In radapprox cooling","cooling method not recognised",ival=od_method)
    return
 end select

!    Tfloor is from input parameters and is background heating
!    Stellar heating
 if (isink_star > 0 .and. Lstar > 0.d0) then
    Tmini4 = Tfloor**4d0 + exp(-coldensi*kappaBari)*(Lstar*solarl/(16d0*pi*steboltz*ri2*udist*udist))
 else
    Tmini4 = Tfloor**4d0
 endif

 call getintenerg_opdep(Tmini4**(1.0/4.0),rhoi_cgs,umini)
 umini = umini/unit_ergg

 opaci = (coldensi**2d0)*kappaBari + (1.d0/kappaParti) ! physical units
 tau_store(i) = coldensi * kappaBari
 dudti_rad = 4.d0*steboltz*(Tmini4 - Ti**4.d0)/opaci/unit_ergg*utime! code units

 if (doFLD) then
    du_tot = duhydro + du_FLDi
 else
    du_tot = duhydro
 endif

 du_store(i) = dudti_rad
 Teqi = du_tot * opaci*unit_ergg/utime ! physical units
 Teqi = Teqi/4.d0/steboltz
 Teqi = Teqi + Tmini4
 du_tot = du_tot + dudti_rad
 !Check if we need to use the temperature floor
 if (Teqi < Tmini4) then
    Teqi = Tmini4**(1.0/4.0)
 else
    Teqi = Teqi**(1.0/4.0)
 endif

 if (Teqi > 9e5) then
    print *,"i=",i, "duhydro=", duhydro, "duradi=", dudti_rad, "Ti=", Ti, &
            "Tmini=", Tmini4**(1.0/4.0),du_tot,Hcomb, "r=",sqrt(ri2), "ui=", ui
 elseif (Teqi < epsilon(Teqi)) then
    print *,  "Teqi=0.0", "Tmini4=", Tmini4, "coldensi=", coldensi, "Tfloor=",Tfloor,&
            "Ti=", Ti, "poti=",poti, "rhoi=", rhoi
 elseif (Teqi < Tfloor) then
    print *,  "Teqi=",Teqi, "Tmini4=", Tmini4, "coldensi=", coldensi, "Tfloor=",Tfloor,&
            "Ti=", Ti, "poti=",poti, "rhoi=", rhoi
 endif

 rhoi_cgs = rhoi*unit_density
 call getintenerg_opdep(Teqi,rhoi_cgs,ueqi)
 ueqi = ueqi/unit_ergg
 ueqi_store(i) = ueqi
 ! calculate thermalization timescale
 if ((du_tot) == 0.d0) then
    tthermi = 0d0
 else
    tthermi = abs((ueqi - ui)/(du_tot))
 endif

 ttherm_store(i) = tthermi


 if (isnan(tthermi) .or. isnan(ueqi)) then
    !    print *, "kappaBari=",kappaBari, "kappaParti=",kappaParti
    print *, "rhoi=",rhoi*unit_density, "Ti=", Ti, "Teqi=", Teqi
    print *, "Tmini=",Tmini4**(1.0/4.0), "ri=", ri2**(0.5)
    print *, "opaci=",opaci,"coldensi=",coldensi,"duhydro",duhydro
    print *,  "tthermi=", tthermi,"umini=", umini
    print *, "dudti_rad=", dudti_rad ,"dudt_fld=",du_fldi,"ueqi=",ueqi,"ui=",ui
    call warning("In Stamatellos cooling","energ=NaN or 0. ui=",val=ui)
    stop
 endif

end subroutine radcool_update_du


subroutine write_options_cooling_radapprox(iunit)
 use infile_utils, only:write_inopt
 use eos_stamatellos, only: eos_file
 integer, intent(in) :: iunit

 !N.B. Tfloor handled in cooling.F90
 call write_inopt(eos_file,'EOS_file','File containing tabulated EOS values',iunit)
 call write_inopt(od_method,'OD method',&
      'Method for estimating optical depth:(1)Stamatellos (2)Lombardi (3)combined (4)modified Lombardi',iunit)
 call write_inopt(Lstar,'Lstar','Luminosity of host star for calculating Tmin (Lsun)',iunit)
 call write_inopt(FLD_opt,'do FLD','Do FLD? (1) yes (0) no',iunit)

end subroutine write_options_cooling_radapprox

subroutine read_options_cooling_radapprox(name,valstring,imatch,igotallstam,ierr)
 use io, only:warning,fatal
 use eos_stamatellos, only: eos_file,doFLD
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotallstam
 integer,          intent(out) :: ierr
 integer       :: ieosread
 integer, save :: ngot = 0

 imatch  = .true.
 igotallstam = .false. ! cooling options are compulsory
 select case(trim(name))
 case('Lstar')
    read(valstring,*,iostat=ierr) Lstar
    if (Lstar < 0.) call fatal('Lstar','Luminosity cannot be negative')
    ngot = ngot + 1
 case('OD method')
    read(valstring,*,iostat=ierr) od_method
    if (od_method < 1 .or. od_method > 4) then
       call fatal('cooling options','od_method must be 1, 2, 3 or 4',var='od_method',ival=od_method)
    endif
    ngot = ngot + 1
 case('EOS_file')
    read(valstring,*,iostat=ierr) eos_file
    ngot = ngot + 1
 case('do FLD')
    read(valstring,*,iostat=ierr) FLD_opt
    if (FLD_opt < 0) call fatal('FLD_opt','FLD option out of range')
    if (FLD_opt == 0) then
       doFLD = .false.
    elseif (FLD_opt == 1) then
       doFLD = .true.
    endif
    ngot = ngot + 1
 case('ieos')
    read(valstring,*,iostat=ierr) ieosread
    if (ieosread /= 24) call fatal('ieosread','For icooling=9, you need ieos=24')
 case default
    imatch = .false.
 end select

 if (ngot >= 4) igotallstam = .true.

end subroutine read_options_cooling_radapprox

end module cooling_radapprox

