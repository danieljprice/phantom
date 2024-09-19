!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!   - EOS_file : File containing tabulated EOS values
!   - Lstar    : Luminosity of host star for calculating Tmin (Lsun)
!
! :Dependencies: eos_stamatellos, infile_utils, io, part, physcon, units
!

 implicit none
 real  :: Lstar = 0d0 ! in units of L_sun
 real,parameter :: dtcool_crit = 0.0001 ! critical dt_rad/dt_hydro for not applying cooling
 integer :: isink_star ! index of sink to use as illuminating star
 integer :: od_method = 4 ! default = Modified Lombardi method (Young et al. 2024)
 integer :: fld_opt = 1 ! by default FLD is switched on
 public :: radcool_update_energ,write_options_cooling_radapprox,read_options_cooling_radapprox
 public :: init_star, radcool_update_energ_loop

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


!
! Do cooling calculation
!
! update energy to return evolved energy array. Called from substep
subroutine radcool_update_energ(i,xi,yi,zi,rhoi,ui,Tfloor,dt,dudti_cool)
 use io,       only:warning
 use physcon,  only:steboltz,pi,solarl,Rg,kb_on_mh,piontwo,rpiontwo
 use units,    only:umass,udist,unit_density,unit_ergg,utime,unit_pressure
 use eos_stamatellos, only:getopac_opdep,getintenerg_opdep,gradP_cool,Gpot_cool,&
          duFLD,doFLD,ttherm_store,teqi_store,opac_store,duSPH
 use part,       only:xyzmh_ptmass,igas
 integer,intent(in) :: i
 real,intent(in) :: xi,yi,zi,rhoi,Tfloor
 real,intent(in) :: ui,dt
 real,intent(out)::dudti_cool
 real            :: coldensi,kappaBari,kappaParti,ri2
 real            :: gmwi,Tmini4,Ti,dudti_rad,Teqi,Hstam,HLom,du_tot
 real            :: cs2,Om2,Hmod2
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
 call getopac_opdep(ui*unit_ergg,rhoi*unit_density,kappaBari,kappaParti,&
           Ti,gmwi)
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

 call getintenerg_opdep(Tmini4**(1.0/4.0),rhoi*unit_density,umini)
 umini = umini/unit_ergg

 opaci = (coldensi**2d0)*kappaBari + (1.d0/kappaParti) ! physical units
 opac_store(i) = opaci
 dudti_rad = 4.d0*steboltz*(Tmini4 - Ti**4.d0)/opaci/unit_ergg*utime! code units

 if (doFLD) then
    du_tot = duSPH(i) + du_FLDi
 else
    du_tot = duSPH(i)
 endif

 ! If radiative cooling is negligible compared to hydrodynamical heating
 ! don't use this method to update energy, just use hydro du/dt. Does it conserve u alright?

 if (abs(du_tot) > epsilon(du_tot) .and.  abs(dudti_rad/du_tot) < dtcool_crit) then
    !      print *, "not cooling/heating for r=",sqrt(ri2),".", dudti_rad,&
    !                 dusph(i)
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

 if (Teqi > 9e5) then
    print *,"i=",i, "duSPH(i)=", duSPH(i), "duradi=", dudti_rad, "Ti=", Ti, &
            "Tmini=", Tmini4**(1.0/4.0),du_tot,Hcomb, "r=",sqrt(ri2), "ui=", ui, &
            "dudt_sph * dti=", dusph(i)*dt
 elseif (Teqi < epsilon(Teqi)) then
    print *,  "Teqi=0.0", "Tmini4=", Tmini4, "coldensi=", coldensi, "Tfloor=",Tfloor,&
            "Ti=", Ti, "poti=",poti, "rhoi=", rhoi
 elseif (Teqi < Tfloor) then
    print *,  "Teqi=",Teqi, "Tmini4=", Tmini4, "coldensi=", coldensi, "Tfloor=",Tfloor,&
            "Ti=", Ti, "poti=",poti, "rhoi=", rhoi
 endif

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
 elseif ( (dt/tthermi) < TINY(ui) ) then
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
    print *, "dudti_rad=", dudti_rad ,"dudt_fld=",du_fldi,"ueqi=",ueqi,"ui=",ui
    call warning("In Stamatellos cooling","energ=NaN or 0. ui=",val=ui)
    stop
 endif

end subroutine radcool_update_energ


!
! Do cooling calculation
!
! update energy to return evolved energy array. Called from evolve.F90
subroutine radcool_update_energ_loop(dtsph,npart,xyzh,energ,dudt_sph,Tfloor)
 use io,       only:warning
 use physcon,  only:steboltz,pi,solarl,Rg,kb_on_mh,piontwo,rpiontwo
 use units,    only:umass,udist,unit_density,unit_ergg,utime,unit_pressure
 use eos_stamatellos, only:getopac_opdep,getintenerg_opdep,gradP_cool,Gpot_cool,&
          duFLD,doFLD,ttherm_store,teqi_store,opac_store
 use part,       only:xyzmh_ptmass,rhoh,massoftype,igas,iactive,isdead_or_accreted
 use part,       only:iphase,maxphase,maxp,iamtype,ibin
 use timestep_ind, only:get_dt
 integer,intent(in) :: npart
 real,intent(in) :: xyzh(:,:),dtsph,Tfloor
 real,intent(inout) :: energ(:),dudt_sph(:)
 real            :: ui,rhoi,coldensi,kappaBari,kappaParti,ri2,dti
 real            :: gmwi,Tmini4,Ti,dudti_rad,Teqi,Hstam,HLom,du_tot
 real            :: cs2,Om2,Hmod2
 real            :: opaci,ueqi,umini,tthermi,poti,presi,Hcomb,du_FLDi
 integer         :: i,ratefile,n_uevo

 coldensi = huge(coldensi)
! write (temp,'(E5.2)') dt
 print *, "radcool min/maxGpot", minval(Gpot_cool),maxval(Gpot_cool)
 print *, "radcool min/max", minval(gradP_cool),maxval(gradP_cool)
 n_uevo = 0
 !$omp parallel do default(none) schedule(runtime) &
 !$omp shared(npart,duFLD,xyzh,energ,massoftype,xyzmh_ptmass,unit_density,Gpot_cool) &
 !$omp shared(isink_star,doFLD,ttherm_store,teqi_store,od_method,unit_pressure,ratefile) &
 !$omp shared(opac_store,Tfloor,dtsph,dudt_sph,utime,udist,umass,unit_ergg,gradP_cool,Lstar) &
 !$omp private(i,poti,du_FLDi,ui,rhoi,ri2,coldensi,kappaBari,Ti,iphase) &
 !$omp private(kappaParti,gmwi,Tmini4,dudti_rad,Teqi,Hstam,HLom,du_tot) &
 !$omp private(cs2,Om2,Hmod2,opaci,ueqi,umini,tthermi,presi,Hcomb,dti) &
 !$omp shared(maxp,maxphase,ibin) reduction(+:n_uevo)

 overpart: do i=1,npart
    if (maxphase==maxp) then
       if (iamtype(iphase(i)) /= igas) cycle
       if (isdead_or_accreted(xyzh(4,i))) cycle
       if (.not. iactive(iphase(i)) ) then
          n_uevo = n_uevo + 1
          cycle
       endif
    endif

    dti = get_dt(dtsph,ibin(i))
    poti = Gpot_cool(i)
    du_FLDi = duFLD(i)
    ui = energ(i)
    if (abs(ui) < epsilon(ui)) print *, "ui zero", i
    rhoi =  rhoh(xyzh(4,i),massoftype(igas))

    if (isink_star > 0) then
       ri2 = (xyzh(1,i)-xyzmh_ptmass(1,isink_star))**2d0 &
            + (xyzh(2,i)-xyzmh_ptmass(2,isink_star))**2d0 &
            + (xyzh(3,i)-xyzmh_ptmass(3,isink_star))**2d0
    else
       ri2 = xyzh(1,i)**2d0 + xyzh(2,i)**2d0 + xyzh(3,i)**2d0
    endif

    ! get opacities & Ti for ui
    call getopac_opdep(ui*unit_ergg,rhoi*unit_density,kappaBari,kappaParti,&
           Ti,gmwi)
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
       print *, "no case!"
       stop
    end select

!    Tfloor is from input parameters and is background heating
!    Stellar heating
    if (isink_star > 0 .and. Lstar > 0.d0) then
       Tmini4 = Tfloor**4d0 + exp(-coldensi*kappaBari)*(Lstar*solarl/(16d0*pi*steboltz*ri2*udist*udist))
    else
       Tmini4 = Tfloor**4d0
    endif

    opaci = (coldensi**2d0)*kappaBari + (1.d0/kappaParti) ! physical units
    opac_store(i) = opaci
    dudti_rad = 4.d0*steboltz*(Tmini4 - Ti**4.d0)/opaci/unit_ergg*utime! code units

    if (doFLD) then
       du_tot = dudt_sph(i) + du_FLDi
    else
       du_tot = dudt_sph(i)
    endif
    !If radiative cooling is negligible compared to hydrodynamical heating
    ! don't use this method to update energy, just use hydro du/dt
    if (abs(dudti_rad/du_tot) < dtcool_crit) then
       !      print *, "not cooling/heating for r=",sqrt(ri2),".", dudti_rad,&
       !                 dudt_sph(i)
       energ(i) = ui + du_tot*dti
       cycle
    endif

    Teqi = du_tot * opaci*unit_ergg/utime ! physical units
    du_tot = du_tot + dudti_rad
    Teqi = Teqi/4.d0/steboltz
    Teqi = Teqi + Tmini4
    if (Teqi < Tmini4) then
       Teqi = Tmini4**(1.0/4.0)
    else
       Teqi = Teqi**(1.0/4.0)
    endif
    teqi_store(i) = Teqi

    if (Teqi > 9e5) then
       print *,"i=",i, "dudt_sph(i)=", dudt_sph(i), "duradi=", dudti_rad, "Ti=", Ti, &
            "Tmini=", Tmini4**(1.0/4.0),du_tot,Hcomb, "r=",sqrt(ri2), "ui=", ui, &
            "dudt_sph * dti=", dudt_sph(i)*dti
    elseif (Teqi < epsilon(Teqi)) then
       print *,  "Teqi=0.0", "Tmini4=", Tmini4, "coldensi=", coldensi, "Tfloor=",Tfloor,&
            "Ti=", Ti, "poti=",poti, "rhoi=", rhoi
    endif

    call getintenerg_opdep(Teqi,rhoi*unit_density,ueqi)
    ueqi = ueqi/unit_ergg

    call getintenerg_opdep(Tmini4**(1.0/4.0),rhoi*unit_density,umini)
    umini = umini/unit_ergg

    ! calculate thermalization timescale
    if ((du_tot) == 0.d0) then
       tthermi = 0d0
    else
       tthermi = abs((ueqi - ui)/(du_tot))
    endif

    ttherm_store(i) = tthermi

    ! evolve energy
    if (tthermi == 0d0) then
       energ(i) = ui ! condition if denominator above is zero
    elseif ( (dti/tthermi) < TINY(ui) ) then
       energ(i) = ui
    else
       energ(i) = ui*exp(-dti/tthermi) + ueqi*(1.d0-exp(-dti/tthermi))  !code units
    endif

    if (isnan(energ(i)) .or. energ(i) < epsilon(ui)) then
       !    print *, "kappaBari=",kappaBari, "kappaParti=",kappaParti
       print *, "rhoi=",rhoi*unit_density, "Ti=", Ti, "Teqi=", Teqi
       print *, "Tmini=",Tmini4**(1.0/4.0), "ri=", ri2**(0.5)
       print *, "opaci=",opaci,"coldensi=",coldensi,"dudt_sphi",dudt_sph(i)
       print *,  "dt=",dti,"tthermi=", tthermi,"umini=", umini
       print *, "dudti_rad=", dudti_rad ,"dudt_fld=",du_fldi,"ueqi=",ueqi,"ui=",ui
       call warning("In Stamatellos cooling","energ=NaN or 0. ui",val=ui)
       stop
    endif

 enddo overpart
 !$omp end parallel do

 print *, "radcool min/max u():", minval(energ(1:npart)), maxval(energ(1:npart))
 print *, "radcool min/max Teqi():", minval(Teqi_store(1:npart)), maxval(Teqi_store(1:npart))
end subroutine radcool_update_energ_loop


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
    if (ieosread /= 23) call fatal('ieosread','For icooling=9, you need ieos=23')
 case default
    imatch = .false.
 end select

 if (ngot >= 4) igotallstam = .true.

end subroutine read_options_cooling_radapprox

end module cooling_radapprox

