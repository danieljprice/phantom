!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module porosity
!
! Contains routine for porosity evolution (growth, bouncing,
! fragmentation, compaction, disruption)
!
! :References:
!   Okuzumi et al. (1997), ApJ 752, 106
!   Garcia, Gonzalez (2020), MNRAS 493, 1788
!   Tatsuuma et Kataoka (2021), ApJ 913, 132
!   Michoulier & Gonzalez (2022), MNRAS 517, 3064
!
! :Owner: Stephane Michoulier
!
! :Runtime parameters:
!   - gammaft     : *Force to torque efficient of gas flow on dust*
!   - ibounce     : *bouncing (0=Off,1=On)*
!   - icompact    : *Compaction during fragmentation (ifrag > 0) (0=off,1=on)*
!   - idisrupt    : *disruption (0=Off,1=On)*
!   - iporosity   : *porosity (0=Off,1=On)*
!   - smonocgs    : *Monomer size in cm (smaller or equal to 1.e-4 cm)*
!   - surfenergSI : *Monomer surface energy in J/m**2*
!   - youngmodSI  : *Monomer young modulus in Pa*
!
! :Dependencies: dim, dust, eos, growth, infile_utils, io, options, part,
!   physcon, random, units, viscosity
!
 use units,        only:umass,udist,unit_energ,unit_pressure,unit_density
 use physcon,      only:Ro,pi,fourpi,roottwo
 implicit none

 !--Default values

 integer, public        :: iporosity    = 0              !--0=Off  1=On    (-1=On for checkup, filfac is initialized but does not evolve)
 integer, public        :: icompact     = 1              !--0=off  1=on    (Compaction of dust grain during fragmentation)
 integer, public        :: ibounce      = 0              !--0=off  1=on    (Allow dust grains to bounce)
 integer, public        :: idisrupt     = 0              !--0=off  1=on    (Rotational disruption by gas flow: Tatsuuma et al. 2021)
 real, public           :: smonocgs     = 1e-4           !--monomer size in cm
 real, public           :: surfenergSI  = 0.20          !--surface energy of monomers in SI: J/m**2 (here for Si: Kimura et al. 2020)
 real, public           :: youngmodSI   = 72e9           !--young modulus of monomers in SI: Pa (here for Si: Yamamoto et al. 2014)
 real, public           :: gammaft      = 0.1            !--force-to-torque efficiency (Tatsuuma et al. 2021)

 real, parameter        :: cratio       = -0.5801454844  !--common ratio for a power
 real, parameter        :: b_oku        = 0.15           !--parameter b (Okuzumi et al. 2012)
 real, parameter        :: maxpacking   = 0.74048        !--max sphere packing for hexagonal close packing

 real, public           :: smono                         !--monomer size
 real, public           :: mmono                         !--monomer mass
 real, public           :: surfenerg
 real, public           :: youngmod
 real                   :: eroll                         !--rolling
 real                   :: grainmassminlog
 real                   :: Yd0                           !test for compaction
 real                   :: Ydpow                         !test for compaction

 public                 :: get_filfac,init_filfac,get_disruption,get_probastick
 public                 :: init_porosity,print_porosity_info,write_options_porosity,read_options_porosity
 public                 :: write_porosity_setup_options,read_porosity_setup_options

 private

contains

!------------------------------------------------
!+
!  Initialise variables for computing porosity
!+
!------------------------------------------------
subroutine init_porosity(ierr)
 use io,                only:error
 use dust,              only:idrag,grainsizecgs,graindenscgs
 integer, intent(out)      :: ierr

 ierr = 0

 !--initialise variables in code units
 smono         = smonocgs / udist
 mmono         = fourpi/3 * (graindenscgs / unit_density) * smono**3
 surfenerg     = surfenergSI * udist * udist * 1000 / unit_energ
 youngmod      = youngmodSI * 10 / unit_pressure
 eroll         = 302.455974078*(surfenerg**5 * smono**4 / youngmod**2)**(1./3.)

 Yd0 = 9.5e6 *10/unit_pressure ! for water+silicate; 9.8e6 for water only
 Ydpow = 6.4    !for silicate+water, 4 for water only

 grainmassminlog = log10(50.*mmono)

 if (smono <= 0.) then
    call error('init_porosity','smonocgs <= 0',var='smonocgs',val=smonocgs)
    ierr = 1
 endif

 if (grainsizecgs < smonocgs) then
    call error('init_porosity','grainsizecgs < smonocgs',var='smonocgs',val=smonocgs)
    ierr = 1
 endif

 if (surfenerg <= 0.) then
    call error('init_porosity','surfenerg <= 0',var='surfenerg',val=surfenerg)
    ierr = 2
 endif

 if (youngmod <= 0.) then
    call error('init_porosity','youngmod <= 0',var='youngmod',val=youngmod)
    ierr = 3
 endif

 if (idrag /= 1) then
    call error('init_porosity','idrag = 1 should be used for porosity',var='idrag',val=real(idrag))
    ierr = 4
 endif

end subroutine init_porosity

!-----------------------------------------------------------------------
!+
!  Compute the initial filling factor
!+
!-----------------------------------------------------------------------
subroutine init_filfac(npart,xyzh,vxyzu)
 use options,           only:use_dustfrac
 use viscosity,         only:shearparam
 use part,              only:idust,igas,iamtype,iphase,massoftype,&
                             rhoh,dustfrac,dustprop,filfac,Omega_k
 use dust,              only:get_viscmol_nu!,grainsizecgs
 use eos,               only:gamma,get_spsound
 integer, intent(in)       :: npart
 real, intent(in)          :: xyzh(:,:)
 real, intent(inout)       :: vxyzu(:,:)

 integer                   :: i,iam
 real                      :: rho,rhogas,cs,cparam,coeff_gei,nu
 real                      :: sfrac,s1,s2,s3,filfacmax
! real                      :: mfrac,m1,m2,m3


 select case (iporosity)   ! add other case for other models here
 case (1)

    !--initialize filling factor (Garcia & Gonzalez 2020, Suyama et al. 2008, Okuzumi et al. 2012)

    if (all(filfac(:) == 0.)) then   ! check if filfac(i) was already initialize by init_filfac or not
       coeff_gei = sqrt(8./(pi*gamma))
       do i=1,npart
          iam = iamtype(iphase(i))
          if (iam == idust .or. (iam == igas .and. use_dustfrac)) then
             sfrac = (dustprop(1,i)/mmono)**(1./3.)
             if (sfrac > 1.) then      ! if grainsize > monomer size, compute filling factor
                !- compute rho, rhogas and cs
                if (iam == igas .and. use_dustfrac) then
                   rho = rhoh(xyzh(4,i),massoftype(igas))
                   rhogas = rho*(1-dustfrac(1,i))
                   cs = get_spsound(3,xyzh(:,i),rhogas,vxyzu(:,i))
                else
                   rhogas = rhoh(xyzh(4,i),massoftype(igas))
                   cs = get_spsound(3,xyzh(:,i),rhogas,vxyzu(:,i))
                   rho = rhogas + rhoh(xyzh(4,i),massoftype(idust))
                endif

                !- molecular viscosity
                nu = get_viscmol_nu(cs,rhogas)

                !- shared parameter for the following filling factors
                cparam = (243.*pi*roottwo/15625.)*(Ro*shearparam*smono**4*dustprop(2,i)*dustprop(2,i)*cs &
                             *Omega_k(i))/(rho*b_oku*eroll)

                !--transition masses m1/mmono and m2/mmono between hit&stick and Epstein/Stokes regimes with St < 1
                s1 = (cparam/(2.*(2.**0.075 - 1.)*coeff_gei))**((1.-cratio)/(1.+8.*cratio))
                s2 = (cparam*cs*smono/(9.*nu*(2.**0.2 - 1.)))**((1.-cratio)/(9.*cratio))

                !--we assume St < 1 here (grainsizecgs < 100-1000 cm)
                if (s1 < s2) then
                   if (sfrac < s1) then      ! filling factor: hit&stick regime
                      filfac(i) = sfrac**(3.*cratio/(1.-cratio))
                   else
                      !- transition masses m3/mmono between Epstein and Stokes regimes with St < 1
                      s3 = s1**((1.+8.*cratio)/(1.-cratio)) / s2**(9.*cratio/(1.-cratio))

                      if (sfrac < s3) then   ! filling factor: Epstein regime - St<1
                         filfac(i) = s1**((1.+8.*cratio)/(3.-3.*cratio))/sfrac**(1./3.)
                      else                   ! filling factor: Stokes regime - St<1
                         filfac(i) = s2**(3.*cratio/(1.-cratio))
                      endif
                   endif
                else
                   if (sfrac < s2) then      ! filling factor: hit&stick regime
                      filfac(i) = sfrac**(3.*cratio/(1.-cratio))
                   else                      ! filling factor: Stokes regime - St<1
                      filfac(i) = s2**(3.*cratio/(1.-cratio))
                   endif
                endif

                !- max value of filfac is maxpacking == max compaction
                filfacmax = 0.5*maxpacking *(1+ sqrt(1 + 4*(1.-maxpacking)/maxpacking/maxpacking*sfrac**(-3)))
                if (filfac(i) > filfacmax) filfac(i) = filfacmax

                !- Compute grain mass of the grain using grain size and filfac
                dustprop(1,i) = filfac(i) * dustprop(1,i)
             else
                filfac(i) = 1.
                dustprop(1,i) = mmono
             endif
          endif
       enddo
    endif
 case (-1)
    !--initialize filling factor for compact grains
    if (all(filfac(:) == 0.)) then   ! check if filfac(i) was already initialize by init_filfac or not
       do i=1,npart
          iam = iamtype(iphase(i))
          if (iam == idust .or. (iam == igas .and. use_dustfrac)) then
             sfrac = (dustprop(1,i)/mmono)**(1./3.)
             if (sfrac > 1.) then      ! if grainsize > monomer size, compute filling factor
                filfac(i) = 1.
             else
                filfac(i) = 1.
                dustprop(1,i) = mmono
             endif
          endif
       enddo
    endif
 end select

end subroutine init_filfac

!----------------------------------------------------------
!+
!  print information about porosity
!+
!----------------------------------------------------------
subroutine print_porosity_info(iprint)
 integer, intent(in) :: iprint

 if (iporosity == 1) then
    write(iprint,"(a)")  ' Using porosity    '
    if (icompact == 1) then
       write(iprint,"(a)")  ' Using compaction during fragmentation    '
    endif
    write(iprint,"(2(a,1pg10.3),a)")' Monomer size = ',smonocgs,' cm = ',smono,' (code units)'
    write(iprint,"(2(a,1pg10.3),a)")' Surface energy = ',surfenergSI,' J/m**2 = ',surfenerg,' (code units)'
    write(iprint,"(2(a,1pg10.3),a)")' Young modulus = ',youngmodSI,' Pa = ',youngmod,' (code units)'
 endif

end subroutine print_porosity_info

!-----------------------------------------------------------------------
!+
!  Compute the final filling factor
!+
!-----------------------------------------------------------------------
subroutine get_filfac(npart,xyzh,mprev,filfac,dustprop,dt)
 use dim,               only:use_dustgrowth
 use options,           only:use_dustfrac
 use part,              only:rhoh,idust,igas,iamtype,iphase,isdead_or_accreted,&
                             massoftype,dustfrac,dustgasprop,VrelVf,probastick
 integer, intent(in)       :: npart
 real,    intent(in)       :: dt
 real,    intent(inout)    :: filfac(:),dustprop(:,:)
 real,    intent(in)       :: xyzh(:,:),mprev(:)
 integer                   :: i,iam
 real                      :: filfacevol,filfacmin,filfacmax
 real                      :: rho,rhod

 select case (iporosity)   ! add other cases for other models here
 case (1)
    !$omp parallel do default(none) &
    !$omp shared(xyzh,npart,iphase,massoftype,use_dustfrac,dustfrac,icompact) &
    !$omp shared(mprev,filfac,dustprop,dustgasprop,VrelVf,probastick,mmono,dt,ibounce) &
    !$omp private(i,iam,rho,rhod,filfacevol,filfacmin,filfacmax)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          iam = iamtype(iphase(i))

          if (iam == idust .or. (iam == igas .and. use_dustfrac)) then
             if (dustprop(1,i) > mmono) then
                !- compute rho = rho_gas + rho_dust

                if (use_dustfrac .and. iam == igas) then
                   rho = rhoh(xyzh(4,i),massoftype(igas))
                   rhod = rho*dustfrac(1,i)
                else
                   rhod = rhoh(xyzh(4,i),massoftype(idust))
                   rho = dustgasprop(2,i) + rhod
                endif

                call get_filfac_min(i,rho,dustprop(1,i)/mmono,dustprop(2,i),dustgasprop(:,i),filfacmin)
                !--if new mass > previous mass, compute the new filling factor due to growth
                if (dustprop(1,i) > mprev(i)) then
                   call get_filfac_growth(mprev(i),dustprop(1,i),filfac(i),dustgasprop(:,i),filfacevol)
                   if (ibounce == 1) call get_filfac_bounce(mprev(i),dustprop(2,i),filfac(i),&
                                                    dustgasprop(:,i),probastick(i),rhod,dt,filfacevol,filfacmin)
                   !--if new mass < previous mass, compute the new filling factor due to fragmentation
                else
                   call get_filfac_frag(mprev(i),dustprop(:,i),filfac(i),dustgasprop(:,i),rhod,VrelVf(i),dt,filfacevol)
                endif
                filfac(i) = filfacevol

                !--check if the filling factor is smaller than the minimum filling factor
                filfac(i) = max(filfac(i),filfacmin)
                !-- max value of filfac is maxpacking == max compaction
                filfacmax = maxpacking + (1.-maxpacking)*mmono/dustprop(1,i)
                filfac(i) = min(filfac(i),filfacmax)
             else
                filfac(i) = 1.
                dustprop(1,i) = mmono
             endif
          endif
       else
          filfac(i) = 0.
       endif
    enddo
    !$omp end parallel do
 end select

end subroutine get_filfac

!-----------------------------------------------------------------------
!+
!  Compute the filling factor during growth
!+
!-----------------------------------------------------------------------
subroutine get_filfac_growth(mprev,mass,filfac,dustgasprop,filfacgrowth)
 use viscosity,         only:shearparam
 use growth,            only:vrelative
 real, intent(in)          :: mprev,mass,filfac
 real, intent(in)          :: dustgasprop(:)
 real, intent(out)         :: filfacgrowth
 real                      :: ekincdt,vrel,vt
 real                      :: j              ! Power of the filling factor dependency in mass

 vt = sqrt(roottwo*Ro*shearparam)*dustgasprop(1)
 vrel = vrelative(dustgasprop,vt)

 !- kinetic energy condition Ekin/(3*b_oku/eroll)
 ekincdt = mprev*vrel*vrel/(12.*b_oku*eroll)

 !-choose power according to the value of ekincdt
 if (ekincdt <= 1.) then
    j = cratio
 else
    j = -0.2
 endif

 !- filling factor due to growth
 filfacgrowth = filfac*(mass/mprev)**j

end subroutine get_filfac_growth

!-----------------------------------------------------------------------
!+
!  Compute the filling factor during bounce
!+
!-----------------------------------------------------------------------
subroutine get_filfac_bounce(mprev,graindens,filfac,dustgasprop,probastick,rhod,dt,filfacevol,filfacmin)
 use viscosity,         only:shearparam
 use growth,            only:vrelative,get_size
 use physcon,           only:fourpi
 real, intent(in)          :: mprev,graindens,filfac,probastick,rhod,dt
 real, intent(in)          :: dustgasprop(:),filfacmin
 real, intent(inout)       :: filfacevol
 real                      :: sdust,vrel,ncoll,vol,deltavol
 real                      :: ekin,pdyn,coeffrest,filfacbnc
 real                      :: vstick,vyield,vend,vt

 if (probastick < 1.) then
    vt = sqrt(roottwo*Ro*shearparam)*dustgasprop(1)
    vrel = vrelative(dustgasprop,vt)
    sdust = get_size(mprev,graindens,filfac)
    vstick = compute_vstick(mprev,sdust)                   !-compute vstick, i.e. max velocity before bouncing appears

    if (vrel >= vstick) then                               !-if vrel>=vstick -> bouncing
       vyield = compute_vyield(vstick)                    !-compute vyield, i.e. max velocity before inelastic collisions appear
       vend = compute_vend(vstick)                        !-compute vend, i.e. max velocity before there is only bouncing => no growth

       if (vrel < vyield) then                            !-elastic collision, no compaction
          filfacbnc = filfac
       else                                               !-inelastic collision, compaction
          vol = fourpi/3. * sdust**3
          ncoll = fourpi*sdust**2*rhod*vrel*dt/mprev     !-number of collision in dt
          ekin = mprev*vrel*vrel/4.
          coeffrest = get_coeffrest(vstick/vrel,vyield/vrel)                !-coefficient of restitution
          !pdyn = eroll * (filfac/(maxpacking - filfac)/smono)**3
          pdyn = eroll /((1./filfac - 1./maxpacking)*smono)**3
          deltavol = (1.-coeffrest*coeffrest)*ekin/pdyn
          if (deltavol > vol) deltavol = vol

          filfacbnc = filfac *(1./(1.-0.5*(deltavol/vol)))**ncoll
          if (filfacbnc > maxpacking) filfacbnc = maxpacking
       endif

       if (vrel < vend) then                              !-final filfac is a combination of filfac due to growth + bouncing
          if (filfacevol < filfacmin) filfacevol = filfacmin
          filfacevol = filfacevol*probastick + (1-probastick)*filfacbnc
       else
          filfacevol = filfacbnc
       endif
    endif
 endif

end subroutine get_filfac_bounce

!-----------------------------------------------------------------------
!+
!  Compute the filling factor during fragmentation
!+
!-----------------------------------------------------------------------
subroutine get_filfac_frag(mprev,dustprop,filfac,dustgasprop,rhod,VrelVf,dt,filfacfrag)
 use viscosity,         only:shearparam
 use growth,            only:vrelative,get_size
 use physcon,           only:fourpi
 real, intent(in)          :: mprev,filfac,rhod,VrelVf,dt
 real, intent(in)          :: dustprop(:),dustgasprop(:)
 real, intent(out)         :: filfacfrag
 real                      :: sdust,vrel,ncoll,vol,deltavol!,compfactor
 real                      :: ekin,pdyn,vt

 select case (icompact)
 case (1)
    ! model Garcia + Kataoka mod
    sdust = get_size(mprev,dustprop(2),filfac)
    vol = fourpi/3. * sdust**3
    vt = sqrt(roottwo*Ro*shearparam)*dustgasprop(1)
    vrel = vrelative(dustgasprop,vt)
    ncoll = fourpi*sdust**2*rhod*vrel*dt/mprev                                 !number of collisions in dt

    ekin = mprev*vrel*vrel/4. - (2.*mprev - dustprop(1))*0.85697283*eroll/mmono       !0.856973 = 3* 1.8 * 48/302.46
    pdyn = eroll /((1./filfac - 1./maxpacking)*smono)**3
    deltavol = ekin/pdyn                                                       !-ekin is kinetic energy - all energy needed to break monomers

    if (deltavol < 0) deltavol = 0.
    if (deltavol > vol) deltavol = vol

    filfacfrag = filfac *(1./(1.-0.5*exp(1-VrelVf**2.)*deltavol/vol))**ncoll
 case default ! (0)
    ! Fragmentation at constant filling factor
    filfacfrag = filfac
 end select

end subroutine get_filfac_frag

!-----------------------------------------------------------------------
!+
!  Compute the filling factor in the collisional compression regime
!+
!-----------------------------------------------------------------------
subroutine get_filfac_col(i,rho,mfrac,graindens,dustgasprop,filfaccol)
 use part,              only:Omega_k,dragreg
 use viscosity,         only:shearparam
 use dust,              only:get_viscmol_nu
 use eos,               only:gamma
 integer, intent(in)       :: i
 real,    intent(in)       :: rho,mfrac,graindens
 real,    intent(in)       :: dustgasprop(:)
 real,    intent(out)      :: filfaccol
 real                      :: cparam,coeff_gei,nu,kwok
 real                      :: m1,m2,m3,m4,m5

 !--compute filling factor due to collisions (Garcia & Gonzalez 2020, Suyama et al. 2008, Okuzumi et al. 2012)

 !- shared parameter for the following filling factors
 cparam = (243.*pi*roottwo/15625.)*(Ro*shearparam*smono**4*graindens*graindens*dustgasprop(1) &
          *Omega_k(i))/(rho*b_oku*eroll)

 coeff_gei = sqrt(8./(pi*gamma))

 !- molecular viscosity
 nu = get_viscmol_nu(dustgasprop(1),dustgasprop(2))

 !- Kwok (1975) correction for supersonic drag is important
 if (dragreg(i) == 2) then
    kwok = sqrt(1.+9.*pi/128.*dustgasprop(4)*dustgasprop(4)/(dustgasprop(1)*dustgasprop(1)))
 else
    kwok = 1.
 endif

 !--transition sizes m1/mmono and m2/mmono between hit&stick and Epstein/Stokes regimes with St < 1
 m1 = (cparam/(2.*(2.**0.075 - 1.)*coeff_gei*kwok))**(0.375/(cratio+0.125))
 m2 = (cparam*dustgasprop(1)*smono/(9.*nu*(2.**0.2 - 1.)))**(1./(3.*cratio))

 if (dustgasprop(3) <= 1) then      !- Stokes < 1
    if (m1 < m2) then
       if (mfrac < m1) then      !- filling factor: hit&stick regime
          filfaccol = mfrac**cratio
       else
          !- transition masses m3/mmono between Epstein and Stokes regimes with St < 1
          m3 = m1**(8.*cratio+1.) / m2**(8*cratio)
          if (mfrac < m3) then  !- filling factor: Epstein regime - St<1
             filfaccol = m1**(cratio+0.125)/mfrac**(0.125)
          else                  !- filling factor: Stokes regime - St<1
             filfaccol = m2**cratio
          endif
       endif
    else
       if (mfrac < m2) then      !- filling factor: hit&stick regime
          filfaccol = mfrac**cratio
       else                      !- filling factor: Stokes regime - St<1
          filfaccol = m2**cratio
       endif
    endif
 else                              !- Stokes > 1
    !--transition masses m4/mmono and m5/mmono between hit&stick and Epstein/Stokes regimes with St > 1
    m4 = (rho*coeff_gei*kwok*dustgasprop(1)/(graindens*smono*Omega_k(i)))**4 / m1**((cratio+0.125)/0.375)
    m5 = (9.*nu*rho/(2.*graindens*smono**2*Omega_k(i)))**1.5 / m2**(0.5*cratio)

    if (m4 < m5) then            !- filling factor: Epstein regime - St>1
       filfaccol = m1**(cratio+0.125) * m4**0.075 / mfrac**0.2
    else                         !- filling factor: Stokes regime - St>1
       filfaccol = m2**cratio * (m5/mfrac)**0.2
    endif
 endif

end subroutine get_filfac_col

!-----------------------------------------------------------------------
!+
!  Compute the minimum filling factor
!+
!-----------------------------------------------------------------------
subroutine get_filfac_min(i,rho,mfrac,graindens,dustgasprop,filfacmin)
 use part,              only:Omega_k
 integer, intent(in)       :: i
 real,    intent(in)       :: rho,mfrac,graindens
 real,    intent(in)       :: dustgasprop(:)
 real,    intent(out)      :: filfacmin
 real                      :: filfaccol,filfacgas,filfacgrav

 call get_filfac_col(i,rho,mfrac,graindens,dustgasprop,filfaccol)

 !--compute filling factor due to gas drag compression (Garcia & Gonzalez 2020, Kataoka et al. 2013a)
 filfacgas = ((mmono*smono*dustgasprop(4)*Omega_k(i))/(pi*eroll*dustgasprop(3)))**(3./7.) * mfrac**(1./7.)

 !--compute filling factor due to self-gravity (Garcia & Gonzalez 2020, Kataoka et al. 2013b)
 filfacgrav = (mmono*mmono/(pi*smono*eroll))**0.6 * mfrac**0.4

 !--return the maximum filling factor between filfaccol, filfacgas and filfacgrav
 filfacmin = max(filfaccol,filfacgas,filfacgrav)

end subroutine get_filfac_min


subroutine get_disruption(npart,xyzh,filfac,dustprop,dustgasprop)
 use options,           only:use_dustfrac
 use part,              only:idust,igas,iamtype,iphase,massoftype,isdead_or_accreted,rhoh
 use growth,            only:check_dustprop,get_size
 use random,            only:ran2
 integer, intent(in)    :: npart
 real, intent(in)       :: xyzh(:,:),dustgasprop(:,:)
 real, intent(inout)    :: dustprop(:,:),filfac(:)
 integer                :: i,iam,seed
 real                   :: stress,strength,filfacmin,rho
 real                   :: grainmasscurlog,grainmassmaxlog,randmass

 select case (idisrupt)
 case(1)
    !$omp parallel do default(none) &
    !$omp shared(xyzh,npart,massoftype,iphase,use_dustfrac) &
    !$omp shared(filfac,dustprop,dustgasprop,mmono,smono,grainmassminlog,surfenerg,gammaft) &
    !$omp private(grainmasscurlog,grainmassmaxlog,randmass,seed) &
    !$omp private(i,iam,rho,filfacmin,stress,strength)
    do i=1, npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          iam = iamtype(iphase(i))
          if (iam == idust .or. (iam == igas .and. use_dustfrac)) then

             stress = 25./36. * dustprop(2,i) * filfac(i) * gammaft**2 * dustgasprop(4,i)**2
             strength = 0.6*filfac(i)**(1.8)*surfenerg/smono
             seed = int(stress)

             if (stress >= strength) then   !-grain is rotationnaly disrupted
                !-compute rho to compute filfacmin
                if (use_dustfrac .and. iam == igas) then
                   rho = rhoh(xyzh(4,i),massoftype(igas))
                else
                   rho = dustgasprop(2,i) + rhoh(xyzh(4,i),massoftype(idust))
                endif

                !-compute current, current/2 and min mass in log10
                grainmasscurlog = log10(dustprop(1,i))
                grainmassmaxlog = log10(dustprop(1,i)/(2.))

                !--call random number between 2 float values to assign a random mass to dustprop(1)
                if (grainmassmaxlog > grainmassminlog) then
                   randmass = (grainmassmaxlog - grainmassminlog) * ran2(seed) + grainmassminlog
                else
                   if (grainmasscurlog > grainmassminlog) then
                      randmass = grainmassminlog
                   else
                      randmass = grainmasscurlog
                   endif
                endif

                dustprop(1,i) = 10.**randmass

                !-compute filfacmin and compare it to filfac(i)
                call get_filfac_min(i,rho,dustprop(1,i)/mmono,dustprop(2,i),dustgasprop(:,i),filfacmin)
                filfac(i) = max(filfac(i),filfacmin)
             endif
          endif
       endif
    enddo
    !$omp end parallel do
 end select

end subroutine get_disruption

!-----------------------------------------------------------------------
!+
!  Compute the probability of bounce and associated growth rate
!+
!-----------------------------------------------------------------------

subroutine get_probastick(npart,xyzh,dmdt,dustprop,dustgasprop,filfac)
 use options,           only:use_dustfrac
 use part,              only:idust,igas,iamtype,iphase,isdead_or_accreted,rhoh,probastick
 use viscosity,         only:shearparam
 use growth,            only:vrelative,get_size
 integer, intent(in)       :: npart
 real, intent(in)          :: filfac(:)
 real, intent(in)          :: xyzh(:,:),dustprop(:,:),dustgasprop(:,:)
 real, intent(inout)       :: dmdt(:)
 integer                   :: i,iam
 real                      :: vrel,vstick,vend,sdust,vt

 if (ibounce == 1) then
    !$omp parallel do default(none) &
    !$omp shared(xyzh,npart,iphase,use_dustfrac) &
    !$omp shared(filfac,dmdt,dustprop,dustgasprop,probastick,shearparam) &
    !$omp private(i,iam,vrel,vstick,vend,sdust,vt)
    do i=1, npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          iam = iamtype(iphase(i))
          if ((iam == idust .or. (iam == igas .and. use_dustfrac))) then
             if (filfac(i) >= 0.3 .and. dmdt(i) >= 0.) then
                vt = sqrt(roottwo*Ro*shearparam)*dustgasprop(1,i)
                vrel = vrelative(dustgasprop(:,i),vt)
                sdust = get_size(dustprop(1,i),dustprop(2,i),filfac(i))
                vstick = compute_vstick(dustprop(1,i),sdust)
                vend = compute_vend(vstick)

                !compute the probability of bounce depending on the velocity
                if (vrel >= vstick) then
                   if (vrel < vend) then
                      probastick(i) = (log(vrel)-log(vend))/(log(vstick)-log(vend))
                   else
                      probastick(i) = 0.      !full bounce -> no growth
                   endif
                else
                   probastick(i) = 1.
                endif
             else
                probastick(i) = 1.
             endif
             !compute new growth rate
             dmdt(i) = dmdt(i)*probastick(i)
          endif
       endif
    enddo
    !$omp end parallel do
 endif

end subroutine get_probastick

!-----------------------------------------------------------------------
!+
!  Write porosity options in the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_porosity(iunit)
 use infile_utils,      only:write_inopt
 integer, intent(in)       :: iunit

 write(iunit,"(/,a)") '# options controlling porosity (require idrag=1)'
 call write_inopt(iporosity,'iporosity','porosity (0=off,1=on) ',iunit)
 if (iporosity == 1 .or. iporosity == -1) then
    call write_inopt(icompact, 'icompact', 'Compaction during fragmentation (ifrag > 0) (0=off,1=on)', iunit)
    call write_inopt(ibounce, 'ibounce', 'Dust bouncing (0=off,1=on)', iunit)
    call write_inopt(idisrupt, 'idisrupt', 'Rotational disruption (0=off,1=on)', iunit)
    call write_inopt(smonocgs,'smonocgs','Monomer size in cm (smaller or equal to 1.e-4 cm)',iunit)
    call write_inopt(surfenergSI,'surfenergSI','Monomer surface energy in J/m**2',iunit)
    call write_inopt(youngmodSI,'youngmodSI','Monomer young modulus in Pa',iunit)
    call write_inopt(gammaft,'gammaft','Force to torque efficient of gas flow on dust',iunit)
 endif

end subroutine write_options_porosity

!-----------------------------------------------------------------------
!+
!  Read porosity options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_porosity(name,valstring,imatch,igotall,ierr)
 use options,                 only: use_porosity
 character(len=*), intent(in)    :: name,valstring
 logical, intent(out)            :: imatch,igotall
 integer, intent(out)            :: ierr

 integer, save                   :: ngot = 0

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('iporosity')
    read(valstring,*,iostat=ierr) iporosity
    ngot = ngot + 1
    if (iporosity == 1 .or. iporosity == -1) use_porosity = .true.
 case('icompact')
    read(valstring,*,iostat=ierr) icompact
    ngot = ngot + 1
 case('ibounce')
    read(valstring,*,iostat=ierr) ibounce
    ngot = ngot + 1
 case('idisrupt')
    read(valstring,*,iostat=ierr) idisrupt
    ngot = ngot + 1
 case('smonocgs')
    read(valstring,*,iostat=ierr) smonocgs
    ngot = ngot + 1
 case('surfenergSI')
    read(valstring,*,iostat=ierr) surfenergSI
    ngot = ngot + 1
 case('youngmodSI')
    read(valstring,*,iostat=ierr) youngmodSI
    ngot = ngot + 1
 case('gammaft')
    read(valstring,*,iostat=ierr) gammaft
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 if ((iporosity == 0) .and. ngot == 1) igotall = .true.
 if ((iporosity /= 0) .and. ngot == 8) igotall = .true.

end subroutine read_options_porosity

!-----------------------------------------------------------------------
!+
!  Write porosity options to the .setup file
!+
!-----------------------------------------------------------------------
subroutine write_porosity_setup_options(iunit)
 use infile_utils,    only:write_inopt
 integer, intent(in)     :: iunit

 write(iunit,"(/,a)") '# options for porosity'
 call write_inopt(iporosity,'iporosity','porosity (0=Off,1=On)',iunit)
 call write_inopt(ibounce,'ibounce','bouncing (0=Off,1=On)',iunit)
 call write_inopt(idisrupt,'idisrupt','disruption (0=Off,1=On)',iunit)

end subroutine write_porosity_setup_options

!-----------------------------------------------------------------------
!+
!  Read growth options from the .setup file
!+
!-----------------------------------------------------------------------
subroutine read_porosity_setup_options(db,nerr)
 use options,         only:use_porosity
 use infile_utils,    only:read_inopt,inopts
 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(inout)                   :: nerr

 call read_inopt(iporosity,'iporosity',db,min=-1,max=1,errcount=nerr)
 if (iporosity == 1 .or. iporosity == -1) use_porosity = .true.
 call read_inopt(ibounce,'ibounce',db,min=0,max=1,errcount=nerr)
 call read_inopt(idisrupt,'idisrupt',db,min=0,max=1,errcount=nerr)

end subroutine read_porosity_setup_options

real function get_coeffrest(vstickvrel,vyieldvrel)
 real, intent(in)   :: vstickvrel,vyieldvrel

 if (vyieldvrel >= 1.) then
    get_coeffrest = sqrt(1.-vstickvrel*vstickvrel)
 else
    get_coeffrest =  sqrt(1.2*sqrt(3.)*(1.-(vyieldvrel*vyieldvrel/6.))*&
                      sqrt(1./(1.+2.*sqrt((1.2/(vyieldvrel*vyieldvrel))-0.2)))-(vstickvrel*vstickvrel))
 endif

end function get_coeffrest

!--velocity limit between full sticking regime and partial sticking + bouncing regime
real function compute_vstick(mass,size)
 real, intent(in) ::mass,size
 compute_vstick = 8.76*((surfenerg**5 * size**4)/(mass**3*youngmod**2))**(1./6.)
end function compute_vstick

!--velocity limit between elastic and inelastic bouncing regime
real function compute_vyield(vstick)
 real, intent(in) ::vstick
 compute_vyield = 10.*vstick
end function compute_vyield

!--velocity limit between partial sticking + bouncing regime and full bouncing regime
real function compute_vend(vstick)
 real, intent(in) ::vstick
 compute_vend = 24343220.*vstick
end function compute_vend

end module porosity
