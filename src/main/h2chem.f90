!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: chem
!
!  DESCRIPTION:
!   Contains routines for Hydrogen Chemistry
!   Routines are originally by Clare Dobbs
!   Modified extensively by Alex Pettitt
!   Translated to Fortran 90 and adapted
!   for use in Phantom by Daniel Price (2011)
!
!  REFERENCES:
!   Nelson & Langer (1997)
!   Glover et al. (2010)
!   Bergin et al. (2004)
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: h2cooling, part, physcon, units
!+
!--------------------------------------------------------------------------

module chem
 implicit none
 public :: init_chem,energ_h2cooling,update_abundances
!
!--some variables needed for CO chemistry, Nelson+Langer97
!
 real, parameter :: k0=5.0d-16
 real, parameter :: k1=5.0d-10

 ! parameter for CO dissociation, Nelson+Langer97 has 1.0d-10
 real, parameter :: phrates=1.56d-10
!
!--some variables needed for H2 chemistry
!
 !  Rate of formation on grains, 6E-18*sqrt(T) in Bergin04
 real, parameter :: Rconst=3.d-1 * 2.2d-18

 ! Rate of dissociation, Bergin04
 real, parameter :: rate_diss0=4.17d-11

 real(kind=8) :: udens,uergg
 private

contains

!----------------------------------------------------------------
!+
!  Routine to initialise a few variables at code startup
!+
!----------------------------------------------------------------
subroutine init_chem()
 use units, only:utime,udist,umass
!
!--some definitions at code startup
!
 udens = umass/udist**3       !Code units of density
 uergg = udist**2/utime**2    !Code units of energy (specific)

end subroutine init_chem

!----------------------------------------------------------------
!+
!  Compute cooling term in energy equation (du/dt)
!  Also updates abundances of various species
!+
!----------------------------------------------------------------
subroutine energ_h2cooling(ui,dudti,rhoi,chemarrays,nchem,dt,xi,yi,zi,divv,idudtcool,ichem)
 use h2cooling, only:nabn,nrates,cool_func,dlq
 use units,     only:utime,udist
 real,    intent(in)        :: ui,rhoi,dt,xi,yi,zi
 real(kind=4), intent(in)   :: divv
 real,    intent(inout) :: dudti
 integer, intent(in)    :: nchem,idudtcool,ichem
 real,    intent(inout) :: chemarrays(nchem)
 real :: abund(nabn)
 real :: ratesq(nrates)
 real                    :: tempiso,ylamq,np1,divv_cgs

!--ARP: changed below to add flags
!--Takes idudtcool flag, determining whether cooling is required
!--By default the routine will call abundance update

 if (dt /= 0.d0) then
!
! Get updated abundances of all species, updates 'chemarrays',
! 'abund' is for internal use only.
!
!--Always call this routine, only returns tempsio and np1 if ichem=0
    call update_abundances(ui,rhoi,chemarrays,nchem,abund,tempiso,np1,dt,xi,yi,zi,ichem)
!
! Call cooling routine, requiring total density, some distance measure and
! abundances in the 'abund' format
!
    if (idudtcool==1) then
!   Dont call cooling if only updating abundances
!
       divv_cgs = divv / utime
       call cool_func(tempiso,np1,dlq,divv_cgs,abund,ylamq,ratesq)
!
!   compute change in u from 'ylamq' above.
!
       dudti = dudti + (-1.d0*ylamq/(rhoi*udens))*utime**3/udist**2
    endif
 endif

end subroutine energ_h2cooling

!----------------------------------------------------------------
!+
!  Internal routine to update abundances of various species
!+
!----------------------------------------------------------------
subroutine update_abundances(ui,rhoi,chemarrays,nchem,abund,tempiso,np1,dt,xi,yi,zi,ichem)
 use part,    only:ih2ratio,iHI,iproton,ielectron,iCO
 use units,   only:utime,udist
 use physcon, only:mp=>mass_proton_cgs,Rg
 use h2cooling, only:nabn,nrates,dchem,cosmic_ray_ion_rate,&
                   abundc,abundo,abundsi,abunde,AV_conversion_factor,hchem,dphotflag,dphot0
 real,    intent(in)    :: ui,rhoi,dt,xi,yi,zi
 integer, intent(in)    :: nchem,ichem
 real,    intent(inout) :: chemarrays(nchem)
 real,    intent(out)   :: abund(nabn)
 real,    intent(out)   :: tempiso,np1
 real :: dtclare,gmwvar
 real :: h2ratio,abHIq,abeq,abco,abhpq
 real :: nh1,th2,h2mol,nh1t,nmol
 real :: cdens,rate,exprate
 real :: dttest,abcp,phrates2,gamma_chx,gamma_co
 real :: tstep,tstep2,beta,coeff1,coeff2,sqrttempiso,tsteptest1 !,tsteptest
 real :: np1tempRconst,k0_np1sq,third,k1ab,gamh2np1,dnp1,gammaco_np1
 real :: tstep10,totH2rate
 real :: rhodisk,ad,bd,MdMo,r2,dphot,bit1
 integer :: i,j,nstep,nstep2


!--If we want all chemistry:
 if (ichem==1) then



!--Is dphot set or varying by radial distance?
    select case(dphotflag)
    case(1)
       !--Vary from centre in accordance with MN stellar disk.
       MdMo   = 8.56*1.0e10  !Mass of MN disk
       ad     = 53.    !5.3kpc  in code units
       bd     = 2.5    !0.25kpc in code units
       r2     = xi*xi+yi*yi
       bit1   =sqrt(zi*zi+bd*bd)
       rhodisk= (bd*bd*MdMo/(4.*3.14159))*(ad*r2+(ad+3.*bit1)*(ad+bit1)**2) &
    / ( (r2+(ad+bit1)**2)**(5./2.)*bit1**(3./2.) )
       dphot  = udist*(0.0001*rhodisk)**(-1./3.)
    case default
       !--Use constant value from cooling module.
       dphot  =dphot0
    end select

!---------------------------------------------------------------------
! Setup chemistry, read in ab., calulate temp, densities and constants
!---------------------------------------------------------------------
    h2ratio = chemarrays(ih2ratio)
    abHIq   = chemarrays(iHI)
    abhpq   = chemarrays(iproton)
    abeq    = chemarrays(ielectron)
    abco    = chemarrays(iCO)

!--Convert physical timestep into seconds
    dtclare=dt*utime
    third = 1.d0/3.d0


!--Calc mean molecular calculated taking into account the molecular gas.
!--This is 1.27 when there is nomolecular hydrogen:
    gmwvar=(2.0d0*h2ratio+(1.d0-2.d0*h2ratio)+0.4d0)/ &
        (0.1d0+h2ratio+(1.d0-2.d0*h2ratio))

!--Calculate an isothermal temperature use gas law
    tempiso = 2.d0/3.d0*ui/(Rg/gmwvar/uergg)
    sqrttempiso = sqrt(tempiso)

!
! This first section is to update the H2 fraction. But this requires using a very small timestep and
! subcycling. First an estimate of the timestep to update the H2 is determined.
!
! np1 =total number density inclusive of protons, CONSTANT
! nh1 =number density of HI inclusive of protons
! nh21=number density of H2

    np1=(rhoi*udens/mp)*5.d0/7.d0     ! n = (5/7)*(rho/mp), gamma=7/5?
    dnp1  = 1.d0/np1                !Inverse for calculations


! Total column density, CONSTANT
!--N_tot = n_tot *ds
    cdens=np1*dphot

!--Calc f_dust absorption rate for H2 photo-diss, and grain constant Bergin04
!  f_dust = exp [-AVconvfac * N_tot]
    rate=-3.74d0*cdens*AV_conversion_factor   !Should that 1.0 be 3.74?
    exprate = exp(rate)
    np1tempRconst = np1*sqrttempiso*Rconst


!--Calc the rates needed for CO generation/loss, Nelson+Langer97/Glover10
!--Go=1.0 simplifies things:
!  Gamma_CO  =phrate*Go*exp[-2.5*AVconvfac*N_tot],   was -0.267*
!  Gamma_CHx =5*Gamma_CO
    phrates2    = phrates*exp(-2.5d0*AV_conversion_factor*cdens)
    gamma_chx   = 5.0d0*phrates2
    gamma_co    = phrates2
    gammaco_np1 = gamma_co*np1
    k0_np1sq    = k0*np1*np1

!---------------------------------------------------------------------
!H2 timsetpping set-up for formation/destruction
!---------------------------------------------------------------------
    th2=10000.d0    !Timestep for H2 initially
    nstep = 5000
    totH2rate = 0.d0
    !--From Bergin et al. 2004, d(n(h2))/dt=(formation-destruction)
    call H2fd_rate(np1,h2ratio,dphot,third,exprate,rate_diss0,np1tempRconst,cosmic_ray_ion_rate,totH2rate,nh1)
    !--Here we are looking to dt such that n(H2) goes to zero, i.e dt=0-n(H2)/(formation-destruction)
    !--First just check that dt is not 0. - but actually this should not occur anyway.
    !--if the fomation-destruction = 0 obviously can't 1/0 to find timestepping as above.
    if (totH2rate==0.d0) then
       !------Not really happening...
       nstep=5000
    else

       !--Calculate dt according to dt=0-n(H2)/(formation-destruction)
       dttest=-h2ratio*np1/totH2rate

! If dttest>0 then h2 is being photodissociated. Set timestep to correspond to 1/10th of the time to
! completely photodissociate to 0 (this step is basically ensuring that the amount of H2 does not
! become negative).
       if (dttest > 0.d0) th2=0.1d0*dttest/dtclare

! Set number of steps for subcycling. In the event that the h2 is 0 (i.e. for the first timestep) or h2
! being created, the number of timesteps is 201 (chosen by running the simulation initially and finding
! a sensible value).
       if (th2 > 0.d0) nstep=MAX(INT(1.d0/th2)+1,201)
    endif
!--Final timesteps:
    tstep=dtclare/nstep
    tstep10 = 10.d0*tstep
!---------------------------------------------------------------------
!Loop to actually calculate H2/CO formation/destruction
!---------------------------------------------------------------------
    do i=1,nstep

!--Determine number density of molecular hydrogen, must be at least 0.:
!  Calc:   n_H2[i+1] = dt*(dn_H2/dt) + n_H2[i]
       totH2rate = 0.d0
       call H2fd_rate(np1,h2ratio,dphot,third,exprate,rate_diss0,np1tempRconst,cosmic_ray_ion_rate,totH2rate,nh1)
       nmol=max( totH2rate*tstep + h2ratio*np1, 0.)

! There should also be a timestep criterion so that H2 abundances do not exceed 1.
! But in my calculations so far, the ratios don't reach 1 (there is no sharp cut off
! in H2 fraction) indicating that the minimum number of timesteps of 201 prevents this.
! But for higher density calculations, a criterion should be included.

!--Calc the number density of H2:
       h2mol=min(np1*0.5,nmol)
!--and ratio of H2 to HI+H2 (i.e. n(H2)/n(HI+H2)
       h2ratio = h2mol*dnp1 ! optimisation: h2mol/np1 -> h2mol*dnp1


!---------------------------------------------------------------------
!--CO timsetpping set-up for formation/destruction
!---------------------------------------------------------------------

!--Subcycle the CO formation inside H2 formation.
!--Assume all CII therefore: ab_CII(t) => ab_CII(0) - ab_CO(t)
       abcp = max(abundc - abco,0.)
!--Assume all OI therefore:  ab_OI(t)  => ab_OI(0)  - ab_CO(t)
       k1ab     = k1*(abundo - abco)
       gamh2np1 = gamma_chx/h2mol !gamma_chx/(h2ratio*np1)
       beta     = k1ab/(k1ab + gamh2np1)
       if (h2ratio <= 0.d0) then
          !As below, if no H2 then no CO, stops some div0's
          !If we have no H2, shouldn't have any CO.
          nstep2=1
          tstep2=tstep
       elseif (abco <= 0.d0) then
          nstep2=INT(rhoi*1000.d0)
          tstep2=tstep/nstep2
       else
!       tsteptest=-abco/(k0*abcp*beta*np1*np1 - gamma_co*abco*np1)
          tsteptest1=-(k0_np1sq*abcp*beta - gammaco_np1*abco)/abco
!           nstep2=max(INT(tstep*10./tsteptest),INT(trho(ipart)*10.),1)
          nstep2=max(INT(tstep10*tsteptest1),1)
          tstep2=tstep/nstep2       !Step at some whole fraction of H2 time step.
       endif

!---------------------------------------------------------------------
!--Sub-Loop to actually calculate CO formation/destruction
!---------------------------------------------------------------------

       do j=1,nstep2
          abcp= max(abundc - abco,0.)
          !If we have no H2, shouldn't have any CO.
          if (h2ratio <= 0.d0) then
             beta = 0.d0
             abco = 0.d0
          else
             k1ab = k1*(abundo-abco)
             beta = k1ab/(k1ab + gamh2np1)
             abco = max(abco+(k0_np1sq*abcp*beta - gammaco_np1*abco)*tstep2,0.)
          endif
       enddo


    enddo

!------------------------------------------------------------------------------------
! End of updating H2/CO ratio. Now to update HI/HII/e- ratio.
!------------------------------------------------------------------------------------
!--If were not including H2, could set h2ratio to a small value (e.g. 1.e-7) and just
!--have this part to calculate heating and cooloing (need nh1 and np1 though).
!
! column density of HI excluding protons
!
    nh1t=nh1*abHIq
!
! update abundances of HI, protons and electrons requires total number density
! and column number density of HI
!
    call hchem(tempiso,np1,nh1t*dchem,abeq,abhpq,coeff1,coeff2,sqrttempiso)
!
! New number density of HI
!
    nh1t=nh1t+(coeff1-coeff2*nh1t)*dtclare
!
! Calculate new abundance of HI
!
!--Need some built in case for when all the HI is in H2, stops a div0 occurring
!--from:  nh1=np1*(1.-2.*h2ratio) when h2ratio=0.5
    if (nh1t==0.0d0 .AND. nh1==0.0d0) then
       abHIq=1.0d0
    else
       abHIq=nh1t/nh1
    endif
!
! Make sure abHI not >1 - should be taken care of by Simon's chemistry.
!
!--Also make sure that the abundance dosen't drop below 0, and above 1
!--by fixing HI we also fix HII and e- problem.
    if (abHIq > 1.d0) abHIq = 1.d0
    if (abHIq < 0.d0) abHIq = 0.d0

!--Abundance of protons
    abhpq = 1.0d0 - abHIq
!--and of electrons
    abeq  = abhpq + abunde
!------------------------------------------------------------------------------------
! Set abundances of H2, HI, electrons and protons for cooling routine only.
! We assume that CI and SiI are not present, and that SiII does not vary from its
! initial value.
!------------------------------------------------------------------------------------
    abund(1)  = h2ratio                             ! H2
    abund(2)  = (1.d0-2.d0*h2ratio)*abHIq           ! HI
    abund(3)  = (1.d0-2.d0*h2ratio)*abhpq + abunde  ! e-
    abund(4)  = (1.d0-2.d0*h2ratio)*abhpq           ! p+ (HII)
    abund(5)  = max(0., abundo - abco)              ! OI
    abund(6)  = 0.                                  ! CI
    abund(7)  = max(0., abundc - abco)              ! CII
    abund(8)  = 0.                                  ! SiI
    abund(9)  = abundsi                             ! SiII
    abund(10) = abco                                ! CO
!------------------------------------------------------------------------------------
! Copy updated abundances back to global array, what is needed for main output.
! (inc. H2, HI, HII, e-, CO)
!------------------------------------------------------------------------------------
    chemarrays(ih2ratio)  = h2ratio                ! H2
    chemarrays(iHI)       = abHIq                  ! HI
    chemarrays(iproton)   = abhpq                  ! p+ (HII)
    chemarrays(ielectron) = abeq                   ! e-
    chemarrays(iCO)       = abco                   ! CO


    !--ARP added:
    !--If we only want to initialise parameters for the cooling:
 elseif(ichem==0) then

    !--read in chemistry of particle
    h2ratio = chemarrays(ih2ratio)
    abHIq   = chemarrays(iHI)
    abhpq   = chemarrays(iproton)
    abeq    = chemarrays(ielectron)
    abco    = chemarrays(iCO)

    !--assign to cool_func format array
    abund(1)  = h2ratio                             ! H2
    abund(2)  = (1.d0-2.d0*h2ratio)*abHIq           ! HI
    abund(3)  = (1.d0-2.d0*h2ratio)*abhpq + abunde  ! e-
    abund(4)  = (1.d0-2.d0*h2ratio)*abhpq           ! p+ (HII)
    abund(5)  = max(0., abundo - abco)              ! OI
    abund(6)  = 0.                                  ! CI
    abund(7)  = max(0., abundc - abco)              ! CII
    abund(8)  = 0.                                  ! SiI
    abund(9)  = abundsi                             ! SiII
    abund(10) = abco                                ! CO

    !--Calc cgs estimates of density and temperature needed by cool_func
    gmwvar= (2.0d0*h2ratio+(1.d0-2.d0*h2ratio)+0.4d0)/ &
        (0.1d0+h2ratio+(1.d0-2.d0*h2ratio))
    tempiso = 2.d0/3.d0*ui/(Rg/gmwvar/uergg)
    np1   = (rhoi*udens/mp)*5.d0/7.d0     ! n = (5/7)*(rho/mp), gamma=7/5?

 else
    !--nothing is done
    print*,'WARNING: Why call the h2chem routine if you don''t want cooling or chemistry?'
 endif

end subroutine update_abundances

!----------------------------------------------------------------
!+
!  Calc the rate for destruction/formation rate of H2 used above
!   --Bergin et al. 2004 formation/destruction equation.
!   --Uses Draine & Bertoldi 1996 photodissociation rate.
!+
!----------------------------------------------------------------
pure subroutine H2fd_rate(np1,h2ratio,dphot,third,exprate,rate_diss0,grainform,cr,totrate,nh1)
 real              :: nh21,cdensH2,cdens2,sqrtcdens2,nH2,rate_diss
 real, intent(in)  :: np1,h2ratio,third,exprate,rate_diss0,grainform,cr
 real, intent(in)  :: dphot
 real, intent(out) :: totrate,nh1

!--Repeat calc at setup for densities of n_HI and n_H2 using new ratios
! nh1 =number density of HI inclusive of protons
! nh21=number density of H2
 nh1 =np1*(1.d0-2.d0*h2ratio)
 nh21=np1*h2ratio

! Column density of H2
!--N_H2 = n_H2 *ds
 cdensH2=nh21*dphot

!--Calc f_shield for H2 rates from Draine96, b=3km/s:
!  x=N_H2/5.4E14
!  f_shield = [0.965/(1+x/b) + 0.035/sqrt(1+x)]*exp[-8.5E-4*sqrt(1+x)]
 cdens2 = cdensH2*2.0d-15      !this is the same as cdensH2/5.e14 - but multiply instead of divide
 sqrtcdens2 = sqrt(1.d0+cdens2)
 nH2=0.965d0/(1.d0+cdens2*third)**2 + 0.035d0*exp(-8.5d-4*sqrtcdens2)/sqrtcdens2

!--Calc total dissossiation rate fot H2: eta_diss = f_shield*f_dust*eta_diss(0)
 rate_diss=nH2*exprate*rate_diss0

!--Add the dissociation rate with cosmic ray destruction rate, and calc
!--difference between this and grain formation rate
 totrate = nh1*grainform  -  (rate_diss+cr)*nh21

end subroutine H2fd_rate

end module chem
