c This version is from http://www.ioffe.ru/astro/EIP/eos22.f on 2026 nov 21
** Equation of state for fully ionized electron-ion plasmas (EOS EIP)
* A.Y.Potekhin & G.Chabrier, Astron. Astrophys., 550, A43 (2013),
*       and references therein,
*     + updates following Baiko & Chugunov, MNRAS, 510, 2628 (2022).
* Please communicate comments/suggestions to Alexander Potekhin:
*                                            palex@astro.ioffe.ru
* Previously distributed versions (obsolete):
*   eos2000, eos2002, eos2004, eos2006, eos2007, eos2009, eos10, eos11,
*     eos13, eos14, and eos17.
* Last update: 09.06.22. All updates since 2008 are listed below.
**   L I S T   O F   S U B R O U T I N E S :
*  MAIN (normally commented-out) - example driving routine.
*  MELANGE9 - for arbitrary ionic mixture, renders total (ion+electron)
*          pressure, internal energy, entropy, heat capacity (all
*          normalized to the ionic ideal-gas values), logarithmic
*          derivatives of pressure over temperature and density.
*  EOSFI22 - nonideal (ion-ion + ion-electron + electron-electron)
*          contributions to the free and internal energies, pressure,
*          entropy, heat capacity, derivatives of pressure over
*          logarithm of temperature and over logarithm of density (all
*          normalized to the ionic ideal-gas values) for one ionic
*          component in a mixture.
*  FITION9 - ion-ion interaction contributions to the free and internal
*          energies, pressure, entropy, heat capacity, derivatives of
*          pressure over logarithms of temperature and density.
*  LIQUBC - quantum corrections in liquid to thermodynamic functions.
*  FSCRliq8 - ion-electron (screening) contributions to the free and
*          internal energies, pressure, entropy, heat capacity,
*          derivatives of pressure over logarithms of temperature and
*          density in the liquid phase for one ionic component in a
*          mixture.
*  FSCRsol8 - ion-electron (screening) contributions to the free and
*          internal energies, pressure, entropy, heat capacity,
*          derivatives of pressure over logarithms of temperature and
*          density for monoionic solid.
*  FHARM12 - harmonic (including static-lattice and zero-point)
*          contributions to the free and internal energies, pressure,
*          entropy, heat capacity, derivatives of pressure over
*          logarithms of temperature and density for solid OCP.
*  HLfit12 - the same as FHARM12, but only for thermal contributions
*  ANHBC - anharmonic contributions to the free and internal energies,
*          pressure, entropy, heat capacity, derivatives of pressure
*          over logarithms of temperature and density for solid OCP.
*  CORMIX - correction to the linear mixing rule for the Coulomb
*          contributions to the thermodynamic functions in the liquid.
*  ELECT11 - for an ideal electron gas of arbitrary degeneracy and
*          relativity at given temperature and electron chemical
*          potential, renders number density (in atomic units), free
*          energy, pressure, internal energy, entropy, heat capacity
*          (normalized to the electron ideal-gas values), logarithmic
*          derivatives of pressure over temperature and density.
*          It uses ELECT11a, which realizes fitting expressions based
*          on BLIN9 (listed below), and ELECT11b, which employs the
*          Sommerfeld expansions for very strong electron degeneracy.
*  SOMMERF provides the Sommerfeld expansion up to order T^2 for the
*          relativistic Fermi-Dirac integrals of orders 1/2, 3/2, 5/2,
*          together with the corresponding approximations for their
*          first, second, and some third derivatives.
*  EXCOR7 - electron-electron (exchange-correlation) contributions to
*          the free and internal energies, pressure, entropy, heat
*          capacity, derivatives of pressure over logarithm of
*          temperature and over logarithm of density (all normalized
*          to the classical electron ideal-gas values).
*  FERINV7 - inverse non-relativistic Fermi integrals of orders -1/2,
*          1/2, 3/2, 5/2, and their first and second derivatives.
*  BLIN9 - relativistic Fermi-Dirac integrals of orders 1/2, 3/2, 5/2,
*          and their first, second, and some third derivatives.
*          It uses BLIN9a, BLIN9b, and BLIN9c - three segments of the
*          fit by Blinnikov, Dunina-Barkovskaya, & Nadyozhin (1996,
*          ApJS 106, 171; erratum: ibid., 118, 603 (1998)).
*  CHEMFIT7 - electron chemical potential at given density and
*          temperature, and its first derivatives over density and
*          temperature and the second derivative over temperature.
** N O T E S   O N   T H E   U N I T S :
* Some variables in the subroutines MELANGE9,  ELECT11(a,b), and CHEMFIT
*   are in atomic units (a.u.), that is, m_e = \hbar = e = k = 1 (here,
*   m_e is the electron mass, \hbar is the reduced Planck constant,
*   e is the elementary charg, and k is the Boltzmann constant).
* Some variables in the subroutines SOMMERF, SUBFERMJ, BLIN9(a,b,c), and
*   CHEMFIT7 are in relativistic units, that is, m_e = \hbar = c = k = 1
*   (c is the speed of light in vacuum).
* The conversions between the different systems of units is based on the
*   CODATA Internationally recommended 2018 values of the fundamental
*   constants, URL https://physics.nist.gov/cuu/Constants/index.html :
*    - m_e = 9.1093837e-31 kg = 9.1093837e-28 g,
*    - \hbar = h/2\pi = 1.05457181765e-34 J s = 1.05457181765e-27 erg s,
*      where the Planck constant h = 6.62607015e-34 J s (exact),
*    - e = 1.602176634e-19 C (exact) = 4.80320471257e-10 statC,
*    - k = 1.380649e-23 J/K = 1.380649e-16 erg/K (exact),
*    - c = 2.99792458e8 m/s = 2.99792458e8 cm/s (exact).
** I M P R O V E M E N T S   S I N C E   2 0 0 8 :
* 2009:
*  - FHARM8 uses a fit HLfit8 to the thermal free energy of the harmonic
*    Coulomb lattice, which is more accurate than its predecessor FHARM7.
*    Resulting corrections amount up to 20% for the ion heat capacity.
*    Accordingly, S/R D3fit and FthCHA7 deleted (not used anymore).
*  - BLIN7 upgraded to BLIN8:
*       - cleaned (a never-reached if-else branch deleted);
*       - Sommerfeld (high-\chi) expansion improved;
*       - some third derivatives added.
*  - CORMIX added (and MELANGE7 upgraded to MELANGE8 accordingly).
*  - ANHARM7 upgraded to ANHARM8, more consistent with Wigner-Kirkwood.
*  - Since the T- and rho-dependences of Z values in a mixture
*    are not considered, the corresponding inputs (AYLR, AYLT) are
*    excluded from MELANGE8 (and EOSFI7 changed to EOSFI8 accordingly).
*  - ELECT7 upgraded to ELECT9 (high-degeneracy behaviour is improved)
* 2010: - ELECT9 upgraded (smooth match of two fits at chi >> 1);
*  - BLIN8 replaced by BLIN9 - smooth fit interfaces at chi=0.6 and 14;
*  - MELANGE8 replaced by MELANGE9 - slightly modified input/output.
* 08.08.11 - corrected mistake (unlikely to have an effect) in CHEMFIT7
* 16.11.11 - ELECT9 upgraded to ELECT11 (additional output)
* 20.04.12 - FHARM8 and HLfit8 upgraded to FHARM12 and HLfit12:
*   output of HLfit12 does not include zero-point vibr., but provides U1
* 22.12.12 - MELANGE9 now includes a correction to the linear mixing
*   rule (LMR) for the Madelung energy in the random bcc multi-ion
*   lattice.
* 14.05.13 - an accidental error in programming the newly introduced
*   correction to the LMR is fixed.
* 20.05.13 - calculation of the Wigner-Kirkwood quantum diffraction term
*   for the liquid plasma is moved from EOSFI8 into MELANGE9.
* 10.12.14 - slight cleaning of the text (no effect on the results)
* 28.05.15 - an accidental error in Wigner-Kirkwood entropy correction
*   is fixed (it was in the line "Stot=Stot+FWK*DENSI" since 20.05.13)
* 29.08.15 - eliminated underflow of exp(-THETA) in CHEMFIT7
* 10.08.16 - modified criteria to avoid accuracy loss (round-off errors)
* 07.02.17 - included possibility to switch off the WK (Wigner) terms
* 27.05.17 - safeguard against Zion < 1 is added in FSCRsol8;
*   safeguard against huge (-CHI) values is added in ELECT11.
* 27.01.19 - safeguard against X1=0 in CORMIX.
* 18.04.20 - corrected Wigner-Kirkwood term for heat capacity.
* 04.03.21 - corrected SUBFERMJ: defined parameter EPS (was undefined).
* 09.06.22 - implemented results of Baiko & Chugunov (2022):
*   - quantum corrections in liquid LIQUBC are included in the OCP,
*     thus upgraded EOSFI8 is named EOSFI22;
*   - anharmonic corrections in crystal (ANHARM8) are upgraded to ANHBC.
*   - Additionally, quantum melting is excluded, because thermodynamics
*     of quantum liquids at corresponding densities is poorly known.
************************************************************************
*                           MAIN program:               Version 02.06.09
* This driving routine allows one to compile and run this code "as is".
* In practice, however, one usually needs to link subroutines from this
* file to another (external) code, therefore the MAIN program is
* normally commented-out.
C%C      implicit double precision (A-H), double precision (O-Z)
C%C      parameter(MAXY=10,UN_T6=.315775,EPS=1.d-7)
C%C      dimension AY(MAXY),AZion(MAXY),ACMI(MAXY)
C%C      write(*,'('' Introduce the chemical composition (up to'',I3,
C%C     *  '' ion species):''/
C%C     *  '' charge number Z_i, atomic weight A_i,'',
C%C     *  '' partial number density x_i, derivatives d x_i / d ln T'',
C%C     *  '' and d x_i / d ln rho''/
C%C     /   '' (non-positive Z, A, or x=1 terminates the input)'')') MAXY
C%C      NMIX=0
C%C    3 continue
C%C      XSUM=0.
C%C      do IX=1,MAXY
C%C         write(*,'(''Z, A ('',I2,''): ''$)') IX
C%C         read*,AZion(IX),ACMI(IX)
C%C        if (AZion(IX).le.0..or.ACMI(IX).le.0.) goto 2
C%C         write(*,'(''x ('',I2,''): ''$)') IX
C%C         read*,AY(IX)
C%C         XSUM=XSUM+AY(IX)
C%C        if (AY(IX).le.0.) goto 2
C%C         NMIX=IX
C%C        if (dabs(XSUM-1.d0).lt.EPS) goto 2
C%C      enddo
C%C    2 continue
C%C      if (NMIX.eq.0) then
C%C      print*,'There must be at least one set of positive (x,Z,A).'
C%C        goto 3
C%C      endif
C%C      write(*,114)
C%C      do IX=1,NMIX
C%C         write(*,113) IX,AZion(IX),ACMI(IX),AY(IX)
C%C      enddo
C%C    9 continue
C%C      write(*,'('' Input T (K) (<0 to stop): ''$)')
C%C      read*,T
C%C      if (T.le.0.) stop
C%C   10 continue
C%C      write(*,'('' Input RHO [g/cc] (<0 to new T): ''$)')
C%C      read*,RHO
C%C      if (RHO.le.0.) goto 9
C%C      RHOlg=dlog10(RHO)
C%C      Tlg=dlog10(T)
C%C      T6=10.d0**(Tlg-6.d0)
C%C      RHO=10.d0**RHOlg
C%C      write(*,112)
C%C    1 continue
C%C      TEMP=T6/UN_T6 ! T [au]
C%C      call MELANGE9(NMIX,AY,AZion,ACMI,RHO,TEMP, ! input
C%C     *   PRADnkT, ! additional output - radiative pressure
C%C     *   DENS,Zmean,CMImean,Z2mean,GAMI,CHI,TPT,LIQSOL, ! output param.
C%C     *   PnkT,UNkT,SNk,CV,CHIR,CHIT) ! output dimensionless TD functions
C%C      Tnk=8.31447d13/CMImean*RHO*T6 ! n_i kT [erg/cc]
C%C      P=PnkT*Tnk/1.d12 ! P [Mbar]
C%C      TEGRAD=CHIT/(CHIT**2+CHIR*CV/PnkT) ! from Maxwell relat.
*   --------------------   OUTPUT   --------------------------------   *
* Here in the output we have:
* RHO - mass density in g/cc
* P - total pressure in Mbar (i.e. in 1.e12 dyn/cm^2)
* PnkT=P/nkT, where n is the number density of ions, T temperature
* CV - heat capacity at constant volume, divided by number of ions, /
* CHIT - logarithmic derivative of pressure \chi_T
* CHIR - logarithmic derivative of pressure \chi_\rho
* UNkT - internal energy divided by NkT, N being the number of ions
* SNk - entropy divided by number of ions, /k
* GAMI - ionic Coulomb coupling parameter
* TPT=T_p/T, where T_p is the ion plasma temperature
* CHI - electron chemical potential, divided by kT
* LIQSOL = 0 in the liquid state, = 1 in the solid state
C%C      write(*,111) RHO,T6,P,PnkT,CV,CHIT,CHIR,UNkT,SNk,GAMI,TPT,CHI,
C%C     *  LIQSOL
C%C      goto 10
C%C  112 format(/
C%C     *  ' rho [g/cc]     T6 [K]      P [Mbar]   P/(n_i kT)  Cv/(N k)',
C%C     *  '     chi_T       chi_r      U/(N k T)    S/(N k)    Gamma_i',
C%C     *  '      T_p/T    chi_e liq/sol')
C%C  111 format(1P,12E12.3,I2)
C%C  113 format(I3,2F8.3,1PE12.4)
C%C  114 format('       Z      CMI        x_j')
C%C      end

      subroutine MELANGE9(NMIX,AY,AZion,ACMI,RHO,TEMP,
     *   PRADnkT,DENS,Zmean,CMImean,Z2mean,GAMImean,CHI,TPT,LIQSOL,
     *   PnkT,UNkT,SNk,CV,CHIR,CHIT)
*                                                       Version 09.06.22
* Difference from v.18.04.20: excluded WK corrections and upgraded
*   EOSFI8 to EOSFI22 by including quantum corrections in liquid and new
*   anharmonic corrections in classical and quantum crystal;
*   excluded quantum melting because of theoretical uncertainties there.
* Limitations:
* (a) inapplicable in the regimes of
*      (1) bound-state formation,
*      (2) presence of particles other than electrons and nuclei;
* (b) for the case of a composition gradually depending on RHO or TEMP,
*  second-order functions (CV,CHIR,CHIT in output) should not be trusted
* Choice of the liquid or solid regime - criterion GAMI [because the
*     choice based on comparison of total (non-OCP) free energies can be
*     sometimes dangerous because of the fit uncertainties ("Local field
*     correction" in solid and quantum effects in liquid are unknown)].
* Input: NMIX - number of different elements;
*        AY - their partial number densities,
*        AZion and ACMI - their charge and mass numbers,
*        RHO - total mass density [g/cc]
*        TEMP - temperature [in a.u. = Hartree / k = 315775 K].
* NB: instead of RHO, a true input is CHI, defined below
*     Hence, disagreement between RHO and DENS is the fit error (<0.4%)
* Output:
*         AY - rescaled so that to sum up to 1 and resorted (by AZion)
*         AZion - resorted in ascending order
*         ACMI - resorted in agreement with AZion
*         DENS - electron number density [in a.u.=6.7483345e24 cm^{-3}]
*         Zmean=<Z>, CMImean=<A> - mean ion charge and mass numbers,
*         Z2mean=<Z^2> - mean-square ion charge number
*         GAMImean - effective ion-ion Coulomb coupling constant
*         CHI = mu_e/kT, where mu_e is the electron chem.potential
*         TPT - effective ionic quantum parameter (T_p/T)
*         LIQSOL=0/1 for liquid/solid
*         SNk - entropy per 1 ion
*         UNkT - internal energy per kT per ion
*         PnkT - pressure / n_i kT, where n_i is the ion number density
*         PRADnkT - radiative pressure / n_i kT
*         CV - heat capacity per ion, div. by Boltzmann const.
*         CHIR - inverse compressibility -(d ln P / d ln V)_T ("\chi_r")
*         CHIT = (d ln P / d ln T)_V ("\chi_T")
      implicit double precision (A-H), double precision (O-Z)
      character CHWK
      save
      parameter(TINY=1.d-7)
      dimension AY(*),AZion(*),ACMI(*)
      parameter (PI=3.141592653d0,C53=5.d0/3.d0,C13=1.d0/3.d0,
     *   AUM=1822.88848d0, ! a.m.u./m_e
     *   GAMIMELT=175., ! OCP value of Gamma_i for melting
     *   RAD=2.55684543477d-7) ! Radiation constant=4\sigma/c [a.u.]
      if (RHO.lt.1.e-19.or.RHO.gt.1.e15) stop'MELANGE: RHO out of range'
      Y=0.
      do IX=1,NMIX
         Y=Y+AY(IX)
      enddo
      if (dabs(Y-1.d0).gt.TINY) then
        do IX=1,NMIX
           AY(IX)=AY(IX)/Y
        enddo
         print*,'MELANGE9: partial densities (and derivatives)',
     *    ' are rescaled by factor',1./Y
      endif
* Sort the elements in ascending order in Z_j:
      KSORT=0
      do I=2,NMIX
         J=I
         Z=AZion(J)
         CMI=ACMI(J)
         Y=AY(J)
    1   if (J.le.1.or.AZion(J-1).le.Z) goto 2
           AZion(J)=AZion(J-1)
           ACMI(J)=ACMI(J-1)
           AY(J)=AY(J-1)
           J=J-1
           KSORT=1
        goto 1
    2    AZion(J)=Z
         ACMI(J)=CMI
         AY(J)=Y
      enddo
      if (KSORT.eq.1) write(*,'('' Ions are resorted as follows:''/
     *  '' i    Z_i       A_i       x_i''/(0P,I3,'':'',1P,3E10.3))')
     *  (J,AZion(J),ACMI(J),AY(J),J=1,NMIX)
* Calculation of average values:
      Zmean=0.
      Z2mean=0.
      Z52=0.
      Z53=0.
      Z73=0.
      Z321=0. ! corr.26.12.09
      CMImean=0.
      do IX=1,NMIX
         Zmean=Zmean+AY(IX)*AZion(IX)
         Z2mean=Z2mean+AY(IX)*AZion(IX)**2
         Z13=AZion(IX)**C13
         Z53=Z53+AY(IX)*Z13**5
         Z73=Z73+AY(IX)*Z13**7
         Z52=Z52+AY(IX)*dsqrt(AZion(IX))**5
         Z321=Z321+AY(IX)*AZion(IX)*dsqrt(AZion(IX)+1.d0)**3 ! 26.12.09
         CMImean=CMImean+AY(IX)*ACMI(IX)
      enddo
* (0) Photons:
      UINTRAD=RAD*TEMP**4
      PRESSRAD=UINTRAD/3.
C      CVRAD=4.*UINTRAD/TEMP
* (1) ideal electron gas (including relativity and degeneracy)  -----  *
      DENS=RHO/11.20587*Zmean/CMImean ! number density of electrons [au]
      call CHEMFIT(DENS,TEMP,CHI)
* NB: CHI can be used as true input instead of RHO or DENS
      call ELECT11(TEMP,CHI,
     *  DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,
     *  DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
* NB: at this point DENS is redefined (the difference can be ~0.1%)
      DTE=DENS*TEMP
      PRESSE=PEid*DTE ! P_e [a.u.]
      UINTE=UEid*DTE ! U_e / V [a.u.]
* (2) non-ideal Coulomb EIP  ----------------------------------------  *
      RS=(.75d0/PI/DENS)**C13 ! r_s - electron density parameter
      RSI=RS*CMImean*Z73*AUM ! R_S - ion density parameter
      GAME=1.d0/RS/TEMP ! electron Coulomb parameter Gamma_e
      GAMImean=Z53*GAME   ! effective Gamma_i - ion Coulomb parameter
      if (GAMImean.lt.GAMIMELT) then
         LIQSOL=0 ! liquid regime
      else
         LIQSOL=1 ! solid regime
      endif
* Calculate partial thermodynamic quantities and combine them together:
      UINT=UINTE
      PRESS=PRESSE
      CVtot=CVE*DENS
      Stot=SEid*DENS
      PDLT=PRESSE*CHITE ! d P_e[a.u.] / d ln T
      PDLR=PRESSE*CHIRE ! d P_e[a.u.] / d ln\rho
      DENSI=DENS/Zmean ! number density of all ions
      PRESSI=DENSI*TEMP ! ideal-ions total pressure (normalization)
      TPT2=0.
      CTP=4.d0*PI/AUM/TEMP**2 ! common coefficient for TPT2.10.12.14
* Add Coulomb+xc nonideal contributions, and ideal free energy:
      do IX=1,NMIX
        if (AY(IX).lt.TINY) goto 10 ! skip this species
         Zion=AZion(IX)
         CMI=ACMI(IX)
         GAMI=Zion**C53*GAME ! Gamma_i for given ion species
         DNI=DENSI*AY(IX) ! number density of ions of given type
         PRI=DNI*TEMP ! = ideal-ions partial pressure (normalization)
         call EOSFI22(LIQSOL,CMI,Zion,RS,GAMI,
     *     FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *     FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
* First-order TD functions:
         UINT=UINT+UC2*PRI ! internal energy density (e+i+Coul.)
         Stot=Stot+DNI*(SC2-dlog(AY(IX))) !entropy per unit volume[a.u.]
         PRESS=PRESS+PC2*PRI ! pressure (e+i+Coul.) [a.u.]
* Second-order functions (they take into account compositional changes):
         CVtot=CVtot+DNI*CV2 ! C_V (e+i+Coul.)/ V (optim.10.12.14)
         PDLT=PDLT+PRI*PDT2 ! d P / d ln T
         PDLR=PDLR+PRI*PDR2 ! d P / d ln\rho
         TPT2=TPT2+CTP*DNI/ACMI(IX)*AZion(IX)**2 ! opt.10.12.14
   10   continue
      enddo ! next IX
* Wigner-Kirkwood perturbative correction for liquid:
      TPT=dsqrt(TPT2) ! effective T_p/T - ion quantum parameter
* (in the case of a mixture, this estimate is crude)
* Corrections to the linear mixing rule:
      if (LIQSOL.eq.0) then ! liquid phase
         call CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321,
     *     FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
      else ! solid phase (only Madelung contribution) [22.12.12]
         FMIX=0.
        do I=1,NMIX
          do J=I+1,NMIX
             RZ=AZion(J)/AZion(I)
             X2=AY(J)/(AY(I)+AY(J))
             X1=dim(1.d0,X2)
            if (X1.lt.TINY) goto 11 ! 27.01.19
            if (X2.lt.TINY) goto 11
             X=X2/RZ+(1.d0-1.d0/RZ)*X2**RZ
             GAMI=AZion(I)**C53*GAME ! Gamma_i corrected 14.05.13
             DeltaG=.012*(1.d0-1.d0/RZ**2)*(X1+X2*RZ**C53)
             DeltaG=DeltaG*X/X2*dim(1.d0,X)/X1
             FMIX=FMIX+AY(I)*AY(J)*GAMI*DeltaG
   11       continue
          enddo
        enddo
         UMIX=FMIX
         PMIX=FMIX/3.d0
         CVMIX=0.
         PDTMIX=0.
         PDRMIX=FMIX/2.25d0
      endif
      UINT=UINT+UMIX*PRESSI
      Stot=Stot+DENSI*(UMIX-FMIX)
      PRESS=PRESS+PMIX*PRESSI
      CVtot=CVtot+DENSI*CVMIX
      PDLT=PDLT+PRESSI*PDTMIX
      PDLR=PDLR+PRESSI*PDRMIX
* First-order:
      PRADnkT=PRESSRAD/PRESSI ! radiative pressure / n_i k T
C      CVtot=CVtot+CVRAD
C      Stot=Stot+CVRAD/3.
      PnkT=PRESS/PRESSI ! P / n_i k T
      UNkT=UINT/PRESSI ! U / N_i k T
C      UNkT=UNkT+UINTRAD/PRESSI
      SNk=Stot/DENSI ! S / N_i k
* Second-order:
      CV=CVtot/DENSI ! C_V per ion
      CHIR=PDLR/PRESS ! d ln P / d ln\rho
      CHIT=PDLT/PRESS ! d ln P / d ln T
C      CHIT=CHIT+4.*PRESSRAD/PRESS ! d ln P / d ln T
      return
      end

      subroutine EOSFI22(LIQSOL,CMI,Zion,RS,GAMI,
     *  FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *  FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
*                                                       Version 09.06.22
* Non-ideal parts of thermodynamic functions in the fully ionized plasma
* Stems from EOSFI8 v.10.12.14.
* Diff.: quantum corrections in liquid LIQUBC and new anharmonic cterms
*   in solid ANHBC according to Baiko & Chugunov, MNRAS 510, 2628 (2022)
* Input: LIQSOL=0/1(liquid/solid),
*        Zion,CMI - ion charge and mass numbers,
*        RS=r_s (electronic density parameter),
*        GAMI=Gamma_i (ion coupling),
* Output: FC1 and UC1 - non-ideal "ii+ie+ee" contribution to the
*         free and internal energies (per ion per kT),
*         PC1 - analogous contribution to pressure divided by (n_i kT),
*         CV1 - "ii+ie+ee" heat capacity per ion [units of k]
*         PDT1=(1/n_i kT)*(d P_C/d ln T)_V
*         PDR1=(1/n_i kT)*(d P_C/d ln\rho)_T
* FC2,UC2,PC2,SC2,CV2 - analogous to FC1,UC1,PC1,SC1,CV1, but including
* the part corresponding to the ideal ion gas. This is useful for
* preventing accuracy loss in some cases (e.g., when SC2 << SC1).
* FC2 does not take into account the entropy of mixing S_{mix}: in a
* mixture, S_{mix}/(N_i k) has to be added externally (see MELANGE9).
* FC2 does not take into account the ion spin degeneracy either.
* When needed, the spin term must be added to the entropy externally.
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(C53=5.d0/3.d0,C76=7.d0/6.d0) ! TINY excl.10.12.14
      parameter (AUM=1822.88848d0) ! a.m.u/m_e
      if (LIQSOL.ne.1.and.LIQSOL.ne.0) stop'EOSFI22: invalid LIQSOL'
      if (CMI.le..1) stop'EOSFI22: too small CMI'
      if (Zion.le..1) stop'EOSFI22: too small Zion'
      if (RS.le..0) stop'EOSFI22: invalid RS'
      if (GAMI.le..0) stop'EOSFI22: invalid GAMI'
      GAME=GAMI/Zion**C53
      call EXCOR7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC) ! "ee"("xc")
* Calculate "ii" part:
      COTPT=dsqrt(3.d0/AUM/CMI)/Zion**C76 ! auxiliary coefficient
      TPT=GAMI/dsqrt(RS)*COTPT            ! = T_p/T in the OCP
      FidION=1.5*dlog(TPT**2/GAMI)-1.323515
* 1.3235=1+0.5*ln(6/pi); FidION = F_{id.ion gas}/(N_i kT), but without
* the term x_i ln x_i = -S_{mix}/(N_i k).
      if (LIQSOL.eq.0) then               ! liquid
         call FITION9(GAMI,
     *     FION,UION,PION,CVii,PDTii,PDRii)
         RSI=3.d0*(GAMI/TPT)**2 ! R_S - ion density parameter
         call LIQUBC(RSI,TPT,FQL,UQL,PQL,CVQL,PDTQL,PDRQL)
         FItot=FION+FidION+FQL
         UItot=UION+1.5d0+UQL
         PItot=PION+1.d0+PQL
         CVItot=CVii+1.5d0+CVQL
         SCItot=UItot-FItot
         PDTi=PDTii+1.d0+PDTQL
         PDRi=PDRii+1.d0+PDTQL
      else                                  ! solid
         call FHARM12(GAMI,TPT,
     *     Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm) ! harm."ii"
         call ANHBC(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah) ! anharm.
         FItot=Fharm+Fah
         FION=FItot-FidION
         UItot=Uharm+Uah
         UION=UItot-1.5d0 ! minus 1.5=ideal-gas, in order to get "ii"
         PItot=Pharm+Pah
         PION=PItot-1.d0 ! minus 1=ideal-gas
         PDTi=PDTharm+PDTah
         PDRi=PDRharm+PDRah
         PDTii=PDTi-1.d0 ! minus 1=ideal-gas
         PDRii=PDRi-1.d0 ! minus 1=ideal-gas
         CVItot=CVharm+CVah
         SCItot=Sharm+Uah-Fah
         CVii=CVItot-1.5d0 ! minus 1.5=ideal-gas
      endif
* Calculate "ie" part:
      if (LIQSOL.eq.1) then
         call FSCRsol8(RS,GAMI,Zion,TPT,
     *     FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR)
      else
         call FSCRliq8(RS,GAME,Zion,
     *     FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR)
         S_SCR=USCR-FSCR
      endif
* Total excess quantities ("ii"+"ie"+"ee", per ion):
      FC0=FSCR+Zion*FXC
      UC0=USCR+Zion*UXC
      PC0=PSCR+Zion*PXC
      SC0=S_SCR+Zion*SXC
      CV0=CVSCR+Zion*CVXC
      PDT0=PDTSCR+Zion*PDTXC
      PDR0=PDRSCR+Zion*PDRXC
      FC1=FION+FC0
      UC1=UION+UC0
      PC1=PION+PC0
      SC1=(UION-FION)+SC0
      CV1=CVii+CV0
      PDT1=PDTii+PDT0
      PDR1=PDRii+PDR0
* Total excess + ideal-ion quantities
      FC2=FItot+FC0
      UC2=UItot+UC0
      PC2=PItot+PC0
      SC2=SCItot+SC0
      CV2=CVItot+CV0
      PDT2=PDTi+PDT0
      PDR2=PDRi+PDR0
      return
      end

* ==================  ELECTRON-ION COULOMB LIQUID  =================== *
      subroutine FITION9(GAMI,FION,UION,PION,CVii,PDTii,PDRii)
*                                                       Version 11.09.08
* Dummy argument Zion is deleted in 2009.
* Non-ideal contributions to thermodynamic functions of classical OCP.
*   Stems from FITION00 v.24.05.00.
* Input: GAMI - ion coupling parameter
* Output: FION - ii free energy / N_i kT
*         UION - ii internal energy / N_i kT
*         PION - ii pressure / n_i kT
*         CVii - ii heat capacity / N_i k
*         PDTii = PION + d(PION)/d ln T = (1/N_i kT)*(d P_{ii}/d ln T)
*         PDRii = PION + d(PION)/d ln\rho
*   Parameters adjusted to Caillol (1999).
      implicit double precision (A-H),double precision (O-Z)
      save
      parameter (A1=-.907347d0,A2=.62849d0,C1=.004500d0,G1=170.0,
     *  C2=-8.4d-5,G2=.0037,SQ32=.8660254038d0) ! SQ32=sqrt(3)/2
      A3=-SQ32-A1/dsqrt(A2)
      F0=A1*(dsqrt(GAMI*(A2+GAMI))-
     -     A2*dlog(dsqrt(GAMI/A2)+dsqrt(1.+GAMI/A2)))+
     +     2.*A3*(dsqrt(GAMI)-datan(dsqrt(GAMI)))
      U0=dsqrt(GAMI)**3*(A1/dsqrt(A2+GAMI)+A3/(1.d0+GAMI))
*   This is the zeroth approximation. Correction:
      UION=U0+C1*GAMI**2/(G1+GAMI)+C2*GAMI**2/(G2+GAMI**2)
      FION=F0+C1*(GAMI-G1*dlog(1.d0+GAMI/G1))+
     +   C2/2.*dlog(1.d0+GAMI**2/G2)
      CVii=-0.5*dsqrt(GAMI)**3*(A1*A2/dsqrt(A2+GAMI)**3+
     +  A3*(1.d0-GAMI)/(1.d0+GAMI)**2) -
     -  GAMI**2*(C1*G1/(G1+GAMI)**2+C2*(G2-GAMI**2)/(G2+GAMI**2)**2)
      PION=UION/3.
      PDRii=(4.*UION-CVii)/9. ! p_{ii} + d p_{ii} / d ln\rho
      PDTii=CVii/3. ! p_{ii} + d p_{ii} / d ln T
      return
      end

      subroutine LIQUBC(RSI,TPT,FQL,UQL,PQL,CVQL,PDTQL,PDRQL)
*                                                       Version 04.06.22
* Quantum additions to classical OCP liquid.
* Input: RSI - ion density parameter = ion-sphere radius / ion Bohr rad.
*        TPT - ion quantum parameter = T_pi/T (T_pi - ion plasma temp.)
* Output - quantum-liquid contributions (QL):
*         FQL - QL free energy / N_i kT
*         UQL - QL internal energy / N_i kT
*         PQL - QL pressure / n_i kT
*         CVQL - QL heat capacity / N_i k
*         PDTQL = PQL + d(PQL)/d ln T = (1/N_i kT)*(d PQL/d ln T)
*         PDRQL = PQL + d(PQL)/d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(TINY=1.d-19,EPS3=2.d-3,BIG=33.d0)
      parameter(P1=0.351d0,P2=0.294d0,P3=90.d0)
      if (TPT.le.0.) stop'LIQBC: non-positive TPT'
      if (RSI.le.0.) stop'LIQBC: non-positive RSI'
      FQL=0.
      UQL=0.
      PQL=0.
      CVQL=0.
      PDTQL=0.
      PDRQL=0.
      if (TPT.lt.TINY) return
      CI1=P1*RSI/(P3+RSI)
      DI1=P3/(P3+RSI)
      do I=1,3
        if (I.eq.1) then
           CI=CI1
           DI=DI1
        elseif (I.eq.2) then
           CI=P2
           DI=0.
        else
           CI=dsqrt(1.d0-CI1**2-P2**2)
           DI=-(CI1/CI)**2*DI1
        endif
         TPTI=CI*TPT
        if (TPTI.lt.EPS3) then
           FQLI=TPTI**2/24.d0
           UQLI=2.d0*FQLI
           CVQLI=-UQLI
           PDTQLI=-FQLI
           PDRQLI=UQLI
        else
          if (TPTI.gt.BIG) then ! asymptote "exp(-TPTI)=0"
             FQLI=TPTI*.5d0-dlog(TPTI)
             UQLI=TPTI*.5d0-1.d0
             CVQLI=-1.d0
          else
             ET=dexp(-TPTI)
             ET1=1.d0-ET
             FQLI=dlog(ET1/(ET*TPTI))-TPTI/2.d0 ! BC'22(34)
             UQLI=TPTI/ET1-TPTI/2.d0-1.d0 ! BC'22(31)
             CVQLI=TPTI**2*ET/ET1**2-1.d0 ! BC'22(35)
          endif
           PDTQLI=.5d0*CVQLI-DI*CVQLI/3.d0 ! BC'22(38)
           PDRQLI=.75d0*UQLI-.25d0*CVQLI-.5d0*DI*(UQLI-CVQLI/3.d0)-
     -       RSI/(P3+RSI)*DI/9.d0*UQLI ! BC'22(39)
        endif
         FQL=FQL+FQLI
         UQL=UQL+UQLI
         PQL=PQL+(.5d0-DI/3.d0)*UQLI
         CVQL=CVQL+CVQLI
         PDTQL=PDTQL+PDTQLI
         PDRQL=PDRQL+PDRQLI
      enddo
      return
      end

      subroutine FSCRliq8(RS,GAME,Zion,
     *     FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR) ! fit to the el.-ion scr.
*                                                       Version 11.09.08
*                                                       cleaned 16.06.09
* Stems from FSCRliq7 v. 09.06.07. Included a check for RS=0.
*   INPUT: RS - density parameter, GAME - electron Coulomb parameter,
*          Zion - ion charge number,
*   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
*           USCR - internal energy per kT per 1 ion (screen.contrib.)
*           PSCR - pressure divided by (n_i kT) (screen.contrib.)
*           CVSCR - heat capacity per 1 ion (screen.contrib.)
*           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      parameter(XRS=.0140047,TINY=1.d-19)
      if (RS.lt.0.) stop'FSCRliq8: RS < 0'
      if (RS.lt.TINY) then
         FSCR=0.
         USCR=0.
         PSCR=0.
         CVSCR=0.
         PDTSCR=0.
         PDRSCR=0.
         return
      endif
      SQG=sqrt(GAME)
      SQR=sqrt(RS)
      SQZ1=dsqrt(1.+Zion)
      SQZ=dsqrt(Zion)
      CDH0=Zion/1.73205 ! 1.73205=sqrt(3.)
      CDH=CDH0*(SQZ1**3-SQZ**3-1.)
      SQG=sqrt(GAME)
      ZLN=dlog(Zion)
      Z13=exp(ZLN/3.) ! Zion**(1./3.)
      X=XRS/RS ! relativity parameter
      CTF=Zion**2*.2513*(Z13-1.+.2/sqrt(Z13))
* Thomas-Fermi constant; .2513=(18/175)(12/\pi)^{2/3}
      P01=1.11*exp(.475*ZLN)
      P03=0.2+0.078*ZLN**2
      PTX=1.16+.08*ZLN
      TX=GAME**PTX
      TXDG=PTX*TX/GAME
      TXDGG=(PTX-1.)*TXDG/GAME
      TY1=1./(1.d-3*Zion**2+2.*GAME)
      TY1DG=-2.*TY1**2
      TY1DGG=-4.*TY1*TY1DG
      TY2=1.+6.*RS**2
      TY2DX=-12.*RS**2/X
      TY2DXX=-3.*TY2DX/X
      TY=RS**3/TY2*(1.+TY1)
      TYX=3./X+TY2DX/TY2
      TYDX=-TY*TYX
      TYDG=RS**3*TY1DG/TY2
      P1=(Zion-1.)/9.
      COR1=1.+P1*TY
      COR1DX=P1*TYDX
      COR1DG=P1*TYDG
      COR1DXX=P1*(TY*(3./X**2+(TY2DX/TY2)**2-TY2DXX/TY2)-TYDX*TYX)
      COR1DGG=P1*RS**3*TY1DGG/TY2
      COR1DXG=-P1*TYDG*TYX
      U0=.78*sqrt(GAME/Zion)*RS**3
      U0DX=-3.*U0/X
      U0DG=.5*U0/GAME
      U0DXX=-4.*U0DX/X
      U0DGG=-.5*U0DG/GAME
      U0DXG=-3.*U0DG/X
      D0DG=Zion**3
      D0=GAME*D0DG+21.*RS**3
      D0DX=-63.*RS**3/X
      D0DXX=252.*RS**3/X**2
      COR0=1.+U0/D0
      COR0DX=(U0DX-U0*D0DX/D0)/D0
      COR0DG=(U0DG-U0*D0DG/D0)/D0
      COR0DXX=(U0DXX-(2.*U0DX*D0DX+U0*D0DXX)/D0+2.*(D0DX/D0)**2)/D0
      COR0DGG=(U0DGG-2.*U0DG*D0DG/D0+2.*U0*(D0DG/D0)**2)/D0
      COR0DXG=(U0DXG-(U0DX*D0DG+U0DG*D0DX)/D0+2.*U0*D0DX*D0DG/D0**2)/D0
* Relativism:
      RELE=dsqrt(1.d0+X**2)
      Q1=.18/dsqrt(dsqrt(Zion))
      Q2=.2+.37/dsqrt(Zion)
      H1U=1.+X**2/5.
      H1D=1.+Q1*X+Q2*X**2
      H1=H1U/H1D
      H1X=.4*X/H1U-(Q1+2.*Q2*X)/H1D
      H1DX=H1*H1X
      H1DXX=H1DX*H1X+
     +  H1*(.4/H1U-(.4*X/H1U)**2-2.*Q2/H1D+((Q1+2.*Q2*X)/H1D)**2)
      UP=CDH*SQG+P01*CTF*TX*COR0*H1
      UPDX=P01*CTF*TX*(COR0DX*H1+COR0*H1DX)
      UPDG=.5*CDH/SQG+P01*CTF*(TXDG*COR0+TX*COR0DG)*H1
      UPDXX=P01*CTF*TX*(COR0DXX*H1+2.*COR0DX*H1DX+COR0*H1DXX)
      UPDGG=-.25*CDH/(SQG*GAME)+
     +  P01*CTF*(TXDGG*COR0+2.*TXDG*COR0DG+TX*COR0DGG)*H1
      UPDXG=P01*CTF*(TXDG*(COR0DX*H1+COR0*H1DX)+
     +  TX*(COR0DXG*H1+COR0DG*H1DX))
      DN1=P03*SQG+P01/RS*TX*COR1
      DN1DX=P01*TX*(COR1/XRS+COR1DX/RS)
      DN1DG=.5*P03/SQG+P01/RS*(TXDG*COR1+TX*COR1DG)
      DN1DXX=P01*TX/XRS*(2.*COR1DX+X*COR1DXX)
      DN1DGG=-.25*P03/(GAME*SQG)+
     +  P01/RS*(TXDGG*COR1+2.*TXDG*COR1DG+TX*COR1DGG)
      DN1DXG=P01*(TXDG*(COR1/XRS+COR1DX/RS)+TX*(COR1DG/XRS+COR1DXG/RS))
      DN=1.+DN1/RELE
      DNDX=DN1DX/RELE-X*DN1/RELE**3
      DNDXX=(DN1DXX-((2.*X*DN1DX+DN1)-3.*X**2*DN1/RELE**2)/RELE**2)/RELE
      DNDG=DN1DG/RELE
      DNDGG=DN1DGG/RELE
      DNDXG=DN1DXG/RELE-X*DN1DG/RELE**3
      FSCR=-UP/DN*GAME
      FX=(UP*DNDX/DN-UPDX)/DN
      FXDG=((UPDG*DNDX+UPDX*DNDG+UP*DNDXG-2.*UP*DNDX*DNDG/DN)/DN-
     -  UPDXG)/DN
      FDX=FX*GAME
      FG=(UP*DNDG/DN-UPDG)/DN
      FDG=FG*GAME-UP/DN
      FDGDH=SQG*DNDG/DN**2 ! d FDG / d CDH
      FDXX=((UP*DNDXX+2.*(UPDX*DNDX-UP*DNDX**2/DN))/DN-UPDXX)/DN*GAME
      FDGG=2.*FG+GAME*((2.*DNDG*(UPDG-UP*DNDG/DN)+UP*DNDGG)/DN-UPDGG)/DN
      FDXG=FX+GAME*FXDG
      USCR=GAME*FDG
      CVSCR=-GAME**2*FDGG
      PSCR=(X*FDX+GAME*FDG)/3.
      PDTSCR=-GAME**2*(X*FXDG+FDGG)/3.
      PDRSCR=(12.*PSCR+X**2*FDXX+2.*X*GAME*FDXG+GAME**2*FDGG)/9.
      return
      end

* ==============   SUBROUTINES FOR THE SOLID STATE   ================= *
      subroutine FSCRsol8(RS,GAMI,ZNUCL,TPT,
     *     FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR)
*                                                       Version 28.05.08
*                    undefined zero variable Q1DXG is wiped out 21.06.10
*                                 accuracy-loss safeguard added 10.08.16
*                              safequard against Zion < 1 added 27.05.17
* Fit to the el.-ion screening in bcc or fcc Coulomb solid
* Stems from FSCRsol8 v.09.06.07. Included a check for RS=0.
*   INPUT: RS - el. density parameter, GAMI - ion coupling parameter,
*          ZNUCL - ion charge, TPT=T_p/T - ion quantum parameter
*   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
*           USCR - internal energy per kT per 1 ion (screen.contrib.)
*           PSCR - pressure divided by (n_i kT) (screen.contrib.)
*           S_SCR - screening entropy contribution / (N_i k)
*           CVSCR - heat capacity per 1 ion (screen.contrib.)
*           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      dimension AP(4) ! parameters of the fit
      parameter (C13=1.d0/3.d0,ENAT=2.7182818285d0,TINY=1.d-19)
      data AP/1.1866,.684,17.9,41.5/,PX/.205/ ! for bcc lattice
cc      data AP/1.1857,.663,17.1,40./,PX/.212/ ! for fcc lattice
      if (RS.lt.0.) stop'FSCRliq8: RS < 0'
      if (RS.lt.TINY) then
         FSCR=0.
         USCR=0.
         PSCR=0.
         S_SCR=0.
         CVSCR=0.
         PDTSCR=0.
         PDRSCR=0.
         return
      endif
      Zion=ZNUCL
      if (Zion.lt.1.d0) then ! 27.05.17
         print*,'FSCRsol8 WARNING: Z =',Zion,' < 1.'
         Zion=1.d0
      endif
      XSR=.0140047/RS ! relativity parameter
      Z13=Zion**C13
      P1=.00352*(1.-AP(1)/Zion**.267+.27/Zion)
      P2=1.d0+2.25/Z13*
     *(1.+AP(2)*Zion**5+.222*Zion**6)/(1.+.222*Zion**6)
      ZLN=dlog(Zion)
      Finf=sqrt(P2/XSR**2+1.)*Z13**2*P1 ! The TF limit
      FinfX=-P2/((P2+XSR**2)*XSR)
      FinfDX=Finf*FinfX
      FinfDXX=FinfDX*FinfX-FinfDX*(P2+3.*XSR**2)/((P2+XSR**2)*XSR)
      R1=AP(4)/(1.+ZLN)
      R2=.395*ZLN+.347/Zion/sqrt(Zion)
      R3=1.d0/(1.d0+ZLN*sqrt(ZLN)*.01+.097/Zion**2)
      Q1U=R1+AP(3)*XSR**2
      Q1D=1.d0+R2*XSR**2
      Q1=Q1U/Q1D
      Q1X=2.*XSR*(AP(3)/Q1U-R2/Q1D)
      Q1XDX=Q1X/XSR+4.*XSR**2*((R2/Q1D)**2-(AP(3)/Q1U)**2)
      Q1DX=Q1*Q1X
      Q1DXX=Q1DX*Q1X+Q1*Q1XDX
* New quantum factor, in order to suppress CVSCR at TPT >> 1
      if (TPT.lt.6./PX) then
         Y0=(PX*TPT)**2
         Y0DX=Y0/XSR
         Y0DG=2.*Y0/GAMI
         Y0DXX=0.
         Y0DGG=Y0DG/GAMI
         Y0DXG=Y0DG/XSR
         Y1=dexp(Y0)
         Y1DX=Y1*Y0DX
         Y1DG=Y1*Y0DG
         Y1DXX=Y1*(Y0DX**2+Y0DXX)
         Y1DGG=Y1*(Y0DG**2+Y0DGG)
         Y1DXG=Y1*(Y0DX*Y0DG+Y0DXG)
         SA=1.d0+Y1
         SUPA=dlog(SA)
         SUPADX=Y1DX/SA
         SUPADG=Y1DG/SA
         SUPADXX=(Y1DXX-Y1DX**2/SA)/SA
         SUPADGG=(Y1DGG-Y1DG**2/SA)/SA
         SUPADXG=(Y1DXG-Y1DX*Y1DG/SA)/SA
         EM2=ENAT-2.d0
         SB=ENAT-EM2/Y1
         SUPB=dlog(SB)
         EM2Y1=EM2/(Y1**2*SB)
         SUPBDX=EM2Y1*Y1DX
         SUPBDG=EM2Y1*Y1DG
         SUPBDXX=EM2Y1*(Y1DXX-2.d0*Y1DX**2/Y1-Y1DX*SUPBDX)
         SUPBDGG=EM2Y1*(Y1DGG-2.d0*Y1DG**2/Y1-Y1DG*SUPBDG)
         SUPBDXG=EM2Y1*(Y1DXG-2.d0*Y1DX*Y1DG/Y1-Y1DG*SUPBDX)
         SUP=dsqrt(SUPA/SUPB)
         SUPX=.5d0*(SUPADX/SUPA-SUPBDX/SUPB)
         SUPDX=SUP*SUPX
         SUPG=.5d0*(SUPADG/SUPA-SUPBDG/SUPB)
         SUPDG=SUP*SUPG
         SUPDXX=SUPDX*SUPX+
     +     SUP*.5d0*(SUPADXX/SUPA-(SUPADX/SUPA)**2-
     -             SUPBDXX/SUPB+(SUPBDX/SUPB)**2)
         SUPDGG=SUPDG*SUPG+
     +     SUP*.5d0*(SUPADGG/SUPA-(SUPADG/SUPA)**2-
     -             SUPBDGG/SUPB+(SUPBDG/SUPB)**2)
         SUPDXG=SUPDX*SUPG+
     +     SUP*.5d0*((SUPADXG-SUPADX*SUPADG/SUPA)/SUPA-
     -             (SUPBDXG-SUPBDX*SUPBDG/SUPB)/SUPB)
      else
         SUP=PX*TPT
         SUPDX=.5d0*PX*TPT/XSR
         SUPDG=PX*TPT/GAMI
         SUPDXX=-.5d0*SUPDX/XSR
         SUPDGG=0.
         SUPDXG=SUPDX/GAMI
      endif
      GR3=(GAMI/SUP)**R3
      GR3X=-R3*SUPDX/SUP
      GR3DX=GR3*GR3X
      GR3DXX=GR3DX*GR3X-R3*GR3*(SUPDXX/SUP-(SUPDX/SUP)**2)
      GR3G=R3*(1.d0/GAMI-SUPDG/SUP)
      GR3DG=GR3*GR3G
      GR3DGG=GR3DG*GR3G+GR3*R3*((SUPDG/SUP)**2-SUPDGG/SUP-1.d0/GAMI**2)
      GR3DXG=GR3DG*GR3X+GR3*R3*(SUPDX*SUPDG/SUP**2-SUPDXG/SUP)
      W=1.d0+Q1/GR3
      WDX=Q1DX/GR3-Q1*GR3DX/GR3**2
      WDG=-Q1*GR3DG/GR3**2
      WDXX=Q1DXX/GR3-
     -   (2.d0*Q1DX*GR3DX+Q1*(GR3DXX-2.d0*GR3DX**2/GR3))/GR3**2
      WDGG=Q1*(2.d0*GR3DG**2/GR3-GR3DGG)/GR3**2
      WDXG=-(Q1DX*GR3DG+Q1*(GR3DXG-2.d0*GR3DX*GR3DG/GR3))/GR3**2
      FSCR=-GAMI*Finf*W
      FDX=-GAMI*(FinfDX*W+Finf*WDX)
      FDXX=-GAMI*(FinfDXX*W+2.d0*FinfDX*WDX+Finf*WDXX)
      FDG=-Finf*W-GAMI*Finf*WDG
      FDGG=-2.d0*Finf*WDG-GAMI*Finf*WDGG
      if (dabs(FDGG).lt.TINY) FDGG=0. ! 10.08.16: roundoff err.safeguard
      FDXG=-FinfDX*W-Finf*WDX-GAMI*(FinfDX*WDG+Finf*WDXG)
      S_SCR=-GAMI**2*Finf*WDG
      USCR=S_SCR+FSCR
      CVSCR=-GAMI**2*FDGG
      PSCR=(XSR*FDX+GAMI*FDG)/3.d0
      PDTSCR=GAMI**2*(XSR*Finf*(FinfX*WDG+WDXG)-FDGG)/3.d0
      PDRSCR=(12.d0*PSCR+XSR**2*FDXX+2.d0*XSR*GAMI*FDXG+
     +  GAMI**2*FDGG)/9.d0
      return
      end

      subroutine FHARM12(GAMI,TPT,
     *   Fharm,Uharm,Pharm,CVth,Sth,PDTharm,PDRharm)
* Thermodynamic functions of a harmonic crystal, incl.stat.Coul.lattice
*
*                                                       Version 27.04.12
* Stems from FHARM8 v.15.02.08
* Replaced HLfit8 with HLfit12: rearranged output.
* Input: GAMI - ionic Gamma, TPT=T_{p,i}/T
* Output: Fharm=F/(N_i T), Uharm=U/(N_i T), Pharm=P/(n_i T),
* CVth=C_V/N_i, Sharm=S/N_i
* PDTharm = Pharm + d Pharm / d ln T, PDRharm = Pharm + d Pharm/d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(CM=.895929256d0) ! Madelung constant for BCC lattice
      call HLfit12(TPT,F,U,CVth,Sth,U1,CW,1)
      U0=-CM*GAMI ! perfect lattice
      E0=1.5d0*U1*TPT ! zero-point energy
      Uth=U+E0
      Fth=F+E0
      Uharm=U0+Uth
      Fharm=U0+Fth
      Pharm=U0/3.d0+Uth/2.d0
      PDTharm=.5d0*CVth
      PDRharm=U0/2.25d0+.75d0*Uth-.25d0*CVth
      return
      end

      subroutine HLfit12(eta,F,U,CV,S,U1,CW,LATTICE)
*                                                       Version 24.04.12
* Stems from HLfit8 v.03.12.08;
*   differences: E0 excluded from  U and F;
*   U1 and d(CV)/d\ln(T) are added on the output.
* Fit to thermal part of the thermodynamic functions.
* Baiko, Potekhin, & Yakovlev (2001).
* Zero-point lattice quantum energy 1.5u_1\eta EXCLUDED (unlike HLfit8).
* Input: eta=Tp/T, LATTICE=1 for bcc, 2 for fcc
* Output: F and U (normalized to NkT) - due to phonon excitations,
*   CV and S (normalized to Nk) in the HL model,
*   U1 - the 1st phonon moment,
*   CW=d(CV)/d\ln(T)
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(EPS=1.d-5,TINY=1.d-99)
      if (LATTICE.eq.1) then ! bcc lattice
         CLM=-2.49389d0 ! 3*ln<\omega/\omega_p>
         U1=.5113875d0
         ALPHA=.265764d0
         BETA=.334547d0
         GAMMA=.932446d0
         A1=.1839d0
         A2=.593586d0
         A3=.0054814d0
         A4=5.01813d-4
         A6=3.9247d-7
         A8=5.8356d-11
         B0=261.66d0
         B2=7.07997d0
         B4=.0409484d0
         B5=.000397355d0
         B6=5.11148d-5
         B7=2.19749d-6
         C9=.004757014d0
         C11=.0047770935d0
      elseif (LATTICE.eq.2) then ! fcc lattice
         CLM=-2.45373d0
         U1=.513194d0
         ALPHA=.257591d0
         BETA=.365284d0
         GAMMA=.9167070d0
         A1=.0
         A2=.532535d0
         A3=.0
         A4=3.76545d-4
         A6=2.63013d-7
         A8=6.6318d-11
         B0=303.20d0
         B2=7.7255d0
         B4=.0439597d0
         B5=.000114295d0
         B6=5.63434d-5
         B7=1.36488d-6
         C9=.00492387d0
         C11=.00437506d0
      else
         stop'HLfit: unknown lattice type'
      endif
      if (eta.gt.1./EPS) then ! asymptote of Eq.(13) of BPY'01
         U=3./(C11*eta**3)
         F=-U/3.
         CV=4.*U
        goto 50
      elseif (eta.lt.EPS) then ! Eq.(17) of BPY'01
        if (eta.lt.TINY) stop'HLfit: eta is too small'
         F=3.*dlog(eta)+CLM-1.5*U1*eta+eta**2/24.
         U=3.-1.5*U1*eta+eta**2/12.
         CV=3.-eta**2/12.
         goto 50
      endif
      eta2=eta**2
      eta3=eta2*eta
      eta4=eta3*eta
      eta5=eta4*eta
      eta6=eta5*eta
      eta7=eta6*eta
      eta8=eta7*eta
      B9=A6*C9
      B11=A8*C11
      UP=1.+A1*eta+A2*eta2+A3*eta3+A4*eta4+A6*eta6+A8*eta8
      DN=B0+B2*eta2+B4*eta4+B5*eta5+B6*eta6+
     +  B7*eta7+eta8*(B9*eta+B11*eta3)
      EA=dexp(-ALPHA*eta)
      EB=dexp(-BETA*eta)
      EG=dexp(-GAMMA*eta)
      F=dlog(1.d0-EA)+dlog(1.d0-EB)+dlog(1.-EG)-UP/DN ! F_{thermal}/NT
      UP1=A1+
     + 2.*A2*eta+3.*A3*eta2+4.*A4*eta3+6.*A6*eta5+8.*A8*eta7
      UP2=2.*A2+6.*A3*eta+12.*A4*eta2+30.*A6*eta4+56.*A8*eta6
      UP3=6.*A3+24.*A4*eta+120.*A6*eta3+336*A8*eta5
      DN1=2.*B2*eta+4.*B4*eta3+5.*B5*eta4+6.*B6*eta5+
     +  7.*B7*eta6+eta8*(9.*B9+11.*B11*eta2)
      DN2=2.*B2+12.*B4*eta2+20.*B5*eta3+30.*B6*eta4+
     +  42.*B7*eta5+72.*B9*eta7+110.*B11*eta8*eta
      DN3=24.*B4*eta+60.*B5*eta2+120.*B6*eta3+
     +  210.*B7*eta4+504.*B9*eta6+990.*B11*eta8
      DF1=ALPHA*EA/(1.d0-EA)+BETA*EB/(1.d0-EB)+GAMMA*EG/(1.d0-EG)-
     -  (UP1*DN-DN1*UP)/DN**2 ! int.en./NT/eta = df/d\eta
      DF2=ALPHA**2*EA/(1.d0-EA)**2+BETA**2*EB/(1.d0-EB)**2+
     +  GAMMA**2*EG/(1.d0-EG)**2+
     +  ((UP2*DN-DN2*UP)*DN-2.*(UP1*DN-DN1*UP)*DN1)/DN**3 ! -d2f/d\eta^2
      U=DF1*eta
      CV=DF2*eta2
      DF3=-ALPHA**3*EA/(1.d0-EA)**3*(1.+EA)-
     -  BETA**3*EB/(1.d0-EB)**3*(1.+EB)-
     -  GAMMA**3*EG/(1.d0-EG)**3*(1.+EG)+
     +  UP3/DN-(3.*UP2*DN1+3.*UP1*DN2+UP*DN3)/DN**2+
     +  6.*DN1*(UP1*DN1+UP*DN2)/DN**3-6.*UP*DN1**3/DN**4 ! -d3f/d\eta^3
      CW=-2.*CV-eta3*DF3
   50 continue
      S=U-F
      return
      end

      subroutine ANHBC(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah)
*                                                       Version 09.06.22
* Anharmonic free energy in OCP crystal according to
*   Baiko & Chugunov, MNRAS 510, 2628 (2022).
* Input: GAMI - ionic Gamma, TPT=Tp/T - ionic quantum parameter
* Output: anharm.free en. Fah=F_{AH}/(N_i kT), internal energy Uah,
*   pressure Pah=P_{AH}/(n_i kT), specific heat CVah = C_{V,AH}/(N_i k),
*   PDTah = Pah + d Pah / d ln T, PDRah = Pah + d Pah / d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(TINY=1.d-19,A1CL=10.2d0,A2CL=248.d0,A3CL=2.03d5,
     *  A1Q=-.62d0/6.d0,A2Q=-.56d0,A3Q=2.35d0,A11=-10.d0,A12=6.d-3)
      data KINIT/0/
      if (TPT.le.0.) stop'LIQBC: non-positive TPT'
      if (GAMI.le.0.) stop'LIQBC: non-positive GAMI'
      Fah=0.
      Uah=0.
      Pah=0.
      CVah=0.
      PDTah=0.
      PDRah=0.
      if (TPT.lt.TINY) return
      if (KINIT.eq.0) then
         A13=-A1CL-A11
         A14=(A1Q-A11*A12)/A13
         A21=(-2.d0*A2Q/A2CL)**(4.d0/3.d0)
         KINIT=1
      endif
      GAMQ=GAMI/TPT
      TPT2=TPT**2
      TPT4=TPT2**2
      GAMQ2=GAMQ**2
      GAMQ3=GAMQ2*GAMQ
      A1U=-A11/TPT2*(3.d0*A12*TPT2+1.d0)/(1.d0+A12*TPT2)**2-
     -  A13/TPT2*(3.d0*A14*TPT2+1.d0)/(1.d0+A14*TPT2)**2+A1Q
      A2U=(A2CL/TPT**3)*(2.d0-A21*TPT4)*.5d0/dsqrt(dsqrt(1.d0+A21*TPT4))
      A3U=A3CL/TPT4+A3Q
      UahBC=A1U/GAMQ2+A2U/GAMQ3+A3U/GAMQ2**2 ! BC'22(5)
      Uah=UahBC*GAMI
      DN12=1.d0+A12*TPT2
      A1F1=A11/TPT2/DN12
      A1F1D=A1F1*(-2.d0-2.d0*A12*TPT2/DN12)
      DN14=1.d0+A14*TPT2
      A1F2=A13/TPT2/DN14
      A1F2D=A1F2*(-2.d0-2.d0*A14*TPT2/DN14)
      A1F=A1F1+A1F2+A1Q
      A1FD=A1F1D+A1F2D
      DN2F=1.d0+A21*TPT4
      A2F=-A2CL*.5d0*(dsqrt(dsqrt(DN2F))/TPT)**3
      A2FD=A2F*(-3.d0+3.d0*A21*TPT4/DN2F)
      A3F=-A3CL/3.d0/TPT4+A3Q
      A3FD=A3CL/.75d0/TPT4
      FahBC=A1F/GAMQ**2+A2F/GAMQ**3+A3F/GAMQ**4 ! BC'22(15)
      Fah=FahBC*GAMI
      A1S1=-2.d0*A11/TPT*(2.d0*A12*TPT2+1.d0)/DN12**2
      A1S2=-2.d0*A13/TPT*(2.d0*A14*TPT2+1.d0)/DN14**2
      A1S=A1S1+A1S2
      A1S1D=A1S1*(-1.d0+4.d0*A12*
     *  (TPT2/(2.d0*A12*TPT2+1.d0)-TPT2/DN12))
      A1S2D=A1S2*(-1.d0+4.d0*A14*
     *  (TPT2/(2.d0*A14*TPT2+1.d0)-TPT2/DN14))
      A1SD=A1S1D+A1S2D
      A2S=1.5d0*A2CL/TPT2/dsqrt(dsqrt(DN2F))
      A2SD=A2S*(-2.d0-A21*TPT4/(DN2F))
      A3S=A3CL/(.75d0*TPT**3)
      A3SD=-3.d0*A3S
      Sah=A1S/GAMQ+A2S/GAMQ2+A3S/GAMQ3
      SahDR=(A1S/6.d0+.5d0*A1SD)/GAMQ+
     +  (A2S/3.d0+.5d0*A2SD)/GAMQ2+
     +  (A3S+A3SD)*.5d0/GAMQ3
      SahDT=-A1SD/GAMQ-A2SD/GAMQ2-A3SD/GAMQ3
      Pah=Sah*.5d0+TPT*(A1F/GAMQ/1.5d0+A2F/1.2d0/GAMQ2+A3F/GAMQ3)
      PahDR=SahDR*.5d0+TPT*((A1F/2.25d0+A1FD/3.d0)/GAMQ+
     +  (A2FD/2.4d0+A2F/1.44d0)/GAMQ2+(A3FD*.5d0+A3F)/GAMQ3)
      PahDT=SahDT*.5d0-TPT*((A1FD+A1F)/1.5d0/GAMQ+
     +  (A2FD+A2F)/1.2d0/GAMQ2+(A3FD+A3F)/GAMQ3)
      PDRah=Pah+PahDR
      PDTah=Pah+PahDT
      CVah=SahDT
      return
      end

      subroutine CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321,
     *  FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
*                                                       Version 02.07.09
* Correction to the linear mixing rule for moderate to small Gamma
* Input: RS=r_s (if RS=0, then OCP, otherwise EIP)
*        GAME=\Gamma_e
*        Zmean=<Z> (average Z of all ions, without electrons)
*        Z2mean=<Z^2>, Z52=<Z^2.5>, Z53=<Z^{5/3}>, Z321=<Z(Z+1)^1.5>
* Output: FMIX=\Delta f - corr.to the reduced free energy f=F/N_{ion}kT
*         UMIX=\Delta u - corr.to the reduced internal energy u
*         PMIX=\Delta u - corr.to the reduced pressure P=P/n_{ion}kT
*         CVMIX=\Delta c - corr.to the reduced heat capacity c_V
*         PDTMIX=(1/n_{ion}kT)d\Delta P / d ln T
*               = \Delta p +  d \Delta p / d ln T
*         PDRMIX=(1/n_{ion}kT)d\Delta P / d ln n_e
* (composition is assumed fixed: Zmean,Z2mean,Z52,Z53=constant)
      implicit double precision (A-H), double precision (O-Z)
      parameter (TINY=1.d-9)
      GAMImean=GAME*Z53
      if (RS.lt.TINY) then ! OCP
         Dif0=Z52-dsqrt(Z2mean**3/Zmean)
      else
         Dif0=Z321-dsqrt((Z2mean+Zmean)**3/Zmean)
      endif
      DifR=Dif0/Z52
      DifFDH=Dif0*GAME*sqrt(GAME/3.) ! F_DH - F_LM(DH)
      D=Z2mean/Zmean**2
      if (dabs(D-1.d0).lt.TINY) then ! no correction
         FMIX=0.
         UMIX=0.
         PMIX=0.
         CVMIX=0.
         PDTMIX=0.
         PDRMIX=0.
         return
      endif
      P3=D**(-0.2)
      D0=(2.6*DifR+14.*DifR**3)/(1.d0-P3)
      GP=D0*GAMImean**P3
      FMIX0=DifFDH/(1.+GP)
      Q=D**2*.0117
      R=1.5/P3-1.
      GQ=Q*GP
      FMIX=FMIX0/(1.+GQ)**R
      G=1.5-P3*GP/(1.+GP)-R*P3*GQ/(1.+GQ)
      UMIX=FMIX*G
      PMIX=UMIX/3.d0
      GDG=-P3**2*(GP/(1.d0+GP)**2+R*GQ/(1.d0+GQ)**2) ! d G /d ln Gamma
      UDG=UMIX*G+FMIX*GDG ! d u_mix /d ln Gamma
      CVMIX=UMIX-UDG
      PDTMIX=PMIX-UDG/3.
      PDRMIX=PMIX+UDG/9.
      return
      end

* ===================  IDEAL ELECTRON GAS  =========================== *
      subroutine ELECT11(TEMP,CHI,
     *  DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,
     *  DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
*                                                       Version 17.11.11
*                 safeguard against huge (-CHI) values is added 27.05.17
* ELECT9 v.04.03.09 + smooth match of two fits at chi >> 1 + add.outputs
* Compared to ELECTRON v.06.07.00, this S/R is completely rewritten:
*        numerical differentiation is avoided now.
* Compared to ELECT7 v.06.06.07,
*    - call BLIN7 is changed to call BLIN9,
*    - Sommerfeld expansion is used at chi >~ 28 i.o. 1.e4
*    - Sommerfeld expansion is corrected: introduced DeltaEF, D1 and D2.
* Ideal electron-gas EOS.
* Input: TEMP - T [a.u.], CHI=\mu/T
* Output: DENS - electron number density n_e [a.u.],
*         FEid - free energy / N_e kT, UEid - internal energy / N_e kT,
*         PEid - pressure (P_e) / n_e kT, SEid - entropy / N_e k,
*         CVE - heat capacity / N_e k,
*         CHITE=(d ln P_e/d ln T)_V, CHIRE=(d ln P_e/d ln n_e)_T
*         DlnDH=(d ln n_e/d CHI)_T = (T/n_e) (d n_e/d\mu)_T
*         DlnDT=(d ln n_e/d ln T)_CHI
*         DlnDHH=(d^2 ln n_e/d CHI^2)_T
*         DlnDTT=(d^2 ln n_e/d (ln T)^2)_CHI
*         DlnDHT=d^2 ln n_e/d (ln T) d CHI
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (CHI2=28.d0,XMAX=20.d0)
      parameter (DCHI2=CHI2-1.d0)
      parameter (XSCAL2=XMAX/DCHI2)
      if (CHI.lt.-1.d2) stop'ELECT11: too large negative CHI' ! 27.05.17
      X2=(CHI-CHI2)*XSCAL2
      if (X2.lt.-XMAX) then
         call ELECT11a(TEMP,CHI,
     *     DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,
     *     DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
      elseif (X2.gt.XMAX) then
         call ELECT11b(TEMP,CHI,
     *     DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,
     *     DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
      else
         call FERMI10(X2,XMAX,FP,FM)
         call ELECT11a(TEMP,CHI,
     *     DENSa,FEida,PEida,UEida,SEida,CVEa,CHITEa,CHIREa,
     *     DlnDHa,DlnDTa,DlnDHHa,DlnDTTa,DlnDHTa)
         call ELECT11b(TEMP,CHI,
     *     DENSb,FEidb,PEidb,UEidb,SEidb,CVEb,CHITEb,CHIREb,
     *     DlnDHb,DlnDTb,DlnDHHb,DlnDTTb,DlnDHTb)
         DENS=DENSa*FP+DENSb*FM
         FEid=FEida*FP+FEidb*FM
         PEid=PEida*FP+PEidb*FM
         UEid=UEida*FP+UEidb*FM
         SEid=SEida*FP+SEidb*FM
         CVE=CVEa*FP+CVEb*FM
         CHITE=CHITEa*FP+CHITEb*FM
         CHIRE=CHIREa*FP+CHIREb*FM
         DlnDH=DlnDHa*FP+DlnDHb*FM
         DlnDT=DlnDTa*FP+DlnDTb*FM
         DlnDHH=DlnDHHa*FP+DlnDHHb*FM
         DlnDHT=DlnDHTa*FP+DlnDHTb*FM
         DlnDTT=DlnDTTa*FP+DlnDTTb*FM
      endif
      return
      end

      subroutine ELECT11a(TEMP,CHI,
     *  DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,
     *  DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
*                                                       Version 16.11.11
* This is THE FIRST PART of ELECT9 v.04.03.09.
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (BOHR=137.03599908d0,PI=3.141592653d0)
      parameter (PI2=PI**2,BOHR2=BOHR**2,BOHR3=BOHR2*BOHR) !cleaned 15/6
      TEMR=TEMP/BOHR2 ! T in rel.units (=T/mc^2)
      call BLIN9(TEMR,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
      TPI=TEMR*dsqrt(2.d0*TEMR)/PI2 ! common pre-factor
      DENR=TPI*(W1*TEMR+W0)
      PR=TEMR*TPI/3.*(W2*TEMR+2.*W1)
      U=TEMR*TPI*(W2*TEMR+W1)
* (these are density, pressure, and internal energy in the rel.units)
      PEid=PR/(DENR*TEMR)
      UEid=U/(DENR*TEMR)
      FEid=CHI-PEid
      DENS=DENR*BOHR3 ! converts from rel.units to a.u.
      SEid=UEid-FEid
* derivatives over T at constant chi:
      dndT=TPI*(1.5*W0/TEMR+2.5*W1+W0DT+TEMR*W1DT) ! (d n_e/dT)_\chi
      dPdT=TPI/3.*(5.*W1+2.*TEMR*W1DT+3.5*TEMR*W2+TEMR**2*W2DT)!dP/dT
      dUdT=TPI*(2.5*W1+TEMR*W1DT+3.5*TEMR*W2+TEMR**2*W2DT)!dU/dT_\chi
* derivatives over chi at constant T and second derivatives:
      dndH=TPI*(W0DX+TEMR*W1DX) ! (d n_e/d\chi)_T
      dndHH=TPI*(W0DXX+TEMR*W1DXX) ! (d^2 n_e/d\chi)_T
      dndTT=TPI*(.75*W0/TEMR**2+3.*W0DT/TEMR+W0DTT+
     +  3.75*W1/TEMR+5.*W1DT+TEMR*W1DTT)
      dndHT=TPI*(1.5*W0DX/TEMR+W0DXT+2.5*W1DX+TEMR*W1DXT)
      DlnDH=dndH/DENR ! (d ln n_e/d\chi)_T
      DlnDT=dndT*TEMR/DENR ! (d ln n_e/d ln T)_\chi
      DlnDHH=dndHH/DENR-DlnDH**2 ! (d^2 ln n_e/d\chi^2)_T
      DlnDTT=TEMR**2/DENR*dndTT+DlnDT-DlnDT**2 ! d^2 ln n_e/d ln T^2
      DlnDHT=TEMR/DENR*(dndHT-dndT*DlnDH) ! d^2 ln n_e/d\chi d ln T
      dPdH=TPI/3.*TEMR*(2.*W1DX+TEMR*W2DX) ! (d P_e/d\chi)_T
      dUdH=TPI*TEMR*(W1DX+TEMR*W2DX) ! (d U_e/d\chi)_T
      CVE=(dUdT-dUdH*dndT/dndH)/DENR
      CHITE=TEMR/PR*(dPdT-dPdH*dndT/dndH)
      CHIRE=DENR/PR*dPdH/dndH ! (dndH*TEMR*PEid) ! DENS/PRE*dPdH/dndH
      return
      end

      subroutine ELECT11b(TEMP,CHI,
     *  DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,
     *  DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
*                                                       Version 17.11.11
* Stems from ELECT9b v.19.01.10, Diff. - additional output.
* Sommerfeld expansion at very large CHI.
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (BOHR=137.03599908d0,PI=3.141592653d0)
      parameter (PI2=PI**2,BOHR2=BOHR**2,BOHR3=BOHR2*BOHR) !cleaned 15/6
      TEMR=TEMP/BOHR2 ! T in rel.units (=T/mc^2)
      EF=CHI*TEMR ! Fermi energy in mc^2 - zeroth aprox. = CMU1
      DeltaEF=PI2*TEMR**2/6.d0*(1.d0+2.d0*EF*(2.d0+EF))/
     /  (EF*(1.d0+EF)*(2.d0+EF)) ! corr. [p.125, equiv.Eq.(6) of PC'10]
      EF=EF+DeltaEF ! corrected Fermi energy (14.02.09)
      G=1.d0+EF ! electron Lorentz-factor
      if (EF.gt.1.d-5) then ! relativistic expansion (Yak.&Shal.'89)
        PF=dsqrt(G**2-1.d0) ! Fermi momentum [rel.un.=mc]
        F=(PF*(1.+2.d0*PF**2)*G-PF**3/.375d0-dlog(PF+G))/8.d0/PI2!F/V
        DF=-TEMR**2*PF*G/6.d0 ! thermal correction to F/V
        P=(PF*G*(PF**2/1.5d0-1.d0)+dlog(PF+G))/8.d0/PI2 ! P(T=0)
        DP=TEMR**2*PF*(PF**2+2.d0)/G/18.d0 ! thermal correction to P
        CVE=PI2*TEMR*G/PF**2
      else ! nonrelativistic limit
        PF=dsqrt(2.d0*EF)
        F=PF**5*0.1d0/PI2
        DF=-TEMR**2*PF/6.d0
        P=F/1.5d0
        DP=TEMR**2*PF/9.d0
        CVE=PI2*TEMR/EF/2.d0
      endif
      F=F+DF
      P=P+DP
      S=-2.d0*DF ! entropy per unit volume [rel.un.]
      U=F+S
      CHIRE=PF**5/(9.d0*PI2*P*G)
      CHITE=2.d0*DP/P
      DENR=PF**3/3.d0/PI2 ! n_e [rel.un.=\Compton^{-3}]
      DENS=DENR*BOHR3 ! conversion to a.u.(=\Bohr_radius^{-3})
* derivatives over chi at constant T and T at constant chi:
      TPI=TEMR*dsqrt(2.d0*TEMR)/PI2 ! common pre-factor
      call SOMMERF(TEMR,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
      dndH=TPI*(W0DX+TEMR*W1DX) ! (d n_e/d\chi)_T
      dndT=TPI*(1.5*W0/TEMR+2.5*W1+W0DT+TEMR*W1DT) ! (d n_e/dT)_\chi
      dndHH=TPI*(W0DXX+TEMR*W1DXX) ! (d^2 n_e/d\chi)_T
      dndTT=TPI*(.75*W0/TEMR**2+3.*W0DT/TEMR+W0DTT+
     +  3.75*W1/TEMR+5.*W1DT+TEMR*W1DTT)
      dndHT=TPI*(1.5*W0DX/TEMR+W0DXT+2.5*W1DX+TEMR*W1DXT)
      DlnDH=dndH/DENR ! (d ln n_e/d\chi)_T
      DlnDT=dndT*TEMR/DENR ! (d ln n_e/d ln T)_\chi
      DlnDHH=dndHH/DENR-DlnDH**2 ! (d^2 ln n_e/d\chi^2)_T
      DlnDTT=TEMR**2/DENR*dndTT+DlnDT-DlnDT**2 ! d^2 ln n_e/d ln T^2
      DlnDHT=TEMR/DENR*(dndHT-dndT*DlnDH) ! d^2 ln n_e/d\chi d ln T
      DT=DENR*TEMR
      PEid=P/DT
      UEid=U/DT
      FEid=F/DT
      SEid=S/DT
* Empirical corrections of 16.02.09:
      D1=DeltaEF/EF
      D2=D1*(4.d0-2.d0*(PF/G))
      CVE=CVE/(1.d0+D2)
      SEid=SEid/(1.d0+D1)
      CHITE=CHITE/(1.d0+D2)
      return
      end

      subroutine SOMMERF(TEMR,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
*                                                       Version 17.11.11
* Sommerfeld expansion for the Fermi-Dirac integrals
* Input: TEMR=T/mc^2; CHI=(\mu-mc^2)/T
* Output: Wk - Fermi-Dirac integral of the order k+1/2
*         WkDX=dWk/dCHI, WkDT = dWk/dT, WkDXX=d^2 Wk / d CHI^2,
*         WkDTT=d^2 Wk / d T^2, WkDXT=d^2 Wk /dCHIdT,
*         W0XXX=d^3 W0 / d CHI^3, W0XTT=d^3 W0 /(d CHI d^2 T),
*         W0XXT=d^3 W0 /dCHI^2 dT
* [Draft source: yellow book pages 124-127]
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(PI=3.141592653d0)
      parameter(PI2=PI**2)
      if (CHI.lt..5d0) stop'SOMMERF: non-degenerate (small CHI)'
      if (TEMR.le.0.d0) stop'SOMMERF: T < 0'
      CMU1=CHI*TEMR ! chemical potential in rel.units
      CMU=1.d0+CMU1
      call SUBFERMJ(CMU1,
     *  CJ00,CJ10,CJ20,
     *  CJ01,CJ11,CJ21,
     *  CJ02,CJ12,CJ22,
     *  CJ03,CJ13,CJ23,
     *  CJ04,CJ14,CJ24,CJ05)
      PIT26=(PI*TEMR)**2/6.d0
CCC      PITAU4=PIT26**2*0.7d0
      CN0=dsqrt(.5d0/TEMR)/TEMR
      CN1=CN0/TEMR
      CN2=CN1/TEMR
      W0=CN0*(CJ00+PIT26*CJ02) ! +CN0*PITAU4*CJ04
      W1=CN1*(CJ10+PIT26*CJ12) ! +CN1*PITAU4*CJ14
      W2=CN2*(CJ20+PIT26*CJ22) ! +CN2*PITAU4*CJ24
      W0DX=CN0*TEMR*(CJ01+PIT26*CJ03) ! +CN0*PITAU4*CJ05
      W1DX=CN0*(CJ11+PIT26*CJ13)
      W2DX=CN1*(CJ21+PIT26*CJ23)
      W0DT=CN1*(CMU1*CJ01-1.5d0*CJ00+PIT26*(CMU1*CJ03+.5d0*CJ02))
CCC     +  CN1*PITAU4*(CMU1*CJ05+2.5d0*CJ04)
      W1DT=CN2*(CMU1*CJ11-2.5d0*CJ10+PIT26*(CMU1*CJ13-.5d0*CJ12))
      W2DT=CN2/TEMR*(CMU1*CJ21-3.5d0*CJ20+PIT26*(CMU1*CJ23-1.5d0*CJ22))
      W0DXX=CN0*TEMR**2*(CJ02+PIT26*CJ04)
      W1DXX=CN0*TEMR*(CJ12+PIT26*CJ14)
      W2DXX=CN0*(CJ22+PIT26*CJ24)
      W0DXT=CN0*(CMU1*CJ02-.5d0*CJ01+PIT26*(CMU1*CJ04+1.5d0*CJ03))
      W1DXT=CN1*(CMU1*CJ12-1.5d0*CJ11+PIT26*(CMU1*CJ14+.5d0*CJ13))
      W2DXT=CN2*(CMU1*CJ22-2.5d0*CJ21+PIT26*(CMU1*CJ24-.5d0*CJ23))
      W0DTT=CN2*(3.75d0*CJ00-3.d0*CMU1*CJ01+CMU1**2*CJ02+
     +  PIT26*(-.25d0*CJ02+CMU1*CJ03+CMU1**2*CJ04))
      W1DTT=CN2/TEMR*(8.75d0*CJ10-5.d0*CMU1*CJ11+CMU1**2*CJ12+
     +  PIT26*(.75d0*CJ12-CMU1*CJ13+CMU1**2*CJ14))
      W2DTT=CN2/TEMR**2*(15.75d0*CJ20-7.d0*CMU1*CJ21+CMU1**2*CJ22+
     +  PIT26*(3.75d0*CJ22-3.d0*CMU1*CJ23+CMU1**2*CJ24))
      W0XXX=CN0*TEMR**3*(CJ03+PIT26*CJ05)
      W0XXT=CN0*TEMR*(CMU1*CJ03+.5d0*CJ02+PIT26*(CMU1*CJ05+2.5d0*CJ04))
      W0XTT=CN1*(.75d0*CJ01-CMU1*CJ02+CMU1**2*CJ03+
     +  PIT26*(.75d0*CJ03+3.d0*CMU1*CJ04+CMU1**2*CJ05))
      return
      end

      subroutine SUBFERMJ(CMU1,
     *  CJ00,CJ10,CJ20,
     *  CJ01,CJ11,CJ21,
     *  CJ02,CJ12,CJ22,
     *  CJ03,CJ13,CJ23,
     *  CJ04,CJ14,CJ24,CJ05)
*                                                       Version 17.11.11
*                                                     corrected 04.03.21
* Supplement to SOMMERF
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(EPS=1.d-4) ! inserted 04.03.21
      if (CMU1.le.0.d0) stop'SUBFERMJ: small CHI'
      CMU=1.d0+CMU1
      X0=dsqrt(CMU1*(2.d0+CMU1))
      X3=X0**3
      X5=X0**5
      if (X0.lt.EPS) then
         CJ00=X3/3.d0
         CJ10=.1d0*X5
         CJ20=X0**7/28.d0
      else
         CL=dlog(X0+CMU)
         CJ00=.5d0*(X0*CMU-CL) ! J_{1/2}^0
         CJ10=X3/3.d0-CJ00 ! J_{3/2}^0
         CJ20=(.75d0*CMU-2.d0)/3.d0*X3+1.25d0*CJ00 ! J_{5/2}^0
      endif
      CJ01=X0 ! J_{1/2}^1
      CJ11=CJ01*CMU1 ! J_{3/2}^1
      CJ21=CJ11*CMU1 ! J_{5/2}^1
      CJ02=CMU/X0 ! J_{1/2}^2
      CJ12=CMU1/X0*(3.d0+2.d0*CMU1) ! J_{3/2}^2
      CJ22=CMU1**2/X0*(5.d0+3.d0*CMU1) ! J_{5/2}^2
      CJ03=-1.d0/X3 ! J_{1/2}^3
      CJ13=CMU1/X3*(2.d0*CMU1**2+6.d0*CMU1+3.d0)
      CJ23=CMU1**2/X3*(6.d0*CMU1**2+2.d1*CMU1+1.5d1)
      CJ04=3.d0*CMU/X5
      CJ14=-3.d0*CMU1/X5
      CJ24=CMU1**2/X5*(6.d0*CMU1**3+3.d1*CMU1**2+45.d0*CMU1+15.d0)
      CJ05=(-12.d0*CMU1**2-24.d0*CMU1-15.d0)/(X5*X0**2)
      return
      end

* ==============  ELECTRON EXCHANGE AND CORRELATION   ================ *
      subroutine EXCOR7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC)
*                                                       Version 09.06.07
*                             Accuracy-loss cut-off modified on 10.08.16
* Exchange-correlation contribution for the electron gas
* Stems from TANAKA1 v.03.03.96. Added derivatives.
* Input: RS - electron density parameter =electron-sphere radius in a.u.
*        GAME - electron Coulomb coupling parameter
* Output: FXC - excess free energy of e-liquid per kT per one electron
*             according to Tanaka & Ichimaru 85-87 and Ichimaru 93
*         UXC - internal energy contr.[per 1 electron, kT]
*         PXC - pressure contribution divided by (n_e kT)
*         CVXC - heat capacity divided by N_e k
*         SXC - entropy divided by N_e k
*         PDTXC,PDRXC = PXC + d PXC / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      parameter(EPS=1.d-8) ! 10.08.16
      THETA=.543*RS/GAME ! non-relativistic degeneracy parameter
      SQTH=dsqrt(THETA)
      THETA2=THETA**2
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
      if (THETA.gt..005) then
         CHT1=dcosh(1.d0/THETA)
         SHT1=dsinh(1.d0/THETA)
         CHT2=dcosh(1.d0/SQTH)
         SHT2=dsinh(1.d0/SQTH)
         T1=SHT1/CHT1 ! dtanh(1.d0/THETA)
         T2=SHT2/CHT2 ! dtanh(1./dsqrt(THETA))
         T1DH=-1./(THETA*CHT1)**2 ! d T1 / d\theta
         T1DHH=2./(THETA*CHT1)**3*(CHT1-SHT1/THETA)
         T2DH=-.5*SQTH/(THETA*CHT2)**2
         T2DHH=(.75*SQTH*CHT2-.5*SHT2)/(THETA*CHT2)**3
      else
         T1=1.
         T2=1.
         T1DH=0.
         T2DH=0.
         T1DHH=0.
         T2DHH=0.
      endif
      A0=.75+3.04363*THETA2-.09227*THETA3+1.7035*THETA4
      A0DH=6.08726*THETA-.27681*THETA2+6.814*THETA3
      A0DHH=6.08726-.55362*THETA+20.442*THETA2
      A1=1.+8.31051*THETA2+5.1105*THETA4
      A1DH=16.62102*THETA+20.442*THETA3
      A1DHH=16.62102+61.326*THETA2
      A=.610887*A0/A1*T1 ! HF fit of Perrot and Dharma-wardana
      AH=A0DH/A0-A1DH/A1+T1DH/T1
      ADH=A*AH
      ADHH=ADH*AH+A*(A0DHH/A0-(A0DH/A0)**2-A1DHH/A1+(A1DH/A1)**2+
     +  T1DHH/T1-(T1DH/T1)**2)
      B0=.341308+12.070873d0*THETA2+1.148889d0*THETA4
      B0DH=24.141746d0*THETA+4.595556d0*THETA3
      B0DHH=24.141746d0+13.786668d0*THETA2
      B1=1.+10.495346d0*THETA2+1.326623*THETA4
      B1DH=20.990692d0*THETA+5.306492*THETA3
      B1DHH=20.990692d0+15.919476d0*THETA2
      B=SQTH*T2*B0/B1
      BH=.5/THETA+T2DH/T2+B0DH/B0-B1DH/B1
      BDH=B*BH
      BDHH=BDH*BH+B*(-.5/THETA2+T2DHH/T2-(T2DH/T2)**2+
     +  B0DHH/B0-(B0DH/B0)**2-B1DHH/B1+(B1DH/B1)**2)
      D0=.614925+16.996055d0*THETA2+1.489056*THETA4
      D0DH=33.99211d0*THETA+5.956224d0*THETA3
      D0DHH=33.99211d0+17.868672d0*THETA2
      D1=1.+10.10935*THETA2+1.22184*THETA4
      D1DH=20.2187*THETA+4.88736*THETA3
      D1DHH=20.2187+14.66208*THETA2
      D=SQTH*T2*D0/D1
      DH=.5/THETA+T2DH/T2+D0DH/D0-D1DH/D1
      DDH=D*DH
      DDHH=DDH*DH+D*(-.5/THETA2+T2DHH/T2-(T2DH/T2)**2+
     +  D0DHH/D0-(D0DH/D0)**2-D1DHH/D1+(D1DH/D1)**2)
      E0=.539409+2.522206*THETA2+.178484*THETA4
      E0DH=5.044412*THETA+.713936*THETA3
      E0DHH=5.044412+2.141808*THETA2
      E1=1.+2.555501*THETA2+.146319*THETA4
      E1DH=5.111002*THETA+.585276*THETA3
      E1DHH=5.111002+1.755828*THETA2
      E=THETA*T1*E0/E1
      EH=1./THETA+T1DH/T1+E0DH/E0-E1DH/E1
      EDH=E*EH
      EDHH=EDH*EH+E*(T1DHH/T1-(T1DH/T1)**2+E0DHH/E0-(E0DH/E0)**2-
     -  E1DHH/E1+(E1DH/E1)**2-1./THETA2)
      EXP1TH=dexp(-1./THETA)
      C=(.872496+.025248*EXP1TH)*E
      CDH=.025248*EXP1TH/THETA2*E+C*EDH/E
      CDHH=.025248*EXP1TH/THETA2*(EDH+(1.-2.*THETA)/THETA2*E)+
     +  CDH*EDH/E+C*EDHH/E-C*(EDH/E)**2
      DISCR=dsqrt(4.*E-D**2)
      DIDH=.5/DISCR*(4.*EDH-2.*D*DDH)
      DIDHH=(-((2.*EDH-D*DDH)/DISCR)**2+2.*EDHH-DDH**2-D*DDHH)/DISCR
      S1=-C/E*GAME
      S1H=CDH/C-EDH/E
      S1DH=S1*S1H
      S1DHH=S1DH*S1H+S1*(CDHH/C-(CDH/C)**2-EDHH/E+(EDH/E)**2)
      S1DG=-C/E ! => S1DGG=0
      S1DHG=S1DG*(CDH/C-EDH/E)
      B2=B-C*D/E
      B2DH=BDH-(CDH*D+C*DDH)/E+C*D*EDH/E**2
      B2DHH=BDHH-(CDHH*D+2.*CDH*DDH+C*DDHH)/E+
     +  (2.*(CDH*D+C*DDH-C*D*EDH/E)*EDH+C*D*EDHH)/E**2
      SQGE=dsqrt(GAME)
      S2=-2./E*B2*SQGE
      S2H=B2DH/B2-EDH/E
      S2DH=S2*S2H
      S2DHH=S2DH*S2H+S2*(B2DHH/B2-(B2DH/B2)**2-EDHH/E+(EDH/E)**2)
      S2DG=.5*S2/GAME
      S2DGG=-.5*S2DG/GAME
      S2DHG=.5*S2DH/GAME
      R3=E*GAME+D*SQGE+1.
      R3DH=EDH*GAME+DDH*SQGE
      R3DHH=EDHH*GAME+DDHH*SQGE
      R3DG=E+.5*D/SQGE
      R3DGG=-.25*D/(GAME*SQGE)
      R3DHG=EDH+.5*DDH/SQGE
      B3=A-C/E
      B3DH=ADH-CDH/E+C*EDH/E**2
      B3DHH=ADHH-CDHH/E+(2.*CDH*EDH+C*EDHH)/E**2-2.*C*EDH**2/E**3
      C3=(D/E*B2-B3)/E ! =D*B2/E**2-B3/E
      C3DH=(DDH*B2+D*B2DH+B3*EDH)/E**2-2.*D*B2*EDH/E**3-B3DH/E
      C3DHH=(-B3DHH+
     +  (DDHH*B2+2.*DDH*B2DH+D*B2DHH+B3DH*EDH+B3*EDHH+B3DH*EDH)/E-
     -  2.*((DDH*B2+D*B2DH+B3*EDH+DDH*B2+D*B2DH)*EDH+D*B2*EDHH)/E**2+
     +  6.*D*B2*EDH**2/E**3)/E
      S3=C3*dlog(R3)
      S3DH=S3*C3DH/C3+C3*R3DH/R3
      S3DHH=(S3DH*C3DH+S3*C3DHH)/C3-S3*(C3DH/C3)**2+
     +  (C3DH*R3DH+C3*R3DHH)/R3-C3*(R3DH/R3)**2
      S3DG=C3*R3DG/R3
      S3DGG=C3*(R3DGG/R3-(R3DG/R3)**2)
      S3DHG=(C3DH*R3DG+C3*R3DHG)/R3-C3*R3DG*R3DH/R3**2
      B4=2.-D**2/E
      B4DH=EDH*(D/E)**2-2.*D*DDH/E
      B4DHH=EDHH*(D/E)**2+2.*EDH*(D/E)**2*(DDH/D-EDH/E)-
     -  2.*(DDH**2+D*DDHH)/E+2.*D*DDH*EDH/E**2
      C4=2.*E*SQGE+D
      C4DH=2.*EDH*SQGE+DDH
      C4DHH=2.*EDHH*SQGE+DDHH
      C4DG=E/SQGE
      C4DGG=-.5*E/(GAME*SQGE)
      C4DHG=EDH/SQGE
      S4A=2./E/DISCR
      S4AH=EDH/E+DIDH/DISCR
      S4ADH=-S4A*S4AH
      S4ADHH=-S4ADH*S4AH-
     -  S4A*(EDHH/E-(EDH/E)**2+DIDHH/DISCR-(DIDH/DISCR)**2)
      S4B=D*B3+B4*B2
      S4BDH=DDH*B3+D*B3DH+B4DH*B2+B4*B2DH
      S4BDHH=DDHH*B3+2.*DDH*B3DH+D*B3DHH+B4DHH*B2+2.*B4DH*B2DH+B4*B2DHH
      S4C=datan(C4/DISCR)-datan(D/DISCR)
      UP1=C4DH*DISCR-C4*DIDH
      DN1=DISCR**2+C4**2
      UP2=DDH*DISCR-D*DIDH
      DN2=DISCR**2+D**2
      S4CDH=UP1/DN1-UP2/DN2
      S4CDHH=(C4DHH*DISCR-C4*DIDHH)/DN1-
     -  UP1*2.*(DISCR*DIDH+C4*C4DH)/DN1**2-
     -  (DDHH*DISCR-D*DIDHH)/DN2+UP2*2.*(DISCR*DIDH+D*DDH)/DN2**2
      S4CDG=C4DG*DISCR/DN1
      S4CDGG=C4DGG*DISCR/DN1-2.*C4*DISCR*(C4DG/DN1)**2
      S4CDHG=(C4DHG*DISCR+C4DG*DIDH-
     -  C4DG*DISCR/DN1*2.*(DISCR*DIDH+C4*C4DH))/DN1
      S4=S4A*S4B*S4C
      S4DH=S4ADH*S4B*S4C+S4A*S4BDH*S4C+S4A*S4B*S4CDH
      S4DHH=S4ADHH*S4B*S4C+S4A*S4BDHH*S4C+S4A*S4B*S4CDHH+
     +  2.*(S4ADH*S4BDH*S4C+S4ADH*S4B*S4CDH+S4A*S4BDH*S4CDH)
      S4DG=S4A*S4B*S4CDG
      S4DGG=S4A*S4B*S4CDGG
      S4DHG=S4A*S4B*S4CDHG+S4CDG*(S4ADH*S4B+S4A*S4BDH)
      FXC=S1+S2+S3+S4
      FXCDH=S1DH+S2DH+S3DH+S4DH
      FXCDG=S1DG+S2DG+S3DG+S4DG
      FXCDHH=S1DHH+S2DHH+S3DHH+S4DHH
      FXCDGG=S2DGG+S3DGG+S4DGG
      FXCDHG=S1DHG+S2DHG+S3DHG+S4DHG
      PXC=(GAME*FXCDG-2.*THETA*FXCDH)/3.
      UXC=GAME*FXCDG-THETA*FXCDH
      SXC=(GAME*S2DG-S2+GAME*S3DG-S3+S4A*S4B*(GAME*S4CDG-S4C))-
     -  THETA*FXCDH
      if (dabs(SXC).lt.EPS*dabs(THETA*FXCDH)) SXC=0. ! accuracy loss
      CVXC=2.*THETA*(GAME*FXCDHG-FXCDH)-THETA**2*FXCDHH-GAME**2*FXCDGG
      if (dabs(CVXC).lt.EPS*dabs(GAME**2*FXCDGG)) CVXC=0. ! accuracy
      PDLH=THETA*(GAME*FXCDHG-2.*FXCDH-2.*THETA*FXCDHH)/3.
      PDLG=GAME*(FXCDG+GAME*FXCDGG-2.*THETA*FXCDHG)/3.
      PDRXC=PXC+(PDLG-2.*PDLH)/3.
      PDTXC=GAME*(THETA*FXCDHG-GAME*FXCDGG/3.)-
     -  THETA*(FXCDH/.75+THETA*FXCDHH/1.5)
      return
      end

* ======================  AUXILIARY SUBROUTINES   ==================== *
      subroutine FERINV7(F,N,X,XDF,XDFF) ! Inverse Fermi intergals
*                                                       Version 24.05.07
* X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84, 101
* q=N-1/2=-1/2,1/2,3/2,5/2 (N=0,1,2,3)
* Input: F - argument, N=q+1/2
* Output: X=X_q, XDF=dX/df, XDFF=d^2 X / df^2
* Relative error: N = 0     1      2      3
*        for X:    3.e-9, 4.2e-9, 2.3e-9, 6.2e-9
* jump at f=4:
*         for XDF: 6.e-7, 5.4e-7, 9.6e-8, 3.1e-7
*       for XDFF: 4.7e-5, 4.8e-5, 2.3e-6, 1.5e-6
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension A(0:5,0:3),B(0:6,0:3),C(0:6,0:3),D(0:6,0:3),
     *  LA(0:3),LB(0:3),LD(0:3)
      data A/-1.570044577033d4,1.001958278442d4,-2.805343454951d3,
     *            4.121170498099d2,-3.174780572961d1,1.d0, ! X_{-1/2}
     *    1.999266880833d4,5.702479099336d3,6.610132843877d2,
     *            3.818838129486d1,1.d0,0., ! X_{1/2}
     *    1.715627994191d2,1.125926232897d2,2.056296753055d1,1.d0,0.,0.,
     *    2.138969250409d2,3.539903493971d1,1.d0,0.,0.,0./, ! X_{5/2}
     *  B/-2.782831558471d4,2.886114034012d4,-1.274243093149d4,
     *            3.063252215963d3,-4.225615045074d2,3.168918168284d1,
     *            -1.008561571363d0, ! X_{-1/2}
     *    1.771804140488d4,-2.014785161019d3,9.130355392717d1,
     *            -1.670718177489d0,0.,0.,0., ! X_{1/2}
     *    2.280653583157d2,1.193456203021d2,1.16774311354d1,
     *            -3.226808804038d-1,3.519268762788d-3,0.,0., ! X_{3/2}
     *    7.10854551271d2,9.873746988121d1,1.067755522895d0,
     *            -1.182798726503d-2,0.,0.,0./, ! X_{5/2}
     *  C/2.206779160034d-8,-1.437701234283d-6,6.103116850636d-5,
     *     -1.169411057416d-3,1.814141021608d-2,-9.588603457639d-2,1.d0,
     *  -1.277060388085d-2,7.187946804945d-2,-4.262314235106d-1,
     *      4.997559426872d-1,-1.285579118012d0,-3.930805454272d-1,1.d0,
     *  -6.321828169799d-3,-2.183147266896d-2,-1.05756279932d-1,
     *       -4.657944387545d-1,-5.951932864088d-1,3.6844711771d-1,1.d0,
     *  -3.312041011227d-2,1.315763372315d-1,-4.820942898296d-1,
     *       5.099038074944d-1,5.49561349863d-1,-1.498867562255d0,1.d0/,
     *  D/8.827116613576d-8,-5.750804196059d-6,2.429627688357d-4,
     *       -4.601959491394d-3,6.932122275919d-2,-3.217372489776d-1,
     *       3.124344749296d0, ! X_{-1/2}
     *    -9.745794806288d-3,5.485432756838d-2,-3.29946624326d-1,
     *       4.077841975923d-1,-1.145531476975d0,-6.067091689181d-2,0.,
     *    -4.381942605018d-3,-1.5132365041d-2,-7.850001283886d-2,
     *     -3.407561772612d-1,-5.074812565486d-1,-1.387107009074d-1,0.,
     *    -2.315515517515d-2,9.198776585252d-2,-3.835879295548d-1,
     *       5.415026856351d-1,-3.847241692193d-1,3.739781456585d-2,
     *       -3.008504449098d-2/, ! X_{5/2}
     *  LA/5,4,3,2/,LB/6,3,4,3/,LD/6,5,5,6/
      if (N.lt.0.or.N.gt.3) stop'FERINV7: Invalid subscript'
      if (F.le.0.) stop'FERINV7: Non-positive argument'
      if (F.lt.4.) then
         T=F
         UP=0.
         UP1=0.
         UP2=0.
         DOWN=0.
         DOWN1=0.
         DOWN2=0.
         do I=LA(N),0,-1
            UP=UP*T+A(I,N)
           if (I.ge.1) UP1=UP1*T+A(I,N)*I
           if (I.ge.2) UP2=UP2*T+A(I,N)*I*(I-1)
         enddo
         do I=LB(N),0,-1
            DOWN=DOWN*T+B(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+B(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+B(I,N)*I*(I-1)
         enddo
         X=dlog(T*UP/DOWN)
         XDF=1.d0/T+UP1/UP-DOWN1/DOWN
         XDFF=-1.d0/T**2+UP2/UP-(UP1/UP)**2-DOWN2/DOWN+(DOWN1/DOWN)**2
      else
         P=-1./(.5+N) ! = -1/(1+\nu) = power index
         T=F**P ! t - argument of the rational fraction
         T1=P*T/F ! dt/df
         T2=P*(P-1.)*T/F**2 ! d^2 t / df^2
         UP=0.
         UP1=0.
         UP2=0.
         DOWN=0.
         DOWN1=0.
         DOWN2=0.
         do I=6,0,-1
            UP=UP*T+C(I,N)
           if (I.ge.1) UP1=UP1*T+C(I,N)*I
           if (I.ge.2) UP2=UP2*T+C(I,N)*I*(I-1)
         enddo
         do I=LD(N),0,-1
            DOWN=DOWN*T+D(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+D(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+D(I,N)*I*(I-1)
         enddo
         R=UP/DOWN
         R1=(UP1-UP*DOWN1/DOWN)/DOWN ! dR/dt
         R2=(UP2-(2.*UP1*DOWN1+UP*DOWN2)/DOWN+2.*UP*(DOWN1/DOWN)**2)/
     /     DOWN
         X=R/T
         RT=(R1-R/T)/T
         XDF=T1*RT
         XDFF=T2*RT+T1**2*(R2-2.*RT)/T
      endif
      return
      end

      subroutine BLIN9(TEMP,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
*                                                       Version 21.01.10
* Stems from BLIN8 v.24.12.08
* Difference - smooth matching of different CHI ranges
* Input: TEMP=T/mc^2; CHI=(\mu-mc^2)/T
* Output: Wk - Fermi-Dirac integral of the order k+1/2
*         WkDX=dWk/dCHI, WkDT = dWk/dT, WkDXX=d^2 Wk / d CHI^2,
*         WkDTT=d^2 Wk / d T^2, WkDXT=d^2 Wk /dCHIdT,
*         W0XXX=d^3 W0 / d CHI^3, W0XTT=d^3 W0 /(d CHI d^2 T),
*         W0XXT=d^3 W0 /dCHI^2 dT
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (CHI1=0.6d0,CHI2=14.d0,XMAX=30.d0)
      parameter (DCHI1=.1d0,DCHI2=CHI2-CHI1-DCHI1)
      parameter (XSCAL1=XMAX/DCHI1,XSCAL2=XMAX/DCHI2)
      X1=(CHI-CHI1)*XSCAL1
      X2=(CHI-CHI2)*XSCAL2
      if (X1.lt.-XMAX) then
         call BLIN9a(TEMP,CHI,
     *     W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *     W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *     W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *     W0XXX,W0XTT,W0XXT)
      elseif (X2.lt.XMAX) then ! match two fits
        if (X1.lt.XMAX) then ! match fits "a" and "b"
           call FERMI10(X1,XMAX,FP,FM)
           call BLIN9a(TEMP,CHI,
     *       W0a,W0DXa,W0DTa,W0DXXa,W0DTTa,W0DXTa,
     *       W1a,W1DXa,W1DTa,W1DXXa,W1DTTa,W1DXTa,
     *       W2a,W2DXa,W2DTa,W2DXXa,W2DTTa,W2DXTa,
     *       W0XXXa,W0XTTa,W0XXTa)
           call BLIN9b(TEMP,CHI,
     *       W0b,W0DXb,W0DTb,W0DXXb,W0DTTb,W0DXTb,
     *       W1b,W1DXb,W1DTb,W1DXXb,W1DTTb,W1DXTb,
     *       W2b,W2DXb,W2DTb,W2DXXb,W2DTTb,W2DXTb,
     *       W0XXXb,W0XTTb,W0XXTb)
        else ! match fits "b" and "c"
           call FERMI10(X2,XMAX,FP,FM)
           call BLIN9b(TEMP,CHI,
     *       W0a,W0DXa,W0DTa,W0DXXa,W0DTTa,W0DXTa,
     *       W1a,W1DXa,W1DTa,W1DXXa,W1DTTa,W1DXTa,
     *       W2a,W2DXa,W2DTa,W2DXXa,W2DTTa,W2DXTa,
     *       W0XXXa,W0XTTa,W0XXTa)
           call BLIN9c(TEMP,CHI,
     *       W0b,W0DXb,W0DTb,W0DXXb,W0DTTb,W0DXTb,
     *       W1b,W1DXb,W1DTb,W1DXXb,W1DTTb,W1DXTb,
     *       W2b,W2DXb,W2DTb,W2DXXb,W2DTTb,W2DXTb,
     *       W0XXXb,W0XTTb,W0XXTb)
        endif
         W0=W0a*FP+W0b*FM
         W0DX=W0DXa*FP+W0DXb*FM !! +(W0a-W0b)*F1
         W0DT=W0DTa*FP+W0DTb*FM
         W0DXX=W0DXXa*FP+W0DXXb*FM !! +2.d0*(W0DXa-W0DXb)*F1+(W0a-W0b)*F2
         W0DTT=W0DTTa*FP+W0DTTb*FM
         W0DXT=W0DXTa*FP+W0DXTb*FM !! +(W0DTa-W0DTb)*F1
         W0XXX=W0XXXa*FP+W0XXXb*FM !! +3.d0*(W0DXXa-W0DXXb)*F1+3.d0*(W0DXa-W0DXb)*F2+(W0a-W0b)*F3
         W0XTT=W0XTTa*FP+W0XTTb*FM !! +(W0DTTa-W0DTTb)*F1
         W0XXT=W0XXTa*FP+W0XXTb*FM !! +2.d0*(W0DXTa-W0DXTb)*F1+(W0DTa-W0DTb)*F2
         W1=W1a*FP+W1b*FM
         W1DX=W1DXa*FP+W1DXb*FM !! +(W1a-W1b)*F1
         W1DT=W1DTa*FP+W1DTb*FM
         W1DXX=W1DXXa*FP+W1DXXb*FM !! +2.d0*(W1DXa-W1DXb)*F1+(W1a-W1b)*F2
         W1DTT=W1DTTa*FP+W1DTTb*FM
         W1DXT=W1DXTa*FP+W1DXTb*FM !! +(W1DTa-W1DTb)*F1
         W2=W2a*FP+W2b*FM
         W2DX=W2DXa*FP+W2DXb*FM !! +(W2a-W2b)*F1
         W2DT=W2DTa*FP+W2DTb*FM
         W2DXX=W2DXXa*FP+W2DXXb*FM !! +2.d0*(W2DXa-W2DXb)*F1+(W2a-W2b)*F2
         W2DTT=W2DTTa*FP+W2DTTb*FM
         W2DXT=W2DXTa*FP+W2DXTb*FM !!
      else
         call BLIN9c(TEMP,CHI,
     *     W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *     W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *     W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *     W0XXX,W0XTT,W0XXT)
      endif
      return
      end

      subroutine BLIN9a(TEMP,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
*                                                       Version 19.01.10
* First part of BILN9: small CHI. Stems from BLIN9 v.24.12.08
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension AC(5,0:2),AU(5,0:2),AA(5,0:2)
      data AC/.37045057 d0, .41258437 d0,
     &        9.777982 d-2, 5.3734153 d-3, 3.8746281 d-5, ! c_i^0
     &        .39603109 d0, .69468795 d0,
     &        .22322760 d0, 1.5262934 d-2, 1.3081939 d-4, ! c_i^1
     &        .76934619 d0, 1.7891437 d0,
     &        .70754974 d0, 5.6755672 d-2, 5.5571480 d-4/ ! c_i^2
      data AU/.43139881 d0, 1.7597537 d0,
     &        4.1044654 d0, 7.7467038 d0, 13.457678 d0, ! \chi_i^0
     &        .81763176 d0, 2.4723339 d0,
     &        5.1160061 d0, 9.0441465 d0, 15.049882 d0, ! \chi_i^1
     &        1.2558461 d0, 3.2070406 d0,
     &        6.1239082 d0, 10.316126 d0, 16.597079 d0/ ! \chi_i^2
      data KRUN/0/
      KRUN=KRUN+1
      if (KRUN.eq.1) then ! initialize
        do J=0,2
        do I=1,5
           AA(I,J)=dexp(-AU(I,J))
        enddo
        enddo
      endif
        do K=0,2
           W=0.
           WDX=0.
           WDT=0.
           WDXX=0.
           WDTT=0.
           WDXT=0.
           WDXXX=0.
           WDXTT=0.
           WDXXT=0.
           ECHI=dexp(-CHI)
          do I=1,5
             SQ=dsqrt(1.d0+AU(I,K)*TEMP/2.)
             DN=AA(I,K)+ECHI ! e^{-\chi_i}+e^{-\chi})
             W=W+AC(I,K)*SQ/DN
             WDX=WDX+AC(I,K)*SQ/DN**2
             WDT=WDT+AC(I,K)*AU(I,K)/(SQ*DN)
             WDXX=WDXX+AC(I,K)*SQ*(ECHI-AA(I,K))/DN**3
             WDTT=WDTT-AC(I,K)*AU(I,K)**2/(DN*SQ**3)
             WDXT=WDXT+AC(I,K)*AU(I,K)/(SQ*DN**2)
             WDXXX=WDXXX+AC(I,K)*SQ*
     *         (ECHI**2-4.*ECHI*AA(I,K)+AA(I,K)**2)/DN**4
             WDXTT=WDXTT-AC(I,K)*AU(I,K)**2/(DN**2*SQ**3)
             WDXXT=WDXXT+AC(I,K)*AU(I,K)*(ECHI-AA(I,K))/(SQ*DN**3)
          enddo
           WDX=WDX*ECHI
           WDT=WDT/4.
           WDXX=WDXX*ECHI
           WDTT=WDTT/16.
           WDXT=WDXT/4.*ECHI
           WDXXX=WDXXX*ECHI
           WDXTT=WDXTT*ECHI/16.
           WDXXT=WDXXT/4.*ECHI
          if (K.eq.0) then
             W0=W
             W0DX=WDX
             W0DT=WDT
             W0DXX=WDXX
             W0DTT=WDTT
             W0DXT=WDXT
             W0XXX=WDXXX
             W0XTT=WDXTT
             W0XXT=WDXXT
          elseif (K.eq.1) then
             W1=W
             W1DX=WDX
             W1DT=WDT
             W1DXX=WDXX
             W1DTT=WDTT
             W1DXT=WDXT
          else
             W2=W
             W2DX=WDX
             W2DT=WDT
             W2DXX=WDXX
             W2DTT=WDTT
             W2DXT=WDXT
          endif
        enddo ! next K
      return
      end

      subroutine BLIN9b(TEMP,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
*                                                       Version 19.01.10
*                                              Small syntax fix 15.03.13
* Second part of BILN9: intermediate CHI. Stems from BLIN8 v.24.12.08
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension AX(5),AXI(5),AH(5),AV(5)
      parameter (EPS=1.d-3)
      data AX/7.265351 d-2, .2694608 d0,
     &        .533122 d0, .7868801 d0, .9569313 d0/ ! x_i
      data AXI/.26356032 d0, 1.4134031 d0,
     &         3.5964258 d0, 7.0858100 d0, 12.640801 d0/ ! \xi_i
      data AH/3.818735 d-2, .1256732 d0,
     &        .1986308 d0, .1976334 d0, .1065420 d0/ ! H_i
      data AV/.29505869 d0, .32064856 d0, 7.3915570 d-2,
     &        3.6087389 d-3, 2.3369894 d-5/ ! \bar{V}_i
      if (CHI.lt.EPS) stop'BLIN9b: CHI is too small'
        do K=0,2
           W=0.
           WDX=0.
           WDT=0.
           WDXX=0.
           WDTT=0.
           WDXT=0.
           WDXXX=0.
           WDXTT=0.
           WDXXT=0.
             SQCHI=dsqrt(CHI)
            do I=1,5
               CE=AX(I)-1.d0
               ECHI=dexp(CE*CHI)
               DE=1.d0+ECHI
               D=1.d0+AX(I)*CHI*TEMP/2.
               H=CHI**(K+1)*SQCHI*dsqrt(D)/DE
               HX=(K+1.5)/CHI+.25*AX(I)*TEMP/D-ECHI*CE/DE
               HDX=H*HX
               HXX=(K+1.5)/CHI**2+.125*(AX(I)*TEMP/D)**2+ECHI*(CE/DE)**2
               HDXX=HDX*HX-H*HXX
               HT=.25*AX(I)*CHI/D
               HDT=H*HT
               HDTT=-H*HT**2
               HTX=1./CHI-.5*AX(I)*TEMP/D
               HDXT=HDX*HT+HDT*HTX
               HDXXT=HDXX*HT+HDX*HT*HTX+HDXT*HTX+
     +           HDT*(.25*(AX(I)*TEMP/D)**2-1./CHI**2)
               HDXTT=HDXT*HT-HDX*.125*(AX(I)*CHI/D)**2+HDTT*HTX+
     +           HDT*.5*AX(I)*(TEMP*.5*AX(I)*CHI/D**2-1./D)
               HXXX=(2*K+3)/CHI**3+.125*(AX(I)*TEMP/D)**3-
     -           ECHI*(1.d0-ECHI)*(CE/DE)**3
               HDXXX=HDXX*HX-2.*HDX*HXX+H*HXXX
               XICHI=AXI(I)+CHI
               DXI=1.d0+XICHI*TEMP/2.
               V=XICHI**K*dsqrt(XICHI*DXI)
               VX=(K+.5)/XICHI+.25*TEMP/DXI
               VDX=V*VX
               VT=.25*XICHI/DXI
               VDT=V*VT
               VXX=(K+.5)/XICHI**2+.125*(TEMP/DXI)**2
               VDXX=VDX*VX-V*VXX
               VDXXX=VDXX*VX-2.*VDX*VXX+
     +           V*((2*K+1)/XICHI**3+.125*(TEMP/DXI)**3)
               VXXT=(1.-.5*TEMP*XICHI/DXI)/DXI
               VDTT=-V*VT**2
               VXT=1./XICHI-.5*TEMP/DXI
               VDXT=VDT*VXT+VDX*VT
               VDXXT=VDXT*VX+VDX*.25*VXXT-VDT*VXX-V*.25*TEMP/DXI*VXXT
               VDXTT=VDTT*VXT-VDT*.5*VXXT+VDXT*VT-
     -           VDX*.125*(XICHI/DXI)**2
               W=W+AH(I)*AX(I)**K*H+AV(I)*V
               WDX=WDX+AH(I)*AX(I)**K*HDX+AV(I)*VDX
               WDT=WDT+AH(I)*AX(I)**K*HDT+AV(I)*VDT
               WDXX=WDXX+AH(I)*AX(I)**K*HDXX+AV(I)*VDXX
               WDTT=WDTT+AH(I)*AX(I)**K*HDTT+AV(I)*VDTT
               WDXT=WDXT+AH(I)*AX(I)**K*HDXT+AV(I)*VDXT
               WDXXX=WDXXX+AH(I)*AX(I)**K*HDXXX+AV(I)*VDXXX
               WDXTT=WDXTT+AH(I)*AX(I)**K*HDXTT+AV(I)*VDXTT
               WDXXT=WDXXT+AH(I)*AX(I)**K*HDXXT+AV(I)*VDXXT
            enddo
          if (K.eq.0) then
             W0=W
             W0DX=WDX
             W0DT=WDT
             W0DXX=WDXX
             W0DTT=WDTT
             W0DXT=WDXT
             W0XXX=WDXXX
             W0XTT=WDXTT
             W0XXT=WDXXT
          elseif (K.eq.1) then
             W1=W
             W1DX=WDX
             W1DT=WDT
             W1DXX=WDXX
             W1DTT=WDTT
             W1DXT=WDXT
          else
             W2=W
             W2DX=WDX
             W2DT=WDT
             W2DXX=WDXX
             W2DTT=WDTT
             W2DXT=WDXT
          endif
        enddo ! next K
      return
      end

      subroutine BLIN9c(TEMP,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
*                                                       Version 19.01.10
* Third part of BILN9: large CHI. Stems from BLIN8 v.24.12.08
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (PI=3.141592653d0,PI26=PI*PI/6.)
      dimension AM(0:2),AMDX(0:2),AMDT(0:2),
     *  AMDXX(0:2),AMDTT(0:2),AMDXT(0:2)
      if (CHI*TEMP.lt..1) then
        do K=0,2
           W=0.
           WDX=0.
           WDT=0.
           WDXX=0.
           WDTT=0.
           WDXT=0.
           WDXXX=0.
           WDXTT=0.
           WDXXT=0.
          do J=0,4 ! for nonrel.Fermi integrals from k+1/2 to k+4.5
             CNU=K+J+.5 ! nonrelativistic Fermi integral index \nu
             CHINU=CHI**(K+J)*dsqrt(CHI) ! \chi^\nu
             F=CHINU*(CHI/(CNU+1.)+PI26*CNU/CHI+ ! nonrel.Fermi
     +         .7*PI26**2*CNU*(CNU-1.)*(CNU-2.)/CHI**3)
             FDX=CHINU*(1.+PI26*CNU*(CNU-1.)/CHI**2+
     +         .7*PI26**2*CNU*(CNU-1.)*(CNU-2.)*(CNU-3.)/CHI**4)
             FDXX=CHINU/CHI*CNU*(1.+PI26*(CNU-1.)*(CNU-2.)/CHI**2+
     +         .7*PI26**2*(CNU-1.)*(CNU-2.)*(CNU-3.)*(CNU-4.)/CHI**4)
             FDXXX=CHINU/CHI**2*CNU*(CNU-1.)*
     *         (1.+PI26*(CNU-2.)*(CNU-3.)/CHI**2+
     +         .7*PI26**2*(CNU-2.)*(CNU-3.)*(CNU-4.)*(CNU-5.)/CHI**4)
            if (J.eq.0) then
               W=F
               WDX=FDX
               WDXX=FDXX
               WDXXX=FDXXX
            elseif (J.eq.1) then
               C=.25*TEMP
               W=W+C*F ! Fermi-Dirac, expressed through Fermi
               WDX=WDX+C*FDX
               WDXX=WDXX+C*FDXX
               WDT=F/4.
               WDXT=FDX/4.
               WDTT=0.
               WDXXX=WDXXX+C*FDXXX
               WDXXT=FDXX/4.
               WDXTT=0.
            else
               C=-C/J*(2*J-3)/4.*TEMP
               W=W+C*F
               WDX=WDX+C*FDX
               WDT=WDT+C*J/TEMP*F
               WDXX=WDXX+C*FDXX
               WDTT=WDTT+C*J*(J-1)/TEMP**2*F
               WDXT=WDXT+C*J/TEMP*FDX
               WDXXX=WDXXX+C*FDXXX
               WDXTT=WDXTT+C*J*(J-1)/TEMP**2*FDX
               WDXXT=WDXXT+C*J/TEMP*FDXX
            endif
          enddo ! next J
          if (K.eq.0) then
             W0=W
             W0DX=WDX
             W0DT=WDT
             W0DXX=WDXX
             W0DTT=WDTT
             W0DXT=WDXT
             W0XXX=WDXXX
             W0XTT=WDXTT
             W0XXT=WDXXT
          elseif (K.eq.1) then
             W1=W
             W1DX=WDX
             W1DT=WDT
             W1DXX=WDXX
             W1DTT=WDTT
             W1DXT=WDXT
          else
             W2=W
             W2DX=WDX
             W2DT=WDT
             W2DXX=WDXX
             W2DTT=WDTT
             W2DXT=WDXT
          endif
        enddo ! next K
*   ----------------------------------------------------------------   *
      else ! CHI > 14, CHI*TEMP > 0.1: general high-\chi expansion
         D=1.d0+CHI*TEMP/2.d0
         R=dsqrt(CHI*D)
         RX=.5d0/CHI+.25d0*TEMP/D
         RDX=R*RX
         RDT=.25d0*CHI**2/R
         RXX=-.5d0/CHI**2-.125d0*(TEMP/D)**2
         RDXX=RDX*RX+R*RXX
         RDTT=-.25d0*RDT*CHI/D
         RXT=.25d0/D-.125d0*CHI*TEMP/D**2
         RDXT=RDT*RX+R*RXT
         RXXX=1.d0/CHI**3+.125d0*(TEMP/D)**3
         RDXXX=RDXX*RX+2.d0*RDX*RXX+R*RXXX
         RXTT=-.25d0/D**2*CHI+.125d0*CHI**2*TEMP/D**3
         RDXTT=RDTT*RX+2.d0*RDT*RXT+R*RXTT
         RXXT=-RXT*TEMP/D
         RDXXT=RDXT*RX+RDX*RXT+RDT*RXX+R*RXXT
        do K=0,2
           DM=K+.5d0+(K+1.d0)*CHI*TEMP/2.d0
           AM(K)=CHI**K*DM/R
           FMX1=.5d0*(K+1.)*TEMP/DM
           FMX2=.25d0*TEMP/D
           FMX=(K-.5d0)/CHI+FMX1-FMX2
           AMDX(K)=AM(K)*FMX
           CKM=.5d0*(K+1.d0)/DM
           FMT1=CKM*CHI
           FMT2=.25d0*CHI/D
           FMT=FMT1-FMT2
           AMDT(K)=AM(K)*FMT
           FMXX=-(K-.5d0)/CHI**2-FMX1**2+2.d0*FMX2**2
           AMDXX(K)=AMDX(K)*FMX+AM(K)*FMXX
           FMTT=2.d0*FMT2**2-FMT1**2
           AMDTT(K)=AMDT(K)*FMT+AM(K)*FMTT
           AMDXT(K)=AMDX(K)*FMT+AM(K)*(CKM*(1.d0-CKM*CHI*TEMP)-
     -       .25d0/D+.125d0*CHI*TEMP/D**2)
          if (K.eq.0) then
             FMXXX=(2*K-1)/CHI**3+2.d0*FMX1**3-8.d0*FMX2**3
             AMDXXX=AMDXX(K)*FMX+2.d0*AMDX(K)*FMXX+AM(K)*FMXXX
             FMT1DX=CKM-TEMP*CHI*CKM**2
             FMT2DX=(.25d0-CHI*TEMP*.125d0/D)/D
             FMXT=FMT1DX-FMT2DX
             FMTTX=4.d0*FMT2*FMT2DX-2.d0*FMT1*FMT1DX
             AMDXTT=AMDXT(K)*FMT+AMDT(K)*FMXT+AMDX(K)*FMTT+AM(K)*FMTTX
             FMX1DT=CKM-CHI*TEMP*CKM**2
             FMX2DT=.25d0/D*(1.d0-.5d0*CHI*TEMP/D)
             FMXXT=4.d0*FMX2*FMX2DT-2.d0*FMX1*FMX1DT
             AMDXXT=AMDXT(K)*FMX+AMDX(K)*FMXT+AMDT(K)*FMXX+AM(K)*FMXXT
          endif
        enddo
         SQ2T=dsqrt(2.d0*TEMP)
           A=1.d0+CHI*TEMP+SQ2T*R
           ADX=TEMP+SQ2T*RDX
           ADT=CHI+R/SQ2T+SQ2T*RDT
           ADXX=SQ2T*RDXX
           ADTT=-R/SQ2T**3+2.d0/SQ2T*RDT+SQ2T*RDTT
           ADXT=1.d0+RDX/SQ2T+SQ2T*RDXT
           ADXTT=-RDX/SQ2T**3+2.d0/SQ2T*RDXT+SQ2T*RDXTT
           ADXXT=RDXX/SQ2T+SQ2T*RDXXT
           XT1=CHI+1.d0/TEMP
           Aln=dlog(A)
           FJ0=.5d0*XT1*R-Aln/SQ2T**3
           ASQ3=A*SQ2T**3
           ASQ3DX=ADX*SQ2T**3
           FJ0DX=.5d0*(R+XT1*RDX)-ADX/ASQ3
           FJ0DT=.5d0*(XT1*RDT-R/TEMP**2)-ADT/ASQ3+
     +       .75d0/(TEMP**2*SQ2T)*Aln
           FJ0DXX=RDX+.5d0*XT1*RDXX+(ADX/A)**2/SQ2T**3-ADXX/ASQ3
           FJ0DTT=R/TEMP**3-RDT/TEMP**2+.5d0*XT1*RDTT+
     +       3.d0/(ASQ3*TEMP)*ADT+
     +     (ADT/A)**2/SQ2T**3-ADTT/ASQ3-1.875d0/(TEMP**3*SQ2T)*Aln
           BXT=1.5d0/TEMP*ADX+ADX*ADT/A-ADXT
           BXXT=1.5d0/TEMP*ADXX+(ADXX*ADT+ADX*ADXT)/A-
     -       (ADX/A)**2*ADT-ADXXT
           FJ0DXT=.5d0*(RDT-RDX/TEMP**2+XT1*RDXT)+BXT/ASQ3
           FJ0XXX=RDXX*1.5d0+.5d0*XT1*RDXXX+
     +      (2.d0*ADX*(ADXX/A-(ADX/A)**2)-
     -      SQ2T*RDXXX+ADXX/ASQ3*ASQ3DX)/ASQ3
           FJ0XTT=RDX/TEMP**3-RDXT/TEMP**2+.5d0*(RDTT+XT1*RDXTT)+
     +      3.d0/TEMP*(ADXT-ADT/ASQ3*ASQ3DX)/ASQ3+
     +      (2.d0*ADT*(ADXT/A-ADT*ADX/A**2)-
     -      ADXTT+ADTT*ASQ3DX/ASQ3)/ASQ3-1.875d0/(TEMP**3*SQ2T)*ADX/A
           FJ0XXT=.5d0*(RDXT-RDXX/TEMP**2+RDXT+XT1*RDXXT)+
     +      (BXXT-BXT*ASQ3DX/ASQ3)/ASQ3
         W0=FJ0+PI26*AM(0)
         W0DX=FJ0DX+PI26*AMDX(0)
         W0DT=FJ0DT+PI26*AMDT(0)
         W0DXX=FJ0DXX+PI26*AMDXX(0)
         W0DTT=FJ0DTT+PI26*AMDTT(0)
         W0DXT=FJ0DXT+PI26*AMDXT(0)
         W0XXX=FJ0XXX+PI26*AMDXXX
         W0XTT=FJ0XTT+PI26*AMDXTT
         W0XXT=FJ0XXT+PI26*AMDXXT
           FJ1=(R**3/1.5d0-FJ0)/TEMP
           FJ1DX=(2.d0*R**2*RDX-FJ0DX)/TEMP
           FJ1DT=(2.d0*R**2*RDT-FJ0DT-FJ1)/TEMP
           FJ1DXX=(4.d0*R*RDX**2+2.d0*R**2*RDXX-FJ0DXX)/TEMP
           FJ1DTT=(4.d0*R*RDT**2+2.d0*R**2*RDTT-FJ0DTT-2.d0*FJ1DT)/TEMP
           FJ1DXT=(4.d0*R*RDX*RDT+2.d0*R**2*RDXT-FJ0DXT-FJ1DX)/TEMP
         W1=FJ1+PI26*AM(1)
         W1DX=FJ1DX+PI26*AMDX(1)
         W1DT=FJ1DT+PI26*AMDT(1)
         W1DXX=FJ1DXX+PI26*AMDXX(1)
         W1DTT=FJ1DTT+PI26*AMDTT(1)
         W1DXT=FJ1DXT+PI26*AMDXT(1)
           FJ2=(.5d0*CHI*R**3-1.25d0*FJ1)/TEMP
           FJ2DX=(.5d0*R**3+1.5d0*CHI*R**2*RDX-1.25d0*FJ1DX)/TEMP
           FJ2DT=(1.5d0*CHI*R**2*RDT-1.25d0*FJ1DT-FJ2)/TEMP
           FJ2DXX=(3.d0*R*RDX*(R+CHI*RDX)+1.5d0*CHI*R**2*RDXX-
     -       1.25d0*FJ1DXX)/TEMP
          FJ2DTT=(3.d0*CHI*R*(RDT**2+.5d0*R*RDTT)-
     -      1.25d0*FJ1DTT-2.d0*FJ2DT)/TEMP
           FJ2DXT=(1.5d0*R*RDT*(R+2.d0*CHI*RDX)+1.5d0*CHI*R**2*RDXT-
     -       1.25d0*FJ1DXT-FJ2DX)/TEMP
         W2=FJ2+PI26*AM(2)
         W2DX=FJ2DX+PI26*AMDX(2)
         W2DT=FJ2DT+PI26*AMDT(2)
         W2DXX=FJ2DXX+PI26*AMDXX(2)
         W2DTT=FJ2DTT+PI26*AMDTT(2)
         W2DXT=FJ2DXT+PI26*AMDXT(2)
      endif
      return
      end

      subroutine CHEMFIT7(DENR,TEMR,CHI,CMU1,KDERIV,
     *  CMUDENR,CMUDT,CMUDTT)
*                                                       Version 29.08.15
* Fit to the chemical potential of free electron gas described in:
*     G.Chabrier & A.Y.Potekhin, Phys.Rev.E 58, 4941 (1998)
* Stems from CHEMFIT v.10.10.96. The main difference - derivatives.
*  All quantities are by default in relativistic units
* Input:  DENR - electron density, TEMR - temperature
*         KDERIV=0 if the derivatives are not required
* Output: CHI=CMU1/TEMR, where CMU1 = \mu-1 - chem.pot.w/o rest-energy
*         CMUDENR = (d\mu/d n_e)_T
*         CMUDT = (d\mu/dT)_V
*         CMUDTT = (d^2\mu/dT^2)_V
* CMUDENR,CMUDT, and CMUDTT =0 on output, if KREDIV=0
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (C13=1.d0/3.d0,PARA=1.612d0,PARB=6.192d0,PARC=.0944d0,
     *  PARF=5.535d0,PARG=.698d0)
      parameter(XEPST=228.d0) ! the largest argument of e^{-X}
      PF0=(29.6088132d0*DENR)**C13 ! Classical Fermi momentum
      if (PF0.gt.1.d-4) then
         TF=dsqrt(1.d0+PF0**2)-1.d0 ! Fermi temperature
      else
         TF=.5d0*PF0**2
      endif
      THETA=TEMR/TF
      THETA32=THETA*dsqrt(THETA)
      Q2=12.d0+8.d0/THETA32
      T1=0.
      if (THETA.lt.XEPST) T1=dexp(-THETA)
      U3=T1**2+PARA
      THETAC=THETA**PARC
      THETAG=THETA**PARG
      D3=PARB*THETAC*T1**2+PARF*THETAG
      Q3=1.365568127d0-U3/D3 ! 1.365...=2/\pi^{1/3}
      if (THETA.gt.1.d-5) then
         Q1=1.5d0*T1/(1.d0-T1)
      else
         Q1=1.5d0/THETA
      endif
      SQT=dsqrt(TEMR)
      G=(1.d0+Q2*TEMR*Q3+Q1*SQT)*TEMR
      H=(1.d0+.5d0*TEMR/THETA)*(1.d0+Q2*TEMR)
      CT=1.d0+G/H
      F=2.d0*C13/THETA32
      call FERINV7(F,1,X,XDF,XDFF)
      CHI=X ! non-relativistic result
     -  -    1.5d0*dlog(CT) ! Relativistic fit
      CMU1=TEMR*CHI ! Fit to chemical potential w/o mc^2
      if (KDERIV.eq.0) then ! DISMISS DERIVATIVES
         CMUDENR=0.
         CMUDT=0.
         CMUDTT=0.
         return
      endif
* CALCULATE DERIVATIVES:
* 1: derivatives of CHI over THETA and T
* (a): Non-relativistic result:
      THETA52=THETA32*THETA
      CHIDY=-XDF/THETA52 ! d\chi/d\theta
      CHIDYY=(XDFF/THETA**4-2.5d0*CHIDY)/THETA ! d^2\chi/d\theta^2
* (b): Relativistic corrections:
      if (THETA.gt.1.d-5) then
         Q1D=-Q1/(1.d0-T1)
         Q1DD=-Q1D*(1.d0+T1)/(1.d0-T1)
      else
         Q1D=-1.5d0/THETA**2
         Q1DD=-2.d0*Q1D/THETA ! sign corrected 08.08.11
      endif
      Q2D=-12.d0/THETA52 ! d q_2 / d \theta
      Q2DD=30.d0/(THETA52*THETA) ! d^2 q_2 / d \theta^2
      U3D=-2.d0*T1**2
      D3D=PARF*PARG*THETAG/THETA+PARB*T1**2*THETAC*(PARC/THETA-2.d0)
      D3DD=PARF*PARG*(PARG-1.d0)*THETAG/THETA**2+
     +PARB*T1**2*THETAC*(PARC*(PARC-1.d0)/THETA**2-4.d0*PARC/THETA+4.d0)
      Q3D=(D3D*U3/D3-U3D)/D3
      Q3DD=(2.d0*U3D+(2.d0*U3D*D3D+U3*D3DD)/D3-2.d0*U3*(D3D/D3)**2)/D3
      GDY=TEMR*(Q1D*SQT+(Q2D*Q3+Q2*Q3D)*TEMR) ! dG/d\theta
      GDT=1.d0+1.5d0*Q1*SQT+2.d0*Q2*Q3*TEMR
      GDYY=TEMR*(Q1DD*SQT+(Q2DD*Q3+2.d0*Q2D*Q3D+Q2*Q3DD)*TEMR)
      GDTT=.75d0*Q1/SQT+2.d0*Q2*Q3
      GDYT=1.5d0*Q1D*SQT+2.d0*(Q2D*Q3+Q2*Q3D)*TEMR
      HDY=(-.5d0/THETA**2+Q2D+.5d0*(Q2D-Q2/THETA)/THETA*TEMR)*TEMR
      HDT=(.5d0+Q2*TEMR)/THETA+Q2
      HDYY=TEMR/THETA**3+Q2DD*TEMR+
     +  TEMR**2*(.5d0*Q2DD-Q2D/THETA+Q2/THETA**2)/THETA
      HDTT=Q2/THETA
      HDYT=Q2D*(1.d0+TEMR/THETA)-(.5d0+Q2*TEMR)/THETA**2
      CTY=GDY/G-HDY/H
      CTT=GDT/G-HDT/H
      GH=G/H
      CTDY=GH*CTY
      CTDT=GH*CTT
      CTDYY=CTDY*CTY+GH*(GDYY/G-(GDY/G)**2-HDYY/H+(HDY/H)**2)
      CTDTT=CTDT*CTT+GH*(GDTT/G-(GDT/G)**2-HDTT/H+(HDT/H)**2)
      CTDYT=CTDT*CTY+GH*(GDYT/G-GDY*GDT/G**2-HDYT/H+HDY*HDT/H**2)
      CHIDY=CHIDY-1.5d0*CTDY/CT
      CHIDT=-1.5d0*CTDT/CT
      CHIDYY=CHIDYY+1.5d0*((CTDY/CT)**2-CTDYY/CT)
      CHIDTT=1.5d0*((CTDT/CT)**2-CTDTT/CT)
      CHIDYT=1.5d0*(CTDY*CTDT/CT**2-CTDYT/CT)
      CMUDENR=-(THETA*PF0)**2/(3.d0*DENR*(1.d0+TF))*CHIDY
      CMUDT=CHI+THETA*CHIDY+TEMR*CHIDT
      CMUDTT=2.d0*(CHIDY/TF+CHIDT+THETA*CHIDYT)+
     +  THETA/TF*CHIDYY+TEMR*CHIDTT
      return
      end

      subroutine CHEMFIT(DENS,TEMP,CHI)
*                                                       Version 07.06.07
* This is merely an interface to CHEMFIT7 for compatibility purposes.
* Input:  DENS - electron density [a.u.=6.7483346e24 cm^{-3}],
* TEMP - temperature [a.u.=2Ryd=3.15775e5 K]
* Output: CHI=\mu/TEMP, where \mu - electron chem.pot.w/o rest-energy
      implicit double precision (A-H), double precision (O-Z)
      save
      DENR=DENS/2.5733806d6 ! n_e in rel.un.=\lambda_{Compton}^{-3}
      TEMR=TEMP/1.8778865d4 ! T in rel.un.=(mc^2/k)=5.93e9 K
      call CHEMFIT7(DENR,TEMR,CHI,CMU1,0,CMUDENR,CMUDT,CMUDTT)
      return
      end

      subroutine FERMI10(X,XMAX,FP,FM)
*                                                       Version 20.01.10
* Fermi distribution function
* Input: X - argument f(x)
*        XMAX - max|X| where it is assumed that 0 < f(x) < 1.
* Output: FP = f(x)
*         FM = 1-f(x)
      implicit double precision (A-H), double precision (O-Z)
      save
      if (XMAX.lt.3.d0) stop'FERMI10: XMAX'
      if (X.gt.XMAX) then
         FP=0.d0
         FM=1.d0
      elseif (X.lt.-XMAX) then
         FP=1.d0
         FM=0.d0
      else
         FP=1.d0/(dexp(X)+1.d0)
         FM=1.d0-FP
      endif
      return
      end subroutine FERMI10
