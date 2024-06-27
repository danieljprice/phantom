module chemistry
 implicit none
       real :: wind_CO_ratio = 0.34 !Based on Ferrarotti and Gail 2002 0.34

 ! Indices for elements and molecules:
 integer, parameter :: nElements = 11
 integer, parameter :: iH=1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10, iMg=11 

 integer, parameter :: nMolecules = 69 
 integer, parameter :: iH2=1, iOH=2, iH2O=3, iCO=4, iCO2=5, iCH4=6, iC2H=7, iC2H2=8, iN2=9, iNH3=10, iCN=11, &
       iHCN=12, iSi2=13, iSi3=14, iSiO=15, iSi2C=16, iSiH4=17, iS2=18, iHS=19, iH2S=20, iSiS=21, &
       iSiH=22, iTiO=23, iTiO2=24, iC2=25, iO2=26, iCH=27,iCOH=28, iC2O=29, iCH2=30, iH2CO=31, iCH3=32, &
       iC2H4=33, iNH=34, iNO=35, iNCO=36, iHCNO=37, iC2N=38, iC2N2=39, iHNO=40, iHNO2=41, iHNO3=42, &
       iNH2=43, iNO2=44, INO3=45, iN2O=46, iN2O4=47, iMgH=48, iMgOH=49, iMgO2H2=50, iMgN=51, iMgO=52, &
       iSiC=53, iSiH2=54, iSiH3=55, iSiN=56, iSiO2=57, iFeO=58, iFeO2H2=59, iCOS=60, iCS=61, &
       iCS2=62, iFeS=63, iH2SO4=64, iMgS=65, iSN=66, iSO=67, iSO2=68, iSO3=69, iTiS=70
 
 real, parameter :: coefs(5,nMolecules) = reshape([&
       4.25321d+05, -1.07123d+05, 2.69980d+01, 5.48280d-04, -3.81498d-08, & !H2-
       4.15670d+05, -1.05260d+05, 2.54985d+01, 4.78020d-04, -2.82416d-08, & !OH-
       8.66184d+05, -2.27851d+05, 5.61473d+01, 7.62548d-04, -4.95254d-08, & !H2O
       3.30340d+05, -2.59792d+05, 3.23662d+01, 3.33339d-04, -1.69521d-08, & !CO
       5.34072d+05, -3.88536d+05, 6.95497d+01, 0.00000d+00,  0.00000d+00, & !CO2
       1.51591d+06, -4.09711d+05, 1.19952d+02, 5.71281d-04, -4.77554d-08, & !CH4
       6.31402d+05, -2.85395d+05, 5.97427d+01, 1.13748d-04, -2.04586d-08, & !C2H
       8.11613d+05, -3.99041d+05, 9.15250d+01, 2.53889d-05, -2.03949d-08, & !C2H2
       3.80144d+05, -2.28698d+05, 3.09554d+01, 2.71276d-04, -6.21348d-09, & !N2
       1.23626d+06, -2.89588d+05, 8.52686d+01, 5.66638d-04, -3.66862d-08, & !NH3
       4.51527d+05, -1.83350d+05, 2.97771d+01, 1.58267d-04, -1.72510d-08, & !CN
       6.02063d+05, -3.08724d+05, 6.00139d+01, 1.66592d-04, -1.13252d-08, & !HCN
      -1.10204d+05, -7.41179d+04, 2.57470d+01, 1.37542d-04, -2.83506d-09, & !Si2
       6.06066d+04, -1.71746d+05, 5.75688d+01, 1.18547d-04, -6.23584d-08, & !Si3
       1.82780d+05, -1.92839d+05, 3.03804d+01, 4.16079d-04, -1.91035d-08, & !SiO
       2.34924d+05, -2.60546d+05, 6.29047d+01, 2.75378d-04, -3.64630d-08, & !Si2C
       9.83348d+05, -3.16438d+05, 1.13885d+02, 2.90172d-04, -2.37569d-08, & !SiH4
       2.24963d+05, -1.03636d+05, 2.85814d+01, 1.77872d-04, -8.41752d-09, & !S2
       3.65656d+05, -8.77172d+04, 2.41557d+01, 3.40917d-04, -2.09217d-08, & !HS
       7.41911d+05, -1.81063d+05, 5.35114d+01, 4.94594d-04, -3.39231d-08, & !H2S
       1.44507d+05, -1.49916d+05, 2.87381d+01, 3.96186d-04, -2.18002d-08, & !SiS
       2.45408d+05, -7.17752d+04, 2.29449d+01, 4.52999d-04, -2.69915d-08, & !SiH
      -4.38897d+05, -1.58111d+05, 2.49224d+01, 1.08714d-03, -5.62504d-08, & !TiO
      -3.32351d+05, -3.04694d+05, 5.86984d+01, 1.17096d-03, -5.06729d-08, & !TiO2
       2.26786d+05, -1.43775d+05, 2.92429d+01, 1.69434d-04, -1.79867d-08, & !C2
       3.17263d+05, -1.21503d+05, 3.11753d+01, 2.53855d-04, -1.73743d-08, & !O2
       4.08366d+05, -8.44491d+04, 2.54152d+01, 1.07036d-04, -1.04264d-08, & !CH
       6.07224d+05, -2.77376d+05, 5.66759d+01, 4.80619d-04, -2.82627d-08, & !COH
       3.56583d+05, -3.36691d+05, 6.30120d+01, -2.15913d-04, 1.12031d-08, & !C2O
       6.99399d+05, -1.88493d+05, 5.33630d+01, 5.15801d-04, -3.28189d-08, & !CH2
       9.91212d+05, -3.70812d+05, 9.03305d+01, 4.39620d-04, -3.01383d-08, & !H2CO
       1.06713d+06, -3.01022d+05, 8.47773d+01, 6.14095d-04, -4.34934d-08, & !CH3
       1.57756d+06, -5.51401d+05, 1.51644d+02, 2.93007d-04, -3.54396d-08, & !C2H4
       4.11181d+05, -7.79746d+04, 2.42428d+01, 4.09611d-04, -2.55866d-08, & !NH
       3.18598d+05, -1.53355d+05, 2.79315d+01, 2.76970d-04, -9.14848d-09, & !NO
       3.22530d+05, -3.08692d+05, 6.20446d+01, 5.46329d-05,  0.00000d+00, & !NCO
       7.64992d+05, -4.26302d+05, 9.19995d+01, 2.86974d-04, -1.52951d-08, & !HCNO
       3.05118d+05, -3.24444d+05, 5.98669d+01, 0.00000d+00,  0.00000d+00, & !C2N
       1.77312d+05, -4.97555d+05, 9.62430d+01, -4.71719d-05, 0.00000d+00, & !C2N2
       6.24780d+05, -2.05700d+05, 5.64594d+01, 4.64490d-04, -2.07625d-08, & !HNO
       6.98937d+05, -3.08668d+05, 8.96184d+01, 2.67463d-04, -1.22202d-08, & !HNO2 (trans)
       8.02373d+05, -3.82689d+05, 1.26355d+02, 1.65050d-05, -2.02960d-10, & !HNO3
       8.79993d+05, -1.77983d+05, 5.28632d+01, 6.37403d-04, -4.63748d-08, & !NH2
       4.70857d+05, -2.28020d+05, 6.18243d+01, 3.24505d-04, -9.54732d-09, & !NO2
       4.42258d+05, -2.78648d+05, 9.82856d+01, 6.50496d-06,  7.35167d-09, & !NO3
       4.31224d+05, -2.69421d+05, 6.44178d+01, 1.28733d-05,  7.77099d-09, & !N2O
       5.92417d+05, -4.67573d+05, 1.63906d+02, -5.25057d-04, 3.57897d-08, & !N2O4
       2.90391d+05, -4.91449d+04, 1.96635d+01, 1.52956d-04, -6.37175d-09, & !MgH
       4.47672d+05, -1.89626d+05, 5.17801d+01, 1.15404d-04,  3.34857d-10, & !MgOH
       7.28023d+05, -4.00948d+05, 1.11212d+02, 5.21735d-05,  1.76309d-11, & !MgO2H2
       2.50267d+05, -8.09495d+04, 2.05587d+01, 4.19383d-05,  1.08947d-08, & !MgN
      -3.68825d+05, -7.91991d+04, 2.24741d+01, -3.33176d-04, 3.46191d-08, & !MgO
       8.92245d+03, -1.06824d+05, 2.65205d+01, 1.67830d-04, -7.50467d-09, & !SiC
       2.46640d+01, -6.88730d+00, 8.37100d-02, -1.00580d-02, 4.92910d-04, & !SiH2
       3.63290d+01, -1.05560d+01, 8.09450d-02, -8.62120d-03, 3.98640d-04, & !SiH3
       3.30323d+05, -1.34105d+05, 2.82397d+01, -4.05990d-05, 5.40508d-09, & !SiN
       2.61948d+05, -3.02036d+05, 6.59588d+01, 1.07405d-04,  0.00000d+00, & !SiO2
       2.18545d+05, -1.00927d+05, 2.69229d+01, 2.53772d-04,  4.37490d-09, & !FeO
       4.18570d+05, -4.05885d+05, 1.14517d+02, -2.48115d-04, 3.53182d-08, & !FeO2H2
       3.84946d+05, -3.33394d+05, 6.59981d+01, -6.67359d-05, 0.00000d+00, & !COS
       2.54293d+05, -1.72618d+05, 3.05596d+01, 2.67035d-04, -1.59172d-08, & !CS
       2.41613d+05, -2.78285d+05, 6.52867d+01, -1.07214d-04, 0.00000d+00, & !CS2
       4.56757d+04, -7.77456d+04, 2.39279d+01, 3.02718d-04,  6.29763d-09, & !FeS
       9.36433d+05, -5.92570d+05, 1.90572d+02, -2.99786d-04, 1.02593d-08, & !H2SO4
      -3.64897d+05, -6.47278d+04, 2.04150d+01, -1.93425d-04, 2.73835d-08, & !MgS
       2.71848d+05, -1.18293d+05, 2.65189d+01, 1.88122d-04, -5.20791d-09, & !SN
       3.21629d+05, -1.27035d+05, 2.89402d+01, 1.25603d-04, -9.85969d-09, & !SO
       4.16261d+05, -2.59787d+05, 6.30161d+01, 2.30250d-04, -1.21022d-08, & !SO2
       4.54902d+05, -3.43672d+05, 1.00999d+02, 5.09991d-05,  0.00000d+00], shape(coefs)) !SO3

       real, parameter :: atomic_mass_unit = 1.66e-24
       real, parameter :: kboltz = 1.3807e-16
       real, parameter :: patm = 1.013250d6
       integer, parameter :: ncols = nElements + nMolecules -2 !Not including Ne and Ti

       character(len=*), parameter :: columns(ncols) =&
                 (/'nH2              ', &
                   'nOH              ', &
                   'nH2O             ', &
                   'nCO              ', &
                   'nCO2             ', &
                   'nCH4             ', &
                   'nC2H             ', &
                   'nC2H2            ', &
                   'nN2              ', &
                   'nNH3             ', &
                   'nCN              ', &
                   'nHCN             ', &
                   'nSi2             ', &
                   'nSi3             ', &
                   'nSiO             ', &
                   'nSi2C            ', &
                   'nSiH4            ', &
                   'nS2              ', &
                   'nHS              ', &
                   'nH2S             ', &
                   'nSiS             ', &
                   'nSiH             ', &
                   'nTiO             ', &
                   'nTiO2            ', &
                   'nC2              ', &
                   'nO2              ', &
                   'nCH              ', &
                   'nCOH             ', &
                   'nC20             ', &
                   'nCH2             ', &
                   'nH2CO            ', &
                   'nCH3             ', &
                   'nC2H4            ', &
                   'nNH              ', &
                   'nNO              ', &
                   'nNCO             ', &
                   'nHCNO            ', &
                   'nC2N             ', &
                   'nC2N2            ', &
                   'nHNO             ', &
                   'nHNO2            ', &
                   'nHNO3            ', &
                   'nNH2             ', &
                   'nNO2             ', &
                   'nNO3             ', &
                   'nN20             ', &
                   'nN204            ', &
                   'nMgH             ', &
                   'nMgOH            ', &
                   'nMgO2H2          ', &
                   'nMgN             ', &
                   'nMgO             ', &
                   'nSiC             ', &
                   'nSiH2            ', &
                   'nSiH3            ', &
                   'nSiN             ', &
                   'nSiO2            ', &
                   'nFeO             ', &
                   'nFeO2H2          ', &
                   'nCOS             ', &
                   'nCS              ', &
                   'nCS2             ', &
                   'nFeS             ', &
                   'nH2SO4           ', &
                   'nMgS             ', &
                   'nSN              ', &
                   'nSO              ', &
                   'nSO2             ', &
                   'nSO3             ', &
                   'nH               ', &
                   'nHe              ', &
                   'nC               ', &
                   'nO               ', &
                   'nN               ', &
                   'nSi              ', &
                   'nS               ', &
                   'nFe              ', &
                   'nMg              '/)

 public

contains

subroutine network(T,rho_cgs,mu,gamma,abundance,pressure_cgs)
 real, intent(in)     :: rho_cgs
 real, intent(out)    :: mu,gamma
 real, intent(out)    :: abundance(ncols)
 real, intent(inout)  :: T
 real, intent(in), optional :: pressure_cgs

 real :: mass_per_H, eps(nElements)
 real :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867, 24.305]
!       real :: epsC, nH_tot 
!       real :: v1 
 real :: pH_tot, rho_cgs_tmp !, pgas !LUIS pH, pH_tot
 real :: KH2 !LUIS, pH2
 integer :: num, j 
 real, dimension(nMolecules)    :: pmol  !, nmol=0.
 real                           :: pelm(nElements)

       pelm = 0.
       pmol = 0.
       num =0

               ! all quantities in cgs
               eps(iH)  = 1.0
               eps(iHe) = 1.009d-1 !1.04d-1
               eps(iOx) = 7.211d-4 !6.87d-4
               eps(iN)  = 2.106d-4 !2.52d-4
               eps(iMg) = 3.859d-5 !3.85d-5
               eps(iNe) = 1.18d-4  !1.17d-4
               eps(iSi) = 3.561d-5 !3.58d-5
               eps(iS)  = 1.863d-5 !1.85d-5
               eps(iFe) = 3.241d-5 !3.24d-5
               eps(iTi) = 8.621d-8 !8.6d-8
               eps(iC)  = eps(iOx) * wind_CO_ratio 
               mass_per_H = atomic_mass_unit*dot_product(Aw,eps)

             !epsC = eps(iC)

      
 !      do i=46, 1000
 !      T = i*10.
      if (present(pressure_cgs) .and. T > 100 .and. T < 1.d5) then
         KH2      = calc_Kd(coefs(:,iH2), T)
         pelm(iH) = newton_method(0.,0.,1.,(1.+eps(iHe))/((1.+2.*eps(iHe))*KH2), &
                    -pressure_cgs/patm/((1.+2.*eps(iHe))*KH2),pressure_cgs/patm)
         pmol(iH2) = KH2*pelm(iH)**2
         pH_tot = pelm(iH)+2.*pmol(iH2)
         rho_cgs_tmp = pH_tot*patm*mass_per_H/(T*kboltz) 
      else
         rho_cgs_tmp = rho_cgs
         pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
      endif
       
      if (T>1.d5) then
         mu = (1.+4.*eps(iHe))/(1.+eps(iHe))
         gamma = 5./3.
         pelm(iH) = pH_tot
         pelm(iHe) = eps(iHe)*pH_tot
      elseif (T>100.) then !(T>600)
         KH2 = calc_Kd(coefs(:,iH2), T)
         pelm(iH) = newton_method(0.,0.,2.*KH2,1.,-pH_tot,pH_tot) !pgas = 1.d-4
         pmol(iH2) = KH2*pelm(iH)**2 !pH2 =   change pH2 for pmol(1)
         pelm(iHe) = eps(iHe)*pH_tot
         call init_muGamma(rho_cgs_tmp,mass_per_H,eps,T, mu, gamma)
         call chemical_equilibrium_light(rho_cgs_tmp, T, pmol, pelm, eps, mass_per_H, mu, gamma, nMolecules)
      endif               

       do j=1,ncols-9
       !if (j.le.69) then
       abundance(j) = pmol(j)*patm/(kboltz*T)  !*patm
       if (abundance(j).lt.1.0E-99) then
               abundance(j) = 0.0E+00
       endif
       enddo

       abundance(70) = pelm(iH)*patm/(kboltz*T) !*patm
       abundance(71) = pelm(iHe)*patm/(kboltz*T) !*patm
       abundance(72) = pelm(iC)*patm/(kboltz*T) !*patm
       abundance(73) = pelm(iOx)*patm/(kboltz*T) !*patm
       abundance(74) = pelm(iN)*patm/(kboltz*T) !*patm
       abundance(75) = pelm(iSi)*patm/(kboltz*T) !*patm
       abundance(76) = pelm(iS)*patm/(kboltz*T) !*patm
       abundance(77) = pelm(iFe)*patm/(kboltz*T) !*patm
       abundance(78) = pelm(iMg)*patm/(kboltz*T) !*patm


end subroutine network

!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance in the gas phase
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_light(rho_cgs,T,pmol,pelm,eps,&
                                      mass_per_H,mu,gamma,nMolecules) 
! all quantities are in cgs
 real, intent(in)    :: rho_cgs !,epsC
 integer, intent(in) :: nMolecules !LUIS epsSi
 real, intent(inout) :: T, mu, gamma, pmol(nMolecules), pelm(nElements)
 real, intent(in)    :: eps(nElements), mass_per_H
 !real, intent(out)   :: pC, pC2, pC2H, pC2H2
 !real, intent(out), optional :: nH, nH2, nHe, nCO, nH2O, nOH,nO, nN, nN2, nCN
 real    :: Kd(nMolecules+1), err(nElements), pH_tot !, a, b, c, d !LUIS err
 !real    :: pH, pCO, pO, pSi, pS, pTi, pN, pNewton !LUIS pMg, pFe
 !real    :: pC_old, pO_old, pSi_old, pS_old, pTi_old, cst !LUIS pFe_old, pMg_old
 real, dimension(nMolecules)    :: pmol_old
 real    :: pelm_old(nElements)
 integer :: i, nit, rndnmbr

 pelm_old = 0.
 pmol_old = 0.

 call calc_muGamma(rho_cgs,mass_per_H,eps,T, mu, gamma, pelm(iH), pH_tot) !pH,pH_tot


 if (T > 1.d4) then
    pelm = eps*pH_tot
    pmol = 0.
        return
 endif

! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 Kd(iTiS) = calc_Kd_TiS(T)

 err      = 1.
 nit      = 0


 pmol(iCO)  = eps(iC)*pH_tot
 pmol(iSiO) = eps(iSi)*pH_tot

 rndnmbr = 0

 do while (maxval(err) > 1.e-6)

 pelm_old(:) = pelm(:)

 if (rndnmbr == 0) then
         pelm(iOx) = newton_method(4.*(pelm(iN)**2*Kd(iN2O4)+pelm(iH)**2*pelm(iS)*Kd(iH2SO4)), &
                         3.*(pelm(iH)*pelm(iN)*Kd(iHNO3)+pelm(iN)*Kd(iNO3)+pelm(iS)*Kd(iSO3)), &
                         2.*(Kd(iO2)+pelm(iC)*Kd(iCO2)+pelm(iH)*pelm(iN)*Kd(iHNO2) &
                             +pelm(iN)*Kd(iNO2)+pelm(iH)**2*pelm(iMg)*Kd(iMgO2H2)+pelm(iSi)*Kd(iSiO2) &
                             +pelm(iFe)*pelm(iH)**2*Kd(iFeO2H2)+pelm(iS)*Kd(iSO2)), &
                         1.+pelm(iH)*Kd(iOH)+pelm(iH)**2*Kd(iH2O)+pelm(iC)*pelm(iH)*Kd(iCOH) &
                         +pelm(iC)**2*Kd(iC2O)+pelm(iC)*pelm(iH)**2*Kd(iH2CO)+pelm(iN)*Kd(iNO) &
                         +pelm(iC)*pelm(iN)*Kd(iNCO)+pelm(iC)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iH)*pelm(iN)*Kd(iHNO) &
                         +pelm(iN)**2*Kd(iN2O)+pelm(iH)*pelm(iMg)*Kd(iMgOH)+pelm(iMg)*Kd(iMgO)+pelm(iFe)*Kd(iFeO) &
                         +pelm(iC)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iSO), &
                         (eps(iC)+eps(iSi)-eps(iOx))*pH_tot,eps(iOx)*pH_tot) !pmol(iCO)  = eps(iC)*pH_tot, pmol(iSiO) = eps(iSi)*pH_tot
                 rndnmbr = rndnmbr + 1
         else
 pelm(iOx) = newton_method(4.*(pelm(iN)**2*Kd(iN2O4)+pelm(iH)**2*pelm(iS)*Kd(iH2SO4)), &
                         3.*(pelm(iH)*pelm(iN)*Kd(iHNO3)+pelm(iN)*Kd(iNO3)+pelm(iS)*Kd(iSO3)), &
                         2.*(Kd(iO2)+pelm(iC)*Kd(iCO2)+pelm(iH)*pelm(iN)*Kd(iHNO2) &
                             +pelm(iN)*Kd(iNO2)+pelm(iH)**2*pelm(iMg)*Kd(iMgO2H2)+pelm(iSi)*Kd(iSiO2) &
                             +pelm(iFe)*pelm(iH)**2*Kd(iFeO2H2)+pelm(iS)*Kd(iSO2)), &
                         1.+pelm(iC)*Kd(iCO)+pelm(iSi)*Kd(iSiO)+pelm(iH)*Kd(iOH)+pelm(iH)**2*Kd(iH2O)+pelm(iC)*pelm(iH)*Kd(iCOH) &
                         +pelm(iC)**2*Kd(iC2O)+pelm(iC)*pelm(iH)**2*Kd(iH2CO)+pelm(iN)*Kd(iNO) &
                         +pelm(iC)*pelm(iN)*Kd(iNCO)+pelm(iC)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iH)*pelm(iN)*Kd(iHNO) & 
                         +pelm(iN)**2*Kd(iN2O)+pelm(iH)*pelm(iMg)*Kd(iMgOH)+pelm(iMg)*Kd(iMgO)+pelm(iFe)*Kd(iFeO) &
                         +pelm(iC)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iSO), &
                         -eps(iOx)*pH_tot,eps(iOx)*pH_tot)
         endif
 pelm(iC) = newton_method(0., &
                          0., &
                          2.*(Kd(iC2)+pelm(iH)*Kd(iC2H)+pelm(iOx)*Kd(iC2O)+pelm(iH)**2*Kd(iC2H2)+pelm(iH)**4*Kd(iC2H4) &
                             +pelm(iN)*Kd(iC2N)+pelm(iN)**2*Kd(iC2N2)), &
                          1.+pelm(iH)*Kd(iCH)+pelm(iOx)*Kd(iCO)+pelm(iOx)*pelm(iH)*Kd(iCOH)+pelm(iOx)**2*Kd(iCO2) &
                          +pelm(iH)**2*Kd(iCH2)+pelm(iOx)*pelm(iH)**2*Kd(iH2CO)+pelm(iH)**3*Kd(iCH3)+pelm(iH)**4*Kd(iCH4) &
                          +pelm(iN)*Kd(iCN)+pelm(iN)*pelm(iH)*Kd(iHCN)+pelm(iOx)*pelm(iN)*Kd(iNCO) &
                          +pelm(iOx)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iSi)*Kd(iSiC)+pelm(iSi)**2*Kd(iSi2C) & 
                          +pelm(iOx)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iCS)+pelm(iS)**2*Kd(iCS2), &
                          -eps(iC)*pH_tot,eps(iC)*pH_tot)

                  pelm(iN) = newton_method(0., 0., &
                                           2.*(Kd(iN2)+pelm(iC)**2*Kd(iC2N2) & 
                                           +pelm(iOx)*Kd(iN2O) &
                                           +pelm(iOx)**4*Kd(iN2O4)), &
                                           1.+pelm(iH)*Kd(iNH)+pelm(iOx)*Kd(iNO)+pelm(iC)*Kd(iCN)+pelm(iC)*pelm(iH)*Kd(iHCN) &
                                           +pelm(iOx)*pelm(iC)*Kd(iNCO)+pelm(iC)*pelm(iOx)*pelm(iH)*Kd(iHCNO)+pelm(iC)**2*Kd(iC2N) &
                                           +pelm(iOx)*pelm(iH)*Kd(iHNO)+pelm(iOx)**2*pelm(iH)*Kd(iHNO2) &
                                           +pelm(iOx)**3*pelm(iH)*Kd(iHNO3)+pelm(iH)**2*Kd(iNH2)+pelm(iH)**3*Kd(iNH3) &
                                           +pelm(iOx)**2*Kd(iNO2)+pelm(iOx)**3*Kd(iNO3)+pelm(iMg)*Kd(iMgN) &
                                           +pelm(iSi)*Kd(iSiN) &
                                           +pelm(iS)*Kd(iSN), &
                                           -eps(iN)*pH_tot,eps(iN)*pH_tot)

                 pelm(iMg) = & 
                 eps(iMg)*pH_tot/(1.+pelm(iH)*Kd(iMgH)+pelm(iOx)*pelm(iH)*Kd(iMgOH) &
                 +pelm(iOx)**2*pelm(iH)**2*Kd(iMgO2H2)+pelm(iN)*Kd(iMgN) &
                 +pelm(iOx)*Kd(iMgO)+pelm(iS)*Kd(iMgS))

      
                pelm(iSi) = newton_method(0., 0., &
                                           2.*(Kd(iSi2)+pelm(iC)*Kd(iSi2C)), &
                                           1.+pelm(iC)*Kd(iSiC)+pelm(iH)*Kd(iSiH)+pelm(iH)**2*Kd(iSiH2) &
                                           +pelm(iH)**3*Kd(iSiH3)+pelm(iH)**4*Kd(iSiH4)+pelm(iN)*Kd(iSiN)+pelm(iOx)*Kd(iSiO) &
                                           +pelm(iOx)**2*Kd(iSiO2) &
                                           +pelm(iS)*Kd(iSiS), &
                                           -eps(iSi)*pH_tot,eps(iSi)*pH_tot)


                pelm(iFe) = &
                eps(iFe)*pH_tot/(1.+pelm(iOx)*Kd(iFeO)+pelm(iOx)**2*pelm(iH)**2*Kd(iFeO2H2) &
                                 +pelm(iS)*Kd(iFeS))

           pelm(iS) = newton_method(0., 0., 2.*(Kd(iS2)+pelm(iC)*Kd(iCS2)), &
                             1.+pelm(iOx)*pelm(iC)*Kd(iCOS)+pelm(iC)*Kd(iCS)+pelm(iFe)*Kd(iFeS)+pelm(iH)*Kd(iHS) &
                             +pelm(iH)**2*Kd(iH2S)+pelm(iOx)**4*pelm(iH)**2*Kd(iH2SO4)+pelm(iMg)*Kd(iMgS)+pelm(iN)*Kd(iSN) &
                             +pelm(iOx)*Kd(iSO)+pelm(iOx)**2*Kd(iSO2)+pelm(iOx)**3*Kd(iSO3)+pelm(iSi)*Kd(iSiS), &
                             -eps(iS)*pH_tot,eps(iS)*pH_tot)

   err = abs((pelm-pelm_old)/(pelm_old+1.0e-30))

  pmol(iH2)  = Kd(iH2)*pelm(iH)**2
  pmol(iOH)  = Kd(iOH)*pelm(iOx)*pelm(iH)
  pmol(iH2O) = Kd(iH2O)*pelm(iOx)*pelm(iH)**2
  pmol(iCO)  = Kd(iCO)*pelm(iOx)*pelm(iC)
  pmol(iCO2)  = Kd(iCO2)*pelm(iC)*pelm(iOx)**2
  pmol(iCH4)  = Kd(iCH4)*pelm(iC)*pelm(iH)**4
  pmol(iC2H)  = Kd(iC2H)*pelm(iH)*pelm(iC)**2
  pmol(iC2H2)  = Kd(iC2H2)*pelm(iC)**2*pelm(iH)**2
  pmol(iN2)  = Kd(iN2)*pelm(iN)**2
  pmol(iNH3)  = Kd(iNH3)*pelm(iN)*pelm(iH)**3
  pmol(iCN)  = Kd(iCN)*pelm(iC)*pelm(iN)
  pmol(iHCN)  = Kd(iHCN)*pelm(iH)*pelm(iC)*pelm(iN)
  pmol(iSi2)  = Kd(iSi2)*pelm(iSi)**2
  pmol(iSi3)  = Kd(iSi3)*pelm(iSi)**3
  pmol(iSiO)  = Kd(iSiO)*pelm(iSi)*pelm(iOx)
  pmol(iSi2C)  = Kd(iSi2C)*pelm(iSi)**2*pelm(iC)
  pmol(iSiH4)  = Kd(iSiH4)*pelm(iH)**4*pelm(iSi)
  pmol(iS2)  = Kd(iS2)*pelm(iS)**2
  pmol(iHS)  = Kd(iHS)*pelm(iH)*pelm(iS)
  pmol(iH2S)  = Kd(iH2S)*pelm(iH)**2*pelm(iS)
  pmol(iSiS)  = Kd(iSiS)*pelm(iSi)*pelm(iS)
  pmol(iSiH)  = Kd(iSiH)*pelm(iSi)*pelm(iH)
  pmol(iC2)  = Kd(iC2)*pelm(iC)**2
  pmol(iO2)  = Kd(iO2)*pelm(iOx)**2
  pmol(iCH)  = Kd(iCH)*pelm(iH)*pelm(iC)
  pmol(iCOH)  = Kd(iCOH)*pelm(iOx)*pelm(iC)*pelm(iH)
  pmol(iC2O)  = Kd(iC2O)*pelm(iOx)*pelm(iC)**2
  pmol(iCH2)  = Kd(iCH2)*pelm(iH)**2*pelm(iC)
  pmol(iH2CO)  = Kd(iH2CO)*pelm(iOx)*pelm(iC)*pelm(iH)**2
  pmol(iCH3)  = Kd(iCH3)*pelm(iH)**3*pelm(iC)
  pmol(iC2H4)  = Kd(iC2H4)*pelm(iH)**4*pelm(iC)**2
  pmol(iNH)  = Kd(iNH)*pelm(iN)*pelm(iH)
  pmol(iNO)  = Kd(iNO)*pelm(iOx)*pelm(iN)
  pmol(iNCO)  = Kd(iNCO)*pelm(iOx)*pelm(iC)*pelm(iN)
  pmol(iHCNO)  = Kd(iHCNO)*pelm(iOx)*pelm(iC)*pelm(iH)*pelm(iN)
  pmol(iC2N)  = Kd(iC2N)*pelm(iN)*pelm(iC)**2
  pmol(iC2N2)  = Kd(iC2N2)*pelm(iN)**2*pelm(iC)**2
  pmol(iHNO)  = Kd(iHNO)*pelm(iOx)*pelm(iH)*pelm(iN)
  pmol(iHNO2)  = Kd(iHNO2)*pelm(iOx)**2*pelm(iH)*pelm(iN)
  pmol(iHNO3)  = Kd(iHNO3)*pelm(iOx)**3*pelm(iH)*pelm(iN)
  pmol(iNH2)  = Kd(iNH2)*pelm(iN)*pelm(iH)**2
  pmol(iNO2)  = Kd(iNO2)*pelm(iN)*pelm(iOx)**2
  pmol(iNO3)  = Kd(iNO3)*pelm(iN)*pelm(iOx)**3
  pmol(iN2O)  = Kd(iN2O)*pelm(iN)**2*pelm(iOx)
  pmol(iN2O4)  = Kd(iN2O4)*pelm(iN)**2*pelm(iOx)**4
  pmol(iMgH)  = Kd(iMgH)*pelm(iH)*pelm(iMg)
  pmol(iMgOH)  = Kd(iMgOH)*pelm(iH)*pelm(iMg)*pelm(iOx)
  pmol(iMgO2H2)  = Kd(iMgO2H2)*pelm(iH)**2*pelm(iMg)*pelm(iOx)**2
  pmol(iMgN)  = Kd(iMgN)*pelm(iN)*pelm(iMg)
  pmol(iMgO)  = Kd(iMgO)*pelm(iOx)*pelm(iMg)
  pmol(iSiC)  = Kd(iSiC)*pelm(iC)*pelm(iSi)
  pmol(iSiH2)  = Kd(iSiH2)*pelm(iH)**2*pelm(iSi)
  pmol(iSiH3)  = Kd(iSiH3)*pelm(iH)**3*pelm(iSi)
  pmol(iSiN)  = Kd(iSiN)*pelm(iN)*pelm(iSi)
  pmol(iSiO2)  = Kd(iSiO2)*pelm(iOx)**2*pelm(iSi)
  pmol(iFeO)  = Kd(iFeO)*pelm(iOx)*pelm(iFe)
  pmol(iFeO2H2)  = Kd(iFeO2H2)*pelm(iH)**2*pelm(iFe)*pelm(iOx)**2
  pmol(iCOS)  = Kd(iCOS)*pelm(iOx)*pelm(iC)*pelm(iS)
  pmol(iCS)  = Kd(iCS)*pelm(iS)*pelm(iC)
  pmol(iCS2)  = Kd(iCS2)*pelm(iC)*pelm(iS)**2
  pmol(iFeS)  = Kd(iFeS)*pelm(iS)*pelm(iFe)
  pmol(iH2SO4)  = Kd(iH2SO4)*pelm(iH)**2*pelm(iS)*pelm(iOx)**4
  pmol(iMgS)  = Kd(iMgS)*pelm(iS)*pelm(iMg)
  pmol(iSN)  = Kd(iSN)*pelm(iN)*pelm(iS)
  pmol(iSO)  = Kd(iSO)*pelm(iOx)*pelm(iS)
  pmol(iSO2)  = Kd(iSO2)*pelm(iOx)**2*pelm(iS)
  pmol(iSO3)  = Kd(iSO3)*pelm(iOx)**3*pelm(iS)

    nit = nit + 1
    if (nit == 200) exit

 enddo
end subroutine chemical_equilibrium_light


!pure real function solve_q(a, b, c)
! real, intent(in) :: a, b, c
! real :: delta
! if (-4.*a*c/b**2 > epsilon(0.)) then
!    delta = max(b**2-4.*a*c, 0.)
!    solve_q = (-b+sqrt(delta))/(2.*a)
! else
!    solve_q = -c/b
! endif
! solve_q = max(solve_q,1e-50)
!end function solve_q

pure real function calc_Kd(coefs, T)
! all quantities are in cgs
 real, intent(in) :: coefs(5), T
 real, parameter :: R = 1.987165
 real :: G, d
 G = coefs(1)/T + coefs(2) + (coefs(3)+(coefs(4)+coefs(5)*T)*T)*T
 d = min(-G/(R*T),222.)
 calc_Kd = exp(d)
end function calc_Kd

pure real function calc_Kd_TiS(T)
! all quantities are in cgs
 real, intent(in) :: T
 real, parameter :: a = 1.3316d1, b = -6.2216, c = 4.5829d-1, d = -6.4903d-2, e = 3.2788d-3
 real :: theta, logKd
 theta = 5040./T
 logKd = a+(b+(c+(d+e*theta)*theta)*theta)*theta
 calc_Kd_TiS = 10.**(-logKd)*patm
end function calc_Kd_TiS

!----------------------------------------
!
!  Calculate mean molecular weight, gamma
!
!----------------------------------------
subroutine calc_muGamma(rho_cgs,mass_per_H,eps,T, mu, gamma, pH, pH_tot)
! all quantities are in cgs
 !use io, only:fatal

 real, intent(in)    :: rho_cgs,mass_per_H,eps(nElements)
 real, intent(inout) :: mu, gamma, T
 real, intent(out)   :: pH, pH_tot
 real :: KH2, pH2
 real :: T_old, mu_old, gamma_old, tol
 logical :: converged
 integer :: i,isolve
 integer, parameter :: itermax = 100
 character(len=30), parameter :: label = 'calc_muGamma'
 
 if (T > 1.d5) then
    mu     = (1.+4.*eps(iHe))/(1.+eps(iHe))
    gamma  = 5./3.
    pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
    pH     = pH_tot 
 elseif (T > 450.) then
! iterate to get consistently pH, T, mu and gamma
    tol       = 1.d-3
    converged = .false.
    isolve    = 0
    pH_tot    = rho_cgs*T*kboltz/(patm*mass_per_H) ! to avoid compiler warning
    pH        = pH_tot ! arbitrary value, overwritten below, to avoid compiler warning
    i = 0
    do while (.not. converged .and. i < itermax)
       i = i+1
       pH_tot    = rho_cgs*T*kboltz/(patm*mass_per_H)
       KH2       = calc_Kd(coefs(:,iH2), T)
       pH        = newton_method(0.,0.,2.*KH2,1.,-pH_tot,pH_tot) !solve_q(2.*KH2, 1., -pH_tot)
       pH2       = KH2*pH**2
       mu_old    = mu
       mu        = (1.+4.*eps(iHe))*pH_tot/(pH+pH2+eps(iHe)*pH_tot)
       gamma_old = gamma
       gamma     = (5.*pH+5.*eps(iHe)*pH_tot+7.*pH2)/(3.*pH+3.*eps(iHe)*pH_tot+5.*pH2)
       T_old     = T
       T         = T_old*mu*(gamma-1.)/(mu_old*(gamma_old-1.))
       !T        = T_old    !uncomment this line to cancel iterations
       converged = (abs(T-T_old)/T_old) < tol
       !print *,i,T_old,T,gamma_old,gamma,mu_old,mu,abs(T-T_old)/T_old
       if (i>=itermax .and. .not.converged) then
          if (isolve==0) then
             isolve = isolve+1
             i      = 0
             tol    = 1.d-2
             print *,'[dust_formation] cannot converge on T(mu,gamma). Trying with lower tolerance'
          else
             print *,'Told=',T_old,',T=',T,',gamma_old=',gamma_old,',gamma=',gamma,',mu_old=',&
                  mu_old,',mu=',mu,',dT/T=',abs(T-T_old)/T_old
            ! call fatal(label,'cannot converge on T(mu,gamma)')
          endif
       endif
    enddo
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
    pH2    = pH_tot/2.
    pH     = 0.
    mu     = (1.+4.*eps(iHe))/(0.5+eps(iHe))
    gamma  = (5.*eps(iHe)+3.5)/(3.*eps(iHe)+2.5)
 endif
end subroutine calc_muGamma



!Fortran subroutine that implements Newton's method for
!root-finding with the given polynomial function
!f(x) = a*x^4 + b*x^3 + c*x^2 + d*x + e

pure real function newton_method(aa, bb, cc, dd, ee, zz)
    real(kind=8), intent(in) :: aa, bb, cc, dd, ee, zz
    real(kind=16) :: a,b,c,d,e,z
    ! Parameters
    integer, parameter :: max_iter = 350
    real(kind=16), parameter :: tolerance = 1.0d-50
    real(kind=16), parameter :: tolerance_rel = 1.0d-50

    ! Local variables
    real(kind=16) :: fx, dfx, x
    integer :: iter

    a = aa
    b = bb
    c = cc
    d = dd
    e = ee
    z = zz

    x = z

    ! Newton's method loop
    do iter = 1, max_iter
        ! Compute function value and its derivative
        fx = a*x**4 + b*x**3 + c*x**2 + d*x + e
        dfx = 4*a*x**3 + 3*b*x**2 + 2*c*x + d

        ! Check for convergence
        if (abs(fx) < tolerance .or. abs(fx) < abs(x) * tolerance_rel) then
                exit
        end if

        ! Update x using Newton's method
        x = x - fx / dfx

    end do
    newton_method = x
      
end function newton_method

!This subroutine takes input variables a, b, c, d, e,
!and the initial guess for the root z. It then
!iteratively applies Newton's method until convergence
!or reaching the maximum number of iterations.


subroutine write_time_file(name_in, cols, time, data_in, ncols, num)
 !outputs a file over a series of dumps
 character(len=*), intent(in) :: name_in
 integer, intent(in)          :: ncols, num
 character(len=*), dimension(ncols), intent(in) :: cols
 character(len=20), dimension(ncols) :: columns
 character(len=40)             :: data_formatter, column_formatter
 character(len(name_in)+9)    :: file_name
 real, intent(in)             :: time
 real, dimension(ncols), intent(in) :: data_in
 integer                      :: i, unitnum

 write(column_formatter, "(a,I2.2,a)") "('#',2x,", ncols+1, "('[',a15,']',3x))"
 write(data_formatter, "(a,I2.2,a)") "(", ncols+1, "(2x,es18.11e2))"
 write(file_name,"(2a,i3.3,a)") name_in, '.ev'
 
 if (num == 0) then
    unitnum = 1000 
 
    open(unit=unitnum, file=file_name, status='replace')
    do i=1,ncols
       write(columns(i), "(I2,a)") i+1, cols(i)
    enddo
 
    !set column headings
    write(unitnum, column_formatter) '1  Temperature', columns(:)
    close(unit=unitnum)
 endif

 unitnum=1001+num

 open(unit=unitnum, file=file_name, position='append')

 write(unitnum,data_formatter) time, data_in(:ncols)

 close(unit=unitnum)

end subroutine write_time_file

!--------------------------------------------
!
!  Initialise mean molecular weight and gamma
!
!--------------------------------------------
subroutine init_muGamma(rho_cgs,mass_per_H,eps,T,mu,gamma,ppH,ppH2)
! all quantities are in cgs
 real, intent(in)              :: rho_cgs,mass_per_H,eps(nElements)
 real, intent(inout)           :: T
 real, intent(out)             :: mu, gamma
 real, intent(out), optional   :: ppH, ppH2
 real :: KH2, pH_tot, pH, pH2

 pH_tot = rho_cgs*kboltz*T/(patm*mass_per_H)
 if (T > 1.d5) then
    pH2 = 0.
    pH  = pH_tot
 elseif (T > 450.) then
    KH2 = calc_Kd(coefs(:,iH2), T)
    pH  = newton_method(0.,0.,2.*KH2,1.,-pH_tot,pH_tot) !solve_q(2.*KH2, 1., -pH_tot)
    pH2 = KH2*pH**2
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH2 = pH_tot/2.
    pH  = 0.
 endif
 mu    = (1.+4.*eps(iHe))*pH_tot/(pH+pH2+eps(iHe)*pH_tot)
 gamma = (5.*pH+5.*eps(iHe)*pH_tot+7.*pH2)/(3.*pH+3.*eps(iHe)*pH_tot+5.*pH2)
 call calc_muGamma(rho_cgs,mass_per_H,eps,T, mu, gamma, pH, pH_tot)
 if (present(ppH))  ppH = pH
 if (present(ppH2)) ppH2 = pH2

end subroutine init_muGamma

end module chemistry
