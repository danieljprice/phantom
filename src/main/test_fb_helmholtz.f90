 
 !! forward backward test for eos_helmholtz_compute_pres_sound, Ali, give a temperature and density, get the energy, 
 !! then use that energy to get the temperature back, and check if we get the same temperature back 
   !  if (firstcall) then
   !    T_input = 77143143.d0  ! T and rho taken from MESA profile
   !    cgsrhotest = 21746381.d0 ! cgs units
   !    write(*,*) 'Ali l 431 (eos), equationofstate: forward backward test for eos_helmholtz_compute_pres_sound: T_input = ',&
   !     T_input, ' cgsrhotest = ', cgsrhotest
   !    call eos_helmholtz_compute_pres_sound(T_input, cgsrhotest, cgsprestest,&
   !     cgsspsoundtest, cgsen_eosforward, cgsdendttest)
   !    !! above should be in cgs
   !    T_recovered = 80143142.d0 ! initial guess
   !    rhotest= cgsrhotest/unit_density
   !    u_target_code_units = cgsen_eosforward / unit_ergg ! u is in cgs, convert to code units for input into eos_helmholtz_pres_sound
   !    !! following is in code units, so convert to cgs
   !    !! backwards test: give the energy from the forward test, and see if we get the same temperature back
   !    write(*,*) 'Ali l 446 (eos),  backward inputs: T_recovered = ',&
   !     T_recovered, ' cgsrhotest = ', cgsrhotest, ' u_target_code_units = ', u_target_code_units
   !    call eos_helmholtz_pres_sound(T_recovered, rhotest, ponrhocodeunits, cs_codeunits, u_target_code_units) ! u here is an input, so use the energy from the forward test

   !       ! ---------------------------
   !       ! STEP 3: report error
   !       ! ---------------------------
   !    write(*,*) 'Ali INPUT T        = ', T_input
   !    write(*,*) 'Ali RECOVERED T    = ', T_recovered
   !    write(*,*) 'Ali REL ERROR      = ', abs(T_recovered - T_input)/T_input
   !    firstcall = .false.
   !  endif

