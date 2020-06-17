!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!         Example programme to test various parameters of NICIL        !
!                                                                      !
!                 Copyright (c) 2015-2019 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
! This is a test programme to calculate grain charge, electron number
! densities and the non-ideal MHD coefficients for a given range of
! densities.  Results from this programme can be directly compared to
! the author's output found in the data folder.
!+
!----------------------------------------------------------------------!
program nicil_ex_eta
 use nicil,  only: nicil_initialise,nicil_get_ion_n,nicil_get_eta,nicil_translate_error
 use nicil,  only: warn_verbose,n_data_out,zeta_of_rho,zeta_cgs
 use nicil,  only: use_fdg_in,zeta_cgs
 use etasup, only: Bconst,which_Bfield,get_Bfield_code,get_temperature
 use etasup, only: write_data_header,write_data_to_file,fatal
 use etasup, only: write_phase_header,write_phase_to_file
 use etasup, only: iprint,iprintdat,iprintwarn
 implicit none
 !--Input parameters
 integer, parameter  :: nlogT           = 10000           ! Number of temperature points to test
 integer, parameter  :: nlogd           = 10000           ! Number of density points to test
 integer, parameter  :: nlogz           = 10000           ! Number of zeta_cr points to test
 integer, parameter  :: nlogp           = 100             ! Number of points to test for phase-space plot
 real,    parameter  :: rho_in          = 1.0d-13         ! density to test [g/cm^3]
 real,    parameter  :: rho_in_low      = 7.43d-18        ! low density to test [g/cm^3] for zeta test
 real,    parameter  :: logrho_min      = -22.0           ! min(log(rho)) to test [g/cm^3]
 real,    parameter  :: logrho_max      =   0.5           ! max(log(rho)) to test [g/cm^3]
 real,    parameter  :: temp_in         =  30.0           ! Temperature of gas [K]
 real,    parameter  :: temp_13         =  13.5           ! Temperature of gas [K] for the zeta test
 real,    parameter  :: temp_min        = 1.0d1           ! minimum temperature to test [K]
 real,    parameter  :: temp_max        = 2.0d5           ! maximum temperature to test [K]
 logical             :: use_input_B     = .false.         ! use user's B (true) or use pre-defined function for B (false)
 real,    parameter  :: B_in            = 1.5657d-3       ! user input magnetic field [G] (set for zeta; previous default = 1.583d-4G)
 real,    parameter  :: logB_min        =  -7.0           ! min(log(B)) to test [G]
 real,    parameter  :: logB_max        =   5.0           ! max(log(B)) to test [G]
 real,    parameter  :: logz_min        = -32.0           ! min(log(z)) to test [s^-1]
 real,    parameter  :: logz_max        =  -8.0           ! max(log(z)) to test [s^-1]
 real,    parameter  :: fdg             = 0.01            ! dust fraction if use_fdg_in = .true.; must be between 0 & 1
 !--Input parameters (other)
 real,    parameter  :: mu0             =  2.38095236     ! Mean molecular mass; calculated with default values in NICIL
 !--Physical Constants
 real,    parameter  :: fourpi          =  12.566370614d0 ! 4pi
 real,    parameter  :: kboltz          = 1.38066d-16     ! Boltzmann constant  [erg/K]
 real,    parameter  :: mass_proton_cgs = 1.67262158d-24  ! Proton mass [g]
 real,    parameter  :: cgsmu0          = fourpi          ! Vacuum permeability [cm/s]
 !--Tests to be run
 logical             :: test_vs_rho     = .true.          ! Run the test using constant temperature
 logical             :: test_vs_T       = .true.          ! Run the test using constant density
 logical             :: test_baro       = .true.          ! Run the test using barotropic EOS
 logical             :: test_colI       = .false.         ! Run the test using T & B from ideal MHD collapse to stellar densities
 logical             :: test_colN       = .false.         ! Run the test using T & B from non-ideal MHD collapse to stellar densities
 logical             :: test_zeta       = .false.         ! Run the test varying zeta_cr
 logical             :: test_phase      = .false.         ! Generate rho-B phase-space
 integer, parameter  :: ntestmax        = 7               ! number of tests
 integer, parameter  :: irho            = 1               ! Index of the test using constant temperature
 integer, parameter  :: iT              = 2               ! Index of the test using constant density
 integer, parameter  :: ibaro           = 3               ! Index of the test using barotropic EOS
 integer, parameter  :: icolI           = 4               ! Index of the test using T & B from ideal MHD collapse to stellar densities
 integer, parameter  :: icolN           = 5               ! Index of the test using T & B from non-ideal MHD collapse to stellar densities
 integer, parameter  :: izeta           = 6               ! Index of the test varying zeta_cr
 integer, parameter  :: iphase          = 7               ! Index for generating rho_B phase space
 !--Local variables
 integer             :: i,j,imax,jmax,nicil_loop,ierr
 integer             :: c0,c1,count_rate,count_max
 real                :: utime,udist,umass
 real                :: unit_density,unit_ndensity,unit_charge,unit_Bfield,unit_eta,mump1
 real                :: T,rho,rho_cgs,Bfield,temp,fBdust,zeta0
 real                :: dlogrho,dlogT,dlogB,dlogz
 real                :: eta_ohm,eta_hall,eta_ambi
 real                :: n_R(4),n_electronT
 real                :: data_out(n_data_out)
 logical             :: do_calculations
 character(len=  1)  :: Btype
 character(len= 10)  :: arg_in
 character(len=128)  :: descrip
 !
 !--Overwrite what tests to use if input given
 call getarg(1,arg_in)
 if ( trim(arg_in) /= "" ) then
    test_vs_rho = .false.
    test_vs_T   = .false.
    test_baro   = .false.
    test_colI   = .false.
    test_colN   = .false.
    test_zeta   = .false.
    test_phase  = .false.
    j = 1
    do while (trim(arg_in) /= "" )
       do i = 1,ntestmax
          if (trim(arg_in)=='all' .or. trim(arg_in)=='vsRho') test_vs_rho = .true.
          if (trim(arg_in)=='all' .or. trim(arg_in)=='vsT')   test_vs_T   = .true.
          if (trim(arg_in)=='all' .or. trim(arg_in)=='baro')  test_baro   = .true.
          if (trim(arg_in)=='all' .or. trim(arg_in)=='colI')  test_colI   = .true.
          if (trim(arg_in)=='all' .or. trim(arg_in)=='colN')  test_colN   = .true.
          if (trim(arg_in)=='all' .or. trim(arg_in)=='zeta')  test_zeta   = .true.
          if (trim(arg_in)=='all' .or. trim(arg_in)=='phase') test_phase  = .true.
       enddo
       j = j + 1
       call getarg(j,arg_in)
    enddo
 endif
 !
 !--Initialise parameters
 udist         = 1.000d16                                 ! = 1 code length unit
 umass         = 1.989d33                                 ! = 1 M_sun = 1 code mass unit
 utime         = 8.681d10                                 ! = 1 code time unit (Chosen such that G=1)
 unit_density  = umass/udist**3                           ! = 1 code density unit
 unit_ndensity = 1.0/udist**3                             ! = 1 code number density unit
 unit_charge   = sqrt(umass*udist/cgsmu0)                 ! = 1 code charge unit
 unit_Bfield   = umass/(utime*unit_charge)                ! = 1 code magentic field unit *may be defined differently in user's code*
 unit_eta      = udist**2/utime                           ! = 1 code eta unit
 mump1         = 1.0/(mu0*mass_proton_cgs)                ! inverse mean mass [CGS]
 warn_verbose  = .true.                                   ! print out all warnings for testing purposes
 fBdust        = 0.0                                      ! Fraction of magnetised dust (not always returned)
 zeta0         = zeta_cgs                                 ! Archive original cosmic ray ionisation rate
 !  Rename the input magnetic field for use in supplementary routines
 Bconst        = B_in
 !  Use input dust-to-gas ratio
 use_fdg_in    = .false.
 !
 !--Call the initialisation routine.
 !  Required only once to calculate and save all the values required for this library
 call nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint,iprintwarn)
 if (ierr/=0) call fatal(ierr)                            ! Abort programme if there are errors in the setup
 !
 !--Open output file of warning long
 open(unit=iprintwarn,file="data/eta_warning.log")
 write(iprintwarn,'(a)') "NICIL: WARNINGS LOG"
 write(iprint,'(a)')     "NICIL: "
 !
 !--This will calculate the coefficients, grain charge and electron number density
 !  We repeat this calculation several times, using different parameters
 do nicil_loop = 1,ntestmax
    n_R             = 0.0           ! (re-)initialise the array
    n_electronT     = 0.0           ! (re-)initialise the array
    imax            = nlogd
    jmax            = 1
    do_calculations = .false.
    if (nicil_loop == irho) then
       if (test_vs_rho) then
          ! A range of densities at a fixed temperatures;
          ! magnetic field and sound speed are related to the density by a given prescription
          do_calculations = .true.
          open(unit=iprintdat, file="data/eta_density.dat")
          descrip = "NICIL: ETA TEST: properties vs number density (with constant T)"
          Btype   = which_Bfield(use_input_B,"P")
          temp    = temp_in
       endif
    elseif (nicil_loop == iT) then
       if (test_vs_T) then
          ! A range of temperatures at a fixed density and magnetic field strength
          do_calculations = .true.
          open(unit=iprintdat,file="data/eta_temperature.dat")
          descrip = "NICIL: ETA TEST: properties vs temperature (with constant density)"
          imax    = nlogT
          Btype   = which_Bfield(use_input_B,"P")
          rho_cgs = rho_in
       endif
    elseif (nicil_loop == ibaro) then
       if (test_baro) then
          ! A barotropic equation of state
          do_calculations = .true.
          open(unit=iprintdat,file="data/eta_barotropic.dat")
          descrip = "NICIL: ETA TEST: properties vs {number density, temperature} using barotropic EOS"
          Btype   = which_Bfield(use_input_B,"U")
       endif
    elseif (nicil_loop == icolI) then
       if (test_colI) then
          ! The rho_max-B_max-T_max relation using ideal MHD in Wurster, Bate & Price (2018d)
          do_calculations = .true.
          open(unit=iprintdat,file="data/eta_collapseIdeal.dat")
          descrip = "NICIL: ETA TEST: properties vs {number density, temperature} using ideal MHD collapse properties"
          Btype   = which_Bfield(use_input_B,"D")
       endif
    elseif (nicil_loop == icolN) then
       if (test_colN) then
          ! The rho_max-B_max-T_max relation using ideal MHD in Wurster, Bate & Price (2018d)
          do_calculations = .true.
          open(unit=iprintdat,file="data/eta_collapseNonideal.dat")
          descrip = "NICIL: ETA TEST: properties vs {number density, temperature} using non-ideal MHD collapse properties"
          Btype   = which_Bfield(use_input_B,"N")
       endif
    elseif (nicil_loop == izeta) then
       if (test_zeta) then
          ! Varying zeta_cr while holding all other values constant
          do_calculations = .true.
          open(unit=iprintdat,file="data/eta_zeta.dat")
          descrip = "NICIL: ZETA TEST:"
          temp    = temp_13
          Btype   = which_Bfield(use_input_B,"I")
       endif
    elseif (nicil_loop == iphase) then
       if (test_phase) then
          ! Generate rho-B phase space
          do_calculations = .true.
          open(unit=iprintdat,file="data/eta_phase.dat")
          descrip = "NICIL: PHASE SPACE:"
          imax    = nlogp
          jmax    = nlogp
       endif
    endif

    !  Set pseudo-grid spacing
    dlogrho = (logrho_max  - logrho_min )/imax         ! spacing of density points
    dlogT   = (log10(temp_max) - log10(temp_min))/imax ! spacing of temperature points
    dlogB   = (logB_max   - logB_min  )/jmax           ! spacing of magnetic field points
    dlogz   = (logz_max  - logz_min )/imax             ! spacing of zeta_cr points

    if (do_calculations) then
       write(iprintwarn,'(a)') descrip
       if (test_phase) then
          call write_phase_header(iprintdat)
       else
          call write_data_header(iprintdat)
       endif
       call system_clock(c0,count_rate,count_max)
       do i = 1,imax
          ! Calculate the new data point in rho-B-T phase space
          rho_cgs = 10**(logrho_min +(i-1)*dlogrho )                      ! density [CGS]
          if (nicil_loop==iT) then
             temp    = 10**(log10(temp_min) +(i-1)*dlogT )                ! temperature
             rho_cgs = rho_in                                             ! density [CGS]
          elseif (nicil_loop==ibaro) then
             temp    = get_temperature(rho_cgs,mump1,"M")                 ! temperature
          elseif (nicil_loop==icolI .or. nicil_loop==icolN .or. nicil_loop==iphase) then
             temp    = get_temperature(rho_cgs,mump1,"S")                 ! temperature
          elseif (nicil_loop==izeta) then
             rho_cgs  = rho_in_low
             zeta_cgs = 10**(logz_min +(i-1)*dlogz )
             !--Call the initialisation routine.
             !  Required continually since zeta is assumed to be fixed and hard-coded into the constants
             call nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint,iprintwarn)
          endif
          rho = rho_cgs/unit_density                                      ! density [code units]
          do j = 1,jmax
             if (nicil_loop==iphase) then
                Bfield = 10**(logB_min +(j-1)*dlogB )/unit_Bfield         ! magnetic field [code units]
             else
                Bfield = get_Bfield_code(rho_cgs,mump1,unit_Bfield,Btype) ! magnetic field  [code units]
             endif

             ! Calculate the number densities
             if (use_fdg_in) then
                call nicil_get_ion_n(rho,temp,n_R,n_electronT,ierr,fdg,fBdust)
             else
                call nicil_get_ion_n(rho,temp,n_R,n_electronT,ierr)
             end if
             if (ierr/=0) then
                call nicil_translate_error(ierr)
                call fatal(ierr,rho*unit_density,temp)
             end if

             ! calculate the eta coefficients.
             ! data_out is an optional out variable that we will pass through to track the values
             if (use_fdg_in) then
                call nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,temp,n_R,n_electronT,ierr,data_out,fdg)
             else
                call nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,temp,n_R,n_electronT,ierr,data_out)
             end if
             if (ierr/=0) then
                call nicil_translate_error(ierr)
                call fatal(ierr,rho*unit_density,temp)
             end if

             ! Write values to file for testing purposes
             if (nicil_loop==iphase) then
                call write_phase_to_file(iprintdat,rho_cgs,temp,Bfield,eta_ohm,eta_hall,eta_ambi,unit_eta,unit_Bfield)
             else
                ! Write values to file for testing purposes
                if (nicil_loop==izeta .or. .not. zeta_of_rho) then
                   ! the unit is required to counteract the unit in the write routine
                   data_out(size(data_out)) = zeta_cgs*utime
                endif
                call write_data_to_file(iprintdat,rho_cgs,temp,Bfield,eta_ohm,eta_hall,eta_ambi,fBdust &
                                       ,data_out,unit_eta,unit_Bfield,unit_density,unit_ndensity,utime)
             endif
          enddo
          if (nicil_loop == iphase) write(iprintdat,'(a)') " "
       end do
       call system_clock(c1,count_rate,count_max)
       write(iprint,'(a)') descrip
       write(iprint,'(a,F6.3,a)') "NICIL: Runtime: ",(c1-c0)/float(count_rate)," seconds"
       close(iprintdat)
       if (nicil_loop==izeta) then
          zeta_cgs = zeta0
          !--Call the initialisation routine to reset the default cosmic ray ionisation rate
          call nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint,iprintwarn)
       endif
    endif
 enddo
 close(iprintwarn)
!----------------------------------------------------------------------!
end program nicil_ex_eta
