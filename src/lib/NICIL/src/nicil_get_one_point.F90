!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!   Interactive programme to calculate the results at one data point   !
!                                                                      !
!                 Copyright (c) 2015-2021 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
! This is an interactive programme that will return the non-ideal MHD
! coefficients, number densities and ionisation rates given an input
! density, temperature and magnetic field.  Input data must be in CGS.
!+
!----------------------------------------------------------------------!
program nicil_get_one_point
 use nicil, only:nicil_initialise,nicil_update_nimhd,nicil_translate_error,n_warn
 use nicil, only:n_nden,n_data_out
 implicit none
 !--Physical Constants
 real,    parameter  :: fourpi          =  12.566370614d0 ! 4pi
 real,    parameter  :: cgsmu0          = fourpi          ! Vacuum permeability [cm/s]
 integer, parameter  :: iprint          = 6               ! unit number to write to the screen
 integer             :: k,ierr,ierrlist(n_warn)
 real                :: utime,udist,umass
 real                :: unit_density,unit_ndensity,unit_charge,unit_Bfield,unit_eta
 real                :: rho_in_cgs,B_in_cgs,T_in,rho_in,B_in
 real                :: nions,ngrains,ngrains1
 real                :: eta_ohm,eta_hall,eta_ambi
 real                :: ndens(n_nden),data_out(n_data_out)
 character(len=128)  :: c_rho_in_cgs,c_B_in_cgs,c_T_in

 !--Get input data (rho,T,B), either from command line or by prompting user
 call getarg(1,c_rho_in_cgs)
 if (trim(c_rho_in_cgs)=="") then
    write(6,'(''Enter the density in g cm^-3: '')')
    read(5,'(a)') c_rho_in_cgs
 endif
 read(c_rho_in_cgs,*) rho_in_cgs
 call getarg(2,c_T_in)
 if ( trim(c_T_in) == "" ) then
    write(*,'(''Enter the temperature in K: '')')
    read(5,'(a)') c_T_in
 endif
 read(c_T_in,*) T_in
 call getarg(3,c_B_in_cgs)
 if ( trim(c_B_in_cgs) == "" ) then
    write(*,'(''Enter the magnetic field strength in G: '')')
    read(5,'(a)') c_B_in_cgs
 endif
 read(c_B_in_cgs,*) B_in_cgs

 !--Initialise parameters
 !  Code Units
 udist         = 1.000d16                                 ! = 1 code length unit
 umass         = 1.989d33                                 ! = 1 M_sun = 1 code mass unit
 utime         = 8.681d10                                 ! = 1 code time unit (Chosen such that G=1)
 unit_density  = umass/udist**3                           ! = 1 code density unit
 unit_ndensity = 1.0/udist**3                             ! = 1 code number density unit
 unit_charge   = sqrt(umass*udist/cgsmu0)                 ! = 1 code charge unit
 unit_Bfield   = umass/(utime*unit_charge)                ! = 1 code magentic field unit *may be defined differently in user's code*
 unit_eta      = udist**2/utime                           ! = 1 code eta unit
 ! Input units in code units
 rho_in        = rho_in_cgs/unit_density
 B_in          = B_in_cgs  /unit_Bfield
 ! Initialise variables
 ndens         = 0.
 ierrlist      = 0

 !--Call the NICIL routines
 call nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint)
 call nicil_update_nimhd(0,eta_ohm,eta_hall,eta_ambi,B_in,rho_in,T_in,ndens,ierrlist,data_out)
 if (any(ierrlist > 0) ) call nicil_translate_error(ierrlist,.false.)
 nions = 0.
 do k = 9,19
    nions = nions + data_out(k)
 enddo
 ngrains = data_out(20) + data_out(21) + data_out(22)  ! Only true for single grain population
 if (ngrains > 0.) then
    ngrains1 = 1./ngrains
 else
    ngrains1 = 0.
 endif

 !--Print out useful quantities
 write(iprint,'(a)')          "--------------------------------------------------------------------------------"
 write(iprint,'(a)')          "The input parameters are "
 write(iprint,'(a,Es18.11,a)') " rho_0 = ",rho_in_cgs," g cm^-3"
 write(iprint,'(a,Es18.11,a)') " T_0   = ",T_in,      " K"
 write(iprint,'(a,Es18.11,a)') " B_0   = ",B_in_cgs,  " G"
 write(iprint,'(a)')          "The output parameters are "
 write(iprint,'(a,Es18.11,a)') " eta_ohm    = ",eta_ohm      *unit_eta,     " cm^2 s^-1"
 write(iprint,'(a,Es18.11,a)') " eta_hall   = ",eta_hall     *unit_eta,     " cm^2 s^-1"
 write(iprint,'(a,Es18.11,a)') " eta_ambi   = ",eta_ambi     *unit_eta,     " cm^2 s^-1"
 write(iprint,'(a,Es18.11,a)') " n_electron = ",data_out(8)  *unit_ndensity," cm^-3"
 write(iprint,'(a,Es18.11,a)') " n_e_therm  = ",ndens(n_nden)*unit_ndensity," cm^-3"
 write(iprint,'(a,Es18.11,a)') " n_ions     = ",nions        *unit_ndensity," cm^-3"
 write(iprint,'(a,Es18.11,a)') " n_neutral  = ",data_out(5)  *unit_ndensity," cm^-3"
 write(iprint,'(a,Es18.11  )') " negative grain fraction: n_g(Z=-1)/n_g = ",data_out(22)*ngrains1
 write(iprint,'(a,Es18.11  )') " neutral  grain fraction: n_g(Z= 0)/n_g = ",data_out(21)*ngrains1
 write(iprint,'(a,Es18.11  )') " positive grain fraction: n_g(Z=+1)/n_g = ",data_out(20)*ngrains1
 write(iprint,'(a,Es18.11  )') " ionisation fraction:     n_i/(n_i+n_n) = ",nions/(nions+data_out(7))
 write(iprint,'(a)')          "--------------------------------------------------------------------------------"

end program nicil_get_one_point
