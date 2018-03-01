!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!         Example programme to test various parameters of NICIL        !
!                                                                      !
!                 Copyright (c) 2015-2017 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
! This is a test programme to calculate grain charge, electron number
! densities and the non-ideal MHD coefficients for a given range of
! cosmic ray coefficients, given a fixed density and magnetic field.
!+
!----------------------------------------------------------------------!
program nicil_ex_zeta
 use nicil,  only:nicil_initialise,nicil_get_ion_n,nicil_get_eta,nicil_translate_error
 use nicil,  only:nelements_max,nelements,nlevels,zeta_cgs,n_data_out
 implicit none
 !--Input parameters
 integer, parameter  :: nlogz           = 10000           ! Number of points to test
 real,    parameter  :: logz_min        = -32.0           ! min(log(z)) to test [s^-1]
 real,    parameter  :: logz_max        = -13.0           ! max(log(z)) to test [s^-1]
 real,    parameter  :: temperature     =  13.5           ! Temperture of gas [K]
 real,    parameter  :: rho_in          = 7.43d-18        ! density to test [g/cm^3]
 real,    parameter  :: B_input         = 1.5657d-3       ! user input magnetic field [G]
 !--Input parameters (other)
 real,    parameter  :: mu0             =  2.38095236     ! Mean molecular mass; calculated with default values in NICIL
 !--Physical Constants
 real,    parameter  :: fourpi          =  12.566370614d0 ! 4pi
 real,    parameter  :: kboltz          = 1.38066d-16     ! Boltzmann constant  [erg/K]
 real,    parameter  :: mass_proton_cgs = 1.67262158d-24  ! Proton mass [g]
 real,    parameter  :: cgsmu0          = fourpi          ! Vacuum permeability [cm/s]
 !--Local variables
 integer             :: i,j,ierr
 real                :: utime,udist,umass
 real                :: unit_density,unit_ndensity,unit_charge,unit_Bfield,unit_eta,mump1
 real                :: T,rho,Bfield,temp
 real                :: dlogz,n_ion,n_total
 real                :: eta_ohm,eta_hall,eta_ambi
 real                :: n_R(4,nlogz),n_electronT(nlogz)
 real                :: data_out(n_data_out)
 character(len=1)    :: Btype
 !
 !
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
 mump1         = 1.0/(mu0*mass_proton_cgs)                ! inverse mean mass [CGS]
 !  Zero Arrays
 n_R           = 0.0
 n_electronT   = 0.0
 !  Set pseudo-grid spacing
 dlogz         = (logz_max  - logz_min )/nlogz
 !  Convert constants into code units
 Bfield        = B_input/unit_Bfield
 rho           = rho_in/unit_density
 !
 !--Open output files & write header
 open(unit=18, file="data/eta_zeta.dat")
 write(18,"('#',16(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'zeta',       &
        2,'density',    &
        3,'temp',       &
        4,'B',          &
        5,'eta_ohm',    &
        6,'eta_Hall',   &
        7,'eta_ambi',   &
        8,'n_e/(nn+ni)',&
        9,'n_i/(nn+ni)',&
       10,'n_electron', &
       11,'n_neutral',  &
       12,'n_ionR_H',   &
       13,'n_ionR_M',   &
       14,'n_grainR_n', &
       15,'n_grainR_0', &
       16,'n_grainR_p'
 write(*,'(a)')     "NICIL: ETA TEST"
 !
 !--This will calculate the coefficients, grain charge and electron number density
 !  for a range of zeta
 do i = 1,nlogz
    zeta_cgs = 10**(logz_min +(i-1)*dlogz )
    !--Call the initialisation routine.
    !  Required continually since zeta is assumed to be fixed and hard-coded into the constants
    call nicil_initialise(utime,udist,umass,unit_Bfield,ierr)
    if (ierr/=0) call fatal(ierr)                            ! Abort programme if there are errors in the setup
    !
    ! Call NICIL to get the eta coefficients.
    ! data_out is an optional out variable that we will pass through to track the values
    call nicil_get_ion_n(rho,temperature,n_R(1:4,i),n_electronT(i),ierr)
    if (ierr/=0) then
       call nicil_translate_error(ierr)
       call fatal(ierr)
    endif
    call nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,temperature,n_R(1:4,i),n_electronT(i),ierr,data_out)
    if (ierr/=0) then
       call nicil_translate_error(ierr)
       call fatal(ierr)
    endif
    !
    ! Write values to file for testing purposes
    n_ion   = data_out(8) + data_out(9) + data_out(10) + data_out(11)
    n_total = n_ion + data_out(7)
    write(18,'(16(1pe18.3,1x))') zeta_cgs,rho_in,temperature,B_input, &
   eta_ohm*unit_eta,eta_hall*unit_eta,eta_ambi*unit_eta, &
   data_out(6)/n_total, n_ion/n_total,data_out(6:9)*unit_ndensity,data_out(12:14)*unit_ndensity
 enddo
 write(*,'(a)') "NICIL: ETA TEST: properties vs number density (with constant T) written to data/eta_density.dat"
 !
 close(18)
 !
!----------------------------------------------------------------------!
end program nicil_ex_zeta
!----------------------------------------------------------------------!
!+
! Terminates the eta test programme, or prints a warning
!+
!----------------------------------------------------------------------!
subroutine fatal(ierr)
 integer,           intent(in) :: ierr
 integer                       :: i
 !
 if (ierr > 0) then
    ! Fatal error encountered.  Aborting.
    write(iprint,'(a)') "NICIL: ETA TEST: fatal error encountered in NICIL.  Aborting!"
    close(18)
    stop
 endif
 !
end subroutine
!----------------------------------------------------------------------!
