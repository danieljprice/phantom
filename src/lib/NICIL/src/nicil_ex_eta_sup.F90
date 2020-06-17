!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!         Example programme to test various parameters of NICIL        !
!                       (supplementary routines)                       !
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
module etasup
 implicit none
 !
 ! unit number for writing to the file/screen
 integer, public, parameter  :: iprint     =  6              ! Unit to write to the screen
 integer, public, parameter  :: iprintdat  = 17              ! Unit to write data file
 integer, public, parameter  :: iprintwarn = 18              ! Unit to write error file
 !
 real,    public             :: Bconst
 !
 public                      :: which_Bfield,get_Bfield_code,get_temperature
 public                      :: write_data_header,write_data_to_file,fatal
 public                      :: write_phase_header,write_phase_to_file
 private
 !
contains
!+
!----------------------------------------------------------------------!
!+
! Determine which magnetic field geometry to use:
! constant from user or functional from functional input
!+
!----------------------------------------------------------------------!
character(len=1) function which_Bfield(use_input_B,Bopt)
 logical, intent(in) :: use_input_B
 character(len=*)    :: Bopt
 !
 if (use_input_B) then
    which_Bfield = "I"
 else
    which_Bfield = Bopt
 endif
 !
end function which_Bfield
!----------------------------------------------------------------------!
!+
! Given an input number density, this will calculate the magnetic field
! strength, as used in
! P: Wardle & Ng (1999)        (i.e. a step function)
! U: Nakan, Nishi & Umebayashi (i.e. a uniform function)
! D: B_max from     ideal MHD curve in fig 1 of Wurster, Bate & Price (2018d)
! N: B_cen from non-ideal MHD curve in fig 1 of Wurster, Bate & Price (2018d)
!+
!----------------------------------------------------------------------!
real function get_Bfield_code(rho_cgs,mump1,unit_Bfield,version)
 real,             intent(in) :: rho_cgs,mump1,unit_Bfield
 character(len=*), intent(in) :: version
 real                         :: Bfield,n_total,B0
 real                         :: nA,nB,nC,nD,nE,nF
 real                         :: gammaA,gammaB,gammaC,gammaD,gammaE,gammaF
 !
 ! Common parameters for versions D & N
 B0      = 1.146e-7 ! [G]
 nA      = 1.e-22   ! [g cm^-3]
 gammaA  = 0.65
 ! Convert mass density to number density
 n_total = rho_cgs*mump1
 if (version=="P") then
    if (n_total < 1.0d6) then
       Bfield = (n_total*1.0d-6)**0.50                ! Magnetic field [mG] for low density
    else
       Bfield = (n_total*1.0d-6)**0.25                ! Magnetic field [mG] for high density
    endif
    Bfield = Bfield * 1.0d-3
 else if (version=="U") then
    Bfield = 1.43d-7*sqrt(n_total)
 else if (version=="I") then
    Bfield = Bconst                                   ! Magnetic field [G] (user's constant value)
 else if (version=="D") then
    ! Piecewise that has been decommissioned for a smooth function:
    !     if (rho_cgs < 4.29e-15) Bfield = 10.**( -3.67 + ( 0.649 )*( log10(rho_cgs) + 16.96 ))
    ! elseif (rho_cgs < 2.88e-13) Bfield = 10.**( -2.00 + ( 0.239 )*( log10(rho_cgs) + 14.42 ))
    ! elseif (rho_cgs < 5.21e-10) Bfield = 10.**( -0.963+ ( 0.783 )*( log10(rho_cgs) + 11.79 ))
    ! elseif (rho_cgs < 9.87e-10) Bfield = 10.
    ! else                        Bfield = 10.**(  1.13 + ( 0.614 )*( log10(rho_cgs) +  8.794))
    nB     =  4.29e-15; nC     = 2.88e-13; nD      = 5.21e-10; nE     =  9.87e-10
    gammaB = -0.42;     gammaC = 0.5;       gammaD = 0.2;      gammaE = -0.32
    Bfield = B0 *         (rho_cgs/nA)  **gammaA &
                * ( 1.0 + (rho_cgs/nB) )**gammaB &
                * ( 1.0 + (rho_cgs/nC) )**gammaC &
                * ( 1.0 + (rho_cgs/nD) )**gammaD &
                * ( 1.0 + (rho_cgs/nE) )**gammaE
 else if (version=="N") then
    ! Piecewise that has been decommissioned for a smooth function:
    !     if (rho_cgs < 2.67e-15) Bfield = 10.**( -3.67 + ( 0.649 )*( log10(rho_cgs) + 16.96 ))
    ! elseif (rho_cgs < 4.43e-13) Bfield = 10.**( -2.07 + ( 0.229 )*( log10(rho_cgs) + 14.35 ))
    ! elseif (rho_cgs < 1.69e-12) Bfield = 10.**( -1.41 + ( 0.907 )*( log10(rho_cgs) + 12.13 ))
    ! elseif (rho_cgs < 4.08e-10) Bfield = 10.**( -1.08 + ( 0.392 )*( log10(rho_cgs) + 11.76 ))
    ! elseif (rho_cgs < 3.45e-09) Bfield = 0.707
    ! else                        Bfield = 10.**( -0.0683+( 0.537 )*( log10(rho_cgs) +  8.312))
    nB     =  2.67e-15; nC     = 4.43e-13; nD     =  1.69e-12; nE     =  6.0e-10; nF     = 3.45e-09
    gammaB = -0.42;     gammaC = 0.50;     gammaD = -0.25;     gammaE = -1.0;     gammaF =  1.1
    Bfield = B0 *         (rho_cgs/nA)  **gammaA &
                * ( 1.0 + (rho_cgs/nB) )**gammaB &
                * ( 1.0 + (rho_cgs/nC) )**gammaC &
                * ( 1.0 + (rho_cgs/nD) )**gammaD &
                * ( 1.0 + (rho_cgs/nE) )**gammaE &
                * ( 1.0 + (rho_cgs/nF) )**gammaF
 else
    write(iprint,'(a)') "That is an invalid magnetic field function"
    call fatal(1)
 endif
 get_Bfield_code = Bfield / unit_Bfield               ! Magnetic field [code units]
 !
end function get_Bfield_code
!----------------------------------------------------------------------!
!+
! Given an input density, this will calculate the temperature, using
! M: Machida et al. (2006) and used in Marchand et al. (2016)
! S: fig 2 of Wurster, Bate & Price (2018a)
!+
!----------------------------------------------------------------------!
real function get_temperature(rho_cgs,mump1,version)
 real,             intent(in) :: rho_cgs,mump1
 character(len=*), intent(in) :: version
 real                         :: T0                    ! [K]
 real                         :: nA,nB,nC              ! [cm^-3]
 real                         :: gammaA,gammaB,gammaC
 real                         :: n_total
 !
 if (version=="M") then
    T0     = 10.0
    nA     = 1.0d11; nB     =  1.0d16; nC     = 1.0d21
    gammaA = 0.4;    gammaB = -0.3;    gammaC = 0.56667
 elseif (version=="S") then
    ! Piecewise that has been decommissioned for a smooth function:
    ! if     (rho_cgs < 7.58d-14) get_temperature = 13.96
    ! elseif (rho_cgs < 1.657d-9) get_temperature = 10.**( 1.85 + ( 0.478 )*( log10(rho_cgs) + 11.65 ))
    ! elseif (rho_cgs < 1.545d-4) get_temperature = 10.**( 3.33 + ( 0.0992)*( log10(rho_cgs) +  7.66 ))
    ! else                        get_temperature = 10.**( 3.99 + ( 0.421 )*( log10(rho_cgs) +  3.15 ))
    T0     = 14.0
    nA     = 2.0e10; nB     =  2.5e14; nC     = 1.0e20
    gammaA = 0.5;    gammaB = -0.4;    gammaC = 0.37
 else
    write(iprint,'(a)') "That is an invalid temperature function"
    call fatal(1)
 endif
 n_total = rho_cgs*mump1
 get_temperature = T0 * sqrt( 1.0 + (n_total/nA)**(2.0*gammaA))         &
                      *     ( 1.0 + (n_total/nB)              )**gammaB &
                      *     ( 1.0 + (n_total/nC)              )**gammaC
 !
end function get_temperature
!----------------------------------------------------------------------!
!+
! Subroutine to write the header to a data file
!+
!----------------------------------------------------------------------!
subroutine write_data_header(i)
 integer, intent (in) :: i
 !
 write(i,"('#',36(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'density',    &
        2,'temp',       &
        3,'B',          &
        4,'eta_ohm',    &
        5,'eta_Hall',   &
        6,'eta_ambi',   &
        7,'sigma_O',    &
        8,'sigma_H',    &
        9,'sigma_P',    &
       10,'rho_n',      &
       11,'rho_ion',    &
       12,'n_electron', &
       13,'n_neutral',  &
       14,'n_ionR_H',   &
       15,'n_ionR_M',   &
       16,'n_ionT_1st', &
       17,'n_ionT_2nd', &
       18,'n_grainR_n', &
       29,'n_grainR_0', &
       20,'n_grainR_p', &
       21,'n_H2',       &
       22,'n_H',        &
       23,'n_H2+',      &
       24,'n_H+',       &
       25,'n_He+',      &
       26,'n_Na+',      &
       27,'n_Mg+',      &
       28,'n_K+',       &
       29,'n_He++',     &
       30,'n_Na++',     &
       31,'n_Mg++',     &
       32,'n_K++',      &
       33,'zeta',       &
       34,'ng(-,+)/ng', &
       35,'n_ion',      &
       36,'n_total'
 !
end subroutine write_data_header
!----------------------------------------------------------------------!
!+
! Subroutine to write the header to a phase-phase data file
!+
!----------------------------------------------------------------------!
subroutine write_phase_header(i)
 integer, intent (in) :: i

 write(i,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'density',  &
        2,'B',        &
        3,'temp',     &
        4,'eta_ohm',  &
        5,'eta_Hall', &
        6,'eta_ambi', &
        7,'max term'

end subroutine write_phase_header
!----------------------------------------------------------------------!
!+
! Subroutine to write data to file
!+
!----------------------------------------------------------------------!
subroutine write_data_to_file(iunit,rho,T,Bfield,eta_ohm,eta_hall,eta_ambi,fBdust, &
                              data_out,unit_eta,unit_Bfield,unit_density,unit_ndensity,utime)
 integer, intent(in)    :: iunit
 real,    intent(in)    :: rho,T,Bfield,eta_ohm,eta_hall,eta_ambi,fBdust
 real,    intent(in)    :: unit_eta,unit_Bfield,unit_density,unit_ndensity,utime
 real,    intent(inout) :: data_out(:)
 integer                :: j,idata
 real                   :: n_ion,n_total
 real,    parameter     :: density_thresh  = 1.0d-60   ! parameter to prevent overflow errors when writing mass densities
 real,    parameter     :: ndensity_thresh = 1.0d-60   ! parameter to prevent overflow errors when writing number densities
 !
 !--Convert mass & number density to cgs units here; done to prevent numerical underflow
 do j = 4,5
    if (data_out(j) < density_thresh) then
       data_out(j) = 0.0
    else
       data_out(j) = data_out(j)*unit_density
    endif
 enddo
 idata = size(data_out)
 do j = 6,idata-1
    if (data_out(j) < ndensity_thresh) then
       data_out(j) = 0.0
    else
       data_out(j) = data_out(j)*unit_ndensity
    endif
 enddo
 data_out(idata) = data_out(idata)/utime
 ! these quantities are required for the zeta test (Wurster, Bate & Price 2018b)
 n_ion           = data_out(8) + data_out(9) + data_out(10) + data_out(11)
 n_total         = n_ion + data_out(7)

 !-- Write values to file for testing purposes; note: sigma = data_out(1:3) is already in cgs
 write(iunit,'(36(1pe18.3,1x))')rho,T,Bfield*unit_Bfield &
                               ,eta_ohm*unit_eta,eta_hall*unit_eta,eta_ambi*unit_eta,data_out,fBdust,n_ion,n_total

end subroutine write_data_to_file
!----------------------------------------------------------------------!
!+
! Subroutine to write phase space data to file
!+
!----------------------------------------------------------------------!
subroutine write_phase_to_file(iunit,rho,temp,Bfield,eta_ohm,eta_hall,eta_ambi,unit_eta,unit_Bfield)
 integer, intent(in)    :: iunit
 real,    intent(in)    :: rho,temp,Bfield,eta_ohm,eta_hall,eta_ambi
 real,    intent(in)    :: unit_eta,unit_Bfield
 integer                :: emax

 emax = 0
 if (    eta_ohm   > abs(eta_hall) .and.     eta_ohm   >     eta_ambi ) emax =  1
 if (abs(eta_hall) >     eta_ohm   .and. abs(eta_hall) >     eta_ambi ) emax =  2
 if (    eta_ambi  >     eta_ohm   .and.     eta_ambi  > abs(eta_hall)) emax =  3
 if (emax == 2                     .and.     eta_hall  < 0.0          ) emax = -2
 write(iunit,'(6(1pe18.3,1x),I18,1x)') rho,Bfield*unit_Bfield, temp, &
                                       eta_ohm*unit_eta,eta_hall*unit_eta,eta_ambi*unit_eta,emax

end subroutine write_phase_to_file
!----------------------------------------------------------------------!
!+
! Terminates the eta test programme, or prints a warning
!+
!----------------------------------------------------------------------!
subroutine fatal(ierr,rho,temperature)
 integer,           intent(in) :: ierr
 real,    optional, intent(in) :: rho,temperature
 !
 ! Warning or fatal error encountered.  Print to file.
 if (present(rho) .and. present(temperature)) then
    write(iprintwarn,'(2(a,Es10.3),a)') "For the above error, rho = ",rho," g/cm^3 and T = ",temperature," K"
 endif
 !
 if (ierr > 0) then
    ! Fatal error encountered.  Aborting.
    write(iprint,'(a)') "NICIL: ETA TEST: fatal error encountered in NICIL.  Aborting!"
    close(iprintdat)
    close(iprintwarn)
    stop
 endif
 !
end subroutine
!----------------------------------------------------------------------!
end module etasup
