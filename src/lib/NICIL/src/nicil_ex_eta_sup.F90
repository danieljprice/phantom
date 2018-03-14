!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!         Example programme to test various parameters of NICIL        !
!                       (supplementary routines)                       !
!                                                                      !
!                 Copyright (c) 2015-2017 James Wurster                !
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
 integer, public, parameter  :: iprint          =  6              ! Unit to write to the screen
 integer, public, parameter  :: iprintrho       = 17              ! Unit to write evolving density (with constant T) to file
 integer, public, parameter  :: iprintbaro      = iprintrho+1     ! Unit to write evolving values using barotropic EOS to file
 integer, public, parameter  :: iprinttemp      = iprintbaro+1    ! Unit to write evolving temperature (with constant n) to file
 integer, public, parameter  :: iprintwarn      = iprinttemp+1    ! Unit to write warnings to file
 !
 real,    public             :: Bconst
 !
 public                      :: which_Bfield,get_Bfield_code,get_temperature
 public                      :: write_data_header,write_data_to_file,fatal
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
! Wardle & Ng (1999)        if version=P (i.e. a step function)
! Nakan, Nishi & Umebayashi if version=U (i.e. a uniform function)
!+
!----------------------------------------------------------------------!
real function get_Bfield_code(n_total,unit_Bfield,version)
 real,             intent(in) :: n_total,unit_Bfield
 character(len=*), intent(in) :: version
 real                         :: Bfield
 !
 if (version=="P") then
    if (n_total < 1.0d6) then
       Bfield = (n_total*1.0d-6)**0.50                ! Magnetic field [mG] for low density
    else
       Bfield = (n_total*1.0d-6)**0.25                ! Magnetic field [mG] for high density
    endif
    Bfield = Bfield * 1.0d-3                          ! Magnetic field [G]
 else if (version=="U") then
    Bfield = 1.43d-7*sqrt(n_total)                    ! Magnetic field [G]
 else if (version=="I") then
    Bfield = Bconst                                   ! Magnetic field [G] (user's constant value)
 else
    write(iprint,'(a)') "That is an invalid magnetic field function"
    call fatal(1)
 endif
 get_Bfield_code = Bfield / unit_Bfield               ! Magnetic field [code units]
 !
end function get_Bfield_code
!----------------------------------------------------------------------!
!+
! Given an input number density, this will calculate the temperature,
! based upon Machida et al. (2006) and used in Marchand et al. (2016)
!+
!----------------------------------------------------------------------!
real function get_temperature(n_total)
 real,    intent(in)    :: n_total
 real,    parameter     :: T0     = 10.0            ! [K]
 real,    parameter     :: nA     = 1.0d11          ! [cm^-3]
 real,    parameter     :: nB     = 1.0d16          ! [cm^-3]
 real,    parameter     :: nC     = 1.0d21          ! [cm^-3]
 real,    parameter     :: gammaA =  0.4
 real,    parameter     :: gammaB = -0.3
 real,    parameter     :: gammaC =  0.56667
 !
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
 write(i,"('#',34(1x,'[',i2.2,1x,a11,']',2x))") &
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
       34,'ng(-,+)/ng'
 !
end subroutine write_data_header
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
 !
 !-- Write values to file for testing purposes; note: sigma = data_out(1:3) is already in cgs
 write(iunit,'(34(1pe18.3,1x))')rho,T,Bfield*unit_Bfield &
                               ,eta_ohm*unit_eta,eta_hall*unit_eta,eta_ambi*unit_eta,data_out,fBdust
 !
end subroutine write_data_to_file
!----------------------------------------------------------------------!
!+
! Terminates the eta test programme, or prints a warning
!+
!----------------------------------------------------------------------!
subroutine fatal(ierr,rho,temperature)
 integer,           intent(in) :: ierr
 real,    optional, intent(in) :: rho,temperature
 integer                       :: i
 !
 ! Warning or fatal error encountered.  Print to file.
 if (present(rho) .and. present(temperature)) then
    write(iprintwarn,'(2(a,Es10.3),a)') "For the above error, rho = ",rho," g/cm^3 and T = ",temperature," K"
 endif
 !
 if (ierr > 0) then
    ! Fatal error encountered.  Aborting.
    write(iprint,'(a)') "NICIL: ETA TEST: fatal error encountered in NICIL.  Aborting!"
    do i = iprintrho,iprintwarn
       close(i)
    enddo
    stop
 endif
 !
end subroutine
!----------------------------------------------------------------------!
end module etasup
