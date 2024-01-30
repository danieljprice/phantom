!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_molecular
!
! Handles molecular cooling of binary wind
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: datafiles, eos, infile_utils, physcon, units
!

 implicit none

 real                                :: CO_abun  = 1.e-4
 real                                :: H2O_abun = 5.e-5
 real                                :: HCN_abun = 3.e-7
 real, dimension(40, 40, 40, 6)      :: coolingTable
 real, dimension(36, 102, 6, 8, 5)   :: cdTable
 logical                             :: do_molecular_cooling = .false.
 real                                :: fit_rho_power, fit_rho_inner, fit_vel, r_compOrb
 real                                :: Tfloor = 0. ![K]  to prevent u < 0 (this is independent of Tfloor in cooling.F90)

contains

!-----------------------------------------------------------------------
!+
!  write the molecular cooling options for the parameter card
!+
!-----------------------------------------------------------------------
subroutine write_options_molecularcooling(iunit)

 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(CO_abun, 'CO_abun', 'set to value>0 to activate CO radiative cooling &
 & (typical value O-rich AGB star=1e-4)',iunit)
 call write_inopt(HCN_abun,'HCN_abun','set to value>0 to activate HCN radiative cooling &
 & (typical value O-rich AGB star=1e-7)',iunit)
 call write_inopt(H2O_abun,'H2O_abun','set to value>0 to activate H2O radiative cooling &
 & (typical value O-rich AGB star=5e-5)',iunit)

end subroutine write_options_molecularcooling

!-----------------------------------------------------------------------
!+
!  read the molecular cooling options from the parameter card for restart
!+
!-----------------------------------------------------------------------
subroutine read_options_molecular_cooling(name,valstring,imatch,igotall,ierr)

 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer :: ngot = 0

 imatch  = .true.
 igotall = .true. ! none of the cooling options are compulsory

 select case(trim(name))
 case('CO_abun')
    read(valstring,*,iostat=ierr) CO_abun
    if (CO_abun < 0.) CO_abun = 0.
    ngot = ngot + 1
 case('HCN_abun')
    read(valstring,*,iostat=ierr) HCN_abun
    if (HCN_abun < 0.) HCN_abun = 0.
    ngot = ngot + 1
 case('H2O_abun')
    read(valstring,*,iostat=ierr) H2O_abun
    if (H2O_abun < 0.) H2O_abun = 0.
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 do_molecular_cooling = ngot > 0 .and. CO_abun+H2O_abun+HCN_abun > 0.

end subroutine read_options_molecular_cooling

!-----------------------------------------------------------------------
!+
!  Initialise molecular cooling tables
!+
!-----------------------------------------------------------------------
subroutine init_cooling_molec

 call loadCoolingTable(coolingTable)
 call loadCDTable(cdTable)

end subroutine init_cooling_molec

!-----------------------------------------------------------------------
!+
!  Calculate the molecular cooling
!+
!-----------------------------------------------------------------------
subroutine calc_cool_molecular( T, r, rho_sph, Q, dlnQdlnT)

 use physcon, only:atomic_mass_unit,kboltz,mass_proton_cgs,au
 use eos,     only:gmw
 use units,   only:udist

! Data dictionary: Arguments
 real, intent(out)  :: Q, dlnQdlnT           ! In CGS and linear scale
 real, intent(in)   :: T, rho_sph            ! In CGS
 real, intent(in)   :: r                     ! code units

! Data dictionary: Additional parameters for calculations
 integer                                     :: i
 real                                        :: rho_H, n_H, Temp, Lambda, r_au
 real                                        :: fit_n_inner,T_log, n_H_log, N_hydrogen, N_coolant_log
 real                                        :: abundance, widthLine_molecule
 real, dimension(3)                          :: lambda_log, params_cool, widthLine,Qi
 real, dimension(4)                          :: params_cd
 real, dimension(3)                          :: molecular_abun
 real, dimension(3), parameter               :: mass_molecules = [28.01, 18.01528, 27.0253]*atomic_mass_unit
 character(len=3), dimension(3), parameter   :: moleculeNames = ["CO ", "H2O", "HCN"]
 character(len=3)                            :: moleculeName

! Initialise variables
 r_au = r*udist/au
 if (r_au <= r_compOrb) then
    rho_H = fit_rho_inner
 else
    rho_H = fit_rho_inner*(r_compOrb/r_au)**fit_rho_power
 endif

 Temp=T
 if (T<=Tfloor) Temp=Tfloor

 n_H                 = rho_H/(gmw*mass_proton_cgs)    ! convert mass density to number density
 fit_n_inner         = fit_rho_inner/(gmw*mass_proton_cgs)
 i                   = -1
 T_log               = log10(Temp)
 n_H_log             = log10(n_H)
 widthLine           = sqrt(2. * kboltz * Temp / mass_molecules) / fit_vel
 N_coolant_log       = -999.
 params_cool         = -999.
 params_cd           = 0.
 lambda_log          = -999.
 molecular_abun      = [CO_abun, H2O_abun, HCN_abun]

! Calculate column density
 do i = 1, 3
    Lambda             = 0.
    widthLine_molecule = widthLine(i)
    abundance          = molecular_abun(i)
    moleculeName       = moleculeNames(i)
    params_cd          = [r_au, widthLine_molecule, fit_rho_power, r_compOrb]
    call ColumnDensity(cdTable, params_cd, N_hydrogen)

    if (N_hydrogen /= 0.) then
       N_coolant_log      = log10(abundance * fit_n_inner * N_hydrogen)
       ! Calculate cooling rate
       params_cool = [T_log, n_H_log, N_coolant_log]

       call CoolingRate(coolingTable, params_cool, moleculeName, lambda_log(i))
       Lambda = Lambda + 10.**lambda_log(i)
    endif
    Qi(i) = -Lambda*abundance*rho_sph/(gmw*mass_proton_cgs)     ! in erg/sec
 enddo

 Q = Qi(1)+Qi(2)+Qi(3)

 if (Q /= 0.) then
    call lambdaGradT(coolingTable, params_cool, Q, dlnQdlnT)
 else
    dlnQdlnT = 0
 endif

end subroutine calc_cool_molecular

!-----------------------------------------------------------------------
!+
!  Load radiative cooling table
!+
!-----------------------------------------------------------------------
subroutine loadCoolingTable(data_array)
 use datafiles,  only:find_phantom_datafile

 real, dimension(40, 40, 40, 6), intent(out) :: data_array

 ! Data dictionary: Read radiative cooling file
 integer            :: i, j, k, o, iunit, istat
 character(len=80)  :: imsg
 integer, parameter :: headerLines = 5
 real               :: T, n_H, N_coolant, lambda_CO, lambda_H2O, lambda_HCN
 character(len=120) :: filename


 ! Initialise variables
 i          = 0
 j          = 0
 k          = 0
 o          = 0
 T          = 0.
 n_H        = 0.
 N_coolant  = 0.
 lambda_CO  = -999.
 lambda_H2O = -999.
 lambda_HCN = -999.
 data_array = -999.

 iunit = 1
 filename = find_phantom_datafile('radcool_all.dat','cooling')
 OPEN(unit=iunit, file=trim(filename), STATUS="OLD", ACTION="read", &
            iostat=istat, IOMSG=imsg)

 ! Begin loading in data
 openif: if (istat == 0) then
    !!! Skip header
    rewind(unit=iunit)
    do o = 1, headerLines
       read(iunit, *, iostat=istat, IOMSG = imsg)
    enddo

    ! Read data
    skipheaderif: if ((istat == 0)) then
       readdo: do
          read(iunit, *, iostat=istat) i, j, k, T, n_H, N_coolant, lambda_CO, lambda_H2O, lambda_HCN
          if (istat /= 0) exit
          data_array(i, j, k, :) = [T, n_H, N_coolant, lambda_CO, lambda_H2O, lambda_HCN]

       enddo readdo

       if (istat > 0) write(*, *) "Error at line ", i, j, k, " during loading data into array."

    else
       write(*, 100) headerLines
100    format("Error: Header consists of more than ", I2, " lines.")
       write(*, *) trim(imsg)
    endif skipheaderif

 else
    write(*, *) "Error: Radiative cooling table ", trim(filename) ," does not exist."
 endif openif
end subroutine loadCoolingTable

!-----------------------------------------------------------------------
!+
!  Load column density table
!+
!-----------------------------------------------------------------------
subroutine loadCDTable(data_array)
 use datafiles,  only:find_phantom_datafile
 real, dimension(36, 102, 6, 8, 5), intent(out) :: data_array

 ! Data dictionary: Read radiative cooling file
 integer            :: i, j, k, l, o, iunit, istat
 character(len=80)  :: imsg
 integer, parameter :: headerLines = 8
 real               :: r_part, widthLine, m_exp, r_sep, N_H
 character(len=120) :: filename


 ! Initialise variables
 i           = 0
 j           = 0
 k           = 0
 l           = 0
 o           = 0
 r_part      = 0.
 widthLine   = 0.
 m_exp       = 0.
 r_sep       = 0.
 N_H         = 0
 data_array  = -999.

 iunit = 1
 filename = find_phantom_datafile('table_cd.dat','cooling')
 open(unit=iunit, file=filename, STATUS="OLD", iostat=istat, IOMSG=imsg)

 ! Begin loading in data
 openif: if (istat == 0) then
    !!! Skip header
    rewind(unit=iunit)
    do o = 1, headerLines
       read(iunit, *, iostat=istat, IOMSG = imsg)
    enddo

    !!! Read data
    skipheaderif: if ((istat == 0)) then
       readdo: do
          read(iunit, *, iostat=istat) i, j, k, l, r_part, widthLine, m_exp, r_sep, N_H
          if (istat /= 0) exit
          data_array(i, j, k, l, :) = [r_part, widthLine, m_exp, r_sep, N_H]

       enddo readdo

       if (istat > 0) write(*, *) "Error at line ", i, j, k, l, " during loading data into array."

    else
       write(*, 100) headerLines
100    format("Error: Header consists of more than ", I2, " lines.")
       write(*, *) trim(imsg)
    endif skipheaderif

 else
    write(*, *) "Error: Column density table ", trim(filename)," does not exist."
 endif openif
end subroutine loadCDTable

!-----------------------------------------------------------------------
!+
!  Calculate the gradient of the cooling rate with varying temperature
!+
!-----------------------------------------------------------------------
subroutine lambdaGradT(data_array, params, dudt, dlnQdlnT)

 ! Data dictionary: Arguments
 real, dimension(40, 40, 40, 6), intent(in)  :: data_array
 real, dimension(3), intent(in)              :: params
 real, intent(in)                            :: dudt
 real, intent(out)                           :: dlnQdlnT

 ! Data dictionary: Additional arguments for calculations
 real                    :: T_lower, T_upper, dT, T, dQ
 integer, dimension(3)   :: index_lower_bound
 real, dimension(3)      :: dQ_lower, dQ_upper
 integer                 :: index_T_lower, index_T_upper, index_nH, index_N
 integer                 :: index_molecule

 ! Initialisation
 index_molecule = 0
 dlnQdlnT       = 0.

 call findLower_cool(data_array, params, index_lower_bound)
 index_T_lower = index_lower_bound(1)
 index_T_upper = index_lower_bound(1) + 1
 index_nH      = index_lower_bound(2)
 index_N       = index_lower_bound(3)

 T            = 10.**params(1)
 T_lower      = 10.**data_array(index_T_lower, index_nH, index_N, 1)
 T_upper      = 10.**data_array(index_T_upper, index_nH, index_N, 1)
 dT           = T_upper - T_lower
 dQ_lower     = 10.**data_array(index_T_lower, index_nH, index_N, 4:6)
 dQ_upper     = 10.**data_array(index_T_upper, index_nH, index_N, 4:6)
 dQ           = sum(dQ_upper) - sum(dQ_lower)
 dlnQdlnT     = T / dudt * dQ / dT
end subroutine lambdaGradT

!-----------------------------------------------------------------------
!+
!  Find closest and strictly smaller indices of parameters in cooling table
!+
!-----------------------------------------------------------------------
subroutine findLower_cool(data_array, params, index_lower_bound)

 ! Data dictionary: Arguments
 real, dimension(3), intent(in)              :: params
 integer, dimension(3), intent(out)          :: index_lower_bound
 real, dimension(40, 40, 40, 6), intent(in)  :: data_array
 real, parameter                             :: tol = 1.E-6

 ! Data dictionary: Find index values
 real, dimension(3) :: params0, dparams
 real, dimension(3) :: real_lower_bound

 params0(:) = data_array(1, 1, 1, 1:3)
 dparams(:) = data_array(2, 2, 2, 1:3) - params0(:)

 ! Calculate index of the parameters in the cooling table
 real_lower_bound  = (params - params0)/dparams
 index_lower_bound = floor(real_lower_bound) + 1

end subroutine findLower_cool

!-----------------------------------------------------------------------
!+
!  Find closest and strictly smaller indices of parameters in column density table
!+
!-----------------------------------------------------------------------
subroutine findLower_cd(data_array, params, index_lower_bound)

 ! Data dictionary: Arguments
 real, dimension(4), intent(in)                  :: params
 integer, dimension(4), intent(out)              :: index_lower_bound
 real, dimension(36, 102, 6, 8, 5), intent(in)   :: data_array

 ! Data dictionary: Find index values
 integer, parameter  :: N_compZone = 15, N_r_part_sample = 36
 real                :: dr_part, r_sep
 real, parameter     :: dv = 0.02, dm = 0.3, widthLine_min = 0.001, widthLine_max = 2., perturbation = 0.0001
 real, parameter     :: r_part_min = 1.1, r_part_max = 250.
 integer             :: i, j, k, l
 real, dimension(4)  :: min_array, max_array

 i = -1
 j = -1
 k = -1
 l = -1
 index_lower_bound = -1

 ! The upper boundary of widthLine is limitless but is minimally equal to zero
 min_array = [10.**data_array(1, 1, 1, 1, 1), 0., data_array(1, 1, 1, 1, 3), data_array(1, 1, 1, 1, 4)]
 max_array = [10.**data_array(36, 102, 6, 8, 1), params(2) + 1, data_array(36, 102, 6, 8, 3), data_array(36, 102, 6, 8, 4)]

 if (all(params <= max_array) .AND. all(min_array <= params)) then
    ! Index r_sep
    l   = nint( params(4) - min_array(4) ) + 1
    r_sep = params(4)

    ! Index m_exp
    k = floor((params(3) - min_array(3)) / dm) + 1

    ! Index widthLine
    outervif: if (params(2) >= widthLine_min) then
       innervif: if (params(2) >= widthLine_max) then
          j = 102
       elseif (params(2) >= 1.999) then
          j = 101
       else
          j = floor((params(2))/dv) + 1
       endif innervif
    endif outervif

    ! Index r_part
    outCompif: if (params(1) >= r_sep) then
       innerr_partif: if (params(1) <= r_sep + perturbation) then
          i = N_compZone
       else
          !dr_part = (log10(r_part_max) - log10(r_sep)) / (N_r_part_sample - N_compZone)
          dr_part = log10(r_part_max/r_sep) / (N_r_part_sample - N_compZone)
          i    = floor(log10(params(1)/r_sep) / dr_part ) + N_compZone
       endif innerr_partif

    else
       dr_part = (r_sep - r_part_min) / N_compZone
       i    = floor((params(1) - r_part_min) / dr_part ) + 1
       if (params(1) >= r_sep / 2.) i = i + 1
    endif outCompif

    index_lower_bound = [i, j, k, l]
 endif
end subroutine findLower_cd

!-----------------------------------------------------------------------
!+
!  Determine the corresponding radiative cooling rate for given parameters of the local gas
!+
!-----------------------------------------------------------------------
subroutine CoolingRate(data_array, params, molecule, lambda)

 ! Data dictionary: Arguments
 real, intent(out)                           :: lambda
 real, dimension(3), intent(in)              :: params
 real, dimension(40, 40, 40, 6), intent(in)  :: data_array
 character(len=3), intent(in)                :: molecule

 ! Data dictionary: Find lambda value
 integer                 :: i
 logical, dimension(3)   :: checkBoundaries
 integer, dimension(3)   :: index_lower_bound

 ! Data dictionary: Interpolation
 real                     :: x_d, y_d, z_d
 real, dimension(3)       :: xyz
 real, dimension(2, 2, 2) :: c_grid_CO, c_grid_H2O, c_grid_HCN
 integer                  :: i0, i1, j0, j1, k0, k1

 ! Initialise variables
 lambda = -999.
 i = -1
 index_lower_bound = -1
 checkBoundaries = .false.
 x_d = -1.
 y_d = -1.
 z_d = -1.
 xyz = -1.
 c_grid_CO = -999.
 c_grid_H2O = -999.
 c_grid_HCN = -999.
 i0 = -1
 i1 = -1
 j0 = -1
 j1 = -1
 k0 = -1
 k1 = -1

 ! Find T1, n_H1, N_coolant1 (which are the lower boundaries for interpolation)
 call findLower_cool(data_array, params, index_lower_bound)
 i0 = index_lower_bound(1)
 j0 = index_lower_bound(2)
 k0 = index_lower_bound(3)

 ! Check two cases: either the indices are in bounds of the data or not
 checkBoundaries = 1 <= index_lower_bound .and. index_lower_bound <= 40

 findlowerif: if (all(checkBoundaries .eqv. .true.)) then
    ! We now have to check whether i, j or k are not located at the upper boundary of the data

    ! Index i
    indexiif: if (i0 == 40) then
       i1  = i0
       x_d = 0.
    else
       i1  = i0 + 1
       x_d = (params(1) - data_array(i0, 1, 1, 1)) / (data_array(i1, 1, 1, 1) - data_array(i0, 1, 1, 1))
    endif indexiif

    ! Index j
    indexjif: if (j0 == 40) then
       j1  = j0
       y_d = 0.
    else
       j1  = j0 + 1
       y_d = (params(2) - data_array(1, j0, 1, 2)) / (data_array(1, j1, 1, 2) - data_array(1, j0, 1, 2))
    endif indexjif

    ! Index k
    indexkif: if (k0 == 40) then
       k1  = k0
       z_d = 0.
    else
       k1  = k0 + 1
       z_d = (params(3) - data_array(1, 1, k0, 3)) / (data_array(1, 1, k1, 3) - data_array(1, 1, k0, 3))
    endif indexkif

    xyz = [x_d, y_d, z_d]

    moleculesif: if (trim(molecule) == "CO") then
       c_grid_CO  = data_array(i0:i1, j0:j1, k0:k1, 4)
       call interpolate(xyz, c_grid_CO, lambda)

    elseif (molecule == "H2O") then
       c_grid_H2O = data_array(i0:i1, j0:j1, k0:k1, 5)
       call interpolate(xyz, c_grid_H2O, lambda)

    elseif (molecule == "HCN") then
       c_grid_HCN = data_array(i0:i1, j0:j1, k0:k1, 6)
       call interpolate(xyz, c_grid_HCN, lambda)
    endif moleculesif

 endif findlowerif

end subroutine CoolingRate

!-----------------------------------------------------------------------
!+
!  Determine the column density based on fit parameters and location of particle
!+
!-----------------------------------------------------------------------
subroutine ColumnDensity(data_array, params, N_hydrogen)

 ! Data dictionary: Arguments
 real, intent(out)                               :: N_hydrogen
 real, dimension(4), intent(in)                  :: params
 real, dimension(36, 102, 6, 8, 5), intent(in)   :: data_array

 ! Data dictionary: Find column density
 integer                 :: i
 logical, dimension(4)   :: checkBoundaries
 integer, dimension(4)   :: index_lower_bound, max_boundary

 ! Data dictionary: Interpolation
 real                     :: x_d, y_d, z_d
 real, dimension(3)       :: xyz
 real, dimension(2, 2, 2) :: c_grid
 integer                  :: i0, i1, j0, j1, k0, k1, l

 ! Initialise variables
 i = -1
 index_lower_bound = -1
 checkBoundaries = .false.
 x_d = -1.
 y_d = -1.
 z_d = -1.
 xyz = -1.
 c_grid = 0.
 i0 = -1
 i1 = -1
 j0 = -1
 j1 = -1
 k0 = -1
 k1 = -1
 l  = -1
 max_boundary = [36, 102, 6, 8]
 N_hydrogen = 0.

 ! Find T1, n_H1, N_coolant1 (which are the lower boundaries for interpolation)
 call findLower_cd(data_array, params, index_lower_bound)
 i0 = index_lower_bound(1)
 j0 = index_lower_bound(2)
 k0 = index_lower_bound(3)
 l  = index_lower_bound(4)

 ! Check two cases: either the indices are in bounds of the data or not
 checkBoundaries = 1 <= index_lower_bound .and. index_lower_bound <= max_boundary
 findlowerif: if (all(checkBoundaries .eqv. .true.)) then
    ! We now have to check whether i, j or k are not located at the upper boundary of the data

    ! Index i
    indexiif: if (i0 == 36) then
       i1  = i0
       x_d = 0.
    else
       i1  = i0 + 1
       x_d = (log10(params(1)) - data_array(i0, 1, 1, l, 1)) / (data_array(i1, 1, 1, l, 1) - data_array(i0, 1, 1, l, 1))
    endif indexiif

    ! Index j
    indexjif: if (j0 == 102) then
       j1  = j0
       y_d = 0.
    else
       j1  = j0 + 1
       y_d = (params(2) - data_array(1, j0, 1, l, 2)) / (data_array(1, j1, 1, l, 2) - data_array(1, j0, 1, l, 2))
    endif indexjif

    ! Index k
    indexkif: if (k0 == 6) then
       k1  = k0
       z_d = 0.
    else
       k1  = k0 + 1
       z_d = (params(3) - data_array(1, 1, k0, l, 3)) / (data_array(1, 1, k1, l, 3) - data_array(1, 1, k0, l, 3))
    endif indexkif

    c_grid  = data_array(i0:i1, j0:j1, k0:k1, l, 5)
    xyz = [x_d, y_d, z_d]
    call interpolate(xyz, c_grid, N_hydrogen)
 endif findlowerif

end subroutine ColumnDensity

!-----------------------------------------------------------------------
!+
!  Three-dimensional interpolation
!+
!-----------------------------------------------------------------------
subroutine interpolate(xyz, c_grid, c_interpol)

 ! Data dictionary: Arguments
 real, dimension(3), intent(in)  :: xyz
 real, dimension(2, 2, 2)        :: c_grid
 real, intent(out)               :: c_interpol

 ! Data dictionary: Additional parameters for interpolation
 real    :: c00, c01, c10, c11, c0, c1, x_d, y_d, z_d
 integer :: i0, i1, j0, j1, k0, k1

 x_d = xyz(1)
 y_d = xyz(2)
 z_d = xyz(3)

 c00 = 0.
 c01 = 0.
 c10 = 0.
 c11 = 0.

 i0 = 1
 i1 = 2
 j0 = 1
 j1 = 2
 k0 = 1
 k1 = 2

 c00 = c_grid(i0, j0, k0) * (1. - y_d) + c_grid(i0, j1, k0) * y_d
 c01 = c_grid(i0, j0, k1) * (1. - y_d) + c_grid(i0, j1, k1) * y_d
 c10 = c_grid(i1, j0, k0) * (1. - y_d) + c_grid(i1, j1, k0) * y_d
 c11 = c_grid(i1, j0, k1) * (1. - y_d) + c_grid(i1, j1, k1) * y_d

 c0 = c00 * (1. - z_d) + c01 * z_d
 c1 = c10 * (1. - z_d) + c11 * z_d

 c_interpol = c0 * (1. - x_d) + c1 * x_d

end subroutine interpolate

end module cooling_molecular
