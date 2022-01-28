module grids_for_setup

 implicit none
 public
 real, dimension(:,:), allocatable :: dataecc, datasigma

 contains
 
 subroutine init_grids_setup()
    use load_from_file, only:load_data_file
 
    call load_data_file('sigma_grid.dat',datasigma,nhead=1)
    call load_data_file('ecc_grid.dat',dataecc,nhead=1)

 end subroutine init_grids_setup

end module grids_for_setup
