module grids_for_setup

 use fileutils, only:load_data_file
 use table_utils, only: differentiate
 use io,       only:warning,error

 implicit none
 public   init_grid_sigma,init_grid_ecc,deallocate_sigma,deallocate_ecc
 real, dimension(:,:), allocatable :: dataecc, datasigma
 real, dimension(:), allocatable :: dsigmadx, deda, ddeda !second derivative
 logical :: ecc_initialised=.false.,sigma_initialised=.false.

 contains
 
 subroutine init_grid_sigma(Rin,Rout)
 !--initialisation occurs in setup_disc.f90, within routine surface_density_profile()
 !--but needs to be reinitialised for every disc as R_in R_out need to be rescaled
    real, intent(in) :: Rin,Rout

    sigma_initialised=.true.
    print*,'init grid sigma with Rin, Rout:',Rin,Rout
    call load_data_file('sigma_grid.dat',datasigma,nhead=1)
    call rescale(Rin,Rout,datasigma)
    call differentiate(datasigma(:,2),datasigma(:,1),dsigmadx)
 end subroutine init_grid_sigma

 subroutine init_grid_ecc(Rin,Rout)
 !--initialisation occurs in set_disc.f90 within the routine set_disc_positions() 
    real, intent(in) :: Rin,Rout

    ecc_initialised=.true.
    print*,'init grid ecc with Rin, Rout:',Rin,Rout
    call load_data_file('ecc_grid.dat',dataecc,nhead=1)
    call rescale(Rin,Rout,dataecc)
    call differentiate(dataecc(:,2),dataecc(:,1),deda)
    call differentiate(deda,dataecc(:,1),ddeda)
 end subroutine init_grid_ecc

 subroutine rescale(Rin,Rout,dataset)
    real, intent(in) :: Rin,Rout
    real, dimension(:,:), intent(inout) :: dataset
    real :: x(size(dataset(:,1))),xin,xout
    integer :: Nsize 

    Nsize=size(dataset(:,1))
    x(:)=dataset(:,1)
    xin=x(1)
    xout=x(Nsize)
    dataset(:,1)=(x(:)-xin)/(xout-xin)*(Rout-Rin)+Rin
    !print*,dataset(:,1)

 end subroutine rescale    

 subroutine deallocate_sigma()
    if(sigma_initialised) then
        deallocate(datasigma)
        deallocate(dsigmadx)
        sigma_initialised=.false.
        write(*,*) '(grids_for_setup) Deallocating: datasigma, dsigmadx'
    else
       call error('grids_for_setup','Trying to deallocate datasigma without having initialised it')
    endif
 end subroutine deallocate_sigma

  subroutine deallocate_ecc()
    if(ecc_initialised) then
         deallocate(dataecc)
         deallocate(deda)
         deallocate(ddeda)
         ecc_initialised=.false.
         write(*,*) '(grids_for_setup) Deallocating: dataecc, deda, ddeda'
    else
      call error('grids_for_setup','Trying to deallocate dataecc without having initialised it')
    endif
 end subroutine deallocate_ecc


end module grids_for_setup
