!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module ptmass_radiation
!
! Implementation of radiation from sink particles
!   Contains routines to compute dust temperature assuming radiative equilibrium
!   Also routine to compute radiative acceleration based on sink particle luminosity
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - alpha_rad       : *fraction of the gravitational acceleration imparted to the gas*
!   - iget_tdust      : *dust temperature (0:Tdust=Tgas 1:T(r) 2:Flux dilution 3:Attenuation 4:Lucy)*
!   - isink_radiation : *sink radiation pressure method (0=off,1=alpha,2=dust,3=alpha+dust)*
!   - tdust_exp       : *exponent of the dust temperature profile*
!
! :Dependencies: dim, dust_formation, infile_utils, io, part, raytracer,
!   units
!


 implicit none
 integer, public  :: isink_radiation = 0
 integer, public  :: iget_tdust      = 0
 integer, public  :: iray_resolution = -1
 real,    public  :: tdust_exp       = 0.5
 real,    public  :: alpha_rad       = 0.

 public :: get_rad_accel_from_ptmass
 public :: read_options_ptmass_radiation,write_options_ptmass_radiation
 public :: get_dust_temperature
 public :: init_radiation_ptmass

 private

contains
!-----------------------------------------------------------------------
!+
!  Initialisation (if needed)
!+
!-----------------------------------------------------------------------
subroutine init_radiation_ptmass(ierr)
 integer, intent(out) :: ierr

 ierr = 0


end subroutine init_radiation_ptmass

!-----------------------------------------------------------------------
!+
!  compute radiative acceleration from ALL sink particles
!+
!-----------------------------------------------------------------------
subroutine get_rad_accel_from_ptmass (nptmass,npart,i,xi,yi,zi,xyzmh_ptmass,fextx,fexty,fextz,tau,fsink_old,extrapfac)
 use part,    only:ilum
 use units,   only:umass,unit_luminosity
 integer,        intent(in)    :: nptmass,npart,i
 real,           intent(in)    :: xi,yi,zi
 real,           intent(in)    :: xyzmh_ptmass(:,:)
 real, optional, intent(in)    :: tau(:)
 real,           intent(inout) :: fextx,fexty,fextz
 real, optional, intent(in)    :: fsink_old(:,:)
 real, optional, intent(in)    :: extrapfac
 real                    :: dx,dy,dz,Mstar_cgs,Lstar_cgs
 integer                 :: j
 logical                 :: extrap

 if (present(fsink_old)) then
    extrap = .true.
 else
    extrap = .false.
 endif

 do j=1,nptmass
    if (xyzmh_ptmass(4,j) < 0.) cycle
    Mstar_cgs  = xyzmh_ptmass(4,j)*umass
    Lstar_cgs  = xyzmh_ptmass(ilum,j)*unit_luminosity
    !compute radiative acceleration if sink particle is assigned a non-zero luminosity
    if (Lstar_cgs > 0.d0) then
       if (extrap) then
          dx = xi - xyzmh_ptmass(1,j) + extrapfac*fsink_old(1,j)
          dy = yi - xyzmh_ptmass(2,j) + extrapfac*fsink_old(2,j)
          dz = zi - xyzmh_ptmass(3,j) + extrapfac*fsink_old(3,j)
       else
          dx = xi - xyzmh_ptmass(1,j)
          dy = yi - xyzmh_ptmass(2,j)
          dz = zi - xyzmh_ptmass(3,j)
       endif
       call calc_rad_accel_from_ptmass(npart,i,dx,dy,dz,Lstar_cgs,Mstar_cgs,fextx,fexty,fextz,tau)
    endif
 enddo

end subroutine get_rad_accel_from_ptmass

!-----------------------------------------------------------------------
!+
!  compute radiative acceleration on all particles
!+
!-----------------------------------------------------------------------
subroutine calc_rad_accel_from_ptmass(npart,i,dx,dy,dz,Lstar_cgs,Mstar_cgs,fextx,fexty,fextz,tau)
 use part,  only:isdead_or_accreted,dust_temp,nucleation,idkappa,idalpha
 use dim,   only:do_nucleation,itau_alloc
 use dust_formation, only:calc_kappa_bowen
 integer,           intent(in)    :: npart,i
 real, optional,    intent(in)    :: tau(:)
 real,              intent(in)    :: dx,dy,dz,Lstar_cgs,Mstar_cgs
 real,              intent(inout) :: fextx,fexty,fextz
 real                             :: r,ax,ay,az,alpha,kappa


 r = sqrt(dx**2 + dy**2 + dz**2)
 if (do_nucleation) then
    if (itau_alloc == 1) then
       call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
               nucleation(idkappa,i),ax,ay,az,nucleation(idalpha,i),tau(i))
    else
       call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
               nucleation(idkappa,i),ax,ay,az,nucleation(idalpha,i))
    endif
 else
    kappa = calc_kappa_bowen(dust_temp(i))
    if (itau_alloc == 1) then
       call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
               kappa,ax,ay,az,alpha,tau(i))
    else
       call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
               kappa,ax,ay,az,alpha)
    endif
 endif
 fextx = fextx + ax
 fexty = fexty + ay
 fextz = fextz + az

end subroutine calc_rad_accel_from_ptmass


!-----------------------------------------------------------------------
!+
!  compute radiative acceleration from a SINGLE sink particle
!  based on sink particle luminosity and computed opacities / column depth
!+
!-----------------------------------------------------------------------
subroutine get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
     kappa,ax,ay,az,alpha,tau_in)
 use units,          only:umass
 use dust_formation, only:calc_Eddington_factor
 real, intent(in)            :: r,dx,dy,dz,Mstar_cgs,Lstar_cgs,kappa
 real, intent(in), optional  :: tau_in
 real, intent(out)           :: ax,ay,az,alpha
 real :: fac,tau

 if (present(tau_in)) then
    tau = tau_in
 else
    tau = 0.
 endif
 select case (isink_radiation)
 case (1)
    ! alpha wind
    alpha = alpha_rad
 case (2)
    ! radiation pressure on dust
    alpha = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, kappa, tau)
 case (3)
    ! radiation pressure on dust + alpha_rad (=1+2)
    alpha = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, kappa, tau) + alpha_rad
 case default
    ! no radiation pressure
    alpha = 0.
 end select
 fac = alpha*Mstar_cgs/(umass*r**3)
 ax = fac*dx
 ay = fac*dy
 az = fac*dz
end subroutine get_radiative_acceleration_from_star

!-----------------------------------------------------------------------
!+
!  compute dust temperature by performing simplistic ray tracing from sink particles
!  through the gas, or by simpler approximations
!+
!-----------------------------------------------------------------------
subroutine get_dust_temperature(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp)
 use part,      only:tau,tau_lucy,ikappa,nucleation
 use raytracer, only:get_all_tau
 use dust_formation, only:calc_kappa_bowen,idust_opacity
 use dim,       only:itau_alloc
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:),xyzmh_ptmass(:,:),eos_vars(:,:)
 real,     intent(out)   :: dust_temp(:)

 !
 ! compute dust temperature based on previous value of tau or tau_lucy
 !
 if (iget_tdust == 4) then
    ! calculate the dust temperature using the value of tau_Lucy from the last timestep
    call get_dust_temperature_from_ptmass(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp,tau_lucy=tau_lucy)
 elseif (iget_tdust == 3) then
    ! calculate the dust temperature using attenuation of stellar flux (exp(-tau)) with the "standard" tau from the last timestep.
    call get_dust_temperature_from_ptmass(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp,tau=tau)
 else
    ! other case : T(r) relation, Flux dilution or Tdust = Tgas
    call get_dust_temperature_from_ptmass(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp)
 endif
 !
 ! do ray tracing to get optical depth : calculate new tau, tau_lucy
 !
 if (iget_tdust == 4) then
    ! update tau_Lucy
    if (idust_opacity == 2) then
       call get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, nucleation(:,ikappa), iray_resolution, tau_lucy)
    else
       call get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, calc_kappa_bowen(dust_temp(1:npart)), iray_resolution, tau_lucy)
    endif
 elseif (itau_alloc == 1) then
    ! update tau
    if (idust_opacity == 2) then
       call get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, nucleation(:,ikappa), iray_resolution, tau)
    else
       call get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, calc_kappa_bowen(dust_temp(1:npart)), iray_resolution, tau)
    endif
 endif
 !
 ! update Tdust with new optical depth. This step gives more consistency but may not be needed. To be checked
 !
 if (iget_tdust == 4) then
    ! update dust temperature with new tau_Lucy
    call get_dust_temperature_from_ptmass(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp,tau_lucy=tau_lucy)
 elseif (iget_tdust == 3) then
    ! update dust temperature using new tau
    call get_dust_temperature_from_ptmass(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp,tau=tau)
 endif

end subroutine get_dust_temperature

!-----------------------------------------------------------------------
!+
!  compute dust temperature once the optical depth is known
!+
!-----------------------------------------------------------------------
subroutine get_dust_temperature_from_ptmass(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp,tau,tau_lucy)
 use part,    only:isdead_or_accreted,iLum,iTeff,iReff,itemp
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:),xyzmh_ptmass(:,:),eos_vars(:,:)
 real,     intent(inout), optional :: tau(:), tau_lucy(:)
 real,     intent(out)   :: dust_temp(:)
 real                    :: r,L_star,T_star,R_star,xa,ya,za
 integer                 :: i,j

 !
 ! sanity check, return zero if no sink particles or dust flag is off
 !
 if (nptmass < 1) then
    dust_temp = 0.
    return
 endif

 !property of the sink particle
 j = 1
 T_star = xyzmh_ptmass(iTeff,j)
 L_star = xyzmh_ptmass(iLum,j)
 R_star = xyzmh_ptmass(iReff,j) !sqrt(L_star/(4.*pi*steboltz*utime**3/umass*R_star**4))
 xa = xyzmh_ptmass(1,j)
 ya = xyzmh_ptmass(2,j)
 za = xyzmh_ptmass(3,j)
 select case (iget_tdust)
 case (1)
    ! simple T(r) relation
    !$omp parallel  do default(none) &
    !$omp shared(npart,xa,ya,za,R_star,T_star,xyzh,dust_temp,tdust_exp) &
    !$omp private(i,r)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          r = sqrt((xyzh(1,i)-xa)**2 + (xyzh(2,i)-ya)**2 + (xyzh(3,i)-za)**2)
          dust_temp(i) = T_star*(R_star/r)**tdust_exp
       endif
    enddo
    !$omp end parallel do
 case(2)
    ! Flux dilution without attenuation
    !$omp parallel  do default(none) &
    !$omp shared(npart,xa,ya,za,R_star,T_star,xyzh,dust_temp,tdust_exp) &
    !$omp private(i,r)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          r = sqrt((xyzh(1,i)-xa)**2 + (xyzh(2,i)-ya)**2 + (xyzh(3,i)-za)**2)
          if (r > R_star) dust_temp(i) = T_star * (.5*(1.-sqrt(1.-(R_star/r)**2)))**(1./4.)
       endif
    enddo
    !$omp end parallel do
 case(3)
    ! Flux dilution with attenuation (exp(-tau))
    !$omp parallel  do default(none) &
    !$omp shared(npart,xa,ya,za,R_star,T_star,xyzh,dust_temp,tdust_exp,tau) &
    !$omp private(i,r)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          r = sqrt((xyzh(1,i)-xa)**2 + (xyzh(2,i)-ya)**2 + (xyzh(3,i)-za)**2)
          if (r > R_star) dust_temp(i) = T_star * (.5*(1.-sqrt(1.-(R_star/r)**2))*exp(-tau(i)))**(1./4.)
       endif
    enddo
    !$omp end parallel do
 case(4)
    ! Lucy approximation for Tdust
    !$omp parallel  do default(none) &
    !$omp shared(npart,xa,ya,za,R_star,T_star,xyzh,dust_temp,tdust_exp,tau_lucy) &
    !$omp private(i,r)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          r = sqrt((xyzh(1,i)-xa)**2 + (xyzh(2,i)-ya)**2 + (xyzh(3,i)-za)**2)
          if (r  <  R_star) r = R_star
          if (isnan(tau_lucy(i))) tau_lucy(i) = 2./3.
          dust_temp(i) = T_star * (.5*(1.-sqrt(1.-(R_star/r)**2)+3./2.*tau_lucy(i)))**(1./4.)
       endif
    enddo
    !$omp end parallel do
 case default
    ! sets Tdust = Tgas
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,eos_vars,dust_temp) &
    !$omp private(i)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          dust_temp(i) = eos_vars(itemp,i)
       endif
    enddo
    !$omp end parallel do
 end select

end subroutine get_dust_temperature_from_ptmass

!-----------------------------------------------------------------------
!+
!  write options to input file
!+
!-----------------------------------------------------------------------
subroutine write_options_ptmass_radiation(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(isink_radiation,'isink_radiation','sink radiation pressure method (0=off,1=alpha,2=dust,3=alpha+dust)',iunit)
 if (isink_radiation == 1 .or. isink_radiation == 3) then
    call write_inopt(alpha_rad,'alpha_rad','fraction of the gravitational acceleration imparted to the gas',iunit)
 endif
 if (isink_radiation >= 2) then
    call write_inopt(iget_tdust,'iget_tdust','dust temperature (0:Tdust=Tgas 1:T(r) 2:Flux dilution 3:Attenuation 4:Lucy)',iunit)
    if (iget_tdust /= 2) call write_inopt(iray_resolution,&
                                   'iray_resolution','set the number of rays to 12*4**iray_resolution (deactivated if <0)',iunit)
 endif
 if (iget_tdust == 1) then
    call write_inopt(tdust_exp,'tdust_exp','exponent of the dust temperature profile',iunit)
 endif

end subroutine write_options_ptmass_radiation

!-----------------------------------------------------------------------
!+
!  read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ptmass_radiation(name,valstring,imatch,igotall,ierr)
 use io,  only:fatal
 use dim, only:itau_alloc,itauL_alloc
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 integer :: ni
 character(len=30), parameter :: label = 'read_options_ptmass_radiation'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('alpha_rad')
    read(valstring,*,iostat=ierr) alpha_rad
    ngot = ngot + 1
    if (alpha_rad < 0.) call fatal(label,'invalid setting for alpha_rad (must be >= 0)')
 case('isink_radiation')
    read(valstring,*,iostat=ierr) isink_radiation
    ngot = ngot + 1
    if (isink_radiation < 0 .or. isink_radiation > 3) call fatal(label,'invalid setting for isink_radiation ([0,3])')
 case('iget_tdust')
    read(valstring,*,iostat=ierr) iget_tdust
    ngot = ngot + 1
    if (iget_tdust == 3) itau_alloc  = 1
    if (iget_tdust == 4) itauL_alloc = 1
    if (iget_tdust < 0 .or. iget_tdust > 4) call fatal(label,'invalid setting for iget_tdust ([0,4])')
 case('tdust_exp')
    read(valstring,*,iostat=ierr) tdust_exp
    ngot = ngot + 1
 case('iray_resolution')
    read(valstring,*,iostat=ierr) iray_resolution
    if (iray_resolution >= 0) itau_alloc = 1
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 ni = 1
 if (isink_radiation > 0) then
    ni = ni+1
 endif
 igotall = (ngot >= ni)
 !when Lucy is activated, no need to calculate optical depth
 if (iget_tdust == 4) itau_alloc = 0

end subroutine read_options_ptmass_radiation

end module ptmass_radiation
