!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!   - Lstar           : *Stellar luminosity (for radiation pressure, in Lsun)*
!   - Mstar           : *Stellar mass (in Msun)*
!   - alpha_rad       : *fraction of the gravitational acceleration imparted to the gas*
!   - iget_tdust      : *dust temperature (0:Tdust=Tgas 1:T(r) 2:Lucy (devel))*
!   - isink_radiation : *sink radiation pressure method (0=off,1=alpha,2=dust,3=alpha+dust)*
!   - tdust_exp       : *exponent of the dust temperature profile*
!
! :Dependencies: dim, dust_formation, eos, infile_utils, io, kernel, part,
!   physcon, units
!


 implicit none
 integer, public  :: isink_radiation = 0
 integer, public  :: iget_tdust      = 0
 real,    public  :: tdust_exp       = 0.5
 real,    public  :: alpha_rad       = 0.

 public :: get_rad_accel_from_ptmass,read_options_ptmass_radiation,write_options_ptmass_radiation
 public :: get_dust_temperature_from_ptmass
 public :: init_radiation_ptmass

 integer, parameter :: N = 1024
 real, parameter :: theta = 0., phi = 0.
 real, parameter :: u(3) = (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)

 private

 real  :: Lstar_lsun      = 5000.
 real  :: Mstar_msun      = 1.

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
subroutine get_rad_accel_from_ptmass(nptmass,npart,xyzh,xyzmh_ptmass,fext)
 use part,    only:ilum
 use units,   only:umass,unit_energ,utime
 use dim,     only:star_radiation
 use physcon, only:solarl,solarm
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:)
 real,     intent(in)    :: xyzmh_ptmass(:,:)
 real,     intent(inout) :: fext(:,:)
 real                    :: xa,ya,za,Mstar_cgs,Lstar_cgs
 integer                 :: j

 if (star_radiation) then
    Lstar_cgs  = Lstar_lsun*solarl
    Mstar_cgs  = Mstar_msun*solarm
    if (Lstar_cgs > 0.d0) then
       xa = xyzmh_ptmass(1,1)
       ya = xyzmh_ptmass(2,1)
       za = xyzmh_ptmass(3,1)
       call calc_rad_accel_from_ptmass(npart,xa,ya,za,Lstar_cgs,Mstar_cgs,xyzh,fext)
    endif
 else
    do j=1,nptmass
       if (xyzmh_ptmass(4,j) < 0.) cycle
       Mstar_cgs  = xyzmh_ptmass(4,j)*umass
       Lstar_cgs  = xyzmh_ptmass(ilum,j)*unit_energ/utime
       !compute radiative acceleration if sink particle is assigned a non-zero luminosity
       if (Lstar_cgs > 0.d0) then
          xa = xyzmh_ptmass(1,j)
          ya = xyzmh_ptmass(2,j)
          za = xyzmh_ptmass(3,j)
          call calc_rad_accel_from_ptmass(npart,xa,ya,za,Lstar_cgs,Mstar_cgs,xyzh,fext)
       endif
    enddo
 endif

end subroutine get_rad_accel_from_ptmass

!-----------------------------------------------------------------------
!+
!  compute radiative acceleration on all particles
!+
!-----------------------------------------------------------------------
subroutine calc_rad_accel_from_ptmass(npart,xa,ya,za,Lstar_cgs,Mstar_cgs,xyzh,fext)
 use part,  only:isdead_or_accreted,dust_temp,nucleation,idkappa,idalpha
 use dim,   only:do_nucleation
 use dust_formation, only:calc_kappa_bowen
 integer,  intent(in)    :: npart
 real,     intent(in)    :: xyzh(:,:)
 real,     intent(in)    :: xa,ya,za,Lstar_cgs,Mstar_cgs
 real,     intent(inout) :: fext(:,:)
 real                    :: dx,dy,dz,r,ax,ay,az,alpha,kappa
 integer                 :: i

 !$omp parallel  do default(none) &
 !$omp shared(nucleation,do_nucleation)&
 !$omp shared(dust_temp) &
 !$omp shared(npart,xa,ya,za,Mstar_cgs,Lstar_cgs,xyzh,fext) &
 !$omp private(i,dx,dy,dz,ax,ay,az,r,alpha,kappa)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       dx = xyzh(1,i) - xa
       dy = xyzh(2,i) - ya
       dz = xyzh(3,i) - za
       r = sqrt(dx**2 + dy**2 + dz**2)
       if (do_nucleation) then
          call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
               nucleation(idkappa,i),ax,ay,az,nucleation(idalpha,i))
       else
          kappa = calc_kappa_bowen(dust_temp(i))
          call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
               kappa,ax,ay,az,alpha)
       endif
       fext(1,i) = fext(1,i) + ax
       fext(2,i) = fext(2,i) + ay
       fext(3,i) = fext(3,i) + az
    endif
 enddo
 !$omp end parallel do
end subroutine calc_rad_accel_from_ptmass


!-----------------------------------------------------------------------
!+
!  compute radiative acceleration from a SINGLE sink particle
!  based on sink particle luminosity and computed opacities / column depth
!+
!-----------------------------------------------------------------------
subroutine get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
     kappa,ax,ay,az,alpha)
 use units,          only:umass
 use dust_formation, only:calc_Eddington_factor
 real, intent(in)    :: r,dx,dy,dz,Mstar_cgs,Lstar_cgs,kappa
 real, intent(out)   :: ax,ay,az,alpha
 real :: fac

 select case (isink_radiation)
 case (1)
    ! alpha wind
    alpha = alpha_rad
 case (2)
    ! radiation pressure on dust
    alpha = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, kappa)
 case (3)
    ! radiation pressure on dust + alpha_rad (=1+2)
    alpha = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, kappa) + alpha_rad
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
subroutine get_dust_temperature_from_ptmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,dust_temp)
 use part,    only:isdead_or_accreted,iLum,iTeff,iReff,rhoh,massoftype,igas,nucleation,idmu,idgamma
 use part,    only:eos_vars,itemp
 use eos,     only:ieos,get_temperature
 use dim,     only:do_nucleation
 use io,      only:fatal
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:),xyzmh_ptmass(:,:),vxyzu(:,:)
 real,     intent(out)   :: dust_temp(:)
 real                    :: r,L_star,T_star,R_star,xa,ya,za,pmassi,vxyzui(4)
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
    !Lucy approximation for Tdust
    print *,'Not implemented yet'
    call fatal('ptmass_radiation','Lucy approximation not implemented')
 case default
    ! sets Tdust = Tgas
    pmassi = massoftype(igas)
    !$omp parallel  do default(none) &
    !$omp shared(npart,ieos,xyzh,vxyzu,eos_vars,pmassi,dust_temp) &
    !$omp shared(nucleation,do_nucleation) &
    !$omp private(i,vxyzui)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          vxyzui= vxyzu(:,i)
          if (do_nucleation) then
             dust_temp(i) = get_temperature(ieos,xyzh(:,i),rhoh(xyzh(4,i),pmassi),vxyzui,&
                  gammai=nucleation(idgamma,i),mui=nucleation(idmu,i))

          else
             dust_temp(i) = eos_vars(itemp,i)
          endif
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
 use infile_utils, only: write_inopt
 use dim,          only: star_radiation!,store_dust_temperature
 integer, intent(in) :: iunit

 call write_inopt(isink_radiation,'isink_radiation','sink radiation pressure method (0=off,1=alpha,2=dust,3=alpha+dust)',iunit)
 if (isink_radiation == 1 .or. isink_radiation == 3) then
    call write_inopt(alpha_rad,'alpha_rad','fraction of the gravitational acceleration imparted to the gas',iunit)
 endif
 if (isink_radiation >= 2) then
    if (star_radiation) then
       call write_inopt(Lstar_lsun,'Lstar','Stellar luminosity (for radiation pressure, in Lsun)',iunit)
       call write_inopt(Mstar_msun,'Mstar','Stellar mass (in Msun)',iunit)
    endif
    call write_inopt(iget_tdust,'iget_tdust','dust temperature (0:Tdust=Tgas 1:T(r) 2:Lucy (devel))',iunit)
 endif
 if (iget_tdust == 1 ) then
    call write_inopt(tdust_exp,'tdust_exp','exponent of the dust temperature profile',iunit)
 endif
 !if (iget_tdust > 0) store_dust_temperature = .true.

end subroutine write_options_ptmass_radiation

!-----------------------------------------------------------------------
!+
!  read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ptmass_radiation(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use dim,     only:star_radiation
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 integer :: ni
 character(len=30), parameter :: label = 'read_options_ptmass_radiation'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('Lstar')
    read(valstring,*,iostat=ierr) Lstar_lsun
    ngot = ngot + 1
    if (Lstar_lsun < 0.) call fatal(label,'invalid setting for Lstar (must be >= 0)')
 case('Mstar')
    read(valstring,*,iostat=ierr) Mstar_msun
    ngot = ngot + 1
    if (Mstar_msun < 0.) call fatal(label,'invalid setting for Mstar (must be >= 0)')
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
    if (iget_tdust < 0 .or. iget_tdust > 2) call fatal(label,'invalid setting for iget_tdust ([0,3])')
 case('tdust_exp')
    read(valstring,*,iostat=ierr) tdust_exp
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 ni = 1
 if (isink_radiation > 0) then
    ni = ni+1
    if (star_radiation)  ni = ni+2
 endif
 igotall = (ngot >= ni)

end subroutine read_options_ptmass_radiation

end module ptmass_radiation
