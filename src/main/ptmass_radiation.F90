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
! :Dependencies: dim, dust_formation, eos, infile_utils, io, kernel,
!   options, part, physcon, units
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
 use options, only:ieos
 use eos,     only:get_temperature
 use dim,     only:do_nucleation
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
    call get_Teq_from_Lucy(npart,xyzh,xa,ya,za,R_star,T_star,dust_temp)
 case default
    ! sets Tdust = Tgas
    pmassi         = massoftype(igas)
    !$omp parallel  do default(none) &
    !$omp shared(npart,ieos,xyzh,vxyzu,pmassi,dust_temp) &
    !$omp shared(nucleation,do_nucleation) &
    !$omp private(i,vxyzui)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          vxyzui= vxyzu(:,i)
          if (do_nucleation) then
             dust_temp(i) = get_temperature(ieos,xyzh(:,i),rhoh(xyzh(4,i),pmassi),vxyzui,&
                  gammai=nucleation(idgamma,i),mui=nucleation(idmu,i))
          else
             dust_temp(i) = get_temperature(ieos,xyzh(:,i),rhoh(xyzh(4,i),pmassi),vxyzui)
          endif
       endif
    enddo
    !$omp end parallel do
 end select

end subroutine get_dust_temperature_from_ptmass

!-------------------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature using the Lucy approximation
!  Performs ray-tracing along 1 direction (could be generalized to include other directions)
!+
!-------------------------------------------------------------------------------
!          UNDER CONSTRUCTION!!!!!!!!!!!!
!-------------------------------------------------------------------------------
subroutine get_Teq_from_Lucy(npart,xyzh,xa,ya,za,R_star,T_star,dust_temp)
 use part,  only:isdead_or_accreted
 use part,  only:nucleation,idK3
 use dim,   only:do_nucleation
 integer,  intent(in)    :: npart
 real,     intent(in)    :: xyzh(:,:),xa,ya,za,R_star,T_star
 real,     intent(out)   :: dust_temp(:)
 real     :: r(3),r0(3),d,dmin,dmax,d2_axis,OR(N),Teq(N),rho_over_r2(2*N+1),rho(N)
 integer  :: i,idx_axis(npart),naxis
#ifdef DUST_NUCLEATION
 real     :: K3(N)
#endif
 !.. find particles that lie within 2 smoothing lengths of the ray axis
 r0(1:3) = (/xa, ya, za/)
 dmin = huge(dmin)
 dmax = 0
 naxis = 0
!$omp parallel do default(none) &
!$omp shared(npart,xyzh,r0,naxis,idx_axis) &
!$omp private(i,r,d,d2_axis) &
!$omp reduction(min:dmin) &
!$omp reduction(max:dmax)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       r = xyzh(1:3,i)-r0
       !d = r(1)**2+r(2)**2+r(3)**2
       d = dot_product(r,r)
       dmin = min(d,dmin)
       dmax = max(d,dmax)
       !distance to the axis
       !d2_axis = sq_distance_to_z(r)
       d2_axis = sq_distance_to_line(r,u)
       if (d2_axis < 4.*xyzh(4,i)*xyzh(4,i)) then
          !$omp critical (naxis_add)
          naxis = naxis+1
          idx_axis(naxis) = i
          !$omp end critical (naxis_add)
       endif
    endif
 enddo
!$omp end parallel do
 dmin = sqrt(dmin)
 dmax = sqrt(dmax)


#ifdef DUST_NUCLEATION
 print *,'WARNING : CHECK THAT idust_opacity > 0'
 if (do_nucleation) then
    call density_along_line(npart, xyzh, r0, naxis, idx_axis, -dmax, dmax, R_star, N, rho, &
         rho_over_r2, dust_temp, Teq, nucleation(idK3,:), K3)
    call calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq, K3)
 else
    call density_along_line(npart, xyzh, r0, naxis, idx_axis, -dmax, dmax, R_star, N, rho, &
         rho_over_r2, dust_temp, Teq)
    call calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq)
 endif
#endif
 call interpolate_on_particles(npart, N, dmax, r0, Teq, dust_temp, xyzh)

end subroutine get_Teq_from_Lucy

!--------------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature along the ray direction
!+
!--------------------------------------------------------------------------
subroutine calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq, K3)
 use dim,            only:do_nucleation
 use dust_formation, only:calc_kappa_dust,calc_kappa_bowen
 integer, intent(in)  :: N
 real,    intent(in)  :: dmax, R_star, T_star, rho(N), rho_over_r2(2*N+1)
 real,    optional, intent(in) :: K3(N)
 real,    intent(out) :: Teq(N)

 real :: OR(N),tau_prime(N),vTeq(N),kappa(N),dTeq,rho_m
 real :: dr, fact
 real, parameter :: tol = 1.e-2, kap_gas = 2.e-4
 integer :: i,istart,iter


 tau_prime = 0.
 iter = 0
 vTeq = 0.
 dTeq = 1.
 dr = dmax/N
 forall(i=1:N) OR(i) = i*dr
 OR(N) = dmax
 fact = dr/2. * R_star**2
 do i = 1,N
    if (OR(i) > R_star) exit
 enddo
 istart = i-1
 if (istart > 0) Teq(1:istart) = T_star
 Teq(istart+1:N) =  T_star*(0.5*(1.-sqrt(1.-(R_star/OR(istart+1:N))**2)))
 vTeq = Teq

 do while (dTeq > tol .and. iter < 20)
    if (iter == 0) dTeq = 0.
    iter = iter+1
    do i=N-1,istart+1,-1
       if (do_nucleation) then
          if (rho(i) > 0.) then
             kappa(i) = calc_kappa_dust(K3(i),Teq(i),rho(i))
          else
             kappa(i) = 0.d0
          endif
       else
          kappa(i) = calc_kappa_bowen(Teq(i))
       endif
       rho_m = (rho_over_r2(N-i)+rho_over_r2(N-i+1)+rho_over_r2(N+i+1)+rho_over_r2(N+i+2))
       tau_prime(i) = tau_prime(i+1) + fact*(kappa(i)+kap_gas)*rho_m
       Teq(i) = T_star*(0.5*(1.-sqrt(1.-(R_star/OR(i))**2)) + 0.75*tau_prime(i))**0.25
       dTeq = max(dTeq,abs(1.-Teq(i)/(1.e-5+vTeq(i))))
       vTeq(i) = Teq(i)
    enddo
 enddo

end subroutine calculate_Teq

!-----------------------------------------------------------------------
!+
!  compute the mean properties along the ray
!+
!-----------------------------------------------------------------------
subroutine density_along_line(npart, xyzh, r0, npart_axis, idx_axis, rmin, rmax, r_star, N, &
     rho_cgs, rho_over_r2, T, Teq, K3, K3i)
 use kernel, only:cnormk,wkern
 use part,   only:massoftype,igas,rhoh
 use units,  only:unit_density
 integer, intent(in)  :: npart,N
 real,    intent(in)  :: xyzh(:,:), T(:), r0(3)
 real, optional, intent(in)  :: K3(:)
 integer, intent(in)  ::  npart_axis, idx_axis(npart)
 real,    intent(in)  :: rmin, rmax, R_star
 real,    intent(out) :: rho_over_r2(2*N+1), Teq(N), rho_cgs(N)
 real, optional, intent(out) :: K3i(N)
 real :: rhoi(2*N+1), OR(2*N+1), Ti(2*N+1), Ki(2*N+1), xnorm(2*N+1)
 real :: OH, d2_axis, HR, q2, q, fact0, fact, h, h2, part_mass
 real :: delta_r, rmin_o, rmin_p, rmax_p, dr, r(3), xfact, rhoinv
 integer :: i, np, j, j_min, j_max, Nr

! Discretization of the line of sight in N segments
 Nr = 2*N+1
 dr = (rmax-rmin)/(Nr-1)
 rmin_o = rmin - dr
 do i=1,Nr
    OR(i) = dr*i+rmin_o
 enddo

 rhoi(:) = 0.
 Teq(:) = 0.
 K3i(:) = 0.
 Ki(:) = 0.
 Ti(:)  = 0.
 xnorm(:) = 0.
 part_mass = massoftype(igas)
 fact0 =  part_mass*cnormk
 do i = 1, npart_axis
    np = idx_axis(i)
    r  = xyzh(1:3,np)-r0(3)
    !distance to z-axis
    !OH = r(3)
    !d2_axis = sq_distance_to_z(r)
    OH = dot_product(r,u)
    d2_axis = sq_distance_to_line(r,u)
    h = xyzh(4,np)
    h2 = h*h
    delta_r = sqrt(4.*h2 - d2_axis)
    ! rmin_p and rmax_p are the positions on the line of the two intersections between the line and the interaction sphere
    rmin_p = OH-delta_r
    rmax_p = OH+delta_r
    j_min = ceiling((rmin_p-rmin_o)/dr)
    j_max = floor((rmax_p-rmin_o)/dr)
    j_min = max(1, j_min)
    j_max = min(Nr, j_max)
    ! Adds the contribution of particle np to density at all the discretized locations in the interaction sphere
    fact = fact0/h**3
    rhoinv = 1./rhoh(h,part_mass)
    do j=j_min, j_max
       HR = OR(j) - OH
       q2 = (d2_axis+HR**2)/h2
       q  = sqrt(q2)
       xfact = fact*wkern(q2,q)
       rhoi(j)  = rhoi(j) + xfact
       xnorm(j) = xnorm(j)+xfact*rhoinv
       Ti(j)    = Ti(j)  + xfact*rhoinv*T(np)
       if (present(K3)) Ki(j) = Ki(j) + xfact*rhoinv*K3(np)
    enddo
 enddo

 ! rho_over_r2 = 0 inside the star so that we do not divide by zero!
 do j=1,Nr
    if (xnorm(j) > 0.) then
       Ti(j)  = Ti(j)/xnorm(j)
       if (present(K3)) Ki(j)  = Ki(j) /xnorm(j)
    endif
    if (abs(OR(j))  <  r_star) then
       rho_over_r2(j) = 0.
    else
       rho_over_r2(j) = rhoi(j)/OR(j)**2
    endif
 enddo
 do j=1,N
    rho_cgs(N+1-j) = (rhoi(j)+rhoi(2*N-j+2))*unit_density/2.
    Teq(N+1-j) = (Ti(j)+Ti(2*N-j+2))/2.
    if (present(K3)) K3i(N+1-j) = (Ki(j)+Ki(2*N-j+2))/2.
 enddo

end subroutine density_along_line

!-----------------------------------------------------------------------
!+
!  Interpolates a quantity computed on the discretized line of sight for all SPH particles
!  (spherical symmetry assumed)
!+
!-----------------------------------------------------------------------
subroutine interpolate_on_particles(npart, N, dmax, r0, Teq, dust_temp, xyzh)
 use part,    only:isdead_or_accreted
 integer, intent(in)  :: npart, N
 real,    intent(in)  :: dmax, r0(3), Teq(N), xyzh(:,:)
 real,    intent(out) :: dust_temp(:)

 real :: r(3), d, dr, d2
 integer :: i, j

 dr = dmax / N
 !should start at nwall
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       r = xyzh(1:3,i) - r0
       d2 = dot_product(r,r)
       d = sqrt(d2)
       j = min(int(d/dr),N-1)
       dust_temp(i) = (d-dr*j)*(Teq(j+1)-Teq(j))/dr + Teq(j)
    endif
 enddo
end subroutine interpolate_on_particles

real function sq_distance_to_z(r)
 real, intent(in) :: r(3)
 sq_distance_to_z = r(1)*r(1)+r(2)*r(2)
end function sq_distance_to_z

real function sq_distance_to_line(r,u)
 real, intent(in) :: r(3),u(3)
 real :: p,d(3)
 p = dot_product(r,u)
 d = r-p*u
 sq_distance_to_line = dot_product(d,d)
end function sq_distance_to_line

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
