!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
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
!   - alpha_rad       : *fraction of the gravitational acceleration imparted to the gas*
!   - iget_tdust      : *method for computing dust temperature (0=none 1=T(r) 2=Lucy 3=MCFOST)*
!   - isink_radiation : *sink radiation pressure method (0=off,1=alpha,2=dust)*
!
! :Dependencies: dust_formation, infile_utils, io, kernel, part, units
!


 implicit none
 integer, public  :: isink_radiation = 0
 integer, public  :: iget_tdust = 0
 real,    public  :: alpha_rad = 0.

 public :: get_rad_accel_from_ptmass,read_options_ptmass_radiation,write_options_ptmass_radiation
 public :: get_dust_temperature_from_ptmass
 private
 integer, parameter :: N = 1024
 real, parameter :: theta = 0., phi = 0.
 real, parameter :: u(3) = (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)

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
 use part,  only:dust_temp,isdead_or_accreted,ilum
 use units, only:umass,unit_energ,utime
#ifdef NUCLEATION
 use part,  only:nucleation
#endif
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:)
 real,     intent(in)    :: xyzmh_ptmass(:,:)
 real,     intent(inout) :: fext(:,:)
 real                    :: dx,dy,dz,xa,ya,za,r,Mstar_cgs,Lstar_cgs,ax,ay,az
 integer                 :: i,j

 do j=1,nptmass
    Mstar_cgs  = xyzmh_ptmass(4,j)*umass
    Lstar_cgs  = xyzmh_ptmass(ilum,j)*unit_energ/utime
    !compute radiative acceleration if sink particle is assigned a non-zero luminosity
    if (Lstar_cgs > 0.d0) then
       xa = xyzmh_ptmass(1,j)
       ya = xyzmh_ptmass(2,j)
       za = xyzmh_ptmass(3,j)
       !$omp parallel  do default(none) &
#ifdef NUCLEATION
       !$omp shared(nucleation)&
#endif
       !$omp shared(dust_temp) &
       !$omp shared(npart,xa,ya,za,Mstar_cgs,Lstar_cgs,xyzh,fext) &
       !$omp private(i,dx,dy,dz,ax,ay,az,r)
       do i=1,npart
          if (.not.isdead_or_accreted(xyzh(4,i))) then
             dx = xyzh(1,i) - xa
             dy = xyzh(2,i) - ya
             dz = xyzh(3,i) - za
             r = sqrt(dx**2 + dy**2 + dz**2)
#ifdef NUCLEATION
             call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,ax,ay,az,kappa=nucleation(8,i))
#else
             call get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,ax,ay,az,Tdust=dust_temp(i))
#endif
             fext(1,i) = fext(1,i) + ax
             fext(2,i) = fext(2,i) + ay
             fext(3,i) = fext(3,i) + az
          endif
       enddo
       !$omp end parallel do
    endif
 enddo

end subroutine get_rad_accel_from_ptmass

!-----------------------------------------------------------------------
!+
!  compute radiative acceleration from a SINGLE sink particle
!  based on sink particle luminosity and computed opacities / column depth
!+
!-----------------------------------------------------------------------
subroutine get_radiative_acceleration_from_star(r,dx,dy,dz,Mstar_cgs,Lstar_cgs,&
     ax,ay,az,Tdust,kappa)
 use units,   only: umass
#ifdef NUCLEATION
 use dust_formation, only:calc_alpha_dust
#else
 use dust_formation, only:calc_alpha_bowen
#endif
 real, intent(in)    :: r,dx,dy,dz,Mstar_cgs,Lstar_cgs
 real, optional, intent(in) :: kappa,Tdust
 real, intent(out)   :: ax,ay,az
 real :: fac,alpha_dust

 select case (isink_radiation)
 case (1)
    ! alpha wind
    fac = alpha_rad*Mstar_cgs/(umass*r**3)
 case (2)
    ! radiation pressure on dust
#ifdef NUCLEATION
    call calc_alpha_dust(Mstar_cgs,Lstar_cgs,kappa,alpha_dust)
#else
    call calc_alpha_bowen(Mstar_cgs,Lstar_cgs,Tdust,alpha_dust)
#endif
    fac = alpha_dust*Mstar_cgs/(umass*r**3)
 case default
    ! no radiation pressure
    fac = 0.
 end select
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
subroutine get_dust_temperature_from_ptmass(npart,xyzh,nptmass,xyzmh_ptmass,dust_temp)
 use part,    only:isdead_or_accreted,iLum,iTeff,iReff
 integer,  intent(in)    :: nptmass,npart
 real,     intent(in)    :: xyzh(:,:),xyzmh_ptmass(:,:)
 real,     intent(out)   :: dust_temp(:)
 real                    :: r,L_star,T_star,R_star,xa,ya,za
 integer                 :: i,j

 !
 ! sanity check, return zero if no sink particles
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
    ! simple T(r) relation
 case (1)
    !$omp parallel  do default(none) &
    !$omp shared(npart,xa,ya,za,R_star,T_star,xyzh,dust_temp) &
    !$omp private(i,r)
    do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          r = sqrt((xyzh(1,i)-xa)**2 + (xyzh(2,i)-ya)**2 + (xyzh(3,i)-za)**2)
          dust_temp(i) = T_star*sqrt(R_star/r)
       endif
    enddo
    !$omp end parallel do
 case(2)
    call get_Teq_from_Lucy(npart,xyzh,xa,ya,za,R_star,T_star,dust_temp)
 end select

end subroutine get_dust_temperature_from_ptmass

!-------------------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature using the Lucy approximation
!  Performs ray-tracing along 1 direction (could be generalized to include other directions)
!+
!-------------------------------------------------------------------------------
subroutine get_Teq_from_Lucy(npart,xyzh,xa,ya,za,R_star,T_star,dust_temp)
 use part,  only:isdead_or_accreted
#ifdef NUCLEATION
 use part,  only:nucleation
#endif
 integer,  intent(in)    :: npart
 real,     intent(in)    :: xyzh(:,:),xa,ya,za,R_star,T_star
 real,     intent(out)   :: dust_temp(:)
 real     :: r(3),r0(3),d,dmin,dmax,d2_axis,OR(N),Teq(N),K3(N),rho_over_r2(2*N+1),rho(N)
 integer  :: i,idx_axis(npart),naxis

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
          !$omp critical
          naxis = naxis+1
          idx_axis(naxis) = i
          !$omp end critical
       endif
    endif
 enddo
!$omp end parallel do
 dmin = sqrt(dmin)
 dmax = sqrt(dmax)


#ifdef NUCLEATION
 call density_along_line(npart, xyzh, r0, naxis, idx_axis, -dmax, dmax, R_star, N, rho, &
      rho_over_r2, dust_temp, Teq, nucleation(5,:), K3)
 call calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq, K3)
#else
 call density_along_line(npart, xyzh, r0, naxis, idx_axis, -dmax, dmax, R_star, N, rho, &
      rho_over_r2, dust_temp, Teq)
 call calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq)
#endif
 call interpolate_on_particles(npart, N, dmax, r0, Teq, dust_temp, xyzh)

end subroutine get_Teq_from_Lucy

!--------------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature along the ray direction
!+
!--------------------------------------------------------------------------
subroutine calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq, K3)
 use dust_formation, only : calc_kappa_dust,kappa_dust_bowen
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
#ifdef NUCLEATION
       if (rho(i) > 0.) then
          call calc_kappa_dust(K3(i),Teq(i),rho(i),kappa(i))
       else
          kappa(i) = 0.d0
       endif
#else
       kappa(i) = kappa_dust_bowen(Teq(i))
#endif
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
 integer, intent(in) :: iunit

 call write_inopt(isink_radiation,'isink_radiation','sink radiation pressure method (0=off,1=alpha,2=dust)',iunit)
 if (isink_radiation == 1) then
    call write_inopt(alpha_rad,'alpha_rad','fraction of the gravitational acceleration imparted to the gas',iunit)
 elseif (isink_radiation == 2) then
    call write_inopt(iget_tdust,'iget_tdust','method for computing dust temperature (0=none 1=T(r) 2=Lucy 3=MCFOST)',iunit)
 endif

end subroutine write_options_ptmass_radiation

!-----------------------------------------------------------------------
!+
!  read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ptmass_radiation(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
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
    if (isink_radiation < 0 .or. isink_radiation > 2) call fatal(label,'invalid setting for isink_radiation ([0,2])')
 case('iget_tdust')
    read(valstring,*,iostat=ierr) iget_tdust
    ngot = ngot + 1
    if (iget_tdust < 0 .or. iget_tdust > 2) call fatal(label,'invalid setting for iget_tdust ([0,2])')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)
 if (isink_radiation > 0) igotall = (ngot >= 2)

end subroutine read_options_ptmass_radiation

end module ptmass_radiation
