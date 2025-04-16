!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for dusty wind testing
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: dim, dust_formation, kernel, part, units
!

 implicit none
 character(len=20), parameter, public :: analysistype = 'dustywind'

 public                               :: do_analysis

 private
 integer, parameter :: N = 1024 !32
 real, parameter :: theta = 0., phi = 0.
 real, parameter :: u(3) = (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use part,  only: nptmass,xyzmh_ptmass,vxyz_ptmass,iLum,iTeff,iReff
 use part,  only: dust_temp,isdead_or_accreted,nucleation
 use dust_formation, only: set_abundances

 !general variables
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 real    :: L_star,T_star,R_star,xa,ya,za
 integer :: j

 call set_abundances
 !property of the sink particle
 j = 1
 T_star = xyzmh_ptmass(iTeff,j)
 L_star = xyzmh_ptmass(iLum,j)
 R_star = xyzmh_ptmass(iReff,j) !sqrt(L_star/(4.*pi*steboltz*utime**3/umass*R_star**4))
 xa = xyzmh_ptmass(1,j)
 ya = xyzmh_ptmass(2,j)
 za = xyzmh_ptmass(3,j)
 call get_Teq_from_Lucy(npart,xyzh,xa,ya,za,R_star,T_star,dust_temp)


end subroutine do_analysis

!-------------------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature using the Lucy approximation
!  Performs ray-tracing along 1 direction (could be generalized to include other directions)
!+
!-------------------------------------------------------------------------------
subroutine get_Teq_from_Lucy(npart,xyzh,xa,ya,za,R_star,T_star,dust_temp)
 use part,  only:isdead_or_accreted,nucleation,idK3
 use dim,   only:do_nucleation
 integer,  intent(in)    :: npart
 real,     intent(in)    :: xyzh(:,:),xa,ya,za,R_star,T_star
 real,     intent(out)   :: dust_temp(:)
 real     :: r(3),r0(3),d,dmin,dmax,d2_axis,OR(N),Teq(N),K3(N),rho_over_r2(2*N+1),rho(N)
 integer  :: i,idx_axis(npart),naxis

 !.. find particles that lie within 2 smoothing lengths of the ray axis
 r0(1:3) = (/xa, ya, za/)
 dmin = 1.d99
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
          !$omp critical (crit_naxis_add)
          naxis = naxis+1
          idx_axis(naxis) = i
          !$omp end critical (crit_naxis_add)
       endif
    endif
 enddo
!$omp end parallel do
 dmin = sqrt(dmin)
 dmax = sqrt(dmax)


 if (do_nucleation) then
    call density_along_line(npart, xyzh, r0, naxis, idx_axis, -dmax, dmax, R_star, N, rho, &
         rho_over_r2, dust_temp, Teq, nucleation(idK3,:), K3)
    call calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq, K3)
 else
    call density_along_line(npart, xyzh, r0, naxis, idx_axis, -dmax, dmax, R_star, N, rho, &
         rho_over_r2, dust_temp, Teq)
    call calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq)
 endif
 call interpolate_on_particles(npart, N, dmax, r0, Teq, dust_temp, xyzh)

end subroutine get_Teq_from_Lucy

!--------------------------------------------------------------------------
!+
!  Calculates the radiative equilibrium temperature along the ray direction
!+
!--------------------------------------------------------------------------
subroutine calculate_Teq(N, dmax, R_star, T_star, rho, rho_over_r2, OR, Teq, K3)
 use dust_formation, only : calc_kappa_dust,calc_kappa_bowen,idust_opacity
 integer, intent(in)  :: N
 real,    intent(in)  :: dmax, R_star, T_star, rho(N), rho_over_r2(2*N+1)
 real,    optional, intent(in) :: K3(N)
 real,    intent(out) :: Teq(N)

 real :: OR(N),tau_prime(N),vTeq(N),kappa(N),dTeq,pTeq(N)
 real :: dr, fact, rho_on_r2(N)
 real, parameter :: tol = 1.d-2, kap_gas = 2.d-4
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
 pTeq= Teq
 rho_on_r2 = 0.
 kappa = 0.

 do while (dTeq > tol .and. iter < 20)
    if (iter == 0) dTeq = 0.
    iter = iter+1
    do i=N-1,istart+1,-1
       if (idust_opacity == 2) then
          if (rho(i) > 0.) then
             kappa(i) = calc_kappa_dust(K3(i),Teq(i),rho(i))
          else
             kappa(i) = 0.d0
          endif
       elseif (idust_opacity == 1) then
          kappa(i) = calc_kappa_bowen(Teq(i))
       endif
       rho_on_r2(i) = rho_over_r2(N-i)+rho_over_r2(N-i+1)+rho_over_r2(N+i+1)+rho_over_r2(N+i+2)
       !if (iter >= 1) print *,'teq loop',i,K3(i),Teq(i),kappa(i),rho_on_r2(i)
       tau_prime(i) = tau_prime(i+1) + fact*(kappa(i)+kap_gas) *rho_on_r2(i)

       Teq(i) = T_star*(0.5*(1.-sqrt(1.-(R_star/OR(i))**2)) + 0.75*tau_prime(i))**(1./4.)
       dTeq = max(dTeq,abs(1.-Teq(i)/(1.d-5+vTeq(i))))
       vTeq(i) = Teq(i)
    enddo
    print *,iter,dTeq
 enddo
 print *,iter
 open(unit=220,file='Teq.dat')
 write(220,*) '# ng z vTeq Teq tau kappa rho_on_r2'
 do i = 1,N
    write(220,*) i,OR(i),pTeq(i),Teq(i),tau_prime(i),kappa(i),rho_on_r2(i)
 enddo
 close(220)

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
    print *,i,OR(i),R_star
 enddo
 print *,'*******',rmax,rmin,r_star

 open(unit=220,file='allpart.dat')
 write(220,*) '# ng x y z rho T K'
 do i = 1, npart
    write (220,*) np,xyzh(1:3,i)-r0(3),rhoh(xyzh(4,i),part_mass),T(i),K3(i)
 enddo
 close(220)
 rhoi(:) = 0.
 Teq(:) = 0.
 K3i(:) = 0.
 Ki(:) = 0.
 Ti(:)  = 0.
 xnorm(:) = 0.
 part_mass = massoftype(igas)
 fact0 =  part_mass*cnormk
 open(unit=221,file='part_axis.dat')
 write(221,*) '# ng x y z rho T K'
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
    write (221,*) np,r,rhoh(h,part_mass),T(np),K3(np)
    do j=j_min, j_max
       HR = OR(j) - OH
       q2 = (d2_axis+HR**2)/h2
       q  = sqrt(q2)
       xfact = fact*wkern(q2,q)
       rhoi(j)  = rhoi(j) + xfact
       xnorm(j) = xnorm(j)+xfact*rhoinv
       Ti(j)    = Ti(j)  + xfact*rhoinv*T(np)
       if (present(K3)) Ki(j) = Ki(j) + xfact*rhoinv*K3(np)
       !print *,j,Ti(j),T(np),part_mass/(rhoh(h,part_mass)*h**3)!rhoh(h,part_mass),part_mass,q,fact,wkern(q2,q)
    enddo
 enddo
 close (221)
! rho_over_r2 = 0 inside the star so that we do not divide by zero!
 open(unit=222,file='ray.dat')
 write(222,*) '# ng z rho T K xnorm rho_over_r2'
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
    print *,j,rho_over_r2(j)
    write (222,*) j,OR(j),rhoi(j),Ti(j),Ki(j),xnorm(j),rho_over_r2(j)
 enddo
 close(222)
 do j=1,N
    rho_cgs(N+1-j) = (rhoi(j)+rhoi(2*N-j+2))*unit_density/2.
    Teq(N+1-j) = (Ti(j)+Ti(2*N-j+2))/2.
    if (present(K3)) K3i(N+1-j) = (Ki(j)+Ki(2*N-j+2))/2.
!    print *,'k3i',j,k3i(j)
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
 open(unit=220,file='all_final.dat')
 write(220,*) '# ng x y z T'
 do i = 1, npart
    write (220,*) i,xyzh(1:3,i)-r0(3),dust_temp(i)
 enddo
 close(220)
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

end module analysis
