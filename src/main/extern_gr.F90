!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_gr
!
! None
!
! :References: None
!
! :Owner: Spencer Magnall
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, metric_tools, part, physcon, timestep, utils_gr
!
 implicit none

 public :: get_grforce, get_grforce_all, update_grforce_leapfrog, get_tmunu_all, get_tmunu_all_exact, get_tmunu

 private

contains

!---------------------------------------------------------------
!+
!  Wrapper subroutine for computing the force due to spacetime curvature
!  (This may be useful in the future if there is something that indicates
!   whether a particle is gas or test particle.)
!+
!---------------------------------------------------------------
subroutine get_grforce(xyzhi,metrici,metricderivsi,veli,densi,ui,pi,fexti,dtf)
 use io, only:iprint,fatal,error
 real, intent(in)  :: xyzhi(4),metrici(:,:,:),metricderivsi(0:3,0:3,3),veli(3),densi,ui,pi
 real, intent(out) :: fexti(3)
 real, intent(out), optional :: dtf
 integer :: ierr

 call forcegr(xyzhi(1:3),metrici,metricderivsi,veli,densi,ui,pi,fexti,ierr)
 if (ierr > 0) then
    write(iprint,*) 'x,y,z = ',xyzhi(1:3)
    call error('get_u0 in extern_gr','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
    call fatal('get_grforce','could not compute forcegr at r = ',val=sqrt(dot_product(xyzhi(1:3),xyzhi(1:3))) )
 endif

 if (present(dtf)) call dt_grforce(xyzhi,fexti,dtf)

end subroutine get_grforce

subroutine get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,dens,fext,dtexternal)
 use timestep, only:C_force
 use eos,      only:ieos,get_pressure
 use part,     only:isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:), metrics(:,:,:,:), metricderivs(:,:,:,:), dens(:)
 real, intent(inout) :: vxyzu(:,:)
 real, intent(out)   :: fext(:,:), dtexternal
 integer :: i
 real    :: dtf,pi

 dtexternal = huge(dtexternal)

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,metrics,metricderivs,vxyzu,dens,fext,ieos,C_force) &
 !$omp private(i,dtf,pi) &
 !$omp reduction(min:dtexternal)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       pi = get_pressure(ieos,xyzh(:,i),dens(i),vxyzu(:,i))
       call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i),pi,fext(1:3,i),dtf)
       dtexternal = min(dtexternal,C_force*dtf)
    endif
 enddo
 !$omp end parallel do

end subroutine get_grforce_all

!--- Subroutine to calculate the timestep constraint from the 'external force'
!    this is multiplied by the safety factor C_force elsewhere
subroutine dt_grforce(xyzh,fext,dtf)
 use physcon,      only:pi
 use metric_tools, only:imetric,imet_schwarzschild,imet_kerr
 real, intent(in)  :: xyzh(4),fext(3)
 real, intent(out) :: dtf
 real :: r,r2,dtf1,dtf2,f2i
 integer, parameter :: steps_per_orbit = 100

 f2i = fext(1)*fext(1) + fext(2)*fext(2) + fext(3)*fext(3)
 if (f2i > 0.) then
    dtf1 = sqrt(xyzh(4)/sqrt(f2i)) ! This is not really accurate since fi is a component of dp/dt, not da/dt
 else
    dtf1 = huge(dtf1)
 endif

 select case (imetric)
 case (imet_schwarzschild,imet_kerr)
    r2   = xyzh(1)*xyzh(1) + xyzh(2)*xyzh(2) + xyzh(3)*xyzh(3)
    r    = sqrt(r2)
    dtf2 = (2.*pi*sqrt(r*r2))/steps_per_orbit
 case default
    dtf2 = huge(dtf2)
 end select

 dtf = min(dtf1,dtf2)

end subroutine dt_grforce


!----------------------------------------------------------------
!+
!  Compute the source terms required on the right hand side of
!  the relativistic momentum equation. These are of the form:
!   T^\mu\nu dg_\mu\nu/dx^i
!+
!----------------------------------------------------------------
pure subroutine forcegr(x,metrici,metricderivsi,v,dens,u,p,fterm,ierr)
 use metric_tools, only:unpack_metric
 use utils_gr,     only:get_u0
 real,    intent(in)  :: x(3),metrici(:,:,:),metricderivsi(0:3,0:3,3),v(3),dens,u,p
 real,    intent(out) :: fterm(3)
 integer, intent(out) :: ierr
 real    :: gcov(0:3,0:3), gcon(0:3,0:3)
 real    :: v4(0:3), term(0:3,0:3)
 real    :: enth, uzero
 integer :: i

 call unpack_metric(metrici,gcov=gcov,gcon=gcon)

 enth = 1. + u + p/dens

 ! lower-case 4-velocity
 v4(0) = 1.
 v4(1:3) = v(:)

 ! first component of the upper-case 4-velocity
 call get_u0(gcov,v,uzero,ierr)

 ! energy-momentum tensor times sqrtg on 2rho*
 do i=0,3
    term(0:3,i) =  0.5*(enth*uzero*v4(0:3)*v4(i) + P*gcon(0:3,i)/(dens*uzero))
 enddo

 ! source term
 fterm = 0.
 do i=0,3
    fterm(1) =  fterm(1) + dot_product(term(:,i),metricderivsi(:,i,1))
    fterm(2) =  fterm(2) + dot_product(term(:,i),metricderivsi(:,i,2))
    fterm(3) =  fterm(3) + dot_product(term(:,i),metricderivsi(:,i,3))
 enddo

end subroutine forcegr


!-------- I don't think this is actually being used at the moment....
subroutine update_grforce_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,densi,ui,pi)
 use io,             only:fatal
 real, intent(in)    :: dt,xi,yi,zi
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(inout) :: fexti(3)
 real, intent(in)    :: densi,ui,pi
!  real                :: fextv(3)
!  real                :: v1x, v1y, v1z, v1xold, v1yold, v1zold, vhalf2, erri, dton2
!  logical             :: converged
!  integer             :: its, itsmax
!  integer, parameter  :: maxitsext = 50 ! maximum number of iterations on external force
!  real, parameter :: tolv = 1.e-2
!  real, parameter :: tolv2 = tolv*tolv
!  real,dimension(3) :: pos,vel
!  real :: dtf
!
!  itsmax = maxitsext
!  its = 0
!  converged = .false.
!  dton2 = 0.5*dt
!
!  v1x = vhalfx
!  v1y = vhalfy
!  v1z = vhalfz
!  vhalf2 = vhalfx*vhalfx + vhalfy*vhalfy + vhalfz*vhalfz
!  fextv = 0. ! to avoid compiler warning
!
!  iterations : do while (its < itsmax .and. .not.converged)
!     its = its + 1
!     erri = 0.
!     v1xold = v1x
!     v1yold = v1y
!     v1zold = v1z
!     pos = (/xi,yi,zi/)
!     vel = (/v1x,v1y,v1z/)
!     call get_grforce(pos,vel,densi,ui,pi,fextv,dtf)
! !    xi = pos(1)
! !    yi = pos(2)
! !    zi = pos(3)
!     v1x = vel(1)
!     v1y = vel(2)
!     v1z = vel(3)
!
!     v1x = vhalfx + dton2*(fxi + fextv(1))
!     v1y = vhalfy + dton2*(fyi + fextv(2))
!     v1z = vhalfz + dton2*(fzi + fextv(3))
!
!     erri = (v1x - v1xold)**2 + (v1y - v1yold)**2 + (v1z - v1zold)**2
!     erri = erri / vhalf2
!     converged = (erri < tolv2)
!
!  enddo iterations
!
!  if (its >= maxitsext) call fatal('update_grforce_leapfrog','VELOCITY ITERATIONS ON EXTERNAL FORCE NOT CONVERGED!!')
!
!  fexti(1) = fextv(1)
!  fexti(2) = fextv(2)
!  fexti(3) = fextv(3)
!
!  fxi = fxi + fexti(1)
!  fyi = fyi + fexti(2)
!  fzi = fzi + fexti(3)

end subroutine update_grforce_leapfrog

subroutine get_tmunu_all(npart,xyzh,metrics,vxyzu,metricderivs,dens,tmunus)
 use eos,         only:ieos,get_pressure
 use part,        only:isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:), metrics(:,:,:,:), metricderivs(:,:,:,:), dens(:)
 real, intent(inout) :: vxyzu(:,:),tmunus(:,:,:)
 real                :: pi
 integer             :: i
 logical             :: verbose

 verbose = .false.
 ! TODO write openmp parallel code
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,metrics,vxyzu,dens,ieos,tmunus) &
 !$omp private(i,pi,verbose)
 do i=1, npart
    !print*, "i: ", i
    if (i==1) then
       verbose = .true.
    else
       verbose = .false.
    endif
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       pi = get_pressure(ieos,xyzh(:,i),dens(i),vxyzu(:,i))
       call get_tmunu(xyzh(:,i),metrics(:,:,:,i),&
         vxyzu(1:3,i),dens(i),vxyzu(4,i),pi,tmunus(:,:,i),verbose)
    endif
 enddo
 !$omp end parallel do
 !print*, "tmunu calc val is: ", tmunus(0,0,5)
end subroutine get_tmunu_all

subroutine get_tmunu_all_exact(npart,xyzh,metrics,vxyzu,metricderivs,dens,tmunus)
 use eos,         only:ieos,get_pressure
 use part,        only:isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:), metrics(:,:,:,:), metricderivs(:,:,:,:), dens(:)
 real, intent(inout) :: vxyzu(:,:),tmunus(:,:,:)
 real                :: pi
 integer             :: i
 logical             :: firstpart
 real                :: tmunu(4,4)
 !print*, "entered get tmunu_all_exact"
 tmunu = 0.
 firstpart = .true.
 ! TODO write openmp parallel code
 do i=1, npart
    if (.not.isdead_or_accreted(xyzh(4,i)) .and. firstpart) then
       pi = get_pressure(ieos,xyzh(:,i),dens(i),vxyzu(:,i))
       call get_tmunu_exact(xyzh(:,i),metrics(:,:,:,i), metricderivs(:,:,:,i), &
         vxyzu(1:3,i),dens(i),vxyzu(4,i),pi,tmunus(:,:,i))
       !print*, "finished get_tmunu call!"
       firstpart = .false.
       !print*, "tmunu: ", tmunu
       !print*, "tmunus: ", tmunus(:,:,i)
       tmunu(:,:) = tmunus(:,:,i)
       !print*, "Got tmunu val: ", tmunu
       !stop
    else
       !print*, "setting tmunu for part: ", i
       tmunus(:,:,i) = tmunu(:,:)
    endif

 enddo
 !print*, "tmunu calc val is: ", tmunus(0,0,5)
end subroutine get_tmunu_all_exact


! Subroutine to calculate the covariant form of the stress energy tensor
! For a particle at position p
subroutine get_tmunu(x,metrici,v,dens,u,p,tmunu,verbose)
 use metric_tools,     only:unpack_metric
 use utils_gr,     only:get_u0
 real,    intent(in)  :: x(3),metrici(:,:,:),v(3),dens,u,p
 real,    intent(out) :: tmunu(0:3,0:3)
 logical, optional, intent(in) :: verbose
 real                 :: w,v4(0:3),uzero,u_upper(0:3),u_lower(0:3)
 real                 :: gcov(0:3,0:3), gcon(0:3,0:3)
 real                 :: gammaijdown(1:3,1:3),betadown(3),alpha
 integer              :: ierr,mu,nu

 ! Reference for all the variables used in this routine:
 ! w - the enthalpy
 ! gcov - the covariant form of the metric tensor
 ! gcon - the contravariant form of the metric tensor
 ! gammaijdown - the covariant form of the spatial metric
 ! alpha - the lapse
 ! betadown - the covariant component of the shift
 ! v4 - the uppercase 4 velocity in covariant form
 ! v - the fluid velocity v^x
 ! vcov - the covariant form of big V_i
 ! bigV - the uppercase contravariant V^i

 ! Calculate the enthalpy
 w = 1 + u + p/dens

 ! Get cov and con versions of the metric + spatial metric and lapse and shift
 ! Not entirely convinced that the lapse and shift calculations are acccurate for the general case!!
 !print*, "Before unpack metric "
 call unpack_metric(metrici,gcov=gcov,gcon=gcon,gammaijdown=gammaijdown,alpha=alpha,betadown=betadown)
 !print*, "After unpack metric"

!  if (present(verbose) .and. verbose) then
!     ! Do we get sensible values
!     print*, "Unpacked metric quantities..."
!     print*, "gcov: ", gcov
!     print*, "gcon: ", gcon
!     print*, "gammaijdown: ", gammaijdown
!     print* , "alpha: ", alpha
!     print*, "betadown: ", betadown
!     print*, "v4: ", v4
!  endif

 ! ! Need to change Betadown to betaup
 ! ! Won't matter at this point as it is allways zero
 ! ! get big V
 ! bigV(:) = (v(:) + betadown)/alpha

 ! ! We need the covariant version of the 3 velocity
 ! ! gamma_ij v^j = v_i where gamma_ij is the spatial metric
 ! do i=1, 3
 !    vcov(i) = gammaijdown(i,1)*bigv(1) + gammaijdown(i,2)*bigv(2) + gammaijdown(i,3)*bigv(3)
 ! enddo


 ! ! Calculate the lorentz factor
 ! lorentz = (1. -  (vcov(1)*bigv(1) + vcov(2)*bigv(2) + vcov(3)*bigv(3)))**(-0.5)

 ! ! Calculate the 4-velocity
 ! velshiftterm = vcov(1)*betadown(1) + vcov(2)*betadown(2) + vcov(3)*betadown(3)
 ! v4(0) = lorentz*(-alpha + velshiftterm)
 ! ! This should be vcov not v
 ! v4(1:3) = lorentz*vcov(1:3)


 ! We are going to use the same Tmunu calc as force GR
 ! And then lower it using the metric
 ! i.e calc T^{\mu\nu} and then lower it using the metric
 ! tensor
 ! lower-case 4-velocity (contravariant)
 v4(0) = 1.
 v4(1:3) = v(:)


 ! first component of the upper-case 4-velocity (contravariant)
 call get_u0(gcov,v,uzero,ierr)

 u_upper = uzero*v4
 do mu=0,3
    u_lower(mu) = gcov(mu,0)*u_upper(0) + gcov(mu,1)*u_upper(1) &
                   + gcov(mu,2)*u_upper(2) + gcov(mu,3)*u_upper(3)
 enddo

 ! Stress energy tensor in contravariant form
 do nu=0,3
    do mu=0,3
       tmunu(mu,nu) = w*dens*u_lower(mu)*u_lower(nu) + p*gcov(mu,nu)
    enddo
 enddo


!  if (present(verbose) .and. verbose) then
!     ! Do we get sensible values
!     print*, "Unpacked metric quantities..."
!     print*, "gcov: ", gcov
!     print*, "gcon: ", gcon
!     print*, "gammaijdown: ", gammaijdown
!     print* , "alpha: ", alpha
!     print*, "betadown: ", betadown
!     print*, "v4: ", v4
!  endif

!  if (verbose) then
!     print*, "tmunu part: ", tmunu
!     print*, "dens: ", dens
!     print*, "w: ", w
!     print*, "p: ", p
!     print*, "gcov: ", gcov
!  endif

 ! print*, "tmunu part: ", tmunu
 ! print*, "dens: ", dens
 ! print*, "w: ", w
 ! print*, "p: ", p
 ! print*, "gcov: ", gcov
 ! stop
end subroutine get_tmunu

subroutine get_tmunu_exact(x,metrici,metricderivsi,v,dens,u,p,tmunu)
 use metric_tools,     only:unpack_metric
 use utils_gr,         only:get_sqrtg
 real,    intent(in)  :: x(3),metrici(:,:,:),metricderivsi(0:3,0:3,3),v(3),dens,u,p
 real,    intent(out) :: tmunu(0:3,0:3)
 real                 :: w,v4(0:3),vcov(3),lorentz
 real                 :: gcov(0:3,0:3), gcon(0:3,0:3)
 real                 :: gammaijdown(1:3,1:3),betadown(3),alpha
 real                 :: velshiftterm
 real                 :: rhostar,rhoprim,negsqrtg
 integer              :: i,j

 ! Calculate the enthalpy
 ! enthalpy should be 1 as we have zero pressure
 ! or should have zero pressure
 w = 1
 ! Calculate the exact value of density from conserved density

 call unpack_metric(metrici,gcov=gcov,gcon=gcon,gammaijdown=gammaijdown,alpha=alpha,betadown=betadown)
 ! We need the covariant version of the 3 velocity
 ! gamma_ij v^j = v_i where gamma_ij is the spatial metric
 do i=1, 3
    vcov(i) = gammaijdown(i,1)*v(1) + gammaijdown(i,2)*v(2) + gammaijdown(i,3)*v(3)
 enddo

 ! Calculate the lorentz factor
 lorentz = (1. -  (vcov(1)*v(1) + vcov(2)*v(2) + vcov(3)*v(3)))**(-0.5)

 ! Calculate the 4-velocity
 velshiftterm = vcov(1)*betadown(1) + vcov(2)*betadown(2) + vcov(3)*betadown(3)
 v4(0) = lorentz*(-alpha + velshiftterm)
 v4(1:3) = lorentz*v(1:3)

 rhostar = 13.294563008157013D0
 call get_sqrtg(gcov,negsqrtg)
 ! Set/Calculate primitive density using rhostar exactly
 rhoprim = rhostar/(negsqrtg/alpha)


 ! Stress energy tensor
 do j=0,3
    do i=0,3
       tmunu(i,j) = rhoprim*w*v4(i)*v4(j) ! + p*gcov(i,j) neglect the pressure term as we don't care
    enddo
 enddo



end subroutine get_tmunu_exact

end module extern_gr
