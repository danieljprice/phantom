!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: radiation_utils
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Sergei Biriukov
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, eos, io, part, physcon, units
!+
!--------------------------------------------------------------------------
module radiation_utils
 implicit none
 public :: update_radenergy!,set_radfluxesandregions
 public :: set_radiation_and_gas_temperature_equal
 public :: get_rad_R

 private

contains
!-------------------------------------------------
!+
!  get R factor needed for flux limited diffusion
!+
!-------------------------------------------------
pure real function get_rad_R(rho,xi,flux,kappa) result(radR)
 real, intent(in) :: rho,xi,flux(3),kappa

 if (abs(xi) > epsilon(xi)) then
    radR = sqrt(dot_product(flux,flux))/(kappa*rho*rho*xi)
 else
    radR = 0.
 endif

end function get_rad_R

!-------------------------------------------------
!+
!  set equal gas and radiation temperatures
!+
!-------------------------------------------------
subroutine set_radiation_and_gas_temperature_equal(npart,gamma,xyzh,vxyzu,massoftype,radiation,opacity)
 use physcon,   only:Rg,steboltz,c
 use units,     only:unit_opacity,unit_ergg,unit_density
 use part,      only:rhoh,igas,iradxi,ikappa
 use eos,       only:gmw

 integer, intent(in) :: npart
 real, intent(in)    :: gamma,xyzh(:,:),vxyzu(:,:),massoftype(:)
 real, intent(inout) :: radiation(:,:)
 real, intent(in), optional :: opacity
 real                :: kappa,kappa_code,Tgas,rhoi,pmassi
 integer             :: i

 kappa = 1e5
 if (present(opacity)) then
    if (opacity > 0.) kappa = opacity
 endif
 kappa_code = kappa/unit_opacity
 pmassi = massoftype(igas)

 do i=1,npart
    rhoi = rhoh(xyzh(4,i),pmassi)
    Tgas = gmw*((gamma-1.)*vxyzu(4,i)*unit_ergg)/Rg
    radiation(iradxi,i) = (4.0*steboltz*Tgas**4.0/c/(rhoi*unit_density))/unit_ergg
    radiation(ikappa,i) = kappa_code
    !print*,i,' Tgas = ',Tgas,'rad=',radiation(iradxi,i),radiation(ikappa,i)
 enddo

end subroutine set_radiation_and_gas_temperature_equal

!--------------------------------------------------------------------
!+
!  integrate radiation energy exchange terms over a time interval dt
!+
!--------------------------------------------------------------------
subroutine update_radenergy(npart,xyzh,fxyzu,vxyzu,radiation,dt)
 use part,         only:rhoh,igas,massoftype,ikappa,iradxi,iphase,iamtype,ithick
 use eos,          only:gmw,gamma
 use units,        only:get_steboltz_code,get_c_code,unit_velocity
 use physcon,      only:Rg
 use io,           only:warning
 use dim,          only:maxphase,maxp
 real, intent(in)    :: dt,xyzh(:,:),fxyzu(:,:)
 real, intent(inout) :: vxyzu(:,:),radiation(:,:)
 integer, intent(in) :: npart
 real :: ui,pmassi,rhoi,xii
 real :: ack,a,cv1,kappa,dudt,etot,unew
 real :: c_code, steboltz_code
 integer :: i

 pmassi        = massoftype(igas)
 steboltz_code = get_steboltz_code()
 c_code        = get_c_code()

 a   = 4.*steboltz_code/c_code
 cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2

 !$omp parallel do default(none)&
 !$omp private(kappa,ack,rhoi,ui)&
 !$omp private(dudt,xii,etot,unew)&
 !$omp shared(radiation,xyzh,vxyzu)&
 !$omp shared(fxyzu,pmassi,maxphase,maxp)&
 !$omp shared(iphase,npart)&
 !$omp shared(dt,cv1,a,steboltz_code)
 do i = 1,npart
    if (maxphase==maxp) then
       if (iamtype(iphase(i)) /= igas) cycle
    endif
    kappa = radiation(ikappa,i)
    ack = 4.*steboltz_code*kappa

    rhoi = rhoh(xyzh(4,i),pmassi)
    ui   = vxyzu(4,i)
    dudt = fxyzu(4,i)
    xii  = radiation(iradxi,i)
    etot = ui + xii
    unew = ui
    if (xii < 0.) then
       call warning('radiation','radiation energy is negative before exchange', i)
    endif
!     if (i==584) then
!        print*, 'Before:  ', 'T_gas=',unew*cv1,'T_rad=',(rhoi*(etot-unew)/a)**(1./4.)
!     endif
    call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt,i)
    ! call solve_internal_energy_implicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_explicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt,di)
    vxyzu(4,i) = unew
    radiation(iradxi,i) = etot - unew
!   if (i==584) then
!      print*, 'After:   ', 'T_gas=',unew*cv1,'T_rad=',unew,etot,(rhoi*(etot-unew)/a)**(1./4.)
!         read*
!     endif
    if (radiation(iradxi,i) < 0.) then
       call warning('radiation','radiation energy negative after exchange', i,var='xi',val=radiation(iradxi,i))
       radiation(iradxi,i) = 0.
    endif
    if (vxyzu(4,i) < 0.) then
       call warning('radiation','thermal energy negative after exchange', i,var='u',val=vxyzu(4,i))
       vxyzu(4,i) = 0.
    endif
 enddo
 !$omp end parallel do
end subroutine update_radenergy

!--------------------------------------------------------------------
!+
!  update internal energy using implicit substeps
!+
!--------------------------------------------------------------------
subroutine solve_internal_energy_implicit_substeps(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1

 real     :: fu,dfu,eps,dts,dunew,unewp,uip
 integer  :: iter,i,level

 unew = ui
 unewp = 2.*ui
 eps = 1e-8
 dunew = (unewp-unew)/unew
 level = 1
 do while((abs(dunew) > eps).and.(level <= 2**10))
    unewp = unew
    dts   = dt/level
    uip   = ui
    do i=1,level
       iter  = 0
       fu    = huge(1.)
       do while((abs(fu) > eps).and.(iter < 10))
          iter = iter + 1
          fu   = unew/dts - uip/dts - dudt - ack*(rho*(etot-unew)/a - (unew*cv1)**4)
          dfu  = 1./dts + ack*(rho/a + 4.*(unew**3*cv1**4))
          unew = unew - fu/dfu
       enddo
       uip = unew
    enddo
    dunew = (unewp-unew)/unew
    level = level*2
 enddo
end subroutine solve_internal_energy_implicit_substeps

!--------------------------------------------------------------------
!+
!  update internal energy using implicit backwards Euler method
!+
!--------------------------------------------------------------------
subroutine solve_internal_energy_implicit(unew,u0,rho,etot,dudt,ack,a,cv1,dt,i)
 real, intent(out) :: unew
 real, intent(in)  :: u0, etot, dudt, dt, rho, ack, a, cv1
 integer, intent(in) :: i
 real     :: fu,dfu,eps,uold
 integer  :: iter

 unew = u0
 uold = 2*u0

 iter = 0
 eps = 1e-16
 do while ((abs(unew-uold) > eps).and.(iter < 10))
    uold = unew
    iter = iter + 1
    fu   = unew/dt - u0/dt - dudt - ack*(rho*(etot-unew)/a - (unew*cv1)**4)
    dfu  = 1./dt + ack*(rho/a + 4.*(unew**3*cv1**4))
    unew = unew - fu/dfu
 enddo

end subroutine solve_internal_energy_implicit

!--------------------------------------------------------------------
!+
!  update internal energy using explicit Euler method
!+
!--------------------------------------------------------------------
subroutine solve_internal_energy_explicit(unew,ui,rho,etot,dudt,ack,a,cv1,dt,di)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1
 integer, intent(in) :: di

 unew = ui + dt*(dudt + ack*(rho*(etot-ui)/a - (ui*cv1)**4))

end subroutine solve_internal_energy_explicit

!--------------------------------------------------------------------
!+
!  update internal energy using series of explicit Euler steps
!+
!--------------------------------------------------------------------
subroutine solve_internal_energy_explicit_substeps(unew,ui,rho,etot,dudt,ack,a,cv1,dt,di)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1
 integer, intent(in) :: di
 real     :: du,eps,dts,unews,uis
 integer  :: i,level

 level = 1
 eps   = 1e-8
 dts   = dt

 unew = ui + dts*(dudt + ack*(rho*(etot-ui)/a - (ui*cv1)**4))
 unews = 2*unew
 du = (unew-unews)/unew
 do while((abs(du) > eps).and.(level <= 2**10))
    unews = unew
    level = level*2
    dts  = dt/level
    unew = 0.
    uis  = ui
    do i=1,level
       unew = uis + dts*(dudt + ack*(rho*(etot-uis)/a - (uis*cv1)**4))
       uis  = unew
    enddo
    du = (unew-unews)/unew
 enddo

end subroutine solve_internal_energy_explicit_substeps

! subroutine set_radfluxesandregions(npart,radiation,xyzh,vxyzu)
!   use part,    only: igas,massoftype,rhoh,ifluxx,ifluxy,ifluxz,ithick,iradxi,ikappa
!   use eos,     only: get_spsound
!   use options, only: ieos
!   use physcon, only:c
!   use units,   only:unit_velocity
!
!   real, intent(inout)    :: radiation(:,:),vxyzu(:,:)
!   real, intent(in)       :: xyzh(:,:)
!   integer, intent(in)    :: npart
!
!   integer :: i
!   real :: pmassi,H,cs,rhoi,r,c_code,lambdai,prevfrac
!
!   pmassi = massoftype(igas)
!
!   prevfrac = 100.*count(radiation(ithick,:)==1)/real(size(radiation(ithick,:)))
!   sch = sch * (1. + 0.005 * (80. - prevfrac))
!   print*, "-}+{- RADIATION Scale Height Multiplier: ", sch
!
!   radiation(ithick,:) = 1
!
!   c_code = c/unit_velocity
!
!   do i = 1,npart
!     rhoi = rhoh(xyzh(4,i),pmassi)
!     ! if (rhoi < 2e-4) then
!     !   if (xyzh(1,i) < 0.) then
!     !     radiation(ifluxy:ifluxz,i) = 0.
!     !     radiation(ifluxx,i)        = -rhoi*abs(radiation(iradxi,i))*0.5
!     !     radiation(ithick,i) = 0
!     !   elseif (xyzh(1,i) > 0.) then
!     !     radiation(ifluxy:ifluxz,i) = 0.
!     !     radiation(ifluxx,i)        =  rhoi*abs(radiation(iradxi,i))*0.5
!     !     radiation(ithick,i) = 0
!     !   endif
!     ! endif
!     cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
!     r  = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
!     H  = cs*sqrt(r**3)
!     if (abs(xyzh(3,i)) > sch*H) then
!        radiation(ithick,i) = 0
!        ! if (xyzh(3,i) <= 0.) then
!        !   radiation(ifluxx:ifluxy,i) = 0.
!        !   radiation(ifluxz,i)        =  rhoi*abs(radiation(iradxi,i))
!        !   radiation(ithick,i) = 0
!        ! elseif (xyzh(3,i) > 0.) then
!        !   radiation(ifluxx:ifluxy,i) = 0.
!        !   radiation(ifluxz,i)        =  -rhoi*abs(radiation(iradxi,i))
!        !   radiation(ithick,i) = 0
!        ! endif
!     endif
!     ! Ri = sqrt(&
!     !   dot_product(radiation(ifluxx:ifluxz,i),radiation(ifluxx:ifluxz,i)))&
!     !   /(radiation(ikappa,i)*rhoi*rhoi*radiation(iradxi,i))
!     ! lambda = (2. + Ri)/(6. + 3*Ri + Ri*Ri)
!     ! lambdai = 1./3
!     ! ! print*, xyzh(4,i), lambdai/radiation(ikappa,i)/rhoi*c_code/cs
!     ! ! read*
!     ! if (xyzh(4,i) > lambdai/radiation(ikappa,i)/rhoi*c_code/cs) then
!     !    radiation(ithick,i) = 0
!     ! endif
!   enddo
! end subroutine set_radfluxesandregions
! subroutine mcfost_do_analysis()
!
! end subroutine mcfost_do_analysis
!
end module radiation_utils
