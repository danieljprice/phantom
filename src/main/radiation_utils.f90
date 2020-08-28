!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module radiation_utils
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, io, part, physcon, units
!
 implicit none
 public :: update_radenergy!,set_radfluxesandregions
 public :: set_radiation_and_gas_temperature_equal
 public :: radiation_and_gas_temperature_equal
 public :: get_rad_R
 public :: radiation_equation_of_state
 public :: T_from_Etot
 public :: radE_from_Trad
 public :: Trad_from_radE
 public :: ugas_from_Tgas
 public :: Tgas_from_ugas

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

!-------------------------------------------------------------
!+
!  set equal gas and radiation temperatures for all particles
!+
!-------------------------------------------------------------
subroutine set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)
 use part,      only:rhoh,igas,iradxi
 use eos,       only:gmw,gamma
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 real, intent(out)   :: rad(:,:)
 real                :: rhoi,pmassi
 integer             :: i

 pmassi = massoftype(igas)
 do i=1,npart
    rhoi = rhoh(xyzh(4,i),pmassi)
    rad(iradxi,i) = radiation_and_gas_temperature_equal(rhoi,vxyzu(4,i),gamma,gmw)
    !print*,i,' Tgas = ',Tgas,'rad=',rad(iradxi,i)
 enddo

end subroutine set_radiation_and_gas_temperature_equal

!-------------------------------------------------
!+
!  set equal gas and radiation temperature
!+
!-------------------------------------------------
real function radiation_and_gas_temperature_equal(rho,u_gas,gamma,gmw) result(xi)
 use physcon,   only:Rg,steboltz,c
 use units,     only:unit_ergg,unit_density,get_steboltz_code,get_c_code
 real, intent(in) :: rho,u_gas,gamma,gmw
 real :: temp,cv1,a,Erad,steboltz_code,c_code

 steboltz_code = get_steboltz_code()
 c_code        = get_c_code()

 a   = 4.*steboltz_code/c_code
 cv1 = (gamma-1.)*gmw/Rg*unit_ergg

 temp = u_gas*cv1
 Erad = temp**4*a
 xi   = Erad /rho

end function radiation_and_gas_temperature_equal

!---------------------------------------------------------
!+
!  solve for the temperature for which Etot=Erad+ugas is
!  satisfied assuming Tgas=Trad
!+
!---------------------------------------------------------
real function T_from_Etot(rho,etot,gamma,gmw) result(temp)
 use physcon,   only:Rg
 use units,     only:unit_ergg,unit_density,get_steboltz_code,get_c_code
 real, intent(in)    :: rho,etot,gamma,gmw
 real                :: steboltz_code,c_code,a,cv1
 real                :: numerator,denominator,correction
 real, parameter     :: tolerance = 1d-15

 steboltz_code = get_steboltz_code()
 c_code        = get_c_code()

 a   = 4.*steboltz_code/c_code
 cv1 = (gamma-1.)*gmw/Rg*unit_ergg

 if (temp <= 0.) then
    temp = etot*cv1  ! Take gas temperature as initial guess
 endif

 correction = huge(0.)
 do while (abs(correction) > tolerance*temp)
    numerator   = etot*rho - rho*temp/cv1 - a*temp**4
    denominator =  - rho/cv1 - 4.*a*temp**3
    correction  = numerator/denominator
    temp        = temp - correction
 enddo
end function T_from_Etot

!---------------------------------------------------------
!+
!  get the radiation energy from the raditaion temperature
!+
!---------------------------------------------------------
real function radE_from_Trad(Trad) result(radE)
 use units,     only:get_steboltz_code,get_c_code
 real, intent(in)  :: Trad
 real              :: a,steboltz_code,c_code

 steboltz_code = get_steboltz_code()
 c_code        = get_c_code()

 a = 4. * steboltz_code/c_code

 radE = Trad**4*a
end function radE_from_Trad

!---------------------------------------------------------
!+
!  get the radiation temperature from the radiation energy
!+
!---------------------------------------------------------
real function Trad_from_radE(radE) result(Trad)
 use units,     only:get_steboltz_code,get_c_code
 real, intent(in)    :: radE
 real                :: a,steboltz_code,c_code

 steboltz_code = get_steboltz_code()
 c_code        = get_c_code()

 a = 4. * steboltz_code/c_code

 Trad = (radE/a)**0.25
end function Trad_from_radE

!---------------------------------------------------------
!+
!  get the internal energy from the gas temperature
!+
!---------------------------------------------------------
real function ugas_from_Tgas(Tgas,gamma,gmw) result(ugas)
 use physcon,   only:Rg
 use units,     only:unit_ergg
 real, intent(in)  :: Tgas,gamma,gmw
 real              :: cv1

 cv1 = (gamma-1.)*gmw/Rg*unit_ergg

 ugas = Tgas/cv1
end function ugas_from_Tgas

!---------------------------------------------------------
!+
!  get gas temperature from the internal energy
!+
!---------------------------------------------------------
real function Tgas_from_ugas(ugas,gamma,gmw) result(Tgas)
 use physcon,   only:Rg
 use units,     only:unit_ergg
 real, intent(in)    :: ugas,gamma,gmw
 real                :: cv1

 cv1 = (gamma-1.)*gmw/Rg*unit_ergg

 Tgas = ugas*cv1
end function Tgas_from_ugas

!--------------------------------------------------------------------
!+
!  integrate radiation energy exchange terms over a time interval dt
!+
!--------------------------------------------------------------------
subroutine update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,dt)
 use part,         only:rhoh,igas,massoftype,ikappa,iradxi,iphase,iamtype,ithick
 use eos,          only:gmw,gamma
 use units,        only:get_steboltz_code,get_c_code,unit_velocity
 use physcon,      only:Rg
 use io,           only:warning
 use dim,          only:maxphase,maxp
 real, intent(in)    :: dt,xyzh(:,:),fxyzu(:,:),radprop(:,:)
 real, intent(inout) :: vxyzu(:,:),rad(:,:)
 integer, intent(in) :: npart
 real :: ui,pmassi,rhoi,xii
 real :: ack,a,cv1,kappa,dudt,etot,unew
 real :: c_code,steboltz_code
 integer :: i

 pmassi        = massoftype(igas)
 steboltz_code = get_steboltz_code()
 c_code        = get_c_code()

 a   = 4.*steboltz_code/c_code
 cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2

 !$omp parallel do default(none)&
 !$omp private(kappa,ack,rhoi,ui)&
 !$omp private(dudt,xii,etot,unew)&
 !$omp shared(rad,radprop,xyzh,vxyzu)&
 !$omp shared(fxyzu,pmassi,maxphase,maxp)&
 !$omp shared(iphase,npart)&
 !$omp shared(dt,cv1,a,steboltz_code)
 do i = 1,npart
    if (maxphase==maxp) then
       if (iamtype(iphase(i)) /= igas) cycle
    endif
    kappa = radprop(ikappa,i)
    ack = 4.*steboltz_code*kappa

    rhoi = rhoh(xyzh(4,i),pmassi)
    ui   = vxyzu(4,i)
    dudt = fxyzu(4,i)
    xii  = rad(iradxi,i)
    etot = ui + xii
    unew = ui
    if (xii < -epsilon(0.)) then
       call warning('radiation','radiation energy is negative before exchange', i)
    endif
!     if (i==584) then
!        print*, 'Before:  ', 'T_gas=',unew*cv1,'T_rad=',(rhoi*(etot-unew)/a)**(1./4.)
!     endif
    call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt,i)
    ! call solve_internal_energy_implicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_explicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt,di)
    vxyzu(4,i) = unew
    rad(iradxi,i) = etot - unew
!   if (i==584) then
!      print*, 'After:   ', 'T_gas=',unew*cv1,'T_rad=',unew,etot,(rhoi*(etot-unew)/a)**(1./4.)
!         read*
!     endif
    if (rad(iradxi,i) < 0.) then
       call warning('radiation','radiation energy negative after exchange', i,var='xi',val=rad(iradxi,i))
       rad(iradxi,i) = 0.
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
          fu   = unew - uip - dts*dudt - dts*ack*(rho*(etot-unew)/a - (unew*cv1)**4)
          dfu  = 1. + dts*ack*(rho/a + 4.*(unew**3*cv1**4))
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
 real     :: fu,dfu,uold
 integer  :: iter

 unew = u0
 uold = 2*u0
 iter = 0
 do while ((abs(unew-uold) > epsilon(unew)).and.(iter < 10))
    uold = unew
    iter = iter + 1
    fu   = unew - u0 - 0.*dt*dudt - dt*ack*(rho*(etot-unew)/a - (unew*cv1)**4)
    dfu  = 1. + dt*ack*(rho/a + 4.*(unew**3*cv1**4))
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

!--------------------------------------------------------------------
!+
!  calculate radiation Pressure from radiation Energy
!+
!--------------------------------------------------------------------
subroutine radiation_equation_of_state(radPi, Xii, rhoi)
 real, intent(out) :: radPi
 real, intent(in) :: Xii, rhoi

 radPi = 1. / 3. * Xii * rhoi

end subroutine radiation_equation_of_state

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
