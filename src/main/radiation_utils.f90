!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module radiation_utils
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - X                  : *hydrogen mass fraction for MESA opacity table*
!   - Z                  : *metallicity for MESA opacity table*
!   - cv_type            : *how to get cv and mean mol weight (0=constant,1=mesa)*
!   - flux_limiter       : *limit radiation flux*
!   - implicit_radiation : *use implicit integration (Whitehouse, Bate & Monaghan 2005)*
!   - iopacity_type      : *opacity method (0=inf,1=mesa,2=constant,-1=preserve)*
!   - itsmax_rad         : *max number of iterations for radiation implicit solve*
!   - kappa_cgs          : *constant opacity value in cm2/g*
!   - tol_rad            : *tolerance on backwards Euler implicit solve of dxi/dt*
!
! :Dependencies: dim, eos, infile_utils, io, mesa_microphysics, part,
!   physcon, units
!
 implicit none
 public :: update_radenergy!,set_radfluxesandregions
 public :: set_radiation_and_gas_temperature_equal
 public :: radiation_and_gas_temperature_equal
 public :: get_rad_R
 public :: radiation_equation_of_state
 public :: radxi_from_Trad
 public :: Trad_from_radxi
 public :: ugas_from_Tgas
 public :: Tgas_from_ugas
 public :: get_opacity
 public :: get_kappa

 ! options for the input file, with default values
 real, public       :: tol_rad = 1.e-6
 integer, public    :: itsmax_rad = 250
 integer, public    :: cv_type = 0
 real, public       :: kappa_cgs = 0.3

 ! following declared public to avoid compiler warnings
 public :: solve_internal_energy_implicit_substeps
 public :: solve_internal_energy_explicit
 public :: solve_internal_energy_explicit_substeps

 ! radiation
 logical, public :: exchange_radiation_energy,limit_radiation_flux,implicit_radiation
 logical, public :: implicit_radiation_store_drad

 public :: set_defaults_radiation
 public :: write_options_radiation
 public :: read_options_radiation

 private

contains

!---------------------------------------------------------
!+
!  set default values for radiation options
!+
!---------------------------------------------------------
subroutine set_defaults_radiation()
 use dim, only:do_radiation
 use eos, only:iopacity_type

 ! radiation
 if (do_radiation) then
    exchange_radiation_energy = .true.
    limit_radiation_flux = .true.
    iopacity_type = 1
    implicit_radiation = .false.
 else
    exchange_radiation_energy = .false.
    limit_radiation_flux = .false.
    iopacity_type = 0
    implicit_radiation = .false.
 endif
 implicit_radiation_store_drad = .false.

end subroutine set_defaults_radiation

!---------------------------------------------------------
!+
!  write options to input file
!+
!---------------------------------------------------------
subroutine write_options_radiation(iunit)
 use infile_utils, only:write_inopt
 use eos,          only:X_in,Z_in,iopacity_type,ieos
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for radiation'
 call write_inopt(implicit_radiation,'implicit_radiation','use implicit integration (Whitehouse, Bate & Monaghan 2005)',iunit)
 call write_inopt(exchange_radiation_energy,'gas-rad_exchange','exchange energy between gas and radiation',iunit)
 call write_inopt(limit_radiation_flux,'flux_limiter','limit radiation flux',iunit)
 call write_inopt(iopacity_type,'iopacity_type','opacity method (0=inf,1=mesa,2=constant,-1=preserve)',iunit)
 if ((iopacity_type == 1) .and. (ieos /= 20)) then  ! for ieos=20, X, Z are already under EoS options
    call write_inopt(X_in,'X','hydrogen mass fraction for MESA opacity table',iunit)
    call write_inopt(Z_in,'Z','metallicity for MESA opacity table',iunit)
 elseif (iopacity_type == 2) then
    call write_inopt(kappa_cgs,'kappa_cgs','constant opacity value in cm2/g',iunit)
 endif
 if (implicit_radiation) then
    call write_inopt(tol_rad,'tol_rad','tolerance on backwards Euler implicit solve of dxi/dt',iunit)
    call write_inopt(itsmax_rad,'itsmax_rad','max number of iterations for radiation implicit solve',iunit)
    call write_inopt(cv_type,'cv_type','how to get cv and mean mol weight (0=constant,1=mesa)',iunit)
 endif

end subroutine write_options_radiation

!---------------------------------------------------------
!+
!  read options from input file
!+
!---------------------------------------------------------
subroutine read_options_radiation(db,nerr)
 use dim,          only:store_dust_temperature
 use eos,          only:iopacity_type
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(implicit_radiation,'implicit_radiation',db,errcount=nerr)
 if (implicit_radiation) store_dust_temperature = .true.
 call read_inopt(exchange_radiation_energy,'gas-rad_exchange',db,errcount=nerr,default=exchange_radiation_energy)
 call read_inopt(limit_radiation_flux,'flux_limiter',db,errcount=nerr,default=limit_radiation_flux)
 call read_inopt(iopacity_type,'iopacity_type',db,errcount=nerr,min=-1,max=2,default=iopacity_type)
 if (iopacity_type == 2) call read_inopt(kappa_cgs,'kappa_cgs',db,errcount=nerr,min=0.)
 if (implicit_radiation) then
    call read_inopt(cv_type,'cv_type',db,errcount=nerr,min=0,max=20,default=cv_type)
    call read_inopt(tol_rad,'tol_rad',db,errcount=nerr,min=epsilon(tol_rad),default=tol_rad)
    call read_inopt(itsmax_rad,'itsmax_rad',db,errcount=nerr,min=1,default=itsmax_rad)
 endif

end subroutine read_options_radiation

!-------------------------------------------------
!+
!  get R factor needed for flux limited diffusion
!+
!-------------------------------------------------
pure real function get_rad_R(rho,xi,flux,kappa) result(radR)
 real, intent(in) :: rho,xi,flux(3),kappa

 ! Note that flux is supposed to be grad(E) = grad(rho*xi)
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
subroutine set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,&
            rad,mu_local,npin)
 use part,      only:rhoh,igas,iradxi
 use eos,       only:gmw,gamma
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 real, intent(out)   :: rad(:,:)
 real,    intent(in), optional :: mu_local(:)
 integer, intent(in), optional :: npin
 real                :: rhoi,pmassi,mu
 integer             :: i,i1

 i1 = 0
 if (present(npin)) i1 = npin

 pmassi = massoftype(igas)
 mu = gmw
 do i=i1+1,npart
    rhoi = rhoh(xyzh(4,i),pmassi)
    if (present(mu_local)) mu = mu_local(i)
    rad(iradxi,i) = radiation_and_gas_temperature_equal(rhoi,vxyzu(4,i),gamma,mu)
 enddo

end subroutine set_radiation_and_gas_temperature_equal

!-------------------------------------------------
!+
!  set equal gas and radiation temperature
!+
!-------------------------------------------------
real function radiation_and_gas_temperature_equal(rho,u_gas,gamma,gmw) result(xi)
 use physcon,   only:Rg
 use units,     only:unit_ergg,get_radconst_code
 real, intent(in) :: rho,u_gas,gamma,gmw
 real :: temp,cv1,Erad

 cv1 = (gamma-1.)*gmw/Rg*unit_ergg

 temp = u_gas*cv1
 Erad = temp**4*get_radconst_code()
 xi   = Erad /rho

end function radiation_and_gas_temperature_equal

!---------------------------------------------------------
!+
!  get specific radiation energy from radiation temperature
!+
!---------------------------------------------------------
real function radxi_from_Trad(rho,Trad) result(radxi)
 use units, only:get_radconst_code
 real, intent(in) :: rho,Trad

 radxi = Trad**4*get_radconst_code()/rho

end function radxi_from_Trad

!---------------------------------------------------------
!+
!  get radiation temperature from the specific radiation energy
!+
!---------------------------------------------------------
real function Trad_from_radxi(rho,radxi) result(Trad)
 use units, only:get_radconst_code
 real, intent(in) :: rho,radxi

 Trad = (rho*radxi/get_radconst_code())**0.25

end function Trad_from_radxi

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
subroutine update_radenergy(npart,xyzh,fxyzu,vxyzu,rad,radprop,dt,mu_local)
 use part,         only:rhoh,igas,massoftype,ikappa,iradxi,iphase,iamtype,ithick
 use eos,          only:gmw,gamma
 use units,        only:get_radconst_code,get_c_code,unit_velocity
 use physcon,      only:Rg
 use io,           only:warning
 use dim,          only:maxphase,maxp
 real, intent(in)    :: dt,xyzh(:,:),fxyzu(:,:),radprop(:,:)
 real, intent(in), optional :: mu_local(:)
 real, intent(inout) :: vxyzu(:,:),rad(:,:)
 integer, intent(in) :: npart
 real :: ui,pmassi,rhoi,xii
 real :: ack,a,cv1,kappa,dudt,etot,unew
 integer :: i

 pmassi        = massoftype(igas)

 a   = get_radconst_code()
 cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2

 !$omp parallel do default(none)&
 !$omp private(kappa,ack,rhoi,ui)&
 !$omp private(dudt,xii,etot,unew)&
 !$omp firstprivate(cv1)&
 !$omp shared(rad,radprop,xyzh,vxyzu,mu_local,gamma)&
 !$omp shared(fxyzu,pmassi,maxphase,maxp)&
 !$omp shared(iphase,npart)&
 !$omp shared(dt,a,unit_velocity)
 do i = 1,npart
    if (maxphase==maxp) then
       if (iamtype(iphase(i)) /= igas) cycle
    endif
    kappa = radprop(ikappa,i)
    ack = get_radconst_code()*get_c_code()*kappa

    rhoi = rhoh(xyzh(4,i),pmassi)
    ui   = vxyzu(4,i)
    dudt = fxyzu(4,i)
    xii  = rad(iradxi,i)
    etot = ui + xii
    unew = ui
    if (xii < -epsilon(0.)) then
       call warning('radiation','radiation energy is negative before exchange', i)
    endif
    if (present(mu_local)) cv1 = (gamma-1.)*mu_local(i)/Rg*unit_velocity**2
    call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt,i)
    ! call solve_internal_energy_implicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_explicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt,di)
    vxyzu(4,i) = unew
    rad(iradxi,i) = etot - unew
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

!--------------------------------------------------------------------
!+
!  get opacity from u and rho in code units (and precalculated cv)
!+
!--------------------------------------------------------------------
real function get_kappa(opacity_type,u,cv,rho) result(kappa)
 integer, intent(in) :: opacity_type
 real, intent(in)    :: u,cv,rho
 real                :: temp

 temp = u/cv
 call get_opacity(opacity_type,rho,temp,kappa)

end function get_kappa

!--------------------------------------------------------------------
!+
!  calculate opacities
!+
!--------------------------------------------------------------------
subroutine get_opacity(opacity_type,density,temperature,kappa)
 use mesa_microphysics, only:get_kappa_mesa
 use units,             only:unit_density,unit_opacity
 real, intent(in)  :: density, temperature
 real, intent(out) :: kappa
 integer, intent(in) :: opacity_type
 real :: kapt,kapr,rho_cgs

 select case(opacity_type)
 case(1)
    !
    ! calculate opacity from the MESA tables
    !
    rho_cgs = density*unit_density
    call get_kappa_mesa(rho_cgs,temperature,kappa,kapt,kapr)
    kappa = kappa/unit_opacity

 case(2)
    !
    ! constant opacity
    !
    kappa = kappa_cgs/unit_opacity

 case default
    !
    ! infinite opacity
    !
    kappa = huge(1.)

 end select

end subroutine get_opacity

! subroutine set_radfluxesandregions(npart,radiation,xyzh,vxyzu)
!   use part,    only: igas,massoftype,rhoh,ifluxx,ifluxy,ifluxz,ithick,iradxi,ikappa
!   use part,    only: eos_vars,ics
!   use options, only:ieos
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
!     cs = eos_vars(ics,i)
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
