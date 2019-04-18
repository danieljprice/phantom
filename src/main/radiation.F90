#ifdef RADIATION
module radiation
 implicit none
 public :: update_radenergy,set_radfluxesandregions

 private

contains

subroutine update_radenergy(npart,xyzh,fxyzu,vxyzu,radenevol,radkappa,dt)
 use part,         only:rhoh,igas,massoftype,radenflux
 use eos,          only:gmw,gamma
 use units,        only:udist,utime,&
                        unit_energ,unit_velocity
 use physcon,      only:Rg,steboltz,c
 use io,           only:fatal

  real, intent(in)    ::dt,xyzh(:,:),fxyzu(:,:),radkappa(:)
  real, intent(inout) ::vxyzu(:,:),radenevol(:)
  integer, intent(in) ::npart

  real :: ui,pmassi,rhoi,xii
  real :: ack,a,cv1,kappa,dudt,etot,unew
  real :: c_code, steboltz_code
  integer :: i

  pmassi = massoftype(igas)
  steboltz_code = steboltz/(unit_energ/(udist**2*utime))
  c_code        = c/unit_velocity

  a   = 4.*steboltz_code/c_code
  cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2

  !$omp parallel do default(none)&
  !$omp private(kappa,ack,rhoi,ui)&
  !$omp private(dudt,xii,etot,unew)&
  !$omp shared(radkappa,xyzh,vxyzu)&
  !$omp shared(fxyzu,pmassi,radenevol)&
  !$omp shared(radenflux)&
  !$omp shared(dt,cv1,a,steboltz_code)
  do i = 1,npart
    ! if (maxphase==maxp) then
    !    if (iamtype(iphase(i)) /= igas) cycle ! skip non-gas
    ! endif
    kappa = radkappa(i)
    ack = 4.*steboltz_code*kappa

    rhoi = rhoh(xyzh(4,i),pmassi)
    ui   = vxyzu(4,i)
    dudt = fxyzu(4,i)
    xii = radenevol(i)
    etot = ui + xii
    unew = ui
    call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    call solve_internal_energy_implicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_explicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    vxyzu(4,i) = unew
    radenevol(i) = etot - unew
    ! print*, 'T_gas=',unew*cv1,'T_rad=',(rhoi*(etot-unew)/a)**(1./4.)
    if (radenevol(i) <= 0) call fatal('radiation','radenevol is negative', i)
  end do
  !$omp end parallel do
end subroutine update_radenergy

subroutine solve_internal_energy_implicit_substeps(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1

 real     :: fu,dfu,eps,dts,dunew,unewp,uip
 integer  :: iter,i,level

 unew = ui
 unewp = huge(1.)
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
     end do
     uip = unew
   end do
   dunew = (unewp-unew)/unew
   level = level*2
 end do
end subroutine solve_internal_energy_implicit_substeps

subroutine solve_internal_energy_implicit(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1
 real     :: fu,dfu,eps
 integer  :: iter

 unew = ui
 fu = huge(1.)
 iter = 0
 eps = 1e-14
 do while((abs(fu) > eps).and.(iter < 10))
   iter = iter + 1
   fu  = unew/dt - ui/dt - dudt - ack*(rho*(etot-unew)/a - (unew*cv1)**4)
   dfu = 1./dt + ack*(rho/a + 4.*(unew**3*cv1**4))
   unew = unew - fu/dfu
 end do
end subroutine solve_internal_energy_implicit

subroutine solve_internal_energy_explicit(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1
 real     :: du,eps,dts,unews,uis
 integer  :: i,level

 du    = huge(1.)
 level = 1
 eps   = 1e-8
 dts   = dt

 unew = ui + dts*(dudt + ack*(rho*(etot-ui)/a - (ui*cv1)**4))

 do while((abs(du) > eps*unew).and.(level <= 2**10))
   unews = unew
   level = level*2
   dts  = dt/level
   unew = 0.
   uis  = ui
   do i=1,level
     unew = uis + dts*(dudt + ack*(rho*(etot-ui)/a - (ui*cv1)**4))
     uis  = unew
   end do
   du = (unew-unews)/unew
 end do
end subroutine solve_internal_energy_explicit

subroutine set_radfluxesandregions(npart,radthick,radenevol,radenflux,xyzh,vxyzu)
  use part,    only: igas,massoftype,rhoh
  use eos,     only: get_spsound
  use options, only: ieos

  real, intent(inout)    :: radenflux(:,:),vxyzu(:,:)
  real, intent(in)       :: xyzh(:,:),radenevol(:)
  logical, intent(inout) :: radthick(:)
  integer, intent(in)    :: npart

  integer :: i
  real :: pmassi,H,cs,rhoi,r

  pmassi = massoftype(igas)
  radthick(:) = .true.

  do i = 1,npart
    rhoi = rhoh(xyzh(4,i),pmassi)
    cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
    r  = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    H  = 2*cs*sqrt(r**3)
    if (abs(xyzh(3,i)) > H) then
      if (xyzh(3,i) <= 0.) then
        radenflux(1:2,i) = 0.
        radenflux(3,i)   =  rhoi*abs(radenevol(i))
        radthick(i) = .false.
      else if (xyzh(3,i) > 0.) then
        radenflux(1:2,i) =  0.
        radenflux(3,i)   = -rhoi*abs(radenevol(i))
        radthick(i) = .false.
      end if
    end if
  end do
end subroutine set_radfluxesandregions
end module radiation
#endif
