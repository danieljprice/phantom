#ifdef RADIATION
module radiation
 implicit none
 public :: update_energy

 private

contains

! print for individual timesteps where radiation is limiting condition for courant dt

subroutine update_energy(npart,xyzh,fxyzu,vxyzu,radenergy,radkappa,dt)
  use part,     only:rhoh,igas,massoftype
  use eos,      only:gmw,gamma
  use units,    only:udist,utime,&
                     unit_energ,unit_velocity
  use physcon,  only:Rg,steboltz,c
  use io,       only:fatal

  real, intent(in)    ::dt,xyzh(:,:),fxyzu(:,:),radkappa(:)
  real, intent(inout) ::vxyzu(:,:),radenergy(:)
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
  !$omp shared(fxyzu,pmassi,radenergy)&
  !$omp shared(dt,cv1,a,steboltz_code)
  do i = 1,npart
    kappa = radkappa(i)
    ack = 4.*steboltz_code*kappa

    rhoi = rhoh(xyzh(4,i),pmassi)
    ui   = vxyzu(4,i)
    dudt = fxyzu(4,i)
    xii = radenergy(i)
    etot = ui + xii
    unew = ui
    call solve_internal_energy(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    vxyzu(4,i) = unew
    radenergy(i) = etot - unew
    if (radenergy(i) < 0) call fatal('radiation','radenergy is negative', i)
  end do
  !$omp end parallel do
end subroutine update_energy

subroutine solve_internal_energy(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in) :: ui, etot, dudt, dt, rho, ack, a, cv1
 real     :: fu,dfu,eps
 integer  :: iter

 unew = ui
 fu = huge(1.)
 iter = 0
 eps = 1e-8
 do while((abs(fu) > eps).and.(dt > 0.).and.(iter < 10))
   iter = iter + 1
   fu  = unew/dt - ui/dt - dudt - ack*(rho*(etot-unew)/a - (unew*cv1)**4)
   dfu = 1./dt + ack*(rho/a + 4.*(unew**3*cv1**4))
   unew = unew - fu/dfu
 end do
end subroutine solve_internal_energy

end module radiation
#endif
