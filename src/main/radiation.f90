#ifdef RADIATION
module radiation
 implicit none
 public :: update_energy

 private

contains

subroutine update_energy(npart,xyzh,fxyzu,dt,vxyzu,radenergy)
  use part,     only:rhoh,igas,massoftype
  use eos,      only:gmw,gamma
  use units,    only:udist,utime,umass,&
                     unit_energ,unit_velocity
  use physcon,  only:Rg,steboltz,c

  real, intent(in)    ::dt,xyzh(:,:),fxyzu(:,:)
  real, intent(inout) ::vxyzu(:,:),radenergy(:)
  integer, intent(in) ::npart

  real :: ui,pmassi,rhoi
  real :: ack,a,cv1,kappa,dudt,etot,unew
  real :: c_code, steboltz_code, kappa_code
  integer :: i

  pmassi = massoftype(igas)
  kappa = 4.e-5
  steboltz_code = steboltz/(unit_energ/(udist**2*utime))
  kappa_code    = kappa/(udist**2/umass)
  c_code        = c/unit_velocity

  ack = 4.*steboltz_code*kappa_code
  a   = 4.*steboltz_code/c_code
  cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2

  do i = 1,npart
    rhoi = rhoh(xyzh(4,i),pmassi)
    ui   = vxyzu(4,i)
    dudt = fxyzu(4,i)
    etot = ui + radenergy(i)
    unew = ui
    call solve_internal_energy(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    vxyzu(4,i) = unew
    radenergy(i) = etot - unew
  end do

end subroutine update_energy

subroutine solve_internal_energy(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in) :: ui, etot, dudt, dt, rho, ack, a, cv1
 real :: fu,dfu

 unew = ui
 fu = huge(1.)
 do while((abs(fu) > 1e-6).and.(dt > 0.))
   fu  = unew/dt - ui/dt - dudt - ack*(rho*(etot-unew)/a - (unew*cv1)**4)
   dfu = 1./dt + ack*(rho/a + 4.*(unew**3*cv1**4))
   unew = unew - fu/dfu
 end do
end subroutine

end module radiation
#endif
