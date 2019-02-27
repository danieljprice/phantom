module radiation
 implicit none
 public :: update_Trad

 private

contains

subroutine update_Trad()
  use part,     only:npart,xyzh,vxyzu,radenergy,rhoh,igas,massoftype
  use eos,      only:gmw,gamma
  use options,  only:ieos
  use units,    only:unit_density,unit_pressure
  use physcon,  only:Rg,steboltz,c

  real :: xi,yi,zi,hi,ui,pmassi,rhoi,ponrho,pri,spsoundi
  real :: Tgas,Trad, Taim
  integer :: i

  pmassi = massoftype(igas)

  do i = 1,npart
    xi     = xyzh(1,i)
    yi     = xyzh(2,i)
    zi     = xyzh(3,i)
    hi     = xyzh(4,i)
    rhoi   = rhoh(hi,pmassi)
    ui     = vxyzu(4,i)
    ! call equationofstate(ieos,ponrho,spsoundi,rhoi,xi,yi,zi,)
    pri    = ui*rhoi*(gamma-1.0)

    Tgas = gmw*(pri*unit_pressure)/(rhoi*unit_density)/Rg
    Trad = (radenergy(i)*(rhoi*unit_density)*c/4.0/steboltz)**(1./4.)
    Taim = (Tgas+Trad)/2

    vxyzu(4,i) = Taim*Rg*(rhoi*unit_density)/gmw/unit_pressure/(gamma-1.0)/rhoi
    radenergy(i) = 4.0*steboltz*Taim**4.0/c/(rhoi*unit_density)
  end do
end subroutine update_Trad

end module radiation
