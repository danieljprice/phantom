module gravwaveutils

implicit none

public :: calculate_strain

contains

subroutine calculate_strain(hx,hp,hxx,hpp,xyzh,vxyz,axyz,pmass,npart)
!  use prompting,      only:prompt
  use units,          only:umass,udist,utime
  use physcon,        only:gg,c
  real, intent(out) :: hx,hp,hxx,hpp
  real, intent(in)  :: xyzh(:,:), vxyz(:,:), axyz(:,:), pmass
  integer, intent(in) :: npart
  real :: q(6), ddq(6)
  !logical, save       :: firstcall = .true.
  integer :: i
  real, save          :: distan

  distan=0.03*3.0e24

  ! initialise quadrupole to zero
    q(:)=0

  ! calculate the components of the traceless quadrupole--not necessary but maybe useful
    do i=1,npart
     if(xyzh(4,i)>tiny(xyzh)) then  !if not accreted
      q(1)=q(1)+pmass*(xyzh(1,i)*xyzh(1,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qxx
      q(2)=q(2)+pmass*(xyzh(1,i)*xyzh(2,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qxy
      q(3)=q(3)+pmass*(xyzh(1,i)*xyzh(3,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qxz
      q(4)=q(4)+pmass*(xyzh(2,i)*xyzh(2,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qyy
      q(5)=q(5)+pmass*(xyzh(2,i)*xyzh(3,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qyz
      q(6)=q(6)+pmass*(xyzh(3,i)*xyzh(3,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qzz
     end if
    enddo

  ! initialise the second time derivative of the quadrupole to zero
    ddq(:)=0

  ! calculate the second time derivative of the traceless quadrupole
   do i=1,npart
     if(xyzh(4,i)>tiny(xyzh)) then !if not accreted
      ddq(1)=ddq(1)+pmass*(2*vxyz(1,i)*vxyz(1,i)+xyzh(1,i)*axyz(1,i)+xyzh(1,i)*axyz(1,i)) !ddqxx
      ddq(2)=ddq(2)+pmass*(2*vxyz(1,i)*vxyz(2,i)+xyzh(1,i)*axyz(2,i)+xyzh(2,i)*axyz(1,i)) !ddqxy
      ddq(3)=ddq(3)+pmass*(2*vxyz(1,i)*vxyz(3,i)+xyzh(1,i)*axyz(3,i)+xyzh(3,i)*axyz(1,i)) !ddqxz
      ddq(4)=ddq(4)+pmass*(2*vxyz(2,i)*vxyz(2,i)+xyzh(2,i)*axyz(2,i)+xyzh(2,i)*axyz(2,i)) !ddqyy
      ddq(5)=ddq(5)+pmass*(2*vxyz(2,i)*vxyz(3,i)+xyzh(2,i)*axyz(3,i)+xyzh(3,i)*axyz(2,i)) !ddqyz
      ddq(6)=ddq(6)+pmass*(2*vxyz(3,i)*vxyz(3,i)+xyzh(3,i)*axyz(3,i)+xyzh(3,i)*axyz(3,i)) !ddqzz
     end if
   enddo

   ! gw strain in the direction perpendicular to the orbit, theta=0
     hp = gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*(ddq(1)-ddq(4))
     hx = 2.*gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*ddq(2)
   ! gw strain in the plane of the orbit, theta=pi/2
     hpp = gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*(ddq(1)-ddq(6))
     hxx = -2.*gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*ddq(3)

    print*, 'ciao'

end subroutine calculate_strain


end module gravwaveutils
