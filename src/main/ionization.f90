!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module ionization_mod
!
! ionization_mod
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none

 real, allocatable, public, dimension(:)  :: eion
 real, allocatable, private, dimension(:) :: logeion,arec,brec,crec,drec,arec1c,brec1c
 real, private                            :: frec,edge,tanh_c,dtanh_c
 real, parameter, private                 :: tanh_edge=3.6467

 interface rapid_tanh
  module procedure rapid_tanhs,rapid_tanhv
 end interface rapid_tanh
 interface rapid_dtanh
  module procedure rapid_dtanhs,rapid_dtanhv
 end interface rapid_dtanh

contains

!-----------------------------------------------------------------------
!+
!  Rapid hyperbolic tangent function for scalars
!+
!-----------------------------------------------------------------------
function rapid_tanhs(x)
 real, intent(in) :: x
 real             :: a,b,x2,rapid_tanhs

 if (abs(x)>=tanh_edge) then
    rapid_tanhs = sign(1.,x)-tanh_c/x
 else
    x2 = x*x
    a = ((x2+105.)*x2+945.)*x
    b = ((x2+28.)*x2+63.)*15.
    rapid_tanhs = a/b
 endif
end function rapid_tanhs

!-----------------------------------------------------------------------
!+
!  Rapid hyperbolic tangent function for vectors
!+
!-----------------------------------------------------------------------
function rapid_tanhv(x)
 real, intent(in)         :: x(:)
 real, dimension(size(x)) :: rapid_tanhv
 real                     :: a,b,x2
 integer                  :: i

 do i = 1,size(x)
    if (abs(x(i)) >= tanh_edge) then
       rapid_tanhv(i) = sign(1.,x(i))-tanh_c/x(i)
    else
       x2 = x(i)*x(i)
       a = ((x2+105.)*x2+945.)*x(i)
       b = ((x2+28.)*x2+63.)*15.
       rapid_tanhv(i) = a/b
    endif
 enddo
end function rapid_tanhv

!-----------------------------------------------------------------------
!+
!  Rapid hyperbolic tangent derivative (1/cosh^2) for scalars
!+
!-----------------------------------------------------------------------
function rapid_dtanhs(x)
 real, intent(in) :: x
 real             :: a,b,x2,rapid_dtanhs

 if (abs(x) >= tanh_edge) then
    rapid_dtanhs = dtanh_c/x**2
 else
    x2=x**2
    a = (((x2-21.)*x2+420.)*x2-6615.)*x2+59535.
    b = (x2+28.)*x2+63.
    rapid_dtanhs = a/(15.*b**2)
 endif
end function rapid_dtanhs

!-----------------------------------------------------------------------
!+
!  Rapid hyperbolic tangent derivative (1/cosh^2) for vectors
!+
!-----------------------------------------------------------------------
function rapid_dtanhv(x)
 implicit none
 real, intent(in)         :: x(:)
 real, dimension(size(x)) :: rapid_dtanhv
 real                     :: a,b,x2
 integer                  :: i

 do i = 1,size(x)
    if (abs(x(i)) >= tanh_edge) then
       rapid_dtanhv(i) = dtanh_c/(x(i)*x(i))
    else
       x2=x(i)**2
       a = (((x2-21.)*x2+420.)*x2-6615.)*x2+59535.
       b = (x2+28.)*x2+63.
       rapid_dtanhv(i) = a/(15.*b**2)
    endif
 enddo
end function rapid_dtanhv

!-----------------------------------------------------------------------
!+
!  Set up all fitting coefficients
!+
!-----------------------------------------------------------------------
subroutine ionization_setup
 use physcon, only:Rg
 real :: x

 allocate(eion(1:4),arec(2:4),brec(2:4),crec(1:4),drec(1:4),&
          arec1c(1:2),brec1c(1:2),logeion(1:4))

 eion(1) = 4.36e12   ! H2   [erg/mol]
 eion(2) = 1.312e13  ! HI   [erg/mol]
 eion(3) = 2.3723e13 ! HeI  [erg/mol]
 eion(4) = 5.2505e13 ! HeII [erg/mol]
 logeion(1:4) = log10(eion(1:4)/Rg)

 ! These fitting parameters are tuned for eosDT in the range X=0.6-0.8
 frec = 0.005
 arec(2:4) = (/ 0.821, 0.829, 0.846 /)
 brec(2:4) = (/ 0.055, 0.055, 0.055 /)
 crec(1:4) = (/ 0.02, 0.025, 0.015, 0.015 /)
 drec(1:4) = (/ 0.05, 0.05, 0.05, 0.05 /)

 arec1c(1:2) = (/ log10(3500.)/logeion(1), 0.753 /)
 brec1c(1:2) = (/ 0., 0.055 /)
 edge = -9.4

 ! Parameter for rapid tanh function
 x=tanh_edge
 tanh_c = x*(1.-((x*x+105.)*x*x+945.)*x/(((x*x+28.)*x*x+63.)*15.))
 dtanh_c= x*x*((((x*x-21.)*x*x+420.)*x*x-6615.)*x*x+59535.)&
               / (15.*((x*x+28.)*x*x+63.)**2)

 return
end subroutine ionization_setup

!-----------------------------------------------------------------------
!+
!  Molecular hydrogen fit coefficients
!+
!-----------------------------------------------------------------------
real function arec1(x)
 real, intent(in):: x ! logd

 if (x<edge) then
    arec1 = arec1c(1)
 else
    arec1 = arec1c(2)
 endif
end function arec1

real function brec1(x)
 real, intent(in):: x

 if (x<edge) then
    brec1 = brec1c(1)
 else
    brec1 = brec1c(2)
 endif
end function brec1

!-----------------------------------------------------------------------
!+
!  Get ionization fractions (and dxdT) given rho and T
!+
!-----------------------------------------------------------------------
subroutine get_xion(logd,T,X,Y,xion,dxion)
 real, intent(in)            :: logd,T,X,Y
 real, intent(out)           :: xion(1:4)
 real, intent(out), optional :: dxion(1:4)
 real                        :: logQ,logT,Yfac
 real, dimension(1:4)        :: Ttra,width,arg

 logT = log10(T)
 logQ = max(-14.,logd)-2.*logT+12.

 Yfac = 1.-frec*Y

 Ttra (1)   = arec1(logd)*logeion(1)+brec1(logd)*logQ
 Ttra (2:4) = Yfac*arec(2:4)*logeion(2:4)+brec(2:4)*logQ
 width(1:4) = Ttra(1:4)*crec(1:4)*(1.+drec(1:4)*logQ)
 arg  (1:4) = (logT-Ttra(1:4))/width(1:4)

 xion(1:4) = 0.5*(rapid_tanh(arg(1:4))+1.)

 if (present(dxion)) then
    dxion(1) = ( width(1)*(1.+2.*brec1(logd)) &
                + 2.*crec(1)*(logT-Ttra(1))&
                  *(brec1(logd)*(1.+drec(1)*logQ)+drec(1)*Ttra(1)) )&
              / (2.*T*width(1)*width(1)) * rapid_dtanh(arg(1))
    dxion(2:4) = ( width(2:4)*(1.+2.*brec(2:4)) &
                 + 2.*crec(2:4)*(logT-Ttra(2:4))&
                   *(brec(2:4)*(1.+drec(2:4)*logQ)+drec(2:4)*Ttra(2:4)) )&
                / (2.*T*width(2:4)*width(2:4)) * rapid_dtanh(arg(2:4))
 endif

end subroutine get_xion

!-----------------------------------------------------------------------
!+
!  Get recombination energy and mean molecular weight given rho and T
!+
!-----------------------------------------------------------------------
subroutine get_erec_imurec(logd,T,X,Y,erec,imurec,derecdT,dimurecdT)
 real, intent(in)            :: logd,T,X,Y
 real, intent(out)           :: erec,imurec
 real, intent(out), optional :: derecdT,dimurecdT
 real, dimension(1:4)        :: e,xi,zi

! CAUTION: This is only a poor man's way of implementing recombination energy.
!          It only should be used for -3.5<logQ<-6 where logQ=logrho-2logT+12.

 e(1) = eion(1)*X*0.5
 e(2) = eion(2)*X
 e(3) = eion(3)*Y*0.25
 e(4) = eion(4)*Y*0.25

 if (present(derecdT).or.present(dimurecdT)) then
    call get_xion(logd,T,X,Y,xi,zi)
 else
    call get_xion(logd,T,X,Y,xi)
 endif

 erec = sum(e(1:4)*xi(1:4))
 if (present(derecdT)) then
    derecdT = sum(e(1:4)*zi(1:4))
 endif

 imurec = (0.5*xi(1)+xi(2))*X+0.25*(xi(3)+xi(4)-1.)*Y+0.5
 if (present(dimurecdT)) then
    dimurecdT = (0.5*zi(1)+zi(2))*X+0.25*(zi(3)+zi(4))*Y
 endif

 return
end subroutine get_erec_imurec

!-----------------------------------------------------------------------
!+
!  Get the mean molecular weight for partially ionised plasma
!+
!-----------------------------------------------------------------------
subroutine get_imurec(logd,T,X,Y,imurec,dimurecdT)
 real, intent(in)            :: logd,T,X,Y
 real, intent(out)           :: imurec
 real, intent(out), optional :: dimurecdT
 real, dimension(1:4)        :: xi,zi

! CAUTION: This is only a poor man's way of implementing recombination energy.
!          It only should be used for -3.5<logQ<-6 where logQ=logrho-2logT+12.

 if (present(dimurecdT)) then
    call get_xion(logd,T,X,Y,xi,zi)
 else
    call get_xion(logd,T,X,Y,xi)
 endif

 imurec = (0.5*xi(1)+xi(2))*X+0.25*(xi(3)+xi(4)-1.)*Y+0.5
 if (present(dimurecdT)) then
    dimurecdT = (0.5*zi(1)+zi(2))*X+0.25*(zi(3)+zi(4))*Y
 endif

 return
end subroutine get_imurec

!-----------------------------------------------------------------------
!+
!  Get recombination energy given rho and T
!+
!-----------------------------------------------------------------------
real function get_erec(logd,T,X,Y)
 real, intent(in)     :: logd,T,X,Y
 real, dimension(1:4) :: e, xi

! CAUTION: This is only a poor man's way of implementing recombination energy.
!          It only should be used for -3.5<logQ<-6 where logQ=logrho-2logT+12.

 e(1) = eion(1)*X*0.5
 e(2) = eion(2)*X
 e(3) = eion(3)*Y*0.25
 e(4) = eion(4)*Y*0.25

 call get_xion(logd,T,X,Y,xi)

 get_erec = sum(e(1:4)*xi(1:4))

 return
end function get_erec

end module ionization_mod
