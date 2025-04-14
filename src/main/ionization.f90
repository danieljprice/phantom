!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
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
! :Dependencies: dim, eos_idealplusrad, io, part, physcon, units,
!   vectorutils
!
 implicit none

 real, allocatable, public, dimension(:)  :: eion
 logical, public                          :: done_ion_setup = .false.
 real, allocatable, private, dimension(:) :: logeion,arec,brec,crec,drec,arec1c,brec1c
 real, private                            :: frec,edge,tanh_c,dtanh_c
 real, parameter, private                 :: tanh_edge = 3.64673859532966

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

 if (.not. done_ion_setup) then
    allocate(eion(1:4),arec(2:4),brec(2:4),crec(1:4),drec(1:4),&
             arec1c(1:2),brec1c(1:2),logeion(1:4))
 endif

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

 done_ion_setup = .true.
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

!-----------------------------------------------------------------------
!+
!  Get recombination energy components
!+
!-----------------------------------------------------------------------
subroutine get_erec_components(logd,T,X,Y,erec)
 real, intent(in)     :: logd,T,X,Y
 real, intent(out)    :: erec(4)
 real, dimension(1:4) :: e,xi

 e(1) = eion(1)*X*0.5
 e(2) = eion(2)*X
 e(3) = eion(3)*Y*0.25
 e(4) = eion(4)*Y*0.25

 call get_xion(logd,T,X,Y,xi)
 erec = e*xi

end subroutine get_erec_components

!----------------------------------------------------------------
!+
!  Calculate thermal (gas + radiation internal energy) energy of a
!  gas particle. Inputs and outputs in code units
!+
!----------------------------------------------------------------
subroutine calc_thermal_energy(particlemass,ieos,xyzh,vxyzu,presi,tempi,ethi,rad)
 use dim,              only:do_radiation
 use part,             only:rhoh,iradxi
 use eos_idealplusrad, only:get_idealgasplusrad_tempfrompres,get_idealplusrad_enfromtemp
 use physcon,          only:radconst,Rg
 use units,            only:unit_density,unit_pressure,unit_ergg,unit_pressure
 integer, intent(in) :: ieos
 real, intent(in)    :: particlemass,presi,tempi,xyzh(4),vxyzu(4)
 real, intent(in), optional :: rad(:)
 real, intent(out)   :: ethi
 real                :: densi_cgs,mui

 select case (ieos)
 case(10,20) ! calculate just gas + radiation thermal energy
    densi_cgs = rhoh(xyzh(4),particlemass)*unit_density
    mui = densi_cgs * Rg * tempi / (presi*unit_pressure - radconst * tempi**4 / 3.) ! Get mu from pres and temp
    call get_idealplusrad_enfromtemp(densi_cgs,tempi,mui,ethi)
    ethi = particlemass * ethi / unit_ergg
 case default ! assuming internal energy = thermal energy
    ethi = particlemass * vxyzu(4)
    if (do_radiation) ethi  = ethi + particlemass*rad(iradxi)
 end select

end subroutine calc_thermal_energy


!----------------------------------------------------------------
!+
!  Solves three Saha equations simultaneously to return ion
!  fractions of hydrogen and helium. Assumes inputs in cgs units
!+
!----------------------------------------------------------------
subroutine ionisation_fraction(dens,temp,X,Y,xh0,xh1,xhe0,xhe1,xhe2)
 use physcon,     only:twopi,kboltz,eV,planckh,mass_electron_cgs,mass_proton_cgs
 use vectorutils, only:matrixinvert3D
 use io,          only:fatal
 real, intent(in)     :: dens,temp,X,Y
 real, intent(out)    :: xh0,xh1,xhe0,xhe1,xhe2
 real                 :: n,nh,nhe
 real                 :: A,B,C,const
 real                 :: xh1g,xhe1g,xhe2g
 real                 :: f,g,h
 real, parameter      :: chih0=13.598,chihe0=24.587,chihe1=54.418
 real, dimension(3,3) :: M,M_inv
 real, dimension(3)   :: dx
 integer              :: i,ierr

 nh = X * dens / mass_proton_cgs
 nhe = Y * dens / (4. * mass_proton_cgs)
 n = nh + nhe

 const = (sqrt(twopi * mass_electron_cgs * kboltz) / planckh)**3 / n

 A = 1. * const * temp**(1.5) * exp(-chih0 * eV / (kboltz * temp))
 B = 4. * const * temp**(1.5) * exp(-chihe0 * eV / (kboltz * temp))
 C = 1. * const * temp**(1.5) * exp(-chihe1 * eV / (kboltz * temp))

 xh1g = 0.4
 xhe1g = 0.3
 xhe2g = 0.2

 do i=1,50
    f = xh1g * (xh1g + xhe1g + 2*xhe2g) - A * ((nh/n) - xh1g)
    g = xhe1g * (xh1g + xhe1g + 2*xhe2g) - B * ((nhe/n) - xhe1g - xhe2g)
    h = xhe2g * (xh1g + xhe1g + 2*xhe2g) - C * xhe1g

    M(1,:) = (/ 2*xh1g + xhe1g + 2*xhe2g + A, xh1g, 2*xh1g /)
    M(2,:) = (/ xhe1g, xh1g + 2*xhe1g + 2*xhe2g + B, 2*xhe1g + B /)
    M(3,:) = (/ xhe2g, xhe2g - C, xh1g + xhe1g + 4*xhe2g /)

    call matrixinvert3D(M,M_inv,ierr)
    if (ierr /= 0) call fatal("ionisation_fraction","Error inverting matrix")

    dx = matmul(M_inv, (/ -f, -g, -h/))

    xh1g = xh1g + dx(1)
    xhe1g = xhe1g + dx(2)
    xhe2g = xhe2g + dx(3)
 enddo

 xh1 = max(xh1g * n / nh,1.e-99)
 xhe1 = max(xhe1g * n / nhe,1.e-99)
 xhe2 = max(xhe2g * n / nhe,1.e-99)
 xh0 = ((nh/n) - xh1g) * n / nh
 xhe0 = ((nhe/n) - xhe1g - xhe2g) * n / nhe

end subroutine ionisation_fraction

end module ionization_mod
