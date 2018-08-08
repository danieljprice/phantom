!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: viscosity
!
!  DESCRIPTION:
!   Routines related to physical viscosity
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, eos, part, timestep
!+
!--------------------------------------------------------------------------
module viscosity
 implicit none
 integer, public :: irealvisc
 real, public :: shearparam, bulkvisc, HoverR

 public :: shearfunc,set_defaults_viscosity,dt_viscosity,viscinfo

 private

contains

!----------------------------------------------------------------
!+
!  set default values for the viscosity coefficients
!+
!----------------------------------------------------------------
subroutine set_defaults_viscosity

 irealvisc = 0   ! Physical viscosity
 shearparam = 0.1  ! alphadisc (if irealvisc=2) or nu if irealvisc=1
 bulkvisc = 0.0   ! bulk viscosity parameter in code units

 return
end subroutine set_defaults_viscosity

!----------------------------------------------------------------
!+
!  function which returns shear parameter \nu given
!  the alpha disc parameter, positions and sound speed
!  ie. \nu = alpha*cs*H
!
!  for irealvisc = 2 alpha is specified in the input file
!  for irealvisc = 1 \nu is specified directly in the input file
!+
!----------------------------------------------------------------
real function shearfunc(xi,yi,zi,spsoundi)
 use part,   only:xyzmh_ptmass
 use eos,    only:polyk,qfacdisc
 real, intent(in) :: xi,yi,zi,spsoundi
 real :: rsph2,omega1,H,r1,r2

 select case(irealvisc)
 case(2)
!
!--Shakura & Sunyaev disc viscosity
!  nu = alpha*cs*H
!
!  but H = cs/Omega, so nu = alpha*cs**2/Omega
!
!  (Omega is assumed to be Keplerian ie. r**-3/2)
!
    rsph2 = xi*xi + yi*yi + zi*zi
    omega1 = rsph2**0.75  ! omega1 = 1/omega
    shearfunc = shearparam*omega1*spsoundi**2

 case(1)
!
!--constant shear viscosity
!
    shearfunc = shearparam

 case(0)
!
!--default is zero
!
    shearfunc = 0.

 case(3)
!
!--Shakura & Sunyaev for ieos=14
!

    r1=sqrt((xi-xyzmh_ptmass(1,1))**2+(yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2)
    r2=sqrt((xi-xyzmh_ptmass(1,2))**2+(yi-xyzmh_ptmass(2,2))**2 + (zi-xyzmh_ptmass(3,2))**2)

    if(r2<r1) then
       H=sqrt(polyk)*r2**(-qfacdisc+1.5)
    else
       H=sqrt(polyk)*r1**(-qfacdisc+1.5)
    endif

    shearfunc=shearparam*spsoundi*H


 case default

    stop 'invalid choice for physical viscosity'

 end select

end function shearfunc

!----------------------------------------------------------------
!+
!  set explicit timestep for physical viscosity
!+
!----------------------------------------------------------------
real function dt_viscosity(xi,yi,zi,hi,spsoundi)
 use timestep,only:C_force
 real, intent(in) :: xi,yi,zi,hi,spsoundi
 real :: viscnu

 viscnu = shearfunc(xi,yi,zi,spsoundi)

 if (viscnu > tiny(viscnu)) then
    dt_viscosity = 0.4*C_force*hi*hi/viscnu
 else
    dt_viscosity = huge(dt_viscosity)
 endif

end function dt_viscosity

!----------------------------------------------------------------
!+
!  prints info about physical viscosity settings into the header
!+
!----------------------------------------------------------------
subroutine viscinfo(ivisc,iprint)
 use dim, only:maxp,maxdvdx
 integer, intent(in) :: ivisc,iprint

 select case(ivisc)
 case(3)
    write(iprint,"(a,es10.3)") ' Shakura-Sunyaev viscosity for binary systems, alpha_SS = ',shearparam
 case(2)
    write(iprint,"(a,es10.3)") ' Shakura-Sunyaev viscosity, alpha_SS = ',shearparam
 case(1)
    write(iprint,"(a,es10.3)") ' Constant physical viscosity, nu = ',shearparam
 case(0)
 case default
    write(iprint,"(a,es10.3)") ' Unknown setting for physical viscosity, nu = ',shearparam
 end select
 if (ivisc /= 0) then
    if (maxdvdx==maxp) then
       write(iprint,"(a,/)") ' (computed using two first derivatives)'
    else
       write(iprint,"(a,/)") ' (computed using direct second derivatives)'
    endif
 endif

end subroutine viscinfo

end module viscosity
