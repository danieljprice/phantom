!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_mesa
!
! None
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: mesa_microphysics, physcon
!

 use mesa_microphysics

 implicit none
 logical,private :: mesa_initialised = .false.

contains

!----------------------------------------------------------------
!+
!  subroutine initialises the mesa eos tables
!+
!----------------------------------------------------------------
subroutine init_eos_mesa(x,z,ierr)
 real, intent(in) :: x,z
 integer, intent(out) :: ierr

 ierr=0
 if ((x + z > 1.) .or. (x < 0.) .or. (z < 0.)) then
    ierr=-1
    return
 endif
 !to not initialise mesa multiple times. Otherwise could call finish_eos_mesa instead of return
 if (mesa_initialised) return
 mesa_initialised = .true.

 call get_environment_variable('MESA_DATA_DIR',mesa_eos_dir)
 mesa_eos_prefix="output_DE_"
 call get_environment_variable('MESA_DATA_DIR',mesa_opacs_dir)
 mesa_opacs_suffix=""
 mesa_eos_full_output=.true.
 mesa_binary_data=.true.

!Initialize the mesa_microphysics module
 mesa_eos_z1=0.d0   ! lowest Z
 mesa_eos_h1=0.d0   ! lowest H
 mesa_eos_dz=0.02d0 ! delta Z
 mesa_eos_dh=0.2d0  ! delta H in EOS tables

 call get_eos_constants_mesa(ierr)
 if (ierr /= 0) return

 call read_eos_mesa(x,z,ierr)
 call get_opacity_constants_mesa
 call read_opacity_mesa(x,z)

end subroutine init_eos_mesa

!----------------------------------------------------------------
!+
!  subroutine closes the mesa eos tables
!+
!----------------------------------------------------------------
subroutine finish_eos_mesa

 call deallocate_arrays_mesa

end subroutine finish_eos_mesa

!----------------------------------------------------------------
!+
!  subroutine returns pressure, temp, gamma1 as a function of
!  density and internal energy
!+
!----------------------------------------------------------------
subroutine get_eos_pressure_temp_gamma1_mesa(den,eint,pres,temp,gam1,ierr)
 real, intent(in) :: den,eint
 real, intent(out) :: pres,temp,gam1
 integer, intent(out) :: ierr

 call getvalue_mesa(den,eint,2,pres,ierr)
 call getvalue_mesa(den,eint,4,temp)
 call getvalue_mesa(den,eint,11,gam1)

end subroutine get_eos_pressure_temp_gamma1_mesa

!----------------------------------------------------------------
!+
!  subroutine returns kappa as a function of
!  density, temperature and composition
!+
!----------------------------------------------------------------
subroutine get_eos_kappa_mesa(den,temp,kappa,kappat,kappar)
 real, intent(in)  :: den, temp
 real, intent(out) :: kappa, kappat, kappar

 call get_kappa_mesa(den,temp,kappa,kappat,kappar)

end subroutine get_eos_kappa_mesa

!----------------------------------------------------------------
!+
!  subroutine returns kappa as a function of
!  density, temperature and composition
!+
!----------------------------------------------------------------
real function get_eos_1overmu_mesa(den,u) result(rmu)
 real, intent(in) :: den,u

 rmu = get_1overmu_mesa(den,u)

end function get_eos_1overmu_mesa

!----------------------------------------------------------------
!+
!  subroutine returns pressure and temperature as
!  a function of density/internal energy
!+
!----------------------------------------------------------------
subroutine get_eos_pressure_temp_mesa(den,eint,pres,temp)
 real, intent(in) :: den, eint
 real, intent(out) :: pres, temp

 call getvalue_mesa(den,eint,2,pres)
 call getvalue_mesa(den,eint,4,temp)

end subroutine get_eos_pressure_temp_mesa

!----------------------------------------------------------------
!+
!  subroutine returns internal energy and temperature from
!  density and internal energy using bisection method
!+
!----------------------------------------------------------------
pure subroutine get_eos_eT_from_rhop_mesa(rho,pres,eint,temp,guesseint)
 real, intent(in)           :: rho,pres
 real, intent(out)          :: eint,temp
 real, intent(in), optional :: guesseint
 real                       :: err,eintguess,eint1,eint2,&
                               eint3,pres1,pres2,pres3,left,right,mid
 real, parameter            :: tolerance = 1d-15
 integer                    :: ierr

 if (present(guesseint)) then
    eintguess = guesseint
    eint1 = 1.005 * eintguess  ! Tight lower bound
    eint2 = 0.995 * eintguess  ! Tight upper bound
 else
    eintguess = 1.5*pres/rho
    eint1 = 10. * eintguess  ! Guess lower bound
    eint2 = 0.1 * eintguess  ! Guess upper bound
 endif

 call getvalue_mesa(rho,eint1,2,pres1,ierr)
 call getvalue_mesa(rho,eint2,2,pres2,ierr)
 left  = pres - pres1
 right = pres - pres2

 ! If lower and upper bounds do not contain roots, extend them until they do
 do while (left*right > 0.)
    eint1 = 0.99 * eint1
    eint2 = 1.01 * eint2
    call getvalue_mesa(rho,eint1,2,pres1,ierr)
    call getvalue_mesa(rho,eint2,2,pres2,ierr)
    left  = pres - pres1
    right = pres - pres2
 enddo

 ! Start bisecting
 err = huge(1.)
 do while (abs(err) > tolerance)
    call getvalue_mesa(rho,eint1,2,pres1,ierr)
    call getvalue_mesa(rho,eint2,2,pres2,ierr)
    left  = pres - pres1
    right = pres - pres2
    eint3 = 0.5*(eint2+eint1)
    call getvalue_mesa(rho,eint3,2,pres3,ierr)
    mid = pres - pres3

    if (left*mid < 0.) then
       eint2 = eint3
    elseif (right*mid < 0.) then
       eint1 = eint3
    elseif (mid == 0.) then
       eint = eint3
       exit
    endif

    eint = eint3
    err = (eint2 - eint1)/eint1
 enddo

 call getvalue_mesa(rho,eint,4,temp,ierr)

end subroutine get_eos_eT_from_rhop_mesa


!----------------------------------------------------------------
!+
!  subroutine returns internal energy from density and internal
!  energy using bisection method. Assumes cgs units
!
!  Note: Needs unit testing (test_eos)
!+
!----------------------------------------------------------------
pure subroutine get_eos_u_from_rhoT_mesa(rho,temp,eint,guesseint)
 use physcon, only:kb_on_mh
 real, intent(in)           :: rho,temp
 real, intent(out)          :: eint
 real, intent(in), optional :: guesseint
 real                       :: err,eintguess,eint1,eint2,&
                               eint3,temp1,temp2,temp3,left,right,mid
 real, parameter            :: tolerance = 1d-15
 integer                    :: ierr

 if (present(guesseint)) then
    eintguess = guesseint
    eint1 = 1.005 * eintguess  ! Tight lower bound
    eint2 = 0.995 * eintguess  ! Tight upper bound
 else
    eintguess = 1.5*kb_on_mh*temp
    eint1 = 10. * eintguess  ! Guess lower bound
    eint2 = 0.1 * eintguess  ! Guess upper bound
 endif

 call getvalue_mesa(rho,eint1,4,temp1,ierr)
 call getvalue_mesa(rho,eint2,4,temp2,ierr)
 left  = temp - temp1
 right = temp - temp2

 ! If lower and upper bounds do not contain roots, extend them until they do
 do while (left*right > 0.)
    eint1 = 0.99 * eint1
    eint2 = 1.01 * eint2
    call getvalue_mesa(rho,eint1,4,temp1,ierr)
    call getvalue_mesa(rho,eint2,4,temp2,ierr)
    left  = temp - temp1
    right = temp - temp2
 enddo

 ! Start bisecting
 err = huge(1.)
 do while (abs(err) > tolerance)
    call getvalue_mesa(rho,eint1,4,temp1,ierr)
    call getvalue_mesa(rho,eint2,4,temp2,ierr)
    left  = temp - temp1
    right = temp - temp2
    eint3 = 0.5*(eint1+eint2)
    call getvalue_mesa(rho,eint3,4,temp3,ierr)
    mid = temp - temp3

    if (left*mid < 0.) then
       eint2 = eint3
    elseif (right*mid < 0.) then
       eint1 = eint3
    elseif (mid == 0.) then
       eint = eint3
       exit
    endif

    eint = eint3
    err = (eint2 - eint1)/eint1
 enddo

end subroutine get_eos_u_from_rhoT_mesa

!----------------------------------------------------------------
!+
!  subroutine returns various quantities as
!  a function of density/internal energy
!+
!----------------------------------------------------------------
subroutine get_eos_various_mesa(den,eint,pres,proint,peint,temp,troint,teint,entrop,abad,gamma1,gam)
 real, intent(in) :: den, eint
 real, intent(out) :: pres, temp, proint, peint, troint, teint, entrop, abad, gamma1, gam

 call getvalue_mesa(den,eint,2,pres)
 call getvalue_mesa(den,eint,4,temp)
 call getvalue_mesa(den,eint,5,proint)
 call getvalue_mesa(den,eint,6,peint)
 call getvalue_mesa(den,eint,7,troint)
 call getvalue_mesa(den,eint,8,teint)
 call getvalue_mesa(den,eint,9,entrop)
 call getvalue_mesa(den,eint,10,abad)
 call getvalue_mesa(den,eint,11,gamma1)
 call getvalue_mesa(den,eint,12,gam)

end subroutine get_eos_various_mesa

end module eos_mesa
