!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setfixedentropysurf
!
! This module replaces the steep surface entropy profile, due to surface
! recombination, with a gentle, linear extrapolation of the entropy
! interior to this region. We solve this equations of hydrostatic
! equilibrium and mass continuity to obtain the density and pressure profile
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: eos, kernel, physcon, table_utils
!
 implicit none

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calculates the constant entropy softened surface
!  profile. Assumes input profile quantities are sorted from centre to
!  surface. 
!+
!-----------------------------------------------------------------------
subroutine set_fixedS_surface(mcore,m,rho,r,pres,ene,temp,ierr)
 use eos,         only:calc_temp_and_ene
 use physcon,     only:pi,gg,solarm,solarr,kb_on_mh
 use table_utils, only:interpolator,flip_array
 real, intent(in)    :: mcore
 real, intent(inout) :: m(:),rho(:),r(:),pres(:),ene(:),temp(:)
 real                :: eneguess,mcore_cm
 real                :: msurf,dm(-1:0),dS_by_dm
 real, allocatable   :: S(:),m_alloc(:),r_alloc(:),rho_alloc(:),pres_alloc(:),testentropy(:),entropyorig(:)
 integer             :: j,i,surfidx,Nmax,ierr
 logical             :: isort_decreasing,iexclude_core_mass
 ! msurf: Mass coordinate beyond which we perform the profile replacement

 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .false.     ! Needs to be true if to be read by Phantom
 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .false.   ! Needs to be true if to be read by Phantom

 msurf = 9.8 * solarm !11.35 * solarm ! Mass coordinate beyond which we flatten the entropy
 call interpolator(m,msurf,surfidx) ! Find index closest to msurf
 Nmax = size(r) - surfidx

 ! Only need dm and S three steps behind the surface to calculate dS/dm there
 ! Weird indexing here so that the index 0 is the surfidx
 dm(-1) = m(surfidx-1) - m(surfidx-2)
 dm(0)  = m(surfidx) - m(surfidx-1)
 allocate(S(-2:0)) 
 do j = -2,0
    S(j) = entropy(rho(surfidx+j),pres(surfidx+j),ierr)
 enddo

 ! Find entropy gradient at m(surfidx-1)
 dS_by_dm = ( dm(-1)**2 * S(0) &
            + ( dm(0)**2 - dm(-1)**2 ) * S(-1) &
            - dm(0)**2 * S(-2) ) &
            / ( dm(-1) * dm(0) * (dm(-1) + dm(0)) )

 ! Fill up entropy grid using the entropy gradient   
 deallocate(S)     
 allocate(S(0:Nmax)) ! Repurpose variable name "S" for full surface entropy grid
 S(0) = entropy(rho(surfidx),pres(surfidx),ierr)
 do i = 1,Nmax
    S(i) = S(0) + dS_by_dm * (m(surfidx+i)-m(surfidx))
 enddo


 ! TEST: OUTPUT ENTROPY TO DATA FILE
 allocate(testentropy(size(r)),entropyorig(size(r)))
 do i = 1,size(r)
    entropyorig(i) = entropy(rho(i),pres(i),ierr)
 enddo
 testentropy(1:surfidx-1) = entropyorig(1:surfidx-1)
 testentropy(surfidx:size(r)) = S(0:Nmax)

 ! TEST: SUPPLY ORIGINAL ENTROPY TO CHECK ORIGINAL RHO AND PRES RECOVERED
 !S = entropyorig(surfidx:size(r))
 !testentropy(surfidx:size(r)) = entropyorig(surfidx:size(r))
 ! TEST: SUPPLY ORIGINAL ENTROPY TO CHECK ORIGINAL RHO AND PRES RECOVERED

 call write_entropy('./profiles/12_0_RSG/entropy_mesa.dat', testentropy, entropyorig)
 !stop
 ! TEST: OUTPUT ENTROPY TO DATA FILE


 ! Make allocatable copies (see instructions of solve_fixedS_surface_rhoP)
 allocate(m_alloc(-1:Nmax))
 m_alloc(-1:Nmax) = m(surfidx-1:size(m))
 allocate(r_alloc(-1:Nmax))
 r_alloc(-1:Nmax) = r(surfidx-1:size(r))
 allocate(rho_alloc(0:Nmax))
 rho_alloc(0:Nmax) = rho(surfidx:size(rho))
 allocate(pres_alloc(-1:Nmax))
 pres_alloc(-1:Nmax) = pres(surfidx-1:size(pres))

 call solve_fixedS_surface_rhoP(m_alloc,S,r_alloc,rho_alloc,pres_alloc,ierr)
 rho(surfidx:size(rho))     = rho_alloc(0:Nmax)
 pres(surfidx-1:size(pres)) = pres_alloc(-1:Nmax)
 r(surfidx:size(r))         = r_alloc(0:Nmax)

 ! Recalculate temperature and energy profiles
 call calc_temp_and_ene(rho(1),pres(1),ene(1),temp(1),ierr)
 do i = 2,size(rho)-1
    eneguess = ene(i-1)
    call calc_temp_and_ene(rho(i),pres(i),ene(i),temp(i),ierr,eneguess)
 enddo
 ene(size(rho))  = 0. ! Zero surface internal energy
 temp(size(rho)) = 0. ! Zero surface temperature

 ! Reverse arrays so that data is sorted from stellar surface to stellar centre.
 if (isort_decreasing) then
    call flip_array(m)
    call flip_array(pres)
    call flip_array(temp)
    call flip_array(r)
    call flip_array(rho)
    call flip_array(ene)
 endif

 if (iexclude_core_mass) then
    mcore_cm = mcore * solarm
    m = m - mcore_cm
 endif

end subroutine set_fixedS_surface


!-----------------------------------------------------------------------
!+
!  Solve for constant entropy surface density and pressure profiles
!+
!-----------------------------------------------------------------------
subroutine solve_fixedS_surface_rhoP(m,S,r,rho,pres,ierr)
 use physcon, only:gg
 real, allocatable, dimension(:), intent(in)    :: m,S
 real, allocatable, dimension(:), intent(inout) :: rho,pres,r
 integer, intent(out)                           :: ierr
 integer                                        :: Nmax,i
 real, allocatable, dimension(:)                :: dm
 real, parameter                                :: one_on_4pi=0.07957747154594767d0

 ! INSTRUCTIONS

 ! Input variables should be given in the following format. msurf refers to the mass
 ! coordinate beyond which the stellar structure is to be replaced with the fixed
 ! entropy profile

 ! m(-1:Nmax):    Array of mass grid to be softened, satisfying m(0)=msurf and m(Nmax)=M
 ! S(0:Nmax):     Entropy grid to find the matching density and pressure profiles for
 ! r(-1,Nmax):    Give r(-1)=(r at m(-1)) as input. Outputs radial grid.
 ! rho(0:Nmax):   Give rho(0)=(rho at msurf) as input. Outputs density profile.
 ! pres(-1:Nmax): Give pres(-1:0)=(p at m(-1:0)) as input. Outputs pressure profile.

 ! ierr: 

 ierr = 0
 Nmax = size(rho)-1

 ! Calculate dm grid
 allocate(dm(0:Nmax))
 do i = 0,Nmax
    dm(i) = m(i)-m(i-1)
 enddo

 do i = 0,Nmax-1
    ! Step forward in pres using hydrostatic equilibrium condition
    pres(i+1) = ( dm(i+1)**2 * pres(i-1) &
                - dm(i) * dm(i+1) * sum(dm(i:i+1)) &
                * gg * m(i) * one_on_4pi / r(i)**4 &
                - ( dm(i+1)**2 - dm(i)**2 ) * pres(i) ) &
                / dm(i)**2
    !pres(i+1) = pres(i) - 0.5*sum(dm(i:i+1)) * gg*m(i)*one_on_4pi/r(i)**4 ! 1st-order finite differencing used in MESA
    
    if (pres(i+1) < 0) then
       print*,'PRES<0 at i = ',i
       stop
    endif

    !call get_rho_from_p_s(pres(i+1),S(i+1),rho(i+1),rho(i),.false.) ! Newton-Raphson solver
    !call get_rho_from_p_s_bisection(pres(i+1),S(i+1),rho(i+1),rho(i),.false.) ! Bisection method solving for rho^0.5
    call get_rho_from_p_s_bisection_v2(pres(i+1),S(i+1),rho(i+1),rho(i),.false.) ! Bisection method solving for rho
 
    if (rho(i+1) < 0) then
      print*,'RHO<0 at i = ',i
      stop
    endif

    ! Step in r using continuity equation
    r(i+1) =    ( dm(i+1)**2 * r(i-1) &
                + dm(i) * dm(i+1) * sum(dm(i:i+1)) &
                * one_on_4pi / (r(i)**2 * rho(i)) &
                - ( dm(i+1)**2 - dm(i)**2 ) * r(i) ) &
                / dm(i)**2    
    !r(i+1) = (r(i)**3 + 3.*one_on_4pi * dm(i+1) / rho(i+1))**(1./3.) ! 1st-order finite differencing used in MESA

    if (r(i+1) < 0) then
      print*,'r<0 at i = ',i
      stop
   endif

 enddo

end subroutine solve_fixedS_surface_rhoP


!-----------------------------------------------------------------------
!+
!  Calculates specific entropy (gas + radiation + recombination)
!  up to an additive integration constant, from density and pressure.
!+
!-----------------------------------------------------------------------
function entropy(rho,pres,ierr)
 use physcon,           only:radconst,kb_on_mh
 use eos,               only:gmw,ieos
 use eos_idealplusrad,  only:get_idealgasplusrad_tempfrompres
 use eos_mesa,          only:get_eos_eT_from_rhop_mesa
 use mesa_microphysics, only:getvalue_mesa
 real, intent(in)               :: rho,pres
 integer, intent(out), optional :: ierr
 real                           :: inv_mu,entropy,logentropy,temp,eint

 if (present(ierr)) ierr=0
 inv_mu = 1/gmw

 select case(ieos)
 case(2) ! Include only gas entropy for adiabatic EoS
     temp = pres * gmw / (rho * kb_on_mh)
     entropy = kb_on_mh * inv_mu * log(temp**1.5/rho)

 case(12) ! Include both gas and radiation entropy gas plus rad. EoS
     temp = pres * gmw / (rho * kb_on_mh) ! Guess for temp
     call get_idealgasplusrad_tempfrompres(pres,rho,gmw,temp) ! First solve for temp from rho and pres
     entropy = kb_on_mh * inv_mu * log(temp**1.5/rho) + 4.*radconst*temp**3 / (3.*rho)
 
 case(10) ! MESA EoS
     call get_eos_eT_from_rhop_mesa(rho,pres,eint,temp)

     ! Get entropy from rho and eint from MESA tables
     if (present(ierr)) then
        call getvalue_mesa(rho,eint,9,logentropy,ierr)
     else
        call getvalue_mesa(rho,eint,9,logentropy)
     endif
     entropy = 10.d0**logentropy

 end select
  
end function entropy


!-----------------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using Newton-Raphson
!  method
!+
!-----------------------------------------------------------------------
subroutine get_rho_from_p_s(pres,S,rho,rhoguess,debug)
 use physcon, only:kb_on_mh
 real, intent(in)  :: pres,S,rhoguess
 real, intent(inout) :: rho
 real(kind=8)              :: corr,dSdsrho,S_plus_dS,srho_plus_dsrho,srho
 real, parameter   :: eoserr=1d-9,dfac=1d-12
 logical, intent(in)           :: debug
 ! We apply the Newton-Raphson method directly to rho^1/2 ("srho") instead
 ! of rho since S(rho) cannot take a negative argument.
 srho = sqrt(rhoguess) ! Initial guess
 corr = huge(corr);
 do while (abs(corr) > eoserr*abs(srho))
    ! First calculate dS/dsrho
    srho_plus_dsrho = srho * (1. + dfac)
    S_plus_dS = entropy(srho_plus_dsrho**2, pres)
    dSdsrho = (S_plus_dS - entropy(srho**2,pres)) / (srho_plus_dsrho - srho)
    corr = ( entropy(srho**2,pres) - S ) / dSdsrho
    srho = srho - corr
 enddo
 rho = srho**2
 return

end subroutine get_rho_from_p_s


!----------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using bisection method
!+
!----------------------------------------------------------------
subroutine get_rho_from_p_s_bisection(pres,S,rho,rhoguess,debug)
 use eos,     only:gmw
 use physcon, only:kb_on_mh
 real, intent(in)           :: pres,S
 real, intent(in), optional :: rhoguess
 real, intent(inout)        :: rho
 logical, intent(in)        :: debug
 real                       :: err,srho1,srho2,srhoguess,&
                               srho3,S1,S2,S3,left,right,mid
 real, parameter            :: tolerance = 1d-15

 if (present(rhoguess)) then
    srho1 = 1.005 * sqrt(rhoguess)  ! Tight lower bound
    srho2 = 0.995 * sqrt(rhoguess)  ! Tight upper bound
 else
    srhoguess = pres**0.2 * (gmw/kb_on_mh)**0.3 * exp(-0.5*gmw/kb_on_mh * S)
    srho1 = 10. * srhoguess  ! Guess lower bound
    srho2 = 0.1 * srhoguess  ! Guess upper bound
 endif

 S1 = entropy(srho1**2, pres)
 S2 = entropy(srho2**2, pres)
 left  = S - S1
 right = S - S2

 ! If lower and upper bounds do not contain roots, extend them until they do
 do while (left*right > 0.)
    srho1 = 0.99 * srho1
    srho2 = 1.01 * srho2
    S1 = entropy(srho1**2, pres)
    S2 = entropy(srho2**2, pres)
    left  = S - S1
    right = S - S2
 enddo

 ! Start bisecting
 err = huge(1.)
 do while (abs(err) > tolerance)
    S1 = entropy(srho1**2, pres)
    S2 = entropy(srho2**2, pres)
    left  = S - S1
    right = S - S2
    srho3 = 0.5*(srho1+srho2)
    S3 = entropy(srho3**2, pres)
    mid = S - S3

    if (left*mid < 0.) then
       srho2 = srho3
    elseif (right*mid < 0.) then
       srho1 = srho3
    elseif (mid == 0.) then
       rho = srho3**2
       exit
    endif

    rho = srho3**2
    err = (srho2 - srho1)/srho1
 enddo

end subroutine get_rho_from_p_s_bisection

!----------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using bisection method
!+
!----------------------------------------------------------------
subroutine get_rho_from_p_s_bisection_v2(pres,S,rho,rhoguess,debug)
 real, intent(in)           :: pres,S
 real, intent(in), optional :: rhoguess
 real, intent(inout)        :: rho
 logical, intent(in)        :: debug
 real                       :: err,rho1,rho2,&
                               rho3,S1,S2,S3,left,right,mid
 real, parameter            :: tolerance = 1d-15

 rho1 = 1.005 * rhoguess  ! Tight lower bound
 rho2 = 0.995 * rhoguess  ! Tight upper bound

 S1 = entropy(rho1, pres)
 S2 = entropy(rho2, pres)
 left  = S - S1
 right = S - S2

 ! If lower and upper bounds do not contain roots, extend them until they do
 do while (left*right > 0.)
    rho1 = 0.99 * rho1
    rho2 = 1.01 * rho2
    S1 = entropy(rho1, pres)
    S2 = entropy(rho2, pres)
    left  = S - S1
    right = S - S2
 enddo

 ! Start bisecting
 err = huge(1.)
 do while (abs(err) > tolerance)
    S1 = entropy(rho1, pres)
    S2 = entropy(rho2, pres)
    left  = S - S1
    right = S - S2
    rho3 = 0.5*(rho1+rho2)
    S3 = entropy(rho3, pres)
    mid = S - S3

    if (left*mid < 0.) then
       rho2 = rho3
    elseif (right*mid < 0.) then
       rho1 = rho3
    elseif (mid == 0.) then
       rho = rho3
       exit
    endif

    rho = rho3
    err = (rho2 - rho1)/rho1
 enddo

end subroutine get_rho_from_p_s_bisection_v2


!----------------------------------------------------------------
!+
!  Write out modified and original entropy profile for debugging
!+
!----------------------------------------------------------------
subroutine write_entropy(outputpath, S, Sorig)
 real, intent(in)                :: S(:),Sorig(:)
 character(len=*), intent(in)    :: outputpath
 integer                         :: i

 open(1, file = outputpath, status = 'replace')

 write(1,'(a)') '[  Entropy  ]  [EntropyOri ]'
 write(1,101) (S(i),Sorig(i),i=1,size(S))
 101 format (es13.6,2x,es13.6)
 close(1, status = 'keep')

end subroutine write_entropy

end module setfixedentropysurf