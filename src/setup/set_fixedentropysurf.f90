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
 integer :: ientropy

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calculates the constant entropy softened surface
!  profile. Assumes input profile quantities are sorted from centre to
!  surface. 
!+
!-----------------------------------------------------------------------
subroutine set_fixedS_surface(mcore,m,rho,r,pres,ene,temp,ierr)
 use io,          only:fatal
 use eos,         only:calc_temp_and_ene,entropy,ieos
 use physcon,     only:pi,gg,solarm,solarr,kb_on_mh
 use table_utils, only:interpolator,flip_array
 real, intent(in)    :: mcore
 real, intent(inout) :: m(:),rho(:),r(:),pres(:),ene(:),temp(:)
 real                :: eneguess,mcore_cm!,dum
 real                :: msurf,dm(-1:0),dS_by_dm
 real, allocatable   :: S(:),m_alloc(:),r_alloc(:),rho_alloc(:),pres_alloc(:),&
                        recalc_S(:),testentropy(:),entropyorig(:),idealplusradtemp(:),mesatemp(:)
 integer             :: j,i,surfidx,Nmax,ierr
 logical             :: isort_decreasing,iexclude_core_mass
 ! msurf: Mass coordinate beyond which we perform the profile replacement

 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .true.     ! Needs to be true if to be read by Phantom
 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .true.   ! Needs to be true if to be read by Phantom

 msurf = 9.8 * solarm !11.35 * solarm ! Mass coordinate beyond which we flatten the entropy
 call interpolator(m,msurf,surfidx) ! Find index closest to msurf
 Nmax = size(r) - surfidx

 select case(ieos)
 case(2)
   ientropy = 1
 case(12)
   ientropy = 2
 case(10)
   ientropy = 2
 case default
   call fatal('setfixedentropycore','ieos not one of 2 (adiabatic), 12 (ideal plus rad.), or 10 (MESA)')
 end select

 ! Only need dm and S three steps behind the surface to calculate dS/dm there
 ! Weird indexing here so that the index 0 is the surfidx
 dm(-1) = m(surfidx-1) - m(surfidx-2)
 dm(0)  = m(surfidx) - m(surfidx-1)
 allocate(S(-2:0)) 
 do j = -2,0
    S(j) = entropy(rho(surfidx+j),pres(surfidx+j),ientropy)
 enddo

 ! Find entropy gradient at m(surfidx-1)
 dS_by_dm = ( dm(-1)**2 * S(0) &
            + ( dm(0)**2 - dm(-1)**2 ) * S(-1) &
            - dm(0)**2 * S(-2) ) &
            / ( dm(-1) * dm(0) * (dm(-1) + dm(0)) )

 ! Fill up entropy grid using the entropy gradient   
!  deallocate(S)     
!  allocate(S(0:Nmax)) ! Repurpose variable name "S" for full surface entropy grid
!  S(0) = entropy(rho(surfidx),pres(surfidx),ientropy)
!  do i = 1,Nmax
!     S(i) = S(0) + dS_by_dm * (m(surfidx+i)-m(surfidx))
!  enddo
 
 ! TEST: OUTPUT ENTROPY TO DATA FILE
!  allocate(testentropy(size(r)),entropyorig(size(r)))
!  do i = 1,size(r)
!     entropyorig(i) = entropy(rho(i),pres(i),ientropy)
!  enddo
!  testentropy(1:surfidx-1) = entropyorig(1:surfidx-1)
!  testentropy(surfidx:size(r)) = S(0:Nmax)

 ! TEST: SUPPLY ORIGINAL ENTROPY TO CHECK ORIGINAL RHO AND PRES RECOVERED
 !S = entropyorig(surfidx:size(r))
 !testentropy(surfidx:size(r)) = entropyorig(surfidx:size(r))
 ! TEST: SUPPLY ORIGINAL ENTROPY TO CHECK ORIGINAL RHO AND PRES RECOVERED
!  allocate(recalc_S(size(r)))
!  do i = 1,size(r)
!    recalc_S(i) = entropy(rho(i),pres(i),ientropy)
!  enddo
 !call write_entropy('./entropy_mesa.dat', m, testentropy, entropyorig, testentropy)!recalc_S)

 ! TEST: OUTPUT ENTROPY TO DATA FILE
 ! TEST

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

 ! TEST: PRINT RECALCULATED TEMPERATURES
!  allocate(idealplusradtemp(size(r)),mesatemp(size(r)))
!  do i = 1,size(r)
!     idealplusradtemp(i) = pres(i) * gmw / (rho(i) * kb_on_mh) ! Guess for temp
!     call get_idealgasplusrad_tempfrompres(pres(i),rho(i),gmw,idealplusradtemp(i))
!     call get_eos_eT_from_rhop_mesa(rho(i),pres(i),dum,mesatemp(i))
!  enddo
!  call write_temp('./temp_comparison.dat', m, idealplusradtemp, mesatemp)
 ! TEST: PRINT RECALCULATED TEMPERATURES

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
 use io, only:fatal
 use physcon, only:gg
 use eos, only:get_rho_from_p_s
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
    if (pres(i+1) < 0) call fatal('setfixedentropysurf','pres < 0 found')

    call get_rho_from_p_s(pres(i+1),S(i+1),rho(i+1),rho(i),ientropy) ! Newton-Raphson solver
    if (rho(i+1) < 0) call fatal('setfixedentropysurf','rho < 0 found')

    ! Step in r using continuity equation
   !  r(i+1) =    ( dm(i+1)**2 * r(i-1) &
   !              + dm(i) * dm(i+1) * sum(dm(i:i+1)) &
   !              * one_on_4pi / (r(i)**2 * rho(i)) &
   !              - ( dm(i+1)**2 - dm(i)**2 ) * r(i) ) &
   !              / dm(i)**2
    r(i+1) = exp( log(r(i)**3 + 3.*one_on_4pi * dm(i+1) / rho(i+1) ) / 3.) ! 1st-order finite differencing used in MESA
    if (r(i+1) < 0) call fatal('setfixedentropysurf','r < 0 found')

 enddo

end subroutine solve_fixedS_surface_rhoP


!----------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using bisection method
!+
!----------------------------------------------------------------
! subroutine get_rho_from_p_s_bisection(pres,S,rho,rhoguess,debug)
!  use eos,     only:gmw,entropy
!  use physcon, only:kb_on_mh
!  real, intent(in)           :: pres,S
!  real, intent(in), optional :: rhoguess
!  real, intent(inout)        :: rho
!  logical, intent(in)        :: debug
!  real                       :: err,srho1,srho2,srhoguess,&
!                                srho3,S1,S2,S3,left,right,mid
!  real, parameter            :: tolerance = 1d-15

!  if (present(rhoguess)) then
!     srho1 = 1.005 * sqrt(rhoguess)  ! Tight lower bound
!     srho2 = 0.995 * sqrt(rhoguess)  ! Tight upper bound
!  else
!     srhoguess = pres**0.2 * (gmw/kb_on_mh)**0.3 * exp(-0.5*gmw/kb_on_mh * S)
!     srho1 = 10. * srhoguess  ! Guess lower bound
!     srho2 = 0.1 * srhoguess  ! Guess upper bound
!  endif

!  S1 = entropy(srho1**2, pres, ientropy)
!  S2 = entropy(srho2**2, pres, ientropy)
!  left  = S - S1
!  right = S - S2

!  ! If lower and upper bounds do not contain roots, extend them until they do
!  do while (left*right > 0.)
!     srho1 = 0.99 * srho1
!     srho2 = 1.01 * srho2
!     S1 = entropy(srho1**2, pres, ientropy)
!     S2 = entropy(srho2**2, pres, ientropy)
!     left  = S - S1
!     right = S - S2
!  enddo

!  ! Start bisecting
!  err = huge(1.)
!  do while (abs(err) > tolerance)
!     S1 = entropy(srho1**2, pres, ientropy)
!     S2 = entropy(srho2**2, pres, ientropy)
!     left  = S - S1
!     right = S - S2
!     srho3 = 0.5*(srho1+srho2)
!     S3 = entropy(srho3**2, pres, ientropy)
!     mid = S - S3

!     if (left*mid < 0.) then
!        srho2 = srho3
!     elseif (right*mid < 0.) then
!        srho1 = srho3
!     elseif (mid == 0.) then
!        rho = srho3**2
!        exit
!     endif

!     rho = srho3**2
!     err = (srho2 - srho1)/srho1
!  enddo

! end subroutine get_rho_from_p_s_bisection

!----------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using bisection method
!+
!----------------------------------------------------------------
! subroutine get_rho_from_p_s_bisection_v2(pres,S,rho,rhoguess,debug)
!  use eos, only:entropy
!  real, intent(in)           :: pres,S
!  real, intent(in), optional :: rhoguess
!  real, intent(inout)        :: rho
!  logical, intent(in)        :: debug
!  real                       :: err,rho1,rho2,&
!                                rho3,S1,S2,S3,left,right,mid
!  real, parameter            :: tolerance = 1d-15

!  rho1 = 1.005 * rhoguess  ! Tight lower bound
!  rho2 = 0.995 * rhoguess  ! Tight upper bound

!  S1 = entropy(rho1, pres, ientropy)
!  S2 = entropy(rho2, pres, ientropy)
!  left  = S - S1
!  right = S - S2

!  ! If lower and upper bounds do not contain roots, extend them until they do
!  do while (left*right > 0.)
!     rho1 = 0.99 * rho1
!     rho2 = 1.01 * rho2
!     S1 = entropy(rho1, pres, ientropy)
!     S2 = entropy(rho2, pres, ientropy)
!     left  = S - S1
!     right = S - S2
!  enddo

!  ! Start bisecting
!  err = huge(1.)
!  do while (abs(err) > tolerance)
!     S1 = entropy(rho1, pres, ientropy)
!     S2 = entropy(rho2, pres, ientropy)
!     left  = S - S1
!     right = S - S2
!     rho3 = 0.5*(rho1+rho2)
!     S3 = entropy(rho3, pres, ientropy)
!     mid = S - S3

!     if (left*mid < 0.) then
!        rho2 = rho3
!     elseif (right*mid < 0.) then
!        rho1 = rho3
!     elseif (mid == 0.) then
!        rho = rho3
!        exit
!     endif

!     rho = rho3
!     err = (rho2 - rho1)/rho1
!  enddo

! end subroutine get_rho_from_p_s_bisection_v2


!----------------------------------------------------------------
!+
!  Write out modified and original entropy profile for debugging
!+
!----------------------------------------------------------------
subroutine write_entropy(outputpath,m,S,Sorig,recalc_S)
 real, intent(in)                :: S(:),Sorig(:),m(:),recalc_S(:)
 character(len=*), intent(in)    :: outputpath
 integer                         :: i

 open(1, file = outputpath, status = 'replace')

 write(1,'(a)') '[   Mass    ]  [  Entropy  ]  [EntropyOri ]  [RecalcS    ]'
 write(1,101) (m(i),S(i),Sorig(i),recalc_S(i),i=1,size(S))
 101 format (es13.6,2x,es13.6,2x,es13.6,2x,es13.6)
 close(1, status = 'keep')

end subroutine write_entropy

!----------------------------------------------------------------
!+
!  Write out modified and original temperature profiles for debugging
!+
!----------------------------------------------------------------
subroutine write_temp(outputpath,m,idealplusradtemp,mesatemp)
   real, intent(in)                :: idealplusradtemp(:),mesatemp(:),m(:)
   character(len=*), intent(in)    :: outputpath
   integer                         :: i
  
   open(1, file = outputpath, status = 'replace')
  
   write(1,'(a)') '[   Mass    ]  [  temp1    ]  [  temp2    ]'
   write(1,101) (m(i),idealplusradtemp(i),mesatemp(i),i=1,size(mesatemp))
   101 format (es13.6,2x,es13.6,2x,es13.6)
   close(1, status = 'keep')
  
end subroutine write_temp

end module setfixedentropysurf