!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: timestep_sts
!
!  DESCRIPTION: This module contains the subroutines to control super-
!               timestepping; the primary control routine, step_STS, is
!               in step_supertimestep.f.  The control routine is in a
!               different file due to dependencies and hence compiling
!               order.
!               For independent timesteps, all particles requiring
!               super-timestepping will be in the largest bin.
!
!  REFERENCES: Alexiades V., Amiez G., Gremaud P.A., 1996, Commun. Numer. Meth. Eng., 12, 31
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, part, timestep_ind
!+
!--------------------------------------------------------------------------
module timestep_sts
 use dim, only: maxsts
 implicit none

 !--Control Variables (Hardcode values if not using STS)
#ifdef STS_TIMESTEPS
 logical,         public    :: use_sts      = .true.
 logical,         public    :: sts_it_n
#else
 logical,         public    :: use_sts      = .false.
 logical,         public    :: sts_it_n     = .true.
#endif
 !--Paramters
#ifdef STS_TIMESTEPS
 integer,         parameter :: nnu          =   12    ! The maximum number of supersteps allowed (stable, but slightly slower if 20)
#else
 integer,         parameter :: nnu          =    2    ! Default to avoid compiler warnings
#endif
 integer,         parameter :: dtcoef_max   =   18    ! The maximum number of dtdiff coefficients (retest if changed)
 integer,         parameter :: ndtau_max    =  256    ! The maximum number Nsts*Nmega allowed
 real,            private   :: nu_min       =  1.0d-4 ! The minimum allowed nu (retest if changed)
 real,            private   :: sts_max_rat  = 20.0    ! The maximum ratio of dtau(1) to dtdiff (retest if changed)
 real,            public    :: bigdt        =  1.0d29 ! =bignumber; duplicated here to avoid continually passing it
 !
 integer,         parameter :: iNosts       =  0      ! No super-timestepping required
 integer,         parameter :: iNsts        =  1      ! can use N \propto sqrt(dt/dtdiff)
 integer,         parameter :: iNmegasts    =  2      ! increased N using Nmegatseps
 integer,         parameter :: iNostsSml    =  3      ! using Nreal ~ dt/didiff (since Nsupersteps = Nreal for small Nreal)
 integer,         parameter :: iNostsBig    =  4      ! using Nreal since Nsts > nnu, or Nsts*Nmega > Nreal
 !
 integer(kind=1), parameter :: isactive_no  =  0      ! for individual timesteps: particle is inactive
 integer(kind=1), parameter :: isactive_yes =  1      ! for individual timesteps: particle is active as per normal
 integer(kind=1), parameter :: isactive_sts =  2      ! for individual timesteps: particle is active and requires STS
 !
 integer(kind=1), parameter :: ino          =  0      ! States that this dtau is in the middle of the sequence; needed if megasteps used
 integer(kind=1), parameter :: iyes         =  1      ! States that this dtau is at the beginning of the sequence; needed if megasteps used
 !
 logical,         private   :: print_nu_to_file = .false.  ! to allow nu to be printed for testing purposes
 !
 !--Variables
 integer(kind=1), private   :: istsactive(maxsts)
 integer(kind=1), public    :: ibin_sts(maxsts),isfirstdtau(ndtau_max)
 real,            public    :: dtau(ndtau_max),nu(nnu,dtcoef_max),dtdiffcoef(dtcoef_max)
 integer,         public    :: ipart_rhomax_sts,nbinmaxsts,Nmegasts_done,Nmegasts_now,Nmegasts_next
 integer,         public    :: Nreal,Nsts,icase_sts
 !
 !--Subroutines

 public                     :: sts_initialise,sts_init_step
 public                     :: sts_get_dtau_next,sts_get_dtau_array
 public                     :: sts_initialise_activity,sts_set_active_particles
 private                    :: sts_init_nu,sts_get_Ndtdiff,sts_get_dtdiff,sts_get_dtau
 private                    :: sts_update_i_nmega

contains

!----------------------------------------------------------------
!+
!  This will initialise the values for sts_check_dt.
!+
!----------------------------------------------------------------
subroutine sts_initialise(ierr,dtdiff)
 use io,             only: iuniteos
 integer, intent(out)   :: ierr
 real,    intent(inout) :: dtdiff
 integer                :: i
 !
 !--Initialise values
 ierr          = 0
 sts_it_n      = .true.
 dtdiff        = bigdt
 nbinmaxsts    = 0
 ibin_sts(:)   = 0
 istsactive(:) = isactive_yes

 do i = 1,dtcoef_max
    dtdiffcoef(i) = 0.9 - (i-1)*0.05
    call sts_init_nu(nu(1:nnu,i),dtdiffcoef(i),ierr)
    if (ierr > 0) return
 enddo

 ! Print N-nu tables to file, if requested
 if (print_nu_to_file) then
    open(unit=iuniteos,file="Nnu.dat",form='formatted',status='replace')
    write(iuniteos,'("# N vs nu values")')
    write(iuniteos,"('#',20(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'N',         &
        2,'coef=0.90', &
        3,'coef=0.85', &
        4,'coef=0.80', &
        5,'coef=0.75', &
        6,'coef=0.70', &
        7,'coef=0.65', &
        8,'coef=0.60', &
        9,'coef=0.55', &
       10,'coef=0.50', &
       11,'coef=0.45', &
       12,'coef=0.40', &
       13,'coef=0.35', &
       14,'coef=0.30', &
       15,'coef=0.25', &
       16,'coef=0.20', &
       17,'coef=0.15', &
       18,'coef=0.10', &
       19,'coef=0.05'
    do i = 1,nnu
       write(iuniteos,'(I18,1x,19(1pe18.10,1x))') i,nu(i,:)
       write(iuniteos+1,'(I18,1x,19(1pe18.10,1x))') i, 0.5*i/sqrt(nu(i,:))* &
                         ((1.0+sqrt(nu(i,:)))**(2*i)-(1.0-sqrt(nu(i,:)))**(2*i))/ &
                         ((1.0+sqrt(nu(i,:)))**(2*i)+(1.0-sqrt(nu(i,:)))**(2*i))
    enddo
    close(iuniteos)
 endif

end subroutine sts_initialise
!----------------------------------------------------------------
!+
!  This will precalculate nu for N = 1,nnu for a given dt & dtdiff,
!  where dt & dtdiff are related via
!     N = int(sqrt(dt/(dtdiffcoef*dtdiff)))+1
!     dtdiff_used = dt/(dtdiffcoef*N**2)
!
!  To get the values to govern the super-timestepping, we should solve
!     f(x) = dt - dtdiff*N/(2.0*sqrt(nu)) * (A-B)/(A+B)
!  where
!      A = (1.0+sqrt(nu))**2N
!      B = (1.0-sqrt(nu))**2N
!      dt,dtdiff are related as above
!  Since we are only interested in x such that f(x) = 0,
!  for stability, we will use
!      f(x) = 2*nu**1.5*(A+B)*dt - N*nu*dtdiff * (A-B)
!  We will use Newton's method to pre-calculate nu.
!+
!----------------------------------------------------------------
pure subroutine sts_init_nu(nu_col,dtdiffcoef_in,ierr)
 integer, intent(out) :: ierr
 real,    intent(in)  :: dtdiffcoef_in
 real,    intent(out) :: nu_col(:)
 integer, parameter   :: ctrmax       = 10000
 integer, parameter   :: iter_tol_max = 4
 real,    parameter   :: tol0         = 1.0e-14
 integer              :: i,ctr,N,twoN,iter_tol
 real                 :: f,A,B,dfdnu,dAdnu,dBdnu
 real                 :: nuold,nunew,nurat,tol,realN
 logical              :: getnu

 ! Initialise global values
 nu_col   = 1.0
 ierr     = 0

 ! Iterate to obtain the remaining nu-values
 do i = 2,nnu
    iter_tol = 0
    tol      = tol0
    getnu    = .true.
    do while ( getnu )
       ! Initialise values
       ctr   = 0
       N     = i
       twoN  = 2*N
       realN = real(N)
       if (nu_col(i-1) > 0.0) then
          nuold = nu_col(i-1)*0.999 ! refined first guess
       else
          nuold = nu_col(i)
       endif
       nurat = tol*2.0
       ! Iterate
       do while (ctr < ctrmax .and. nurat > tol)
          A     =                    (1.0d0+sqrt(nuold))**twoN
          B     =                    (1.0d0-sqrt(nuold))**twoN
          dAdnu =  realN/sqrt(nuold)*(1.0d0+sqrt(nuold))**(twoN-1)
          dBdnu = -realN/sqrt(nuold)*(1.0d0-sqrt(nuold))**(twoN-1)
          f     = 2.0d0*nuold**1.5*(A+B)*realN*dtdiffcoef_in - nuold*(A-B)
          dfdnu = realN*sqrt(nuold)*dtdiffcoef_in*(3.0d0*(A+B)+2.0d0*nuold*(dAdnu+dBdnu))   &
                - (A-B+nuold*(dAdnu-dBdnu))
          nunew = nuold - f/dfdnu
          nurat = abs( 1.0d0 - nunew/nuold )
          nuold = nunew
          ctr   = ctr + 1
          if (nunew > 1.0 .or. nunew < epsilon(nunew)) then
             ! if a solution does not exist within 0 < nu < 1, then set to 0
             nuold = 0.0
             nurat = 0.0
             getnu = .false.
          endif
       enddo
       if (ctr >= ctrmax .and. iter_tol < iter_tol_max) then
          tol      = tol * 10.0
          iter_tol = iter_tol + 1
       else
          getnu = .false.
       endif
    enddo
    if (nuold > 1.0 .or. nuold < 0.0 .or. iter_tol == iter_tol_max) then
       ierr = 1                       ! call error
       return
    elseif (nuold > nu_min) then
       nu_col(i) = nuold              ! save value
    else
       nu_col(i) = 0.0                ! solution is outside of 0 < nu < 1; set to 0 to be ignored
    endif
 enddo

end subroutine sts_init_nu
!------------------------------------------------------------
!+
!  After twas is initialised from init_step in step_leapfrog.F90,
!  overwrite twas for the particle what will undergo super-timestepping
!+
!------------------------------------------------------------
subroutine sts_init_step(npart,timei,dtmax,dtau_in)
#ifdef IND_TIMESTEPS
 use timestep_ind, only: get_dt
 use part,         only: ibin,twas
#endif
 integer,         intent(in)  :: npart
 real,            intent(in)  :: timei,dtau_in,dtmax
#ifdef IND_TIMESTEPS
 integer                      :: i

!$omp parallel default(none) &
!$omp shared(npart,ibin,ibin_sts,twas,timei,dtau_in,dtmax) &
!$omp private(i)
!$omp do schedule(runtime)
 do i=1,npart
    if (ibin_sts(i) > ibin(i)) then
       twas(i) = timei + 0.5*dtau_in
    endif
 enddo
!$omp enddo
!$omp end parallel
#endif

end subroutine sts_init_step
!----------------------------------------------------------------
!+
!  Calculate the next dt value, dtau_next
!+
!----------------------------------------------------------------
subroutine sts_get_dtau_next(dtau_next,dt_in,dtmax,dtdiff_in,nbinmax)
 integer(kind=1), intent(in),  optional :: nbinmax
 real,            intent(out)           :: dtau_next
 real,            intent(in)            :: dt_in,dtmax,dtdiff_in
 real                                   :: dt_next,dt_remain
 integer                                :: i,Nmega_tmp,Nms_p1

 ! For individual dt, nbinmax may change during our superstepping, thus
 ! always compare against this.
 if (present(nbinmax)) then
    dt_next = dtmax/2**nbinmax    ! individual timesteps
 else
    dt_next = dt_in               ! global timesteps
 endif

 if ( sts_it_n ) then
    ! We are currently on the final super-timestep, thus need to predict
    ! what will happen next loop so that we can calculate its dtau_next
    Nmegasts_done = 0 ! Since we are resetting the dtau array
    call sts_get_dtau_array(Nmegasts_next,dt_next,dtdiff_in)
    dtau_next = dtau(1)
 else
    Nms_p1 = Nmegasts_done+1
    dtau_next = dtau(Nms_p1)
    if (isfirstdtau(Nms_p1)==iyes .and. dtau_next > dt_next) then
       write(*,'(a)') "Super-timestep: sts_get_dtau_next: dtau_next > dt_next thus modifying mid-mega-step"
       dt_remain = 0.
       do i = Nms_p1,Nmegasts_now
          dt_remain = dt_remain + dtau(i)
       enddo
       Nmega_tmp = 0
       do while (dtau_next > dt_next)
          Nmega_tmp = Nmega_tmp + 1
          call sts_get_dtau_array(Nmegasts_now,dt_remain,dtdiff_in,Nmega_tmp)
          dtau_next = dtau(Nms_p1)
       enddo
    endif
 endif

end subroutine sts_get_dtau_next
!----------------------------------------------------------------
!+
!  Calculates all the dtau's require for this Nsts*Nmega; can
!  revise on the fly if necessary
!+
!----------------------------------------------------------------
subroutine sts_get_dtau_array(Nmegasts,dt_next,dtdiff_in,Nmega_in)
 use io, only: fatal
 real,    intent(in)           :: dt_next,dtdiff_in
 integer, intent(out)          :: Nmegasts
 integer, intent(in), optional :: Nmega_in
 real                          :: dtdiff_used,nu
 integer                       :: i,j,k,Nmega
 logical                       :: calc_dtau

 nu          = 0.2       ! to avoid compiler warnings
 dtdiff_used = dtdiff_in ! to avoid compiler warnings

 ! Determine the number of real steps required;
 ! if Nmegasts_done > 0, then this is the real steps remaining
 Nreal = int(dt_next/dtdiff_in) + 1

 ! Calculate the number of super and mega steps
 if (present(Nmega_in)) then
    call sts_get_Ndtdiff(dt_next/float(Nmega_in),dtdiff_in,dtdiff_used,Nsts,Nmega,nu,Nreal,icase_sts)
    Nmega = Nmega * Nmega_in
 else
    call sts_get_Ndtdiff(dt_next,dtdiff_in,dtdiff_used,Nsts,Nmega,nu,Nreal,icase_sts)
 endif
 Nmegasts = Nmegasts_done + Nsts*Nmega

 ! set cases
 if ( Nmegasts > ndtau_max) then
    call fatal('sts_get_dtau_array','dtau array is too small',var="Nmegasts",ival=Nmegasts)
 endif
 if (icase_sts==iNosts) then
    dtau(1)        = dt_next
    isfirstdtau(1) = iyes
    return
 elseif (icase_sts==iNsts .or. icase_sts==iNmegasts) then
    calc_dtau = .true.
 else ! includes iNostsSml, iNostsBig & iNosts (where Nreal = 1)
    calc_dtau = .false.
 endif

 ! set the time array and note which times are the first in the set of supersteps
 k = Nmegasts_done
 do j = 1,Nmega
    do i = 1,Nsts
       k = k + 1
       if (calc_dtau) then
          dtau(k) = sts_get_dtau(i,Nsts,nu,dtdiff_used)
          if (i==1) then
             isfirstdtau(k) = iyes
          else
             isfirstdtau(k) = ino
          endif
       else
          dtau(k)        = dt_next/real(Nreal)
          isfirstdtau(k) = iyes
       endif
    enddo
 enddo

end subroutine sts_get_dtau_array
!----------------------------------------------------------------
!+
!  Calculates N and the revised dtdiff
!+
!----------------------------------------------------------------
subroutine sts_get_Ndtdiff(dt,dtdiff_in,dtdiff_out,Nsts,Nmega,nu_local,Nreal,icase)
 use io, only: fatal
 real,    intent(in)    :: dt,dtdiff_in
 real,    intent(out)   :: dtdiff_out,nu_local
 integer, intent(in)    :: Nreal
 integer, intent(out)   :: Nsts,icase,Nmega
 integer                :: i
 real                   :: dtau_local
 logical                :: find_dtdiff

 ! Determine values for super-timestepping
 if (dt > dtdiff_in .and. dtdiff_in > tiny(dtdiff_in) .and. dtdiff_in < bigdt) then
    find_dtdiff = .true.
    i           = 1
    Nmega       = 1
    do while (find_dtdiff)
       Nsts = int( sqrt(dt/(dtdiffcoef(i)*dtdiff_in*Nmega)) ) + 1
       if (Nsts*Nmega >= Nreal) find_dtdiff = .false.
       if (find_dtdiff) then
          dtdiff_out = sts_get_dtdiff(i,dt/float(Nmega),Nsts)
          if (Nsts < nnu) then
             nu_local   = nu(Nsts,i)
             dtau_local = sts_get_dtau(1,Nsts,nu_local,dtdiff_out)
             ! if ratio is too big, or nu is too small, try again with different values
             if (dtau_local > dtdiff_in*sts_max_rat .or. nu_local < nu_min) then ! tests of an extreme case showed this is required
                call sts_update_i_nmega(i,Nmega)
             else
                find_dtdiff = .false.
             endif
          else
             call sts_update_i_nmega(i,Nmega)
          endif
       endif
    enddo

    ! Set the icase number & modify Nsts & Nmega as required
    if (Nmega==1) then
       if (Nsts == 1) then
          icase = iNosts
       else if (Nsts==Nreal) then
          icase = iNostsSml
       else if (Nsts > nnu .or. Nsts > Nreal) then
          icase = iNostsBig
          Nsts  = Nreal
       else
          icase = iNsts
       endif
    else if (Nmega > 1) then
       if (Nsts == 1 .or. Nsts*Nmega >= Nreal) then
          icase = iNostsBig
          Nmega = 1
          Nsts  = Nreal
       else
          icase = iNmegasts
       endif
    endif
 else
    Nmega      = 1
    Nsts       = 1
    dtdiff_out = bigdt
    icase      = iNosts
 endif

end subroutine sts_get_Ndtdiff
!
subroutine sts_update_i_nmega(i,Nmega)
 integer, intent(inout) :: i,Nmega
 i = i + 1
 if (i > dtcoef_max) then
    Nmega = Nmega + 1
    i     = 1
 endif
end subroutine sts_update_i_nmega
!----------------------------------------------------------------
!+
!  Calculate the revised diffusive timestep only
!+
!----------------------------------------------------------------
pure function sts_get_dtdiff(i,dt,N)
 real,    intent(in)  :: dt
 integer, intent(in)  :: i,N
 real                 :: sts_get_dtdiff

 sts_get_dtdiff =  dt/(dtdiffcoef(i)*real(N)**2)

end function sts_get_dtdiff
!----------------------------------------------------------------
!+
!  Calculate the timestep
!+
!----------------------------------------------------------------
pure function sts_get_dtau(j,N,nu0,dtdiff_in)
 integer, intent(in)  :: j,N
 real,    intent(in)  :: nu0,dtdiff_in
 real                 :: sts_get_dtau,pibytwo

 pibytwo      = 1.5707963268d0
 sts_get_dtau = dtdiff_in /((nu0-1.0d0)*cos(pibytwo*real(2*j-1)/real(N)) + 1.0d0+nu0)

end function sts_get_dtau
!----------------------------------------------------------------
!+
!  Set active particle list to distinguish between
!  inactive particles, (normal) active particle, and active
!  particles requiring sts; this will prevent confusion of
!  particles switching ibin's during while super-timestepping
!+
!----------------------------------------------------------------
subroutine sts_initialise_activity(nactive_sts,npart,ibin,iphase)
 integer,         intent(in)  :: npart
 integer(kind=1), intent(in)  :: ibin(:),iphase(:)
 integer        , intent(out) :: nactive_sts
 integer                      :: i

 nactive_sts = 0
!$omp parallel default(none) &
!$omp shared(npart,ibin,ibin_sts,iphase,istsactive) &
!$omp private(i) &
!$omp reduction(+:nactive_sts)
!$omp do
 do i = 1,npart
    if (iphase(i) > 0) then
       ! particle is active
       if (ibin_sts(i) > ibin(i)) then
          istsactive(i) = isactive_sts
          nactive_sts   = nactive_sts + 1
       else
          istsactive(i) = isactive_yes
       endif
    else
       istsactive(i) = isactive_no
    endif
 enddo
!$omp enddo
!$omp end parallel

end subroutine sts_initialise_activity
!----------------------------------------------------------------
!+
!  routine to set iactive flag determining whether a particle's
!  forces are to be evaluated on the current timestep.
!  A modified version of set_active_particles
!+
!----------------------------------------------------------------
subroutine sts_set_active_particles(npart,nactive,all_active)
 use part, only: isdead_or_accreted,iamtype,isetphase,iphase,xyzh
 integer,         intent(in)    :: npart
 integer,         intent(out)   :: nactive
 logical,         intent(in)    :: all_active
 integer                        :: i,itype
 integer(kind=1)                :: isactive_opt2
 logical                        :: iactivei

 nactive = 0
 if (all_active) then
    isactive_opt2 = isactive_yes
 else
    isactive_opt2 = isactive_sts
 endif
!$omp parallel default(none) &
!$omp shared(npart,iphase,xyzh,istsactive,isactive_opt2) &
!$omp private(i,itype,iactivei) &
!$omp reduction(+:nactive)
!$omp do
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (istsactive(i)==isactive_sts .or. istsactive(i)==isactive_opt2) then
          iactivei = .true.
          nactive  = nactive + 1
       else
          iactivei = .false.
       endif
       itype     = iamtype(iphase(i))
       iphase(i) = isetphase(itype,iactivei)
    endif
 enddo
!$omp enddo
!$omp end parallel

end subroutine sts_set_active_particles
!----------------------------------------------------------------
end module timestep_sts
