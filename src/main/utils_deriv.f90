!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module derivutils
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: io, kernel, mpiutils, options, timing
!
 use timing, only:timers,itimer_dens,itimer_force,itimer_tree,itimer_balance,itimer_cons2prim,&
                   itimer_radiation,itimer_rad_save,itimer_rad_neighlist,itimer_rad_arrays,itimer_rad_its,&
                   itimer_rad_flux,itimer_rad_diff,itimer_rad_update,itimer_rad_store,itimer_HII

 implicit none

 private

 public :: do_timing,get_grkern_corrected

contains

!-------------------------------------------------------------
!+
!  interface to timing routines
!+
!-------------------------------------------------------------
subroutine do_timing(label,tlast,tcpulast,start,lunit)
 use io,       only:iverbose,id,master
 use mpiutils, only:reduce_mpi
 use timing,   only:increment_timer,log_timing,get_timings
 character(len=*), intent(in)    :: label
 real(kind=4),     intent(inout) :: tlast,tcpulast
 logical,          intent(in), optional :: start
 integer,          intent(in), optional :: lunit
 integer                                :: itimer
 real(kind=4) :: t2,tcpu2,tcpu,twall

 call get_timings(t2,tcpu2)
 twall = reduce_mpi('+',t2-tlast)
 tcpu  = reduce_mpi('+',tcpu2-tcpulast)

 select case (label)
 case ('dens')
    itimer = itimer_dens
 case ('force')
    itimer = itimer_force
 case ('tree')
    itimer = itimer_tree
 case ('cons2prim')
    itimer = itimer_cons2prim
 case ('radiation')
    itimer = itimer_radiation
 case ('radsave')
    itimer = itimer_rad_save
 case ('radneighlist')
    itimer = itimer_rad_neighlist
 case ('radarrays')
    itimer = itimer_rad_arrays
 case ('radits')
    itimer = itimer_rad_its
 case ('radflux')
    itimer = itimer_rad_flux
 case ('raddiff')
    itimer = itimer_rad_diff
 case ('radupdate')
    itimer = itimer_rad_update
 case ('radstore')
    itimer = itimer_rad_store
 case ('HII_region')
    itimer = itimer_HII
 case default
    itimer = 1
 end select

 call increment_timer(itimer,t2-tlast,tcpu2-tcpulast)

 if (iverbose >= 2 .and. id==master) then
    if (present(start)) then
       call log_timing(label,t2-tlast,tcpu,start=.true.)
    elseif (present(lunit)) then
       call log_timing(label,t2-tlast,tcpu,iunit=lunit)
    else
       call log_timing(label,t2-tlast,tcpu)
    endif
 endif
 tlast = t2
 tcpulast = tcpu2

end subroutine do_timing

!-------------------------------------------------------------
!+
!  Kernel gradient with the smoothing-length correction for h set from
!  number density: grad W + zeta/(Omega~ m_j) grad W~ (with grad W~ = grad W
!  unless two_kernel), so that m_j*grkern_corr is the term in d(rho_i)/dr_i.
!  h41 is 1/h^4; zetaomega is zeta_i/Omega~_i; pmassj is the neighbour mass.
!+
!-------------------------------------------------------------
pure subroutine get_grkern_corrected(q2,q,h41,zetaomega,pmassj,grkern_corr,grkern_tilde)
 use kernel,  only:grkern,cnormk,cnormk_tilde,get_kernel_tilde
 use options, only:two_kernel
 real, intent(in)  :: q2,q,h41,zetaomega,pmassj
 real, intent(out) :: grkern_corr,grkern_tilde
 real :: wtilde

 grkern_corr = grkern(q2,q)*h41*cnormk
 if (two_kernel) then
    call get_kernel_tilde(q2,q,wtilde,grkern_tilde)
    grkern_tilde = grkern_tilde*h41*cnormk_tilde
 else
    grkern_tilde = grkern_corr
 endif
 if (pmassj > tiny(pmassj)) grkern_corr = grkern_corr + (zetaomega/pmassj)*grkern_tilde

end subroutine get_grkern_corrected

end module derivutils
