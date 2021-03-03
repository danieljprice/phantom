!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cons2prim
!
! None
!
! :References: None
!
! :Owner: Elisabeth Borchert
!
! :Runtime parameters: None
!
! :Dependencies: cons2primsolver, cullendehnen, dim, eos, io, nicil,
!   options, part, radiation_utils, units, utils_gr
!
 use cons2primsolver, only:ien_entropy
 implicit none

 public :: cons2primall,cons2prim_everything
 public :: prim2consall,prim2consi

 private

contains

!-------------------------------------
!
!  Primitive to conservative routines
!
!-------------------------------------

subroutine prim2consall(npart,xyzh,metrics,vxyzu,dens,pxyzu,use_dens)
 use part,         only:isdead_or_accreted
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:),metrics(:,:,:,:),vxyzu(:,:)
 real,    intent(inout) :: dens(:)
 real,    intent(out) :: pxyzu(:,:)
 logical, intent(in), optional :: use_dens
 logical :: usedens
 integer :: i

!  By default, use the smoothing length to compute primitive density, and then compute the conserved variables.
!  (Alternatively, use the provided primitive density to compute conserved variables.
!   Depends whether you have prim dens prior or not.)
 if (present(use_dens)) then
    usedens = use_dens
 else
    usedens = .false.
 endif

!$omp parallel do default (none) &
!$omp shared(xyzh,metrics,vxyzu,dens,pxyzu,npart,usedens) &
!$omp private(i)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call prim2consi(xyzh(:,i),metrics(:,:,:,i),vxyzu(:,i),dens(i),pxyzu(:,i),usedens)
    endif
 enddo
!$omp end parallel do

end subroutine prim2consall

subroutine prim2consi(xyzhi,metrici,vxyzui,dens_i,pxyzui,use_dens)
 use cons2primsolver, only:primitive2conservative
 use utils_gr,        only:h2dens
 use eos,             only:equationofstate,ieos,gamma
 real, dimension(4), intent(in)  :: xyzhi, vxyzui
 real,               intent(in)  :: metrici(:,:,:)
 real, intent(inout)             :: dens_i
 real, dimension(4), intent(out) :: pxyzui
 logical, intent(in), optional   :: use_dens
 logical :: usedens
 real :: rhoi,Pi,ui,xyzi(1:3),vi(1:3),pondensi,spsoundi,densi

 !  By default, use the smoothing length to compute primitive density, and then compute the conserved variables.
 !  (Alternatively, use the provided primitive density to compute conserved variables.
 !   Depends whether you have prim dens prior or not.)
 if (present(use_dens)) then
    usedens = use_dens
 else
    usedens = .false.
 endif

 xyzi = xyzhi(1:3)
 vi   = vxyzui(1:3)
 ui   = vxyzui(4)
 if (usedens) then
    densi = dens_i
 else
    call h2dens(densi,xyzhi,metrici,vi) ! Compute dens from h
    dens_i = densi                      ! Feed the newly computed dens back out of the routine
 endif
 call equationofstate(ieos,pondensi,spsoundi,densi,xyzi(1),xyzi(2),xyzi(3),ui)
 pi = pondensi*densi
 call primitive2conservative(xyzi,metrici,vi,densi,ui,Pi,rhoi,pxyzui(1:3),pxyzui(4),ien_entropy,gamma)

end subroutine prim2consi

!-------------------------------------
!
!  Conservative to primitive routines
!
!-------------------------------------

subroutine cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
 use cons2primsolver, only:conservative2primitive
 use part,            only:isdead_or_accreted,massoftype,igas,rhoh,igasP,ics
 use io,              only:fatal
 use eos,             only:equationofstate,ieos,gamma,done_init_eos,init_eos
 integer, intent(in)    :: npart
 real,    intent(in)    :: pxyzu(:,:),xyzh(:,:),metrics(:,:,:,:)
 real,    intent(inout) :: vxyzu(:,:),dens(:)
 real,    intent(out)   :: eos_vars(:,:)
 integer :: i, ierr
 real    :: p_guess,rhoi,pondens,spsound

 if (.not.done_init_eos) call init_eos(ieos,ierr)

!$omp parallel do default (none) &
!$omp shared(xyzh,metrics,vxyzu,dens,pxyzu,npart,massoftype) &
!$omp shared(ieos,gamma,eos_vars) &
!$omp private(i,ierr,spsound,pondens,p_guess,rhoi)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       ! Construct a guess for pressure (dens is already passed in and is also a guess coming in, but correct value gets passed out)
       call equationofstate(ieos,pondens,spsound,dens(i),xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))
       p_guess = pondens*dens(i)
       rhoi    = rhoh(xyzh(4,i),massoftype(igas))
       call conservative2primitive(xyzh(1:3,i),metrics(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i), &
                                  p_guess,rhoi,pxyzu(1:3,i),pxyzu(4,i),ierr,ien_entropy,gamma)
       eos_vars(igasP,i)     = p_guess
       eos_vars(ics,i)       = spsound
       if (ierr > 0) then
          print*,' pmom =',pxyzu(1:3,i)
          print*,' rho* =',rhoh(xyzh(4,i),massoftype(igas))
          print*,' en   =',pxyzu(4,i)
          call fatal('cons2prim','could not solve rootfinding',i)
       endif
    endif
 enddo
!$omp end parallel do

end subroutine cons2primall

!-------------------------------------
!
!  Primitive variables from conservative variables
!
!-------------------------------------

subroutine cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,&
                                gamma_chem,Bevol,Bxyz,dustevol,dustfrac,alphaind)
 use part,              only:isdead_or_accreted,massoftype,igas,rhoh,igasP,iradP,iradxi,ics,&
                             iohm,ihall,n_R,n_electronT,eta_nimhd,iambi,get_partinfo,iphase,this_is_a_test,&
                             ndustsmall,itemp,ikappa
 use eos,               only:equationofstate,ieos,gamma,get_temperature,done_init_eos,init_eos
 use radiation_utils,   only:radiation_equation_of_state,get_opacity
 use dim,               only:store_temperature,store_gamma,mhd,maxvxyzu,maxphase,maxp,use_dustgrowth,&
                             do_radiation,nalpha,mhd_nonideal
 use nicil,             only:nicil_get_ion_n,nicil_get_eta,nicil_translate_error
 use io,                only:fatal,real4
 use cullendehnen,      only:get_alphaloc,xi_limiter
 use options,           only:alpha,alphamax,use_dustfrac,iopacity_type
 use units,             only:unit_density,unit_opacity

 integer,      intent(in)    :: npart
 real,         intent(in)    :: xyzh(:,:),rad(:,:),gamma_chem(:),Bevol(:,:),dustevol(:,:)
 real(kind=4), intent(in)    :: dvdx(:,:)
 real,         intent(inout) :: vxyzu(:,:)
 real(kind=4), intent(inout) :: alphaind(:,:)
 real,         intent(out)   :: eos_vars(:,:),radprop(:,:),Bxyz(:,:),dustfrac(:,:)
 integer      :: i,ierr
 real         :: rhoi,pondens,spsound,p_on_rhogas,rhogas,gasfrac,pmassi
 real         :: Bxi,Byi,Bzi,psii,xi_limiteri,Bi,temperaturei
 real         :: xi,yi,zi,hi
 integer      :: iamtypei
 logical      :: iactivei,iamgasi,iamdusti

 iactivei = .true.
 iamtypei = igas
 iamgasi  = .true.
 iamdusti = .false.
 if (.not.done_init_eos) then
    call init_eos(ieos,ierr)
    if (ierr /= 0) call fatal('eos','could not initialise equation of state')
 endif

!$omp parallel do default (none) &
!$omp shared(xyzh,vxyzu,npart,rad,eos_vars,radprop,Bevol,Bxyz) &
!$omp shared(ieos,gamma,gamma_chem,n_R,n_electronT,eta_nimhd) &
!$omp shared(alpha,alphamax,iphase,maxphase,maxp,massoftype) &
!$omp shared(use_dustfrac,dustfrac,dustevol,this_is_a_test,ndustsmall,alphaind,dvdx) &
!$omp shared(unit_density,unit_opacity,iopacity_type) &
!$omp private(i,spsound,pondens,rhoi,p_on_rhogas,rhogas,gasfrac) &
!$omp private(Bxi,Byi,Bzi,psii,xi_limiteri,Bi,temperaturei,ierr,pmassi) &
!$omp private(xi,yi,zi,hi) &
!$omp firstprivate(iactivei,iamtypei,iamgasi,iamdusti)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       !
       !--get pressure (actually pr/dens) and sound speed from equation of state
       !
       xi      = xyzh(1,i)
       yi      = xyzh(2,i)
       zi      = xyzh(3,i)
       hi      = xyzh(4,i)

       if (maxphase==maxp) call get_partinfo(iphase(i),iactivei,iamgasi,iamdusti,iamtypei)

       pmassi  = massoftype(iamtypei)
       rhoi    = rhoh(hi,pmassi)
       !
       !--Convert dust variable to dustfrac
       !
       if (use_dustfrac) then
          !--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
          if (.not.(use_dustgrowth .and. this_is_a_test)) &
             dustfrac(1:ndustsmall,i) = dustevol(:,i)**2/(1.+dustevol(:,i)**2)
          gasfrac = (1. - sum(dustfrac(1:ndustsmall,i)))  ! rhogas/rho
          rhogas  = rhoi*gasfrac       ! rhogas = (1-eps)*rho
       else
          rhogas  = rhoi
       endif
       if (.not. iamgasi) cycle  !stop here if not a gas particle

       !
       !--Calling Equation of state
       !
       if (maxvxyzu >= 4) then
          if (store_gamma) then
             call equationofstate(ieos,p_on_rhogas,spsound,rhogas,xi,yi,zi,eni=vxyzu(4,i),gamma_local=gamma_chem(i),&
                                  tempi=temperaturei)
          else
             call equationofstate(ieos,p_on_rhogas,spsound,rhogas,xi,yi,zi,eni=vxyzu(4,i),tempi=temperaturei)
          endif
       else
          !isothermal
          call equationofstate(ieos,p_on_rhogas,spsound,rhogas,xi,yi,zi,tempi=temperaturei)
       endif

       eos_vars(igasP,i)  = p_on_rhogas*rhogas
       eos_vars(ics,i)    = spsound
       eos_vars(itemp,i)  = temperaturei

       if (do_radiation) then
          !
          ! Get the opacity from the density and temperature if required
          !
          if (iopacity_type > 0) call get_opacity(iopacity_type,rhogas,temperaturei,radprop(ikappa,i))
          !
          ! Get radiation pressure from the radiation energy, i.e. P = 1/3 E if optically thick
          !
          call radiation_equation_of_state(radprop(iradP,i),rad(iradxi,i),rhogas)
       endif
       !
       ! Cullen & Dehnen (2010) viscosity switch, set alphaloc
       !
       if (nalpha >= 2) then
          xi_limiteri = xi_limiter(dvdx(:,i))
          alphaind(2,i) = real4(get_alphaloc(real(alphaind(3,i)),spsound,hi,xi_limiteri,alpha,alphamax))
       endif

       if (mhd) then
          ! construct B from B/rho (conservative to primitive)
          Bxi = Bevol(1,i) * rhoi
          Byi = Bevol(2,i) * rhoi
          Bzi = Bevol(3,i) * rhoi
          psii = Bevol(4,i)

          ! store primitive variables
          Bxyz(1,i) = Bxi
          Bxyz(2,i) = Byi
          Bxyz(3,i) = Bzi
          !
          !--calculate species number densities & non-ideal MHD coefficients
          !
          if (mhd_nonideal .and. iactivei) then
             Bi = sqrt(Bxi*Bxi + Byi*Byi + Bzi*Bzi)
             call nicil_get_ion_n(rhoi,temperaturei,n_R(:,i),n_electronT(i),ierr)
             if (ierr/=0) then
                call nicil_translate_error(ierr)
                if (ierr > 0) call fatal('densityiterate','error in Nicil in calculating number densities')
             endif
             call nicil_get_eta(eta_nimhd(iohm,i),eta_nimhd(ihall,i),eta_nimhd(iambi,i),Bi, &
                            rhoi,temperaturei,n_R(:,i),n_electronT(i),ierr)
             if (ierr/=0) then ! ierr is reset in the above subroutine
                call nicil_translate_error(ierr)
                if (ierr > 0) call fatal('densityiterate','error in Nicil in calculating eta')
             endif
          endif
       endif
    endif
 enddo
!$omp end parallel do

end subroutine cons2prim_everything

end module cons2prim
