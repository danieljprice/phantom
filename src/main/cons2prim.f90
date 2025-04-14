!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cons2prim
!
! Subroutines to swap between primitive variables needed on RHS of fluid
! equations (density,velocity,internal energy) and conserved/evolved
! variables on the LHS of the fluid equations (rho*,momentum,entropy)
!
! This is complicated in the GR code but also useful to structure
! things this way in the non-GR code, e.g. B/rho is the evolved variable
! while B is the "primitive" variable for magnetic field. Similarly
! sqrt(rho_d) is the evolved variable for dust species but the primitive
! variable is the dust mass fraction.
!
! :References:
!   Liptai & Price (2019), MNRAS 485, 819-842
!   Ballabio et al. (2018), MNRAS 477, 2766-2771
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: cons2primsolver, cullendehnen, dim, eos, io, nicil,
!   options, part, radiation_utils, utils_gr
!
 implicit none

 public :: cons2primall,cons2prim_everything,cons2primall_sink
 public :: prim2consall,prim2consi

 private

contains

!----------------------------------------------------------------------
!+
!  Primitive to conservative transform (for GR):
!  Construct conserved variables (rho*,momentum,entropy)
!  from the primitive/fluid rest frame variables
!  (density,velocity,internal energy), for ALL particles
!+
!----------------------------------------------------------------------
subroutine prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,use_dens,dens,use_sink)
 use part, only:isdead_or_accreted,ien_type,eos_vars,igasP,igamma,itemp
 use eos,  only:gamma,ieos
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:),metrics(:,:,:,:),vxyzu(:,:)
 real,    intent(out) :: pxyzu(:,:)
 real,    intent(inout), optional :: dens(:)
 logical, intent(in),    optional :: use_dens, use_sink
 logical :: usedens
 integer :: i
 real    :: pri,tempi,xyzhi(4),vxyzui(4),densi

!  By default, use the smoothing length to compute primitive density, and then compute the conserved variables.
!  (Alternatively, use the provided primitive density to compute conserved variables.
!   Depends whether you have prim dens prior or not.)
 if (present(use_dens)) then
    usedens = use_dens
 else
    usedens = .false.
 endif

 !$omp parallel do default (none) &
 !$omp shared(xyzh,metrics,vxyzu,dens,pxyzu,npart,usedens,ien_type,eos_vars,gamma,ieos,use_sink,use_dens) &
 !$omp private(i,pri,tempi,xyzhi,vxyzui,densi)
 do i=1,npart

    if (present(use_sink)) then
       xyzhi(1:3) = xyzh(1:3,i) ! save positions
       xyzhi(4) = xyzh(5,i) ! save smoothing length, h
       vxyzui(1:3) = vxyzu(1:3,i)
       vxyzui(4) = 0. ! assume energy as 0. for sink
       densi = 1.
       call prim2consi(xyzhi,metrics(:,:,:,i),vxyzui,pri,tempi,pxyzu(:,i),ien_type,&
                   use_sink=use_sink,dens_i=densi) ! this returns temperature and pressure as 0.
    else
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          call prim2consi(xyzh(:,i),metrics(:,:,:,i),vxyzu(:,i),pri,tempi,pxyzu(:,i),ien_type,&
                   use_dens=usedens,dens_i=dens(i))

          ! save eos vars for later use
          eos_vars(igasP,i)  = pri
          eos_vars(itemp,i)  = tempi
          if (vxyzu(4,i) > 0. .and. ieos == 12) then
             eos_vars(igamma,i) = 1. + pri/(dens(i)*vxyzu(4,i))
          else
             ! prevent getting NaN or Infinity when u = 0
             eos_vars(igamma,i) = gamma
          endif
       endif
    endif
 enddo
 !$omp end parallel do

end subroutine prim2consall

!----------------------------------------------------------------------
!+
!  Primitive to conservative transform (for GR):
!  for a single SPH particle
!+
!----------------------------------------------------------------------
subroutine prim2consi(xyzhi,metrici,vxyzui,pri,tempi,pxyzui,ien_type,use_dens,use_sink,dens_i)
 use cons2primsolver, only:primitive2conservative
 use utils_gr,        only:h2dens
 use eos,             only:equationofstate,ieos
 real, dimension(4), intent(in)  :: xyzhi, vxyzui
 real,               intent(in)  :: metrici(:,:,:)
 real, intent(inout)             :: pri,tempi
 integer,            intent(in)  :: ien_type
 real, dimension(4), intent(out) :: pxyzui
 logical, intent(in), optional   :: use_dens,use_sink
 real, intent(inout), optional   :: dens_i
 logical :: usedens
 real    :: rhoi,ui,xyzi(1:3),vi(1:3),pondensi,spsoundi,densi

 pondensi = 0.
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
    if (present(use_sink)) then
       densi    = 1.    ! using a value of 0. results in NaN values for the pxyzui array.
       pondensi = 0.
    else
       call h2dens(densi,xyzhi,metrici,vi) ! Compute dens from h
       dens_i = densi ! Feed the newly computed dens back out of the routine
       call equationofstate(ieos,pondensi,spsoundi,densi,xyzi(1),xyzi(2),xyzi(3),tempi,ui)
    endif
 endif

 pri = pondensi*densi
 call primitive2conservative(xyzi,metrici,vi,densi,ui,pri,rhoi,pxyzui(1:3),pxyzui(4),ien_type)

end subroutine prim2consi

!----------------------------------------------------------------------
!+
!  Conservative to primitive routines (for GR):
!  Solve for primitive variables (density,velocity,internal energy)
!  from the evolved/conservative variables (rho*,momentum,entropy)
!+
!----------------------------------------------------------------------
subroutine cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
 use cons2primsolver, only:conservative2primitive
 use part,            only:isdead_or_accreted,massoftype,igas,rhoh,igasP,ics,ien_type,&
                           itemp,igamma,aprmassoftype,apr_level
 use io,              only:fatal
 use eos,             only:ieos,done_init_eos,init_eos,get_spsound
 use dim,             only:use_apr
 use eos,             only:ieos,done_init_eos,init_eos,get_spsound
 integer, intent(in)    :: npart
 real,    intent(in)    :: pxyzu(:,:),xyzh(:,:),metrics(:,:,:,:)
 real,    intent(inout) :: vxyzu(:,:),dens(:)
 real,    intent(out)   :: eos_vars(:,:)
 integer :: i, ierr
 real    :: p_guess,rhoi,tempi,gammai,pmassi

 if (.not.done_init_eos) call init_eos(ieos,ierr)

!$omp parallel do default (none) &
!$omp shared(xyzh,metrics,vxyzu,dens,pxyzu,npart,massoftype,aprmassoftype) &
!$omp shared(ieos,eos_vars,ien_type,apr_level) &
!$omp private(i,ierr,p_guess,rhoi,tempi,gammai,pmassi)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       ! get pressure, temperature and gamma from previous step as the initial guess
       p_guess = eos_vars(igasP,i)
       tempi   = eos_vars(itemp,i)
       gammai  = eos_vars(igamma,i)
       if (use_apr) then
          pmassi = aprmassoftype(igas,apr_level(i))
       else
          pmassi = massoftype(igas)
       endif
       rhoi    = rhoh(xyzh(4,i),pmassi)

       call conservative2primitive(xyzh(1:3,i),metrics(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i), &
                                  p_guess,tempi,gammai,rhoi,pxyzu(1:3,i),pxyzu(4,i),ierr,ien_type)
       ! store results
       eos_vars(igasP,i)  = p_guess
       eos_vars(ics,i)    = get_spsound(ieos,xyzh(1:3,i),dens(i),vxyzu(:,i),gammai)
       eos_vars(itemp,i)  = tempi
       eos_vars(igamma,i) = gammai
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

!----------------------------------------------------------------------
!+
!  Conservative to primitive routines (for GR sink particles):
!  Solve for primitive variables (density,velocity,internal energy)
!  from the evolved/conservative variables (rho*,momentum,entropy)
!+
!----------------------------------------------------------------------
subroutine cons2primall_sink(nptmass,xyzmh_ptmass,metrics_ptmass,pxyzu_ptmass,vxyz_ptmass,eos_vars)
 use cons2primsolver, only:conservative2primitive
 use io,              only:fatal
 use part,            only:ien_type
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: pxyzu_ptmass(:,:),xyzmh_ptmass(:,:),metrics_ptmass(:,:,:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(out), optional   :: eos_vars(:,:)
 integer :: i, ierr
 real    :: p_guess,rhoi,tempi,gammai,eni,densi

!$omp parallel do default (none) &
!$omp shared(xyzmh_ptmass,metrics_ptmass,vxyz_ptmass,pxyzu_ptmass,nptmass,ien_type) &
!$omp private(i,ierr,p_guess,rhoi,tempi,gammai,eni,densi)
 do i=1,nptmass
    p_guess = 0.
    tempi   = 0.
    gammai  = 0.
    rhoi    = 1.
    densi   = 1.
    ! conservative 2 primitive
    call conservative2primitive(xyzmh_ptmass(1:3,i),metrics_ptmass(:,:,:,i),vxyz_ptmass(1:3,i),densi,eni, &
                              p_guess,tempi,gammai,rhoi,pxyzu_ptmass(1:3,i),pxyzu_ptmass(4,i),ierr,ien_type)

    if (ierr > 0) then
       print*,' pmom =',pxyzu_ptmass(1:3,i)
       print*,' rho* =',rhoi
       print*,' en   =',eni
       call fatal('cons2prim','could not solve rootfinding',i)
    endif

 enddo
!$omp end parallel do

end subroutine cons2primall_sink

!-----------------------------------------------------------------------------
!+
!  Solve for primitive variables (v,u,P,B,dustfrac) from evolved variables
!  (v,energy variable,B/rho,sqrt(rho*eps)) in the non-relativistic code
!
!  In this case no "solver" is required, but we do need to call the
!  equation of state to get the pressure
!+
!-----------------------------------------------------------------------------
subroutine cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,&
                                Bevol,Bxyz,dustevol,dustfrac,alphaind)
 use part,              only:isdead_or_accreted,massoftype,igas,rhoh,igasP,iradP,iradxi,ics,imu,iX,iZ,&
                             iohm,ihall,nden_nimhd,eta_nimhd,iambi,get_partinfo,iphase,this_is_a_test,&
                             ndustsmall,itemp,ikappa,idmu,idgamma,icv,aprmassoftype,apr_level,isionised
 use part,              only:nucleation,igamma
 use eos,               only:equationofstate,ieos,eos_outputs_mu,done_init_eos,init_eos,gmw,X_in,Z_in,gamma
 use radiation_utils,   only:radiation_equation_of_state,get_opacity
 use dim,               only:mhd,maxvxyzu,maxphase,maxp,use_dustgrowth,&
                             do_radiation,nalpha,mhd_nonideal,do_nucleation,use_krome,update_muGamma,use_apr
 use nicil,             only:nicil_update_nimhd,nicil_translate_error,n_warn
 use io,                only:fatal,real4,warning
 use cullendehnen,      only:get_alphaloc,xi_limiter
 use options,           only:alpha,alphamax,use_dustfrac,iopacity_type,use_var_comp,implicit_radiation
 integer,      intent(in)    :: npart
 real,         intent(in)    :: xyzh(:,:),rad(:,:),Bevol(:,:),dustevol(:,:)
 real(kind=4), intent(in)    :: dvdx(:,:)
 real,         intent(inout) :: vxyzu(:,:)
 real(kind=4), intent(inout) :: alphaind(:,:)
 real,         intent(out)   :: eos_vars(:,:),radprop(:,:),Bxyz(:,:),dustfrac(:,:)
 integer      :: i,iamtypei,ierr
 integer      :: ierrlist(n_warn)
 real         :: rhoi,spsound,p_on_rhogas,rhogas,gasfrac,pmassi,uui
 real         :: Bxi,Byi,Bzi,psii,xi_limiteri,Bi,temperaturei,mui,X_i,Z_i,gammai
 real         :: xi,yi,zi,hi
 logical      :: iactivei,iamgasi,iamdusti

 iactivei = .true.
 iamtypei = igas
 iamgasi  = .true.
 iamdusti = .false.
 ierrlist = 0
 if (.not.done_init_eos) then
    call init_eos(ieos,ierr)
    if (ierr /= 0) call fatal('eos','could not initialise equation of state')
 endif
 gammai = gamma
 mui    = gmw
 X_i    = X_in
 Z_i    = Z_in

!$omp parallel do default (none) &
!$omp shared(xyzh,vxyzu,npart,rad,eos_vars,radprop,Bevol,Bxyz,apr_level) &
!$omp shared(ieos,nucleation,nden_nimhd,eta_nimhd) &
!$omp shared(alpha,alphamax,iphase,maxphase,maxp,massoftype,aprmassoftype,isionised) &
!$omp shared(use_dustfrac,dustfrac,dustevol,this_is_a_test,ndustsmall,alphaind,dvdx) &
!$omp shared(iopacity_type,use_var_comp,do_nucleation,update_muGamma,implicit_radiation) &
!$omp private(i,spsound,rhoi,p_on_rhogas,rhogas,gasfrac,uui) &
!$omp private(Bxi,Byi,Bzi,psii,xi_limiteri,Bi,temperaturei,ierr,pmassi) &
!$omp private(xi,yi,zi,hi) &
!$omp firstprivate(iactivei,iamtypei,iamgasi,iamdusti,mui,gammai,X_i,Z_i) &
!$omp reduction(+:ierrlist)
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

       if (use_apr) then
          pmassi = aprmassoftype(iamtypei,apr_level(i))
       else
          pmassi  = massoftype(iamtypei)
       endif
       rhoi    = rhoh(hi,pmassi)
       !
       !--Convert dust variable to dustfrac
       !
       if (use_dustfrac) then
          !--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
          if (.not.(use_dustgrowth .and. this_is_a_test)) &
             dustfrac(1:ndustsmall,i) = dustevol(1:ndustsmall,i)**2/(1.+dustevol(1:ndustsmall,i)**2)
          gasfrac = (1. - sum(dustfrac(1:ndustsmall,i)))  ! rhogas/rho
          rhogas  = rhoi*gasfrac       ! rhogas = (1-eps)*rho
          if (gasfrac < 0.) call warning('cons2prim','total dust fraction > 1: try limiting dust flux',i,'eps',1.-gasfrac)
       else
          rhogas  = rhoi
       endif
       if (.not. iamgasi) cycle  !stop here if not a gas particle

       !
       !--Calling Equation of state
       !
       temperaturei = eos_vars(itemp,i) ! needed for initial guess for idealplusrad
       if (use_var_comp) then
          mui = eos_vars(imu,i)
          X_i = eos_vars(iX,i)
          Z_i = eos_vars(iZ,i)
       endif
       if (do_nucleation) then
          mui    = nucleation(idmu,i)
          gammai = nucleation(idgamma,i)
       endif
       if (update_muGamma) then
          mui    = eos_vars(imu,i)
          gammai = eos_vars(igamma,i)
       endif
       if (use_krome) gammai = eos_vars(igamma,i)
       if (maxvxyzu >= 4) then
          uui = vxyzu(4,i)
          if (uui < 0. .and. .not. ieos==23) then
             call warning('cons2prim','Internal energy < 0',i,'u',uui)
          endif
          call equationofstate(ieos,p_on_rhogas,spsound,rhogas,xi,yi,zi,temperaturei,eni=uui,&
                               gamma_local=gammai,mu_local=mui,Xlocal=X_i,Zlocal=Z_i,isionised=isionised(i))
       else
          !isothermal
          call equationofstate(ieos,p_on_rhogas,spsound,rhogas,xi,yi,zi,temperaturei,mu_local=mui, &
                               isionised=isionised(i))
       endif

       eos_vars(igasP,i)  = p_on_rhogas*rhogas
       eos_vars(ics,i)    = spsound
       eos_vars(itemp,i)  = temperaturei
       if (use_var_comp .or. eos_outputs_mu(ieos) .or. do_nucleation .or. update_muGamma) eos_vars(imu,i) = mui

       if (do_radiation) then
          if (temperaturei > tiny(0.)) then
             radprop(icv,i) = vxyzu(4,i)/temperaturei
          else
             radprop(icv,i) = 1. ! arbitrary, but should give zero for u when u=cv*T
          endif
          if (.not. implicit_radiation) then
             !
             ! Get the opacity from the density and temperature if required
             !
             if (iopacity_type > 0) call get_opacity(iopacity_type,rhogas,temperaturei,radprop(ikappa,i))
          endif
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
          if (mhd_nonideal) then
             Bi = sqrt(Bxi*Bxi + Byi*Byi + Bzi*Bzi)
             ! sanity check the temperature
             if (temperaturei < 1.) call warning('cons2prim',&
                'T < 1K in non-ideal MHD library',i,'T',temperaturei)

             call nicil_update_nimhd(0,eta_nimhd(iohm,i),eta_nimhd(ihall,i),eta_nimhd(iambi,i), &
                                     Bi,rhoi,temperaturei,nden_nimhd(:,i),ierrlist)
          endif
       endif
    endif
 enddo
!$omp end parallel do

 if (mhd_nonideal) then
    ! look for fatal errors in nicil and kill if necessary
    if ( any(ierrlist > 0) ) then
       call nicil_translate_error(ierrlist,.true.)
       call fatal('cons2prim_everything','error in Nicil')
    endif
 endif

end subroutine cons2prim_everything

end module cons2prim
