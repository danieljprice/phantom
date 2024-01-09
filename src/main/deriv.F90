!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module deriv
!
! this module is a wrapper for the main derivative evaluation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: cons2prim, densityforce, derivutils, dim, externalforces,
!   forces, forcing, growth, io, linklist, metric_tools, options, part,
!   photoevap, ptmass, ptmass_radiation, radiation_implicit, timestep,
!   timestep_ind, timing
!
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 public :: derivs, get_derivs_global
 real, private :: stressmax

 private

contains

!-------------------------------------------------------------
!+
!  calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!+
!-------------------------------------------------------------
subroutine derivs(icall,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                  Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
                  dustevol,ddustevol,dustfrac,eos_vars,time,dt,dtnew,pxyzu,dens,metrics)
 use dim,            only:maxvxyzu,mhd,fast_divcurlB,gr,periodic,do_radiation,&
                          sink_radiation,use_dustgrowth
 use io,             only:iprint,fatal,error
 use linklist,       only:set_linklist
 use densityforce,   only:densityiterate
 use ptmass,         only:ipart_rhomax,ptmass_calc_enclosed_mass,ptmass_boundary_crossing
 use externalforces, only:externalforce
 use part,           only:dustgasprop,dvdx,Bxyz,set_boundaries_to_active,&
                          nptmass,xyzmh_ptmass,sinks_have_heating,dust_temp,VrelVf,fxyz_drag
#ifdef IND_TIMESTEPS
 use timestep_ind,   only:nbinmax
#else
 use timestep,       only:dtcourant,dtforce,dtrad
#endif
 use timestep,       only:dtmax
#ifdef DRIVING
 use forcing,        only:forceit
#endif
#ifdef PHOTO
 use photoevap,      only:find_ionfront,photo_ionize
 use part,           only:massoftype
#endif
 use growth,         only:get_growth_rate
 use ptmass_radiation, only:get_dust_temperature
 use timing,         only:get_timings
 use forces,         only:force
 use part,           only:mhd,gradh,alphaind,igas,iradxi,ifluxx,ifluxy,ifluxz,ithick
 use derivutils,     only:do_timing
 use cons2prim,      only:cons2primall,cons2prim_everything,prim2consall
 use metric_tools,   only:init_metric
 use radiation_implicit, only:do_radiation_implicit,ierr_failed_to_converge
 use options,        only:implicit_radiation,implicit_radiation_store_drad

 integer,      intent(in)    :: icall
 integer,      intent(inout) :: npart
 integer,      intent(in)    :: nactive
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(inout) :: vxyzu(:,:)
 real,         intent(inout) :: fxyzu(:,:)
 real,         intent(in)    :: fext(:,:)
 real(kind=4), intent(out)   :: divcurlv(:,:)
 real(kind=4), intent(out)   :: divcurlB(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real,         intent(out)   :: dBevol(:,:)
 real,         intent(inout) :: rad(:,:)
 real,         intent(out)   :: eos_vars(:,:)
 real,         intent(out)   :: drad(:,:)
 real,         intent(inout) :: radprop(:,:)
 real,         intent(in)    :: dustevol(:,:)
 real,         intent(inout) :: dustprop(:,:)
 real,         intent(out)   :: dustfrac(:,:)
 real,         intent(out)   :: ddustevol(:,:),ddustprop(:,:)
 real,         intent(in)    :: time,dt
 real,         intent(out)   :: dtnew
 real,         intent(inout) :: pxyzu(:,:), dens(:)
 real,         intent(inout) :: metrics(:,:,:,:)
 integer                     :: ierr,i
 real(kind=4)                :: t1,tcpu1,tlast,tcpulast

 t1    = 0.
 tcpu1 = 0.
 call get_timings(t1,tcpu1)
 tlast    = t1
 tcpulast = tcpu1
!
!--check for errors in input options
!
 if (icall < 0 .or. icall > 2) call fatal('deriv','invalid icall on input')
!
! icall is a flag to say whether or not positions have changed
! since the last call to derivs.
!
! icall = 1 is the "standard" call to derivs: calculates all derivatives
! icall = 2 does not remake the link list and does not recalculate density
!           (ie. only re-evaluates the SPH force term using updated values
!            of the input variables)
!
! call link list to find neighbours
!
 if (icall==1 .or. icall==0) then
    call set_linklist(npart,nactive,xyzh,vxyzu)

    if (gr) then
       ! Recalculate the metric after moving particles to their new tasks
       call init_metric(npart,xyzh,metrics)
       call prim2consall(npart,xyzh,metrics,vxyzu,dens,pxyzu,use_dens=.false.)
    endif

    if (nptmass > 0 .and. periodic) call ptmass_boundary_crossing(nptmass,xyzmh_ptmass)
 endif

 call do_timing('link',tlast,tcpulast,start=.true.)

#ifdef PHOTO
 !
 ! update location of particles on grid and calculate the location of the ionization front
 !
 call find_ionfront(time,npart,xyzh,massoftype(igas))
 !
 ! update the temperatures of the particles depending on whether ionized or not
 !
 call photo_ionize(vxyzu,npart)
#endif
!
! calculate density by direct summation
!
 if (icall==1) then
    call densityiterate(1,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                        stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx)
    if (.not. fast_divcurlB) then
       ! Repeat the call to calculate all the non-density-related quantities in densityiterate.
       ! This needs to be separate for an accurate calculation of divcurlB which requires an up-to-date rho.
       ! if fast_divcurlB = .false., then all additional quantities are calculated during the previous call
       call densityiterate(3,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                           stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx)
    endif
    set_boundaries_to_active = .false.     ! boundary particles are no longer treated as active
    call do_timing('dens',tlast,tcpulast)
 endif
 
 if (gr) then
    call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
 else
    call cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,Bevol,Bxyz,dustevol,dustfrac,alphaind)
 endif
 call do_timing('cons2prim',tlast,tcpulast)

 !
 ! implicit radiation update
 !
 if (do_radiation .and. implicit_radiation .and. dt > 0.) then
    call do_radiation_implicit(dt,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    if (ierr /= 0 .and. ierr /= ierr_failed_to_converge) call fatal('radiation','Failed in radiation')
 endif
 
!
! compute forces
!
#ifdef DRIVING
 ! forced turbulence -- call driving routine
 call forceit(time,npart,xyzh,vxyzu,fxyzu)
 call do_timing('driving',tlast,tcpulast)
#endif
 stressmax = 0.
 if (sinks_have_heating(nptmass,xyzmh_ptmass)) call ptmass_calc_enclosed_mass(nptmass,npart,xyzh)
 call force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&
            rad,drad,radprop,dustprop,dustgasprop,dustfrac,ddustevol,fext,fxyz_drag,&
            ipart_rhomax,dt,stressmax,eos_vars,dens,metrics)
 call do_timing('force',tlast,tcpulast)

 if (use_dustgrowth) then ! compute growth rate of dust particles
    call get_growth_rate(npart,xyzh,vxyzu,dustgasprop,VrelVf,dustprop,ddustprop(1,:))!--we only get ds/dt (i.e 1st dimension of ddustprop)
 endif

!
! compute dust temperature
!
 if (sink_radiation .and. maxvxyzu == 4) then
    call get_dust_temperature(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp)
 endif

 if (do_radiation .and. implicit_radiation .and. .not.implicit_radiation_store_drad) then
    !$omp parallel do shared(drad,fxyzu,npart) private(i)
    do i=1,npart
       drad(:,i) = 0.
       fxyzu(4,i) = 0.
    enddo
    !$omp end parallel do
 endif
!
! set new timestep from Courant/forces condition
!
#ifdef IND_TIMESTEPS
 dtnew = dtmax/2**nbinmax  ! minimum timestep over all particles
#else
 dtnew = min(dtforce,dtcourant,dtrad,dtmax)
#endif

 call do_timing('total',t1,tcpu1,lunit=iprint)

end subroutine derivs

!--------------------------------------
!+
!  wrapper for the call to derivs
!  so only one line needs changing
!  if interface changes
!
!  this should NOT be called during timestepping, it is useful
!  for when one requires just a single call to evaluate derivatives
!  and store them in the global shared arrays
!+
!--------------------------------------
subroutine get_derivs_global(tused,dt_new,dt)
 use part,   only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
                dustfrac,ddustevol,eos_vars,pxyzu,dens,metrics,dustevol
 use timing, only:printused,getused
 use io,     only:id,master
 real(kind=4), intent(out), optional :: tused
 real,         intent(out), optional :: dt_new
 real,         intent(in), optional  :: dt  ! optional argument needed to test implicit radiation routine
 real(kind=4) :: t1,t2
 real :: dtnew
 real :: time,dti

 time = 0.
 dti = 0.
 if (present(dt)) dti = dt
 call getused(t1)
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
             rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,dustfrac,eos_vars,&
             time,dti,dtnew,pxyzu,dens,metrics)
 call getused(t2)
 if (id==master .and. present(tused)) call printused(t1)
 if (present(tused)) tused = t2 - t1
 if (present(dt_new)) dt_new = dtnew

end subroutine get_derivs_global

subroutine calc_lambda(npart,xyzh,vxyzu,rho,lambda)
  use io,    only:error
  use dim,   only:maxneigh
  use part,  only:get_partinfo,iamgas,iboundary,igas,maxphase,massoftype,iphase,gradh
#ifdef PERIODIC
  use boundary,  only:dxbound,dybound,dzbound
#endif
  use kernel,   only:radkern2,wkern,grkern,cnormk,get_kernel
  use linklist, only:get_neighbour_list
  use units,    only:unit_density,unit_ergg,unit_opacity
  use eos_stamatellos, only:getopac_opdep

  integer,         intent(in)     :: npart 
  real,            intent(in)     :: xyzh(:,:)
  real,            intent(in)     :: vxyzu(:,:)
  real,            intent(in)     :: rho(:)
  real,            intent(inout)  :: lambda(:)

  integer, parameter  :: maxcellcache = 10000
  logical :: iactivei,iamdusti,iamgasi,iactivej,iamgasj,iamdustj
  integer :: iamtypei,i,j,n,iamtypej,mylistneigh(maxneigh),nneigh
  integer(kind=1) :: iphasei,iphasej
  real                :: rhoi,rhoj,uradi,xi,yi,zi
  real                :: dradi,Ti,Tj,kappaBarj,kappaPartj,gmwj,gammaj,wkerni,grkerni
  real                :: kappaBari,kappaParti,gmwi,dx,dy,dz,rij2,rij,rij1,Wi,dWi
  real                :: dradxi,dradyi,dradzi,runix,runiy,runiz,R_rad,dT4
  real                :: pmassi,pmassj,hi,hi21,hi1,q,q2,q2i,xj,yj,zj,qi,q2j
  real                :: hj21,hi4i,hj,hj1,hi41,hi31,gradhi
  real                :: xyzcache(maxcellcache,4)
  logical             :: added_self,ignoreself
  integer             :: success_count
  ignoreself = .true.
  
  print *, "calculating lambda in deriv"
  success_count = 0
! loop over parts                                                                     
!$omp parallel do schedule(runtime) default(none) &
!$omp shared(npart,xyzh,iphase,ignoreself,unit_opacity) &
!$omp shared(rho,lambda,massoftype,vxyzu,unit_density,unit_ergg,gradh) &
!$omp private(i,j,n,mylistneigh,nneigh,xyzcache,iphasei,iactivei) &
!$omp private(iamdusti,iamgasi,iactivej,iamgasj,iamdustj,iamtypei,iamtypej) &
!$omp private(rhoi,rhoj,uradi,xi,yi,zi,dradi,Ti,Tj,kappaBarj,kappaPartj,gmwj) &
!$omp private(gammaj,wkerni,grkerni,kappaBari,kappaParti,gmwi,dx,dy,dz,rij2) &
!$omp private(rij,rij1,Wi,dWi,dradxi,dradyi,dradzi,runix,runiy,runiz,R_rad,dT4) &
!$omp private(pmassi,pmassj,hi,hi21,hi1,q,q2,q2i,xj,yj,zj,qi,q2j,hj21,hi4i,hj,hj1) &
!$omp private(hi41,hi31,gradhi,added_self,success_count,iphasej)

  over_parts: do i = 1,npart
    !                                                                              
    !---loop of neighbours to calculate radiation energy density
    !                                                    
     !check active and gas particle                         
     call get_neighbour_list(i,mylistneigh,nneigh,xyzh,xyzcache,maxcellcache, &
                           getj=.true.) 
     iphasei = iphase(i)
     call get_partinfo(iphasei,iactivei,iamgasi,iamdusti,iamtypei)
     print *, i, iactivei,iamgasi,rho(i),"nneigh:",nneigh
     if (.not. iactivei) cycle over_parts
     if (iamtypei == iboundary) cycle over_parts
     if (.not. iamgasi) cycle over_parts
    
     uradi = 0.
     dradi = 0.
     pmassi = massoftype(iamtypei)
     rhoi = rho(i)
     added_self = .false.

     call getopac_opdep(vxyzu(4,i)*unit_ergg,rhoi*unit_density,kappabari, &
          kappaparti,ti,gmwi)
     hi    = xyzh(4,i)
     hi1   = 1./hi
     hi21  = hi1*hi1
     hi31  = hi1*hi21
     hi41  = hi21*hi21
     xi = xyzh(1,i)
     yi = xyzh(2,i)
     zi = xyzh(3,i)
     !                                                                  
     !--compute density and related quantities from the smoothing length
     !                                                                  
     pmassi = massoftype(iamtypei)
     dradxi = 0.0
     dradyi = 0.0
     dradzi = 0.0

 loop_over_neighbours: do n = 1,nneigh
        j = mylistneigh(n)
        if (j < 0) cycle loop_over_neighbours
        iphasej = iphase(j) 
        call get_partinfo(iphasej,iactivej,iamgasj,iamdustj,iamtypej)
        if (.not. iactivej) cycle loop_over_neighbours
        if (iamtypej == iboundary) cycle loop_over_neighbours
        if (.not. iamgasj) cycle loop_over_neighbours
        if ((ignoreself) .and. (i==j)) cycle loop_over_neighbours
        if (xyzh(4,j) < 0d0) cycle loop_over_neighbours
        print *, "n=", n
        if (n <= maxcellcache) then
           ! positions from cache are already mod boundary                              
           xj = xyzcache(n,1)
           yj = xyzcache(n,2)
           zj = xyzcache(n,3)
        else
           xj = xyzh(1,j)
           yj = xyzh(2,j)
           zj = xyzh(3,j)
        endif

        dx = xi - xj
        dy = yi - yj
        dz = zi - zj
#ifdef PERIODIC
        if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
        if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
        if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
        rij2 = dx*dx + dy*dy + dz*dz + TINY(0.)
        rij = SQRT(rij2)
        q2i = rij2*hi21
        qi = SQRT(q2i)

       !--hj is in the cell cache but not in the neighbour cache
       !  as not accessed during the density summation          
        hj1 = 1./xyzh(4,j)
        hj = 1./hj1
        hj21 = hj1*hj1
        q2j  = rij2*hj21
        is_sph_neighbour: if (q2i < radkern2 .or. q2j < radkern2) then
           call get_partinfo(iphase(j),iactivej,iamgasj,iamdustj,iamtypej)
           if (.not.iamgasj) cycle loop_over_neighbours
           if (.not.iactivej) cycle loop_over_neighbours
           !Get kernel quantities                  
           if (gradh(1,i) > 0.) then
              gradhi = gradh(1,i)
           else
  !call error('force','stored gradh is zero, resetting to 1')                
              gradhi = 1.
           endif
           !print *, "q2i,q2j", q2i,q2j
           call get_kernel(q2,q,wkerni,grkerni)
           Wi = wkerni*cnormk*hi21*hi1
           dWi = grkerni*cnormk*hi21*hi21*gradh(1,i)
           pmassj = massoftype(iamtypej)
           rhoj = rho(j)
           if (rhoj < tiny(rhoj)) cycle loop_over_neighbours
          call getopac_opdep(vxyzu(4,j)*unit_ergg,rhoj*unit_density,kappaBarj, &
               kappaPartj,Tj,gmwj)
!          uradi = uradi + arad*pmassj*Tj**4.0d0*Wi/(rhoj)
          if (uradi > 1e40) then
             print *, "Urad huge- i,j,pmassj,Tj,Wi,hi,rhoj:", i,j,pmassj,Tj,Wi,1/hi1,rhoj,xj,yj,zj 
             print *, "wkerni,cnormk,hi21,hi1", wkerni,cnormk,hi21,hi1,tiny(rhoj)
             print *, "max min rho", maxval(rho), minval(rho)
             stop
          endif
          ! unit vector components    
          runix = dx/rij
          runiy = dy/rij
          runiz = dz/rij

          dT4 = Ti**4d0 - Tj**4d0
!          dradxi = dradxi + arad*pmassj*dT4*dWi*runix/rhoj
!          dradyi = dradyi + arad*pmassj*dT4*dWi*runiy/rhoj
!          dradzi = dradzi + arad*pmassj*dT4*dWi*runiz/rhoj
       endif is_sph_neighbour

    enddo loop_over_neighbours

    if (.not. added_self) then
       !       print *, "Has not added self in lambda hybrid
!       uradi = uradi + cnormk*hi1*hi21*pmassj*arad*Ti**4d0/rhoi ! add self contribution
    endif
    !    urad  = cnormk*uradi + arad*Ti**4.0d0
    dradi = cnormk*SQRT(dradxi**2.0d0 + dradyi**2.0d0 + dradzi**2.0d0)
    !Now calculate flux limiter coefficients                          
    !Calculate in cgs (converted to code units in forcei)

    if ((dradi.eq.0.0d0).or.(uradi.eq.0.0d0)) then
       R_rad = 0.0d0
    else
       R_rad = dradi/(uradi*rhoi*kappaParti/unit_opacity)
    endif

    lambda(i) = (2.0d0+R_rad)/(6.0d0+3.0d0*R_rad+R_rad*R_rad)
    if (isnan(lambda(i))) then
       print *, "lambda isnan when calculated. i, R_Rad, uradi,dradi,rhoi,kappaParti, Ti", &
            i,R_Rad,uradi,dradi,rhoi,kappaParti,Ti
    elseif (lambda(i) < tiny(lambda(i))) then
       print *, "lambda is 0,i, rhoi, uradi", i, rhoi, uradi
       stop
    else
!$omp critical
       success_count = success_count + 1
!$omp end critical
    endif
 enddo over_parts
!$omp end parallel do
 print *, "success =", success_count
 do i=1,npart
    if (lambda(i) < tiny(lambda(i))) then
       print *, "lambda is 0,i, rhoi, uradi", i, rhoi, uradi,kappaParti
    endif
 enddo
end subroutine calc_lambda

end module deriv
