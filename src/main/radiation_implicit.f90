!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module radiation_implicit
!
! Implicit scheme for radiative transfer in flux limited diffusion
! approximation
!
! :References:
!   Whitehouse, Bate & Monaghan (2005), MNRAS 364, 1367-1377
!   Whitehouse & Bate (2004), MNRAS 353, 1078-1094
!   Bate & Keto (2015), MNRAS 449, 2643-2667
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, eos, io, kdtree, kernel, linklist, options,
!   part, physcon, quartic, radiation_utils, units
!
 use part,            only:ikappa,ilambda,iedd,idkappa,iradxi,icv,ifluxx,ifluxy,ifluxz,igas,rhoh,massoftype,imu
 use eos,             only:iopacity_type
 use radiation_utils, only:get_kappa,get_1overmu
 use eos,             only:get_cv
 implicit none
 integer, parameter :: ierr_failed_to_converge = 1,&
                       ierr_negative_opacity = 2, &
                       ierr_neighbourlist_empty = 3
 integer, parameter :: gas_dust_collisional_term_type = 0

 ! options for Bate & Keto ISM radiative transfer (not working yet)
 logical, parameter :: dustRT = .false.
 logical, parameter :: H2formation_heating = .false.
 logical, parameter :: use_cosmic_ray_heating = .false.
 logical, parameter :: use_photoelectric_heating = .false.
 real, parameter    :: Tdust_threshold = 100.

 ! options for the input file, with default values
 real, public       :: tol_rad = 1.e-4
 integer, public    :: itsmax_rad = 250
 integer, public    :: cv_type = 0

 character(len=*), parameter :: label = 'radiation_implicit'

 private
 public :: do_radiation_implicit,ierr_failed_to_converge

contains

!---------------------------------------------------------
!+
!  main subroutine
!+
!---------------------------------------------------------
subroutine do_radiation_implicit(dt,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
 use io, only:fatal,warning
 integer, intent(in)  :: npart
 real, intent(in)     :: dt,xyzh(:,:)
 real, intent(inout)  :: radprop(:,:),rad(:,:),vxyzu(:,:),drad(:,:)
 integer, intent(out) :: ierr
 integer              :: nsubsteps,i,nit
 logical              :: failed,moresweep
 real                 :: dtsub,errorE,errorU
 real, allocatable    :: origEU(:,:),EU0(:,:)

 ierr = 0

 allocate(origEU(2,npart),EU0(2,npart),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','could not allocate memory to origEU and EU0')

 call save_radiation_energies(npart,rad,xyzh,vxyzu,radprop,drad,origEU,.false.)

 nsubsteps = 1
 moresweep = .true.
 do while (moresweep)
    moresweep = .false.
    dtsub = dt/nsubsteps
    over_substeps: do i = 1,nsubsteps
       call do_radiation_onestep(dtsub,rad,xyzh,vxyzu,radprop,origEU,EU0,failed,nit,errorE,errorU,moresweep,ierr)
       if (failed .or. moresweep) then
          ierr = ierr_failed_to_converge
          call warning('radiation_implicit','integration failed - using U and E values anyway')
          moresweep = .false.
          !exit over_substeps
       endif
       if (i /= nsubsteps) call save_radiation_energies(npart,rad,xyzh,vxyzu,radprop,drad,origEU,.true.)
    enddo over_substeps

    !if (moresweep) then
    !    call restore_radiation_energies(npart,rad,xyzh,vxyzu,radprop,origEU,.true.)
    !    nsubsteps = nsubsteps*2
    !endif
 enddo

 deallocate(origEU,EU0)

end subroutine do_radiation_implicit


!---------------------------------------------------------
!+
!  save values of E, U
!+
!---------------------------------------------------------
subroutine save_radiation_energies(npart,rad,xyzh,vxyzu,radprop,drad,origEU,save_cv)
 integer, intent(in) :: npart
 real, intent(in)    :: rad(:,:),xyzh(:,:),vxyzu(:,:)
 logical, intent(in) :: save_cv
 real, intent(inout) :: radprop(:,:)
 real, intent(out)   :: origEU(:,:),drad(:,:)
 integer             :: i
 real                :: rhoi

 !$omp parallel do schedule(static) default(none) &
 !$omp shared(rad,origeu,vxyzu,xyzh,npart,radprop,save_cv,iopacity_type,massoftype,drad,cv_type) &
 !$omp private(i,rhoi)
 do i = 1,npart
    origEU(1,i) = rad(iradxi,i)
    origEU(2,i) = vxyzu(4,i)
    if (save_cv) then
       rhoi = rhoh(xyzh(4,i),massoftype(igas))
       radprop(icv,i) = get_cv(rhoi,vxyzu(4,i),cv_type)
       radprop(ikappa,i) = get_kappa(iopacity_type,vxyzu(4,i),radprop(icv,i),rhoi)
    endif
    drad(:,i) = 0.  ! Set dxi/dt = 0 for implicit scheme
 enddo
 !$omp end parallel do

end subroutine save_radiation_energies

!---------------------------------------------------------
!+
!  perform single iteration
!+
!---------------------------------------------------------
subroutine do_radiation_onestep(dt,rad,xyzh,vxyzu,radprop,origEU,EU0,failed,nit,errorE,errorU,moresweep,ierr)
 use io,      only:fatal,error,iverbose,warning
 use part,    only:hfact
 use physcon, only:pi
 use kernel,  only:radkern
 real, intent(in)     :: dt,xyzh(:,:),origEU(:,:)
 real, intent(inout)  :: radprop(:,:),rad(:,:),vxyzu(:,:)
 logical, intent(out) :: failed,moresweep
 integer, intent(out) :: nit,ierr
 real, intent(out)    :: errorE,errorU,EU0(:,:)
 integer, allocatable :: ivar(:,:),ijvar(:)
 integer              :: ncompact,ncompactlocal,npart,icompactmax,nneigh_average,its
 real, allocatable    :: vari(:,:),varij(:,:),varij2(:,:),varinew(:,:)
 real :: maxerrE2,maxerrU2,maxerrE2last,maxerrU2last
 logical :: converged

 failed = .false.
 errorE = 0.
 errorU = 0.
 ierr = 0

 npart = size(xyzh(1,:))
 nneigh_average = int(4./3.*pi*(radkern*hfact)**3) + 1
 icompactmax = int(1.2*nneigh_average*npart)
 allocate(ivar(3,npart),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for ivar')
 allocate(ijvar(icompactmax),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for ijvar')
 allocate(vari(2,npart),varij(4,icompactmax),varij2(3,icompactmax),varinew(3,npart),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for vari, varij, varij2, varinew')

 !dtimax = dt/imaxstep
 call get_compacted_neighbour_list(xyzh,ivar,ijvar,ncompact,ncompactlocal)
 ! check for errors
 if (ncompact <= 0 .or. ncompactlocal <= 0) then
    call error('radiation_implicit','empty neighbour list - need to call set_linklist first?')
    ierr = ierr_neighbourlist_empty
    return
 endif
 call fill_arrays(ncompact,ncompactlocal,npart,icompactmax,dt,&
                  xyzh,vxyzu,ivar,ijvar,radprop,rad,vari,varij,varij2,EU0)

 maxerrE2last = huge(0.)
 maxerrU2last = huge(0.)

 iterations: do its=1,itsmax_rad
    call compute_flux(ivar,ijvar,ncompact,npart,icompactmax,varij,varij2,vari,EU0,varinew,radprop)
    call calc_lambda_and_eddington(ivar,ncompactlocal,npart,vari,EU0,radprop,ierr)
    call calc_diffusion_term(ivar,ijvar,varij,ncompact,npart,icompactmax,radprop,vari,EU0,varinew,ierr)
    call update_gas_radiation_energy(ivar,ijvar,vari,ncompact,npart,ncompactlocal,&
                                     vxyzu,radprop,rad,origEU,varinew,EU0,moresweep,maxerrE2,maxerrU2)

    if (iverbose >= 2) print*,'iteration: ',its,' error = ',maxerrE2,maxerrU2
    converged = (maxerrE2 <= tol_rad .and. maxerrU2 <= tol_rad)
    if (converged) exit iterations

    maxerrU2last = maxerrU2
 enddo iterations

 if (converged) then
    if (iverbose >= 0) print "(1x,a,i4,a,es10.3,a,es10.3)", &
          trim(label)//': succeeded with ',its,' iterations: xi err:',maxerrE2,' u err:',maxerrU2
 else
    call warning('radiation_implicit','maximum iterations reached')
    moresweep = .true.
 endif

 call store_radiation_results(ncompactlocal,npart,ivar,EU0,rad,vxyzu)

end subroutine do_radiation_onestep


!---------------------------------------------------------
!+
!  get compacted neighbour list
!+
!---------------------------------------------------------
subroutine get_compacted_neighbour_list(xyzh,ivar,ijvar,ncompact,ncompactlocal)
 use dim,      only:periodic,maxphase,maxp
 use linklist, only:ncells,get_neighbour_list,listneigh,ifirstincell
 use kdtree,   only:inodeparts,inoderange
 use boundary, only:dxbound,dybound,dzbound
 use part,     only:iphase,igas,get_partinfo
 use kernel,   only:radkern2
 use io,       only:fatal
 real, intent(in)                  :: xyzh(:,:)
 integer, intent(out)              :: ivar(:,:),ijvar(:)
 integer, intent(out)              :: ncompact,ncompactlocal
 integer                           :: icell,nneigh,i,j,k,n,ip,ncompact_private
 integer                           :: icompact_private,icompact,icompactmax,iamtypei,nneigh_trial
 integer, parameter                :: maxcellcache = 10000
 integer, save, allocatable        :: neighlist(:)
 real                              :: dx,dy,dz,hi21,rij2,q2i
 real, save, allocatable           :: xyzcache(:,:)
 !real, save                        :: xyzcache(maxcellcache,3)
 !$omp threadprivate(xyzcache,neighlist)
 logical                           :: iactivei,iamdusti,iamgasi

 if (.not. allocated(neighlist)) then
    !$omp parallel
    allocate(neighlist(size(xyzh(1,:))),xyzcache(maxcellcache,3))
    !$omp end parallel
 endif

 ncompact = 0
 ncompactlocal = 0
 icompact = 0
 icompactmax = size(ijvar)
 !$omp parallel do schedule(runtime)&
 !$omp shared(ncells,xyzh,inodeparts,inoderange,iphase,dxbound,dybound,dzbound,ifirstincell)&
 !$omp shared(ncompact,icompact,icompactmax)&
 !$omp private(icell,i,j,k,n,ip,iactivei,iamgasi,iamdusti,iamtypei,dx,dy,dz,rij2,q2i)&
 !$omp private(hi21,ncompact_private,icompact_private,nneigh_trial,nneigh)

 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)

    !--skip empty cells AND inactive cells
    if (i <= 0) cycle over_cells

    !
    !--get the neighbour list and fill the cell cache
    !
    call get_neighbour_list(icell,listneigh,nneigh_trial,xyzh,xyzcache,maxcellcache)

    over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
       i = inodeparts(ip)

       if (maxphase==maxp) then
          call get_partinfo(iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
       else
          iactivei = .true.
          iamtypei = igas
          iamdusti = .false.
          iamgasi  = .true.
       endif

       if (.not.iactivei .or. .not.iamgasi) then ! skip if particle is inactive or not gas
          cycle over_parts
       endif
       nneigh = 0
       hi21 = 1./xyzh(4,i)**2

       loop_over_neigh: do n = 1,nneigh_trial

          j = listneigh(n)
          !--do self contribution separately to avoid problems with 1/sqrt(0.)
          if (j==i) cycle loop_over_neigh

          if (n <= maxcellcache) then
             ! positions from cache are already mod boundary
             dx = xyzh(1,i) - xyzcache(n,1)
             dy = xyzh(2,i) - xyzcache(n,2)
             dz = xyzh(3,i) - xyzcache(n,3)
          else
             dx = xyzh(1,i) - xyzh(1,j)
             dy = xyzh(2,i) - xyzh(2,j)
             dz = xyzh(3,i) - xyzh(3,j)
          endif
          if (periodic) then
             if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
             if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
             if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
          endif
          rij2 = dx*dx + dy*dy + dz*dz
          q2i = rij2*hi21
          !
          !--do interaction if r/h < compact support size
          !
          is_sph_neighbour: if (q2i < radkern2) then ! .or. q2j < radkern2) then
             nneigh = nneigh + 1
             neighlist(nneigh) = j
          endif is_sph_neighbour

       enddo loop_over_neigh

!$omp critical(listcompact)
       ncompact = ncompact + 1
       ncompact_private = ncompact
       icompact_private = icompact
       icompact = icompact + nneigh
!$omp end critical (listcompact)
       if (icompact_private+nneigh > icompactmax) then
          print*,'i=',i,'nneigh=',nneigh,'desired size=',icompact_private+nneigh,' actual size=',icompactmax
          call fatal('radiation-implicit','not enough memory allocated for neighbour list', &
                     var='icompactmax',ival=icompactmax)
       endif
       ivar(1,ncompact_private) = nneigh
       ivar(2,ncompact_private) = icompact_private
       ivar(3,ncompact_private) = i

       do k = 1, nneigh
          j = neighlist(k)
          ijvar(icompact_private + k) = j
       enddo
    enddo over_parts
 enddo over_cells
 !$omp end parallel do

 ncompactlocal = ncompact

end subroutine get_compacted_neighbour_list


!---------------------------------------------------------
!+
!  fill arrays
!+
!---------------------------------------------------------
subroutine fill_arrays(ncompact,ncompactlocal,npart,icompactmax,dt,xyzh,vxyzu,ivar,ijvar,radprop,rad,vari,varij,varij2,EU0)
 use dim,             only:periodic
 use boundary,        only:dxbound,dybound,dzbound
 use part,            only:dust_temp,nucleation,gradh,dvdx
 use units,           only:get_c_code
 use kernel,          only:grkern,cnormk
 integer, intent(in) :: ncompact,ncompactlocal,icompactmax,npart
 integer, intent(in) :: ivar(:,:),ijvar(:)
 real, intent(in)    :: dt,xyzh(:,:),vxyzu(:,:),rad(:,:)
 real, intent(inout) :: radprop(:,:)
 real, intent(out)   :: vari(:,:),EU0(2,npart),varij(4,icompactmax),varij2(3,icompactmax)
 integer             :: n,i,j,k,icompact
 real :: cv_effective,pmi,hi,hi21,hi41,rhoi,dx,dy,dz,rij2,rij,rij1,dr,dti,&
         dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi,&
         pmj,rhoj,hj,hj21,hj41,v2i,vi,v2j,vj,dWi,dWj,dvx,dvy,dvz,rhomean,&
         dvdotdr,dv,vmu,dvdWimj,dvdWimi,dvdWjmj,c_code,&
         dWidrlightrhorhom,dWiidrlightrhorhom,dWjdrlightrhorhom,&
         pmjdWrijrhoi,pmjdWrunix,pmjdWruniy,pmjdWruniz,&
         dust_kappai,dust_cooling,heatingISRi,dust_gas

 c_code = get_c_code()
 dti = dt
 !$omp parallel do default(none) &
 !$omp shared(EU0,radprop,rad,xyzh,vxyzu,c_code,vari,ivar,ijvar,varij,varij2,dvdx,dxbound,dybound,dzbound) &
 !$omp shared(dust_temp,ncompactlocal,ncompact,massoftype,iopacity_type,nucleation,dt,gradh,cv_type) &
 !$omp firstprivate(dti) &
 !$omp private(n,i,j,k,rhoi,icompact,pmi,dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi) &
 !$omp private(dvzdxi,dvzdyi,dvzdzi,dx,dy,dz,rij2,rij,rij1,dr,pmj,rhoj,hi,hj,hi21,hj21,hi41,hj41) &
 !$omp private(v2i,vi,v2j,vj,dWi,dWj,dvx,dvy,dvz,rhomean,dvdotdr,dv,vmu,dvdWimj,dvdWimi,dvdWjmj) &
 !$omp private(dWidrlightrhorhom,pmjdWrijrhoi,dWjdrlightrhorhom,dWiidrlightrhorhom,cv_effective) &
 !$omp private(pmjdWrunix,pmjdWruniy,pmjdWruniz,dust_kappai,dust_cooling,heatingISRi,dust_gas)

 do n = 1,ncompact
    i = ivar(3,n)
    !  if (iphase(i) == 0) then
    EU0(1,i) = rad(iradxi,i)
    EU0(2,i) = vxyzu(4,i)
    rhoi = rhoh(xyzh(4,i), massoftype(igas))
    radprop(icv,i) = get_cv(rhoi,vxyzu(4,i),cv_type)
    radprop(ikappa,i) = get_kappa(iopacity_type,vxyzu(4,i),radprop(icv,i),rhoi)
    !
    !--Diffuse ISM: Set dust temperature and opacity
    !
    if (dustRT .and. n<=ncompactlocal) then
       dust_temp(i) = dust_temperature(rad(iradxi,i),vxyzu(4,i),rhoi,dust_kappai,dust_cooling,heatingISRi,dust_gas)
       nucleation(idkappa,i) = dust_kappai
    endif
    !
    !--Note that CV and Kappa have already been done in ASS
    !
    cv_effective = radprop(icv,i)/get_1overmu(rhoi,vxyzu(4,i),cv_type)
    dvxdxi = 0.
    dvxdyi = 0.
    dvxdzi = 0.
    dvydxi = 0.
    dvydyi = 0.
    dvydzi = 0.
    dvzdxi = 0.
    dvzdyi = 0.
    dvzdzi = 0.

    pmi = massoftype(igas)
    hi = xyzh(4,i)
    hi21 = 1./(hi*hi)
    hi41 = hi21*hi21
    rhoi = rhoh(xyzh(4,i), massoftype(igas))

    do k = 1,ivar(1,n) ! Looping from 1 to nneigh
       icompact = ivar(2,n) + k
       j = ijvar(icompact)
       !
       !--Need to make sure that E and U values are loaded for non-active neighbours
       !
       EU0(1,j) = rad(iradxi,j)
       EU0(2,j) = vxyzu(4,j)
       !
       !--Note that CV and Kappa have already been done in ASS
       !
       rhoj = rhoh(xyzh(4,j), massoftype(igas))
       cv_effective = radprop(icv,j)/get_1overmu(rhoj,vxyzu(4,j),cv_type)
       !dti = dt
       !
       !--Calculate other quantities
       !
       dx = xyzh(1,i) - xyzh(1,j)
       dy = xyzh(2,i) - xyzh(2,j)
       dz = xyzh(3,i) - xyzh(3,j)
       if (periodic) then
          if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
          if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
          if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
       endif
       rij2 = dx*dx + dy*dy + dz*dz + tiny(0.)
       rij = sqrt(rij2)
       rij1 = 1./rij
       dr = rij

       pmj = massoftype(igas)

       hj = xyzh(4,j)
       hj21 = 1./(hj*hj)
       hj41 = hj21*hj21

       v2i = rij2*hi21
       vi = rij/hi

       v2j = rij2*hj21
       vj = rij/hj

       dWi = grkern(v2i,vi)*hi41*cnormk*gradh(1,i)
       dWj = grkern(v2j,vj)*hj41*cnormk*gradh(1,j)

       dvx = vxyzu(1,i) - vxyzu(1,j)
       dvy = vxyzu(2,i) - vxyzu(2,j)
       dvz = vxyzu(3,i) - vxyzu(3,j)

       dvdotdr = dvx*dx + dvy*dy + dvz*dz
       dv = dvdotdr/dr

       if (dvdotdr > 0.) then
          vmu = 0.
       else
          vmu = dv
       endif

       ! Coefficients in radiative flux term in radiation energy density equation (e.g. eq 22 & 25, Whitehouse & Bate 2004)
       dvdWimj = pmj*dv*dWi
       dvdWimi = pmi*dv*dWi
       dvdWjmj = pmj*dv*dWj

       ! Coefficients for p(div(v))/rho term in gas energy equation (e.g. eq 26, Whitehouse & Bate 2004)
       dWidrlightrhorhom = c_code*dWi/dr*pmj/(rhoi*rhoj)
       dWiidrlightrhorhom = c_code*dWi/dr*pmi/(rhoi*rhoj)
       dWjdrlightrhorhom = c_code*dWj/dr*pmj/(rhoi*rhoj)

       pmjdWrijrhoi = pmj*dWi*rij1/rhoi
       pmjdWrunix = pmjdWrijrhoi*dx
       pmjdWruniy = pmjdWrijrhoi*dy
       pmjdWruniz = pmjdWrijrhoi*dz
       !
       !--Calculates density(i) times the gradient of velocity
       !
       dvxdxi = dvxdxi - dvx*pmjdWrunix
       dvxdyi = dvxdyi - dvx*pmjdWruniy
       dvxdzi = dvxdzi - dvx*pmjdWruniz
       dvydxi = dvydxi - dvy*pmjdWrunix
       dvydyi = dvydyi - dvy*pmjdWruniy
       dvydzi = dvydzi - dvy*pmjdWruniz
       dvzdxi = dvzdxi - dvz*pmjdWrunix
       dvzdyi = dvzdyi - dvz*pmjdWruniy
       dvzdzi = dvzdzi - dvz*pmjdWruniz

       varij(1,icompact) = rhoj
       varij(2,icompact) = dWiidrlightrhorhom
       varij(3,icompact) = dWidrlightrhorhom
       varij(4,icompact) = dWjdrlightrhorhom

       varij2(1,icompact) = pmjdWrunix
       varij2(2,icompact) = pmjdWruniy
       varij2(3,icompact) = pmjdWruniz
    enddo
    dvdx(1,i) = real(dvxdxi,kind=kind(dvdx)) ! convert to real*4 explicitly to avoid warnings
    dvdx(2,i) = real(dvxdyi,kind=kind(dvdx))
    dvdx(3,i) = real(dvxdzi,kind=kind(dvdx))
    dvdx(4,i) = real(dvydxi,kind=kind(dvdx))
    dvdx(5,i) = real(dvydyi,kind=kind(dvdx))
    dvdx(6,i) = real(dvydzi,kind=kind(dvdx))
    dvdx(7,i) = real(dvzdxi,kind=kind(dvdx))
    dvdx(8,i) = real(dvzdyi,kind=kind(dvdx))
    dvdx(9,i) = real(dvzdzi,kind=kind(dvdx))

    vari(1,n) = dti
    vari(2,n) = rhoi
    ! endif
 enddo
 !$omp end parallel do

end subroutine fill_arrays


!---------------------------------------------------------
!+
!  compute radiative flux
!+
!---------------------------------------------------------
subroutine compute_flux(ivar,ijvar,ncompact,npart,icompactmax,varij,varij2,vari,EU0,varinew,radprop)
 integer, intent(in) :: ivar(:,:),ijvar(:),ncompact,npart,icompactmax
 real, intent(in)    :: varij(4,icompactmax),varij2(3,icompactmax),vari(2,npart),EU0(2,npart)
 real, intent(inout) :: radprop(:,:)
 real, intent(out)   :: varinew(3,npart)  ! we use this parallel loop to set varinew to zero
 integer             :: i,j,k,n,icompact
 real                :: rhoi,rhoj,pmjdWrunix,pmjdWruniy,pmjdWruniz,dedxi,dedyi,dedzi,dradenij

 !$omp parallel do default(none)&
 !$omp shared(vari,ivar,EU0,varij2,ijvar,varij,ncompact,radprop,varinew)&
 !$omp private(i,j,k,n,dedxi,dedyi,dedzi,rhoi,rhoj,icompact)&
 !$omp private(pmjdWrunix,pmjdWruniy,pmjdWruniz,dradenij)

 do n = 1,ncompact
    i = ivar(3,n)
    varinew(1,i) = 0.
    varinew(2,i) = 0.
    !varinew(3,i) = 0.
    !          if (iphase(i)==0) then
    dedxi = 0.
    dedyi = 0.
    dedzi = 0.

    rhoi = vari(2,n)

    do k = 1,ivar(1,n)
       icompact = ivar(2,n) + k
       j = ijvar(icompact)
       rhoj = varij(1,icompact)
       pmjdWrunix = varij2(1,icompact)
       pmjdWruniy = varij2(2,icompact)
       pmjdWruniz = varij2(3,icompact)

       ! Calculates the gradient of E (where E=rho*e, and e is xi)
       dradenij = rhoj*EU0(1,j) - rhoi*EU0(1,i)
       dedxi = dedxi + dradenij*pmjdWrunix
       dedyi = dedyi + dradenij*pmjdWruniy
       dedzi = dedzi + dradenij*pmjdWruniz
    enddo

    radprop(ifluxx,i) = dedxi
    radprop(ifluxy,i) = dedyi
    radprop(ifluxz,i) = dedzi
 enddo
 !$omp end parallel do

end subroutine compute_flux


!---------------------------------------------------------
!+
!  calculate flux limiter (lambda) and eddington factor
!+
!---------------------------------------------------------
subroutine calc_lambda_and_eddington(ivar,ncompactlocal,npart,vari,EU0,radprop,ierr)
 use io,              only:error
 use part,            only:dust_temp,nucleation
 use radiation_utils, only:get_rad_R
 use options,         only:limit_radiation_flux
 integer, intent(in) :: ivar(:,:),ncompactlocal,npart
 real, intent(in)    :: vari(:,:),EU0(2,npart)
 real, intent(inout) :: radprop(:,:)
 integer             :: n,i,ierr
 real                :: rhoi,gradE1i,opacity,radRi

 ierr = 0
 !$omp parallel do default(none)&
 !$omp shared(vari,ivar,radprop,ncompactlocal,EU0,dust_temp,nucleation,limit_radiation_flux) &
 !$omp private(i,n,rhoi,gradE1i,opacity,radRi)&
 !$omp reduction(max:ierr)
 do n = 1,ncompactlocal
    i = ivar(3,n)
    rhoi = vari(2,n)
    !
    ! If using diffuse ISM, use Rosseland mean opacity from the frequency
    ! dependent opacity when dust temperatures are cold (T_d<100 K).
    ! Otherwise use the tabulated grey dust opacities.
    !
    opacity = radprop(ikappa,i)
    if (dustRT) then
       if (dust_temp(i) < Tdust_threshold) opacity = nucleation(idkappa,i)
    endif
    if (opacity < 0.) then
       ierr = max(ierr,ierr_negative_opacity)
       call error(label,'Negative opacity',val=opacity)
    endif

    if (limit_radiation_flux) then
       radRi = get_rad_R(rhoi,EU0(1,i),radprop(ifluxx:ifluxz,i),opacity)
    else
       radRi = 0.
    endif
    radprop(ilambda,i) = (2. + radRi ) / (6. + 3.*radRi + radRi**2)  ! Levermore & Pomraning's flux limiter (e.g. eq 12, Whitehouse & Bate 2004)
    radprop(iedd,i) = radprop(ilambda,i) + radprop(ilambda,i)**2 * radRi**2  ! e.g., eq 11, Whitehouse & Bate (2004)
 enddo
 !$omp end parallel do

end subroutine calc_lambda_and_eddington


!---------------------------------------------------------
!+
!  calculate diffusion coefficients
!+
!---------------------------------------------------------
subroutine calc_diffusion_term(ivar,ijvar,varij,ncompact,npart,icompactmax, &
                               radprop,vari,EU0,varinew,ierr)
 use io,   only:error
 use part, only:dust_temp,nucleation
 integer, intent(in)  :: ivar(:,:),ijvar(:),ncompact,npart,icompactmax
 real, intent(in)     :: vari(:,:),varij(4,icompactmax),EU0(2,npart),radprop(:,:)
 integer, intent(out) :: ierr
 real, intent(inout)  :: varinew(3,npart)
 integer              :: n,i,j,k,icompact
 real                 :: rhoi,rhoj,opacityi,opacityj,bi,bj,b1
 real                 :: dWiidrlightrhorhom,dWidrlightrhorhom,dWjdrlightrhorhom
 real                 :: diffusion_numerator,diffusion_denominator,tempval1,tempval2

 ierr = 0
 !$omp parallel do default(none)&
 !$omp shared(vari,varij,ivar,ijvar,radprop,ncompact,EU0,varinew,dust_temp,nucleation)&
 !$omp private(i,j,k,n,rhoi,rhoj,opacityi,opacityj,bi,bj,b1,diffusion_numerator,diffusion_denominator)&
 !$omp private(dWiidrlightrhorhom,dWidrlightrhorhom,dWjdrlightrhorhom,tempval1,tempval2,icompact)&
 !$omp reduction(max:ierr)
 do n = 1,ncompact
    i = ivar(3,n)
    !  if (iphase(i) == 0) then
    rhoi = vari(2,n)
    !
    !--NOTE: Needs to do this loop even for boundaryparticles because active
    !     boundary particles will need to contribute to the varinew()
    !     quantities (i.e. diffusion terms) of particle j due to the way that
    !     particle j only finds neighbours inside h_j or non-active particles
    !     inside h_i.  The varinew() quantities of a boundaryparticle are
    !     not used, but its contributions to j are.
    !
    !--Initialising counters to zero for this particle
    !
    diffusion_numerator = 0.
    diffusion_denominator = 0.
    !
    !--All the neighbours loop
    !
    do k = 1,ivar(1,n)
       icompact = ivar(2,n) + k
       j = ijvar(icompact)
       rhoj = varij(1,icompact)
       dWiidrlightrhorhom = varij(2,icompact)
       dWidrlightrhorhom = varij(3,icompact)
       dWjdrlightrhorhom = varij(4,icompact)
       !
       !--Set c*lambda/kappa*rho term (radiative diffusion coefficient) for current quantities
       !
       opacityi = radprop(ikappa,i)
       opacityj = radprop(ikappa,j)
       if (dustRT) then
          if (dust_temp(i) < Tdust_threshold) opacityi = nucleation(idkappa,i)
          if (dust_temp(j) < Tdust_threshold) opacityj = nucleation(idkappa,j)
       endif
       if ((opacityi <= 0.) .or. (opacityj <= 0.)) then
          ierr = max(ierr,ierr_negative_opacity)
          call error(label,'Negative or zero opacity',val=min(opacityi,opacityj))
       endif
       bi = radprop(ilambda,i)/(opacityi*rhoi)
       bj = radprop(ilambda,j)/(opacityj*rhoj)
       !
       !--Choose the 'average' diffusion value.  The (bi+bj) quantity biased in
       !     favour of the particle with the lowest opacity.  The other average
       !     is that original recommended in Cleary & Monaghan for heat diffusion.
       b1 = bi + bj
       !
       !--Diffusion numerator and denominator
       !
       diffusion_numerator = diffusion_numerator - 0.5*dWidrlightrhorhom*b1*EU0(1,j)*rhoj
       diffusion_denominator = diffusion_denominator + 0.5*dWidrlightrhorhom*b1*rhoi
       !
       !--If particle j is active, need to add contribution due to i for hj
       !
       !  if (iactive(iphase(j))) then  !
       tempval1 = 0.5*dWiidrlightrhorhom*b1
       tempval2 = tempval1*rhoj
       tempval1 = tempval1*EU0(1,i)*rhoi
       !$omp atomic
       varinew(1,j) = varinew(1,j) - tempval1
       !$omp atomic
       varinew(2,j) = varinew(2,j) + tempval2
       !  else
       !     diffusion_numerator = diffusion_numerator - 0.5*dWjdrlightrhorhom*b1*EU0(1,J)*rhoj
       !     diffusion_denominator = diffusion_denominator + 0.5*dWjdrlightrhorhom*b1*rhoi
       !  endif
       ! ENDIF         !--iphase(j)==0
    enddo
    !$omp atomic
    varinew(1,i) = varinew(1,i) + diffusion_numerator
    !$omp atomic
    varinew(2,i) = varinew(2,i) + diffusion_denominator
 enddo
 !$omp end parallel do

end subroutine calc_diffusion_term

!---------------------------------------------------------
!+
!  update gas and radiation energy
!+
!---------------------------------------------------------
subroutine update_gas_radiation_energy(ivar,ijvar,vari,ncompact,npart,ncompactlocal,&
                                       vxyzu,radprop,rad,origEU,varinew,EU0,moresweep,maxerrE2,maxerrU2)
 use io,      only:fatal,error
 use part,    only:pdvvisc=>luminosity,dvdx,nucleation,dust_temp,eos_vars,drad,iradxi,fxyzu
 use units,   only:get_radconst_code,get_c_code,unit_density
 use physcon, only:mass_proton_cgs
 use eos,     only:metallicity=>Z_in
 use options, only:implicit_radiation_store_drad
 integer, intent(in) :: ivar(:,:),ijvar(:),ncompact,npart,ncompactlocal
 real, intent(in)    :: vari(:,:),varinew(3,npart),rad(:,:),origEU(:,:),vxyzu(:,:)
 real, intent(inout) :: radprop(:,:),EU0(2,npart)
 real, intent(out)   :: maxerrE2,maxerrU2
 logical, intent(out):: moresweep
 integer             :: i,j,n,ieqtype,ierr
 logical             :: moresweep2,skip_quartic
 real                :: dti,rhoi,diffusion_numerator,diffusion_denominator,gradEi2,gradvPi,rpdiag,rpall
 real                :: radpresdenom,stellarradiation,gas_temp,xnH2,betaval,gammaval,tfour,betaval_d,chival
 real                :: gas_dust_val,dust_tempi,dust_kappai,a_code,c_code,dustgammaval,gas_dust_cooling
 real                :: cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term,e_planetesimali
 real                :: u4term,u1term,u0term,pcoleni,dust_cooling,heatingISRi,dust_gas
 real                :: pres_numerator,pres_denominator,mui,U1i,E1i,Tgas,dUcomb,dEcomb
 real                :: residualE,residualU,xchange,maxerrU2old,Tgas4,Trad4,ck,ack

 a_code = get_radconst_code()
 c_code = get_c_code()
 moresweep = .false.
 maxerrE2 = 0.
 maxerrU2 = 0.

 !$omp parallel do default(none)&
 !$omp shared(vari,ivar,ijvar,radprop,rad,ncompact,ncompactlocal,EU0,varinew) &
 !$omp shared(dvdx,origEU,nucleation,dust_temp,eos_vars,implicit_radiation_store_drad,cv_type) &
 !$omp shared(moresweep,pdvvisc,metallicity,vxyzu,iopacity_type,a_code,c_code,massoftype,drad,fxyzu) &
 !$omp private(i,j,n,rhoi,dti,diffusion_numerator,diffusion_denominator,U1i,skip_quartic,Tgas,E1i,dUcomb,dEcomb) &
 !$omp private(gradEi2,gradvPi,rpdiag,rpall,radpresdenom,stellarradiation,dust_tempi,dust_kappai,xnH2) &
 !$omp private(dust_cooling,heatingISRi,dust_gas,gas_dust_val,dustgammaval,gas_dust_cooling,cosmic_ray) &
 !$omp private(cooling_line,photoelectric,h2form,dust_heating,dust_term,betaval,chival,gammaval,betaval_d,tfour) &
 !$omp private(e_planetesimali,u4term,u1term,u0term,pcoleni,pres_numerator,pres_denominator,moresweep2,mui,ierr) &
 !$omp private(residualE,residualU,xchange,maxerrU2old,gas_temp,ieqtype,unit_density,Tgas4,Trad4,ck,ack) &
 !$omp reduction(max:maxerrE2,maxerrU2)
 main_loop: do n = 1,ncompactlocal
    i = ivar(3,n)

    ! if (iphase(i)==0) then
    dti = vari(1,n)
    rhoi = vari(2,n)
    !  if (.NOT.boundaryparticle(i,xyzmh,rhoi)) then
    diffusion_numerator = varinew(1,i)
    diffusion_denominator = varinew(2,i)
    pres_numerator = pdvvisc(i)/massoftype(igas) ! in phantom pdvvisc->luminosity which is m*du/dt not du/dt
    pres_denominator = 0.
    !
    !--Radiation pressure...
    !
    gradEi2 = dot_product(radprop(ifluxx:ifluxz,i),radprop(ifluxx:ifluxz,i))

    if (gradEi2 < tiny(0.)) then
       gradvPi = 0.
    else
       rpdiag = 0.5*(1.-radprop(iedd,i))  ! Diagonal component of Eddington tensor (eq 10, Whitehouse & Bate 2004)
       rpall = 0.5*(3.*radprop(iedd,i)-1.)/gradEi2  ! n,n-component of Eddington tensor, where n is the direction of grad(E) (or -ve flux)
       gradvPi = (((rpdiag+rpall*radprop(ifluxx,i)**2)*dvdx(1,i))+ &
                 ((rpall*radprop(ifluxx,i)*radprop(ifluxy,i))*dvdx(2,i))+ &
                 ((rpall*radprop(ifluxx,i)*radprop(ifluxz,i))*dvdx(3,i))+ &
                 ((rpall*radprop(ifluxy,i)*radprop(ifluxx,i))*dvdx(4,i))+ &
                 ((rpdiag+rpall*radprop(ifluxy,i)**2)*dvdx(5,i))+ &
                 ((rpall*radprop(ifluxy,i)*radprop(ifluxz,i))*dvdx(6,i))+ &
                 ((rpall*radprop(ifluxz,i)*radprop(ifluxx,i))*dvdx(7,i))+ &
                 ((rpall*radprop(ifluxz,i)*radprop(ifluxy,i))*dvdx(8,i))+ &
                 ((rpdiag+rpall*radprop(ifluxz,i)**2)*dvdx(9,i)))  ! e.g. eq 23, Whitehouse & Bate (2004)
    endif

    radpresdenom = gradvPi * EU0(1,i)

    stellarradiation = 0. ! set to zero
    e_planetesimali = 0.
    pcoleni = 0. ! specific collision energy from planetesimals

    !
    !--For idustRT>0, replace the gas-dust coupling term in the u equation
    !     by the dust-radiation term, but keep the T_d rather than T_g
    !
    if (dustRT) then
       radprop(ikappa,i) = get_kappa(iopacity_type,EU0(2,i),radprop(icv,i),rhoi)
       dust_tempi = dust_temperature(rad(iradxi,i),EU0(2,i),rhoi,dust_kappai,dust_cooling,heatingISRi,dust_gas)
       gas_temp = EU0(2,i)/radprop(icv,i)
       mui = eos_vars(imu,i)
       xnH2 = rhoi*unit_density/(mui*mass_proton_cgs)  ! Mike: Check units
    endif

    skip_quartic = .false.
    if (dustRT .and. &
        ( abs(dust_tempi-(rhoi*EU0(1,i)/a_code)**0.25) > 1. .and. xnH2 < 1.e11/metallicity) ) then
       !--For low densities and temperatures, use the form of the equations that
       !     includes the gas-dust coupling term.  Also explicitly set the gas
       !     opacity to zero (i.e. betaval=gammaval=tfour=0).  This works well
       !     until the gas-dust coupling term gets very large, where upon the
       !     convergence fails because even a small difference (<1 K) in the
       !     gas and dust temperatures results in a large term.  The gas-dust
       !     coupling term potentially becomes large for high densities, but
       !     a density criterion alone is not sufficient.  What we really want
       !     to see is that the matter is well coupled with the radiation field
       !     already.  So if we have high densities AND the difference between
       !     the temperatures of the dust and local radiation field is small
       !     (<1 K), then we abandon the gas-dust coupling term.
       !
       call set_heating_cooling_low_rhoT(i,EU0(1,i),EU0(2,i),origEU(1,i),origEU(2,i),&
                                         radprop(icv,i),dti,diffusion_denominator,&
                                         pres_numerator,radpresdenom,rhoi,xnH2,heatingISRi,e_planetesimali,&
                                         metallicity,gas_temp,ieqtype,betaval,betaval_d,gammaval,&
                                         chival,tfour,dust_tempi,gas_dust_val,dustgammaval,gas_dust_cooling,&
                                         cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term,skip_quartic,U1i,ierr)
       if (ierr > 0) then
          !$omp critical (moresweepset)
          moresweep = .true.
          !$omp end critical (moresweepset)
          cycle main_loop
       endif

    elseif (dustRT) then
       !--Replaces the gas-dust coupling term in the u equation by the dust-radiation
       !     term and assumes that T_g=T_d so that the dust-radiation term is
       !     actually the gas-radiation term (i.e. uses kappa from opacity tables
       !     which is based on the gas temperature, and uses (u/cv) rather than T_d.)
       call set_heating_cooling(i,EU0(2,i),radprop(icv,i),rhoi,mui,heatingISRi,metallicity,ieqtype,dust_tempi,&
                                gas_dust_val,dustgammaval,gas_dust_cooling,&
                                cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term)

    else
       !--Else, this is the original version of radiative transfer
       !     (Whitehouse & Bate 2006), which does not include the
       !     diffuse ISM model of Bate (2015).
       call turn_heating_cooling_off(ieqtype,dust_tempi,gas_dust_val,dustgammaval,gas_dust_cooling,&
                                     cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term)
    endif
    !
    !--Now solve those equations... (these are eqns 22 in Whitehouse, Bate & Monaghan 2005)
    !
    Tgas4 = (EU0(2,i)/radprop(icv,i))**4
    Trad4 = rhoi*EU0(1,i)/a_code
    ck  = c_code*radprop(ikappa,i)
    ack = a_code*ck

    betaval  = ck*rhoi*dti
    chival   = dti*(diffusion_denominator-radpresdenom/EU0(1,i))-betaval
    gammaval = ack/radprop(icv,i)**4
    tfour    = ack*(Trad4 - Tgas4)
    u4term   = gammaval*dti*(dti*(diffusion_denominator-radpresdenom/EU0(1,i)) - 1.)
    u1term   = (chival-1.)*(1.-dti*pres_denominator + dti*gas_dust_val/radprop(icv,i)) &
               - betaval*dti*gas_dust_val/radprop(icv,i)
    u0term   = betaval*(origEU(1,i) - dti*gas_dust_val*dust_tempi + dti*dust_heating) + &
               (chival-1.)*(-origEU(2,i) - dti*pres_numerator - dti*e_planetesimali &
               - dti*gas_dust_val*dust_tempi - dti*cosmic_ray + dti*cooling_line - dti*photoelectric &
               - dti*h2form + dti*dust_term) + dti*diffusion_numerator*betaval &
               + stellarradiation*betaval - (chival-1.)*pcoleni

    if (u1term > 0. .and. u0term > 0. .or. u1term < 0. .and. u0term < 0.) then
       !$omp critical(quart)
       print *,"ngs ",u4term,u1term,u0term,betaval,chival,gammaval
       print *,"    ",radprop(ikappa,i),rhoi,dti
       print *,"    ",diffusion_denominator,diffusion_numerator
       print *,"    ",pres_denominator,pres_numerator !,uradconst
       print *,"    ",radpresdenom,EU0(1,i),EU0(2,i) !,ekcle(3,i)
       print *,"    ",c_code,origEU(1,i),origEU(2,i)
       !$omp end critical(quart)
       !$omp critical (moresweepset)
       moresweep = .true.
       !$omp end critical (moresweepset)
       cycle main_loop
    endif

    if (.not. skip_quartic) then
       u1term = u1term/u4term
       u0term = u0term/u4term
       moresweep2 = .false.
       call solve_quartic(u1term,u0term,EU0(2,i),U1i,moresweep2,ierr)  ! U1i is the quartic solution
       if (ierr /= 0) then
          print*,'Error in solve_quartic'
          print*,'i=',i,'u1term=',u1term,'u0term=',u0term,'EU0(2,i)=',EU0(2,i),'U1i=',U1i,'moresweep=',moresweep
          print*,"info: ",EU0(2,i)/radprop(icv,i)
          print*,"info2: ",u0term,u1term,u4term,gammaval,radprop(ikappa,i),radprop(icv,i)
          print*,"info3: ",chival,betaval,dti,rhoi
          print*,"info4: ",pres_denominator,origEU(1,i),pres_numerator
          print*,"info5: ",diffusion_numerator,stellarradiation,diffusion_denominator
          print*,"info6: ",radpresdenom,EU0(1,i)
          print*,"Tgas: ",EU0(2,i)/radprop(icv,i)," Trad:",(rhoi*EU0(1,i)/a_code)**0.25,' ack*(Tgas^4 - Trad^4): ',tfour

          call fatal('solve_quartic','Fail to solve')
       endif

       if (moresweep2) then
!$omp critical (moresweepset)
          moresweep = .true.
          print*,"info: ",EU0(2,i)/radprop(icv,i)
          print*,"info2: ",u0term,u1term,u4term,gammaval,radprop(ikappa,i),radprop(icv,i)
          print*,"info3: ",chival,betaval,dti
          print*,"info4: ",pres_denominator,origeu(1,i),pres_numerator
          print*,"info5: ",diffusion_numerator,stellarradiation,diffusion_denominator
          print*,"info6: ",radpresdenom,EU0(1,i)
          print*,"info7: ",cosmic_ray,heatingisri
          print*,"info8: ",cooling_line,photoelectric,h2form
!$omp end critical (moresweepset)
          cycle main_loop
       endif
    endif

    E1i = (origEU(1,i) + dti*diffusion_numerator &
                       + gammaval*dti*U1i**4 &
                       + dustgammaval*dti &
                       + dti*gas_dust_val*(U1i/radprop(icv,i) - dust_tempi) &
                       + dti*dust_heating &
                       + stellarradiation)/(1.-chival)
    dUcomb = pres_numerator + pres_denominator*EU0(2,i) + tfour &
           - gas_dust_cooling + cosmic_ray - cooling_line &
           + photoelectric + h2form + e_planetesimali + pcoleni
    dEcomb = diffusion_numerator + diffusion_denominator * EU0(1,i) &
             - tfour - radpresdenom + stellarradiation + dust_heating + gas_dust_cooling
    !
    !--Tests for negativity
    !
    if (U1i <= 0.) then
!$omp critical (moresweepset)
       print*, "radiation_implicit: u has gone negative ",i,u1term,u0term,u4term,EU0(2,i),U1i,moresweep,ierr
       moresweep=.true.
       print*, "radiation_implicit: u has gone negative ",i,U1i
       print*,'Error in solve_quartic'
       print*,'i=',i,'u1term=',u1term,'u0term=',u0term,'EU0(2,i)=',EU0(2,i),'U1i=',U1i,'moresweep=',moresweep
       print*,"info: ",EU0(2,i)/radprop(icv,i)
       print*,"info2: ",u0term,u1term,u4term,gammaval,radprop(ikappa,i),radprop(icv,i)
       print*,"info3: ",chival,betaval,dti
       print*,"info4: ",pres_denominator,origEU(1,i),pres_numerator
       print*,"info5: ",diffusion_numerator,stellarradiation,diffusion_denominator
       print*,"info6: ",radpresdenom,EU0(1,i)
!$omp end critical (moresweepset)
    endif
    if (E1i <= 0.) then
!$omp critical (moresweepset)
       moresweep=.true.
       call error(label,'e has gone negative',i)
!$omp end critical (moresweepset)
    endif
    !
    ! And the error is...
    !
    Tgas = EU0(2,i)/radprop(icv,i)
    if (Tgas >= 0.) then
       maxerrE2 = max(maxerrE2, 1.*abs((EU0(1,i) - E1i)/E1i))
       residualE = 0.
    else
       xchange = abs((origEU(1,i) + (dEcomb)*dti - E1i)/E1i)
       maxerrE2 = max(maxerre2,xchange)
       residualE = origEU(1,i) + (dEcomb)*dti - E1i
    endif

    if (Tgas >= 2000.) then
       maxerrU2 = max(maxerrU2, 1.*abs((EU0(2,i) - U1i)/U1i))
       residualU = 0.
    else
       maxerrU2old = maxerrU2
       maxerrU2 = max(maxerrU2, abs((origEU(2,i)+(dUcomb)* dti - U1i)/U1i))
       residualU = origEU(2,i)+(dUcomb)*dti - U1i
    endif
    !
    !--Copy values
    !
    EU0(1,i) = E1i
    EU0(2,i) = U1i
    radprop(icv,i) = get_cv(rhoi,EU0(2,i),cv_type)
    radprop(ikappa,i) = get_kappa(iopacity_type,EU0(2,i),radprop(icv,i),rhoi)

    if (implicit_radiation_store_drad) then  ! use this for testing
       drad(iradxi,i) = (E1i - origEU(1,i))/dti  ! dxi/dt
       fxyzu(4,i) = (U1i - origEU(2,i))/dti      ! du/dt
    endif

    if (dustRT) then
       dust_temp(i) = dust_temperature(rad(iradxi,i),EU0(2,i),rhoi,dust_kappai,&
                                       dust_cooling,heatingISRi,dust_gas)
       nucleation(idkappa,i) = dust_kappai
    endif

 enddo main_loop
!$omp end parallel do

end subroutine update_gas_radiation_energy


subroutine set_heating_cooling_low_rhoT(i,eradi,ugasi,orig_eradi,orig_ugasi,cvi,dti,&
           diffusion_denominator,pres_numerator,radpresdenom,rhoi,xnH2,heatingISRi,&
           e_planetesimali,metallicity,gas_temp,ieqtype,betaval,betaval_d,gammaval,&
           chival,tfour,dust_tempi,gas_dust_val,dustgammaval,gas_dust_cooling,&
           cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term,skip_quartic,U1i,ierr)
 use units,   only:unit_pressure,unit_ergg,utime
 use physcon, only:eV
 use eos,     only:get_u_from_rhoT,ieos
 integer, intent(in)  :: i
 real, intent(in)     :: eradi,ugasi,orig_eradi,orig_ugasi,xnH2,metallicity,gas_temp,rhoi
 real, intent(in)     :: heatingISRi,dti,diffusion_denominator,radpresdenom,pres_numerator,e_planetesimali
 real, intent(inout)  :: cvi
 logical, intent(out) :: skip_quartic
 integer, intent(out) :: ieqtype,ierr
 real, intent(inout)  :: dust_tempi
 real, intent(out)    :: betaval,betaval_d,gammaval,tfour,gas_dust_val,dustgammaval,chival,&
                         gas_dust_cooling,cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term,U1i
 integer              :: itry,iterationloop
 real                 :: u_found,t_found,t_orig,u_last,t_last,t_plus,u_plus
 real                 :: tdiff,cooling_line1,cooling_line2,dline_du,h2form1,h2form2,dh2form_dt
 real                 :: photoelectric1,photoelectric2,dphoto_du,func,derivative,func_old,gas_dust_dustT,cv1,h2fraci

 skip_quartic = .false.
 ierr = 0
 U1i = 0.
 ieqtype = 3
 gammaval = 0.
 tfour = 0.
 betaval = 0.
 betaval_d = 0.
 chival = dti*(diffusion_denominator - radpresdenom/eradi) -betaval - betaval_d  ! eradi = EU0(1,i)
 gas_dust_val = 0.
 gas_dust_val = gas_dust_collisional_term(xnH2,metallicity,gas_temp) / rhoi / unit_pressure*utime
 !
 !--Implement crude effect of dust sublimation on gas-dust thermal coupling
 !
 if (gas_temp > 1000.) gas_dust_val = max(gas_dust_val*(1500.-gas_temp)/500., 0.)
 gas_dust_cooling = gas_dust_val*(gas_temp - dust_tempi)
 gas_dust_dustT = gas_dust_val*dust_tempi
 cosmic_ray = cosmic_ray_heating(xnH2) / rhoi / unit_pressure*utime
 !
 !--This adds the dust heating term into the equations.
 !
 dust_heating = heatingISRi/unit_ergg*utime
 cooling_line = cooling_line_rate(i,gas_temp,xnH2,metallicity) / rhoi / unit_pressure*utime
 !
 !--H_2 formation heating
 !
 if (H2formation_heating) then
    h2form = h2_formation(h2fraci,gas_temp,dust_tempi,xnH2,metallicity)*4.48*eV / rhoi / unit_pressure*utime
    !
    !--Potentially add heating from UV destruction and pumping of H_2
    !     (see Bate 2015) for effect of this.
    !
    !     &                 + h2_destruction(i,xnH2,.FALSE.) *
    !     &            (0.4 + 2.0/(1.0+criticaln(i,gas_temp)/(2.0*xnH2)))*
    !     &                 eleccharge/rhoi/uergcc*utime
 else
    h2form = 0.
 endif
 dustgammaval = 0.
 dust_term = 0.
 photoelectric = photoelectric_heating(i,gas_temp,xnH2,metallicity) / rhoi / unit_pressure*utime

 !--In the low-density regime, the line cooling term (and sometimes
 !     other terms) can be very large, so need to solve for U1i using
 !     Newton-Raphson so that it doesn't go negative
 if (.true.) then
    try_loop: do itry = 1,2
       u_found = ugasi  ! ugasi = EU0(2,i)
       t_found = ugasi/get_cv(rhoi,u_found,cv_type)
       t_orig = t_found
       do iterationloop = 1,100
          u_last = u_found
          t_last = t_found

          u_found = get_u_from_rhoT(rhoi,t_found,ieos)
          cv1 = get_cv(rhoi,u_found,cv_type)
          t_found = u_found/cv1
          !
          !--For calculating numerical derivative with gas temperature,
          !     set t_plus to be 1 K higher, or 1% higher if temperature is small
          !
          if (t_found > 10.) then
             t_plus = t_found + 1.
          else
             t_plus = t_found * 1.01
          endif
          tdiff = t_plus - t_found
          u_plus = get_u_from_rhoT(rhoi,t_plus,ieos)
          !--Molecular line cooling.  NOTE that because cooling_line_rate() sets the
          !     abundances of carbon (in the chemistry() array) it it important that
          !     the t_found call is done after the t_plus call otherwise the carbon
          !     abundance that is stored in the chemistry() array is produced using
          !     the wrong temperature!
          cooling_line2 = cooling_line_rate(i,t_plus,xnH2,metallicity)/rhoi/unit_pressure*utime
          cooling_line1 = cooling_line_rate(i,t_found,xnH2,metallicity)/rhoi/unit_pressure*utime

          dline_du = (cooling_line2-cooling_line1)/tdiff
          !
          !--H_2 formation heating
          !
          if (H2formation_heating) then
             h2form1 = h2_formation(h2fraci,t_found,dust_tempi,xnH2,metallicity) * 4.48*eV/rhoi/unit_pressure*utime
             !
             !--Potentially add heating from UV destruction and pumping of H_2
             !     (see Bate 2015) for effect of this.
             !
             !     &                          + h2_destruction(i,xnH2,.FALSE.) *
             !     &               (0.4 + 2.0/(1.0+criticaln(i,t_found)/(2.0*xnH2)))*
             !     &                          eleccharge/rhoi/uergcc*utime
             h2form2 = h2_formation(h2fraci,t_plus,dust_tempi,xnH2,metallicity) * 4.48*eV/rhoi/unit_pressure*utime
             !
             !--Potentially add heating from UV destruction and pumping of H_2
             !     (see Bate 2015) for effect of this.
             !
             !     &                          + h2_destruction(i,xnH2,.FALSE.) *
             !     &               (0.4 + 2.0/(1.0+criticaln(i,t_plus)/(2.0*xnH2)))*
             !     &                          eleccharge/rhoi/uergcc*utime
             dh2form_dt = (h2form2-h2form1)/tdiff
          else
             h2form1 = 0.
             h2form2 = 0.
             dh2form_dt = 0.
          endif
          !
          !--Photoelectric heating
          !
          photoelectric1 = photoelectric_heating(i,t_found,xnH2,metallicity)/rhoi/unit_pressure*utime
          photoelectric2 = photoelectric_heating(i,t_plus,xnH2,metallicity)/rhoi/unit_pressure*utime
          dphoto_du = (photoelectric2-photoelectric1)/tdiff
          !
          !--Now perform Newton-Raphson iteration
          !
          func = u_found + dti*gas_dust_val*(u_found/cv1 - dust_tempi) - orig_ugasi - dti*pres_numerator &  ! orig_ugasi = origEU(2,i)
                           - dti*e_planetesimali - dti*cosmic_ray + dti*cooling_line1 - dti*photoelectric1 &
                           + dti*dust_term - dti*h2form1

          derivative = (u_plus - u_found)/tdiff + dti*gas_dust_val + dti*dline_du - dti*dphoto_du - dti*dh2form_dt
          !
          !--If failed, do it again, but this time print out diagnostics
          !
          if (itry == 2) then
!$omp critical(nrprint)
             print*,'N-R',iterationloop,i,t_orig,rhoi
             print*,'t_f,t_p',t_found,t_plus
             print*,'u_f,u_l,u_p',u_found,u_last,u_plus,cv1
             print*,dti*gas_dust_val*(u_found/cv1 - dust_tempi),orig_ugasi,-dti*cosmic_ray, & ! orig_ugasi = origEU(2,i)
                    dti*cooling_line1,-dti*photoelectric1,dust_tempi,-dti*pres_numerator,dust_term,func,derivative
             print*,u_plus,u_found,tdiff,(u_plus - u_found)/tdiff,dti*gas_dust_val,dti*dline_du, &
                    -dti*dphoto_du,photoelectric2,photoelectric1
!$omp end critical(nrprint)
          endif
          !
          !--Limit the change to be no larger than a certain fraction each iteration
          !
          func_old = func
          if (func/derivative / t_found > 0.3) then
             func = 0.3*derivative * t_found
          elseif (func/derivative / t_found < -0.3) then
             func = -0.3*derivative * t_found
          endif
          t_found = t_found - func/derivative

          if (t_found < 1.) t_found = 1.

          u_found = get_u_from_rhoT(rhoi,t_found,ieos)
          !
          !--Test for success
          !
          if (abs((t_found - t_last)/t_orig) < 1.e-3) then
             U1i = get_u_from_rhoT(rhoi,t_found,ieos)
             cvi = cv1
             photoelectric = photoelectric1
             cooling_line = cooling_line1
             h2form = h2form1
             skip_quartic = .true.
             return
          endif
       enddo
!$omp critical(quart)
       print *,"N-R failed for ieqtype=3, try again"
       !  print *,"ngs ",u4term,u1term,u0term,betaval,chival,gammaval
       !  print *,"    ",rhoi,dti
       !  print *,"    ",diffusion_denominator,diffusion_numerator
       !  print *,"    ",pres_denominator,pres_numerator,uradconst
       !  print *,"    ",radpresdenom,eradi,egasi,cvi
       print *,"    ",orig_eradi,orig_ugasi
       !  print *,"    ",dti*gas_dust_val*dust_tempi
       !  print *,"    ",dti*dust_heating,dti*gas_dust_val*dust_tempi
       !  print *,"    ",dti*cosmic_ray,dti*cooling_line
       !  print *,"    ",dti*photoelectric,dti*dust_term
       !  print *,"    ",ieqtype,u_found,u_last,u_plus,cv1
       !  print *,"    ",cooling_line1,cooling_line2,dline_du,func
       !  print *,"    ",derivative,h2form1,h2form2,dh2form_dt
!$omp end critical(quart)
       if (itry == 2) ierr = 1  ! if itry=1, try again, otherwise quit
    enddo try_loop
 endif  !-turn on/off n-r solve

end subroutine set_heating_cooling_low_rhoT


subroutine set_heating_cooling(i,ugasi,cvi,rhoi,mui,heatingISRi,metallicity,ieqtype, &
           dust_tempi,gas_dust_val,dustgammaval,gas_dust_cooling, &
           cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term)
 use units, only:unit_ergg,utime,unit_pressure,unit_density
 use physcon, only:eV,mass_proton_cgs
 integer, intent(in)  :: i
 real, intent(in)     :: ugasi,rhoi,cvi,heatingISRi,metallicity,mui
 integer, intent(out) :: ieqtype
 real, intent(inout)  :: dust_tempi
 real, intent(out)    :: gas_dust_val,dustgammaval,gas_dust_cooling,cosmic_ray,cooling_line,&
                         photoelectric,h2form,dust_heating,dust_term
 real                 :: gas_temp,xnH2,h2fraci

 ieqtype = 2
 gas_temp = ugasi/cvi
 xnH2 = rhoi*unit_density/(mui*mass_proton_cgs)
 !--Use existing values of betaval, gammaval, chival, tfour
 !     because the values from getkappa include both
 !     the gas and dust opacities already.
 gas_dust_val = 0.
 dustgammaval = 0.
 gas_dust_cooling = 0.
 dust_heating = 0.
 dust_term = 0.

 cosmic_ray = cosmic_ray_heating(xnH2) / rhoi / unit_pressure*utime

 !--This adds the dust heating term into the equations.  However, because the
 !     above assumption is that T_g=T_d, this makes the low-density gas
 !     as hot as the dust whereas in fact it should be a lot cooler.
 cosmic_ray = cosmic_ray + heatingISRi/unit_ergg*utime
 cooling_line = cooling_line_rate(i,gas_temp,xnH2,metallicity) / rhoi / unit_pressure*utime
 !
 !--H_2 formation heating
 !
 if (H2formation_heating) then
    h2form = h2_formation(h2fraci,gas_temp,dust_tempi,xnH2,metallicity)*4.48*eV / rhoi / unit_pressure*utime
    !
    !--Potentially add heating from UV destruction and pumping of H_2
    !     (see Bate 2015) for effect of this.
    !c
    !c     &                 + h2_destruction(i,xnH2,.FALSE.) *
    !c     &              (0.4 + 2.0/(1.0+criticaln(i,gas_temp)/(2.0*xnH2)))*
    !c     &                 eleccharge/rhoi/uergcc*utime
 else
    h2form = 0.
 endif

 photoelectric = photoelectric_heating(i,gas_temp,xnH2,metallicity) / rhoi / unit_pressure*utime

end subroutine set_heating_cooling


subroutine turn_heating_cooling_off(ieqtype,dust_tempi,gas_dust_val,dustgammaval,gas_dust_cooling,&
                                    cosmic_ray,cooling_line,photoelectric,h2form,dust_heating,dust_term)
 integer, intent(out) :: ieqtype
 real, intent(out)    :: dust_tempi,gas_dust_val,dustgammaval,gas_dust_cooling,cosmic_ray,cooling_line,&
                         photoelectric,h2form,dust_heating,dust_term

 ieqtype = 0

 dust_tempi = 0.
 gas_dust_val = 0.
 dustgammaval = 0.
 gas_dust_cooling = 0.

 cosmic_ray = 0.
 cooling_line = 0.
 photoelectric = 0.
 h2form = 0.

 dust_heating = 0.
 dust_term = 0.

end subroutine turn_heating_cooling_off


subroutine store_radiation_results(ncompactlocal,npart,ivar,EU0,rad,vxyzu)
 integer, intent(in) :: ncompactlocal,npart,ivar(:,:)
 real, intent(in)    :: EU0(2,npart)
 real, intent(out)   :: rad(:,:),vxyzu(:,:)
 integer :: i,n

 !$omp parallel do default(none) &
 !$omp shared(ncompactlocal,ivar,vxyzu,EU0,rad) &
 !$omp private(n,i)
 do n = 1,ncompactlocal
    i = ivar(3,n)
    rad(iradxi,i) = EU0(1,i)
    vxyzu(4,i) = EU0(2,i)
 enddo
!$omp end parallel do

end subroutine store_radiation_results


real function dust_temperature(xi,u,rho,dust_kappa,dust_cooling,heatingISR,dust_gas)
 real, intent(in)    :: xi,u,rho
 real, intent(out)   :: dust_kappa,dust_cooling,heatingISR,dust_gas

 dust_temperature = 0.
 dust_cooling = 0.
 heatingISR = 0.
 dust_gas = 0.
 dust_kappa = 0.

end function dust_temperature

!---------------------------------------------------------
!+
!  Following Goldsmith et al. (2001), we set the cosmic ray heating rate to be
!  Gamma_cr = 1.0E-27 * n_H2  in erg/cm^3/s.
!+
!---------------------------------------------------------
real function cosmic_ray_heating(xnH2)
 real, intent(in) :: xnH2
 if (use_cosmic_ray_heating) then
    cosmic_ray_heating = 1.e-27*xnH2
 else
    cosmic_ray_heating = 0.
 endif
end function cosmic_ray_heating

real function cooling_line_rate(ipart,gas_temp,xnH2,metallicity)
 integer, intent(in) :: ipart
 real, intent(in)    :: gas_temp,xnH2,metallicity

 cooling_line_rate = 0.

end function cooling_line_rate

!---------------------------------------------------------
!+
!  Based on Glover et al. 2010 (rate 165 in appendix)
!+
!---------------------------------------------------------
real function h2_formation(h2fraci,gas_temp,dust_temp,xnH2,metallicity)
 real, intent(in) :: h2fraci,gas_temp,dust_temp,xnH2,metallicity
 real             :: grain_formation_rate,G0,fA,fB,xnH

 !--nH is twice nH2, where nH is actually the number density of protons from
 !     hydrogen (i.e. it does not depend on the fractions of H and H_2 )
 xnH = 2.*xnH2
 G0 = 1.
 !
 !--Partial rate from Glover et al. 2010 (rate 165 in appendix)
 !
 if (dust_temp > 20.) then
    fA = 1./(1. + 1.e4*exp(-600./dust_temp))
 else
    fA = 1.
 endif
 fB = 1/(1 + 0.04*sqrt(gas_temp + dust_temp) + 0.002*gas_temp + 8e-6*gas_temp**2)
 grain_formation_rate = 3.E-18*sqrt(gas_temp) * fA * fB * ((1. - 2.*h2fraci)*xnH)**2 * metallicity
 h2_formation = grain_formation_rate
 !-Energy released is 4.48eV times the formation rate (eV in erg),
 ! with H_2 formation *heating* the gas

end function h2_formation

!---------------------------------------------------------
!+
!  Following Young et al. (2004) we set the photoelectric heating rate
!  of the gas by electrons liberated from grains as
!
!     Gamma_pe = 1E-24 * 0.05 * G(r) * n_H
!
!  and we take n_H = 2*n_H2 + 4*n_He = 2.66 nH2
!
!  The units are erg/cm^3/s and the function G(r) is the attenuation
!  factor for high energy photons (=1 for no extinction).
!+
!---------------------------------------------------------
real function photoelectric_heating(ipart,gas_temp,xnH2,metallicity)
 integer, intent(in) :: ipart
 real, intent(in)    :: xnH2,gas_temp,metallicity
!  real                :: xne,phiPAH,xnH,G0,ep,electron_fraction

 photoelectric_heating = 0. ! Mike: Set to zero for now
!  xnH = 2.heatingISR(2,ipart)  ! Mike: Where does heatingISR come from?
! ! G0 = 1.

!  ep = 0.049/(1.+0.004*(G0*sqrt(gas_temp)/xne/phiPAH)**0.73) +&
!      0.037*(gas_temp/1.e+4)**0.7/&
!      (1. + 2.e-4*(G0*sqrt(gas_temp)/xne/phiPAH))
! ! ep = 0.05

!  if (use_photoelectric_heating) then
!     photoelectric_heating = 1.3e-24*ep*g0*xnh*metallicity
!  else
!     photoelectric_heating = 0.
!  endif*xnH2
!  xne = electron_fraction(xnH,phiPAH)
!  G0 = hea

end function photoelectric_heating

!---------------------------------------------------------
!+
!  Need to approximate electron density (approximation from Fig 10 of
!  Wolfire et al. 2003).  Note that this needs n(H) rather than
!  n(H2).
!+
!---------------------------------------------------------
real function electron_fraction(xnH,phiPAH_opt)
 real, intent(in)           :: xnH
 real, intent(in), optional :: phiPAH_opt
 real                       :: xne_over_nH,phiPAH

 if (present(phiPAH_opt)) then
    phiPAH = phiPAH_opt
 else
    !--Need to set phiPAH (Wolfire et al. use 0.5, but I use 0.55 to get
    !     closer to 10^4 K at very low ISM densities (<1 cm^-3)
    phiPAH = 0.55
 endif

 xne_over_nH = 0.008/xnH
 if (xne_over_nH < 1.e-4) xne_over_nH = 1.e-4
 if (xne_over_nH > 1.) xne_over_nH = 1.
 electron_fraction = xne_over_nH*xnH

end function electron_fraction

!---------------------------------------------------------
!+
!  Returns the first bit of the gas-dust collisional heating/cooling
!  term used by Keto & Field (2005), which is
!  Lambda_gd = 1.0E-33 * n_H2**2 * sqrt(T_gas) * (T_gas - T_dust)
!  or the coupling term used by Glover & Clark (2012), which is
!  Lambda_gd = 1.5E-32 * n_H2**2 * sqrt(T_gas) * (T_gas - T_dust) *
!              (1.0 + 0.8*exp(-75.0/T_gas)
!  in erg/cm^3/s.  This function also includes a metallicity term
!  which essentially assumes that the dust number density depends
!  linearly on the metallicity.
!
!  The value is calculated in cgs units and excludes the
!  (T_gas-T_dust) term because of the need for derivatives w.r.t.
!  T_gas in solving the system of equations.
!+
!---------------------------------------------------------
real function gas_dust_collisional_term(xnH2,metallicity,gas_temp)
 real, intent(in) :: xnH2,metallicity,gas_temp
 real             :: expfac

 if (gas_dust_collisional_term_type == 1) then
    gas_dust_collisional_term = 1.e-33*(xnH2**2)*sqrt(gas_temp)*metallicity
 elseif (gas_dust_collisional_term_type == 2) then
    if (gas_temp > 4.) then
       expfac = 1.-0.8*exp(-75./gas_temp)
    else
       expfac = 1.
    endif
    gas_dust_collisional_term = 1.5e-32*(xnH2**2)*sqrt(gas_temp)*expfac*metallicity
 else
    gas_dust_collisional_term = 0.
 endif

end function gas_dust_collisional_term

!---------------------------------------------------------
!+
!  solve quartic (eq 22, Whitehouse et al. 2005)
!+
!---------------------------------------------------------
subroutine solve_quartic(u1term,u0term,uold,soln,moresweep,ierr)
 use quartic, only:quarticsolve
 real, intent(in) :: u1term,u0term,uold
 integer, intent(out)   :: ierr
 logical, intent(inout) :: moresweep
 real,    intent(out)   :: soln
 real :: a(0:3)

 ! Between eq 22 & 23
 ! a4 = 1. but this is already assumed in the solver
 a(3) = 0.
 a(2) = 0.
 a(1) = u1term
 a(0) = u0term

 call quarticsolve(a,uold,soln,moresweep,ierr)

end subroutine solve_quartic

end module radiation_implicit
