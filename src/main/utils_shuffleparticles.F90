!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_shuffleparticles
!
! Using only particle positions, this will shuffle those particles so
! that their density profile will match an input density profile.
! Current permitted profiles:
!   A uniform density
!   sphere-in-box (i.e. uniform density sphere in a uniform density medium)
!   BE sphere with a uniform background from an input array
!   an input array [UNTESTED]
!   using a different particle set
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, boundary, densityforce, dim, io, kdtree,
!   kernel, linklist, mpidomain, part
!
 implicit none
 public :: shuffleparticles

 ! identifiers for exact profile
 integer, parameter :: iuniform     = 1
 integer, parameter :: isphereinbox = 2
 integer, parameter :: iarray       = 3
 integer, parameter :: ireference   = 4
 ! additional parameters
 integer, parameter :: maxcellcache = 50000
 integer            :: ncall        =     1  ! to track the number of calls to this subroutine for debugging
 logical            :: firstcall    = .true. ! logical to control debugging print statements
 private

contains
!-------------------------------------------------------------------
!
! Shuffle the particles to match an input profile
! Options and their specific inputs
!  Uniform density:
!     duniform: density of the uniform medium
!  Sphere-in-box:
!     rsphere: radius of the sphere
!     dsphere: density of the sphere
!     dmedium: density of the medium
!  Array with uniform background (designed for a BE sphere in a background medium)
!     ntab: array lengths
!     rtab: array radial positions
!     dtab: array with density positions
!     dmedium: density of the medium
!     dcontrast: density contrast between medium and edge of sphere
!  Array without a background
!     ntab: array lengths
!     rtab: array radial positions
!     dtab: array with density positions
!  Particle distribution:
!     xyzh_ref:  array containing the position of the reference particles
!     pmass_ref: mass of the reference particle
!     n_ref:     number of reference particles
!
!-------------------------------------------------------------------
subroutine shuffleparticles(iprint,npart,xyzh,pmass,duniform,rsphere,dsphere,dmedium,ntab,rtab,dtab,dcontrast, &
                            xyzh_ref,pmass_ref,n_ref,is_setup,prefix)
 use io,           only:id,master,fatal
 use dim,          only:maxneigh,maxp_hard
 use part,         only:vxyzu,divcurlv,divcurlB,Bevol,fxyzu,fext,alphaind,iphase,igas
 use part,         only:gradh,rad,radprop,dvdx,rhoh,hrho,apr_level
 use densityforce, only:densityiterate
 use linklist,     only:ncells,ifirstincell,set_linklist,get_neighbour_list,allocate_linklist,listneigh
 use kernel,       only:cnormk,wkern,grkern,radkern2
#ifdef PERIODIC
 use boundary,     only:dxbound,dybound,dzbound,cross_boundary
 use mpidomain,    only:isperiodic
#endif
 use kdtree,       only:inodeparts,inoderange
 use allocutils,   only:allocate_array
 integer,           intent(in)    :: iprint
 integer,           intent(inout) :: npart
 real,              intent(in)    :: pmass
 real,              intent(inout) :: xyzh(:,:)
 integer, optional, intent(in)    :: ntab,n_ref
 real,    optional, intent(in)    :: duniform,rsphere,dsphere,dmedium,dcontrast,pmass_ref
 real,    optional, intent(in)    :: rtab(:),dtab(:)
 real,    optional, intent(inout) :: xyzh_ref(:,:)
 logical, optional, intent(in)    :: is_setup
 character(len=*) , optional, intent(in)    :: prefix
 integer      :: i,j,k,jm1,p,ip,icell,ineigh,idebug,ishift,nshiftmax,iprofile,nparterr,nneigh,ncross,nlink,n_part
 real         :: stressmax,redge,dedge,dmed
 real         :: max_shift2,max_shift_thresh,max_shift_thresh2,link_shift,link_shift_thresh,radkern12
 real         :: xi,yi,zi,hi,hi12,hi14,radi,radinew,coefi,rhoi,rhoi1,rij2,qi2,qj2,denom,rhoe,drhoe,err
 real         :: xj,yj,zj,hj,hj1,termi
 real         :: maggradi,maggrade,magshift,rinner,router
 real         :: dx_shift(3,npart),rij(3),runi(3),grrhoonrhoe(3),grrhoonrhoi(3)
 real         :: errmin(3,2),errmax(3,2),errave(3,2),stddev(2,2),rthree(3),totalshift(3,npart)
 real         :: kernsum,wkerni
#ifdef SPLITTING
 real         :: gradhi,gradhj,grrhoonrhoe_ref(3,maxp_hard)
 real,save    :: xyzcache_ref(maxcellcache,4)
#endif
 real,save    :: xyzcache(maxcellcache,4)
 logical      :: shuffle,at_interface,use_ref_h,call_linklist,is_ref
 character(len=128) :: prefix0,fmt1,fmt2,fmt3
 !$omp threadprivate(xyzcache)
#ifdef SPLITTING
 !$omp threadprivate(xyzcache_ref)
#endif

 !--Initialise free parameters
 idebug            =    1  ! 0 = off; 1=errors; 2=initial & final distribution + (1); 3=print every step + (2)
 nshiftmax         =  200 !600 ! maximum number of shuffles/iterations
 max_shift_thresh  = 0. !4.d-3 ! will stop shuffling once (maximum shift)/h is less than this value
 link_shift_thresh = 0.01  ! will recalculate the link list when the cumulative maximum relative shift surpasses this limit (=0 will call every loop)
 !--Initialise remaining parameters
 rthree            = 0.
 use_ref_h         = .true. ! to prevent compiler warnings
 call_linklist     = .true.
 max_shift_thresh2 = max_shift_thresh*max_shift_thresh
 n_part            = npart
 is_ref            = .false.
 radkern12         = 1.0/radkern2
 ncross            = 0

 !--Determine what has been passed in
 if (idebug > 0 .and. id==master) then
    write(fmt1,'(a)') '(1x,a)'
    write(fmt2,'(a)') '(1x,a,I8)'
    write(fmt3,'(a)') '(1x,a,Es18.6)'
    if (present(ntab))      write(iprint,fmt2) 'Shuffling: optional variable is present: ntab: ', ntab
    if (present(duniform))  write(iprint,fmt1) 'Shuffling: optional variable is present: duniform'
    if (present(rsphere))   write(iprint,fmt1) 'Shuffling: optional variable is present: rsphere'
    if (present(dsphere))   write(iprint,fmt1) 'Shuffling: optional variable is present: dsphere'
    if (present(dmedium))   write(iprint,fmt1) 'Shuffling: optional variable is present: dmedium'
    if (present(dcontrast)) write(iprint,fmt1) 'Shuffling: optional variable is present: dcontrast'
    if (present(rtab))      write(iprint,fmt2) 'Shuffling: optional variable is present: rtab: ', size(rtab)
    if (present(dtab))      write(iprint,fmt2) 'Shuffling: optional variable is present: dtab: ', size(rtab)
    if (present(xyzh_ref))  write(iprint,fmt1) 'Shuffling: optional variable is present: xyzh_ref'
    if (present(pmass_ref)) write(iprint,fmt3) 'Shuffling: optional variable is present: pmass_ref: ',pmass_ref
    if (present(n_ref))     write(iprint,fmt2) 'Shuffling: optional variable is present: n_ref: ',n_ref
    if (present(is_setup))  write(iprint,fmt1) 'Shuffling: optional variable is present: is_setup'
    write(iprint,fmt2) 'Shuffling: non-optional variable: npart = ',npart
    write(iprint,fmt3) 'Shuffling: non-optional variable: pmass = ',pmass
 endif

 !--Determine the exact profile
 iprofile = 0
 redge    = 0.
 rinner   = 0.
 router   = 0.
 dmed     = 0.
 dedge    = 0.
 if (present(duniform)) then
    if (id==master) write(iprint,'(1x,a)') 'Shuffling: Shuffling to a uniform density medium'
    iprofile = iuniform
    dmed     = duniform
 elseif (present(rsphere) .and. present(dsphere) .and. present(dmedium)) then
    if (id==master) write(iprint,'(1x,a)') 'Shuffling: Shuffling to a uniform sphere-in-box profile'
    iprofile = isphereinbox
    dedge    = dsphere
    dmed     = dmedium
    redge    = rsphere
 elseif (present(ntab) .and. present(rtab) .and. present(dtab) .and. present(dmedium) .and. present(dcontrast)) then
    if (id==master) then
       write(iprint,'(1x,2a)') 'Shuffling: Shuffling to non-uniform sphere-in-box profile',&
                               ' (e.g., a BE sphere in a background medium)'
    endif
    iprofile = iarray
    dedge    = dmedium*dcontrast
    dmed     = dmedium
    redge    = rtab(ntab)
 elseif (present(ntab) .and. present(rtab) .and. present(dtab)) then
    if (id==master) then
       write(iprint,*) 'Shuffling to a non-uniform sphere'
       write(iprint,*) 'This is untested without a background, so aborting until properly tested'
    endif
    return
    iprofile = iarray
    redge    = rtab(ntab)
 elseif (present(xyzh_ref) .and. present(pmass_ref) .and. present(n_ref)) then
    if (id==master) write(iprint,'(1x,a)') 'Shuffling: Shuffling to a reference array'
    if (id==master) write(iprint,'(1x,a)') 'Shuffling: This option is currently not available for public use'
    if (id==master) write(iprint,'(1x,a)') 'Shuffling: Key features have been removed.  Aborting.'
    return
    iprofile = ireference
    if (pmass > pmass_ref) then
       use_ref_h = .false.
    else
       use_ref_h = .true.
    endif
#ifdef SPLITTING
    !--Ensure maxp_hard is large enough to linklist the primary & reference simultansously
    n_part = npart + n_ref
    if (n_part > maxp_hard) call fatal('shuffling','npart + n_ref > maxp_hard',var='n_part',ival=n_part)

#else
    call fatal('shuffling','SPLITTING==yes is required for iprofile == ireference')
#endif
    !--Append the following arrays with the reference particles for linking
    xyzh(:,npart+1:n_part) = xyzh_ref(:,1:n_ref)
    iphase(npart+1:n_part) = igas                ! don't really care what this value is, just can't be zero
 else
    if (id==master) write(iprint,'(1x,a)') 'Shuffling: Density profile undefined.  Aborting.'
    return
 endif
 rinner = redge - 4.*hrho(dedge,pmass)
 router = redge + 4.*hrho(dedge,pmass)

 !--Open debugging files and print the initial particle placements
 if (idebug > 0 .and. id==master) then
    if (present(prefix)) then
       prefix0 = trim(prefix)//'_'
    else
       prefix0 = ''
    endif
    if (firstcall) then
       if (idebug > 1) then
          open(unit=333,file=trim(prefix0)//'shuffling_partposition.dat')
          write(333,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
             1,'ipart',  &
             2,'x',      &
             3,'y',      &
             4,'z',      &
             5,'h',      &
             6,'rho',    &
             7,'Sum dx', &
             8,'Sum dy', &
             9,'Sum dz', &
            10,'dx',     &
            11,'dy',     &
            12,'dz',     &
            13,'ncall'
       endif
       open(unit=334,file=trim(prefix0)//'shuffling_error_all.dat')
       if (iprofile/=ireference) then
          open(unit=335,file=trim(prefix0)//'shuffling_error_notinterface.dat')
          p = 335
       else
          p = 334
       endif
       do j = 334,p
          write(j,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
             1,'ishift',&
             2,'min rho err', &
             3,'ave rho err', &
             4,'max rho err', &
             5,'mn dd/d err', &
             6,'av dd/d err', &
             7,'mx dd/d err', &
             8,'min shift/h', &
             9,'ave shift/h', &
            10,'max shift/h', &   ! This is the criteria used to exit the shuffling loop
            11,'stddev_rho',  &
            12,'stddev_grho', &
            13,'ncall'
       enddo
       firstcall = .false.
    else
       if (idebug > 1) then
          open(unit=333,file=trim(prefix0)//'shuffling_partposition.dat',position='append')
          write(333,'(a)') ' '
       endif
       open(unit=334,file=trim(prefix0)//'shuffling_error_all.dat',position='append')
       write(334,'(a)') ' '
       if (iprofile==isphereinbox .or. iprofile==iarray) then
          open(unit=335,file=trim(prefix0)//'shuffling_error_notinterface.dat',position='append')
          write(335,'(a)') ' '
       endif
    endif
 endif
 if (idebug > 1 .and. id==master) then
    do i = 1,npart
       write(333,'(I18,1x,11(es18.10,1x),I18)') i,xyzh(1:4,i),rhoh(xyzh(4,i),pmass),rthree,rthree,ncall
    enddo
    write(333,'(a)') ' '
 endif

 !--initialise memory for linklist
 if (present(is_setup)) then
    if (is_setup) call allocate_linklist()
 endif

 !--Shuffle particles
 shuffle    = .true.
 ishift     = 0
 totalshift = 0.
 nlink      = 0
 link_shift = 0.

 do while (shuffle .and. ishift < nshiftmax)

    ! update densities
    if (call_linklist .or. iprofile==ireference) then
       call set_linklist(npart,npart,xyzh,vxyzu)
       nlink      = nlink + 1
       link_shift = 0.
    endif
    call densityiterate(2,npart,npart,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                               fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level)
    if (iprofile==ireference) then
       call set_linklist(n_part,n_part,xyzh,vxyzu)
    endif

    ! initialise variables for this loop
    dx_shift   = 0.
    errmin     = huge(errmin)
    errmax     = 0.
    errave     = 0.
    stddev     = 0.
    nparterr   = 0
    max_shift2 = 0.
    radi       = 0.
    radinew    = 0.

    ! determine how much to shift by
!$omp parallel default (none) &
!$omp shared(xyzh,npart,pmass,dx_shift,gradh,max_shift_thresh2,redge,iprofile,dedge,dmed,idebug) &
!$omp shared(use_ref_h,n_ref,xyzh_ref,totalshift,pmass_ref,radkern12,ishift) &
!$omp shared(rtab,dtab,ntab,rinner,router,inodeparts,inoderange,ifirstincell,ncells,ncall) &
#ifdef PERIODIC
!$omp shared(dxbound,dybound,dzbound) &
#endif
#ifdef SPLITTING
!$omp shared(grrhoonrhoe_ref) &
#endif
!$omp private(i,j,k,ip,ineigh,icell,jm1,xi,yi,zi,hi,hi12,rij,rij2,runi,qi2,qj2,hi14,rhoi,rhoi1,coefi) &
!$omp private(xj,yj,zj,hj,hj1,termi) &
!$omp private(grrhoonrhoe,grrhoonrhoi,rhoe,drhoe,denom,nneigh) &
!$omp private(err,maggradi,maggrade,magshift,at_interface) &
!$omp private(kernsum,wkerni) &
!$omp firstprivate(is_ref,radi,radinew) &
!$omp reduction(min: errmin) &
!$omp reduction(max: errmax,max_shift2) &
!$omp reduction(+:   errave,stddev,nparterr)
!$omp do schedule(runtime)
    over_cells: do icell=1,int(ncells)
       k = ifirstincell(icell)

       ! Skip empty cells AND inactive cells
       if (k <= 0) cycle over_cells

       ! Get the neighbour list and fill the cell cache
       call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)

       ! Loop over particles in the cell
       over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
          i  = inodeparts(ip)
          if (i > npart) cycle over_parts ! skip reference particles; guaranteed to be false if not splitting
          xi = xyzh(1,i)
          yi = xyzh(2,i)
          zi = xyzh(3,i)
          hi = xyzh(4,i)
          hi12   = 1.0/(hi*hi)
          hi14   = hi12*hi12
          rhoi   = rhoh(hi,pmass)
          rhoi1  = 1.0/rhoi
          termi  = cnormk*gradh(1,i)*hi14*rhoi1
          coefi  = 0.25*hi*hi
          grrhoonrhoi = 0.0
          grrhoonrhoe = 0.0
          drhoe       = 0.0
          rhoe        = 0.0
          kernsum     = 0.0

          ! calculate (grad rho) / rho for particle i
          over_neighbours: do ineigh = 1,nneigh
             j = abs(listneigh(ineigh))
             if (i==j) cycle over_neighbours
             if (iprofile==ireference) then
                ! determine if primary or reference particle
                if (j > npart) then
                   is_ref = .true.
                else
                   is_ref = .false.
                endif
             endif

             if (ineigh <= maxcellcache) then
                ! positions from cache are already mod boundary
                xj  = xyzcache(ineigh,1)
                yj  = xyzcache(ineigh,2)
                zj  = xyzcache(ineigh,3)
                hj1 = xyzcache(ineigh,4)
                hj  = 1./hj1
             else
                xj  = xyzh(1,j)
                yj  = xyzh(2,j)
                zj  = xyzh(3,j)
                hj  = xyzh(4,j)
                hj1 = 1./hj
             endif

             rij(1) = xj - xi
             rij(2) = yj - yi
             rij(3) = zj - zi
#ifdef PERIODIC
             if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
             if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
             if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
#endif
             rij2 = dot_product(rij,rij)
             if (rij2 < epsilon(rij2)) then
                if (idebug > 0 .and. .not.is_ref) print*, 'pairing may be occurring between ', i, j
                cycle
             endif

             if (is_ref) then
                ! This is for splitting only!
#ifdef SPLITTING

                ! Determine the kernel-weighted average grad rho / rho for the reference distribution
                if (use_ref_h) then  ! Always use the larger smoothing length for pseudo-neighbour finding.
                   qi2 = rij2*hj1**2
                else
                   qi2 = rij2*hi12
                endif

#endif
             else
                ! Calculate (grad rho)/rho of the primary particle type
                qi2 = rij2*hi12
                qj2 = rij2*hj1**2
                if (qi2 < radkern2 .or. qj2 < radkern2) then
                   runi = rij/sqrt(rij2)
                   if (qi2 < radkern2) then
                      grrhoonrhoi = grrhoonrhoi - runi*termi*grkern(qi2,sqrt(qi2))
                   endif
                   if (qj2 < radkern2) then
                      denom = hj**4 * rhoh(hj,pmass)
                      grrhoonrhoi = grrhoonrhoi - runi*cnormk*gradh(1,j)*grkern(qj2,sqrt(qj2))/denom
                   endif
                endif
             endif
          enddo over_neighbours

          ! calculate the exact (grad rho) / rho
          radi = 0.
          if (iprofile==ireference) then
             ! Finalise calculation for splitting
             if (kernsum > 0.) then
                grrhoonrhoe = grrhoonrhoe/kernsum
             else
                print*, 'Shuffling: WARNING! No neighbours when determining reference value!'
             endif
          elseif (iprofile==iuniform) then
             rhoe = dmed
          else
             radi = sqrt(xi*xi + yi*yi + zi*zi)
             if (radi < redge) then
                if (iprofile==isphereinbox) then
                   rhoe = dedge
                elseif (iprofile==iarray) then
                   j = 1
                   do while (radi > rtab(j) .and. j < ntab)
                      j = j + 1
                   enddo
                   if (j==1) then
                      jm1 = j + 1
                   else
                      jm1 = j - 1
                   endif
                   drhoe = (dtab(j) - dtab(jm1))/(rtab(j) - rtab(jm1))
                   rhoe  = dtab(jm1) +  drhoe*(radi - rtab(jm1))
                   grrhoonrhoe(1) = drhoe*xi/(rhoe*radi)
                   grrhoonrhoe(2) = drhoe*yi/(rhoe*radi)
                   grrhoonrhoe(3) = drhoe*zi/(rhoe*radi)
                endif
             else
                rhoe  = dmed
             endif
          endif

          ! shift the particles
          grrhoonrhoi  = coefi*pmass*grrhoonrhoi   ! multiply in the correct coefficients for the i'th particle
          grrhoonrhoe  = coefi      *grrhoonrhoe   ! multiply in the correct coefficients for the exact value
          dx_shift(:,i) = grrhoonrhoi - grrhoonrhoe

          ! limit the particle shift to 0.5h
          magshift   = dot_product(dx_shift(1:3,i),dx_shift(1:3,i))
          if (magshift > 0.25*hi*hi) then
             dx_shift(:,i) = 0.5*hi*dx_shift(:,i)/sqrt(magshift)
             magshift      = dot_product(dx_shift(1:3,i),dx_shift(1:3,i))
          endif

          ! for spheres-in-boxes, prevent particles from moving from one medium to the other
          if (iprofile==isphereinbox .or. iprofile==iarray) then
             radinew = sqrt( (xi-dx_shift(1,i))**2 + (yi-dx_shift(2,i))**2 + (zi-dx_shift(3,i))**2)
             if (radi < redge .and. radinew > redge) then
                ! prevent particle from moving outside of the sphere
                dx_shift(:,i) = dx_shift(:,i)*(redge-radi)/sqrt(magshift) - dx_shift(:,i)*(radinew-redge)/sqrt(magshift)
             elseif (radi > redge .and. radinew < redge) then
                ! prevent particle from moving inside of the sphere
                dx_shift(:,i) = dx_shift(:,i)*(radinew-redge)/sqrt(magshift) - dx_shift(:,i)*(redge-radi)/sqrt(magshift)
             endif
             magshift = dot_product(dx_shift(1:3,i),dx_shift(1:3,i))
          endif
          max_shift2 = max(max_shift2,magshift*hi12) ! metric to track for exiting loop

          ! debugging statements
          if (idebug > 0) then
             if (idebug > 1) totalshift(:,i) = totalshift(:,i) - dx_shift(:,i)
             if (radi < rinner .or. radi > router .or. iprofile==ireference) then
                at_interface = .false.
                nparterr     = nparterr + 1
             else
                at_interface = .true.
             endif
             ! calculate error in density
             err       = abs(1.0-rhoe*rhoi1)
             errmin(1,1) = min(errmin(1,1),err)
             errmax(1,1) = max(errmax(1,1),err)
             errave(1,1) = errave(1,1) + err
             if (.not. at_interface) then
                errmin(1,2) = min(errmin(1,2),err)
                errmax(1,2) = max(errmax(1,2),err)
                errave(1,2) = errave(1,2) + err
             endif
             ! calculate error in grad rho / rho
             maggrade  = dot_product(grrhoonrhoe,grrhoonrhoe)
             maggradi  = dot_product(grrhoonrhoi,grrhoonrhoi)
             if (maggrade > 0.) then
                err    = abs(1.0-sqrt(maggradi/maggrade))
             else
                err    = sqrt(maggradi)
             endif
             errmin(2,1) = min(errmin(2,1),err)
             errmax(2,1) = max(errmax(2,1),err)
             errave(2,1) = errave(2,1) + err
             if (.not. at_interface) then
                errmin(2,2) = min(errmin(2,2),err)
                errmax(2,2) = max(errmax(2,2),err)
                errave(2,2) = errave(2,2) + err
             endif
             ! calculate the magnitude of the shift
             err         = sqrt(magshift)/hi
             errmin(3,1) = min(errmin(3,1),err)
             errmax(3,1) = max(errmax(3,1),err)
             errave(3,1) = errave(3,1) + err
             if (.not. at_interface) then
                errmin(3,2) = min(errmin(3,2),err)
                errmax(3,2) = max(errmax(3,2),err)
                errave(3,2) = errave(3,2) + err
             endif
             ! calculate the standard deviation compared to the exact values
             stddev(1,1) = stddev(1,1) + (rhoe - rhoi)**2
             stddev(2,1) = stddev(2,1) + (maggrade - maggradi)**2
             if (.not. at_interface) then
                stddev(1,2) = stddev(1,2) + (rhoe - rhoi)**2
                stddev(2,2) = stddev(2,2) + (maggrade - maggradi)**2
             endif
          endif
       enddo over_parts
    enddo over_cells
!$omp enddo
!$omp end parallel

    ! shift particles
    xyzh(1:3,1:npart) = xyzh(1:3,1:npart) - dx_shift(:,1:npart)

    ! write diagnostics to file
    if (idebug > 0 .and. id==master) then
       stddev(:,1) = sqrt(stddev(:,1)/npart)
       stddev(:,2) = sqrt(stddev(:,1)/nparterr)

       write(334,'(I18,1x,11(es18.10,1x),I18)') ishift &
                          ,errmin(1,1),errave(1,1)/npart,errmax(1,1) &
                          ,errmin(2,1),errave(2,1)/npart,errmax(2,1) &
                          ,errmin(3,1),errave(3,1)/npart,errmax(3,1) &
                          ,stddev(1:2,1),ncall
       if (iprofile==isphereinbox .or. iprofile==iarray) then
          write(335,'(I18,1x,11(es18.10,1x),I18)') ishift &
                          ,errmin(1,2),errave(1,2)/nparterr,errmax(1,2) &
                          ,errmin(2,2),errave(2,2)/nparterr,errmax(2,2) &
                          ,errmin(3,2),errave(3,2)/nparterr,errmax(3,2) &
                          ,stddev(1:2,2),ncall
       endif
       if (idebug==3) then
          do i = 1,npart
             write(333,'(I18,1x,11(es18.10,1x),I18)') i,xyzh(1:4,i),rhoh(xyzh(4,i),pmass),&
             totalshift(:,i),dx_shift(:,i),ncall
          enddo
          write(333,'(a)') ' '
       endif
    endif

    ! Update the maximum shift; this is summing max shifts per loop, so we're not necessarily using contributions from the same particle
    link_shift = link_shift + sqrt(max_shift2)

    ! Determine if particles are shuffling less than max_shift_thresh of their hi
    if (max_shift2 < max_shift_thresh2) then
       if (call_linklist) shuffle = .false.  ! can only end on an iteration where the linklist was called
       call_linklist = .true.
    else
       if (link_shift > link_shift_thresh) then
          call_linklist = .true.
       else
          call_linklist = .false.
       endif
    endif

    ! update counter; ensure linklist will be called on final loop
    ishift = ishift + 1
    if (ishift == nshiftmax-1) call_linklist = .true.

    if (idebug > 0 .and. id==master) then
       write(iprint,'(1x,a,I5,a,Es18.6)') 'Shuffling: iteration ',ishift,' completed with max shift of ',sqrt(max_shift2)
    endif
 enddo

 ! re-adjust particles to ensure none crossed the boundary
#ifdef PERIODIC
 ncross = 0
!$omp parallel do default (none) &
!$omp shared(npart,xyzh,isperiodic) &
!$omp private(i) &
!$omp reduction(+:ncross)
 do i = 1,npart
    call cross_boundary(isperiodic,xyzh(:,i),ncross)
 enddo
!$omp end parallel do
 if (idebug > 0 .and. ncross > 0 .and. id==master) then
    write(iprint,'(1x,a,I6)') 'Shuffling: number of particles that cross the boundary: ',ncross
 endif
#endif

 ! final debugging print-statements
 if (idebug > 0 .and. id==master) then
    if (idebug > 1) then
       do i = 1,npart
          write(333,'(I18,1x,11(Es18.10,1x),I18)') i,xyzh(1:4,i),rhoh(xyzh(4,i),pmass),totalshift(1:3,i),rthree,ncall
       enddo
    endif
    close(333)
    close(334)
    if (iprofile==isphereinbox .or. iprofile==iarray) close(335)
 endif
 if (id==master) then
    write(iprint,'(1x,3(a,I6),a,Es18.10)') 'Shuffling: completed with ',ishift,' iterations on call number ',ncall, &
                                           ' with ',nlink,' calls to linklist and max shift of ',sqrt(max_shift2)
 endif
 ncall = ncall + 1

end subroutine shuffleparticles
!-------------------------------------------------------------------
end module utils_shuffleparticles
