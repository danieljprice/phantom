!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! This a clump-finding algorithm
!  On the first call, user will be prompted for simulation-dependent values; these will
!  be written to a parameter, clumpfindWB23.in, for subsequent exectutions of the routine.
!  This follows a similar boundness criteria as per Bate (2018) & Wurster, Bate & Price (2019)
!  The clump's mass & energy is summed from all the bound particles
!  We set the spherical radius as \bar{r}+2*sigma to remove the distant, less bound particles;
!  We also fit an ellipse using the moment of inertia method of Rocha, Velho & Carvalho (2002)
!  To do this, we fit three ellipses with the data projected onto the three Cardinal planes. We then rotate
!  the data based upon the ellipse fitting routine and repeat.  When the data is no longer being rotated, we
!  take the final values to determine the lengths of the three semi-axes.
!  WARNING! this currently cannot track across dumps
!  WARNING! there is a bug fix on line labeled 'Used in WB23'; the patch is untested
!
! :References: Wurster & Bonnell (2023) MNRAS, 522, 891-911.
!
! :Owner: James Wurster
!
! :Runtime parameters:
!   - dx_grid    : *grid resolution to determine the clump shape (pc)*
!   - ekin_coef  : *a clump will be bound if coef*ekin + epot < 0*
!   - nclumpmax  : *the maximum number of clumps (recommend 100000 on clusters)*
!   - npartmax   : *maximum number of particles per clump (recommend 100000 on clusters)*
!   - rhominbkg  : *density (cgs) above which a particle can be a member of a clump*
!   - rhominpeak : *density (cgs) above which a particle can be the lead particle in a clump*
!
! :Dependencies: boundary, dim, infile_utils, io, kernel, part, physcon,
!   prompting, sortutils, timing, units
!
 use dim,        only:maxptmass,maxvxyzu,mhd,maxp_hard
 use part,       only:Bxyz,xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,periodic,iorig
 use kernel,     only:radkern,kernel_softening
 use sortutils,  only:indexx
 implicit none
 character(len=20), parameter, public :: analysistype = 'clumpfindWB23'
 public  :: do_analysis

 ! The following parameter/variables need to be set to the specific project
 integer, parameter :: ngrid          =    1000  ! the number of grid points (in one direction) when fitting ellipses to the clumps
 integer            :: nclumpmax      =   10000  ! the maximum number of unique clumps (recommend 100000 is using a computer cluster)
 integer            :: npartmax       =   50000  ! the maximum number of particles per clump clumps (recommend 100000 is using a computer cluster)
 real               :: rhominpeak_cgs = 1.2d-22  ! particles with rho < rhominpeak will not be considered for the lead particle in a clumps (default for WB23)
 real               :: rhominbkg_cgs  = 1.0d-23  ! particles with rho < rhominbkg  will not be considered for membership in a clumps (default for WB23)
 real               :: dxgrid0_pc     = 0.1      ! The default grid spacing (pc) when fitting ellipses to the clumps (default for WB23)
 real               :: ekin_coef      = 1.       ! Will be a clump if ekin_coef*ekin + epot < 0; ekin_coef can be read in (default for WB23)
 real, parameter    :: racc_fac       = 4.0      ! Sinks of distance racc_fac*racc will be checked for boundness during the merger process
 logical            :: debug          = .false.  ! will print useful statements to the log file
 logical            :: print_part     = .true.   ! will create a file for each clump containing all the particle positions
 logical, parameter :: soft_potential = .true.   ! use the kernel softened gravitational potential

 ! The following are common variable not to be modified
 !integer(kind=8)   :: iorigold(maxp_hard) ! will be required when updated to permit tracking across dumps
 integer            :: idclumpold(maxp_hard),idclump(maxp_hard)
 integer, allocatable, dimension(:,:) :: idlistpart(:,:),idlistsink(:,:)
 integer            :: idlistsinkold(maxptmass)
 real               :: dxgrid0
 logical            :: firstcall = .true.   ! required logical; do not change

 type sphclump
    integer           :: nmember,nmemberpv,nptmass,nptmasspv,ID
    real,dimension(3) :: r,v
    real,dimension(6) :: ellipse
    real              :: mass,size
    real              :: kinetic,potential,thermal,magnetic,rhomax,rhoaveln,rhoavelg
 end type sphclump
 type(sphclump), allocatable :: clump(:)

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only:rhoh,massoftype,igas
 use physcon,    only:pc
 use units,      only:udist,umass
 use prompting,  only:prompt
#ifdef PERIODIC
 use boundary,   only:dxbound,dybound,dzbound
#endif
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i,j,j1,j2,k,p,q,n1,nclump,idx,izero,ndense,nclumpcandidate,nclumpprev,ictr,imin,imax,ierr
 integer                      :: nnotcached,nmembermax,lctr
 integer                      :: lst_rho(npart),idense(npart),igrid(ngrid,ngrid)
 real                         :: rhomin_peak,rhomin_bkg
 real                         :: xi,yi,zi,hi,hmax,hmax1,twoh,vxi,vyi,vzi,rhoi,pmassi
 real                         :: dx,dy,dz,dxc,dyc,dzc,dr,dv2,x1,x2,dx0,xc,yc
 real                         :: rad2,size2,sep,ekin,epot,etot,dr_ave,stddev,qp,psoft,fsoft,walltime
 real                         :: rho(npart),rho_dense(npart),rad_pos(npart),rtmp(3),vtmp(3),rold(3)
 real                         :: ella(3),ellb(3),ellp(3),phi(3)
 real                         :: x1grid(ngrid),x2grid(ngrid)
 real, allocatable, dimension(:,:) :: xyz(:,:),eclumpcandidate(:,:)
 real, allocatable, dimension(:)   :: hsmooth(:)
 logical                      :: iexist,connected,connected_global,use_npart,new_config,force_merge,has_changed
 character(len=128)           :: filename,prefix

 !--Look for input file & read.  Else prompt for answers
 if (firstcall) then
    inquire(file='clumpfindWB23.in',exist=iexist)
    if (iexist) call read_clumpparams('clumpfindWB23.in',ierr)
    if (.not. iexist .or. ierr > 0) then
       ! prompt user for initial parameters that are relevant to simulation
       ! defaults are those of Wurster & Bonnell (2023)
       call prompt('Enter the maximum number of clumps (recommend 100000 on clusters) ', nclumpmax, 0)
       call prompt('Enter the maximum number of particles per clump (recommend 100000 on clusters) ', npartmax, 0)
       call prompt('Enter the density (cgs) above which a particle can be the lead particle in a clump', rhominpeak_cgs, 0.0)
       call prompt('Enter the density (cgs) above which a particle can be a member of a clump', rhominbkg_cgs, 0.0)
       call prompt('Enter coef, where a clump will be bound if coef*ekin + epot < 0', ekin_coef, 0.0)
       call prompt('Enter the grid resolution to determine the clump shape (pc) (selected ~ size/100) ', dxgrid0_pc, 0.0)
       call write_clumpparams('clumpfindWB23.in')
    endif
 endif

 !--Allocate arrays
 allocate(clump(nclumpmax))
 allocate(idlistpart(nclumpmax,npartmax))
 allocate(idlistsink(nclumpmax,maxptmass))
 allocate(eclumpcandidate(3,nclumpmax))
 !--Initialise values
 izero       = 0
 nclump      = 0
 nclumpprev  = -1
 new_config  = .true.
 rhomin_peak = rhominpeak_cgs*udist**3/umass
 rhomin_bkg  = rhominbkg_cgs*udist**3/umass
 dxgrid0     = dxgrid0_pc*pc/udist
 idclump     = 0
 idlistsink  = 0
 idlistpart  = 0
 pmassi      = massoftype(igas)
 idx         = index(dumpfile,'_')
 prefix      = dumpfile(1:idx-1)
 walltime    = 0.
 if (firstcall) then
    idclumpold    = 0
    idlistsinkold = 0
    !iorigold     = iorig
 endif
 call update_time(walltime)

 !--Print free parameters
 write(*,'(a)')             ' '
 write(*,'(a)')             '==============================================================='
 write(*,'(a)')             '  Free parameters for clumpfindWB23'
 write(*,'(a)')             ' '
 write(*,'(a,I10)')         '  The maximum number of unique clumps:                                              ',nclumpmax
 write(*,'(a,I10)')         '  The maximum number of particles per clump before losing computational efficiency: ',npartmax
 write(*,'(a,Es10.3,a)')    '  The minimum density required to be the lead particle in a clump: ',rhominpeak_cgs,' g/cm^3'
 write(*,'(a,Es10.3,a)')    '  The minimum density required to a member of a clump:             ',rhominbkg_cgs, ' g/cm^3'
 write(*,'(a,f7.2,a)')      '  Groups of particles are considered a clump if ',ekin_coef,'*ekin + epot < 0'
 if (dxgrid0_pc < 0.01 .or. dxgrid0_pc > 1000.) then
    write(*,'(a,Es10.3,a)') '  The resolution of the ellipse-fitting routine to determine the clump size: ',dxgrid0_pc,' pc'
 else
    write(*,'(a,f7.2,a)')   '  The resolution of the ellipse-fitting routine to determine the clump size: ',dxgrid0_pc,' pc'
 endif
 write(*,'(a)')             '==============================================================='
 write(*,'(a)')             ' WARNING! there is a bug fix on line labeled "Used in WB23"; the patch is untested'
 write(*,'(a)')             '==============================================================='
 write(*,'(a)')             ' '

 !--calculate all densities
 rho    = 0.0
 idense = 0
 ndense = 0
 do i = 1,npart
    hi = xyzh(4,i)
    if (hi > 0.0) then
       rho(i) = rhoh(hi,particlemass)
       if (rho(i) > rhomin_bkg) then
          ndense            = ndense + 1
          idense(ndense)    = i
          rho_dense(ndense) = rho(i)
       endif
    endif
 enddo
 print*, 'Dense gas located: ',ndense,' particles'

 !--Sort the dense particles by position and density
 call indexx(ndense,rho_dense,lst_rho)
 print*, 'array sorting complete'

 !--Initialise all clumps
 do i = 1,nclumpmax
    call reset_clump(i)
 enddo

 !--Initialise the clumps around sinks
 do i = 1,nptmass
    nclump = nclump + 1
    if (nclump > nclumpmax) then
       print*, 'recompile with larger nclumpmax.  aborting'
       return
    endif
    if (xyzmh_ptmass(4,nclump) > 0.) then ! accounts for dead sinks from mergers
       ! initialise clump
       call reset_clump(nclump)
       clump(nclump)%ID      = nclump
       clump(nclump)%r(1:3)  = xyzmh_ptmass(1:3,i)
       clump(nclump)%v(1:3)  = vxyz_ptmass(1:3,i)
       clump(nclump)%mass    = xyzmh_ptmass(4,i)
       clump(nclump)%size    = xyzmh_ptmass(ihacc,i)*racc_fac
       clump(nclump)%nptmass = 1
       idlistsink(nclump,1)  = i
    else
       nclump = nclump - 1
    endif
 enddo
 print*, 'initialised clumps around sinks'
 call update_time(walltime)

 !--Find the clumps around dense gas particles
 ictr = 1
 do while (nclump /= nclumpprev .or. new_config)
    print*, ' '
    print*, 'Finding clumps: Loop ',ictr
    ictr = ictr + 1
    nclumpprev = nclump
    idclumpold = idclump

    do k = ndense,1,-1
       i = idense(lst_rho(k))
       if (idclump(i) > 0) cycle
       xi   = xyzh(1,i)
       yi   = xyzh(2,i)
       zi   = xyzh(3,i)
       hi   = xyzh(4,i)
       twoh = radkern*hi
       vxi  = vxyzu(1,i)
       vyi  = vxyzu(2,i)
       vzi  = vxyzu(3,i)
       rhoi = rho(i)
       eclumpcandidate = 0.
       nclumpcandidate = 0
       sep  = 0.
       epot = 0.
       ekin = 0.
       connected_global = .false.

       ! search through all clumps; if we are close enough, verify energies & connectedness
!$omp parallel default(none) &
!$omp shared(npart,nclump,xyzh,xyzmh_ptmass,xi,yi,zi,hi,twoh,vxi,vyi,vzi,eclumpcandidate,idlistpart,clump,ekin_coef,k) &
!$omp shared(pmassi,connected_global,idlistsink) &
!$omp private(j,p,q,rad2,connected,ekin,epot,sep,dv2,dr,dx,dy,dz,dxc,dyc,dzc,n1,use_npart,qp,psoft,fsoft,hmax,hmax1) &
!$omp reduction(+:nclumpcandidate)
!$omp do schedule(runtime)
       do j = 1,nclump
          rad2 = (xi-clump(j)%r(1))**2 + (yi-clump(j)%r(2))**2 + (zi-clump(j)%r(3))**2
          if (rad2 < (clump(j)%size + twoh)**2) then
             epot = 0.
             sep  = rad2
             dv2  = (vxi-clump(j)%v(1))**2 + (vyi-clump(j)%v(2))**2 + (vzi-clump(j)%v(3))**2
             connected = .false.
             call n_in_clump(j,npart,n1,use_npart)
             do q = 1,n1
                p = p_val(j,q,use_npart)
                if (p > 0) then
                   dx   = xyzh(1,p) - xi
                   dy   = xyzh(2,p) - yi
                   dz   = xyzh(3,p) - zi
                   hmax = max(xyzh(4,p),hi)
                   dxc  = xyzh(1,p) - clump(j)%r(1)
                   dyc  = xyzh(2,p) - clump(j)%r(2)
                   dzc  = xyzh(3,p) - clump(j)%r(3)
                   dr   = sqrt(dx*dx + dy*dy + dz*dz)
                   sep  = max(sep,dxc*dxc + dyc*dyc + dzc*dzc) ! this is to approximate the size of the cloud for sampling reasons
                   if (dr > radkern*hmax .or. .not.soft_potential) then
                      if (dr > 0.) epot = epot - pmassi/dr
                   else
                      hmax1 = 1.0/hmax
                      qp    = dr*hmax1
                      call kernel_softening(qp*qp,qp,psoft,fsoft)
                      epot = epot + pmassi*psoft*hmax1
                   endif
                   if (dr < radkern*hi) connected = .true.
                endif
             enddo
             if (clump(j)%nptmass > 0) then
                do q = 1,clump(j)%nptmass
                   p    = idlistsink(j,q)
                   dx   = xyzmh_ptmass(1,p) - xi
                   dy   = xyzmh_ptmass(2,p) - yi
                   dz   = xyzmh_ptmass(3,p) - zi
                   hmax = max(0.5*xyzmh_ptmass(ihacc,p),hi)
                   dxc  = xyzmh_ptmass(1,p) - clump(j)%r(1)
                   dyc  = xyzmh_ptmass(2,p) - clump(j)%r(2)
                   dzc  = xyzmh_ptmass(3,p) - clump(j)%r(3)
                   dr   = sqrt(dx*dx + dy*dy + dz*dz)
                   sep  = max(sep,dxc*dxc + dyc*dyc + dzc*dzc)
                   if (dr > radkern*hmax .or. .not.soft_potential) then
                      if (dr > 0.) epot = epot - xyzmh_ptmass(4,p)/dr
                   else
                      hmax1 = 1.0/hmax
                      qp    = dr*hmax1
                      call kernel_softening(qp*qp,qp,psoft,fsoft)
                      epot = epot + xyzmh_ptmass(4,p)*psoft*hmax1
                   endif
                   if (dr < radkern*hi .or. dr < xyzmh_ptmass(ihacc,p)*racc_fac) connected = .true.
                enddo
             endif
             !ekin = 0.5*(pmassi*clump(j)%mass)/(pmassi+clump(j)%mass)*dv2  ! Used in WB23
             ekin = clump(j)%kinetic   + 0.5*(pmassi*clump(j)%mass)/(pmassi+clump(j)%mass)*dv2  ! likely the correct version (untested)
             epot = clump(j)%potential + epot*pmassi
             if (connected .and. ekin_coef*ekin + epot < 0.) then
                eclumpcandidate(1,j) = ekin
                eclumpcandidate(2,j) = epot
                eclumpcandidate(3,j) = sep
                nclumpcandidate      = nclumpcandidate + 1
             endif
             if (connected) connected_global = .true.
          endif
       enddo
!$omp enddo
!$omp end parallel

       if (nclumpcandidate > 1) then
          ! there are multiple clumps; see if the two that particle i is most bound to are bound.
          etot = 0.
          j1   = 1.
          do j = 1,nclump
             ekin = eclumpcandidate(1,j)
             epot = eclumpcandidate(2,j)
             if (ekin_coef*ekin + epot < etot) then
                j1   = j
                etot = ekin_coef*ekin + epot
             endif
          enddo
          etot = 0.
          j2   = nclump
          do j = 1,nclump
             ekin = eclumpcandidate(1,j)
             epot = eclumpcandidate(2,j)
             if (ekin_coef*ekin + epot < etot .and. j /= j1) then
                j2   = j
                etot = ekin_coef*ekin + epot
             endif
          enddo
          ! the two clumps are j1 & j2
          call are_clumps_joined(j1,j2,npart,xyzh,vxyzu,rtmp,vtmp,sep,pmassi,ekin,epot,i,xi,yi,zi,hi,vxi,vyi,vzi)
          if ( ekin_coef*ekin + epot < 0. ) then
             if (debug) write(*,'(3(a,I10))') 'merging clumps ',j1,' & ',j2,' and particle ', i
             call merge_clumps(j1,j2,npart,nclump,vxyzu,Bxyz,rho,rtmp,vtmp,sep,pmassi,ekin,epot,i)
          else
             ! associate particle with only one clump
             nclumpcandidate = 1
             do j = 1,nclump
                if (j /= j1) eclumpcandidate(:,j) = 0
             enddo
          endif
       endif

       if (nclumpcandidate==0 .and. .not.connected_global .and. rhoi > rhomin_peak) then
          ! create a new clump
          nclump = nclump + 1
          if (nclump > nclumpmax) then
             print*, 'recompile with larger nclumpmax.  aborting'
             return
          endif
          call reset_clump(nclump)
          clump(nclump)%ID      = nclump
          clump(nclump)%r(1:3)  = xyzh(1:3,i)
          clump(nclump)%v(1:3)  = vxyzu(1:3,i)
          clump(nclump)%mass    = pmassi
          clump(nclump)%size    = radkern*hi
          clump(nclump)%rhomax  = rhoi
          clump(nclump)%nmember = 1
          idlistpart(nclump,1)  = i
          idclump(i)            = nclump
          if (maxvxyzu>=4) clump(nclump)%thermal  = vxyzu(4,i)*pmassi
          if (mhd)         clump(nclump)%magnetic = 0.5*pmassi*dot_product(Bxyz(1:3,i),Bxyz(1:3,i))/rhoi
       elseif (nclumpcandidate==1) then
          ! add particle to clump
          do j = 1,nclump
             if (eclumpcandidate(2,j) < 0.) then
                rold = clump(j)%r(:)
                do p=1,3
                   clump(j)%r(p)  = (clump(j)%mass*clump(j)%r(p) + pmassi*xyzh( p,i))/(clump(j)%mass + pmassi)
                   clump(j)%v(p)  = (clump(j)%mass*clump(j)%v(p) + pmassi*vxyzu(p,i))/(clump(j)%mass + pmassi)
                enddo
                dx = rold(1) - clump(j)%r(1)
                dy = rold(2) - clump(j)%r(2)
                dz = rold(3) - clump(j)%r(3)
                clump(j)%mass      = clump(j)%mass + pmassi
                clump(j)%size      = sqrt(eclumpcandidate(3,j)) + sqrt(dx*dx+dy*dy+dz*dz)
                clump(j)%nmember   = clump(j)%nmember + 1
                clump(j)%kinetic   = eclumpcandidate(1,j)
                clump(j)%potential = eclumpcandidate(2,j)
                clump(j)%rhomax    = max(clump(j)%rhomax,rhoi)
                idclump(i)         = j
                if (clump(j)%nmember <= npartmax) idlistpart(j,clump(j)%nmember) = i
                if (maxvxyzu>=4) clump(j)%thermal  = clump(j)%thermal  + vxyzu(4,i)*pmassi
                if (mhd)         clump(j)%magnetic = clump(j)%magnetic + 0.5*pmassi*dot_product(Bxyz(1:3,i),Bxyz(1:3,i))/rhoi
             endif
          enddo
       endif
    enddo  ! end of k-loop
    print*, 'Loop ',ictr, ': done adding particles to existing clumps. nclump = ',nclump
    call update_time(walltime)


    !--Merge clumps
    !  the order we do this should be irrelevant since if multiple cores are bound, then they should all merge, independent of order
    !  if not, it should be picked up on subsequent loops
    j1 = 1
    do while (j1 <= nclump-1)
       j2 = j1 + 1
       xi = clump(j1)%r(1)
       yi = clump(j1)%r(2)
       zi = clump(j1)%r(3)
       has_changed = (clump(j1)%nmember /= clump(j1)%nmemberpv .or. clump(j1)%nptmass /= clump(j1)%nptmasspv)
       do while (j2 <= nclump)
          if (has_changed .or. clump(j2)%nmember /= clump(j2)%nmemberpv .or. clump(j2)%nptmass /= clump(j2)%nptmasspv) then
             dx    = xi-clump(j2)%r(1)
             dy    = yi-clump(j2)%r(2)
             dz    = zi-clump(j2)%r(3)
             rad2  = dx*dx + dy*dy + dz*dz
             size2 = 0.0
             if (clump(j1)%nptmass == 1) size2 = (xyzmh_ptmass(ihacc,idlistsink(j1,1))*racc_fac)**2
             if (clump(j2)%nptmass == 1) size2 = (xyzmh_ptmass(ihacc,idlistsink(j2,1))*racc_fac)**2
             force_merge = .false.
             !  Force the merger if the centre of mass of a sinkless clump is within the clump with only one sink
             if (clump(j1)%nptmass + clump(j2)%nptmass == 1 .and. rad2 < size2) force_merge = .true.

             ! determine if the clumps are bound, and merge if they are
             if (rad2 < (clump(j1)%size+clump(j2)%size)**2 .or. force_merge) then
                call are_clumps_joined(j1,j2,npart,xyzh,vxyzu,rtmp,vtmp,sep,pmassi,ekin,epot)
                if (debug) write(*,'(2(a,I6),a,Es18.6)') 'Clumps ',j1,' & ',j2,' are connected: etot = ',ekin_coef*ekin + epot
                if ( ekin_coef*ekin + epot < 0. .or. force_merge) then
                   if (debug) then
                      if (force_merge) then
                         write(*,'(2(a,I10),a,Es18.6)') 'forced merging clumps ',j1,' & ',j2,': E = ',ekin_coef*ekin + epot
                      else
                         write(*,'(2(a,I10))') 'merging clumps ',j1,' & ',j2
                      endif
                   endif
                   call merge_clumps(j1,j2,npart,nclump,vxyzu,Bxyz,rho,rtmp,vtmp,sep,pmassi,ekin,epot)
                endif
             endif ! if within the same sizes
          endif ! if we have not previously tested this pair
          j2 = j2 + 1
       enddo
       j1 = j1 + 1
    enddo
    print*, 'Loop ',ictr, ': done merging clumps. nclump = ',nclump
    call update_time(walltime)

    new_config = .false.
!$omp parallel default(none) &
!$omp shared(npart,idclump,idclumpold,new_config,ictr) &
!$omp private(i)
!$omp do schedule(runtime)
    do i = 1,npart
       if (idclump(i)/= idclumpold(i)) then
          new_config = .true.
          if (ictr > 35) print*, 'NC',i, idclump(i), idclumpold(i)
       endif
    enddo
!$omp enddo
!$omp end parallel

    nnotcached = 0
!$omp parallel default(none) &
!$omp shared(nclump,clump,npartmax) &
!$omp private(i) &
!$omp reduction(+:nnotcached)
!$omp do schedule(runtime)
    do i = 1,nclump
       if (clump(i)%nmember > npartmax) nnotcached = nnotcached + 1
    enddo
!$omp enddo
!$omp end parallel

    clump(:)%nmemberpv = clump(:)%nmember ! so we can avoid repeated calculations to see is the same pair are bound
    clump(:)%nptmasspv = clump(:)%nptmass ! so we can avoid repeated calculations to see is the same pair are bound
    print*, 'nclump, nclumpprev = ',nclump, nclumpprev
    print*, 'nclumps with nmember > npartmax = ',nnotcached
    !--Repeat the entire process until the number of clumps stabilises
 enddo
 call update_time(walltime)

 !--Remove all clumps with 57 or fewer particles
 !  Do this here since poorly-resolved clumps may grow & become resolved on subsequent iterations through the density list
 j = 1
 do while (j <= nclump)
    if (clump(j)%nmember < 58 .and. clump(j)%nptmass==0) then
       if (debug) print*, 'removing clump ',j,' which has ',clump(j)%nmember,' neighbours'
       ! remove clump
       do i = 1,npart
          if (idclump(i)==j) idclump(i) = 0
       enddo
       if (j < nclump) then
          ! move clump nclump to clump j
          do i = 1,npart
             if (idclump(i)==nclump) idclump(i) = j
          enddo
          clump(j)        = clump(nclump)
          clump(j)%ID     = j
          idlistpart(j,:) = idlistpart(nclump,:)
          idlistsink(j,:) = idlistsink(nclump,:)
       endif
       nclump = nclump - 1
    else
       j = j + 1
    endif
 enddo
 print*, 'Done removing unresolved clumps'
 call update_time(walltime)

 nmembermax = 0
 do j = 1,nclump
    nmembermax = max(nmembermax,clump(j)%nmember+clump(j)%nptmass)
 enddo
 print*, 'The maximum number of clump members is ',nmembermax
 allocate(xyz(3,nmembermax))
 allocate(hsmooth(nmembermax))
 rad_pos = 0.

! Calculate a reasonable radius & average density
! The spherical radius is done by selecting the radius of \bar{r} + 2sigma
! The elliptical parameters are set rotating the data until the ellipse lies along the cardinal axes,
! then using three 2D ellipse fitting routine of Rocha, Velho & Carvalho (2002) who use Moments of Inertia to fit the ellipse.
! Not all particles are fit within the ellipse, and visual inspection suggests these are good fits.
!$omp parallel default(none) &
!$omp shared(npart,nclump,pmassi,clump,xyzh,xyzmh_ptmass,idlistsink,idclump,rho,rad_pos,debug) &
#ifdef PERIODIC
!$omp shared(dxbound,dybound,dzbound) &
#endif
!$omp private(j,q,p,dx,dy,dz,n1,use_npart,ictr,dr_ave,stddev,xyz,ella,ellb,ellp,phi,x1,x2,imin,imax) &
!$omp private(lctr,dx0,x1grid,x2grid,hsmooth,xc,yc,igrid)
!$omp do schedule(runtime)
 do j = 1,nclump
    clump(j)%size = 0.
    call n_in_clump(j,npart,n1,use_npart)
    if (debug) then
       write(*,'(3(a,I10),a,L)') 'start: calculating radius for clump ',j,' of ',nclump,'. n1 = ',n1,' & use_npart = ',use_npart
    endif
    ictr = 0
    if (n1 > 0) then
       dr_ave = 0.
       do q = 1,n1
          p = p_val(j,q,use_npart)
          if (p > 0) then
             dx = xyzh(1,p) - clump(j)%r(1)
             dy = xyzh(2,p) - clump(j)%r(2)
             dz = xyzh(3,p) - clump(j)%r(3)
#ifdef PERIODIC
             if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
             if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
             if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
             rad_pos(p) = sqrt(dx*dx + dy*dy + dz*dz)
             dr_ave     = dr_ave + rad_pos(p)
             clump(j)%rhoaveln = clump(j)%rhoaveln + rho(p)
             clump(j)%rhoavelg = clump(j)%rhoavelg + log10(rho(p))
             ictr        = ictr + 1
             xyz(1,ictr) = dx
             xyz(2,ictr) = dy
             xyz(3,ictr) = dz
             hsmooth(ictr) = xyzh(4,p)
          endif
       enddo
       if (ictr > 0) then
          dr_ave = dr_ave/ictr
          clump(j)%rhoaveln = clump(j)%rhoaveln/ictr
          clump(j)%rhoavelg = 10**(clump(j)%rhoavelg/ictr)
       else
          print*, 'This should not be possible!'
       endif

       stddev = 0.
       do q = 1,n1
          p = p_val(j,q,use_npart)
          if (p > 0) then
             stddev = stddev + (rad_pos(p) - dr_ave)**2
          endif
       enddo
       stddev = sqrt(stddev/ictr)

       ! find the separations of the particles within the radial criteria
       do q = 1,n1
          p = p_val(j,q,use_npart)
          if (p > 0) then
             if (rad_pos(p) < dr_ave+2.0*stddev) then
                clump(j)%size = max(clump(j)%size,rad_pos(p)**2)
                idclump(p) =  abs(idclump(p))
             else
                idclump(p) = -abs(idclump(p))
             endif
          endif
       enddo
    endif

    ! Determine size constraints based upon the sink particle
    do q = 1,clump(j)%nptmass
       p  = idlistsink(j,q)
       dx = xyzmh_ptmass(1,p) - clump(j)%r(1)
       dy = xyzmh_ptmass(2,p) - clump(j)%r(2)
       dz = xyzmh_ptmass(3,p) - clump(j)%r(3)
#ifdef PERIODIC
       if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
       if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
       if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
       clump(j)%size  = max(clump(j)%size,dx*dx + dy*dy + dz*dz)

       ictr = ictr + 1
       xyz(1,ictr) = dx
       xyz(2,ictr) = dy
       xyz(3,ictr) = dz
       hsmooth(ictr) = xyzmh_ptmass(ihacc,p)
    enddo
    clump(j)%size = sqrt(clump(j)%size)

    ! Fit the ellipse if we have particles
    if (n1 > 0) then

       ellp = 1.
       phi  = 0.
       lctr = 0
       do while (lctr < 1000 .and. (abs(ellp(1)) > 1e-3 .or. abs(ellp(2)) > 1e-3 .or. abs(ellp(3)) > 1e-3))
          lctr = lctr + 1

          call set_grid(ictr,1,2,xyz,hsmooth,dx0,x1grid,x2grid)
          call fill_igrid(ictr,1,2,dx0,xyz,hsmooth,x1grid,x2grid,igrid)
          call calc_moment(x1grid,x2grid,igrid,ella(1),ellb(1),ellp(1),xc,yc)
          call rotate_data(ictr,1,2,xc,yc,ellp(1),xyz)

          call set_grid(ictr,1,3,xyz,hsmooth,dx0,x1grid,x2grid)
          call fill_igrid(ictr,1,3,dx0,xyz,hsmooth,x1grid,x2grid,igrid)
          call calc_moment(x1grid,x2grid,igrid,ella(2),ellb(2),ellp(2),xc,yc)
          call rotate_data(ictr,1,3,xc,yc,ellp(2),xyz)

          call set_grid(ictr,2,3,xyz,hsmooth,dx0,x1grid,x2grid)
          call fill_igrid(ictr,2,3,dx0,xyz,hsmooth,x1grid,x2grid,igrid)
          call calc_moment(x1grid,x2grid,igrid,ella(3),ellb(3),ellp(3),xc,yc)
          call rotate_data(ictr,2,3,xc,yc,ellp(3),xyz)

          phi = phi + ellp
       enddo


       clump(j)%ellipse(4:6) = phi

       ! Although the ellipse are rotated to be on the Cartesian axes before the final ellipse fit, the gridding means
       ! that the projections are not perfect and we do not have three independent length.  We pair them off, and take the
       ! larger value for each axis.

       imin = minloc(ella,1)
       imax = maxloc(ella,1)
       clump(j)%ellipse(1) = ella(imax)*0.5  ! The semi-major axis
       clump(j)%ellipse(2) = ella(imin)*0.5  ! The semi-middle axis
       imin = minloc(ellb,1)
       imax = maxloc(ellb,1)
       ellb(imin) = ellb(imax)
       imin = minloc(ellb,1)
       clump(j)%ellipse(3) = ellb(imin)*0.5  ! The semi-minor axis

       if (debug) write(*,'(a,9Es18.6)') 'ellipse:',ella,ellb,clump(j)%ellipse(1:3)
    endif

    ! Ensure that the ellipse is at least as big at the sink radius
    ! This should be redundant, but leave it in anyway since it won't hurt anything
    if (clump(j)%nptmass > 0) then
       p = idlistsink(j,1)
       clump(j)%size       = max(clump(j)%size,      xyzmh_ptmass(ihacc,p))
       clump(j)%ellipse(1) = max(clump(j)%ellipse(1),xyzmh_ptmass(ihacc,p))
       clump(j)%ellipse(2) = max(clump(j)%ellipse(2),xyzmh_ptmass(ihacc,p))
       clump(j)%ellipse(3) = max(clump(j)%ellipse(3),xyzmh_ptmass(ihacc,p))
    endif

    if (debug) write(*,'(a,I10)') 'done:  calculating radius for clump ',j
 enddo
!$omp enddo
!$omp end parallel
 deallocate(xyz)
 deallocate(hsmooth)
 print*, 'Done recalculating size'
 call update_time(walltime)

 !--Write results to file (both all the clumps at the current time, and to the file for each clump)
 write(filename,'(2a)') trim(dumpfile),'clumps'
 open(unit=iunit,file=trim(filename))
 write(iunit,'(a,I6,a,Es18.6)') '#Nclumps = ',nclump,'; Time = ',time
 write(iunit,"('#',24(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'clump ID',    &
       2,'N gas',       &
       3,'N sink',      &
       4,'x',           &
       5,'y',           &
       6,'z',           &
       7,'vx',          &
       8,'vy',          &
       9,'vz',          &
      10,'mass',        &
      11,'rho_max',     &
      12,'rho_ave,lin', &
      13,'rho_ave,log', &
      14,'E_kinetic',   &
      15,'E_potential', &
      16,'E_thermal',   &
      17,'E_magnetic',  &
      18,'dr sphere',   &
      19,'ell a',       &
      20,'ell b',       &
      21,'ell c',       &
      22,'ell phi_xy',  &
      23,'ell phi_xz',  &
      24,'ell phi_yz'

 ictr = 0
 do i = 1,nclump
    if (clump(i)%nmember==0) ictr = ictr + 1
    write(iunit,'(3(I18,1x),21(Es18.6,1x))') clump(i)%ID,clump(i)%nmember,clump(i)%nptmass,                          &
                                             clump(i)%r(:),clump(i)%v(:),clump(i)%mass,                              &
                                             clump(i)%rhomax,clump(i)%rhoaveln,clump(i)%rhoavelg,                    &
                                             clump(i)%kinetic,clump(i)%potential,clump(i)%thermal,clump(i)%magnetic, &
                                             clump(i)%size,clump(i)%ellipse
 enddo
 close(iunit)

 ! Print characteristics for each dump.  This file is not currently useful since we cannot track clumps across dumps
 do i = 1,nclump
    if (.false.) then
       write(filename,'(2a,I5.5)') trim(prefix),'clumpID',clump(i)%ID
       inquire(file=trim(filename),exist=iexist)
       if (iexist) then
          open(iunit,file=trim(filename), position='append')
       else
          open(iunit,file=trim(filename))
          write(iunit,"('#',26(1x,'[',i2.2,1x,a11,']',2x))") &
             1,'time',        &
             2,'idump',       &
             3,'clump ID',    &
             4,'N gas',       &
             5,'N sink',      &
             6,'x',           &
             7,'y',           &
             8,'z',           &
             9,'vx',          &
            10,'vy',          &
            11,'vz',          &
            12,'mass',        &
            13,'rho_max',     &
            14,'rho_ave,lin', &
            15,'rho_ave,log', &
            16,'E_kinetic',   &
            17,'E_potential', &
            18,'E_thermal',   &
            19,'E_magnetic',  &
            20,'dr sphere',   &
            21,'ell a',       &
            22,'ell b',       &
            23,'ell c',       &
            24,'ell phi_xy',  &
            25,'ell phi_xz',  &
            26,'ell phi_yz'
       endif
       write(iunit,'(Es18.6,1x,4(I18,1x),21(Es18.6,1x))') time,num,clump(i)%ID,clump(i)%nmember,clump(i)%nptmass,   &
                                                          clump(i)%r(:),clump(i)%v(:),clump(i)%mass,                &
                                                          clump(i)%rhomax,clump(i)%rhoaveln,clump(i)%rhoavelg,      &
                                                          clump(i)%kinetic,clump(i)%potential,                      &
                                                          clump(i)%thermal,clump(i)%magnetic,                       &
                                                          clump(i)%size,clump(i)%ellipse
       close(iunit)
    endif

    if (print_part) then
       ! print the particles in each clump to file
       ! column 1: >0: gas; < 0: sink
       ! column 6: >0 means within the defined radius; < 0 means outside the defined radius
       write(filename,'(2a,I5.5)') trim(prefix),'clump',i
       inquire(file=trim(filename),exist=iexist)
       open(iunit,file=trim(filename))

       call n_in_clump(i,npart,n1,use_npart)
       write(iunit,'(2(a,I8))')   '#N = ',clump(i)%nmember,' ; nsink = ',clump(i)%nptmass
       write(iunit,'(a, Es18.6)') '#spherical radius = ',clump(i)%size
       write(iunit,'(a,6Es18.6)') '#ellipse(a,b,c,phi_xy,phy_xz,phi_yz) = ',clump(i)%ellipse
       write(iunit,'(a,3Es18.6)') '#origin     = ',clump(i)%r(:)
       write(iunit,'(a)') ' '
       do j = 1,n1
          p = p_val(i,j,use_npart)
          if (p > 0) write(iunit,'(I12,5Es18.6,I8)') p, xyzh(1:4,p),rho(p),idclump(p)
       enddo
       do j = 1,clump(i)%nptmass
          p = idlistsink(i,j)
          write(iunit,'(I12,5Es18.6,I8)') -p, xyzmh_ptmass(1:3,p),xyzmh_ptmass(ihacc,p),xyzmh_ptmass(4,p),clump(i)%ID
       enddo
       close(iunit)
    endif
 enddo
 print*, 'There are ',nclump,', of which ',ictr,' have only sink particles'

 !--Save values for next dumpfile; we should do this in a stand-alone loop here
 !idclumpold     = idclump
 !idlistsinkold = idlistsink
 !iorigold       = iorig
 firstcall      = .false.
 call update_time(walltime)

 deallocate(clump)
 deallocate(idlistpart)
 deallocate(idlistsink)
 deallocate(eclumpcandidate)

end subroutine do_analysis
!--------------------------------------------------------------------------
!--Zero the clump
subroutine reset_clump(i)
 integer, intent(in) :: i

 clump(i)%r(:)       = 0.0
 clump(i)%v(:)       = 0.0
 clump(i)%ellipse(:) = 0.0
 clump(i)%ID         = 0
 clump(i)%mass       = 0.0
 clump(i)%size       = 0.0
 clump(i)%nmember    = 0
 clump(i)%nmemberpv  = 0
 clump(i)%kinetic    = 0.0
 clump(i)%potential  = 0.0
 clump(i)%thermal    = 0.0
 clump(i)%magnetic   = 0.0
 clump(i)%rhomax     = 0.0
 clump(i)%rhoaveln   = 0.0
 clump(i)%rhoavelg   = 0.0
 clump(i)%nptmass    = 0
 clump(i)%nptmasspv  = 0

end subroutine reset_clump
!--------------------------------------------------------------------------
! This will determine if two clumps should be merged, either due to themselves or a joining particle
subroutine are_clumps_joined(j1,j2,npart,xyzh,vxyzu,rtmp,vtmp,sep,pmassi,ekin,epot,i,xi,yi,zi,hi,vxi,vyi,vzi)
 integer, intent(in)           :: j1,j2,npart
 integer, intent(in), optional :: i
 real,    intent(in)           :: xyzh(:,:),vxyzu(:,:),pmassi
 real,    intent(out)          :: ekin,epot,sep,rtmp(:),vtmp(:)
 real,    intent(in), optional :: xi,yi,zi,hi,vxi,vyi,vzi
 integer                       :: ii,j,q,p1,p2,zero_one,n1,n2
 real                          :: pmassii,dx,dy,dz,dr,dvx,dvy,dvz,xii,yii,zii,hii,vxii,vyii,vzii,ekini,epoti,sepi
 real                          :: qp,psoft,fsoft,hmax,hmax1,pmassi2
 logical                       :: connected,use_npart1,use_npart2

 if (present(xi) .and. present(yi) .and. present(zi) .and. present(vxi) .and. present(vyi) .and. present(vzi)) then
    ! We are merging two clumps & a joining particle
    pmassii  = pmassi
    xii      = xi
    yii      = yi
    zii      = zi
    hii      = hi
    vxii     = vxi
    vyii     = vyi
    vzii     = vzi
    ii       = i
    zero_one = 1
 else
    ! We are merging two clumps only
    pmassii  = 0.
    xii      = 0.
    yii      = 0.
    zii      = 0.
    hii      = 0.
    vxii     = 0.
    vyii     = 0.
    vzii     = 0.
    ii       = 1 ! just a dummy value, but will always be multiplied by zero
    zero_one = 0
 endif
 pmassi2 = pmassi*pmassi

 !--determine how we will be looping through the particles based upon clump membership
 call n_in_clump(j1,npart,n1,use_npart1)
 call n_in_clump(j2,npart,n2,use_npart2)

 !--calculate the centre of mass & velocity
 do q = 1,3
    rtmp(q) = (clump(j1)%mass*clump(j1)%r(q) + clump(j2)%mass*clump(j2)%r(q) + pmassii*xyzh(q,ii)) &
                    / (clump(j1)%mass + clump(j2)%mass + pmassii)

    vtmp(q) = (clump(j1)%mass*clump(j1)%v(q) + clump(j2)%mass*clump(j2)%v(q) + pmassii*vxyzu(q,ii)) &
            / (clump(j1)%mass                + clump(j2)%mass                + pmassii)
 enddo

 !--calculate potential energy
 connected = .false.
 epot = 0.
 ekin = 0.
!$omp parallel default(none) &
!$omp shared(j1,j2,n1,n2,clump,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,idlistsink,connected) &
!$omp shared(xii,yii,zii,hii,vxii,vyii,vzii,pmassii,pmassi,pmassi2,zero_one,use_npart1,use_npart2) &
!$omp private(q,p1,p2,dx,dy,dz,dr,dvx,dvy,dvz,qp,psoft,fsoft,hmax,hmax1) &
!$omp reduction(+:epot,ekin)
!$omp do schedule(runtime)
 jloop: do j = 1,n1
    p1 = p_val(j1,j,use_npart1)
    if (p1 == 0) cycle jloop
    qloop: do q = 1,n2
       p2 = p_val(j2,q,use_npart2)
       if (p2 == 0) cycle qloop
       dx   = xyzh(1,p1) - xyzh(1,p2)
       dy   = xyzh(2,p1) - xyzh(2,p2)
       dz   = xyzh(3,p1) - xyzh(3,p2)
       hmax = max(xyzh(4,p1),xyzh(4,p2))
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epot = epot - pmassi2/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epot  = epot + pmassi2*psoft*hmax1
       endif
       if (dr < radkern*max(xyzh(4,p1),xyzh(4,p2))) connected = .true.  ! looking for neighbours, not overlapping h's
       dvx  = vxyzu(1,p1) - vxyzu(1,p2)
       dvy  = vxyzu(2,p1) - vxyzu(2,p2)
       dvz  = vxyzu(3,p1) - vxyzu(3,p2)
       ekin = ekin + 0.25*pmassi*(dvx*dvx + dvy*dvy + dvz*dvz) ! since using reduced mass: (pmassi*pmassi/(pmassi+pmassi)) = 0.5*pmassi
    enddo qloop
    qloopsink: do q = 1,clump(j2)%nptmass
       p2   = idlistsink(j2,q)
       dx   = xyzh(1,p1) - xyzmh_ptmass(1,p2)
       dy   = xyzh(2,p1) - xyzmh_ptmass(2,p2)
       dz   = xyzh(3,p1) - xyzmh_ptmass(3,p2)
       hmax = max(xyzh(4,p1),0.5*xyzmh_ptmass(ihacc,p2))
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epot = epot - pmassi*xyzmh_ptmass(4,p2)/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epot  = epot + pmassi*xyzmh_ptmass(4,p2)*psoft*hmax1
       endif
       if (dr < max(radkern*xyzh(4,p1),xyzmh_ptmass(ihacc,p2)*racc_fac)) connected = .true.
       dvx  = vxyzu(1,p1) - vxyz_ptmass(1,p2)
       dvy  = vxyzu(2,p1) - vxyz_ptmass(2,p2)
       dvz  = vxyzu(3,p1) - vxyz_ptmass(3,p2)
       ekin = ekin + 0.5*(dvx*dvx + dvy*dvy + dvz*dvz)*(xyzmh_ptmass(4,p2)*pmassi)/(xyzmh_ptmass(4,p2)+pmassi)
    enddo qloopsink
    if (zero_one==1) then
       dx   = xyzh(1,p1) - xii
       dy   = xyzh(2,p1) - yii
       dz   = xyzh(3,p1) - zii
       hmax = max(xyzh(4,p1),hii)
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epot = epot - pmassii*pmassi/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epot  = epot + pmassii*pmassi*psoft*hmax1
       endif
       if (dr < radkern*xyzh(4,p1)) connected = .true. ! always just testing to see if candidate is a neighbour of a clump member
       dvx  = vxyzu(1,p1) - vxii
       dvy  = vxyzu(2,p1) - vyii
       dvz  = vxyzu(3,p1) - vzii
       ekin = ekin + 0.25*pmassii*(dvx*dvx + dvy*dvy + dvz*dvz)
    endif
 enddo jloop
!$omp enddo
!$omp end parallel
 epot  = epot + clump(j1)%potential + clump(j2)%potential   ! self-potential energy of the two clumps
 ekini = 0.
 epoti = 0.
!$omp parallel default(none) &
!$omp shared(j1,j2,n2,clump,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,idlistsink,connected) &
!$omp shared(xii,yii,zii,hii,vxii,vyii,vzii,pmassii,pmassi,use_npart2,zero_one) &
!$omp private(j,q,p1,p2,dx,dy,dz,dr,dvx,dvy,dvz,qp,psoft,fsoft,hmax,hmax1) &
!$omp reduction(+:epoti,ekini)
!$omp do schedule(runtime)
 qiloop: do q = 1,n2
    p2 = p_val(j2,q,use_npart2)
    if (p2 == 0) cycle qiloop
    jloopsink: do j = 1,clump(j1)%nptmass
       p1   = idlistsink(j1,j)
       dx   = xyzh(1,p2) - xyzmh_ptmass(1,p1)
       dy   = xyzh(2,p2) - xyzmh_ptmass(2,p1)
       dz   = xyzh(3,p2) - xyzmh_ptmass(3,p1)
       hmax = max(xyzh(4,p2),0.5*xyzmh_ptmass(ihacc,p1))
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epoti = epoti - pmassi*xyzmh_ptmass(4,p1)/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epoti = epoti + pmassi*xyzmh_ptmass(4,p1)*psoft*hmax1
       endif
       if (dr < max(radkern*xyzh(4,p2),xyzmh_ptmass(ihacc,p1)*racc_fac)) connected = .true.
       dvx   = vxyzu(1,p2) - vxyz_ptmass(1,p1)
       dvy   = vxyzu(2,p2) - vxyz_ptmass(2,p1)
       dvz   = vxyzu(3,p2) - vxyz_ptmass(3,p1)
       ekini = ekini + 0.5*(dvx*dvx + dvy*dvy + dvz*dvz)*(xyzmh_ptmass(4,p1)*pmassi)/(xyzmh_ptmass(4,p1)+pmassi)
    enddo jloopsink
    if (zero_one==1) then
       dx   = xyzh(1,p2) - xii
       dy   = xyzh(2,p2) - yii
       dz   = xyzh(3,p2) - zii
       hmax = max(xyzh(4,p2),hii)
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epoti = epoti - pmassii*pmassi/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epoti = epoti + pmassii*pmassi*psoft*hmax1
       endif
       if (dr < radkern*xyzh(4,p2)) connected = .true.
       dvx   = vxyzu(1,p2) - vxii
       dvy   = vxyzu(2,p2) - vyii
       dvz   = vxyzu(3,p2) - vzii
       ekini = ekini + 0.25*pmassii*(dvx*dvx + dvy*dvy + dvz*dvz)
    endif
 enddo qiloop
!$omp enddo
!$omp end parallel
 epot = epot + epoti
 ekin = ekin + ekini

 do j = 1,clump(j1)%nptmass
    p1   = idlistsink(j1,j)
    do q = 1,clump(j2)%nptmass
       p2   = idlistsink(j2,q)
       dx   = xyzmh_ptmass(1,p1) - xyzmh_ptmass(1,p2)
       dy   = xyzmh_ptmass(2,p1) - xyzmh_ptmass(2,p2)
       dz   = xyzmh_ptmass(3,p1) - xyzmh_ptmass(3,p2)
       hmax = max(0.5*xyzmh_ptmass(ihacc,p1),0.5*xyzmh_ptmass(ihacc,p2))
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epot = epot - xyzmh_ptmass(4,p1)*xyzmh_ptmass(4,p2)/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epot  = epot + xyzmh_ptmass(4,p1)*xyzmh_ptmass(4,p2)*psoft*hmax1
       endif
       if (dr < max(xyzmh_ptmass(ihacc,p1)*racc_fac,xyzmh_ptmass(ihacc,p2)*racc_fac)) connected = .true.
       dvx  = vxyz_ptmass(1,p1) - vxyz_ptmass(1,p2)
       dvy  = vxyz_ptmass(2,p1) - vxyz_ptmass(2,p2)
       dvz  = vxyz_ptmass(3,p1) - vxyz_ptmass(3,p2)
       ekin = ekin + 0.5*(dvx*dvx + dvy*dvy + dvz*dvz) &
                   *(xyzmh_ptmass(4,p1)*xyzmh_ptmass(4,p2))/(xyzmh_ptmass(4,p1)+xyzmh_ptmass(4,p2))
    enddo
    if (zero_one==1) then
       dx   = xyzmh_ptmass(1,p1) - xii
       dy   = xyzmh_ptmass(2,p1) - yii
       dz   = xyzmh_ptmass(3,p1) - zii
       hmax = max(0.5*xyzmh_ptmass(ihacc,p1),hii)
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epot = epot - pmassii*xyzmh_ptmass(4,p1)/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epot  = epot + pmassii*xyzmh_ptmass(4,p1)*psoft*hmax1
       endif
       if (dr < xyzmh_ptmass(ihacc,p1)*racc_fac) connected = .true.
       dvx  = vxyz_ptmass(1,p1) - vxii
       dvy  = vxyz_ptmass(2,p1) - vyii
       dvz  = vxyz_ptmass(3,p1) - vzii
       ekin = ekin + 0.5*(dvx*dvx + dvy*dvy + dvz*dvz)*(xyzmh_ptmass(4,p1)*pmassi)/(xyzmh_ptmass(4,p1)+pmassii)
    endif
 enddo
 if (zero_one==1) then
    do q = 1,clump(j2)%nptmass
       p2   = idlistsink(j2,q)
       dx   = xyzmh_ptmass(1,p2) - xii
       dy   = xyzmh_ptmass(2,p2) - yii
       dz   = xyzmh_ptmass(3,p2) - zii
       hmax = max(0.5*xyzmh_ptmass(ihacc,p2),hii)
       dr   = sqrt(dx*dx + dy*dy + dz*dz)
       if (dr > radkern*hmax .or. .not.soft_potential) then
          if (dr > 0.) epot = epot - pmassii*xyzmh_ptmass(4,p2)/dr
       else
          hmax1 = 1.0/hmax
          qp    = dr*hmax1
          call kernel_softening(qp*qp,qp,psoft,fsoft)
          epot  = epot + pmassii*xyzmh_ptmass(4,p2)*psoft*hmax1
       endif
       if (dr < xyzmh_ptmass(ihacc,p2)*racc_fac) connected = .true.
       dvx  = vxyz_ptmass(1,p2) - vxii
       dvy  = vxyz_ptmass(2,p2) - vyii
       dvz  = vxyz_ptmass(3,p2) - vzii
       ekin = ekin + 0.5*(dvx*dvx + dvy*dvy + dvz*dvz)*(xyzmh_ptmass(4,p2)*pmassi)/(xyzmh_ptmass(4,p2)+pmassii)
    enddo
 endif

 ! calculate separation
 sep = 0.
 if (zero_one==1) then
    dx  = xii - rtmp(1)
    dy  = yii - rtmp(3)
    dz  = zii - rtmp(3)
    dr  = dx*dx + dy*dy + dz*dz
    sep = max(sep,dr)
 endif

 sepi = 0.
!$omp parallel default(none) &
!$omp shared(j1,n1,xyzh,rtmp,use_npart1) &
!$omp private(j,p1,dx,dy,dz,dr) &
!$omp reduction(max:sepi)
!$omp do schedule(runtime)
 do j = 1,n1
    p1 = p_val(j1,j,use_npart1)
    if (p1 > 0) then
       dx   = xyzh(1,p1) - rtmp(1)
       dy   = xyzh(2,p1) - rtmp(3)
       dz   = xyzh(3,p1) - rtmp(3)
       dr   = dx*dx + dy*dy + dz*dz
       sepi = max(sepi,dr)
    endif
 enddo
!$omp enddo
!$omp end parallel
 sep = max(sep,sepi)

 sepi = 0.
!$omp parallel default(none) &
!$omp shared(j2,n2,xyzh,rtmp,use_npart2) &
!$omp private(q,p2,dx,dy,dz,dr) &
!$omp reduction(max:sepi)
!$omp do schedule(runtime)
 do q = 1,n2
    p2 = p_val(j2,q,use_npart2)
    if (p2 > 0) then
       dx   = xyzh(1,p2) - rtmp(1)
       dy   = xyzh(2,p2) - rtmp(2)
       dz   = xyzh(3,p2) - rtmp(3)
       dr   = dx*dx + dy*dy + dz*dz
       sepi = max(sepi,dr)
    endif
 enddo
!$omp enddo
!$omp end parallel
 sep = max(sep,sepi)

 do j = 1,clump(j1)%nptmass
    p1  = idlistsink(j1,j)
    dx  = xyzmh_ptmass(1,p1) - rtmp(1)
    dy  = xyzmh_ptmass(2,p1) - rtmp(2)
    dz  = xyzmh_ptmass(3,p1) - rtmp(3)
    dr  = dx*dx + dy*dy + dz*dz
    sep = max(sep,dr)
 enddo

 do q = 1,clump(j2)%nptmass
    p2  = idlistsink(j2,q)
    dx  = xyzmh_ptmass(1,p2) - rtmp(1)
    dy  = xyzmh_ptmass(2,p2) - rtmp(2)
    dz  = xyzmh_ptmass(3,p2) - rtmp(3)
    dr  = dx*dx + dy*dy + dz*dz
    sep = max(sep,dr)
 enddo

 ! This will ensure positive energy and not join clumps
 ! We want nearby clumps with sinks to be bound even if they are not overlapping since it is not always practical for a sink to overlap with gas
 if (zero_one==0 .and. .not. connected .and. max(clump(j1)%nptmass,clump(j2)%nptmass)==0) epot = abs(epot)

end subroutine are_clumps_joined
!--------------------------------------------------------------------------
! merge two clumps, and possibly the joining particle
subroutine merge_clumps(j1,j2,npart,nclump,vxyzu,Bxyz,rho,rtmp,vtmp,sep,pmassi,ekin,epot,i)
 integer, intent(in)           :: j1,j2,npart
 integer, intent(inout)        :: nclump
 integer, intent(in), optional :: i
 real,    intent(in)           :: vxyzu(:,:),Bxyz(:,:),rho(:),rtmp(:),vtmp(:),pmassi,sep,ekin,epot
 integer                       :: p,q,ii,zero_one,npartnew
 real                          :: pmassii,rhoii

 if (present(i)) then
    pmassii  = pmassi
    rhoii    = rho(i)
    ii       = i
    zero_one = 1
 else
    pmassii  = 0.
    rhoii    = 0.
    ii       = 1 ! just a dummy value, but will always be multiplied by zero
    zero_one = 0
 endif

 do q = 1,3
    clump(j1)%r(q)  = rtmp(q)
    clump(j1)%v(q)  = vtmp(q)
 enddo
 clump(j1)%mass      = clump(j1)%mass + clump(j2)%mass + pmassii
 clump(j1)%size      = sep
 clump(j1)%kinetic   = ekin
 clump(j1)%potential = epot
 clump(j1)%rhomax    = max(clump(j1)%rhomax,clump(j2)%rhomax,rhoii)
 idlistsink(j1,clump(j1)%nptmass+1:clump(j1)%nptmass+clump(j2)%nptmass) = idlistsink(j2,1:clump(j2)%nptmass)  ! particle i cannot be a sink
 if (present(i)) then
    ! add to clump j1
    idclump(i) = j1
    npartnew   = clump(j1)%nmember + 1
    clump(j1)%nmember = npartnew
    if (npartnew <= npartmax) then
       idlistpart(j1,npartnew) = i
    endif
 endif
 npartnew = clump(j1)%nmember + clump(j2)%nmember
 if (npartnew <= npartmax) then
    idlistpart(j1,clump(j1)%nmember+1:clump(j1)%nmember+clump(j2)%nmember) = idlistpart(j2,1:clump(j2)%nmember)
 endif
 clump(j1)%nmember   = npartnew
 clump(j1)%nptmass   = clump(j1)%nptmass + clump(j2)%nptmass
 if (maxvxyzu>=4) clump(j1)%thermal  = clump(j1)%thermal  + clump(j2)%thermal  + vxyzu(4,ii)*pmassii
 if (mhd)         clump(j1)%magnetic = clump(j1)%magnetic + clump(j2)%magnetic &
                                                          + 0.5*pmassii*dot_product(Bxyz(1:3,ii),Bxyz(1:3,ii))/rho(ii)
 do p = 1,npart
    if (idclump(p)==j2) idclump(p) = j1
 enddo
 do p = 1,npart
    if (idclump(p)==nclump) idclump(p) = j2
 enddo
 clump(j2)        = clump(nclump)
 clump(j2)%ID     = j2
 idlistpart(j2,:) = idlistpart(nclump,:)
 idlistsink(j2,:) = idlistsink(nclump,:)
 nclump           = nclump - 1

end subroutine merge_clumps
!--------------------------------------------------------------------------
! Determine the loop based upon the number of members
subroutine n_in_clump(j,npart,nout,use_npart)
 integer, intent(in)  :: j,npart
 integer, intent(out) :: nout
 logical, intent(out) :: use_npart

 if (clump(j)%nmember > npartmax) then
    nout      = npart
    use_npart = .true.
 else
    nout      = clump(j)%nmember
    use_npart = .false.
 endif

end subroutine n_in_clump

!--------------------------------------------------------------------------
! Get the particle index; return 0 if we need to skip
integer function p_val(jclump,jpart,use_npart)
 integer, intent(in) :: jclump,jpart
 logical, intent(in) :: use_npart

 if (use_npart) then
    p_val = jpart
    if (idclump(jpart) /= jclump) p_val = 0
 else
    p_val = idlistpart(jclump,jpart)
 endif

end function p_val
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The next set of subroutine are required for ellipse fitting.
! First, we define a grid, then project the particles onto the grid
! Next, we use the grid to determine the ellipses, and rotate until the
! ellipses lies along the Cartesian axes.  This is necessary, since using
! the SPH particles directly is effectively a density-weighted fit and
! represents the core rather than the full cloud
!----------------------------------------------------------------------
! Define the grid for ellipse fitting
subroutine set_grid(n,i1,i2,xyz,hsmooth,dx0,x1grid,x2grid)
 integer, intent(in)  :: n,i1,i2
 real,    intent(in)  :: xyz(:,:),hsmooth(:)
 real,    intent(out) :: dx0,x1grid(ngrid),x2grid(ngrid)
 integer              :: i,nx1,nx2
 real                 :: xmin(3),xmax(3)

 ! find extreme values
 xmin =  huge(xmin(1))
 xmax = -huge(xmax(1))
 do i = 1,n
    xmin(1) = min(xmin(1),xyz(i1,i))
    xmax(1) = max(xmax(1),xyz(i1,i))
    xmin(2) = min(xmin(2),xyz(i2,i))
    xmax(2) = max(xmax(2),xyz(i2,i))
    xmin(3) = min(xmin(3),hsmooth(i))
    xmax(3) = max(xmax(3),hsmooth(i))
 enddo

 ! Verify the grid size is large enough to encompass the cloud; otherwise reset dx
 nx1 = int((xmax(1) - xmin(1) + 4.*xmax(3))/dxgrid0) + 2
 nx2 = int((xmax(2) - xmin(2) + 4.*xmax(3))/dxgrid0) + 2
 if (nx1 > ngrid .or. nx2 > ngrid) then
    if (nx1 > nx2) then
       dx0 = (xmax(1) - xmin(1) + 4.*xmax(3))/(ngrid*0.9)
    else
       dx0 = (xmax(2) - xmin(2) + 4.*xmax(3))/(ngrid*0.9)
    endif
    print*, 'Cloud too big for grid.  Resetting dx = ',dx0, nx1,nx2,ngrid
 else
    dx0 = dxgrid0
 endif

 ! define the grid
 x1grid(ngrid/2) = 0.5*(xmax(1) + xmin(1))
 x2grid(ngrid/2) = 0.5*(xmax(2) + xmin(2))
 do i = ngrid/2+1,ngrid
    x1grid(i) = x1grid(i-1) + dx0
    x2grid(i) = x2grid(i-1) + dx0
 enddo
 do i = ngrid/2-1,1,-1
    x1grid(i) = x1grid(i+1) - dx0
    x2grid(i) = x2grid(i+1) - dx0
 enddo

end subroutine set_grid
!----------------------------------------------------------------------
! Create a binary grid where the entry is one if a particle overlaps a cell
subroutine fill_igrid(n,i1,i2,dx,xyz,h,x1grid,x2grid,igrid)
 integer, intent(in)  :: n,i1,i2
 integer, intent(out) :: igrid(:,:)
 real,    intent(in)  :: dx,xyz(:,:),h(:),x1grid(:),x2grid(:)
 integer              :: i,j,p,imin,imax,jmin,jmax,hondx
 real                 :: h2,r2,x10,x20

 igrid = 0
 do p = 1,n
    i = 1
    do while (xyz(i1,p) > x1grid(i) .and. i < ngrid)
       i = i + 1
    enddo
    j = 1
    do while (xyz(i2,p) > x2grid(j) .and. j < ngrid)
       j = j + 1
    enddo
    if (i==1 .or. i==ngrid) print*, 'WARNING! SPH PARICLE CENTRE IS OUTSIDE OF THE GRID', i,x1grid(1),xyz(i1,p),x1grid(ngrid)
    if (j==1 .or. j==ngrid) print*, 'WARNING! SPH PARICLE CENTRE IS OUTSIDE OF THE GRID', j,x2grid(1),xyz(i2,p),x2grid(ngrid)
    igrid(i,j) = 1
    ! tag the surrounding cells to account for smoothing length
    hondx = int(h(p)/dx)
    imin = max(1,    i-hondx - 1)
    imax = min(ngrid,i+hondx + 1)
    jmin = max(1,    j-hondx - 1)
    jmax = min(ngrid,j+hondx + 1)
    x10  = x1grid(i)
    x20  = x2grid(j)
    h2   = h(p)*h(p)
    do i = imin,imax
       do j = jmin,jmax
          r2 = (x1grid(i)-x10)**2 + (x2grid(j)-x20)**2
          if (h2 > r2) igrid(i,j) = 1
       enddo
    enddo
 enddo

end subroutine fill_igrid

!----------------------------------------------------------------------
subroutine rotate_data(n,i1,i2,xc,yc,theta,xyz)
 integer, intent(in)    :: n,i1,i2
 real,    intent(in)    :: xc,yc,theta
 real,    intent(inout) :: xyz(:,:)
 integer                :: i
 real                   :: x1,x2

 do i = 1,n
    x1 = xyz(i1,i) - xc
    x2 = xyz(i2,i) - yc
    xyz(i1,i) =  x1*cos(theta) + x2*sin(theta) + xc
    xyz(i2,i) = -x1*sin(theta) + x2*cos(theta) + yc
 enddo

end subroutine rotate_data
!--------------------------------------------------------------------------
! 2D ellipse fitting routine of Rocha,Velho & Carvalho (2002)
subroutine calc_moment(x1,x2,igrid,l,w,phi,xc,yc)
 real,    intent(in)  :: x1(:),x2(:)
 integer, intent(in)  :: igrid(:,:)
 real,    intent(out) :: l,w,phi,xc,yc
 integer              :: i,j
 real                 :: m00,m01,m10,m11,m20,m02,a,b,c

 ! the moments
 m00 = 0.
 m01 = 0.
 m10 = 0.
 m11 = 0.
 m02 = 0.
 m20 = 0.
 do i = 1,ngrid
    do j = 1,ngrid
       m00 = m00 +             igrid(i,j)
       m10 = m10 + x1(i)      *igrid(i,j)
       m01 = m01 + x2(j)      *igrid(i,j)
       m11 = m11 + x1(i)*x2(j)*igrid(i,j)
       m20 = m20 + x1(i)**2   *igrid(i,j)
       m02 = m02 + x2(j)**2   *igrid(i,j)
    enddo
 enddo

 ! the centre of the ellipse
 xc = m10/m00
 yc = m01/m00

 ! intermediate calculations
 a = m20/m00 - xc**2
 b = 2.0*(m11/m00 - xc*yc)
 c = m02/m00 - yc**2

 ! the results: rotation angle, minor axis & major axis
 phi = 0.5*atan2(b,a-c)
 w   = sqrt(6.*(a+c - sqrt(b*b+(a-c)**2)))
 l   = sqrt(6.*(a+c + sqrt(b*b+(a-c)**2)))

end subroutine calc_moment
!--------------------------------------------------------------------------
! To time the subroutine
subroutine update_time(walltime)
 use timing,     only:wallclock
 real, intent(inout) :: walltime

 walltime = walltime + wallclock()
 if (walltime < 120.) then
    write(*,'(a,F8.3,a)') 'The elapsed time is ',walltime,'s'
 elseif (walltime < 7200.) then
    write(*,'(a,F8.3,a)') 'The elapsed time is ',walltime/60.,'m'
 else
    write(*,'(a,F8.3,a)') 'The elapsed time is ',walltime/3600.,'h'
 endif

end subroutine update_time
!--------------------------------------------------------------------------
! Write information to params file
subroutine write_clumpparams(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 print "(a)",' writing analysis options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a,/)") '# options when performing clumpfinding analysis'
 call write_inopt(nclumpmax,     'nclumpmax', 'the maximum number of clumps (recommend 100000 on clusters)',iunit)
 call write_inopt(npartmax,      'npartmax',  'maximum number of particles per clump (recommend 100000 on clusters)',iunit)
 call write_inopt(rhominpeak_cgs,'rhominpeak','density (cgs) above which a particle can be the lead particle in a clump',iunit)
 call write_inopt(rhominbkg_cgs, 'rhominbkg', 'density (cgs) above which a particle can be a member of a clump',iunit)
 call write_inopt(ekin_coef,     'ekin_coef', 'a clump will be bound if coef*ekin + epot < 0',iunit)
 call write_inopt(dxgrid0_pc,    'dx_grid',   'grid resolution to determine the clump shape (pc)',iunit)

 close(iunit)

end subroutine write_clumpparams
!--------------------------------------------------------------------------
! Read information from params file
subroutine read_clumpparams(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer,          parameter   :: iunit = 21
 integer                       :: nerr
 type(inopts),     allocatable :: db(:)

 print "(a)",'reading analysis options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(nclumpmax,      'nclumpmax', db,errcount=nerr)
 call read_inopt(npartmax,       'npartmax',  db,errcount=nerr)
 call read_inopt(rhominpeak_cgs, 'rhominpeak',db,errcount=nerr)
 call read_inopt(rhominbkg_cgs,  'rhominbkg', db,errcount=nerr)
 call read_inopt(ekin_coef,      'ekin_coef', db,errcount=nerr)
 call read_inopt(dxgrid0_pc,     'dx_grid',   db,errcount=nerr)

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of params file: re-writing...'
    ierr = nerr
 endif

end subroutine read_clumpparams
!--------------------------------------------------------------------------
end module analysis
