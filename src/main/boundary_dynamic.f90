!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module boundary_dyn
!
! This module contains variables and subroutines relating to boundaries,
! including dynamically adjusting periodic boundaries
!
! :References: Wurster & Bonnell (2023) MNRAS, 522, 891-911 for dynamic boundaries
!
! :Owner: James Wurster
!
! :Runtime parameters:
!   - dynamic_bdy    : *turn on/off dynamic boundaries*
!   - n_dtmax        : *particles must not reach v*n_dtmax*dtmax of the boundary*
!   - rho_thresh_bdy : *threshold density separating dense gas from background gas*
!   - vbdyx          : *velocity of the x-boundary*
!   - vbdyy          : *velocity of the y-boundary*
!   - vbdyz          : *velocity of the z-boundary*
!   - width_bkg_nx   : *width of the boundary in the -x direction*
!   - width_bkg_ny   : *width of the boundary in the -y direction*
!   - width_bkg_nz   : *width of the boundary in the -z direction*
!   - width_bkg_px   : *width of the boundary in the +x direction*
!   - width_bkg_py   : *width of the boundary in the +y direction*
!   - width_bkg_pz   : *width of the boundary in the +z direction*
!
! :Dependencies: boundary, dim, infile_utils, io, kernel, mpidomain, part
!

 use dim, only: maxvxyzu
 use io,  only: fatal
 use boundary, only: xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound,totvol,cross_boundary
 implicit none

 logical, public :: dynamic_bdy    = .false.
 integer, public :: irho_bkg_ini   = 0
 real,    public :: rho_thresh_bdy = 0.
 real,    public :: rho_bkg_ini    = 0.
 real,    public :: n_dtmax        = 3.         ! interesting particles must be > ndtmax*dtmax*v away from the boundary
 real,    public :: dxyz           = 0.
 real,    public :: vbdyx          = 0.
 real,    public :: vbdyy          = 0.
 real,    public :: vbdyz          = 0.
 real,    public :: width_bkg(3,2) = 0.         ! (:,1) == min side; (:,2) == max side
 real,    public :: xyz_n(3),xyz_x(3)           ! extreme positions of the dense particles
 real,    public :: vdtxyz_n(3),vdtxyz_x(3)     ! extreme positions of where particles will end up
 real,    public :: v_bkg(maxvxyzu),B_bkg(3)
 real,    public :: rho_bkg_ini1   = 0.
 integer, public :: ibkg,n_bkg
 logical, public :: vbdyx_not_0 = .false.
 logical, public :: vbdyy_not_0 = .false.
 logical, public :: vbdyz_not_0 = .false.
 logical         :: remove_accreted = .true.      ! remove accreted particles from the list
 real,    public, parameter :: bdy_vthresh = 0.05 ! tolerance on velocity
 ! if a particle's velocity is within this tolerance, then it is a
 ! 'boring' particle that may be used to determine the background.
 ! Splash dumps suggest 1% is too low and 10% is too high.

 public :: adjust_particles_dynamic_boundary,find_dynamic_boundaries,update_xyzminmax
 public :: in_domain,update_boundaries,init_dynamic_bdy
 public :: write_options_boundary,read_options_boundary

 private

contains

!----------------------------------------------------------------
!+
!  Initialise dynamic boundaries
!  This must be done twice, once before derivs to ensure
!  boundaries are in the correct spot, and once after to ensure
!  the paramters are correct.
!+
!---------------------------------------------------------------
subroutine init_dynamic_bdy(icall,npart,nptmass,dtmax)
 use part, only: xyzh,rhoh,iorig,massoftype,igas
 integer, intent(in)    :: icall,nptmass
 integer, intent(inout) :: npart
 real,    intent(in)    :: dtmax
 real    :: xyz_n_all(3),xyz_x_all(3)
 integer :: ndummy1,ndummy2,ierr
 logical :: abortrun

 if (icall==1) then
    ! Update the background medium, if required.  Do this prior
    ! to calling derivs so the new particles are properly initialised
    call set_dynamic_bdy_width()
    call find_dynamic_boundaries(npart,nptmass,dtmax,xyz_n_all,xyz_x_all,ierr)
    call update_boundaries(ndummy1,ndummy2,npart,abortrun)
 elseif (icall==2) then
    ! Reset the background density for dynamic boundaries, if necessary;
    ! this is to ensure a consistent density
    if (irho_bkg_ini > 0.) then
       print*, 'original rho_bkg_ini = ',rho_bkg_ini
       rho_bkg_ini  = rhoh(xyzh(4,iorig(irho_bkg_ini)),massoftype(igas))
       rho_bkg_ini1 = 1.0/rho_bkg_ini
       print*, 'revised rho_bkg_ini = ',rho_bkg_ini
       irho_bkg_ini = 0
    endif
    call set_dynamic_bdy_width()
 endif

end subroutine init_dynamic_bdy

!----------------------------------------------------------------
!+
!  Determine if the boundary is moving and in what direction
!+
!---------------------------------------------------------------
subroutine set_dynamic_bdy_width()
 use part,   only: massoftype,igas,hfact
 use kernel, only: radkern
 integer :: i,j
 real    :: minborder

 if (abs(vbdyx) > epsilon(vbdyx)) vbdyx_not_0 = .true.
 if (abs(vbdyy) > epsilon(vbdyy)) vbdyy_not_0 = .true.
 if (abs(vbdyz) > epsilon(vbdyz)) vbdyz_not_0 = .true.

 minborder = 10.0*radkern*hfact*(massoftype(igas)*rho_bkg_ini1)**(1.0/3.0)
 do j = 1,2
    do i = 1,3
       width_bkg(i,j) = max(width_bkg(i,j),minborder)
    enddo
 enddo

end subroutine set_dynamic_bdy_width

!----------------------------------------------------------------
!+
!  Adjust particles to ensure that they are in the domain
!+
!---------------------------------------------------------------
subroutine adjust_particles_dynamic_boundary(npart,xyzh)
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:)
 integer                :: i

 do i=1,npart
    if (xyzh(1,i) < xmin) xyzh(1,i) = xyzh(1,i) + dxbound
    if (xyzh(1,i) > xmax) xyzh(1,i) = xyzh(1,i) - dxbound
    if (xyzh(2,i) < ymin) xyzh(2,i) = xyzh(2,i) + dybound
    if (xyzh(2,i) > ymax) xyzh(2,i) = xyzh(2,i) - dybound
    if (xyzh(3,i) < zmin) xyzh(3,i) = xyzh(3,i) + dzbound
    if (xyzh(3,i) > zmax) xyzh(3,i) = xyzh(3,i) - dzbound
 enddo

end subroutine adjust_particles_dynamic_boundary

!----------------------------------------------------------------
!+
!  For dynamic boundaries, will determine if the particle is within
!  the gas domain
!+
!---------------------------------------------------------------
subroutine update_xyzminmax(dt)
 real, intent(in) :: dt

 xmin = xmin + dt*vbdyx
 xmax = xmax + dt*vbdyx
 ymin = ymin + dt*vbdyy
 ymax = ymax + dt*vbdyy
 zmin = zmin + dt*vbdyz
 zmax = zmax + dt*vbdyz

end subroutine update_xyzminmax

!----------------------------------------------------------------
!+
!  For dynamic boundaries, will determine if the particle is within
!  the gas domain
!+
!---------------------------------------------------------------
logical function in_domain(xyz,domain)
 real, intent(in) :: xyz(3),domain(3,2)

 if ( xyz(1) < domain(1,1) .or. xyz(1) > domain(1,2) .or. &
      xyz(2) < domain(2,1) .or. xyz(2) > domain(2,2) .or. &
      xyz(3) < domain(3,1) .or. xyz(3) > domain(3,2) ) then
    in_domain = .false.
 else
    in_domain = .true.
 endif

end function in_domain

!----------------------------------------------------------------
!+
!  This will calculate the location of the dynamic boundaries
!+
!---------------------------------------------------------------
subroutine find_dynamic_boundaries(npart,nptmass,dtmax,xyz_n_all,xyz_x_all,ierr)
 use io,     only:id,master
 use part,   only: maxp,maxphase,mhd,massoftype,igas,ics,isdead_or_accreted,rhoh,iamtype
 use part,   only: xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,Bevol,eos_vars,iphase
 use kernel, only: radkern
 integer, intent(in)  :: npart,nptmass
 integer, intent(out) :: ierr
 real,    intent(in)  :: dtmax
 real,    intent(out) :: xyz_n_all(3),xyz_x_all(3)
 integer              :: i,itype,ibkg_thread,n_bkg
 real                 :: xi,yi,zi,hi,pmassi,rhoi,rho1i,vxi,vyi,vzi,v2i,vi1,Bxi,Byi,Bzi,B2i,valfven2i,spsoundi,spsound2i
 real                 :: fourh,ndtmax,vmsxi,vmsyi,vmszi,va2cs2,rho1cs2,rho_bkg_thread,rho_bkg,v_bkg(3),B_bkg(3)
 logical              :: high_density_gas,bdy_is_interesting

 itype     = igas
 pmassi    = massoftype(igas)
 ierr      = 0
 xyz_x     = -huge(xyz_x(1))
 xyz_n     =  huge(xyz_n(1))
 vdtxyz_x  = -huge(vdtxyz_x(1))
 vdtxyz_n  =  huge(vdtxyz_n(1))
 xyz_x_all = -huge(xyz_x_all(1))
 xyz_n_all =  huge(xyz_n_all(1))
 rho_bkg   = huge(rho_bkg)
 v_bkg     = 0.
 B_bkg     = 0.
 n_bkg     = 0
 ndtmax    = n_dtmax*dtmax
 high_density_gas = .false.

!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(npart,xyzh,vxyzu,Bevol,eos_vars,iphase,massoftype,id) &
!$omp shared(nptmass,xyzmh_ptmass,vxyz_ptmass) &
!$omp shared(rho_thresh_bdy,rho_bkg,ibkg,high_density_gas,rho_bkg_ini1) &
!$omp shared(xmin,xmax,ymin,zmin,ymax,zmax) &
!$omp shared(vbdyx,vbdyy,vbdyz,vbdyx_not_0,vbdyy_not_0,vbdyz_not_0,ndtmax) &
!$omp private(i,xi,yi,zi,hi,rhoi,rho1i,vxi,vyi,vzi,Bxi,Byi,Bzi,B2i,v2i,vi1) &
!$omp private(spsoundi,spsound2i,va2cs2,rho1cs2,valfven2i,fourh) &
!$omp private(rho_bkg_thread,ibkg_thread,vmsxi,vmsyi,vmszi,bdy_is_interesting) &
!$omp firstprivate(itype,pmassi) &
!$omp reduction(min:xyz_n,xyz_n_all,vdtxyz_n) &
!$omp reduction(max:xyz_x,xyz_x_all,vdtxyz_x) &
!$omp reduction(+:v_bkg,B_bkg,n_bkg)
 ibkg_thread    = 0
 rho_bkg_thread = huge(rho_bkg_thread)
!$omp do
 do i = 1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (itype <= 0) call fatal('energies','particle type <= 0')
          pmassi = massoftype(itype)
       endif
       rhoi = rhoh(hi,pmassi)

       ! determine the particle whose density is closest to the original background
       ! all new particles will initially be copied from this particle
       if ( abs(rhoi*rho_bkg_ini1 - 1.0) < rho_bkg_thread) then
          rho_bkg_thread = abs(rhoi*rho_bkg_ini1 - 1.0)
          ibkg_thread    = i
       endif
       ! The actual boundary of the particles
       xyz_n_all(1) = min(xyz_n_all(1),xi)
       xyz_n_all(2) = min(xyz_n_all(2),yi)
       xyz_n_all(3) = min(xyz_n_all(3),zi)
       xyz_x_all(1) = max(xyz_x_all(1),xi)
       xyz_x_all(2) = max(xyz_x_all(2),yi)
       xyz_x_all(3) = max(xyz_x_all(3),zi)
       ! The boundary of the dense particles
       if (rhoi > rho_thresh_bdy) then
          xyz_n(1) = min(xyz_n(1),xi)
          xyz_n(2) = min(xyz_n(2),yi)
          xyz_n(3) = min(xyz_n(3),zi)
          xyz_x(1) = max(xyz_x(1),xi)
          xyz_x(2) = max(xyz_x(2),yi)
          xyz_x(3) = max(xyz_x(3),zi)
          high_density_gas = .true. ! technically a race condition since shared, but can only become true, and want true if even only particle is true
       endif

       ! The boundary of the particles satisfying the velocity requirement
       bdy_is_interesting = .false.
       vxi = vxyzu(1,i)
       vzi = vxyzu(2,i)
       vyi = vxyzu(3,i)
       v2i = vxi*vxi + vyi*vyi + vzi*vzi
       vi1   = 1.0/sqrt(v2i)
       fourh = 2.0*radkern*hi
       ! x-velocity
       if (vbdyx_not_0) then
          ! tag if differs by more than bdy_vthresh from the preset background velocity
          if (abs(1.0 - abs(vbdyx)*vi1) > bdy_vthresh) bdy_is_interesting = .true.
       else
          ! tag if velocity component is more than bdy_vthresh of the total velocity
          if (abs(vxi)*vi1 > bdy_vthresh) bdy_is_interesting = .true.
       endif
       ! y-velocity
       if (vbdyy_not_0) then
          if (abs(1.0 - abs(vbdyy)*vi1) > bdy_vthresh) bdy_is_interesting = .true.
       else
          if (abs(vyi)*vi1 > bdy_vthresh) bdy_is_interesting = .true.
       endif
       ! z-velocity
       if (vbdyz_not_0) then
          if (abs(1.0 - abs(vbdyz)*vi1) > bdy_vthresh) bdy_is_interesting = .true.
       else
          if (abs(vzi)*vi1 > bdy_vthresh) bdy_is_interesting = .true.
       endif

       ! particle satisfies velocity criteria, and we need to ensure that it will not promptly reach the boundary
       if (bdy_is_interesting) then
          ! ignore particles within 4h of the boundary since they may be affected by non-periodic gravity
          if (xi-fourh < xmin .or. xi+fourh > xmax) bdy_is_interesting = .false.
          if (yi-fourh < ymin .or. yi+fourh > ymax) bdy_is_interesting = .false.
          if (zi-fourh < zmin .or. zi+fourh > zmax) bdy_is_interesting = .false.
          if (bdy_is_interesting) then
             spsoundi   = eos_vars(ics,i)
             ! Ensure interesting particles will stay v*ndtmax away from the boundary
             if (mhd .and. itype==igas) then
                Bxi       = Bevol(1,i)*rhoi
                Byi       = Bevol(2,i)*rhoi
                Bzi       = Bevol(3,i)*rhoi
                B2i       = Bxi*Bxi + Byi*Byi + Bzi*Bzi
                rho1i     = 1./rhoi
                valfven2i = B2i*rho1i

                ! velocity dependent boundaries, assuming MHD; use v \pm fast magnetosonic wave
                spsound2i  = spsoundi*spsoundi
                va2cs2     = valfven2i+spsound2i
                rho1cs2    = 4.0*spsound2i*rho1i

                vmsxi = sqrt(0.5*(va2cs2 + sqrt(va2cs2**2 - Bxi*Bxi*rho1cs2)))
                vdtxyz_n(1) = min(vdtxyz_n(1),xi+ndtmax*(vxi+vmsxi-vbdyx),xi+ndtmax*(vxi-vmsxi-vbdyx))
                vdtxyz_x(1) = max(vdtxyz_x(1),xi+ndtmax*(vxi+vmsxi-vbdyx),xi+ndtmax*(vxi-vmsxi-vbdyx))

                vmsyi = sqrt(0.5*(va2cs2 + sqrt(va2cs2**2 - Byi*Byi*rho1cs2)))
                vdtxyz_n(2) = min(vdtxyz_n(2),yi+ndtmax*(vyi+vmsyi-vbdyy),yi+ndtmax*(vyi-vmsyi-vbdyy))
                vdtxyz_x(2) = max(vdtxyz_x(2),yi+ndtmax*(vyi+vmsyi-vbdyy),yi+ndtmax*(vyi-vmsyi-vbdyy))

                vmszi = sqrt(0.5*(va2cs2 + sqrt(va2cs2**2 - Bzi*Bzi*rho1cs2)))
                vdtxyz_n(3) = min(vdtxyz_n(3),zi+ndtmax*(vzi+vmszi-vbdyz),zi+ndtmax*(vzi-vmszi-vbdyz))
                vdtxyz_x(3) = max(vdtxyz_x(3),zi+ndtmax*(vzi+vmszi-vbdyz),zi+ndtmax*(vzi-vmszi-vbdyz))
             else
                ! velocity dependent boundaries, assuming pure hydrodynamics
                vdtxyz_n(1) = min(vdtxyz_n(1),xi+ndtmax*(vxi-spsoundi-vbdyx))
                vdtxyz_n(2) = min(vdtxyz_n(2),yi+ndtmax*(vyi-spsoundi-vbdyy))
                vdtxyz_n(3) = min(vdtxyz_n(3),zi+ndtmax*(vzi-spsoundi-vbdyz))
                vdtxyz_x(1) = max(vdtxyz_x(1),xi+ndtmax*(vxi+spsoundi-vbdyx))
                vdtxyz_x(2) = max(vdtxyz_x(2),yi+ndtmax*(vyi+spsoundi-vbdyy))
                vdtxyz_x(3) = max(vdtxyz_x(3),zi+ndtmax*(vzi+spsoundi-vbdyz))
             endif
          endif
       endif
       ! find the average velocity & magnetic field strength of the background particles not near the boundaries
       if (.not.bdy_is_interesting) then
          ! remove particles near the boundaries
          if (xi-fourh < xmin .or. xi+fourh > xmax) bdy_is_interesting = .true.
          if (yi-fourh < ymin .or. yi+fourh > ymax) bdy_is_interesting = .true.
          if (zi-fourh < zmin .or. zi+fourh > zmax) bdy_is_interesting = .true.
          ! ensure we're in a narrow density range
          if ( abs(rhoi*rho_bkg_ini1 - 1.0) > bdy_vthresh) bdy_is_interesting = .true.
          ! add uninteresting particles to the averages
          if (.not.bdy_is_interesting) then
             n_bkg = n_bkg + 1
             v_bkg = v_bkg + vxyzu(1:3,i)
             if (mhd) B_bkg = B_bkg + Bevol(:,i)*rhoi
          endif
       endif
    endif
 enddo
 if (id==master) then
    !$omp do
    do i=1,nptmass
       xi     = xyzmh_ptmass(1,i)
       yi     = xyzmh_ptmass(2,i)
       zi     = xyzmh_ptmass(3,i)
       pmassi = xyzmh_ptmass(4,i)
       if (pmassi < 0.) cycle

       vxi    = vxyz_ptmass(1,i)
       vyi    = vxyz_ptmass(2,i)
       vzi    = vxyz_ptmass(3,i)

       ! boundary of all particles or a sink
       xyz_n_all(1) = min(xyz_n_all(1),xi)
       xyz_n_all(2) = min(xyz_n_all(2),yi)
       xyz_n_all(3) = min(xyz_n_all(3),zi)
       xyz_x_all(1) = max(xyz_x_all(1),xi)
       xyz_x_all(2) = max(xyz_x_all(2),yi)
       xyz_x_all(3) = max(xyz_x_all(3),zi)
       ! boundary of dense particles or a sink
       xyz_n(1) = min(xyz_n(1),xi)
       xyz_n(2) = min(xyz_n(2),yi)
       xyz_n(3) = min(xyz_n(3),zi)
       xyz_x(1) = max(xyz_x(1),xi)
       xyz_x(2) = max(xyz_x(2),yi)
       xyz_x(3) = max(xyz_x(3),zi)
       high_density_gas = .true. ! or a sink particle
       ! boundary of sinks when accounting for sink velocity
       vdtxyz_n(1) = min(vdtxyz_n(1),xi+ndtmax*(vxi-vbdyx))
       vdtxyz_n(2) = min(vdtxyz_n(2),yi+ndtmax*(vyi-vbdyy))
       vdtxyz_n(3) = min(vdtxyz_n(3),zi+ndtmax*(vzi-vbdyz))
       vdtxyz_x(1) = max(vdtxyz_x(1),xi+ndtmax*(vxi-vbdyx))
       vdtxyz_x(2) = max(vdtxyz_x(2),yi+ndtmax*(vyi-vbdyy))
       vdtxyz_x(3) = max(vdtxyz_x(3),zi+ndtmax*(vzi-vbdyz))
    enddo
    !$omp enddo
 endif

!$omp critical(collatedata)
 if (rho_bkg_thread < rho_bkg) then
    rho_bkg = rho_bkg_thread
    ibkg    = ibkg_thread
 endif
!$omp end critical(collatedata)
!$omp end parallel

 if (n_bkg > 1) then
    v_bkg = v_bkg/n_bkg
    B_bkg = B_bkg/n_bkg
 endif

 if (.not. high_density_gas) ierr = 1

end subroutine find_dynamic_boundaries

!----------------------------------------------------------------
!+
!  Update the dynamic boundaries
!+
!---------------------------------------------------------------
subroutine update_boundaries(nactive,nalive,npart,abortrun_bdy)
 use dim,       only: maxp_hard,mhd
 use mpidomain, only: isperiodic
 use io,        only: iprint
 use part,      only: set_particle_type,copy_particle_all,shuffle_part,kill_particle,&
                      isdead_or_accreted,npartoftype,xyzh,igas,vxyzu,Bxyz
 integer, intent(inout) :: nactive,nalive,npart
 logical, intent(out)   :: abortrun_bdy
 integer                :: i,npart0,ndie,nadd,ncross,naccreted
 real                   :: dx,dy,dz
 real                   :: border(3,2)
 logical                :: updated_bdy

 !--Initialise variables
 npart0       = npart
 nadd         = 0
 ndie         = 0
 naccreted    = 0
 ncross       = 0
 abortrun_bdy = .false.
 updated_bdy  = .false.

 !--Dynamically calculate the border based upon how far a 'fast' particle can move within v*n_dtmax*dtmax (calculated in energies)
 border(:,1) = min(vdtxyz_n(:),xyz_n(:) - width_bkg(:,1))  ! -x, -y, -z boundaries
 border(:,2) = max(vdtxyz_x(:),xyz_x(:) + width_bkg(:,2))  ! +x, +y, +z boundaries

 if (.false.) then
    !--Determine if we need to expand the buffer zone using the location of dense clumps if using density threshold only
    !  this should be obsolete given the velocity criteria
    if (xyz_n(1) - xmin < 0.5*width_bkg(1,1)) then
       width_bkg(1,1) = 2.0*width_bkg(1,1)
       write(iprint,*) 'Updating Boundaries: increasing -x background medium width to ',width_bkg(1,1)
    endif
    if (xyz_n(2) - ymin < 0.5*width_bkg(2,1)) then
       width_bkg(2,1) = 2.0*width_bkg(2,1)
       write(iprint,*) 'Updating Boundaries: increasing -y background medium width to ',width_bkg(2,1)
    endif
    if (xyz_n(3) - zmin < 0.5*width_bkg(3,1)) then
       width_bkg(3,1) = 2.0*width_bkg(3,1)
       write(iprint,*) 'Updating Boundaries: increasing -z background medium width to ',width_bkg(3,1)
    endif
    if (xmax - xyz_x(1) < 0.5*width_bkg(1,2)) then
       width_bkg(1,2) = 2.0*width_bkg(1,2)
       write(iprint,*) 'Updating Boundaries: increasing +x background medium width to ',width_bkg(1,2)
    endif
    if (ymax - xyz_x(2) < 0.5*width_bkg(2,2)) then
       width_bkg(2,2) = 2.0*width_bkg(2,2)
       write(iprint,*) 'Updating Boundaries: increasing +y background medium width to ',width_bkg(1,2)
    endif
    if (zmax - xyz_x(3) < 0.5*width_bkg(3,2)) then
       width_bkg(3,2) = 2.0*width_bkg(3,2)
       write(iprint,*) 'Updating Boundaries: increasing +z background medium width to ',width_bkg(3,2)
    endif
    !--Reset the boundaries
    border(:,1) = xyz_n(:) - width_bkg(:,1)  ! -x, -y, -z boundaries
    border(:,2) = xyz_x(:) + width_bkg(:,2)  ! +x, +y, +z boundaries
 endif

 !--Ensure that there are no particles outside the boundaries
 !$omp parallel default(none) &
 !$omp shared(npart,xyzh) &
 !$omp private(i) &
 !$omp shared(isperiodic) &
 !$omp reduction(+:ncross)
 !$omp do
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call cross_boundary(isperiodic,xyzh(:,i),ncross)
    endif
 enddo
 !$omp enddo
 !$omp end parallel
 if (ncross > 0) print*, 'There were ',ncross,' particles on the wrong side of the boundary'

 !--Add particles on faces to which the interesting particles are advancing
 !  Must be serial to avoid conflicts when adding new particles

 ! add particles in the -x direction
 do while (xmin - dxyz > border(1,1))
    updated_bdy = .true.
    call write_bdy_update(iprint,'increasing','-x','xmin',xmin,xmin-dxyz)
    !write(iprint,*) 'Updating Boundaries: increasing -x: xmin_old, xmin_new = ',xmin,xmin-dxyz
    xmin = xmin - dxyz
    dx = xmin + 0.5*dxyz
    dz = zmin + 0.5*dxyz
    do while (dz < zmax)
       dy = ymin + 0.5*dxyz
       do while (dy < ymax)
          npart = min(npart + 1,maxp_hard)
          call copy_particle_all(ibkg,npart,.true.)
          xyzh(1,npart) = dx
          xyzh(2,npart) = dy
          xyzh(3,npart) = dz
          dy            = dy + dxyz
          nadd          = nadd + 1
       enddo
       dz = dz + dxyz
    enddo
 enddo

 ! add particles in the -y direction
 do while (ymin - dxyz > border(2,1))
    updated_bdy = .true.
    call write_bdy_update(iprint,'increasing','-y','ymin',ymin,ymin-dxyz)
    !write(iprint,*) 'Updating Boundaries: increasing -y: ymin_old, ymin_new = ',ymin,ymin-dxyz
    ymin = ymin - dxyz
    dy = ymin + 0.5*dxyz
    dz = zmin + 0.5*dxyz
    do while (dz < zmax)
       dx = xmin + 0.5*dxyz
       do while (dx < xmax)
          npart = min(npart + 1,maxp_hard)
          call copy_particle_all(ibkg,npart,.true.)
          xyzh(1,npart) = dx
          xyzh(2,npart) = dy
          xyzh(3,npart) = dz
          dx            = dx + dxyz
          nadd          = nadd + 1
       enddo
       dz = dz + dxyz
    enddo
 enddo

 ! add particles in the -z direction
 do while (zmin - dxyz > border(3,1))
    updated_bdy = .true.
    call write_bdy_update(iprint,'increasing','-z','zmin',zmin,zmin-dxyz)
    !write(iprint,*) 'Updating Boundaries: increasing -z: zmin_old, zmin_new = ',zmin,zmin-dxyz
    zmin = zmin - dxyz
    dy = ymin + 0.5*dxyz
    dz = zmin + 0.5*dxyz
    do while (dy < ymax)
       dx = xmin + 0.5*dxyz
       do while (dx < xmax)
          npart = min(npart + 1,maxp_hard)
          call copy_particle_all(ibkg,npart,.true.)
          xyzh(1,npart) = dx
          xyzh(2,npart) = dy
          xyzh(3,npart) = dz
          dx            = dx + dxyz
          nadd          = nadd + 1
       enddo
       dy = dy + dxyz
    enddo
 enddo

 ! add particles in the +x direction
 do while (xmax + dxyz < border(1,2))
    updated_bdy = .true.
    call write_bdy_update(iprint,'increasing','+x','xmax',xmax,xmax+dxyz)
    !write(iprint,*) 'Updating Boundaries: increasing +x: xmax_old, xmax_new = ',xmax,xmax+dxyz
    xmax = xmax + dxyz
    dx = xmax - 0.5*dxyz
    dz = zmin + 0.5*dxyz
    do while (dz < zmax)
       dy = ymin + 0.5*dxyz
       do while (dy < ymax)
          npart = min(npart + 1,maxp_hard)
          call copy_particle_all(ibkg,npart,.true.)
          xyzh(1,npart) = dx
          xyzh(2,npart) = dy
          xyzh(3,npart) = dz
          dy            = dy + dxyz
          nadd          = nadd + 1
       enddo
       dz = dz + dxyz
    enddo
 enddo

 ! add particles in the +y direction
 do while (ymax + dxyz < border(2,2))
    updated_bdy = .true.
    call write_bdy_update(iprint,'increasing','+y','ymax',ymax,ymax+dxyz)
    !write(iprint,*) 'Updating Boundaries: increasing +y: ymax_old, ymax_new = ',ymax,ymax+dxyz
    ymax = ymax + dxyz
    dy = ymax - 0.5*dxyz
    dz = zmin + 0.5*dxyz
    do while (dz < zmax)
       dx = xmin + 0.5*dxyz
       do while (dx < xmax)
          npart = min(npart + 1,maxp_hard)
          call copy_particle_all(ibkg,npart,.true.)
          xyzh(1,npart) = dx
          xyzh(2,npart) = dy
          xyzh(3,npart) = dz
          dx            = dx + dxyz
          nadd          = nadd + 1
       enddo
       dz = dz + dxyz
    enddo
 enddo

 ! add particles in the +z direction
 do while (zmax + dxyz < border(3,2))
    updated_bdy = .true.
    call write_bdy_update(iprint,'increasing','+z','zmax',zmax,zmax+dxyz)
    !write(iprint,*) 'Updating Boundaries: increasing +z: zmax_old, zmax_new = ',zmax,zmax+dxyz
    zmax = zmax + dxyz
    dz = zmax - 0.5*dxyz
    dy = ymin + 0.5*dxyz
    do while (dy < ymax)
       dx = xmin + 0.5*dxyz
       do while (dx < xmax)
          npart = min(npart + 1,maxp_hard)
          call copy_particle_all(ibkg,npart,.true.)
          xyzh(1,npart) = dx
          xyzh(2,npart) = dy
          xyzh(3,npart) = dz
          dx            = dx + dxyz
          nadd          = nadd + 1
       enddo
       dy = dy + dxyz
    enddo
 enddo
 npartoftype(igas) = npartoftype(igas) + nadd

 !--Failsafe
 if (npart==maxp_hard) call fatal('update_boundary','npart >=maxp_hard.  Recompile with larger maxp and rerun')

 !--Reset boundaries to remove particles
 do while (xmin + dxyz < border(1,1))
    call write_bdy_update(iprint,'decreasing','-x','xmin',xmin,xmin+dxyz)
    !write(iprint,*) 'Updating Boundaries: decreasing -x: xmin_old, xmin_new = ',xmin,xmin+dxyz
    xmin = xmin + dxyz
 enddo
 do while (ymin + dxyz < border(2,1))
    call write_bdy_update(iprint,'decreasing','-y','ymin',ymin,ymin+dxyz)
    !write(iprint,*) 'Updating Boundaries: decreasing -y: ymin_old, ymin_new = ',ymin,ymin+dxyz
    ymin = ymin + dxyz
 enddo
 do while (zmin + dxyz < border(3,1))
    call write_bdy_update(iprint,'decreasing','-z','zmin',zmin,zmin+dxyz)
    !write(iprint,*) 'Updating Boundaries: decreasing -z: zmin_old, zmin_new = ',zmin,zmin+dxyz
    zmin = zmin + dxyz
 enddo
 do while (xmax - dxyz > border(1,2))
    call write_bdy_update(iprint,'decreasing','+x','xmax',xmax,xmax-dxyz)
    !write(iprint,*) 'Updating Boundaries: decreasing +x: xmax_old, xmax_new = ',xmax,xmax-dxyz
    xmax = xmax - dxyz
 enddo
 do while (ymax - dxyz > border(2,2))
    call write_bdy_update(iprint,'decreasing','+y','ymax',ymax,ymax-dxyz)
    !write(iprint,*) 'Updating Boundaries: decreasing +y: ymax_old, ymax_new = ',ymax,ymax-dxyz
    ymax = ymax - dxyz
 enddo
 do while (zmax - dxyz > border(3,2))
    call write_bdy_update(iprint,'decreasing','+z','zmax',zmax,zmax-dxyz)
    !write(iprint,*) 'Updating Boundaries: decreasing +z: zmax_old, zmax_new = ',zmax,zmax-dxyz
    zmax = zmax - dxyz
 enddo

 !--For each new particle, replace its velocity and magnetic field with the average values
 if (n_bkg > 0) then
    !$omp parallel default(none) &
    !$omp shared(npart0,npart,v_bkg,B_bkg,vxyzu,Bxyz) &
    !$omp private(i)
    !$omp do
    do i = npart0+1,npart
       vxyzu(:,i) = v_bkg
       if (mhd) Bxyz(:,i) = B_bkg
    enddo
    !$omp enddo
    !$omp end parallel
 endif

 !--Kill unwanted particles
 border(1,1) = xmin
 border(2,1) = ymin
 border(3,1) = zmin
 border(1,2) = xmax
 border(2,2) = ymax
 border(3,2) = zmax
 do i = 1,npart
    if (.not. in_domain(xyzh(1:3,i),border)) then
       ndie = ndie + 1
       call kill_particle(i,npartoftype) ! This subroutine is not thread-safe!
    elseif (remove_accreted .and. isdead_or_accreted(xyzh(4,i))) then
       naccreted = naccreted + 1
       call kill_particle(i,npartoftype) ! This subroutine is not thread-safe!
    endif
 enddo

 !--Shuffle particles to remove the dead particles from the list
 if (ndie > 0 .or. naccreted > 0) call shuffle_part(npart)

 !--Reset logical to calculate properties on boundaries once & all particle numbers
 nactive = npart
 nalive  = nactive
 dxbound = xmax - xmin
 dybound = ymax - ymin
 dzbound = zmax - zmin
 totvol  = dxbound*dybound*dzbound

 !--Cleanly end at next full dump if we predict to go over maxp_hard next time we add particles
 if (nadd+npart > maxp_hard) abortrun_bdy = .true.

 !--Final print-statements
 if (nadd > 0 .or. ndie > 0) then
    write(iprint,'(a,3I12)') 'Updated  boundaries: initial npart, final npart, particles added, particles killed: ', &
    npart0,npart,nadd,ndie
 endif
 if (naccreted > 0) then
    write(iprint,'(a,I12,a)') 'Updated  boundaries: removed ',naccreted,' accreted particles from the list'
 endif

end subroutine update_boundaries

!-----------------------------------------------------------------------
!+
!  Short subroutine to cleanly write information regarding boundary updates
!+
!-----------------------------------------------------------------------
subroutine write_bdy_update(iprint,inde_crease,cval,cminmax,val_old,val_new)
 integer,           intent(in) :: iprint
 real,              intent(in) :: val_old,val_new
 character(len= *), intent(in) :: inde_crease,cval,cminmax
 character(len=12)             :: fmt

 if (abs(val_old) < 1.d5) then
    fmt = '(9a,2f12.5)'
 else
    fmt = '(9a,2Es14.6)'
 endif
 write(iprint,fmt) 'Updating Boundaries: ',trim(inde_crease),':',cval,': ',cminmax,'_old, ',cminmax,'_new = ',val_old,val_new

end subroutine write_bdy_update

!-----------------------------------------------------------------------
!+
!  writes boundary options to the input file (for dynamic boundaries only)
!+
!-----------------------------------------------------------------------
subroutine write_options_boundary(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 if (dynamic_bdy) then
    write(iunit,"(/,a)") '# options controlling dynamic boundaries particles [all values in code units]'
    call write_inopt(dynamic_bdy,'dynamic_bdy','turn on/off dynamic boundaries',iunit)
    call write_inopt(rho_thresh_bdy,'rho_thresh_bdy','threshold density separating dense gas from background gas',iunit)
    call write_inopt(width_bkg(1,1),'width_bkg_nx','width of the boundary in the -x direction',iunit)
    call write_inopt(width_bkg(2,1),'width_bkg_ny','width of the boundary in the -y direction',iunit)
    call write_inopt(width_bkg(3,1),'width_bkg_nz','width of the boundary in the -z direction',iunit)
    call write_inopt(width_bkg(1,2),'width_bkg_px','width of the boundary in the +x direction',iunit)
    call write_inopt(width_bkg(2,2),'width_bkg_py','width of the boundary in the +y direction',iunit)
    call write_inopt(width_bkg(3,2),'width_bkg_pz','width of the boundary in the +z direction',iunit)
    call write_inopt(vbdyx,'vbdyx','velocity of the x-boundary', iunit)
    call write_inopt(vbdyy,'vbdyy','velocity of the y-boundary', iunit)
    call write_inopt(vbdyz,'vbdyz','velocity of the z-boundary', iunit)
    call write_inopt(n_dtmax,'n_dtmax','particles must not reach v*n_dtmax*dtmax of the boundary', iunit)
 endif

end subroutine write_options_boundary

!-----------------------------------------------------------------------
!+
!  reads boundary options from the input file (for dynamic boundaries only)
!+
!-----------------------------------------------------------------------
subroutine read_options_boundary(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter  :: label = 'read_options_boundary'

 imatch  = .true.
 select case(trim(name))
 case('dynamic_bdy')
    read(valstring,*,iostat=ierr) dynamic_bdy
    ngot = ngot + 1
 case('rho_thresh_bdy')
    read(valstring,*,iostat=ierr) rho_thresh_bdy
    if (rho_thresh_bdy < 0.) call fatal(label,'rho_thresh_bdy < 0')
    ngot = ngot + 1
 case('width_bkg_nx')
    read(valstring,*,iostat=ierr) width_bkg(1,1)
    ngot = ngot + 1
 case('width_bkg_ny')
    read(valstring,*,iostat=ierr) width_bkg(2,1)
    ngot = ngot + 1
 case('width_bkg_nz')
    read(valstring,*,iostat=ierr) width_bkg(3,1)
    ngot = ngot + 1
 case('width_bkg_px')
    read(valstring,*,iostat=ierr) width_bkg(1,2)
    ngot = ngot + 1
 case('width_bkg_py')
    read(valstring,*,iostat=ierr) width_bkg(2,2)
    ngot = ngot + 1
 case('width_bkg_pz')
    read(valstring,*,iostat=ierr) width_bkg(3,2)
    ngot = ngot + 1
 case('vbdyx')
    read(valstring,*,iostat=ierr) vbdyx
    ngot = ngot + 1
 case('vbdyy')
    read(valstring,*,iostat=ierr) vbdyy
    ngot = ngot + 1
 case('vbdyz')
    read(valstring,*,iostat=ierr) vbdyz
    ngot = ngot + 1
 case('n_dtmax')
    read(valstring,*,iostat=ierr) n_dtmax
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
 if (dynamic_bdy) then
    igotall = (ngot == 12)
 else
    igotall = .true.
 endif

end subroutine read_options_boundary

!-----------------------------------------------------------------------
end module boundary_dyn
