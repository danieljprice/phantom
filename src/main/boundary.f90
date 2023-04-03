!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module boundary
!
! This module contains variables and subroutines relating to boundaries,
! including dynamically adjusting periodic boundaries
!
! :References: Wurster & Bonnell (2023) for dynamic boundaries
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
! :Dependencies: dim, infile_utils, io, kernel, mpidomain, part
!

 use dim, only: maxvxyzu
 implicit none
 real,    public :: xmin,xmax,ymin,ymax,zmin,zmax
 real,    public :: dxbound,dybound,dzbound
 real,    public :: totvol

 ! Dynamic Boundaries: Initialise values for the .in & dump files
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
 real,    public :: rho_bkg_ini1,ndtmaxdtmax
 integer, public :: ibkg,n_bkg
 logical, public :: vbdyx_not_0 = .false.
 logical, public :: vbdyy_not_0 = .false.
 logical, public :: vbdyz_not_0 = .false.
 logical         :: remove_accreted = .true.      ! remove accreted particles from the list
 real,    public, parameter :: bdy_vthresh = 0.05 ! tolerance on velocity
 ! if a particle's velocity is within this tolerance, then it is a
 ! 'boring' particle that may be used to determine the background.
 ! Splash dumps suggest 1% is too low and 10% is too high.

 public :: set_boundary
 public :: cross_boundary
 public :: in_domain,update_boundaries,init_dynamic_bdy
 public :: write_options_boundary,read_options_boundary

 private

contains

!---------------------------------------------------------------
!+
!  Routine to set the (cartesian) boundaries for the code
!  when using periodic boundary conditions
!
!  This sets the module variables xmin,xmax,ymin,ymax,zmin,zmax
!  as well as the subsidiary variables dxbound,dybound,dzbound
!  and totvol.
!
!  Can be called with no arguments, which gives the defaults:
!
!   call set_boundary()
!
!  Or with each boundary set individually:
!
!   call set_boundary(0.,1.,0.,1.,0.,1.)
!
!  Or with an array of values array=(/xmin,xmax,ymin,ymax,zmin,zmax/)
!
!   call set_boundary(pos=array)
!+
!---------------------------------------------------------------
subroutine set_boundary(x_min,x_max,y_min,y_max,z_min,z_max,pos)
 real, intent(in), optional :: x_min, x_max, y_min, y_max, z_min, z_max
 real, intent(in), optional :: pos(6)
 !
 ! give default values if no settings given
 !
 xmin = -0.5
 xmax = 0.5
 ymin = -0.5
 ymax = 0.5
 zmin = -0.5
 zmax = 0.5
 !
 ! set the min and max position in each coordinate
 !
 if (present(pos)) then
    xmin = pos(1)
    xmax = pos(2)
    ymin = pos(3)
    ymax = pos(4)
    zmin = pos(5)
    zmax = pos(6)
 else
    if (present(x_min)) xmin = x_min
    if (present(x_max)) xmax = x_max
    if (present(y_min)) ymin = y_min
    if (present(y_max)) ymax = y_max
    if (present(z_min)) zmin = z_min
    if (present(z_max)) zmax = z_max
 endif
 !
 ! set subsidiary quantities that depend on above settings
 ! these are to save computation in periodicity calculations
 !
 dxbound = xmax - xmin
 dybound = ymax - ymin
 dzbound = zmax - zmin

 totvol = dxbound*dybound*dzbound

end subroutine set_boundary

!----------------------------------------------------------------
!+
!  This subroutine determines whether particles should cross
!  the boundary
!+
!---------------------------------------------------------------
subroutine cross_boundary(isperiodic,xyz,ncross)
 logical, intent(in)    :: isperiodic(3)
 real,    intent(inout) :: xyz(:)
 integer, intent(inout) :: ncross
!
!--check for crossings of periodic boundaries
!
 if (isperiodic(1)) then
    if (xyz(1) < xmin) then
       xyz(1) = xyz(1) + dxbound
       ncross = ncross + 1
    elseif (xyz(1) > xmax) then
       xyz(1) = xyz(1) - dxbound
       ncross = ncross + 1
    endif
 endif

 if (isperiodic(2)) then
    if (xyz(2) < ymin) then
       xyz(2) = xyz(2) + dybound
       ncross = ncross + 1
    elseif (xyz(2) > ymax) then
       xyz(2) = xyz(2) - dybound
       ncross = ncross + 1
    endif
 endif

 if (isperiodic(3)) then
    if (xyz(3) < zmin) then
       xyz(3) = xyz(3) + dzbound
       ncross = ncross + 1
    elseif (xyz(3) > zmax) then
       xyz(3) = xyz(3) - dzbound
       ncross = ncross + 1
    endif
 endif

 return
end subroutine cross_boundary

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
!  Determine if the boundary is moving and in what direction
!+
!---------------------------------------------------------------
subroutine init_dynamic_bdy()
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

end subroutine init_dynamic_bdy
!----------------------------------------------------------------
!+
!  Update the dynamic boundaries
!+
!---------------------------------------------------------------
subroutine update_boundaries(nactive,nalive,npart,abortrun_bdy)
 use dim,       only: maxp_hard,mhd
 use mpidomain, only: isperiodic
 use io,        only: fatal,iprint
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
 use io,         only:fatal
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
end module boundary
