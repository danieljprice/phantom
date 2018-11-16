module inject
 use physcon, only:pi
 implicit none

 character(len=*), parameter, public :: inject_type = 'grbondi'

 !-- Public subroutines in module
 public :: init_inject,          &
           inject_particles,     &
           write_options_inject, &
           read_options_inject

 !-- Runtime variables read from input file
 real,    public :: rin           = 18.1
 real,    public :: drdp          = 1.
 integer, public :: iboundspheres = 3

 !-- Variables calculated from the previous parameters
 real,    public  :: dtsphere
 real,    private :: masssphere,neighdist,vin
 integer, private :: npsphere,iwindres

 private

 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 logical, parameter :: wind_verbose = .false.

contains

! Wrapper
subroutine get_solution(rho,v,u,r,gamma)
 use bondiexact, only:get_bondi_solution
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,gamma
 real, parameter   :: mass1 = 1.

 !-- n.b. GR code is written such that central mass is always 1 (i.e. distance, time, and BH spin are in units of M always)
 call get_bondi_solution(rho,v,u,r,mass1,gamma)

end subroutine get_solution

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use injectutils, only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 use icosahedron, only:compute_matrices,compute_corners
 use part,        only:massoftype,igas,iboundary
 use eos,         only:gamma
 use io,          only:iverbose
 integer, intent(out) :: ierr
 real :: mVonMdotR,masspart
 real :: dr,dp,masspart1,mdot,speedin
 real :: rhoin,uthermin

!-- return without error
 ierr = 0

 call get_solution(rhoin,vin,uthermin,rin,gamma)
 speedin = abs(vin)
 mdot  = 4.*pi*rin**2*rhoin*speedin

!-- compute the dimensionless resolution factor m V / (Mdot R)
!   where m = particle mass and V, Mdot and R are wind parameters
 masspart              = massoftype(igas)
 massoftype(iboundary) = masspart
 mVonMdotR             = masspart*speedin/(mdot*rin)

!-- solve for the integer resolution of the geodesic spheres
!   gives number of particles on the sphere via N = 20*(2*q*(q - 1)) + 12
 iwindres  = get_sphere_resolution(drdp,mVonMdotR)
 npsphere  = get_parts_per_sphere(iwindres)
 neighdist = get_neighb_distance(iwindres)

 masssphere = masspart*npsphere
 dtsphere   = masssphere/mdot

 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 if (iverbose >= 1) then
    print*,'mass of particle = ',massoftype(igas)
    masspart1 = drdp*get_neighb_distance(4)*rin*mdot/(get_parts_per_sphere(4)*speedin)
    print*,'require mass of particle = ',masspart1,' to get 492 particles per sphere'
    print*,'Mdot*R/(m*V) ',1./mVonMdotR
    print*,'wind_resolution ',iwindres
    print*,'npsphere ',npsphere
    print*,'neighdist ',neighdist
    dp = neighdist*rin
    dr = speedin*masspart*npsphere/mdot
    print*,'particle separation on spheres = ',dp
    print*,'distance between spheres = ',dr
    print*,'got dr/dp = ',dr/dp,' compared to desired dr on dp = ',drdp
    print*,'masspart ',masspart
    print*,'masssphere ',masssphere
    print*,'dtsphere ',dtsphere
 endif

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,dtinject)
 use io,          only:iprint,warning
 use eos,         only:gamma
 use part,        only:igas,iboundary
 use injectutils, only:inject_geodesic_sphere
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real,    parameter          :: pi3   = pi/3. !-- irrational number close to one
 real,    parameter          :: shift = 0.
 integer                     :: outer_sphere, inner_sphere, inner_boundary_sphere, i, ierr, itype, ipartbegin
 real                        :: tlocal, GM, r, v, u, rho, e, x0(3), v0(3)
 logical, save               :: first_run = .true.
 character(len=*), parameter :: label = 'inject_particles'

 if (first_run) then
    call init_inject(ierr)
    first_run = .false.
 endif

 x0 = 0.
 v0 = 0.
 GM = 0.

 outer_sphere = floor((time-dtlast)/dtsphere) + 1
 inner_sphere = floor(time/dtsphere)
 inner_boundary_sphere = inner_sphere + iboundspheres

 if (inner_sphere-outer_sphere > iboundspheres) call warning(label,'problem with boundary spheres, timestep likely too large!')
! cs2max = 0.

 do i=inner_sphere+iboundspheres,outer_sphere,-1
    tlocal = time - (i-shift) * dtsphere
    call compute_sphere_properties(time,tlocal,gamma,GM,r,v,u,rho,e,i,inner_sphere,inner_boundary_sphere)

    if (wind_verbose) then
       write(iprint,*) '   Sphere            : ',i,(i-shift)
       write(iprint,*) '   Local Time        : ',tlocal,time,(i-shift)*dtsphere
       write(iprint,*) '   Radius            : ',r
       write(iprint,*) '   Expansion velocity: ',v
       write(iprint,*) '   Density           : ',rho
       write(iprint,*) ''
    endif

    if (i > inner_sphere) then
       !-- boundary sphere
       ipartbegin = (iboundspheres-i+inner_sphere)*npsphere+1
    else
       !-- injected sphere
       ipartbegin = npart + 1
    endif

    itype = igas

    call inject_geodesic_sphere(i,ipartbegin,iwindres,r,v,u,rho,geodesic_R,geodesic_V,npart,npartoftype,xyzh,vxyzu,itype,x0,v0)
 enddo

!-- Return timestep constraint to ensure that time between sphere
!   injections is adequately resolved
 dtinject = .5*pi3*dtsphere

end subroutine inject_particles


!-----------------------------------------------------------------------
!+
!  Compute the radius and velocity of a sphere at the current local time
!+
!-----------------------------------------------------------------------
subroutine compute_sphere_properties(time,tlocal,gamma,GM,r,v,u,rho,e,isphere,inner_sphere,inner_boundary_sphere)
 integer, intent(in)  :: isphere, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,tlocal,gamma,GM
 real,    intent(out) :: r, v, u, rho, e

 call integrate_solution(tlocal,r,v,u,rho,gamma)
 e = 0. !???

end subroutine compute_sphere_properties

subroutine integrate_solution(tlocal,r,v,u,rho,gamma)
 real, intent(in)   :: tlocal,gamma
 real, intent(out)  :: r,v,u,rho
 integer, parameter :: N = 10000
 integer :: i
 real    :: dt,v1,v2,v3,v4

 dt = tlocal/N
 r  = rin
 v  = vin

!--- n.b. I don't know if RK4 is necessary here, perhaps simple Euler would do, but speed isn't really affected it seems.
 ! iterations
 do i=1,N
    call get_solution(rho,v1,u,r,gamma)
    call get_solution(rho,v2,u,r+dt/2.*v1,gamma)
    call get_solution(rho,v3,u,r+dt/2.*v2,gamma)
    call get_solution(rho,v4,u,r+dt*v3,gamma)
    r = r + dt/6. * (v1 + 2.*v2 + 2.*v3 + v4)
    call get_solution(rho,v,u,r,gamma)
 enddo

end subroutine integrate_solution

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(rin          ,'rin'          ,'radius of injection of the wind'                           ,iunit)
 call write_inopt(drdp         ,'drdp'         ,'desired ratio of sphere spacing to particle spacing'       ,iunit)
 call write_inopt(iboundspheres,'iboundspheres','number of boundary spheres (integer)'                      ,iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save                 :: ngot = 0
 integer                       :: noptions
 character(len=30), parameter  :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('rin')
    read(valstring,*,iostat=ierr) rin
    ngot = ngot + 1
    if (rin < 0.) call fatal(label,'invalid setting for rin (<0)')
 case('drdp')
    read(valstring,*,iostat=ierr) drdp
    ngot = ngot + 1
    if (drdp <= 0.)        call fatal(label,'drdp must be >=0')
 case('iboundspheres')
    read(valstring,*,iostat=ierr) iboundspheres
    ngot = ngot + 1
    if (iboundspheres <= 0)     call fatal(label,'iboundspheres must be > 0')
 case default
    imatch = .false.
 end select

 noptions = 3
 igotall  = (ngot >= noptions)

end subroutine read_options_inject

end module inject
