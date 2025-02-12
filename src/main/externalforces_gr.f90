!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module externalforces
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - accradius1      : *soft accretion radius of black hole*
!   - accradius1_hard : *hard accretion radius of black hole*
!
! :Dependencies: dump_utils, infile_utils, io, metric, metric_tools, part,
!   units
!
 use metric, only:mass1
 implicit none

 private
 public :: externalforce,externalforce_vdependent
 public :: accrete_particles,was_accreted
 public :: write_options_externalforces,read_options_externalforces
 public :: initialise_externalforces,is_velocity_dependent
 public :: update_vdependent_extforce
 public :: update_externalforce
 public :: write_headeropts_extern,read_headeropts_extern

 !
 ! enumerated list of external forces
 !
 integer, parameter, public :: iext_gr = 1

 public :: mass1  ! exported from metric module
 real, public :: accradius1 = 0.
 real, public :: accradius1_hard = 0.
 real, public :: accretedmass1 = 0.
 real, public :: accretedmass2 = 0.

 logical, public :: extract_iextern_from_hdr = .false.

 ! (the following for compatibility with non-relativistic code)
 integer, parameter, public :: iext_lensethirring = -1
 integer, parameter, public :: iext_einsteinprec = -2
 integer, parameter, public :: iext_binary = -3
 integer, parameter, public :: iext_spiral = -4
 integer, parameter, public :: iext_star = -5
 integer, parameter, public :: iext_corotate = -6
 integer, parameter, public :: iext_corot_binary = -7
 integer, parameter, public :: iext_gwinspiral = -8
 integer, parameter, public :: iext_densprofile = -9
 real, public :: omega_corotate = 0.

 !
 ! Human-readable labels for these
 !
 integer, parameter, public  :: iexternalforce_max = 1
 character(len=*), parameter, public :: externalforcetype(iexternalforce_max) = (/'gr'/)

contains
!-----------------------------------------------------------------------
!+
!  Computes external (body) forces on a particle given its co-ordinates
!  This doesn't doesn't actually get used in gr...
!+
!-----------------------------------------------------------------------
subroutine externalforce(iexternalforce,xi,yi,zi,hi,ti,fextxi,fextyi,fextzi,phi,dtf,ii)
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: xi,yi,zi,hi,ti
 real,    intent(out) :: fextxi,fextyi,fextzi,phi
 real,    intent(out), optional :: dtf
 integer, intent(in),  optional :: ii ! NOTE: index-base physics can be dangerous;
 !
 !  This doesn't doesn't actually get used in gr...
 !
end subroutine externalforce

!-----------------------------------------------------------------------
!+
!  Query function to determine whether or not an external force
!  has velocity dependent parts
!+
!-----------------------------------------------------------------------
logical function is_velocity_dependent(iexternalforce)
 integer, intent(in) :: iexternalforce

 select case(iexternalforce)
 case default
    is_velocity_dependent = .true.
 end select

end function is_velocity_dependent

!-----------------------------------------------------------------------
!+
!  Returns the part of the external force that is velocity dependent
!  (these are treated implicitly in the leapfrog integrator)
!  This routine returns an explicit evaluation
!
! This doesn't doesn't actually get used in gr...
!
!+
!-----------------------------------------------------------------------
subroutine externalforce_vdependent(iexternalforce,xyzi,veli,fexti,poti,densi,ui)
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: xyzi(3),veli(3)
 real,    intent(out) :: fexti(3)
 real,    intent(inout) :: poti
 real,    intent(in),    optional :: densi
 real,    intent(inout), optional :: ui
 !
 ! This doesn't doesn't actually get used in gr...
 !
end subroutine externalforce_vdependent

!-----------------------------------------------------------------------
!+
!  Solves for velocity-dependent part of external force via an
!  implicit inversion of v^1 = v^1/2 + 1/2*dt*f0 + 1/2*dt*f1(x^1,v^1)
!  necessary for using v-dependent forces in leapfrog
!+
!-----------------------------------------------------------------------
subroutine update_vdependent_extforce(iexternalforce, &
           vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,densi,ui)
 integer, intent(in)    :: iexternalforce
 real,    intent(in)    :: dt,xi,yi,zi
 real,    intent(in)    :: vhalfx,vhalfy,vhalfz
 real,    intent(inout) :: fxi,fyi,fzi
 real,    intent(out)   :: fexti(3)
 real,    intent(in),    optional :: densi
 real,    intent(inout), optional :: ui
 !
 ! This doesn't doesn't actually get used in gr...
 !
end subroutine update_vdependent_extforce

!-----------------------------------------------------------------------
!+
!  update time-dependent external forces where necessary
!+
!-----------------------------------------------------------------------
subroutine update_externalforce(iexternalforce,ti,dmdt)
 integer, intent(in) :: iexternalforce
 real,    intent(in) :: ti,dmdt

end subroutine update_externalforce

!-----------------------------------------------------------------------
!+
!  Performs checks to see whether or not a particle should be
!  accreted by crossing a boundary of the external potential
!  (at the moment this is just a hard boundary, but could in principle
!   add checks to see if particle is bound etc. here)
!+
!-----------------------------------------------------------------------
subroutine accrete_particles(iexternalforce,xi,yi,zi,hi,mi,ti,accreted,i)
 use metric_tools, only:imet_minkowski,imet_schwarzschild,imet_kerr,imetric
 use part,         only:set_particle_type,iboundary,maxphase,maxp,igas,npartoftype
 integer, intent(in)    :: iexternalforce
 real,    intent(in)    :: xi,yi,zi,mi,ti
 real,    intent(inout) :: hi
 logical, intent(out)   :: accreted
 integer, intent(in), optional :: i
 logical, save :: first = .true.
 real :: r2

 accreted = .false.
 select case(imetric)
 case(imet_minkowski)
    if (first) print*,"WARNING: Accrete particles: but Metric = Minkowski"

 case(imet_schwarzschild,imet_kerr)
    r2 = xi*xi + yi*yi + zi*zi
    if (accradius1>accradius1_hard .and. r2 < accradius1**2 .and. maxphase==maxp .and. present(i)) then
       call set_particle_type(i,iboundary)
       npartoftype(igas) = npartoftype(igas) - 1
       npartoftype(iboundary) = npartoftype(iboundary) + 1
    endif
    if (r2 < (accradius1_hard)**2) accreted = .true.

 end select

 if (accreted) then
    hi = -abs(hi)
 endif

 first = .false.

end subroutine accrete_particles

!-----------------------------------------------------------------------
!+
!  query function for particles that have already been accreted
!+
!-----------------------------------------------------------------------
pure logical function was_accreted(iexternalforce,hi)
 use metric_tools, only:imet_minkowski,imetric
 integer, intent(in) :: iexternalforce
 real,    intent(in) :: hi

 ! An accreted particle is indicated by h < 0.
 ! Note less than, but not equal.
 ! (h=0 indicates dead MPI particle)
 was_accreted = (hi < 0.)

 if (imetric==imet_minkowski) was_accreted = .false.

end function was_accreted

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_externalforces(iunit,iexternalforce)
 use metric_tools, only:imet_minkowski,imetric
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit,iexternalforce

 if (imetric /= imet_minkowski) then
    write(iunit,"(/,a)") '# options relating to GR external forces'
    if (accradius1_hard < tiny(0.)) accradius1_hard = accradius1
    call write_inopt(accradius1,'accradius1','soft accretion radius of black hole',iunit)
    call write_inopt(accradius1_hard,'accradius1_hard','hard accretion radius of black hole',iunit)
 endif

end subroutine write_options_externalforces

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_extern(iexternalforce,hdr,time,ierr)
 use dump_utils, only:dump_h
 integer,      intent(in)    :: iexternalforce
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time
 integer,      intent(out)   :: ierr

end subroutine write_headeropts_extern

!-----------------------------------------------------------------------
!+
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_extern(iexternalforce,hdr,nptmass,ierr)
 use dump_utils, only:dump_h
 integer,      intent(in)  :: iexternalforce,nptmass
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr = 0

end subroutine read_headeropts_extern

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_externalforces(name,valstring,imatch,igotall,ierr,iexternalforce)
 use io,           only:fatal,warn
 use metric_tools, only:imet_minkowski,imetric
 character(len=*), intent(in)    :: name,valstring
 logical,          intent(out)   :: imatch,igotall
 integer,          intent(out)   :: ierr
 integer,          intent(inout) :: iexternalforce
 integer, save :: ngot = 0
 character(len=*), parameter :: tag = 'externalforces_gr'

 imatch            = .true.
 igotall           = .false.

 if (imetric /= imet_minkowski) then

    select case(trim(name))
    case('accradius1')
       read(valstring,*,iostat=ierr) accradius1
       if (accradius1 < 0.)    call fatal(tag,'negative accretion radius')
       if (imetric == imet_minkowski) call warn(tag,'Minkowski metric: ignoring accradius1 value')
       ngot = ngot + 1
    case('accradius1_hard')
       read(valstring,*,iostat=ierr) accradius1_hard
       if (accradius1_hard > accradius1) call fatal(tag,'hard accretion boundary must be within soft accretion boundary')
       if (imetric == imet_minkowski) call warn(tag,'Minkowski metric: ignoring accradius1_hard value')
       ngot = ngot + 1
    case default
       imatch = .false.
    end select

    igotall = (ngot >= 2)

 else

    igotall = .true.
    imatch  = .false.

    if (accradius1 > 0.)      call fatal(tag,'accradius1 > 0 when metric = Minkowski')
    if (accradius1_hard > 0.) call fatal(tag,'accradius1_hard > 0 when metric = Minkowski')

 endif

end subroutine read_options_externalforces

!-----------------------------------------------------------------------
!+
!  interface to initialisation/checking routines for external forces
!+
!-----------------------------------------------------------------------
subroutine initialise_externalforces(iexternalforce,ierr)
 use io,    only:error
 use units, only:G_is_unity,c_is_unity,get_G_code,get_c_code
 integer, intent(in)  :: iexternalforce
 integer, intent(out) :: ierr

 ierr = 0

 !
 !--check that G=1 in code units
 !
 if (.not.G_is_unity()) then
    call error('units',trim(externalforcetype(iexternalforce))//&
               ' external force assumes G=1 in code units but we have',var='G',val=real(get_G_code()))
    ierr = ierr + 1
 endif

 !
 !--check that c=1 in code units
 !
 if (.not.c_is_unity()) then
    call error('units',trim(externalforcetype(iexternalforce))//&
                  ' external force assumes c=1 in code units but we have',var='c',val=real(get_c_code()))
    ierr = ierr + 1
 endif

end subroutine initialise_externalforces

end module externalforces
