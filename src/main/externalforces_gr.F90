module externalforces
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 private
 public :: externalforce,externalforce_vdependent
 public :: accrete_particles,was_accreted
 public :: write_options_externalforces,read_options_externalforces
 public :: initialise_externalforces,is_velocity_dependent
 public :: update_vdependent_extforce_leapfrog
 public :: update_externalforce
 public :: write_headeropts_extern,read_headeropts_extern

 !
 ! enumerated list of external forces
 !
 integer, parameter, public :: iext_gr = 1

 real, public :: accradius1 = 0.
 real, public :: accradius1_hard = 0.


 ! (the following for compatibility with non-relativistic code)
 integer, parameter, public :: iext_lensethirring = -1
 integer, parameter, public :: iext_einsteinprec = -1
 integer, parameter, public :: iext_binary = -1
 integer, parameter, public :: iext_spiral = -1
 integer, parameter, public :: iext_star = -1
 integer, parameter, public :: iext_corotate = -1

 !
 ! Human-readable labels for these
 !
 integer, parameter, public  :: iexternalforce_max = 1
 character(len=*), parameter, public :: externalforcetype(iexternalforce_max) = (/'gr'/)

contains
!-----------------------------------------------------------------------
!+
!  Computes external (body) forces on a particle given its co-ordinates
!+
!-----------------------------------------------------------------------
subroutine externalforce(iexternalforce,xi,yi,zi,hi,ti,fextxi,fextyi,fextzi,phi,dtf,ii)
#ifdef FINVSQRT
 use fastmath,         only:finvsqrt
#endif
 use io,                 only:fatal
 !use part,               only:rhoh,massoftype,igas
 ! use force_gr,           only:get_forcegr
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: xi,yi,zi,hi,ti
 real,    intent(out) :: fextxi,fextyi,fextzi,phi
 real,    intent(out), optional :: dtf
 integer, intent(in),  optional :: ii ! NOTE: index-base physics can be dangerous;
 real :: f2i, dtf1, dtf2

!-----------------------------------------------------------------------
!
!--set external force to zero
!
 fextxi = 0.
 fextyi = 0.
 fextzi = 0.
 phi    = 0.

 ! call get_forcegr(x,v,dens,u,p,fterm)

!
!--return a timestep based only on the external force
!  so that we can do substeps with only the external force call
!
 if (present(dtf)) then
    f2i = fextxi*fextxi + fextyi*fextyi + fextzi*fextzi
    if (abs(f2i) > epsilon(f2i)) then
       !
       !--external force timestep based on sqrt(h/accel)
       !
       if (hi > epsilon(hi)) then
#ifdef FINVSQSRT
          dtf1 = sqrt(hi*finvsqrt(f2i))
#else
          dtf1 = sqrt(hi/sqrt(f2i))
#endif
       else
          dtf1 = huge(dtf1)
       endif
       !
       !--external force timestep based on sqrt(phi)/accel
       !
       if (abs(phi) > epsilon(phi)) then
          dtf2 = sqrt(abs(phi)/f2i)
       else
          dtf2 = huge(dtf2)
       endif
       dtf  = min(dtf1,dtf2)
       !if (dtf2 < dtf1) print*,' phi timestep = ',dtf2,' h/a = ',dtf1, ' ratio = ',dtf2/dtf1
    else
       dtf = huge(dtf)
    endif
 endif

 return
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
!+
!-----------------------------------------------------------------------
subroutine externalforce_vdependent(iexternalforce,xyzi,veli,fexti,poti,densi,ui)
 use extern_gr, only:get_grforce
 use eos,       only:equationofstate,ieos
 use io,        only:fatal
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: xyzi(3),veli(3)
 real,    intent(out) :: fexti(3)
 real,    intent(inout) :: poti
 real,    intent(in),    optional :: densi
 real,    intent(inout), optional :: ui
 real :: pi,pondensi,spsoundi,dtf

 if (.not. present(densi) .or. .not. present(ui)) call fatal('externalforce_vdependent','densi and ui not present')
 call equationofstate(ieos,pondensi,spsoundi,densi,xyzi(1),xyzi(2),xyzi(3),ui)
 pi = pondensi*densi

 call get_grforce(xyzi,veli,densi,ui,pi,fexti,dtf)

end subroutine externalforce_vdependent

!-----------------------------------------------------------------------
!+
!  Solves for velocity-dependent part of external force via an
!  implicit inversion of v^1 = v^1/2 + 1/2*dt*f0 + 1/2*dt*f1(x^1,v^1)
!  necessary for using v-dependent forces in leapfrog
!+
!-----------------------------------------------------------------------
subroutine update_vdependent_extforce_leapfrog(iexternalforce, &
           vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,densi,ui)
 use extern_gr, only:update_grforce_leapfrog
 use eos,       only:equationofstate,ieos
 use io,        only:fatal
 integer, intent(in)    :: iexternalforce
 real,    intent(in)    :: dt,xi,yi,zi
 real,    intent(in)    :: vhalfx,vhalfy,vhalfz
 real,    intent(inout) :: fxi,fyi,fzi
 real,    intent(out)   :: fexti(3)
 real,    intent(in),    optional :: densi
 real,    intent(inout), optional :: ui
 real :: pi,pondensi,spsoundi

 if (.not. present(densi) .or. .not. present(ui)) call fatal('update_vdependent_extforce_leapfrog','densi and ui not present')
 call equationofstate(ieos,pondensi,spsoundi,densi,xi,yi,zi,ui)
 pi = pondensi*densi

 call update_grforce_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,densi,ui,pi)

end subroutine update_vdependent_extforce_leapfrog

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
 use metric,       only:imetric
 use metric_tools, only:imet_minkowski,imet_schwarzschild,imet_kerr
 use part,         only:set_particle_type,iboundary,maxphase,maxp,igas,npartoftype
 integer, intent(in)    :: iexternalforce
 real,    intent(in)    :: xi,yi,zi,mi,ti
 real,    intent(inout) :: hi
 logical, intent(out)   :: accreted
 integer, intent(in), optional :: i
 integer, save :: ifirst = 0
 real :: r2

 accreted = .false.
 select case(imetric)
 case(imet_minkowski)
    if(ifirst==0) print*,"WARNING: Accrete particles: but Metric = Minkowski"

 case(imet_schwarzschild,imet_kerr)
    r2 = xi*xi + yi*yi + zi*zi
    if (r2 < accradius1**2 .and. maxphase==maxp .and. present(i)) then
       call set_particle_type(i,iboundary)
       npartoftype(igas) = npartoftype(igas) - 1
       npartoftype(iboundary) = npartoftype(iboundary) + 1
    endif
    if (r2 < (accradius1_hard)**2) accreted = .true.

 end select

 if (accreted) then
    hi = -abs(hi)
 endif

 ifirst = ifirst + 1

end subroutine accrete_particles

!-----------------------------------------------------------------------
!+
!  query function for particles that have already been accreted
!+
!-----------------------------------------------------------------------
pure logical function was_accreted(iexternalforce,hi)
 integer, intent(in) :: iexternalforce
 real,    intent(in) :: hi

 ! An accreted particle is indicated by h < 0.
 ! Note less than, but not equal.
 ! (h=0 indicates dead MPI particle)
 was_accreted = (hi < 0.)

end function was_accreted

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_externalforces(iunit,iexternalforce)
 use metric,       only:imetric
 use metric_tools, only:imet_minkowski
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
 use dump_utils,        only:dump_h,add_to_rheader
 integer,      intent(in)    :: iexternalforce
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time
 integer,      intent(out)   :: ierr

 ! call write_headeropts_gwinspiral(hdr,ierr)

end subroutine write_headeropts_extern

!-----------------------------------------------------------------------
!+
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_extern(iexternalforce,hdr,ierr)
 use dump_utils,        only:dump_h,extract
 integer,      intent(in)  :: iexternalforce
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr = 0
 !call read_headeropts_gwinspiral(hdr,ierr)

end subroutine read_headeropts_extern

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_externalforces(name,valstring,imatch,igotall,ierr,iexternalforce)
 use io,           only:fatal,warn
 use metric,       only:imetric
 use metric_tools, only:imet_minkowski
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
 use io,                   only:error
 use units,                only:umass,utime,udist
 use physcon,              only:gg,c
 integer, intent(in)  :: iexternalforce
 integer, intent(out) :: ierr
 real(kind=8) :: gcode, ccode

 ierr = 0

 !
 !--check that G=1 in code units
 !
 gcode = gg*umass*utime**2/udist**3
 if (abs(gcode-1.) > 1.e-10) then
    call error('units',trim(externalforcetype(iexternalforce))//&
               ' external force assumes G=1 in code units but we have',var='G',val=real(gcode))
    ierr = ierr + 1
 endif

 !
 !--check that c=1 in code units
 !
 ccode = c*utime/udist
 if (abs(ccode-1.) > 1.e-10) then
    call error('units',trim(externalforcetype(iexternalforce))//&
                  ' external force assumes c=1 in code units but we have',var='c',val=real(ccode))
    ierr = ierr + 1
 endif

end subroutine initialise_externalforces

end module externalforces
