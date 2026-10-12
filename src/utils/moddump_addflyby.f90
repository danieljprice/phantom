!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! moddump to add a single flyby with general orbital parameters
!
! :References: None
!
! :Owner: Arcelia Hermosillo Ruiz
!
! :Runtime parameters:
!   - accr2   : *accretion radius of secondary*
!   - deltat  : *output interval as fraction of binary period*
!   - m2      : *mass of secondary (in code units)*
!   - norbits : *maximum number of binary orbits*
!
! :Dependencies: centreofmass, dim, io, moddump_utils, part,
!   prompting, setorbit, timestep
!

 use setorbit,      only:orbit_t
 implicit none

 type(orbit_t) :: orbit
 real    :: deltat,accr2,m2,m1
 integer :: norbits

 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:nsinkproperties
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc
 use centreofmass,      only:reset_centreofmass
 use setorbit,          only:set_orbit,get_orbital_time
 use io,             only:id,master,fatal
 use moddump_utils,  only:get_moddump_options
 use timestep,       only:tmax,dtmax

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: hacc1,period
 integer :: nptmass_in,ierr
 real :: xyzmh_ptmass_in(nsinkproperties,2),vxyz_ptmass_in(3,2)

 if (nptmass > 0) then
    m1 = xyzmh_ptmass(4,1)
    hacc1 = xyzmh_ptmass(ihacc,1)
    if (id==master) print "(a,es10.3)", ' Using first sink particle as primary with M = ',m1
 else
    call fatal('moddump','no sink particles in dump file')
 endif

 !--read parameters from the .moddump file, prompting if it does not exist
 call set_defaults()
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 xyzmh_ptmass_in = 0.
 vxyz_ptmass_in = 0.
 nptmass_in = 0
 call set_orbit(orbit,m1,m2,hacc1,accr2,xyzmh_ptmass_in,vxyz_ptmass_in,&
                        nptmass_in,(id==master),ierr)

 nptmass = nptmass + 1
 if (nptmass > size(xyzmh_ptmass, 2)) then
    call fatal('moddump', 'Not enough space in xyzmh_ptmass for another sink particle')
 endif

! translate the secondary so that it is correctly set relative to the primary
 xyzmh_ptmass(:,nptmass) = xyzmh_ptmass_in(:,2)

 xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,nptmass) &
                          + xyzmh_ptmass(1:3,1) &
                          - xyzmh_ptmass_in(1:3,1)

 vxyz_ptmass(:,nptmass) = vxyz_ptmass_in(:,2) &
                       + vxyz_ptmass(:,1) &
                       - vxyz_ptmass_in(:,1)

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 ! set .in parameters below

 !--time of flyby
 call get_orbital_time(orbit,m1,m2,period)

 if (period > 0.) then
    if (deltat > 0.) dtmax = deltat*period
    if (norbits >= 0) tmax = norbits*period
 endif

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Do you want to specify the flyby orbit as bound (elliptic) or'// &
                  ' unbound (parabolic/hyperbolic) or as observed dx,dv?'//new_line('A')// &
                  ' 0=bound'//new_line('A')//' 1=unbound'//new_line('A')// &
                  ' 2=orbit reconstructor'//new_line('A')// '3=observed dx,dv'//new_line('A'),orbit%input_type,0,4)
 select case (orbit%input_type)
 case (0)
    !--bound
    m2       = 0.2
    orbit%elems%a = '10.'
    accr2    = 1.0
 case default
    !--unbound (flyby)
    m2       = 1.
    accr2    = 1.
    !
    ! the following is only if we want to override defaults in set_defaults_orbit
    ! so for input_type=2 or 3 we will just get those defaults
    !
    if (orbit%input_type >= 1) then
       orbit%flyby%rp = '200.'
       orbit%flyby%d = '2000.'
    endif
    orbit%e = 2.0
 end select
end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 use setorbit,      only:write_options_orbit
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# perturber parameters'
 call write_inopt(m2,'m2','mass of secondary (in code units)',iunit)
 call write_inopt(accr2,'accr2','accretion radius of secondary',iunit)
 call write_options_orbit(orbit,iunit)

 write(iunit,"(/,a)") '# timestepping'
 call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
 call write_inopt(deltat,'deltat','output interval as fraction of binary period',iunit)

end subroutine write_options_moddump

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 use setorbit,      only:read_options_orbit
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(m2,'m2',db,errcount=nerr,min=0.)
 call read_inopt(accr2,'accr2',db,errcount=nerr,min=0.)
 call read_options_orbit(orbit,m1,m2,db,nerr)
 call read_inopt(norbits,'norbits',db,errcount=nerr)
 call read_inopt(deltat,'deltat',db,errcount=nerr,default=0.1)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  set default values of moddump parameters
!+
!----------------------------------------------------------------
subroutine set_defaults()
 use setorbit, only:set_defaults_orbit

 call set_defaults_orbit(orbit)
 orbit%input_type = 1
 m2 = 0.1
 accr2 = 1.0
 norbits = 100
 deltat = 0.1

end subroutine set_defaults

end module moddump

