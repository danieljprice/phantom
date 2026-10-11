!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Moddump to convert sink particles into resolved gaseous spheres or stars
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, moddump_utils, mpidomain, part,
!   setstar
!
 use eos,     only:ieos,gamma,X_in,Z_in,use_var_comp,polyk
 use setstar, only:set_stars,shift_stars,set_defaults_stars,star_t,write_options_stars,read_options_stars
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''
 type(star_t), allocatable :: stars(:)
 logical :: relax = .true., write_rho_to_file = .false.

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:nptmass,xyzmh_ptmass,vxyz_ptmass,eos_vars,rad,hfact
 use io,            only:id,master
 use mpidomain,     only:i_belong
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, allocatable :: xyzmh_ptmass_in(:,:),vxyz_ptmass_in(:,:)
 integer(kind=8) :: npart_total
 integer :: ierr,nstars
 real    :: rhozero

 call set_defaults()
 call get_moddump_options(read_options_moddump,write_options_moddump)

 nstars = nptmass
 nptmass = 0
 !
 !--allocate temporary arrays and copy existing ptmass arrays into them
 !
 allocate(xyzmh_ptmass_in, source=xyzmh_ptmass)
 allocate(vxyz_ptmass_in, source=vxyz_ptmass)
 !
 !--setup and relax stars as needed
 !
 polyk = 0.
 call set_stars(id,master,nstars,stars,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,&
                X_in,Z_in,relax,use_var_comp,write_rho_to_file,&
                rhozero,npart_total,i_belong,ierr)
 !
 !--place stars into orbit, or keep as sink particles if iprofile=0
 !
 call shift_stars(nstars,stars,xyzmh_ptmass_in,vxyz_ptmass_in,xyzh,vxyzu,&
                  xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,nptmass)

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 integer, intent(in) :: iunit

 write(iunit,"(a)") '# parameter file for sink2gas moddump'
 call write_options_stars(stars,relax,write_rho_to_file,ieos,iunit,nstar=size(stars))

end subroutine write_options_moddump

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_options_stars(stars,ieos,relax,write_rho_to_file,db,nerr)

end subroutine read_options_moddump

!---------------------------------------------------------------------
!+
!  set default options for each star from sink particles in the dump
!+
!---------------------------------------------------------------------
subroutine set_defaults()
 use part, only:nptmass,xyzmh_ptmass
 use io, only:fatal
 integer :: i

 if (nptmass <= 0) then
    call fatal('moddump','no sink particles present in file')
 endif
 !
 ! allocate blank options templates for each sink particle
 !
 allocate(stars(nptmass))
 call set_defaults_stars(stars)
 !
 ! fill in the mass and accretion radius for each body from sink
 ! particles already present. Also set the default option to iprofile=0
 ! which just preserves the body as a sink particle
 !
 stars(:)%iprofile = 0
 do i=1,nptmass
    print*,'sink ',i,'m = ',xyzmh_ptmass(4,i),' h = ',xyzmh_ptmass(5,i)
    write(stars(i)%m,"(es20.10)") xyzmh_ptmass(4,i)
    write(stars(i)%hacc,"(es20.10)") xyzmh_ptmass(5,i)
 enddo
end subroutine set_defaults

end module moddump
