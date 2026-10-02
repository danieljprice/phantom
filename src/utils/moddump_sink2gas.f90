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
! :Dependencies: eos, infile_utils, io, moddump_utils, mpidomain, part,
!   setstar
!
 use eos,     only:ieos,gamma,X_in,Z_in,use_var_comp,polyk
 use setstar, only:set_stars,shift_stars,set_defaults_stars,star_t,write_options_stars,read_options_stars
 use moddump_utils, only:prompt_for_params,write_moddump_header
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''
 type(star_t), allocatable :: stars(:)
 logical :: relax = .true., write_rho_to_file = .false.

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .false.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,eos_vars,rad,hfact
 use io,           only:id,master
 use mpidomain,    only:i_belong
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, allocatable :: xyzmh_ptmass_in(:,:),vxyz_ptmass_in(:,:)
 integer(kind=8) :: npart_total
 integer :: ierr,nstars
 real    :: rhozero
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

!
!---Read/write moddump file------------------------------------------------
!
subroutine write_moddump(filename)

 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(a)") '# parameter file for sink2gas moddump'
 call write_options_stars(stars,relax,write_rho_to_file,ieos,iunit,nstar=size(stars))
 close(iunit)

end subroutine write_moddump

subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_options_stars(stars,ieos,relax,write_rho_to_file,db,nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine init_moddump()
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
end subroutine init_moddump

end module moddump
