!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Dependencies: eos, infile_utils, io, mpidomain, part, setstar
!
 use eos,     only:ieos,gamma,X_in,Z_in,use_var_comp,polyk
 use setstar, only:set_stars,shift_stars,set_defaults_stars,star_t,write_options_stars,read_options_stars
 implicit none
 type(star_t), allocatable :: stars(:)
 logical :: relax = .true., write_rho_to_file = .false.

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,eos_vars,rad,hfact
 use io,        only:fatal,id,master,error
 use mpidomain, only:i_belong
 integer,  intent(inout) :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: massoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=120)      :: filename
 real, allocatable :: xyzmh_ptmass_in(:,:),vxyz_ptmass_in(:,:)
 integer(kind=8) :: npart_total
 integer :: ierr,nstars,i
 logical :: iexist
 real    :: rhozero
 !
 ! check there are sink particles present
 !
 if (nptmass <= 0) then
    call fatal('moddump','no sink particles present in file')
    return
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
 !
 ! read the parameter file if it exists
 !
 filename = 'sink2gas.moddump'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,nptmass,ieos,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename,nptmass)
       print*,' Edit '//trim(filename)//' and rerun phantommoddump'
    endif
    stop
 endif

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
 call shift_stars(nstars,stars,xyzmh_ptmass_in(1:3,:),vxyz_ptmass_in(1:3,:),xyzh,vxyzu,&
                  xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,nptmass)

end subroutine modify_dump

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename,nstars)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer,          intent(in) :: nstars
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# parameter file for sink2star moddump'
 call write_options_stars(stars,relax,write_rho_to_file,ieos,iunit,nstar=nstars)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,nstars,ieos,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 integer,          intent(in)    :: nstars
 integer,          intent(inout) :: ieos
 character(len=*), intent(in)    :: filename
 integer,          intent(out)   :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_stars(stars,ieos,relax,write_rho_to_file,db,nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of moddump parameter file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module moddump
