!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup procedure for stars
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, externalforces, infile_utils, io, kernel,
!   mpidomain, options, part, physcon, setstar, setunits, setup_params,
!   timestep
!
 use io,             only:fatal,error,warning,master
 use part,           only:gravity,gr
 use physcon,        only:solarm,solarr,km,pi,c,radconst
 use options,        only:nfulldump,iexternalforce,calc_erot,use_var_comp
 use timestep,       only:tmax,dtmax
 use eos,                only:ieos
 use externalforces,     only:iext_densprofile
 use setstar,            only:star_t
 use setunits,           only:dist_unit,mass_unit
 implicit none
 !
 ! Input parameters
 !
 real               :: maxvxyzu
 logical            :: iexist
 logical            :: relax_star_in_setup,write_rho_to_file
 type(star_t)       :: star(1)

 public             :: setpart
 private

contains

!-----------------------------------------------------------------------
!+
!  Setup routine for stars / spherical collapse calculations
!+
!-----------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use kernel,          only:hfact_default
 use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,eos_vars,rad
 use eos,             only:X_in,Z_in
 use mpidomain,       only:i_belong
 use setup_params,    only:rhozero,npart_total
 use setstar,         only:set_defaults_stars,set_stars,shift_stars,ibpwpoly,ievrard
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer                          :: ierr
 logical                          :: setexists
 character(len=120)               :: setupfile,inname
 real                             :: x0(3,1),v0(3,1)
 !
 ! Initialise parameters, including those that will not be included in *.setup
 !
 time         = 0.
 polyk        = 1.0
 gamma        = 5./3.
 hfact        = hfact_default
 maxvxyzu     = size(vxyzu(:,1))
 !
 ! set default options
 !
 dist_unit   = 'solarr'
 mass_unit   = 'solarm'
 call set_defaults_stars(star)
 !
 ! determine if the .in file exists
 !
 inname=trim(fileprefix)//'.in'
 inquire(file=inname,exist=iexist)
 if (.not. iexist) then
    tmax  = 100.
    dtmax = 1.0
    ieos  = 2
 endif
 !
 ! determine if the .setup file exists
 !
 setupfile = trim(fileprefix)//'.setup'
 inquire(file=setupfile,exist=setexists)
 if (setexists) then
    call read_setupfile(setupfile,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(setupfile)
       stop 'please rerun phantomsetup with revised .setup file'
    endif
    !--Prompt to get inputs and write to file
 elseif (id==master) then
    print "(a,/)",trim(setupfile)//' not found: using interactive setup'
    call setup_interactive(ieos)
    call write_setupfile(setupfile)
    stop 'please check and edit .setup file and rerun phantomsetup'
 endif

 !
 ! Verify correct pre-processor commands
 !
 if (.not.gravity) then
    iexternalforce = iext_densprofile
 endif
 write_rho_to_file = .true.
 !
 ! set up particles
 !
 npartoftype(:) = 0
 npart          = 0
 nptmass        = 0
 call set_stars(id,master,1,star,xyzh,vxyzu,eos_vars,rad,&
                npart,npartoftype,massoftype,hfact,&
                xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
                relax_star_in_setup,use_var_comp,write_rho_to_file,&
                rhozero,npart_total,i_belong,ierr)
 !
 ! put the star at the origin with zero velocity,
 ! or replace with sink particle
 !
 x0 = 0.
 v0 = 0.
 call shift_stars(1,star,x0,v0,xyzh,vxyzu,&
                  xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,nptmass)
 !
 ! override some default settings in the .in file for some cases
 !
 select case(star(1)%iprofile)
 case(ibpwpoly) ! piecewise polytrope
    calc_erot = .true.
 case(ievrard)  ! Evrard Collapse
    if (.not.iexist) then
       tmax      = 3.0
       dtmax     = 0.1
       nfulldump = 1
    endif
 end select

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Ask questions of the user to determine which setup to use
!+
!-----------------------------------------------------------------------
subroutine setup_interactive(ieos)
 use setstar,       only:set_stars_interactive
 use setunits,      only:set_units_interactive
 integer, intent(inout) :: ieos

 ! units
 call set_units_interactive(gr)

 ! star
 call set_stars_interactive(star,ieos,relax_star_in_setup)

end subroutine setup_interactive

!-----------------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use dim,           only:tagline
 use setstar,       only:write_options_stars
 use setunits,      only:write_options_units
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(a)") '# input file for Phantom star setup'

 call write_options_units(iunit,gr)
 call write_options_stars(star,relax_star_in_setup,write_rho_to_file,ieos,iunit)
 close(iunit)

end subroutine write_setupfile
!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,  only:open_db_from_file,inopts,close_db
 use setstar,       only:read_options_stars
 use setunits,      only:read_options_and_set_units
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return

 nerr = 0

 ! units
 call read_options_and_set_units(db,nerr,gr)

 ! star options
 call read_options_stars(star,ieos,relax_star_in_setup,write_rho_to_file,db,nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_star: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile

end module setup
