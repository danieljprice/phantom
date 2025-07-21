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
! :Dependencies: apr, apr_region, dim, eos, externalforces, infile_utils,
!   io, kernel, mpidomain, options, part, physcon, setstar, setunits,
!   setup_params, timestep
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
 use apr,             only:use_apr
 use infile_utils,    only:get_options,infile_exists
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
 if (.not. infile_exists(fileprefix)) then
    tmax  = 100.
    dtmax = 1.0
    ieos  = 2
 endif
 !
 ! read/write from .setup file
 !
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 !
 ! Verify correct pre-processor commands
 !
 if (.not.gravity) then
    iexternalforce = iext_densprofile
 endif
 write_rho_to_file = .true.

 !
 ! if apr is being used, read in the parameters here or assume the defaults
 ! note that this needs to be done before the particles are set up
 !
 if (use_apr .and. relax_star_in_setup) then
    if (infile_exists(fileprefix)) then
       call read_aprsetupfile(trim(fileprefix)//'.in',ierr)
    else
       call warning('setup_star','apr options needed for relaxation not found; making you a .in file, update and try again')
       relax_star_in_setup = .false.
    endif
 endif

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
subroutine setup_interactive()
 use setstar,  only:set_stars_interactive
 use setunits, only:set_units_interactive

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

!-----------------------------------------------------------------------
!+
!  Read apr parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_aprsetupfile(filename,ierr)
 use infile_utils,  only:open_db_from_file,inopts,close_db,read_inopt
 use setstar,       only:read_options_stars
 use setunits,      only:read_options_and_set_units
 use apr_region,           only:apr_max_in,ref_dir,apr_type,apr_rad,apr_drad
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return

 nerr = 0

 ! apr options
 call read_inopt(apr_max_in,'apr_max',db,nerr)
 call read_inopt(ref_dir,'ref_dir',db,nerr)
 call read_inopt(apr_type,'apr_type',db,nerr)
 call read_inopt(apr_rad,'apr_rad',db,nerr)
 call read_inopt(apr_drad,'apr_drad',db,nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_star: ',nerr,' error(s) during read of apr options from .in file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_aprsetupfile

end module setup
