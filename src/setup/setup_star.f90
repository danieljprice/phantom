!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Runtime parameters:
!   - EOSopt            : *EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)*
!   - X                 : *hydrogen mass fraction*
!   - gamma             : *Adiabatic index*
!   - ieos              : *1=isothermal,2=adiabatic,10=MESA,12=idealplusrad*
!   - initialtemp       : *initial temperature of the star*
!   - irecomb           : *Species to include in recombination (0: H2+H+He, 1:H+He, 2:He*
!   - metallicity       : *metallicity*
!   - mu                : *mean molecular weight*
!   - polyk             : *polytropic constant (cs^2 if isothermal)*
!   - relax_star        : *relax star(s) automatically during setup*
!   - use_var_comp      : *Use variable composition (X, Z, mu)*
!   - write_rho_to_file : *write density profile(s) to file*
!
! :Dependencies: apr, dim, eos, eos_gasradrec, eos_piecewise,
!   extern_densprofile, externalforces, infile_utils, io, kernel,
!   mpidomain, mpiutils, options, part, physcon, prompting, relaxstar,
!   setstar, setunits, setup_params, timestep, units
!
 use io,             only:fatal,error,warning,master
 use part,           only:gravity,gr
 use physcon,        only:solarm,solarr,km,pi,c,kb_on_mh,radconst
 use options,        only:nfulldump,iexternalforce,calc_erot,use_var_comp
 use timestep,       only:tmax,dtmax
 use eos,                only:ieos
 use externalforces,     only:iext_densprofile
 use extern_densprofile, only:nrhotab
 use setstar,            only:ibpwpoly,ievrard,imesa,star_t,need_polyk
 use setunits,           only:dist_unit,mass_unit
 implicit none
 !
 ! Input parameters
 !
 integer            :: EOSopt
 integer            :: need_iso
 real               :: maxvxyzu
 logical            :: iexist
 logical            :: relax_star_in_setup,write_rho_to_file
 type(star_t)       :: star

 public             :: setpart
 private

contains

!-----------------------------------------------------------------------
!+
!  Setup routine for stars / spherical collapse calculations
!+
!-----------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use units,           only:set_units,select_unit
 use kernel,          only:hfact_default
 use eos,             only:init_eos,finish_eos,gmw,X_in,Z_in
 use eos_piecewise,   only:init_eos_piecewise_preset
 use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,eos_vars,rad
 use mpiutils,        only:reduceall_mpi
 use mpidomain,       only:i_belong
 use setup_params,    only:rhozero,npart_total
 use setstar,         only:set_star
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
 !
 ! Initialise parameters, including those that will not be included in *.setup
 !
 time         = 0.
 polyk        = 1.0
 gamma        = 5./3.
 hfact        = hfact_default
 maxvxyzu     = size(vxyzu(:,1))
 relax_star_in_setup = .false.
 write_rho_to_file = .false.

 !
 ! set default options
 !
 dist_unit   = 'solarr'
 mass_unit   = 'solarm'
 EOSopt      = 1
 gmw         = 0.5988
 X_in        = 0.74
 Z_in        = 0.02
 use_var_comp = .false.
 !
 ! defaults needed for error checking
 !
 need_iso = 0 ! -1 = no; 0 = doesn't matter; 1 = yes
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
    call read_setupfile(setupfile,gamma,polyk,need_iso,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(setupfile,gamma,polyk)
       stop 'please rerun phantomsetup with revised .setup file'
    endif
    !--Prompt to get inputs and write to file
 elseif (id==master) then
    print "(a,/)",trim(setupfile)//' not found: using interactive setup'
    call setup_interactive(polyk,gamma,iexist,id,master,ierr)
    call write_setupfile(setupfile,gamma,polyk)
    stop 'please check and edit .setup file and rerun phantomsetup'
 endif

 !
 ! Verify correct pre-processor commands
 !
 if (.not.gravity) then
    iexternalforce = iext_densprofile
    write_rho_to_file = .true.
 endif

 if (maxvxyzu > 3  .and. need_iso == 1) call fatal('setup','require ISOTHERMAL=yes')
 if (maxvxyzu < 4  .and. need_iso ==-1) call fatal('setup','require ISOTHERMAL=no')
 !
 ! initialise the equation of state
 !
 if (ieos==9) call init_eos_piecewise_preset(EOSopt)
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('setup','could not initialise equation of state')
 !
 ! set up particles
 !
 npartoftype(:) = 0
 npart          = 0
 nptmass        = 0
 vxyzu          = 0.0
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
               massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
               X_in,Z_in,relax_star_in_setup,use_var_comp,write_rho_to_file,&
               rhozero,npart_total,i_belong,ierr)
 !
 ! finish/deallocate equation of state tables
 !
 call finish_eos(ieos,ierr)

 !
 ! override some default settings in the .in file for some cases
 !
 select case(star%iprofile)
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
subroutine setup_interactive(polyk,gamma,iexist,id,master,ierr)
 use prompting,     only:prompt
 use units,         only:select_unit
 use eos,           only:X_in,Z_in,gmw
 use eos_gasradrec, only:irecomb
 use setstar,       only:set_star_interactive
 use setunits,      only:set_units_interactive
 real, intent(out)    :: polyk,gamma
 logical, intent(in)  :: iexist
 integer, intent(in)  :: id,master
 integer, intent(out) :: ierr

 ierr = 0

 ! units
 call set_units_interactive(gr)

 ! star
 call set_star_interactive(id,master,star,need_iso,use_var_comp,ieos,polyk)

 ! equation of state
 call prompt('Enter the desired EoS (1=isothermal,2=adiabatic,10=MESA,12=idealplusrad)',ieos)
 select case(ieos)
 case(15) ! Helmholtz
    call prompt('Enter temperature',star%initialtemp,1.0e3,1.0e11)
 case(9)
    write(*,'(a)') 'EOS options: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)'
    call prompt('Enter equation of state type',EOSopt,1,4)
 case(2)
    call prompt('Enter gamma (adiabatic index)',gamma,1.,7.)
 case(20)
    call prompt('Enter irecomb (0: H2+H+He, 1:H+He, 2:He)',irecomb,0)
 end select

 if (need_polyk(star%iprofile)) then
    call prompt('Enter polytropic constant (cs^2 if isothermal)',polyk,0.)
 endif

 if ((.not. use_var_comp) .and. (star%isoftcore<=0)) then
    if ( (ieos==12) .or. (ieos==2) ) call prompt('Enter mean molecular weight',gmw,0.)
    if ( (ieos==10) .or. (ieos==20) ) then
       call prompt('Enter hydrogen mass fraction (X)',X_in,0.,1.)
       call prompt('Enter metals mass fraction (Z)',Z_in,0.,1.)
    endif
 endif

 call prompt('Relax star automatically during setup?',relax_star_in_setup)

end subroutine setup_interactive

!-----------------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename,gamma,polyk)
 use infile_utils,  only:write_inopt
 use dim,           only:tagline,use_apr
 use relaxstar,     only:write_options_relax
 use eos,           only:X_in,Z_in,gmw
 use eos_gasradrec, only:irecomb
 use setstar,       only:write_options_star,need_polyk
 use setunits,      only:write_options_units
 use apr,           only:write_options_apr
 real,             intent(in) :: gamma,polyk
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(a)") '# input file for Phantom star setup'

 call write_options_units(iunit,gr)
 call write_options_star(star,iunit)

 write(iunit,"(/,a)") '# equation of state'
 call write_inopt(ieos,'ieos','1=isothermal,2=adiabatic,10=MESA,12=idealplusrad',iunit)

 if (star%iprofile==imesa) then
    call write_inopt(use_var_comp,'use_var_comp','Use variable composition (X, Z, mu)',iunit)
 endif

 select case(ieos)
 case(15) ! Helmholtz
    call write_inopt(star%initialtemp,'initialtemp','initial temperature of the star',iunit)
 case(9)
    write(iunit,"(/,a)") '# Piecewise Polytrope default options'
    call write_inopt(EOSopt,'EOSopt','EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)',iunit)
 case(2)
    call write_inopt(gamma,'gamma','Adiabatic index',iunit)
    if ((star%isoftcore<=0) .and. (.not. use_var_comp)) call write_inopt(gmw,'mu','mean molecular weight',iunit)
 case(10,20)
    if (ieos==20) call write_inopt(irecomb,'irecomb','Species to include in recombination (0: H2+H+He, 1:H+He, 2:He',iunit)
    if ( (.not. use_var_comp) .and. (star%isoftcore <= 0) ) then
       call write_inopt(X_in,'X','hydrogen mass fraction',iunit)
       call write_inopt(Z_in,'Z','metallicity',iunit)
    endif
 case(12)
    call write_inopt(gamma,'gamma','Adiabatic index',iunit)
    if ((star%isoftcore<=0) .and. (.not. use_var_comp)) call write_inopt(gmw,'mu','mean molecular weight',iunit)
 end select

 if (need_polyk(star%iprofile)) call write_inopt(polyk,'polyk','polytropic constant (cs^2 if isothermal)',iunit)

 write(iunit,"(/,a)") '# relaxation options'
 call write_inopt(relax_star_in_setup,'relax_star','relax star(s) automatically during setup',iunit)
 if (relax_star_in_setup) call write_options_relax(iunit)

 call write_inopt(write_rho_to_file,'write_rho_to_file','write density profile(s) to file',iunit)

 if (use_apr) call write_options_apr(iunit)

 close(iunit)

end subroutine write_setupfile
!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,gamma,polyk,need_iso,ierr)
 use infile_utils,  only:open_db_from_file,inopts,close_db,read_inopt
 use io,            only:error
 use units,         only:select_unit
 use relaxstar,     only:read_options_relax
 use eos,           only:X_in,Z_in,gmw
 use eos_gasradrec, only:irecomb
 use setstar,       only:read_options_star
 use setunits,      only:read_options_and_set_units
 use apr,           only:apr_max_in,ref_dir,apr_type,apr_rad,apr_drad
 use dim,           only:use_apr
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: need_iso,ierr
 real,             intent(out) :: gamma,polyk
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return

 nerr = 0

 ! units
 call read_options_and_set_units(db,nerr,gr)

 ! star options
 call read_options_star(star,need_iso,ieos,polyk,db,nerr)

 ! equation of state
 call read_inopt(ieos,'ieos',db,errcount=nerr)
 if (star%iprofile==imesa) call read_inopt(use_var_comp,'use_var_comp',db,errcount=nerr)

 select case(ieos)
 case(15) ! Helmholtz
    call read_inopt(star%initialtemp,'initialtemp',db,errcount=nerr)
 case(9)
    call read_inopt(EOSopt,'EOSopt',db,errcount=nerr)
 case(2)
    call read_inopt(gamma,'gamma',db,errcount=nerr)
    if ( (.not. use_var_comp) .and. (star%isoftcore <= 0)) call read_inopt(gmw,'mu',db,errcount=nerr)
 case(10,20)
    if (ieos==20) call read_inopt(irecomb,'irecomb',db,errcount=nerr)
    ! if softening stellar core, composition is automatically determined at R/2
    if ( (.not. use_var_comp) .and. (star%isoftcore <= 0)) then
       call read_inopt(X_in,'X',db,errcount=nerr)
       call read_inopt(Z_in,'Z',db,errcount=nerr)
    endif
 case(12)
    ! if softening stellar core, mu is automatically determined at R/2
    call read_inopt(gamma,'gamma',db,errcount=nerr)
    if ( (.not. use_var_comp) .and. (star%isoftcore <= 0)) call read_inopt(gmw,'mu',db,errcount=nerr)
 end select

 if (need_polyk(star%iprofile)) call read_inopt(polyk,'polyk',db,errcount=nerr)

 ! relax star options
 call read_inopt(relax_star_in_setup,'relax_star',db,errcount=nerr)
 if (relax_star_in_setup) call read_options_relax(db,nerr)
 if (nerr /= 0) ierr = ierr + 1

 ! option to write density profile to file
 call read_inopt(write_rho_to_file,'write_rho_to_file',db)

 if (use_apr) then
    call read_inopt(apr_max_in,'apr_max',db,errcount=nerr)
    call read_inopt(ref_dir,'ref_dir',db,errcount=nerr)
    call read_inopt(apr_type,'apr_type',db,errcount=nerr)
    call read_inopt(apr_rad,'apr_rad',db,errcount=nerr)
    call read_inopt(apr_drad,'apr_drad',db,errcount=nerr)
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_star: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile

end module setup
