!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of two stars or sink particles in a binary going into RLOF
!
! :References: Juarez-Garcia et al 2025 , Jackson et al 2017 ApJ 835,145
!
! :Owner: Ana Lourdes Juarez
!
! :Runtime parameters:
!   - a       : *semi-major axis*
!   - gamma   : *adiabatic index*
!   - gastemp : *surface temperature of the donor star in K*
!   - hacc    : *accretion radius of the companion star*
!   - hdon   : *accretion radius of the donor star*
!   - macc    : *mass of the companion star*
!   - mdon    : *mass of the donor star*
!   - mdot    : *mass transfer rate given by MESA in solar mass / yr*
!   - pmass   : *particle mass in code units*
!
! :Dependencies: centreofmass, eos, extern_corotate, externalforces,
!   infile_utils, inject, io, options, part, partinject, physcon,
!   setbinary, setunits, timestep, units
!

 use inject, only:init_inject,lattice_type,wind_radius,wind_injection_x,&
                  mach,v_inf,mdot_msun_yr,filemesa,use_mesa_file

 implicit none
 public :: setpart
 
 private 
 real    :: a,mdon,hdon,macc,hacc,mdot,pmass
 integer :: sink_off
 real :: gastemp = 3000.

contains

!----------------------------------------------------------------
!+
!  setup for binary star simulations (with or without gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,            only:ihsoft,igas,nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc
 use setbinary,       only:set_binary,get_period_from_a
 use centreofmass,    only:reset_centreofmass
 use options,         only:iexternalforce
 use units,           only:set_units,umass,utime
 use externalforces,  only:iext_corotate,omega_corotate
 use extern_corotate, only:icompanion_grav,companion_xpos,companion_mass,hsoft,&
                           primarycore_xpos,primarycore_mass,primarycore_hsoft
 use physcon,         only:solarm,solarr,pi,gg,years
 use io,              only:master,fatal
 use eos,             only:ieos, gmw
 use setunits,        only:mass_unit,dist_unit
 use timestep,        only:tmax
! use readwrite_mesa,  only:read_masstransferrate
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer :: ierr,np
 logical :: iexist
 real    :: period,ecc,mass_ratio
 real    :: XL1,rad_inj,rho_l1,vel_l1,mach_l1,mdot_code
!
!--general parameters
!
 dist_unit = 'solarr'
 mass_unit = 'solarm'
 time = 0.
 polyk = 0.
 gamma = 5./3.
!
!--space available for injected gas particles
!  in case only sink particles are used
!

 iexternalforce = iext_corotate
 sink_off = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0
 a    = 266.34
 mdon = 6.97
 macc = 1.41
 mdot = 1.e-4
 hdon = 100.
 hacc = 1.
 ieos = 2
 gmw  = 0.6
 ecc  = 0.
 pmass = 1.e-8

 filemesa = 'test_data.txt'

 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the shooting particles at a star setup'

 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif

 period = get_period_from_a(mdon,macc,a)
 print*,' period is ',period*utime/years,' yrs'
 tmax = 10.*period
 !dtmax = tmax/200.  !avoid resetting

 ! default value for particle mass based on default mdot
 mdot_code = mdot*(solarm/years)/(umass/utime)
 print*,' suggested pmass for 10,000 particles at end of simulation = ',mdot_code*tmax/10000.
 !
 !--now setup orbit using fake sink particles
 !
 call set_binary(mdon,macc,a,ecc,hdon,hacc,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                  verbose=(id==master))

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)


 if (ierr /= 0) call fatal ('setup_binary','error in call to set_binary')

 massoftype(igas) = pmass

 if (ierr /= 0) call fatal('windtunnel','errors in setup parameters')

 call L1(xyzmh_ptmass,vxyz_ptmass,mdot_code,rad_inj,XL1,rho_l1,vel_l1,mach_l1)

 if (sink_off == 0) then
    icompanion_grav = 0
 elseif (sink_off == 1) then
    icompanion_grav = 2
    companion_mass = xyzmh_ptmass(4,2)
    companion_xpos = xyzmh_ptmass(1,2)
    primarycore_mass = xyzmh_ptmass(4,1)
    primarycore_xpos = xyzmh_ptmass(1,1)
    mass_ratio = mdon / macc
    primarycore_hsoft = hdon !0.1 *  a * 0.49 * mass_ratio**(2./3.) / (0.6*mass_ratio**(2./3.) + &
                  !log( 1. + mass_ratio**(1./3.) ) )
    hsoft = hacc !0.1 *  a * 0.49 * mass_ratio**(-2./3.) / (0.6*mass_ratio**(-2./3.) + &
                 ! log( 1. + mass_ratio**(-1./3.) ) )
    nptmass = 0 !--delete sinks
 endif

 ! Wind parameters (see inject_windtunnel module)
 v_inf = vel_l1
 mach = mach_l1

 ! Wind injection settings
 wind_radius = rad_inj ! in code units
 wind_injection_x = XL1    ! in code units

 npart = 0
 np=0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.


end subroutine setpart

!----------------------------------------------------------------
!+
!  Roche lobe properties
!+
!----------------------------------------------------------------
subroutine L1(xyzmh_ptmass,vxyz_ptmass,mdot_l1,rad_l1,XL1,rho_l1,vel_l1,mach_l1)
 use physcon,    only:pi,twopi,solarm,years,gg,kboltz,mass_proton_cgs
 use units,      only:unit_velocity
 use partinject, only:add_or_update_particle
 use setbinary,  only:L1_point
 use eos,        only:gmw
 use part,       only:igas
 use io,         only:fatal
 real, intent (in)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),mdot_l1
 real, intent (out) :: XL1,rho_l1,vel_l1,mach_l1,rad_l1

 real :: m1,m2,q,radL1,A,Porb,r12
 real :: xyzL1(3),dr(3),x1(3),x2(3),v1(3),v2(3)
 real :: lsutime,cs,u_part,b1,b2,dy,dz,omega,mtot

 x1 = xyzmh_ptmass(1:3,1)
 x2 = xyzmh_ptmass(1:3,2)
 v1 = vxyz_ptmass(1:3,1)
 v2 = vxyz_ptmass(1:3,2)
 m1 = xyzmh_ptmass(4,1)
 m2 = xyzmh_ptmass(4,2)
 q  = m2/m1
 radL1 = L1_point(1./q)   ! find L1 point given binary mass ratio
 dr = x2 - x1
 r12 = sqrt(dot_product(dr, dr))
 xyzL1(1:3) = radL1*dr(:)
 mtot = m1 + m2
 lsutime = sqrt(r12**3/mtot)
 Porb    = twopi * lsutime
 omega = sqrt(mtot/r12**3)

 cs = sqrt(gastemp*kboltz/(gmw*mass_proton_cgs))/unit_velocity !Isothermal sound speed in code units
 u_part = 1.5*cs**2

 b1 = 2.*3.**(2./3.)
 b2 = 0.25*b1 - 2.
 A = 4. + b1/(b2+q**(1./3.)+q**(-1./3.)) !Eq (10) in Jackson 2017
 dy = (sqrt(2.)*cs)/(sqrt(A-1.)*omega)   !Eq (8) in Jackson 2017
 dz = (sqrt(2.)*cs)/(sqrt(A)*omega)      !Eq (9) in Jackson 2017
 rad_l1 = sqrt(dy*dz)

 print*, 'dy =', dy, 'dz =', dz, 'rad_l1 =', rad_l1, 'cs/omega =', cs/omega, 'sqrt(A) =', sqrt(A), 'period =', 2.*pi/omega

 xyzL1(1:3) = radL1*dr(:) + x1
 XL1 = xyzL1(1)

 mach_l1 = 0.1
 vel_l1 = mach_l1*cs
 rho_l1 = mdot_l1/(pi*rad_l1**2*vel_l1)

end subroutine L1


!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setunits,     only:write_options_units
 use eos,          only:gamma
 use setunits,      only:write_options_units
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'

 call write_options_units(iunit)

 write(iunit,"(/,a)") '# orbit settings'
 call write_inopt(a,'a','semi-major axis',iunit)
 call write_inopt(mdon,'mdon','mass of the donor star',iunit)
 call write_inopt(hdon,'hdon','accretion radius of the donor star',iunit)
 call write_inopt(macc,'macc','mass of the companion star',iunit)
 call write_inopt(hacc,'hacc','accretion radius of the companion star',iunit)

 write(iunit,"(/,a)") '# Replace sink particles with fixed gravitational potentials'
 call write_inopt(sink_off,'sink_off','0 = both stars are sink particles, 1 = both stars are fixed gravitational potentials',iunit)

 write(iunit,"(/,a)") '# mass resolution'
 call write_inopt(pmass,'pmass','particle mass in code units',iunit)

 write(iunit,"(/,a)") '# mass injection settings'
 call write_inopt(use_mesa_file,'use_mesa_file','use_mesa_file',iunit)
 if (use_mesa_file) then
    call write_inopt(filemesa,'filemesa','mesa file path',iunit)
 else
    call write_inopt(mdot_msun_yr,'mdot','mass transfer rate in solar mass / yr',iunit)
 endif
 call write_inopt(gastemp,'gastemp','surface temperature of the donor star in K',iunit)
 call write_inopt(gamma,'gamma','adiabatic index',iunit)

 close(iunit)

end subroutine write_setupfile


!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,  only:open_db_from_file,inopts,read_inopt,close_db
 use io,            only:error,fatal
 use setunits,      only:read_options_and_set_units
 use eos,           only:gamma
 character(len=*), intent(in) :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0

 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr)

 call read_inopt(a,'a',db,errcount=nerr)
 call read_inopt(mdon,'mdon',db,errcount=nerr)
 call read_inopt(hdon,'hdon',db,errcount=nerr)
 call read_inopt(macc,'macc',db,errcount=nerr)
 call read_inopt(hacc,'hacc',db,errcount=nerr)

 call read_inopt(sink_off,'sink_off',db,errcount=nerr)

 call read_inopt(pmass,'pmass',db,errcount=nerr)
 call read_inopt(gastemp,'gastemp',db,errcount=nerr)
 call read_inopt(gamma,'gamma',db,errcount=nerr)
 call read_inopt(use_mesa_file,'use_mesa_file',db,errcount=nerr)

 if (use_mesa_file) then
    call read_inopt(filemesa,'filemesa',db,errcount=nerr)
 else
    call read_inopt(mdot_msun_yr,'mdot',db,errcount=nerr)
 endif

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

 call close_db(db)
end subroutine read_setupfile

end module setup
