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
! :References: None
!
!=======
! :Owner: Ana Lourdes Juarez

!
! :Runtime parameters:
!   - a    : *semi-major axis*
!   - hacc : *accretion radius of the companion star*
!   - macc : *mass of the companion star*
!   - mdon : *mass of the donor star*
!
! :Dependencies: centreofmass, eos, extern_corotate, externalforces,
!   infile_utils, io, options, part, setbinary, setunits, timestep
!

use inject, only:init_inject,nstar,Rstar,lattice_type,handled_layers,&
                  wind_radius,wind_injection_x,wind_length,&
                  rho_inf,pres_inf,v_inf,windonly

 implicit none
 public :: setpart
 real    :: a,mdon,macc,hacc,mdot
 integer :: nstar_in

 real, private :: gastemp = 3000.

 private

contains

!----------------------------------------------------------------
!+
!  setup for binary star simulations (with or without gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,        only:ihsoft,igas
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc
 use setbinary,      only:set_binary,get_period_from_a
 use centreofmass,    only:reset_centreofmass
 use options,        only:iexternalforce
 use units,       only:udist,umass,utime,set_units,unit_velocity,unit_density,unit_pressure
 use externalforces, only:iext_corotate,omega_corotate
 use extern_corotate, only:icompanion_grav,companion_xpos,companion_mass,hsoft
 use physcon,     only:solarm,solarr,pi,gg
 use io,             only:master,fatal
 use eos,            only:ieos, gmw
 use setunits,       only:mass_unit,dist_unit
 use timestep,       only:tmax,dtmax
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
 real    :: period,ecc,hdon,mass_ratio
 real    :: rhocentre,pmass
 real    :: XL1,rad_inj,rho_l1,vel_l1,pr_l1
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
 icompanion_grav = 1
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0
 a    = 266.34
 mdon = 6.97
 macc = 1.41
 Mdot = 1.e-2
 hacc = 1.
 ieos = 2
 gmw  = 0.6
 ecc  = 0.
 hdon = 1.
 nstar_in = 1000

 period = get_period_from_a(mdon,macc,a)
 tmax = 10.*period
 dtmax = tmax/200.

 !
 !--now setup orbit using fake sink particles
 !
 call set_binary(mdon,macc,a,ecc,hdon,hacc,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                  verbose=(id==master))

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)


 if (ierr /= 0) call fatal ('setup_binary','error in call to set_binary')

 pmass = Mdot / real(nstar_in)
 massoftype(igas) = pmass
 call check_setup(pmass,ierr)

 if (ierr /= 0) call fatal('windtunnel','errors in setup parameters')

 call L1(xyzmh_ptmass,vxyz_ptmass,Mdot,pmass,nstar_in,XL1,rad_inj,rho_l1,vel_l1,pr_l1)

 ! Wind parameters (see inject_windtunnel module)
 v_inf    = vel_l1 !/ unit_velocity
 rho_inf  = rho_l1 !/ unit_density
 pres_inf = pr_l1  !/ unit_pressure

 ! Wind injection settings
 lattice_type = 1
 handled_layers = 4
 wind_radius = rad_inj/solarr ! in units of Rstar
 wind_injection_x = XL1/solarr     ! in units of Rstar
 wind_length = rad_inj/solarr    


 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the shooting at a star setup'

 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ieos,polyk,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif
 !
 !
 !--if a is negative or is given time units, interpret this as a period
 !


 ! Initialise particle injection
 call init_inject(ierr)
 npart = 0
 np=0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
! massoftype = 0.

 companion_mass = mdon
 companion_xpos = xyzmh_ptmass(1,1)
 mass_ratio = mdon / macc
 hsoft = 0.1 * 0.49 * mass_ratio**(2./3.) / (0.6*mass_ratio**(2./3.) + &
               log( 1. + mass_ratio**(1./3.) ) ) * a
 !
 !--delete donor sink
 !
 !nptmass=1
 !xyzmh_ptmass(:,1) = xyzmh_ptmass(:,2)
 !vxyz_ptmass(1:3,1) = 0.

 !--restore options
 !


end subroutine setpart

!----------------------------------------------------------------
!+
!  Roche lobe properties
!+
!----------------------------------------------------------------     
subroutine L1(xyzmh_ptmass,vxyz_ptmass,mdot,pmass,nstar_in,XL1,rad_inj,rho_l1,vel_l1,pr_l1)
 use physcon,  only:pi,twopi,solarm,years,gg,kboltz,mass_proton_cgs
 use units,    only:utime,udist,umass
 use partinject,only:add_or_update_particle
 use setbinary, only:L1_point
 use eos,        only:gamma,gmw
 use part,       only:igas
 use io,         only:fatal
 real, intent (in)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),mdot,pmass
 integer, intent (in) :: nstar_in
 real, intent (out) :: XL1,rad_inj,rho_l1,vel_l1,pr_l1

 real :: m1,m2,q,radL1,theta_s,A,mu,Porb,r12,r2L1,smag,vol_l1
 real :: xyzL1(3),dr(3),x1(3),x2(3),x0(3),vxyzL1(3),v1(3),v2(3),xyzinj(3),s(3)
 real :: eps,spd_inject,lm12,lm1,lm32,sw_chi,sw_gamma,U1,lsutime,cs,u_part

 x1 = xyzmh_ptmass(1:3,1)
 x2 = xyzmh_ptmass(1:3,2)
 v1 = vxyz_ptmass(1:3,1)
 v2 = vxyz_ptmass(1:3,2)
 dr = x2 - x1
 r12 = dist(x2,x1)
 r2L1 = dist(xyzL1, x2) 
 m1 = xyzmh_ptmass(4,1)
 m2 = xyzmh_ptmass(4,2)
 q  = m2/m1
 mu = 1./(1 + q)
 radL1      = L1_point(1./q)                     ! find L1 point given binary mass ratio
 XL1        = radL1-(1-mu)
 lsutime = ((m1+m2)/(r12**3))**(-1./2)*utime

 Porb      = twopi * sqrt( (r12*udist)**3 / (gg*(m1+m2)*umass) )
 eps = Porb/(twopi*r12*udist) * (gastemp*kboltz/(gmw*mass_proton_cgs))**0.5
 A  = mu / abs(XL1 - 1. + mu)**3 + (1. - mu)/abs(XL1 + mu)**3! See Lubow & Shu 1975, eq 13
 theta_s = -acos( -4./(3.*A)+(1-8./(9.*A))**0.5)/2.              ! See Lubow & Shu 1975, eq 24
 s =  (/1.,0.,0./)*r2L1*eps/200.0 !(/cos(theta_s),sin(theta_s),0.0/)*r2L1*eps/200.0   ! last factor still a "magic number". Fix. ANA:WHAT??
 smag = sqrt(dot_product(s,s))
 xyzL1(1:3) = XL1*dr(:) 
 xyzinj(1:3) = xyzL1 + s
 vxyzL1 = v1*dist(xyzL1,x0)/dist(x0, x1) ! orbital motion of L1 point
 U1 = 3*A*sin(2*theta_s)/(4*lsutime/utime)
 vel_l1 = abs(U1*dist(xyzinj,xyzL1))  !L&S75 eq 23b
 lm12 = (A - 2. + sqrt(A*(9.*A - 8.)))/2.                        ! See Heerlein+99, eq A8
 lm32 = lm12 - A + 2.                                            ! See Heerlein+99, eq A15
 lm1  = sqrt(A)
 sw_chi = (1e8/udist)**(-2)
 sw_gamma = (1e8/udist)**(-2)
 rad_inj = sw_chi**(-0.5)
 vol_l1 = (4./3.)*pi*rad_inj**3
 rho_l1 = pmass*real(nstar_in)/vol_l1
 u_part = 3.*(kboltz*gastemp/(mu*mass_proton_cgs))/2.
 cs = (u_part*gamma*(gamma-1.))
 pr_l1 = rho_l1*cs**2
end subroutine L1

!
! Function to get the distance between two points
!
real function dist(x1,x2)
 real, intent(in)  :: x1(3), x2(3)
 real :: dr(3)

 dr = x1 - x2
 dist = sqrt(dot_product(dr, dr))
 return
end function dist


!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setunits,     only:write_options_units
 use eos,          only:write_options_eos,gamma
 use units,         only:unit_density,unit_pressure,unit_velocity
 use setunits,      only:write_options_units
 use units,         only:unit_density
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'

 call write_options_units(iunit)
  call write_options_eos(iunit)

 write(iunit,"(/,a)") '# orbit settings'
 call write_inopt(a,'a','semi-major axis',iunit)
 call write_inopt(mdon,'mdon','mass of the donor star',iunit)
 call write_inopt(macc,'macc','mass of the companion star',iunit)
 call write_inopt(hacc,'hacc','accretion radius of the companion star',iunit)

 write(iunit,"(/,a)") '# sphere settings'
 call write_inopt(nstar_in,'nstar','number of particles resolving gas sphere',iunit)  ! note: this is an estimate, actual no. of particles is npart outputted from set_sphere

 write(iunit,"(/,a)") '# wind settings'
 call write_inopt(v_inf*unit_velocity/1.e5,'v_inf','wind speed / km s^-1',iunit)
 call write_inopt(rho_inf*unit_density,'rho_inf','wind density / g cm^-3',iunit)
 call write_inopt(pres_inf*unit_pressure,'pres_inf','wind pressure / dyn cm^2',iunit)
 call write_inopt(gamma,'gamma','adiabatic index',iunit)

 write(iunit,"(/,a)") '# wind injection settings'
 call write_inopt(lattice_type,'lattice_type','0: cubic, 1: close-packed cubic',iunit)
 call write_inopt(handled_layers,'handled_layers','number of handled layers',iunit)
 call write_inopt(wind_radius,'wind_radius','injection radius in units of Rstar',iunit)
 call write_inopt(wind_injection_x,'wind_injection_x','injection x in units of Rstar',iunit)
 call write_inopt(wind_length,'wind_length','wind length in units of Rstar',iunit)

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ieos,polyk,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error,fatal
 use units,         only:select_unit,unit_density,unit_pressure,unit_velocity
 use setunits,      only:read_options_and_set_units
 use eos,           only:gamma
 character(len=*), intent(in) :: filename
 integer,          intent(inout) :: ieos
 real,             intent(inout) :: polyk
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0

 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr)

 call read_inopt(ieos,'ieos',db,errcount=nerr) ! equation of state

 call read_inopt(a,'a',db,errcount=nerr)
 call read_inopt(mdon,'mdon',db,errcount=nerr)
 call read_inopt(macc,'macc',db,errcount=nerr)
 call read_inopt(hacc,'hacc',db,errcount=nerr)

 call read_inopt(nstar_in,'nstar',db,errcount=nerr)

 call read_inopt(v_inf,'v_inf',db,errcount=nerr)
 call read_inopt(rho_inf,'rho_inf',db,errcount=nerr)
 call read_inopt(pres_inf,'pres_inf',db,errcount=nerr)
 call read_inopt(gamma,'gamma',db,errcount=nerr)

 ! Convert wind quantities to code units
 v_inf = v_inf / unit_velocity * 1.e5
 rho_inf = rho_inf / unit_density
 pres_inf = pres_inf / unit_pressure

 call read_inopt(lattice_type,'lattice_type',db,errcount=nerr)
 call read_inopt(handled_layers,'handled_layers',db,errcount=nerr)
 call read_inopt(wind_radius,'wind_radius',db,errcount=nerr)
 call read_inopt(wind_injection_x,'wind_injection_x',db,errcount=nerr)
 call read_inopt(wind_length,'wind_length',db,errcount=nerr)


 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

call close_db(db)
end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Check that setup is sensible
!+
!-----------------------------------------------------------------------
subroutine check_setup(pmass,ierr)
 real, intent(in)     :: pmass
 integer, intent(out) :: ierr
 real                 :: min_layer_sep

 ierr = 0

 min_layer_sep = (pmass / rho_inf)**(1./3.)

 if ( abs(wind_injection_x - 1.)*Rstar < real(handled_layers)*min_layer_sep ) then
    print*,'error: Handled layers overlap with sphere. Try decreasing wind_injection_x or handled_layers'
    ierr = 1
 endif
 if (wind_radius < 1.) then
    print*,'error: Wind cross-section should not be smaller than the sphere'
    ierr = 1
 endif
 if ( wind_injection_x + wind_length < 1. ) then
    print*,'error: Wind not long enough to cover initial sphere position. Try increasing wind_injection_x or wind_length'
    ierr = 1
 endif

end subroutine check_setup

end module setup

