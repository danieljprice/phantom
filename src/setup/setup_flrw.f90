!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for a constant density + petrubtations FLRW universe 
! as described in Magnall et al. 2023 
!
! :References: None
!
! :Owner: Spencer Magnall
!
! :Runtime parameters:
!   - cs0                 : *initial sound speed in code units*
!   - dist_unit           : *distance unit (e.g. au)*
!   - ilattice            : *lattice type (1=cubic, 2=closepacked)*
!   - mass_unit           : *mass unit (e.g. solarm)*
!   - nx                  : *number of particles in x direction*
!   - radiation_dominated : *Radiation dominated universe (yes/no)*
!   - rhozero             : *initial density in code units*
!
! :Dependencies: boundary, dim, infile_utils, io, mpidomain, mpiutils,
!   options, part, physcon, prompting, setup_params, stretchmap, unifdis,
!   units, utils_gr
!
 use dim,          only:use_dust
 use setup_params, only:rhozero
 use physcon, only:radconst
 implicit none
 public :: setpart

 integer           :: npartx,ilattice
 real              :: cs0,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,ampl,phaseoffset
 character(len=20) :: dist_unit,mass_unit,perturb_direction,perturb,radiation_dominated
 real              :: perturb_wavelength
 real              :: rho_matter
 real(kind=8)      :: udist,umass

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu,gr
 use setup_params, only:npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis,rho_func!,mass_func
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use part,         only:periodic
 use physcon,      only:years,pc,solarm,pi
 use units,        only:set_units
 use mpidomain,    only:i_belong
 use stretchmap,   only:set_density_profile
 use utils_gr, only:perturb_metric, get_u0, get_sqrtg
 !use cons2primsolver, only:primative2conservative

 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(inout) :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=40) :: filename,lattice
 real    :: totmass,deltax
 integer :: i,ierr
 logical :: iexist
 real    :: kwave,denom,length, c1,c3,lambda
 real    :: xval
 real    :: Vup(0:3),phi,sqrtg,gcov(0:3,0:3),alpha,hub
 real    :: last_scattering_temp
 procedure(rho_func), pointer :: density_func

 density_func => rhofunc  ! desired density function
 

 !
 !--general parameters
 !
 perturb_wavelength = 1.0
 time = 0.
 if (maxvxyzu < 4) then
    gamma = 1.
 else
    ! 4/3 for radiation dominated case
    ! irrelevant for
    gamma = 4./3.
 endif
 
 !
 ! default units
 !
 mass_unit = 'solarm'
 dist_unit = 'mpc'
 !
 ! set boundaries to default values
 !
 xmini = xmin; xmaxi = xmax
 ymini = ymin; ymaxi = ymax
 zmini = zmin; zmaxi = zmax
 !
 ! set default values for input parameters
 !
 npartx   = 64
 ilattice = 1
 perturb  = '"no"'
 perturb_direction = '"none"'
 radiation_dominated = '"no"'

 ! Ideally this should read the values of the box length
 ! and initial Hubble parameter from the par file.
 ! Then it should be set using the Friedmann equation:
 !!!!!! rhozero = (3H^2)/(8*pi*a*a)

 hub =  10.553495658357338!/10.
 !hub = 23.588901903912664
 !hub = 0.06472086375185665
 rhozero = 3. * hub**2 / (8. * pi)
 phaseoffset = 0.
 ampl = 0.

 ! Approx Temp of the CMB in Kelvins
 !last_scattering_temp = 3000
 !last_scattering_temp = (rhozero/radconst)**(1./4.)*0.999999999999999d0
 last_scattering_temp = 0.

 ! Define some parameters for Linear pertubations
 ! We assume ainit = 1, but this may not always be the case
 c1 = 1./(4.*pi*rhozero)
 !c2 = We set g(x^i) = 0 as we only want to extract the growing mode
 c3 =  - sqrt(1./(6.*pi*rhozero)) 
 !c3 = hub/(4.d0*PI*rhozero)


 if (gr) then
    ! 0 Because dust?
    cs0 = 0.
 else
    cs0 = 1.
 endif

 ! get disc setup parameters from file or interactive setup
 !
 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    !--read from setup file
    call read_setupfile(filename,ierr)
    if (id==master) call write_setupfile(filename)
    if (ierr /= 0) then
       stop
    endif
 elseif (id==master) then
    call setup_interactive(id,polyk)
    call write_setupfile(filename)
    stop 'rerun phantomsetup after editing .setup file'
 else
    stop
 endif
 !
 ! set units and boundaries
 !
 if (gr) then
    call set_units(dist=udist,c=1.d0,G=1.d0)
 else
    call set_units(dist=udist,mass=umass,G=1.d0)
 endif
 call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
 !
 ! setup particles
 !

 npart = 0
 npart_total = 0
 length = xmaxi - xmini
 deltax = length/npartx
!
! general parameters
!
! time should be read in from the par file
 !time   = 0.08478563386065302 
 time = 0.18951066686763596 ! z~1000
 lambda = perturb_wavelength*length
 kwave  = (2.d0*pi)/lambda
 denom = length - ampl/kwave*(cos(kwave*length)-1.0)
 rhozero = 3.d0 * hub**2 / (8.d0 * pi)
 print*, rhozero

 select case(radiation_dominated)
 case('"yes"')

    ! Set a value of rho_matter 
    rho_matter = 1.e-40
    !rhozero = rhozero - radconst*last_scattering_temp**4
    ! Solve for temperature 
    last_scattering_temp = ((rhozero-rho_matter)/radconst)**(1./4.)
    rhozero = rho_matter
 end select

 xval = density_func(0.75)
 xval = density_func(0.5)

 select case(ilattice)
 case(2)
    lattice = 'closepacked'
 case default
    if (ilattice /= 1) print*,' error: chosen lattice not available, using cubic'
    lattice = 'cubic'
 end select

 select case(perturb)
 case('"yes"')
    select case(perturb_direction)
       !TODO Z AND Y LINEAR PERTURBATIONS
    case('"x"')
       call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,&
                     npart,xyzh,periodic,nptot=npart_total,mask=i_belong,rhofunc=density_func)!,massfunc=mass_function)
    case('"y"')
       call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,&
         npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
       call set_density_profile(npart,xyzh,min=ymin,max=ymax,rhofunc=density_func,&
               geom=1,coord=2)
    case('"all"')
       call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,&
                     npart,xyzh,periodic,nptot=npart_total,mask=i_belong,rhofunc=density_func)
       call set_density_profile(npart,xyzh,min=ymin,max=ymax,rhofunc=density_func,&
               geom=1,coord=2)
       call set_density_profile(npart,xyzh,min=zmin,max=zmax,rhofunc=density_func,&
               geom=1,coord=3)
    end select
 case('"no"')
    call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,&
      npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 end select

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total


 totmass = rhozero*dxbound*dybound*dzbound
 massoftype(1) = totmass/npart_total
 if (id==master) print*,' particle mass = ',massoftype(1)
 if (id==master) print*,' initial sound speed = ',cs0,' pressure = ',cs0**2/gamma



 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs0**2
 else
    polyk = 0.
 endif
 do i=1,npart

    select case(perturb_direction)
    case ('"x"')
       ! should not be zero, for a perturbed wave
       vxyzu(1,i)  = kwave*c3*ampl*cos((2.*pi*xyzh(1,i))/lambda - phaseoffset)
       phi = ampl*sin(kwave*xyzh(1,i)-phaseoffset)
       Vup(1)  = kwave*c3*ampl*cos(2.*pi*xyzh(1,i) - phaseoffset)
       Vup(2:3) = 0.
       call perturb_metric(phi,gcov)
       call get_sqrtg(gcov,sqrtg)

       alpha = sqrt(-gcov(0,0))
       vxyzu(1,i) = Vup(1)*alpha
       vxyzu(2:3,i) = 0.
    case ('"y"')
       vxyzu(2,i)  = kwave*c3*ampl*cos((2.*pi*xyzh(2,i))/lambda - phaseoffset)
       phi = ampl*sin(kwave*xyzh(2,i)-phaseoffset)
       Vup = 0.
       Vup(2)  = kwave*c3*ampl*cos(2.*pi*xyzh(2,i) - phaseoffset)

       call perturb_metric(phi,gcov)
       call get_sqrtg(gcov,sqrtg)

       alpha = sqrt(-gcov(0,0))
       vxyzu(:,i) = 0.
       vxyzu(2,i) = Vup(2)*alpha

    case ('"all"')
       phi = ampl*(sin(kwave*xyzh(1,i)-phaseoffset) - sin(kwave*xyzh(2,i)-phaseoffset) - sin(kwave*xyzh(3,i)-phaseoffset))
       Vup(1)  = kwave*c3*ampl*cos((2.*pi*xyzh(1,i))/lambda - phaseoffset)
       Vup(2)  = kwave*c3*ampl*cos((2.*pi*xyzh(2,i))/lambda - phaseoffset)
       Vup(3)  = kwave*c3*ampl*cos((2.*pi*xyzh(3,i))/lambda - phaseoffset)

       call perturb_metric(phi,gcov)
       call get_sqrtg(gcov,sqrtg)

       alpha = sqrt(-gcov(0,0))

       ! perturb the y and z velocities
       vxyzu(1,i)  = Vup(1)*alpha
       vxyzu(2,i)  = Vup(2)*alpha
       vxyzu(3,i)  = Vup(3)*alpha
    end select
    ! Setup the intial internal energy here?
    ! This should be u = aT^4/\rho
    ! Choose an initial temp of the cmb ~ 3000K
    ! Set a=1 for now
    ! Asssuming that this is constant density/pressure for now so I'm making sure that
    ! Note that rhozero != rho
    ! rhozero = rho + rho*u as this is the energy density
    select case(radiation_dominated)
    case('"yes"')
       if (maxvxyzu >= 4 .and. gamma > 1.) vxyzu(4,i) = (radconst*(last_scattering_temp**4))/rhozero !vxyzu(4,i) = cs0**2/(gamma*(gamma-1.))
       ! Check that the pressure is correct
       print*, "Pressure: ", (gamma-1)*rhozero*vxyzu(4,i)
       print*, "Pressure from energy density: ", 3. * hub**2 / (8. * pi)/3.
       print*, "Pressure 1/3 \rho u: ",radconst*(last_scattering_temp**4)/3.
       print*, "particle mass: ", massoftype
    end select
 enddo


contains
!----------------------------------------------------
!+
!  callback function giving desired density profile
!+
!----------------------------------------------------
real function rhofunc(x)
 use utils_gr, only:perturb_metric, get_u0, get_sqrtg
 real, intent(in) :: x
 real :: const, phi, rhoprim, gcov(0:3,0:3), sqrtg,u0,v(3),Vup(3)
 real :: alpha
 integer :: ierr

 !rhofunc = 1.d0 + ampl*sin(kwave*(x-xmin))
 !rhofunc = ampl*sin(kwave*(x-xmin))
 ! Eq 28. in Macpherson+ 2017
 ! Although it is missing a negative sign
 const = -kwave*kwave*c1 - 2.
 phi = ampl*sin(kwave*x-phaseoffset)
 !rhofunc = rhozero*(1.d0 + const*ampl*sin(kwave*x))
 ! Get the primative density from the linear perb
 rhoprim = rhozero*(1.d0+const*phi)

 ! Get the perturbed 4-metric
 call perturb_metric(phi,gcov)
 ! Get sqrt(-det(g))
 call get_sqrtg(gcov,sqrtg)
 ! Define the 3 velocities to calculate u0
 ! Three velocity will need to be converted from big V to small v
 !
 Vup(1) = kwave*c3*ampl*cos((2.*pi*x)/lambda-phaseoffset)
 Vup(2:3) = 0.
 alpha = sqrt(-gcov(0,0))
 v(1) = Vup(1)*alpha
 v(2:3) = 0.
 ! calculate u0
 ! TODO Should probably handle this error at some point
 call get_u0(gcov,v,u0,ierr)
 ! Perform a prim2cons
 rhofunc = rhoprim*u0*sqrtg

end function rhofunc

real function massfunc(x,xmin)
 use utils_gr, only:perturb_metric, get_u0, get_sqrtg,dot_product_gr
 real, intent(in) :: x,xmin
 real :: const, gcov(0:3,0:3), sqrtg,u0,v(3),Vup(3)
 real :: massprimx,massprimmin,massprim
 real :: lorrentz, bigv2

 ! The value inside the bracket
 const = -kwave*kwave*c1 - 2.
 phi = ampl*sin(kwave*x-phaseoffset)
 !expr = ampl*(-(1./kwave))*cos(phaseoffset - (2.d0*pi*x)/lambda)
 !exprmin = ampl*(-(1./kwave))*cos(phaseoffset - (2.d0*pi*xmin)/lambda)
 massprimx = (x+deltaint(x))
 massprimmin = (xmin+deltaint(xmin))
 ! Evalutation of the integral
 ! rho0[x-Acos(kx)]^x_0
 massprim = rhozero*(massprimx - massprimmin)
 print*, massprim

 ! Get the perturbed 4-metric
 call perturb_metric(phi,gcov)
 ! Get sqrt(-det(g))
 call get_sqrtg(gcov,sqrtg)
 ! Define the 3 velocities to calculate u0
 ! Three velocity will need to be converted from big V to small v
 !
 Vup(1) = kwave*c3*ampl*cos((2.*pi*x)/lambda-phaseoffset)
 Vup(2:3) = 0.
 alpha = sqrt(-gcov(0,0))
 !v(0) = 1
 v(1) = Vup(1)*alpha
 v(2:3) = 0.
 bigv2 = dot_product_gr(Vup,Vup,gcov)
 lorrentz =  1./sqrt(1.-bigv2)
 call get_u0(gcov,v,u0,ierr)
 massfunc = (massprim)!*lorrentz
 massfunc = massprim!*sqrtg*u0

end function massfunc

real function deltaint(x)
   real, intent(in) :: x

   deltaint = (1./kwave)*(kwave*kwave*c1 - 2)*ampl*cos(2*pi*x/lambda)

end function deltaint

end subroutine setpart

!------------------------------------------------------------------------
!
! interactive setup
!
!------------------------------------------------------------------------
subroutine setup_interactive(id,polyk)
 use io,        only:master
 use mpiutils,  only:bcast_mpi
 use dim,       only:maxp,maxvxyzu
 use prompting, only:prompt
 use units,     only:select_unit
 integer, intent(in)  :: id
 real,    intent(out) :: polyk
 integer              :: ierr

 if (id==master) then
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo

    call prompt('enter xmin boundary',xmini)
    call prompt('enter xmax boundary',xmaxi,xmini)
    call prompt('enter ymin boundary',ymini)
    call prompt('enter ymax boundary',ymaxi,ymini)
    call prompt('enter zmin boundary',zmini)
    call prompt('enter zmax boundary',zmaxi,zmini)
 endif
 !
 ! number of particles
 !
 if (id==master) then
    print*,' uniform setup... (max = ',nint((maxp)**(1/3.)),')'
    call prompt('enter number of particles in x direction ',npartx,1)
 endif
 call bcast_mpi(npartx)
 !
 ! mean density
 !
 if (id==master) call prompt(' enter density (gives particle mass)',rhozero,0.)
 call bcast_mpi(rhozero)
 !
 ! sound speed in code units
 !
 if (id==master) then
    call prompt(' enter sound speed in code units (sets polyk)',cs0,0.)
 endif
 call bcast_mpi(cs0)
 !
 ! type of lattice
 !
 if (id==master) then
    call prompt(' select lattice type (1=cubic, 2=closepacked)',ilattice,1)
 endif
 call bcast_mpi(ilattice)
end subroutine setup_interactive

!------------------------------------------------------------------------
!
! write setup file
!
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for uniform setup routine'

 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 !
 ! boundaries
 !
 write(iunit,"(/,a)") '# boundaries'
 call write_inopt(xmini,'CoordBase::xmin','xmin boundary',iunit)
 call write_inopt(xmaxi,'CoordBase::xmax','xmax boundary',iunit)
 call write_inopt(ymini,'CoordBase::ymin','ymin boundary',iunit)
 call write_inopt(ymaxi,'CoordBase::ymax','ymax boundary',iunit)
 call write_inopt(zmini,'CoordBase::zmin','zmin boundary',iunit)
 call write_inopt(zmaxi,'CoordBase::zmax','zmax boundary',iunit)



 !
 ! other parameters
 !
 write(iunit,"(/,a)") '# setup'
 call write_inopt(npartx,'nx','number of particles in x direction',iunit)
 call write_inopt(rhozero,'rhozero','initial density in code units',iunit)
 call write_inopt(cs0,'cs0','initial sound speed in code units',iunit)
 call write_inopt(perturb,'FLRWSolver::FLRW_perturb','Pertrubations of FLRW?',iunit)
 call write_inopt(ampl,'FLRWSolver::phi_amplitude','Pertubation amplitude',iunit)
 call write_inopt(phaseoffset,'FLRWSolver::phi_phase_offset','Pertubation phase offset',iunit)
 call write_inopt(perturb_direction, 'FLRWSolver::FLRW_perturb_direction','Pertubation direction',iunit)
 call write_inopt(radiation_dominated, 'radiation_dominated','Radiation dominated universe (yes/no)',iunit)
 call write_inopt(perturb_wavelength,'FLRWSolver::single_perturb_wavelength','Perturbation wavelength',iunit)
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!
! read setup file
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use units,        only:select_unit
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 ! units
 !
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr)
 call read_inopt(dist_unit,'dist_unit',db,errcount=nerr)
 !
 ! boundaries
 !
 call read_inopt(xmini,'CoordBase::xmin',db,errcount=nerr)
 call read_inopt(xmaxi,'CoordBase::xmax',db,min=xmini,errcount=nerr)
 call read_inopt(ymini,'CoordBase::ymin',db,errcount=nerr)
 call read_inopt(ymaxi,'CoordBase::ymax',db,min=ymini,errcount=nerr)
 call read_inopt(zmini,'CoordBase::zmin',db,errcount=nerr)
 call read_inopt(zmaxi,'CoordBase::zmax',db,min=zmini,errcount=nerr)
 !
 ! other parameters
 !
 call read_inopt(npartx,'nx',db,min=8,errcount=nerr)
 call read_inopt(rhozero,'rhozero',db,min=0.,errcount=nerr)
 call read_inopt(cs0,'cs0',db,min=0.,errcount=nerr)

 call read_inopt(perturb_direction,'FLRWSolver::FLRW_perturb_direction',db,errcount=nerr)
 call read_inopt(ampl, 'FLRWSolver::phi_amplitude',db,errcount=nerr)
 call read_inopt(phaseoffset,'FLRWSolver::phi_phase_offset',db,errcount=nerr)
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
 ! TODO Work out why this doesn't read in correctly
 call read_inopt(perturb,'FLRWSolver::FLRW_perturb',db,errcount=nerr)
 call read_inopt(radiation_dominated,'radiation_dominated',db,errcount=nerr)
 call read_inopt(perturb_wavelength,'FLRWSolver::single_perturb_wavelength',db,errcount=nerr)
 !print*, db
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_unifdis','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_unifdis','length unit not recognised')
    ierr = ierr + 1
 endif


end subroutine read_setupfile

end module setup
