!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION: Read particle file from SpEC BH-NS merger simulation and set up SPH particles
!
!  REFERENCES: None
!
!  OWNER: Siva Darbha
!
!  RUNTIME PARAMETERS:
!    mhole         -- mass of black hole (solar mass)
!    particle_file_name -- name of SpEC particle file
!
!  DEPENDENCIES: eos, externalforces, infile_utils, io, metric, part,
!    physcon, rho_profile, timestep, units, vectorutils
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 real    :: mhole
 character(len=120) :: particle_file_name

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,ihacc,ihsoft,igas,set_particle_type,rhoh,gravity
 use units,     only:set_units,umass,udist,unit_velocity,unit_energ,unit_density,unit_ergg
 use physcon,   only:solarm,pi,solarr,gg,c,eV,radconst,kboltz
 use io,        only:master,fatal,warning
 use timestep,  only:tmax,dtmax
 use metric,    only:mass1
 use eos,       only:ieos
 use externalforces,only:accradius1,accradius1_hard
 use allocutils, only:allocate_array
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 !----- xyzh(:,i) is the array (x,y,z,h) for the ith SPH particle, where (x,y,z) are the coordinates and h is the smoothing length.
 !----- In Price et al PASA (2018) and in Liptai & Price (2019), h is given by the symbols h and h_a.
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 !----- vxyzu(:,i) is the array (vx,vy,vz,u) for the ith SPH particle, where (vx,vy,vz) are the coordinate velocities and u is the specific internal energy in the rest frame.
 !----- In Price et al PASA (2018) and in Liptai & Price (2019), u is given the symbol u.
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer, parameter :: ntab=5000
 integer :: ierr,i
 logical :: iexist
 real    :: x,y,z
 real, allocatable  :: m_array(:), r_array(:), theta_array(:), phi_array(:), e_array(:), l_array(:), u_r_array(:), u_theta_array(:), y_e_array(:), s_array(:), T_array(:), rho_array(:)
 real    :: m_ejecta
 real    :: r, theta, phi, e, l, u_r, u_theta, dt_dtau, dr_dtau, dtheta_dtau, dphi_dtau, dr_dt, dtheta_dt, dphi_dt
 real    :: vx, vy, vz
 real    :: m, rho, delta
 real    :: T, y_e

!
!-- general parameters
!
 time  = 0.
 ieos  = 2 !----- flag for EOS type. The value ieos=2 chooses an adiabatic EOS
 gamma = 4./3. !----- gamma in adiabatic gamma-law EOS
 polyk = 0.
! if (.not.gravity) call fatal('setup','recompile with GRAVITY=yes')

!
!-- space available for injected gas particles
!
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.
 nptmass        = 0

!
!-- Default runtime parameters
!
!
 mhole         = 0.   ! (solar masses)
 particle_file_name = 'particle_file_name.dat'

!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' BH-NS merger ejecta in GR'
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

!
!-- Convert to code untis
!
 mhole = mhole*solarm
 call set_units(mass=mhole,c=1.,G=1.) !--Set central mass to M=1 in code units

 accradius1_hard = 2.02*mass1
 accradius1      = accradius1_hard

 call get_num_particles(npart)

 allocate(m_array(npart))
 allocate(r_array(npart))
 allocate(theta_array(npart))
 allocate(phi_array(npart))
 allocate(e_array(npart))
 allocate(l_array(npart))
 allocate(u_r_array(npart))
 allocate(u_theta_array(npart))
 allocate(y_e_array(npart))
 allocate(s_array(npart))
 allocate(T_array(npart))
 allocate(rho_array(npart))

 call read_particle_data(npart, m_array, r_array, theta_array, phi_array, e_array, l_array, u_r_array, u_theta_array, y_e_array, s_array, T_array, rho_array)

 m_ejecta = 0
 do i=1,npart
    !----- Convert to code units, i.e. first convert to cgs then normalize to code units
    !----- Note: the particle file is in units of G = c = M_Sun = 1 (for most of the variables), whereas in set_units above we set the code units to be G = c = M_bh = 1
    !----- Note: the particle file gives temperature in units of MeV, which we have to convert here
    m_array(i) = m_array(i) * solarm / umass
    r_array(i) = r_array(i) * (gg*solarm/(c**2.)) / udist
    e_array(i) = e_array(i) * (c**2.) / (unit_velocity**2.)
    l_array(i) = l_array(i) * (gg*solarm/c) / (udist*unit_velocity)
    u_r_array(i) = u_r_array(i) * c / unit_velocity
    u_theta_array(i) = u_theta_array(i) * (gg*solarm/c) / (udist*unit_velocity)
    ! y_e_array(i) = y_e_array(i) !----- Y_e is dimensionless, so we don't need to do any unit conversion
    ! s_array(i) = s_array(i)* !----- I'm not sure what the units are for the entropy in the particle file, I have to look into this
    T_array(i) = T_array(i) * (1.e6)*eV / unit_energ
    rho_array(i) = rho_array(i) * ((c**6.)/((gg**3.)*(solarm**2.))) / unit_density

    m_ejecta = m_ejecta + m_array(i)
 enddo

 npartoftype(igas) = npart
 !----------- The SPH gas particles must all have the same mass, since massoftype(igas) can only take one value -------------
 massoftype(igas)  = m_array(1) !----- set the mass of the SPH gas particles
 do i=1,npart
    call set_particle_type(i,igas)
    r = r_array(i)
    theta = theta_array(i)
    phi = phi_array(i)
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
    e = e_array(i)
    u_r = u_r_array(i)
    u_theta = u_theta_array(i)
    l = l_array(i)
    dt_dtau = e / (1. - (2.*mass1/r))
    dr_dtau = (1. - (2.*mass1/r)) * u_r
    dtheta_dtau = u_theta / (r**2.)
    dphi_dtau = l / ((r**2.)*(sin(theta)**2.))
    dr_dt = dr_dtau / dt_dtau
    dtheta_dt = dtheta_dtau / dt_dtau
    dphi_dt = dphi_dtau / dt_dtau
    vx = sin(theta)*cos(phi)*dr_dt + r*cos(theta)*cos(phi)*dtheta_dt - r*sin(theta)*sin(phi)*dphi_dt
    vy = sin(theta)*sin(phi)*dr_dt + r*cos(theta)*sin(phi)*dtheta_dt + r*sin(theta)*cos(phi)*dphi_dt
    vz = cos(theta)*dr_dt - r*sin(theta)*dtheta_dt
    !----- xyzh(:,i) is the array (x,y,z,h) for the ith SPH particle, where (x,y,z) are the coordinates and h is the smoothing length.
    !----- In Price et al PASA (2018) and in Liptai & Price (2019), h is given by the symbols h and h_a.
    m = m_array(i)
    rho = rho_array(i)
    delta = (m/rho)**(1./3.) !----- local interparticle spacing
    xyzh(4,i)    = hfact*delta !----- smoothing length
    xyzh(1:3,i)  = (/x,y,z/)
    !----- vxyzu(:,i) is the array (vx,vy,vz,u) for the ith SPH particle, where (vx,vy,vz) are the coordinate velocities and u is the specific internal energy in the rest frame.
    !----- In Price et al PASA (2018) and in Liptai & Price (2019), u is given the symbol u.
    T = T_array(i)
    y_e = y_e_array(i)
    vxyzu(4,i)   = ( radconst * ((T*unit_energ/kboltz)**4.) / (rho*unit_density) ) / unit_ergg !----- the specific internal energy in the rest frame, u = a * T^4 / rho
    vxyzu(1:3,i) = (/vx,vy,vz/)
 enddo

 deallocate(m_array)
 deallocate(r_array)
 deallocate(theta_array)
 deallocate(phi_array)
 deallocate(e_array)
 deallocate(l_array)
 deallocate(u_r_array)
 deallocate(u_theta_array)
 deallocate(y_e_array)
 deallocate(s_array)
 deallocate(T_array)
 deallocate(rho_array)

 if (id==master) then
    print "(/,a)",       ' EJECTA SETUP:'
    print "(a,1f10.3)"  ,' EOS gamma = ', gamma
 endif

 if (id==master) print "(/,a,i10,/)",' Number of particles setup = ',npart

 if (npart == 0)   call fatal('setup','no particles setup')
 if (ierr /= 0)    call fatal('setup','ERROR during setup')

end subroutine setpart

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'
 call write_inopt(mhole,        'mhole',        'mass of black hole (solar mass)',            iunit)
 call write_inopt(particle_file_name,'particle_file_name','name of SpEC particle file',                  iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mhole,        'mhole',        db,min=0.,errcount=nerr)
 call read_inopt(particle_file_name,'particle_file_name',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!
!---Get number of particles in particle data--------------------------------------------------
!
subroutine get_num_particles(npart)
 integer :: iostatus, i
 integer, intent(out) :: npart
 real :: m, r, theta, phi, e, l, u_r, u_theta, y_e, s, T, rho

 print *, 'getting number of particles from '//trim(particle_file_name)

 open(1, file=particle_file_name, status='old', action='read')

 do i=1,13
    read(1,*)
 enddo

 i = 0
 do
    read(1,*,IOSTAT=iostatus) m, r, theta, phi, e, l, u_r, u_theta, y_e, s, T, rho
    ! read failed
    if (iostatus > 0) then
        print *, 'read failed, input data not in correct format'
        exit
    ! reached end of file
    else if (iostatus < 0) then
        exit
    ! read successful
    else
        i = i + 1
    endif
 enddo
 npart = i

 close(1)

end subroutine get_num_particles

!
!---Read particle data--------------------------------------------------
!
subroutine read_particle_data(npart, m_array, r_array, theta_array, phi_array, e_array, l_array, u_r_array, u_theta_array, y_e_array, s_array, T_array, rho_array)
 integer :: iostatus, i
 integer, intent(in) :: npart
 real, intent(out)    :: m_array(:), r_array(:), theta_array(:), phi_array(:), e_array(:), l_array(:), u_r_array(:), u_theta_array(:), y_e_array(:), s_array(:), T_array(:), rho_array(:)

 print *, 'reading particle data from '//trim(particle_file_name)

 open(1, file=particle_file_name, status='old', action='read')

 do i=1,13
    read(1,*)
 enddo

 i = 1
 do
    read(1,*,IOSTAT=iostatus) m_array(i), r_array(i), theta_array(i), phi_array(i), e_array(i), l_array(i), u_r_array(i), u_theta_array(i), y_e_array(i), s_array(i), T_array(i), rho_array(i)
    ! read failed
    if (iostatus > 0) then
        print *, 'read failed, input data not in correct format'
        exit
    ! reached end of file
    else if (iostatus < 0) then
        exit
    ! read successful
    else
        i = i + 1
    endif
 enddo

 close(1)

end subroutine read_particle_data

end module setup