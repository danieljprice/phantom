!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for BH-NS merger simulations
!
!  REFERENCES: None
!
!  OWNER: Siva Darbha
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'bhnsejecta'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use metric,    only:mass1
 use units,     only:umass,udist,utime,unit_velocity,unit_energ,unit_density,unit_ergg
 use physcon,   only:solarm,pi,solarr,gg,c,eV,radconst,kboltz
 use part,      only:hfact,pxyzu,dens,metrics,metricderivs
 use metric_tools, only:init_metric
 use utils_gr,  only:h2dens
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: numfile,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: pmass,time
 integer :: iostatus, i
 real    :: x,y,z,h
 real    :: r, theta, phi, e, l, u_r, u_theta, dt_dtau, dr_dtau, dtheta_dtau, dphi_dtau, dr_dt, dtheta_dt, dphi_dt
 real    :: vx, vy, vz, u
 real    :: m, densi
 real    :: T, s, y_e
 real    :: q
 real    :: time_output
 character(len=120) :: particle_file_name, info_file_name



 particle_file_name = trim(dumpfile)//'.dat'

 open(1, file=particle_file_name, action='write', status='replace')

 write(1,'(A)',IOSTAT=iostatus) '# Particle data'
 write(1,'(A)',IOSTAT=iostatus) '# [1] = Mass'
 write(1,'(A)',IOSTAT=iostatus) '# [2] = Radius'
 write(1,'(A)',IOSTAT=iostatus) '# [3] = Theta'
 write(1,'(A)',IOSTAT=iostatus) '# [4] = Phi'
 write(1,'(A)',IOSTAT=iostatus) '# [5] = e (-u_t)'
 write(1,'(A)',IOSTAT=iostatus) '# [6] = l (u_phi)'
 write(1,'(A)',IOSTAT=iostatus) '# [7] = u_r'
 write(1,'(A)',IOSTAT=iostatus) '# [8] = u_theta'
 write(1,'(A)',IOSTAT=iostatus) '# [9] = Y_e'
 write(1,'(A)',IOSTAT=iostatus) '# [10]= s'
 write(1,'(A)',IOSTAT=iostatus) '# [11]= T'
 write(1,'(A)',IOSTAT=iostatus) '# [12]= rho'

 call init_metric(npart,xyzh,metrics,metricderivs)

 do i=1,npart

    m = pmass

    !----- xyzh(:,i) is the array (x,y,z,h) for the ith SPH particle, where (x,y,z) are the coordinates and h is the smoothing length.
    !----- In Price et al PASA (2018) and in Liptai & Price (2019), h is given by the symbols h and h_a.
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    h = xyzh(4,i) !----- the smoothing length

    !----- vxyzu(:,i) is the array (vx,vy,vz,u) for the ith SPH particle, where (vx,vy,vz) are the coordinate velocities and u is the specific internal energy in the rest frame.
    !----- In Price et al PASA (2018) and in Liptai & Price (2019), u is given the symbol u.
    vx = vxyzu(1,i) !----- vx here is defined as vx = dx/dt
    vy = vxyzu(2,i) !----- vy here is defined as vy = dy/dt
    vz = vxyzu(3,i) !----- vz here is defined as vz = dz/dt
    u = vxyzu(4,i) !----- the specific internal energy in the rest frame

    r = (x**2. + y**2. + z**2.)**(1./2.)
    q = (x**2. + y**2.)**(1./2.)
    theta = atan2(q,z)
    phi = atan2(y,x)
    if (phi < 0) then
       phi = phi + 2.*pi
    endif

    dr_dt = (x/r)*vx + (y/r)*vy + (z/r)*vz
    dtheta_dt = (z/(r**2.))*(x/q)*vx + (z/(r**2.))*(y/q)*vy - (q/(r**2.))*vz
    dphi_dt = (-y/(q**2.))*vx + (x/(q**2.))*vy

    dt_dtau = ( -1. * ( -1.*(1. - (2.*mass1/r)) + ((1. - (2.*mass1/r))**(-1.)) * (dr_dt**2.) + (r**2.) * (dtheta_dt**2.) + ((r**2.)*(sin(theta)**2.)) * (dphi_dt**2.) ) )**(-1./2.)

    dr_dtau = dr_dt * dt_dtau
    dtheta_dtau = dtheta_dt * dt_dtau
    dphi_dtau = dphi_dt * dt_dtau

    e = dt_dtau * (1. - (2.*mass1/r))
    u_r = dr_dtau / (1. - (2.*mass1/r))
    u_theta = dtheta_dtau * (r**2.)
    l = dphi_dtau * ((r**2.)*(sin(theta)**2.))

    !----- dens(:) is the array of rest mass densities of the SPH particles, where dens(i) is the rest mass density of the ith SPH particle
    !----- The dens(:) array resides in the part module in the file part.F90, and is updated when the code calls the conservative to primitive variable transformation function
    !----- The routines in GR Phantom use the following naming convention:
    !----- The variable name 'dens' is the rest mass density. In Liptai & Price (2019), the rest mass density is given the symbol rho
    !----- The variable name 'rho' is the "relativistic rest-mass density"/"conserved density". In Liptai & Price (2019), the "relativistic rest-mass density" is given the symbol rho^*
    !call h2dens(densi, xyzh(:,i), metrics(:,:,:,i), vxyzu(1:3,i))
    densi = dens(i)

    T = ( (u*unit_ergg) * (densi*unit_density) / radconst )**(1./4.) !----- temperature in K

    y_e = 0. !----- The SPH gas particles do not store the electron fraction Y_e since the code doesn't evolve the composition, so I just set it to zero here

    s = 0. !----- I'm not sure what the units are for the entropy in the particle file, so I just set it to zero here



    !----- Convert from code units to units of particle file, i.e. first convert to cgs then normalize to units of particle file
    !----- Note: the particle file is in units of G = c = M_Sun = 1 (for most of the variables), whereas in set_units in setup_grbhnsejecta.f90 we set the code units to be G = c = M_bh = 1
    !----- Note: the particle file gives temperature in units of MeV, which we have to convert here
    m = m * umass / solarm
    r = r * udist / (gg*solarm/(c**2.))
    e = e * (unit_velocity**2.) / (c**2.)
    l = l * (udist*unit_velocity) / (gg*solarm/c)
    u_r = u_r * unit_velocity / c
    u_theta = u_theta * (udist*unit_velocity) / (gg*solarm/c)
    T = kboltz * T / ( (1.e6)*eV ) !----- temperature in MeV
    densi = densi * unit_density / ( solarm / (gg*solarm/(c**2.))**3. )



    write(1,'(ES11.4,1X)',advance='no',IOSTAT=iostatus) m
    write(1,'(ES12.5,1X)',advance='no',IOSTAT=iostatus) r
    write(1,'(F9.5,1X)',advance='no',IOSTAT=iostatus) theta
    write(1,'(F9.5,1X)',advance='no',IOSTAT=iostatus) phi
    write(1,'(F10.5,1X)',advance='no',IOSTAT=iostatus) e
    write(1,'(ES12.5,1X)',advance='no',IOSTAT=iostatus) l
    write(1,'(ES12.5,1X)',advance='no',IOSTAT=iostatus) u_r
    write(1,'(ES12.5,1X)',advance='no',IOSTAT=iostatus) u_theta
    write(1,'(F10.7,1X)',advance='no',IOSTAT=iostatus) y_e
    write(1,'(F10.7,1X)',advance='no',IOSTAT=iostatus) s
    write(1,'(ES12.4,1X)',advance='no',IOSTAT=iostatus) T
    write(1,'(ES12.5)',IOSTAT=iostatus) densi

 enddo

 close(1)



 !----- Convert time from code units to units of particle file, i.e. first convert to cgs then normalize to units of particle file
 !----- Note: the particle file is in units of G = c = M_Sun = 1 (for most of the variables), whereas in set_units in setup_grbhnsejecta.f90 we set the code units to be G = c = M_bh = 1
 time_output = time * utime / (gg*solarm/(c**3.))

 info_file_name = trim(dumpfile)//'.info'

 open(2, file=info_file_name, action='write', status='replace')

 write(2,'(A,1X,ES10.3,1X,A)',IOSTAT=iostatus) 'time =', time_output, '[G = c = M_Sun = 1]'

 close(2)



end subroutine do_analysis

end module