!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up N converging flows in a box.  The flows are
!   elliptical, and their properties can be individually set.
!
! :References: Wurster & Bonnell (2023)
!
! :Owner: James Wurster
!
! :Runtime parameters:
!   - Bzero              : *Bzero in G*
!   - Ncloud             : *number of clouds*
!   - Temperature_medium : *Temperature of the background*
!   - density_contrast   : *density contrast*
!   - dynamic_bdy        : *use dynamic boundary conditions*
!   - h_acc              : *accretion radius (code units)*
!   - h_soft_sinksink    : *sink-sink softening radius (code units)*
!   - icreate_sinks      : *1: create sinks.  0: do not create sinks*
!   - np                 : *actual number of particles in cloud 1*
!   - plasmaB            : *plasma beta*
!   - r_crit             : *critical radius (code units)*
!   - rho_crit_cgs       : *sink formation density (cgs)*
!
! :Dependencies: boundary, boundary_dyn, cooling, datafiles, dim, eos,
!   infile_utils, io, kernel, mpidomain, options, part, physcon, prompting,
!   ptmass, setunits, setup_params, spherical, timestep, unifdis, units,
!   velfield
!
 use part,         only:mhd
 use dim,          only:maxvxyzu,maxp_hard
 use boundary_dyn, only:dynamic_bdy
 implicit none
 public :: setpart

 private
 !--private module variables
 integer, parameter :: Ncloud_max = 2 ! < 10 (might need to modify the code for N > 2)
 integer            :: Ncloud,np
 integer            :: icreate_sinks_setup
 real               :: xmini(3), xmaxi(3)
 real               :: r_cloud(4,Ncloud_max),v_cloud(3,Ncloud_max)
 real               :: cen_cloud(3,Ncloud_max),ang_cloud(3,Ncloud_max)
 real               :: mass_cloud(Ncloud_max),Temp_cloud(Ncloud_max)
 real               :: cs_cloud(Ncloud_max)
 real               :: rms_mach(Ncloud_max),density_contrast,T_bkg,plasmaB,Bzero,angB(3)
 real               :: r_crit_setup,h_acc_setup,h_soft_sinksink_setup,rho_crit_cgs_setup
 logical            :: input_plasmaB
 character(len= 1), parameter :: labelx(4) = (/'x','y','z','r'/)

contains

!----------------------------------------------------------------
!+
!  setup for a ellipses-in-a-box to model converging flows
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,hours,years,au,pc,kboltz,mass_proton_cgs
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:master,fatal,warning
 use unifdis,      only:set_unifdis,get_xyzmin_xyzmax_exact
 use spherical,    only:set_ellipse
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use boundary_dyn, only:width_bkg,rho_thresh_bdy,rho_bkg_ini,dxyz,vbdyx,vbdyy,vbdyz,in_domain,irho_bkg_ini
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_Bfield,unit_velocity,umass,udist
 use eos,          only:polyk2,gmw
 use part,         only:Bxyz,Bextx,Bexty,Bextz,igas,idust,set_particle_type,periodic
 use timestep,     only:dtmax,tmax
 use options,      only:nfulldump,icooling
 use kernel,       only:hfact_default,radkern
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use mpidomain,    only:i_belong
 use ptmass,       only:icreate_sinks,r_crit,h_acc,h_soft_sinksink,rho_crit_cgs
 use cooling,      only:Tfloor
 use setunits,     only:dist_unit,mass_unit,set_units_interactive
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=20), parameter     :: filevx = 'cube_v1.dat'
 character(len=20), parameter     :: filevy = 'cube_v2.dat'
 character(len=20), parameter     :: filevz = 'cube_v3.dat'
 real(kind=8)       :: h_acc_in
 real               :: dy,y0,ang,ang0
 real               :: dens_medium,ndens_medium_cgs,cs_medium,mass_medium,tmax0,dtmax0
 real               :: vol_cloud(Ncloud_max),dens_cloud(Ncloud_max),ndens_cloud_cgs(Ncloud_max)
 real               :: przero(Ncloud_max),xtmp(3),Trho,psep0
 real               :: r_cloud_21(4,Ncloud_max),r_cloud_mid21(4,Ncloud_max),r_cloud_out21(4,Ncloud_max)
 real               :: totmass,vol_box,psep,psep_box,v2,v2max,rsys_max
 real               :: xj,yj,zj,r2,twoh2,sinA,cosA,sinB,cosB,rell2,xij,yij,zij,r0,epot,ekin
 real               :: rmsmach,v2i,turbfac,vrat,angrat,rad2,dm
 real               :: v_cloud_rot(3,Ncloud_max),v_bkg(3),angtmp(3),xyz_nx(3,2)
 integer            :: i,j,k,idb,np0,np_in,nps,npe,npc,npmax,npbkg,ierr,ndead
 integer            :: Np_cloud(0:Ncloud_max)
 logical            :: iexist,in_cloud,make_sinks
 logical            :: moving_bkg   = .true.   ! For each component, will set the background velocity to that of the clouds, if the clouds are the same
 character(len=120) :: filex,filey,filez
 character(len=100) :: filename
 character(len= 40) :: fmt,lattice
 character(len= 10) :: h_acc_char
!
!--Initialise variables required for setup
!
 gmw    = 2.381 ! mean molecular mass
 gamma  = 5./3.
 Np_in  = 0    ! to prevent compiler warning
 dy     = 0.
 ang    = 0.
 y0     = 0.
 ang0   = 0.
 v2max  = 0.
 Bzero  = 0.
 plasmaB = 0.
 npmax   = maxp_hard
 input_plasmaB = .false.
 make_sinks    = .true.
 dynamic_bdy   = .true.
 dist_unit = 'pc'
 mass_unit = 'solarm'

 filex  = find_phantom_datafile(filevx,'velfield')
 filey  = find_phantom_datafile(filevy,'velfield')
 filez  = find_phantom_datafile(filevz,'velfield')
!
!--Read setup file, else prep prompt user for inputs
!
 filename=trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,a),/,1x,63('-'),/)",'  Converging Flows-in-box setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    np_in = np
    if (ierr /= 0) then
       print*, 'Error in reading all required values from setup file.  Aborting.'
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    !
    ! units
    !
    call set_units_interactive()
    !
    ! prompt user for settings
    !
    Ncloud = 2
    npmax        = int(2.0/3.0*size(xyzh(1,:)))/(2*Ncloud) ! approx max number allowed in sphere given size(xyzh(1,:))
    np           = npmax
    np           = 10000
    v_cloud      =   0.0       ! velocity in km/s
    r_cloud(1,:) =  50.0       ! semi-major-axis (x) in pc
    r_cloud(2,:) =  12.5       ! semi-minor-axis (y) in pc
    r_cloud(3,:) =  12.5       ! semi-minor-axis (z) in pc
    v_cloud(1,:) =  21.75      ! velocity along the x-axis in km/s
    r0           =  50.        ! Distance that the leading edge of the cloud travels before collision [pc]
    cs_cloud     =   0.061522  ! sound speed in code units
    cen_cloud    =   0.        ! centre of the cloud in pc
    ang_cloud    =   0.        ! angle of the cloud with respect to the origin
    plasmaB      =   1.0d2     ! plasma beta in the first cloud
    Bzero        =   2.65d-6   ! initial magnetic field strength in G
    angB         =   0.        ! angle of the magnetic field with respect
    angB(1)      = 270.        ! angle of the magnetic field with respect [default is +Bx]
    rms_mach     =  -1.        ! the rms Mach number of the cloud (<0 for a fraction of Epot)
    Temp_cloud   =  10.        ! place-holder temperature in K
    mass_cloud   =   1.0       ! place-holder mass in M_sun
    density_contrast = 100.    ! density contrast between the first cloud and background
    Trho             = 5000.   ! place-holder value
    ndens_cloud_cgs  = 3.      ! initial number density of the cloud (Could also be 30 for a higher density contrast)

    call prompt('Enter the number of clouds to include',Ncloud,0,Ncloud_max)
    do i = 1,Ncloud
       if (i > 1) then
          !--Reset initial conditons all subsequent clouds have the same default as cloud 1
          r_cloud( :,i) = r_cloud( :,1)
          v_cloud( :,i) = v_cloud( :,1)
          ndens_cloud_cgs(i) = ndens_cloud_cgs(1)
          mass_cloud(i) = mass_cloud(1)
          Temp_cloud(i) = Temp_cloud(1)
          cs_cloud(i)   = cs_cloud(1)
          rms_mach(i)   = rms_mach(1)
       endif
       write(*,*) ''
       write(*,'(a,I3,a)') '*** Enter properties for cloud ',i,' currently aligned along xyz ***'
       call prompt('Enter semi-major axis (x) in units of '//dist_unit,r_cloud(1,i),0.)
       call prompt('Enter semi-minor axis (y) in units of '//dist_unit,r_cloud(2,i),0.)
       call prompt('Enter semi-minor axis (z) in units of '//dist_unit,r_cloud(3,i),0.)
       call prompt('Enter the number density of the cloud in units of cm^-3',ndens_cloud_cgs(i))
       ! Mass based upon size & density from cooling curve [calculated here since this is historically a .setup value]
       mass_cloud(i) = 4.0/3.0*pi*r_cloud(1,i)*r_cloud(2,i)*r_cloud(3,i)*(ndens_cloud_cgs(i)*gmw*mass_proton_cgs)*udist**3/umass
       if (i==1) then
          call KIcoolingcurve(ndens_cloud_cgs(i),Temp_cloud(i),Trho,.true.,ierr)
          if (ierr==1) call fatal('setup','Iterated to negative temperature.  Aborting.')
          if (ierr==2) call fatal('setup','Could not converge to a temperature since T < 0.  Aborting.')
          if (ierr==3) call fatal('setup','Could not converge to a temperature since T >> 1.  Aborting.')
          print*, 'Predicted temperature is based upon given number density and the KIO3 cooling curve'
       endif
       call prompt('Enter temperature of the cloud in units of K',Temp_cloud(i),0.)
       call prompt('Enter velocity in x-direction in units of km/s',v_cloud(1,i))
       call prompt('Enter velocity in y-direction in units of km/s',v_cloud(2,i))
       call prompt('Enter velocity in z-direction in units of km/s',v_cloud(3,i))
       v2       = v_cloud(1,i)**2 + v_cloud(2,i)**2 + v_cloud(3,i)**2
       v2max    = max(v2,v2max)
       call prompt('Enter the Mach number of the cloud turbulence (<0 for a fraction of Epot)',rms_mach(i))
       if (i==1) then
          dy   = r0 + r_cloud(1,1)
          ang  = 90.
          y0   = -0.5*(Ncloud-1)*dy
          ang0 = -0.5*(Ncloud-1)*ang
       endif
       if (Ncloud > 1) then
          ang_cloud(1,i) = ang0 + (i-1)*ang
          call prompt('Enter angle of the cloud wrt x in degrees',ang_cloud(1,i))
          call prompt('Enter angle of the cloud wrt z in degrees',ang_cloud(3,i))
          cen_cloud(2,i) = -dy*sin(ang_cloud(1,i)*pi/180.)
          call prompt('Enter x-centre in units of '//dist_unit,cen_cloud(1,i))
          call prompt('Enter y-centre in units of '//dist_unit,cen_cloud(2,i))
          call prompt('Enter z-centre in units of '//dist_unit,cen_cloud(3,i))
       endif
    enddo

    write(*,*) ''
    write(*,'(a)') '*** Enter global properties based upon cloud 1 ***'
    call prompt('Enter the approximate number of particles in cloud 1',np,0,npmax)
    call KIcoolingcurve(ndens_medium_cgs,T_bkg,Trho,.false.,ierr)
    if (ierr==0) then
       density_contrast = ndens_cloud_cgs(1)/ndens_medium_cgs
       print*, 'Using the KI02 cooling curve, for pressure equilibrium, we suggest', &
                ndens_medium_cgs,'cm^-3, which is a density contrast of ',density_contrast
    else
       if (ierr==1) call warning('setup','Iterated to negative temperature.  Aborting.')
       if (ierr==2) call warning('setup','Could not converge to a temperature since T < 0.  Aborting.')
       if (ierr==3) call warning('setup','Could not converge to a temperature since T >> 1.  Aborting.')
       print*, 'Could not predict a background density in pressure equilibrium'
       density_contrast = ndens_cloud_cgs(1)
    endif
    call prompt('Enter density contrast between cloud 1 and the box ',density_contrast,1.)
    ndens_medium_cgs = ndens_cloud_cgs(1)/density_contrast
    call KIcoolingcurve(ndens_medium_cgs,T_bkg,Trho,.true.,ierr)
    print*, 'Using the KI02 cooling curve, the equilibrium temperature is T_bkg = ',T_bkg,'K'
    if (ierr==1) call fatal('setup','Iterated to negative temperature.  Aborting.')
    if (ierr==2) call fatal('setup','Could not converge to a temperature since T < 0.  Aborting.')
    if (ierr==3) call fatal('setup','Could not converge to a temperature since T >> 1.  Aborting.')
    call prompt('Enter temperature of the medium ',T_bkg,1.)
    if (mhd) then
       call prompt('Use Plasma beta to determine the magnetic field strength? else input B0',input_plasmaB)
       if (input_plasmaB) then
          call prompt('Enter the plasma beta in cloud 1 (to determine magnetic field strength)',plasmaB,0.)
       else
          call prompt('Enter the magnetic field strength in G',Bzero,0.)
       endif
       call prompt('Enter the angle between B and the x-axis in degrees ',angB(1))
       call prompt('Enter the angle between B and the z-axis in degrees ',angB(3))
    endif
    np_in    = np

    ! ask about sink particle details; these will not be saved to the .setup file since they exist in the .in file
    write(*,*) ''
    write(*,'(a)') '*** Enter sink properties ***'
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)
    if (make_sinks) then
       rho_crit_cgs_setup = 1.0d-20
       call prompt('Enter the sink formation density (cgs) ',rho_crit_cgs_setup)
       h_acc_char  = '0.25pc'
       call prompt('Enter the accretion radius of the sink (with units; e.g. au,pc,kpc,0.1pc) ',h_acc_char)
       call select_unit(h_acc_char,h_acc_in,ierr)
       h_acc_setup = h_acc_in
       if (ierr==0 ) h_acc_setup = h_acc_setup/udist
       r_crit_setup        = h_acc_setup
       icreate_sinks_setup = 1
       h_soft_sinksink_setup = 2.0*h_acc_setup
    else
       icreate_sinks_setup = 0
    endif

    write(*,*) ''
    write(*,'(a)') '*** Enter box properties ***'
    call prompt('Do you want to use dynamic boundary conditions? ',dynamic_bdy)
    if (.not. dynamic_bdy) then
       xmini(1) =  -2.0*r_cloud(1,1)
       xmaxi(1) =  36.0*r_cloud(1,1)
       xmini(2) = -10.0*r_cloud(1,1)
       xmaxi(2) = -xmini(2)
       xmini(3) =  -5.0*r_cloud(1,1)
       xmaxi(3) = -xmini(3)
       call prompt('Enter minimum x coordinate [code units] ',xmini(1))
       call prompt('Enter maximum x coordinate [code units] ',xmaxi(1))
       call prompt('Enter minimum y coordinate [code units] ',xmini(2))
       call prompt('Enter maximum y coordinate [code units] ',xmaxi(2))
       call prompt('Enter minimum z coordinate [code units] ',xmini(3))
       call prompt('Enter maximum z coordinate [code units] ',xmaxi(3))
    endif
    if (id==master) call write_setupfile(filename)
    stop 'please edit .setup file and rerun phantomsetup'
 else
    stop ! MPI, stop on other threads, interactive on master
 endif
 !
 !--general parameters
 !
 time        = 0.
 hfact       = hfact_default
 tmax0       = 25.*(1.d6*years)/utime
 dtmax0      = tmax0/50.
 !
 !--general parameters of the cloud
 !
 do i = 1,Ncloud
    v2max = max(v2max,v_cloud(1,i)**2 + v_cloud(2,i)**2 +v_cloud(3,i)**2)
 enddo
 v_cloud    = v_cloud/(unit_velocity*1.d-5)        ! from km/s -> code units
 v2max      = v2max  /(unit_velocity*1.d-5)**2     ! from km/s -> code units
 vol_cloud  = 4./3.*pi*r_cloud(1,:)*r_cloud(2,:)*r_cloud(3,:)
 dens_cloud = mass_cloud / vol_cloud
 rhozero    = dens_cloud(1)
 cs_cloud   = sqrt(kboltz*Temp_cloud/(gmw*mass_proton_cgs*(gamma-1.0)))/unit_velocity
 polyk      = cs_cloud(1)**2
 przero     = cs_cloud**2*dens_cloud
 ang_cloud  = ang_cloud*pi/180.
 angB       = angB*pi/180.
 ndens_cloud_cgs = dens_cloud*unit_density/(gmw*mass_proton_cgs)
 do i = 1,Ncloud
    r_cloud(4,i) = max(r_cloud(1,i),r_cloud(2,i),r_cloud(3,i))**2
 enddo
 !
 !--magnetic field parameters
 !
 Bzero       = Bzero  /unit_Bfield                ! from G    -> code units
 if (mhd .and. (plasmaB > 0. .or. Bzero > 0.)) then
    if (input_plasmaB) then
       Bzero = sqrt(2.0*przero(1)/plasmaB)
    else
       plasmaB = 2.0*przero(1)/Bzero**2
    endif
 else
    Bzero   = 0.
    plasmaB = 0.
 endif
 Bextx  = Bzero
 Bexty  = 0.
 Bextz  = 0.
 !
 !--Determine boundary type
 if (dynamic_bdy) then
    if (.not. periodic) call fatal('setup','this iteration requires periodic=true for dynamic boundaries')
    rho_thresh_bdy =  10.*rhozero/density_contrast
    rho_bkg_ini    =  rhozero/density_contrast
    if (density_contrast < 20.) rho_thresh_bdy =  0.5*rhozero*(1. + 1./density_contrast)
    width_bkg      = 100.*pc/udist      ! default in Wurster & Bonnell (2023)
    width_bkg(1,2) = 200.*pc/udist
    width_bkg      = 0.1*r_cloud(1,1)   ! preferred (unimportant value as per Appendix A of Wurster & Bonnell (2023)
    lattice        = 'cubic'
 else
    lattice        = 'closepacked'
 endif
 !
 !--setup particles in the ellipsoids; use this routine to get N_ellipsoid as close to np as possible
 !
 Np_cloud = 0
 xyz_nx(:,1) =  huge(xyz_nx(1,1))
 xyz_nx(:,2) = -huge(xyz_nx(1,1))
 do i = 1,Ncloud
    if (i==1) then
       np0 = np
    else
       np0 = int(mass_cloud(i)/mass_cloud(1)*Np_cloud(1))
    endif

    psep = (r_cloud(1,i)*r_cloud(2,i)*r_cloud(3,i)/np)**(1./3.)
    call set_ellipse(trim(lattice),id,master,r_cloud(1:3,i),psep,hfact,xyzh,npart,npart_total,np0,i_belong)
    psep0       = psep
    Np_cloud(i) = npart
    print "(a,es10.3)",' Particle separation (x-direction) in ellipsoid = ',psep
    print "(a,I8)",' Number of particles in ellipsoid = ', Np_cloud(i)-Np_cloud(i-1)
    if (i==1 .and. np_in/=npart) np = npart
    nps = Np_cloud(i-1) + 1
    npe = Np_cloud(i)
    npc = npe - nps + 1

    ! Calculate turbulence
    rsys_max = max(r_cloud(1,i),r_cloud(2,i),r_cloud(3,i))
    call set_velfield_from_cubes(xyzh(:,nps:npe),vxyzu(:,nps:npe),npc,filex,filey,filez,1.,rsys_max,.false.,ierr)
    if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')
    rmsmach = 0.0
    if (rms_mach(i) > 0.) then
       !--scale the turbulent velocity to the input RMS mach number
       print*, 'Turbulence being set by user'
       do j = nps,npe
          v2i     = dot_product(vxyzu(1:3,j),vxyzu(1:3,j))
          rmsmach = rmsmach + v2i/cs_cloud(i)**2
       enddo
       rmsmach = sqrt(rmsmach/npc)
       if (rmsmach > 0.) then
          turbfac = rms_mach(i)/rmsmach ! normalise the energy to the desired mach number
       else
          turbfac = 0.
       endif
       rms_mach(i) = 0.
       do j = nps,npe
          v2i         = dot_product(vxyzu(1:3,j)*turbfac,vxyzu(1:3,j)*turbfac)
          rms_mach(i) = rms_mach(i) + v2i/cs_cloud(i)**2
       enddo
       rms_mach(i) = sqrt(rms_mach(i)/npc)
    elseif (rms_mach(i) < 0.) then
       print*, 'Turbulence being set by fraction of Epot'
       !--scale the turbulent velocity to a fraction of the potential energy
       epot = 0.
       ekin = 0.
       dm   = mass_cloud(1)/Np_cloud(1)
       !$omp parallel do default(none) &
       !$omp shared(nps,npe,xyzh,vxyzu,dm) &
       !$omp private(j,k,xj,yj,zj,rad2,v2i) &
       !$omp reduction(+:epot,ekin)
       do j = nps,npe
          xj  = xyzh(1,j)
          yj  = xyzh(2,j)
          zj  = xyzh(3,j)
          v2i = dot_product(vxyzu(1:3,j),vxyzu(1:3,j))
          ekin = ekin + 0.5*dm*v2i
          do k = j+1,npe
             rad2 = (xj-xyzh(1,k))**2 + (yj-xyzh(2,k))**2 + (zj-xyzh(3,k))**2
             if (rad2 > 0.) epot = epot + dm*dm/sqrt(rad2)
          enddo
       enddo
       !$omp end parallel do
       turbfac = sqrt(-rms_mach(i)*epot/ekin)
       !--calculate the actual Mach number of the turbulence
       rms_mach(i) = 0.
       ekin        = 0.
       do j = nps,npe
          v2i         = dot_product(vxyzu(1:3,j)*turbfac,vxyzu(1:3,j)*turbfac)
          rms_mach(i) = rms_mach(i) + v2i/cs_cloud(i)**2
          ekin        = ekin + 0.5*dm*v2i
       enddo
       rms_mach(i) = sqrt(rms_mach(i)/npc)
    else
       turbfac = 0.
    endif

    ! Rotate clouds and velocity
    sinA = sin(ang_cloud(1,i))
    cosA = cos(ang_cloud(1,i))
    sinB = sin(ang_cloud(3,i))
    cosB = cos(ang_cloud(3,i))
    xtmp = v_cloud(1:3,i)
    v_cloud_rot(1,i) = xtmp(1)*cosA - xtmp(2)*sinA*cosB + xtmp(3)*sinA*sinB
    v_cloud_rot(2,i) = xtmp(1)*sinA + xtmp(2)*cosA*cosB - xtmp(3)*sinA*sinB
    v_cloud_rot(3,i) =                xtmp(2)*     sinB + xtmp(3)*     cosB
    do j = nps,npe
       xtmp         = xyzh(1:3,j)
       xyzh(1,j)    = xtmp(1)*cosA - xtmp(2)*sinA*cosB + xtmp(3)*sinA*sinB + cen_cloud(1,i)
       xyzh(2,j)    = xtmp(1)*sinA + xtmp(2)*cosA*cosB - xtmp(3)*sinA*sinB + cen_cloud(2,i)
       xyzh(3,j)    =                xtmp(2)*     sinB + xtmp(3)*     cosB + cen_cloud(3,i)
       if (maxvxyzu >= 4) then
          vxyzu(4,j) = cs_cloud(i)**2
       endif
       call set_particle_type(j,igas)
       do k = 1,3
          xyz_nx(k,1) = min(xyz_nx(k,1),xyzh(k,j))
          xyz_nx(k,2) = max(xyz_nx(k,2),xyzh(k,j))
       enddo
    enddo

    ! Superimpose translational, turbulent velocity & radial infall
    do j = nps,npe
       vxyzu(1:3,j) = vxyzu(1:3,j)*turbfac + v_cloud_rot(1:3,i)
    enddo
 enddo
 massoftype(igas) = mass_cloud(1)/Np_cloud(1)
 !
 !--Set background velocity
 !
 if (moving_bkg) then
    v_bkg  = v_cloud_rot(1:3,1)
    angtmp = abs(ang_cloud(1:3,1))
    do i = 2,Ncloud
       do j = 1,3
          vrat = v_cloud_rot(j,i)/(v_bkg(j)+epsilon(v_bkg(j)))
          if (vrat > 1.0000001 .or. vrat < 0.99999) v_bkg(j) = 0.
          angrat = abs(ang_cloud(j,i))/(angtmp(j)+epsilon(angtmp(j)))
          if (angrat > 1.0000001 .or. angrat < 0.99999) angtmp(j) = 0.
       enddo
    enddo
 else
    v_bkg = 0.
 endif
 vbdyx = v_bkg(1)
 vbdyy = v_bkg(2)
 vbdyz = v_bkg(3)
 !
 !--general parameters of the box (i.e. the background medium)
 !
 psep_box = psep*(density_contrast)**(1./3.)             ! calculate psep in box (this is x-direction)
 if (dynamic_bdy) then
    dxyz = psep_box
    do k = 1,3
       do j = 1,2
          idb            = nint(width_bkg(k,j)/dxyz)
          width_bkg(k,j) = idb*dxyz
       enddo
    enddo
    do k = 1,3
       xmini(k) = xyz_nx(k,1) - width_bkg(k,1)
       xmaxi(k) = xyz_nx(k,2) + width_bkg(k,2)
    enddo
 else
    if (trim(lattice)=='closepacked') then
       dy = psep_box*sqrt(3./4.)
       xmaxi(2) = xmini(2) + int((xmaxi(2) - xmini(2))/dy)*dy
       dy = psep_box*sqrt(6.)/3.
       xmaxi(3) = xmini(3) + int((xmaxi(3) - xmini(3))/dy)*dy
    endif
 endif
 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))
 cs_medium   = sqrt(kboltz*T_bkg/(gmw*mass_proton_cgs*(gamma-1.0)))/unit_velocity
 vol_box     = dxbound*dybound*dzbound
 dens_medium = rhozero/density_contrast
 polyk2      = cs_medium**2
 mass_medium = vol_box
 totmass     = 0.
 do i = 1,Ncloud
    mass_medium = mass_medium - vol_cloud(i)
    totmass     = totmass + mass_cloud(i)
 enddo
 mass_medium = mass_medium*dens_medium
 totmass     = totmass + mass_medium
 !
 !--setup surrounding low density medium (and then remove overlapping particles)
 !
 call get_xyzmin_xyzmax_exact(trim(lattice),xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3),ierr,psep_box)
 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3)) ! since the boundary is slightly adjusted so dx=dy=dz
 if (trim(lattice)=='random') then
    npbkg = nint(mass_medium/massoftype(igas))
    call set_unifdis(trim(lattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box,hfact,npart,xyzh,.true., &
    npnew_in=npbkg,nptot=npart_total)
 else
    call set_unifdis(trim(lattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box,hfact,npart,xyzh,.true.,nptot=npart_total)
 endif
 vxyzu(:,Np_cloud(Ncloud)+1:npart) = 0.

 dxyz          = psep_box
 ndead         = 0
 r_cloud_21    = 1./(    r_cloud)**2
 r_cloud_mid21 = 1./(2.0*r_cloud)**2
 r_cloud_out21 = 1./(4.0*r_cloud)**2
 twoh2         = hfact*(massoftype(igas)/dens_cloud(1))**(1./3.)
!$omp parallel default(none) &
!$omp shared(Np_cloud,npart,xyzh,vxyzu,Ncloud,cs_medium,twoh2,cen_cloud,r_cloud,r_cloud_21,r_cloud_mid21,r_cloud_out21) &
!$omp shared(v_bkg,v_cloud_rot,ang_cloud,v_cloud,turbfac) &
!$omp private(j,i,xj,yj,zj,r2,in_cloud,sinA,cosA,sinB,cosB,rell2,xij,yij,zij,xtmp,rad2) &
!$omp reduction(+:ndead)
!$omp do schedule(runtime)
 do j = Np_cloud(Ncloud)+1,npart
    if (maxvxyzu >= 4) then
       vxyzu(4,j) = cs_medium**2
    endif
    xj = xyzh(1,j)
    yj = xyzh(2,j)
    zj = xyzh(3,j)
    in_cloud = .false.
    do i = 1,Ncloud
       sinA = sin(-ang_cloud(1,i))
       cosA = cos(-ang_cloud(1,i))
       sinB = sin(-ang_cloud(3,i))
       cosB = cos(-ang_cloud(3,i))
       xij   = xj-cen_cloud(1,i)
       yij   = yj-cen_cloud(2,i)
       zij   = zj-cen_cloud(3,i)
       xtmp(1) = xij*cosA - yij*sinA*cosB + zij*sinA*sinB
       xtmp(2) = xij*sinA + yij*cosA*cosB - zij*sinA*sinB
       xtmp(3) =            yij*     sinB + zij*     cosB
       xij   = xtmp(1)
       yij   = xtmp(2)
       zij   = xtmp(3)
       rell2 = xij*xij*r_cloud_21(1,i) + yij*yij*r_cloud_21(2,i) + zij*zij*r_cloud_21(3,i)
       if (rell2 <= 1.) in_cloud = .true.
    enddo
    if (in_cloud) then
       ndead = ndead + 1
       xyzh(4,j) = -abs(xyzh(4,j))
    endif
    vxyzu(1:3,j) = vxyzu(1:3,j) + v_bkg
 enddo
!$omp enddo
!$omp end parallel
 print*, 'There are ',ndead, ' medium particles within the clouds to remove'
 j = Np_cloud(Ncloud)
 k = j
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 do while (k <= npart)
    do while (xyzh(4,k) < 0. .and. k <=npart)
       k     = k + 1
       npartoftype(igas) = npartoftype(igas) - 1
    enddo
    if (k <= npart) then
       xyzh(:,j)  = xyzh(:,k)
       vxyzu(:,j) = vxyzu(:,k)
    endif
    j = j + 1
    k = k + 1
 enddo
 print "(a,i10,a)",' There remain ',npart-Np_cloud(Ncloud),' particles in low-density medium'
 print*, ""
 !
 !--choose the particle to represent the background density
 irho_bkg_ini = 0
 i = 1
 do while (irho_bkg_ini == 0 .and. i <= npart)
    if (xyzh(1,i)+2.*radkern*xyzh(4,i) > xmax) irho_bkg_ini = i
    i = i + 1
 enddo
 print*, 'irho_bkg_ini = ',irho_bkg_ini
 !
 !--set particle properties (again), plus the properties of the medium
 !
 npart = npartoftype(igas)
 do i = 1,npart
    call set_particle_type(i,igas)
 enddo
 !
 !--add the magnetic field
 !
 if (mhd) then
    ihavesetupB = .true.
    Bxyz = 0.
    do j = 1,npart
       Bxyz(1,j) = -Bzero*sin(angB(1))*cos(angB(3))
       Bxyz(2,j) =  Bzero*cos(angB(1))*cos(angB(3))
       Bxyz(3,j) =  Bzero*             sin(angB(3))
    enddo
 else
    ihavesetupB = .false.
 endif
 !
 !--set default runtime parameters if .in file does not exist
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax      = tmax0
    dtmax     = dtmax0
    nfulldump = 1
    icooling  = 5
    Tfloor    = 3.
    icreate_sinks   = icreate_sinks_setup
    r_crit          = r_crit_setup
    h_acc           = h_acc_setup
    h_soft_sinksink = h_soft_sinksink_setup
    rho_crit_cgs    = rho_crit_cgs_setup
 endif
!
!--summary of cloud properties
!
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,totmass*umass,' g'
 do i = 1,Ncloud
    print"(a,I2)", ' Cloud ',i
    print fmt,'   Mass           : ',mass_cloud(i),mass_cloud(i)*umass,' g'
    print fmt,'   Density        : ',dens_cloud(i),dens_cloud(i)*unit_density,' g/cm^3'
    print fmt,'   Number Density : ',ndens_cloud_cgs(i)/unit_density,ndens_cloud_cgs(i),' cm^-3'
    print fmt,'   T              : ',Temp_cloud(i),Temp_cloud(i),' K'
    print fmt,'   cs             : ',cs_cloud(i),cs_cloud(i)*unit_velocity*1.d-5,' km/s'
    if (mhd) &
    print fmt,'   Alfven speed   : ',Bzero/sqrt(dens_cloud(i)),Bzero/sqrt(dens_cloud(i))*unit_velocity*1.d-5,' km/s'
    print fmt,'   v_{bulk,x}     : ',v_cloud_rot(1,i),v_cloud_rot(1,i)*unit_velocity*1.d-5,' km/s'
    print fmt,'   v_{bulk,y}     : ',v_cloud_rot(2,i),v_cloud_rot(2,i)*unit_velocity*1.d-5,' km/s'
    print fmt,'   v_{bulk,z}     : ',v_cloud_rot(3,i),v_cloud_rot(3,i)*unit_velocity*1.d-5,' km/s'
    print fmt,'   Mach#          : ',rms_mach(i),rms_mach(i)
    print fmt,'   N              : ',Np_cloud(i)-Np_cloud(i-1)
 enddo
 print fmt,' Max cloud vel.   : ',sqrt(v2max),sqrt(v2max)*unit_velocity*1.d-5,' km/s'
 print*, 'Background '
 print fmt,'   Mass           : ',mass_medium,mass_medium*umass,' g'
 print fmt,'   Density        : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 print fmt,'   Number Density : ',dens_medium*umass/(gmw*mass_proton_cgs), &
                                  dens_medium*unit_density/(gmw*mass_proton_cgs),' cm^-3'
 print fmt,'   T              : ',T_bkg,T_bkg, ' K'
 print fmt,'   cs             : ',cs_medium,cs_medium*udist/utime*1.d-5,' km/s'
 if (mhd) &
 print fmt,'   Alfven speed   : ',Bzero/sqrt(dens_medium),Bzero/sqrt(dens_medium)*unit_velocity*1.d-5,' km/s'
 print fmt,'   v_x            : ',v_bkg(1),v_bkg(1)*unit_velocity*1.d-5,' km/s'
 print fmt,'   v_y            : ',v_bkg(2),v_bkg(3)*unit_velocity*1.d-5,' km/s'
 print fmt,'   v_z            : ',v_bkg(3),v_bkg(3)*unit_velocity*1.d-5,' km/s'
 print fmt,'   N              : ',npart-Np_cloud(Ncloud)

 if (dynamic_bdy) then
    print fmt,' dynamic_bdy      : ',dynamic_bdy
    print fmt,' dx               : ',dxyz
    print fmt,' rho_thresh       : ',rho_thresh_bdy,rho_thresh_bdy*unit_density,' g/cm^3'
    print fmt,' medium width     : ',width_bkg(1,1),width_bkg(1,1)*udist/pc,' pc'
 endif

 if (mhd) then
    print fmt,' B field          : ',Bzero,Bzero*unit_Bfield*1.d6,' micro-G'
    print fmt,' plasma beta      : ',plasmaB
 endif

 print fmt,' particle mass    : ',massoftype(igas),massoftype(igas)*umass,' g'
 print fmt,' N_total          : ',npart
 print "(1x,50('-'))"

 print*, 'set object rect from ',xmini(1),',',xmini(2),' to ',xmaxi(1),',',xmaxi(2)

end subroutine setpart

!----------------------------------------------------------------
!+
!  Calculate values using the Koyama & Inutsuka (2002) cooling curve
!  get_T = .true.: Given a density, calculate a temperature
!  get_T = .false.: Given a pressure, find a lower density at the same pressure
!+
!----------------------------------------------------------------
subroutine KIcoolingcurve(nrho,T,pres,get_T,ierr)
 integer, intent(out)   :: ierr
 real,    intent(inout) :: nrho,T,pres
 logical, intent(in)    :: get_T
 integer                :: ctr
 real                   :: Gamma,dGamma,Lambda,dLambda,fatT,fatTdT,Trat,Tnew,dnrho
 logical                :: iterate

 !--initialise variables
 Gamma   = 2.0d-26
 dGamma  = 0.0
 iterate = .true.
 ctr     = 0.
 ierr    = 0
 T       = 10000.
 if (get_T) dnrho = 0.

 !--iterate
 do while ( iterate )
    if (.not. get_T) then
       nrho  =  pres/T
       dnrho = -pres/T**2
    endif
    Lambda  = 1.d7*exp(-1.184d5/(T+1.d3)) + 0.014*sqrt(T)*exp(-92./T) ! This is actually Lamda / Gamma
    dLambda = 0.007*exp(-92./T)*(T+184.)*T**(-1.5) + 1.184d12*exp(-1.184d5/(T+1.d3))*(T+1.d3)**(-2)
    fatT    =  Lambda*Gamma*nrho                      -  Gamma
    fatTdT  = dLambda*Gamma*nrho + Lambda*Gamma*dnrho - dGamma
    Tnew    = abs(T - fatT/fatTdT)
    Trat    = abs( 1.0 - T/Tnew )
    T       = Tnew
    ctr     = ctr + 1
    if (Trat < 1d-6) then
       ! success!
       iterate = .false.
    endif
    if (T < 0.) then
       ierr    = 1
       iterate = .false.
    endif
    if (ctr > 2000) then
       ierr    = 2
       iterate = .false.
    endif
    if (T > 1.d10) then
       ierr    = 3
       iterate = .false.
    endif
 enddo

 if (get_T) then
    pres = T*nrho
 else
    nrho = pres/T
 endif

end subroutine KIcoolingcurve

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 use setunits,     only: write_options_units
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20
 integer                      :: i,j
 character(len=128)           :: label

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for converging flows-in-box setup routines'
 call write_options_units(iunit)

 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','actual number of particles in cloud 1',iunit)

 write(iunit,"(/,a)") '# options for box'
 call write_inopt(dynamic_bdy,'dynamic_bdy','use dynamic boundary conditions',iunit)
 if (.not. dynamic_bdy) then
    do j=1,3
       call write_inopt(xmini(j),labelx(j)//'min',labelx(j)//' min',iunit)
       call write_inopt(xmaxi(j),labelx(j)//'max',labelx(j)//' max',iunit)
    enddo
 endif

 write(iunit,"(/,a)") '# Number of clouds'
 call write_inopt(Ncloud,'Ncloud','number of clouds',iunit)
 do i = 1,Ncloud
    write(iunit,"(/,a,I2)") '# options for cloud', i
    do j = 1,3
       write(label,'(a,I1,a)') 'r_cloud',i,labelx(j)
       call write_inopt(r_cloud(j,i),trim(label),labelx(j)//'-semi-axis of the cloud in code units',iunit)
    enddo
    do j = 1,3
       write(label,'(a,I1,a)') 'v_cloud',i,labelx(j)
       call write_inopt(v_cloud(j,i),trim(label),labelx(j)//'-velocity of the cloud in km/s',iunit)
    enddo
    do j = 1,3
       write(label,'(a,I1,a)') 'cen_cloud',i,labelx(j)
       call write_inopt(cen_cloud(j,i),trim(label),labelx(j)//'-centre of the cloud in code units',iunit)
    enddo
    do j = 1,3
       write(label,'(a,I1,a)') 'ang_cloud',i,labelx(j)
       call write_inopt(ang_cloud(j,i),trim(label),labelx(j)//'-angle of the cloud in degrees',iunit)
    enddo
    write(label,'(a,I1)') 'mass_cloud',i
    call write_inopt(mass_cloud(i),trim(label),'mass of the cloud in code units',iunit)
    write(label,'(a,I1)') 'Temp_cloud',i
    call write_inopt(Temp_cloud(i),trim(label),'Temperature of the cloud in K',iunit)
    write(label,'(a,I1)') 'rms_mach',i
    call write_inopt(rms_mach(i),trim(label),'RMS Mach number of cloud turbulence (<0 for a fraction of Epot)',iunit)
 enddo

 write(iunit,"(/,a)") '# additional properties'
 call write_inopt(density_contrast,'density_contrast','density contrast',iunit)
 call write_inopt(T_bkg,'Temperature_medium','Temperature of the background',iunit)
 if (mhd) then
    if (input_plasmaB) then
       call write_inopt(plasmaB,'plasmaB','plasma beta',iunit)
    else
       call write_inopt(Bzero,'Bzero','Bzero in G',iunit)
    endif
    do j = 1,3
       write(label,'(a,a)') 'angB',labelx(j)
       call write_inopt(angB(j),trim(label),labelx(j)//'-angle of the magnetic field',iunit)
    enddo
 endif

 write(iunit,"(/,a)") '# Sink properties (values in .in file, if present, will take precedence)'
 call write_inopt(icreate_sinks_setup,'icreate_sinks','1: create sinks.  0: do not create sinks',iunit)
 if (icreate_sinks_setup==1) then
    call write_inopt(rho_crit_cgs_setup,'rho_crit_cgs','sink formation density (cgs)',iunit)
    call write_inopt(h_acc_setup,'h_acc','accretion radius (code units)',iunit)
    call write_inopt(r_crit_setup,'r_crit','critical radius (code units)',iunit)
    call write_inopt(h_soft_sinksink_setup,'h_soft_sinksink','sink-sink softening radius (code units)',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use setunits,     only: read_options_and_set_units
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: i,j,nerr,jerr
 character(len=128)            :: label
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr)

 call read_inopt(np,'np',db,ierr)
 call read_inopt(dynamic_bdy,'dynamic_bdy',db,ierr)
 if (.not. dynamic_bdy) then
    do j=1,3
       call read_inopt(xmini(j),labelx(j)//'min',db,ierr)
       call read_inopt(xmaxi(j),labelx(j)//'max',db,ierr)
    enddo
 endif
 call read_inopt(Ncloud,'Ncloud',db,ierr)
 do i = 1,Ncloud
    do j = 1,3
       write(label,'(a,I1,a)') 'r_cloud',i,labelx(j)
       call read_inopt(r_cloud(j,i),label,db,ierr)
    enddo
    do j = 1,3
       write(label,'(a,I1,a)') 'v_cloud',i,labelx(j)
       call read_inopt(v_cloud(j,i),label,db,ierr)
    enddo
    do j = 1,3
       write(label,'(a,I1,a)') 'cen_cloud',i,labelx(j)
       call read_inopt(cen_cloud(j,i),label,db,ierr)
    enddo
    do j = 1,3
       write(label,'(a,I1,a)') 'ang_cloud',i,labelx(j)
       call read_inopt(ang_cloud(j,i),label,db,ierr)
    enddo
    write(label,'(a,I1)') 'mass_cloud',i
    call read_inopt(mass_cloud(i),label,db,ierr)
    write(label,'(a,I1)') 'Temp_cloud',i
    call read_inopt(Temp_cloud(i),label,db,ierr)
    write(label,'(a,I1)') 'rms_mach',i
    call read_inopt(rms_mach(i),label,db,ierr)
 enddo
 call read_inopt(density_contrast,'density_contrast',db,ierr)
 call read_inopt(T_bkg,'Temperature_medium',db,ierr)
 if (mhd) then
    call read_inopt(plasmaB,'plasmaB',db,jerr)
    if (jerr == 0) then
       input_plasmaB = .true.
    else
       call read_inopt(Bzero,'Bzero',db,jerr)
       if (jerr == 0) then
          input_plasmaB = .false.
       else
          ierr = jerr
       endif
    endif
    do j = 1,3
       write(label,'(a,a)') 'angB',labelx(j)
       call read_inopt(angB(j),label,db,ierr)
    enddo
 endif
 call read_inopt(icreate_sinks_setup,'icreate_sinks',db,ierr)
 if (icreate_sinks_setup==1) then
    call read_inopt(rho_crit_cgs_setup,'rho_crit_cgs',db,ierr)
    call read_inopt(h_acc_setup,'h_acc',db,ierr)
    call read_inopt(r_crit_setup,'r_crit',db,ierr)
    call read_inopt(h_soft_sinksink_setup,'h_soft_sinksink',db,ierr)
 endif
 call close_db(db)

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_convergingflows: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup
