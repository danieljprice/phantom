!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! this module does setup
!
! :References: None
!
! :Owner: Taj Jankovič & Aleksej Jurca
!
! :Runtime parameters: None
!
! :Dependencies: physcon, units
!
 implicit none
 public :: setpart

 private
 real :: L = 10. !-- Length of the disk (x-direction)
 real :: W = 10. !-- Width of the disk (y-direction)
 real :: H = 3. !-- Height of the disk (z-direction)
 integer :: Nparts = 500000 !-- Number of particles
! real :: Rstar = 1 !-- radius of the star in 1Rsun
 real :: Mintercept = 1.2e-7 !-- Intercepted mass of the disk in Msun for 1Rsun stars
 real :: disk_std = 1. !-- Standard deviation for the disk density profile along z-direction
 real :: u_0 = 1.e-5 !-- Internal energy
 character(len=120) :: disk_profile = 'uniform' !-- Vertical density profile of the disk (uniform, gauss)
 integer,private :: iseed = -123456789
contains

!----------------------------------------------------------------
!+
!  empty setup for driven simulation
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use units,       only:set_units,umass,unit_density,udist
 use physcon,     only:au,solarm,pi,solarr,c
 use io,          only:master,fatal
 use part,        only:igas,maxp
 use prompting,   only:prompt
 use random,      only:ran2,gauss_random
 use boundary,    only:set_boundary
 use timestep,    only:tmax,dtmax

 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)

 character(len=120) :: filename,infile

 real :: xmin,xmax,ymin,ymax,zmin,zmax
 real :: rho,u,z_p,mpart,total_m,rho_0,r_unit,v_unit, sigma
 integer :: i,ierr
 logical :: iexist

 infile = trim(fileprefix)//'.in'
 iexist = .false.
 inquire(file=trim(infile),exist=iexist)
 !disk_profile = 'gauss' ! uniform or gauss
 !
 !-- Read runtime parameters from setup file
 !
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Cubic lattice'
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

 ! units
 r_unit = solarr
 v_unit = 0.1*c ! typical velocity in QPEs
 print *,'-----Setting units for r_unit = ', r_unit,', and v_unit = ',v_unit
 print *,'-----Time unit = ', r_unit/v_unit

 !set new boundary
 xmin = -L
 xmax = L
 ymin = -W
 ymax = W
 zmin = -H
 zmax = H

 call set_boundary(xmin,xmax,ymin,ymax,zmin,zmax)
 call set_units(dist=r_unit,time=r_unit/v_unit,G=1.)
 time = 0.
 polyk = 0.
 gamma = 4./3.
 hfact = 1.2

 !get the number of particles
 if (Nparts > maxp) then
   call fatal('set_localdisk','Nparts exceeds maxp; reduce Nparts or increase maxp/recompile')
 endif
 npart = Nparts

 ! sigma = column density (integrated over z) in code units
 sigma = Mintercept * solarm / (pi * solarr**2)/umass*udist**2

 ! total mass in the box
 total_m = sigma * 4 * W * L

 mpart = total_m / npart
 npartoftype(igas) = npart
 massoftype = mpart

 ! parameters
 if (disk_profile=='uniform') then
   rho_0 = total_m / (8*L* W * H)
 else
   !rho_0 = total_m / (4*L* W * sqrt(2*pi)*disk_std)
   !rho_0 = sigma/disk_std ! following Huang et al. 2025 (neglect sqrt(2*pi)
   rho_0 = sigma / (sqrt(2.0*pi) * disk_std * erf( H / (sqrt(2.0) * disk_std) ))
 endif
 u = u_0
 tmax      = 1000.
 dtmax     = 10.
 rho = rho_0

 print *,'-----Number of particles:', npart
 print *,'-----Size of the box (x,y,z):', L, W, H
 if (disk_profile=='gauss') then
   print *,'-----Gaussian vertical density profile with std=', disk_std
 else
   print *,'-----Uniform vertical density profile'
 endif
 print *,'-----Mass of a particle:', mpart, 'Total mass:', total_m
 print *,'-----Central density is (code units):', rho_0, 'Central density is (cgs):', rho_0*unit_density
 print *,'-----Surface density is (code units):', sigma, 'Surface density is (cgs):', sigma*unit_density*udist


 do i=1,npart
    !get random postion
    if (disk_profile=='uniform') then
       rho = rho_0
       xyzh(3,i) = (2.*ran2(iseed) - 1.)*H
    elseif (disk_profile=='gauss') then
       ! generate z_p randomly until inside the disk
       do
         z_p = gauss_random(iseed)
         if (abs(z_p * disk_std) <= H) exit
      end do
      !u = 3 * p_c / rho ! ~constant pressure
      xyzh(3,i) = z_p * disk_std
      rho = rho_0 * exp(-0.5 * xyzh(3,i)**2/disk_std**2)
    endif

    xyzh(1,i) = (2.*ran2(iseed) -1.)*L
    xyzh(2,i) = (2.*ran2(iseed) -1.)*W
    xyzh(4,i) = hfact*(mpart/rho)**(1./3.)

    vxyzu(1,i) = 0.
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    vxyzu(4,i) = u

 enddo
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
! write(iunit,"(a)") '# Half size of the cube (stream height)'
 call write_inopt(L,'L','Length of the disk (x-direction)',iunit)
 call write_inopt(W,'W','Width of the disk (y-direction)',iunit)
 call write_inopt(H,'H','Height of the disk (z-direction)',iunit)
 call write_inopt(Nparts,'Nparts','Number of particles',iunit)
! call write_inopt(Rstar,'Rstar','Radius of the star in 1Rsun',iunit)
 call write_inopt(Mintercept,'Mintercept','Intercepted mass of the disk in Msun for 1Rsun star',iunit)
 call write_inopt(disk_std,'disk_std','Standard deviation for the disk density profile along z-direction',iunit)
 call write_inopt(u_0,'u_0','Internal energy',iunit)
 call write_inopt(disk_profile,'disk_profile','Vertical density profile of the disk (uniform, gauss)',iunit)
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
 call read_inopt(L,'L',db,errcount=nerr)
 call read_inopt(W,'W',db,errcount=nerr)
 call read_inopt(H,'H',db,errcount=nerr)
 call read_inopt(Nparts,'Nparts',db,errcount=nerr)
! call read_inopt(Rstar,'Rstar',db,errcount=nerr)
 call read_inopt(Mintercept,'Mintercept',db,errcount=nerr)
 call read_inopt(disk_std,'disk_std',db,errcount=nerr)
 call read_inopt(u_0,'u_0',db,errcount=nerr)
 call read_inopt(disk_profile,'disk_profile',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile


end module setup
