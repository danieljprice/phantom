!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! setup for Quebec collaboration (ask T. Tricco for what this does)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - np : *number of particles in the sphere*
!
! :Dependencies: infile_utils, io, part, physcon, prompting, setup_params,
!   spherical, units
!
 implicit none
 public :: setpart

 real, parameter :: rho_crit = 10.0     ! g/cm^3
 real, parameter :: temp_crit = 1.0e7   ! K
 integer :: np = 100000

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,&
                   vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:maxvxyzu,igas
 use spherical, only:set_sphere,rho_func
 use io,        only:master
 use setup_params,   only:rhozero,npart_total
 use units,          only:set_units,udist,umass,utime
 use physcon,        only:pi,solarr,solarm,Rg,gg
 use infile_utils,   only:get_options
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: a
 real :: rmin,rmax,psep,totvol,total_mass,r
 integer :: i,nx,ierr
 real :: Rg_codeunits
 procedure(rho_func), pointer :: density_func

 call set_units(dist=solarr, mass=solarm, G=1.0d0)

 !
 ! parameters for this setup
 !
 ! cgs units
 a = sqrt((Rg * temp_crit) / (pi * gg * rho_crit))
 total_mass = sqrt(16.0 * pi / rho_crit) * (Rg * temp_crit / gg)**(1.5)

 ! convert to code units
 a = a / udist
 total_mass = total_mass / umass
 !
 ! general parameters
 !
 time = 0.
 hfact = 1.2
 gamma = 7./3.
 rmin = 0.0
 rmax = pi*a
 !
 ! setup
 !
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 polyk = 0.

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                 read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 totvol = 4./3.*pi*rmax**3
 nx = int((2.2/3.0*np)**(1./3.))
 psep = totvol**(1./3.)/real(nx)
 print *, ' rho_crit = ', rho_crit, ' in g/cm^3'
 print *, ' radius = ', rmax, ' in solar radius'
 print *, ' total volume = ',totvol,' particle separation = ',psep
 print *, ' totalvol / partsep**3 = ', totvol/psep**3
 print *, ''

 npart = 0
 npart_total = 0

 density_func => rhofunc
 call set_sphere('closepacked',id,master,rmin,rmax,psep,&
      hfact,npart,xyzh,rhofunc=density_func,nptot=npart_total)
!
!--set particle properties
!
 npart_total = npart
 rhozero = total_mass / totvol
 print *, ' total mass = ',total_mass,' mean density = ',rhozero

 npartoftype(:) = 0
 npartoftype(igas) = npart
 print *, ' npart = ',npart_total

 massoftype(igas) = total_mass / npart_total
 print *, ' particle mass = ',massoftype(igas)
 print *, ''

 ! convert Rg to code units
 Rg_codeunits = Rg * utime * utime / (udist * udist)

 do i=1,npart
    vxyzu(1:3,i) = 0.
    r = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))
    if (maxvxyzu >= 4) then
       vxyzu(4,i) = 1.5 * Rg_codeunits * temp_crit * a * sin(r/a) / r
    endif
 enddo

end subroutine setpart

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for Quebec setup routine'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','number of particles in the sphere',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(np,'np',db,ierr,min=10)
 call close_db(db)

end subroutine read_setupfile

!----------------------------------------------------------------
!+
!  Interactive setup routine
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use part,      only:maxp

 call prompt('Enter the approximate number of particles in the sphere ',np,0,maxp)

end subroutine setup_interactive

!--------------------------------------
!+
!  functional form of density profile
!+
!--------------------------------------
real function rhofunc(r)
 use units,          only:udist,umass
 use physcon,        only:pi,Rg,gg,solarm,solarr
 real, intent(in) :: r
 real :: a

 ! cgs units
 a = sqrt((Rg * temp_crit) / (pi * gg * rho_crit))

 ! convert to code units
 a = a / udist

 ! cgs units
 if (r > 0.) then
    rhofunc = rho_crit * a * sin(r/a) / r
 else
    rhofunc = huge(0.)
 endif

 ! convert to code units
 rhofunc = rhofunc / (umass / udist**3)

end function rhofunc

end module setup
