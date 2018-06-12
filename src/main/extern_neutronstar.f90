!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_neutronstar
!
!  DESCRIPTION:
!    This module contains routines relating to the computation
!    of the background gravitational force/potential for a
!    (neutron) star, from the equilibrium density profile
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: datafiles, io, physcon, units
!+
!--------------------------------------------------------------------------
module extern_neutronstar
 implicit none

 real, public, allocatable :: r2tab(:), ftab(:)
 real, public                            :: r2surf, fsurf
 integer, public                         :: ntab

 ! *** Add option to .in file to specify density profile / mass enclosed filename? ***
 character(15), public, parameter        :: rhotabfile = "ns-rdensity.tab"
 integer, public, parameter              :: nrhotab = 10001  ! maximum allowed size of r rho tabulated arrays

 public :: neutronstar_force, load_extern_neutronstar, read_rhotab, write_rhotab, calc_menc
 public ::  read_rhotab_wrapper

 private

contains

!----------------------------------------------
!+
!  compute the force on a given particle
!+
!----------------------------------------------
subroutine neutronstar_force(xi,yi,zi,fxi,fyi,fzi,phi)
 real, intent(in)  :: xi, yi, zi
 real, intent(out) :: fxi, fyi, fzi, phi

 real :: ri2, f

 ! check tables have been initialised
 if (.not.(allocated(r2tab) .and. allocated(ftab))) then
    fxi = 0.
    fyi = 0.
    fzi = 0.
    return
 endif
 ri2 = xi*xi + yi*yi + zi*zi

 if (ri2 >= r2surf) then
    f = fsurf * (r2surf / ri2)**1.5
 else
    f = yinterp(ri2, r2tab, ftab, ntab)
 endif

 fxi = f * xi
 fyi = f * yi
 fzi = f * zi
 phi = f * ri2

end subroutine neutronstar_force


!----------------------------------------------
!+
!  load_extern_neutronstar
!
!  Read tabulated r, rho for NS and set
!  up tables of values to interpolate
!  over in extern_neutronstar
!+
!----------------------------------------------
subroutine load_extern_neutronstar(ierr)
 use units,   only: umass, utime, udist
 use physcon, only: pi, gg
 use io,      only: error
 integer, intent(out) :: ierr

 real :: rtab(nrhotab), rhotab(nrhotab), menctab(nrhotab)
 real    :: polyk, gamma, rhoc, gcode
 integer :: i

 call read_rhotab(rhotabfile,nrhotab,rtab,rhotab,ntab,polyk,gamma,rhoc,ierr)

 if (ierr == 0) then
    ! Calculate r^2, F arrays to store
    allocate(r2tab(ntab), ftab(ntab), stat=ierr)
    if (ierr /= 0) then
       call error('extern_neutronstar','Error allocating in load_extern_neutronstar')
       return
    endif
    call calc_menc(ntab, rtab, rhotab, menctab)
    ! Store r^2 and (- G Menc / r^3) values
    r2tab = rtab(1:ntab)**2
    gcode = gg * umass * utime**2 / udist**3
    ftab(1) = 0.0    ! get NaN if dividing by zero at r = 0
    do i = 2, ntab
       ftab(i) = - gcode * menctab(i) / rtab(i)**3
    enddo
    ftab(1) = 0.0
    r2surf = r2tab(ntab)
    fsurf = ftab(ntab)
 else
    ! Set ierr to 1, so that it gives a uniform non-zero error for any failure
    ! (needed to pass external force tests)
    ierr = 1
 endif

end subroutine load_extern_neutronstar


! Read tabulated r, rho from file
subroutine read_rhotab(filename, rsize, rtab, rhotab, nread, polyk, gamma, rhoc, ierr)
 use io, only:error,id,master
 character(len=*), intent(in)  :: filename
 integer,          intent(in)  :: rsize
 integer,          intent(out) :: nread
 integer,          intent(out) :: ierr
 real,             intent(out) :: rtab(rsize), rhotab(rsize)
 real,             intent(out) :: polyk, gamma, rhoc

 integer           :: i, iunit
 character(len=1)  :: hash   ! for reading in leading '#' character

 iunit = 15
 ierr = 0
 open(iunit,file=filename,access='sequential',action='read',status='old',iostat=ierr)
 if (ierr /= 0) then
    if (id==master) call error('extern_neutronstar','Error opening '//trim(filename))
    return
 endif

 ! First line: # K gamma rhoc
 read(iunit, *, iostat=ierr) hash,polyk, gamma, rhoc
 if (ierr /= 0) then
    call error('extern_neutronstar','Error reading first line of header from '//trim(filename))
    return
 endif
 ! Second line: # nentries  (number of r density entries in file)
 read(iunit,*, iostat=ierr) hash,nread
 if (ierr /= 0) then
    call error('extern_neutronstar','Error reading second line of header from '//trim(filename))
    return
 endif
 if (nread > nrhotab) then
    call error('extern_neutronstar','Error with too many entries in density profile file')
    ierr = 1
    return
 endif
 ! Loop over 'n' lines: r and density separated by space
 do i = 1,nread
    read(iunit,*, iostat=ierr) rtab(i), rhotab(i)
    if (ierr /= 0) then
       call error('extern_neutronstar','Error reading data from '//trim(filename))
       return
    endif
 enddo

 close(iunit)

 if (id==master) print *, 'Read density profile...  ierr = ', ierr

end subroutine read_rhotab

! Write tabulated r, rho to file
subroutine write_rhotab(rtab, rhotab, ntab, polyk, gamma, rhoc, ierr)
 real,    intent(in)    :: rtab(:), rhotab(:)
 integer, intent(in)    :: ntab
 real,    intent(in)    :: polyk, gamma, rhoc
 integer, intent(inout) :: ierr

 integer           :: i, iunit

 iunit = 15
 ierr = 0

 open(iunit, file=rhotabfile, access='sequential', action='write', status='replace')

 ! First line: # K gamma rhoc
 write(iunit, *) '# ', polyk, gamma, rhoc

 ! Second line: # nentries  (number of r density entries in file)
 write(iunit,*) '# ', ntab

 if (ntab > nrhotab) then
    print *, 'Error importing density profile: nrhotab < number of entries to read'
    print *, 'nrhotab = ', nrhotab, ' ntab = ', ntab
    ierr = 1
    return
 endif

 ! Loop over 'n' lines: r and density separated by space
 do i = 1,ntab
    write(iunit,*) rtab(i), rhotab(i)
 enddo
 close(iunit)

end subroutine write_rhotab

! Integrate to find table of m_enc values
subroutine calc_menc(n, r, rho, menc_out, totmass)
 use physcon, only:pi
 integer,           intent(in)  :: n
 real,              intent(in)  :: r(:), rho(:)
 real,    optional, intent(out) :: menc_out(n)
 real,    optional, intent(out) :: totmass
 integer                        :: i
 real                           :: r2(n), r2rho(n),menc(n),totalmass

 r2 = r(1:n)**2
 r2rho = r2(1:n) * rho(1:n)

 if (.false.) then
    ! NB: Ensure that this mass calculation is correct if it is to be used (J. Wurster)
    ! Use trapezoid, Simpson's and Simpson's 3/8 for first entries then Simpson's for remaining
    ! (trapezoidal term has largest order error: avoid using it as part of sum for later terms)
    menc(1) = 0.
    menc(2) = (r(2)-r(1)) * (r2rho(1) + r2(2) * rho(2)) / 2.
    menc(3) = (r(3)-r(1)) * (r2rho(1) + 4.*r2rho(2) + r2rho(3)) / 6.
    menc(4) = (r(4)-r(1)) * (r2rho(1) + 3.*r2rho(2) + 3.*r2rho(3) + r2rho(4)) / 8.
    totalmass = menc(1) + menc(2) + menc(3) + menc(4)
    do i = 5, n
       menc(i) = menc(i-2) + (r2rho(i-2) + 4.*r2rho(i-1) + r2rho(i)) * (r(i) - r(i-2)) / 6.
       totalmass = totalmass + menc(i)
    enddo
    menc(:)   = 4.0 * pi * menc(:)
    totalmass = 4.0 * pi * totalmass
 else
    menc(1)   = 4.0/3.0*pi*r(1)**3 * rho(1)
    do i = 2,n
       menc(i)  = menc(i-1) + 4.0/3.0*pi*(r(i)**3 - r(i-1)**3) * rho(i)
    enddo
    totalmass = menc(n)
 endif
 if (present(menc_out)) menc_out = menc
 if (present(totmass))  totmass  = totalmass

end subroutine calc_menc

! Linear 1D interpolation
real function yinterp(x, xtab, ytab, ntab)
 real,    intent(in) :: x
 real,    intent(in) :: ytab(:),xtab(:)
 integer, intent(in) :: ntab

 integer :: ibelow, iabove
 real    :: slope

 yinterp = 0.
 if (x <= xtab(1)) then
    yinterp = ytab(1)
    return
 else if (x >= xtab(ntab)) then
    yinterp = ytab(ntab)
    return
 endif

 ibelow = indexbelow(x, xtab, ntab)
 iabove = ibelow + 1

 slope = (ytab(iabove) - ytab(ibelow)) / (xtab(iabove) - xtab(ibelow))
 yinterp = ytab(ibelow) + (x - xtab(ibelow)) * slope

end function yinterp

! Find index below value x in monotomic array xtab
pure integer function indexbelow(x, xtab, ntab)
 real,    intent(in) :: x
 real,    intent(in) :: xtab(:)
 integer, intent(in) :: ntab

 integer :: ibelow, imid, iabove

 ibelow = 1
 iabove = ntab
 do while ((iabove - ibelow) > 1)
    imid = (iabove + ibelow) / 2
    if ((x > xtab(imid)) .eqv. (xtab(ntab) > xtab(1))) then
       ibelow = imid
    else
       iabove = imid
    endif
 enddo

 indexbelow = ibelow

end function indexbelow

!----------------------------------------------
!+
! Wrapper to get the density profile (J.Wurster)
!+
!----------------------------------------------
subroutine read_rhotab_wrapper(densityfile,ng,r,den,npts,&
                               polyk,gamma,rhoc,Mstar,iexist,ierr)
 use datafiles, only:find_phantom_datafile
 integer,          intent(in)  :: ng
 integer,          intent(out) :: npts,ierr
 real,             intent(out) :: r(ng),den(ng),polyk,gamma,rhoc,Mstar
 logical,          intent(out) :: iexist
 character(len=*), intent(in)  :: densityfile
 character(len=200)            :: fulldensityfile

 ierr = 0
 !
 ! get filename and determine if it exists
 !
 fulldensityfile = find_phantom_datafile(densityfile,'neutronstar')
 inquire(file=trim(fulldensityfile),exist=iexist)
 if (.not.iexist) return
 !
 ! Read the density file
 !
 call read_rhotab(trim(fulldensityfile),ng,r,den,npts,polyk,gamma,rhoc,ierr)
 !
 ! Get the total mass
 !
 call calc_menc(npts,r,den,totmass=Mstar)

end subroutine read_rhotab_wrapper

end module extern_neutronstar
