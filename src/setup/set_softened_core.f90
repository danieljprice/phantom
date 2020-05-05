!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setsoftenedcore
!
!  DESCRIPTION:
!   This module softens the core of a MESA stellar profile with a cubic
!   density profile, given a softening length and core mass, in preparation
!   for adding a sink particle core. CAUTION: This module does not output
!   self-consistent internal energy and temperature profiles, and just
!   returns the input MESA data.
!
!  REFERENCES: None
!
!  OWNER: Mike Lau
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon
!+
!--------------------------------------------------------------------------
module setsoftenedcore
 use physcon,   only:pi,gg
 implicit none
 real(kind=8)  :: hsoft,msoft,mcore
 integer       :: hidx

 ! hsoft: Softening length of core particle
 ! mcore: Mass of core particle
 ! msoft: Softened mass (mass at softening length minus mass of core particle) 

  contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that reads an input MESA file, calls the subroutines
!  to calculate the cubic core profile, and writes an output Phantom-
!  readable file.
!+
!-----------------------------------------------------------------------
subroutine set_softened_core(filepath,outputpath,mcore,hsoft)
 use physcon,    only:solarm,solarr
 implicit none
 real(kind=8), allocatable, dimension(:) :: m(:),rho(:),pres(:),r(:),phi(:),ene(:),temp(:)
 character(len=120), intent(in)          :: filepath, outputpath
 real(kind=8), intent(in)                :: mcore,hsoft
 real(kind=8)                            :: mc,h
 integer                                 :: i
 logical                                 :: sortedDecreasing, excludeCoreMass

 call read_mesa(rho, r, pres, m, ene, temp, filepath)
 h = hsoft * solarr ! Convert to cm
 call interpolator(r, h, hidx) ! Find index in r closest to h
 mc = mcore * solarm ! Convert to g
 msoft = m(hidx) - mc
 if (msoft < 0) then
     stop 'ERROR: mcore cannot exceed m(r = h)'
 endif
 !
 ! Output data to be sorted from stellar surface to interior?
 sortedDecreasing = .true. ! Needs to be true if to be read by Phantom
 !
 ! Exclude core mass in output mass coordinate?
 excludeCoreMass = .true. ! Needs to be true if to be read by Phantom
 !
 ! Calculate density profile inside softening length and the mass profile inside softening length
 ! (Note: This mass includes the contribution by the core particle)
 !
 call calc_rho_and_m(rho, m, r, mc, h)
 !
 ! Calculate gravitational potential
 call calc_phi(r, m-mc, phi, mc, h)
 !
 ! Calculate pressure profile inside softening length
 call calc_pres(r, rho, phi, pres)
 !
 ! Write data
 ! Note: The temperature and internal energy are fake, since they are from the original mesa file
 !
 ! Reverse arrays so that data is sorted from stellar surface to stellar centre.
 !
 if (sortedDecreasing) then
     call flip_array(m)
     call flip_array(pres)
     call flip_array(temp)
     call flip_array(r)
     call flip_array(rho)
     call flip_array(ene)
     call flip_array(phi)
 end if

 if (excludeCoreMass) then
     m = m - mc
 endif

 open(1, file = outputpath, status = 'new')  
 write(1,'(a)') '[    Mass   ]  [  Pressure ]  [Temperature]  [   Radius  ]  [  Density  ]  [   E_int   ]'
 write(1,42) (m(i), pres(i), temp(i), r(i), rho(i), ene(i), i = 1, size(r))
 42 format (es13.7, 2x, es13.7, 2x, es13.7, 2x, es13.7, 2x, es13.7, 2x, es13.7)
 close(1,status='keep')
end subroutine set_softened_core


subroutine calc_rho_and_m(rho,m,r,mc,h)
 implicit none
 real(kind=8) :: a, b, d, drhodr_h
 real(kind=8), intent(in) :: mc, h
 real(kind=8), dimension(1:hidx+1), intent(in) :: r
 real(kind=8), dimension(1:hidx+1), intent(inout) :: rho, m
            
 ! a, b, d: Coefficients of cubic density profile defined by rho(r) = ar**3 + br**2 + d
 drhodr_h = (rho(hidx+1) - rho(hidx)) / (r(hidx+1) - r(hidx)) ! drho/dr at r = h
 a = 2d0/h**2d0 * drhodr_h - 1d1/h**3d0 * rho(hidx) + 7.5d0/pi/h**6d0 * msoft
 b = 0.5d0*drhodr_h/h - 1.5d0*a*h
 d = rho(hidx) - 0.5d0*h*drhodr_h + 0.5d0*a*h**3d0

 rho(1:hidx) = a*r(1:hidx)**3d0 + b*r(1:hidx)**2d0 + d

 ! Mass is then given by m(r) = mcore + 4*pi (1/6 a r^6 + 1/5 b r^5 + 1/3 d r^3)
 m(1:hidx) = mc + 4d0*pi * (1d0/6d0 * a * r(1:hidx)**6d0 + 0.2d0 * b * r(1:hidx)**5d0 + &
                            1d0/3d0 * d * r(1:hidx)**3d0)
end subroutine calc_rho_and_m


subroutine calc_phi(r,mgas,phi,mc,h)
 implicit none
 real(kind=8), dimension(:), intent(in) :: r(:), mgas(:)
 real(kind=8), dimension(:), allocatable :: q
 real(kind=8), intent(in) :: mc, h
 real(kind=8), dimension(:), allocatable, intent(out) :: phi
 real(kind=8), dimension(size(r)) :: phi_core, phi_gas
 integer :: idx2h, i
 ! The gravitational potential is needed to integrate the pressure profile using the
 ! equation of hydrostatic equilibrium. First calculate gravitational potential due
 ! to point mass core, according to the cubic spline kernel in Price & Monaghan (2006)
 ! to be consistent with Phantom, then calculate the gravitational potential due to the
 ! softened gas.
 allocate(phi(size(r)))

 ! (i) Gravitational potential due to core particle (cubic spline softening)
 ! For 0 <= r/h < 1
 allocate(q(1:hidx-1))
 q = r(1:hidx-1) / h
 phi_core(1:hidx-1) = gg*mc/h * (2d0/3d0*q**2d0 - 0.3d0*q**4d0 + 0.1d0*q**5d0 &
                                   - 7d0/5d0)
 deallocate(q)

 ! For 1 <= r/h < 2
 call interpolator(r, 2*h, idx2h) ! Find index corresponding to r = 2h
 allocate(q(hidx:idx2h-1))
 q = r(hidx:idx2h-1) / h
 phi_core(hidx:idx2h-1) = gg*mc/h * (4d0/3d0*q**2d0 - q**3d0 + 0.3d0*q**4d0 &
                                             - 1d0/30d0*q**5d0 - 1.6d0 + 1/15d0/q)
 deallocate(q)

 ! For 2 <= r/h
 phi_core(idx2h:size(r)) = - gg * mc / r(idx2h:size(r))

 ! (ii) Gravitational potential due to softened gas
 phi_gas(size(r)) = - gg * mgas(size(r)) / r(size(r)) ! Surface boundary condition for phi
 do i = 1, size(r) - 1
  phi_gas(size(r)-i) = phi_gas(size(r)-i+1) - gg * mgas(size(r)-i) / r(size(r)-i)**2d0 &
                                              * (r(size(r)-i+1) - r(size(r)-i))
 end do
      
 ! (iii) Add the potentials 
 phi = phi_gas + phi_core 
end subroutine calc_phi


subroutine calc_pres(r, rho, phi,pres)
 ! Calculates pressure by integrating the equation of hydrostatic equilibrium
 ! given the gravitational potential and the density profile
 implicit none
 real(kind=8), dimension(:), intent(in) :: rho(:), phi(:), r(:)
 real(kind=8), dimension(size(rho)), intent(out) :: pres(1:size(rho))
 integer :: i

 pres(size(r)) = 0 ! Set boundary condition of zero pressure at stellar surface
 do i = 1, size(r)-1
  ! Reverse Euler
  pres(size(r) - i) = pres(size(r)-i+1) + rho(size(r)-i+1) * (phi(size(r)-i+1) - phi(size(r)-i))
 end do
end subroutine calc_pres
    

subroutine interpolator(array, value, valueidx)
 implicit none
 real(kind=8), dimension(:), intent(in) :: array(:) 
 real(kind=8), intent(in) :: value
 integer, intent(out) :: valueidx
 ! A subroutine to interpolate an array given a value, returning the index closest to the
 ! required value. Only works if array is ordered.
 valueidx = minloc(abs(array - value), dim = 1)
end subroutine interpolator
  
  
subroutine flip_array(array)
 implicit none
 real(kind=8), dimension(:), intent(inout) :: array(:)
 real(kind=8), dimension(:), allocatable :: flipped_array(:)
 integer :: i
 ! A subroutine that reverses the elements of a 1-d array
 allocate(flipped_array(size(array)))
 do i = 1, size(array)
  flipped_array(i) = array(size(array) - i + 1)
 end do
 array = flipped_array
end subroutine flip_array


subroutine read_mesa(rho,r,pres,m,ene,temp,filepath)
 implicit none
 integer                                           :: lines,rows=0,i
 character(len=120), intent(in)                    :: filepath
 character(len=10000)                              :: dumc
 character(len=24),allocatable                     :: header(:),dum(:)
 real(kind=8),allocatable,dimension(:,:)           :: dat
 real(kind=8),allocatable,dimension(:),intent(out) :: rho,r,pres,m,ene,temp
    
 ! reading data from datafile ! -----------------------------------------------
 open(unit=40,file=filepath,status='old')
 read(40,'()')
 read(40,'()')
 read(40,*) lines, lines
 read(40,'()')
 read(40,'()')
 read(40,'(a)') dumc! counting rows
 allocate(dum(500)) ; dum = 'aaa'
 read(dumc,*,end=101) dum
 101 do i = 1, 500
  if (dum(i)=='aaa') then
   rows = i-1
  exit
  end if
 end do

 allocate(header(1:rows),dat(1:lines,1:rows))
 header(1:rows) = dum(1:rows)
 deallocate(dum)
        
 do i = 1, lines
  read(40,*) dat(lines-i+1,1:rows)
 end do
        
 allocate(m(1:lines),r(1:lines),pres(1:lines),rho(1:lines),ene(1:lines), &
          temp(1:lines))
                    
 do i = 1, rows
!   if(trim(header(i))=='[    Mass   ]') m(1:lines) = dat(1:lines,i)
!   if(trim(header(i))=='[  Density  ]') rho(1:lines) = dat(1:lines,i)
!   if(trim(header(i))=='[   E_int   ]') ene(1:lines) = dat(1:lines,i)
!   if(trim(header(i))=='[   Radius  ]') r(1:lines) = dat(1:lines,i)
!   if(trim(header(i))=='[  Pressure ]') pres(1:lines) = dat(1:lines,i)
!   if(trim(header(i))=='[Temperature]') temp(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='mass_grams') m(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='rho') rho(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='cell_specific_IE') ene(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='radius_cm') r(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='pressure') pres(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='temperature') temp(1:lines) = dat(1:lines,i)
 end do
end subroutine read_mesa

end module setsoftenedcore