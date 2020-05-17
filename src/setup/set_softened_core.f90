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
 use physcon,          only:pi,gg
 use eos_idealplusrad, only:get_idealgasplusrad_tempfrompres,&
                            get_idealplusrad_enfromtemp
 implicit none
 real(kind=8)  :: hsoft,msoft,mcore
 integer       :: hidx

 ! hsoft: Softening length of core particle
 ! mcore: Mass of core particle
 ! msoft: Softened mass (mass at softening length minus mass of core
 !        particle)
 ! hphi:  Softening length for the point particle potential, defined in
 !        Price & Monaghan (2006). Set to be 0.5*hsoft. 

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
 real(kind=8), allocatable, dimension(:) :: m0(:),rho0(:),pres0(:),&
                                            m(:),rho(:),pres(:),r(:),phi(:),&
                                            ene(:),temp(:),drho(:)
 character(len=120), intent(in)          :: filepath,outputpath
 real(kind=8), intent(inout)             :: mcore,hsoft
 real(kind=8)                            :: mc,mh,h,hphi,tolerance,tempi,eni
 logical                                 :: isort_decreasing,iexclude_core_mass
 integer                                 :: i

 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .true. ! Needs to be true if to be read by Phantom
 !
 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .true. ! Needs to be true if to be read by Phantom

 call read_mesa(rho0, r, pres0, m0, ene, temp, filepath)
 tolerance = 1.5 ! How much we allow the softened density to exceed the original profile by

 ! We iteratively look for a value of mcore for given hsoft that produces a nice softened
 ! density profile, and vice-versa. We also check for errors.
 !
 ! Case 1 (default): hsoft specified only
 !
 if ((hsoft > 0) .and. (mcore < 0)) then ! Default case
     h = hsoft * solarr ! Convert to cm
     call interpolator(r, h, hidx) ! Find index in r closest to h
     mc = 0.7*m0(hidx) ! Initialise profile to have very large softened mass
     do
         msoft = m0(hidx) - mc
         rho = rho0 ! Reset density
         m = m0 ! Reset mass
         call calc_rho_and_m(rho, m, r, mc, h)
         call diff(rho, drho)
         if (all(rho/rho0 < tolerance) .and. all(drho(1:hidx) < 0)) exit
         if (mc > 0.98*m0(hidx)) then
             stop 'ERROR: Cannot find mcore that produces nice profile (mcore/m(h) > 0.98 reached)'
         endif
         mc = mc + 0.01*m0(hidx) ! Increase mcore/m(h) by 1 percent
     enddo
     ! Write out mcore
     mcore = mc / solarm
 !
 ! Case 2: mcore specified only
 !
 elseif ((hsoft < 0) .and. (mcore > 0)) then
     mc = mcore * solarm ! Convert to g
     mh = mc / 0.7 ! Initialise h such that m(h) to be much larger than mcore
     do
         call interpolator(m0, mh, hidx) 
         h = r(hidx) 
         msoft = mh - mc
         rho = rho0 ! Reset density
         m = m0 ! Reset mass
         call calc_rho_and_m(rho, m, r, mc, h)
         call diff(rho, drho)
         if (all(rho/rho0 < tolerance) .and. all(drho(1:hidx) < 0)) exit
         if (mc > 0.98*m0(hidx)) then
             stop 'ERROR: Cannot find softening length that produces nice profile (h/r(mcore) < 1.02 reached)'
         endif
         call interpolator(m0, 1.3*mc, hidx)
         mh = 1./(1./mh + 0.01/mc) ! Increase mcore/m(h) by 1 percent
     enddo
     ! Write out hsoft
     hsoft = h / solarr
 !
 ! Case 3: Both hsoft and mcore specified
 !
 elseif ((hsoft > 0) .and. (mcore > 0)) then
     h = hsoft * solarr ! Convert to cm
     mc = mcore * solarm ! Convert to g
     call interpolator(r, h, hidx) ! Find index in r closest to h
     msoft = m0(hidx) - mc
     ! Check for sensible choices
     if (msoft < 0) then
         print*,'mcore = ',mcore,', m(r=h) = ',m0(hidx)/solarm
         stop 'ERROR: mcore cannot exceed m(r=h)'
     endif
     rho = rho0
     m = m0
     call calc_rho_and_m(rho, m, r, mc, h)
     ! Test if profile is sensible
     if (any(rho/rho0 > tolerance)) then 
         print*,'Warning: softenedrho/rho > tolerance'
     endif
     call diff(rho, drho)
     if (any(drho(1:hidx) > 0)) then
         stop 'ERROR: drho/dr > 0 found in softened profile'
     endif

 else
     stop 'ERROR: Neither hsoft nor mcore were specified.'
 endif
 !
 ! Calculate gravitational potential
 !
 hphi = 0.5*h
 call calc_phi(r, m-mc, phi, mc, hphi)
 !
 ! Calculate pressure profile inside softening length
 !
 pres = pres0
 call calc_pres(r, rho, phi, pres)
 !
 ! Calculate temperature profile from pressure and density
 ! (only implemented for ideal gas plus radiation pressure
 ! EoS for now)
 !
 temp(size(r)) = 0.
 do i = 1,size(r)-1
     tempi = temp(i)
     call get_idealgasplusrad_tempfrompres(pres(i), rho(i), tempi)
     temp(i) = tempi
 enddo
 !
 ! Calculate internal energy per unit mass
 ! (only implemented for ideal gas plus radiation pressure
 ! EoS for now)
 !
 do i = 1,size(r)
    call get_idealplusrad_enfromtemp(rho(i), temp(i), eni)
    ene(i) = eni
 enddo
 !
 ! Write data
 ! Note: The temperature and internal energy are fake, since they are from the original mesa file
 !
 ! Reverse arrays so that data is sorted from stellar surface to stellar centre.
 !
 if (isort_decreasing) then
     call flip_array(m)
     call flip_array(pres)
     call flip_array(temp)
     call flip_array(r)
     call flip_array(rho)
     call flip_array(ene)
     call flip_array(phi)
 end if

 if (iexclude_core_mass) then
     m = m - mc
 endif

 call write_softened_profile(outputpath, m, pres, temp, r, rho, ene)
end subroutine set_softened_core


subroutine calc_rho_and_m(rho,m,r,mc,h)
 implicit none
 real(kind=8) :: a, b, d, drhodr_h
 real(kind=8), intent(in) :: mc, h
 real(kind=8), dimension(1:hidx+1), intent(in) :: r
 real(kind=8), dimension(1:hidx+1), intent(inout) :: rho, m
 
! a, b, d: Coefficients of cubic density profile defined by rho(r) = ar**3 + br**2 + d
 drhodr_h = (rho(hidx+1) - rho(hidx)) / (r(hidx+1) - r(hidx)) ! drho/dr at r = h
 a = 2./h**2. * drhodr_h - 10./h**3. * rho(hidx) + 7.5/pi/h**6. * msoft
 b = 0.5*drhodr_h/h - 1.5*a*h
 d = rho(hidx) - 0.5*h*drhodr_h + 0.5*a*h**3.

 rho(1:hidx) = a*r(1:hidx)**3. + b*r(1:hidx)**2. + d

 ! Mass is then given by m(r) = mcore + 4*pi (1/6 a r^6 + 1/5 b r^5 + 1/3 d r^3)
 m(1:hidx) = mc + 4.*pi * (1./6. * a * r(1:hidx)**6. + 0.2 * b * r(1:hidx)**5. + &
                           1./3. * d * r(1:hidx)**3.)
end subroutine calc_rho_and_m


subroutine calc_phi(r,mgas,phi,mc,hphi)
 implicit none
 real(kind=8), dimension(:), intent(in) :: r(:), mgas(:)
 real(kind=8), dimension(:), allocatable :: q
 real(kind=8), intent(in) :: mc, hphi
 real(kind=8), dimension(:), allocatable, intent(out) :: phi
 real(kind=8), dimension(size(r)) :: phi_core, phi_gas
 integer :: idx2hphi, idxhphi, i
 ! The gravitational potential is needed to integrate the pressure profile using the
 ! equation of hydrostatic equilibrium. First calculate gravitational potential due
 ! to point mass core, according to the cubic spline kernel in Price & Monaghan (2006)
 ! to be consistent with Phantom, then calculate the gravitational potential due to the
 ! softened gas.
 allocate(phi(size(r)))
 ! (i) Gravitational potential due to core particle (cubic spline softening)
 ! For 0 <= r/hphi < 1
 call interpolator(r, hphi, idxhphi) ! Find index corresponding to r = 2*hphi
 allocate(q(1:idxhphi-1))
 q = r(1:idxhphi-1) / hphi
 phi_core(1:idxhphi-1) = gg*mc/hphi * (2./3.*q**2. - 0.3*q**4. + 0.1*q**5. &
                                   - 7./5.)
 deallocate(q)

 ! For 1 <= r/hphi < 2
 call interpolator(r, 2*hphi, idx2hphi) ! Find index corresponding to r = 2*hphi
 allocate(q(idxhphi:idx2hphi-1))
 q = r(idxhphi:idx2hphi-1) / hphi
 phi_core(idxhphi:idx2hphi-1) = gg*mc/hphi * (4./3.*q**2. - q**3. + 0.3*q**4. &
                                        - 1./30.*q**5. - 1.6 + 1./15./q)
 deallocate(q)

 ! For 2 <= r/hphi
 phi_core(idx2hphi:size(r)) = - gg * mc / r(idx2hphi:size(r))

 ! (ii) Gravitational potential due to softened gas
 phi_gas(size(r)) = - gg * mgas(size(r)) / r(size(r)) ! Surface boundary condition for phi
 do i = 1, size(r) - 1
     phi_gas(size(r)-i) = phi_gas(size(r)-i+1) - gg * mgas(size(r)-i) / r(size(r)-i)**2. &
                                               * (r(size(r)-i+1) - r(size(r)-i))
 end do
      
 ! (iii) Add the potentials 
 phi = phi_gas + phi_core 
end subroutine calc_phi


subroutine calc_pres(r, rho, phi, pres)
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


subroutine diff(array, darray)
 implicit none
 real(kind=8), dimension(:), intent(in)               :: array(:)
 real(kind=8), dimension(:), allocatable, intent(out) :: darray(:)
 integer                                              :: i
 ! A subroutine that takes a n-dim array and produces a (n-1)-dim array
 ! with the ith element is the (i+1)th element minus the ith element of
 ! the original array.
 allocate(darray(size(array)-1))
 do i = 1, size(array)-1
    darray(i) = array(i+1) - array(i)
 end do
end subroutine diff


subroutine write_softened_profile(outputpath, m, pres, temp, r, rho, ene)
 implicit none
 real(kind=8), allocatable, dimension(:) :: m(:),rho(:),pres(:),r(:),ene(:),temp(:)
 character(len=120), intent(in)          :: outputpath
 integer                                 :: i
 open(1, file = outputpath, status = 'new')  
 write(1,'(a)') '[    Mass   ]  [  Pressure ]  [Temperature]  [   Radius  ]  [  Density  ]  [   E_int   ]'
 write(1,42) (m(i), pres(i), temp(i), r(i), rho(i), ene(i), i = 1, size(r))
 42 format (es13.7, 2x, es13.7, 2x, es13.7, 2x, es13.7, 2x, es13.7, 2x, es13.7)
 close(1, status = 'keep')
end subroutine write_softened_profile


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