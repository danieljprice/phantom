!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: photoevap
!
!  DESCRIPTION: This module contains all the subroutines necessary for the
!     photoevaporation switch
!
!  REFERENCES: Alexander, Clarke & Pringle (2006), MNRAS 369, 216-228
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    ionflux_cgs    -- Stellar EUV flux rate
!    mu_cgs         -- Mean molecular weight
!    recombrate_cgs -- Recombination rate (alpha)
!
!  DEPENDENCIES: allocutils, dim, eos, externalforces, infile_utils,
!    physcon, units
!+
!--------------------------------------------------------------------------
module photoevap

 implicit none

 !--# of grid nodes for photoevaporation grid
 integer, parameter :: Nr = 400
 integer, parameter :: Nphi = 400
 integer, parameter :: Ntheta = 400

 !--Index to identify which cell particles belong
 integer, allocatable :: Rnum(:)
 integer, allocatable :: Thetanum(:)
 integer, allocatable :: Phinum(:)

 !--# of particles per cell and ray.
 integer :: Cellpartnum(Nr-1,Ntheta-1,Nphi-1)
 integer :: Raypartnum(Ntheta-1,Nphi-1)

 !--Location of ionization front and # of ionized particles in ray
 integer :: Ionfront(Ntheta-1,Nphi-1)
 integer :: Nion(Ntheta-1,Nphi-1)

 !--Reciprical of cell volume
 real    :: rCellvol(Nr-1,Nphi-1)

 !--Fraction of ions to neutrals in boundary cell
 real    :: Ionfrac(Ntheta-1,Nphi-1)

 !--Change in # of ions per unit time due to stellar flux (constant in time)
 real    :: dN_ion(Ntheta-1)

 !--Photoevaporation grid minimums and increment values
 real    :: Rgrid_min,  Thetagrid_min,  Phigrid_min
 real    :: dr_grid,    dphi_grid,      dtheta_grid

 !--Recombination rate, ionization flux, mean Mol. weight, temperature of ions
 real    :: recombrate
 real    :: recombrate_cgs = 2.6d-13
 real    :: ionflux
 real    :: ionflux_cgs = 1.d41
 real    :: mu
 real    :: mu_cgs = 1.26
 real    :: mH
 real    :: temp_ion
 real    :: energy_to_temperature_ratio

 real    :: prev_time

 public  :: allocate_photoevap
 public  :: deallocate_photoevap
 public  :: write_options_photoevap
 public  :: read_options_photoevap
 public  :: photo_ionize
 public  :: find_ionfront
 public  :: set_photoevap_grid

 private

contains

!***************************************************************************************
!***************************************************************************************

subroutine allocate_photoevap
 use dim, only:maxp
 use allocutils, only:allocate_array

 call allocate_array('Rnum', Rnum, maxp)
 call allocate_array('Thetanum', Thetanum, maxp)
 call allocate_array('Phinum', Phinum, maxp)
end subroutine allocate_photoevap

subroutine deallocate_photoevap
 deallocate(Rnum)
 deallocate(Thetanum)
 deallocate(Phinum)
end subroutine deallocate_photoevap

!----------------------------------------------------------------
!+
!  This subroutine makes a spherical grid for photoevaporation
!  Note: this routine is ment to be called only once at the
!        beginning of the simulation to get the grid spacing.
!+
!----------------------------------------------------------------
subroutine set_photoevap_grid
 use units,          only:udist,umass,utime
 use physcon,        only:pi,atomic_mass_unit,mass_proton_cgs,kboltz,Rg
 use externalforces, only:accradius1
 use eos,            only:gamma

 integer :: i,j

 !--Inner and outer radius of grid
 real    :: R_in
 real    :: R_out

 !--Photoevaporation grid min and max in each direction
 real    :: Rgrid_max
 real    :: Thetagrid_max
 real    :: Phigrid_max

 !--photoevaporation grid in r,theta directions
 real    :: r_grid(Nr)
 real    :: theta_grid(Ntheta)

 !--TODO: try to read these in from setup_photoevap
 R_in  = accradius1
 R_out = 10.0

 !--Set the temperature of ions to 10,000 K
 temp_ion = 1.d4

 !--Constant that converts specific energy density (u) to temperature in K.
 !  Note: the (udist/utime)**2 comes from adding physical units back onto u.
 !energy_to_temperature_ratio = Rg/(mu*(gamma-1.))/(udist/utime)**2
 energy_to_temperature_ratio = kboltz/((gamma-1.)*mu*atomic_mass_unit)*(utime/udist)**2

 !--Mass of hydrogen gas molecule in code units
 mH = mu*mass_proton_cgs/umass

 !--Time at previous time step (Initialized to zero here)
 prev_time = 0.

 !--photoevaporation grid's inner/outer radius and other dimensions
 Rgrid_min     = R_in
 Rgrid_max     = 1.3*R_out
 Thetagrid_min = 0.
 Thetagrid_max = pi
 Phigrid_min   = -pi
 Phigrid_max   = pi

 !--photoevaporation grid spacing in each direction (with Nr x Nphi x Ntheta nodes)
 dr_grid     = (Rgrid_max     - Rgrid_min    )/(Nr-1)
 dtheta_grid = (Thetagrid_max - Thetagrid_min)/(Ntheta-1)
 dphi_grid   = (Phigrid_max   - Phigrid_min  )/(Nphi-1)

 !--photoevaporation grid in r,theta directions
 r_grid      = (/( Rgrid_min     + dr_grid*i    , i = 0, Nr-1     )/)
 theta_grid  = (/( Thetagrid_min + dtheta_grid*i, i = 0, Ntheta-1 )/)

 do i = 1,Ntheta-1
    do j = 1,Nr-1
       !--Calculate the volume of each grid cell (symmetrical in phi so only need two dimensional array)
       rCellvol(j,i) = (0.5*(r_grid(j)+r_grid(j+1)))**2*sin(0.5*(theta_grid(i)+theta_grid(i+1)))*dr_grid*dtheta_grid*dphi_grid
    enddo

    !--Calculate the ionization rate in each ray (symmetrical in phi so only need one dimensional array)
    dN_ion(i) = (0.25/pi*sin(0.5*(theta_grid(j)+theta_grid(j+1)))*dtheta_grid*dphi_grid)*ionflux
 enddo

 !--Because the reciprical of Cellvol is only ever used: rCellvol = 1/Cellvol
 rCellvol = 1./rCellvol

 return
end subroutine set_photoevap_grid

!-----------------------------------------------------------------------
!+
!  Subroutine to identify in which grid cell particles reside, solve the
!  ionization/recombination balance for each ray, and finally find the
!  location for the ionization front.
!+
!-----------------------------------------------------------------------
subroutine find_ionfront(timei,npart,xyzh,pmassi)

 integer, intent(in) :: npart
 real,    intent(in) :: timei
 real,    intent(in) :: pmassi
 real,    intent(in) :: xyzh(:,:)

 integer :: i,j,k

 !--Cumulative sum of particles along ray just below current cell
 integer :: curr_ray_count

 !--Change in # of ions per unit time due to recombination rate
 real    :: dN_recomb

 !--Position of particle in spherical coordinates
 real    :: r_pos,theta_pos,phi_pos

 real    :: pmass_on_mH
 real    :: dt

 !--Find how much time has elapsed since the last call
 dt = timei - prev_time

 if ( dt == 0. ) then
    print*,'WARNING! find_ionfront was called needlessly!'
 else

    !--Gives the number of hydrogen gas molecules per SPH particle
    pmass_on_mH = pmassi/mH

!$omp parallel do private(i,r_pos,theta_pos,phi_pos) schedule(static)
    do i = 1,npart
       r_pos     = sqrt(sum(xyzh(1:3,i)**2))
       theta_pos = acos(xyzh(3,i)/r_pos)
       phi_pos   = atan2(xyzh(2,i),xyzh(1,i))

       ! Find the (*INTEGER*) grid node just below the particle in each direction
       Rnum(i)     = int((r_pos-Rgrid_min)/dr_grid)+1
       Thetanum(i) = int((theta_pos-Thetagrid_min)/dtheta_grid)+1
       Phinum(i)   = int((phi_pos-Phigrid_min)/dphi_grid)+1
    enddo
!$omp end parallel do

    !--Re-initialize/re-calculate Cellpartnum every time step
    Cellpartnum = 0
!$omp parallel do private(i) schedule(static)
    do i = 1,npart
       Cellpartnum(Rnum(i),Thetanum(i),Phinum(i)) = Cellpartnum(Rnum(i),Thetanum(i),Phinum(i)) + 1
    enddo
!$omp end parallel do

    !--Find the total number of particles along each ray (used to speed up loop below if ray is empty)
    Raypartnum(:,:) = sum(Cellpartnum,1)

    !
    !--Solve for ionization/recombination balance and update Nion, ionization
    !  front, and fraction of ionized particles at the front
    !
!$omp parallel do default(none)   &
!$omp private(i,j,k,dN_recomb,curr_ray_count) &
!$omp shared(Ionfrac,Ionfront,Raypartnum,Cellpartnum,rCellvol,recombrate,dN_ion,Nion,dt,pmass_on_mH) &
!$omp schedule(dynamic)
    do i = 1,Nphi-1
       do j = 1,Ntheta-1
          !--Save radial location of ionfront for current ray in k
          k = Ionfront(j,i)

          !--Find the change in Nion due to recombination
          if ( k == 1 ) then
             dN_recomb = Ionfrac(j,i)*Cellpartnum(k,j,i)**2*rCellvol(k,j)
          else
             dN_recomb = sum(Cellpartnum(1:k-1,j,i)**2*rCellvol(1:k-1,j)) + Ionfrac(j,i)*Cellpartnum(k,j,i)**2*rCellvol(k,j)
          endif
          dN_recomb = recombrate*dN_recomb*pmass_on_mH

          !--Update the # of ionized particles in each radial column
          if ( Raypartnum(j,i) > 0 ) then
             Nion(j,i) = Nion(j,i) + nint(dt*(dN_ion(j)/pmass_on_mH-dN_recomb))

             !--Make sure that flux doesn't "build-up" in the fully ionized columns
             !  (i.e. the excess light escapes the system)
             if ( Nion(j,i) > Raypartnum(j,i) ) then
                Nion(j,i) = Raypartnum(j,i)
             endif

!!!!!!
!       print*,dN_ion(j)/pmass_on_mH,dN_recomb,Nion(j,i)
!!!!!!

             if ( Nion(j,i) < 0  ) then
                print*,'Warning! Negative ion number!',Nion(j,i),j,i
                Nion(j,i) = 0
             endif
          else
             Nion(j,i) = 0 !--If no particles, then Nion must be reset
          endif

          !--Now that we have the # of ions in each column, integrate from the star
          !  out to Nion to find where the ionization front is located
          k = 1
          curr_ray_count = Cellpartnum(k,j,i)
          do while ( curr_ray_count < Nion(j,i) )
             if ( k < Nr-1 ) then
                k = k+1
             else
                exit
             endif
             curr_ray_count = curr_ray_count + Cellpartnum(k,j,i)
          enddo

          !--Save the new ionization front radial cell # for the next iteration
          Ionfront(j,i) = k

          !--Find the fraction of ions to neutrals in the ionization front
          !  This only needs to be done for cells with more particles than Nion
          if ( Raypartnum(j,i) <= Nion(j,i) ) then
             Ionfrac(j,i) = 1.
          else
             if ( Cellpartnum(k,j,i) == 0 ) then
                Ionfrac(j,i) = 1.
             else
                Ionfrac(j,i) = (Nion(j,i)-(curr_ray_count-Cellpartnum(k,j,i)))/real(Cellpartnum(k,j,i))
             endif
          endif

          if ( Ionfrac(j,i) < 0 ) then
             print*,'Ionfrac is less than zero!'
             stop
          elseif ( Ionfrac(j,i) > 1 ) then
             print*,'Ionfrac is greater than 1!'
             stop
          endif

       enddo
    enddo
!$omp end parallel do

    prev_time = timei
 endif

 return
end subroutine find_ionfront

!-----------------------------------------------------------------------
!+
!  Update the temperatures of the particles (Ionized,Boundary,Neutral)
!+
!-----------------------------------------------------------------------
subroutine photo_ionize(vxyzu,npart)
 integer, intent(in)    :: npart
 real,    intent(inout) :: vxyzu(:,:)

 integer  :: ipart

 !--Temperature of the particle in Kelvin
 real     :: temperature,tempi

 if ( size(vxyzu,1) /= 4 ) then
    print*,'Fatal error! vxyzu does not have an energy column!'
    stop
 endif

!$omp parallel do default(none)          &
!$omp private(ipart,tempi,temperature)   &
!$omp shared(npart,vxyzu,energy_to_temperature_ratio,Rnum,Thetanum,Phinum,Ionfrac,Ionfront,temp_ion)
 do ipart = 1,npart
    tempi = vxyzu(4,ipart)/energy_to_temperature_ratio

    if ( Rnum(ipart) < Ionfront(Thetanum(ipart),Phinum(ipart)).or.          &
        Rnum(ipart) >= Nr                                          ) then
       ! Above ionization front (ionized particles)
       temperature = temp_ion

    elseif ( Rnum(ipart) == Ionfront(Thetanum(ipart),Phinum(ipart)) ) then
       ! Ionization front (fractionally ionized)
       temperature =     Ionfrac(Thetanum(ipart),Phinum(ipart))*temp_ion +   &
                   (1.-Ionfrac(Thetanum(ipart),Phinum(ipart)))*tempi

    elseif ( Rnum(ipart) > Ionfront(Thetanum(ipart),Phinum(ipart)) ) then
       ! Below ionization front (neutral particles)
       temperature = tempi

    endif

    vxyzu(4,ipart)     = temperature*energy_to_temperature_ratio

!   if ( vxyzu(4,ipart) < 0 ) then
!     print*,'vxyzu is negative!', vxyzu(4,ipart),ipart
!     stop
!   endif

 enddo
!$omp end parallel do

! do ipart = 1,npart
!    print*,vxyzu(4,ipart)/energy_to_temperature_ratio
! enddo

 return
end subroutine photo_ionize

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_photoevap(iunit)
 use infile_utils, only:write_inopt

 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling photoevaporation'
 call write_inopt(mu_cgs,'mu_cgs','Mean molecular weight',iunit)
 call write_inopt(recombrate_cgs,'recombrate_cgs','Recombination rate (alpha)',iunit)
 call write_inopt(ionflux_cgs,'ionflux_cgs','Stellar EUV flux rate',iunit)

end subroutine write_options_photoevap

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_photoevap(name,valstring,imatch,igotall,ierr)
 use units,    only:udist,utime

 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot = 0

 imatch = .false.
 igotall = .true.
 ierr = 0

 select case(trim(name))
 case('mu_cgs')
    read(valstring,*,iostat=ierr) mu_cgs
    mu = mu_cgs
    ngot = ngot + 1
 case('recombrate_cgs')
    read(valstring,*,iostat=ierr) recombrate_cgs
    recombrate = recombrate_cgs*utime/udist**3
    ngot = ngot + 1
 case('ionflux_cgs')
    read(valstring,*,iostat=ierr) ionflux_cgs
    ionflux = ionflux_cgs*utime
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 igotall = ( ngot >= 3 )

end subroutine read_options_photoevap

end module photoevap
