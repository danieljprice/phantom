!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: coolfunc
!
!  DESCRIPTION:
!   Implements cooling defined using a tabulated cooling table
!   produced e.g. by CLOUDY. We implement the "exact cooling"
!   method by Townsend to avoid the timestep constraint due
!   to cooling
!
!  REFERENCES:
!   Townsend (2009), ApJS 181, 391-397
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    cooltable  -- data file containing cooling function
!    habund     -- Hydrogen abundance assumed in cooling function
!    temp_floor -- Minimum allowed temperature in K
!
!  DEPENDENCIES: datafiles, eos, infile_utils, io, physcon, units
!+
!--------------------------------------------------------------------------
module coolfunc
 implicit none
 integer, parameter :: maxt = 1000
 real    :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt)
 integer :: nt
 !
 ! set default values for input parameters
 !
 real :: habund     = 0.7
 real :: temp_floor = 1.e4
 character(len=120) :: cooltable = 'cooltable.dat'

 public :: init_coolfunc,write_options_coolfunc,read_options_coolfunc
 public :: energ_coolfunc

 public :: find_in_table

 private

contains

!-----------------------------------------------------------------------
!+
!  read cooling table from file and initialise arrays
!+
!-----------------------------------------------------------------------
subroutine init_coolfunc(ierr)
 use io,        only:fatal
 use datafiles, only:find_phantom_datafile
 integer, intent(out) :: ierr
 integer, parameter :: iu = 127
 integer :: i
 character(len=120) :: filepath

 !
 ! read the cooling table from file
 !
 filepath=find_phantom_datafile(cooltable,'cooling')
 open(unit=iu,file=filepath,status='old',iostat=ierr)
 if (ierr /= 0) call fatal('coolfunc','error opening cooling table')
 i = 0
 do while(ierr==0 .and. i < maxt)
    i = i + 1
    read(iu,*,iostat=ierr) temper(i),lambda(i)
 enddo
 nt = i-1
 if (nt==maxt) call fatal('coolfunc','size of cooling table exceeds array size')
 if (nt < 2) call fatal('coolfunc','size of cooling table is too small',ival=nt,var='nt')
 !
 ! calculate the slope of the cooling function
 !
 do i=1,nt-1
    slope(i) = log(lambda(i+1)/lambda(i))/log(temper(i+1)/temper(i))
 enddo
 slope(nt) = slope(nt-1)

 !
 ! initialise the functions required for Townsend method
 !
 yfunc(nt) = 0.
 do i=nt-1,1,-1
    if (abs(slope(i)-1.) < tiny(0.)) then
       yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/(lambda(i)*temper(nt))*log(temper(i)/temper(i+1))
    else
       yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/((1. - slope(i))*lambda(i)*temper(nt))&
                 *(1.- (temper(i)/temper(i+1))**(slope(i) - 1.))
    endif
 enddo

end subroutine init_coolfunc

!-----------------------------------------------------------------------
!+
!  implement du/dt equation
!+
!-----------------------------------------------------------------------
subroutine energ_coolfunc(uu,rho,dt,dudt)
 use eos,     only:gamma,gmw
 use physcon, only:atomic_mass_unit,kboltz,Rg
 use units,   only:unit_density,unit_ergg,utime
 real, intent(in)    :: rho,dt
 real, intent(inout) :: uu,dudt
 real    :: gam1,density_cgs,dt_cgs,amue,amuh,dtemp,durad
 real    :: sloperef,slopek,temp,temp1,tref,yfunx,yinv0
 integer :: k

 gam1 = gamma - 1.
 temp = gam1*uu/Rg*gmw*unit_ergg

 tref     = temper(nt)
 sloperef = slope(nt)

 if (temp < temp_floor) then
    temp1 = temp_floor
 else
    amue = 2.*atomic_mass_unit/(1. + habund)
    amuh = atomic_mass_unit/habund
    density_cgs = rho*unit_density
    dt_cgs      = dt*utime

!original version
!    dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))* &
!         sloperef/tref*dt_cgs
    print *,'check coolfunc'
!Lionel Siess : I think there is an error. dtemp should write (sloperef <-> lambda(nt)
     dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))* &
         lambda(nt)/tref*dt_cgs

    k = find_in_table(nt,temper,temp)

    slopek = slope(k)
    if (abs(slopek - 1.) < tiny(0.)) then
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt))*log(temper(k)/temp)
    else
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt)*(1. - slopek)) &
                          *(1. - (temper(k)/temp)**(slopek-1.))
    endif
    yfunx = yfunx + dtemp

    if (abs(slopek - 1.) < tiny(0.)) then
       temp1 = max(temper(k)*exp(-lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))),temp_floor)
    else
       yinv0 = 1. - (1. - slopek)*lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))
       if (yinv0 > 0.) then
          temp1 = max(temper(k)*yinv0**(1./(1. - slopek)),temp_floor)
       else
          temp1 = temp_floor
       endif
    endif
 endif

 durad = (temp1 - temp)*Rg/(gam1*gmw*unit_ergg)
 !dudt = dudt + durad/dt
 uu = uu + durad

end subroutine energ_coolfunc

!-----------------------------------------------------------------------
!+
!  utility to find the index of closest value in a table
!+
!-----------------------------------------------------------------------
pure integer function find_in_table(n,table,val) result(i)
 integer, intent(in) :: n
 real,    intent(in) :: table(n), val
 integer :: i0,i1

 i0 = 0
 i1 = n + 1
 do while (i1 - i0 > 1)
    i = (i0 + i1)/2
    if ((table(n) >= table(1)).eqv.(val >= table(i))) then
       i0 = i
    else
       i1 = i
    endif
 enddo
 if (abs(val-table(1)) < tiny(0.)) then
    i = 1
 elseif (abs(val-table(n)) < tiny(0.)) then
    i = n-1
 else
    i = i0
 endif

end function find_in_table

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_coolfunc(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(cooltable,'cooltable','data file containing cooling function',iunit)
 call write_inopt(habund,'habund','Hydrogen abundance assumed in cooling function',iunit)
 call write_inopt(temp_floor,'temp_floor','Minimum allowed temperature in K',iunit)

end subroutine write_options_coolfunc

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_coolfunc(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .false.  ! cooling options are compulsory
 select case(trim(name))
 case('cooltable')
    read(valstring,*,iostat=ierr) cooltable
    ngot = ngot + 1
 case('habund')
    read(valstring,*,iostat=ierr) habund
    ngot = ngot + 1
 case('temp_floor')
    read(valstring,*,iostat=ierr) temp_floor
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if (ngot >= 3) igotall = .true.

end subroutine read_options_coolfunc

end module coolfunc
