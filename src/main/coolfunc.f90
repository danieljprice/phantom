!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: coolfunc
!
!  DESCRIPTION:
!   Implements cooling defined using a tabulated cooling table
!   produced e.g. by CLOUDY
!
!  REFERENCES: 
!   Townsend (2009), ApJS 181, 391-397
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    C_cool    -- factor controlling cooling timestep
!    beta_cool -- beta factor in Gammie (2001) cooling
!    icooling  -- cooling function (0=off, 1=Gammie cooling 2=SD93)
!
!  DEPENDENCIES: h2cooling, infile_utils, io, options, part, timestep
!+
!--------------------------------------------------------------------------
module coolfunc
 implicit none
 integer, parameter :: maxt = 256
 real    :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt)
 integer :: nt
 !
 ! set default values for input parameters
 !
 real :: habund = 0.7
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
subroutine init_coolfunc()
 use io,        only:fatal
 use datafiles, only:find_phantom_datafile
 integer, parameter :: iu = 127
 integer :: ierr,i
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
 if (i==maxt) call fatal('coolfunc','size of cooling table exceeds array size')
 nt = i-1
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
subroutine energ_coolfunc(uu,rho,dt)
 use eos,     only:gamma,gmw
 use physcon, only:atomic_mass_unit,kboltz,Rg
 use units,   only:unit_density,unit_ergg,utime
 real, intent(in)    :: rho,dt
 real, intent(inout) :: uu
 real    :: gam1,density_cgs,dt_cgs,amue,amuh,dtemp,durad
 real    :: sloperef,slopek,temp,temp1,temp_min,tref,yfunx,yinv0
 integer :: k
 
 gam1 = gamma - 1.
 temp = gam1*uu/Rg*gmw*unit_ergg
 temp_min = 10. !gam1*thermal/Rg*gmw*unit_ergg
 
 tref = temper(nt)
 sloperef = slope(nt)
 
 if (temp < temp_min) then
    temp = temp_min
 else
    amue = 2.*atomic_mass_unit/(1. + habund)
    amuh = atomic_mass_unit/habund
    density_cgs = rho*unit_density
    dt_cgs      = dt*utime
    
    dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))* &
            sloperef/tref*dt_cgs
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
       temp1 = max(temper(k)*exp(-lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))),temp_min)    
    else
       yinv0 = 1. - (1. - slopek)*lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))
       if (yinv0 > 0.) then
          temp1 = max(temper(k)*yinv0**(1./(1. - slopek)),temp_min)
       else
          temp1 = temp_min
       endif
    endif
 endif
 
 durad = (temp1 - temp)*Rg/(gam1*gmw*unit_ergg)
 
 uu = uu + durad

end subroutine energ_coolfunc

!-----------------------------------------------------------------------
!+
!  utility to find the index of closest value in a table
!+
!-----------------------------------------------------------------------
integer function find_in_table(n,table,val) result(i)
 integer, intent(in) :: n
 real,    intent(in) :: table(n), val
 integer :: i0,i1
 
 i0 = 0
 i1 = n + 1
 print*,'val=',val,' table = ',table(1),table(nt),' i0,i1=',i0,i1
 do while (i1 - i0 > 1)
    i = (i0 + i1)/2
    if ((table(nt) >= table(1)).eqv.(val >= table(i))) then
       i0 = i
    else
       i1 = i
    endif
    print*,'lims=',i0,i1,' i = ',i
 enddo
 if (abs(val-table(1)) < tiny(0.)) then
    i = 1
 elseif (abs(val-table(nt)) < tiny(0.)) then
    i = n-1
 else
    i = i0
 endif
 read*

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
 case default
    imatch = .false. 
 end select
 if (ngot >= 1) igotall = .true.

end subroutine read_options_coolfunc

end module coolfunc
