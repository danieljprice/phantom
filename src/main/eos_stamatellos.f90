!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_stamatellos
!
! eos_stamatellos
!
! :References: None
!
! :Owner: Alison Young
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, datafiles, dim, io, physcon, units
!

 implicit none
 real,allocatable,public :: optable(:,:,:)
 real,allocatable,public :: Gpot_cool(:),duFLD(:),gradP_cool(:),lambda_FLD(:),urad_FLD(:) !gradP_cool=gradP/rho
 real,allocatable,public :: ttherm_store(:),ueqi_store(:),duSPH(:)
 real,allocatable,public :: du_store(:),tau_store(:) ! Only saved to write to dumps
 character(len=25), public :: eos_file= 'eos_lom.dat' !default name of tabulated EOS file
 logical,public :: doFLD = .True., floor_energy = .False.
 integer,public :: iunitst=19
 integer,save :: nx,ny ! dimensions of optable read in

 public :: read_optab,getopac_opdep,init_coolra,getintenerg_opdep,finish_coolra
 public :: get_k_fld

contains

subroutine init_coolra()
 use dim, only:maxp,maxp_alloc
 use allocutils, only:allocate_array
 integer :: np
 if (maxp == 0) then
    np = int(maxp_alloc)
 else
    np = maxp
 endif
 print *, "Allocating cooling arrays for np=",np
 call allocate_array('gradP_cool',gradP_cool,np)
 call allocate_array('Gpot_cool',Gpot_cool,np)
 call allocate_array('duFLD',duFLD,np)
 call allocate_array('lambda_fld',lambda_fld,np)
 call allocate_array('urad_FLD',urad_FLD,np)
 call allocate_array('duSPH',duSPH,np)
 if (.not. allocated(ttherm_store)) then
    call allocate_array('ttherm_store',ttherm_store,np)
    call allocate_array('ueqi_store',ueqi_store,np)
    call allocate_array('tau_store',tau_store,np)
    call allocate_array('du_store',du_store,np)
 endif

 Gpot_cool(:) = 0d0
 gradP_cool(:) = 0d0
 urad_FLD(:) = 0d0
 duFLD(:) = 0d0
 ueqi_store(:) = 0d0
 ttherm_store(:) = 0d0
 tau_store(:) = 0d0
 du_store(:) = 0d0
 duSPH(:) = 0d0
 !open(unit=iunitst,file='EOSinfo.dat',status='replace')
 if (doFLD) then
    print *, "Using Forgan+ 2009 hybrid cooling method (FLD)"
 else
    print *, "NOT using FLD. Using cooling only"
 endif
end subroutine init_coolra

subroutine finish_coolra()
 if (allocated(optable)) deallocate(optable)
 if (allocated(gradP_cool)) deallocate(gradP_cool)
 if (allocated(Gpot_cool)) deallocate(Gpot_cool)
 if (allocated(duFLD)) deallocate(duFLD)
 if (allocated(lambda_fld)) deallocate(lambda_fld)
 if (allocated(urad_FLD)) deallocate(urad_FLD)
 if (allocated(ttherm_store)) deallocate(ttherm_store)
 if (allocated(ueqi_store)) deallocate(ueqi_store)
 if (allocated(tau_store)) deallocate(tau_store)
 if (allocated(du_store)) deallocate(du_store)
 if (allocated(duSPH)) deallocate(duSPH)
! close(iunitst)
end subroutine finish_coolra

subroutine read_optab(eos_file,ierr)
 use datafiles, only:find_phantom_datafile
 character(len=*),intent(in) :: eos_file
 integer, intent(out) :: ierr
 integer :: i,j,errread
 character(len=120) :: filepath,junk

 ! read in data file for interpolation
 filepath=find_phantom_datafile(eos_file,'eos/lombardi')
 print *,"EOS file: FILEPATH:",filepath
 open(10,file=filepath,form="formatted",status="old",iostat=ierr)
 if (ierr > 0) return
 do
    read(10,'(A120)') junk
!    print *, junk
    if (len(trim(adjustl(junk))) == 0) cycle ! blank line
    if ((index(adjustl(junk),"::") == 0) .and. (index(adjustl(junk),"#")  /=  1 )) then !ignore comment lines
       junk = adjustl(junk)
       read(junk, *,iostat=errread) nx, ny
       exit
    endif
 enddo
! allocate array for storing opacity tables
 allocate(optable(nx,ny,6))
 do i = 1,nx
    do j = 1,ny
       read(10,*) OPTABLE(i,j,1),OPTABLE(i,j,2),OPTABLE(i,j,3),&
              OPTABLE(i,j,4),OPTABLE(i,j,5),OPTABLE(i,j,6)
    enddo
 enddo
end subroutine read_optab

!
! Main subroutine for interpolating tables to get EOS values
!
subroutine getopac_opdep(ui,rhoi,kappaBar,kappaPart,Ti,gmwi)
 use io, only:warning
 real, intent(in)  :: ui,rhoi
 real, intent(out) :: kappaBar,kappaPart,Ti,gmwi

 integer :: irho,iu
 real :: m,c
 real :: kbar1,kbar2
 real :: kappa1,kappa2
 real :: Tpart1,Tpart2
 real :: gmw1,gmw2
 real :: rhomin,umin

 rhomin = OPTABLE(1,1,1)
 umin = OPTABLE(1,1,3)
 ! interpolate through OPTABLE to find corresponding kappaBar, kappaPart and T
 ! check values are in range of tables
 if (rhoi > OPTABLE(nx,1,1) ) then
    call warning('getopac_opdep','rhoi above range of EOS table.',var='rhoi',val=rhoi)
 elseif  (rhoi < rhomin) then    
    call warning('getopac_opdep','rhoi below range of EOS table.',var='rhoi',val=rhoi)
 elseif (ui > OPTABLE(1,ny,3) .or. ui < umin) then
    call warning('getopac_opdep','ui out of range',var='ui',val=ui)
 endif

 ! Find index of rhoi in table such that array(ind) < rhoi < array(ind+1)
! print *, "search for rhoi"
 irho = search_table(optable(:,1,1),nx,rhoi)
 iu = search_table(optable(irho,:,3),ny,ui)

 m = (optable(irho,iu,5)-optable(irho,iu+1,5))/(optable(irho,iu,3)-optable(irho,iu+1,3))
 c = optable(irho,iu+1,5) - m*optable(irho,iu+1,3)
 kbar1 = m*ui + c

 m = (optable(irho,iu,6)-optable(irho,iu+1,6))/(optable(irho,iu,3)-optable(irho,iu+1,3))
 c = optable(irho,iu+1,6) - m*optable(irho,iu+1,3)
 kappa1 = m*ui + c

 m = (optable(irho,iu,2) - optable(irho,iu+1,2))/(optable(irho,iu,3)-optable(irho,iu+1,3))
 c = optable(irho,iu+1,2) - m*optable(irho,iu+1,3)
 Tpart1 = m*ui + c

 m = (OPTABLE(irho,iu,4) - OPTABLE(irho,iu+1,4))/(OPTABLE(irho,iu,3)-OPTABLE(irho,iu+1,3))
 c = OPTABLE(irho,iu+1,4) - m*OPTABLE(irho,iu+1,3)
 gmw1 = m*ui + c

 ! Search for ui in irho+1 list for interpolation
 iu = search_table(optable(irho+1,:,3),ny,ui)

 m = (OPTABLE(irho+1,iu,5) - OPTABLE(irho+1,iu+1,5))/(OPTABLE(irho+1,iu,3) - OPTABLE(irho+1,iu+1,3))
 c = OPTABLE(irho+1,iu+1,5) - m*OPTABLE(irho+1,iu+1,3)
 kbar2 = m*ui + c

 m = (OPTABLE(irho+1,iu,6) - OPTABLE(irho+1,iu+1,6))/(OPTABLE(irho+1,iu,3) - OPTABLE(irho+1,iu+1,3))
 c = OPTABLE(irho+1,iu+1,6) - m*OPTABLE(irho+1,iu+1,3)
 kappa2 = m*ui + c

 m = (OPTABLE(irho+1,iu,2) - OPTABLE(irho+1,iu+1,2))/(OPTABLE(irho+1,iu,3) - OPTABLE(irho+1,iu+1,3))
 c = OPTABLE(irho+1,iu+1,2) - m*OPTABLE(irho+1,iu+1,3)
 Tpart2 = m*ui + c
 
 m = (OPTABLE(irho+1,iu,4) - OPTABLE(irho+1,iu+1,4))/(OPTABLE(irho+1,iu,3) - OPTABLE(irho+1,iu+1,3))
 c = OPTABLE(irho+1,iu+1,4) - m*OPTABLE(irho+1,iu+1,3)
 gmw2 = m*ui + c

 m = (kappa2 - kappa1)/(OPTABLE(irho+1,1,1)-OPTABLE(irho,1,1))
 c = kappa2 - m*OPTABLE(irho+1,1,1)
 kappaPart = m*rhoi + c

 m = (kbar2 - kbar1)/(OPTABLE(irho+1,1,1)-OPTABLE(irho,1,1))
 c = kbar2 - m*OPTABLE(irho+1,1,1)
 kappaBar = m*rhoi + c

 m = (Tpart2 - Tpart1)/(OPTABLE(irho+1,1,1)-OPTABLE(irho,1,1))
 c = Tpart2 - m*OPTABLE(irho+1,1,1)
 Ti = m*rhoi + c
 
 m = (gmw2 - gmw1)/(OPTABLE(irho+1,1,1)-OPTABLE(irho,1,1))
 c = gmw2 - m*OPTABLE(irho+1,1,1)
 gmwi = m*rhoi + c
 
end subroutine getopac_opdep

subroutine getintenerg_opdep(Teqi, rhoi, ueqi)
 use io, only:warning
 real, intent(out) :: ueqi
 real, intent(in)    :: Teqi,rhoi

 real :: u1, u2
 real :: m, c
 integer :: irho,itemp

 if (rhoi > OPTABLE(nx,1,1) .or. rhoi < OPTABLE(1,1,1)) then
    call warning('getintenerg_opdep','rhoi out of range',var='rhoi',val=rhoi)
 elseif (Teqi > OPTABLE(1,ny,2) .or. Teqi < OPTABLE(1,1,2)) then
    call warning('getintenerg_opdep','Ti out of range',var='Ti',val=Teqi)
 endif


 ! interpolate through OPTABLE to obtain equilibrium internal energy
 irho = search_table(optable(:,1,1),nx,rhoi)
 itemp = search_table(optable(irho,:,2),ny,Teqi)
 
 m = (OPTABLE(irho,itemp,3) - OPTABLE(irho,itemp+1,3))/(OPTABLE(irho,itemp,2) - OPTABLE(irho,itemp+1,2))
 c = OPTABLE(irho,itemp+1,3) - m*OPTABLE(irho,itemp+1,2)
 u1 = m*Teqi + c

 itemp = search_table(optable(irho+1,:,2),ny,Teqi)
 
 m = (OPTABLE(irho+1,itemp,3) - OPTABLE(irho+1,itemp+1,3))/&
      (OPTABLE(irho+1,itemp,2) - OPTABLE(irho+1,itemp+1,2))
 c = OPTABLE(irho+1,itemp+1,3) - m*OPTABLE(irho+1,itemp+1,2)
 u2 = m*Teqi + c

 m = (u2 - u1)/(OPTABLE(irho+1,1,1)-OPTABLE(irho,1,1))
 c = u2 - m*OPTABLE(irho+1,1,1)

 ueqi = m*rhoi + c
end subroutine getintenerg_opdep

!
! Binary search given array
!
integer function search_table(array,arrlen,invalue) result(outind)
  real,intent(in)    :: array(:),invalue
  integer,intent(in) :: arrlen
  integer            :: leftind,rightind,midind

  leftind = 1; rightind = arrlen
  do
     if (rightind - leftind == 1 .or. rightind == leftind) then
        outind = leftind
        return
     endif
     midind = floor((rightind - leftind) / 2.) + leftind
     if (invalue == array(midind) ) then
        outind = midind
        return
     endif
     if (invalue < array(midind)) then
        rightind = midind
     else
        leftind = midind
     endif
  enddo

end function search_table

!
! Calculate factor for FLD
!
subroutine get_k_fld(rhoi,eni,i,ki,Ti)
 use physcon,  only:c,fourpi
 use units,    only:unit_density,unit_ergg,unit_opacity,get_radconst_code
 real,intent(in)    :: rhoi,eni
 integer,intent(in) :: i
 real               :: kappaBar,gmwi,kappaPart,eni_ergg,rhoi_g
 real,intent(out)   :: ki,Ti

 if (lambda_FLD(i) == 0d0) then
    ki = 0.
 else
    eni_ergg = eni*unit_ergg
    rhoi_g = rhoi*unit_density
    call getopac_opdep(eni_ergg,rhoi_g,kappaBar,kappaPart,Ti,gmwi)
    kappaPart = kappaPart/unit_opacity
    ! steboltz constant = 4pi/c * arad
    ki = 16d0*(fourpi/c)*get_radconst_code()*lambda_FLD(i)*Ti**3 /rhoi/kappaPart
    if (isnan(ki)) then
       print *, "WARNING k isnan, lambda_FLDi,Ti,rhoi,kappaPart", &
             lambda_FLD(i), Ti, rhoi,kappaPart
    endif
 endif
end subroutine get_k_fld

end module eos_stamatellos


