!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Dependencies: allocutils, datafiles, dim, io
!

 implicit none
 real,allocatable,public :: optable(:,:,:)
 real,allocatable,public :: gradP_cool(:)!gradP_cool=gradP/rho
 real,allocatable,public :: ttherm_store(:),ueqi_store(:),opac_store(:),duSPH(:)
 character(len=25), public :: eos_file= 'eos_lom.dat' !default name of tabulated EOS file
 logical,public :: floor_energy = .False.
 integer,public :: iunitst=19
 integer,save :: nx,ny ! dimensions of optable read in

 public :: read_optab,getopac_opdep,init_coolra,getintenerg_opdep,finish_coolra

contains


subroutine init_coolra()
 use dim, only:maxp
 use allocutils, only:allocate_array

 print *, "Allocating cooling arrays for maxp=",maxp
 call allocate_array('gradP_cool',gradP_cool,maxp)
 call allocate_array('ttherm_store',ttherm_store,maxp)
 call allocate_array('ueqi_store',ueqi_store,maxp)
 call allocate_array('opac_store',opac_store,maxp)
 call allocate_array('duSPH',duSPH,maxp)

 gradP_cool(:) = 0.
 ueqi_store(:) = 0.
 ttherm_store(:) = 0.
 opac_store(:) = 0.
 duSPH(:) = 0.

 print *, "NOT using FLD. Using cooling only"

end subroutine init_coolra

subroutine finish_coolra()

 if (allocated(optable)) deallocate(optable)
 if (allocated(gradP_cool)) deallocate(gradP_cool)
 if (allocated(ttherm_store)) deallocate(ttherm_store)
 if (allocated(ueqi_store)) deallocate(ueqi_store)
 if (allocated(opac_store)) deallocate(opac_store)
 if (allocated(duSPH)) deallocate(duSPH)

end subroutine finish_coolra




subroutine read_optab(eos_file,ierr)
 use datafiles, only:find_phantom_datafile
 character(len=*),intent(in) :: eos_file
 integer, intent(out) :: ierr
 integer i,j,errread
 character(len=120) :: filepath,junk

 ! read in EOS and opacity data file for interpolation
 filepath=find_phantom_datafile(eos_file,'eos/lombardi')
 print *,"EOS file: FILEPATH:",filepath
 open(10,file=filepath,form="formatted",status="old",iostat=ierr)
 if (ierr > 0) return
 do
    read(10,'(A120)') junk
    print *, junk
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
       read(10,*) optable(i,j,1),optable(i,j,2),optable(i,j,3),&
              optable(i,j,4),optable(i,j,5),optable(i,j,6)
    enddo
 enddo
end subroutine read_optab

!
! Main subroutine for interpolating tables to get EOS values
!
subroutine getopac_opdep(ui,rhoi,kappaBar,kappaPart,Ti,gmwi)
 use io, only:fatal
 real, intent(in)  :: ui,rhoi
 real, intent(out) :: kappaBar,kappaPart,Ti,gmwi

 integer i,j
 real m,c
 real kbar1,kbar2
 real kappa1,kappa2
 real Tpart1,Tpart2
 real gmw1,gmw2
 real ui_, rhoi_,rhomin,umin

 rhomin = optable(1,1,1)
 umin = optable(1,1,3)
 ! interpolate through optable to find corresponding kappaBar, kappaPart and T

 ! check values are in range of tables
 if (rhoi > optable(nx,1,1) .or. rhoi < optable(1,1,1)) then
    print *, "optable rho min =", rhomin
    call fatal('getopac_opdep','rhoi out of range. Collapsing clump?',var='rhoi',val=rhoi)
 elseif (ui > optable(1,ny,3) .or. ui < optable(1,1,3)) then
    call fatal('getopac_opdep','ui out of range',var='ui',val=ui)
 endif

 if (rhoi <  rhomin) then
    rhoi_ = rhomin
 else
    rhoi_ = rhoi
 endif

 i = 2
 do while((optable(i,1,1) <= rhoi_).and.(i < nx))
    i = i + 1
 enddo

 if (ui < umin) then
    ui_ = umin
 else
    ui_ = ui
 endif

 j = 2
 do while ((optable(i-1,j,3) <= ui_).and.(j < ny))
    j = j + 1
 enddo

 m = (optable(i-1,j-1,5) - optable(i-1,j,5))/(optable(i-1,j-1,3) - optable(i-1,j,3))
 c = optable(i-1,j,5) - m*optable(i-1,j,3)

 kbar1 = m*ui_ + c

 m = (optable(i-1,j-1,6) - optable(i-1,j,6))/(optable(i-1,j-1,3) - optable(i-1,j,3))
 c = optable(i-1,j,6) - m*optable(i-1,j,3)

 kappa1 = m*ui_ + c

 m = (optable(i-1,j-1,2) - optable(i-1,j,2))/(optable(i-1,j-1,3) - optable(i-1,j,3))
 c = optable(i-1,j,2) - m*optable(i-1,j,3)

 Tpart1 = m*ui_ + c

 m = (optable(i-1,j-1,4) - optable(i-1,j,4))/(optable(i-1,j-1,3) - optable(i-1,j,3))
 c = optable(i-1,j,4) - m*optable(i-1,j,3)

 gmw1 = m*ui_ + c

 j = 2
 do while ((optable(i,j,3) <= ui).and.(j < ny))
    j = j + 1
 enddo

 m = (optable(i,j-1,5) - optable(i,j,5))/(optable(i,j-1,3) - optable(i,j,3))
 c = optable(i,j,5) - m*optable(i,j,3)

 kbar2 = m*ui_ + c

 m = (optable(i,j-1,6) - optable(i,j,6))/(optable(i,j-1,3) - optable(i,j,3))
 c = optable(i,j,6) - m*optable(i,j,3)

 kappa2 = m*ui_ + c

 m = (optable(i,j-1,2) - optable(i,j,2))/(optable(i,j-1,3) - optable(i,j,3))
 c = optable(i,j,2) - m*optable(i,j,3)

 Tpart2 = m*ui_ + c

 m = (optable(i,j-1,4) - optable(i,j,4))/(optable(i,j-1,3) - optable(i,j,3))
 c = optable(i,j,4) - m*optable(i,j,3)

 gmw2 = m*ui_ + c

 m = (kappa2 - kappa1)/(optable(i,1,1)-optable(i-1,1,1))
 c = kappa2 - m*optable(i,1,1)

 kappaPart = m*rhoi_ + c

 m = (kbar2 - kbar1)/(optable(i,1,1)-optable(i-1,1,1))
 c = kbar2 - m*optable(i,1,1)

 kappaBar = m*rhoi_ + c

 m = (Tpart2 - Tpart1)/(optable(i,1,1)-optable(i-1,1,1))
 c = Tpart2 - m*optable(i,1,1)

 Ti = m*rhoi_ + c

 m = (gmw2 - gmw1)/(optable(i,1,1)-optable(i-1,1,1))
 c = gmw2 - m*optable(i,1,1)

 gmwi = m*rhoi_ + c
end subroutine getopac_opdep

subroutine getintenerg_opdep(Teqi, rhoi, ueqi)
 use io, only:warning
 real, intent(out) :: ueqi
 real, intent(in)    :: Teqi,rhoi

 real u1, u2
 real m, c
 integer i, j
 real rhoi_

 if (rhoi > optable(nx,1,1) .or. rhoi < optable(1,1,1)) then
    call warning('getintenerg_opdep','rhoi out of range',var='rhoi',val=rhoi)
 elseif (Teqi > optable(1,ny,2) .or. Teqi < optable(1,1,2)) then
    call warning('getintenerg_opdep','Ti out of range',var='Ti',val=Teqi)
 endif


 ! interpolate through optable to obtain equilibrium internal energy

 if (rhoi < 1.0e-24) then
    rhoi_ = 1.0e-24
 else
    rhoi_ = rhoi
 endif

 i = 2
 do while((optable(i,1,1) <= rhoi_).and.(i < nx))
    i = i + 1
 enddo

 j = 2
 do while ((optable(i-1,j,2) <= Teqi).and.(j < ny))
    j = j + 1
 enddo


 m = (optable(i-1,j-1,3) - optable(i-1,j,3))/(optable(i-1,j-1,2) - optable(i-1,j,2))
 c = optable(i-1,j,3) - m*optable(i-1,j,2)

 u1 = m*Teqi + c

 j = 2
 do while ((optable(i,j,2) <= Teqi).and.(j < ny))
    j = j + 1
 enddo

 m = (optable(i,j-1,3) - optable(i,j,3))/(optable(i,j-1,2) - optable(i,j,2))
 c = optable(i,j,3) - m*optable(i,j,2)

 u2 = m*Teqi + c

 m = (u2 - u1)/(optable(i,1,1)-optable(i-1,1,1))
 c = u2 - m*optable(i,1,1)

 ueqi = m*rhoi_ + c
end subroutine getintenerg_opdep

end module eos_stamatellos


