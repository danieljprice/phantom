!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Dependencies: datafiles, part
!

 implicit none
 real,allocatable,public :: optable(:,:,:)
 real,allocatable,public :: Gpot_cool(:),duFLD(:),gradP_cool(:),lambda_FLD(:),urad_FLD(:) !gradP_cool=gradP/rho
 real,allocatable,public :: ttherm_store(:),teqi_store(:),opac_store(:),duSPH(:)
 character(len=25), public :: eos_file= 'eos_lom.dat' !default name of tabulated EOS file
 logical,public :: doFLD = .True., floor_energy = .False.
 integer,public :: iunitst=19
 integer,save :: nx,ny ! dimensions of optable read in

 public :: read_optab,getopac_opdep,init_S07cool,getintenerg_opdep,finish_S07cool
 public :: get_k_fld

contains

subroutine init_S07cool()
 use part, only:npart,maxradprop
 print *, "Allocating cooling arrays"
 allocate(gradP_cool(npart))
 allocate(Gpot_cool(npart))
 allocate(duFLD(npart))
 allocate(lambda_fld(npart))
 allocate(urad_FLD(npart))
 allocate(ttherm_store(npart))
 allocate(teqi_store(npart))
 allocate(opac_store(npart))
 allocate(duSPH(npart))
 Gpot_cool(:) = 0d0
 gradP_cool(:) = 0d0
 urad_FLD(:) = 0d0
 duFLD(:) = 0d0
 teqi_store(:) = 0d0
 ttherm_store(:) = 0d0
 opac_store(:) = 0d0
 duSPH(:) = 0d0
 open (unit=iunitst,file='EOSinfo.dat',status='replace')
 if (doFLD) then
    print *, "Using Forgan+ 2009 hybrid cooling method (FLD)"
 else
    print *, "NOT using FLD. Using cooling only"
 endif
end subroutine init_S07cool

subroutine finish_S07cool()
 deallocate(optable)
 if (allocated(gradP_cool)) deallocate(gradP_cool)
 if (allocated(Gpot_cool)) deallocate(Gpot_cool)
 if (allocated(duFLD)) deallocate(duFLD)
 if (allocated(lambda_fld)) deallocate(lambda_fld)
 if (allocated(urad_FLD)) deallocate(urad_FLD)
 if (allocated(ttherm_store)) deallocate(ttherm_store)
 if (allocated(teqi_store)) deallocate(teqi_store)
 if (allocated(opac_store)) deallocate(opac_store)
 if (allocated(duSPH)) deallocate(duSPH)
 close(iunitst)
end subroutine finish_S07cool

subroutine read_optab(eos_file,ierr)
 use datafiles, only:find_phantom_datafile
 character(len=*),intent(in) :: eos_file
 integer, intent(out) :: ierr
 integer i,j,errread
 character(len=120) :: filepath,junk

 ! read in data file for interpolation
 filepath=find_phantom_datafile(eos_file,'eos/lombardi')
 print *,"EOS file: FILEPATH:",filepath
 open(10, file=filepath, form="formatted", status="old",iostat=ierr)
 if (ierr > 0) return
 do
    read(10,'(A120)') junk
    if (len(trim(adjustl(junk))) == 0) cycle ! blank line
    if ((index(adjustl(junk),"::") == 0) .and. (index(adjustl(junk),"#") .ne. 1 )) then !ignore comment lines
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

 rhomin = OPTABLE(1,1,1)
 umin = OPTABLE(1,1,3)
 ! interpolate through OPTABLE to find corresponding kappaBar, kappaPart and T

 ! check values are in range of tables
 if (rhoi > OPTABLE(nx,1,1) .or. rhoi < OPTABLE(1,1,1)) then
    call fatal('getopac_opdep','rhoi out of range. Collapsing clump?',var='rhoi',val=rhoi)
 elseif (ui > OPTABLE(1,ny,3) .or. ui < OPTABLE(1,1,3)) then
    call fatal('getopac_opdep','ui out of range',var='ui',val=ui)
 endif

 if (rhoi <  rhomin) then
    rhoi_ = rhomin
 else
    rhoi_ = rhoi
 endif

 i = 2
 do while((OPTABLE(i,1,1) <= rhoi_).and.(i < nx))
    i = i + 1
 enddo

 if (ui < umin) then
    ui_ = umin
 else
    ui_ = ui
 endif

 j = 2
 do while ((OPTABLE(i-1,j,3) <= ui_).and.(j < ny))
    j = j + 1
 enddo

 m = (OPTABLE(i-1,j-1,5) - OPTABLE(i-1,j,5))/(OPTABLE(i-1,j-1,3) - OPTABLE(i-1,j,3))
 c = OPTABLE(i-1,j,5) - m*OPTABLE(i-1,j,3)

 kbar1 = m*ui_ + c

 m = (OPTABLE(i-1,j-1,6) - OPTABLE(i-1,j,6))/(OPTABLE(i-1,j-1,3) - OPTABLE(i-1,j,3))
 c = OPTABLE(i-1,j,6) - m*OPTABLE(i-1,j,3)

 kappa1 = m*ui_ + c

 m = (OPTABLE(i-1,j-1,2) - OPTABLE(i-1,j,2))/(OPTABLE(i-1,j-1,3) - OPTABLE(i-1,j,3))
 c = OPTABLE(i-1,j,2) - m*OPTABLE(i-1,j,3)

 Tpart1 = m*ui_ + c

 m = (OPTABLE(i-1,j-1,4) - OPTABLE(i-1,j,4))/(OPTABLE(i-1,j-1,3) - OPTABLE(i-1,j,3))
 c = OPTABLE(i-1,j,4) - m*OPTABLE(i-1,j,3)

 gmw1 = m*ui_ + c

 j = 2
 do while ((OPTABLE(i,j,3) <= ui).and.(j < ny))
    j = j + 1
 enddo

 m = (OPTABLE(i,j-1,5) - OPTABLE(i,j,5))/(OPTABLE(i,j-1,3) - OPTABLE(i,j,3))
 c = OPTABLE(i,j,5) - m*OPTABLE(i,j,3)

 kbar2 = m*ui_ + c

 m = (OPTABLE(i,j-1,6) - OPTABLE(i,j,6))/(OPTABLE(i,j-1,3) - OPTABLE(i,j,3))
 c = OPTABLE(i,j,6) - m*OPTABLE(i,j,3)

 kappa2 = m*ui_ + c

 m = (OPTABLE(i,j-1,2) - OPTABLE(i,j,2))/(OPTABLE(i,j-1,3) - OPTABLE(i,j,3))
 c = OPTABLE(i,j,2) - m*OPTABLE(i,j,3)

 Tpart2 = m*ui_ + c

 m = (OPTABLE(i,j-1,4) - OPTABLE(i,j,4))/(OPTABLE(i,j-1,3) - OPTABLE(i,j,3))
 c = OPTABLE(i,j,4) - m*OPTABLE(i,j,3)

 gmw2 = m*ui_ + c

 m = (kappa2 - kappa1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = kappa2 - m*OPTABLE(i,1,1)

 kappaPart = m*rhoi_ + c

 m = (kbar2 - kbar1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = kbar2 - m*OPTABLE(i,1,1)

 kappaBar = m*rhoi_ + c

 m = (Tpart2 - Tpart1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = Tpart2 - m*OPTABLE(i,1,1)

 Ti = m*rhoi_ + c

 m = (gmw2 - gmw1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = gmw2 - m*OPTABLE(i,1,1)

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

 if (rhoi > OPTABLE(nx,1,1) .or. rhoi < OPTABLE(1,1,1)) then
    call warning('getintenerg_opdep','rhoi out of range',var='rhoi',val=rhoi)
 elseif (Teqi > OPTABLE(1,ny,2) .or. Teqi < OPTABLE(1,1,2)) then
    call warning('getintenerg_opdep','Ti out of range',var='Ti',val=Teqi)
 endif


 ! interpolate through OPTABLE to obtain equilibrium internal energy

 if (rhoi < 1.0e-24) then
    rhoi_ = 1.0e-24
 else
    rhoi_ = rhoi
 endif

 i = 2
 do while((OPTABLE(i,1,1) <= rhoi_).and.(i < nx))
    i = i + 1
 enddo

 j = 2
 do while ((OPTABLE(i-1,j,2) <= Teqi).and.(j < ny))
    j = j + 1
 enddo


 m = (OPTABLE(i-1,j-1,3) - OPTABLE(i-1,j,3))/(OPTABLE(i-1,j-1,2) - OPTABLE(i-1,j,2))
 c = OPTABLE(i-1,j,3) - m*OPTABLE(i-1,j,2)

 u1 = m*Teqi + c

 j = 2
 do while ((OPTABLE(i,j,2) <= Teqi).and.(j < ny))
    j = j + 1
 enddo

 m = (OPTABLE(i,j-1,3) - OPTABLE(i,j,3))/(OPTABLE(i,j-1,2) - OPTABLE(i,j,2))
 c = OPTABLE(i,j,3) - m*OPTABLE(i,j,2)

 u2 = m*Teqi + c

 m = (u2 - u1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = u2 - m*OPTABLE(i,1,1)

 ueqi = m*rhoi_ + c
end subroutine getintenerg_opdep

subroutine get_k_fld(rhoi,eni,i,ki,Ti)
 use physcon,  only:c,fourpi
 use units,    only:unit_density,unit_ergg,unit_opacity,get_radconst_code
 real,intent(in)    :: rhoi,eni
 integer,intent(in) :: i
 real               :: kappaBar,gmwi,kappaPart
 real,intent(out)   :: ki,Ti

 if (lambda_FLD(i) == 0d0) then
    ki = 0.
 else
    call getopac_opdep(eni*unit_ergg,rhoi*unit_density,kappaBar,kappaPart,Ti,gmwi)
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


