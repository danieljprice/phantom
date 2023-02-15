!--------------------------------------------------------------------------
! Script to interpolate the opacity table myeos.dat
!--------------------------------------------------------------------------

module eos_stamatellos
 implicit none
 real, public :: optable(260,1001,6)
 public :: read_optab,getopac_opdep
contains

subroutine read_optab(ierr)
 use datafiles, only:find_phantom_datafile

 integer, intent(out) :: ierr
 integer i,j,nx,ny
 character(len=120) :: filepath

 ! read in data file for interpolation
 filepath=find_phantom_datafile('myeos.dat','cooling')
 print *,"FILEPATH:",filepath 
 open(10, file=filepath, form="formatted", status="old",iostat=ierr)
 if (ierr > 0) return
 read(10, *) nx, ny
 do i = 1,nx
  do j = 1,ny
   read(10,*) OPTABLE(i,j,1),OPTABLE(i,j,2),OPTABLE(i,j,3),&
              OPTABLE(i,j,4),OPTABLE(i,j,5),OPTABLE(i,j,6)
  enddo
 enddo

end subroutine read_optab

subroutine getopac_opdep(ui,rhoi,kappaBar,kappaPart,Ti,gmwi,gammai)
 real, intent(in)  :: ui,rhoi
 real, intent(out) :: kappaBar,kappaPart,Ti,gmwi,gammai

 integer i,j,nx,ny
 real m,c
 real kbar1,kbar2
 real kappa1,kappa2
 real Tpart1,Tpart2
 real gmw1,gmw2
 real cv
 real ui_, rhoi_

 ! interpolate through OPTABLE to find corresponding kappaBar, kappaPart and T

 if (rhoi.lt.1.0e-24) then
  rhoi_ = 1.0e-24
 else
  rhoi_ = rhoi
 endif  

 i = 1
 do while((OPTABLE(i,1,1).le.rhoi_).and.(i.lt.260))
  i = i + 1
 enddo

 if (ui.lt.0.5302E8) then
  ui_ = 0.5302E8 
 else
  ui_ = ui
 endif

 j = 1
 do while ((OPTABLE(i-1,j,3).le.ui_).and.(j.lt.1000))
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
 
 j = 1
 do while ((OPTABLE(i,j,3).le.ui).and.(j.lt.1000))
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
 !kappaPart = kappaPart*kappa_corr

 m = (kbar2 - kbar1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = kbar2 - m*OPTABLE(i,1,1)

 kappaBar = m*rhoi_ + c
 !kappaBar = kappaBar*kappa_corr

 m = (Tpart2 - Tpart1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = Tpart2 - m*OPTABLE(i,1,1)
 
 Ti = m*rhoi_ + c

 m = (gmw2 - gmw1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = gmw2 - m*OPTABLE(i,1,1)

 gmwi = m*rhoi_ + c

 cv = ui_/Ti
 gammai = 1.0d0 + 1.38d-16/1.67d-24/gmwi/cv
end subroutine getopac_opdep

subroutine getintenerg_opdep(Teqi, rhoi, ueqi, OPTABLE)
 real, intent(out) :: ueqi
 real, intent(in)    :: Teqi,rhoi
 !real, intent(inout)  :: OPTABLE(260,1001,6)

 real u1, u2 
 real m, c
 integer i, j
 real rhoi_
 real, intent(in)  :: OPTABLE(260,1001,6)

 ! interpolate through OPTABLE to obtain equilibrium internal energy

 if (rhoi.lt.1.0e-24) then
  rhoi_ = 1.0e-24 
 else
  rhoi_ = rhoi
 endif

 i = 1
 do while((OPTABLE(i,1,1).le.rhoi_).and.(i.lt.260))
  i = i + 1
 enddo

 j = 1
 do while ((OPTABLE(i-1,j,2).le.Teqi).and.(j.lt.1000))
  j = j + 1
 enddo

 m = (OPTABLE(i-1,j-1,3) - OPTABLE(i-1,j,3))/(OPTABLE(i-1,j-1,2) - OPTABLE(i-1,j,2))
 c = OPTABLE(i-1,j,3) - m*OPTABLE(i-1,j,2)

 u1 = m*Teqi + c

 j = 1
 do while ((OPTABLE(i,j,2).le.Teqi).and.(j.lt.1000))
  j = j + 1
 enddo

 m = (OPTABLE(i,j-1,3) - OPTABLE(i,j,3))/(OPTABLE(i,j-1,2) - OPTABLE(i,j,2))
 c = OPTABLE(i,j,3) - m*OPTABLE(i,j,2)

 u2 = m*Teqi + c

 m = (u2 - u1)/(OPTABLE(i,1,1)-OPTABLE(i-1,1,1))
 c = u2 - m*OPTABLE(i,1,1)

 ueqi = m*rhoi_ + c
end subroutine getintenerg_opdep

end module eos_stamatellos


