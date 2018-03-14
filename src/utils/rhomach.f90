!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: rhomach
!
!  DESCRIPTION:
!  Analysis routine to compute density variances and Mach number
!  (Used for Price, Federrath and Brunt 2010)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module rhomach
 implicit none
 integer, public :: lunit_rhomach = 87
 public :: get_rhomach_grid,write_rhomach

contains

!------------------------------------------------------------------------
!+
!  Routine to calculate the linear and logarithmic density variance
!  as well as the RMS velocity (=RMS Mach number if isothermal and cs=1)
!  from gridded data on a uniform mesh
!
!  Input: datgrid, a (4,nx,ny,nz) dimensional array
!                  containing density,vx,vy,vz
!
!  Also returns a one dimensional array of the density field (rhogrid)
!+
!------------------------------------------------------------------------
subroutine get_rhomach_grid(datgrid,nx,ny,nz,smin,rhomean,rhovar,smean,svar,rmsv,rhogrid)
 real,    intent(in)  :: datgrid(:,:,:,:)
 integer, intent(in)  :: nx,ny,nz
 real,    intent(in)  :: smin
 real,    intent(out) :: rhomean,rhovar
 real,    intent(out) :: smean,svar,rmsv
 real,    intent(out), allocatable :: rhogrid(:)
 integer :: i,j,k,n,isize
 real    :: rhoi,si,vx2,vy2,vz2
 real    :: rhomin,rhomax,totvol,totmass,dn
 logical :: densityonly

 if (nx*ny*nz <= 0) then
    print*,'get_rhomach_grid: ERROR: nx,ny,nz = ',nx,ny,nz
    return
 endif
 isize = size(datgrid(:,1,1,1))
 densityonly = .false.
 if (isize==1) then
    print*,'get_rhomach_grid: WARNING: datgrid has only one array, assuming density only'
    densityonly = .true.
 elseif (isize < 4) then
    print*,'get_rhomach_grid: ERROR: datgrid does not contain enough data (rho,vx,vy,vz)'
    return
 endif
!
!--calculate rhomach quantities from the gridded data
!
 rhomin = huge(rhomin)
 rhomax  = 0.
 rhomean = 0.
 smean   = 0.
 rmsv    = 0.
 do k=1,nz
    do j=1,ny
       do i=1,nx
          rhoi = datgrid(1,i,j,k)
          if (rhoi > 0.) then
             si = log(rhoi)
          else
             si = smin
          endif
          if (.not.densityonly) then
             vx2 = datgrid(2,i,j,k)**2
             vy2 = datgrid(3,i,j,k)**2
             vz2 = datgrid(4,i,j,k)**2
             rmsv    = rmsv + (vx2 + vy2 + vz2)
          endif
          rhomin  = min(rhomin,rhoi)
          rhomax  = max(rhomax,rhoi)
          rhomean = rhomean + rhoi
          smean   = smean + si
       enddo
    enddo
 enddo

 dn = 1./real(nx*ny*nz)
 rmsv    = sqrt(rmsv*dn)
 rhomean = rhomean*dn
 smean   = smean*dn
 totvol  = 1.
 totmass = totvol*rhomean
 print*,'On grid: Total Mass = ',totmass,' mean dens = ',rhomean
 print*,'On grid: max. dens = ',rhomax,' min. dens = ',rhomin
 print*,'On grid: rms v = ',rmsv,' smean = ',smean
 !
 !--calculate variances from second pass
 !
 if (.not.allocated(rhogrid)) allocate(rhogrid(nx*ny*nz))
 n = 0
 svar   = 0.
 rhovar = 0.
 do k=1,nz
    do j=1,ny
       do i=1,nx
          n = n + 1
          rhoi     = datgrid(1,i,j,k)
          if (rhoi > 0.) then
             si = log(rhoi)
          else
             si = smin
          endif
          rhogrid(n) = si
          rhovar = rhovar + (rhoi - rhomean)**2
          svar   = svar   + (si - smean)**2
       enddo
    enddo
 enddo

 rhovar = rhovar*dn
 svar   = svar*dn
 print*,' var(rho) = ',rhovar,' var(lnrho) = ',svar

end subroutine get_rhomach_grid

!----------------------------------------------------------------
!+
!  Routine to write the volume and mass weighted variances,
!  standard deviations and RMS Mach numbers to a file
!  (I usually call the file "rhomach.out")
!+
!----------------------------------------------------------------
subroutine write_rhomach(filename,time,rhomeanvw,rhomeanmw,rhovarvw,rhovarmw,&
                         rmsv,rmsvmw,smeanvw,smeanmw,svarvw,svarmw)
 character(len=*), intent(in) :: filename
 real,             intent(in) :: time,rhomeanvw,rhomeanmw,rhovarvw,rhovarmw
 real,             intent(in) :: rmsv,rmsvmw,smeanvw,smeanmw,svarvw,svarmw
 character(len=20) :: fmtstring
 real              :: bval,bvalmw
 integer           :: ierr

 print "(a)",' writing rhomach data to '//trim(filename)
 if (rmsv > 0.) then
    bval = sqrt(svarvw)/rmsv
 else
    bval = 0.
 endif
 if (rmsvmw > 0.) then
    bvalmw = sqrt(svarmw)/rmsvmw
 else
    bvalmw = 0.
 endif

 open(unit=lunit_rhomach,file=trim(filename),status='replace',form='formatted')
 write(fmtstring,"('(',i2,'(es18.10,1x))')",iostat=ierr) 17
 write(lunit_rhomach,fmtstring) time,rhomeanvw,rhomeanmw,rhovarvw,rhovarmw,sqrt(rhovarvw),sqrt(rhovarmw),&
                           rmsv,rmsvmw,bval,bvalmw,smeanvw,smeanmw,svarvw,svarmw,sqrt(svarvw),sqrt(svarmw)
 close(unit=lunit_rhomach)

end subroutine write_rhomach

end module rhomach
