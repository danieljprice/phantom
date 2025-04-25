!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! This will calculate the total masses above given threshholds
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, part, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'ConvergStats'
 integer, parameter :: nden = 4
 real               :: dthresh_cgs(nden),dthresh(nden)
 logical, private   :: firstcall = .true.

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,          only: maxvxyzu
 use part,         only: nptmass,xyzmh_ptmass,rhoh,isdead_or_accreted
 use units,        only: unit_density
 use eos,          only: ieos,get_spsound
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i,j,iparts(nden)
 real                         :: rhoi,v2i,csi,rmsmachi,masses(nden+1),rmsmach(nden),vxyzui(maxvxyzu)
 character(len=200)           :: fileout

 !
 ! Initialise values & Open file
 !
 write(fileout,'(2a)') trim(dumpfile(1:index(dumpfile,'_')-1)),'_MassEvolution.dat'
 if ( firstcall ) then
    firstcall = .false.
    dthresh_cgs(1) = 1.0d-23
    dthresh_cgs(2) = 1.0d-22
    dthresh_cgs(3) = 1.0d-21
    dthresh_cgs(4) = 1.0d-20
    dthresh = dthresh_cgs/unit_density

    ! open files
    open(iunit,file=fileout,status='replace')
    write(iunit,  "('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'idump',   &
          2,'time',    &
          3,'mass > rho1', &
          4,'mass > rho2', &
          5,'mass > rho3', &
          6,'mass > rho4', &
          7,'msink'      , &
          8,'rmsmch>rho1', &
          9,'rmsmch>rho2', &
         10,'rmsmch>rho3', &
         11,'rmsmch>rho4'
 else
    open(iunit,file=fileout,position='append')
 endif

 masses  = 0.
 rmsmach = 0.
 iparts  = 0
 ! get gas masses
 print*, 'printing gas distributions for dthresh=',dthresh_cgs,dthresh
!$omp parallel do default(none) &
!$omp shared(ieos,npart,xyzh,vxyzu,particlemass,dthresh) &
!$omp private(i,j,rhoi,v2i,csi,vxyzui,rmsmachi) &
!$omp reduction(+:masses,rmsmach,iparts)
 do i = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       vxyzui = vxyzu(:,i)
       rhoi = rhoh(xyzh(4,i),particlemass)
       v2i  = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       csi  = get_spsound(ieos,xyzh(1:3,i),rhoi,vxyzui)
       rmsmachi = v2i/csi**2
       do j = 1,4
          if (rhoi > dthresh(j)) then
             masses(j)  = masses(j)  + particlemass
             rmsmach(j) = rmsmach(j) + rmsmachi
             iparts(j)  = iparts(j)  + 1
          endif
       enddo
    endif
 enddo
!omp end parallel do
 rmsmach = sqrt(rmsmach/iparts)

 ! get total sink mass
 print*, 'printing sink stats'
 do i = 1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) masses(5) = masses(5) + xyzmh_ptmass(4,i)
 enddo

 print*, 'printing to file'
 write(iunit,'((I18,1x),10(es18.10,1x))') num,time,masses,rmsmach
 close(iunit)

end subroutine do_analysis

end module analysis
