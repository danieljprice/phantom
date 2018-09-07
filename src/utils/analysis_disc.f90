!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for discs
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: discanalysisutils, infile_utils, io, part, physcon
!+
!--------------------------------------------------------------------------
module analysis
 use discanalysisutils, only:disc_analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'CJN'
 public :: do_analysis

 integer, parameter :: nr = 300
 real,dimension(nr) :: twist,twistprev

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use part,    only:xyzmh_ptmass,vxyz_ptmass,nptmass
 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=9) :: output
 character(len=20) :: filename
 character(len=20) :: discprefix
 integer :: i,ierr,iline
 real :: R_in,R_out,H_R,p_index,q_index,M_star
 real :: G,rmin,rmax
 real :: tilt(nr),Lx(nr),Ly(nr),Lz(nr)
 real :: rad(nr),h_smooth(nr),sigma(nr),H(nr)
 real :: unitlx(nr),unitly(nr),unitlz(nr),ecc(nr)
 real :: psi(nr),tilt_acc(nr)
 integer :: ninbin(nr)
 logical :: assume_Ltot_is_same_as_zaxis,iexist

 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 logical :: do_precession,ifile

 do_precession = .false.

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'angm',numfile
 write(*,'("Output file name is ",A)') output

! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G==1")')
 G = 1.0

 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 inquire(file=trim(discprefix)//'.discparams', exist=ifile)
 if (ifile) then
    call read_discparams(trim(discprefix)//'.discparams',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read .discparams file')
 else
    call read_discparams('discparams.list',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')
 endif

! Print out the parameters
 write(*,*)
 write(*,'("Parameters are:")')
 write(*,*) 'R_in    = ',R_in
 write(*,*) 'R_out   = ',R_out
 write(*,*) 'H/R_ref = ',H_R
 write(*,*) 'p_index = ',p_index
 write(*,*) 'q_index = ',q_index
 write(*,*) 'M_star  = ',M_star
 write(*,*)
 write(*,*)

! Setup rmin and rmax for the analysis
 rmin = R_in
 rmax = R_out

! This variable should be set to false for any discs that use sink particles to set
! the potential or any discs that have a warp
! For any setup that uses iexternalforce and assumes that the vast majority of the angular
! momentum is held by the central potential, this should be set to true

 assume_Ltot_is_same_as_zaxis = .false.

 call disc_analysis(xyzh,vxyz,npart,pmass,time,nr,rmin,rmax,H_R,G,M_star,q_index,&
                     tilt,tilt_acc,twist,twistprev,psi,H,rad,h_smooth,sigma,unitlx,unitly,unitlz,&
                     Lx,Ly,Lz,ecc,ninbin,assume_Ltot_is_same_as_zaxis,xyzmh_ptmass,vxyz_ptmass,nptmass)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'radius', &
       2,'sigma', &
       3,'<h>/H', &
       4,'lx', &
       5,'ly', &
       6,'lz', &
       7,'tilt', &
       8,'twist', &
       9,'psi', &
       10,'H/R', &
       11,'|e|'

 do i = 1,nr
    if (ninbin(i) > 0) then
       write(iunit,'(13(es18.10,1X))') rad(i),sigma(i),h_smooth(i),unitlx(i),unitly(i),unitlz(i),&
                                         tilt(i),twistprev(i),psi(i),H(i)/rad(i),ecc(i)
    endif



! Printing time and twist for each radius bin
    if (do_precession) then
       write(filename,"(a,i3.3)")"precess",i
       inquire(file=filename,exist=iexist)
       if (.not.iexist .or. numfile==0) then
          open(unit=iprec,file=filename,status="replace")
          write(iprec,'("# tilt and twist with time for r = ",es18.10)') rad(i)
          write(iprec,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'rad', &
               2,'time', &
               3,'tilt', &
               4,'twist', &
               5,'tot twist', &
               6,'|e|'
       else
          open(unit=iprec,file=filename,status="old",position="append")
       endif
       write(iprec,'(6(es18.10,1X))') rad(i),time,tilt(i),twist(i),twistprev(i),ecc(i)
       close(unit=iprec)
    endif

 enddo

 close(iunit)

end subroutine do_analysis

!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------
subroutine read_discparams(filename,R_in,R_out,H_R,p_index,q_index,M_star,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,H_R,p_index,q_index,M_star
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)

! Read in parameters from the file discparams.list
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(R_in,'R_in',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_out,'R_out',db,ierr)
 if (ierr /= 0) return
 call read_inopt(H_R,'H/R_ref',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_discparams

end module analysis
