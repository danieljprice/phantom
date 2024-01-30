!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Compute gravitational wave strain in post-processing
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: deriv, externalforces, gravwaveutils, infile_utils,
!   initial, part, prompting, readwrite_infile, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'gws'

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunitone)
 use externalforces,   only:initialise_externalforces,update_externalforce,externalforce,externalforce_vdependent
 use deriv,            only:derivs
 use initial,          only:initialise,startrun,endrun
 use readwrite_infile, only:read_infile
 use prompting,        only:prompt
 use units,            only:utime,umass,udist,set_units
 use part,             only:isdead_or_accreted,fxyzu,fext,nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass
 use gravwaveutils,    only:calculate_strain
 use infile_utils,     only:open_db_from_file,inopts,read_inopt,close_db
 integer           :: iunit
 type(inopts), allocatable :: db(:)
 character(len=*),     intent(in)        :: dumpfile
 real,                 intent(inout)     :: xyzh(:,:),vxyzu(:,:)
 real,                 intent(inout)     :: pmass,time
 integer,              intent(in)        :: iunitone,numfile
 integer,              intent(inout)     :: npart
 character(len=120)    :: infile,logfile,evfile,dfile
 integer               :: ierr,i
 real                  :: dtextforce,x0(3),v0(3),a0(3),q(6),r2,x,y,z
 integer, parameter    :: iu = 1993, iuu=1994, iuuu=1995
 real,parameter        :: onethird=1/3.
 logical, save         :: first = .true., firstdump=.true.,firstdumpa=.true.
 real                  :: hp(4),hx(4),ddq_xy(3,3)
 real                  :: utime_tmp,umass_tmp,udist_tmp,dtmax

!
!--store units, otherwise initialise() put them to 1
!
 utime_tmp = utime
 umass_tmp = umass
 udist_tmp = udist
!
!--initialize COM position
!
 x0(:)=0.
 v0(:)=0.
 a0(:)=0.

! so if there are very small numbers the code does not complain
 dtextforce=huge(dtextforce)
 print*,'pmass',pmass
!
!--read dtmax from the infile
!
 infile = dumpfile(1:index(dumpfile,'_')-1)//'.in'
 call initialise()
 call set_units(udist_tmp,umass_tmp,utime_tmp)
 call read_infile(infile,logfile,evfile,dfile)
 call open_db_from_file(db,infile,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(dtmax,'dtmax',db,ierr)
 print*, "dtmax read from infile", dtmax
 call startrun(infile,logfile,evfile,dfile,noread=.true.)
 call calculate_strain(hx,hp,pmass,ddq_xy,x0,v0,a0,npart,xyzh,vxyzu(1:3,:),fxyzu(1:3,:),fext(1:3,:),&
                       nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)
 call endrun()

! Write a file where I append all the values of the strain wrt time
 if (first) then
    first = .false.
    open(unit=iu, file='strain.gw',status='replace')
    write(iu,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
          1, 'time',  &
          2, 'hx_0',  &
          3, 'hp_0',  &
          4, 'hx_{30}', &
          5, 'hp_{30}', &
          6, 'hx_{60}', &
          7, 'hx_{60}', &
          8, 'hx_{90}', &
          9, 'hx_{90}'
 else
    open(unit=iu, file='strain.gw',position='append')
 endif

 print*, 'time', time
 write(iu,'(9(es18.10,1X))') time,hx(1),hp(1),hx(2),hp(2),hx(3),hp(3),hx(4),hp(4)

 q(:)=0
 do i=1,npart
    if (xyzh(4,i) > tiny(xyzh)) then  !if not accreted
       x   = xyzh(1,i) - x0(1)
       y   = xyzh(2,i) - x0(2)
       z   = xyzh(3,i) - x0(3)
       r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
! calculate the components of the traceless quadrupole
       q(1) = q(1) + pmass*(x*x-onethird*r2) !qxx
       q(2) = q(2) + pmass*(x*y-onethird*r2) !qxy
       q(3) = q(3) + pmass*(x*z-onethird*r2) !qxz
       q(4) = q(4) + pmass*(y*y-onethird*r2) !qyy
       q(5) = q(5) + pmass*(y*z-onethird*r2) !qyz
       q(6) = q(6) + pmass*(z*z-onethird*r2) !qzz
    endif
 enddo

! Write a file where I append all the values of the strain wrt time
 if (firstdump) then
    firstdump = .false.
    open(unit=iuu, file='quadrupole.txt',status='replace')
    write(iuu,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
          1, 'q11', &
          2, 'q12', &
          3, 'q13', &
          4, 'q22', &
          5, 'q23', &
          6, 'q33'
 else
    open(unit=iuu, file='quadrupole.txt',position='append')
 endif

 write(iuu,'(6(es18.10,1X))') q(1), q(2), q(3), q(4), q(5), q(6)

 if (firstdumpa) then
    firstdumpa = .false.
    open(unit=iuuu, file='second_time_quadrupole.txt',status='replace')
    write(iuuu,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
          1, 'time',&
          2, 'ddq_xy(1,1)', &
          3, 'ddq_xy(1,2)', &
          4, 'ddq_xy(1,3)', &
          5, 'ddq_xy(2,2)', &
          6, 'ddq_xy(2,3)', &
          7, 'ddq_xy(3,3)'
 else
    open(unit=iuuu, file='second_time_quadrupole.txt',position='append')
 endif

 write(iuuu,'(7(es18.10,1X))') time,ddq_xy(1,1),ddq_xy(1,2),ddq_xy(1,3),ddq_xy(2,2),ddq_xy(2,3),ddq_xy(3,3)

 close(iu)
 close(iuu)
 close(iuuu)
 print*, 'aussie crocos'

end subroutine do_analysis

end module analysis
