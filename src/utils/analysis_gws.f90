module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'gws'

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunitone)
 use externalforces, only:initialise_externalforces,update_externalforce,externalforce,externalforce_vdependent
 use timestep,       only:C_force
 use prompting,      only:prompt
 use units,          only:utime
 use options,        only:iexternalforce
 use part,           only:isdead_or_accreted
 use io,             only:fatal
 use gravwaveutils,  only:calculate_strain
 character(len=*),   intent(in) :: dumpfile
 real,               intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(inout) :: pmass,time
 integer,            intent(in) :: npart,iunitone,numfile
 real,               dimension(3,npart):: fext
 integer             :: ierr,i
 real                :: dtextforce,dtf,poti,fextv(3)
 integer, parameter  :: iu = 1993
 logical, save       :: first = .true.
 logical, save       :: firstcall = .true.
 logical, save       :: firstdump = .true.
 real                :: hp,hx,hpp,hxx
 real,save           :: d

! read the particle accelerations
 dtextforce=huge(dtextforce)

 if (firstdump) then
    firstdump=.false.
    call prompt('Write the external force [1:16]:', iexternalforce)
    print*, 'WARNING: you have to check that the external force used is reasonable for a physical point of view.'
 endif

 call initialise_externalforces(iexternalforce,ierr)
 call update_externalforce(iexternalforce,time,0.)
 if (ierr /= 0) call fatal('initial','error in external force settings/initialisation')
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,fext,time,iexternalforce,C_force) &
 !$omp private(i,poti,dtf,fextv) &
 !$omp reduction(min:dtextforce)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call externalforce(iexternalforce,xyzh(1,i),xyzh(2,i),xyzh(3,i), &
      xyzh(4,i),time,fext(1,i),fext(2,i),fext(3,i),poti,dtf,i)
       dtextforce = min(dtextforce,C_force*dtf)
       ! add velocity-dependent part
       call externalforce_vdependent(iexternalforce,xyzh(1:3,i),vxyzu(1:3,i),fextv,poti)
       fext(1:3,i) = fext(1:3,i) + fextv
    endif
 enddo
 !$omp end parallel do

 ! Returns the gravitational wave strain at 1Mpc
 call calculate_strain(hx,hp,hxx,hpp,xyzh,vxyzu(1:3,:),fext(1:3,:),pmass,npart)

 if (firstcall) then
    firstcall=.false.
    d = 1.
    call prompt('Write the distance from the source (Mpc):',d)
 endif

! Write a file where I append all the values of the strain wrt time
 if (first) then
    first = .false.
    open(unit=iu, file='strain.gw',status='replace')
    write(iu,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
    1, 'time', &
    2, 'hp', &
    3, 'hx', &
    4, 'hpp', &
    5, 'hxx'
 else
    open(unit=iu, file='strain.gw',position='append')
 endif

 ! Rescale the strains by the desired distance
 write(iu,'(5(es18.10,1X))') time*utime,hp/d,hx/d,hpp/d,hxx/d

 close(iu)

end subroutine do_analysis

end module analysis
