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
 use units,          only: umass,udist,utime
 use physcon,        only:gg,c
 use options,        only: iexternalforce
 use part,           only: isdead_or_accreted
 use io,             only:fatal
 character(len=*),   intent(in) :: dumpfile
 real,               intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(inout) :: pmass,time
 integer,            intent(in) :: npart,iunitone,numfile
 real,               dimension(3,npart):: fext
 integer             :: ierr,i
 real                ::dtextforce,dtf,poti,fextv(3)
 integer, parameter  :: iu = 1993
 logical, save       :: first = .true.
 logical             :: firstcall = .true.
 real                :: hp,hx,q(6),ddq(6),hp_,hx_,d
 real, save          :: distan

! read the particle accelerations
  dtextforce=huge(dtextforce)

  iexternalforce=14 !-->check

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

! initialise quadrupole to zero
  q(:)=0

! calculate the components of the traceless quadrupole--not necessary but maybe useful
  do i=1,npart
   if(xyzh(4,i)>tiny(xyzh)) then  !if not accreted
    q(1)=q(1)+pmass*(xyzh(1,i)*xyzh(1,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qxx
    q(2)=q(2)+pmass*(xyzh(1,i)*xyzh(2,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qxy
    q(3)=q(3)+pmass*(xyzh(1,i)*xyzh(3,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qxz
    q(4)=q(4)+pmass*(xyzh(2,i)*xyzh(2,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qyy
    q(5)=q(5)+pmass*(xyzh(2,i)*xyzh(3,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qyz
    q(6)=q(6)+pmass*(xyzh(3,i)*xyzh(3,i)-0.3*(xyzh(1,i)**2.+xyzh(2,i)**2.+xyzh(3,i)**2.)) !qzz
   end if
  enddo

! initialise the second time derivative of the quadrupole to zero
  ddq(:)=0
! calculate the second time derivative of the traceless quadrupole
 do i=1,npart
   if(xyzh(4,i)>tiny(xyzh)) then !if not accreted
    ddq(1)=ddq(1)+pmass*(2*vxyzu(1,i)*vxyzu(1,i)+xyzh(1,i)*fext(1,i)+xyzh(1,i)*fext(1,i)) !ddqxx
    ddq(2)=ddq(2)+pmass*(2*vxyzu(1,i)*vxyzu(2,i)+xyzh(1,i)*fext(2,i)+xyzh(2,i)*fext(1,i)) !ddqxy
    ddq(3)=ddq(3)+pmass*(2*vxyzu(1,i)*vxyzu(3,i)+xyzh(1,i)*fext(3,i)+xyzh(3,i)*fext(1,i)) !ddqxz
    ddq(4)=ddq(4)+pmass*(2*vxyzu(2,i)*vxyzu(2,i)+xyzh(2,i)*fext(2,i)+xyzh(2,i)*fext(2,i)) !ddqyy
    ddq(5)=ddq(5)+pmass*(2*vxyzu(2,i)*vxyzu(3,i)+xyzh(2,i)*fext(3,i)+xyzh(3,i)*fext(2,i)) !ddqyz
    ddq(6)=ddq(6)+pmass*(2*vxyzu(3,i)*vxyzu(3,i)+xyzh(3,i)*fext(3,i)+xyzh(3,i)*fext(3,i)) !ddqzz
   end if
 enddo


  if (firstcall) then
   firstcall=.false.
   call prompt('Write the distance from the source (Mpc):',d)
   distan=d*3.10e24!--convert Mpc in cm
  endif

! gw strain in the direction perpendicular to the orbit
  hx=gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*(ddq(1)-ddq(4))
  hp=gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*(ddq(2))

! gw strain in the plane of the orbit
  hx_=gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*(ddq(6)-ddq(4))
  hp_=-gg*c**(-4.)*distan**(-1.)*umass*udist**(2.)*utime**(-2.)*(ddq(5))


! Write a file where I append all the values of the strain wrt time
 if (first) then
  first = .false.
  open(unit=iu, file='strain.gw',status='replace')
   write(iu,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
    1, 'time in cu', &
    2, 'hp', &
    3, 'hx', &
    4, 'hp_', &
    5, 'hx_'
  else
   open(unit=iu, file='strain.gw',position='append')
   endif
  write(iu,'(5(es18.10,1X))') time,hp,hx,hp_,hx_
  close(iu)

  end subroutine do_analysis

end module analysis
