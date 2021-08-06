module analysis
implicit None
character(len=20), parameter, public :: analysistype = 'test_analysis'
public :: do_analysis

logical, private :: firstcall = .true.
private

contains
  subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
    use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass
    use centreofmass, only: get_centreofmass
    character(len=*), intent(in) :: dumpfile
    integer,          intent(in) :: num,npart,iunit
    real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
    real,             intent(in) :: particlemass,time
    logical                      :: iexist
    character(len=200)           :: fileout

    fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_com.dat'
    inquire(file=fileout,exist=iexist)
    if ( .not.iexist .or. firstcall ) then
       firstcall = .false.
       open(iunit,file=fileout,status='replace')
       write(iunit,*)'hello world!'
     else
        open(iunit,file=fileout,position='append')

     endif
     write(iunit,*) trim(dumpfile)
     close(iunit)
  end subroutine do_analysis
end module analysis
