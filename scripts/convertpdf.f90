program convertpdf
 implicit none
 integer :: iarg,nargs,nheaderlines,inum,ierr,i,ires,lstr
 integer, parameter :: iunit = 1,iout=2
 character(len=128) :: filename,fileout
 character(len=6) :: cdum
 character(len=2) :: c1
 character(len=30) :: var
 real :: vallnrho,dum,valpdf
 
 nargs = command_argument_count()
 if (nargs.lt.1) stop 'Usage: convertpdf file(s)'
 
 do iarg=1,nargs
    call get_command_argument(iarg,filename)
    
    read(filename,"(a2,i3,a1,i4,a)") c1,ires,cdum,inum
    open(unit=iunit,file=filename,status='old',form='formatted',action='read')
    nheaderlines = 11
    do i=1,nheaderlines
       read(iunit,*)
    enddo
    
    var = 'density'
    lstr = len_trim(filename)
    if (lstr.ge.6) then
       if (filename(lstr-5:lstr)=='NOITER') then
          var = 'density_tracersNOITER'
       elseif (filename(lstr-3:lstr)=='ITER') then
          var = 'density_tracersITER'
       endif
    endif
!
!--open output file
!
    if (ires.eq.512) then
       write(fileout,"(a,i3.3,a)") 'turb',inum*2,'_pdf_'//trim(var)//'.dat'
    else
       write(fileout,"(a,i3.3,a)") 'turb',inum+1,'_pdf_'//trim(var)//'.dat'
    endif
    print "(/,a)",trim(filename)//' -> '//trim(fileout)
    open(unit=iout,file=fileout,status='replace')

    i = 0
    ierr = 0
    do while (ierr.eq.0)
       i = i + 1
       read(iunit,*,iostat=ierr) vallnrho,dum,valpdf
       if (ierr.eq.0) write(iout,*) exp(vallnrho),valpdf
    enddo
    print*,' last bin = ',i-1
    close(iunit)
    close(iout)
 enddo
 
end program convertpdf
