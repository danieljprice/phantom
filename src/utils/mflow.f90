!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: mflow
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: mflow[int time] file01.mf file02.mf ...
!
!  DEPENDENCIES: mf_write
!+
!--------------------------------------------------------------------------
program mflow
 use mf_write, only: nradi,ncolsi
 implicit none

 integer ::lu=11,iout=15,nthoutput=313
 integer ::ncols, ndump=0, nrad,inttime,nthcount
 integer ::nargs,istart,imf,nlines,ierr,i
 character(len=120) ::filename,outfile,formathead,formatout,formatint,num,nthfile,nth
 real, allocatable  :: dat(:)
 real, allocatable  :: datprev(:),datflow(:),rad(:)
 integer, allocatable :: intcol(:)
 real               :: tprev, pmass
 logical            :: nthtime

 !
 !---If more radii needed adjust it in mf_write.f90
 !

 ncols=ncolsi
 nrad=nradi
 write(num,*) (nrad)
 formathead="(2(a18,1x),"//trim(num)//"(es18.10,1x))"
 write(num,*) nrad
 formatout="(1(1x,es18.10,1x),(I18,1x),"//trim(num)//"(es18.10,1x))"
 write(num,*) (ncols-2)
 formatint="('#',"//trim(num)//"(I18,1x))"
 pmass=1.1000132002E-08

 allocate(datprev(nrad),dat(ncols),datflow(nrad),rad(nrad),intcol((ncols-2)))

 do i=1,ncols-2
    intcol(i)=i
 enddo


 nargs = command_argument_count()
 call get_command_argument(0,filename)

 if (nargs <= 0) then
    print "(a)",' Usage: '//trim(filename)//'[int time] file01.mf file02.mf ...'
    stop
 endif

 call get_command_argument(1,filename)
 read(filename,'(I7.7)',iostat=ierr) inttime
 if (ierr /= 0) then
    inttime = 0
    istart = 1
    nthtime=.false.
 else
    print*, 'Creating also flow at',inttime,'-th time '
    istart = 2
    nthtime=.true.
    write(nth,'(I5.5)')inttime

 endif



 over_args: do i=istart,nargs
    call get_command_argument(i,filename)
    open(unit=lu,file=trim(filename),status='old',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR opening '//trim(filename)//' ...skipping'
    else
       !
       !--- open output file
       !
       imf = index(filename,'.mf') - 1

       if (imf <= 0) imf = len_trim(filename)
       outfile = trim(filename(1:imf))//'.mflow'
       print*,trim(filename)//' --> '//trim(outfile)
       nthfile=trim(filename(1:imf))//'-'//trim(nth)//'.mflow'

       if(nthtime) print*,trim(outfile)//'--->'//trim(nthfile)
       open(unit=iout,file=trim(outfile),status='replace',form='formatted')

       ierr=1
       nlines = 0
       rad(:)=0
       dat(:)=0

       do while(nlines < 4 )
          nlines = nlines + 1
          read(lu,*)
       enddo

       read(lu,*,iostat=ierr) rad(1:nrad)
       write(iout,formatint) intcol(1:(ncols-2))
       write(iout,formathead)'t','ndump',rad(:)

       !
       !---skip header lines
       !
       read(lu,*)
       nlines = nlines + 1
       read(lu,*,iostat=ierr) dat(1:ncols)

       tprev=dat(1)
       datprev(:)=dat(5:ncols)

       do while(ierr==0)

          read(lu,*,iostat=ierr) dat(1:ncols)

          datflow(:)=(dat(5:ncols)-datprev(1:nrad))/(dat(1)-tprev)!pmass

          datprev(:)=dat(5:ncols)
          tprev=dat(1)
          if(ierr==0) write(iout,formatout) dat(1),ndump,datflow

          if(nthtime) then
             if(ndump==inttime)then
                open(unit=nthoutput,file=trim(nthfile),status='replace',form='formatted')
                write(nthoutput,"('#',(1x,a17),es18.10)")"time:",dat(1)
                do nthcount=1,nrad
                   write(nthoutput,"(2(1x,es18.10,1x))") rad(nthcount), datflow(nthcount)
                enddo
                close(nthoutput)
             endif
          endif

          ndump=ndump+1
       enddo
       close(unit=iout)
    endif
    close(unit=lu)

 enddo over_args

end program mflow
