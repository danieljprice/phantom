!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom2grid
!
!  DESCRIPTION: Interpolate Phantom data to a grid
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom2grid ncells gridformat dumpfilename(s)
!
!  DEPENDENCIES: boundary, dim, hdf5utils, interpolations3D, io,
!    iso_c_binding, part, readwrite_dumps
!+
!--------------------------------------------------------------------------
program phantom2grid
 use part,             only:rhoh,xyzh,massoftype,hfact,vxyzu,npart,Bxyz
 use io,               only:set_io_unit_numbers,iprint,idisk1,iverbose,real4,formatreal,cstring
 use boundary,         only:xmin,ymin,zmin,xmax,ymax,zmax
 use readwrite_dumps,  only:read_dump,write_gadgetdump
 use interpolations3D, only:interpolate3D
#ifdef HDF5
 use hdf5utils,        only:write_grid_hdf5
 use iso_c_binding,    only:c_float
#endif
 implicit none
 real, allocatable :: datgrid(:,:,:,:) !rho, vx, vy, vz, Bx, By, Bz
 real               :: weight,time,pixwidth
 real :: partval(8),partmin(8),partmax(8),partmean(8)
 real :: gridmin(8),gridmax(8),gridmean(8)
! real, dimension(maxp) :: rho
 real               :: offset,rhomin,dati
! real               :: trans(6),datmax
! real, dimension(npixx,npixy) :: coldens
#ifdef HDF5
 real(kind=c_float), allocatable :: dattemp(:)
 character(len=120) :: filenameout
 integer            :: inum,n
#endif
 integer            :: maxp,maxmhd
 integer            :: npixx,npixy,npixz
 integer            :: nargs,iarg
 integer            :: ifile,k,j,i,ierr,mx,my,mz,mv,irec
 integer            :: ioutformat,nx,ilendat
 character(len=20)  :: string
 character(len=120) :: dumpfile,fileout(8),fileprefix
 integer, parameter :: maxoutformats = 6
 character(len=7), dimension(maxoutformats), parameter :: label = &
  (/'spyros ','spyros ','spyros ','kitp   ','gridbin','hdf5   '/)
 namelist /dim/ mx,my,mz,mv,offset

 call set_io_unit_numbers
 iprint = 6
 iverbose = 2
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 3) then
    call print_usage
    stop
 endif
!
!--get number of grid cells from the command line
!
 npixx = 0
 call get_command_argument(1,string)
 read(string,*,iostat=ierr) npixx
 if (ierr /= 0 .or. npixx <= 0 .or. npixx > 100000) then
    call print_usage
    print "(/,a,i5,a,/)",' *** Error in number of pixels from command line (got ',npixx,')'
    stop
 endif
 npixy = npixx
 npixz = npixx
!
!--get the output format from the command line,
!  as either integer or as a string matching the label
!
 ioutformat = 0
 call get_command_argument(2,string)
 read(string,*,iostat=ierr) ioutformat
 do while (ioutformat <= 0 .or. ioutformat > maxoutformats)
    do i=maxoutformats,1,-1
       if (trim(label(i))==trim(adjustl(string(1:len(label(i)))))) then
          ioutformat = i
          print "(a,i2,a)",' got output format '//trim(label(i))//'(iout = ',ioutformat,') from command line'
       endif
    enddo
    if (ioutformat <= 0 .or. ioutformat >= maxoutformats) then
       call print_usage
       print "(/,a,i8,a,/)",' *** Error reading output format from command line (got ',ioutformat,')'
       stop
    endif
 enddo

 print "(/,a,/)",' Phantom2grid: we are pleased to do you service'

 if (IACHAR(TRANSFER(1,"a")) == 0) then
    print "(a,/)",' native endian is BIG'
 else
    print "(a,/)",' native endian is LITTLE'
 endif

 over_args: do iarg=3,nargs

    call get_command_argument(iarg,dumpfile)
!
!--read particle setup from dumpfile
!
    call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) then
       if (allocated(datgrid)) deallocate(datgrid)
       stop 'error reading dumpfile'
    endif

    !
    !--warn if grid is not square
    !
    if ((abs(abs(xmax-xmin)-abs(ymax-ymin)) > tiny(xmin)) .or. &
     (abs(abs(xmax-xmin)-abs(zmax-zmin)) > tiny(xmin)) .or. &
     (abs(abs(ymax-ymin)-abs(zmax-zmin)) > tiny(xmin))) then
       print "(3(6x,'[ ',f8.2,' < ',a,' < ',f8.2,' ]',/))",xmin,'x',xmax,ymin,'y',ymax,zmin,'z',zmax
       print "(a)",' *** ERROR: boundaries indicate that grid is not square'
       print "(a,/)",' *** INTERPOLATION FOR NON-SQUARE GRIDS IS NOT CURRENTLY IMPLEMENTED'
       stop
    endif

! do i=1,npart
!    rho(i) = rhoh(xyzh(4,i),massoftype(1))
! enddo
! call write_gadgetdump(trim(fileprefix),time,xyzh,particlemass,vxyzu,rho,1.5*polyk,npart)
! cycle over_args

    !
    !--allocate memory for the grid
    !
    if (.not.allocated(datgrid)) then
       write(*,"(a,i5,2(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixx,npixy,npixz,' grid ...'
       maxp = size(xyzh,2)   ! do this to avoid clash with namelist
       maxmhd = size(Bxyz,2)
       ilendat = 4 + 4*(maxmhd/maxp)
       allocate(datgrid(ilendat,npixx,npixy,npixz),stat=ierr)
       if (ierr /= 0) then
          write(*,*) 'FAILED: NOT ENOUGH MEMORY'
          stop
       else
          write(*,*) 'OK'
       endif
    endif
!
!--interpolate to the 3D grid
!
    pixwidth = (xmax-xmin)/npixx
    weight = 1.d0/hfact**3
    print*,'using hfact = ',hfact,' weight =',weight,massoftype(1)/(rhoh(xyzh(4,1),massoftype(1))*xyzh(4,1)**3)
    print*,' xmin,xmax = ',xmin,xmax
    print*,' ymin,ymax = ',ymin,ymax
    print*,' zmin,zmax = ',zmin,zmax
    if (maxmhd==maxp) then
       print "(a)",' interpolating hydro + MHD quantities to grid...'
       call interpolate3D(xyzh,weight,massoftype(1),vxyzu,npart,xmin,ymin,zmin,datgrid,npixx,npixy,npixz,&
                          pixwidth,.true.,4,Bxyz)
    else
       print "(a)",' interpolating hydro quantities to grid...'
       call interpolate3D(xyzh,weight,massoftype(1),vxyzu,npart,xmin,ymin,zmin,datgrid,npixx,npixy,npixz,&
                       pixwidth,.true.,0)
    endif

    if (maxmhd==maxp) then
       print "(22x,7(a10,1x))",'    rho   ','    vx    ','    vy    ','    vz    ',&
                                         '    Bx    ','    By    ','    Bz    '
    else
       print "(22x,4(a10,1x))",'    rho   ','    vx    ','    vy    ','    vz    '
    endif
!
!--calculate max and min values on grid
!
    gridmax(:) = -huge(gridmax)
    gridmin(:) = huge(gridmin)
    gridmean(:) = 0.
    rhomin = huge(rhomin)
    !$omp parallel do schedule(static) &
    !$omp reduction(min:gridmin,rhomin) &
    !$omp reduction(max:gridmax) &
    !$omp reduction(+:gridmean) &
    !$omp private(k,j,i,ifile,dati)
    do k=1,npixz
       do j=1,npixy
          do i=1,npixx
             do ifile=1,ilendat
                dati = datgrid(ifile,i,j,k)
                gridmax(ifile) = max(gridmax(ifile),dati)
                gridmin(ifile) = min(gridmin(ifile),dati)
                gridmean(ifile) = gridmean(ifile) + dati
                if (ifile==1 .and. dati > 0.) rhomin = min(dati,rhomin)
             enddo
          enddo
       enddo
    enddo
    gridmean(:) = gridmean(:)/(npixx*npixy*npixz)
!
!--calculate max and min values on particles
!
    partmax(:) = -huge(partmax)
    partmin(:) = huge(partmin)
    partmean(:) = 0.0
    !$omp parallel do reduction(min:partmin) &
    !$omp reduction(max:partmax) &
    !$omp reduction(+:partmean) &
    !$omp private(j,partval)
    do i=1,npart
       partval(1) = xyzh(4,i)
       partval(2) = vxyzu(1,i)
       partval(3) = vxyzu(2,i)
       partval(4) = vxyzu(3,i)
       if (maxmhd==maxp) then
          partval(5) = Bxyz(1,i)
          partval(6) = Bxyz(2,i)
          partval(7) = Bxyz(3,i)
       endif
       do j=1,ilendat
          partmin(j) = min(partmin(j),partval(j))
          partmax(j) = max(partmax(j),partval(j))
          partmean(j) = partmean(j) + partval(j)
       enddo
    enddo
    partmean(:) = partmean(:)/npart

    print "(a,7(1pe10.3,1x))",'max (on particles)  = ',rhoh(partmin(1),massoftype(1)),partmax(2:ilendat)
    print "(a,7(1pe10.3,1x))",'max      (on grid)  = ',gridmax(1:ilendat)
    print "(a,7(1pe10.3,1x))",'min (on particles)  = ',rhoh(partmax(1),massoftype(1)),partmin(2:ilendat)
    print "(a,7(1pe10.3,1x))",'min      (on grid)  = ',gridmin(1:ilendat)
    print "(a,7(1pe10.3,1x))",'mean (on particles) = ',massoftype(1)*npart/((xmax-xmin)*(ymax-ymin)*(zmax-zmin)),&
                                                    partmean(2:ilendat)
    print "(a,7(1pe10.3,1x))",'mean      (on grid) = ',gridmean(1:ilendat)

    print*,' total mass = ',massoftype(1)*npart
    print*,' rho mean (mass/vol) = ',massoftype(1)*npart/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))
    print*,' rho mean  (on grid) = ',gridmean(1)
!
!--set minimum density on the grid
!
    print*,'enforcing minimum rho on grid = ',rhomin

    !$omp parallel do private(k) schedule(static)
    do k=1,npixz
       where (datgrid(1,:,:,k) <= tiny(datgrid))
          datgrid(1,:,:,k) = rhomin
       end where
    enddo
!
!--calculate and plot column density on the grid as a check
!
! do j=1,npixy
!   do i=1,npixx
!      coldens(i,j) = datgrid(1,i,j,npixz/2) !sum(datgrid(1,i,j,:))/(zmax-zmin)
!      if (any(datgrid(1,i,j,:) <= tiny(datgrid))) then
!         do k=1,npixz
!           if (datgrid(1,i,j,k) <= tiny(datgrid)) print*,i,j,k,' dat = ',datgrid(1,i,j,k)
!         enddo
!     endif
!   enddo
! enddo
! call pgopen('/xw')
! call pgenv(xmin,xmax,ymin,ymax,1,0)
! trans = 0.
! trans(1) = xmin - 0.5*pixwidth
! trans(2) = pixwidth
! trans(4) = ymin - 0.5*pixwidth
! trans(6) = pixwidth
! datmax = maxval(coldens(1:npixx,1:npixy))
! call pgimag(coldens,npixx,npixy,1,npixx,1,npixy,0.0,datmax,trans)
! call pgend
!
!--write to the output files, using Spyros' naming convention
!
    select case(label(ioutformat))
    case('spyros')
!
!--construct output filenames, using Spyros' naming convention
!  (t=2 is 020, t=12 is 120)
!
       ifile = nint(10.*time)
       write(fileprefix,"(a6,i3.3,'_',i3.3)") dumpfile(1:2)//'_ic_',npixx,ifile
       print*,' filename = ',trim(fileprefix)
       fileout(1) = trim(adjustl(fileprefix))//'_dens.dat'
       fileout(2) = trim(adjustl(fileprefix))//'_xvel.dat'
       fileout(3) = trim(adjustl(fileprefix))//'_yvel.dat'
       fileout(4) = trim(adjustl(fileprefix))//'_zvel.dat'
       !--binary files
       fileout(5) = 'BDENS'//trim(adjustl(fileprefix))
       fileout(6) = 'BXVEL'//trim(adjustl(fileprefix))
       fileout(7) = 'BYVEL'//trim(adjustl(fileprefix))
       fileout(8) = 'BZVEL'//trim(adjustl(fileprefix))

       do ifile=1,8
          if (ifile <= 4 .and. (ioutformat==1 .or. ioutformat==3)) then
             print "(a)",' writing ascii output to file '//trim(fileout(ifile))
             print*,' max =',gridmax(ifile)
             open(unit=16+ifile,file=trim(fileout(ifile)),form='formatted')
             !write(16+ifile) real4(time),real4(gridmax(ifile)),real4(npixx)
             do k=1,npixz
                do j=1,npixy
                   do i=1,npixx
                      write(16+ifile,*) real4(datgrid(ifile,i,j,k))
                   enddo
                enddo
             enddo
             close(unit=16+ifile)
          elseif (ifile > 4 .and. (ioutformat==1 .or. ioutformat==2)) then
             print "(a)",' writing binary output to file '//trim(fileout(ifile))
             print*,' max =',gridmax(ifile-4)
             print*,' first few numbers = ',(real4(datgrid(ifile-4,i,1,1)),i=1,5)
             open(unit=16+ifile,file=trim(fileout(ifile)),form='unformatted')
             if (ifile==5) then
                write(16+ifile) real4(time),real4(gridmax(ifile-4)),real4(npixx)
             else
                write(16+ifile) real4(time),real4(gridmax(ifile-4)),real4(npixx),real4(-1.0)
             endif
             write(16+ifile) (((real4(datgrid(ifile-4,i,j,k)),i=1,npixx),j=1,npixy),k=1,npixz)
             close(unit=16+ifile)
          endif
       enddo
    case('kitp')
       if (npart > 1e9) then
          nx = 1024
       elseif (npart > 130e6) then
          nx = 512
       elseif (npart > 16e6) then
          nx = 256
       else
          nx = 128
       endif
       if (nx > 999) then
          write(fileprefix,"(a,i4.4,a,i3.3,a,f4.2)") 'price',nx,'grid',npixx,'_t',time
       else
          write(fileprefix,"(a,i3.3,a,i3.3,a,f4.2)") 'price',nx,'grid',npixx,'_t',time
       endif
       fileout(1) = trim(adjustl(fileprefix))//'.dim'
       fileout(2) = trim(adjustl(fileprefix))//'.dat'
       print "(a)",' output files= ',trim(fileout(1)),trim(fileout(2))
       open(10,file=trim(fileout(1)),form='formatted',status='unknown')
       mx = npixx
       my = npixy
       mz = npixz
       if (maxmhd==maxp) then
          mv = 7
       else
          mv = 4
       endif
       offset = 0.
       write(10,NML=dim)
       close(10)

       open(unit=16,file=trim(fileout(2)),form='unformatted',access='direct',status='unknown',recl=npixx*npixy*npixz) !*4)
       do irec=1,mv
          write(16,rec=irec) (((real4(datgrid(irec,i,j,k)),i=1,npixx),j=1,npixy),k=1,npixz)
       enddo
       close(unit=16)
    case('mine','gridbin')
       !if (npart > 1e9) then
       !   nx = 1024
       !elseif (npart > 130e6) then
       !   nx = 512
       !elseif (npart > 16e6) then
       !   nx = 256
       !else
       !   nx = 128
       !endif
       !if (nx > 999) then
       !   write(fileprefix,"(a,i4.4,a,i3.3,a)") 'price',nx,'grid',npixx
       !else
       !   write(fileprefix,"(a,i3.3,a,i3.3,a)") 'price',nx,'grid',npixx
       !endif
       !call formatreal(time,string,precision=1.e-2)
       !fileout(1) = trim(adjustl(fileprefix))//'_t'//trim(string)//'.dat'

       write(fileout(1),"(a,i3.3,a)") trim(dumpfile)//'_grid',npixx,'.grid'
       print "(a)",' output file= '//trim(fileout(1))

       if (maxmhd==maxp) then
          mv = 7
       else
          mv = 4
       endif
       open(unit=16,file=trim(fileout(1)),form='unformatted',status='replace')
       write(16) npixx,npixy,npixz,mv
       do irec=1,mv
          write(16) (((real4(datgrid(irec,i,j,k)),i=1,npixx),j=1,npixy),k=1,npixz)
       enddo
       close(unit=16)
#ifdef HDF5
    case('hdf5')
       if (npart > 1e9) then
          nx = 1024
       elseif (npart > 130e6) then
          nx = 512
       elseif (npart > 16e6) then
          nx = 256
       else
          nx = 128
       endif
       inum = nint(time/0.005)
       if (nx > 999) then
          write(fileprefix,"(a,i4.4,a,i4.4,a)") 'HD',nx,'_',inum,'_'
!       write(fileprefix,"(a,i4.4,a,i3.3,a)") 'price',nx,'grid',npixx
       else
          write(fileprefix,"(a,i3.3,a,i4.4,a)") 'HD',nx,'_',inum,'_'
!       write(fileprefix,"(a,i3.3,a,i3.3,a)") 'price',nx,'grid',npixx
       endif
!    call formatreal(time,string,precision=1.e-2)
       if (maxmhd==maxp) then
          mv = 7
       else
          mv = 4
       endif
       fileout(1) = 'dens'
       fileout(2) = 'velx'
       fileout(3) = 'vely'
       fileout(4) = 'velz'
       fileout(5) = 'Bx'
       fileout(6) = 'By'
       fileout(7) = 'Bz'

       allocate(dattemp(npixx*npixy*npixz))
       do irec=1,mv
!       filenameout = trim(adjustl(fileprefix))//'_t'//trim(string)//'_'//trim(fileout(irec))//'.dat'
          filenameout = trim(adjustl(fileprefix))//trim(fileout(irec))
          print "(a)",' output file= '//trim(filenameout)//' writing dataset '//trim(fileout(irec))
          n = 0
          do k=1,npixz
             do j=1,npixy
                do i=1,npixx
                   n = n + 1
                   dattemp(n) = datgrid(irec,i,j,k)
                enddo
             enddo
          enddo
          call write_grid_hdf5(cstring(filenameout),cstring(fileout(irec)),dattemp,npixx,npixy,npixz,ierr)
          if (ierr /= 0) then
             print*,'*** ERRORS writing to '//trim(filenameout)
             cycle over_args
          endif
       enddo
       deallocate(dattemp)

#endif
    end select

 enddo over_args

 if (allocated(datgrid)) deallocate(datgrid)

 print "(/,a,/)",' Phantom2grid: we thank you for your kind custom'

contains

subroutine print_usage
 use dim, only:tagline

 print "(a,/)",trim(tagline)
 print "(a)",' Usage: phantom2grid ncells gridformat dumpfilename(s)'

 print "(/,2x,a,f6.2,a,f6.2,a)",'ncells     : number of grid cells along x axis (',xmin,' < x < ',xmax,')'
 print "(2x,a,5(/,14x,a))",'gridformat : 1 = Potsdam comparison format: ascii & binary ', &
            ' 2 = Potsdam comparison format: binary only', &
            ' 3 = Potsdam comparison format: ascii only', &
            ' 4 = KITP comparison format ', &
            ' 5 = gridbinary format (KITP without recl)', &
#ifdef HDF5
            ' 6 = FLASH-style HDF5 files'
#else
 ' 6 = FLASH-style HDF5 files (NOT AVAILABLE: compile with -DHDF5)'
#endif

end subroutine print_usage

end program phantom2grid

