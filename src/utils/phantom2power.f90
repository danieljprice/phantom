!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom2power
!
!  DESCRIPTION: This program is a utility for calculating power spectra
!   direct from phantom dumps
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom2power dumpfilename [npixx]
!
!  DEPENDENCIES: boundary, dim, fft3d, infile_utils, interpolations3D, io,
!    part, physcon, power, powerspec, readwrite_dumps
!+
!--------------------------------------------------------------------------
program phantom2power
 use dim,  only:tagline,periodic,maxp_hard
 use part, only:xyzh,vxyzu,npart,massoftype,rhoh,hfact
 use io, only:set_io_unit_numbers,iprint,idisk1,real4
 use boundary, only:xmin,ymin,zmin,dxbound,dybound,dzbound
 use readwrite_dumps, only:read_dump
 use interpolations3D, only:interpolate3D
 use powerspec, only:sphfft3D,sphfft3D_fast,sphpow3D,power_part,write_power_ascii,open_power_ascii
 use fft3d, only:fft3df
 use power, only:open_power,write_power,power3d,mlabel,power_unit
 use physcon, only:pi
 implicit none
 integer :: npixx,npixy,npixz,numk
 integer :: nargs,ierr,n_out,i,j,k
 integer, allocatable :: nk(:)
 integer, parameter :: iunit=12
 real, allocatable :: pk(:),xk(:)
 real, allocatable :: datgrid(:,:,:,:) !rho, vx, vy, vz
 real, allocatable :: dati(:,:,:),datx(:,:,:),daty(:,:,:)
 real :: rho(maxp_hard)
 real :: comp(12),xminpart(3),xmaxpart(3),xminp(3),dxmaxp(3)
 real :: ptot,dxmin,compensate,totvol,totmass,totrho2,rhoi,rhopow,rhopowprev,rhomin,rhomax
 real :: weight,time,pixwidth,ekin,ekinrho,ekinrho23,ekx,eky,ekz,dvol,vx2,vy2,vz2
 logical, parameter :: interpolate = .true.
 logical, parameter :: use_rhopow_from_particles = .false.
 logical :: do_average, read_input_file
 character(len=120) :: dumpfile,fileout
 character(len=mlabel) :: variable
 character(len=20) :: string
 character(len=1), parameter :: labelx(3) = (/'x','y','z'/)
 character(len=20), parameter :: filename = 'phantom2power.in'

 call set_io_unit_numbers
 iprint = 6
 npixx = 256
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1 .or. nargs > 2) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantom2power dumpfilename [npixx]'
    print "(a)",'        (default npixx=256)'
    stop
 endif
 call get_command_argument(1,dumpfile)

 print "(/,a,/)",' Phantom2power: we are pleased to do you service'

 ! set boundaries from
 xminp(:)  = (/xmin,ymin,zmin/)
 dxmaxp(:) = (/dxbound,dybound,dzbound/)

 read_input_file = .false.
 call read_analysis_options(xminp,dxmaxp,filename,iunit,ierr)
 if (ierr == 0) read_input_file = .true.

 if (.not.read_input_file .and. periodic) then
    call write_analysis_options(xminp,xminp+dxmaxp,filename,iunit+1)
    stop
 endif

 if (nargs >= 2) then
    call get_command_argument(2,string)
    read(string,*,iostat=ierr) npixx
    if (ierr /= 0) stop 'error reading npixx from command line'
 endif
 if (mod(npixx,2**(nint(log(real(npixx))/log(2.)))) /= 0) stop 'error: npixx must be power of 2 for fft'

 pixwidth = dxmaxp(1)/npixx
 npixy = nint(dxmaxp(2)/pixwidth)
 npixz = nint(dxmaxp(3)/pixwidth)
! npixy = npixx; npixz = npixx

 if (interpolate) then
    print "(a,i4,a)",' (allocating memory for ',npixx,'^3 grid)'
    allocate(datgrid(4,npixx,npixy,npixz),dati(npixx,npixy,npixz))
 endif

 numk = npixx/2
!
!--allocate memory for k-space arrays
!
 print "(a,i3,a)",' (allocating memory for ',numk,' k-values)'
 allocate(nk(numk),pk(numk),xk(numk))
!
!--read particle setup from dumpfile
!
 call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading dumpfile'

 n_out = 1
 compensate = 0. ! this does compensation *inside* power3d - makes very little difference
 ! to the final power spectrum so we don't use this.
!
! compensation factors -
! these are just used to write extra columns in the ascii powerspectrum file
!
 comp(1) = 1./3.
 comp(2) = 1./2.
 comp(3) = 2./3.
 comp(4) = 1.
 comp(5) = 4./3.
 comp(6) = 1.5
 comp(7) = 5./3.
 comp(8) = 2.
 comp(9) = 7./3.
 comp(10) = 2.5
 comp(11) = 8./3.
 comp(12) = 3.
 comp(:) = comp(:) - compensate
 do_average = .true.
!
! calculate Ekin on the particles to check against the values on the grid
!
 ekin = 0.
 ekinrho = 0.
 ekinrho23 = 0.
 ekx = 0.
 eky = 0.
 ekz = 0.
 totrho2 = 0.
 rhomin = huge(rhomin)
 rhomax = 0.
 totvol = 0.
 print*,'hfact = ',hfact
 xminpart = huge(xminpart)
 xmaxpart = -xminpart
 do i=1,npart
    do j=1,3
       xminpart(j) = min(xminpart(j),xyzh(j,i))
       xmaxpart(j) = max(xmaxpart(j),xyzh(j,i))
    enddo
    rho(i) = rhoh(xyzh(4,i),massoftype(1))
    dvol = massoftype(1)/rho(i)
    totrho2 = totrho2 + massoftype(1)*rho(i)
    rhomin = min(rhomin,rho(i))
    rhomax = max(rhomax,rho(i))
    vx2 = vxyzu(1,i)**2
    vy2 = vxyzu(2,i)**2
    vz2 = vxyzu(3,i)**2
    ekin = ekin + dvol*(vx2 + vy2 + vz2)
    ekx = ekx + dvol*vx2
    eky = eky + dvol*vy2
    ekz = ekz + dvol*vz2
    ekinrho = ekinrho + massoftype(1)*(vx2 + vy2 + vz2)
    ekinrho23 = ekinrho23 + massoftype(1)*rho(i)**(-1./3.)*(vx2 + vy2 + vz2)
    totvol = totvol + dvol
 enddo
 print*,'On parts: Total Volume =',totvol
 print*,'On parts: Ekin(vol) = ',Ekin,' Ekin_x = ',ekx,' Ekin_y = ',eky,' Ekin_z = ',ekz
 print*,'On parts: Ekin(rho) = ',Ekinrho,' Total Mass = ',massoftype(1)*npart,' int(rho^2)dV = ',totrho2
 print*,'On parts: Ekin((rho1/3 * v)^2) =',Ekinrho23,' max. dens = ',rhomax,' min. dens = ',rhomin

 if (.not.periodic) then
    print *,' MIN, MAX of particle distribution is:'
    call print_bounds(xminpart,xmaxpart)
    if (.not.read_input_file) then
       call write_analysis_options(xminpart,xmaxpart,filename,iunit+1)
       stop
    endif
    print *,' MIN, MAX of grid is:'
    call print_bounds(xminp,xminp+dxmaxp)
 endif
!---------------------------------------------
 if (interpolate) then
!
!--interpolate to the 3D grid
!
    write(fileout,"(a,i3.3,a)") trim(dumpfile)//'_grid',npixx,'.pow'

    !pixwidth = dxmaxp(1)/npixx
    weight = 1.d0/hfact**3
    print*,'using hfact = ',hfact,' weight = ',weight

    call interpolate3D(xyzh,weight,massoftype(1),vxyzu,npart,xminp(1),xminp(2),xminp(3),datgrid, &
                       npixx,npixy,npixz,pixwidth,.true.,0)
!
!--calculate Ekin on the grid as a check (both against the particles and against Parseval's theorem)
!
    ekin = 0.
    ekinrho = 0.
    ekinrho23 = 0.
    ekx = 0.
    eky = 0.
    ekz = 0.
    totmass = 0.
    totrho2 = 0.
    rhomin = huge(rhomin)
    rhomax = 0.
    do k=1,npixz
       do j=1,npixy
          do i=1,npixx
             rhoi = datgrid(1,i,j,k)
             vx2 = datgrid(2,i,j,k)**2
             vy2 = datgrid(3,i,j,k)**2
             vz2 = datgrid(4,i,j,k)**2
             rhomin = min(rhomin,rhoi)
             rhomax = max(rhomax,rhoi)
             ekx = ekx + vx2
             eky = eky + vy2
             ekz = ekz + vz2
             ekin = ekin + vx2 + vy2 + vz2
             ekinrho = ekinrho + rhoi*(vx2 + vy2 + vz2)
             ekinrho23 = ekinrho23 + rhoi**(2./3.)*(vx2 + vy2 + vz2)
             totmass = totmass + rhoi
             totrho2 = totrho2 + rhoi*rhoi
          enddo
       enddo
    enddo
    totvol = dxbound*dybound*dzbound
    dvol = totvol/float(npixx*npixy*npixz)
    ekx = ekx*dvol
    eky = eky*dvol
    ekz = ekz*dvol
    ekin = ekin*dvol
    ekinrho = ekinrho*dvol
    ekinrho23 = ekinrho23*dvol
    totmass = totmass*dvol
    totrho2 = totrho2*dvol
    print*,'On grid: Ekin(vol) = ',Ekin,' Ekin_x = ',ekx,' Ekin_y = ',eky,' Ekin_z = ',ekz
    print*,'On grid: Ekin(rho) = ',Ekinrho,' Total Mass = ',totmass,' int(rho^2)dV = ',totrho2
    print*,'On grid: Ekin((rho1/3 * v)^2) = ',Ekinrho23,' max. dens = ',rhomax,' min. dens = ',rhomin
    rhopowprev = 0.
    rhopow = 0.
!
!--get the power spectrum for various combinations of rho and v
!  Kinetic energy power spectra (i > 4) require extra storage so these are done last
!
    do i=1,7
       select case(i)
       case(1)
          write(variable,10) 'grid',npixx,'rho'
          string = 'int(rho^2)dV'
       case(2)
          write(variable,10) 'grid',npixx,'vx'
          string = 'Ekin_x'
       case(3)
          write(variable,10) 'grid',npixx,'vy'
          string = 'Ekin_y'
       case(4)
          write(variable,10) 'grid',npixx,'vz'
          string = 'Ekin_z'
       case(5)
          write(variable,10) 'grid',npixx,'ekin'
          string = 'Ekin(vol)'
       case(6)
          rhopow = 1./3.
          write(variable,10) 'grid',npixx,'rho0.33ekin'
          string = 'Ekin((rho1/3 * v)^2)'
       case(7)
          rhopow = 0.5
          write(variable,10) 'grid',npixx,'rho0.50ekin'
          string = 'Ekin(rho)'
       end select
10     format (a,i3.3,a)
       if (i <= 4) then
          print*,'taking fft of '//trim(variable)
          !$omp parallel
          call fft3df(datgrid(i,:,:,:),dati,npixx,npixy,npixz)
          call power3d(dati,npixx,npixy,npixz,pk,xk,nk,numk,ptot,compensate,do_average)
          !$omp end parallel
       elseif (i==5) then
!
!--for Ekin, get fft of each component, then get total power from sum of squares
!
          print*,'(allocating memory for x,y components)'
          allocate(datx(npixx,npixy,npixz),daty(npixx,npixy,npixz))
          print 10,' taking fft of grid',npixx,'vy'
          !$omp parallel
          call fft3df(datgrid(3,:,:,:),daty,npixx,npixy,npixz)
          !$omp end parallel
          print 10,' taking fft of grid',npixx,'vx'
          !$omp parallel
          call fft3df(datgrid(2,:,:,:),datx,npixx,npixy,npixz)
          call power3d(datx,npixx,npixy,npixz,pk,xk,nk,numk,ptot,compensate,do_average,fty=daty,ftz=dati)
          !$omp end parallel
       elseif (i >= 6) then
!
!--these are Kinetic energy spectra, with various powers of rho
!  rho**0.5 gives the mass-weighted Kinetic energy (i.e., the usual SPH one)
!  rho**0.33 gives a 5/3 slope in the spectrum for compressible turbulence
!
          if (use_rhopow_from_particles) then
!
!--here we use rho defined on the particles to get the rho**power
!  by re-interpolating from the particles to the grid
!
             do j=1,npart
                vxyzu(1,j) = rho(j)**(rhopow - rhopowprev)*vxyzu(1,j)
                vxyzu(2,j) = rho(j)**(rhopow - rhopowprev)*vxyzu(2,j)
                vxyzu(3,j) = rho(j)**(rhopow - rhopowprev)*vxyzu(3,j)
             enddo
             call interpolate3D(xyzh,weight,massoftype(1),vxyzu,npart,xmin,ymin,zmin,datgrid, &
                                npixx,npixy,npixz,pixwidth,.true.,0)

             print 10,' taking fft of '//trim(variable)//'_z'
             !$omp parallel
             call fft3df(datgrid(4,:,:,:),dati,npixx,npixy,npixz)
             !$omp end parallel
             print 10,' taking fft of '//trim(variable)//'_y'
             !$omp parallel
             call fft3df(datgrid(3,:,:,:),daty,npixx,npixy,npixz)
             !$omp end parallel
             print 10,' taking fft of '//trim(variable)//'_x'
             !$omp parallel
             call fft3df(datgrid(2,:,:,:),datx,npixx,npixy,npixz)
             call power3d(datx,npixx,npixy,npixz,pk,xk,nk,numk,ptot,compensate,do_average,fty=daty,ftz=dati)
             !$omp end parallel
             rhopowprev = rhopow
          else
!
!--here we use rho defined on the grid to get the rho**power
!  this one seems to give better results, so it is the default
!
             print 10,' taking fft of '//trim(variable)//'_z'
             !$omp parallel
             call fft3df(datgrid(4,:,:,:)*datgrid(1,:,:,:)**rhopow,dati,npixx,npixy,npixz)
             !$omp end parallel
             print 10,' taking fft of '//trim(variable)//'_y'
             !$omp parallel
             call fft3df(datgrid(3,:,:,:)*datgrid(1,:,:,:)**rhopow,daty,npixx,npixy,npixz)
             !$omp end parallel
             print 10,' taking fft of '//trim(variable)//'_x'
             !$omp parallel
             call fft3df(datgrid(2,:,:,:)*datgrid(1,:,:,:)**rhopow,datx,npixx,npixy,npixz)
             call power3d(datx,npixx,npixy,npixz,pk,xk,nk,numk,ptot,compensate,do_average,fty=daty,ftz=dati)
             !$omp end parallel
             rhopowprev = rhopow
          endif
       endif

       print*,'In fourier space: '//trim(string)//' = ',ptot
       if (i==1) call open_power(fileout,dumpfile(1:mlabel),numk,xk,nk,n_out)
       !print*,'writing to ',trim(fileout)
       call write_power(numk,0.,pk,ptot,variable)
       call write_power_ascii(power_unit,trim(dumpfile),trim(variable),numk,xk,pk,ptot,comp=comp)
    enddo
 else
!
!--define output filename
!
    write(fileout,"(a,i3.3,a)") trim(dumpfile)//'_part',npixx,'k.pow'
!
!--allocate memory for k values
!
    print "(a,i4,a)",' (allocating memory for ',npixx,'^3 k values)'
    allocate(dati(npixx,npixy,npixz))
!
!--get the power spectrum for various combinations of rho and v
!  Kinetic energy power spectra (i > 4) require extra storage so these are done last
!
    do i=1,1
       select case(i)
       case(1)
          string = 'int(rho^2)dV'
          variable = 'rho'
          print*,'taking sph-fft of '//trim(variable)
          !call sphpow3D(rho,xyzh,rho,massoftype(1),npart,pk,xk,nk,numk,ptot,npixx,npixy,npixz,dxbound,.false.)
!          call sphpow3D(rho,xyzh,rho,massoftype(1),npart,pk,npixx,npixy,npixz,dxbound)
       case(2)
          variable = 'vx'
          string = 'Ekin_x'
          print*,'taking sph-fft of '//trim(variable)
          call sphfft3D(vxyzu(1,:),xyzh,rho,massoftype(1),npart,dati,npixx,npixy,npixz,dxbound)
       case(3)
          variable = 'vy'
          string = 'Ekin_y'
          print*,'taking sph-fft of '//trim(variable)
          call sphfft3D(vxyzu(2,:),xyzh,rho,massoftype(1),npart,dati,npixx,npixy,npixz,dxbound)
       case(4)
          variable = 'vz'
          string = 'Ekin_z'
          print*,'taking sph-fft of '//trim(variable)
          call sphfft3D(vxyzu(3,:),xyzh,rho,massoftype(1),npart,dati,npixx,npixy,npixz,dxbound)
       case(5)
          variable = 'Ekin'
          string = 'Ekin(vol)'
       case(6)
          rhopow = 1./3.
          variable = 'rho0.33ekin'
          string = 'Ekin((rho1/3 * v)^2)'
       case(7)
          rhopow = 0.5
          variable = 'rho0.50ekin'
          string = 'Ekin(rho)'
       end select
!       !$omp parallel
!       call power3d(dati,npixx,npixy,npixz,pk,xk,nk,numk,ptot,compensate,do_average)
!       !$omp end parallel

       print*,'In fourier space: '//trim(string)//' = ',ptot

       print*,'calculating power of '//trim(variable)
       dxmin = dxbound/256.
       call power_part(dxbound,dxmin,numk,xk,pk,xyzh,rho,rho,massoftype,npart,0)

       if (i==1) call open_power(fileout,dumpfile(1:mlabel),numk,xk,nk,n_out)
       print*,'writing to ',trim(fileout)
       call write_power(numk,0.,pk,ptot,variable)
       call write_power_ascii(power_unit,trim(dumpfile),trim(variable),numk,xk,pk,ptot,comp=comp)
       !close(power_unit)
    enddo
 endif

 if (allocated(datgrid)) deallocate(datgrid)
 if (allocated(dati)) deallocate(dati)
 if (allocated(datx)) deallocate(datx)
 if (allocated(daty)) deallocate(daty)
 if (allocated(pk)) deallocate(pk)
 if (allocated(nk)) deallocate(nk)
 if (allocated(xk)) deallocate(xk)
!---------------------------------------------

 print "(/,a,/)",' Phantom2power: we thank you for your kind custom'

contains
subroutine write_analysis_options(xmin,xmax,filename,iunit)
 use infile_utils, only:write_inopt
 real,             intent(in) :: xmin(:),xmax(:)
 character(len=*), intent(in) :: filename
 integer,          intent(in) :: iunit
 integer :: ierr,i

 open(unit=iunit,file=filename,status='replace',iostat=ierr)
 if (ierr /= 0) then
    print*,' ERROR writing to '//trim(filename)
    return
 endif
 print "(a)",' WRITING PHANTOM2POWER OPTIONS TO '//trim(filename)
 print "(a)",' RERUN PHANTOM2POWER AFTER EDITING THIS FILE'
 write(iunit,"(a)") '# box dimensions'
 do i=1,size(xmin)
    call write_inopt(xmin(i),labelx(i)//'min',' min boundary',iunit)
    call write_inopt(xmax(i),labelx(i)//'max',' max boundary',iunit)
 enddo
 close(iunit)

end subroutine write_analysis_options

subroutine read_analysis_options(xmin,dxmax,filename,iunit,ierr)
 use infile_utils, only:read_inopt,inopts,close_db,open_db_from_file
 real,             intent(out) :: xmin(:),dxmax(:)
 character(len=*), intent(in)  :: filename
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)
 real, dimension(size(xmin)) :: xmax
 integer :: i,nerr

 nerr = 0
 xmin = 0.
 xmax = 0.
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr == 0) then
    do i=1,3
       call read_inopt(xmin(i),labelx(i)//'min',db,errcount=nerr)
       call read_inopt(xmax(i),labelx(i)//'max',db,errcount=nerr)
    enddo
 endif
 if (nerr > 0) ierr = 1
 dxmax = xmax - xmin
 call close_db(db)

end subroutine read_analysis_options

subroutine print_bounds(xmin,xmax)
 real, intent(in) :: xmin(:),xmax(:)
 character(len=1), parameter :: labelx(3) = (/'x','y','z'/)
 integer :: i

 do i=1,size(xmin)
    print "(2(1x,a,1pg10.3))",labelx(i)//'min = ',xmin(i),&
                               labelx(i)//'max = ',xmax(i)
 enddo

end subroutine print_bounds
end program phantom2power
