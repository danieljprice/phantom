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
!  Analysis routine to produce accurate volume-weighted
!  Probability Distribution Functions from SPH particle data
!  by interpolation to an adaptive mesh
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    binspacing --  bin width in ln(rho) for PDF
!    rhologmax  --  max ln(rho) for PDF
!    rhologmin  --  min ln(rho) for PDF
!
!  DEPENDENCIES: adaptivemesh, boundary, dim, infile_utils,
!    interpolations3D_amr, part, pdfs
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'PDF on AMR mesh'
 public :: do_analysis
 real, public :: rhologmax = 15.
 real, public :: rhologmin = -15.
 real, public :: Bsqlogmax = 15.
 real, public :: Bsqlogmin = -20.
 character(len=1), parameter :: labelx(3) = (/'x','y','z'/)
 real, public :: binspacing = 0.1

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use adaptivemesh,         only:build_mesh,nsub,ndim,ifirstlevel
 use boundary,             only:xmin,ymin,zmin,dxbound,dybound,dzbound
 use part,                 only:hfact,rhoh,mhd,Bxyz,isdead_or_accreted
 use interpolations3D_amr, only:interpolate3D_amr
 use pdfs,                 only:pdf_write
 use dim,                  only:periodic,tagline
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
!--
 real, allocatable :: datgrid(:,:,:)
 real, allocatable :: xval(:), pdf(:)
 real :: xminp(3),dxmaxp(3),xminpart(3),xmaxpart(3)
 real               :: weight,rhoi,hi
 real               :: rmsv,rhomin,rhomean,rhomax,smean,totvol,totmass
 real               :: rhomin_box,rhomean_box,rhomax_box
 real               :: rhomeanmw,rhovarmw,svarmw,smeanmw,bvalmw,rmsvmw
 real               :: rhovar,svar,totprob,bval,t2,t1
 integer            :: maxp,maxmhd,np,ninbox
 integer            :: nmesh,i,j,nbins,ibin,ierr,ilendat
 character(len=120) :: tagline1,fileout
 character(len=20)  :: fmtstring,fmt
 character(len=20), parameter :: filename = 'analysis.in'
 logical            :: read_input_file

!
!--build adaptive mesh structure
!
 xminp(1)  = xmin
 xminp(2)  = ymin
 xminp(3)  = zmin
 dxmaxp(1) = dxbound
 dxmaxp(2) = dybound
 dxmaxp(3) = dzbound
 !call print_bounds(xminp,xminp+dxmaxp)

 read_input_file = .false.
 call read_analysis_options(xminp,dxmaxp,filename,iunit,ierr)
 if (ierr == 0) read_input_file = .true.

 if (.not.read_input_file .and. periodic) then
    call write_analysis_options(xminp,xminp+dxmaxp,filename,iunit+1)
    stop
 endif

 rhomax  = 0.
 rhomin  = huge(rhomin)
 if (npart < 1) then
    print*, 'ERROR: no SPH particles... quitting '
    return
 endif
 xminpart(:) = xyzh(1:3,1)
 xmaxpart(:) = xyzh(1:3,1)
 np = 0
 ninbox = 0
 rhomax_box = 0.
 rhomin_box = huge(rhomin_box)
 rhomean_box = 0.
 do i=1,npart
    hi     = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       rhoi   = rhoh(hi,particlemass)
       rhomax = max(rhomax,rhoi)
       rhomin = min(rhomin,rhoi)
       rhomean = rhomean + rhoi
       do j=1,3
          xminpart(j) = min(xyzh(j,i),xminpart(j))
          xmaxpart(j) = max(xyzh(j,i),xmaxpart(j))
       enddo
       if (.not.periodic) then
          if (in_bounds(xyzh(1:3,i),xminp,dxmaxp)) then
             ninbox = ninbox + 1
             rhomax_box = max(rhomax_box,rhoi)
             rhomean_box = rhomean_box + rhoi
             rhomin_box = min(rhomin_box,rhoi)
          endif
       endif
       np = np + 1
    endif
 enddo
 if (np > 0) rhomean = rhomean/real(np)
 if (ninbox > 0) rhomean_box = rhomean_box/real(ninbox)
 if (.not.periodic) then
    print *,' MIN, MAX of particle distribution is:'
    call print_bounds(xminpart,xmaxpart)
    if (.not.read_input_file) then
       call write_analysis_options(xminpart,xmaxpart,filename,iunit+1)
       stop
    endif
 endif

 print "(72('-'))"
 fmt = "(1x,2(a,1pg10.3))"
 print fmt,'On parts: Total Mass =',particlemass*npart,' mean dens (by mass) =',rhomean
 print fmt,'On parts: max. dens  =',rhomax,            ' min. dens =',rhomin
 if (.not.periodic) then
    print fmt,'On parts: Total Mass in box =',particlemass*ninbox,' mean dens in box =',rhomean_box
    print fmt,'On parts: max. dens  in box =',rhomax_box,         ' min. dens in box =',rhomin_box
    print "(3(a,i10),a)",' Got ',ninbox, ' of ',npart,' particles in box (not counting ',npart-np,' accreted)'
 endif
 print fmt,'Total mass accreted =',(npart-np)*particlemass
 print "(72('-'))"

 call cpu_time(t1)
 call build_mesh(xyzh,npart,nmesh,xminp,dxmaxp)
 call cpu_time(t2)
 print*,'built mesh in ',t2-t1,' cpu-s'

 maxp = size(xyzh, 2)
 maxmhd = size(Bxyz, 2)
 ilendat = 4*(maxmhd/maxp)

! allocate memory for the grid data
!
 allocate(datgrid(4+ilendat,nsub**ndim,nmesh))
!
!--interpolate to the 3D grid
!
 weight = 1.d0/hfact**3
 print*,'using hfact = ',hfact,' weight = ',weight

 if (maxmhd > 0) then
    call interpolate3D_amr(xyzh,weight,particlemass,vxyzu,npart, &
                           xminp,datgrid,nmesh,dxmaxp,.true.,ilendat,Bxyz)
 else
    call interpolate3D_amr(xyzh,weight,particlemass,vxyzu,npart, &
                           xminp,datgrid,nmesh,dxmaxp,.true.,0)
 endif

 rhomin  = huge(rhomin)
 rhomean = 0.
 rhomax  = 0.
 rmsv    = 0.
 smean   = 0.
 totvol  = 0.
 totmass = 0.
 call get_rhomach(1,ifirstlevel,datgrid,rmsv,rhomin,rhomean,rhomax,smean,totvol,totmass)
 rmsv    = sqrt(rmsv)
 ! global integrals assume unit box, scale to actual box size
 totvol  = totvol*product(dxmaxp)
 totmass = totmass*product(dxmaxp)

 print "(72('-'))"
 print fmt,'On grid: Total Mass =',totmass,' mean dens (by vol) =',rhomean
 print fmt,'On grid: max. dens  =',rhomax, ' min. dens =',rhomin
 print fmt,'On grid: rms v      =',rmsv,   ' mean(lnrho) =',smean
 print fmt,'On grid: volume     =',totvol
 print "(72('-'))"

 svar    = 0.
 rhovar  = 0.
 call get_variance(1,ifirstlevel,datgrid,rhomean,smean,rhovar,svar)
 bval    = sqrt(svar)/rmsv

 print*,' var(rho) = ',rhovar,' var(lnrho) = ',svar
 !
 !--fill in the blanks
 !
 rhomeanmw = 0.
 rhovarmw  = 0.
 svarmw    = 0.
 smeanmw   = 0.
 bvalmw    = 0.
 rmsvmw    = 0.
 !
 !--write line to rhomach file
 !
 write(fileout,"(a)") 'rhomach_'//trim(dumpfile)//'_amrgrid.out'
 print "(a)",' writing to '//trim(fileout)

 open(unit=iunit,file=trim(fileout),status='replace',form='formatted')
 write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) 17
 write(iunit,fmtstring) time,rhomean,rhomeanmw,rhovar,rhovarmw,sqrt(rhovar),sqrt(rhovarmw),&
                           rmsv,rmsvmw,bval,bvalmw,smean,smeanmw,svar,svarmw,sqrt(svar),sqrt(svarmw)
 close(unit=iunit)
!
! calculate ln(rho/rho_0) PDF
!

 print "(/,a)",' Computing PDF of ln(density)...'

! allocate memory for PDF calculation
 nbins = nint((rhologmax - rhologmin)/binspacing)
 print "(a,i3,a)",' (allocating memory for ',nbins,' PDF bins)'
 if (.not.allocated(xval)) allocate(xval(nbins),pdf(nbins))

! setup PDF bins
 print "(a,es10.3)",' bin width = ',binspacing
 print "(2(a,1pg10.3))",' lnrho min = ',rhologmin,' max = ',rhologmax
 do ibin=1,nbins
    xval(ibin) = rhologmin + (ibin-1)*binspacing
 enddo
 pdf(:) = 0.
 call get_pdf_lnrho(1,ifirstlevel,datgrid,rhologmin,rhologmax,binspacing,nbins,pdf)

! get total area under PDF
 totprob = 0.
 do ibin=1,nbins
    totprob = totprob + binspacing*pdf(ibin)
 enddo
 print*,'normalisation factor = ',totprob,' should be ',binspacing
 pdf(:) = pdf(:)/totprob

! write PDF to file
 xval = exp(xval)
 tagline1 = 'Phantomanalysis: '//trim(analysistype)//', part of '//trim(tagline)
 call pdf_write(nbins,xval,pdf,'lnrho',.true.,trim(dumpfile),trim(tagline1))

 if (allocated(xval)) deallocate(xval)
 if (allocated(pdf)) deallocate(pdf)

!
! calculate log(B^2) PDF
!
 if (maxmhd > 0) then
    print "(/,a)",' Computing PDF of magnetic energy ...'
    nbins = nint((Bsqlogmax - Bsqlogmin)/binspacing)
    print "(a,i3,a)",' (allocating memory for ',nbins,' PDF bins)'
    if (.not.allocated(xval)) allocate(xval(nbins),pdf(nbins))

    print "(a,1pe10.3)",' bin width = ',binspacing
    do ibin=1,nbins
       xval(ibin) = Bsqlogmin + (ibin-1)*binspacing
    enddo
    pdf(:) = 0.
    call get_pdf_logBsq(1,ifirstlevel,datgrid,Bsqlogmin,Bsqlogmax,binspacing,nbins,pdf)

    xval = exp(xval)
    call pdf_write(nbins,xval,pdf,'logBsq',.true.,trim(dumpfile),trim(tagline1))
 endif
!
! cleanup memory
!
 if (allocated(xval)) deallocate(xval)
 if (allocated(pdf)) deallocate(pdf)
 if (allocated(datgrid)) deallocate(datgrid)

end subroutine do_analysis

subroutine print_bounds(xmin,xmax)
 real, intent(in) :: xmin(:),xmax(:)
 character(len=1), parameter :: labelx(3) = (/'x','y','z'/)
 integer :: i

 do i=1,size(xmin)
    print "(2(1x,a,1pg10.3))",labelx(i)//'min = ',xmin(i),&
                              labelx(i)//'max = ',xmax(i)
 enddo

end subroutine print_bounds

!
! function to determine whether a given particle is
! inside the box or not
!
logical function in_bounds(xyzi,xmin,dx)
 real, intent(in) :: xyzi(:),xmin(:),dx(:)
 real :: xmax(size(xmin))
 logical :: in_box(size(xmin))

 xmax = xmin + dx
 in_box = .false.
 where (xyzi > xmin .and. xyzi < xmax)
    in_box = .true.
 end where
 in_bounds = all(in_box)

end function in_bounds

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
 print "(a)",' WRITING ANALYSIS OPTIONS TO '//trim(filename)
 print "(a)",' RERUN THE ANALYSIS AFTER EDITING THIS FILE'
 write(iunit,"(a)") '# box dimensions'
 do i=1,size(xmin)
    call write_inopt(xmin(i),labelx(i)//'min',' min boundary',iunit)
    call write_inopt(xmax(i),labelx(i)//'max',' max boundary',iunit)
 enddo
 write(iunit,"(a)") '# analysis options'
 call write_inopt(rhologmin,'rhologmin',' min ln(rho) for PDF',iunit)
 call write_inopt(rhologmax,'rhologmax',' max ln(rho) for PDF',iunit)
 call write_inopt(binspacing,'binspacing',' bin width in ln(rho) for PDF',iunit)
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
    call read_inopt(rhologmin,'rhologmin',db,errcount=nerr)
    call read_inopt(rhologmax,'rhologmax',db,errcount=nerr)
    call read_inopt(binspacing,'binspacing',db,errcount=nerr)
 endif
 dxmax = xmax - xmin
 call close_db(db)

end subroutine read_analysis_options

recursive subroutine get_rhomach(imesh,level,datgrid,rmsv,rhomin,rhomean,rhomax,smean,totvol,totmass)
 use adaptivemesh, only:nsub,ndim,gridnodes
 real,    intent(inout) :: rmsv,rhomin,rhomean,rhomax,smean,totvol,totmass
 integer, intent(in)    :: imesh,level
 real,    intent(in)    :: datgrid(:,:,:)
 real    :: dn,weightl,rhoi,si,vx2,vy2,vz2
 integer :: isubmesh,icell
 real, parameter :: rhologmin = -15.
!
!--now calculate mean density from the gridded data
!
! weightl = 1./((nsub**ndim)**(level-ifirstlevel))
 dn      = 1./(nsub**ndim)   ! get weight in two steps
 weightl = dn**level         ! to avoid problems if level >= 12

 do icell=1,nsub**ndim
    isubmesh = gridnodes(icell,imesh) !grid(imesh)%daughter(icell)
    if (isubmesh > 0) then
       call get_rhomach(isubmesh,level+1,datgrid,rmsv,rhomin,rhomean,rhomax,smean,totvol,totmass)
    else
       rhoi    = datgrid(1,icell,imesh)

       if (rhoi > 0.) then
          si   = log(rhoi)
       else
          si   = rhologmin
       endif
       vx2     = datgrid(2,icell,imesh)**2
       vy2     = datgrid(3,icell,imesh)**2
       vz2     = datgrid(4,icell,imesh)**2
       rmsv    = rmsv + weightl*(vx2 + vy2 + vz2)
       rhomin  = min(rhomin,rhoi)
       rhomax  = max(rhomax,rhoi)
       rhomean = rhomean + weightl*rhoi
       smean   = smean + weightl*si
       totvol  = totvol + weightl
       totmass = totmass + weightl*rhoi
       !print*,'mesh = ',imesh,' cell = ',icell,' dat = ',datgrid(2,icell,imesh)
    endif
 enddo
end subroutine get_rhomach

recursive subroutine get_variance(imesh,level,datgrid,rhomean,smean,rhovar,svar)
 use adaptivemesh, only:nsub,ndim,gridnodes
 real,    intent(in)    :: rhomean,smean
 real,    intent(inout) :: rhovar,svar
 integer, intent(in)    :: imesh,level
 real,    intent(in)    :: datgrid(:,:,:)
 real    :: dn,weightl,rhoi,si
 integer :: isubmesh,icell
!
!--calculate variance in rho and ln(rho) from the gridded data
!
! weightl = 1./((nsub**ndim)**(level-ifirstlevel))
 dn      = 1./(nsub**ndim)    ! get weight in two steps
 weightl = dn**level          ! to avoid problems if level >= 12

 do icell=1,nsub**ndim
    isubmesh = gridnodes(icell,imesh) !grid(imesh)%daughter(icell)
    if (isubmesh > 0) then
       call get_variance(isubmesh,level+1,datgrid,rhomean,smean,rhovar,svar)
    else
       rhoi   = datgrid(1,icell,imesh)
       if (rhoi > 0.) then
          si  = log(rhoi)
       else
          si  = rhologmin
       endif
       rhovar = rhovar + weightl*(rhoi - rhomean)**2
       svar   = svar   + weightl*(si - smean)**2
       !print*,'mesh = ',imesh,' cell = ',icell,' dat = ',datgrid(2,icell,imesh)
    endif
 enddo
end subroutine get_variance


recursive subroutine get_pdf_lnrho(imesh,level,datgrid,smin,smax,ds,nbins,pdf)
 use adaptivemesh, only:nsub,ndim,ifirstlevel,gridnodes
 real,    intent(in)    :: smin,smax,ds
 integer, intent(in)    :: imesh,level,nbins
 real,    intent(in)    :: datgrid(:,:,:)
 real,    intent(inout) :: pdf(:)
 real    :: dn,weightl,rhoi,si
 integer :: isubmesh,icell,ibin
!
!--calculate variance in rho and ln(rho) from the gridded data
!
! weightl = 1./((nsub**ndim)**(level-ifirstlevel))
 dn      = 1./(nsub**ndim)  ! get weight in two steps
 weightl = dn**level        ! to avoid problems if level >= 12

 over_cells: do icell=1,nsub**ndim
    isubmesh = gridnodes(icell,imesh) !grid(imesh)%daughter(icell)
    if (isubmesh > 0) then
       call get_pdf_lnrho(isubmesh,level+1,datgrid,smin,smax,ds,nbins,pdf)
    else
       rhoi   = datgrid(1,icell,imesh)
       if (rhoi > 0.) then
          si  = log(rhoi)
       else
          cycle over_cells
       endif
       ibin = int((si - smin)/ds) + 1
       if (ibin < 1) cycle over_cells !    ibin = 1
       if (ibin > nbins) cycle over_cells !ibin = nbins

       pdf(ibin) = pdf(ibin) + weightl
    endif
 enddo over_cells

end subroutine get_pdf_lnrho

recursive subroutine get_pdf_logBsq(imesh,level,datgrid,Bsqmin,Bsqmax,binspacing,nbins,pdf)
 use adaptivemesh, only:nsub,ndim,ifirstlevel,gridnodes
 implicit none
 real,    intent(in)    :: Bsqmin,Bsqmax,binspacing
 integer, intent(in)    :: imesh,level,nbins
 real, dimension(:,:,:), intent(in) :: datgrid
 real, dimension(:),  intent(inout) :: pdf
 real    :: dn,weightl,Bxi,Byi,Bzi,Bsqi,logBsqi
 integer :: isubmesh,icell,ibin

! weightl = 1./((nsub**ndim)**(level-ifirstlevel))
 dn      = 1./(nsub**ndim)  ! get weight in two steps
 weightl = dn**level        ! to avoid problems if level >= 12

 do icell=1,nsub**ndim
    isubmesh = gridnodes(icell,imesh) !grid(imesh)%daughter(icell)
    if (isubmesh > 0) then
       call get_pdf_logBsq(isubmesh,level+1,datgrid,Bsqmin,Bsqmax,binspacing,nbins,pdf)
    else
       Bxi  = datgrid(5,icell,imesh)
       Byi  = datgrid(6,icell,imesh)
       Bzi  = datgrid(7,icell,imesh)
       Bsqi = Bxi*Bxi + Byi*Byi + Bzi*Bzi
       if (Bsqi > 0.) then
          logBsqi  = log10(Bsqi)
       else
          logBsqi  = Bsqlogmin
       endif
       ibin = int((logBsqi - Bsqmin)/binspacing) + 1
       if (ibin < 1)     ibin = 1
       if (ibin > nbins) ibin = nbins

       pdf(ibin) = pdf(ibin) + weightl
    endif
 enddo

end subroutine get_pdf_logBsq

end module
