!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_mesa
!
! Utility routines for reading stellar profiles from the MESA code
!
! :References: Paxton et al. (2011), ApJS 192, 3
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: datafiles, fileutils, physcon, units
!
 implicit none

 public :: read_mesa,write_mesa

 private

contains
!-----------------------------------------------------------------------
!+
!  Read quantities from MESA profile or from profile in the format of
!  the P12 star (phantom/data/star_data_files/P12_Phantom_Profile.data)
!+
!-----------------------------------------------------------------------
subroutine read_mesa(filepath,rho,r,pres,m,ene,temp,X_in,Z_in,Xfrac,Yfrac,Mstar,ierr,cgsunits)
 use physcon,   only:solarm,solarr
 use fileutils, only:get_nlines,get_ncolumns,string_delete,lcase,read_column_labels
 use datafiles, only:find_phantom_datafile
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 character(len=*), intent(in)               :: filepath
 integer, intent(out)                       :: ierr
 real,    intent(in)                        :: X_in,Z_in
 real, allocatable,dimension(:),intent(out) :: rho,r,pres,m,ene,temp,Xfrac,Yfrac
 real, intent(out)                          :: Mstar
 logical, intent(in), optional              :: cgsunits
 integer                                    :: lines,i,ncols,nheaderlines,nlabels
 integer                                    :: idir,iu
 character(len=120)                         :: fullfilepath
 character(len=24),allocatable              :: header(:)
 logical                                    :: iexist,usecgs,ismesafile,got_column
 real,allocatable,dimension(:,:)            :: dat

 usecgs = .false.
 if (present(cgsunits)) usecgs = cgsunits
 !
 !--Get path name
 !
 ierr = 0
 fullfilepath = find_phantom_datafile(filepath,'star_data_files')
 inquire(file=trim(fullfilepath),exist=iexist)
 if (.not.iexist) then
    ierr = 1
    return
 endif
 lines = get_nlines(fullfilepath) ! total number of lines in file

 print "(1x,a)",trim(fullfilepath)
 open(newunit=iu,file=fullfilepath,status='old',iostat=ierr)
 if (ierr /= 0) then
    print "(a,/)",' ERROR opening file '//trim(fullfilepath)
    return
 endif

 call get_ncolumns(iu,ncols,nheaderlines)
 if (nheaderlines == 6) then ! Assume file is a MESA profile, and so it has 6 header lines, and (row=3, col=2) = number of zones
    read(iu,'()',iostat=ierr)
    read(iu,'()',iostat=ierr)
    read(iu,*,iostat=ierr) lines,lines
    read(iu,'()',iostat=ierr)
    read(iu,'()',iostat=ierr)
    if (ierr /= 0) then
       print "(a,/)",' ERROR reading MESA file header'
       return
    endif
    ismesafile = .true.
 else
    ismesafile = .false.
    lines = lines - nheaderlines
    do i = 1,nheaderlines-1
       read(iu,'()',iostat=ierr)
    enddo
    if (ierr /= 0) then
       print "(a,/)",' ERROR reading file header [not MESA format]'
       return
    endif
 endif
 if (lines <= 0) then ! file not found
    ierr = 1
    return
 endif

 ! extract column labels from the file header
 allocate(header(ncols),dat(lines,ncols))
 call read_column_labels(iu,nheaderlines,ncols,nlabels,header)
 if (nlabels /= ncols) print*,' WARNING: different number of labels compared to columns'

 allocate(m(lines))
 m = -1.
 allocate(r,pres,rho,ene,temp,Xfrac,Yfrac,source=m)

 over_directions: do idir=1,2   ! try backwards, then forwards
    if (idir==1) then
       ! read MESA file backwards, from surface to centre
       do i = 1,lines
          read(iu,*,iostat=ierr) dat(lines-i+1,1:ncols)
       enddo
    else
       ! read file forwards, from centre to surface
       do i = 1,lines
          read(iu,*,iostat=ierr) dat(i,1:ncols)
       enddo
    endif
    if (ierr /= 0) then
       print "(a,/)",' ERROR reading data from file: reached end of file?'
       return
    endif

    ! Set mass fractions to fixed inputs if not in file
    Xfrac = X_in
    Yfrac = 1. - X_in - Z_in
    do i = 1,ncols
       if (header(i)(1:1) == '#' .and. .not. trim(lcase(header(i)))=='#mass') then
          print '("Detected wrong header entry : ",a," in file ",a)',trim(lcase(header(i))),trim(fullfilepath)
          ierr = 2
          return
       endif
       got_column = .true.
       select case(trim(lcase(header(i))))
       case('mass_grams')
          m = dat(1:lines,i)
       case('mass','#mass','m')
          m = dat(1:lines,i)
          if (ismesafile .or. maxval(m) < 1.e-10*solarm) m = m * solarm  ! If reading MESA profile, 'mass' is in units of Msun
       case('logM','log_mass')
          m = 10**dat(1:lines,i)
          if (ismesafile .or. maxval(m) < 1.e-10*solarm) m = m * solarm  ! If reading MESA profile, 'mass' is in units of Msun
       case('rho','density')
          rho = dat(1:lines,i)
       case('logrho')
          rho = 10**(dat(1:lines,i))
       case('energy','e_int','e_internal')
          ene = dat(1:lines,i)
       case('loge')
          ene = 10**dat(1:lines,i)
       case('radius_cm')
          r = dat(1:lines,i)
       case('radius_km')
          r = dat(1:lines,i) * 1e5
       case('radius','r')
          r = dat(1:lines,i)
          if (ismesafile .or. maxval(r) < 1e-10*solarr) r = r * solarr
       case('logr')
          r = (10**dat(1:lines,i)) * solarr
       case('logr_cm')
          r = 10**dat(1:lines,i)
       case('pressure','p')
          pres = dat(1:lines,i)
       case('logp')
          pres = 10**dat(1:lines,i)
       case('temperature','t')
          temp = dat(1:lines,i)
       case('logt')
          temp = 10**dat(1:lines,i)
       case('x_mass_fraction_h','x','xfrac','h1')
          Xfrac = dat(1:lines,i)
       case('log_x')
          Xfrac = 10**dat(1:lines,i)
       case('y_mass_fraction_he','y','yfrac','he4')
          Yfrac = dat(1:lines,i)
       case('log_y')
          Yfrac = 10**dat(1:lines,i)
       case default
          got_column = .false.
       end select
       if (got_column .and. idir==1) print "(1x,i0,': ',a)",i,trim(header(i))
    enddo
    if (idir==1) print "(a)"

    ! quit the loop over directions if the radius increases
    if (idir==1 .and. r(2) > r(1)) exit over_directions

    ! otherwise rewind and re-skip header
    rewind(iu)
    do i=1,nheaderlines
       read(iu,*,iostat=ierr)
    enddo
 enddo over_directions
 close(iu)

 if(min(minval(pres),minval(rho))<0d0)ierr = 1

 if (ierr /= 0) then
    print "(a,/)",' ERROR reading MESA file [missing required columns]'
    return
 endif

 if (.not. usecgs) then
    m = m / umass
    r = r / udist
    pres = pres / unit_pressure
    rho = rho / unit_density
    ene = ene / unit_ergg
 endif

 Mstar = m(lines)

end subroutine read_mesa

!----------------------------------------------------------------
!+
!  Write stellar profile in format readable by read_mesa;
!  used in star setup to write softened stellar profile.
!+
!----------------------------------------------------------------
subroutine write_mesa(outputpath,m,pres,temp,r,rho,ene,Xfrac,Yfrac,csound,mu)
 real, intent(in)                :: m(:),rho(:),pres(:),r(:),ene(:),temp(:)
 real, intent(in), optional      :: Xfrac(:),Yfrac(:),csound(:),mu(:)
 character(len=120), intent(in)  :: outputpath
 character(len=200)              :: headers(100)
 integer                         :: i,ncols,noptionalcols,j,iu
 real, allocatable               :: optionalcols(:,:)
 character(len=*), parameter     :: fmtstring = "(5(es24.16e3,2x),es24.16e3)"

 ncols = 6
 headers(1:ncols) = (/'       Mass','   Pressure','Temperature','     Radius','    Density','       Eint'/)

 ! Add optional columns
 noptionalcols = 0
 allocate(optionalcols(size(r),10))
 if (present(Xfrac)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = '      Xfrac'
    optionalcols(:,noptionalcols) = Xfrac
 endif
 if (present(Yfrac)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = '      Yfrac'
    optionalcols(:,noptionalcols) = Yfrac
 endif
 if (present(mu)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = '        mu'
    optionalcols(:,noptionalcols) = mu
 endif
 if (present(csound)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = 'Sound speed'
    optionalcols(:,noptionalcols) = csound
 endif

 open(newunit=iu,file=outputpath,status='replace')
 do i = 1,noptionalcols+ncols-1
    write(iu,'(a24,2x)',advance="no") trim(headers(i))
 enddo
 write(iu,'(a24)') trim(headers(noptionalcols+ncols))

 do i=1,size(r)
    if (noptionalcols <= 0) then
       write(iu,fmtstring) m(i),pres(i),temp(i),r(i),rho(i),ene(i)
    else
       write(iu,fmtstring,advance="no") m(i),pres(i),temp(i),r(i),rho(i),ene(i)
       do j=1,noptionalcols
          if (j==noptionalcols) then
             write(iu,'(2x,es24.16e3)') optionalcols(i,j)
          else
             write(iu,'(2x,es24.16e3)',advance="no") optionalcols(i,j)
          endif
       enddo
    endif
 enddo
 close(iu)

end subroutine write_mesa

end module readwrite_mesa
