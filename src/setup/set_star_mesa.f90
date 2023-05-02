!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setstar_mesa
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
 use fileutils, only:get_nlines,get_ncolumns,string_delete,lcase
 use datafiles, only:find_phantom_datafile
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 integer                                    :: lines,rows,i,ncols,nheaderlines,iu
 character(len=*), intent(in)               :: filepath
 logical, intent(in), optional              :: cgsunits
 integer, intent(out)                       :: ierr
 character(len=10000)                       :: dumc
 character(len=120)                         :: fullfilepath
 character(len=24),allocatable              :: header(:),dum(:)
 logical                                    :: iexist,usecgs,ismesafile,got_column
 real,allocatable,dimension(:,:)            :: dat
 real, intent(in)                           :: X_in,Z_in
 real,allocatable,dimension(:),intent(out)  :: rho,r,pres,m,ene,temp,Xfrac,Yfrac
 real, intent(out)                          :: Mstar

 rows = 0
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
 open(newunit=iu,file=fullfilepath,status='old')
 call get_ncolumns(iu,ncols,nheaderlines)
 if (nheaderlines == 6) then ! Assume file is a MESA profile, and so it has 6 header lines, and (row=3, col=2) = number of zones
    read(iu,'()')
    read(iu,'()')
    read(iu,*) lines,lines
    read(iu,'()')
    read(iu,'()')
    ismesafile = .true.
 else
    ismesafile = .false.
    lines = lines - nheaderlines
    do i = 1,nheaderlines-1
       read(iu,'()')
    enddo
 endif
 if (lines <= 0) then ! file not found
    ierr = 1
    return
 endif

 read(iu,'(a)') dumc! counting rows
 call string_delete(dumc,'[')
 call string_delete(dumc,']')
 allocate(dum(500)) ; dum = 'aaa'
 read(dumc,*,end=101) dum
101 continue
 do i = 1,500
    if (dum(i)=='aaa') then
       rows = i-1
       exit
    endif
 enddo
 allocate(header(rows),dat(lines,rows))
 header(1:rows) = dum(1:rows)
 deallocate(dum)
 do i = 1,lines
    read(iu,*) dat(lines-i+1,1:rows)
 enddo

 allocate(m(lines),r(lines),pres(lines),rho(lines),ene(lines), &
             temp(lines),Xfrac(lines),Yfrac(lines))

 close(iu)
 ! Set mass fractions to fixed inputs if not in file
 Xfrac = X_in
 Yfrac = 1. - X_in - Z_in
 do i = 1,rows
    if (header(i)(1:1) == '#' .and. .not. trim(lcase(header(i)))=='#mass') then
       print '("Detected wrong header entry : ",A," in file ",A)',trim(lcase(header(i))),trim(fullfilepath)
       ierr = 2
       return
    endif
    got_column = .true.
    select case(trim(lcase(header(i))))
    case('mass_grams')
       m = dat(1:lines,i)
    case('mass','#mass')
       m = dat(1:lines,i)
       if (ismesafile) m = m * solarm  ! If reading MESA profile, 'mass' is in units of Msun
    case('rho','density')
       rho = dat(1:lines,i)
    case('logrho')
       rho = 10**(dat(1:lines,i))
    case('energy','e_int','e_internal')
       ene = dat(1:lines,i)
    case('radius_cm')
       r = dat(1:lines,i)
    case('radius')
       r = dat(1:lines,i)
       if (ismesafile) r = r * solarr
    case('logr')
       r = (10**dat(1:lines,i)) * solarr
    case('pressure')
       pres = dat(1:lines,i)
    case('temperature')
       temp = dat(1:lines,i)
    case('x_mass_fraction_h','xfrac')
       Xfrac = dat(1:lines,i)
    case('y_mass_fraction_he','yfrac')
       Yfrac = dat(1:lines,i)
    case default
       got_column = .false.
    end select
    if (got_column) print "(1x,i0,': ',a)",i,trim(header(i))
 enddo
 print "(a)"

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
 character(len=200)              :: headers
 integer                         :: i,noptionalcols,j,iu
 real, allocatable               :: optionalcols(:,:)
 character(len=*), parameter     :: fmtstring = "(5(es13.6,2x),es13.6)"

 headers = '[    Mass   ]  [  Pressure ]  [Temperature]  [   Radius  ]  [  Density  ]  [   E_int   ]'

 ! Add optional columns
 allocate(optionalcols(size(r),10))
 noptionalcols = 0
 if (present(Xfrac)) then
    noptionalcols = noptionalcols + 1
    headers = trim(headers) // '  [   Xfrac   ]'
    optionalcols(:,noptionalcols) = Xfrac
 endif
 if (present(Yfrac)) then
    noptionalcols = noptionalcols + 1
    headers = trim(headers) // '  [   Yfrac   ]'
    optionalcols(:,noptionalcols) = Yfrac
 endif
 if (present(mu)) then
    noptionalcols = noptionalcols + 1
    headers = trim(headers) // '  [    mu     ]'
    optionalcols(:,noptionalcols) = mu
 endif
 if (present(csound)) then
    noptionalcols = noptionalcols + 1
    headers = trim(headers) // '  [Sound speed]'
    optionalcols(:,noptionalcols) = csound
 endif

 open(newunit=iu, file = outputpath, status = 'replace')
 write(iu,'(a)') headers

 do i=1,size(r)
    if (noptionalcols <= 0) then
       write(iu,fmtstring) m(i),pres(i),temp(i),r(i),rho(i),ene(i)
    else
       write(iu,fmtstring,advance="no") m(i),pres(i),temp(i),r(i),rho(i),ene(i)
       do j=1,noptionalcols
          if (j==noptionalcols) then
             write(iu,'(2x,es13.6)') optionalcols(i,j)
          else
             write(iu,'(2x,es13.6)',advance="no") optionalcols(i,j)
          endif
       enddo
    endif
 enddo
 close(iu)

end subroutine write_mesa

end module setstar_mesa
