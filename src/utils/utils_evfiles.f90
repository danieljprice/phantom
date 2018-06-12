!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: evutils
!
!  DESCRIPTION:
!   Contains supplementary subroutines useful for dealing with .ev files
!   -> currently used by phantomevcompare and ev2mdot
!
!  REFERENCES:
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    Concise -- write columns that exist in all files (F: write columns that exist in any file)
!
!  DEPENDENCIES: infile_utils
!+
!--------------------------------------------------------------------------
module evutils
 implicit none
 !
 ! Subroutines
 !
 public :: read_evfile
 public :: get_column_labels_from_ev,write_columns_to_file
 public :: read_evin_file,read_evin_filenames,write_evin_file
 public :: find_column,verboseentry
 integer, parameter, public :: max_columns = 150
 integer, parameter :: inumev = max_columns
 !
 ! logical unit number for read/write operations
 !
 integer, parameter :: ianalysis = 25

 private

contains
!----------------------------------------------------------------
!+
!  Read the column labels from the .ev file
!  this is a wrapper to the subroutine below
!+
!----------------------------------------------------------------
subroutine get_column_labels_from_ev(evfile,labels,numcol,ierr)
 integer,             intent(out) :: numcol,ierr
 character(len=*   ), intent(in)  :: evfile
 character(len=*   ), intent(out) :: labels(:)
 logical                          :: iexist
 character(len=2048)              :: line

 ierr = 0
 inquire(file=trim(evfile),exist=iexist)
 if (iexist) then
    open(unit=ianalysis,file=trim(evfile))
    read(ianalysis,'(a)') line
    close(ianalysis)
    call get_column_labels(line,labels,numcol)
 else
    ierr   = 1
    numcol = 0
 endif

end subroutine get_column_labels_from_ev

!----------------------------------------------------------------
!+
!  find the column labels from the corresponding header line
!  Note: This formatting is specific to the formating in evwrite.f90
!  i.e. [01  blah] [02  blah] [03  blah]
!+
!----------------------------------------------------------------
subroutine get_column_labels(line,labels,numcol)
 character(len=*), intent(inout) :: line
 character(len=*), intent(out)   :: labels(:)
 integer,          intent(out)   :: numcol
 integer                         :: i,iopen,iclose

 i = 0
 iopen = 1 ! to get into the loop
 do while ( iopen > 0 .and. i < size(labels))
    iopen  = index(line,'[')
    iclose = index(line,']')
    i = i + 1
    labels(i) = trim(adjustl(line(iopen +3:iclose-1)))
    line      = line(iclose+1:)
 enddo
 numcol = i - 1

end subroutine get_column_labels

!----------------------------------------------------------------
!+
!  find the column number matching a particular label
!+
!----------------------------------------------------------------
integer function find_column(labels,tag)
 character(len=*), intent(in) :: labels(:)
 character(len=*), intent(in) :: tag
 integer :: i

 find_column = 0
 do i=1,size(labels)  ! find first one that matches
    if (trim(adjustl(labels(i)))==trim(adjustl(tag))) then
       find_column = i
       exit
    endif
 enddo

end function find_column

!----------------------------------------------------------------
!+
!  subroutine to read entire contents of .ev file into an array
!  we also allocate the correct amount of memory for dat
!  returns both data and column labels
!+
!----------------------------------------------------------------
subroutine read_evfile(filename,dat,labels,ncols,nsteps,ierr)
 character(len=*),  intent(in)  :: filename
 real, allocatable, intent(out) :: dat(:,:)
 character(len=*),  intent(out) :: labels(:)
 integer,           intent(out) :: ncols,nsteps,ierr
 integer :: i,lu,nheaderlines
 character(len=2048) :: line
 !
 ! open the file
 !
 lu = ianalysis
 open(unit=lu,file=trim(filename),status='old',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    !print "(a,i2)",' ERROR opening '//trim(filename)
    return
 endif
 !
 ! read the first line and extract the column labels
 !
 labels(:) = ' '
 nheaderlines = 1
 read(lu,"(a)",iostat=ierr) line
 call get_column_labels(line,labels,ncols)
 if (ncols <= 0 .or. ierr /= 0) then
    !print "(a,i2)",' ERROR: could not determine number of columns in '//trim(filename)
    return
 endif
 !
 ! count the number of useable lines (timesteps) in the file
 !
 nsteps = 0
 do while(ierr == 0)
    nsteps = nsteps + 1
    read(lu,*,iostat=ierr)
 enddo
 nsteps = nsteps - 1
 !
 ! allocate memory
 !
 allocate(dat(ncols,nsteps))
 dat(:,:) = 0.
 !
 ! rewind the file
 !
 rewind(lu)
 !
 ! skip header lines again
 !
 do i=1,nheaderlines
    read(lu,*,iostat=ierr)
 enddo
 !
 ! read the data
 !
 do i = 1,nsteps
    read(lu,*,iostat=ierr) dat(1:ncols,i)
 enddo
 close(lu)

end subroutine read_evfile

!----------------------------------------------------------------
!+
!  Write the .columns file
!+
!----------------------------------------------------------------
subroutine write_columns_to_file(numcol0,columns0,outputprefix)
 integer,           intent(in) :: numcol0
 character(len=*),  intent(in) :: columns0(numcol0),outputprefix
 integer                       :: i
 character(len=200)            :: columnsfile,label
 !
 ! create and open file
 !
 write(columnsfile,'(2a)') trim(outputprefix),'.columns'
 open(ianalysis,file=trim(columnsfile))
 do i = 1,numcol0
    label = verboseentry(columns0(i))
    write(ianalysis,'(a)') trim(verboseentry(columns0(i)))
 enddo
 close(ianalysis)

end subroutine write_columns_to_file
!----------------------------------------------------------------
!+
!  Read the .evin file
!+
!----------------------------------------------------------------
! Determine non-filename parameters
subroutine read_evin_file(filename,outprefix,single_output,concise,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 character(len=*), intent(out) :: outprefix
 logical,          intent(out) :: single_output,concise
 integer,          intent(out) :: ierr
 integer,          parameter   :: iunit = 21
 type(inopts),     allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(outprefix,     'Output prefix', db,ierr)
 call read_inopt(single_output, 'Num Outputs',   db,ierr)
 call read_inopt(concise,       'Concise',       db,ierr)
 call close_db(db)

end subroutine read_evin_file

!----------------------------------------------------------------
!+
! Determine filenames
!+
!----------------------------------------------------------------
subroutine read_evin_filenames(maxfiles,filename,infilenames,nummodels,ierr)
 integer,            intent(in)    :: maxfiles
 integer,            intent(out)   :: nummodels
 integer,            intent(inout) :: ierr
 character(len=  *), intent(in)    :: filename
 character(len=  *), intent(out)   :: infilenames(maxfiles)
 integer,            parameter     :: iunit = 21
 integer                           :: io,imodel,iequal,iem,ispace
 character(len=512)                :: evinline,infiletmp

 nummodels = 0
 io        = 0
 open(unit=iunit,file=filename)
 do while (io==0)
    read(iunit,'(a)',iostat=io) evinline
    imodel = index(evinline,'Model')
    if (imodel > 0) then
       iequal = index(evinline,'=')
       iem    = index(evinline,'!')
       nummodels = nummodels + 1
       infiletmp = evinline(iequal+1:iem-1)
       ispace = 1
       do while (ispace > 0)
          ispace = index(trim(infiletmp),' ')
          if (ispace==1) then
             infiletmp = infiletmp(2:)
          else if (ispace > 0) then
             infiletmp = trim(infiletmp(:ispace))
          endif
       enddo
       infilenames(nummodels) = trim(infiletmp)
    endif
 enddo
 close(iunit)
 if (nummodels==0) ierr = ierr + 1

end subroutine read_evin_filenames
!----------------------------------------------------------------
!+
!  Write a .evin file
!+
!----------------------------------------------------------------
subroutine write_evin_file(outprefix,infilenames,nummodels,single_output,concise)
 use infile_utils, only:write_inopt
 integer,          intent(in) :: nummodels
 character(len=*), intent(in) :: outprefix,infilenames(:)
 logical,          intent(in) :: single_output,concise
 integer,          parameter  :: iunit = 668
 integer                      :: i
 character(len=120)           :: filename

 filename = trim(outprefix)//".evin"
 print "(a)",' writing input options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(/,a)") '# Output prefix'
 call write_inopt(trim(outprefix),        'Output prefix',   'output file name prefix of all compared files',iunit)

 write(iunit,"(/,a)") '# Files to compare'
 do i = 1,nummodels
    call write_inopt(trim(infilenames(i)),"Model",    'prefix of an .ev file to compare',iunit)
 enddo

 write(iunit,"(/,a)") '# Output Format'
 call write_inopt(single_output,'Num Outputs','output one .ev file per model (F: one output per .ev per model)',iunit)
 call write_inopt(concise,'Concise','write columns that exist in all files (F: write columns that exist in any file)',iunit)

 close(iunit)

end subroutine write_evin_file

!----------------------------------------------------------------
!+
!  convert concise label into verbose label
!+
!----------------------------------------------------------------
character(len=128) function verboseentry(label)
 character(len=*), intent(in) :: label
 character(len=128)           :: l
 !
 ! List of the verbose entries
 !
 l = trim(adjustl(label))
 select case(trim(adjustl(label)))
 case('time');        l='Time'
 case('ekin');        l='Kinetic energy'
 case('etherm');      l='Thermal energy'
 case('emag');        l='Magnetic energy'
 case('epot');        l='Potential energy'
 case('etot');        l='Total energy'
 case('totmom');      l='Total momentum'
 case('angtot');      l='Total angular momentum'
 case('rhomax');      l='max( rho )'
 case('rhomean');     l='mean( rho )'
 case('dt');          l='dt'
 case('totentrop');   l='Total entropy'
 case('rmsmach');     l='RMS Mach number'
 case('rho dust X');  l='max( rho_dust )'
 case('rho dust A');  l='mean( rho_dust )'
 case('rho bdy X');   l='max( rho_boundary )'
 case('rho bdy A');   l='mean( rho_boundary )'
 case('rho star X');  l='max( rho_star )'
 case('rho star A');  l='mean( rho_star )'
 case('rho dm X');    l='max( rho_darkmatter )'
 case('rho dm A');    l='mean( rho_darkmatter )'
 case('rho blg X');   l='max( rho_bulge )'
 case('rho blg A');   l='mean( rho_bulge )'
 case('rho gas X');   l='max( rho_gas )'
 case('rho gas A');   l='mean( rho_gas )'
 case('alpha max');   l='max( alpha )'
 case('alphaB max');  l='max( alphaB )'
 case('divB max');    l='max( divB )'
 case('divB ave');    l='mean( divB )'
 case('hdivB/B max'); l='max( hdivB/B )'
 case('hdivB/B ave'); l='mean( hdivB/B )'
 case('beta max');    l='max( Plasma beta )'
 case('beta ave');    l='mean( Plasma beta )'
 case('beta min');    l='min( Plasma beta )'
 case('temp max');    l='max( Temperature )'
 case('temp ave');    l='mean( Temperature )'
 case('eta_ar max');  l='max( eta_artificial )'
 case('eta_ar ave');  l='mean( eta_artificial )'
 case('eta_o max');   l='max( eta_Ohm )'
 case('eta_o ave');   l='mean( eta_Ohm )'
 case('eta_o min');   l='min( eta_Ohm )'
 case('eta_o/art X'); l='max( eta_Ohm/eta_artificial )'
 case('eta_o/art A'); l='mean( eta_Ohm/eta_artificial )'
 case('eta_o/art N'); l='min( eta_Ohm/eta_artificial )'
 case('eta_h max');   l='max( eta_hall )'
 case('eta_h ave');   l='mean( eta_hall )'
 case('eta_h min');   l='min( eta_hall )'
 case('|eta_h| max'); l='max( |eta_hall| )'
 case('|eta_h| ave'); l='mean( |eta_hall| )'
 case('|eta_h| min'); l='min( |eta_hall| )'
 case('eta_h/art X'); l='max( eta_hall/eta_artificial )'
 case('eta_h/art A'); l='mean( eta_hall/eta_artificial )'
 case('eta_h/art N'); l='min( eta_hall/eta_artificial )'
 case('|e_h|/art X'); l='max( |eta_hall|/eta_artificial )'
 case('|e_h|/art A'); l='mean( |eta_hall|/eta_artificial )'
 case('|e_h|/art N'); l='min( |eta_hall|/eta_artificial )'
 case('eta_a max');   l='max( eta_ambipolar )'
 case('eta_a ave');   l='mean( eta_ambipolar )'
 case('eta_a min');   l='min( eta_ambipolar )'
 case('eta_a/art X'); l='max( eta_ambipolar/eta_artificial )'
 case('eta_a/art A'); l='mean( eta_ambipolar/eta_artificial )'
 case('eta_a/art N'); l='min( eta_ambipolar/eta_artificial )'
 case('n_e/n max');   l='max( n_e/n )'
 case('n_e/n ave');   l='mean( n_e/n )'
 case('n_e max');     l='max( n_electron )'
 case('n_e ave');     l='mean( n_electron )'
 case('n_n max');     l='max( n_neutral )'
 case('n_n ave');     l='mean( n_neutral )'
 case('Z_Grain max'); l='max( Z_grain )'
 case('Z_Grain ave'); l='mean( Z_grain )'
 case('Z_Grain min'); l='min( Z_grain )'
 case('n_ihR max');   l='max( n_ion_lightElements (cosmic rays) )'
 case('n_ihR ave');   l='mean( n_ion_lightElements (cosmic rays) )'
 case('n_imR max');   l='max( n_ion_metallicElements (cosmic rays) )'
 case('n_imR ave');   l='mean( n_ion_metallicElements (cosmic rays) )'
 case('n_gR max');    l='max( n_grain (cosmic rays) )'
 case('n_gR ave');    l='mean( n_grain (cosmic rays) )'
 case('n_iT max');    l='max( n_ion (thermal) )'
 case('n_iT ave');    l='mean( n_ion (thermal) )'
 case('n_H+ max');    l='max( n_Hydrogen+ )'
 case('n_H+ ave');    l='mean( n_Hydrogen+ )'
 case('n_He+ max');   l='max( n_Helium+ )'
 case('n_He+ ave');   l='mean( n_Helium+ )'
 case('n_Na+ max');   l='max( n_Sodium+ )'
 case('n_Na+ ave');   l='mean( n_Sodium+ )'
 case('n_Mg+ max');   l='max( n_Magnesium+ )'
 case('n_Mg+ ave');   l='mean( n_Magnesium+ )'
 case('n_K+ max');    l='max( n_Potassium+ )'
 case('n_K+ ave');    l='mean( n_Potassium+ )'
 case('n_He++ max');  l='max( n_Helium++ )'
 case('n_He++ ave');  l='mean( n_Helium++ )'
 case('n_Na++ max');  l='max( n_Sodium++ )'
 case('n_Na++ ave');  l='mean( n_Sodium++ )'
 case('n_Mg++ max');  l='max( n_Magnesium++ )'
 case('n_Mg++ ave');  l='mean( n_Magnesium++ )'
 case('n_K++ max');   l='max( n_Potassium++ )'
 case('n_K++ ave');   l='mean( n_Potassium++ )'
 case('dust/gas X');  l='max( dust/gas )'
 case('dust/gas A');  l='mean( dust/gas )'
 case('dust/gas N');  l='min( dust/gas )'
 case('t_s ave');     l='mean( t_s )'
 case('t_s min');     l='min( t_s )'
 case('totmomall');   l='Total momentum (incl. accreted)'
 case('angall');      l='Total anglular momentum (incl. accreted)'
 case('Macc sink 1'); l='Accreted mass onto star 1'
 case('Macc sink 2'); l='Accreted mass onto star 2'
 case('accretedmas'); l='Accreted mass'
 case('eacc');        l='Accreted energy'
 case('tot lum');     l='Total luminosity'
 case('erot');        l='Rotational energy'
 case('visc_rat X');  l='max( Viscious ratio )'
 case('visc_rat A');  l='mean( Viscious ratio )'
 case('visc_rat N');  l='min( Viscious ratio )'
 end select
 verboseentry = l

end function verboseentry

end module evutils
