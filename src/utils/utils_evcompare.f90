!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
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
!  DEPENDENCIES: infile_utils, io
!+
!--------------------------------------------------------------------------
module evutils
 implicit none
 !
 ! Subroutines
 !
 public :: get_column_labels_from_ev,write_columns_to_file
 public :: read_evin_file,read_evin_filenames,write_evin_file
 public :: find_column
 !
 ! maximum number of .ev files
 !
 integer, parameter, public  :: inumev = 150
 !
 ! logical unit number for read/write operations
 !
 integer, parameter :: ianalysis = 25

 private

contains
!----------------------------------------------------------------
!+
!  Read the column labels from the .ev file
!  Note: This formatting is specific to the formating in evwrite.f90
!+
!----------------------------------------------------------------
subroutine get_column_labels_from_ev(evfile,labels,numcol,ierr)
 integer,             intent(out) :: numcol,ierr
 character(len=*   ), intent(in)  :: evfile
 character(len=*   ), intent(out) :: labels(:)
 integer                          :: i,iopen,iclose
 logical                          :: iexist
 character(len=2048)              :: line

 ierr = 0
 inquire(file=trim(evfile),exist=iexist)
 if (iexist) then
    open(unit=ianalysis,file=trim(evfile))
    read(ianalysis,'(a)') line
    close(ianalysis)
    i     = 0
    iopen = 1 ! to get into the loop
    do while ( iopen > 0 .and. i < inumev)
       iopen  = index(line,'[')
       iclose = index(line,']')
       i = i + 1
       labels(i) = line(iopen +3:iclose-1)
       line      = line(iclose+1:len(line))
    enddo
    numcol = i - 1
 else
    ierr   = 1
    numcol = 0
 endif

end subroutine get_column_labels_from_ev


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
!  A list of the possible columns, written verbosely
!+
!----------------------------------------------------------------
character(len=128) function verboseentry(label)
 character(len=*), intent(in) :: label
 character(len=12)            :: columnT(inumev)
 character(len=128)           :: columnV(inumev)
 integer                      :: i,j
 logical                      :: keep_searching
 !
 ! List of the verbose entries
 !
 i = 1
 columnT(i) = "        time" ; columnV(i) = "Time"; i=i+1
 columnT(i) = "        ekin" ; columnV(i) = "Kinetic energy"; i=i+1
 columnT(i) = "      etherm" ; columnV(i) = "Thermal energy"; i=i+1
 columnT(i) = "        emag" ; columnV(i) = "Magnetic energy"; i=i+1
 columnT(i) = "        epot" ; columnV(i) = "Potential energy"; i=i+1
 columnT(i) = "        etot" ; columnV(i) = "Total energy"; i=i+1
 columnT(i) = "      totmom" ; columnV(i) = "Total momentum"; i=i+1
 columnT(i) = "      angtot" ; columnV(i) = "Total angular momentum"; i=i+1
 columnT(i) = "      rhomax" ; columnV(i) = "max( rho )"; i=i+1
 columnT(i) = "     rhomean" ; columnV(i) = "mean( rho )"; i=i+1
 columnT(i) = "          dt" ; columnV(i) = "dt"; i=i+1
 columnT(i) = "   totentrop" ; columnV(i) = "Total entropy"; i=i+1
 columnT(i) = "     rmsmach" ; columnV(i) = "RMS Mach number"; i=i+1
 columnT(i) = "  rho dust X" ; columnV(i) = "max( rho_dust )"; i=i+1
 columnT(i) = "  rho dust A" ; columnV(i) = "mean( rho_dust )"; i=i+1
 columnT(i) = "   rho bdy X" ; columnV(i) = "max( rho_boundary )"; i=i+1
 columnT(i) = "   rho bdy A" ; columnV(i) = "mean( rho_boundary )"; i=i+1
 columnT(i) = "  rho star X" ; columnV(i) = "max( rho_star )"; i=i+1
 columnT(i) = "  rho star A" ; columnV(i) = "mean( rho_star )"; i=i+1
 columnT(i) = "    rho dm X" ; columnV(i) = "max( rho_darkmatter )"; i=i+1
 columnT(i) = "    rho dm A" ; columnV(i) = "mean( rho_darkmatter )"; i=i+1
 columnT(i) = "   rho blg X" ; columnV(i) = "max( rho_bulge )"; i=i+1
 columnT(i) = "   rho blg A" ; columnV(i) = "mean( rho_bulge )"; i=i+1
 columnT(i) = "   rho gas X" ; columnV(i) = "max( rho_gas )"; i=i+1
 columnT(i) = "   rho gas A" ; columnV(i) = "mean( rho_gas )"; i=i+1
 columnT(i) = "   alpha max" ; columnV(i) = "max( alpha )"; i=i+1
 columnT(i) = "  alphaB max" ; columnV(i) = "max( alphaB )"; i=i+1
 columnT(i) = "    divB max" ; columnV(i) = "max( divB )"; i=i+1
 columnT(i) = "    divB ave" ; columnV(i) = "mean( divB )"; i=i+1
 columnT(i) = " hdivB/B max" ; columnV(i) = "max( hdivB/B )"; i=i+1
 columnT(i) = " hdivB/B ave" ; columnV(i) = "mean( hdivB/B )"; i=i+1
 columnT(i) = "    beta max" ; columnV(i) = "max( Plasma beta )"; i=i+1
 columnT(i) = "    beta ave" ; columnV(i) = "mean( Plasma beta )"; i=i+1
 columnT(i) = "    beta min" ; columnV(i) = "min( Plasma beta )"; i=i+1
 columnT(i) = "    temp max" ; columnV(i) = "max( Temperature )"; i=i+1
 columnT(i) = "    temp ave" ; columnV(i) = "mean( Temperature )"; i=i+1
 columnT(i) = "  eta_ar max" ; columnV(i) = "max( eta_artificial )"; i=i+1
 columnT(i) = "  eta_ar ave" ; columnV(i) = "mean( eta_artificial )"; i=i+1
 columnT(i) = "   eta_o max" ; columnV(i) = "max( eta_Ohm )"; i=i+1
 columnT(i) = "   eta_o ave" ; columnV(i) = "mean( eta_Ohm )"; i=i+1
 columnT(i) = "   eta_o min" ; columnV(i) = "min( eta_Ohm )"; i=i+1
 columnT(i) = " eta_o/art X" ; columnV(i) = "max( eta_Ohm/eta_artificial )"; i=i+1
 columnT(i) = " eta_o/art A" ; columnV(i) = "mean( eta_Ohm/eta_artificial )"; i=i+1
 columnT(i) = " eta_o/art N" ; columnV(i) = "min( eta_Ohm/eta_artificial )"; i=i+1
 columnT(i) = "   eta_h max" ; columnV(i) = "max( eta_hall )"; i=i+1
 columnT(i) = "   eta_h ave" ; columnV(i) = "mean( eta_hall )"; i=i+1
 columnT(i) = "   eta_h min" ; columnV(i) = "min( eta_hall )"; i=i+1
 columnT(i) = " |eta_h| max" ; columnV(i) = "max( |eta_hall| )"; i=i+1
 columnT(i) = " |eta_h| ave" ; columnV(i) = "mean( |eta_hall| )"; i=i+1
 columnT(i) = " |eta_h| min" ; columnV(i) = "min( |eta_hall| )"; i=i+1
 columnT(i) = " eta_h/art X" ; columnV(i) = "max( eta_hall/eta_artificial )"; i=i+1
 columnT(i) = " eta_h/art A" ; columnV(i) = "mean( eta_hall/eta_artificial )"; i=i+1
 columnT(i) = " eta_h/art N" ; columnV(i) = "min( eta_hall/eta_artificial )"; i=i+1
 columnT(i) = " |e_h|/art X" ; columnV(i) = "max( |eta_hall|/eta_artificial )"; i=i+1
 columnT(i) = " |e_h|/art A" ; columnV(i) = "mean( |eta_hall|/eta_artificial )"; i=i+1
 columnT(i) = " |e_h|/art N" ; columnV(i) = "min( |eta_hall|/eta_artificial )"; i=i+1
 columnT(i) = "   eta_a max" ; columnV(i) = "max( eta_ambipolar )"; i=i+1
 columnT(i) = "   eta_a ave" ; columnV(i) = "mean( eta_ambipolar )"; i=i+1
 columnT(i) = "   eta_a min" ; columnV(i) = "min( eta_ambipolar )"; i=i+1
 columnT(i) = " eta_a/art X" ; columnV(i) = "max( eta_ambipolar/eta_artificial )"; i=i+1
 columnT(i) = " eta_a/art A" ; columnV(i) = "mean( eta_ambipolar/eta_artificial )"; i=i+1
 columnT(i) = " eta_a/art N" ; columnV(i) = "min( eta_ambipolar/eta_artificial )"; i=i+1
 columnT(i) = "   n_e/n max" ; columnV(i) = "max( n_e/n )"; i=i+1
 columnT(i) = "   n_e/n ave" ; columnV(i) = "mean( n_e/n )"; i=i+1
 columnT(i) = "     n_e max" ; columnV(i) = "max( n_electron )"; i=i+1
 columnT(i) = "     n_e ave" ; columnV(i) = "mean( n_electron )"; i=i+1
 columnT(i) = "     n_n max" ; columnV(i) = "max( n_neutral )"; i=i+1
 columnT(i) = "     n_n ave" ; columnV(i) = "mean( n_neutral )"; i=i+1
 columnT(i) = " Z_Grain max" ; columnV(i) = "max( Z_grain )"; i=i+1
 columnT(i) = " Z_Grain ave" ; columnV(i) = "mean( Z_grain )"; i=i+1
 columnT(i) = " Z_Grain min" ; columnV(i) = "min( Z_grain )"; i=i+1
 columnT(i) = "   n_ihR max" ; columnV(i) = "max( n_ion_lightElements (cosmic rays) )"; i=i+1
 columnT(i) = "   n_ihR ave" ; columnV(i) = "mean( n_ion_lightElements (cosmic rays) )"; i=i+1
 columnT(i) = "   n_imR max" ; columnV(i) = "max( n_ion_metallicElements (cosmic rays) )"; i=i+1
 columnT(i) = "   n_imR ave" ; columnV(i) = "mean( n_ion_metallicElements (cosmic rays) )"; i=i+1
 columnT(i) = "    n_gR max" ; columnV(i) = "max( n_grain (cosmic rays) )"; i=i+1
 columnT(i) = "    n_gR ave" ; columnV(i) = "mean( n_grain (cosmic rays) )"; i=i+1
 columnT(i) = "    n_iT max" ; columnV(i) = "max( n_ion (thermal) )"; i=i+1
 columnT(i) = "    n_iT ave" ; columnV(i) = "mean( n_ion (thermal) )"; i=i+1
 columnT(i) = "    n_H+ max" ; columnV(i) = "max( n_Hydrogen+ )"; i=i+1
 columnT(i) = "    n_H+ ave" ; columnV(i) = "mean( n_Hydrogen+ )"; i=i+1
 columnT(i) = "   n_He+ max" ; columnV(i) = "max( n_Helium+ )"; i=i+1
 columnT(i) = "   n_He+ ave" ; columnV(i) = "mean( n_Helium+ )"; i=i+1
 columnT(i) = "   n_Na+ max" ; columnV(i) = "max( n_Sodium+ )"; i=i+1
 columnT(i) = "   n_Na+ ave" ; columnV(i) = "mean( n_Sodium+ )"; i=i+1
 columnT(i) = "   n_Mg+ max" ; columnV(i) = "max( n_Magnesium+ )"; i=i+1
 columnT(i) = "   n_Mg+ ave" ; columnV(i) = "mean( n_Magnesium+ )"; i=i+1
 columnT(i) = "    n_K+ max" ; columnV(i) = "max( n_Potassium+ )"; i=i+1
 columnT(i) = "    n_K+ ave" ; columnV(i) = "mean( n_Potassium+ )"; i=i+1
 columnT(i) = "  n_He++ max" ; columnV(i) = "max( n_Helium++ )"; i=i+1
 columnT(i) = "  n_He++ ave" ; columnV(i) = "mean( n_Helium++ )"; i=i+1
 columnT(i) = "  n_Na++ max" ; columnV(i) = "max( n_Sodium++ )"; i=i+1
 columnT(i) = "  n_Na++ ave" ; columnV(i) = "mean( n_Sodium++ )"; i=i+1
 columnT(i) = "  n_Mg++ max" ; columnV(i) = "max( n_Magnesium++ )"; i=i+1
 columnT(i) = "  n_Mg++ ave" ; columnV(i) = "mean( n_Magnesium++ )"; i=i+1
 columnT(i) = "   n_K++ max" ; columnV(i) = "max( n_Potassium++ )"; i=i+1
 columnT(i) = "   n_K++ ave" ; columnV(i) = "mean( n_Potassium++ )"; i=i+1
 columnT(i) = "  dust/gas X" ; columnV(i) = "max( dust/gas )"; i=i+1
 columnT(i) = "  dust/gas A" ; columnV(i) = "mean( dust/gas )"; i=i+1
 columnT(i) = "  dust/gas N" ; columnV(i) = "min( dust/gas )"; i=i+1
 columnT(i) = "     t_s ave" ; columnV(i) = "mean( t_s )"; i=i+1
 columnT(i) = "     t_s min" ; columnV(i) = "min( t_s )"; i=i+1
 columnT(i) = "   totmomall" ; columnV(i) = "Total momentum (incl. accreted)"; i=i+1
 columnT(i) = "      angall" ; columnV(i) = "Total anglular momentum (incl. accreted)"; i=i+1
 columnT(i) = " Macc sink 1" ; columnV(i) = "Accreted mass onto star 1"; i=i+1
 columnT(i) = " Macc sink 2" ; columnV(i) = "Accreted mass onto star 2"; i=i+1
 columnT(i) = " accretedmas" ; columnV(i) = "Accreted mass"; i=i+1
 columnT(i) = "        eacc" ; columnV(i) = "Accreted energy"; i=i+1
 columnT(i) = "     tot lum" ; columnV(i) = "Total luminosity"; i=i+1
 columnT(i) = "        erot" ; columnV(i) = "Rotational energy"; i=i+1
 columnT(i) = "  visc_rat X" ; columnV(i) = "max( Viscious ratio )"; i=i+1
 columnT(i) = "  visc_rat A" ; columnV(i) = "mean( Viscious ratio )"; i=i+1
 columnT(i) = "  visc_rat N" ; columnV(i) = "min( Viscious ratio )"; i=i+1
 !
 !--Determine the entry to use
 !
 verboseentry = label ! default value
 j = 0
 keep_searching = .true.
 do while (keep_searching .and. j < i - 1)
    j = j + 1
    if (label==columnT(j)) then
       verboseentry   = columnV(j)
       keep_searching = .false.
    endif
 enddo

end function verboseentry
!-----------------------------------------------------------------------
end module evutils
