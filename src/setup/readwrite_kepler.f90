!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_kepler
!
! Utilities for reading in/out stellar profiles from the Kepler stellar
! evolution code
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: datafiles, fileutils, table_utils, units
!
 implicit none

 public :: read_kepler_file,write_kepler_comp

 private

contains

!-----------------------------------------------------------------------
!+
!  Read in datafile from the KEPLER stellar evolution code
!+
!-----------------------------------------------------------------------
subroutine read_kepler_file(filepath,ng_max,n_rows,rtab,rhotab,ptab,mtab,temperature,&
                            enitab,totmass,composition,comp_label,Xfrac,Yfrac,columns_compo,&
                            ierr,mcut,rcut,cgsunits)
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 use datafiles, only:find_phantom_datafile
 use fileutils, only:get_ncolumns,get_nlines,skip_header,get_column_labels

 integer,intent(in)                       :: ng_max
 integer,intent(out)                      :: ierr,n_rows
 real,allocatable,intent(out)             :: rtab(:),rhotab(:),ptab(:),temperature(:),enitab(:),mtab(:),Xfrac(:),Yfrac(:)
 real,intent(out),allocatable             :: composition(:,:)
 real,intent(out)                         :: totmass
 real,intent(out),optional                :: rcut
 real,intent(in),optional                 :: mcut
 character(len=20),allocatable,intent(out):: comp_label(:)
 character(len=20),allocatable            :: all_label(:) !This contains all labels read from KEPLER file.
 character(len=*),intent(in)              :: filepath
 integer,intent(out)                      :: columns_compo
 logical, intent(in), optional            :: cgsunits

 character(len=120)                       :: fullfilepath
 character(len=10000)                     :: line

 integer                                  :: k,aloc,i,column_no,n_cols,n_labels,iu
 integer                                  :: nheaderlines,skip_no
 real,allocatable                         :: stardata(:,:)
 logical                                  :: iexist,n_too_big,composition_available,usecgs

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
 !
 !--Read data from file
 !
 k = 1
 n_rows = 0
 n_cols = 0
 n_too_big = .false.
 composition_available = .false.
 skip_no = 0
 !
 !--Calculate the number of rows, columns and comments in kepler file.
 !
 n_rows = get_nlines(trim(fullfilepath),skip_comments=.true.,n_columns=n_cols,n_headerlines=nheaderlines)
 !
 !--Check if the number of rows is 0 or greater than ng_max.
 !
 if (n_rows < 1 .or. n_cols < 1) then
    ierr = 2
    return
 endif

 if (n_rows >= ng_max) n_too_big = .true.

 if (n_too_big) then
    ierr = 3
    return
 endif

 print "(a,i5)",' number of data rows = ', n_rows
 !
 !--If there is 'nt1' in the column headings, we know the file contains composition.
 !--This is used as a test for saving composition.
 !
 ierr = 0
 open(newunit=iu,file=trim(fullfilepath))
 !The row with the information about column headings is at nheaderlines-1.
 !we skip the first nheaderlines-2 rows and then read the nheaderlines-1 to find the substrings
 call skip_header(iu,nheaderlines-2,ierr)
 read(iu, '(a)',iostat=ierr) line

 !read the column labels and store them in an array.
 allocate(all_label(n_cols))
 call get_column_labels(line,n_labels,all_label)
 close(iu)

 !check which label gives nt1.
 do i = 1,len(all_label)
    if (all_label(i) == 'nt1') then
       skip_no = i - 1
       composition_available = .true.
    endif
 enddo

 !Allocate memory for saving data
 allocate(stardata(n_rows, n_cols))
 allocate(Xfrac(n_rows),Yfrac(n_rows))
 allocate(rtab(n_rows),rhotab(n_rows),ptab(n_rows),temperature(n_rows),enitab(n_rows),mtab(n_rows))
 !
 !--Read the file again and save the data in stardata tensor.
 !
 open(newunit=iu,file=trim(fullfilepath))
 call skip_header(iu,nheaderlines,ierr)
 do k=1,n_rows
    read(iu,*,iostat=ierr) stardata(k,:)
 enddo
 close(iu)
 !
 !--Save the relevant information we require into arrays that can be used later.
 !--convert relevant data from CGS to code units
 !
 !radius
 stardata(1:n_rows,4)  = stardata(1:n_rows,4)
 rtab(1:n_rows)        = stardata(1:n_rows,4)

 !density
 stardata(1:n_rows,6)  = stardata(1:n_rows,6)
 rhotab(1:n_rows)      = stardata(1:n_rows,6)
 !mass
 stardata(1:n_rows,3)  = stardata(1:n_rows,3)
 mtab(1:n_rows)        = stardata(1:n_rows,3)
 totmass               = stardata(n_rows,3)
 !pressure
 stardata(1:n_rows,8)  = stardata(1:n_rows,8)
 ptab(1:n_rows)        = stardata(1:n_rows,8)

 !temperature
 temperature(1:n_rows) = stardata(1:n_rows,7)

 !specific internal energy
 stardata(1:n_rows,9)  = stardata(1:n_rows,9)
 enitab(1:n_rows)      = stardata(1:n_rows,9)

 if (.not. usecgs) then
    mtab = mtab / umass
    rtab = rtab / udist
    ptab = ptab / unit_pressure
    rhotab = rhotab / unit_density
    enitab = enitab / unit_ergg
    totmass = totmass / umass
 endif

 print*, 'Total Mass: ', totmass, 'Max radius: ', rtab(n_rows)

 !if elements were found in the file read, save the composition by allocating an array
 !else set it to 0
 if (composition_available) then
    !saving composition. In a composition file of KEPLER, we have first 13 columns that do not contain composition
    !We skip them and store the rest into a composition tensor.
    print*, 'Kepler file has composition.'
    columns_compo = 0
    allocate(composition(n_rows,n_cols-skip_no))
    allocate(comp_label(n_cols-skip_no))
    comp_label(:) = all_label(skip_no+1:n_cols)
    column_no = skip_no + 1
    do i = 1, n_cols-skip_no
       composition(:,i) = stardata(:,column_no)
       column_no        = column_no + 1
    enddo
    columns_compo = n_cols-skip_no
 else
    allocate(composition(0,0))
    allocate(comp_label(0))
 endif

 Xfrac = composition(:,2) ! save H1 fraction
 Yfrac = composition(:,4) ! save He4 fraction

 print*, shape(composition),'shape of composition array'

 if (present(rcut) .and. present(mcut)) then
    aloc = minloc(abs(stardata(1:n_rows,1) - mcut),1)
    rcut = rtab(aloc)
    print*, 'rcut = ', rcut
 endif
 print*, 'Finished reading KEPLER file'
 print*, '------------------------------------------------------------'

end subroutine read_kepler_file

!-----------------------------------------------------------------------
!+
!  Write kepler.comp which is file that contains interpolated composition
!  for each particle. This works for ikepler only.
!  Interpolating composition for each particle in the star.
!  For now I write a file with the interpolate data
!  This is used in analysis kepler file to bin composition.
!
!+
!-----------------------------------------------------------------------
subroutine write_kepler_comp(filename,composition,comp_label,columns_compo,r,&
                             xyzh,npart,npts,composition_exists,npin)

 use table_utils, only:yinterp
 character(len=*), intent(in)               :: filename
 integer, intent(in)                        :: columns_compo,npart,npts
 real,    intent(in)                        :: xyzh(:,:)
 real, allocatable,intent(in)               :: r(:)
 real, allocatable, intent(in)              :: composition(:,:)
 character(len=20), allocatable,intent(in)  :: comp_label(:)
 logical, intent(out)                       :: composition_exists
 integer, intent(in), optional              :: npin
 real , allocatable                         :: compositioni(:,:)
 real, allocatable                          :: comp(:)
 integer                                    :: i,j,iu,i1,ierr
 real                                       :: ri

 composition_exists = .false.
 i1 = 0
 if (present(npin)) i1=npin

 ! !Check if composition exists. If composition array is non-zero, we use it to interpolate composition for each particle
 ! !in the star.
 if (columns_compo /= 0) then
    composition_exists = .true.
 endif

 if (composition_exists) then
    print*, 'Writing the stellar composition for each particle into '//trim(filename)

    open(newunit=iu,file=trim(filename),status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print*,' ERROR writing to '//trim(filename)
       return
    endif
    write(iu,"('#',50(1x,'[',1x,a7,']',2x))") &
            comp_label
    !Now setting the composition of star if the case used was ikepler
    allocate(compositioni(columns_compo,1))
    allocate(comp(1:npts))
    do i = i1+1,npart
       !Interpolate compositions
       ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))

       do j = 1,columns_compo
          comp(1:npts)      = composition(1:npts,j)
          compositioni(j,1) = yinterp(comp(1:npts),r(1:npts),ri)
       enddo
       write(iu,'(50(es18.10,1X))') &
            (compositioni(j,1),j=1,columns_compo)
    enddo
    close(iu)
 endif

end subroutine write_kepler_comp

end module readwrite_kepler
