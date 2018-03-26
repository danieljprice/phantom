!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantomevcompare
!
!  DESCRIPTION: For the given input .ev files, will rewrite them using
!               common headers
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantomevcompare [no arguments]
!
!  DEPENDENCIES: evutils, prompting
!+
!--------------------------------------------------------------------------
program phantomevcompare
 use prompting, only: prompt
 use evutils,   only: max_columns,get_column_labels_from_ev,read_evin_file,read_evin_filenames, &
                      write_evin_file,write_columns_to_file
 implicit none
 integer, parameter  :: maxfiles  = 100
 logical, parameter  :: write_columns = .true.
 integer             :: i,j,k,ki,kj,iio,idot,numcol,numcol0,nummodels,nargs,ierr
 integer             :: headerorder(max_columns)
 real                :: oldorder(max_columns),neworder(max_columns)
 logical             :: interactive,single_output,concise
 logical             :: enter_filename,add_prefix,add_prefix_verified,get_newfile
 logical             :: initial_entry,keep_searching,iexist,ifailed,readfailed,make_outname
 character(len=  12) :: columns0(max_columns),columns(max_columns)
 character(len=  19) :: columnsEV(max_columns)
 character(len=  64) :: filename,evfile,evinfile,outprefix,modelprefix,ev_fmtD,ev_fmtH
 character(len=  64) :: infilenames(maxfiles),outfilenames(maxfiles),outfilename0a,outfilename0b
 character(len= 256) :: evout_new
 character(len=4096) :: cdummy
 integer, parameter  :: icolA = 26, icolB = 27

 !
 !--If argument exists, read it in
 !
 interactive = .true.
 nargs       = command_argument_count()
 evinfile    = ''
 if (nargs > 0) then
    call get_command_argument(1,evinfile)
    idot = index(evinfile,'.evin')  ! strip .evin if present
    if ( idot > 0 ) evinfile = evinfile(1:idot-1)
    evinfile = trim(evinfile)//'.evin'
    inquire(file=trim(evinfile),exist=iexist)
    if (iexist) then
       ! read in non-filename parameters from .evin
       call read_evin_file(evinfile,outprefix,single_output,concise,ierr)
       ! read in filenames
       call read_evin_filenames(maxfiles,evinfile,infilenames,nummodels,ierr)
       if (ierr > 0) then
          write(*,'(3a)') "Did not successfully read ",trim(evinfile),".  Entering interactive mode"
       else
          interactive = .false.
       endif
    else
       write(*,'(3a)') "WARNING: ",trim(evinfile)," does not exists.  Entering interactive mode"
    endif
 endif
 !
 !--Through an interactive session, determine which files will be compared, and by which method
 !
 if (interactive) then
    !--Ask for details about what columns to compare
    !  (i.e. use all available columns (T) or only columns exisitng in all files (F))
    concise = .true.
    call prompt('Do you wish to include only columns existing in all files?',concise)
    !
    !--Make a new .ev file for each input .ev file (F), or one per model (T)
    single_output = .true.
    call prompt('Do you want only one output file per model?',single_output)
    !
    !--List the models to compare
    !  Note: will check to ensure that 01.ev exists
    enter_filename = .true.
    write(*,'(a)') ' '
    write(*,'(a)') 'Enter the file prefixes you wish to compare; press enter when complete.'
    i = 0
    do while (enter_filename)
       if (i==maxfiles) then
          !--Failsafe to complie with output file name format
          write(*,'(a)') 'WARNING! You have exceeded the maximum number of files; beginning comparison now.'
          enter_filename = .false.
       else
          add_prefix = .true.
          modelprefix = ''
          call prompt('Enter file prefix:',modelprefix)
       endif
       if (trim(modelprefix)=='') then
          !--Completed entering model names
          enter_filename = .false.
       else
          !--Verify 01.ev exists
          filename=trim(modelprefix)//'01.ev'
          inquire(file=filename,exist=iexist)
          if (.not. iexist) then
             write(*,'(3a)') 'WARNING! ',trim(filename),' does not exists.  Not adding to lists of file to compare.'
             add_prefix = .false.
          endif
          !--A candidate to add to the list of models to compare
          if (add_prefix) then
             add_prefix_verified = .true.
             if (i > 0) then
                !--Ensure the model is not already in the list
                do j = 0,i-1
                   if (trim(modelprefix)==trim(infilenames(j+1))) then
                      write(*,'(3a)') 'WARNING! ',trim(modelprefix),' is already in the list!'
                      add_prefix_verified = .false.
                   endif
                enddo
             endif
             if (add_prefix_verified) then
                !--Add model to the list of models to compare
                i = i + 1
                infilenames(i) = trim(modelprefix)
             endif
          endif
       endif       ! endif: (trim(modelprefix)=='')
    enddo          ! enddo:  while (enter_filename)
    nummodels = i   ! the number of models to compare
 else
    !
    !--Ensure that all the 01.ev files from the .evin file exists, and ensure that there are no duplicates
    do i = 1,nummodels
       filename=trim(infilenames(i))//'01.ev'
       inquire(file=filename,exist=iexist)
       if (iexist) then
          do j = i+1,nummodels
             if (trim(infilenames(i))==trim(infilenames(j))) infilenames(j) = "doesnotexists"
          enddo
       else
          infilenames(i) = "doesnotexists"
          write(*,'(3a)') 'WARNING! ',trim(filename),' does not exists. Removing from list of files to compare.'
       endif
    enddo
    !--Check list and remove any non-existent names
    i = 0
    do while (i < nummodels)
       if (trim(infilenames(i))=="doesnotexists") then
          do j = i,nummodels-1
             infilenames(j) = infilenames(j+1)
          enddo
          nummodels = nummodels - 1
       else
          i = i + 1
       endif
    enddo
 endif   ! endif: interactive
 if (nummodels==0) then
    write(*,'(a)') 'There are no files to compare'
    stop
 endif
 !
 !--Convert infilenames to outfilenames.  Only required if comparing files in different folders
 outfilenames = infilenames
 do i = 1,nummodels
    make_outname = .true.
    do while (make_outname)
       idot = index(outfilenames(i),'..')
       if (idot > 0) then
          outfilename0a = outfilenames(i)
          write(outfilename0b,'(2a)') trim(outfilename0a(:idot-1)),trim(outfilename0a(idot+2:))
          outfilenames(i) = outfilename0b
       else
          make_outname = .false.
       endif
    enddo
    make_outname = .true.
    do while (make_outname)
       idot = index(outfilenames(i),'/')
       if (idot > 0) then
          outfilename0a = outfilenames(i)
          write(outfilename0b,'(3a)') trim(outfilename0a(:idot-1)),'_',trim(outfilename0a(idot+1:))
          outfilenames(i) = outfilename0b
       else
          make_outname = .false.
       endif
    enddo
 enddo
 !
 !--List models to compare
 !  This version allows one model to be 'compared' since it may combine multiple .ev files, and print a .columns file
 !
 write(*,'(a)') ' '
 write(*,'(a,I3,a)') 'The following ',nummodels,' simulations will be compared:'
 do i = 1,nummodels
    write(*,'(5x,a)') trim(infilenames(i))
 enddo
 write(*,'(a)') ' '
 !
 !--Read in every .ev file to determine which columns are to be included
 !
 initial_entry = .true.
 do j = 1,nummodels
    i = 0
    ierr = 0
    do while (ierr==0)
       i = i + 1
       write(evfile,'(a,I2.2,a)') trim(infilenames(j)),i,'.ev'
       if (initial_entry) then
          !--Get the columns of the first file and use as a baseline
          call get_column_labels_from_ev(evfile,columns0,numcol0,ierr)
          if (ierr /=0) stop 'evcompare: file does not exist (sanity check: should not happen)'
          initial_entry = .false.
       else
          !--Read the columns from the next file and compare; make modifications to columns0
          call get_column_labels_from_ev(evfile,columns,numcol,ierr)
          if (ierr==0) then ! if /=0, file does not exist, which acceptable/expected for i > 1
             !
             if ( concise ) then
                !--Remove entries from columns0 that do not exists in columns
                ki = 0
                do while ( ki < numcol0 )
                   keep_searching = .true.
                   ki = ki + 1
                   do kj = 1,numcol
                      if (columns(kj)==columns0(ki)) keep_searching = .false.
                   enddo
                   if ( keep_searching ) then
                      !--Delete entry ki from columns0
                      do k = ki,numcol0-1
                         columns0(k) = columns0(k+1)
                      enddo
                      numcol0 = numcol0 - 1
                      ki      = ki - 1       ! the new entry needs to be verified
                   endif
                enddo
             else
                !--Add entries to columns0 that exist in columns
                do kj = 1,numcol
                   ki = 0
                   keep_searching = .true.
                   do while ( ki < numcol0 )
                      ki = ki + 1
                      if (columns(kj)==columns0(ki)) keep_searching = .false.
                   enddo
                   if ( keep_searching ) then
                      !--Add entry kj to columns0
                      numcol0           = numcol0 + 1
                      columns0(numcol0) = columns(kj)
                   endif
                enddo
             endif! end:if (concise)
             !
          endif   ! end: if file exists
       endif      ! end:if (initial_entry)
    enddo         ! end : do while (ierr==0)
 enddo            ! end: do j = 1,nummodels
 !
 !--Determine prefix for new output (if not read in)
 !
 if (evinfile=='') then
    i           = 0
    get_newfile = .true.
    do while (get_newfile)
       i = i + 1
       write(filename,'(a,I2.2,a)')'evCompare',i,'.columns'
       inquire(file=filename,exist=iexist)
       if (.not.iexist .or. i==99) then  ! i=99 is a failsafe for the current format
          write(outprefix,'(a,I2.2)')'evCompare',i
          get_newfile = .false.
       endif
    enddo
 else
    idot      = index(evinfile,'.evin')    ! .evin will be in the filename by construction
    outprefix = evinfile(1:idot-1)
 endif
 !
 !--Write .columns file (verbosely)
 if ( write_columns ) then
    call write_columns_to_file(numcol0,columns0,outprefix)
 endif
 write(*,'(3a,I4,a)') 'The .columns file is ',trim(outprefix),'.columns, with ',numcol0,' entries'
 !
 !--Write and store labels in proper format for .ev files
 !
 do i = 1,numcol0
    write(columnsEV(i),"(1x,'[',i2.2,a12,']',2x)")i,columns0(i)
 enddo
 !
 !--Rewrite all the .ev files using the new columns ordering
 !
 write(ev_fmtH,'(a,I3,a)') '(',numcol0+1,'a)'
 write(ev_fmtD,'(a,I3,a)') '(',numcol0,'(1pe18.10,1x))'
 write(*,'(a)') ' '
 write(*,'(a)') 'The following new .ev files are being created: '
 ifailed    = .false.
 readfailed = .false.
 do j = 1,nummodels
    i    = 0
    ierr = 0
    do while (ierr==0)
       i = i + 1
       write(evfile,'(a,I2.2,a)') trim(infilenames(j)),i,'.ev'
       call get_column_labels_from_ev(evfile,columns,numcol,ierr)
       if (ierr==0) then
          !--Compare the header to the master header
          headerorder = 0
          do kj = 1,numcol
             ki = 0
             keep_searching = .true.
             do while (keep_searching .and. ki < numcol0)
                ki = ki + 1
                if (columns(kj)==columns0(ki)) then
                   keep_searching  = .false.
                   headerorder(kj) = ki
                endif
             enddo
          enddo
          !--Get the output filename
          if ( single_output ) then
             if ( i==1 ) then
                ifailed = .false.
                write(evout_new,'(4a)')trim(outprefix),'.',trim(outfilenames(j)),'00.ev'
                open(unit=icolB,file=evout_new)
                write(icolB,ev_fmtH)'#',columnsEV(1:numcol0)
             endif
          else
             ifailed = .false.
             write(evout_new,'(3a,I2.2,a)')trim(outprefix),'.',trim(outfilenames(j)),i,'.ev'
             open(unit=icolB,file=evout_new)
             write(icolB,ev_fmtH)'#',columnsEV(1:numcol0)
          endif
          !--Rewrite the .ev file in the new order with the new name
          open(icolA,file=trim(evfile))
          read(icolA,*) cdummy ! read and ignore header
          neworder = 0.0
          iio      = 0
          do while (iio==0)
             read(icolA,*,iostat=iio) oldorder(1:numcol)
             if (iio==0) then
                do k = 1,numcol
                   if (headerorder(k) > 0) neworder(headerorder(k)) = oldorder(k)
                enddo
                write(icolB,ev_fmtD) neworder(1:numcol0)
             endif
          enddo
          close(icolA)
          if (iio>0) then
             ifailed    = .true.
             readfailed = .true.
             write(*,'(5x,3a,I3)') "ERROR reading data from ",trim(evfile),".  Failed with error code ",iio
          endif
          if (.not.single_output) then
             close(icolB)
             if ( .not. ifailed ) write(*,'(5x,a)') trim(evout_new)
          endif
       endif
    enddo
    if (single_output) then
       close(icolB)
       if ( .not. ifailed ) write(*,'(5x,a)') trim(evout_new)
    endif
 enddo
 if ( readfailed ) then
    write(*,'(a)') ' '
    write(*,'(a)') "Some files were created with no or incomplete data.  Only trust files listed above."
    write(*,'(a)') "Try recompiling after changing 'oldorder' and 'neworder' to kind=16 arrays."
    write(*,'(a)') "We apologise for any inconvenience that this may have caused."
 endif
 write(*,'(a)') ' '
 !
 !--Write inputs to .evin file
 !
 call write_evin_file(outprefix,infilenames(1:nummodels),nummodels,single_output,concise)

end program phantomevcompare
