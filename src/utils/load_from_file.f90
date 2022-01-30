!--------------------------------------------------------------
!   This module provides routines to load txt files as matrices
!   Usage:
!
!   use load_from_file, only:load_data_file
!   real, dimension(:,:), allocatable :: sigmadat
!   call load_data_file('sigma_density.dat',sigmadat)
!------------------------------------------------------------

module load_from_file

 implicit none

 public :: load_data_file,write_in_file

 interface write_in_file
    module procedure  write_in_file_1d,write_in_file_2d,write_in_file_1dx2
 end interface

 private
 
 contains

 subroutine load_data_file(namefile,datafile,nhead)
 character(len=*), intent(in) :: namefile
 integer, intent(in), optional :: nhead
 integer           :: nrows,mcolumns,nheadlines,iunit,ierr,i
 character :: c
 real, dimension(:,:), allocatable, intent(inout) :: datafile

 ierr=0
 iunit=1
 !--specify number of headlines
 nheadlines=1
 if(present(nhead)) then
    nheadlines=nhead
 endif
 write(*,*) 'Skipping ',nheadlines,' head lines'

 open(unit=iunit,file=namefile,status='old',action='read')

 mcolumns=number_of_columns(iunit,nheadlines)
 nrows=number_of_rows(iunit)
 write(*,*) 'Loading data from file:',namefile
 write(*,*) 'Found nrows, mcolumns:',nrows,mcolumns

 allocate(datafile(nrows-nheadlines,mcolumns))
 do i=1,nheadlines
    read(iunit,*) c
 enddo
 do i=1, nrows-nheadlines
    read(iunit,*) datafile(i,:) 
 enddo
 close(iunit)

 end subroutine load_data_file

 subroutine write_in_file_1d(namefile,arraytowrite)
    character(len=*), intent(in) :: namefile
    real, intent(in), dimension(:) :: arraytowrite
    integer :: iunit,i

    iunit=1
 
    open(unit=iunit,file=namefile,status='replace',action='write')
    do i=1,size(arraytowrite(:))
       write(iunit,*) arraytowrite(i)
    enddo
    close(unit=iunit)

 end subroutine write_in_file_1d    
 

subroutine write_in_file_2d(namefile,arraytowrite)
    character(len=*), intent(in) :: namefile
    real, intent(in), dimension(:,:) :: arraytowrite
    integer :: iunit,i
 
    iunit=1
 
    open(unit=iunit,file=namefile,status='replace',action='write')

    do i=1,size(arraytowrite(:,1))
       write(iunit,*) arraytowrite(i,:)
    enddo

    close(unit=iunit)

 end subroutine write_in_file_2d    

 subroutine write_in_file_1dx2(namefile,arraytowrite1,arraytowrite2)
    character(len=*), intent(in) :: namefile
    real, intent(in), dimension(:) :: arraytowrite1,arraytowrite2
    integer :: iunit,i

    iunit=1 

    open(unit=iunit,file=namefile,status='replace',action='write')
    do i=1,size(arraytowrite1(:))
       write(iunit,*) arraytowrite1(i),arraytowrite2(i)
    enddo
    close(unit=iunit)

 end subroutine write_in_file_1dx2    
 

 integer function number_of_columns(s,nhead)
    !! version: experimental
    !!
    !! determine number of columns
    integer,intent(in) :: s,nhead

    integer :: ios,i
    character :: c
    logical :: lastblank

    rewind(s)
    number_of_columns = 0
    lastblank = .true.
    !--skipping head lines for calculating mcolumns
    do i=1,nhead
       read(s,*) c
    enddo

    do
      read(s, '(a)', advance='no', iostat=ios) c
      if (ios /= 0) exit
      ! basically it counts number of free spaces between between the characters in the first line 
      if (lastblank .and. (.not. is_blank(c))) number_of_columns = number_of_columns + 1
      lastblank = is_blank(c)
    end do
    rewind(s)

 end function number_of_columns


 integer function number_of_rows(s) result(nrows)
    !! version: experimental
    !!
    !! determine number or rows
    integer,intent(in)::s

    integer :: ios
    character  :: r

    rewind(s)
    nrows = 0
    do
      read(s, *, iostat=ios) r
      if (ios /= 0) exit
      nrows = nrows + 1
    end do

    rewind(s)

end function number_of_rows
 
 pure logical function is_blank(c)
        character(len=1), intent(in) :: c !! The character to test.
        integer :: ic
        ic = iachar(c)             ! TAB
        is_blank = (c == ' ') .or. (ic == int(z'09'));
 end function

end module load_from_file

