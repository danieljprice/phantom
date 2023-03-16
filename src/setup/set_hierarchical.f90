!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module sethierarchical
  !
  ! Contains utility routines for setting up generic hierarchical systems
  !
  ! :References:
  !
  ! :Owner: Simone Ceppi
  !
  ! :Runtime parameters: None
  !
  ! :Dependencies: 
  !
  
  implicit none

  public :: set_hierarchical
  public :: set_hierarchical_interactively
  public :: write_hierarchical_setupfile
  public :: read_hierarchical_setupfile
  public :: set_hierarchical_default_options

  public :: get_hierarchical_level_com

  public :: get_hier_level_mass

  private
  
contains
  
  subroutine set_hierarchical(hierarchy,sink_num, sink_labels, hl_labels, hl_num, mass, accr, a, e, inc, O, w, f, &
                              nptmass, xyzmh_ptmass, vxyz_ptmass, ierr)
    use setbinary, only:set_multiple
    character(len=100), intent(in) :: hierarchy
    character(len=10), intent(in) :: sink_labels(:)
    character(len=10), intent(in) :: hl_labels(:)
    integer, intent(out)    :: sink_num, hl_num
    real, intent(in)     :: a(:), e(:), inc(:), O(:), w(:), f(:), mass(:), accr(:)
    real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
    integer, intent(inout) :: nptmass
    integer, intent(out)   :: ierr


    integer :: i
    
    real :: m1, m2, accr1, accr2, binary_a, binary_e, binary_i, binary_O, binary_w, binary_f
    
    integer :: splits, sink_num_temp, hl_index, subst
    character(len=10) :: sink_list(10), split_list(10), hl_temp

    splits = 0
    sink_list = sink_labels
    sink_num_temp = sink_num

    call recursive_splitting(sink_num_temp, sink_list, split_list, splits)
    
    do i=splits+hl_num*0,1,-1
       
       hl_temp = trim(split_list(i))
       
       m1 = get_hier_level_mass(trim(hl_temp)//'1', mass, sink_num, sink_labels)
       m2 = get_hier_level_mass(trim(hl_temp)//'2', mass, sink_num, sink_labels)
       
       if (any(sink_list == trim(hl_temp)//'1')) then
          accr1 = accr(findloc(sink_list,trim(hl_temp)//'1', 1))
       else
          accr1 = 1.
       end if
       
       if (any(sink_list == trim(hl_temp)//'2')) then
          accr2 = accr(findloc(sink_list,trim(hl_temp)//'2', 1))
       else
          accr2 = 1.
       end if
       
       hl_index = findloc(hl_labels, trim(hl_temp), 1)
       
       binary_a =  a(hl_index)
       binary_e = e(hl_index)
       binary_O = O(hl_index)
       binary_w = w(hl_index)
       binary_i = inc(hl_index)
       binary_f = f(hl_index)
       
       read(hl_temp,*,iostat=subst) subst
       
       call set_multiple(m1,m2,semimajoraxis=binary_a,eccentricity=binary_e, &
            posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
            f=binary_f,accretion_radius1=accr1,accretion_radius2=accr2, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr, subst=subst)
       
    end do
  end subroutine set_hierarchical
  
  
  
  subroutine set_hierarchical_interactively(hierarchy, sink_num, sink_labels, hl_labels, hl_num)
    use prompting, only:prompt

    character(len=100), intent(out) :: hierarchy
    character(len=10), intent(out) :: sink_labels(:)
    character(len=10), intent(out) :: hl_labels(:)
    integer, intent(out)    :: sink_num, hl_num


    
    call prompt('What is the hierarchy?',hierarchy)
    
    call process_hierarchy(hierarchy,sink_num, sink_labels, hl_labels, hl_num) 

  end subroutine set_hierarchical_interactively
  
  
  
  subroutine write_hierarchical_setupfile(hierarchy, iunit, sink_num, sink_labels, hl_labels, hl_num, &
                                          mass, accr, a, e, inc, O, w, f)
    use infile_utils, only:write_inopt
    character(len=100), intent(in) :: hierarchy
    integer, intent(in) :: iunit
    integer, intent(out) :: sink_num, hl_num
    character(len=10), intent(out) :: sink_labels(:), hl_labels(:)
    real, intent(in)     :: a(:), e(:), inc(:), O(:), w(:), f(:), mass(:), accr(:)

    integer :: i

    write(iunit,"(/,a)") '# options for hierarchical system'
      
    call write_inopt(hierarchy, 'hier','', iunit)
    
    call process_hierarchy(hierarchy,sink_num, sink_labels, hl_labels, hl_num)
    
    write(iunit,"(/,a)") '### sink properties'
    do i=1,sink_num
       call write_inopt(mass(i), trim(sink_labels(i))//'_mass','', iunit)
       call write_inopt(accr(i), trim(sink_labels(i))//'_accr','', iunit)
    end do
    
    write(iunit,"(/,a)") '### orbit properties'
    
    do i=1,hl_num
       call write_inopt(a(i), trim(hl_labels(i))//'_a','',iunit)
       call write_inopt(e(i), trim(hl_labels(i))//'_e','',iunit)
       call write_inopt(inc(i), trim(hl_labels(i))//'_i','',iunit)
       call write_inopt(O(i), trim(hl_labels(i))//'_O','',iunit)
       call write_inopt(w(i), trim(hl_labels(i))//'_w','',iunit)
       call write_inopt(f(i), trim(hl_labels(i))//'_f','',iunit)
    end do
    
    
  end subroutine write_hierarchical_setupfile
  
  
  subroutine read_hierarchical_setupfile(hierarchy, db, nerr, sink_num, sink_labels, hl_labels, hl_num, &
                                         mass, accr, a, e, inc, O, w, f)
    
    use infile_utils, only:read_inopt, inopts
    character(len=100), intent(out) :: hierarchy
    type(inopts), allocatable, intent(inout) :: db(:)
    integer, intent(inout) :: nerr
    integer, intent(out) :: sink_num, hl_num
    character(len=10), intent(out) :: sink_labels(:), hl_labels(:)
    real, intent(out)     :: a(:), e(:), inc(:), O(:), w(:), f(:), mass(:), accr(:)

    integer :: i
  
    call read_inopt(hierarchy,'hier',db,errcount=nerr)
      
    call process_hierarchy(hierarchy,sink_num, sink_labels, hl_labels, hl_num)
    
    do i=1,sink_num
       call read_inopt(mass(i), trim(sink_labels(i))//'_mass',db,errcount=nerr)
       call read_inopt(accr(i), trim(sink_labels(i))//'_accr',db,errcount=nerr)
    end do
    
    do i=1,hl_num
       call read_inopt(a(i), trim(hl_labels(i))//'_a',db,errcount=nerr)
       call read_inopt(e(i), trim(hl_labels(i))//'_e',db,errcount=nerr)
       call read_inopt(inc(i), trim(hl_labels(i))//'_i',db,errcount=nerr)
       call read_inopt(O(i), trim(hl_labels(i))//'_O',db,errcount=nerr)
       call read_inopt(w(i), trim(hl_labels(i))//'_w',db,errcount=nerr)
       call read_inopt(f(i), trim(hl_labels(i))//'_f',db,errcount=nerr)
    end do
  end subroutine read_hierarchical_setupfile


  subroutine set_hierarchical_default_options(hierarchy, sink_num, sink_labels, hl_labels, hl_num, &
                                              mass, accr, a, e, inc, O, w, f)
    character(len=100), intent(in) :: hierarchy
    character(len=10), intent(out) :: sink_labels(:)
    character(len=10), intent(out) :: hl_labels(:)
    integer, intent(out)    :: sink_num, hl_num
    real, intent(out) :: mass(:), accr(:), a(:), e(:), inc(:), O(:), w(:), f(:)

    integer :: i

    call process_hierarchy(hierarchy,sink_num, sink_labels, hl_labels, hl_num)

    do i=1,sink_num
       mass(i)=1.
       accr(i)=1.
    end do
    
    do i=1,hl_num
       a(i)=10000./10**len(trim(hl_labels(i)))
       e(i)=0.
       inc(i)=0.
       O(i)=0.
       w(i)=0.
       f(i)=0.
    end do

  end subroutine set_hierarchical_default_options
  
  
  !--------------------------------------------------------------------------
  !
  ! Process the .setup hierarchy string to extract information for building the system 
  !
  !--------------------------------------------------------------------------
  subroutine process_hierarchy(hierarchy, sink_num, sink_labels, hl_labels, hl_num)
    character(len=100), intent(in) :: hierarchy
    character(len=10), intent(out) :: sink_labels(:)
    character(len=10), intent(out) :: hl_labels(:)
    integer, intent(out)    :: sink_num, hl_num
    
    integer :: i,j, del_pos, pre_del_pos
    character(len=:), allocatable :: sink, temp_hl
    
    sink_labels = '          '
    hl_labels = '          '
    
    
    pre_del_pos = 0
    sink_num = 0
    
    ! count the number of sinks in the system
    del_pos = scan(hierarchy, ',')
    do while (del_pos > pre_del_pos)
       sink_labels(sink_num+1) = hierarchy(pre_del_pos+1:del_pos-1)
       sink_num=sink_num+1
       pre_del_pos = del_pos
       del_pos = scan(hierarchy(pre_del_pos+1:), ',')+pre_del_pos
    end do
    
    sink_labels(sink_num+1) = hierarchy(pre_del_pos+1:)
    sink_num = sink_num+1
    
    hl_num = 0
    do i=1,sink_num
       sink = trim(sink_labels(i))
       del_pos = len(trim(sink))!-1
       do j=1,del_pos!, 1, -1
          temp_hl = trim(sink(1:j))
          if ( .not. any(hl_labels == trim(temp_hl))) then
             if ( any(sink_labels == trim(temp_hl)) .or. any(sink_labels == trim(temp_hl))) then
                ! add to discs
             else
                hl_num=hl_num+1
                hl_labels(hl_num) = trim(temp_hl)
             end if
          end if
       end do
    end do
    
  end subroutine process_hierarchy
  
  !--------------------------------------------------------------------------
  !
  ! Reverse the splitting process from sink_labels to build the system
  !
  !--------------------------------------------------------------------------
  subroutine recursive_splitting(sink_num, sink_list, split_list, splits)
    character(len=10), intent(in) :: sink_list(:)
    character(len=10), intent(inout) :: split_list(:)
    integer, intent(inout)    :: splits
    integer, intent(in) ::sink_num
    
    integer :: i, longest, longests_len, count
    character(len=10) :: longests(10), new_splits(10), new_sink_list(10), longest_cut, sink_list_temp(10)

    sink_list_temp = sink_list
    
    longest = 0
    if (sink_num>1) then 
       ! Find the longest
       do i=1,sink_num
          if (longest < len(trim(sink_list_temp(i)))) then
             longest = len(trim(sink_list_temp(i)))
          end if
       end do
       
       ! Select the longests and cut them
       longests_len = 0
       do i=1,sink_num
          if (len(trim(sink_list_temp(i))) == longest) then
             longests_len = longests_len+1
             longests(longests_len) = trim(sink_list_temp(i))
             
             sink_list_temp(i) = sink_list_temp(i)(:len(trim(longests(longests_len)))-1)!//' '
          end if
       end do
       
       
       ! Cut the longest and add to split list with no doubles
       count=0
       do i=1,longests_len
          longest_cut = longests(i)
          if (.not. any(new_splits == longest_cut(:longest-1))) then
             count = count + 1
             new_splits(count) = longests(i)(1:longest-1)
          end if
       end do
       
       ! Add new splits to split_list
       do i=splits+1, splits+count
          split_list(i) = new_splits(i-splits)
       end do
       splits = splits + count
       
       ! Clean sink_list_temp from doubles
       
       count=0
       do i=1,sink_num
          if (.not. any(new_sink_list == sink_list_temp(i))) then
             count = count + 1
             new_sink_list(count) = sink_list_temp(i)
          end if
       end do
              
       call recursive_splitting(count, new_sink_list(:count), split_list(:splits), splits)
       
    end if
  end subroutine recursive_splitting
  
  !--------------------------------------------------------------------------
  !
  ! Compute the total mass of an hierarchical level
  !
  !--------------------------------------------------------------------------
  real function get_hier_level_mass(level, mass, sink_num, sink_list)
    character(*), intent(in) :: level
    character(len=10), intent(in) :: sink_list(:)
    real, intent(in) :: mass(:)
    integer, intent(in) :: sink_num
    
    real :: part_mass
    integer :: i
    
    part_mass = 0
    do i=1, sink_num
       if ((len(trim(sink_list(i))) >= len(trim(level))) .and. (sink_list(i)(:len(trim(level))) == trim(level))) then
          part_mass = part_mass + mass(i)
       end if
    end do

    get_hier_level_mass = part_mass
    
  end function get_hier_level_mass


  !--------------------------------------------------------------------------
  !
  ! Compute position and velocity of hierarchical level
  !
  !--------------------------------------------------------------------------
  subroutine get_hierarchical_level_com(level, xorigin, vorigin, xyzmh_ptmass, vxyz_ptmass)
    character(*), intent(in) :: level
    real, intent(out) :: xorigin(:), vorigin(:)
    real, intent(in) :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)

    integer :: int_sinks(10), inner_sinks_num, i
    real :: mass

    call read_HIERARCHY_index(level, int_sinks, inner_sinks_num)

    xorigin = 0.
    vorigin = 0.
    mass = 0.
    
    do i=1,inner_sinks_num
       xorigin(:) = xorigin(:)+xyzmh_ptmass(4,int_sinks(i))*xyzmh_ptmass(1:3,int_sinks(i))
       vorigin(:) = vorigin(:)+xyzmh_ptmass(4,int_sinks(i))*vxyz_ptmass(1:3,int_sinks(i))
       mass = mass + xyzmh_ptmass(4,int_sinks(i))
    end do
    xorigin = xorigin/mass
    vorigin = vorigin/mass


  end subroutine get_hierarchical_level_com

  !--------------------------------------------------------------------------
  !
  ! Retrieve sink position in HIERARCHY file
  !
  !--------------------------------------------------------------------------
  subroutine read_HIERARCHY_index(level, int_sinks, inner_sinks_num)
    character(*), intent(in) :: level
    integer, intent(out) :: inner_sinks_num
    integer, intent(out)                :: int_sinks(10)
    
    real, dimension(24,10) :: data
    integer                :: i, io, lines, ierr, h_index
    logical                :: iexist



    character(len=10)      :: label = '         '

    read(level, *, iostat=h_index) h_index

    inquire(file='HIERARCHY', exist=iexist)
 
    if (iexist) then
       open(1, file = 'HIERARCHY', status = 'old')
       lines=0
       do
          read(1, *, iostat=io) data(lines+1,:)
          if (io/=0) exit
          lines = lines + 1
       enddo
       close(1)
    else
       print "(1x,a)",'ERROR: set_multiple: there is no HIERARCHY file, cannot perform subtitution.'
       ierr = 100!ierr_HIER2
    endif

    inner_sinks_num = 0
    do i=1, lines
       write(label, '(i0)') int(data(i,2))
       if (data(i,1) > 0 .and. (len(trim(label)) >= len(trim(level))) .and. (label(:len(trim(level))) == trim(level))) then
          inner_sinks_num = inner_sinks_num+1
          int_sinks(inner_sinks_num) = int(data(i,1))
       end if
    end do

  end subroutine read_HIERARCHY_index


  
  
end module sethierarchical
