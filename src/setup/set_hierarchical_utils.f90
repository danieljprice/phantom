!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module sethier_utils
!
! Utility routines needed for the hierarchical sink particle setup
!
! :References:
!  - Ceppi et al. (2022) MNRAS 514, 906
!  - Ceppi et al. (2023) MNRAS 520, 5817
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!

 implicit none
 integer, parameter, public :: max_hier_levels=10
 ! one row per sink and per intermediate level of the hierarchy
 integer, parameter, public :: hier_db_size=2*max_hier_levels, hier_db_prop=10
 integer, parameter, public :: lenhierstring = 100
 real, parameter, public :: pi = 4.*atan(1.)

 integer, parameter, public :: ierr_HIER2      = 6, &
                               ierr_subststar  = 7, &
                               ierr_Omegasubst = 8, &
                               ierr_missstar   = 9

 type hs_labels
    integer :: sink_num, hl_num
    character(len=10) :: sink(max_hier_levels), hl(max_hier_levels)
 end type hs_labels

 type hs_orbital_elements
    real :: a,e,O,w,inc,f
 end type hs_orbital_elements

 type hs_sink_properties
    real :: mass, accr
 end type hs_sink_properties

 type hierarchical_system
    type(hs_labels) :: labels
    type(hs_orbital_elements) :: levels(max_hier_levels)
    type(hs_sink_properties) :: sinks(max_hier_levels)
 end type hierarchical_system

 public :: hierarchical_system,process_hierarchy,count_hierarchy_labels
 public :: recursive_splitting, get_hierarchical_level_mass
 public :: load_hierarchy_file,update_hierarchy_file,check_substitution
 public :: find_hierarchy_index,find_hier_level_orb_elem,find_data_index
 public :: gen_rotate
 public :: find_ptmass_index
 public :: findloc_local
 public :: hierarchy_from_brackets,hierarchy_to_brackets,print_hierarchy
 public :: write_hierarchy_file,hierarchy_filename

 interface findloc_local
  module procedure findloc_int
  module procedure findloc_string
  module procedure findloc_logical
 end interface findloc_local

 private

contains

!--------------------------------------------------------------------------
!+
!  Process the .setup hierarchy string to extract information
!  for building the system
!+
!--------------------------------------------------------------------------
pure function process_hierarchy(string) result(mylabels)
 character(len=*), intent(in) :: string
 integer :: i,j,del_pos,pre_del_pos
 type(hs_labels) :: mylabels
 character(len=:), allocatable :: sink,temp_hl

 mylabels%sink = ''
 mylabels%hl = ''

 pre_del_pos = 0
 mylabels%sink_num = 0

 ! count the number of sinks in the system (labels beyond
 ! max_hier_levels are ignored: see count_hierarchy_labels)
 del_pos = scan(string, ',')
 do while (del_pos > pre_del_pos .and. mylabels%sink_num < max_hier_levels-1)
    mylabels%sink(mylabels%sink_num+1) = string(pre_del_pos+1:del_pos-1)
    mylabels%sink_num=mylabels%sink_num+1
    pre_del_pos = del_pos
    del_pos = scan(string(pre_del_pos+1:), ',')+pre_del_pos
 enddo

 mylabels%sink(mylabels%sink_num+1) = string(pre_del_pos+1:)
 mylabels%sink_num = mylabels%sink_num+1

 mylabels%hl_num = 0
 do i=1,mylabels%sink_num
    sink = trim(mylabels%sink(i))
    del_pos = len(trim(sink))
    do j=1,del_pos
       temp_hl = trim(sink(1:j))
       if ( .not. any(mylabels%hl == trim(temp_hl))) then
          if ( any(mylabels%sink == trim(temp_hl)) .or. any(mylabels%sink == trim(temp_hl))) then
             ! add to discs
          else
             mylabels%hl_num=mylabels%hl_num+1
             mylabels%hl(mylabels%hl_num) = trim(temp_hl)
          endif
       endif
    enddo
 enddo

end function process_hierarchy

!--------------------------------------------------------------------------
!+
!  number of labels (comma-separated) in a hierarchy string
!+
!--------------------------------------------------------------------------
pure integer function count_hierarchy_labels(string) result(n)
 character(len=*), intent(in) :: string
 integer :: i

 n = 1
 do i=1,len_trim(string)
    if (string(i:i) == ',') n = n + 1
 enddo

end function count_hierarchy_labels

!--------------------------------------------------------------------------
!+
!  reverse the splitting process from sink_labels to build the system
!+
!--------------------------------------------------------------------------
pure recursive subroutine recursive_splitting(sink_num, sink_list, split_list, splits)
 character(len=10), intent(in)    :: sink_list(:)
 integer,           intent(in)    :: sink_num
 character(len=10), intent(inout) :: split_list(:)
 integer,           intent(inout) :: splits
 integer :: i, longest, longests_len, count
 character(len=10) :: longests(max_hier_levels), new_splits(max_hier_levels)
 character(len=10) :: new_sink_list(max_hier_levels), longest_cut, sink_list_temp(max_hier_levels)

 sink_list_temp(:) = sink_list(:)
 ! these are searched before being filled, so must start empty
 new_splits(:)    = ''
 new_sink_list(:) = ''

 longest = 0
 if (sink_num>1) then
    ! Find the longest
    do i=1,sink_num
       if (longest < len(trim(sink_list_temp(i)))) then
          longest = len(trim(sink_list_temp(i)))
       endif
    enddo

    ! Select the longests and cut them
    longests_len = 0
    do i=1,sink_num
       if (len(trim(sink_list_temp(i))) == longest) then
          longests_len = longests_len+1
          longests(longests_len) = trim(sink_list_temp(i))

          sink_list_temp(i) = sink_list_temp(i)(:len(trim(longests(longests_len)))-1)
       endif
    enddo

    ! Cut the longest and add to split list with no doubles
    count=0
    do i=1,longests_len
       longest_cut = longests(i)
       if (.not. any(new_splits == longest_cut(:longest-1))) then
          count = count + 1
          new_splits(count) = longests(i)(1:longest-1)
       endif
    enddo

    ! Add new splits to split_list
    do i=splits+1, min(splits+count,size(split_list))
       split_list(i) = new_splits(i-splits)
    enddo
    splits = splits + count

    ! Clean sink_list_temp from doubles
    count=0
    do i=1,sink_num
       if (.not. any(new_sink_list == sink_list_temp(i))) then
          count = count + 1
          new_sink_list(count) = sink_list_temp(i)
       endif
    enddo

    call recursive_splitting(count, new_sink_list(:), split_list(:), splits)

 endif

end subroutine recursive_splitting

!--------------------------------------------------------------------------
!+
!  compute the total mass of an hierarchical level
!+
!--------------------------------------------------------------------------
real function get_hierarchical_level_mass(level, hs) result(part_mass)
 character(len=*),          intent(in) :: level
 type(hierarchical_system), intent(in) :: hs
 integer :: i

 part_mass = 0
 do i=1, hs%labels%sink_num
    !print*,'-',trim(hs%labels%sink(i)),'-', '   ', '-',trim(level),'- ', len(trim(level))
    if ((len(trim(hs%labels%sink(i))) >= len(trim(level))) .and. &
         ((hs%labels%sink(i)(:len(trim(level))) == trim(level)))) then
       part_mass = part_mass + hs%sinks(i)%mass
    endif
 enddo

end function get_hierarchical_level_mass

!--------------------------------------------------------------------------
!+
!  read information from the hierarchy file
!+
!--------------------------------------------------------------------------
subroutine load_hierarchy_file(filename, data, lines, ierr)
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr,lines
 real,             intent(out) :: data(hier_db_size,hier_db_prop)
 integer :: io
 logical :: iexist

 ierr = 0
 inquire(file=trim(filename), exist=iexist)

 if (iexist) then
    open(2,file=trim(filename),status='old')
    lines=0
    do while (lines < hier_db_size)
       read(2, *,iostat=io) data(lines+1,:)
       if (io/=0) exit
       lines = lines + 1
    enddo
    close(2)
 else
    lines = 0
    print "(1x,a)",'ERROR: set_multiple: there is no HIERARCHY file, cannot perform subtitution.'
    ierr = ierr_HIER2
 endif
end subroutine load_hierarchy_file

!--------------------------------------------------------------------------
!+
!  write information to the hierarchy file
!+
!--------------------------------------------------------------------------
subroutine update_hierarchy_file(filename, hs, data, lines, hier_prefix, i1, i2, ierr)
 integer,                   intent(in)    :: i1, i2
 character(len=20),         intent(in)    :: hier_prefix
 type(hierarchical_system), intent(in)    :: hs
 real,                      intent(inout) :: data(hier_db_size,hier_db_prop)
 integer,                   intent(inout) :: lines
 integer,                   intent(out)   :: ierr
 character(len=*),          intent(in)    :: filename
 integer :: iu
 real :: mprimary, msecondary, semimajoraxis, eccentricity, incl, arg_peri, posang_ascnode, binary_f
 real :: accr1, accr2, period
 logical :: iexist

 !print*,'updating: ',trim(adjustl(hier_prefix))

 call find_hier_level_orb_elem(hier_prefix, hs, mprimary, msecondary, accr1, accr2, &
       semimajoraxis, eccentricity, incl, posang_ascnode, &
       arg_peri, binary_f)

 !print*,'elements of ', hier_prefix, ' : ', mprimary, msecondary, accr1, accr2, &
 !     semimajoraxis, eccentricity, incl, posang_ascnode, &
 !     arg_peri, binary_f

 period = sqrt(4.*pi**2*semimajoraxis**3/(mprimary+msecondary))

 if (lines > 0) then
    call write_hierarchy_file(filename,data,lines)
    open(newunit=iu,file=trim(filename),status='old',position='append')
 else
    inquire(file=trim(filename), exist=iexist)
    if (iexist) print "(1x,a)",'WARNING: set_multiple: deleting an existing HIERARCHY file.'
    open(newunit=iu,file=trim(filename),status='replace')
 endif
 write(iu,*) i1, trim(hier_prefix)//"1", mprimary, msecondary, semimajoraxis, eccentricity, &
       period, incl, arg_peri, posang_ascnode
 write(iu,*) i2, trim(hier_prefix)//"2", msecondary, mprimary, semimajoraxis, eccentricity, &
       period, incl, arg_peri, posang_ascnode
 close(iu)

 call load_hierarchy_file(filename, data, lines, ierr)

end subroutine update_hierarchy_file

!--------------------------------------------------------------------------
!+
!  name of the hierarchy file: prefix.hierarchy, or HIERARCHY if no
!  prefix is given
!+
!--------------------------------------------------------------------------
function hierarchy_filename(prefix) result(filename)
 character(len=20), intent(in), optional :: prefix
 character(len=30) :: filename ! prefix (len 20) + '.hierarchy'

 if (present(prefix)) then
    filename = trim(prefix)//'.hierarchy'
 else
    filename = 'HIERARCHY'
 endif

end function hierarchy_filename

!--------------------------------------------------------------------------
!+
!  write the whole table of the hierarchy file
!+
!--------------------------------------------------------------------------
subroutine write_hierarchy_file(filename, data, lines)
 character(len=*), intent(in) :: filename
 real,             intent(in) :: data(hier_db_size,hier_db_prop)
 integer,          intent(in) :: lines
 integer :: i,iu

 open(newunit=iu,file=trim(filename),status='replace')
 do i=1,lines
    write(iu,*) int(data(i,1)), int(data(i,2)), data(i,3:)
 enddo
 close(iu)

end subroutine write_hierarchy_file

!--------------------------------------------------------------------------
!+
!  convert a hierarchy given with brackets, where o is a star and
!  (A-B) is a binary of the components A and B, e.g. (o-o)-o, into
!  the numbered labels of the stars, e.g. 111,112,12. The outermost
!  pair needs no brackets. A single star, o, is labelled 1
!+
!--------------------------------------------------------------------------
subroutine hierarchy_from_brackets(brackets,hierarchy,ierr)
 character(len=*), intent(in)  :: brackets
 character(len=*), intent(out) :: hierarchy
 integer,          intent(out) :: ierr
 character(len=len(brackets)+2) :: str
 integer :: i,n,pos

 ! remove spaces
 str = ''
 n = 0
 do i=1,len_trim(brackets)
    if (brackets(i:i) /= ' ') then
       n = n + 1
       str(n:n) = brackets(i:i)
    endif
 enddo

 ! first try the string as a single node, e.g. o or (o-o)
 pos = 1
 hierarchy = ''
 call parse_bracket_node(str(1:n),pos,'1',hierarchy,ierr)
 if (ierr == 0 .and. pos == n+1) return

 ! otherwise the outermost pair has no brackets, e.g. (o-o)-o
 str = '('//str(1:n)//')'
 pos = 1
 hierarchy = ''
 call parse_bracket_node(str(1:n+2),pos,'1',hierarchy,ierr)
 if (ierr == 0 .and. pos /= n+3) ierr = 1

end subroutine hierarchy_from_brackets

!--------------------------------------------------------------------------
!+
!  parse one node of a bracketed hierarchy starting at str(pos:),
!  where the node has the numbered label given: a star, o, is added
!  to the list of labels; a binary (A-B) gives its components the
!  labels label//'1' and label//'2'
!+
!--------------------------------------------------------------------------
recursive subroutine parse_bracket_node(str,pos,label,hierarchy,ierr)
 character(len=*), intent(in)    :: str,label
 integer,          intent(inout) :: pos
 character(len=*), intent(inout) :: hierarchy
 integer,          intent(out)   :: ierr

 ierr = 1
 if (pos > len(str) .or. len(label) > len(hierarchy) .or. len(label) > 10) return
 select case(str(pos:pos))
 case('o')
    if (len_trim(hierarchy) + len(label) + 1 > len(hierarchy)) return
    if (len_trim(hierarchy) > 0) then
       hierarchy = trim(hierarchy)//','//label
    else
       hierarchy = label
    endif
    pos = pos + 1
    ierr = 0
 case('(')
    pos = pos + 1
    call parse_bracket_node(str,pos,label//'1',hierarchy,ierr)
    if (ierr /= 0) return
    ierr = 1
    if (pos > len(str)) return
    if (str(pos:pos) /= '-') return
    pos = pos + 1
    call parse_bracket_node(str,pos,label//'2',hierarchy,ierr)
    if (ierr /= 0) return
    ierr = 1
    if (pos > len(str)) return
    if (str(pos:pos) /= ')') return
    pos = pos + 1
    ierr = 0
 end select

end subroutine parse_bracket_node

!--------------------------------------------------------------------------
!+
!  convert the numbered labels of the stars, e.g. 111,112,12, into a
!  hierarchy given with brackets, e.g. (o-o)-o. The labels must form
!  a complete binary tree: every binary has both components, and no
!  star is also a binary. The brackets list the stars from left to
!  right, which may differ from the order of the numbered labels
!+
!--------------------------------------------------------------------------
subroutine hierarchy_to_brackets(hierarchy,brackets,ierr)
 character(len=*), intent(in)  :: hierarchy
 character(len=*), intent(out) :: brackets
 integer,          intent(out) :: ierr
 type(hs_labels) :: labels
 integer :: i,nfound

 labels = process_hierarchy(hierarchy)
 brackets = ''
 ierr = 1
 ! labels contain only 1 and 2 and start from the top level, 1
 do i=1,labels%sink_num
    if (len_trim(labels%sink(i)) < 1) return
    if (verify(trim(labels%sink(i)),'12') /= 0 .or. labels%sink(i)(1:1) /= '1') return
 enddo

 nfound = 0
 call write_bracket_node('1',labels%sink(1:labels%sink_num),brackets,nfound,ierr)
 if (ierr /= 0) return
 ! every label must have been found exactly once
 if (nfound /= labels%sink_num) then
    ierr = 1
    return
 endif
 ! the outermost pair needs no brackets (blanked and shifted left,
 ! rather than copying a slice of the string onto itself)
 i = len_trim(brackets)
 if (i > 2) then
    brackets(1:1) = ' '
    brackets(i:i) = ' '
    brackets = adjustl(brackets)
 endif

end subroutine hierarchy_to_brackets

!--------------------------------------------------------------------------
!+
!  write the node with the numbered label given, as o if it is a star
!  or as (A-B) if it is a binary
!+
!--------------------------------------------------------------------------
recursive subroutine write_bracket_node(label,sinks,brackets,nfound,ierr)
 character(len=*), intent(in)    :: label
 character(len=*), intent(in)    :: sinks(:)
 character(len=*), intent(inout) :: brackets
 integer,          intent(inout) :: nfound
 integer,          intent(out)   :: ierr
 integer :: i,nstar,nbelow

 ierr = 1
 if (len(label) > 10) return
 nstar = 0
 nbelow = 0
 do i=1,size(sinks)
    if (trim(sinks(i)) == label) then
       nstar = nstar + 1
    elseif (len_trim(sinks(i)) > len(label)) then
       if (sinks(i)(1:len(label)) == label) nbelow = nbelow + 1
    endif
 enddo
 if (nstar == 1 .and. nbelow == 0) then
    ! a star
    brackets = trim(brackets)//'o'
    nfound = nfound + 1
    ierr = 0
 elseif (nstar == 0 .and. nbelow >= 2) then
    ! a binary
    brackets = trim(brackets)//'('
    call write_bracket_node(label//'1',sinks,brackets,nfound,ierr)
    if (ierr /= 0) return
    brackets = trim(brackets)//'-'
    call write_bracket_node(label//'2',sinks,brackets,nfound,ierr)
    if (ierr /= 0) return
    brackets = trim(brackets)//')'
 endif

end subroutine write_bracket_node

!--------------------------------------------------------------------------
!+
!  print a hierarchy given by the numbered labels of the stars, both
!  with brackets, e.g. (o-o)-o, and as a tree drawn from the top level
!  down, e.g. for 11,121,122:
!
!           __ 1 __
!         /         \
!       11           12
!       |           /  \
!       11        121  122
!       1          2    3      <- sink id
!
!  where stars are carried down to the bottom row, below which the sink
!  number of each star (its position in the hierarchy string) is given
!+
!--------------------------------------------------------------------------
subroutine print_hierarchy(hierarchy)
 character(len=*), intent(in) :: hierarchy
 integer, parameter :: margin = 3, lenrow = 200
 type(hs_labels) :: labels
 character(len=lenhierstring) :: brackets
 character(len=10) :: leaves(max_hier_levels),tmp
 character(len=lenrow) :: row
 character(len=12) :: num
 integer :: ierr,i,j,n,nleaves,maxdepth,spacing,depth,x,xl,xr,isink,i1,i2,nul
 integer :: xleaf(max_hier_levels)

 call hierarchy_to_brackets(hierarchy,brackets,ierr)
 if (ierr /= 0) then
    print "(a)",' WARNING: hierarchy '//trim(hierarchy)//' is not a complete binary tree'
    return
 endif
 labels = process_hierarchy(hierarchy)
 print "(a)",' hierarchy (o = star): '//trim(brackets)

 ! stars from left to right: for labels made of 1 and 2 this is alphabetical order
 nleaves = labels%sink_num
 leaves(1:nleaves) = labels%sink(1:nleaves)
 do i=2,nleaves
    tmp = leaves(i)
    j = i - 1
    do while (j >= 1)
       if (leaves(j) <= tmp) exit
       leaves(j+1) = leaves(j)
       j = j - 1
    enddo
    leaves(j+1) = tmp
 enddo

 ! one column per star, wide enough for the longest label
 maxdepth = 0
 do i=1,nleaves
    maxdepth = max(maxdepth,len_trim(leaves(i))-1)
 enddo
 spacing = maxdepth + 4
 do i=1,nleaves
    xleaf(i) = margin + spacing/2 + (i-1)*spacing
 enddo

 print "(a)",''
 do depth=0,maxdepth
    ! labels of the binaries and stars at this depth, with the stars above carried down
    row = ''
    do n=1,labels%hl_num
       if (len_trim(labels%hl(n))-1 == depth) then
          x  = node_column(trim(labels%hl(n)),leaves(1:nleaves),xleaf)
          xl = (x + node_column(trim(labels%hl(n))//'1',leaves(1:nleaves),xleaf))/2
          xr = (x + node_column(trim(labels%hl(n))//'2',leaves(1:nleaves),xleaf))/2
          ! underline from above each branch to one space short of the label, as in __ 1 __
          i1 = x - (len_trim(labels%hl(n))-1)/2
          i2 = i1 + len_trim(labels%hl(n)) - 1
          ! the same length on both sides, as long as the shorter side allows
          nul = min(i1-2-xl,xr-2-i2)
          if (nul > 0) then
             row(i1-1-nul:i1-2) = repeat('_',nul)
             row(i2+2:i2+1+nul) = repeat('_',nul)
          endif
          call put_label(row,trim(labels%hl(n)),x)
       endif
    enddo
    do i=1,nleaves
       if (len_trim(leaves(i))-1 <= depth) call put_label(row,trim(leaves(i)),xleaf(i))
    enddo
    print "(a)",trim(row)
    if (depth == maxdepth) exit

    ! branches down to the next depth
    row = ''
    do n=1,labels%hl_num
       if (len_trim(labels%hl(n))-1 == depth) then
          x  = node_column(trim(labels%hl(n)),leaves(1:nleaves),xleaf)
          xl = (x + node_column(trim(labels%hl(n))//'1',leaves(1:nleaves),xleaf))/2
          xr = (x + node_column(trim(labels%hl(n))//'2',leaves(1:nleaves),xleaf))/2
          row(xl:xl) = '/'
          row(xr:xr) = '\'
       endif
    enddo
    do i=1,nleaves
       if (len_trim(leaves(i))-1 <= depth) row(xleaf(i):xleaf(i)) = '|'
    enddo
    print "(a)",trim(row)
 enddo

 ! sink number of each star, i.e. its position in the hierarchy string
 row = ''
 do i=1,nleaves
    isink = findloc_local(labels%sink(1:labels%sink_num),leaves(i))
    write(num,"(i0)") isink
    call put_label(row,trim(num),xleaf(i))
 enddo
 print "(a)",row(1:xleaf(nleaves)+len_trim(num))//'     <- sink id'
 print "(a)",''

end subroutine print_hierarchy

!--------------------------------------------------------------------------
!+
!  column of a node in the tree: the column of a star, or the middle
!  of the columns of the two components of a binary
!+
!--------------------------------------------------------------------------
recursive integer function node_column(label,leaves,xleaf) result(x)
 character(len=*), intent(in) :: label
 character(len=*), intent(in) :: leaves(:)
 integer,          intent(in) :: xleaf(:)
 integer :: i

 do i=1,size(leaves)
    if (trim(leaves(i)) == label) then
       x = xleaf(i)
       return
    endif
 enddo
 x = (node_column(label//'1',leaves,xleaf) + node_column(label//'2',leaves,xleaf))/2

end function node_column

!--------------------------------------------------------------------------
!+
!  write a label into a row of text, centred on column x
!+
!--------------------------------------------------------------------------
subroutine put_label(row,label,x)
 character(len=*), intent(inout) :: row
 character(len=*), intent(in)    :: label
 integer,          intent(in)    :: x
 integer :: i1,i2

 if (len(label) == 0) return
 i1 = max(1,x - (len(label)-1)/2)
 i2 = min(len(row),i1 + len(label) - 1)
 row(i1:i2) = label(1:i2-i1+1)

end subroutine put_label

!--------------------------------------------------------------------------
!+
!  check that a substitution is valid
!+
!--------------------------------------------------------------------------
subroutine check_substitution(hier_prefix, semimajoraxis, filename, ierr)
 character(len=20), intent(in)  :: hier_prefix
 character(len=*),  intent(in)  :: filename
 real,              intent(in)  :: semimajoraxis
 integer,           intent(out) :: ierr
 integer :: lines, subst_index
 real :: ma, mb, a_comp, q_comp, e_comp, period_ratio, criterion
 real :: data(hier_db_size,hier_db_prop)

 ! Check that star to be substituted exists in HIERARCHY file
 call find_data_index(hier_prefix, subst_index, filename, ierr)

 if (subst_index == 0) then
    print "(1x,a)",'ERROR: set_multiple: star '//trim(hier_prefix)//' not present in HIERARCHY file.'
    ierr = ierr_missstar
    return
 endif

 ! Check that star to be substituted has not already been substituted
 call load_hierarchy_file(filename, data, lines, ierr)

 if (data(subst_index,1) == 0) then
    print "(1x,a)",'ERROR: set_multiple: star '//trim(hier_prefix)//' substituted yet.'
    ierr = ierr_subststar
    return
 endif

 ! Check that substituting binary Omega = 0
 if (data(subst_index,10) /= 0) then
    print "(1x,a)",'ERROR: set_multiple: at the moment phantom can subst only Omega=0 binaries.'
    ierr = ierr_Omegasubst
    return
 endif

 ! Mardling&Aarseth (2001) criterion check
 ma = data(subst_index, 3)
 mb = data(subst_index, 4)
 a_comp = data(subst_index, 5)
 e_comp = data(subst_index, 6)

 q_comp = ma/mb
 if (q_comp>1) q_comp=q_comp**(-1)

 ! Mardling&Aarseth (2001) criterion check

 period_ratio = sqrt((a_comp*a_comp*a_comp)/(ma+mb)/(semimajoraxis*semimajoraxis*semimajoraxis)*(ma)) ! Po/Pi
 criterion = 4.7*(1-e_comp)**(-1.8)*(1+e_comp)**(0.6)*(1+q_comp)**(0.1)

 if (criterion > period_ratio) then
    print "(1x,a)",'WARNING: set_multiple: orbital parameters does not satisfy Mardling and Aarseth stability criterion.'
 endif

end subroutine check_substitution

!--------------------------------------------------------------------------
!+
!  find orbital elements for a particular binary
!+
!--------------------------------------------------------------------------
subroutine find_hier_level_orb_elem(hl_temp, hs, m1, m2, accr1, accr2, &
                                    binary_a, binary_e, binary_i, binary_O, &
                                    binary_w, binary_f)
 character(len=20),         intent(in)  :: hl_temp
 type(hierarchical_system), intent(in)  :: hs
 real,                      intent(out) :: m1, m2, accr1, accr2, binary_a, binary_e, binary_i, binary_O, binary_w, binary_f

 integer :: hl_index, idx

 !print*, 'find ', trim(hl_temp)

 m1 = get_hierarchical_level_mass(trim(adjustl(hl_temp))//'1', hs)
 m2 = get_hierarchical_level_mass(trim(adjustl(hl_temp))//'2', hs)

!print *,'labels passing: ', trim(adjustl(hl_temp))//'1 ', m1,trim(adjustl(hl_temp))//'2 ',m2

 idx = findloc_local(hs%labels%sink, trim(adjustl(hl_temp))//'1')
 if (idx > 0) then
    accr1 = hs%sinks(idx)%accr
 else
    accr1 = 1.
 endif

 idx = findloc_local(hs%labels%sink, trim(adjustl(hl_temp))//'2')
 if (idx > 0) then
    accr2 = hs%sinks(idx)%accr
 else
    accr2 = 1.
 endif

 hl_index = findloc_local(hs%labels%hl, trim(adjustl(hl_temp)))

 if (hl_index > 0) then
    binary_a = hs%levels(hl_index)%a
    binary_e = hs%levels(hl_index)%e
    binary_O = hs%levels(hl_index)%O
    binary_w = hs%levels(hl_index)%w
    binary_i = hs%levels(hl_index)%inc
    binary_f = hs%levels(hl_index)%f
 else
    print*,'ERROR: could not find level ',trim(adjustl(hl_temp)),' in hierarchy'
    print*,'Available levels: ',hs%labels%hl
    binary_a = 100.
    binary_e = 0.
    binary_O = 0.
    binary_w = 0.
    binary_i = 0.
    binary_f = 0.
 endif

end subroutine find_hier_level_orb_elem

!--------------------------------------------------------------------------
!+
!  find index of ptmass corresponding to a particular label
!+
!--------------------------------------------------------------------------
subroutine find_ptmass_index(hier_label, index, filename, ierr)
 integer,           intent(out) :: index, ierr
 character(len=20), intent(in)  :: hier_label
 character(len=*),  intent(in)  :: filename

 real :: data(hier_db_size,hier_db_prop)
 integer :: lines, hier_int, io, loc

 call load_hierarchy_file(filename, data, lines, ierr)

 read(hier_label,*,iostat=io) hier_int

 loc = findloc_local(int(data(1:lines,2)), hier_int)
 if (loc > 0) then
    index = int(data(loc,1))
 else
    index = 0
 endif

end subroutine find_ptmass_index

!--------------------------------------------------------------------------
!+
!  retrieve sink ptmass index or hierarchical level's sink ptmass
!  indexes in HIERARCHY file
!+
!--------------------------------------------------------------------------
subroutine find_hierarchy_index(level, int_sinks, inner_sinks_num, filename)
 character(len=10), intent(in)  :: level
 integer,           intent(out) :: inner_sinks_num
 integer,           intent(out) :: int_sinks(max_hier_levels)
 character(len=*),  intent(in)  :: filename
 real :: data(hier_db_size,hier_db_prop)
 integer                :: i, lines, ierr, h_index, io
 character(len=10)      :: label = '         '

 read(level,*,iostat=io) h_index

 call load_hierarchy_file(filename, data, lines, ierr)

 inner_sinks_num = 0
 do i=1, lines
    write(label, '(i0)') int(data(i,2))
    if (data(i,1) > 0 .and. (len(trim(label)) >= len(trim(level))) .and. (label(:len(trim(level))) == trim(level))) then
       inner_sinks_num = inner_sinks_num+1
       int_sinks(inner_sinks_num) = int(data(i,1))
    endif
 enddo

end subroutine find_hierarchy_index

!--------------------------------------------------------------------------
!+
!  return position of hier_label in hierarchy file
!+
!--------------------------------------------------------------------------
subroutine find_data_index(hier_label, index, filename, ierr)
 integer,           intent(out) :: index, ierr
 character(len=20), intent(in)  :: hier_label
 character(len=*),  intent(in)  :: filename
 real :: data(hier_db_size,hier_db_prop)
 integer :: lines, hier_int, io

 call load_hierarchy_file(filename, data, lines, ierr)

 read(hier_label,*,iostat=io) hier_int
 index = findloc_local(int(data(1:lines,2)), hier_int)

end subroutine find_data_index

!------------------------------------------------------
!+
!  Rotate an (x,y,z) point by theta radians around an
!  axis with alpha, beta and gamma eulerian angles
!+
!------------------------------------------------------
pure subroutine gen_rotate(xyz,alpha,beta,gamma,theta)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: alpha, beta, gamma,theta
 real :: xi,yi,zi,A,B,C,D,E,F,G,H,I,nx,ny,nz

 nx=cos(alpha)
 ny=cos(beta)
 nz=cos(gamma)

 A=cos(theta)+nx**2*(1-cos(theta))
 B=nx*ny*(1-cos(theta))-nz*sin(theta)
 C=nx*nz*(1-cos(theta))+ny*sin(theta)
 D=nx*ny*(1-cos(theta))+nz*sin(theta)
 E=cos(theta)+ny**2*(1-cos(theta))
 F=ny*nz*(1-cos(theta))-nx*sin(theta)
 G=nx*nz*(1-cos(theta))-ny*sin(theta)
 H=ny*nz*(1-cos(theta))+nx*sin(theta)
 I=cos(theta)+nz**2*(1-cos(theta))

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) = A*xi+B*yi+C*zi
 xyz(2) = D*xi+E*yi+F*zi
 xyz(3) = G*xi+H*yi+I*zi

end subroutine gen_rotate

!------------------------------------------------------
!+
!  findloc replacement for non-working findloc
!  intrinsic implemented in aocc compiler
!+
!------------------------------------------------------
function findloc_int(array, value) result(index)
 integer, intent(in) :: array(:)
 integer, intent(in) :: value
 integer :: index, i

 index = 0
 do i=1,size(array)
    if (array(i) == value) then
       index = i
       return
    endif
 enddo

end function findloc_int

function findloc_string(array, value) result(index)
 character(len=*), intent(in) :: array(:)
 character(len=*), intent(in) :: value
 integer :: index, i

 index = 0
 do i=1,size(array)
    if (trim(array(i)) == trim(value)) then
       index = i
       return
    endif
 enddo

end function findloc_string

function findloc_logical(array, value) result(index)
 logical, intent(in) :: array(:)
 logical, intent(in) :: value
 integer :: index, i

 index = 0
 do i=1,size(array)
    if (array(i) .eqv. value) then
       index = i
       return
    endif
 enddo

end function findloc_logical

end module sethier_utils
