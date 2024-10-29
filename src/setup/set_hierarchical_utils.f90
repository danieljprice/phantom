!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module sethier_utils
!
! sethierarchical
!
! :References: None
!
! :Owner: Simone Ceppi
!
! :Runtime parameters: None
!
! :Dependencies: None
!

 implicit none
 integer, parameter, public :: max_hier_levels=10
 integer, parameter, public :: hier_db_size=10, hier_db_prop=10
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

 public :: hierarchical_system,process_hierarchy
 public :: recursive_splitting, get_hierarchical_level_mass
 public :: load_hierarchy_file,update_hierarchy_file,check_substitution
 public :: find_hierarchy_index,find_hier_level_orb_elem,find_data_index
 public :: gen_rotate
 public :: find_ptmass_index

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

 integer :: i,j, del_pos, pre_del_pos
 type(hs_labels)    :: mylabels
 character(len=:), allocatable :: sink, temp_hl

 mylabels%sink = ''
 mylabels%hl = ''

 pre_del_pos = 0
 mylabels%sink_num = 0

 ! count the number of sinks in the system
 del_pos = scan(string, ',')
 do while (del_pos > pre_del_pos)
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

 !print *, 'sink to generate: ', sink_list_temp

 longest = 0
 if (sink_num>1) then
    ! Find the longest
    do i=1,sink_num
       if (longest < len(trim(sink_list_temp(i)))) then
          longest = len(trim(sink_list_temp(i)))
       endif
    enddo

    !print*, 'longest lenght: ', longest

    ! Select the longests and cut them
    longests_len = 0
    do i=1,sink_num
       if (len(trim(sink_list_temp(i))) == longest) then
          longests_len = longests_len+1
          longests(longests_len) = trim(sink_list_temp(i))

          sink_list_temp(i) = sink_list_temp(i)(:len(trim(longests(longests_len)))-1)
       endif
    enddo

    !print*, 'found ', longests_len, ' sinks long ', longest
    !print*, '  they are ', longests

    ! Cut the longest and add to split list with no doubles
    count=0
    do i=1,longests_len
       longest_cut = longests(i)
       if (.not. any(new_splits == longest_cut(:longest-1))) then
          count = count + 1
          new_splits(count) = longests(i)(1:longest-1)
       endif
    enddo

    !print *, 'new splits to generate longests are ', count
    !print *, '    i.e. ', new_splits

    !print *, 'up to now splits are ', split_list

    ! Add new splits to split_list
    do i=splits+1, min(splits+count,size(split_list))
       split_list(i) = new_splits(i-splits)
    enddo
    splits = splits + count

    !print *, 'adding new splits: ', split_list

    ! Clean sink_list_temp from doubles
    count=0
    do i=1,sink_num
       if (.not. any(new_sink_list == sink_list_temp(i))) then
          count = count + 1
          new_sink_list(count) = sink_list_temp(i)
       endif
    enddo

    !print *, 'sink to generate after new splits: ', sink_list_temp

    !print *, '_'
    !print *, new_sink_list(:count)
    !print *, split_list(:splits)
    !print *, splits

    call recursive_splitting(count, new_sink_list(:), split_list(:), splits)

 endif

end subroutine recursive_splitting

!--------------------------------------------------------------------------
!+
!  compute the total mass of an hierarchical level
!+
!--------------------------------------------------------------------------
!pure
real function get_hierarchical_level_mass(level, hs) result(part_mass)
 character(len=*), intent(in) :: level
 type(hierarchical_system), intent(in) :: hs

 integer :: i

 !print*,'get_hierarchical_lvl_mass ', trim(level)

 part_mass = 0
 !print*, 'computing!'
 do i=1, hs%labels%sink_num
    !print*,'-',trim(hs%labels%sink(i)),'-', '   ', '-',trim(level),'- ', len(trim(level))
    if ((len(trim(hs%labels%sink(i))) >= len(trim(level))) .and. &
         ((hs%labels%sink(i)(:len(trim(level))) == trim(level)))) then
       !print*,'inside ', trim(hs%labels%sink(i)), hs%sinks(i)%mass
       part_mass = part_mass + hs%sinks(i)%mass
    endif
 enddo

end function get_hierarchical_level_mass

!--------------------------------------------------------------------------
!+
!  read information from the hierarchy file
!+
!--------------------------------------------------------------------------
subroutine load_hierarchy_file(prefix, data, lines, ierr)
 character(len=20), intent(in), optional :: prefix
 integer,           intent(out) :: ierr, lines
 real, dimension(hier_db_size,hier_db_prop), intent(out) :: data

 character(len=20) :: filename
 integer :: io
 logical :: iexist

 if (present(prefix)) then
    filename = trim(prefix)//'.hierarchy'
 else
    filename = 'HIERARCHY'
 endif

 inquire(file=trim(filename), exist=iexist)

 if (iexist) then
    open(2,file=trim(filename),status='old')
    lines=0
    do
       read(2, *,iostat=io) data(lines+1,:)
       if (io/=0) exit
       lines = lines + 1
    enddo
    close(2)
 else
    print "(1x,a)",'ERROR: set_multiple: there is no HIERARCHY file, cannot perform subtitution.'
    ierr = ierr_HIER2
 endif
end subroutine load_hierarchy_file

!--------------------------------------------------------------------------
!+
!  write information to the hierarchy file
!+
!--------------------------------------------------------------------------
subroutine update_hierarchy_file(prefix, hs, data, lines, hier_prefix, i1, i2, ierr)
 integer, intent(in) :: i1, i2
 character(len=20), intent(in), optional:: prefix
 character(len=20), intent(in) :: hier_prefix
 type(hierarchical_system), intent(in) :: hs
 real, dimension(hier_db_size,hier_db_prop), intent(inout) :: data
 integer,    intent(inout)    :: lines
 integer,    intent(out)    :: ierr

 integer :: i,iu
 real :: mprimary, msecondary, semimajoraxis, eccentricity, incl, arg_peri, posang_ascnode, binary_f
 real :: accr1, accr2
 real :: period
 character(len=20) :: filename
 logical :: iexist

 !print*,'updating: ',trim(adjustl(hier_prefix))

 call find_hier_level_orb_elem(hier_prefix, hs, mprimary, msecondary, accr1, accr2, &
       semimajoraxis, eccentricity, incl, posang_ascnode, &
       arg_peri, binary_f)

 !print*,'elements of ', hier_prefix, ' : ', mprimary, msecondary, accr1, accr2, &
 !     semimajoraxis, eccentricity, incl, posang_ascnode, &
 !     arg_peri, binary_f

 period = sqrt(4.*pi**2*semimajoraxis**3/(mprimary+msecondary))

 if (present(prefix)) then
    filename = trim(prefix)//'.hierarchy'
 else
    filename = 'HIERARCHY'
 endif

 if (lines > 0) then
    open(newunit=iu,file=trim(filename),status='old')
    do i=1,lines
       write(iu,*) int(data(i,1)), int(data(i,2)), data(i,3:)
    enddo
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

 call load_hierarchy_file(prefix, data, lines, ierr)

end subroutine update_hierarchy_file

!--------------------------------------------------------------------------
!+
!  check that a substitution is valid
!+
!--------------------------------------------------------------------------
subroutine check_substitution(hier_prefix, semimajoraxis, prefix, ierr)
 character(len=20), intent(in) :: hier_prefix, prefix
 real, intent(in) :: semimajoraxis
 integer, intent(out) :: ierr

 integer :: lines, subst_index
 real :: ma, mb, a_comp, q_comp, e_comp, period_ratio, criterion
 real, dimension(hier_db_size,hier_db_prop) :: data

 ! Check that star to be substituted exists in HIERARCHY file
 call find_data_index(hier_prefix, subst_index, prefix, ierr)

 if (subst_index == 0) then
    print "(1x,a)",'ERROR: set_multiple: star '//trim(hier_prefix)//' not present in HIERARCHY file.'
    ierr = ierr_missstar
    return
 endif

 ! Check that star to be substituted has not already been substituted
 call load_hierarchy_file(prefix, data, lines, ierr)

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

 character(len=20),         intent(in) :: hl_temp
 type(hierarchical_system), intent(in) :: hs
 real, intent(out) :: m1, m2, accr1, accr2, binary_a, binary_e, binary_i, binary_O, binary_w, binary_f

 integer :: hl_index

 !print*, 'find ', trim(hl_temp)

 m1 = get_hierarchical_level_mass(trim(adjustl(hl_temp))//'1', hs)
 m2 = get_hierarchical_level_mass(trim(adjustl(hl_temp))//'2', hs)

 !print *,'labels passing: ', trim(adjustl(hl_temp))//'1 ', m1,trim(adjustl(hl_temp))//'2 ',m2

 if (any(hs%labels%sink == trim(adjustl(hl_temp))//'1')) then
    accr1 = hs%sinks(findloc(hs%labels%sink,trim(adjustl(hl_temp))//'1', 1))%accr
 else
    accr1 = 1.
 endif

 if (any(hs%labels%sink == trim(adjustl(hl_temp))//'2')) then
    accr2 = hs%sinks(findloc(hs%labels%sink,trim(adjustl(hl_temp))//'2', 1))%accr
 else
    accr2 = 1.
 endif

 hl_index = findloc(hs%labels%hl, trim(adjustl(hl_temp)), 1)

 binary_a = hs%levels(hl_index)%a
 binary_e = hs%levels(hl_index)%e
 binary_O = hs%levels(hl_index)%O
 binary_w = hs%levels(hl_index)%w
 binary_i = hs%levels(hl_index)%inc
 binary_f = hs%levels(hl_index)%f

end subroutine find_hier_level_orb_elem

!--------------------------------------------------------------------------
!+
!  find index of ptmass corresponding to a particular label
!+
!--------------------------------------------------------------------------
subroutine find_ptmass_index(hier_label, index, prefix, ierr)
 integer,    intent(out)    :: index, ierr
 character(len=20), intent(in), optional:: prefix, hier_label

 real, dimension(hier_db_size,hier_db_prop) :: data
 integer :: lines, hier_int, io

 call load_hierarchy_file(prefix, data, lines, ierr)

 read(hier_label,*,iostat=io) hier_int

 index = int(data(findloc(int(data(:,2)), hier_int, 1),1))

end subroutine find_ptmass_index

!--------------------------------------------------------------------------
!+
!  retrieve sink ptmass index or hierarchical level's sink ptmass
!  indexes in HIERARCHY file
!+
!--------------------------------------------------------------------------
subroutine find_hierarchy_index(level, int_sinks, inner_sinks_num, prefix)
 character(len=10), intent(in) :: level
 integer, intent(out) :: inner_sinks_num
 integer, intent(out)                :: int_sinks(max_hier_levels)
 character(len=20), optional, intent(in) :: prefix

 real, dimension(hier_db_size,hier_db_prop) :: data
 integer                :: i, lines, ierr, h_index, io

 character(len=10)      :: label = '         '

 read(level, *,iostat=io) h_index

 call load_hierarchy_file(prefix, data, lines, ierr)

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
subroutine find_data_index(hier_label, index, prefix, ierr)
 integer,    intent(out)    :: index, ierr
 character(len=20), intent(in), optional:: prefix, hier_label

 real, dimension(hier_db_size,hier_db_prop) :: data
 integer :: lines, hier_int, io

 call load_hierarchy_file(prefix, data, lines, ierr)

 read(hier_label,*,iostat=io) hier_int

 index = findloc(int(data(:,2)), hier_int, 1)

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

end module sethier_utils
