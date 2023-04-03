!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module sethierarchical
!
! sethierarchical
!
! :References: None
!
! :Owner: Simone Ceppi
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, prompting, setbinary
!

 implicit none

 public :: set_hierarchical, set_multiple
 public :: set_hierarchical_interactively
 public :: write_hierarchical_setupfile
 public :: read_hierarchical_setupfile
 public :: set_hierarchical_default_options

 public :: get_hierarchical_level_com

 public :: get_hier_level_mass

 integer, parameter :: max_hier_levels=10
 integer, parameter :: hier_db_size=24, hier_db_prop=10

 

 character(len=100) :: hierarchy = '111,112,121,1221,1222'
 integer :: sink_num, hl_num
 character(len=10) :: sink_labels(max_hier_levels), hl_labels(max_hier_levels)
 real :: mass(max_hier_levels), accr(max_hier_levels)
 real :: a(max_hier_levels), e(max_hier_levels), inc(max_hier_levels), O(max_hier_levels), w(max_hier_levels), f(max_hier_levels)

 public :: hierarchy
 public :: sink_num, hl_num
 public :: sink_labels, hl_labels
 public :: mass, accr
 public :: a, e, inc, O, w, f

 private

  real, parameter :: pi = 4.*atan(1.)
 integer, parameter :: &
   ierr_HIER2 = 6, &
   ierr_subststar = 7, &
   ierr_Omegasubst = 8, &
   ierr_missstar = 9

 

contains

subroutine set_hierarchical(prefix, nptmass, xyzmh_ptmass, vxyz_ptmass, ierr)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 integer, intent(out)   :: ierr
 character(len=20), intent(in) :: prefix

 integer :: i
 real :: m1, m2, accr1, accr2, binary_a, binary_e, binary_i, binary_O, binary_w, binary_f

 integer :: splits, sink_num_temp, hl_index, subst
 character(len=10) :: sink_list(max_hier_levels), split_list(max_hier_levels)
 character(len=20) :: hl_temp

 splits = 0
 sink_list = sink_labels
 sink_num_temp = sink_num

 call recursive_splitting(sink_num_temp, sink_list, split_list, splits)

 do i=splits,1,-1

    hl_temp = trim(split_list(i))

    call find_hier_level_orb_elem(hl_temp, m1, m2, accr1, accr2, &
                                  binary_a, binary_e, binary_i, binary_O, &
                                  binary_w, binary_f)

    read(hl_temp,*,iostat=subst) subst
    call set_multiple(m1,m2,semimajoraxis=binary_a,eccentricity=binary_e, &
            posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
            f=binary_f,accretion_radius1=accr1,accretion_radius2=accr2, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,&
            ierr=ierr, subst=subst, prefix=prefix)

 enddo
end subroutine set_hierarchical



subroutine set_hierarchical_interactively()
 use prompting, only:prompt

 call prompt('What is the hierarchy?',hierarchy)

 call process_hierarchy()

end subroutine set_hierarchical_interactively



subroutine write_hierarchical_setupfile(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 integer :: i

 write(iunit,"(/,a)") '# options for hierarchical system'

 call write_inopt(hierarchy, 'hierarchy','', iunit)

 call process_hierarchy()

 write(iunit,"(/,a)") '### sink properties'
 do i=1,sink_num
    call write_inopt(mass(i), trim(sink_labels(i))//'_mass','', iunit)
    call write_inopt(accr(i), trim(sink_labels(i))//'_accr','', iunit)
 enddo

 write(iunit,"(/,a)") '### orbit properties'

 do i=1,hl_num
    call write_inopt(a(i), trim(hl_labels(i))//'_a','',iunit)
    call write_inopt(e(i), trim(hl_labels(i))//'_e','',iunit)
    call write_inopt(inc(i), trim(hl_labels(i))//'_i','',iunit)
    call write_inopt(O(i), trim(hl_labels(i))//'_O','',iunit)
    call write_inopt(w(i), trim(hl_labels(i))//'_w','',iunit)
    call write_inopt(f(i), trim(hl_labels(i))//'_f','',iunit)
 enddo


end subroutine write_hierarchical_setupfile


subroutine read_hierarchical_setupfile(db, nerr)

 use infile_utils, only:read_inopt, inopts
 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(inout) :: nerr

 integer :: i

 call read_inopt(hierarchy,'hierarchy',db,errcount=nerr)

 call process_hierarchy()

 do i=1,sink_num
    call read_inopt(mass(i), trim(sink_labels(i))//'_mass',db,errcount=nerr)
    call read_inopt(accr(i), trim(sink_labels(i))//'_accr',db,errcount=nerr)
 enddo

 do i=1,hl_num
    call read_inopt(a(i), trim(hl_labels(i))//'_a',db,errcount=nerr)
    call read_inopt(e(i), trim(hl_labels(i))//'_e',db,errcount=nerr)
    call read_inopt(inc(i), trim(hl_labels(i))//'_i',db,errcount=nerr)
    call read_inopt(O(i), trim(hl_labels(i))//'_O',db,errcount=nerr)
    call read_inopt(w(i), trim(hl_labels(i))//'_w',db,errcount=nerr)
    call read_inopt(f(i), trim(hl_labels(i))//'_f',db,errcount=nerr)
 enddo
end subroutine read_hierarchical_setupfile


subroutine set_hierarchical_default_options()
 integer :: i

 call process_hierarchy()

 do i=1,sink_num
    mass(i)=1.
    accr(i)=1.
 enddo

 do i=1,hl_num
    a(i)=10000./10**len(trim(hl_labels(i)))
    e(i)=0.
    inc(i)=0.
    O(i)=0.
    w(i)=0.
    f(i)=0.
 enddo

end subroutine set_hierarchical_default_options


 !--------------------------------------------------------------------------
 !
 ! Process the .setup hierarchy string to extract information for building the system
 !
 !--------------------------------------------------------------------------
subroutine process_hierarchy()
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
 enddo

 sink_labels(sink_num+1) = hierarchy(pre_del_pos+1:)
 sink_num = sink_num+1

 hl_num = 0
 do i=1,sink_num
    sink = trim(sink_labels(i))
    del_pos = len(trim(sink))
    do j=1,del_pos
       temp_hl = trim(sink(1:j))
       if ( .not. any(hl_labels == trim(temp_hl))) then
          if ( any(sink_labels == trim(temp_hl)) .or. any(sink_labels == trim(temp_hl))) then
             ! add to discs
          else
             hl_num=hl_num+1
             hl_labels(hl_num) = trim(temp_hl)
          endif
       endif
    enddo
 enddo

end subroutine process_hierarchy


 !--------------------------------------------------------------------------
 !
 ! Reverse the splitting process from sink_labels to build the system
 !
 !--------------------------------------------------------------------------
recursive subroutine recursive_splitting(sink_num, sink_list, split_list, splits)
 character(len=10), intent(in) :: sink_list(:)
 character(len=10), intent(inout) :: split_list(:)
 integer, intent(inout)    :: splits
 integer, intent(in) ::sink_num

 integer :: i, longest, longests_len, count
 character(len=10) :: longests(max_hier_levels), new_splits(max_hier_levels)
 character(len=10) :: new_sink_list(max_hier_levels), longest_cut, sink_list_temp(max_hier_levels)

 sink_list_temp = sink_list

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
    do i=splits+1, splits+count
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

    call recursive_splitting(count, new_sink_list(:count), split_list(:splits), splits)

 endif
end subroutine recursive_splitting

 !--------------------------------------------------------------------------
 !
 ! Compute the total mass of an hierarchical level
 !
 !--------------------------------------------------------------------------
real function get_hier_level_mass(level)
 character(*), intent(in) :: level

 real :: part_mass
 integer :: i

 part_mass = 0
 do i=1, sink_num
    if ((len(trim(sink_labels(i))) >= len(trim(level))) .and. (sink_labels(i)(:len(trim(level))) == trim(level))) then
       part_mass = part_mass + mass(i)
    endif
 enddo

 get_hier_level_mass = part_mass

end function get_hier_level_mass


 !--------------------------------------------------------------------------
 !
 ! Compute position and velocity of hierarchical level
 !
 !--------------------------------------------------------------------------
subroutine get_hierarchical_level_com(level, xorigin, vorigin, xyzmh_ptmass, vxyz_ptmass, prefix)
 character(*), intent(in) :: level
 real, intent(out) :: xorigin(:), vorigin(:)
 real, intent(in) :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 character(len=20), intent(in) :: prefix

 integer :: int_sinks(max_hier_levels), inner_sinks_num, i
 real :: mass

 call find_hierarchy_index(level, int_sinks, inner_sinks_num, prefix)

 xorigin = 0.
 vorigin = 0.
 mass = 0.

 do i=1,inner_sinks_num
    xorigin(:) = xorigin(:)+xyzmh_ptmass(4,int_sinks(i))*xyzmh_ptmass(1:3,int_sinks(i))
    vorigin(:) = vorigin(:)+xyzmh_ptmass(4,int_sinks(i))*vxyz_ptmass(1:3,int_sinks(i))
    mass = mass + xyzmh_ptmass(4,int_sinks(i))
 enddo
 xorigin = xorigin/mass
 vorigin = vorigin/mass


end subroutine get_hierarchical_level_com

 !--------------------------------------------------------------------------
 !
 ! Retrieve sink index or hierarchical level's sink indexes in HIERARCHY file
 !
 !--------------------------------------------------------------------------
subroutine find_hierarchy_index(level, int_sinks, inner_sinks_num, prefix)
 character(*), intent(in) :: level
 integer, intent(out) :: inner_sinks_num
 integer, intent(out)                :: int_sinks(max_hier_levels)
 character(len=20), optional, intent(in) :: prefix

  real, dimension(hier_db_size,hier_db_prop) :: data
 integer                :: i, lines, ierr, h_index

 character(len=10)      :: label = '         '

 read(level, *, iostat=h_index) h_index

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

subroutine find_hier_level_orb_elem(hl_temp, m1, m2, accr1, accr2, &
                                    binary_a, binary_e, binary_i, binary_O, &
                                    binary_w, binary_f)

  character(len=20), intent(in) :: hl_temp
  real, intent(out) :: m1, m2, accr1, accr2, binary_a, binary_e, binary_i, binary_O, binary_w, binary_f

  integer :: hl_index
  
  m1 = get_hier_level_mass(trim(hl_temp)//'1')
  m2 = get_hier_level_mass(trim(hl_temp)//'2')
  
  if (any(sink_labels == trim(hl_temp)//'1')) then
     accr1 = accr(findloc(sink_labels,trim(hl_temp)//'1', 1))
  else
     accr1 = 1.
  endif
  
  if (any(sink_labels == trim(hl_temp)//'2')) then
     accr2 = accr(findloc(sink_labels,trim(hl_temp)//'2', 1))
  else
     accr2 = 1.
  endif
  
  hl_index = findloc(hl_labels, trim(hl_temp), 1)
  
  binary_a =  a(hl_index)
  binary_e = e(hl_index)
  binary_O = O(hl_index)
  binary_w = w(hl_index)
  binary_i = inc(hl_index)
  binary_f = f(hl_index)

end subroutine find_hier_level_orb_elem


!subroutine read_hierarchy_file(hier_label, index, prefix, ierr)
!  integer,    intent(out)    :: index, ierr
!  character(len=20), intent(in), optional:: prefix
!  
!  real, dimension(hier_db_size,hier_db_prop) :: data
!  integer :: lines, hier_int
!  
!  call load_hierarchy_file(prefix, data, lines, ierr)

!  read(hier_label,*,iostat=hier_int) hier_int
  
!  index = findloc(int(data(:,1)), hier_int, 1)

!  print*, index
  
  
!end subroutine read_hierarchy_file


subroutine load_hierarchy_file(prefix, data, lines, ierr)
  integer,    intent(out)    :: ierr, lines
  character(len=20), intent(in), optional:: prefix
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
     open(2, file = trim(filename), status = 'old')
     lines=0
     do
        read(2, *, iostat=io) data(lines+1,:)
        if (io/=0) exit
        lines = lines + 1
     enddo
     close(2)
  else
     print "(1x,a)",'ERROR: set_multiple: there is no HIERARCHY file, cannot perform subtitution.'
     ierr = ierr_HIER2
  endif
end subroutine load_hierarchy_file

subroutine update_hierarchy_file(prefix, data, lines, hier_prefix, i1, i2, ierr)
  integer,    intent(in)    :: lines, i1, i2
  character(len=20), intent(in), optional:: prefix
  real, dimension(hier_db_size,hier_db_prop), intent(in) :: data
  character(len=20), intent(in) :: hier_prefix
  integer,    intent(out)    :: ierr

  integer :: i
  real :: mprimary, msecondary, semimajoraxis, eccentricity, incl, arg_peri, posang_ascnode, binary_f
  real :: accr1, accr2
  real :: period
  character(len=20) :: filename
  logical :: iexist

  call find_hier_level_orb_elem(adjustl(hier_prefix), mprimary, msecondary, accr1, accr2, &
       semimajoraxis, eccentricity, incl, posang_ascnode, &
       arg_peri, binary_f)

  period = sqrt(4.*pi**2*semimajoraxis**3/(mprimary+msecondary))
  
  if (present(prefix)) then
     filename = trim(prefix)//'.hierarchy'
  else
     filename = 'HIERARCHY'
  endif
  
  if (lines>0) then  
     open(2, file = trim(filename), status = 'old')
     do i=1,lines
        write(2,*) int(data(i,1)), int(data(i,2)), data(i,3:)
     enddo
  else
     inquire(file=trim(filename), exist=iexist)
     if (iexist) then
        print "(1x,a)",'WARNING: set_multiple: deleting an existing HIERARCHY file.'
        open(1, file=trim(filename), status='old')
        close(1, status='delete')
     endif
     
     open(2, file = trim(filename), status = 'new')
  end if
  write(2,*) i1, trim(hier_prefix)//"1", mprimary, msecondary, semimajoraxis, eccentricity, &
       period, incl, arg_peri, posang_ascnode
  write(2,*) i2, trim(hier_prefix)//"2", msecondary, mprimary, semimajoraxis, eccentricity, &
       period, incl, arg_peri, posang_ascnode
  close(2)
  
end subroutine update_hierarchy_file

subroutine check_substitution()

end subroutine check_substitution



!----------------------------------------------------------------
!+
!  setup for a multiple, using set_binary
!+
!----------------------------------------------------------------
subroutine set_multiple(m1,m2,semimajoraxis,eccentricity, &
                      accretion_radius1,accretion_radius2, &
                      xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                      posang_ascnode,arg_peri,incl,f,verbose,subst, prefix)
  use setbinary, only: set_binary
 real,    intent(in)    :: m1,m2
 real,    intent(in)    :: semimajoraxis,eccentricity
 real,    intent(in)    :: accretion_radius1,accretion_radius2
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 integer, intent(out)   :: ierr
 real,    intent(in),  optional :: posang_ascnode,arg_peri,incl,f
 integer, intent(in),  optional :: subst
 real,    intent(out), optional :: omega_corotate
 logical, intent(in),  optional :: verbose
 character(len=20), optional, intent(in) :: prefix

 integer :: i1,i2,i,subst_index
 real    :: mtot,period
 real    :: x_subst(3),v_subst(3)
 real    :: omega,inc
 !logical :: do_verbose

 real, dimension(hier_db_size,hier_db_prop) :: data
 character(len=20)      :: filename
 character(len=20)      :: hier_prefix
 logical                :: iexist
 integer                :: io, lines
 real                   :: period_ratio,criterion,q_comp,a_comp,e_comp,m_comp
 real                   :: rel_posang_ascnode=0.,rel_arg_peri=0.,rel_incl=0.
 real                   :: q2,mprimary,msecondary
 real                   :: alpha_y, beta_y, gamma_y, alpha_z, beta_z, gamma_z, sign_alpha, sign_gamma

 ierr = 0
 !do_verbose = .true.
 !if (present(verbose)) do_verbose = verbose

  !--- Load/Create HIERARCHY file: xyzmh_ptmass index | hierarchical index | star mass | companion star mass | semi-major axis | eccentricity | period | inclination | argument of pericenter | ascending node longitude
 if (present(subst) .and. subst>10) then
    call load_hierarchy_file(prefix, data, lines, ierr)
 else
    mtot = m1 + m2
    period = sqrt(4.*pi**2*semimajoraxis**3/mtot)

    hier_prefix = '1'

    call update_hierarchy_file(prefix, data, 0, hier_prefix, 1, 2, ierr) 
 endif
 subst_index = 0
 !--- Checks to avoid bad substitutions
 if (present(subst) .and. subst>10) then
    write(hier_prefix, *) subst

    !call read_hierarchy_file(hier_prefix, index, prefix, ierr)!check_substitution()
    
    io=0
    mtot = 0.
    do i=1,lines
       if (data(i,2)==abs(subst)) then ! Check that star to be substituted exists in HIERARCHY file
          if (data(i,1)==0) then ! Check that star to be substituted has not already been substituted
             print "(1x,a)",'ERROR: set_multiple: star '//trim(hier_prefix)//' substituted yet.'
             ierr = ierr_subststar
          endif
          subst_index = int(data(i,1))
          data(i,1) = 0

          if (subst>0) then
             rel_posang_ascnode = data(i, 10)

             if (rel_posang_ascnode /= 0) then
                print "(1x,a)",'ERROR: set_multiple: at the moment phantom can subst only Omega=0 binaries.'
                ierr = ierr_Omegasubst
             endif

             rel_arg_peri= data(i, 9)
             rel_incl = data(i, 8)
          else
             rel_posang_ascnode = posang_ascnode
             rel_arg_peri = arg_peri
             rel_incl = incl
          endif

          mtot = data(i, 3)
          m_comp = data(i, 4)
          a_comp = data(i, 5)
          e_comp = data(i, 6)

          q_comp = mtot/m_comp
          if (q_comp>1) q_comp=q_comp**(-1)

          ! Mardling&Aarseth (2001) criterion check

          period_ratio = sqrt((a_comp*a_comp*a_comp)/(m_comp+mtot)/(semimajoraxis*semimajoraxis*semimajoraxis)*(mtot)) ! Po/Pi
          criterion = 4.7*(1-e_comp)**(-1.8)*(1+e_comp)**(0.6)*(1+q_comp)**(0.1)

          if (criterion > period_ratio) then
             print "(1x,a)",'WARNING: set_multiple: orbital parameters does not satisfy Mardling and Aarseth stability criterion.'
          endif

          q2=m2/m1
          mprimary = mtot/(1+q2)
          msecondary = mtot*q2/(1+q2)

          io=1
          exit
       endif
    enddo

    if (io == 0) then
       print "(1x,a)",'ERROR: set_multiple: star '//trim(hier_prefix)//' not present in HIERARCHY file.'
       ierr = ierr_missstar
    endif

    if (subst_index > 0 .and. subst_index <= size(xyzmh_ptmass(1,:))) then ! check for seg fault
       x_subst(:)=xyzmh_ptmass(1:3,subst_index)
       v_subst(:)=vxyz_ptmass(:,subst_index)
    endif

    period = sqrt(4.*pi**2*semimajoraxis**3/mtot)
 else
    mprimary = m1
    msecondary = m2

    if (present(posang_ascnode)) rel_posang_ascnode = posang_ascnode
    if (present(arg_peri)) rel_arg_peri= arg_peri
    if (present(incl)) rel_incl = incl

 endif
 !--- Create the binary
 call set_binary(mprimary,msecondary,semimajoraxis=semimajoraxis,eccentricity=eccentricity, &
            posang_ascnode=rel_posang_ascnode,arg_peri=rel_arg_peri,incl=rel_incl, &
            f=f,accretion_radius1=accretion_radius1,accretion_radius2=accretion_radius2, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass, ierr=ierr)

 if (present(subst) .and. subst>10) then
    !--- lower nptmass, copy one of the new sinks to the subst star
    nptmass = nptmass-1
    i1 = subst_index
    i2 = nptmass

    ! positions and accretion radii
    xyzmh_ptmass(1:6,i1) = xyzmh_ptmass(1:6,nptmass+1)

    ! velocities
    vxyz_ptmass(:,i1) = vxyz_ptmass(:,nptmass+1)

    !---
    ! Rotate the substituting binary with orientational parameters
    ! referring to the substituted star's orbital plane
    if (subst>0) then

       omega     = rel_arg_peri*pi/180.
       !big_omega = rel_posang_ascnode*pi/180.! + 0.5*pi
       inc       = rel_incl*pi/180.

       ! Retrieve eulerian angles of the substituted star orbit's semi-major axis (y axis)
       if (omega <= pi/2) then
          beta_y = omega
          sign_alpha=-1
          if (inc <= pi) then
             sign_gamma=1
          else
             sign_gamma=-1
          endif
       else
          beta_y = 2*pi-omega
          sign_alpha=1
          if (inc <= pi) then
             sign_gamma=-1
          else
             sign_gamma=1
          endif
       endif
       gamma_y=acos(sign_gamma*sin(beta_y)*sin(inc))
       alpha_y=acos(sign_alpha*sqrt(abs(sin(beta_y)**2-cos(gamma_y)**2))) ! Needs abs cause float approx for cos

       ! Retrieve eulerian angles of the axis perpendicular to the substituted star orbital plane (z axis)
       beta_z = pi/2.
       gamma_z = inc
       alpha_z = pi/2. - inc
       if (inc <= pi) then
          gamma_z=inc
          if (inc <= pi/2.) then
             alpha_z = pi/2.-inc
          elseif (inc > pi/2.) then
             alpha_z = inc-pi/2.
          endif
       elseif (inc < 2.*pi .and. inc > pi) then
          gamma_z = 2.*pi-inc
          if (inc <= 3.*pi/2.) then
             alpha_z = inc-pi/2
          elseif (inc > 3.*pi/2.) then
             alpha_z = 5.*pi/2.-inc
          endif
       endif

       
       ! Rotate substituting sinks by argument of pericenter around the z axis
       call gen_rotate(xyzmh_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, arg_peri*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, arg_peri*pi/180)
       call gen_rotate(xyzmh_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, arg_peri*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, arg_peri*pi/180)

       ! Rotate substituting sinks by inclination around the y axis
       call gen_rotate(xyzmh_ptmass(1:3,i1),alpha_y,beta_y,gamma_y, incl*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i1),alpha_y,beta_y,gamma_y, incl*pi/180)
       call gen_rotate(xyzmh_ptmass(1:3,i2),alpha_y,beta_y,gamma_y, incl*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i2),alpha_y,beta_y,gamma_y, incl*pi/180)

       ! Rotate substituting sinks by ascending node longitude around the z axis
       call gen_rotate(xyzmh_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i1),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
       call gen_rotate(xyzmh_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
       call gen_rotate(vxyz_ptmass(1:3,i2),alpha_z,beta_z,gamma_z, posang_ascnode*pi/180)
    endif
    
    ! Move the substituting binary's center of mass in the substituted star position
    xyzmh_ptmass(1:3,i1) = xyzmh_ptmass(1:3,i1)+x_subst
    xyzmh_ptmass(1:3,i2) = xyzmh_ptmass(1:3,i2)+x_subst
    ! Set the substituting binary's center of mass velocity
    vxyz_ptmass(:,i1) = vxyz_ptmass(:,i1)+v_subst
    vxyz_ptmass(:,i2) = vxyz_ptmass(:,i2)+v_subst

    ! Write updated HIERARCHY file with the two new stars and the substituted one
    call update_hierarchy_file(prefix, data, lines, hier_prefix, i1, i2, ierr) 
 endif

end subroutine set_multiple

!------------------------------------
! Rotate an (x,y,z) point by theta
! radiants around an axis with alpha,
! beta and gamma eulerian angles
!------------------------------------
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




end module sethierarchical


