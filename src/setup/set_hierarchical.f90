!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module sethierarchical
!
! setup for general hierarchical systems (binaries, triples, quadruples, etc.)
! this is specified in the form of a hierarchy string:
!
!  e.g. hierarchy = '111,112,121,1221,1222'
!
! where each component of the multiple system is given a number
!
! :References:
!  - Ceppi et al. (2022) MNRAS 514, 906
!  - Ceppi et al. (2023) MNRAS 520, 5817
!
! :Owner: Simone Ceppi
!
! :Runtime parameters:
!   - hierarchy : *string definining the hierarchy (e.g. 111,112,121,1221,1222)*
!
! :Dependencies: infile_utils, setbinary, sethier_utils
!
 use sethier_utils, only:process_hierarchy,max_hier_levels,lenhierstring,hierarchical_system,&
      hier_db_size,hier_db_prop,recursive_splitting, gen_rotate, find_hier_level_orb_elem, &
      find_hierarchy_index, find_data_index, find_ptmass_index,&
      load_hierarchy_file,update_hierarchy_file,check_substitution,&
      ierr_HIER2,ierr_subststar,ierr_Omegasubst,ierr_missstar,pi
 implicit none

 public :: set_hierarchical, set_multiple, set_hier_multiple
 public :: write_hierarchical_setupfile
 public :: read_hierarchical_setupfile
 public :: set_hierarchical_default_options
 public :: get_hierarchical_level_com
 public :: get_hier_level_mass
 public :: print_chess_logo, generate_hierarchy_string

 character(len=lenhierstring) :: hierarchy = '111,112,121,1221,1222'

 type(hierarchical_system) :: hs
 public :: hierarchy
 public :: hs

 private

contains

!--------------------------------------------------------------------------
!+
!  routine to actually place the particles based on configuration
!  specified in the hierarchy file
!+
!--------------------------------------------------------------------------
subroutine set_hierarchical(prefix, nptmass, xyzmh_ptmass, vxyz_ptmass, ierr)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 integer, intent(out)   :: ierr
 character(len=20), intent(in) :: prefix

 integer :: i, io
 real :: m1, m2, accr1, accr2, binary_a, binary_e, binary_i, binary_O, binary_w, binary_f

 integer :: splits, sink_num_temp, subst
 character(len=10) :: sink_list(max_hier_levels), split_list(max_hier_levels)
 character(len=20) :: hl_temp

 call print_chess_logo()

 splits = 0
 sink_list = hs%labels%sink
 sink_num_temp = hs%labels%sink_num

 call recursive_splitting(sink_num_temp, sink_list, split_list, splits)

 do i=splits,1,-1

    hl_temp = trim(split_list(i))

    call find_hier_level_orb_elem(hl_temp, hs, m1, m2, accr1, accr2, &
                                  binary_a, binary_e, binary_i, binary_O, &
                                  binary_w, binary_f)
    !print*,'elements of ', hl_temp, ' : ', m1, m2, accr1, accr2, &
    !                              binary_a, binary_e, binary_i, binary_O, &
    !                              binary_w, binary_f

    read(hl_temp,*,iostat=io) subst

    !print*, 'passing subst = ', subst
    call set_hier_multiple(m1,m2,semimajoraxis=binary_a,eccentricity=binary_e, &
            posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
            f=binary_f,accretion_radius1=accr1,accretion_radius2=accr2, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,&
            ierr=ierr, subst=subst, prefix=prefix)

 enddo

end subroutine set_hierarchical

!--------------------------------------------------------------------------
!+
!  write options to .setup file
!+
!--------------------------------------------------------------------------
subroutine write_hierarchical_setupfile(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit
 integer :: i

 write(iunit,"(/,a)") '# options for hierarchical system'
 call write_inopt(hierarchy, 'hierarchy','string definining the hierarchy (e.g. 111,112,121,1221,1222)', iunit)

 hs%labels = process_hierarchy(hierarchy)

 write(iunit,"(/,a)") '# sink properties'
 do i=1,hs%labels%sink_num
    call write_inopt(hs%sinks(i)%mass, trim(hs%labels%sink(i))//'_mass',&
                     'mass of object '//trim(hs%labels%sink(i)), iunit)
 enddo
 do i=1,hs%labels%sink_num
    call write_inopt(hs%sinks(i)%accr, trim(hs%labels%sink(i))//'_accr',&
                     'accretion radius for object '//trim(hs%labels%sink(i)), iunit)
 enddo

 write(iunit,"(/,a)") '# orbit properties'
 do i=1,hs%labels%hl_num
    write(iunit,"(a)") '# binary '//trim(hs%labels%hl(i))
    call write_inopt(hs%levels(i)%a, trim(hs%labels%hl(i))//'_a',&
         'semi-major axis for binary '//trim(hs%labels%hl(i)),iunit)
    call write_inopt(hs%levels(i)%e, trim(hs%labels%hl(i))//'_e',&
         'eccentricity for binary '//trim(hs%labels%hl(i)),iunit)
    call write_inopt(hs%levels(i)%inc, trim(hs%labels%hl(i))//'_i',&
         'i [deg] inclination for binary '//trim(hs%labels%hl(i)),iunit)
    call write_inopt(hs%levels(i)%O, trim(hs%labels%hl(i))//'_O',&
         'Omega [deg] PA of ascending node for binary '//trim(hs%labels%hl(i)),iunit)
    call write_inopt(hs%levels(i)%w, trim(hs%labels%hl(i))//'_w',&
         'w [deg] argument of periapsis for binary '//trim(hs%labels%hl(i)),iunit)
    call write_inopt(hs%levels(i)%f, trim(hs%labels%hl(i))//'_f',&
         'f [deg] true anomaly for binary '//trim(hs%labels%hl(i)),iunit)
 enddo

end subroutine write_hierarchical_setupfile

!--------------------------------------------------------------------------
!+
!  read options from .setup file
!+
!--------------------------------------------------------------------------
subroutine read_hierarchical_setupfile(db, nerr)
 use infile_utils, only:read_inopt, inopts
 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(inout) :: nerr
 integer :: i

 call read_inopt(hierarchy,'hierarchy',db,errcount=nerr)

 hs%labels = process_hierarchy(hierarchy)

 do i=1,hs%labels%sink_num
    call read_inopt(hs%sinks(i)%mass, trim(hs%labels%sink(i))//'_mass',db,errcount=nerr)
    call read_inopt(hs%sinks(i)%accr, trim(hs%labels%sink(i))//'_accr',db,errcount=nerr)
 enddo

 do i=1,hs%labels%hl_num
    call read_inopt(hs%levels(i)%a, trim(hs%labels%hl(i))//'_a',db,errcount=nerr)
    call read_inopt(hs%levels(i)%e, trim(hs%labels%hl(i))//'_e',db,errcount=nerr)
    call read_inopt(hs%levels(i)%inc, trim(hs%labels%hl(i))//'_i',db,errcount=nerr)
    call read_inopt(hs%levels(i)%O, trim(hs%labels%hl(i))//'_O',db,errcount=nerr)
    call read_inopt(hs%levels(i)%w, trim(hs%labels%hl(i))//'_w',db,errcount=nerr)
    call read_inopt(hs%levels(i)%f, trim(hs%labels%hl(i))//'_f',db,errcount=nerr)
 enddo

end subroutine read_hierarchical_setupfile

!--------------------------------------------------------------------------
!+
!  set default options for hierarchical setup
!+
!--------------------------------------------------------------------------
subroutine set_hierarchical_default_options()
 integer :: i

 hs%labels = process_hierarchy(hierarchy)

 do i=1,hs%labels%sink_num
    hs%sinks(i)%mass=1.
    hs%sinks(i)%accr=1.
 enddo

 do i=1,hs%labels%hl_num
    hs%levels(i)%a=10000./10**len(trim(hs%labels%hl(i)))
    hs%levels(i)%e=0.
    hs%levels(i)%inc=0.
    hs%levels(i)%O=0.
    hs%levels(i)%w=0.
    hs%levels(i)%f=0.
 enddo

end subroutine set_hierarchical_default_options

!--------------------------------------------------------------------------
!+
!  compute position and velocity of hierarchical level
!+
!--------------------------------------------------------------------------
subroutine get_hierarchical_level_com(level, xorigin, vorigin, xyzmh_ptmass, vxyz_ptmass, prefix)
 character(len=10), intent(in)  :: level
 real,              intent(out) :: xorigin(:), vorigin(:)
 real,              intent(in)  :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 character(len=20), intent(in)  :: prefix

 integer :: int_sinks(max_hier_levels), inner_sinks_num, i
 real :: mass

 call find_hierarchy_index(level, int_sinks, inner_sinks_num, prefix)
 !print*, 'found ', inner_sinks_num, ' sinks: ', int_sinks(:inner_sinks_num)

 xorigin = 0.
 vorigin = 0.
 mass = 0.

 !print*, xyzmh_ptmass(4,:5)

 do i=1,inner_sinks_num
    xorigin(:) = xorigin(:)+xyzmh_ptmass(4,int_sinks(i))*xyzmh_ptmass(1:3,int_sinks(i))
    vorigin(:) = vorigin(:)+xyzmh_ptmass(4,int_sinks(i))*vxyz_ptmass(1:3,int_sinks(i))
    mass = mass + xyzmh_ptmass(4,int_sinks(i))
 enddo
 !print*, '--> ', mass, xorigin, vorigin
 xorigin = xorigin/mass
 vorigin = vorigin/mass

end subroutine get_hierarchical_level_com

!--------------------------------------------------------------------------
!+
!  compute the total mass of an hierarchical level
!+
!--------------------------------------------------------------------------
real function get_hier_level_mass(level)
 use sethier_utils, only:get_hierarchical_level_mass
 character(len=10), intent(in) :: level

 get_hier_level_mass = get_hierarchical_level_mass(level, hs)

end function get_hier_level_mass

!----------------------------------------------------------------
!+
!  setup for a multiple, using set_binary
!+
!----------------------------------------------------------------
subroutine set_hier_multiple(m1,m2,semimajoraxis,eccentricity, &
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

 integer :: i1,i2,subst_index
 real    :: mtot,period
 real    :: x_subst(3)=0,v_subst(3)=0
 real    :: omega,inc
 !logical :: do_verbose

 real, dimension(hier_db_size,hier_db_prop) :: data
 character(len=20)      :: hier_prefix
 integer                :: lines
 real                   :: rel_posang_ascnode=0.,rel_arg_peri=0.,rel_incl=0.
 real                   :: q2
 real                   :: alpha_y, beta_y, gamma_y, alpha_z, beta_z, gamma_z, sign_alpha, sign_gamma

 ierr = 0
 !--- Load/Create HIERARCHY file: xyzmh_ptmass index | hierarchical index | star mass | companion star mass | semi-major axis | eccentricity | period | inclination | argument of pericenter | ascending node longitude
 if (present(subst) .and. subst > 10) then
    call load_hierarchy_file(prefix, data, lines, ierr)
 else
    mtot = m1 + m2
    period = sqrt(4.*pi**2*semimajoraxis**3/mtot)

    hier_prefix = '1'

    lines = 0
    call update_hierarchy_file(prefix, hs, data, lines, hier_prefix, 1, 2, ierr)
 endif
 subst_index = 0
 !--- Checks to avoid bad substitutions
 if (present(subst) .and. subst > 10) then
    write(hier_prefix, *) subst

    call check_substitution(hier_prefix, semimajoraxis, prefix, ierr)
    !call check_substitution('11                  ', semimajoraxis, prefix, ierr)

    !print*,data(:,1)
    !print*,data(:,2)
    !print*,hier_prefix


    call find_data_index(hier_prefix, subst_index, prefix, ierr) ! QUI VOGLIO DATA INDEX

    data(subst_index,1) = 0

    call find_ptmass_index(hier_prefix, subst_index, prefix, ierr) ! QUI VOGLIO PTMASS INDEX

    q2=m2/m1

    if (subst_index > 0 .and. subst_index <= size(xyzmh_ptmass(1,:))) then ! check for seg fault
       !print*, '!!!!!!!!!!!!!!!!!!! ', subst_index
       x_subst(:)=xyzmh_ptmass(1:3,subst_index)
       v_subst(:)=vxyz_ptmass(:,subst_index)
    endif

    mtot = m1 + m2
    period = sqrt(4.*pi**2*semimajoraxis**3/mtot)
 else

    if (present(posang_ascnode)) rel_posang_ascnode = posang_ascnode
    if (present(arg_peri)) rel_arg_peri= arg_peri
    if (present(incl)) rel_incl = incl

 endif

 !if (subst==111) then
 !   print*, 'quiiiit'
 !   !return
 !endif
 !print*, 'set_binary args: ', mprimary, msecondary
 !--- Create the binary
 call set_binary(m1,m2,semimajoraxis=semimajoraxis,eccentricity=eccentricity, &
            posang_ascnode=rel_posang_ascnode,arg_peri=rel_arg_peri,incl=rel_incl, &
            f=f,accretion_radius1=accretion_radius1,accretion_radius2=accretion_radius2, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass, ierr=ierr)

 if (present(subst) .and. subst>10) then
    !--- lower nptmass, copy one of the new sinks to the subst star
    nptmass = nptmass-1
    i2 = subst_index
    i1 = nptmass
    !print *,'--> ',i1,i2
    ! positions and accretion radii
    xyzmh_ptmass(:,i2) = xyzmh_ptmass(:,nptmass+1)

    ! velocities
    vxyz_ptmass(:,i2) = vxyz_ptmass(:,nptmass+1)

    !
    ! Rotate the substituting binary with orientational parameters
    ! referring to the substituted star's orbital plane
    !
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
    call update_hierarchy_file(prefix, hs, data, lines, hier_prefix, i1, i2, ierr)
    !print*,'---------------'
    !print*,x_subst
    !call get_hierarchical_level_com(hier_prefix, x_subst, v_subst, xyzmh_ptmass, vxyz_ptmass, prefix)
    !print*,x_subst
    !print*,'---------------'
 endif

end subroutine set_hier_multiple

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

 real, dimension(24,10) :: data
 character(len=20)      :: hier_prefix, filename
 logical                :: iexist
 integer                :: io, lines
 real                   :: period_ratio,criterion,q_comp,a_comp,e_comp,m_comp
 real                   :: rel_posang_ascnode=0.,rel_arg_peri=0.,rel_incl=0.
 real                   :: q2,mprimary,msecondary
 real                   :: alpha_y, beta_y, gamma_y, alpha_z, beta_z, gamma_z, sign_alpha, sign_gamma

 ierr = 0
 !do_verbose = .true.
 !if (present(verbose)) do_verbose = verbose

 if (present(prefix)) then
    filename = trim(prefix)//'.hierarchy'
 else
    filename = 'HIERARCHY'
 endif

 !--- Load/Create HIERARCHY file: xyzmh_ptmass index | hierarchical index | star mass | companion star mass | semi-major axis | eccentricity | period | inclination | argument of pericenter | ascending node longitude
 inquire(file=trim(filename), exist=iexist)
 if (present(subst)) then
    if (subst>10) then
       if (iexist) then
          open(1,file=trim(filename),status='old')
          lines=0
          do
             read(1, *,iostat=io) data(lines+1,:)
             if (io/=0) exit
             lines = lines + 1
          enddo
          close(1)
       else
          print "(1x,a)",'ERROR: set_multiple: there is no HIERARCHY file, cannot perform subtitution.'
          ierr = ierr_HIER2
       endif
    endif
 else
    if (iexist) then
       print "(1x,a)",'WARNING: set_multiple: deleting an existing HIERARCHY file.'
       open(1,file=trim(filename),status='old')
       close(1,status='delete')
    endif

    mtot = m1 + m2
    period = sqrt(4.*pi**2*semimajoraxis**3/mtot)

    open(1,file=trim(filename),status='new')
    if (present(incl)) then
       if (present(posang_ascnode) .and. present(arg_peri)) then
          write(1,*) 1, 11, m1, m2, semimajoraxis, eccentricity, period, incl, arg_peri, posang_ascnode
          write(1,*) 2, 12, m2, m1, semimajoraxis, eccentricity, period, incl, arg_peri, posang_ascnode
       else ! set binary at apastron with inclination
          write(1,*) 1, 11, m1, m2, semimajoraxis, eccentricity, period, incl, 0, 0
          write(1,*) 2, 12, m2, m1, semimajoraxis, eccentricity, period, incl, 0, 0
       endif
    else ! set binary at apastron without inclination
       write(1,*) 1, 11, m1, m2, semimajoraxis, eccentricity, period, 0, 0, 0
       write(1,*) 2, 12, m2, m1, semimajoraxis, eccentricity, period, 0, 0, 0
    endif
    close(1)
 endif

 subst_index = 0
 !--- Checks to avoid bad substitutions
 if (present(subst)) then
    if (subst>10) then
       write(hier_prefix, *) subst
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

             ! Mardling & Aarseth (2001) criterion check
             period_ratio = sqrt((a_comp*a_comp*a_comp)/(m_comp+mtot)/&
                                 (semimajoraxis*semimajoraxis*semimajoraxis)*(mtot)) ! Po/Pi
             criterion = 4.7*(1-e_comp)**(-1.8)*(1+e_comp)**(0.6)*(1+q_comp)**(0.1)

             if (criterion > period_ratio) then
                print "(1x,a)",'WARNING: set_multiple: orbital parameters do not satisfy '//&
                               'Mardling & Aarseth stability criterion.'
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
       !i1 = subst_index
       !i2 = nptmass + 1
       !nptmass = nptmass + 1

       period = sqrt(4.*pi**2*semimajoraxis**3/mtot)
    endif
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

 if (present(subst)) then
    if (subst>10) then
       !--- lower nptmass, copy one of the new sinks to the subst star
       nptmass = nptmass-1
       i1 = subst_index
       i2 = nptmass

       ! positions and accretion radii
       xyzmh_ptmass(1:6,i1) = xyzmh_ptmass(1:6,nptmass+1)

       ! test Jolien
!    print "(5(2x,a,g12.3,/),2x,a,g12.3)", &
!    'i1     :',i1, &
!     'mass i1:',xyzmh_ptmass(4,i1), &
!     'i2     :',i2, &
!     'mass i2:',xyzmh_ptmass(4,i2)

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
       open(1,file=trim(filename),status='old')
       do i=1,lines
          write(1,*) int(data(i,1)), int(data(i,2)), data(i,3:)
       enddo
       write(1,*) i1, trim(hier_prefix)//"1", mprimary, msecondary, semimajoraxis, eccentricity, &
         period, incl, arg_peri, posang_ascnode
       write(1,*) i2, trim(hier_prefix)//"2", msecondary, mprimary, semimajoraxis, eccentricity, &
         period, incl, arg_peri, posang_ascnode
       close(1)
    endif
 endif

end subroutine set_multiple

subroutine generate_hierarchy_string(nsinks)
 integer, intent(in) :: nsinks

 integer :: i, pos
 character(len=10) :: label

 hierarchy = '11,12'

 do i=1,nsinks-2
    pos = scan(hierarchy, ',', .true.)

    label = trim(hierarchy(pos+1:))

    hierarchy = trim(hierarchy(:pos-1))//','//trim(label)//'1,'//trim(label)//'2'

    !print*,label
 enddo

end subroutine generate_hierarchy_string


subroutine print_chess_logo()!id)
 !use io,               only:master
 !integer,           intent(in) :: id

! if (id==master) then
 print*,"                                                                      "
 print*,"                                                       _:_            "
 print*,"                                                      '-.-'           "
 print*,"                                             ()      __.'.__          "
 print*,"                                          .-:--:-.  |_______|         "
 print*,"                                   ()      \____/    \=====/          "
 print*,"                                   /\      {====}     )___(           "
 print*,"                        (\=,      //\\      )__(     /_____\          "
 print*,"        __    |'-'-'|  //  .\    (    )    /____\     |   |           "
 print*,"       /  \   |_____| (( \_  \    )__(      |  |      |   |           "
 print*,"       \__/    |===|   ))  `\_)  /____\     |  |      |   |           "
 print*,"      /____\   |   |  (/     \    |  |      |  |      |   |           "
 print*,"       |  |    |   |   | _.-'|    |  |      |  |      |   |           "
 print*,"       |__|    )___(    )___(    /____\    /____\    /_____\          "
 print*,"      (====)  (=====)  (=====)  (======)  (======)  (=======)         "
 print*,"      }===={  }====={  }====={  }======{  }======{  }======={         "
 print*,"     (______)(_______)(_______)(________)(________)(_________)        "
 print*,"                                                                      "
 print*,"          _             _       _    _           _           _        "
 print*,"        /\ \           / /\    / /\ /\ \        / /\        / /\      "
 print*,"       /  \ \         / / /   / / //  \ \      / /  \      / /  \     "
 print*,"      / /\ \ \       / /_/   / / // /\ \ \    / / /\ \__  / / /\ \__  "
 print*,"     / / /\ \ \     / /\ \__/ / // / /\ \_\  / / /\ \___\/ / /\ \___\ "
 print*,"    / / /  \ \_\   / /\ \___\/ // /_/_ \/_/  \ \ \ \/___/\ \ \ \/___/ "
 print*,"   / / /    \/_/  / / /\/___/ // /____/\      \ \ \       \ \ \       "
 print*,"  / / /          / / /   / / // /\____\/  _    \ \ \  _    \ \ \      "
 print*," / / /________  / / /   / / // / /______ /_/\__/ / / /_/\__/ / /      "
 print*,"/ / /_________\/ / /   / / // / /_______\\ \/___/ /  \ \/___/ /       "
 print*,"\/____________/\/_/    \/_/ \/__________/ \_____\/    \_____\/        "
 print*,"                                                                      "

 print "(/,65('-'),1(/,a),/,65('-'),/)",&
         '  Welcome to CHESS (Complete Hierarchical Endless System Setup)'


 !    print "(/,65('-'),1(/,a),/,1(a),/,65('-'),/)",&
 !         '  Welcome to CHESS (Complete Hierarchical Endless System Setup)', &
 !         '        simulate the universe as a hierarchical system'

! endif
end subroutine print_chess_logo


end module sethierarchical
