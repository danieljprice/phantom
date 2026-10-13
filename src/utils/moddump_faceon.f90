!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Realign a disc on the xy plane and move the origin
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - outer_radius : *radius from centre of mass within which to measure L*
!   - system_type  : *1=single 2=binary 3=triple(inner) 4=triple(outer)*
!
! :Dependencies: moddump_utils, part, prompting,
!   vectorutils
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters (written to / read from the prefix.moddump file)
 integer :: system_type  = 2    ! 1=single/isink1, 2=binary, 3=triple (centred on binary), 4=triple (external)
 real    :: outer_radius = 400. ! radius from centre of mass within which to measure L

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:igas,xyzmh_ptmass,vxyz_ptmass,nptmass
 use vectorutils,   only:rotatevec,cross_product3D
 use moddump_utils, only:get_moddump_options

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real    :: radius,pmass
 real    :: Ltot(3),Lunit(3),z_axis(3),axis(3),angle
 real    :: centre_of_mass_sinks(3)

 Ltot = 0.

 pmass = massoftype(igas)

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 select case(system_type)
 case(1)
    centre_of_mass_sinks = xyzmh_ptmass(1:3,1)
 case(2)
    centre_of_mass_sinks = (xyzmh_ptmass(1:3,1)*xyzmh_ptmass(4,1)+xyzmh_ptmass(1:3,2)*xyzmh_ptmass(4,2))&
                               /(xyzmh_ptmass(4, 1)+xyzmh_ptmass(4, 2))
 case(3)
    centre_of_mass_sinks = (xyzmh_ptmass(1:3,2)*xyzmh_ptmass(4,2)+xyzmh_ptmass(1:3,3)*xyzmh_ptmass(4,3))&
                               /(xyzmh_ptmass(4, 2)+xyzmh_ptmass(4, 3))
 case(4)
    centre_of_mass_sinks = xyzmh_ptmass(1:3,1)
 end select

 ! Shift all the sink particles so our chosen centre of mass is at the origin
 do i = 1,nptmass
    xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - centre_of_mass_sinks
 enddo

 do i = 1,npart
    ! Shift the particles so that they're centred on our chosen centre of mass
    xyzh(1:3,i)=xyzh(1:3,i)-centre_of_mass_sinks

    radius = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    if (radius < outer_radius) then
       Ltot(1) = Ltot(1) + pmass*(xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i))
       Ltot(2) = Ltot(2) + pmass*(xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i))
       Ltot(3) = Ltot(3) + pmass*(xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i))
    endif
 enddo

 Lunit = Ltot / sqrt(Ltot(1)**2 + Ltot(2)**2 + Ltot(3)**2)
 z_axis = (/0.,0.,1./)

 call cross_product3D(Lunit,z_axis,axis)
 angle = acos(dot_product(Lunit,z_axis))

 ! Now we rotate everything about this axis
 do i = 1,npart
    call rotatevec(xyzh(1:3,i),axis,angle)
    call rotatevec(vxyzu(1:3,i),axis,angle)
 enddo

 do i = 1,nptmass
    call rotatevec(xyzmh_ptmass(1:3,i),axis,angle)
    call rotatevec(vxyz_ptmass(1:3,i),axis,angle)
 enddo

end subroutine modify_dump

!-----------------------------------------------------------------------
!+
!  set parameters interactively (when no .moddump file is found)
!+
!-----------------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('System (1=single, 2=binary, 3=triple about binary, 4=triple external)',system_type,1,4)
 call prompt('Outer radius within which to measure angular momentum',outer_radius,0.)

end subroutine moddump_interactive

!-----------------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!-----------------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(system_type,'system_type',db,errcount=nerr,min=1,max=4)
 call read_inopt(outer_radius,'outer_radius',db,errcount=nerr,min=0.)

end subroutine read_options_moddump

!-----------------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!-----------------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# face-on realignment parameters'
 call write_inopt(system_type,'system_type','1=single 2=binary 3=triple(inner) 4=triple(outer)',iunit)
 call write_inopt(outer_radius,'outer_radius','radius from centre of mass within which to measure L',iunit)

end subroutine write_options_moddump

end module moddump
