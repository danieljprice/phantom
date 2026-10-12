!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Interactively change sink particle properties
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters:
!   - delete_sink : *delete the sink instead of modifying it*
!   - hsoft       : *new softening length for the sink [code units]*
!   - isinkpart   : *sink particle number to modify (0 = none)*
!   - lnuc_cgs    : *new sink heating luminosity [erg/s]*
!   - mass        : *new mass for the sink [code units]*
!   - newx        : *new x-coordinate for the sink [code units]*
!   - racc        : *new accretion radius for the sink [code units]*
!   - reset_cm    : *reset centre of mass*
!
! :Dependencies: centreofmass, moddump_utils, part,
!   prompting, ptmass_heating, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer :: isinkpart   = 0        ! sink particle number to modify (0 = none)
 logical :: delete_sink = .false.  ! delete the sink instead of modifying it
 logical :: reset_CM    = .false.  ! reset centre of mass
 real    :: mass     = 0.          ! new mass [code units]
 real    :: racc     = 0.          ! new accretion radius [code units]
 real    :: hsoft    = 0.          ! new softening length [code units]
 real    :: newx     = 0.          ! new x-coordinate [code units]
 real    :: Lnuc_cgs = 0.          ! new sink heating luminosity [erg/s]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,ihsoft,ilum
 use centreofmass,   only:reset_centreofmass
 use units,          only:unit_energ,utime
 use moddump_utils,  only:get_moddump_options
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 integer                :: i

 print*,'Sink particles in dump:'
 do i=1,nptmass
    print "(a,1x,i4,a)",'Sink',i,':'
    print "(7(a5,1x,a,1x,es24.16e3,/))",&
             'x','=',xyzmh_ptmass(1,i),&
             'y','=',xyzmh_ptmass(2,i),&
             'z','=',xyzmh_ptmass(3,i),&
             'mass','=',xyzmh_ptmass(4,i),&
             'h','=',xyzmh_ptmass(ihsoft,i),&
             'hacc','=',xyzmh_ptmass(ihacc,i),&
             'Lnuc','=',xyzmh_ptmass(ilum,i)
    if (i > 10) then
       print*, "The rest of the sink particles are not displayed"
       exit
    endif
 enddo
 !
 ! read parameters from the .moddump file, prompting if it does not exist
 ! (edits a single sink; rerun to modify additional sinks)
 !
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 if (isinkpart > 0 .and. isinkpart <= nptmass) then
    if (delete_sink) then
       if (isinkpart==nptmass) then
          xyzmh_ptmass(:,isinkpart) = 0.
          vxyz_ptmass(:,isinkpart) = 0.
       else
          xyzmh_ptmass(:,isinkpart:nptmass-1) = xyzmh_ptmass(:,isinkpart+1:nptmass)
          vxyz_ptmass(:,isinkpart:nptmass-1) = vxyz_ptmass(:,isinkpart+1:nptmass)
       endif
       nptmass = nptmass - 1
    else
       xyzmh_ptmass(4,isinkpart) = mass
       print*,'Mass changed to ',mass
       xyzmh_ptmass(ihacc,isinkpart) = racc
       print*,'Accretion radius changed to ',racc
       xyzmh_ptmass(ihsoft,isinkpart) = hsoft
       print*,'Softening length changed to ',hsoft
       xyzmh_ptmass(1,isinkpart) = newx
       print*,'x-coordinate changed to ',xyzmh_ptmass(1,isinkpart)
       xyzmh_ptmass(ilum,isinkpart) = Lnuc_cgs / unit_energ * utime
       print*,'Luminosity [erg/s] changed to ',xyzmh_ptmass(ilum,isinkpart) * unit_energ / utime
    endif
 endif

 if (reset_CM) call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting,      only:prompt
 use part,           only:xyzmh_ptmass,nptmass,ihacc,ihsoft,ilum
 use ptmass_heating, only:Lnuc
 use units,          only:unit_energ,utime
 real :: mass_old

 isinkpart = 2
 call prompt('Enter the sink particle number to modify (0 to exit):',isinkpart,0,nptmass)
 if (isinkpart <= 0) return

 call prompt('Delete sink?',delete_sink,.false.)
 if (.not. delete_sink) then
    mass = xyzmh_ptmass(4,isinkpart)
    mass_old = mass
    call prompt('Enter new mass for the sink:',mass,0.)

    racc = xyzmh_ptmass(ihacc,isinkpart)
    ! rescaling accretion radius for updated mass
    if (mass_old > 0.) racc = racc * (mass/mass_old)**(1./3)
    call prompt('Enter new accretion radius for the sink:',racc,0.)

    hsoft = xyzmh_ptmass(ihsoft,isinkpart)
    call prompt('Enter new softening length for the sink:',hsoft,0.)

    newx = xyzmh_ptmass(1,isinkpart)
    call prompt('Enter new x-coordinate for the sink in code units:',newx)

    Lnuc = xyzmh_ptmass(ilum,isinkpart)
    Lnuc_cgs = Lnuc * unit_energ / utime
    call prompt('Enter new sink heating luminosity in erg/s:',Lnuc_cgs,0.)
 endif

 call prompt('Reset centre of mass?',reset_CM)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(isinkpart,'isinkpart','sink particle number to modify (0 = none)',iunit)
 call write_inopt(delete_sink,'delete_sink','delete the sink instead of modifying it',iunit)
 call write_inopt(mass,'mass','new mass for the sink [code units]',iunit)
 call write_inopt(racc,'racc','new accretion radius for the sink [code units]',iunit)
 call write_inopt(hsoft,'hsoft','new softening length for the sink [code units]',iunit)
 call write_inopt(newx,'newx','new x-coordinate for the sink [code units]',iunit)
 call write_inopt(Lnuc_cgs,'lnuc_cgs','new sink heating luminosity [erg/s]',iunit)
 call write_inopt(reset_CM,'reset_cm','reset centre of mass',iunit)

end subroutine write_options_moddump

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(isinkpart,'isinkpart',db,min=0,errcount=nerr)
 call read_inopt(delete_sink,'delete_sink',db,errcount=nerr)
 call read_inopt(mass,'mass',db,min=0.,errcount=nerr)
 call read_inopt(racc,'racc',db,min=0.,errcount=nerr)
 call read_inopt(hsoft,'hsoft',db,min=0.,errcount=nerr)
 call read_inopt(newx,'newx',db,errcount=nerr)
 call read_inopt(Lnuc_cgs,'lnuc_cgs',db,min=0.,errcount=nerr)
 call read_inopt(reset_CM,'reset_cm',db,errcount=nerr)

end subroutine read_options_moddump

end module moddump
