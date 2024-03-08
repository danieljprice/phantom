!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part, prompting, ptmass_heating, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,ihsoft,ilum
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass
 use ptmass_heating, only:Lnuc
 use units,          only:unit_energ,utime
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 integer                :: i,isinkpart
 real                   :: racc,hsoft,mass,mass_old,newx,Lnuc_cgs
 logical                :: iresetCM

 print*,'Sink particles in dump:'
 do i=1,nptmass
    print "(a,1x,i4,a)",'Sink',i,':'
    print "(7(a5,1x,a,1x,f13.7,/))",&
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

 isinkpart = 2
 do while (isinkpart /= 0)
    call prompt('Enter the sink particle number to modify (0 to exit):',isinkpart,0,nptmass)
    if (isinkpart <= 0) exit

    mass = xyzmh_ptmass(4,isinkpart)
    mass_old = mass
    call prompt('Enter new mass for the sink:',mass,0.)
    print*,'Mass changed to ',mass
    xyzmh_ptmass(4,isinkpart) = mass

    racc = xyzmh_ptmass(ihacc,isinkpart)
    ! rescaling accretion radius for updated mass
    racc = racc * (mass/mass_old)**(1./3)
    call prompt('Enter new accretion radius for the sink:',racc,0.)
    print*,'Accretion radius changed to ',racc
    xyzmh_ptmass(ihacc,isinkpart) = racc

    hsoft = xyzmh_ptmass(ihsoft,isinkpart)
    call prompt('Enter new softening length for the sink:',hsoft,0.)
    print*,'Softening length changed to ',hsoft
    xyzmh_ptmass(ihsoft,isinkpart) = hsoft

    newx = xyzmh_ptmass(1,isinkpart)
    call prompt('Enter new x-coordinate for the sink in code units:',newx,0.)
    xyzmh_ptmass(1,isinkpart) = newx
    print*,'x-coordinate changed to ',xyzmh_ptmass(1,isinkpart)

    Lnuc = xyzmh_ptmass(ilum,isinkpart)
    Lnuc_cgs = Lnuc * unit_energ / utime
    call prompt('Enter new sink heating luminosity in erg/s:',Lnuc_cgs,0.)
    xyzmh_ptmass(ilum,isinkpart) = Lnuc_cgs / unit_energ * utime
    print*,'Luminosity [erg/s] changed to ',xyzmh_ptmass(ilum,isinkpart) * unit_energ / utime
 enddo

 iresetCM = .false.
 call prompt('Reset centre of mass?',iresetCM)
 if (iresetCM) call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 return
end subroutine modify_dump

end module moddump
