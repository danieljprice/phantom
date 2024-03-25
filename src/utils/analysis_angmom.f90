!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! compute the angular momentum budget of a simulation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, energies, part, physcon, units, vectorutils
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'angmom'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use energies, only:compute_energies,angtot
 use part,     only:nptmass,xyzmh_ptmass,vxyz_ptmass,ispinx,ispinz,igas,massoftype
 use units,    only:unit_angmom,utime
 use physcon,  only:years
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 logical, save      :: first = .true.
 real    :: Lhat(3),Ltot(3),Ltot_mag,inc,rot
 real    :: Ltot_sink(3),Lhat_sink(3),Lsink_mag,inc_sink,rot_sink
 real    :: Lspin(3),Lspin_mag,spini(3),L_total(3),L_total_mag
 real    :: dx(3),sep,mgas,msinks,mtot
 integer :: i,iu

! Print the analysis being done
 write(*,'(1x,"Performing analysis type ",a)') analysistype
 write(*,'(1x,"Input file name is ",a,/)') dumpfile
 !
 ! get angular momentum of the gas and the sink particles
 !
 call get_angmom(npart,xyzh,vxyzu,Lhat,Ltot,Ltot_mag,inc,rot,'gas')
 call get_angmom(nptmass,xyzmh_ptmass,vxyz_ptmass,Lhat_sink,Ltot_sink,Lsink_mag,&
                 inc_sink,rot_sink,'sink')

 Lspin = 0.
 do i=1,nptmass
    spini = xyzmh_ptmass(ispinx:ispinz,i)
    print "(' L_spin (',i2,') = ',g0)" ,i,sqrt(dot_product(spini,spini))
    Lspin(:) = Lspin(:) + spini
 enddo
 Lspin_mag = sqrt(dot_product(Lspin,Lspin))
 !
 ! get total angular momentum from the main code routine
 !
 call compute_energies(time)

 print*
 print*,' L_gas   = ',Ltot_mag,' dir = ',Lhat
 print*,' L_sinks = ',Lsink_mag,' dir = ',Lhat_sink
 print*,' L_spin  = ',Lspin_mag,' dir = ',Lspin/Lspin_mag

 if (nptmass >= 2) then
    dx = xyzmh_ptmass(1:3,2) - xyzmh_ptmass(1:3,1)
    sep = sqrt(dot_product(dx,dx))
    print*,' sink separation = ',sep,' mass = ',xyzmh_ptmass(4,1),xyzmh_ptmass(4,2)
 endif

 L_total = Ltot + Ltot_sink + Lspin
 L_total_mag = sqrt(dot_product(L_total,L_total))
 print*,' L_total (from adding the above) = ',L_total_mag
 print*,' L_total (from energies routine) = ',angtot

 mgas = massoftype(igas)*npart
 msinks = sum(xyzmh_ptmass(4,1:nptmass))
 mtot = mgas + msinks
 print*,' total mass in gas = ',mgas
 print*,' total mass in sinks = ',msinks
 print*,' total mass = ',mtot

 ! Write angular momentum information
 if (first) then
    first = .false.
    open(newunit=iu, file='angmom.ev',status='replace')
    write(iu,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',&
          2,'L_{gas}',  &
          3,'L_{orbit}',  &
          4,'L_{spin}',  &
          5,'L_{total}'
 else
    open(newunit=iu, file='angmom.ev',position='append')
 endif
 write(iu,'(6(es18.10,1X))') time*utime/years,Ltot_mag*unit_angmom,Lsink_mag*unit_angmom,&
       Lspin_mag*unit_angmom,L_total_mag*unit_angmom
 close(iu)

end subroutine do_analysis

!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine get_angmom(n,xyz_arr,vxyz,Lhat,Ltot,Lmag,inc,rot,type)
 use physcon,     only:pi
 use vectorutils, only:cross_product3D
 use part,        only:isdead_or_accreted,massoftype,igas,iamtype,iphase
 use dim,         only:maxphase,maxp
 integer, intent(in) :: n
 real, intent(in)    :: xyz_arr(:,:),vxyz(:,:)
 real, intent(out)   :: Lhat(3),Ltot(3),Lmag,inc,rot
 character(len=*), intent(in) :: type
 integer :: i
 real    :: Li(3),massi

 Ltot = 0.
 select case(type)
 case('sink')
    do i=1,n
       massi = xyz_arr(4,i)  ! xyzmh_ptmass
       call cross_product3D(xyz_arr(1:3,i),vxyz(1:3,i),Li)
       Ltot = Ltot + massi*Li
    enddo
 case default
    massi = massoftype(igas)
    do i=1,n
       if (.not.isdead_or_accreted(xyz_arr(4,i))) then
          if (maxphase==maxp) massi = massoftype(iamtype(iphase(i)))
          call cross_product3D(xyz_arr(1:3,i),vxyz(1:3,i),Li)
          Ltot = Ltot + massi*Li
       endif
    enddo
 end select

 Lmag = sqrt(dot_product(Ltot,Ltot))
 Lhat = Ltot/Lmag
 inc  = acos(dot_product( Lhat,(/0.,0.,1./)))*180./pi ! Angle from +z axis -- should always be 0<inc<180 degrees
 rot  = atan2(Lhat(2),Lhat(1))*180./pi                ! Angle around in xy plane (from +x axis) -- should always -180<rot<180

end subroutine get_angmom

end module analysis
