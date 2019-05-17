!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!   Compute total gas and dust mass when using one fluid dust algorithm
!   Use this to check the conservation of total gas and total dust mass
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, options, part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'dustmass'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,     only:maxp,maxdusttypes
 use part,    only:maxphase,isdead_or_accreted,dustfrac,massoftype,igas,&
                   iphase,iamtype
 use options, only:use_dustfrac
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real          :: Mtot,Mgas,Mdust1,Mdust2,Macc,pmassi
 real          :: dustfraci(maxdusttypes),dustfracisum
 real, save    :: Mtot_in,Mgas_in,Mdust1_in,Mdust2_in
 integer       :: i,itype,lu
 logical, save :: init = .false.

 Mtot   = 0.
 Mgas   = 0.
 Mdust1 = 0. !--one-fluid dust mass
 Mdust2 = 0. !--two-fluid dust mass
 Macc   = 0.
 pmassi = massoftype(igas)
 dustfraci(:) = 0.
 do i=1,npart
    if (maxphase==maxp) then
       itype = iamtype(iphase(i))
       pmassi = massoftype(itype)
    endif
    Mtot = Mtot + pmassi
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (itype==1) then
          if (use_dustfrac) dustfraci(:) = dustfrac(:,i)
          dustfracisum = sum(dustfraci)
          Mgas   = Mgas   + pmassi*(1. - dustfracisum)
          Mdust1 = Mdust1 + pmassi*dustfracisum
       elseif (itype==2) then
          Mdust2 = Mdust2 + pmassi
       endif
    else
       Macc  = Macc + pmassi
    endif
 enddo

 if (.not.init) then
    open(newunit=lu,file='dustmass.ev',status='replace')
    write(lu,"('# ',5('[',i2.2,1x,a12,']',1x))") 1,'time',2,'Mtot',3,'Mgas',4,'Mdust1',5,'Mdust2',6,'Macc'
    init = .true.
    Mgas_in  = Mgas
    Mdust1_in = Mdust1
    Mdust2_in = Mdust2
    Mtot_in  = Mtot
 else
    open(newunit=lu,file='dustmass.ev',status='old',position='append')
 endif

 print*,' Mtot = ',Mtot,' Mgas = ',Mgas,' Mdust1 = ',Mdust1,' Mdust2 = ',Mdust2,' Macc = ',Macc
 print "(4(/,a,2pf6.2,'%'))",' dMtot  = ',(Mtot-Mtot_in)/Mtot_in,&
       ' dMgas  = ',(Mgas-Mgas_in)/Mtot_in,&
       ' dMdust1 = ',(Mdust1-Mdust1_in)/Mtot_in,&
       ' dMdust2 = ',(Mdust2-Mdust2_in)/Mtot_in,&
       ' dMacc  = ',Macc/Mtot_in
 write(lu,*) time,Mtot,Mgas,Mdust1,Mdust2,Macc
 close(lu)
 !print "(a,es10.3,e

end subroutine do_analysis

end module
