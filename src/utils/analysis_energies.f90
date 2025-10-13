!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine computing the energy accounting for accreted
! particles
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: energies, evwrite, metric_tools, options, part
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'energies'
 public :: do_analysis

 logical, private :: first = .true.
 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use energies,     only:compute_energies,track_mass,ekin,emag,etherm,epot,etot,&
                        eacc,etotall,totmom,angtot,angall
 use metric_tools, only:init_metric
 use part,         only:metrics,metricderivs,gr
 use evwrite,      only:init_evfile
 use options,      only:iexternalforce
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 if (gr) then
    call init_metric(npart,xyzh,metrics,metricderivs)
    iexternalforce = 1
 endif
 if (first) then
    call init_evfile(1,'crap.ev',open_file=.false.)
 endif
 track_mass = .true.
 call compute_energies(time)

 if (first) then
    open(unit=1,file='energies.ev',status='new',action='write')
    write(1,"(a)") '# time,ekin,etherm,emag,epot,etot,eacc,etot+eacc,totmom,angtot,etotall,angall'
    first = .false.
 endif
 write(1,*) time,ekin,etherm,emag,epot,etot,eacc,etot+eacc,totmom,angtot,etotall,angall

 print*,'                  TOTAL ENERGY IS: ',etot
 print*,' TOTAL ENERGY INCLUDING ACCRETION: ',etotall

end subroutine do_analysis

end module analysis
