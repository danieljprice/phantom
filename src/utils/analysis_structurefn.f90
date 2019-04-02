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
!  Analysis routine to produce structure functions
!  from SPH particle data
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, io_structurefn, part, structurefn_part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'structure functions'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use boundary,         only:dxbound,dybound,dzbound
 use structurefn_part, only:get_structure_fn
 use io_structurefn,   only:mfile,power_unit,openw_sf,write_sf,write_structfiles,write_cfstruct
 use part,             only:rhoh
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                   :: i,ierr,iunitw
 character(len=mfile)      :: fileout,origin
 integer, parameter        :: norder = 10
 integer, parameter        :: nbins = 1024 !int(npart**(1./3.)) + 1
 integer(kind=8) :: ncount(nbins)
 real(kind=8) :: sf(2,norder,nbins)
 real :: distbins(nbins)
 real :: orders(norder)
 real, parameter           :: rho_power = 0.
 real                      :: distmin,distmax
 !--local memory
 real, allocatable :: rho(:)

 do i=1,norder
    orders(i) = i
 enddo
 iunitw = 54

! max distance in a periodic box is a diagonal across the box
 distmax = sqrt(0.25*dxbound**2 + 0.25*dybound**2 + 0.25*dzbound**2)
 distmin = distmax/real(nbins)
 print*,' Calculating structure functions using ',nbins,' bins between ',distmin,' and ',distmax
!
!--allocate memory for density array
!
 allocate(rho(npart),stat=ierr)
 if (ierr /= 0) then
    print*,' ERROR ALLOCATING MEMORY for rho'
    return
 endif
 !--get density from h
 do i=1,npart
    rho(i) = rhoh(xyzh(4,i),particlemass)
 enddo
!
!--calculate structure function
!
 call get_structure_fn(sf,nbins,norder,distmin,distmax,distbins,ncount, &
                       npart,xyzh,vxyzu,rho,dxbound,dybound,dzbound,.true.,ierr)
 deallocate(rho)
 if (ierr /= 0) return
!
!--write to Aake format .sfn file
!
 fileout = trim(dumpfile)//'.sfn'
 origin = ' '
 origin = trim(dumpfile)
 print "(a)",' writing to '//trim(fileout)
 call openw_sf(fileout,origin,nbins,distbins,norder,1)
 call write_sf(nbins,norder,sf,orders,rho_power)
 close(power_unit)
!
!--write to plain ascii format .struct files
!
 call write_structfiles(trim(dumpfile),nbins,distbins,norder,sf,rho_power,time,iunitw)
!
!--write to Christoph Federrath format ascii files
!
 call write_cfstruct(trim(dumpfile),nbins,distbins,norder,ncount,sf,rho_power,time,iunitw)

end subroutine do_analysis

end module
