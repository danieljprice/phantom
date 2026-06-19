!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Output the maximum quantity and its position in a file
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters: None
!
! :Dependencies: growth, infile_utils, io, part, physcon, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'disc'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use part,    only:xyzmh_ptmass,dustfrac,rhoh!dustprop
 use infile_utils, only:open_db_from_file,read_inopt,close_db,inopts
 use units, only:unit_density
 use growth, only:get_size
 character(len=*), intent(in)    :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in)    :: npart,iunit,numfile
 real :: dens,maxquant,xmax,ymax,zmax
 integer :: imax,i,outunit
 character(len=25) :: filename

 ! select quantity
 maxquant = -1e99
 imax = 0
 do i=1,npart
    dens = rhoh(xyzh(4,i),pmass) * dustfrac(1,i) * unit_density !1-f
    !dens = get_size(dustprop(1,i),dustprop(2,i),1.0) * udist   !2-f
    if (dens>maxquant) then
       maxquant = dens
       imax = i
    endif
 enddo

 ! centering on sink
 xmax = xyzh(1,imax) - xyzmh_ptmass(1,1)
 ymax = xyzh(2,imax) - xyzmh_ptmass(2,1)
 zmax = xyzh(3,imax) - xyzmh_ptmass(3,1)

 outunit = 99
 !write(filename,'("maxgrainsize_",I0,".txt")') numfile
 write(filename,'("maxdustdenss_",I0,".txt")') numfile
 open(unit=outunit,file=filename,status='unknown',action='write')

 write(outunit,*) maxquant
 write(outunit,*) xmax
 write(outunit,*) ymax
 write(outunit,*) zmax

end subroutine do_analysis

end module analysis

