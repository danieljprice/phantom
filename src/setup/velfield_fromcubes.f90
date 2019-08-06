!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: velfield
!
!  DESCRIPTION:
! this module does setup of the velocity field
! from sphNG cubes
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module velfield
 public :: set_velfield_from_cubes

 private

contains

!----------------------------------------
!+
! Trilinear interpolation from the grid
!+
!----------------------------------------
pure real function vfield_interp(vgrid,iposx,iposy,iposz,delx,dely,delz,deli)
 integer,      intent(in) :: iposx,iposy,iposz
 real(kind=4), intent(in) :: vgrid(:,:,:)
 real,         intent(in) :: delx,dely,delz,deli
 real :: velx1,velx2,vely1,vely2

 velx1 = vgrid(iposx,iposy,iposz)   + delx/deli*(vgrid(iposx+1,iposy,iposz)-vgrid(iposx,iposy,iposz))
 velx2 = vgrid(iposx,iposy+1,iposz) + delx/deli*(vgrid(iposx+1,iposy+1,iposz)-vgrid(iposx,iposy+1,iposz))
 vely1 = velx1 + dely/deli*(velx2-velx1)

 velx1 = vgrid(iposx,iposy,iposz+1)   + delx/deli*(vgrid(iposx+1,iposy,iposz+1)-vgrid(iposx,iposy,iposz+1))
 velx2 = vgrid(iposx,iposy+1,iposz+1) + delx/deli*(vgrid(iposx+1,iposy+1,iposz+1)-vgrid(iposx,iposy+1,iposz+1))
 vely2 = velx1 + dely/deli*(velx2-velx1)

 vfield_interp = vely1 + delz/deli*(vely2-vely1)

end function vfield_interp

!----------------------------------------
!+
!  Read data cube from file
!+
!----------------------------------------
subroutine read_cube(filename,vgrid,nspace,iunit,ierr)
 character(len=*), intent(in)  :: filename
 real(kind=4),     intent(out) :: vgrid(:,:,:)
 integer,          intent(in)  :: iunit,nspace
 integer,          intent(out) :: ierr
 integer :: i,j,k

 open(unit=iunit,file=filename,form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'ERROR opening '//trim(filename)
    return
 endif
 read(iunit,iostat=ierr) (((vgrid(i,j,k),i=1,nspace),j=1,nspace),k=1,nspace)
 if (ierr > 0) then
    write(*,"(1x,a,i2,a)") 'ERROR reading '//trim(filename)//' on unit ',iunit,' (wrong endian?): aborting'
 elseif (ierr < 0) then
    write(*,*) 'ERROR reading '//trim(filename)//' (hit end of file/record): aborting'
 endif
 close(iunit)

end subroutine read_cube

!------------------------------------------------------------------------
!+
!  set up velocity field from cubes
!  rmax is the radius of sphere (use rmax = xmax for cartesian geometry)
!+
!------------------------------------------------------------------------
subroutine set_velfield_from_cubes(xyzh,vxyzu,npart,filevx,filevy,filevz,&
                                   amplitude,rmax,falloffnearedge,ierr)
 integer,          intent(in)  :: npart
 real,             intent(in)  :: xyzh(:,:)
 real,             intent(out) :: vxyzu(:,:)
 character(len=*), intent(in)  :: filevx,filevy,filevz
 logical,          intent(in)  :: falloffnearedge
 real,             intent(in)  :: amplitude,rmax
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 43

 integer :: i,nspace,nbytes_per_vel,nbytes_fheader
 integer :: ifilesizex,ifilesizey,ifilesizez
 integer :: ierrx,ierry,ierrz,iposx,iposy,iposz
 logical :: iexistx,iexisty,iexistz
 real(kind=4), allocatable :: velx(:,:,:),vely(:,:,:),velz(:,:,:)
 real :: xi,yi,zi,deli,delx,dely,delz,radnorm,radius,factor

 ierr = 0
 !
 !--set velocities to zero initially
 !
 vxyzu(1:3,1:npart) = 0.
 !
 !--check that the velocity files exist
 !
 inquire(file=filevx,exist=iexistx,size=ifilesizex)
 inquire(file=filevy,exist=iexisty,size=ifilesizey)
 inquire(file=filevz,exist=iexistz,size=ifilesizez)
 if (.not.iexistx) write(*,*) 'ERROR: '//trim(filevx)//' not found'
 if (.not.iexisty) write(*,*) 'ERROR: '//trim(filevy)//' not found'
 if (.not.iexistz) write(*,*) 'ERROR: '//trim(filevz)//' not found'
 if (.not.(iexistx .and. iexisty .and. iexistz)) then
    ierr = 1
    return
 endif
 if (.not.((ifilesizey==ifilesizex) .and. (ifilesizez==ifilesizex))) then
    write(*,*) 'ERROR: velocity files have different sizes (should all be the same): aborting...'
    write(*,*) 'File sizes in Bytes = ',ifilesizex,ifilesizey,ifilesizez
    ierr = 1
    return
 endif
!
!--work out the size of the velocity files from the filesize
!  and allocate memory for the velocity grid arrays
!
 nbytes_per_vel = kind(velx)
 nbytes_fheader = 8   ! Fortran inserts 4 bytes before and after a write statement
 nspace = nint(((ifilesizex - nbytes_fheader)/nbytes_per_vel)**(1./3.))
 write(*,"(1x,a,i4,a)") 'Size of velocity grid (from filesize) = ',nspace,'^3'
 !write (*,*) 'enter size of velocity files (e.g.n=32^3)'

 allocate(velx(nspace,nspace,nspace),stat=ierrx)
 allocate(vely(nspace,nspace,nspace),stat=ierry)
 allocate(velz(nspace,nspace,nspace),stat=ierrz)
 if (ierrx /= 0 .or. ierry /= 0 .or. ierrz /= 0) then
    write(*,*) 'ERROR allocating memory'
    call deallocate_vels
    ierr = 3
    return
 endif
 !
 !--read data cubes
 !
 call read_cube(filevx,velx,nspace,iunit,ierrx)
 call read_cube(filevy,vely,nspace,iunit,ierry)
 call read_cube(filevz,velz,nspace,iunit,ierrz)
 if (ierrx /= 0 .or. ierry /= 0 .or. ierrz /= 0) then
    ierr = 2
    call deallocate_vels
    return
 endif

 radnorm = rmax
 deli = radnorm/real(nspace/2)
 do i = 1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    iposx = int(xi/radnorm*(nspace/2)+(nspace/2)+0.5)
    iposx = min(max(iposx, 1),nspace-1)
    iposy = int(yi/radnorm*(nspace/2)+(nspace/2)+0.5)
    iposy = min(max(iposy, 1),nspace-1)
    iposz = int(zi/radnorm*(nspace/2)+(nspace/2)+0.5)
    iposz = min(max(iposz, 1),nspace-1)
    delx = xi - (iposx-(nspace/2)-0.5)/real(nspace/2)*radnorm
    dely = yi - (iposy-(nspace/2)-0.5)/real(nspace/2)*radnorm
    delz = zi - (iposz-(nspace/2)-0.5)/real(nspace/2)*radnorm
    radius = sqrt(xi**2 + yi**2 + zi**2)
!
!--Find interpolated velocities
!
    vxyzu(1,i) = amplitude*vfield_interp(velx,iposx,iposy,iposz,delx,dely,delz,deli)
    vxyzu(2,i) = amplitude*vfield_interp(vely,iposx,iposy,iposz,delx,dely,delz,deli)
    vxyzu(3,i) = amplitude*vfield_interp(velz,iposx,iposy,iposz,delx,dely,delz,deli)

    if (falloffnearedge) then
       if (radius/radnorm > 0.9) then
          factor = 0.0
       elseif (radius/radnorm < 0.8) then
          factor = 1.0
       else
          factor = (0.9-radius/radnorm)*10.0
       endif
       vxyzu(1,i) = vxyzu(1,i)*factor
       vxyzu(2,i) = vxyzu(2,i)*factor
       vxyzu(3,i) = vxyzu(3,i)*factor
    endif
 enddo

 call deallocate_vels
 return

contains

subroutine deallocate_vels

 if (allocated(velx)) deallocate(velx)
 if (allocated(vely)) deallocate(vely)
 if (allocated(velz)) deallocate(velz)

end subroutine deallocate_vels

end subroutine set_velfield_from_cubes

end module velfield
