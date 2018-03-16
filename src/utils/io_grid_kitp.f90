!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: io_grid
!
!  DESCRIPTION:
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, params
!+
!--------------------------------------------------------------------------

module io_grid
 implicit none
 integer, parameter, private :: lio=4

contains

!***********************************************************************
subroutine openr_grid(lun,file,mx,my,mz,mv,sx,sy,sz)
!
! OPENING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 use params
 character(len=mfile) file
 integer lun,mx,my,mz,mv
 real sx,sy,sz

 if (do_seq_read) then
    call openr_grid_seq(lun,file,mx,my,mz,mv,sx,sy,sz)
 else
    call openr_grid_dir(lun,file,mx,my,mz,mv,sx,sy,sz)
 endif
end subroutine openr_grid

!***********************************************************************
subroutine openw_grid(lun,file,mx,my,mz,mv,sx,sy,sz)
!
! OPENING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 use params
 character(len=mfile) file
 integer lun,mx,my,mz,mv
 real sx,sy,sz

 if (do_seq_write) then
    call openw_grid_seq(lun,file,mx,my,mz,mv,sx,sy,sz)
 else
    call openw_grid_dir(lun,file,mx,my,mz,mv,sx,sy,sz)
 endif
end subroutine openw_grid

!***********************************************************************
subroutine read_grid(lun,mx,my,mz,ivar,v)
!
! READING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 use params
 integer lun,ivar,mx,my,mz,iz
 real :: v(mx,my,mz)

 if (do_seq_read) then
    call read_grid_seq(lun,mx,my,mz,ivar,v)
 else
    call read_grid_dir(lun,mx,my,mz,ivar,v)
 endif
end subroutine read_grid

!***********************************************************************
subroutine write_grid(lun,mx,my,mz,ivar,v)
!
! READING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 use params
 integer lun,ivar,mx,my,mz,iz
 real :: v(mx,my,mz)

 if (do_seq_write) then
    call write_grid_seq(lun,mx,my,mz,ivar,v)
 else
    call write_grid_dir(lun,mx,my,mz,ivar,v)
 endif
end subroutine write_grid

!***********************************************************************
subroutine openr_grid_dir(lun,file,mx,my,mz,mv,sx,sy,sz)
!
! OPENING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 use params
 use io
 character(len=mfile) file
 integer lun,mx,my,mz,mv
 real sx,sy,sz
 namelist /dim/mx,my,mz,mv,sx,sy,sz,offset,x_index,do_exp,do_massflux

 print *,'reading ',file(1:index(file,'.',back=.true.))//'dim'
 open(lun,file=file(1:index(file,'.',back=.true.))//'dim',form='formatted',status='old')
 sx=1.; sy=sx; sz=sx; mv=4; offset=-.5
 do_exp=.false.; do_massflux=.false.
 x_index=(/2,5/)
 read(lun,dim)
 write(*,dim)
 close(lun)

 print *,'reading direct access file: ',trim(file)
 open(lun,file=file,form='unformatted',status='old', &
        access='direct',recl=mx*my*lio)
end subroutine openr_grid_dir

!***********************************************************************
subroutine read_grid_dir(lun,mx,my,mz,ivar,v)
!
! READING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 integer lun,ivar,mx,my,mz,iz
 real :: v(mx,my,mz)

 print '(i3,1x,a)',ivar,'reading direct access record'
 do iz=1,mz
    if (mz >= 1000 .and. mod(iz,100)==0) print'(i7)',iz
    read(lun,rec=iz+(ivar-1)*mz) v(:,:,iz)
 enddo
end subroutine read_grid_dir

!***********************************************************************
subroutine openw_grid_dir(lun,file,mx,my,mz,mv,sx,sy,sz)
!
! OPENING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 use params
 use io
 character(len=mfile) file
 integer lun,mx,my,mz,mv
 real sx,sy,sz
 namelist /dim/mx,my,mz,mv,sx,sy,sz,offset,x_index,do_exp,do_massflux

 print *,'writing ',file(1:index(file,'.',back=.true.))//'dim'
 open(lun,file=file(1:index(file,'.',back=.true.))//'dim',form='formatted',status='unknown')
 write(lun,dim)
 write(*,dim)
 close(lun)

 print *,'opening direct access file '//trim(file)//' for writing'
 open(lun,file=file,form='unformatted',status='unknown', &
        access='direct',recl=mx*my*lio)
end subroutine openw_grid_dir

!***********************************************************************
subroutine write_grid_dir(lun,mx,my,mz,ivar,v)
!
! WRITING DIRECT ACCESS GRID DATA
!
!-----------------------------------------------------------------------
 integer lun,ivar,mx,my,mz,iz
 real :: v(mx,my,mz)

 print '(i3,1x,a)',ivar,'writing direct access record'
 do iz=1,mz
    write(lun,rec=iz+(ivar-1)*mz) v(:,:,iz)
 enddo
end subroutine write_grid_dir

!***********************************************************************
subroutine openr_grid_seq(lun,file,mx,my,mz,mv,sx,sy,sz)
!
! OPEN SEQUENTIAL GRID DATA FILE
!
!-----------------------------------------------------------------------
 use params
 character(len=mfile) file
 integer lun,mx,my,mz,mv
 real sx,sy,sz

 write(*,*) 'reading from sequential grid data file ',trim(file)
 open(lun,file=file,form='unformatted',status='old')
 read(lun) mx,my,mz,mv
 read(lun) sx,sy,sz
 write(*,*) 'mx,my,mz,mv =',mx,my,mz,mv
 write(*,*) 'sx,sy,sz =',sx,sy,sz
end subroutine openr_grid_seq

!***********************************************************************
subroutine read_grid_seq(lun,mx,my,mz,ivar,v)
!
! READING SEQUENTIAL GRID DATA
!
!-----------------------------------------------------------------------
 integer lun,mx,my,mz,iz,ivar
 real :: v(mx,my,mz)

 print '(i3,1x,a)',ivar,' reading sequential access record'
 read(lun) v
end subroutine read_grid_seq

!***********************************************************************
subroutine openw_grid_seq(lun,file,mx,my,mz,mv,sx,sy,sz)
!
! OPEN SEQUENTIAL GRID DATA FILE
!
!-----------------------------------------------------------------------
 use params
 character(len=mfile) file
 integer lun,mx,my,mz,mv
 real sx,sy,sz

 write(*,*) 'opening sequential grid data file ',trim(file)
 open(lun,file=file,form='unformatted',status='unknown')
 write(lun) mx,my,mz,mv
 write(lun) sx,sy,sz
 write(*,*) 'mx,my,mz,mv =',mx,my,mz,mv
 write(*,*) 'sx,sy,sz =',sx,sy,sz
end subroutine openw_grid_seq

!***********************************************************************
subroutine write_grid_seq(lun,mx,my,mz,ivar,v)
!
! WRITING SEQUENTIAL GRID DATA
!
!-----------------------------------------------------------------------
 integer lun,mx,my,mz,iz,ivar
 real :: v(mx,my,mz)

 print '(i3,1x,a)',ivar,'writing sequential access record'
 write(lun) v
end subroutine write_grid_seq

end module io_grid
