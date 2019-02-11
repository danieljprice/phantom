!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: mf_write
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: infile_utils, io, part, physcon
!+
!--------------------------------------------------------------------------
module mf_write

 implicit none

 public mflow_init, mflow_write,vmflow_init,vmflow_write,binpos_init,binpos_write,flow_analysis,&
        createbins,nradi,ncolsi

 !--If need more radial resolution increase nradi
 !--resolution depends only from number of particles

 !--------------
 integer :: nradi=300
 integer :: ncolsi=304 !ncolsi nradi+4
 real    :: gasrad=2 !change to decide gpos radius
 !------------

 character (len=120) :: mflowname
 real             :: maxradius
 integer, private :: imflow,ivmflow,ibinpos,ipartatrad
 character (len=120), private ::evfileprivate

 private

contains
subroutine mflow_write(time,dt)
 use part,             only:npart,xyzh,massoftype
!  use analysis,         only:read_discparams
 use io,               only:fatal
 real, intent(in)                :: time,dt
 real                            :: mass(nradi),ri,rad(nradi),dr,pmass,totalmassrad,totalmass

 integer                         :: i,k
 character (len=120)             :: num, formatout



 pmass=massoftype(1)
 mass(:)=0
 rad(:)=0
 totalmass=0

 dr = (maxradius)/real(nradi-1)
 do i=1,nradi
    rad(i)= real(i-1)*dr
 enddo



 do i = 1,npart

    if (xyzh(4,i) < -tiny(xyzh)) then !ACCRETED
       !!!Mass of accreted put in radius 0
       mass(1)=mass(1)+pmass

    elseif (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri=sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       do k=2,nradi
          if (ri<rad(k)) mass(k)=mass(k)+pmass
       enddo
    endif

    totalmass=totalmass+pmass
 enddo

 totalmassrad=mass(nradi)
 write(num,*) nradi+4
 formatout="("//trim(num)//"(es18.10,1x))"

 !open(unit=47,file=trim(mflowname))

 !---Correct the missing particle when accreted
 do i=2, nradi
    mass(i) = mass(i)+mass(1)
 enddo

 write(imflow,formatout) time,dt,totalmass,totalmassrad,mass

 call flush(imflow)

end subroutine mflow_write


subroutine mflow_init(iflow,evfile,infile)

!  use analysis, only:read_discparams
 use io,       only:fatal
 use part,     only:massoftype

 integer,           intent(in)   :: iflow
 character (len=*), intent(in)   :: evfile,infile
 integer                         :: idot,i,iline
 real, dimension(nradi)          :: rad
 real                            :: dr
 real                            :: R_in,R_out,p_index,q_index,M_star,H_R
 integer                         :: iparams=10,ierr
 character (len=120)             :: num,formatout,discprefix

 evfileprivate=evfile
 imflow=iflow
 idot = index(evfile,'.ev') - 1
 mflowname = evfile(1:idot)//'.mf'  !create .mf
 iline = index(infile,'.')
 discprefix = infile(1:iline-1)


 call read_discparams(trim(discprefix)//'.discparams',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
 if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')

 maxradius=R_out





 open(unit=imflow,file=mflowname,form='formatted',status='replace')

 !---Write header
 dr = (maxradius)/real(nradi-1)
 do i=1,nradi
    rad(i)= real(i-1)*dr
 enddo


 write(num,*) nradi
 formatout="("//trim(num)//"(es18.10,1x))"




 write(imflow,*)"  Mass at radius 0 is accreted mass"
 write(imflow,"('#',(1x,a17),I18)")" Number of radii:",nradi
 write(imflow,"('#',(1x,a17),2(es18.10,1X))")" Maxradius,pmass:", maxradius,massoftype(1)
 write(imflow,"('#',(1x,a17))")"radii:"
 write(imflow,formatout)rad

 write(imflow,"('#',4(1x,'[',i2.2,1x,a11,']',2x),(1x,'[',a4,1x,a9,']',2x))") &
       1,'time', &
       2,'dt', &
       3,'totmass', &
       4,'totmlastr', &
       "5-104","massinrad"

end subroutine mflow_init



subroutine vmflow_write(time,dt)
 use part,             only:npart,xyzh,vxyzu,massoftype,rhoh
 !  use analysis,         only:flow_analysis !contained in analysis_disc_MFlow
 use io,               only:fatal
 use physcon,          only:pi

 real, intent(in)                :: time,dt
 real                            :: flow(nradi),rad(nradi),dr,pmass,totalmass,rmax,rmin
 character (len=120)             :: num, formatout



 pmass=massoftype(1)
 rmax=maxradius
 rmin=0.0
 flow(:)=0
 rad(:)=0
 totalmass=0

 call createbins(rad,nradi,rmax,rmin,dr)
 call flow_analysis(xyzh,vxyzu,pmass,flow,npart,rad,nradi,dr)


 write(num,*) nradi+2
 formatout="("//trim(num)//"(es18.10,1x))"
 write(ivmflow,formatout) time,dt,flow
 call flush(ivmflow)

end subroutine vmflow_write


subroutine vmflow_init(ivflow,evfile,infile)

 ! use analysis, only:read_discparams
 use io,       only:fatal

 integer,           intent(in)   :: ivflow
 character (len=*), intent(in)   :: evfile,infile
 integer                         :: idot,i,iline
 real, dimension(nradi)          :: rad
 real                            :: dr
 real                            :: R_in,R_out,p_index,q_index,M_star,H_R,rmax,rmin
 integer                         :: iparams=10,ierr,numbercol(nradi+2)
 character (len=120)             :: num,formatout,formatcol,discprefix

 evfileprivate=evfile
 ivmflow=ivflow
 idot = index(evfile,'.ev') - 1
 mflowname = evfile(1:idot)//'_v.mflowv'  !create .vmf
 iline = index(infile,'.')
 discprefix = infile(1:iline-1)


 call read_discparams(trim(discprefix)//'.discparams',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)

 if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')

 maxradius=R_out
 rmax=maxradius
 rmin=0.0

 open(unit=ivmflow,file=mflowname,form='formatted',status='replace')

 call createbins(rad,nradi,rmax,rmin,dr)

 do i=1,(nradi+2)
    numbercol(i)=i
 enddo

 !--write format
 write(num,*) nradi
 formatout="(2(2x,'[',i2.2,a4,']',9x),"//trim(num)//"(es18.10,1x))"
 write(num,*) (nradi+2)
 formatcol="('#',"//trim(num)//"(I18,1x))"

!    rad=rad+dr/2

 !---Write header
 write(ivmflow,formatcol)numbercol
 write(ivmflow,formatout) &
       1,'time', &
       2,'xxdt', &
       rad

end subroutine vmflow_init

!----------------------------------------------------

subroutine binpos_write(time,dt)

 use part,               only:xyzmh_ptmass,npart,xyzh
 use physcon,            only:pi
 use io,                 only:igpos

 real, intent(in)                :: time,dt
 real                            :: ri
 integer                         :: i,idot
 character(len=20)               :: num


 write(ibinpos,'(7(es18.10,1x))')time, xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1), &
                                       xyzmh_ptmass(1,2),xyzmh_ptmass(2,2),xyzmh_ptmass(3,2)

 call flush(ibinpos)

!-----follow one gas particle

 if (xyzh(4,ipartatrad)  >  tiny(xyzh)) then ! IF ACTIVE
    write(igpos,'(4(es18.10,1x))')time, xyzh(1,ipartatrad),xyzh(2,ipartatrad),xyzh(3,ipartatrad)

 elseif(xyzh(4,ipartatrad) < -tiny(xyzh)) then !ACCRETED
    close(unit=igpos)
    !---Look for new particle to follow and create new file
    !--choose the gas particle at radius=gasrad

    do i=1,npart
       ri=sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       if(ri<(gasrad+0.1) .and. ri>(gasrad-0.1))  then
          ipartatrad=i
          exit
       elseif(i==npart)  then
          ipartatrad=npart
       endif
    enddo

    write(num,'(I8.8)') ipartatrad

    idot = index(evfileprivate,'.ev') - 1
    mflowname = evfileprivate(1:idot)//'-'//trim(num)//'.gpos'

    open(unit=igpos,file=mflowname,form='formatted',status='replace')

    write(igpos,"('#',4(1x,'[',i2.2,1x,a11,']',2x),(1x,'[',a4,1x,a9,']',2x))") &
        1,'time', &
        2,'xg', &
        3,'yg', &
        4,'zg'

    write(igpos,'(4(es18.10,1x))')time, xyzh(1,ipartatrad),xyzh(2,ipartatrad),xyzh(3,ipartatrad)

 endif

 call flush(igpos)

end subroutine binpos_write


subroutine binpos_init(ibinposi,evfile)

!  use analysis, only:read_discparams
 use io,       only:fatal,igpos
 use part,     only:npart,xyzh

 integer,           intent(in)   :: ibinposi
 character (len=*), intent(in)   :: evfile
 integer                         :: i,idot
 real                            :: ri
 character (len=120)             :: num

 evfileprivate=evfile
 ibinpos=ibinposi
 idot = index(evfile,'.ev') - 1
 mflowname = evfile(1:idot)//'.binpos'


 open(unit=ibinpos,file=mflowname,form='formatted',status='replace')


 write(ibinpos,"('#',7(1x,'[',i2.2,1x,a11,']',2x),(1x,'[',a4,1x,a9,']',2x))") &
       1,'time', &
       2,'x1', &
       3,'y1', &
       4,'z1', &
       5,'x2', &
       6,'y2', &
       7,'z2'

 !--choose the gas particle at radius=gasrad
 do i=1,npart

    ri=sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))

    if(ri<(gasrad+0.1) .and. ri>(gasrad-0.1))  then
       ipartatrad=i
       exit
    elseif(i==npart)  then
       ipartatrad=npart
    endif

 enddo

 write(num,'(I8.8)') ipartatrad

 mflowname = evfile(1:idot)//'-'//trim(num)//'.gpos'

 open(unit=igpos,file=mflowname,form='formatted',status='replace')

 write(igpos,"('#',4(1x,'[',i2.2,1x,a11,']',2x),(1x,'[',a4,1x,a9,']',2x))") &
       1,'time', &
       2,'xg', &
       3,'yg', &
       4,'zg'

end subroutine binpos_init

subroutine createbins(rad,nr,rmax,rmin,dr)
 use io, only: fatal

 real,    intent(inout)   :: dr
 real,    intent(in)      :: rmax,rmin
 real,    intent(inout)   :: rad(:)
 integer, intent(in)      :: nr
 integer                  :: i

 if(size(rad)<nr) call fatal('subroutine createbin','size(rad)<nr')

 dr = (rmax-rmin)/real(nr-1)
 do i=1,nr
    rad(i)=rmin + real(i-1)*dr
 enddo


end subroutine createbins

subroutine flow_analysis(xyzh,vxyz,pmass,flow,npart,rad,nr,dr)

 real, intent(inout) :: flow(:)
 real, intent(in)    :: xyzh(:,:),vxyz(:,:),rad(:),dr,pmass
 integer, intent(in) :: npart,nr
 real                :: rcili,vri
 integer             :: i,ii

 flow(:)=0

 do i=1,npart
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       rcili = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       ii = int((rcili-rad(1))/dr + 1)!ii=1 if rcili>rad(1)

       if (ii > nr) cycle
       if (ii < 1)  cycle



       !compute mass flow trough rad(ii)
       rcili=sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       vri=vxyz(1,i)*xyzh(1,i)/rcili+vxyz(2,i)*xyzh(2,i)/rcili
       flow(ii)=flow(ii)-vri*pmass/dr !normalized on dr
    endif

 enddo

end subroutine flow_analysis

subroutine read_discparams(filename,R_in,R_out,H_R,p_index,q_index,M_star,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,H_R,p_index,q_index,M_star
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)

! Read in parameters from the file discparams.list
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(R_in,'R_in',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_out,'R_out',db,ierr)
 if (ierr /= 0) return
 call read_inopt(H_R,'H_R',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_discparams



end module mf_write
