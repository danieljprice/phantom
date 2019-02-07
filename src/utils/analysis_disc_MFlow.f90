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
!  Analysis routine for discs by MFlow
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: infile_utils, io, part, physcon, sortutils
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'MFlow'
 public :: do_analysis,nr,rmin,rmax,createbins,flow_analysis,read_discparams

 integer, parameter :: nr = 400
 real, parameter :: rmin = 0,  rmax = 15


 interface read_discparams
  module procedure read_discparams, read_discparams2
 end interface

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use part,    only:rhoh
 use sortutils, only:indexx

 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(:,:),vxyz(:,:)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile


 character(len=9) :: output
 character(len=20) :: filename
 integer :: i,ii,ierr
 real :: R_in,R_out,H_R,p_index,q_index,M_star,Sig0
 real :: G,dr,cs0,angx,angy,angz,ri,area,vri,vphi,phi
 real :: angtot,Ltot
 real :: Li(3)
 real :: rad(nr),Lx(nr),Ly(nr),Lz(nr),h_smooth(nr),sigma(nr),cs(nr),H(nr),omega(nr),e1(nr)
 real :: unitlx(nr),unitly(nr),unitlz(nr),sigmavrsini(nr),sigmavphi(nr),sigmavrcosi(nr)

 real :: flow(nr),Hperc(nr),mass(nr),matradi(nr),buf
 real,   allocatable ::z(:,:)
 integer,  allocatable ::indexz(:,:)
 integer :: j(nr),ninbin(nr),gausslimit(nr)
 integer :: k,l

 integer :: idot
 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 logical :: comment= .false.

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile
 idot = index(dumpfile,'_') - 1
 filename = dumpfile(1:idot)  !create filename



 write(output,"(a4,i5.5)") 'angm',numfile
 write(*,'("Output file name is ",A)') output

! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G==1")')
 G = 1.0

 call read_discparams(''//trim(filename)//'.discparams',R_in,R_out,H_R,p_index,q_index,M_star,Sig0,iparams,ierr)
 if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')


! Print out the parameters
 write(*,*)
 write(*,'("Parameters are:")')
 write(*,*) 'R_in    = ',R_in
 write(*,*) 'R_out   = ',R_out
 write(*,*) 'H_R     = ',H_R
 write(*,*) 'p_index = ',p_index
 write(*,*) 'q_index = ',q_index
 write(*,*) 'M_star  = ',M_star
 write(*,*)
 write(*,*)


 call createbins(rad,nr,rmax,rmin,dr)



! Initialise arrays to zero
 ninbin(:)=0
 lx(:)=0.0
 ly(:)=0.0
 lz(:)=0.0
 sigmavrcosi(:)=0.0
 sigmavrsini(:)=0.0
 sigmavphi(:)=0.0


 h_smooth(:)=0.0
 sigma(:)=0.0
 Hperc(:)=0.0
 mass(:)=0.0


! Set up cs0: cs = cs0 * R^-q
 cs0 = H_R * sqrt(G*M_star) * R_in**(q_index-0.5)

! And thus the sound speed array
 do i=1,nr
    cs(i) = cs0 * rad(i)**(-q_index)
    omega(i) = sqrt(G*M_star/rad(i)**3)
 enddo


! and thus the disc scale height
 do i=1,nr
    H(i) = cs(i)/omega(i)
 enddo

 angx = 0.0
 angy = 0.0
 angz = 0.0

! Loop over particles putting properties into the correct bin
 do i = 1,npart


    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       phi= atan(xyzh(2,i)/xyzh(1,i))
       vri=vxyz(1,i)*xyzh(1,i)/ri+vxyz(2,i)*xyzh(2,i)/ri
       vphi=-vxyz(1,i)*xyzh(2,i)/ri+vxyz(2,i)*xyzh(1,i)/ri
       ii = int((ri-rad(1))/dr + 1)
!       print*,vri,vphi

       if (ii > nr) cycle
       if (ii < 1)  cycle

       mass(ii)=mass(ii)+pmass
       if(ii==1)then
          area=pi*(dr/2)**2
       else
          area = (pi*((rad(ii)+dr/2.)**2-(rad(ii)- dr/2.)**2))
       endif
       sigma(ii) = sigma(ii) + pmass/area
       sigmavrcosi(ii)=sigmavrcosi(ii)+vri*cos(phi)!/xyzh(4,i)**3
       sigmavrsini(ii)=sigmavrsini(ii)+vri*sin(phi)!/xyzh(4,i)**3
       sigmavphi(ii)=sigmavphi(ii)+vphi!/xyzh(4,i)**3

       Li(1) = pmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
       Li(2) = pmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
       Li(3) = pmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))

       Lx(ii)=Lx(ii)+Li(1)
       Ly(ii)=Ly(ii)+Li(2)
       Lz(ii)=Lz(ii)+Li(3)

       h_smooth(ii) = h_smooth(ii) + xyzh(4,i)

       ninbin(ii) = ninbin(ii) + 1

!       print*, 'stai sforando ciccio11!',ii,j(ii),ninbin(ii)




    elseif (xyzh(4,i) < -tiny(xyzh)) then !ACCRETED
       angx = angx + pmass*(xyzh(2,i)*vxyz(3,i) - xyzh(3,i)*vxyz(2,i))
       angy = angy + pmass*(xyzh(3,i)*vxyz(1,i) - xyzh(1,i)*vxyz(3,i))
       angz = angz + pmass*(xyzh(1,i)*vxyz(2,i) - xyzh(2,i)*vxyz(1,i))
    endif
 enddo

 do i=1,nr
    if(sigmavphi(i)==0)then
       e1(i)=0
    else
       e1(i)=sqrt(sigmavrcosi(i)**2+sigmavrsini(i)**2)/abs(sigmavphi(i))
    endif
 enddo
 !e1(:)=0
 sigma=sigma/Sig0
 call flow_analysis(xyzh,vxyz,pmass,flow,npart,rad,nr,dr)

!print*,sigmavrcosi

 l=maxval(ninbin)


 allocate(z(nr,l))
 allocate(indexz(nr,l))
 indexz(:,:)=0
 z(:,:)=0

 do i=1,nr

    do k=1,ninbin(i)
       indexz(i,k)=k
    enddo

 enddo

 j(:)=1

 do i=1,npart
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       ii = int((ri-rad(1))/dr + 1)

       if (ii > nr) cycle
       if (ii < 1)  cycle

       z(ii,j(ii))= abs(xyzh(3,i))
       j(ii)=j(ii)+1
       if(j(ii)>ninbin(ii)+1)then
          print*, 'out of array limit (ii,j(ii),ninbin(ii):)',ii,j(ii),ninbin(ii)
       endif


    endif

 enddo


 do i=1,nr
    call  indexx(ninbin(i),z(i,:),indexz(i,:))
 enddo

 gausslimit=int(0.68*ninbin)


!----index at which i found sorted z below which 68% of particles, 1 gaussian sigma

 do i=1,nr
    if (gausslimit(i) < 1)  cycle
    Hperc(i)=z(i,indexz(i,gausslimit(i)))
 enddo

!-----How much mass contained inside a certain radius?

 buf=0.
 do i=1,nr
    buf=buf+mass(i)
    matradi(i)=buf
 enddo



! Print angular momentum of accreted particles
 angtot = sqrt(angx*angx + angy*angy + angz*angz)
 print*,' angular momentum of accreted particles = ',angtot

! Convert total angular momentum into a unit vector, and average h_smooth
 do i = 1,nr
    Ltot = sqrt(Lx(i)*Lx(i) + Ly(i)*Ly(i) + Lz(i)*Lz(i))

    unitlx(i) = Lx(i)/Ltot
    unitly(i) = Ly(i)/Ltot
    unitlz(i) = Lz(i)/Ltot

    if (ninbin(i) > 0) h_smooth(i) = h_smooth(i)/ninbin(i)
 enddo

! Now loop over rings to calculate required quantities
 do i = 1, nr
    if(ninbin(i)<2) then
       unitlx(i)=0.0
       unitly(i)=0.0
       unitlz(i)=0.0
       sigma(i)=0.0
       h_smooth(i) = 0.0
    else
       h_smooth(i) = h_smooth(i)/Hperc(i)
    endif
 enddo

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'radius', &
       2,'sigma', &
       3,'<h>/H', &
       4,'lx', &
       5,'ly', &
       6,'lz', &
       7,'Hperc', &
       8,'flow', &
       9,'e1'

 do i=1,nr

    write(iunit,'(9(es18.10,1X))') rad(i),sigma(i),h_smooth(i),unitlx(i),unitly(i),unitlz(i),Hperc(i),flow(i),e1(i)

 enddo

 close(iunit)

 print*,"Number of particles in each bin:"

 if(comment)then
    do i=1, nr
       print*,"i, ninbin(i),Hperc(i):",i,ninbin(i),Hperc(i)
    enddo
 endif
end subroutine do_analysis


!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------

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

subroutine read_discparams2(filename,R_in,R_out,H_R,p_index,q_index,M_star,Sig0,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,H_R,p_index,q_index,M_star,Sig0
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
 call read_inopt(Sig0,'Sig0',db,ierr)
 if (ierr /= 0) return


 call close_db(db)

end subroutine read_discparams2




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


end module analysis
