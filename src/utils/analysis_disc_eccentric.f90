!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for discs based on semimajor axis of particles
!  for studying eccentric orbits
!
!  REFERENCES: None
!
!  OWNER: Enrico Ragusa
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: infile_utils, io, part, physcon, sortutils
!+
!--------------------------------------------------------------------------

!--------------------------- N.B. ------------------------------ !
! discfrac is not the density, to compute the density one should !
! compute the surface of an eccentric semimajor axis bin.        !
! One should also keep in mind that the eccentricity equations   !
! used here are for purely keplerian orbits. Pressure            !
! corrections to the velocity introduce a spurious eccentricity  !
! of e_spur=beta*(H/R)^2. beta=1.5+p+q                           !
!--------------------------------------------------------------- !

module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'eccentric'
 public :: do_analysis,nr,createbins,read_discparams

 integer, parameter :: nr = 400
 !real, parameter :: rmin = 0,  rmax = 15


 interface read_discparams
  module procedure read_discparams, read_discparams2
 end interface

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use part,    only:rhoh,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use sortutils, only:indexx

 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(:,:),vxyz(:,:)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile


 character(len=9) :: output
 character(len=20) :: filename
 integer :: i,ii,ierr
 real :: R_in,R_out,H_R,p_index,q_index,M_star,Sig0,timesRout
 real :: G,dr,angx,angy,angz,ri,vri,vphi,phi,ai,phasei,ecci,Entoti,vcil2i,ji
 real :: Li(3), xx(3),vv(3),CMpos(3),CMVel(3)
 real :: a(nr),Lx(nr),Ly(nr),Lz(nr),h_smooth(nr),discfrac(nr),ecc(nr),phase(nr),RR(nr)

 real :: Hperc(nr),mass(nr),honH
 real,   allocatable ::z(:,:)
 integer,  allocatable ::indexz(:,:)
 integer :: j(nr),ninbin(nr),gausslimit(nr)
 integer :: k,l

 integer :: idot
 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 logical :: comment= .true.

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
 if (ierr /= 0) call fatal('analysis','could not open/read .discparams file')


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

 timesRout=2.
 call createbins(a,nr,R_out*timesRout,R_in,dr)

! Initialise arrays to zero
 ninbin(:)=1
 lx(:)=0.0
 ly(:)=0.0
 lz(:)=0.0

 ecc(:)=0.0
 phase(:)=0.0


 h_smooth(:)=0.0
 discfrac(:)=0.0
 Hperc(:)=0.0
 mass(:)=0.0
 RR(:)=0.0
 xx(:)=0.0
 vv(:)=0.0

 CMpos=(xyzmh_ptmass(4,1)*xyzmh_ptmass(1:3,1)+xyzmh_ptmass(4,2)*xyzmh_ptmass(1:3,2))/(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2))
 CMVel=(xyzmh_ptmass(4,1)*vxyz_ptmass(1:3,1)+xyzmh_ptmass(4,2)*vxyz_ptmass(1:3,2))/(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2))

 do i = 1,npart


    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE

       !change coordinate system
       xx=xyzh(1:3,i)-CMpos
       vv=vxyz(1:3,i)-CMVel

       ri = sqrt(dot_product(xx(1:2),xx(1:2)))
       phi= atan2(xx(2),xx(1))
       vcil2i=dot_product(vv(1:2),vv(1:2))
       vri=vv(1)*xx(1)/ri+vv(2)*xx(2)/ri
       vphi=-vv(1)*xx(2)/ri+vv(2)*xx(1)/ri
       ji=xx(1)*vv(2)-xx(2)*vv(1)

       Entoti=1./2.*vcil2i-M_star/ri
       ai=-M_star/(2.*Entoti)
       phasei=atan2(-(ji/M_star*vv(1)+xx(2)/ri),(ji/M_star*vv(2)-xx(1)/ri))
       ecci=sqrt(1.-ji**2/(M_star*ai))
       ii = int((ai-a(1))/dr + 1)


       if (ii > nr) cycle
       if (ii < 1)  cycle

       mass(ii)=mass(ii)+pmass

       Li(1) = pmass*(xx(2)*vv(3)-xx(3)*vv(2))
       Li(2) = pmass*(xx(3)*vv(1)-xx(1)*vv(3))
       Li(3) = pmass*(xx(1)*vv(2)-xx(2)*vv(1))

       Lx(ii)=Lx(ii)+Li(1)
       Ly(ii)=Ly(ii)+Li(2)
       Lz(ii)=Lz(ii)+Li(3)

       !Eccentricity related quantities
       ecc(ii)=ecc(ii)+ecci
       phase(ii)=phase(ii)+phasei
       RR(ii)=RR(ii)+ri


       h_smooth(ii) = h_smooth(ii) + xyzh(4,i)

       ninbin(ii) = ninbin(ii) + 1


    elseif (xyzh(4,i) < -tiny(xyzh)) then !ACCRETED
       angx = angx + pmass*(xx(2)*vv(3) - xx(3)*vv(2))
       angy = angy + pmass*(xx(3)*vv(1) - xx(1)*vv(3))
       angz = angz + pmass*(xx(1)*vv(2) - xx(2)*vv(1))
    endif
 enddo

!-------------------------------------
! Compute disc thickness using std dev
!-------------------------------------

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
       !change coordinate system
       xx=xyzh(1:3,i)-CMpos
       vv=vxyz(1:3,i)-CMVel

       ri = sqrt(dot_product(xx(1:2),xx(1:2)))
       vcil2i=dot_product(vv(1:2),vv(1:2))
       Entoti=1./2.*vcil2i-M_star/ri
       ai=-M_star/(2*Entoti)
       ii = int((ai-a(1))/dr + 1)

       if (ii > nr) cycle
       if (ii < 1)  cycle

       z(ii,j(ii))= abs(xx(3))
       j(ii)=j(ii)+1
       if(j(ii)>ninbin(ii)+1)then
          print*, 'out of array limit (ii,j(ii),ninbin(ii):)',ii,j(ii)-1,ninbin(ii)
       endif
    endif
 enddo


 do i=1,nr
    call  indexx(ninbin(i),z(i,:),indexz(i,:))
 enddo

 gausslimit=int(0.68*ninbin)


!----index at which i found sorted z below which 68% of particles (1 gaussian st dev)

 do i=1,nr
    if (gausslimit(i) < 1)  cycle
    Hperc(i)=z(i,indexz(i,gausslimit(i)))
 enddo

 !Average quantities on particle number into the bin
 discfrac(:)=real(ninbin(:))/npart
 ecc(:)=ecc(:)/ninbin(:)
 phase(:)=phase(:)/ninbin(:)
 h_smooth(:)=h_smooth(:)/ninbin(:)
 RR(:)=RR(:)/ninbin(:)


!-------------------------
! writing in analysis file
!-------------------------

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',10(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'a', &
       2,'discfrac', &
       3,'<h>/H', &
       4,'Lx', &
       5,'Ly', &
       6,'Lz', &
       7,'H', &
       8,'ecc',&
       9,'phase',&
       10,'R'

 do i=1,nr
    !if H=0 does not divide
    if(.not. Hperc(i)==0.) then
       honH=h_smooth(i)/Hperc(i)
    else
       honH=0.
    endif

    write(iunit,'(10(es18.10,1X))') a(i),discfrac(i),honH,Lx(i),Ly(i),Lz(i),&
                                                Hperc(i),ecc(i),phase(i)/acos(-1.)*180.,RR(i)
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
 call read_inopt(H_R,'H/R_in',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call read_inopt(Sig0,'sig_in',db,ierr)
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





end module analysis
