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
!  Analysis calculating magnetic field and wave characteristics.
!  Currently only possible when using fixed, non-calculated coefficients
!  for ambipolar diffusion and/or the Hall effect
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
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'bzrms'
 public :: do_analysis

 real,    private   :: rhoin0,Bxin0,amplitude,kwave
 real,    private   :: C_HE,C_AD
 real,    private   :: Bzrmsc,Bzrmsav
 real,    parameter :: meanmolmass = 2.3809523809523809
 logical, private   :: ambitest,halltest
 logical, private   :: firstcall = .true.


 private

contains
!-----------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,      only: id,master,fatal
 use part,    only: Bxyz,rhoh,mhd
 use physcon, only: pi,fourpi,qe,c,mass_proton_cgs
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i
 real                         :: vmax(3),vrms(3),vaave,B2,vA
 real                         :: h0,quada,quadb,quadc,omegaI,omegaR,hoft,pdiffW
 real                         :: termANA,termACT,ratioANA,ratioACT,pdiffR,npart1
 real                         :: Bave(3),Bzrms,Bzrmsnoa
 character(len=200)           :: fileout,filename
 logical                      :: iexist
 !
 !--Open file (appendif exists)
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_bzrms.dat'
 inquire(file=fileout,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    Bzrmsc  = 0.0
    Bzrmsav = 0.0
    write(iunit,"('#',19(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',       &
          2,'Bz_rms',     &
          3,'Bz_ave',     &
          4,'Bz_rmsnoave',&
          5,'h(t)',       &
          6,'% diff_wave',&
          7,'vx_max',     &
          8,'vy_max',     &
          9,'vz_max',     &
         10,'CRMSE',      &
         11,'vx_rms',     &
         12,'vy_rms',     &
         13,'vz_rms',     &
         14,'Bz_timeave', &
         15,'vA_ave',     &
         16,'Bx_ave',     &
         17,'w/kv_ana',   &
         18,'w/kv_act',   &
         19,'% diff_w/kv'
 else
    open(iunit,file=fileout,position='append')
 endif
 !
 !--Read the setup file to get the values of interest
 filename=trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,mhd)
 elseif (id==master) then
    rhoin0     = 1.0
    Bxin0      = 1.0
    amplitude  = 0.01
    kwave      = 2.0
 endif
 !
 !--Get coefficient values from the .in file
 filename=trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'.in'
 inquire(file=filename,exist=iexist)
 C_AD = 0.0
 C_HE = 0.0
 if (iexist) call read_infile(filename)
 if (ambitest .and. C_AD==0.0) call fatal('analysis','cannot perform ambi test without C_AD from the .in file')
 if (halltest .and. C_HE==0.0) call fatal('analysis','cannot perform Hall test without C_HE from the .in file')
 !
 !--Calculate the average & maximums
 Bave  = 0.0
 vmax  = 0.0
 vaave = 0.0
 do i = 1,npart
    B2        = dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
    Bave(1:3) = Bave(1:3) + Bxyz(1:3,i)
    vaave     = vaave     + sqrt( B2 / rhoh(xyzh(4,i),particlemass) )
    vmax(1)   = max(vmax(1),abs(vxyzu(1,i)))
    vmax(2)   = max(vmax(2),abs(vxyzu(2,i)))
    vmax(3)   = max(vmax(3),abs(vxyzu(3,i)))
 enddo
 npart1   = 1.0/npart
 Bave     = Bave *npart1
 vaave    = vaave*npart1
 !
 !--Calculate the rms average
 Bzrms    = 0.0 ! subtracting average
 Bzrmsnoa = 0.0 ! not subtracting average
 vrms     = 0.0
 do i = 1,npart
    Bzrms     = Bzrms     + (Bxyz(3,i)-Bave(3))**2
    Bzrmsnoa  = Bzrmsnoa  + Bxyz(3,i)**2
    vrms(1:3) = vrms(1:3) + vxyzu(1:3,i)**2
 enddo
 Bzrms    = sqrt(Bzrms   *npart1)
 vrms     = sqrt(vrms    *npart1)
 Bzrmsnoa = sqrt(Bzrmsnoa*npart1)
 !
 !--Produce no results if ambitest .and. halltest
 if (ambitest .and. halltest) then
    ! return zeroed values
    pdiffW   = 0.0
    hoft     = 0.0
    ratioANA = 0.0
    ratioACT = 0.0
    pdiffR   = 0.0
 else
    vA      = Bxin0/sqrt(rhoin0)
    h0      = amplitude*Bxin0/sqrt(2.0)
    if (ambitest) then
       !--Solution to the ambipolar diffusion equation (Choi test)
       write(*,*) "Ambipolar diffusion being tested"
       write(*,*) "Calculating analytical results using (B_0,rho,v_amp,C_AD) = " &
                  ,Bxin0,rhoin0,amplitude,C_AD
       quada   = 1.0
       quadb   = (vA*kwave*pi)**2*C_AD
       quadc   = -(vA*kwave*pi)**2
       omegaI  = -0.5*quadb
       omegaR  = 0.5*sqrt( -quadb*quadb - 4.0*quada*quadc )
       hoft    = h0*abs(sin(omegaR*time))*exp(omegaI*time)
       pdiffR  = 0.0
    else if (halltest) then
       !--Solution to the Hall wave equation (Sano & Stone)
       write(*,*) "Hall effect being tested"
       write(*,*) "Calculating analytical results using (B_0,rho,eta_hall,k/pi,v_amp) = " &
                 ,Bxin0,rhoin0,C_HE,kwave,amplitude
       termANA  = C_HE*Bxin0*kwave*pi/(2.0*vA)
       ratioANA = termANA + sqrt(termANA**2 + 1.0)
       termACT  = C_HE*Bave(1)*kwave*pi/(2.0*vaave)
       ratioACT = termACT + sqrt(termACT**2 + 1.0)
       pdiffR   = abs(ratioANA - ratioACT)/ratioANA
       !
       quada    = 1.0
       quadb    = -C_HE*Bxin0*(kwave*pi)**2
       quadc    = -(vA*kwave*pi)**2
       omegaR   = 0.5*( -quadb + sqrt(quadb**2 - 4.0*quada*quadc ) )
       hoft     = h0*abs(sin(omegaR*time))
    else
       !--Solution to the ideal wave equation
       write(*,*) "Ideal MHD being tested"
       write(*,*) "Calculating analytical results using (B_0,rho,v_amp) = " &
                 ,Bxin0,rhoin0,amplitude
       omegaR   = vA*kwave*pi
       hoft     = h0*abs(sin(omegaR*time))
    endif
    if (hoft > 0.0) then
       pdiffW = abs(hoft - Bzrms)/hoft
    else
       pdiffW = 0.0
    endif
    Bzrmsc    = Bzrmsc  + (hoft-Bzrms)**2
    Bzrmsav   = Bzrmsav + hoft
 endif
 !
 !--Write results to file
 write(iunit,'(19(es18.10,1x))') &
    time, Bzrms, Bave(3), Bzrmsnoa, hoft,pdiffW ,vmax(1:3), &
    sqrt(Bzrmsc/float(num)),vrms, Bzrmsav/sqrt(float(num)), &
    vaave,Bave(1),ratioANA,ratioACT,pdiffR
 close(iunit)
 !
end subroutine do_analysis
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,mhd)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in) :: filename
 logical,          intent(in) :: mhd
 integer, parameter           :: iunit = 21
 integer                      :: ierr
 type(inopts), allocatable    :: db(:)
 !
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(rhoin0,'rhoin',db,ierr)
 if (mhd) call read_inopt(Bxin0,'Bxin',db,ierr)
 call read_inopt(amplitude,'amplitude',db,ierr)
 call read_inopt(kwave,'kwave',db,ierr)
 if (mhd) then
    call read_inopt(ambitest,'ambitest',db,ierr)
    call read_inopt(halltest,'halltest',db,ierr)
 endif
 !
 call close_db(db)
 !
end subroutine read_setupfile
!-----------------------------------------------------------------------
subroutine read_infile(filename)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 21
 integer                      :: ierr
 type(inopts), allocatable    :: db(:)
 !
 print "(a)",' reading in options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(C_HE,'C_HE',db,ierr)
 call read_inopt(C_AD,'C_AD',db,ierr)
 call close_db(db)
 !
end subroutine read_infile
!-----------------------------------------------------------------------
end module
