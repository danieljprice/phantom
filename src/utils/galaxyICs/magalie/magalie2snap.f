      program galie_to_snapshot

      !-------------------------------------------------------------------
      !Code to read in magalie data file 
      !Reads in 'm.dat' setup files built from magalie (NEMO)
      !
      !Heavily modified by Alex Pettitt from a similar code included in Gadget2
      !-------------------------------------------------------------------

      implicit none

      !GADGET STUFF
      integer, parameter :: SNAPSHOT = 0       ! number of dump
      integer, parameter :: FILES = 1          ! number of files per snapshot
      character*200, parameter::path=''
      character*200 filename, snumber, fnumber,cack,gasit,ojunk
      integer*4 npart(0:5), nall(0:5)
      real*8    massarr(0:5)
      real*8    a,junk
      real*8    redshift
      integer*4 unused(34)
      integer*4 fn,i,nstart,ns,flag_sfr,flag_feedback
      integer*4 N,Ntot,j,iall,ntypes
      real*4,allocatable    :: pos(:,:),vel(:,:)
      integer*4,allocatable :: id(:),temp(:)
      real*4,allocatable    :: PartPos(:,:),dummy(:,:),unown(:)

      !ARP asciifile stuff:
      real*4,allocatable    ::xyzD(:,:),vxyzD(:,:),iphase(:),mh(:,:)
      real*4,allocatable    ::xyzB(:,:),vxyzB(:,:)
      real*4,allocatable    ::xyzH(:,:),vxyzH(:,:)
      real*4,allocatable    ::xyzT(:,:),vxyzT(:,:)
      real*4,allocatable    ::xyzG(:,:),vxyzG(:,:)
      integer*4             ::n1,nbulge,ndisc,nhalo,nbndry,ngas,nstar
      real*4                ::m1,mbulge,mdisc,mhalo,mbndry,mgas,mstar
      real*4                ::umass,udist,uvel,Mo,gasfac,Dno,To,Ho
      real*4                ::x1,y1,vx1,vy1,gasoff,unit_scale,kms,kpc
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !BEGIN READ OF ASCII DATA
      !Convert NEMO units to Gadget units

      !Re-scale using:
      !         G  =  v^2 * R /  M
      !   where we set G=1
      !Set:
      kms   = 1e3
      kpc   = 3.086e19
      Mo    = 1.99e30      

      !#######IMPORTANT: these values are used to re-sale the disc
      print*,'WARNING:----------------------------------------'
      print*,'The code requires some distance and vel scales'
      print*,'I advise you look at the plotted rotation curves'
      print*,'to ensure this looks like what you want.'
      !The scales for uvel and udist can be changed as we wish.
      uvel  = 200.
      udist = 3.5
      print*,'uvel  = ',uvel
      print*,'udist = ',udist
      print*,'            <HIT RETURN TO COMPLY>'
      print*,'WARNING:----------------------------------------'
      read(*,*)
      umass = ((uvel*kms*uvel*kms) * (udist*kpc)) /Mo/ 6.67e-11
      umass = umass/1e10  !umass for GADGET/GIZMO in 1e10Mo

      print*,' Enter gas initial thermal energy in Kelvin'
      print*,' [if you are galaxying, I advice 100-10000K]'
      read*,To
      Dno    = 0.
      Ho     = 0.
      print*,' Enter gas to stellar disc mass ratio [0.05-0.5 advised]'
      read*,gasfac
      gasoff = 90.            !angular shift of gas wrt stellar disc.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      !Read in the junk from the .in file (not stored in m.dat):
      print*,'Openning paramter file...'
      OPEN(20,file='magalie.in',form='formatted')
      DO i=1,55  !Set no. of lines
        if (i.eq.3) then
          read(20,*)ndisc,cack
        elseif (i.eq.25) then
          read(20,*)nbulge,cack
        elseif (i.eq.37) then
          read(20,*)nhalo,cack
        else
          read(20,*)ojunk
        endif
      END DO
      print*,'Read in complete, closing parameter file'
      CLOSE(20)

      nstar  = 0
      nbndry = 0
      ngas   = ndisc
      print*,'Number of particles:'
      print*,'-ndisc ',ndisc
      print*,'-ngas  ',ngas
      print*,'-nbulge',nbulge
      print*,'-nhalo ',nhalo
      print*,'-nstar ',nstar
      print*,'-nbndry',nbndry
      print*,'-------------------------------------------------'
      read(*,*)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      print*,"Reading in this order (no gas in magalie)"
      print*,'1.disc, 2.bulge 3.halo'

      print*,'Openning data file...'
      OPEN(22,file='m.dat',form='formatted')
      read(22,*)ntot,ntypes,junk
      ntot = ntot + ngas

      allocate(xyzD(3,ndisc),vxyzD(3,ndisc))
      print*,'Reading in',ndisc,'particles from disc...'
      DO i=1,ndisc
        READ(22,*)mdisc,xyzD(1,i),xyzD(2,i),xyzD(3,i),
     &          vxyzD(1,i),vxyzD(2,i),vxyzD(3,i)
        xyzD(1,i) = xyzD(1,i)*udist
        xyzD(2,i) = xyzD(2,i)*udist
        xyzD(3,i) = xyzD(3,i)*udist
        vxyzD(1,i)= vxyzD(1,i)*uvel
        vxyzD(2,i)= vxyzD(2,i)*uvel
        vxyzD(3,i)= vxyzD(3,i)*uvel
      END DO
      allocate(xyzB(3,nbulge),vxyzB(3,nbulge))
      print*,'Reading in',nbulge,'particles from bulge...'
      DO i=1,nbulge
        READ(22,*)mbulge,xyzB(1,i),xyzB(2,i),xyzB(3,i),
     &          vxyzB(1,i),vxyzB(2,i),vxyzB(3,i)
        xyzB(1,i) = xyzB(1,i)*udist
        xyzB(2,i) = xyzB(2,i)*udist
        xyzB(3,i) = xyzB(3,i)*udist
        vxyzB(1,i)= vxyzB(1,i)*uvel
        vxyzB(2,i)= vxyzB(2,i)*uvel
        vxyzB(3,i)= vxyzB(3,i)*uvel
      END DO
      allocate(xyzH(3,nhalo),vxyzH(3,nhalo))
      print*,'Reading in',nhalo,'particles from halo...'
      DO i=1,nhalo
        READ(22,*)mhalo,xyzH(1,i),xyzH(2,i),xyzH(3,i),
     &          vxyzH(1,i),vxyzH(2,i),vxyzH(3,i)
        xyzH(1,i) = xyzH(1,i)*udist
        xyzH(2,i) = xyzH(2,i)*udist
        xyzH(3,i) = xyzH(3,i)*udist
        vxyzH(1,i)= vxyzH(1,i)*uvel
        vxyzH(2,i)= vxyzH(2,i)*uvel
        vxyzH(3,i)= vxyzH(3,i)*uvel
      END DO
      
      print*,'Read in complete, closing file.'
      CLOSE(22)

      
      print*,'Makin some gas.'
      print*,'We are simply taking the stellar disc,'
      print*,'reducing mass and rotating by.',gasoff
      print*,'*ngas must be ndisc for now.'
      gasoff = gasoff  *  3.4159/180.
      i=0
      allocate(xyzG(3,ngas),vxyzG(3,ngas),temp(ngas))
      do i=1,ngas
        xyzG(1,i) =  xyzD(1,i)*COS(gasoff) - xyzD(2,i)*SIN(gasoff)
        xyzG(2,i) =  xyzD(1,i)*SIN(gasoff) + xyzD(2,i)*COS(gasoff)
        xyzG(3,i) =  xyzD(3,i)
        vxyzG(1,i)= vxyzD(1,i)*COS(gasoff)-vxyzD(2,i)*SIN(gasoff)
        vxyzG(2,i)= vxyzD(1,i)*SIN(gasoff)+vxyzD(2,i)*COS(gasoff)
        vxyzG(3,i)= vxyzD(3,i)
      end do
      print*,'Set set, assigning mass.'
      print*,'Re-scaling the disc mass by',gasfac,'for gas.'
      mgas = gasfac*mdisc

      if (ntot .ne. (ngas+nhalo+ndisc+nbulge+nstar+nbndry)) then
        print*, '---OAHHHH YOU MUCKED UP'
        print*, 'particle numbers dont add up'
        stop
      endif
      mgas   = mgas*umass
      mhalo  = mhalo*umass
      mdisc  = mdisc*umass
      mbulge = mbulge*umass
      mstar   = 0.
      mbndry  = 0.

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !BEGIN WRITE OUT OF GADGET2 BINARY FORMAT
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      
      filename = 'newsnap_000'
      print *,'opening...  '//filename(1:len_trim(filename))     
      open (2, file=filename, form='unformatted')
      print*,'Filling header arrays...'
      npart(0) = ngas
      npart(1) = nhalo
      npart(2) = ndisc
      npart(3) = nbulge
      npart(4) = nstar
      npart(5) = nbndry
      massarr(0) = mgas
      massarr(1) = mhalo
      massarr(2) = mdisc
      massarr(3) = mbulge
      massarr(4) = mstar
      massarr(5) = mbndry
      a             = 0.
      redshift      = 0.
      flag_sfr      = 0
      flag_feedback = 0
      nall          = npart
      unused        = 0

      print*,'Writing ascii data into snapshot array'
      write (2)npart,massarr,a,redshift,flag_sfr,flag_feedback,
     & nall,unused
      print*,npart,massarr,a,redshift,flag_sfr,flag_feedback,
     & nall,unused

      allocate(xyzT(3,ntot),vxyzT(3,ntot))
      iall=1
      print*,'Making gas array'
      do i=1,ngas
        xyzT(1,iall)=xyzG(1,i)
       	xyzT(2,iall)=xyzG(2,i)
        xyzT(3,iall)=xyzG(3,i)
        vxyzT(1,iall)=-vxyzG(1,i)
        vxyzT(2,iall)=-vxyzG(2,i)
        vxyzT(3,iall)=vxyzG(3,i)
        iall=iall+1
      end do
      print*,'Making halo array'
      do i=1,nhalo
        xyzT(1,iall)=xyzH(1,i)
       	xyzT(2,iall)=xyzH(2,i)
        xyzT(3,iall)=xyzH(3,i)
        vxyzT(1,iall)=-vxyzH(1,i)
        vxyzT(2,iall)=-vxyzH(2,i)
        vxyzT(3,iall)=vxyzH(3,i)
        iall=iall+1
      end do
      print*,'Making disc array'
      do i=1,ndisc
        xyzT(1,iall)=xyzD(1,i)
        xyzT(2,iall)=xyzD(2,i)
        xyzT(3,iall)=xyzD(3,i)
        vxyzT(1,iall)=-vxyzD(1,i)
        vxyzT(2,iall)=-vxyzD(2,i)
        vxyzT(3,iall)=vxyzD(3,i)
       	iall=iall+1
      end do
      print*,'Making bulge array'
      do i=1,nbulge
        xyzT(1,iall)=xyzB(1,i)
        xyzT(2,iall)=xyzB(2,i)
        xyzT(3,iall)=xyzB(3,i)
        vxyzT(1,iall)=-vxyzB(1,i)
        vxyzT(2,iall)=-vxyzB(2,i)
        vxyzT(3,iall)=vxyzB(3,i)
       	iall=iall+1
      end do

      print*,'Writing all data (x)'
      write (2)xyzT
      print*,'Writing all data (v)'
      write (2)vxyzT

      !--Fill ID's for all components:
      allocate(unown(ntot))
      do i=1,ntot
        unown(i)=i
      end do
      !--Fill hydrodynamics for gas ony:
      do i=1,ngas
        !unown(i)=i
        temp(i)=To
        !dens(i)=Do
        !smh(i) =Ho
      end do

      write (2) unown  !ID's
      write (2) temp   !Temps (energy for gas)
      !write (2) dens   !Density
      !write (2) smh    !Sm. length
      print*,'Assigning values of:'
      print*,'  T=',To
      print*,'rho=',Dno
      print*,'  h=',Ho

      close (2)
      print*,'Closing file'

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !DONE, BELOW IS FOR THE READ IN.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      print*,' <<<<<<<<<<FIN>>>>>>>>>>>'
      stop
      
      end program galie_to_snapshot

