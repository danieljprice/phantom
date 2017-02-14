      program ascii_to_snapshot

      !-------------------------------------------------------------------
      !Code to read in asciifiles 
      !Reads in 'asciifile' setup files built for sphNG runs
      !
      !Heavily modified by Alex Pettitt from a similar code included in Gadget2
      !-------------------------------------------------------------------

      implicit none
      !GADGET STUFF
      integer, parameter :: SNAPSHOT = 0       ! number of dump
      integer, parameter :: FILES = 1          ! number of files per snapshot
      character*200, parameter::path=''
      character*200 filename, snumber, fnumber
      integer*4 npart(0:5), nall(0:5)
      real*8    massarr(0:5)
      real*8    a,Ho,Dno,To
      real*8    redshift
      integer*4 unused(34)
      integer*4 fn,i,nstart,ns,flag_sfr,flag_feedback
      integer*4 N,Ntot,j
      real*4,allocatable    :: pos(:,:),vel(:,:)
      integer*4,allocatable :: id(:),temp(:),dens(:),smh(:)
      real*4,allocatable    :: PartPos(:,:),dummy(:,:),unown(:)
    
      !ARP asciifile stuff:
      real*4,allocatable    ::xyz(:,:),vxyz(:,:),iphase(:),mh(:,:)
      integer*4             ::n1,nbulge,ndisc,nhalo,nbndry,ngas,nstar
      real*4                ::m1,mbulge,mdisc,mhalo,mbndry,mgas,mstar
      real*4                ::umass,udist,uvel,Mo
     
      character(100)        ::galsetupic,sometext
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !--BEGIN READ OF ASCII DATA
      print*,'Reading in metadata from txt file...'
      galsetupic = 'galsetic.txt'
      OPEN(21,file=galsetupic,form='formatted')
      do i=1,5
        if (i.eq.1) then
          read(21,*)sometext,ngas
        elseif (i.eq.2) then
          read(21,*),sometext,ndisc
        elseif (i.eq.3) then
          read(21,*),sometext,nbulge
        elseif (i.eq.4) then
          read(21,*),sometext,nhalo
        endif
      enddo
      nstar  = 0.
      mstar  = 0.
      mbndry = 0.
      nbndry = 0.

      close(21)
      print*,'Read in the IC parameter file: ',galsetupic
      print*,' ngas   =',ngas
      print*,' ndisc  =',ndisc
      print*,' nbulge =',nbulge
      print*,' nhalo  =',nhalo
      print*,' nstar  =',nstar
      print*,' nbndry =',nbndry
      n1     = ngas+ndisc+nbulge+nhalo+nstar+nbndry

      print*,' Enter gas initial thermal energy in Kelvin'
      print*,' [if you are galaxying, I advice 100-10000K]'
      read*,To
      Dno    = 0.
      Ho     = 0.
      Mo     = 1.99e30        !convert to 1e10Mo
      udist  = 3.0856e19      !convert to kpc
      uvel   = 1.e3           !convert to km/s
      umass  = 1.0e10         !mass unit used by gadget
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      
      print*,"NOTE: DATA MUST BE WRITTEN OUT IN THIS ORDER"
      print*,'1.gas, 2.halo, 3.disc, 4.bulge etc.'
      print*,'MESS WITH AT YOUR OWN PERIL'
      
      allocate(xyz(3,n1),vxyz(3,n1),mh(2,n1),iphase(n1),unown(n1)) 
      print*,'If you want to read in rho and T, need to edit source.'
      !allocate(temp(n1),dens(n1),smh(n1))

      print*,'Reading in',n1,'particles from asciifile...'

      if (ngas.ne.0) then
        OPEN(23,file='asciifile_G',form='formatted')
        write(*,*)'Reading gas particles without h'
        do while(i <= ngas)
          READ(23,*)xyz(1,i),xyz(2,i),xyz(3,i),mh(1,i),
     &          vxyz(1,i),vxyz(2,i),vxyz(3,i),iphase(i)
          xyz(1,i) = xyz(1,i)/udist
          xyz(2,i) = xyz(2,i)/udist
          xyz(3,i) = xyz(3,i)/udist
          mgas     = mh(1,i)
          vxyz(1,i)= vxyz(1,i)/uvel
          vxyz(2,i)= vxyz(2,i)/uvel
          vxyz(3,i)= vxyz(3,i)/uvel
          i = i + 1
        END DO
        CLOSE(23)
      else
        mgas = 0.
      endif

      if (nhalo.ne.0) then
        OPEN(24,file='asciifile_H',form='formatted')
        write(*,*)'Reading halo particles'
        do while(i <= ngas+nhalo)
          READ(24,*)xyz(1,i),xyz(2,i),xyz(3,i),mh(1,i),
     &          vxyz(1,i),vxyz(2,i),vxyz(3,i),iphase(i)
          xyz(1,i) = xyz(1,i)/udist
          xyz(2,i) = xyz(2,i)/udist
          xyz(3,i) = xyz(3,i)/udist
          mhalo     = mh(1,i)
          vxyz(1,i)= vxyz(1,i)/uvel
          vxyz(2,i)= vxyz(2,i)/uvel
          vxyz(3,i)= vxyz(3,i)/uvel
          i = i + 1
        END DO
        CLOSE(24)
      else
        mhalo = 0.
      endif

      if (ndisc.ne.0) then
        OPEN(25,file='asciifile_D',form='formatted')
        write(*,*)'Reading disc particles'
        do while(i <= ngas+nhalo+ndisc)
          READ(25,*)xyz(1,i),xyz(2,i),xyz(3,i),mh(1,i),
     &          vxyz(1,i),vxyz(2,i),vxyz(3,i),iphase(i)
          xyz(1,i) = xyz(1,i)/udist
          xyz(2,i) = xyz(2,i)/udist
          xyz(3,i) = xyz(3,i)/udist
          mdisc     = mh(1,i)
          vxyz(1,i)= vxyz(1,i)/uvel
          vxyz(2,i)= vxyz(2,i)/uvel
          vxyz(3,i)= vxyz(3,i)/uvel
          i = i + 1
        END DO
        CLOSE(25)
      else
        mdisc = 0.
      endif

      if (nbulge.ne.0) then
        OPEN(26,file='asciifile_B',form='formatted')
        write(*,*)'Reading gas particles without h'
        do while(i <= ngas+nhalo+ndisc+nbulge)
          READ(26,*)xyz(1,i),xyz(2,i),xyz(3,i),mh(1,i),
     &          vxyz(1,i),vxyz(2,i),vxyz(3,i),iphase(i)
          xyz(1,i) = xyz(1,i)/udist
          xyz(2,i) = xyz(2,i)/udist
          xyz(3,i) = xyz(3,i)/udist
          mbulge     = mh(1,i)
          vxyz(1,i)= vxyz(1,i)/uvel
          vxyz(2,i)= vxyz(2,i)/uvel
          vxyz(3,i)= vxyz(3,i)/uvel
          i = i + 1
        END DO
        CLOSE(26)
      else
        mdisc = 0.
      endif

      print*,'Read in complete, closing asciifiles.'
      print*,'Note: no smoothing length set, dealt with by Gadget2.'

      !--converting to gadget units:
      mgas   = mgas/Mo/umass   
      mhalo  = mhalo/Mo/umass
      mdisc  = mdisc/Mo/umass  
      mbulge = mbulge/Mo/umass 

      allocate(temp(ngas),dens(ngas),smh(ngas))
      
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
      write (2)xyz
      write (2)vxyz

      !--Fill ID's for all components:
      do i=1,n1 
        unown(i)=i
      end do
      !--Fill hydrodynamics for gas ony:
      do i=1,ngas
        temp(i)=To
        !dens(i)=Dno
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

      print*,' <<<<<<<<<<FIN>>>>>>>>>>>'
      stop
      
      end program ascii_to_snapshot

