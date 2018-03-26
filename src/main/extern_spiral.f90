!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_spiral
!
!  DESCRIPTION:
!    This module contains everything related to
!    Clare and Alex's galactic potentials:
!
!  REFERENCES:
!    Cox and Gomez (2002), ApJS 142, 261-267
!    Pichardo et al. (2003), ApJ 582, 230-245
!    Freeman (1970), ApJ 160, 811-830
!    Miyamoto and Nagi (1975), PASJ 27, 533-543
!    Dobbs, Bonnell & Pringle (2006), MNRAS 371, 1663-1674
!    Cladwell and Ostriker (1981), APJ 251, 61-87
!    Allen and Martos (1986), RMAA 13, 137-147
!    Khorperskov et al. (2012), arXiv:1207.5162v1
!    Long and Murali (1992), APJ 397, 44-48
!    Voigt and Letelier (2011), MNRAS 411, 2371-2382
!    and other basic forms (see Binney and Tremaine 1987)
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    NN     -- No of arms in stellar spiral potential
!    a_bar  -- Major axis of galactic bar (in x, kpc)
!    b_bar  -- Minor axis of galactic bar (in y, kpc)
!    c_bar  -- Minor axis of galactic bar (in z, kpc)
!    iarms  -- type of arm potential (1:C&G 2:4 P&M spheroids+linear)
!    ibulg  -- type of bulge potential (1:Plummer 2:Hernquist 3:Hubble)
!    idisk  -- type of disk potential (1:log 2:flattened 3:Freeman w. bar/arm)
!    ihalo  -- type of halo potential (1:C&O 2:Flat 3:A&M 4:K&B 5:NFW)
!    iread  -- Read in potential from file (1=y,0=n)
!    phib   -- Bar(s) potential pattern speed (km/s/kpc)
!    phir   -- Spiral potential pattern speed (km/s/kpc)
!    pitchA -- Pitch angle of spiral arms (deg)
!
!  DEPENDENCIES: infile_utils, io, mathfunc, physcon, units
!+
!--------------------------------------------------------------------------
module extern_spiral
 implicit none
 public :: s_potential,schmidt_potential,pichardo_potential,initialise_spiral,&
           LogDisc,MNDisc,KFDiscSp,PlumBul,HernBul,HubbBul,COhalo,Flathalo,AMhalo,&
           KBhalo,LMXbar,LMTbar,Orthog_basisbar,DehnenBar,VogtSbar,Wadabar,&
           BINReadPot3D,NFWhalo
 public :: write_options_spiral, read_options_spiral
 real(kind=8) :: r0,rS,p0,Hz,Co,Rc,Rtot,rhalomax,HaloCi,dzq
 real(kind=8) :: Cz(3)
 real(kind=8) :: rc2,rc21,p1,sinalpha,cotalpha,t0,Rtotterm,gcode,strength,rcore,rdcore
 real(kind=8) :: barmass,ebar,rHubbMax,zdisk,MHubbbulge,RHubbbulge,betabar,epsBar
 real(kind=8) :: Mdisk,adisk,bdisk,Mbulge,Rbulge,vhalo,rhalo,Mhalo,HubbCi
 real(kind=8) :: a_0,c_0,d_0,e_0,Ri,Rf,Rl,Rsarms,Nshape,Mspiral,a_02,c_02,e_02
 real(kind=8) :: ra1,ra2,ra3,rsK,phi0,eta_Khop,tau_arm,rbars,rho0b,rlim,StrBarW, &
       StrBarD,radiusofbar,VoB,Rob,alphaBar,MNmdisk,MNadisk,MNbdisk,AMmhalo,AMrhalo,rLhalo,&
       potxmax,potymax,potzmax,dxpot,dypot,dzpot,Mnfw,Cnfw,rnfw
 real(kind=8), allocatable :: spiralsum(:)
 real(kind=8), allocatable :: Rspheroids(:,:),shapefn(:,:),den0(:,:)
 integer :: Nt,NNi,j,jj,i,ios
 character(len=100) :: potfilename
 integer :: potlenx, potleny, potlenz
 real(kind=8), allocatable :: newpot3D(:,:,:)


!
!--the following are parameters to be written/read from the input file
!
 integer      :: idisk=  1       !Disk id flag
!--1: B&T log-disk
!--2: M&N flattened-disk
!--3: K&F Freeman-disk with build in bar/arm pretubation
 integer      :: ibulg=  0        !Bulge id flag
!--1: Plummer bulge
!--2: Henrquist bulge
!--3: Hubble bulge
 integer      :: ihalo=  0        !DM halo id flag
!--1: C&O halo
!--2: Flat halo
!--3: A&M halo
!--4: K&B halo
!--5: NFW halo
 integer      :: iarms=  0        !Arm id flag
!--1:C&G TWA spirals
!--2:P&M linear density inside spheroids(a)  linear density from galactic centre(r)
!--3:P&M linear density inside spheroids(a), log density from galactic centre(r)
!--4:P&M inverse density inside spheroids(a),linear density from galactic centre(r)
 integer      :: ibar =  0
!--1:L&M biaxial bar
!--2:L&M triaxial bar
!--3:W&D G2+G3 bar combination
!--4:Dehnen cosine bar
!--5:V&L S-shape bar
!--6:Wada Bar
 integer      :: iread =  0        !Read in potential file?
!--0: Dont read in
!--1: Read in
 real(kind=8) :: NN     = 2.      ! should be integer, but used in real expressions
 real(kind=8) :: phibar = 40.0    !km/s/kpc, for bars from Dehnen99
 real(kind=8) :: phir   = 20.0    !km/s/kpc, for arms
 real(kind=8) :: alpha  = 15.0    !The arm pitch angle in deg
 real(kind=8) :: LMabar = 4.00    !Bar major axis (x1kpc)
 real(kind=8) :: LMbbar = 1.0     !Bar minor axis (x1kpc)
 real(kind=8) :: LMcbar = 1.0     !Bar minor axis (x1kpc), triaxial bar
 integer      :: initialSP = 1   !Flag to track if the above are in code units (1=no, 0=yes)

contains

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_spiral(iunit)
 use infile_utils, only:write_inopt
 use units, only:utime,udist
 use physcon, only:pc,solarm,pi,gg,kpc,km
 integer, intent(in) :: iunit

 !Convert to Galactic units:
 if (initialSP==0) then
    phir  = phir    /(km/kpc * utime)
    phibar= phibar  /(km/kpc * utime)
    alpha = alpha   /(pi/180.0)
    LMabar= LMabar  /(kpc/udist)
    LMbbar= LMbbar  /(kpc/udist)
    LMcbar= LMcbar  /(kpc/udist)
 endif

 write(iunit,"(/,a)") '# options relating to spiral potentials'
 call write_inopt(idisk,'idisk','type of disk potential (1:log 2:flattened 3:Freeman w. bar/arm)',iunit)
 call write_inopt(ibulg,'ibulg','type of bulge potential (1:Plummer 2:Hernquist 3:Hubble)',iunit)
 call write_inopt(ihalo,'ihalo','type of halo potential (1:C&O 2:Flat 3:A&M 4:K&B 5:NFW)',iunit)
 call write_inopt(iarms,'iarms','type of arm potential (1:C&G 2:4 P&M spheroids+linear)',iunit)
 call write_inopt(ibar,'ibar',&
     'type of bar potential (1:biaxial 2:triaxial 3:G2G3 4:Dehnen cos 5:S-shape 6:Wada cos)',iunit)
 call write_inopt(iread,'iread','Read in potential from file (1=y,0=n)',iunit)
 call write_inopt(NN,'NN','No of arms in stellar spiral potential',iunit)
 call write_inopt(alpha,'pitchA','Pitch angle of spiral arms (deg)',iunit)
 call write_inopt(phir,'phir','Spiral potential pattern speed (km/s/kpc)',iunit)
 call write_inopt(phibar,'phib','Bar(s) potential pattern speed (km/s/kpc)',iunit)
 call write_inopt(LMabar,'a_bar','Major axis of galactic bar (in x, kpc)',iunit)
 call write_inopt(LMbbar,'b_bar','Minor axis of galactic bar (in y, kpc)',iunit)
 call write_inopt(LMcbar,'c_bar','Minor axis of galactic bar (in z, kpc)',iunit)

 !Convert back to code units:
 if (initialSP==0) then
    phir  = phir    *(km/kpc * utime)
    phibar= phibar  *(km/kpc * utime)
    alpha = alpha   *(pi/180.0)
    LMabar= LMabar  *(kpc/udist)
    LMbbar= LMbbar  *(kpc/udist)
    LMcbar= LMcbar  *(kpc/udist)
 endif

end subroutine write_options_spiral

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_spiral(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 use units, only:utime,udist
 use physcon, only:pc,solarm,pi,gg,kpc,km
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_spiral'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('idisk')
    read(valstring,*,iostat=ierr) idisk
    if (idisk < 0) call fatal(label,'idisk < 0')
    ngot = ngot + 1
 case('ibulg')
    read(valstring,*,iostat=ierr) ibulg
    if (ibulg < 0) call fatal(label,'ibulg < 0')
 case('ihalo')
    read(valstring,*,iostat=ierr) ihalo
    if (ihalo < 0) call fatal(label,'ihalo < 0')
 case('iarms')
    read(valstring,*,iostat=ierr) iarms
    if (iarms < 0) call fatal(label,'iarms < 0')
 case('ibar')
    read(valstring,*,iostat=ierr) ibar
    if (ibar < 0) call fatal(label,'ibar < 0')
 case('iread')
    read(valstring,*,iostat=ierr) iread
 case('a_bar')
    read(valstring,*,iostat=ierr) LMabar
    if (LMabar < 0) call fatal(label,'a_bar < 0')
 case('b_bar')
    read(valstring,*,iostat=ierr) LMbbar
    if (LMbbar < 0) call fatal(label,'b_bar < 0')
 case('c_bar')
    read(valstring,*,iostat=ierr) LMcbar
    if (LMcbar < 0) call fatal(label,'c_bar < 0')
 case('pitchA')
    read(valstring,*,iostat=ierr) alpha
    if (alpha < 0) call fatal(label,'pitchA < 0')
 case('NN')
    read(valstring,*,iostat=ierr) NN
    if (NN < 0) call fatal(label,'NN < 0')
 case('phib')
    read(valstring,*,iostat=ierr) phibar
    if (phibar < 0) call fatal(label,'phib < 0')
 case('phir')
    read(valstring,*,iostat=ierr) phir
    if (phir < 0) call fatal(label,'phir < 0')
 case default
    imatch = .false.
 end select

 !Convert to code units:
 if (initialSP==0.) then
    phir  = phir    *(km/kpc * utime)
    phibar= phibar  *(km/kpc * utime)
    alpha = alpha   *(pi/180.0)
    LMabar= LMabar  *(kpc/udist)
    LMbbar= LMbbar  *(kpc/udist)
    LMcbar= LMcbar  *(kpc/udist)
    initialSP = 0
 endif

 igotall = (ngot >= 1)  !+1 for each required potential

end subroutine read_options_spiral

!----------------------------------------------------------------
subroutine initialise_spiral(ierr)
 use physcon, only:pc,solarm,pi,gg,kpc,km
 use units,   only:udist,umass,utime
 use io,      only:id,master
 integer, intent(out) :: ierr

 ierr = 0
 ! Spiral potential parameters: initialise the first time only
 ! cgs units are entered raw, divide by udist etc to set to natural units of PHANTOM

 !--Time at initial conditions (100Myr)
 t0=3.153d+15/utime
 gcode= gg*umass*utime**2/(udist**3)
 !--Convert the below from km/s/kpc to code units (rad/Myr) etc.
 phir  = (phir   *km/kpc) * utime
 phibar= (phibar *km/kpc) * utime
 alpha = alpha*pi/180.0
 LMabar= LMabar*kpc/udist
 LMbbar= LMbbar*kpc/udist
 LMcbar= LMcbar*kpc/udist
 initialSP = 0   !integration has started, so now the above are in code units.

 !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 !<><><><><><><><><><><><><>SETUP-PARAMETERS<><><><><><><><><><><><><><><>
 !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

 !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=DISKS=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 !--Dobbs&B&P log-disk potential.
 Rc  = 1.0 *kpc/udist        !1kpc was 0.1kpc, C.Dobbs changed for ring
 Co  = (2.31d14*(utime)**2/udist**2) !Co=(Vc^2)/2 at Ro(=8kpc) Vc=215km/s
 dzq = 1./0.7                        !z potential constant
 !--Miyamoto&Naigi disk
 MNmdisk = 8.56 * 1.0d10*solarm/umass
 MNadisk = 5.3178*kpc/udist
 MNbdisk = 0.2500*kpc/udist
 !--Modified Freeman (Type II) disk
 Mdisk   = 4.0 * 1.0d10*solarm/umass
 zdisk   = 0.1*kpc/udist
 adisk   = 3.0*kpc/udist      !3.5*(3.086d21/udist)

 !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=HALOS=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 !--Cladwell&Ostriker log-atanh halo
 p1   =(6.19d-31*(udist)**3/(umass)) ! rho*4pi*G, in cgs, 1000Mo/kpc3
 rc2  =(2.019d21/udist)   !0.65kpc sacle term for halo, a second rc NOT SQUARED.
 rc21 =1./rc2
 Rtot = 225.*kpc/udist                !Needed in potential, not force, unimportant
 Rtotterm = (1.+(Rtot*rc21)**2)
 !--A simple log halo
 vhalo  = (2.00d07*utime/udist)  !186, 200, 230km/s etc
 rLhalo = 12.0*kpc/udist
 !--Allen&Martos halo
 AMmhalo= 10.7068 * 1.0d10*solarm/umass
 AMrhalo= 12.00*kpc/udist        ! r=7 for Allen halo
 rlim   = 100.0*kpc/udist
 !--Begemen quasi-isothermal halo (from Khoperskov 2012)
 Mhalo   = 6.4 * 1.0d10*solarm/umass
 rhalo   = 6.00*kpc/udist
 rhalomax= 12.0*kpc/udist
 HaloCi  = 1./(rhalomax/rhalo - atan(rhalomax/rhalo ))
 !--NFW halo
 Mnfw= 63.0 * 1.0d10*solarm/umass
 Cnfw= 5.0d0
 rnfw= (122.0d0*3.086d21/Cnfw)/udist

 !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=BULGES-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 !--Plummer bulge        or
 !--Hernequist bulge
 Mbulge  = 1.4 * 1.0d10*solarm/umass
 Rbulge  = 0.39*kpc/udist
 !--Hubble bulge (only the inverse of the constant is used):
 MHubbbulge  = 3.0* 1.0d10*solarm/umass  !Khop says: 0.92*(1.99d43/umass)
 RHubbbulge  = 0.2 *kpc/udist
 rHubbMax    = 12.0*kpc/udist
 HubbCi  = 1./( log(rHubbMax/RHubbbulge + sqrt(1.+(rHubbMax/RHubbbulge)**2)) &
   - rHubbMax/(RHubbbulge*sqrt(1.+(rHubbMax/RHubbbulge)**2))  )

 !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=BARS-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 VoB         = 2.2d7 * (utime/udist)         !220km/s
 RoB         = 8.0  *  kpc/udist
 !--Long&Murali bars
 barmass = 0.625* 1.0d10*solarm/umass  ! From Lee et al. 1999
 ebar    =  0.10     !Only for non-straight bars
 !--Dwek/Ferres bars:
 rbars = kpc/udist
 rho0b = 10.0 * (1.99E33/kpc**3)  *  udist**3/umass
 !--Dehnen bar:
 alphaBar    = 0.01  !Strength paramater for the bar
 radiusofbar = 0.8*VoB/phibar  !was set to 3kpc
 betabar     = 1.0
 rdcore      = radiusofbar/(exp(1.0)-1.0)**(1./betabar)
 StrBarD     = alphaBar*VoB*VoB*(RoB/radiusofbar)**3 / 3.
 !--Wada Bar:
 epsBar      = 0.05
 rcore       = 2.d0*kpc/udist
 StrBarW     = epsBar*VoB*VoB*2.59807621   !sqrt(27/4)

 !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=SPIRALS=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 !--Global spiral parameters, pitch and velocity:
 sinalpha = sin(alpha)
 cotalpha = 1./tan(alpha)
 !--Khoperskov TWA spiral potential:
 rSK  = 3.0*kpc/udist
 ra1  = 3.0*kpc/udist
 ra2  = 7.0*kpc/udist
 ra3  = 7.0*kpc/udist
 phi0 = 0.
 eta_Khop = 0.1
 Tau_arm  = 220.0*1d6*365.0*24.0*3600.0/utime  !220Myrs
 !--Cox&Gomez TWA spiral potential
 r0=(2.47518d22/udist)    !   8kpc, reference radii
 rS=(2.16578d22/udist)    !   7kpc, potential scale length
 Hz=(5.56916d20/udist)    !   0.18kpc, z-scale height
 p0=(2.12889d-24*(udist)**3/(umass))         !ref. density
 Cz=(/ 8.d0/(3.d0*pi), 0.5d0, 8.d0/(15.d0*pi) /) !co-efs of expansion
 !--Pichardo spheroidal potential
 Mspiral=0.0175* 8.56 * 1.0d10*solarm/umass  !0.05,0.03,0.0175  of 8.56E10 Mo
 c_0=0.500*kpc/udist
 d_0=0.135*kpc/udist  !This not mathematical, just a good value (0.1 before)
 a_0=1.000*kpc/udist
 e_0=sqrt( 1 - (c_0/a_0)**2 )
 a_02=a_0*a_0
 c_02=c_0*c_0
 e_02=e_0*e_0
 Ri    =3.30*kpc/udist
 Rsarms=3.00*kpc/udist
 Rf    =12.0*kpc/udist
 Rl    =2.50*kpc/udist
 Nshape=100.0d0
 strength= -12.5663706*gcode*(sqrt(1.0-e_02)/(e_02*e_02))

 select case(iarms)
 case(2,3,4)
    if (id==master) print*,'Setting up spiral arm spheroids...'
    Nt=INT((Rf-Ri)/d_0)
    NNi=int(NN)
    if (id==master) print*,'There are ',Nt,' spheroids in each arm'
    !-Build the arrays for the radii and angle of the S.Shperoids.
    allocate(Rspheroids(NNi,Nt),shapefn(NNi,Nt),den0(NNi,Nt),&
  spiralsum(NNi))
    !-Loop over arms:
    do jj=1,NNi
       spiralsum(jj)=0.0d0
       !-Loop over spheroids
       do j=1,Nt
          Rspheroids(jj,j) = Ri+(DBLE(j)-1.d0)*d_0
          shapefn(jj,j)    = (cotalpha/Nshape) * &
    log(1.d0+(Rspheroids(jj,j)/Rsarms)**Nshape) + jj*2.0d0*pi/DBLE(NNi)
          !print*,jj,j,Rspheroids(jj,j),shapefn(jj,j)
          select case(iarms)
          case(2,4)
             !--For a linear density drop off from galactic centre:
             den0(jj,j)       = (Rf-Rspheroids(jj,j))*3.d0*Mspiral &
      / (DBLE(NNi)*pi*a_0*a_0*c_0)
             spiralsum(jj) = spiralsum(jj) + (Rf-Rspheroids(jj,j))
          case(3)
             !--For a log density drop off from galactic centre:
             den0(jj,j)       = exp((Ri-Rspheroids(jj,j))/Rl)*3.d0*Mspiral &
      / (DBLE(NNi)*pi*a_0*a_0*c_0)
             spiralsum(jj) = spiralsum(jj) + exp((Ri-Rspheroids(jj,j))/Rl)
          end select
       enddo
    enddo
    !-Properly re-scale the density array using the radius sum.
    do jj=1,NNi
       do j=1,Nt
          den0(jj,j)=den0(jj,j)/spiralsum(jj)
       enddo
    enddo
 end select

 if (id==master) then
    select case(iarms)
    case(2)
       print*,'linear density inside spheroids(a)  linear density from galactic centre(r)'
    case(3)
       print*,'linear density inside spheroids(a)  log density from galactic centre(r)'
    case(4)
       print*,'inverse density inside spheroids(a)  linear density from galactic centre(r)'
    case default
       print*,"No spiral arm spheroids set."
    end select
 endif

 select case(iread)
 case(1)
    potfilename = 'pot3D.bin'
    if (id==master) print*,'Reading in potential from an external file (BINARY): ',potfilename
    open (unit =1, file = TRIM(potfilename), status='old', form='UNFORMATTED', access='SEQUENTIAL', iostat=ios)
    if (ios /= 0 .and. id==master) then
       print*, 'Error opening file:', TRIM(potfilename)
    endif
    !Read in the grid lengths if they exist in the header.
    read(1) potlenz,potlenx,potleny
    !Read in the grid dimensions if they exist in the header.
    read(1) potzmax,potxmax,potymax
    !Read in the grid.
    potxmax     = potxmax*kpc/udist
    potymax     = potymax*kpc/udist
    potzmax     = potzmax*kpc/udist
    dxpot       = 2.d0*potxmax/potlenx
    dypot       = 2.d0*potymax/potleny
    dzpot       = 2.d0*potzmax/potlenz
    allocate(newpot3D(potlenz,potlenx,potleny))
    read (1) newpot3D
    close(1)
    if (id==master) then
       print*,'Potential file read in successfully.'
       print*,'Grid lengths [z:x:y],'
       print*, potlenz,potlenx,potleny
       print*,'Grid physical size, in kpc [z:x:y],'
       print*, potzmax*udist/kpc,potxmax*udist/kpc,potymax*udist/kpc
    endif

 case default
    if (id==master) print*,'No potential to be read in.'
 end select

end subroutine initialise_spiral


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><>POTENTIAL-SUBROUTINES><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~DISCS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--Binney&Tremmaine--:Logarithmic disc potential
!----------------------------------------------------------------
subroutine LogDisc(xi,yi,zi,d2,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: d2,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: r2term,dr2term,pot,fxi,fyi,fzi
 !--Not a true disc potential as it has a flat rotation curve, is in
 !--fact a disc+halo.
 r2term   = (Rc)**2 + d2 + (zi*dzq)**2
 dr2term  = 1./r2term
 pot      = +Co*log(r2term)  !Co=0.5*vo^2
 fxi      = -2.*Co*xi*dr2term
 fyi      = -2.*Co*yi*dr2term
 fzi      = -2.*Co*zi*dr2term*dzq**2
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine LogDisc

!--Miyamoto&Naigi--:Two parameter disc potential
!----------------------------------------------------------------
subroutine MNDisc(xi,yi,zi,d2,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: d2,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: diskterm,diskterm2,pot,fxi,fyi,fzi
 !--If adisk=0 gives a Plummer spherical potential.
 !--If bdisk=0 gives a Plummer-Kuzmin infinitely thin disk.
 diskterm = sqrt(zi*zi+MNbdisk*MNbdisk)
 diskterm2= sqrt( (MNadisk+diskterm)*(MNadisk+diskterm)+d2 )
 pot   = - gcode*MNmdisk/diskterm2
 fxi   = - xi*gcode*MNmdisk/(diskterm2*diskterm2*diskterm2)
 fyi   = - yi*gcode*MNmdisk/(diskterm2*diskterm2*diskterm2)
 fzi   = - zi*gcode*MNmdisk*(MNadisk+diskterm) &
    /(diskterm2*diskterm2*diskterm2*diskterm)
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine MNDisc

!--Khoperskov&Freeeman--:A double exp disc + optional TWA spirals
!----------------------------------------------------------------

subroutine KFDiscSp(xi,yi,zi,d2,r,phii,ti,phi,fextxi,fextyi,fextzi)
 use physcon, only:pi
 use mathfunc, only: bessk0_s,bessi0_s,bessk1_s,bessi1_s,poly,IK01A
 real, intent(in)    :: xi,yi,zi,d2,r,phii,ti
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: rratio,phioff,eta1,eta2,eta3,eta_arm,sdens_g,&
   potdis,potarm,potdisdr,potarmdr,potarmdphi,pot,fxi,fyi,fzi,&
   BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1,fphii,fri,d1,kap1,kap2,kap3
 !--Must be added as a perturbation to the disk potential.
 !--Returns to Freeman disc if eta_Khop=0
 d1     = sqrt(d2)
 rratio = d1/(2.d0*adisk)
 phioff = phi0+phir*(t0+ti)
 eta1=d1/ra1
 eta2=d1/ra2
 eta3=d1/ra3
 kap1=eta1*eta1
 kap2=eta2*eta3
 kap3=eta3*eta3
 eta_arm = MIN(eta_Khop,0.1*ti/Tau_arm)
 sdens_g = gcode* Mdisk/(2.*pi*adisk*adisk)

 !--Originally used these bessesl fns but are deemed too inaccurate for double precision,
 !--do not satisfy: force=-GRAD(potential)
 !potdis = (pi*zdisk*log(COSH(zi/zdisk)) - &
 !    (bessi0_s(rratio)*bessk1_s(rratio) - bessi1_s(rratio)*bessk0_s(rratio)) &
 !    *pi*adisk*rratio) * sdens_g
 !potdisdr  =  -(bessi0_s(rratio)*bessk0_s(rratio) - bessi1_s(rratio)*bessk1_s(rratio))* &
 !    pi*rratio  * sdens_g

 !--Disc only component in dr (dz simple and evluated at force update below)
 call   IK01A(rratio,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
 potdis  =  (pi*zdisk*log(COSH(zi/zdisk)) - (BI0*BK1 - BI1*BK0)*pi*adisk*rratio)*sdens_g
 potdisdr= (BI0*BK0 - BI1*BK1)*pi*rratio*sdens_g
 !--Spiral arm purtubtaion.
 potarm=  kap1*cos(2.d0*(phii+phioff                     ))/(1.d0+kap1)**(3./2.)+ &
             kap2*cos(2.d0*(phii+phioff-cotalpha*log(d1/rSK)))/(1.d0+kap2)**(3./2.)+ &
             kap3*cos(4.d0*(phii+phioff-cotalpha*log(d1/rSK)))/(1.d0+kap3)**(3./2.)
 !--Radial derivative.
 potarmdr= -eta1*(kap1-2.d0)*cos(2.d0*(phii+phioff                      ))/(1.d0+kap1)**(5./2.)/ra1 &
              +(2.*cotalpha*eta2*(kap2+1.d0)*sin(2.d0*(phii+phioff-cotalpha*log(d1/rSK)))              &
              -eta2*(kap2-2.d0)*cos(2.d0*(phii+phioff-cotalpha*log(d1/rSK))))/(1.d0+kap2)**(5./2.)/ra2 &
              +(4.*cotalpha*eta3*(kap3+1.d0)*sin(4.d0*(phii+phioff-cotalpha*log(d1/rSK)))              &
              -eta3*(kap3-2.d0)*cos(4.d0*(phii+phioff-cotalpha*log(d1/rSK))))/(1.d0+kap3)**(5./2.)/ra3
 !--Azimuthal derivative.
 potarmdphi= -kap1*2.d0*sin(2.d0*(phii+phioff                     ))/(1.d0+kap1)**(3./2.) &
                -kap2*2.d0*sin(2.d0*(phii+phioff-cotalpha*log(d1/rSK)))/(1.d0+kap2)**(3./2.) &
                -kap3*4.d0*sin(4.d0*(phii+phioff-cotalpha*log(d1/rSK)))/(1.d0+kap3)**(3./2.)
 !--Combine disc and arm terms.
 pot    = + potdis*(1.d0+eta_arm*potarm)
 fri    = -( (1.d0+eta_arm*potarm)*potdisdr + potdis*eta_arm*potarmdr )
 fphii  = - potdis*eta_arm*potarmdphi/d1
 !--Convert forces in cylindrical to Cartesian co-ordinates.
 fxi    = (fri*cos(phii) - fphii*sin(phii))
 fyi    = (fri*sin(phii) + fphii*cos(phii))
 fzi    = -(1.+eta_arm*potarm)*TANH(zi/zdisk)*pi* sdens_g
 !--Update the input forces/potential
 phi    = phi    + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine KFDiscSp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~BULGES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--Plummer--:A Plummer, one parameter, bulge.
!----------------------------------------------------------------
subroutine PlumBul(xi,yi,zi,r,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: r2term,pot,fxi,fyi,fzi
 !--Add a Plummer sphere
 r2term = sqrt(r*r + Rbulge*Rbulge)
 pot = - gcode*Mbulge/r2term
 fxi = - xi*gcode*Mbulge/(r2term*r2term*r2term)
 fyi = - yi*gcode*Mbulge/(r2term*r2term*r2term)
 fzi = - zi*gcode*Mbulge/(r2term*r2term*r2term)
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine PlumBul

!--Hernaquist--:A Hernaquist, one parameter, bulge.
!----------------------------------------------------------------
subroutine HernBul(xi,yi,zi,r,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: r2term,pot,fxi,fyi,fzi
 !--Add Hernequist bulge
 r2term = ( r + Rbulge )
 pot = - gcode*Mbulge/r2term
 fxi = - xi*gcode*Mbulge/(r2term*r2term*r)
 fyi = - yi*gcode*Mbulge/(r2term*r2term*r)
 fzi = - zi*gcode*Mbulge/(r2term*r2term*r)
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine HernBul

!--HubbleBul--:A Hubble luminosity profile bulge.
!----------------------------------------------------------------
subroutine HubbBul(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,dr,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: rratio,r2term,Hubbforce,pot,fxi,fyi,fzi
 !--Add a Hubble luminosity profile bulge
 rratio    = r/RHubbbulge
 r2term    = sqrt(1.0 + rratio*rratio)
 Hubbforce = MHubbbulge*gcode*HubbCi*(dr*dr*log(rratio+r2term) &
     - dr/(RHubbbulge*r2term))
 pot = - dr*HubbCi*gcode*MHubbbulge*log(r/rHubbbulge + r2term)
 fxi = - xi*dr*Hubbforce
 fyi = - yi*dr*Hubbforce
 fzi = - zi*dr*Hubbforce
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine HubbBul

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~HALOS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--Caldwell&Ostriker--:A logarithmic/arctan DM halo.
!----------------------------------------------------------------
subroutine COhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,dr,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: atanterm,pot,fxi,fyi,fzi
 !--NOTE: in phi we write log(sqrt(x)) = 0.5*log(x)
 atanterm = atan(r*rc21)
 pot = + p1*rc2*(1.+0.5*log(Rtotterm/(1.+(r*rc21)**2))-atanterm*(rc2*dr))
 fxi = - p1*(rc2*dr)**2*(r-rc2*atanterm)*xi*dr
 fyi = - p1*(rc2*dr)**2*(r-rc2*atanterm)*yi*dr
 fzi = - p1*(rc2*dr)**2*(r-rc2*atanterm)*zi*dr
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine COhalo

!--A flat logarithmic DM halo.
!----------------------------------------------------------------
subroutine Flathalo(xi,yi,zi,r,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: pot,fxi,fyi,fzi
 pot = + 0.5*vhalo*vhalo*log(r*r + rLhalo*rLhalo)
 fxi = - vhalo*vhalo*xi/(rLhalo*rLhalo + r*r)
 fyi = - vhalo*vhalo*yi/(rLhalo*rLhalo + r*r)
 fzi = - vhalo*vhalo*zi/(rLhalo*rLhalo + r*r)
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine Flathalo

!--Allen&Martos--: A logarithmic/arctan DM halo
!----------------------------------------------------------------
subroutine AMhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,dr,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: pot,fxi,fyi,fzi,haloterm,Massinhalo,haloforce
 haloterm = ((r/AMrhalo)**1.02)
 Massinhalo= AMmhalo*haloterm*(r/AMrhalo)/(1.0+haloterm)
 pot       =  (- gcode*Massinhalo/r - gcode*AMmhalo/(1.02*AMrhalo)  &
    * ( -1.02/(1.0+(rlim/AMrhalo)**1.02) + log(1.0+(rlim/AMrhalo)**1.02) &
       +1.02/(1.0+haloterm)           - log(1.0+haloterm)          ) )
 haloforce  = -gcode*Massinhalo*dr*dr
 fxi = + xi*dr*haloforce
 fyi = + yi*dr*haloforce
 fzi = + zi*dr*haloforce
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine AMhalo

!--Khoperskov/Begeman--: A normalised logarithmic/arctan DM halo
!----------------------------------------------------------------
subroutine KBhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,dr,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: pot,fxi,fyi,fzi,fri,rratio,atanterm
 rratio=r/rhalo
 atanterm = atan(rratio)

 pot = + Mhalo*gcode*HaloCi*( log(rratio) + atanterm/rratio + &
     0.5*log((1.+rratio**2)/(rratio**2)) )/rhalo
 fri = + dr*dr * Mhalo*gcode*HaloCi *(atanterm -rratio)
 fxi = + xi*dr*fri
 fyi = + yi*dr*fri
 fzi = + zi*dr*fri

 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine KBhalo

!--Navarro/Frenk/White--: DM halo fitted to observations
!----------------------------------------------------------------
subroutine NFWhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: r,dr,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: potent1,fxi,fyi,fzi,fri
 !--Evaluate potential
 potent1= -log(1.d0+r/rnfw)*dr*gcode*Mnfw/(log(1.d0+Cnfw)-Cnfw/(1.d0+Cnfw))

 !--Calculate the force in r, and use to calc Cartesian forces.
 fri    =(log(1.d0+r/rnfw)*dr*dr-1.d0/((rnfw*r)*(1.d0+r/rnfw)))*gcode*Mnfw/(log(1.d0+Cnfw)&
     -Cnfw/(1.d0+Cnfw))
 fxi = -fri*xi*dr
 fyi = -fri*yi*dr
 fzi = -fri*zi*dr

 !--Update potential and forces
 phi    = phi + potent1
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine NFWhalo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~ARMS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--Cox&Gomez--:Altered cosine of TWA
!----------------------------------------------------------------
subroutine s_potential(xi,yi,zi,ti,potout,fextxi,fextyi,fextzi)
 use physcon, only:pi
 real, intent(in)    :: xi,yi,zi,ti
 real, intent(inout) :: potout,fextxi,fextyi,fextzi
 real(kind=8) :: pot,fxi,fyi,fzi,rcyl, phi, rcyl2, dgamdx, dgamdy
 real(kind=8) :: sum, gamma, Kn, Bn, Dn, KnHz, KnHz2, KnHz3, sumg, sumk, sumz
 real(kind=8) :: sechterm,sechtermBn,dsechtermdKn,dsechtermdz
 real(kind=8) :: cdkterm,dcdkdKn,costerm,sinterm,Dnterm,dKnterm,kzonb,term
 real(kind=8) :: costerm2,sinterm2,coshterm,sinhterm,tanhterm
 real(kind=8) :: expx,expxm,As,p0t
 real(kind=8) :: cosin(2,3)
 integer :: n

!--Use the Dehnen00 smooth activation function:
 As = softpot(phir,1.d0,ti)
 p0t =  As * p0

!
!--Convert to cylindrical polars:
!
 phi = atan2(yi,xi)
 rcyl2 = xi*xi + yi*yi
 rcyl  = sqrt(rcyl2)
!
!--Calculate the potential
!
! gamma=NN*(phi+pi/2.-log(r/r0)/tan(alpha))
 gamma  =  NN*(phi+phir*(t0+ti)-log(rcyl/r0)*cotalpha)
 dgamdx = -NN*(yi + xi*cotalpha)/rcyl2
 dgamdy =  NN*(xi - yi*cotalpha)/rcyl2

 sum  = 0.
 sumg = 0.
 sumk = 0.
 sumz = 0.

 !--optimisation: precalculate cos(n*gamma) using
 !  double and triple angle formulae up to n=3
 costerm = cos(gamma)
 sinterm = sin(gamma)
 costerm2 = costerm*costerm
 sinterm2 = sinterm*sinterm
 cosin(1,1) = costerm
 cosin(2,1) = sinterm
 cosin(1,2) = 2.*costerm2 - 1.   ! cos(2*gamma)
 cosin(2,2) = 2.*sinterm*costerm ! sin(2*gamma)
 cosin(1,3) = (4.*costerm2 - 3.)*costerm  ! cos(3*gamma)
 cosin(2,3) = (3. - 4.*sinterm2)*sinterm  ! sin(3*gamma)

 do n=1,3

    Kn      =  n*NN/(rcyl*sinalpha)
    dKnterm = -n*NN/(rcyl*rcyl2*sinalpha)

    KnHz    = Kn*Hz
    KnHz2   = KnHz*KnHz
    KnHz3   = KnHz2*KnHz
    Dnterm  = (1. + KnHz + 0.3*KnHz2)
    Dn      = Dnterm/(1. + 0.3*KnHz)
    Bn      = KnHz*(1.+0.4*KnHz)

    kzonb   = (Kn*zi)/Bn
    expx    = exp(kzonb)
    expxm   = 1./expx !exp(-kzonb)
    coshterm     = 0.5*(expx + expxm) !cosh(kzonb)
    sechterm     = 1./coshterm
    sechtermBn   = sechterm**Bn
    sinhterm     = 0.5*(expx - expxm) !sinh(kzonb)
    tanhterm     = sinhterm*sechterm
    dsechtermdKn = sechtermBn*((0.4*KnHz*Hz + Hz*(1. + 0.4*KnHz))*log(sechterm) &
                              + 0.4*KnHz*zi*tanhterm/(1. + 0.4*KnHz))
    dsechtermdz  = -Kn*sechtermBn*sechterm*sinhterm

    cdkterm  =  Cz(n)/(Dn*Kn)
    dcdkdKn  = -Cz(n)*(1. + 2.*KnHz + 1.2*KnHz2 + 0.18*KnHz3)/(Kn*Dnterm)**2

    !gamman   = n*gamma
    !costerm  = cos(gamman)
    !sinterm  = sin(gamman)
    costerm  = cosin(1,n)
    sinterm  = cosin(2,n)

    sum      = sum  + cdkterm*costerm*sechtermBn
    sumg     = sumg - cdkterm*sinterm*n*sechtermBn
    sumk     = sumk + (dcdkdKn*costerm*sechtermBn + cdkterm*costerm*dsechtermdKn)*dKnterm
    sumz     = sumz + cdkterm*costerm*dsechtermdz

 enddo

 term = -4.*pi*Hz*p0t*exp(-(rcyl-r0)/rS)
 pot   = term*sum

 fxi   = + term*xi/(rS*rcyl)*sum &
         - term*sumg*dgamdx      &
         - term*sumk*xi
 fyi   = + term*yi/(rS*rcyl)*sum &
         - term*sumg*dgamdy      &
         - term*sumk*yi
 fzi   = - term*sumz

 !    logpart=Co*log(Rc**2.+r**2.+(z/0.7)**2.)

 !--Update the input forces/potential
 potout    = potout + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine s_potential

!--Pichardo&Martos--:Schmidt spheroids of rho=p0+p1*a
!----------------------------------------------------------------
subroutine pichardo_potential(xi,yi,zi,rcyl2,ti,phii,phi,fextxi,fextyi,fextzi)
 use physcon, only:pc,solarm,pi,gg,kpc
 real, intent(in)    :: xi,yi,zi,rcyl2,ti,phii
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8) :: rcyl, betain(2)
!--Parameters for the spheroids
 real(kind=8) :: p0i,p1i,om,om2,z2,ys,xs,beta,psi,gArm,fomi,fznew,potnew
 real(kind=8) :: sbeta,cbeta,tbeta,t1,t2,t3,t4,pot,fxi,fyi,fzi
 integer :: ii,n

!
!--Convert to cylindrical polars:
!
 rcyl  = sqrt(rcyl2)
 z2    = zi*zi

 pot =0.d0
 fomi=0.d0
 fzi =0.d0
 fxi =0.d0
 fyi =0.d0
!--Loop over an arm:
 do ii=1,NNi
    !--Loop over all the schmidt spheroids in a single arm:
    n=1
    fomi =0.d0
    fznew=0.d0
    do n=1,Nt
       gArm = shapefn(ii,n) - phir*(t0+ti)
       om2=Rspheroids(ii,n)*Rspheroids(ii,n) + rcyl2 - 2.0d0*Rspheroids(ii,n)*rcyl*&
     cos(gArm-phii)
       om=sqrt(om2)

       !--Are we in or outside the spheroid?
       if ((om2/a_02 + z2/c_02)<=1.0d0) then
          !Inside
          beta = asin(e_0)
          sbeta= e_0
       else
          !Outside
          betain= betafn(om2,z2)
          beta  = betain(2)
          sbeta = betain(1)
       endif

       ys  = Rspheroids(ii,n)*sin(gArm)
       xs  = Rspheroids(ii,n)*cos(gArm)
       psi = atan2((yi-ys),(xi-xs))

       !FOR A PICHARDO SPHEROID ARM STRUCTURE
       !Form: rho = p0 + p1*a
       p0i= den0(ii,n)
       p1i= -p0i/a_0
       cbeta=cos(beta)
       tbeta=TAN(beta)
       t1=sqrt(om2*cbeta*cbeta+z2)
       t2=sqrt(om2+z2)
       t3=(t1+zi)*(t2-zi) /(om2*cbeta)
       t4=log( ((t1+zi)*(t2-zi))/((t1-zi)*(t2+zi)) )

       potnew=  strength *(&
     +a_0*beta*e_0*e_02*((1.d0/3.d0)*p1i*a_02+0.5d0*p0i*a_0)&
     -0.5d0*p0i*e_0*(0.5d0*om2*(beta-sbeta*cbeta) +z2*(tbeta-beta))&
     +(1.d0/9.d0)*p1i*(t2*t2*t2-t1*t1*t1)&
     +t2*((1.d0/6.d0)*p1i*z2-(1.d0/3.d0)*p1i*(om2-z2))&
     +t1*((1.d0/3.d0)*p1i*(om2-z2)-(1.d0/6.d0)*p1i*z2/(cbeta*cbeta))&
     -p1i*zi*(0.5d0*om2-(1.d0/3.d0)*z2)*log(t3)         )

       fomi = strength*om*( 0.5d0*e_0*p0i*(beta-sbeta*cbeta) &
      + p1i*(  +(1.0d0/(3.0d0*om2))*(t1*t1*t1-t2*t2*t2)&
               +t2 -t1 +0.5d0*zi*t4     ))

       fznew = strength*zi*(e_0*p0i*(tbeta-beta) &
      + p1i*(  (1.0d0/(2.0d0*cbeta*cbeta) + 1.0d0)*t1 - 1.5d0*t2 &
       + 0.25d0*((om2-2.0d0*z2)/zi)*t4  ))

       !--Sum forces in x and y (NOT omega), and z.
       pot = pot+potnew
       fxi = fxi+fomi*cos(psi)
       fyi = fyi+fomi*sin(psi)
       fzi = fzi+fznew
       !--Next spheroid...
    enddo
!--Next arm...
 enddo
!--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine pichardo_potential

!--Pichardo&Martos--:Schmidt spheroids of rho=p0+p1/a
!----------------------------------------------------------------
subroutine schmidt_potential(xi,yi,zi,rcyl2,ti,phii,phi,fextxi,fextyi,fextzi)
 use physcon, only:pc,solarm,pi,gg,kpc
 real, intent(in)    :: xi,yi,zi,rcyl2,ti,phii
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8) :: rcyl, betain(2)
!--Parameters for the spheroids
 real(kind=8) :: p0i,p1i,om,om2,z2,ys,xs,beta,psi,gArm,fomi,fznew,potnew
 real(kind=8) :: sbeta,cbeta,tbeta,t1,t2,t3,t4,pot,fxi,fyi,fzi
 integer :: ii,n
!
!--Convert to cylindrical polars:
!
 rcyl  = sqrt(rcyl2)
 z2    = zi*zi

 pot =0.0
 fomi=0.0
 fzi =0.0
 fxi =0.0
 fyi =0.0
!--Loop over an arm:
 do ii=1,NNi
    !--Loop over all the schmidt spheroids in a single arm:
    n=1
    fomi =0.0
    fznew=0.0
    do n=1,Nt
       gArm = shapefn(ii,n) - phir*(t0+ti)
       om2=Rspheroids(ii,n)*Rspheroids(ii,n) + rcyl2 - 2.0*Rspheroids(ii,n)*rcyl*&
       cos(gArm-phii)
       om=sqrt(om2)

       !--Are we in or outside the spheroid?
       if ((om2/a_02 + z2/c_02)<=1.0d0) then
          !Inside
          beta = asin(e_0)
          sbeta= e_0
       else
          !Outside
          betain= betafn(om2,z2)
          beta = betain(2)
          sbeta= betain(1)
       endif

       ys  = Rspheroids(ii,n)*sin(gArm)
       xs  = Rspheroids(ii,n)*cos(gArm)
       psi = atan2((yi-ys),(xi-xs))

       !FOR A POINT MASS ARM STRUCTURE
       !pot = pot + gcode*0.05*Mspiral/(om+0.4*a_0)
       !fomi = -gcode*0.05*Mspiral/((om+0.4*a_0)**2)
       !fznew= 0.0

       !FOR A SCHMIDT SPHEROID ARM STRUCTURE
       !Form: rho = p0 + p1/a
       p0i= den0(ii,n)
       p1i= -p0i*a_0
       cbeta=cos(beta)
       tbeta=TAN(beta)
       t1=sqrt(om2*cbeta*cbeta+z2)
       t2=sqrt(om2+z2)
       t3=(t1+zi)*(t2-zi) /(om2*cbeta)
       t4=log(t3)

       potnew= (strength*e_0)*(a_0*e_02*beta*(0.5d0*p0i*a_0+p1i) &
       - 0.5d0*p0i*(0.5d0*om2*(beta-cbeta*sbeta)+z2*(tbeta-beta)) &
       -p1i*e_0*(t2-t1+zi*t4)   )

       fomi  = (strength*e_0)*om*(0.5d0*p0i*(beta-sbeta*cbeta) +p1i*e_0*(t2-t1)/om2 )
       fznew = (strength*e_0)*zi*(p0i*(tbeta-beta) + (p1i*e_0/zi)*t4 )

       !--Sum forces in x and y (NOT omega), and z.
       pot = pot+potnew
       fxi = fxi+fomi*cos(psi)
       fyi = fyi+fomi*sin(psi)
       fzi = fzi+fznew
       !--Next spheroid...
    enddo
!--Next arm...
 enddo

!--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi

 return
end subroutine schmidt_potential

!----------------------------------------------------------------
function betafn(om2,z2)
 !--Calcs beta using quad equation if particle is inside spheroid.
 !--The difference here compared to betafnold is this only evaluates
 !--the realistic outcomes to save computation time.
 use io, only:fatal
 real(kind=8) :: betafn(2)
 real(kind=8), intent(in) :: om2,z2
 real(kind=8) :: aq,bq,cq,sol1,sol0

 aq = -om2
 bq = om2+z2+a_02*e_02
 cq = -(a_02*e_02)
 sol0= sqrt( bq*bq - 4.0d0*aq*cq )
 sol1= (-bq+sol0)/(2.0d0*aq)

 betafn(1)=sqrt(sol1)  !sine beta
 betafn(2)=asin(betafn(1))  !beta

 return
end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~BARS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--Long&Murali--:Biaxial bar, x-aligned
!----------------------------------------------------------------
subroutine LMXbar(offset,ascale,bscale,xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: offset,ti,xi,yi,zi,ascale,bscale
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: pot,fxi,fyi,fzi,xir,yir,zir,Splus,Sminus,fextxr,fextyr,&
  abar,bbar,Ab,barmasst !--Scale the bar from a default size:
 abar = ascale*LMabar
 bbar = bscale*LMbbar
 !--Rotate our reference frame so the bar is still x-allinged.
 xir = xi*cos(ti*phibar+offset)-yi*sin(ti*phibar+offset)
 yir = xi*sin(ti*phibar+offset)+yi*cos(ti*phibar+offset)
 zir = zi

 !--Use the Dehnen00 smooth activation function:
 Ab = softpot(phibar,1.d0,ti)
 barmasst =  Ab * barmass

 Splus   =sqrt(bbar*bbar + yir*yir+zir*zir + (abar+xir)*(abar+xir))
 Sminus  =sqrt(bbar*bbar + yir*yir+zir*zir + (abar-xir)*(abar-xir))

 pot = + (gcode*barmasst/(2.d0*abar))*log( (xir-abar+Sminus)/(xir+abar+Splus) )
 fextxr =  -2.d0*gcode*barmasst*xir/(Splus*Sminus*(Splus+Sminus))
 fextyr =  -gcode*barmasst*yir*(Sminus+Splus - 4.d0*xir*xir/(Sminus+Splus))&
 / ( 2.d0*Splus*Sminus*(yir*yir+zir*zir+bbar*bbar) )
 fzi =  -gcode*barmasst*zir*(Sminus+Splus - 4.d0*xir*xir/(Sminus+Splus))&
 / ( 2.d0*Splus*Sminus*(yir*yir+zir*zir+bbar*bbar) )

 !--Rotate the reference frame of the forces back to those of the observer,
 !--rather than the bar.
 fxi=+ (fextxr*cos(-ti*phibar-offset)-fextyr*sin(-ti*phibar-offset))
 fyi=+ (fextxr*sin(-ti*phibar-offset)+fextyr*cos(-ti*phibar-offset))

 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine LMXbar

!--Long&Murali--:Triaxial bar, x-aligned
!----------------------------------------------------------------
subroutine LMTbar(xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: xi,yi,zi,ti
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: pot,fxi,fyi,fzi,xir,yir,zir,Tplus,Tminus,fextxr,fextyr,&
   barterm,barterm2,barmasst,Ab

 !--Returns same as above as c->0.
 !--Rotate our reference frame so the bar is still x-allinged.
 xir = xi*cos(ti*phibar)-yi*sin(ti*phibar)
 yir = xi*sin(ti*phibar)+yi*cos(ti*phibar)
 zir = zi

 !--Use the Dehnen00 smooth activation function:
 Ab = softpot(phibar,1.d0,ti)
 barmasst =  Ab * barmass

 barterm = (LMbbar+sqrt(LMcbar*LMcbar+zir*zir))
 barterm2 = barterm*barterm
 Tplus = sqrt((LMabar+xir)*(LMabar+xir)+yir*yir+ barterm2)
 Tminus= sqrt((LMabar-xir)*(LMabar-xir)+yir*yir+ barterm2)
 pot = + gcode*barmasst*log( (xir-LMabar+Tminus)/(xir+LMabar+Tplus) )/(2.d0*LMabar)
 fextxr = -2.d0*gcode*barmasst*xir/(Tplus*Tminus*(Tplus+Tminus))
 fextyr = -gcode*barmasst*yir*(Tminus+Tplus - 4.d0*xir*xir/(Tminus+Tplus)) /&
  ( 2.d0*Tplus*Tminus*(yir*yir+barterm2) )
 fzi    = -gcode*barmasst*zir*(Tminus+Tplus - 4.d0*xir*xir/(Tminus+Tplus)) *&
  ( barterm/(barterm-LMbbar) )  / ( 2.d0*Tplus*Tminus*(yir*yir+barterm2) )

 !--Rotate the reference frame of the forces back to those of the observer,
 !--rather than the bar.
 fxi= + (fextxr*cos(-ti*phibar)-fextyr*sin(-ti*phibar))
 fyi= + (fextxr*sin(-ti*phibar)+fextyr*cos(-ti*phibar))

 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine LMTbar

!--Dehnen--:Constant Quadrupole moment bar, Dehnen00
!----------------------------------------------------------------
subroutine DehnenBar(xi,yi,d2,phii,ti,phi,fextxi,fextyi)
 use units,   only:udist
 use physcon, only:kpc,pi
 real, intent(in)    :: ti,xi,yi,d2,phii
 real, intent(inout) :: phi,fextxi,fextyi
 real(kind=8)  :: bartrig,Cbartrig,Sbartrig,rcubed,pot,fxi,fyi,Ab,t1,eta,d1
 !--Basic quadrupole from Duhnen00

 !--The full bar potential has a smooth switch:
 t1       = 4.*pi/phibar
 eta      = 2.*ti/t1 - 1.
 Ab       = StrBarD*min((0.1875d0*eta**5-0.625d0*eta**3+0.9375d0*eta+0.5d0),1.d0)
 bartrig  = 2.d0*(phii   +  phibar*ti )
 Cbartrig = cos(bartrig)
 Sbartrig = sin(bartrig)
 rcubed=radiusofbar*radiusofbar*radiusofbar

 !--Should all be in terms of x,y [d2], not x,y,z [r**2]
 d1    =sqrt(d2)
 if (d1>radiusofbar) then
    pot = + Ab*Cbartrig * (-rcubed/(d2*d1))
    fxi = - Ab*rcubed *(3.*xi*Cbartrig - yi*2.*Sbartrig)/(d2*d2*d1)
    fyi = - Ab*rcubed *(3.*yi*Cbartrig + xi*2.*Sbartrig)/(d2*d2*d1)
 elseif (d1<radiusofbar .and. d1 > 0.05*kpc/udist) then
    pot = + Ab*Cbartrig*(d2*d1/rcubed - 2.)
    fxi = - Ab*(  2.*yi*Sbartrig*(d2*d1/rcubed-2.)/d2 + 3.*xi*d1*Cbartrig/rcubed )
    fyi = - Ab*( -2.*xi*Sbartrig*(d2*d1/rcubed-2.)/d2 + 3.*yi*d1*Cbartrig/rcubed )
 else
    !Avoid div0 errors if very close to centre
    pot = + Ab*Cbartrig*(d2*d1/rcubed - 2.)
    fxi = 0.
    fyi = 0.
 endif
 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 !--Just a perturbation of a disc,
 !--so this bar has no effect on the z-direction...

 return
end subroutine DehnenBar

!--WadaBar--:similar to Dehnen bar
!----------------------------------------------------------------
subroutine Wadabar(xi,yi,d2,phii,ti,hi,phi,fextxi,fextyi)
 real, intent(in)    :: ti,xi,yi,d2,hi,phii
 real, intent(inout) :: phi,fextxi,fextyi
 real(kind=8)  :: d1,Cbartrig,Sbartrig,fxi,fyi,bartrig,&
   pot,rratio,rterm_inv,dpotdr,dpotdtheta,Ab,StrBarWt

 !--Use the Dehnen00 smooth activation function:
 Ab = softpot(phibar,1.d0,ti)
 StrBarWt =  Ab * StrBarW

 d1 =sqrt(d2)
 bartrig   = 2.0*(phii +  phibar*ti )
 Cbartrig  = cos(bartrig)
 Sbartrig  = sin(bartrig)

 rratio = d1/rcore
 rterm_inv  = 1./(1. + rratio*rratio)

 pot = + StrBarWt * Cbartrig * rratio*rratio *rterm_inv*rterm_inv
 dpotdr     = +StrBarWt * Cbartrig *2.*d1/(rcore*rcore*rcore*rcore)*(rcore*rcore-d2)*rterm_inv*rterm_inv*rterm_inv
 dpotdtheta = -2. * StrBarWt * Sbartrig *  rratio*rratio *rterm_inv*rterm_inv

 fxi = - ( dpotdr*xi/d1 - dpotdtheta*yi/d2 )
 fyi = - ( dpotdr*yi/d1 + dpotdtheta*xi/d2 )

 !--Update the input forces/potential
 phi    = phi + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 return
end subroutine Wadabar

!--Dwek&Wang--:The Dwek G2+G3 bar+bulge.
!----------------------------------------------------------------
subroutine Wang_bar(ri,phii,thetai,pot)
 use physcon, only:pi
 use mathfunc, only:legendre_associated,gegenbauer_poly
 real(kind=8), intent(in)  :: ri,phii,thetai
 real(kind=8), intent(out) :: pot
 !--Bulge+bar, x0=(1,49,0.58,0.4):
 real(kind=8), dimension(20) :: Anlm=(/1.509,-0.086,-0.033,-0.020,-2.606,&
 -0.221,-0.001,0.665,0.129,0.006,6.406,1.295,-0.660,-0.140,0.044,&
 -0.012,-5.859,0.984,-0.030,0.001/)

 real(kind=8) :: thisphi,AlmnSum,s,Plm,Gnl
 integer, dimension(20) :: n_wbar=(/0,1,2,3,0,1,2,0,1,2,0,1,0,1,0,1,0,0,0,0/)
 integer, dimension(20) :: l_wbar=(/0,0,0,0,2,2,2,2,2,2,4,4,4,4,4,4,6,6,6,6/)
 integer, dimension(20) :: m_wbar=(/0,0,0,0,0,0,0,2,2,2,0,0,2,2,4,4,0,2,4,6/)
 real(kind=8), allocatable :: GnlA(:),PlmA(:)
 integer :: i,l,m,n

 !Need to solve via finite differences
 pot=0.
 thisphi=0.
 AlmnSum=0.
 s = ri/rbars

 do i=1,20
    n = n_wbar(i)
    l = l_wbar(i)
    m = m_wbar(i)
    allocate(GnlA(n+1))
    call gegenbauer_poly(n,(2.*real(l,kind=8)+1.5),(s-1.)/(s+1.),GnlA)
    Gnl=GnlA(n+1)

    allocate(PlmA(l+1))
    call legendre_associated(l,m,cos(thetai),PlmA)
    Plm=PlmA(l+1)
    thisphi = Anlm(i) * (s**REAL(l))/((1.+s)**(2.*REAL(l)+1.)) * Gnl * Plm * cos(REAL(m)*(phii))
    AlmnSum = AlmnSum + thisphi

    deallocate(GnlA,PlmA)
 enddo

 pot = +gcode*barmass*AlmnSum/rbars

 return
end subroutine Wang_bar

!--OthgBasis--:Wang spherical harmonic basis set, inc. boxy bulge
!----------------------------------------------------------------
subroutine Orthog_basisbar(xi,yi,zi,r,dr,ti,hi,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: ti,xi,yi,zi,hi,r,dr
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8)  :: potent1,potx,poty,potz,thetai,phii,xit,yit,zit,rit,&
   xir,yir,zir,phiit,thetait,fxi,fyi,fzi,fextxr,fextyr,fextzr,dhi

 !--The r terms are for rotation, t terms for finite difference.
 xir = xi*cos(ti*phibar)-yi*sin(ti*phibar)
 yir = xi*sin(ti*phibar)+yi*cos(ti*phibar)
 zir = zi
 dhi = 0.000001*hi
 !--The potential must be evlulated by finite difference,
 !--this is done over d(x,y,z) = h/1000
 potent1=0.
 potx   =0.
 poty   =0.
 potz   =0.
 phii   = atan2(yir,xir)
 thetai = acos(zir*dr)
 call Wang_bar(real(r,kind=8),phii,thetai,potent1)

 xit    = xir+dhi
 rit    = sqrt(xit*xit+yir*yir+zir*zir)
 phiit  = atan2(yir,xit)
 thetait= acos(zir/rit)
 call Wang_bar(rit,phiit,thetait,potx)

 yit    = yir+dhi
 rit    = sqrt(xir*xir+yit*yit+zir*zir)
 phiit  = atan2(yit,xir)
 thetait= acos(zir/rit)
 call Wang_bar(rit,phiit,thetait,poty)

 zit    = zir+dhi
 rit    = sqrt(xir*xir+yir*yir+zit*zit)
 phiit  = atan2(yir,xir)
 thetait= acos(zit/rit)
 call Wang_bar(rit,phiit,thetait,potz)

 fextxr=-(potx-potent1)/dhi
 fextyr=-(poty-potent1)/dhi
 fextzr=-(potz-potent1)/dhi

 !--Rotate the reference frame of the forces back to those of the observer,
 !--rather than the bar.
 fxi = + (fextxr*cos(-ti*phibar)-fextyr*sin(-ti*phibar))
 fyi = + (fextxr*sin(-ti*phibar)+fextyr*cos(-ti*phibar))
 fzi = + fextzr

 !--Update the input forces/potential
 phi    = phi + potent1
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi

 return
end subroutine Orthog_basisbar

!--VogtS--:A smoothed bar with a S shape pertabation
!----------------------------------------------------------------
subroutine VogtSbar(xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
 real, intent(in)    :: ti,xi,yi,zi
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real(kind=8) :: lamb,Rbent1,Rbent2,Rbent13,Rbent23,barterm1,barterm2,&
   barterm3,pot,fxi,fyi,fzi,fextxr,fextyr,&
   xir,yir,zir
 lamb= barmass/(2.*LMabar)
 !--Rotate our reference frame so the bar is still x-allinged.
 xir = xi*cos(ti*phibar)-yi*sin(ti*phibar)
 yir = xi*sin(ti*phibar)+yi*cos(ti*phibar)
 zir = zi

 Rbent1=sqrt(zir*zir+yir*yir+LMbbar*LMbbar+(xir+LMabar)**2)
 Rbent2=sqrt(zir*zir+yir*yir+LMbbar*LMbbar+(xir-LMabar)**2)
 Rbent13=Rbent1*Rbent1*Rbent1
 Rbent23=Rbent2*Rbent2*Rbent2

 barterm1=zir*zir+yir*yir+xir*xir+LMbbar*LMbbar !z^2+y^2+x^2+b^2
 barterm2=zir*zir+yir*yir+LMbbar*LMbbar         !z^2+y^2+b^2
 barterm3=yir*yir-zir*zir-LMbbar*LMbbar         !y^2-z^2-b^2

 pot = -(gcode*lamb*log( (xir-LMabar+Rbent2)/(xir+LMabar+Rbent1) ) &
  + (ebar*gcode*lamb*yir/(LMabar*LMabar*barterm2*Rbent1*Rbent2)) &
  *(3.*xir*barterm2*Rbent1*Rbent2*log((xir-LMabar+Rbent2)/(xir+LMabar+Rbent1)) &
  +LMabar*xir*(5.*zir*zir+5.*yir*yir+5.*LMbbar*LMbbar-xir*xir)*(Rbent1+Rbent2) &
  - (Rbent1-Rbent2)*(LMabar*LMabar*barterm2 &
  +barterm1*(2.*zir*zir+2.*yir*yir+2.*LMbbar*LMbbar-xir*xir))    ))

 fextxr = +gcode*lamb*(Rbent2-Rbent1)/(Rbent1*Rbent2) -&
  3.*ebar*gcode*lamb*yir*log((xir-LMabar+Rbent2)/(xir+LMabar+Rbent1))/(LMabar*LMabar) -&
  (ebar*gcode*lamb*yir/(LMabar*LMabar*(zir*zir+yir*yir+LMbbar*LMbbar)*Rbent13*Rbent23))*&
   ( Rbent13*(3.*xir*barterm1**2 - &
    LMabar*barterm2*(6.*xir*xir+3.*LMabar*xir-4.*LMabar*LMabar) +&
     3.*LMabar*barterm2**2 -&
      3.*LMabar*xir*xir*(3.*xir*xir-3.*LMabar*xir+LMabar*LMabar))-&
   Rbent23*(3.*xir*barterm1**2 +&
    LMabar*barterm2*(6.*xir*xir-3.*LMabar*xir-4.*LMabar*LMabar) -&
     3.*LMabar*barterm2**2 +&
      3.*LMabar*xir*xir*(3.*xir*xir+3.*LMabar*xir+LMabar*LMabar))  )

 fextyr = +gcode*lamb*((xir-LMabar)*Rbent1-(xir+LMabar)*Rbent2)*yir/((zir*zir+yir*yir+LMbbar*LMbbar)*Rbent1*Rbent2) &
  -3.*ebar*gcode*lamb*xir*log((xir-LMabar+Rbent2)/(xir+LMabar+Rbent1))/(LMabar*LMabar) &
  + (ebar*gcode*lamb/(LMabar*LMabar*Rbent13*Rbent23*(zir*zir+yir*yir+LMbbar*LMbbar)**2))&
   *(+Rbent13*( (barterm1**2)*( 2.*barterm2*(zir*zir+2.*yir*yir+LMbbar*LMbbar)+xir*xir*barterm3 ) &
      -LMabar*xir*xir*xir*barterm3*(LMabar*LMabar-3.*LMabar*xir+3.*xir*xir)  &
       +(LMabar**4)*(barterm2**2)  &
        -3.*LMabar*(barterm2**2) * (xir*(3.*zir*zir+5.*yir*yir+3.*LMbbar*LMbbar)-LMabar*(zir*zir+2.*yir*yir+LMbbar*LMbbar)) &
         -LMabar*xir*barterm2*(-3.*LMabar*xir*(4.*zir*zir+7.*yir*yir+4.*LMbbar*LMbbar) &
           +6.*xir*xir*(zir*zir+3.*yir*yir+LMbbar*LMbbar)+LMabar*LMabar*(7.*zir*zir+10.*yir*yir+7.*LMbbar*LMbbar))  ) &
     -Rbent23*( (barterm1**2)*( 2.*barterm2*(zir*zir+2.*yir*yir+LMbbar*LMbbar)+xir*xir*barterm3 ) &
      +LMabar*xir*xir*xir*barterm3*(LMabar*LMabar+3.*LMabar*xir+3.*xir*xir)  &
       +(LMabar**4)*(barterm2**2)  &
        +3.*LMabar*(barterm2**2) * (xir*(3.*zir*zir+5.*yir*yir+3.*LMbbar*LMbbar)+LMabar*(zir*zir+2.*yir*yir+LMbbar*LMbbar)) &
         +LMabar*xir*barterm2*(+3.*LMabar*xir*(4.*zir*zir+7.*yir*yir+4.*LMbbar*LMbbar) &
           +6.*xir*xir*(zir*zir+3.*yir*yir+LMbbar*LMbbar)+LMabar*LMabar*(7.*zir*zir+10.*yir*yir+7.*LMbbar*LMbbar))  ) )

 fzi=+gcode*lamb*((xir-LMabar)*Rbent1-(xir+LMabar)*Rbent2)*zir/((zir*zir+yir*yir+LMbbar*LMbbar)*Rbent1*Rbent2)

 !!--Rotate the reference frame of the forces back to those of the observer, rather than the bar.
 fxi= +(fextxr*cos(-ti*phibar)-fextyr*sin(-ti*phibar))
 fyi= +(fextxr*sin(-ti*phibar)+fextyr*cos(-ti*phibar))

 !--Update the input forces/potential
 phi    = phi    - pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi

 return
end subroutine VogtSbar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~SOFTENER~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function softpot(pspeed,softfac,ti)
 !--Softens the potentials for bars/spirals to avoid "hard" activation.
 !--Taken from Dehnen00.
 !--The value for softpot reduces the density/mass of the potential.
 use io, only:fatal
 real(kind=8) :: softpot
 real(kind=8), intent(in) :: pspeed,softfac
 real,         intent(in) :: ti
 real(kind=8) :: Tperiod,t1,epsact

 !--rotation period keeping consistent code units:
 Tperiod  = 2.*3.14159 / pspeed
 !--potential activation time, D00 used 2.
 t1       = softfac * Tperiod
 epsact  =  2.*ti/t1 - 1.
 if (ti<t1) then
    softpot = (3./16.*epsact**5 - 5./8.*epsact**3 + 15./16.*epsact + 0.5 )
 else
    softpot = 1.
 endif

 return
end function


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~READIN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine BINReadPot3D(xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
 use physcon, only:pi
 real, intent(in)    :: xi,yi,zi,ti
 real, intent(inout) :: phi,fextxi,fextyi,fextzi
 real    :: fxi,fyi,fzi,xir,yir,zir,fextxr,fextyr,fextzr,&
   pot
 integer :: xindex, yindex, zindex

 !--Rotate our reference frame so the arms are still x-allinged.
 xir = xi*cos(ti*phir)-yi*sin(ti*phir)
 yir = xi*sin(ti*phir)+yi*cos(ti*phir)
 zir = zi

 !--Use the array here, i.e. this routine will put the particles in
 !--an appropriate grid cell and simply read in fextxi, fextyi, fextzi.
 !--The +1 grid cell is because there is no zeroth cell.
 xindex = nint(potlenx*(xir+potxmax)/(2.*potxmax) + 1)
 yindex = nint(potleny*(yir+potymax)/(2.*potymax) + 1)
 zindex = nint(potlenz*(zir+potzmax)/(2.*potzmax) + 1)

 !--Look up the potential at the grid cell, and use Finite Diffs to
 !--find the forces, using the size of grid cells as the distance element.
 pot    = newpot3D(zindex,xindex,yindex)
 fextxr = (pot - newpot3D(zindex,xindex+1,yindex))/dxpot
 fextyr = (pot - newpot3D(zindex,xindex,yindex+1))/dypot
 fextzr = (pot - newpot3D(zindex+1,xindex,yindex))/dzpot

 !--Rotate the reference frame of the forces back to those of the observer,
 !--rather than the arms.
 fxi= + (fextxr*cos(-ti*phir)-fextyr*sin(-ti*phir))
 fyi= + (fextxr*sin(-ti*phir)+fextyr*cos(-ti*phir))
 fzi= + fextzr

 !--Update the input forces/potential
 phi    = phi    + pot
 fextxi = fextxi + fxi
 fextyi = fextyi + fyi
 fextzi = fextzi + fzi
 return
end subroutine BINReadPot3D

!----------------------------------------------------------------
end module extern_spiral
!----------------------------------------------------------------
