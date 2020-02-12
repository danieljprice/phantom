!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  transforms dustgrowth dump into multigrain dump for mcfost usage
!
!  REFERENCES: None
!
!  OWNER: Arnaud Vericel
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES:
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust,maxdusttypes,use_dustgrowth
 use part,           only:dustprop,idust,grainsize,graindens,iamtype,iphase,&
                        set_particle_type,isdead_or_accreted,ndusttypes,ndustlarge
 use table_utils,    only:logspace
 use prompting,      only:prompt
 use units,          only:udist
 use io,             only:error
 use timestep,       only:nmax
 use initial_params, only:mdust_in
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                         :: i,j,itype,ndustold,ndustnew,bins_per_dex,nbinmax,nbins,iforce_smax
 real                            :: smaxtmp,smintmp,smin,smax,smax_user
 real                            :: mdustold,mdustnew,code_to_mum,tolm
 real, allocatable, dimension(:) :: grid
 logical                         :: force_smax,file_exists
 character(len=20)               :: infile  = "bin_param.txt", &
                                    outfile = "bin_distrib.dat" 
 
 if ((.not. use_dust) .or. (.not. use_dustgrowth)) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes AND DUSTGROWTH=yes'
    stop
 endif
 
 !- initialise variables
 smaxtmp      = 0.
 smax_user    = 1.
 smintmp      = 1.e26
 bins_per_dex = 5
 ndustold     = 0
 ndustnew     = 0
 mdustold     = 0.
 mdustnew     = 0.
 nbinmax      = 25
 tolm         = 1.e-5
 iforce_smax  = 2
 force_smax   = .false.
 graindens    = maxval(dustprop(2,:)) !- it's actually a constant, we don't care
 code_to_mum  = udist*1.e4
 
 !- set nmax to zero, we just want deriv called once after moddump
 nmax         = 0
 
 !- check if param file exists, created by python script growthtomcfost.py
 inquire(file=infile, exist=file_exists)
 
 !- if not, switch to interactive method
 if (.not.file_exists) then
    call prompt('Set smax manually? (1=yes, 2=no)',iforce_smax,1,2)
    if (iforce_smax == 1) then
       force_smax = .true.
    else
       force_smax = .false.
    endif
    if (force_smax) call prompt('Enter smax in cm',smax_user,0.05)
    call prompt('Enter number of bins per dex',bins_per_dex,1)
 else
    !- file created by phantom/scripts/growthtomcfost.py module
    open (unit=420, file=infile)
    read(420,*) force_smax, smax_user, bins_per_dex
    close(unit=420)
 endif
 
 !- loop over particles, find min and max on non-accreted dust particles
 do i = 1,npart
    itype = iamtype(iphase(i))
    if (itype==idust .and. .not.isdead_or_accreted(xyzh(4,i))) then
       if (dustprop(1,i) < smintmp) smintmp = dustprop(1,i)
       if (dustprop(1,i) > smaxtmp) smaxtmp = dustprop(1,i)
    endif
 enddo
 
 !- force smax if needed
 if (force_smax) then
    smax = smax_user/udist
 else
    smax = smaxtmp
 endif
 smin = smintmp

 !- set ndusttypes based on desired log size spacing
 nbins      = int((log10(smax)-log10(smin))*bins_per_dex + 1.)
 ndusttypes = min(nbins, nbinmax) !- prevent memory allocation errors
 ndustlarge = ndusttypes !- this is written to the header
 
 !- allocate memory for a grid of ndusttypes+1 elements
 allocate(grid(ndusttypes+1))
 
 !- bin sizes in ndusttypes bins
 write(*,"(a,f10.1,a,f10.1,a,i3,a)") "Binning sizes between ",smin*code_to_mum, " (µm) and ",&
                                     smax*code_to_mum," (µm) in ",ndusttypes, " bins"

 call logspace(grid(1:ndusttypes+1),smin,smax)

 !- find representative size for each bin
 do i = 1,ndusttypes
    grainsize(i) = sqrt(grid(i)*grid(i+1))
 enddo

 !- transfer particles from bin to bin depending on their size
 ndustold = npartoftype(idust)
 mdustold = massoftype(idust)*npartoftype(idust) !- initial total mass
 do i=1,npart
    if (iamtype(iphase(i))==idust .and. .not.isdead_or_accreted(xyzh(4,i))) then
       !- figure out which bin
       do j=1,ndusttypes
          if ((dustprop(1,i) >= grid(j)) .and. (dustprop(1,i) < grid(j+1))) then
             if (j > 1) then
                npartoftype(idust+j-1) = npartoftype(idust+j-1) + 1
                npartoftype(idust)     = npartoftype(idust) - 1
                call set_particle_type(i,idust+j-1)
             endif
          endif
          !- if smax has been forced, put larger grains inside last bin
          if ((j==ndusttypes) .and. force_smax .and. (dustprop(1,i) >= grid(j+1))) then
             npartoftype(idust+j-1) = npartoftype(idust+j-1) + 1
             npartoftype(idust)     = npartoftype(idust) - 1
             call set_particle_type(i,idust+j-1)
          endif
       enddo
    endif
 enddo
 
 !- set massoftype for each bin and print info
 open (unit=3693, file=outfile, status="replace")
 write(*,"(a3,a1,a10,a1,a10,a1,a10,a5,a6)") "Bin #","|","s_min","|","s","|","s_max","|==>","npart"
 
 do itype=idust,idust+ndusttypes-1
    write(*,"(i3,a1,f10.1,a1,f10.1,a1,f10.1,a5,i6)") itype-idust+1,"|",grid(itype-idust+1)*code_to_mum,"|", &
                                                     grainsize(itype-idust+1)*code_to_mum, &
                                                     "|",grid(itype-idust+2)*code_to_mum,"|==>",npartoftype(itype)
                                                     
    if (itype > idust) massoftype(itype) = massoftype(idust)
    mdust_in(itype) = massoftype(itype)*npartoftype(itype)
    mdustnew        = mdustnew + mdust_in(itype)
    ndustnew        = ndustnew + npartoftype(itype)
    
    write(3693,*) itype-idust+1,grid(itype-idust+1)*code_to_mum,grainsize(itype-idust+1)*code_to_mum,&
                  grid(itype-idust+2)*code_to_mum,npartoftype(itype)
 enddo
 
 close(unit=3693)
 
 !- sanity check for total number of dust particles
 if (sum(npartoftype(idust:)) /= ndustold) then
    write(*,*) 'ERROR! npartoftype not conserved'
    write(*,*) sum(npartoftype(idust:)), " <-- new vs. old --> ",ndustold
 endif
 
 !- sanity check for total dust mass
 if (abs(mdustold-mdustnew)/mdustold > tolm) then
    write(*,*) 'ERROR! total dust mass not conserved'
    write(*,*) mdustnew, " <-- new vs. old --> ",mdustold
 endif
 
end subroutine modify_dump

end module moddump