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
!  Analysis routine for discs by DF, adapted from a routine by CJN
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, eos, io, kernel, linklist, part, physcon, prompting,
!    units
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'disc_stresses'
 public :: do_analysis, radial_binning, calc_gravitational_forces

 ! Variables for neighbour lists
 integer, parameter :: neighmax = 500
 integer, parameter :: maxcellcache = 50000

 integer, allocatable, dimension(:) :: neighcount, eigenpart
 integer, allocatable, dimension(:,:) :: neighb

 real :: meanneigh, sdneigh, neighcrit
 logical :: neigh_overload

 ! Variables for radial binning
 integer :: nbins
 real :: rin, rout,dr
 integer, allocatable,dimension(:) :: ipartbin

 real,allocatable,dimension(:) :: rad,ninbin,sigma,csbin,vrbin,vphibin, omega
 real,allocatable,dimension(:) :: H,Q, toomre_q,epicyc
 real,allocatable,dimension(:) :: alpha_reyn,alpha_grav,alpha_mag,alpha_art

 real,allocatable,dimension(:) :: rpart,phipart,vrpart,vphipart, gr,gphi,Br,Bphi
 real,allocatable,dimension(:,:) :: gravxyz

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,      only:fatal
 use dim,     only:maxp
 use part,    only:gravity,mhd,Bxyz

 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile

 character(len=9) :: output


 ! Code calculates the following alphas:
 ! Reynolds stress: 2*dvr*dvphi/3*cs^2
 ! Gravitational stress: gr*gphi/(4*pi*cs^2*Grho)
 ! Maxwell stress: 2*-Br*Bphi/3(4*pi*rho*cs^2)
 ! Numerical stress: 0.01*h/H

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'stresses_',numfile
 write(*,'("Output file name is ",A)') output

! Read analysis options
 call read_analysis_options

 if (mhd) print*, 'This is an MHD dump: will calculate Maxwell Stress'

 if (gravity) then
! Calculate gravitational forces for all particles (gradient of the potential)
    call calc_gravitational_forces(dumpfile,npart,xyzh,vxyzu)
 else
    allocate(gravxyz(3,npart))
    gravxyz(:,:) = 0.0
 endif

! Calculate the radial and tangential velocities, forces and fields
 call transform_to_cylindrical(npart,xyzh,vxyzu)

! Bin particles by radius
 call radial_binning(npart,xyzh,vxyzu,pmass)

! Calculate stresses
 call calc_stresses(npart,xyzh,vxyzu,pmass)

! Write out data to file
 call write_radial_data(iunit,output,time)

! End of analysis
 call deallocate_arrays

 print '(a,a)', 'Analysis complete for dump ',dumpfile

end subroutine do_analysis

!-------------------------------------------
!+
! Read options for analysis from file
!+
!-------------------------------------------
subroutine read_analysis_options

 use prompting, only:prompt

 implicit none

 logical :: inputexist
 character(len=21) :: inputfile

! Check for existence of input file
 inputfile = 'disc_stresses.options'
 inquire(file=inputfile, exist=inputexist)

 if (inputexist) then

    print '(a,a,a)', "Parameter file ",inputfile, " found: reading analysis options"

    open(10,file=inputfile, form='formatted')
    read(10,*) nbins
    read(10,*) rin
    read(10,*) rout
    close(10)

 else

    print '(a,a,a)', "Parameter file ",inputfile, " NOT found"

    call prompt('Enter the number of radial bins: ', nbins)
    call prompt('Enter the disc inner radius: ', rin)
    call prompt('Enter the disc outer radius: ', rout)

! Write choices to new inputfile

    open(10,file=inputfile, status='new', form='formatted')
    write(10,*) nbins, "      Number of radial bins"
    write(10,*) rin,  "      Inner Disc Radius"
    write(10,*) rout, "      Outer Disc Radius"
    close(10)
 endif


 print*, 'Inner Disc Radius (code units): ', rin
 print*, 'Outer Disc Radius (code units): ', rout
 print*, 'Number of bins: ', nbins

end subroutine read_analysis_options


!---------------------------------------------------
!+
! Computes gravitational force vector from potential gradients
! Requires neighbours to compute SPH derivatives
!+
!---------------------------------------------------
subroutine calc_gravitational_forces(dumpfile,npart,xyzh,vxyzu)

 use dim, only:gravity,maxp
 use part, only:poten,igas,iphase,maxphase,rhoh,massoftype,get_partinfo
 use kernel, only: get_kernel,get_kernel_grav1,cnormk

 implicit none

 character(len=*),intent(in) :: dumpfile
 real,intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer,intent(in) :: npart


 integer :: j,k,igrav,iamtypei,ipart
 real,dimension(3) :: dr
 real :: rij,rij2, hj1,hj21,hj41,q2i,qi
 real :: rhoi, rhoj, wabi, grkerni, dphidhi, grpmrho1

 logical :: iactivei,iamdusti, iamgasi,existneigh
 character(100) :: neighbourfile

 if (allocated(gravxyz))deallocate(gravxyz)
 allocate(gravxyz(3,npart))

 gravxyz(:,:) = 0.0

! Construct neighbour lists for derivative calculations

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 inquire(file=neighbourfile,exist = existneigh)

 if (existneigh.eqv..true.) then
    print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
    call read_neighbours(neighbourfile,npart)

 else

    ! If there is no neighbour file, generate the list
    print*, 'No neighbour file found: generating'

    call generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile)

 endif

 print*, 'Neighbour lists acquired'
 print*, 'Calculating gravitational force from potential derivative'

 do ipart=1,npart

    rhoi = rhoh(xyzh(4,ipart),massoftype(igas))

    ! Smoothing lengths
    hj1 = 1.0/xyzh(4,ipart)
    hj21 = hj1*hj1
    hj41 = hj21*hj21

    over_neighbours: do k = 1, neighcount(ipart)

       if (k>neighmax) exit over_neighbours
       j = neighb(ipart,k)

       if (maxphase==maxp) then
          call get_partinfo(iphase(j), iactivei,iamdusti,iamtypei)
          iamgasi = (iamtypei ==igas)
       else
          iactivei = .true.
          iamtypei = igas
          iamdusti = .false.
          iamgasi = .true.
       endif

       if (ipart==j) cycle
       if (.not.iamgasi) cycle

       ! Calculate gradient of SPH kernel
       ! Separation of particles

       rij2 = 0.0
       do igrav = 1,3
          dr(igrav) = xyzh(igrav,ipart) - xyzh(igrav,j)
          rij2 = rij2 + dr(igrav)*dr(igrav)
       enddo

       rij = sqrt(rij2)
       dr(:) = dr(:)/rij



       ! r/h
       q2i = rij2*hj21
       qi = rij*hj1

       !--kernel and gradient
       if (gravity) then
          call get_kernel_grav1(q2i,qi,wabi,grkerni,dphidhi)
       else
          call get_kernel(q2i,qi,wabi,grkerni)
       endif

       ! Prefactor to SPH sum = grad*mass/rho
       rhoj = rhoh(xyzh(4,j),massoftype(igas))

       grpmrho1 = grkerni*massoftype(igas)/rhoi

       ! Now calculate each gravitational force vector element
       ! F = - grad phi

       grav_force: do igrav = 1,3
          gravxyz(igrav,ipart) = gravxyz(igrav,ipart) - grpmrho1*hj41*dr(igrav)*(poten(j)-poten(ipart))
       enddo grav_force

       ! End of matrix loop

    enddo over_neighbours
    ! End of loop over nearest neighbours

    gravxyz(:,ipart) = gravxyz(:,ipart)/(massoftype(igas))

 enddo
! End loop over all particles

end subroutine calc_gravitational_forces


!---------------------------------------------------
!+
! Calculates cylindrical transformation for positions, forces and fields
!+
!---------------------------------------------------
subroutine transform_to_cylindrical(npart,xyzh,vxyzu)

 use part, only: mhd,gravity,Bxyz
 implicit none

 integer, intent(in) :: npart
 real,intent(in) ::xyzh(:,:),vxyzu(:,:)
 integer :: ipart

 allocate(rpart(npart))
 allocate(phipart(npart))
 allocate(vrpart(npart))
 allocate(vphipart(npart))
 allocate(gr(npart))
 allocate(gphi(npart))
 allocate(Br(npart))
 allocate(Bphi(npart))

 rpart(:) = 0.0
 phipart(:) = 0.0
 vrpart(:) = 0.0
 vphipart(:) = 0.0
 gr(:) = 0.0
 gphi(:) = 0.0
 Br(:) = 0.0
 Bphi(:) = 0.0

 do ipart=1,npart

    ! Calculate positions first
    rpart(ipart) = sqrt(xyzh(1,ipart)*xyzh(1,ipart) + &
        xyzh(2,ipart)*xyzh(2,ipart))

    phipart(ipart) = atan2(xyzh(2,ipart),xyzh(1,ipart))

    ! Now velocities
    vrpart(ipart) = vxyzu(1,ipart)*cos(phipart(ipart)) + &
        vxyzu(2,ipart)*sin(phipart(ipart))

    vphipart(ipart) = vxyzu(2,ipart)*cos(phipart(ipart)) - &
        vxyzu(1,ipart)*sin(phipart(ipart))

    ! Now gravitational forces
    if (gravity) then
       gr(ipart) = gravxyz(1,ipart)*cos(phipart(ipart)) + &
           gravxyz(2,ipart)*sin(phipart(ipart))

       gphi(ipart) = gravxyz(2,ipart)*cos(phipart(ipart)) - &
           gravxyz(1,ipart)*sin(phipart(ipart))
    endif

    ! Finally, B-Fields
    if (mhd) then
       Br(ipart) = Bxyz(1,ipart)*cos(phipart(ipart)) + &
           Bxyz(2,ipart)*sin(phipart(ipart))
       Bphi(ipart) = Bxyz(2,ipart)*cos(phipart(ipart)) - &
           Bxyz(1,ipart)*sin(phipart(ipart))
    endif

 enddo

end subroutine transform_to_cylindrical

!---------------------------------------------------------------
!+
! Bins particles by radius (and computes some binned quantities)
!+
!---------------------------------------------------------------

subroutine radial_binning(npart,xyzh,vxyzu,pmass)
 use physcon, only:pi
 use eos, only: gamma

 implicit none

 integer,intent(in) :: npart
 real,intent(in) :: pmass
 real,intent(in) :: xyzh(:,:),vxyzu(:,:)

 integer :: ibin,ipart,nbinned
 real :: area

 print '(a,I4)', 'Carrying out radial binning, number of bins: ',nbins

 allocate(ipartbin(npart))
 allocate(rad(nbins))
 allocate(ninbin(nbins))
 allocate(sigma(nbins))
 allocate(csbin(nbins))
 allocate(omega(nbins))
 allocate(vrbin(nbins))
 allocate(vphibin(nbins))

 ipartbin(:) = 0
 ninbin(:) = 0.0
 sigma(:) = 0.0
 csbin(:) = 0.0
 omega(:) = 0.0
 vrbin(:) = 0.0
 vphibin(:) = 0.0

! Set up radial bins

 dr = (rout-rin)/real(nbins-1)
 do ibin=1,nbins
    rad(ibin)=rin + real(ibin-0.5)*dr
 enddo

! Now bin all particles

 nbinned = 0
 do ipart=1,npart

    ! i refers to particle, ii refers to bin
    if (xyzh(4,ipart)  >  tiny(xyzh)) then ! IF ACTIVE

       ibin = int((rpart(ipart)-rad(1))/dr + 1)

       ! Check particle in the binning range - if not skip it
       if (ibin > nbins) ibin=0
       if (ibin < 1)  ibin=0

       if (ibin<=0) cycle

       nbinned = nbinned +1

       ninbin(ibin) = ninbin(ibin) +1
       ipartbin(ipart) = ibin

       csbin(ibin) = csbin(ibin) + sqrt(gamma*(gamma-1)*vxyzu(4,ipart))

       area = pi*((rad(ibin)+0.5*dr)**2-(rad(ibin)- 0.5*dr)**2)
       sigma(ibin) = sigma(ibin) + pmass/area

       vrbin(ibin) = vrbin(ibin) + vrpart(ipart)
       vphibin(ibin) = vphibin(ibin) + vphipart(ipart)
       omega(ibin) = omega(ibin) + vphipart(ipart)/rad(ibin)

    endif

 enddo

 print*, nbinned, ' particles have been binned'

 where(ninbin(:)/=0)
    csbin(:) = csbin(:)/ninbin(:)
    vrbin(:) = vrbin(:)/ninbin(:)
    vphibin(:) = vphibin(:)/ninbin(:)
    omega(:) = omega(:)/ninbin(:)
 end where

 print*, 'Binning Complete'

end subroutine radial_binning


!--------------------------------------------------------------
!+
! Calculates the disc stresses (alphas) and other radial quantities
!+
!--------------------------------------------------------------
subroutine calc_stresses(npart,xyzh,vxyzu,pmass)
 use physcon, only: pi,gg
 use units,   only: print_units, umass,udist,utime,unit_velocity,unit_density,unit_Bfield
 use dim,     only: gravity
 use part,    only: mhd,rhoh,alphaind
 use eos,     only: gamma

 implicit none

 integer, intent(in) :: npart
 real,intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,intent(in) :: pmass

 integer :: ibin,ipart
 real :: cs2, dvr,dvphi,Keplog,rhopart,unit_force

 print*, 'Calculating Disc Stresses'

 allocate(alpha_reyn(nbins))
 allocate(alpha_grav(nbins))
 allocate(alpha_mag(nbins))
 allocate(alpha_art(nbins))
 allocate(toomre_q(nbins))
 allocate(epicyc(nbins))
 allocate(H(nbins))

 alpha_reyn(:) = 0.0
 alpha_grav(:) = 0.0
 alpha_mag(:) = 0.0
 alpha_art(:) = 0.0
 toomre_q(:) = 0.0
 epicyc(:) = 0.0

! Convert data into cgs
 print*, 'Converting data to cgs units'

 call print_units

 sigma(:) = sigma(:)*umass/(udist*udist)
 csbin(:) = csbin(:)*unit_velocity
 omega(:) = omega(:)/utime

 Keplog = 1.5
 unit_force = (unit_velocity/utime)

 if (gravity) then
    gr(:) = gr(:)*unit_force
    gphi(:) = gphi(:)*unit_force
 endif

 if (mhd) then
    Br(:) = Br(:)*unit_Bfield
    Bphi(:) = Bphi(:)*unit_Bfield
 endif

 do ipart=1,npart
    ibin = ipartbin(ipart)

    if (ibin<=0) cycle

    cs2 = gamma*(gamma-1)*vxyzu(4,ipart)*unit_velocity*unit_velocity

    dvr = (vrpart(ipart) - vrbin(ibin))*unit_velocity
    dvphi = (vphipart(ipart) -vphibin(ibin))*unit_velocity
    rhopart = rhoh(xyzh(4,ipart),pmass)*unit_density

    alpha_reyn(ibin) = alpha_reyn(ibin) + dvr*dvphi

    alpha_art(ibin) = alpha_art(ibin) + alphaind(1,ipart)*xyzh(4,ipart)*udist

    if (gravity) alpha_grav(ibin) = alpha_grav(ibin) + gr(ipart)*gphi(ipart)/rhopart

    if (mhd) alpha_mag(ibin) = alpha_mag(ibin) + Br(ipart)*Bphi(ipart)/rhopart

 enddo

! Normalise stresses (by number of particles in each bin)
! Also calculate other properties that rely on binned data
! (e.g. Toomre Q)

 do ibin=1,nbins
    cs2 = csbin(ibin)*csbin(ibin)

    ! calculate epicyclic frequency
    if (ibin>1) then
       epicyc(ibin) = (rad(ibin)*rad(ibin)*omega(ibin)-rad(ibin-1)*rad(ibin-1)*omega(ibin-1))/(rad(ibin)-rad(ibin-1))
    else
       epicyc(ibin) = rad(ibin)*omega(ibin)
    endif

    if (epicyc(ibin)*omega(ibin)>=0.0) then
       epicyc(ibin) = sqrt(2.0*epicyc(ibin)*omega(ibin)/rad(ibin))
    else
       epicyc(ibin) = -sqrt(2.0*abs(epicyc(ibin)*omega(ibin))/rad(ibin))
    endif

    if (sigma(ibin)>0.0) then
       toomre_q(ibin) = csbin(ibin)*omega(ibin)/(pi*gg*sigma(ibin))
    else
       toomre_q(ibin) = 1.0e30
    endif

    if (abs(omega(ibin))>0.0) then
       H(ibin) = csbin(ibin)/omega(ibin)
    else
       H(ibin) = 1.0e30
    endif

    if (ninbin(ibin) >0) then
       alpha_reyn(ibin) =alpha_reyn(ibin)/(Keplog*cs2*ninbin(ibin))
       alpha_art(ibin) = 0.1*alpha_art(ibin)/(ninbin(ibin)*H(ibin))

       if (gravity) alpha_grav(ibin) = alpha_grav(ibin)/(Keplog*4.0*pi*gg*ninbin(ibin)*cs2)
       if (mhd) alpha_mag(ibin) = alpha_mag(ibin)/(Keplog*cs2*ninbin(ibin))

    else
       alpha_reyn(ibin) = 0.0
       alpha_art(ibin) = 0.0
       if (gravity) alpha_grav(ibin) = 0.0
       if (mhd) alpha_mag(ibin) = 0.0
    endif
 enddo

! Convert H back to code units

 H(:) = H(:)/udist

 print*, 'Stresses calculated'
end subroutine calc_stresses

!--------------------------------------------------------------
!+
! Writes radially binned data to file
!+
!--------------------------------------------------------------
subroutine write_radial_data(iunit,output,time)
 implicit none
 integer, intent(in) :: iunit
 real, intent(in) :: time
 character(len=*) :: output
 integer :: ibin

 print '(a,a)', 'Writing to file ',output
 open(iunit,file=output)
 write(iunit,'("# Disc Stress data at t = ",es20.12)') time
 write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'radius (AU)', &
       2,'sigma (cgs)', &
       3,'cs (cgs)', &
       4,'omega (cgs)', &
       5,'epicyc (cgs)',&
       6,'H (AU)',&
       7,'Toomre Q',&
       8,'alpha_reyn',&
       9,'alpha_grav',&
       10,'alpha_mag',&
       11,'alpha_art'

 do ibin=1,nbins
    write(iunit,'(11(es18.10,1X))') rad(ibin),sigma(ibin),csbin(ibin), &
            omega(ibin),epicyc(ibin),H(ibin), abs(toomre_q(ibin)),alpha_reyn(ibin), &
            alpha_grav(ibin),alpha_mag(ibin),alpha_art(ibin)
 enddo

 close(iunit)

 print '(a)', 'File Write complete'

end subroutine write_radial_data

!--------------------------------------------------------
!+
! Deallocate arrays
!+
!-------------------------------------------------------
subroutine deallocate_arrays

 implicit none

 deallocate(gravxyz)
 deallocate(neighcount,neighb)
 deallocate(ipartbin,ninbin,rad)
 deallocate(rpart,phipart,vrpart,vphipart)
 deallocate(gr,gphi,Br,Bphi,vrbin,vphibin)
 deallocate(sigma,csbin,H,toomre_q,omega,epicyc)
 deallocate(alpha_reyn,alpha_grav,alpha_mag,alpha_art)

end subroutine deallocate_arrays

!---------------------------------------------------------
!+
! Generate neighbour lists for all particles
!+
!--------------------------------------------------------
subroutine generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile)

 use dim, only: maxneigh,maxp
 use kernel, only: radkern2
 use linklist, only: ncells, ifirstincell, set_linklist, get_neighbour_list
 use part, only: get_partinfo, igas, iboundary,maxphase, ll, iphase, gravity

 implicit none

 real, intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in) :: npart
 character(len=*), intent(in) :: dumpfile

 real,allocatable,dimension(:,:) :: dumxyzh

 integer :: i, j, iamtypei, icell, ineigh, nneigh
 real :: dx,dy,dz, rij2
 real :: hi1,hj1,hi21,hj21, q2i, q2j

 integer,save :: listneigh(maxneigh)
 real, save:: xyzcache(4,maxcellcache)
 !$omp threadprivate(xyzcache,listneigh)
 real :: fgrav(20)

 logical :: iactivei, iamdusti,iamgasi, ifilledcellcache

 character(100) :: neighbourfile

 !****************************************
 ! 1. Build kdtree and linklist
 ! --> global (shared) neighbour lists for all particles in tree cell
 !****************************************

 print*, 'Building kdtree and linklist: '

 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 print*, '- Done'

 print*, 'Allocating arrays for neighbour storage : '

 allocate(neighcount(npart))
 allocate(neighb(npart,neighmax))

 neighcount(:) = 0
 neighb(:,:) = 0

 print*, '- Done'
 print "(A,I6)", 'Maximum neighbour number allocated:  ', neighmax

 !***************************************
 ! 2. Assign neighbour lists to particles by searching shared list of host cell
 !***************************************

 print*, 'Creating neighbour lists for particles'

 ! Loop over cells

 !$omp parallel default(none) &
 !$omp shared(maxp,maxphase) &
 !$omp shared(ncells,ll,ifirstincell,npart) &
 !$omp shared(xyzh,vxyzu,iphase) &
 !$omp shared(neighcount,neighb) &
 !$omp private(icell,i, j)&
 !$omp private(iamtypei,iamgasi,iamdusti,iactivei) &
 !$omp private(ifilledcellcache,nneigh) &
 !$omp private(hi1,hi21,hj1,hj21,rij2,q2i,q2j) &
 !$omp private(fgrav, dx,dy,dz)
 !$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)

    i = ifirstincell(icell)

    ! Skip empty/inactive cells
    if (i<=0) cycle over_cells

    ! Get neighbour list for the cell

    if (gravity) then
       call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.,f=fgrav)
    else
       call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)
    endif
    ifilledcellcache = .true.

    ! Loop over particles in the cell

    over_parts: do while(i /=0)
       !print*, i, icell, ncells
       if (i<0) then ! i<0 indicates inactive particles
          i = ll(abs(i))
          cycle over_parts
       endif

       if (maxphase==maxp) then
          call get_partinfo(iphase(i), iactivei,iamdusti,iamtypei)
          iamgasi = (iamtypei ==igas)
       else
          iactivei = .true.
          iamtypei = igas
          iamdusti = .false.
          iamgasi = .true.
       endif

       ! Catches case where first particle is inactive
       if (.not.iactivei) then
          i = ll(i)
          cycle over_parts
       endif

       ! do not compute neighbours for boundary particles
       if (iamtypei ==iboundary) cycle over_parts


       ! Fill neighbour list for this particle

       neighcount(i) = 0

       over_neighbours: do ineigh = 1,nneigh
          !print*, i,ineigh, listneigh(ineigh)
          j = abs(listneigh(ineigh))

          ! Skip self-references
          if (i==j) cycle over_neighbours

          dx = xyzh(1,i) - xyzh(1,j)
          dy = xyzh(2,i) - xyzh(2,j)
          dz = xyzh(3,i) - xyzh(3,j)

          rij2 = dx*dx + dy*dy +dz*dz

          hi1 = 1.0/xyzh(4,i)
          hi21 = hi1*hi1

          q2i = rij2*hi21

          hj1 = 1.0/xyzh(4,j)
          hj21 = hj1*hj1
          q2j = rij2*hj21

          is_sph_neighbour: if (q2i < radkern2 .or. q2j < radkern2) then
             !$omp critical
             neighcount(i) = neighcount(i) + 1
             if (neighcount(i) <=neighmax) neighb(i,neighcount(i)) = j
             !$omp end critical
          endif is_sph_neighbour

       enddo over_neighbours
       ! End loop over neighbours

       i = ll(i)
    enddo over_parts
    ! End loop over particles in the cell

 enddo over_cells
 !$omp enddo
 !$omp end parallel

 ! End loop over cells in the kd-tree

 ! Do some simple stats on neighbour numbers

 meanneigh = 0.0
 sdneigh = 0.0
 neighcrit = 0.0

 call neighbours_stats(npart)

 !**************************************
 ! 3. Output neighbour lists to file
 !**************************************

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 call write_neighbours(neighbourfile, npart)

 print*, 'Neighbour finding complete for file ', TRIM(dumpfile)
 deallocate(dumxyzh)
end subroutine generate_neighbour_lists



!--------------------------------------------------------------------
!+
! Calculates the mean and standard deviation of the neighbour number
! Also calculates a 5 sigma deviation from meanneigh
! (This is principally used as a diagnostic aid for structure finding
! algorithms that rely on the nearest neighbours, like CLUMPFIND)
!+
!--------------------------------------------------------------------
subroutine neighbours_stats(npart)

 implicit none
 integer, intent(in) :: npart
 integer :: ipart

 real :: minimum, maximum

 ! Calculate mean and standard deviation of neighbour counts

 maximum = maxval(neighcount)
 minimum = minval(neighcount)
 print*, 'The maximum neighbour count is ', maximum
 print*, 'The minimum neighbour count is ', minimum

 if (maximum > neighmax) then
    print*, 'WARNING! Neighbour count too large for allocated arrays'
 endif

 meanneigh = sum(neighcount)/REAL(npart)
 sdneigh = 0.0

!$omp parallel default(none) &
!$omp shared(neighcount,meanneigh,npart)&
!$omp private(ipart) &
!$omp reduction(+:sdneigh)
!$omp do schedule(runtime)
 do ipart=1,npart
    sdneigh = sdneigh+(neighcount(ipart)-meanneigh)**2
 enddo
 !$omp enddo
 !$omp end parallel

 sdneigh = sqrt(sdneigh/REAL(npart))

 print*, 'Mean neighbour number is ', meanneigh
 print*, 'Standard Deviation: ', sdneigh

 return
end subroutine neighbours_stats

!---------------------------------------
!+
! Reads in a pre-written neighbours file
!+
!---------------------------------------
subroutine read_neighbours(neighbourfile,npart)

 implicit none

 integer, intent(in) :: npart
 character(100), intent(in) ::neighbourfile

 integer :: i,j,neighcheck, tolcheck

 neigh_overload = .false.

 allocate(neighcount(npart))
 allocate(neighb(npart,neighmax))
 neighcount(:) = 0
 neighb(:,:) = 0

 print*, 'Reading neighbour file ', TRIM(neighbourfile)

 open(2, file= neighbourfile,  form = 'UNFORMATTED')

 read(2)  neighcheck, tolcheck, meanneigh,sdneigh,neighcrit

 if (neighcheck/=neighmax) print*, 'WARNING: mismatch in neighmax: ', neighmax, neighcheck

 read(2) (neighcount(i), i=1,npart)
 do i=1,npart

    if (neighcount(i) > neighmax) then
       neigh_overload = .true.
       read(2) (neighb(i,j), j=1,neighmax)
    else
       read(2) (neighb(i,j), j=1,neighcount(i))
    endif

 enddo
 close(2)

 call neighbours_stats(npart)

 if (neigh_overload) then
    print*, 'WARNING! File Read incomplete: neighbour count exceeds array size'
 else
    print*, 'File Read Complete'
 endif

end subroutine read_neighbours



!--------------------------------------------------------------------
!+
! Writes neighbour data to binary file
!+
!--------------------------------------------------------------------
subroutine write_neighbours(neighbourfile,npart)

 implicit none

 integer, intent(in) :: npart

 integer :: i,j
 character(100)::neighbourfile
 ! This is a dummy parameter, used to keep file format similar to other codes
 ! (Will probably delete this later)

 real, parameter :: tolerance = 2.0e0

 neigh_overload = .false.

 neighbourfile = TRIM(neighbourfile)

 print*, 'Writing neighbours to file ', neighbourfile

 OPEN (2, file=neighbourfile, form='unformatted')

 write(2)  neighmax, tolerance, meanneigh,sdneigh,neighcrit
 write(2) (neighcount(i), i=1,npart)
 do i=1,npart
    if (neighcount(i) > neighmax) then
       neigh_overload = .true.
       write(2) (neighb(i,j), j=1,neighmax)
    else
       write(2) (neighb(i,j), j=1,neighcount(i))
    endif
 enddo

 close(2)


 if (neigh_overload) then
    print*, 'WARNING! File write incomplete: neighbour count exceeds array size'
 else
    print*, 'File Write Complete'
 endif

 return
end subroutine write_neighbours



end module analysis
