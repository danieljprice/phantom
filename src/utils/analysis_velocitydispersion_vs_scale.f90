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
!  Analysis routine calculates the velocity dispersion as a function of size scale
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
!    sortutils
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=27), parameter, public :: analysistype = 'velocitydispersion_vs_scale'
 public :: do_analysis

 ! Variables for neighbour lists
 integer, parameter :: neighmax = 500
 integer, parameter :: maxcellcache = 50000

 integer, allocatable, dimension(:) :: neighcount
 integer, allocatable, dimension(:,:) :: neighb

 real :: meanneigh, sdneigh, neighcrit
 logical :: neigh_overload

! Variables for particle list construction
 integer :: ninside
 integer, allocatable,dimension(:) :: particlelist, checked, rhosort

! Variables for mean velocity and velocity dispersion

 real,allocatable,dimension(:,:) :: vmean, vdisp
 real, allocatable,dimension(:) :: rhopart
 real, allocatable,dimension(:) :: ekin, egrav,etherm, emag


 ! Variables for scales
 integer :: nscale
 real :: rscale, volume_scale,rscalemin, rscalemax,drscale

 ! Variables for outputs
 real, dimension(3) :: vdisp_exp, vdisp_var
 real :: ekin_exp, etherm_exp, egrav_exp, emag_exp
 real :: ekin_var, etherm_var, egrav_var, emag_var

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,      only:fatal
 use dim,     only:maxp
 use part,    only:gravity,mhd,Bxyz,rhoh,igas,&
              get_partinfo,maxphase,maxp,iphase,massoftype,poten
 use eos,     only: utherm
 use physcon, only:pi
 use sortutils, only: indexx

 implicit none

 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile

 character(len=11) :: output

 character(100) :: neighbourfile

 integer :: j,k,l, ipart,jpart,iscale,ncalculated,iamtypei
 real :: percent, percentcount,vmeansum,vdispsum,dv
 real :: pmassi, rhoj,rhoj1
 logical :: existneigh, iactivei,iamdusti,iamgasi

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a6,i5.5)") 'scale_',numfile
 write(*,'("Output file name is ",A)') output

 call write_output_header(iunit,output,time)

 ! Read analysis options
 call read_analysis_options

 !***************************************
 !1. Obtain Neighbour Lists
 !****************************************

 print '(a)', 'Obtaining Neighbour Lists'

 ! Check if a neighbour file is present

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 inquire(file=neighbourfile,exist = existneigh)

 if(existneigh.eqv..true.) then
    print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
    call read_neighbours(neighbourfile,npart)

 else

    ! If there is no neighbour file, generate the list
    print*, 'No neighbour file found: generating'

    call generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile)

 endif

 print*, 'Neighbour lists acquired'

 ! Sort particles into descending density order

 allocate(rhopart(npart))
 rhopart(:) = 0.0

 do ipart=1,npart
    if(maxphase==maxp) then
       call get_partinfo(iphase(ipart), iactivei,iamdusti,iamtypei)
       iamgasi = (iamtypei ==igas)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi = .true.
    endif

    if(.not.iamgasi) cycle
    rhopart(ipart) = rhoh(xyzh(4,ipart), massoftype(igas))
 enddo

 allocate(rhosort(npart))
 call indexx(npart,rhopart, rhosort)

 print*, 'Particles sorted into descending density order'

! Allocate array to track particles that have already participated in the calculation
 allocate(checked(npart))
 allocate(particlelist(npart))

 allocate(vmean(3,npart))
 allocate(vdisp(3,npart))

 allocate(ekin(npart))
 allocate(etherm(npart))
 allocate(egrav(npart))
 allocate(emag(npart))

 vmean(:,:) = 0.0
 vdisp(:,:) = 0.0

 ekin(:) = 0.0
 etherm(:) = 0.0
 egrav(:) = 0.0
 emag(:) = 0.0

 pmassi = massoftype(igas)

 ! Loop over size scales
 do iscale=1,nscale

    rscale = rscalemin + (iscale-1)*drscale

    volume_scale = 4.0*pi*rscale*rscale*rscale/3.0

    checked(:) = 0

    write(*,'(a,f5.2)') 'Calculating velocity dispersion at size scale ',rscale
    ! Loop over particles (in density order)

    percentcount = 10.0
    do j=1,npart

       percent = 100.0*real(j)/real(npart)
       if(percent> percentcount) then
          write(*,'(f5.1,a)') percentcount, '% complete'
          percentcount = percentcount + 10.0
       endif

       ipart = rhosort(npart-j+1)

       ! If particle has zero density, skip
       if(rhopart(ipart) < 1.0e-20) cycle

       ! If the size scale is less than 2h, skip
       if(rscale < 2.0*xyzh(4,ipart)) cycle

       ! Find all particles within a distance rscale
       ! That haven't already been found (checked=1)

       !print*, 'Finding particles in range'
       call find_particles_in_range(ipart,npart,xyzh,particlelist,rscale)
       !print*, 'Done'

       ! if number of particles in range less than mean neighbour number, then skip
       if(ninside < meanneigh) cycle

       ! calculate mean velocity in this volume (in x-direction)
       ! Calculate energies - kinetic, gravitational, thermal, magnetic

       vmean(:,ipart) = 0.0
       vdisp(:,ipart) = 0.0
       ekin(ipart) = 0.0
       egrav(ipart) = 0.0
       etherm(ipart) = 0.0
       emag(ipart) = 0.0

       do l = 1,ninside
          jpart = particlelist(l)

          rhoj = rhoh(xyzh(4,jpart),pmassi)
          rhoj1 = 1.0/rhoj

          do k=1,3
             vmean(k,ipart) = vmean(k,ipart) + vxyzu(k,jpart)
          enddo

          ekin(ipart) = ekin(ipart) + 0.5*pmassi*dot_product(vxyzu(:,jpart), vxyzu(:,jpart))
          etherm(ipart) = etherm(ipart) + pmassi*utherm(vxyzu(4,jpart),rhoj)
          if(gravity) egrav(ipart) = egrav(ipart) + pmassi*poten(jpart)
          if(mhd) emag(ipart) = emag(ipart) + pmassi*dot_product(Bxyz(:,jpart), Bxyz(:,jpart))*rhoj1

       enddo


       ! Divide by number inside to get mean
       vmean(:,ipart) = vmean(:,ipart)/real(ninside)

       ! Calculate velocity dispersion
       do l = 1,ninside
          jpart = particlelist(l)

          do k=1,3
             dv = vxyzu(k,ipart)-vmean(k,ipart)
             vdisp(k,ipart)  = vdisp(k,ipart) + dv*dv
          enddo
       enddo

       do k=1,3
          vdisp(k,ipart) = sqrt(vdisp(k,ipart)/real(ninside-1))
       enddo

    enddo

    ! End of loop over particles

    ekin(:) = ekin(:)/volume_scale
    etherm(:) = etherm(:)/volume_scale
    if(gravity) egrav(:) = egrav(:)/volume_scale
    if(mhd) emag(:) = emag(:)/volume_scale

    ! Now calculate expectation value, variance of velocity dispersion
    ! Expectation and variance of energies as well TODO

    vdisp_exp(:) = 0.0
    ekin_exp = 0.0
    etherm_exp = 0.0
    emag_exp = 0.0
    egrav_exp = 0.0


    ncalculated = 0

    do ipart=1,npart

       vmeansum = sum(vmean(:,ipart))
       vdispsum = sum(vdisp(:,ipart))
       if(vmeansum<1.0e-10 .and. vdispsum < 1.0e-10) cycle
       ncalculated = ncalculated + 1

       do k=1,3
          vdisp_exp(k) = vdisp_exp(k) + vdisp(k,ipart)
       enddo

       ekin_exp = ekin_exp + ekin(ipart)
       etherm_exp = etherm_exp + etherm(ipart)
       if(gravity) egrav_exp = egrav_exp + egrav(ipart)
       if(mhd) emag_exp = emag_exp + emag(ipart)

    enddo

    vdisp_exp(:) = vdisp_exp(:)/real(ncalculated)
    ekin_exp = ekin_exp/real(ncalculated)
    etherm_exp = etherm_exp/real(ncalculated)
    if(gravity) egrav_exp = egrav_exp/real(ncalculated)
    if(mhd) emag_exp = emag_exp/real(ncalculated)

    vdisp_var(:) = 0.0
    ekin_var = 0.0
    etherm_var = 0.0
    egrav_var = 0.0
    emag_var = 0.0

    ncalculated = 0

    do ipart=1,npart

       vmeansum = sum(vmean(:,ipart))
       vdispsum = sum(vdisp(:,ipart))
       if(vmeansum<1.0e-10 .and. vdispsum < 1.0e-10) cycle

       ncalculated = ncalculated+1
       do k=1,3
          dv = vdisp_exp(k)-vdisp(k,ipart)
          vdisp_var(k) = vdisp_var(k) + dv*dv
       enddo

       ekin_var = ekin_var + (ekin_exp - ekin(ipart))*(ekin_exp - ekin(ipart))
       etherm_var = etherm_var + (etherm_exp - etherm(ipart))*(etherm_exp - etherm(ipart))
       if(gravity) egrav_var = egrav_var + (egrav_exp - egrav(ipart))*(egrav_exp - egrav(ipart))
       if(mhd) emag_var = emag_var + (emag_exp - emag(ipart))*(emag_exp - emag(ipart))

    enddo

    vdisp_var(:) = sqrt(vdisp_var(:)/real(ncalculated-1))
    ekin_var = sqrt(ekin_var/real(ncalculated-1))
    etherm_var = sqrt(etherm_var/real(ncalculated-1))
    if(gravity) egrav_var = sqrt(egrav_var/real(ncalculated-1))
    if(mhd) emag_var = sqrt(emag_var/real(ncalculated-1))

    ! Write this data to file
    call write_output_data(iunit,output)

 enddo


 close(iunit)
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
 character(len=35) :: inputfile

! Check for existence of input file
 inputfile = 'velocitydispersion_vs_scale.options'
 inquire(file=inputfile, exist=inputexist)

 if(inputexist) then

    print '(a,a,a)', "Parameter file ",inputfile, " found: reading analysis options"

    open(10,file=inputfile, form='formatted')
    read(10,*) nscale
    read(10,*) rscalemin
    read(10,*) rscalemax
    close(10)

    drscale = (rscalemax-rscalemin)/real(nscale)

 else

    print '(a,a,a)', "Parameter file ",inputfile, " NOT found"

    call prompt('Enter the number of scale evaluations: ', nscale)
    call prompt('Enter the minimum scale (code units): ', rscalemin)
    call prompt('Enter the maximum scale (code units): ', rscalemax)

! Write choices to new inputfile

    open(10,file=inputfile, status='new', form='formatted')
    write(10,*) nscale, "      Number of scale evaluations"
    write(10,*) rscalemin,  "      Minimum scale (code units)"
    write(10,*) rscalemax, "      Maximum scale (code units)"
    close(10)
 endif


 print*, 'Minimum scale (code units): ', rscalemin
 print*, 'Maximum scale (code units): ', rscalemax
 print*, 'Number of Evaluations ', nscale

end subroutine read_analysis_options

!--------------------------------------------------------------
!+
! Find all particles within a range rscale of the selected particle
! Uses a linked list search across the neighbour list
!+
!-------------------------------------------------------------

subroutine find_particles_in_range(ipart,npart,xyzh,particlelist,d)
 implicit none

 integer, intent(in) :: ipart,npart
 real, intent(in) :: d
 real, intent(in) :: xyzh(:,:)
 integer,intent(inout) :: particlelist(:)

 real,parameter :: tolerance = 2.0

 integer, allocatable, dimension(:) :: teststack

 integer :: jpart,l,nstack
 real :: sep

! Create a stack of particle IDs to test
 allocate(teststack(npart))

 teststack(:) = 0
 nstack = 0

! List of particles within range rscale
 particlelist(:) = 0
 ninside = 0

! Add neighbours to test stack

 do l=1,neighcount(ipart)
    nstack = nstack+1
    teststack(nstack) = neighb(ipart,l)
 enddo

 particlelist(:) = 0
 ninside = 0

! Continue testing particles in the stack until the end is reached
 jpart = 1
 do while(jpart < nstack .and. nstack < npart)


    if(checked(jpart)==1) cycle ! Skip particles that have already been calculated
    if(rhopart(jpart) < 1.0e-20) cycle ! Skip particles that have zero density

    ! Calculate separation
    sep = (xyzh(1,ipart)-xyzh(1,jpart))*(xyzh(1,ipart)-xyzh(1,jpart)) + &
        (xyzh(2,ipart)-xyzh(2,jpart))*(xyzh(2,ipart)-xyzh(2,jpart)) + &
        (xyzh(3,ipart)-xyzh(3,jpart))*(xyzh(3,ipart)-xyzh(3,jpart))
    sep = sqrt(sep)

    ! If particle closer than rscale, add it to the particle list
    if(sep < d) then
       ninside = ninside +1
       particlelist(ninside) = jpart
    endif

    ! If particle within the tolerance distance, then add its neighbours to the stack for testing
    if(sep < tolerance*d) then

       do l=1,neighcount(jpart)
          if(nstack==npart) exit
          nstack = nstack+1
          teststack(nstack) = neighb(jpart,l)
       enddo
    endif

    jpart = jpart+1

 enddo

 deallocate(teststack)

end subroutine


!--------------------------------------------------------------
!+
! Writes header of output file
!+
!--------------------------------------------------------------
subroutine write_output_header(iunit,output,time)
 implicit none
 integer, intent(in) :: iunit
 real, intent(in) :: time
 character(len=*) :: output


 print '(a,a)', 'Writing to file ',output
 open(iunit,file=output)
 write(iunit,'("# Velocity Dispersion data at t = ",es20.12)') time
 write(iunit,"('#',15(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'size scale (pc)', &
       2,'E(vel disp x.) ', &
       3,'V(vel disp x.)', &
       4,'E(vel disp y.) ', &
       5,'V(vel disp y.)', &
       6,'E(vel disp z.) ', &
       7,'V(vel disp z.)', &
       8,'E(kinetic E)', &
       9,'V(kinetic E)', &
       10,'E(thermal E)', &
       11,'V(thermal E)', &
       12,'E(gravity E)', &
       13,'V(gravity E)', &
       14,'E(magnetic E)', &
       15,'V(magnetic E)'


 print '(a)', 'File Header Written'
 call flush(iunit)

end subroutine write_output_header

!-----------------------------------------------------------------
!+
! Write a line of output to file
!+
!-----------------------------------------------------------------
subroutine write_output_data(iunit,output)
 implicit none
 integer, intent(in) :: iunit
 character(len=*) :: output


 write(iunit,'(15(es18.10,1X))') rscale, &
     vdisp_exp(1), vdisp_var(1), &
     vdisp_exp(2), vdisp_var(2), &
     vdisp_exp(3), vdisp_var(3), &
     ekin_exp, ekin_var, etherm_exp, etherm_var, &
     egrav_exp, egrav_var, emag_exp, emag_var
 call flush(iunit)

end subroutine write_output_data


!--------------------------------------------------------
!+
! Deallocate arrays
!+
!-------------------------------------------------------
subroutine deallocate_arrays

 implicit none

 deallocate(neighcount,neighb)
 deallocate(particlelist,checked,rhosort,rhopart)
 deallocate(vmean,vdisp)
 deallocate(ekin,egrav,etherm,emag)

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
    if(i<=0) cycle over_cells

    ! Get neighbour list for the cell

    if(gravity) then
       call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.,f=fgrav)
    else
       call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)
    endif
    ifilledcellcache = .true.

    ! Loop over particles in the cell

    over_parts: do while(i /=0)
       !print*, i, icell, ncells
       if(i<0) then ! i<0 indicates inactive particles
          i = ll(abs(i))
          cycle over_parts
       endif

       if(maxphase==maxp) then
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
       if(iamtypei ==iboundary) cycle over_parts


       ! Fill neighbour list for this particle

       neighcount(i) = 0

       over_neighbours: do ineigh = 1,nneigh
          !print*, i,ineigh, listneigh(ineigh)
          j = abs(listneigh(ineigh))

          ! Skip self-references
          if(i==j) cycle over_neighbours

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

          is_sph_neighbour: if(q2i < radkern2 .or. q2j < radkern2) then
             !$omp critical
             neighcount(i) = neighcount(i) + 1
             if(neighcount(i) <=neighmax) neighb(i,neighcount(i)) = j
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

 if(maximum > neighmax) then
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
 neighcrit = meanneigh-5.0*sdneigh

 print*, 'Clumps created only if neighbour number greater than ', neighcrit

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

 if(neighcheck/=neighmax) print*, 'WARNING: mismatch in neighmax: ', neighmax, neighcheck

 read(2) (neighcount(i), i=1,npart)
 do i=1,npart

    if(neighcount(i) > neighmax) then
       neigh_overload = .true.
       read(2) (neighb(i,j), j=1,neighmax)
    else
       read(2) (neighb(i,j), j=1,neighcount(i))
    endif

 enddo
 close(2)

 call neighbours_stats(npart)

 if(neigh_overload) then
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
    if(neighcount(i) > neighmax) then
       neigh_overload = .true.
       write(2) (neighb(i,j), j=1,neighmax)
    else
       write(2) (neighb(i,j), j=1,neighcount(i))
    endif
 enddo

 close(2)


 if(neigh_overload) then
    print*, 'WARNING! File write incomplete: neighbour count exceeds array size'
 else
    print*, 'File Write Complete'
 endif

 return
end subroutine write_neighbours



end module analysis
