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
!  Analysis routine which computes the velocity shear tensor for all particles
!
!  Sij = 0.5*(di vj + dj vi)
!
!  and writes the eigenvalues and eigenvectors to file
!  REFERENCES: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, kernel, linklist, part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'velocityshear'

 integer, parameter :: neighmax = 350
 integer, parameter :: maxcellcache = 50000

 integer,parameter :: it_max = 10

 integer, allocatable, dimension(:) :: neighcount, eigenpart
 integer, allocatable, dimension(:,:) :: neighb

 real,allocatable,dimension(:) :: xbin,ybin,zbin
 real,allocatable,dimension(:,:) :: eigenvalues
 real,allocatable,dimension(:,:,:) :: eigenvectors

 real :: meanneigh, sdneigh, neighcrit
 logical :: neigh_overload

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use dim, only: maxp
 use part, only: iphase,maxphase, igas, get_partinfo

 implicit none

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 integer :: i, iwrite, iamtypei
 integer :: it_num, rot_num, ndigits


 real, dimension(3) :: eigen
 real, dimension(3,3) :: eigenvec, tensor

 logical :: iactivei, iamdusti,iamgasi
 logical :: existneigh

 character(100) :: neighbourfile, valuefile, vectorfile
 character(10) :: fmtstring,numstring

 !***************************************
 !1. Obtain Neighbour Lists
 !****************************************

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
 print*, 'Calculating velocity shear tensor and computing eigenvalues'

 allocate(xbin(npart), ybin(npart), zbin(npart))
 allocate(eigenpart(npart))
 allocate(eigenvalues(3,npart))
 allocate(eigenvectors(3,3,npart))

 iwrite = 0
 over_parts: do i=1,npart

    ! Only calculate for gas particles

    if(maxphase==maxp) then
       call get_partinfo(iphase(i), iactivei,iamdusti,iamtypei)
       iamgasi = (iamtypei ==igas)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi = .true.
    endif

    if(.not.iamgasi) cycle over_parts

    iwrite = iwrite +1

    !************************************************
    ! 2. Calculate the velocity shear tensor for all particles
    !***********************************************

    tensor(:,:) = 0.0
    call calc_velocitysheartensor(i,tensor, xyzh,vxyzu)

    !print*, tensor(:,:)
    !***********************************************
    ! 3. Calculate the eigenvalues and eigenvectors
    !***********************************************

    it_num = 0
    rot_num = 0

    call jacobi_eigenvalue(3,tensor,it_max,eigenvec,eigen,it_num,rot_num)
    !print*, eigen


    ! Store data in arrays for later writing to file

    xbin(iwrite) = xyzh(1,i)
    ybin(iwrite) = xyzh(2,i)
    zbin(iwrite) = xyzh(3,i)
    eigenpart(iwrite) = i
    eigenvalues(:,iwrite) = eigen(:)
    eigenvectors(:,:,iwrite) = eigenvec(:,:)

 enddo over_parts
 print*, '- Done'
! End of loop over particles

!*********************************************
! 4. Write to eigenvalue and eigenvector files
!*********************************************

 if(num==0) then
    ndigits = 1
 else
    ndigits = int(log10(real(num)))+1
 endif

 print*, num, ndigits
 write(fmtstring, "('(I',I1,')')") ndigits

 write(numstring, fmtstring) num
 valuefile = 'eig0'//TRIM(numstring)
 vectorfile = 'evc0'//TRIM(numstring)


 call write_eigenfiles(valuefile,vectorfile, iwrite)

 ! Deallocate memory for next dump
 deallocate(neighcount,neighb)

end subroutine do_analysis


!---------------------------------------------------------
!+
! Generate neighbour lists for all particles
!+
!--------------------------------------------------------
subroutine generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile)

 use dim, only: maxneigh,maxp
 use kernel, only: radkern2
 use linklist, only: ncells, ifirstincell, set_linklist, get_neighbour_list
 use part, only: get_partinfo, igas, iboundary,maxphase, ll, iphase

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
 print "(A,I3)", 'Maximum neighbour number allocated:  ', neighmax

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

#ifdef GRAVITY
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.,f=fgrav)
#else
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)
#endif
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


!-----------------------------------------------------
!+
! Calculates the velocity shear tensor for particle ipart
!+
!-----------------------------------------------------

subroutine calc_velocitysheartensor(ipart,tensor, xyzh,vxyzu)

 use dim, only: gravity, maxp
 use kernel, only: get_kernel, get_kernel_grav1
 use part, only: igas, iphase, maxphase, rhoh, massoftype, get_partinfo

 integer, intent(in) :: ipart
 real, intent(in) :: xyzh(:,:), vxyzu(:,:)
 real, dimension(3,3), intent(out) :: tensor

 integer :: j,k, imat, jmat, iamtypei

 real :: rij,rij2, hj1,hj21,hj41,q2i,qi
 real :: rhoj, wabi, grkerni, dphidhi, grpmrho1

 logical :: iactivei,iamdusti, iamgasi

 real, dimension(3) :: dr

 over_neighbours: do k = 1, neighcount(ipart)

    if (k>neighmax) exit

    j = neighb(ipart,k)

    if(maxphase==maxp) then
       call get_partinfo(iphase(j), iactivei,iamdusti,iamtypei)
       iamgasi = (iamtypei ==igas)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi = .true.
    endif

    if(ipart==j) cycle
    if(.not.iamgasi) cycle

    ! Calculate gradient of SPH kernel
    ! Separation of particles

    rij2 = 0.0
    do jmat = 1,3
       dr(jmat) = xyzh(jmat,ipart) - xyzh(jmat,j)
       rij2 = rij2 + dr(jmat)*dr(jmat)
    enddo

    rij = sqrt(rij2)

    ! Smoothing lengths
    hj1 = 1.0/xyzh(4,j)
    hj21 = hj1*hj1
    hj41 = hj21*hj21

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
    grpmrho1 = grkerni*massoftype(igas)/rhoj

    ! Loops for matrix elements i,j
    tensor_calc: do imat = 1,3
       do jmat = 1,3
          tensor(imat,jmat) = tensor(imat,jmat) - &
    grpmrho1*hj41*(dr(imat)*(vxyzu(jmat,j)-vxyzu(jmat,ipart)) + &
    dr(jmat)*(vxyzu(imat,j)-vxyzu(imat,ipart)))
       enddo
    enddo tensor_calc
    ! End of matrix loop

 enddo over_neighbours
 ! End of loop over nearest neighbours

 tensor(:,:) = 0.5*tensor(:,:)

 return
end subroutine calc_velocitysheartensor


!----------------------------------------------------
!+
! Calculate eigenvalues and eigenvectors of a tensor
! using the classical Jacobi rotation method
! (authored by John Burkardt,
! distributed under the GNU LGPL license)
!+
!----------------------------------------------------

subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
!
!    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
!
!    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
!
 implicit none

 integer ( kind = 4 ) n

 real ( kind = 8 ) a(n,n)
 real ( kind = 8 ) bw(n)
 real ( kind = 8 ) c
 real ( kind = 8 ) d(n)
 real ( kind = 8 ) g
 real ( kind = 8 ) gapq
 real ( kind = 8 ) h
 integer ( kind = 4 ) i
 integer ( kind = 4 ) it_max
 integer ( kind = 4 ) it_num
 integer ( kind = 4 ) j
 integer ( kind = 4 ) k
 integer ( kind = 4 ) l
 integer ( kind = 4 ) m
 integer ( kind = 4 ) p
 integer ( kind = 4 ) q
 integer ( kind = 4 ) rot_num
 real ( kind = 8 ) s
 real ( kind = 8 ) t
 real ( kind = 8 ) tau
 real ( kind = 8 ) term
 real ( kind = 8 ) termp
 real ( kind = 8 ) termq
 real ( kind = 8 ) theta
 real ( kind = 8 ) thresh
 real ( kind = 8 ) v(n,n)
 real ( kind = 8 ) w(n)
 real ( kind = 8 ) zw(n)

 do j = 1, n
    do i = 1, n
       v(i,j) = 0.0D+00
    enddo
    v(j,j) = 1.0D+00
 enddo

 do i = 1, n
    d(i) = a(i,i)
 enddo

 bw(1:n) = d(1:n)
 zw(1:n) = 0.0D+00
 it_num = 0
 rot_num = 0

 do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
       do i = 1, j - 1
          thresh = thresh + a(i,j) ** 2
       enddo
    enddo

    thresh = sqrt ( thresh ) / real ( 4 * n, kind = 8 )

    if ( thresh == 0.0D+00 ) then
       exit
    endif

    do p = 1, n
       do q = p + 1, n

          gapq = 10.0D+00 * abs ( a(p,q) )
          termp = gapq + abs ( d(p) )
          termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
          if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

             a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
          else if ( thresh <= abs ( a(p,q) ) ) then

             h = d(q) - d(p)
             term = abs ( h ) + gapq

             if ( term == abs ( h ) ) then
                t = a(p,q) / h
             else
                theta = 0.5D+00 * h / a(p,q)
                t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
                if ( theta < 0.0D+00 ) then
                   t = - t
                endif
             endif

             c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
             s = t * c
             tau = s / ( 1.0D+00 + c )
             h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
             zw(p) = zw(p) - h
             zw(q) = zw(q) + h
             d(p) = d(p) - h
             d(q) = d(q) + h

             a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
             do j = 1, p - 1
                g = a(j,p)
                h = a(j,q)
                a(j,p) = g - s * ( h + g * tau )
                a(j,q) = h + s * ( g - h * tau )
             enddo

             do j = p + 1, q - 1
                g = a(p,j)
                h = a(j,q)
                a(p,j) = g - s * ( h + g * tau )
                a(j,q) = h + s * ( g - h * tau )
             enddo

             do j = q + 1, n
                g = a(p,j)
                h = a(q,j)
                a(p,j) = g - s * ( h + g * tau )
                a(q,j) = h + s * ( g - h * tau )
             enddo
!
!  Accumulate information in the eigenvector matrix.
!
             do j = 1, n
                g = v(j,p)
                h = v(j,q)
                v(j,p) = g - s * ( h + g * tau )
                v(j,q) = h + s * ( g - h * tau )
             enddo

             rot_num = rot_num + 1

          endif

       enddo
    enddo

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

 enddo
!
!  Restore upper triangle of input matrix.
!
 do j = 1, n
    do i = 1, j - 1
       a(i,j) = a(j,i)
    enddo
 enddo
!
!  Ascending sort the eigenvalues and eigenvectors.
!
 do k = 1, n - 1

    m = k

    do l = k + 1, n
       if ( d(l) < d(m) ) then
          m = l
       endif
    enddo

    if ( m /= k ) then

       t    = d(m)
       d(m) = d(k)
       d(k) = t

       w(1:n)   = v(1:n,m)
       v(1:n,m) = v(1:n,k)
       v(1:n,k) = w(1:n)

    endif

 enddo

 return
end subroutine jacobi_eigenvalue



!----------------------------------------------------
!+
! Write eigenvalue and eigenvector data to file
!+
!----------------------------------------------------
subroutine write_eigenfiles(valuefile,vectorfile, ngas)

 integer, intent(in) :: ngas
 integer :: i
 character(100) :: valuefile, vectorfile


! Write eigenvalues to file
 print*, 'Writing eigenvalues to file ', TRIM(valuefile)
 open(27,file=TRIM(valuefile), status='unknown',form='unformatted')
 write(27) ngas
 write(27) (eigenpart(i),i=1,ngas)
 write(27) (xbin(i), i=1,ngas)
 write(27) (ybin(i), i=1,ngas)
 write(27) (zbin(i), i=1,ngas)
 write(27) (eigenvalues(1,i), i=1,ngas)
 write(27) (eigenvalues(2,i), i=1,ngas)
 write(27) (eigenvalues(3,i), i=1,ngas)
 close(27)

! Now write the eigenvectors to file
 print*, 'Writing eigenvectors to file ', TRIM(vectorfile)
 open(27,file=TRIM(vectorfile),status='unknown', form='unformatted')
 write(27) ngas
 write(27) (eigenpart(i),i=1,ngas)
 write(27) (eigenvectors(1,1:3,i),i=1,ngas)
 write(27) (eigenvectors(2,1:3,i),i=1,ngas)
 write(27) (eigenvectors(3,1:3,i),i=1,ngas)


end subroutine write_eigenfiles

end module
