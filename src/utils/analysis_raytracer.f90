!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which computes neighbour lists for all particles
!
! :References: None
!
! :Owner: Mats Esseldeurs
!
! :Runtime parameters: None
!
! :Dependencies: raytracer, part, dump_utils, dust_formation, linklist
!
   use raytracer_copy,   only:get_all_tau_inwards, get_all_tau_adaptive
   use raytracer,        only:get_all_tau
   use part,             only:rhoh,isdead_or_accreted,nsinkproperties,iReff
   use dump_utils,       only:read_array_from_file!, read_header, dump_h, print_header, open_dumpfile_r, lenid, read_block_header
   use getneighbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                              neighcount,neighb,neighmax
   use dust_formation,   only:calc_kappa_bowen
   use physcon,          only:kboltz,mass_proton_cgs,au,solarm
   use linklist, only:set_linklist,allocate_linklist,deallocate_linklist

   implicit none

   character(len=20), parameter, public :: analysistype = 'raytracer'
   real :: gamma = 1.2
   real :: mu = 2.381
   public :: do_analysis

   private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
   use omp_lib

   character(len=*), intent(in) :: dumpfile
   integer,          intent(in) :: num,npart,iunit
   real(kind=8),     intent(in) :: xyzh(:,:),vxyzu(:,:)
   real(kind=8),     intent(in) :: particlemass,time

   logical :: existneigh
   character(100) :: neighbourfile
   character(100)   :: jstring, kstring
   real             :: primsec(4,2), rho(npart), kappa(npart), temp(npart), u(npart), &
        xyzh2(4,npart), vxyzu2(4,npart), xyzmh_ptmass(nsinkproperties,2)
   real, dimension(:), allocatable :: tau
   integer :: i,j,k,ierr,iu1,iu2,iu3,iu4, npart2!,iu
   integer :: start, finish, method, analyses, minOrder, maxOrder, order
   real :: totalTime, timeTau, Rstar, Rcomp
   logical :: SPH, calcInwards

   real, parameter :: udist = au, umass = solarm

   Rstar = 2.37686663
   Rcomp = 0.1
   xyzmh_ptmass = 0.
   xyzmh_ptmass(iReff,1) = Rstar
   xyzmh_ptmass(iReff,2) = Rcomp

   print*,'("Reading kappa from file")'
   call read_array_from_file(123,dumpfile,'kappa',kappa(:),ierr, 1)
   if (ierr/=0) then
      print*,''
      print*,'("WARNING: could not read kappa from file. It will be set to zero")'
      print*,''
      kappa = 0.
   endif

   if (kappa(1) <= 0. .and. kappa(2) <= 0. .and. kappa(2) <= 0.) then
      print*,'("Reading temperature from file")'
      call read_array_from_file(123,dumpfile,'temperature',temp(:),ierr, 1)
      if (temp(1) <= 0. .and. temp(2) <= 0. .and. temp(2) <= 0.) then
         print*,'("Reading internal energy from file")'
         call read_array_from_file(123,dumpfile,'u',u(:),ierr, 1)
         do i=1,npart
            temp(i)=(gamma-1.)*mu*u(i)*mass_proton_cgs*kboltz
         enddo
      endif
      do i=1,npart
         kappa(i)=calc_kappa_bowen(temp(i))
      enddo
   endif

   j=1
   do i=1,npart
      if (.not.isdead_or_accreted(xyzh(4,i))) then
         xyzh2(:,j) = xyzh(:,i)
         vxyzu2(:,j) = vxyzu(:,i)
         kappa(j) = kappa(i)
         j=j+1
      endif
   enddo
   npart2 = j-1
   call set_linklist(npart2,npart2,xyzh2,vxyzu)
   print*,'npart = ',npart2
   allocate(tau(npart2))

   !get position of sink particles (stars)
   call read_array_from_file(123,dumpfile,'x',primsec(1,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'y',primsec(2,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'z',primsec(3,:),ierr, 2)
   call read_array_from_file(123,dumpfile,'h',primsec(4,:),ierr, 2)
   if (primsec(1,1) == xyzh(1,1) .and.  primsec(2,1) == xyzh(2,1) .and. primsec(3,1) == xyzh(3,1)) then
      primsec(:,1) = (/0.,0.,0.,1./)
   endif
   xyzmh_ptmass(1:4,1) = primsec(:,1)
   xyzmh_ptmass(1:4,2) = primsec(:,2)


   print *,'What do you want to an analyses?'
   print *, '(1) Analysis'
   print *, '(2) Method'
   print *, '(3) Calculate tau as in phantom'
   print *, '(3) Preloaded settings'
   read *,analyses
   ! analyses=4

   if (analyses == 1) then
      print *,'Which analysis would you like to run?'
      print *, '(1) Outwards Integration'
      print *, '(2) Adaptive (Outward) Integration'
      print *, '(3) Scaling'
      print *, '(4) Time evoluation for mutiple files'
      read *,method
      if (method == 1) then
         SPH = .true.
         calcInwards = .false.
         print *,'At which order would you like to start?'
         read *,minOrder
         print *,'At which order would you like to stop?'
         read *,maxOrder
      else if (method == 2) then
         SPH = .true.
         calcInwards = .false.
         print *,'At which order would you like to start?'
         read *,minOrder
         print *,'At which order would you like to stop?'
         read *,maxOrder
      else if (method == 3) then

      endif
   else if (analyses == 2) then
      print *,'Which algorithm would you like to run?'
      print *, '(1) Inwards'
      print *, '(2) Outwards'
      print *, '(3) Adaptive'
      read *,method
      if (method == 1) then
         print *,'Do you want to use SPH neighbours? (T/F)'
         read*,SPH
      else if (method == 2) then
         print *,'What order do you want to run?'
         read*,j
         write(jstring,'(i0)') j
      else if (method == 3) then
         print *,'What order do you want to run? (integer below 7)'
         read*,j
         write(jstring,'(i0)') j
         print *,'What refinement level do you want to run? (integer below 7)'
         read*,k
         write(kstring,'(i0)') k
      endif
   endif

   if ((analyses == 1 .and. calcInwards) .or. (analyses == 2 .and. method==1)) then ! get neighbours
      if (SPH) then
         neighbourfile = 'neigh_'//TRIM(dumpfile)
         inquire(file=neighbourfile,exist = existneigh)
         if (existneigh) then
            print*, 'Neighbour file ', TRIM(neighbourfile), ' found'
            call read_neighbours(neighbourfile,npart2)
         else
            ! If there is no neighbour file, generate the list
            print*, 'No neighbour file found: generating'
            call system_clock(start)
            call generate_neighbour_lists(xyzh2,vxyzu2,npart2,dumpfile, .false.)
            call system_clock(finish)
            totalTime = (finish-start)/1000.
            print*,'Time = ',totalTime,' seconds.'
            call write_neighbours(neighbourfile, npart2)
            print*, 'Neighbour finding complete for file ', TRIM(dumpfile)
         endif
      else
         allocate(neighb(npart2+2,100))
         neighb = 0
         inquire(file='neighbors_tess.txt',exist = existneigh)
         if (existneigh) then
            print*, 'Neighbour file neighbors.txt found'
         else
            call execute_command_line('python3 getNeigh.py -f '//'points_'//dumpfile//'.txt')
         endif
         open(newunit=iu4, file='neighbors_tess.txt', status='old', action='read')
         do i=1, npart2+2
               read(iu4,*) neighb(i,:)
         enddo
         close(iu4)
      endif
   endif

   if (analyses == 1) then

      ! OUTWARD INTEGRATION SCHEME
      if (method == 1) then
         if (calcInwards) then
            print*,''
            print*, 'Start calculating optical depth inwards'
            if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
               call system_clock(start)
               call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, kappa, Rstar, tau)
               call system_clock(finish)
            else
               call system_clock(start)
               call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, kappa, Rstar, tau, primsec(:,2),Rcomp)
               call system_clock(finish)
            endif
            timeTau = (finish-start)/1000.
            print*,'Time = ',timeTau,' seconds.'
            open(newunit=iu4, file='times_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) timeTau
            close(iu4)
            totalTime = timeTau
            open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
            do i=1, size(tau)
               write(iu2, *) tau(i)
            enddo
            close(iu2)
            deallocate(neighb)
         else
            open(newunit=iu4, file='times_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) 0.
            close(iu4)
            totalTime=0
         endif

         do j = minOrder, maxOrder
            write(jstring,'(i0)') j
            print*,''
            print*, 'Start calculating optical depth outwards: ', trim(jstring)
            if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
               call system_clock(start)
               call get_all_tau(npart2, 1, xyzmh_ptmass, xyzh2, kappa, j, tau)
               call system_clock(finish)
            else
               call system_clock(start)
               call get_all_tau(npart2, 2, xyzmh_ptmass, xyzh2, kappa, j, tau)
               call system_clock(finish)
            endif
            timeTau = (finish-start)/1000.
            print*,'Time = ',timeTau,' seconds.'
            open(newunit=iu4, file='times_'//dumpfile//'.txt',position='append', status='old', action='write')
            write(iu4, *) timeTau
            close(iu4)
            totalTime = totalTime + timeTau
            open(newunit=iu2, file='taus_'//dumpfile//'_'//trim(jstring)//'.txt', status='replace', action='write')
            do i=1, size(tau)
               write(iu2, *) tau(i)
            enddo
            close(iu2)
         enddo
         print*,''
         print*,'Total time of the calculation = ',totalTime,' seconds.'

      !ADAPTIVE (OUTWARD) INTEGRATION SCHEME
      else if (method == 2) then
         if (calcInwards) then
            print*,''
            print*, 'Start calculating optical depth inwards'
            if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
               call system_clock(start)
               call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, kappa, Rstar, tau)
               call system_clock(finish)
            else
               call system_clock(start)
               call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, kappa, Rstar, tau, primsec(:,2),Rcomp)
               call system_clock(finish)
            endif
            timeTau = (finish-start)/1000.
            print*,'Time = ',timeTau,' seconds.'
            open(newunit=iu4, file='times_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) timeTau
            close(iu4)
            totalTime = timeTau
            open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
            do i=1, size(tau)
               write(iu2, *) tau(i)
            enddo
            close(iu2)
            deallocate(neighb)
         else
            open(newunit=iu4, file='times_opt_'//dumpfile//'.txt', status='replace', action='write')
               write(iu4, *) 0.
            close(iu4)
            totalTime=0
         endif

         do j = minOrder, maxOrder
            write(jstring,'(i0)') j
            do k = minOrder,maxOrder-j
               write(kstring,'(i0)') k
               print*,''
               print*, 'Start calculating optical depth outwards: minOrder = ', trim(jstring),', refineLevel = ', trim(kstring)
               if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
                  call system_clock(start)
                  call get_all_tau_adaptive(primsec(:,1), xyzh2, kappa, Rstar, j, k, tau)
                  call system_clock(finish)
               else
                  call system_clock(start)
                  call get_all_tau_adaptive(primsec(:,1), xyzh2, kappa, Rstar, j, k, tau, primsec(:,2),Rcomp)
                  call system_clock(finish)
               endif
               timeTau = (finish-start)/1000.
               print*,'Time = ',timeTau,' seconds.'
               open(newunit=iu4, file='times_opt_'//dumpfile//'.txt',position='append', status='old', action='write')
               write(iu4, *) timeTau
               close(iu4)
               totalTime = totalTime + timeTau
               open(newunit=iu2, file='taus_'//dumpfile//'_opt_'//trim(jstring)// &
                     '_'//trim(kstring)//'.txt', status='replace', action='write')
               do i=1, size(tau)
                  write(iu2, *) tau(i)
               enddo
               close(iu2)
            enddo
         enddo
         print*,''
         print*,'Total time of the calculation = ',totalTime,' seconds.'
      else if (method == 3) then
         order = 5
         print*,'Start doing scaling analysis with order =',order
         open(newunit=iu4, file='times_'//dumpfile//'_scaling.txt', status='replace', action='write')
         close(iu4)
         do i=1, omp_get_max_threads()
            call omp_set_num_threads(i)
            call deallocate_linklist
            call allocate_linklist
            call set_linklist(npart2,npart2,xyzh2,vxyzu)
            if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
               call system_clock(start)
               call get_all_tau(npart2, 1, xyzmh_ptmass, xyzh2, kappa, order, tau)
               call system_clock(finish)
            else
               call system_clock(start)
               call get_all_tau(npart2, 2, xyzmh_ptmass, xyzh2, kappa, order, tau)
               call system_clock(finish)
            endif
            timeTau = (finish-start)/1000.
            print*,'nthread = ',omp_get_max_threads(),': Time = ',timeTau,' seconds.'
            open(newunit=iu4, file='times_'//dumpfile//'_scaling.txt',position='append', status='old', action='write')
            write(iu4, *) omp_get_max_threads(), timeTau
            close(iu4)
         enddo
      else if (method == 4) then
            order = 5
            print*,'Start doing scaling analysis with order =',order
            if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
               call system_clock(start)
               call get_all_tau(npart2, 1, xyzmh_ptmass, xyzh2, kappa, order, tau)
               call system_clock(finish)
            else
               call system_clock(start)
               call get_all_tau(npart2, 2, xyzmh_ptmass, xyzh2, kappa, order, tau)
               call system_clock(finish)
            endif
            timeTau = (finish-start)/1000.
            print*,'Time = ',timeTau,' seconds.'
            open(newunit=iu1, file='npart_wind.txt',position='append', action='write')
            write(iu1, *) npart2
            close(iu1)
            open(newunit=iu4, file='times_wind.txt',position='append', action='write')
            write(iu4, *) timeTau
            close(iu4)
            totalTime = totalTime + timeTau
            open(newunit=iu2, file='taus_'//dumpfile//'.txt', status='replace', action='write')
            do i=1, size(tau)
               write(iu2, *) tau(i)
            enddo
            close(iu2)
      endif

   else if (analyses == 2) then
      !ADAPTIVE (OUTWARD) INTEGRATION SCHEME
      if (method == 1) then
         print*,''
         print*, 'Start calculating optical depth inwards'
         if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
            call system_clock(start)
            call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, kappa, Rstar, tau)
            call system_clock(finish)
         else
            call system_clock(start)
            call get_all_tau_inwards(primsec(:,1), xyzh2, neighb, kappa, Rstar, tau, primsec(:,2),Rcomp)
            call system_clock(finish)
         endif
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         if (SPH) then
            open(newunit=iu2, file='taus_'//dumpfile//'_inwards.txt', status='replace', action='write')
         else
            open(newunit=iu2, file='taus_'//dumpfile//'_tess_inwards.txt', status='replace', action='write')
         endif
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      else if (method == 2) then
         print*,''
         print*, 'Start calculating optical depth outwards: ', trim(jstring)
         if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
            call system_clock(start)
            call get_all_tau(npart2, 1, xyzmh_ptmass, xyzh2, kappa, j, tau)
            call system_clock(finish)
         else
            call system_clock(start)
            call get_all_tau(npart2, 2, xyzmh_ptmass, xyzh2, kappa, j, tau)
            call system_clock(finish)
         endif
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         open(newunit=iu2, file='taus_'//dumpfile//'_'//trim(jstring)//'.txt', status='replace', action='write')
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      else if (method == 3) then
         print*,''
         print*, 'Start calculating optical depth adaptive: minOrder = ', trim(jstring),', refineLevel = ', trim(kstring)
         if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
            call system_clock(start)
            call get_all_tau_adaptive(primsec(:,1), xyzh2, kappa, Rstar, j, k, tau)
            call system_clock(finish)
         else
            call system_clock(start)
            call get_all_tau_adaptive(primsec(:,1), xyzh2, kappa, Rstar, j, k, tau, primsec(:,2),Rcomp)
            call system_clock(finish)
         endif
         timeTau = (finish-start)/1000.
         print*,'Time = ',timeTau,' seconds.'
         totalTime = totalTime + timeTau
         open(newunit=iu2, file='taus_'//dumpfile//'_opt_'//trim(jstring)// &
               '_'//trim(kstring)//'.txt', status='replace', action='write')
         do i=1, size(tau)
            write(iu2, *) tau(i)
         enddo
         close(iu2)
      endif

   else if (analyses == 3) then
      order = 4
      print*,'Start calculating optical depth'
      if (primsec(1,2) == 0. .and. primsec(2,2) == 0. .and. primsec(3,2) == 0.) then
         call system_clock(start)
         call get_all_tau(npart2, 1, xyzmh_ptmass, xyzh2, kappa, order, tau)
         call system_clock(finish)
      else
         call system_clock(start)
         call get_all_tau(npart2, 2, xyzmh_ptmass, xyzh2, kappa, order, tau)
         call system_clock(finish)
      endif
      timeTau = (finish-start)/1000.
      print*,'Time = ',timeTau,' seconds.'
      open(newunit=iu, file='taus_'//dumpfile//'.txt', status='replace', action='write')
      do i=1, size(tau)
         write(iu, *) tau(i)
      enddo
      close(iu)

   else if (analyses == 4) then
      do i=1,npart
         if (norm2(xyzh2(1:3,i) - (/10.,10.,10./)) < 4.) then
            kappa(i) = 1e10
         endif
      enddo
      ! allocate(neighb(npart2+2,100))
      ! neighb = 0
      ! open(newunit=iu4, file='neighbors_tess.txt', status='old', action='read')
      ! do i=1, npart2+2
      !       read(iu4,*) neighb(i,:)
      ! enddo
      ! close(iu4)
      print*,''
      order = 7
      print*, 'Start calculating optical depth outwards, order=',order
      call system_clock(start)
      call get_all_tau(npart2, 2, xyzmh_ptmass, xyzh2, kappa, order, tau)
      call system_clock(finish)
      timeTau = (finish-start)/1000.
      print*,'Time = ',timeTau,' seconds.'
      totalTime = totalTime + timeTau
      open(newunit=iu2, file='taus_'//dumpfile//'_raypolation_7.txt', status='replace', action='write')
      do i=1, size(tau)
         write(iu2, *) tau(i)
      enddo
      close(iu2)
   endif

   if (.false.) then! Write out the location and density of all the points
      open(newunit=iu1, file='points_'//dumpfile//'.txt', status='replace', action='write')
      do i=1, npart2+2
         write(iu1, *) xyzh2(1:3,i)
      enddo
      close(iu1)

      open(newunit=iu3, file='rho_'//dumpfile//'.txt', status='replace', action='write')
      do i=1,npart2
         rho(i) = rhoh(xyzh2(4,i), particlemass)
         write(iu3, *) rho(i)
      enddo
      close(iu3)
   endif

end subroutine do_analysis
end module analysis