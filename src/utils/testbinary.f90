!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program testbin
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: testbin [no arguments]
!
! :Dependencies: datafiles, fileutils, orbits, prompting, setorbit,
!   testbinary
!
 use testbinary, only:test_binary
 use prompting,  only:prompt
 use datafiles,  only:find_phantom_datafile
 use fileutils,  only:get_ncolumns,read_column_labels,find_column,load_data_file
 use setorbit,   only:orbit_t,set_orbit_elements
 use orbits,     only:orbit_is_parabolic
 implicit none
 integer :: j,k,ierr,ierr1,itex,iunit,iu_output
 integer :: nargs,ncolumns,nheader,nlabels
 integer :: icol_dv(3),icol_dx(3),icol_m1,icol_m2,icol_chisq
 real :: a,e,inc,o,w,f,m1,m2,chisq,rp
 character(len=120) :: filename
 character(len=32), allocatable :: labels(:) ! column labels
 real, allocatable :: array(:,:)
 type(orbit_t) :: orbit

 ! get filename from command line
 nargs = command_argument_count()
 if (nargs > 0) then
    call get_command_argument(1,filename)
 else
    filename = find_phantom_datafile('orbits.txt','orbits')
 endif

 itex = 0
 open(newunit=iunit,file=filename,status='old',iostat=ierr)
 if (ierr /= 0) then
    print*,' ERROR: could not open file ',trim(filename)
    stop
 endif
 call get_ncolumns(iunit,ncolumns,nheader)

 if (ncolumns == 8) then
    ! read orbits from orbits.txt file or similar
    if (ierr==0) then
       open(newunit=itex,file='orbits.tex',status='replace',iostat=ierr1)
       print "(a)",' writing to orbits.tex'
    endif
    j = 0
    do while(ierr==0)
       j = j + 1
       read(iunit,*,iostat=ierr) a,e,inc,o,w,f,m1,m2
       if (ierr==0) call test_binary(m1,m2,a,e,inc,o,w,f,j,itex)
    enddo
 else
    ! process orbits from optimisation runs, convert output from
    ! the Orbit Reconstructor^TM to orbital elements
    allocate(labels(ncolumns))
    icol_dv = 0; icol_dx = 0; icol_m1 = 0; icol_m2 = 0; icol_chisq = 0
    call read_column_labels(iunit,nheader,ncolumns,nlabels,labels)
    icol_dv(1) = find_column(labels,'binary_dvx',verbose=.true.)
    if (icol_dv(1) > 0) then
       orbit%input_type = 2
       orbit%obs%dx(:) = [' 346.','-242.','   0.']
       orbit%obs%dv(:) = [' 0.040','-0.067',' 0.000']
       orbit%flyby%d = '1200.0'
       m1 = 1.0; m2 = 0.8
       chisq = 0.
       icol_dv(2) = find_column(labels,'binary_dvy',verbose=.true.)
       icol_dv(3) = find_column(labels,'binary_dvz',verbose=.true.)
       icol_dx(1) = find_column(labels,'binary_dx',verbose=.true.)
       icol_dx(2) = find_column(labels,'binary_dy',verbose=.true.)
       icol_dx(3) = find_column(labels,'binary_dz',verbose=.true.)
       icol_m1 = find_column(labels,'m1',verbose=.true.)
       icol_m2 = find_column(labels,'m2',verbose=.true.)
       icol_chisq = find_column(labels,'chisq',verbose=.true.)
    endif
    close(iunit)
    call load_data_file(filename,array,nheader)
    print*,' loaded ',size(array,1),' rows of data'
    print*,' first row of data: ',array(1,:)
    if (allocated(array) .and. icol_dv(1) > 0) then
       open(newunit=iu_output,file='orbits.csv',status='new',action='write')
       write(iu_output,"(a)") 'iteration,rp,e,i,O,w,f,m1,m2,chisq'
       do j=1,size(array,1)
          print*,' processing row ',j,' dvx=',array(j,icol_dv(1))
          if (icol_m1 > 0) m1 = array(j,icol_m1)
          if (icol_m2 > 0) m2 = array(j,icol_m2)
          do k=1,3
             if (icol_dv(k) > 0) write(orbit%obs%dv(k),"(g12.4)") array(j,icol_dv(k))
          enddo
          do k=1,3
             if (icol_dx(k) > 0) write(orbit%obs%dx(k),"(g12.4)") array(j,icol_dx(k))
          enddo
          if (icol_chisq > 0) chisq = array(j,icol_chisq)
          call set_orbit_elements(orbit,m1,m2,verbose=.true.)
          if (orbit_is_parabolic(orbit%e)) then
             rp = orbit%a
          else
             rp = orbit%a*(1.-orbit%e)
          endif
          write(iu_output,"(i0,',',8(g0,','),g0)") j,rp,orbit%e,orbit%i,orbit%O,orbit%w,orbit%f,m1,m2,chisq
       enddo
       close(iu_output)
    endif
    ! print*,' got ',nlabels,' labels: ',labels(1:nlabels)
 endif

 close(iunit)
 if (itex > 0) close(itex)

end program testbin
