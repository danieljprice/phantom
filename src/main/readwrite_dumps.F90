!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_dumps
!
! readwrite_dumps
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, io, readwrite_dumps_fortran,
!   readwrite_dumps_hdf5
!

 use readwrite_dumps_fortran
#ifdef HDF5
 use readwrite_dumps_hdf5
#endif

 implicit none

 public :: write_smalldump,write_fulldump,read_smalldump,read_dump,write_gadgetdump

#ifdef AOCC
 logical, pointer, public    :: opened_full_dump
 logical, pointer, public    :: dt_read_in
#else
 logical, pointer, public    :: opened_full_dump => opened_full_dump_fortran      ! for use in analysis files if user wishes to skip small dumps
 logical, pointer, public    :: dt_read_in => dt_read_in_fortran           ! to determine if dt has been read in so that ibin & ibinold can be set on restarts
#endif

 integer, parameter, public :: is_small_dump = 1978
 integer, parameter, public :: is_not_mhd = 1979

 procedure(read_dump_fortran), pointer :: read_dump => read_dump_fortran
 procedure(read_smalldump_fortran), pointer :: read_smalldump => read_smalldump_fortran
 procedure(write_smalldump_fortran), pointer :: write_smalldump => write_smalldump_fortran
 procedure(write_fulldump_fortran), pointer :: write_fulldump => write_fulldump_fortran


contains

subroutine init_readwrite_dumps()
 logical :: lhdf5

#ifdef HDF5
 lhdf5 = .true. ! we can potentially use both file formats if needed

 if (lhdf5) then
    read_dump => read_dump_hdf5
    read_smalldump => read_smalldump_hdf5
    write_smalldump => write_smalldump_hdf5
    write_fulldump => write_fulldump_hdf5

    opened_full_dump => opened_full_dump_hdf5
    dt_read_in => dt_read_in_hdf5
 else
    read_dump => read_dump_fortran
    read_smalldump => read_smalldump_fortran
    write_smalldump => write_smalldump_fortran
    write_fulldump => write_fulldump_fortran

    opened_full_dump => opened_full_dump_fortran
    dt_read_in => dt_read_in_fortran
 endif
#else
 lhdf5 = .false.

 read_dump => read_dump_fortran
 read_smalldump => read_smalldump_fortran
 write_smalldump => write_smalldump_fortran
 write_fulldump => write_fulldump_fortran

 opened_full_dump => opened_full_dump_fortran
 dt_read_in => dt_read_in_fortran
#endif

end subroutine init_readwrite_dumps

 !--------------------------------------------------------------------
 !+
 !  subroutine to write output to full dump file
 !  in GADGET format
 !+
 !-------------------------------------------------------------------
subroutine write_gadgetdump(dumpfile,t,xyzh,particlemass,vxyzu,rho,utherm,npart)
 use io,       only:iprint,idump,real4
#ifdef PERIODIC
 use boundary, only:dxbound
#endif
 real,             intent(in) :: t,particlemass,utherm
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: npart
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: rho(:)

 integer(kind=4) :: particleid(size(rho))
 integer :: npartoftype(6),nall(6),ncrap(6)
 real(kind=8) :: massoftype(6)
 real(kind=8)                          :: time,boxsize
 real(kind=8), parameter               :: dumz = 0.d0
 real(kind=4) :: unused(15)
 integer, parameter :: iflagsfr = 0, iflagfeedback = 0, iflagcool = 0
 integer, parameter :: nfiles = 1
 integer            :: ierr,i,j
 !
 !--open dumpfile
 !
 write(iprint,"(/,/,'-------->   TIME = ',f12.4,"// &
         "': full dump written to file ',a,'   <--------',/)")  t,trim(dumpfile)

 write(iprint,*) 'writing to unit ',idump
 open(unit=idump,file=dumpfile,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'error: can''t create new dumpfile ',trim(dumpfile)
    stop
 endif

 npartoftype(:) = 0
 npartoftype(1) = npart
 nall(:)  = npartoftype(:)
 ncrap(:) = 0
 time     = t
#ifdef PERIODIC
 boxsize = dxbound
#else
 boxsize = 0.
#endif

 massoftype(:) = 0.
 massoftype(1) = particlemass
 unused(:) = 0

 do i=1,npart
    particleid(i) = i
 enddo
 write(idump,iostat=ierr) npartoftype(1:6),massoftype(1:6),time,dumz, &
         iflagsfr,iflagfeedback,nall(1:6),iflagcool,nfiles,boxsize, &
         dumz,dumz,dumz,iflagsfr,iflagsfr,ncrap(1:6),iflagsfr,unused(:)

 write(idump,iostat=ierr) ((real4(xyzh(j,i)),j=1,3),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing positions'
    return
 endif
 write(idump,iostat=ierr) ((real4(vxyzu(j,i)),j=1,3),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing velocities'
    return
 endif
 write(idump,iostat=ierr) (particleid(i),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing particle ID'
    return
 endif
 if (size(vxyzu(:,1)) >= 4) then
    write(idump,iostat=ierr) (real4(vxyzu(4,i)),i=1,npart)
 else
    write(idump,iostat=ierr) (real4(utherm),i=1,npart)
 endif
 if (ierr /= 0) then
    print*,' error writing utherm'
    return
 endif
 write(idump,iostat=ierr) (real4(rho(i)),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing rho'
    return
 endif
 write(idump,iostat=ierr) (real4(xyzh(4,i)),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing h'
    return
 endif
 print*,' finished writing file -- OK'

 return
end subroutine write_gadgetdump

end module readwrite_dumps
