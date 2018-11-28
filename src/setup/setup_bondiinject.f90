module setup
 implicit none
 public :: setpart

 private

!
!-- Defaults for setup-time parameters
!
 logical :: filldomain = .true.
 real    :: pmassi      = 4.e-4

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma_eos,hfact,time,fileprefix)
 use part,           only:igas,gr,xyzmh_ptmass,vxyz_ptmass
 use options,        only:iexternalforce
 use units,          only:set_units
 use inject,         only:init_inject,inject_particles,dtsphere,rin,inject_interactive
 use timestep,       only:tmax
 use io,             only:iprint
 use eos,            only:gamma
 use prompting,      only:prompt
 use metric,         only:imetric
 use metric_tools,   only:imet_schwarzschild
 use externalforces, only:accradius1,accradius1_hard
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma_eos,hfact
 real,              intent(inout) :: time
 character(len=*),  intent(in)    :: fileprefix
 logical :: iexist
 integer :: ierr,nspheres
 real :: dtinject,tinfall

 if (.not.gr) call fatal('setup','This setup only works with GR on')
 if (imetric/=imet_schwarzschild) call fatal('setup','This setup is meant for use with the Schwarzschild metric')
 call set_units(G=1.,c=1.)

 time            = 0.
 polyk           = 0.
 iexternalforce  = 1

 call read_write_setupfile(id,fileprefix)

 !-- Don't overwrite these values if infile exists
 inquire(file=trim(fileprefix)//'.in',exist=iexist)
 if (.not.iexist) then
    tmax            = 360.
    accradius1      = 2.1
    accradius1_hard = accradius1
    write(iprint,'(/,a,/)') trim(fileprefix)//'.in not found'
    write(iprint,*) 'Using interactive setup to set injection parameters:'
    call inject_interactive()
    call prompt('Enter tmax',tmax,0.)
    call prompt('Enter accretion radius',accradius1,0.)
 endif

 npart            = 0
 npartoftype(:)   = 0
 massoftype(igas) = pmassi
 gamma            = 5./3. !Set gamma in module eos since init_inject needs to know about it.
 gamma_eos        = gamma

 call init_inject(ierr)

 if (filldomain) then
!
!--- Inject 'real' spheres into the whole domain
!    (get the number of spheres required self-consistently from the infall time)
!
    call get_tinfall(tinfall,r1=accradius1,r2=rin,gamma=gamma)
    nspheres = int(tinfall/dtsphere) !27!100!20!
    print*,'number of "real" spheres: ',nspheres
    call inject_particles(dtsphere*nspheres,dtsphere*nspheres,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,dtinject)
 endif

end subroutine setpart

!
!-- Basic integration of the velocty to get the infall time between two points
!
subroutine get_tinfall(tinfall,r1,r2,gamma)
 use bondiexact, only:get_bondi_solution
 real, intent(in)  :: r1,r2,gamma
 real, intent(out) :: tinfall
 integer, parameter :: N=10000
 integer :: i
 real :: dr,rhor,vr,ur,r
 dr = abs(r2-r1)/N
 tinfall = 0.

 r = r1 + dr
 do i=1,N
    call get_bondi_solution(rhor,vr,ur,r,mass=1.,gamma_eos=gamma)
    tinfall = tinfall + dr/abs(vr)
    r = r + dr
 enddo

end subroutine get_tinfall

!
!---Read/write setup file--------------------------------------------------
!

subroutine read_write_setupfile(id,fileprefix)
 use io, only:master
 integer,          intent(in) :: id
 character(len=*), intent(in) :: fileprefix
 character(len=120) :: filename
 logical :: iexist
 integer :: ierr

 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Bondi injection in GR'
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
   if (id==master) then
      call setup_interactive()
      call write_setupfile(filename)
      print*,' Edit '//trim(filename)//' and rerun phantomsetup'
   endif
   stop
 endif

end subroutine read_write_setupfile

subroutine setup_interactive()
 use prompting, only:prompt

 call prompt('Enter particle mass',pmassi,0.)
 call prompt('Dou want to prefill the domain with gas?',filldomain)

end subroutine setup_interactive

subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'
 call write_inopt(pmassi,    'pmassi',    'particle mass',                           iunit)
 call write_inopt(filldomain,'filldomain','filldomain to accretion radius (logical)',iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(pmassi,    'pmassi',    db,min=0.,errcount=nerr)
 call read_inopt(filldomain,'filldomain',db,       errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
