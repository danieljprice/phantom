 module setup 
 
 use externalforces,   only:iext_brokenstarp,iext_broken_omegapot
 use extern_broken,    only:rbreak,delta_r,bmass1,bmass2,omega1,nratio
 use io,               only:master,warning,error,fatal
 use kernel,           only:hfact_default
 use options,          only:use_dustfrac,iexternalforce,use_hybrid,use_porosity
 use options,          only:use_mcfost,use_mcfost_stellar_parameters
 use part,             only:xyzmh_ptmass,maxvxyzu,vxyz_ptmass,ihacc,ihsoft,&
                            iJ2,ispinx,ispinz,iReff,igas
 use physcon,          only:au,solarm,jupiterm,earthm,pi,twopi,years,hours,deg_to_rad
 use units,            only:umass,udist,utime

 implicit none

 public  :: setpart

 private

 !-- broken star potential
 logical :: brokendisc
 logical :: constomega

contains 

 !--------------------------------------------------------------------------
!+
! This is the only public subroutine of the module
!+
!--------------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use externalforces, only:iexternalforce
 use units, only:set_units,select_unit
 use setdisc, only: set_disc
 use options, only:ieos 
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk
 real,              intent(out)   :: gamma
 real,              intent(out)   :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 write(*,"(/,65('-'),/,/,5x,a,/,/,65('-'))") 'Lets make some hexagons!'

 brokendisc = .false.
 constomega = .false. 

 !--simulation time
 deltat  = 0.1
 norbits = 100

call select_unit(mass_unit,umass,ierr)
if (ierr /= 0) call error('setup_hexagons','mass unit not recognised')
call select_unit(dist_unit,udist,ierr)
if (ierr /= 0) call error('setup_hexagons','length unit not recognised')
call set_units(dist=udist,mass=umass,G=1.d0)





! Setup velocities 

if (iexternalforce==iext_broken_omegapot) then
  rl = rbreak - 0.5*delta_r
  rr = rbreak + 0.5*delta_r
  if (r <= rl) then
    vphi = r*omega1
  elseif (r >= rr) then
    vphi = r*nratio*omega1
  else
    x = (r-rl)/delta_r        ! 0–1
    vphi = omega1 + (nratio*omega1-omega1)*x     ! linear S(x)=x
  endif
    vr = 0.0
elseif (iexternalforce==iext_brokenstarp) then
    rl = rbreak - 0.5*delta_r
    vphi = sqrt(G*mass_eff(R, bmass1, bmass2, rl, delta_r)*R*R/(R*R + eps2_soft)**1.5)
    vr = 0.0
endif  
 !--remind user to check for warnings and errors
 write(*,20)
20 format(/, &
   "-----------------------------------------------------------------",/, &
   "",/, &
   "     Check output for warnings and errors",/, &
   "",/, &
   "-----------------------------------------------------------------",/)


 return

end subroutine setpart

!--------------------------------------------------------------------------

subroutine setup_potential(fileprefix)
 use externalforces, only:iexternalforce
 use externalforces,       only:mass1,accradius1,bmass1,bmass2,delta_r,rbreak
 use options, only:ieos 

 character(len=20), intent(in)    :: fileprefix
select case (iexternalforce)
    case (1)
       print "(/,a)",' Keplerian potential with discontinuity at rbreak'
       print "(a,g10.3,a)",'   Inner mass:      ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   Outer mass:      ', m2,    trim(mass_unit)
       print "(a,g10.3,a)",'   Break radius:    ', rbreak, trim(dist_unit)
       print "(a,g10.3,a)",'   Delta r:         ', delta_r, trim(dist_unit)
       brokendisc = .true.
       bmass1     = m1
       bmass2     = m2
     case (2)
       print "(/,a)",' Solid body rotation with discontinuity at rbreak'
       print "(a,g10.3,a)",'   Omega1:      ', omega1,    1./trim(time_unit)
       print "(a,g10.3,a)",'   nratio:      ', nratio
       print "(a,g10.3,a)",'   rbreak:      ', rbreak,    trim(dist_unit)
       print "(a,g10.3,a)",'   delta_r:     ', delta_r,   trim(dist_unit)
       brokendisc = .true.
end select

end subroutine setup_potential

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')

 call write_inopt(np,'np','number of gas particles',iunit)

select case (iexternalforce)
    case (1)
       !--keplerian potential with discontinuity at rbreak
       call write_inopt(m1,'m1','inner mass',iunit)
       call write_inopt(m2,'m2','outer mass',iunit)
       call write_inopt(rbreak,'rbreak','break radius',iunit)
       call write_inopt(delta_r,'delta_r','delta r',iunit)
    case (2)
       !--keplerian potential with discontinuity at rbreak
       call write_inopt(omega1,'omega1','inner mass',iunit)
       call write_inopt(nratio,'nratio','outer mass',iunit)
       call write_inopt(rbreak,'rbreak','break radius',iunit)
       call write_inopt(delta_r,'delta_r','delta r',iunit)
end select

 call write_inopt(norbits,'norbits','maximum number of orbits at outer disc',iunit)
 call write_inopt(deltat,'deltat','output interval as fraction of orbital period',iunit)

 close(iunit)

end subroutine write_setupfile
!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)

 !--Read values
 print "(a)",' reading setup options from '//trim(filename)


 call close_db(db)
 ierr = nerr
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
 endif
end subroutine read_setupfile

end module setup
