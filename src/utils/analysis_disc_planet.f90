!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for a disc and planet interaction (set with two ptmass)
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, discanalysisutils, infile_utils, io,
!   options, part, physcon, vectorutils
!
 use discanalysisutils, only:disc_analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'disc_planet'
 public :: do_analysis

 integer, parameter :: nr = 100

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use part,    only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use options, only:iexternalforce
 use centreofmass, only:get_total_angular_momentum,reset_centreofmass
 use vectorutils, only:cross_product3D,rotatevec
 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=9) :: output
 character(len=20) :: filename
 character(len=20) :: discprefix
 integer :: i,ierr,iline,ii
 real :: R_in,R_out,H_R,p_index,q_index,M_star
 real :: G,rmin,rmax,tilt(nr),twist(nr)
 real :: rad(nr),h_smooth(nr),sigma(nr),H(nr),dr
 real :: unitlx(nr),unitly(nr),unitlz(nr),ecc(nr)
 real :: psi(nr),tilt_acc(nr),Lx(nr),Ly(nr),Lz(nr),twistprev(nr)
 real :: L_tot(3),L_p(3),L_inner_mag,L_outer_mag
 real :: L_p_mag,L_ratio_inner,L_ratio_outer,e_sinks,ecc_sinks(3)
 real :: rad_sinks,twist_inner,twist_outer,tilt_inner,tilt_outer,minclin,minclout,mannulus
 real :: m_red,mu,rotate_about_y,rotate_about_z,sinks_mass,pos_sinks(3),vel_sinks(3)
 real :: temp(3),temp_mag,term(3),tilt_sinks,twist_sinks,L_tot_mag
 real :: unitl_in(3),unitl_out(3),eff_tilt,eff_outertilt
 real :: break_radius,break_disc_at,eff_innertilt,dlhs,drhs
 integer :: ninbin(nr),n_count_inner,n_count_outer,nptmassinit
 logical :: assume_Ltot_is_same_as_zaxis

 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 integer, parameter :: isinks = 23
 logical :: do_precession,ifile,iexist

 do_precession = .false.

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'angm',numfile
 write(*,'("Output file name is ",A)') output

! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G==1")')
 G = 1.0

 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 inquire(file=trim(discprefix)//'.discparams', exist=ifile)
 if (ifile) then
    call read_discparams(trim(discprefix)//'.discparams',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read .discparams file')
 else
    call read_discparams('discparams.list',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')
 endif

! Print out the parameters
 write(*,*)
 write(*,'("Parameters are:")')
 write(*,*) 'R_in    = ',R_in
 write(*,*) 'R_out   = ',R_out
 write(*,*) 'H/R_ref = ',H_R
 write(*,*) 'p_index = ',p_index
 write(*,*) 'q_index = ',q_index
 write(*,*) 'M_star  = ',M_star
 write(*,*)
 write(*,*)

! Setup rmin and rmax for the analysis
 rmin = R_in
 rmax = R_out

! This variable should be set to false for any discs that use sink particles to set
! the potential or any discs that have a warp
! For any setup that uses iexternalforce and assumes that the vast majority of the angular
! momentum is held by the central potential, this should be set to true

 assume_Ltot_is_same_as_zaxis = .false.

 call disc_analysis(xyzh,vxyz,npart,pmass,time,nr,rmin,rmax,G,M_star,&
                     tilt,tilt_acc,twist,twistprev,psi,H,rad,h_smooth,sigma,unitlx,unitly,unitlz,&
                     Lx,Ly,Lz,ecc,ninbin,assume_Ltot_is_same_as_zaxis,xyzmh_ptmass,vxyz_ptmass,nptmass)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'radius', &
       2,'sigma', &
       3,'<h>/H', &
       4,'lx', &
       5,'ly', &
       6,'lz', &
       7,'tilt', &
       8,'twist', &
       9,'psi', &
       10,'H/R', &
       11,'|e|'

 do i = 1,nr
    if (ninbin(i) > 0) then
       write(iunit,'(13(es18.10,1X))') rad(i),sigma(i),h_smooth(i),unitlx(i),unitly(i),unitlz(i),&
                                         tilt(i),twist(i),psi(i),H(i)/rad(i),ecc(i)
    endif

! Printing time and twist for each radius bin
    if (do_precession) then
       write(filename,"(a,i3.3)")"precess",i
       inquire(file=filename,exist=iexist)
       if (.not.iexist .or. numfile==0) then
          open(unit=iprec,file=filename,status="replace")
          write(iprec,'("# tilt and twist with time for r = ",es18.10)') rad(i)
          write(iprec,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'rad', &
               2,'time', &
               3,'tilt', &
               4,'twist', &
               5,'tot twist', &
               6,'|e|'
       else
          open(unit=iprec,file=filename,status="old",position="append")
       endif
       write(iprec,'(6(es18.10,1X))') rad(i),time,tilt(i),twist(i),twistprev(i),ecc(i)
       close(unit=iprec)
    endif

 enddo

 close(iunit)

 if (iexternalforce /=0) then
    nptmassinit = 0
 else
    nptmassinit = 1
 endif

 ! Prepare for calculations below
 if (nptmass > 0) then
    call reset_centreofmass(npart,xyzh,vxyz,nptmass,xyzmh_ptmass,vxyz_ptmass)
 endif

 ! Prepare for rotations later on
 call get_total_angular_momentum(xyzh,vxyz,npart,L_tot,xyzmh_ptmass,vxyz_ptmass,nptmass)

 ! Rotate first about z-axis
 temp = (/L_tot(1),L_tot(2),0./)
 temp_mag = sqrt(dot_product(temp,temp))
 rotate_about_z = acos(dot_product((/1.,0.,0./),temp/temp_mag))

 ! Rotate second about y-axis
 L_tot_mag = sqrt(dot_product(L_tot,L_tot))
 rotate_about_y = -acos(dot_product((/0.,0.,1./),L_tot/L_tot_mag))

! Set up the radius array
 dr = (rmax-rmin)/real(nr-1)

! Calculating and printing information for the sinks
 if (nptmass>nptmassinit) then
    do i=nptmassinit+1,nptmass
       write(filename,"(a,i3.3)")"sinks_",i-1
       inquire(file=filename,exist=iexist)
       if (.not.iexist .or. numfile == 0) then
          open(isinks,file=filename,status="replace")
          write(isinks,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'time', &
               2,'rad(sph)', &
               3,'Tilt sink', &
               4,'Twist sink', &
               5,'Tilt ID', &
               6,'Tilt OD', &
               7,'Twist ID', &
               8,'Twist OD', &
               9,'Lin/Lp', &
               10,'Lout/Lp',&
               11,'Rel. misalignment', &
               12,'Sink-ID', &
               13,'Sink-OD'

       else
          open(isinks,file=filename,status="old",position="append")
       endif


       ! Properties of sink
       sinks_mass = xyzmh_ptmass(4,i)
       pos_sinks = xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,1)
       rad_sinks = sqrt(dot_product(pos_sinks,pos_sinks))
       vel_sinks = vxyz_ptmass(1:3,i) - vxyz_ptmass(1:3,1)

       ! Identify if there is a break in the disc based on sigma
       break_radius = -1
       do ii=2,nr-2
          dlhs = (sigma(ii) - sigma(ii-1))/(rad(ii) - rad(ii-1))
          drhs = (sigma(ii+1) - sigma(ii))/(rad(ii+1) - rad(ii))
          if (dlhs < 0. .and. drhs > 0.) then
             dlhs = (sigma(ii-1) - sigma(ii-2))/(rad(ii-1) - rad(ii-2))
             drhs = (sigma(ii+2) - sigma(ii+1))/(rad(ii+2) - rad(ii+1))
             if (dlhs < 0. .and. drhs > 0.) break_radius = rad(ii)
          endif
       enddo

       ! If there is a break in the disc use this to decide the inner and outer disc
       break_disc_at = rad_sinks
       if (break_radius > rad(1)) break_disc_at = break_radius

       tilt_inner = 0.
       tilt_outer = 0.
       twist_inner = 0.
       twist_outer = 0.
       L_inner_mag = 0.
       L_outer_mag = 0.
       n_count_inner = 0
       n_count_outer = 0
       unitl_in = 0.
       unitl_out = 0.
       minclin = 0.
       minclout = 0.

       do ii=1,nr
          if (ninbin(ii) > 0) then
             mannulus = 2.*pi*rad(ii)*sigma(ii)*dr
             if (rad(ii) < break_disc_at) then
                n_count_inner = n_count_inner + ninbin(ii)
                minclin = minclin + mannulus
                tilt_inner = tilt_inner + tilt(ii)*mannulus
                twist_inner = twist_inner + twist(ii)*mannulus
                L_inner_mag = L_inner_mag + sqrt(Lx(ii)**2 + Ly(ii)**2 + Lz(ii)**2)
                unitl_in = unitl_in + (/unitlx(ii),unitly(ii),unitlz(ii)/)*mannulus
             else
                n_count_outer = n_count_outer + ninbin(ii)
                minclout = minclout + mannulus
                tilt_outer = tilt_outer + tilt(ii)*mannulus
                twist_outer = twist_outer + twist(ii)*mannulus
                L_outer_mag = L_outer_mag + sqrt(Lx(ii)**2 + Ly(ii)**2 + Lz(ii)**2)
                unitl_out = unitl_out + (/unitlx(ii),unitly(ii),unitlz(ii)/)*mannulus
             endif
          endif
       enddo

       ! Average the tilt and twist of inner and outer disc
       if (n_count_inner > 0) then
          tilt_inner = tilt_inner/minclin
          twist_inner = twist_inner/minclin
          unitl_in = unitl_in/minclin
       else
          tilt_inner = 0.
          twist_inner = 0.
          unitl_in = 0.
       endif
       if (n_count_outer > 0) then
          tilt_outer = tilt_outer/minclout
          twist_outer = twist_outer/minclout
          unitl_out = unitl_out/minclout
       else
          tilt_outer = 0.
          twist_outer = 0.
          unitl_out = 0.
       endif

       call cross_product3D(pos_sinks,vel_sinks,L_p)
       L_p = L_p*sinks_mass
       L_p_mag = sqrt(dot_product(L_p,L_p))

       if (.not.assume_Ltot_is_same_as_zaxis) then
          ! Rotate sinks vector such that Ltot is parallel to z-axis
          call rotatevec(L_p,(/0.,0.,1.0/),-rotate_about_z)
          call rotatevec(L_p,(/0.,1.0,0./),rotate_about_y)
       endif
       tilt_sinks = acos(L_p(3)/L_p_mag)

       ! Angular momentum ratios
       L_ratio_inner = L_inner_mag/L_p_mag
       L_ratio_outer = L_outer_mag/L_p_mag

       ! Calculate the effective tilt between the inner and outer disc
       eff_tilt = acos(dot_product(unitl_in,unitl_out))
       ! Calculate the effective tilt between discs and binary
       ! NB: here we assume that each sink in the binary has the same
       ! direction of AM
       eff_innertilt = acos(dot_product(unitl_in,L_p/L_p_mag))
       eff_outertilt = acos(dot_product(unitl_out,L_p/L_p_mag))

       ! For now, measure twist of sinks as in disc
       twist_sinks = atan2(L_p(2),L_p(1))
       if (twist_sinks < 0.) twist_sinks = twist_sinks + 2.*pi

       ! Calculate the eccentricity of the point mass
       ! These are calculated using the relative position and velocities,
       ! rotations above do not affect these
       m_red = sinks_mass*xyzmh_ptmass(4,1)/(sinks_mass + xyzmh_ptmass(4,1))
       mu = 1.0*(sinks_mass + xyzmh_ptmass(4,1))  !mu = GM, G=1
       call cross_product3D(vel_sinks,L_p/sinks_mass,term)
       ecc_sinks = term/mu - pos_sinks/rad_sinks

       if (ecc_sinks(1) < -1.0) then
          print*,'Warning: eccentricity is undefined, returning e=1.'
          term = 0.
       endif
       e_sinks = sqrt(dot_product(ecc_sinks,ecc_sinks))

       write(isinks,'(13(es18.10,1X))') time,break_disc_at, &
            tilt_sinks,twist_sinks,tilt_inner,tilt_outer, &
            twist_inner,twist_outer,L_ratio_inner,L_ratio_outer, &
            eff_tilt,eff_innertilt,eff_outertilt

       if (numfile == 0) then
          print*,'For sink',i,'located at',rad_sinks
          print*,'and for a reference radius of',break_disc_at
          print*,'AM RATIOS: Lin/Lp=',L_ratio_inner,' and Lout/Lp=',L_ratio_outer
          print*,'Mass inside',real(n_count_inner)*pmass
          print*,'Mass outside',real(n_count_outer)*pmass
       endif

       close(isinks)
    enddo
 endif

end subroutine do_analysis

!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------
subroutine read_discparams(filename,R_in,R_out,H_R,p_index,q_index,M_star,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,H_R,p_index,q_index,M_star
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)

! Read in parameters from the file discparams.list
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(R_in,'R_in',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_out,'R_out',db,ierr)
 if (ierr /= 0) return
 call read_inopt(H_R,'H/R_ref',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_discparams

end module analysis
