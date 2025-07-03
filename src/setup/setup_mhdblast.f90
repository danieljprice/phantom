!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the MHD blast wave problem
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, infile_utils, io, kernel, mpidomain,
!   options, part, physcon, setup_params, slab, timestep, unifdis
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for MHD blast wave
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:mhd,isothermal
 use setup_params, only:rhozero,ihavesetupB,npart_total
 use unifdis,      only:set_unifdis
 use io,           only:master,fatal
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use physcon,      only:pi
 use timestep,     only:tmax,dtmax
 use options,      only:nfulldump
 use kernel,       only:wkern,cnormk,radkern2,hfact_default
 use part,         only:Bxyz,igas,periodic
 use mpidomain,    only:i_belong
 use infile_utils, only:infile_exists
 use slab,         only:get_options_slab
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real                             :: deltax,totmass,toten
 real                             :: Bx,By,Bz,Pblast,Pmed,Rblast,r2
 real                             :: plasmaB0,pfrac,plasmaB
 integer                          :: i,ierr,nx
 !
 ! quit if not properly compiled
 !
 if (.not.periodic) call fatal('setup','require PERIODIC=yes')
 if (.not.mhd)      call fatal('setup','require MHD=yes')
 if (isothermal)    call fatal('setup','require ISOTHERMAL=no')
 !
 !--general parameters
 !
 time        = 0.
 hfact       = hfact_default
 Bx          = 10./sqrt(2.0)
 By          = 0.0
 Bz          = 10./sqrt(2.0)
 rhozero     = 1.0
 Pblast      = 100.0
 Pmed        = 1.0
 Rblast      = 0.125
 nx          = 64
 gamma       = 1.4
 polyk       = 0.
 plasmaB0    = 2.0*Pblast/(Bx*Bx + By*By + Bz*Bz)
 plasmaB     = plasmaB0
 ihavesetupB = .true.
 if (.not. infile_exists(fileprefix)) then
    tmax      = 0.020
    dtmax     = 0.005
    nfulldump = 1
 endif

 if (id==master) print "(/,1x,63('-'),1(/,1x,a),/,1x,63('-'),/)", 'MHD Blast Wave.'
 call get_options_slab(fileprefix,id,master,nx,rhozero,ierr,plasmab=plasmaB)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 deltax = dxbound/nx
 !
 ! Put particles on grid
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 ! Finalise particle properties
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 totmass           = rhozero*dxbound*dybound*dzbound
 massoftype(igas)  = totmass/npart_total
 if (id==master) print*,' particle mass = ',massoftype(igas)

 ! Reset magnetic field to get the requested plasma beta
 pfrac = sqrt(plasmaB0/plasmaB)
 Bx = Bx*pfrac
 By = By*pfrac
 Bz = Bz*pfrac

 toten = 0.
 do i=1,npart
    vxyzu(:,i) = 0.
    Bxyz(1,i) = Bx
    Bxyz(2,i) = By
    Bxyz(3,i) = Bz
    r2        = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    if (r2 < Rblast**2) then
       vxyzu(4,i) = Pblast/(rhozero*(gamma - 1.0))
    else
       vxyzu(4,i) = Pmed/(rhozero*(gamma - 1.0))
    endif
 enddo

 if (id==master) then
    write(*,'(2x,a,3es11.4)') 'Magnetic field (Bx,By,Bz): ',Bx,By,Bz
    write(*,'(2x,a,2es11.4)') 'Pressure in blast, medium: ',Pblast,Pmed
    write(*,'(2x,a,2es11.4)') 'Plasma beta in blast, medium: ',plasmaB,2.0*Pmed/(Bx*Bx + By*By + Bz*Bz)
    write(*,'(2x,a, es11.4)') 'Initial blast radius: ',Rblast
 endif

end subroutine setpart

end module setup
