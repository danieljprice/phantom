!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_Bfield
!
!  DESCRIPTION:
!   This routine implements external forces related to various external
!   magnetic field configurations
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    Bphi          -- toroidal magnetic field
!    Rtorus        -- radius of Torus containing tokamak B field
!    a_on_R        -- aspect ratio (a/R) of Torus containing tokamak B field
!    currJ0        -- normalisation of magnetic current at R=0
!    itype_externB -- type of external B field (0=none,1=Wesson torus,2=Gaussian torus)
!    nutorus       -- winding number of Torus (for itype_externB = 1)
!
!  DEPENDENCIES: infile_utils, io, physcon
!+
!--------------------------------------------------------------------------
module extern_Bfield
 implicit none
 private

 public :: get_externalB_force, externBfield
 public :: write_options_externB, read_options_externB, check_externB_settings
 public :: vec_rthetaphi_to_xyz, vec_xyz_to_rthetaphi

 ! default values for runtime options
 real, public :: Rtorus      = 1.
 real, public :: a_on_Rtorus = 0.2
 real, public :: currJ0      = 1.

 integer, private :: itype_externB = 0
 integer, private :: nutorus       = 2
 real,    private :: Bphi          = 0.

 real, private :: atorus, da2

contains

!--------------------------------------------------------
!+
!    This subroutine and associated routines handles
!    everything to do with external B fields
!    (with non-zero spatial derivatives)
!
!    Originally written by D. Price & C. Toniolo, 2006
!    Rewritten for Phantom by D. Price 2015
!+
!--------------------------------------------------------

subroutine externBfield(xi,yi,zi,hi,vxi,vyi,vzi,rhoi, &
                        Bintx, Binty, Bintz, &
                        currJintx, currJinty, currJintz, &
                        Bextx, Bexty, Bextz, &
                        fextx, fexty, fextz, &
                        vdotgradBx, vdotgradBy, vdotgradBz, &
                        string)

 use io,  only:warning,error
 real, intent(in) :: xi,yi,zi,hi,vxi,vyi,vzi,rhoi
 real, intent(in) :: Bintx,Binty,Bintz
 real, intent(in) :: currJintx,currJinty,currJintz
 character(len=*), intent(in) :: string
 real, intent(out) :: Bextx,Bexty,Bextz,fextx,fexty,fextz
 real, intent(out) :: vdotgradBx,vdotgradBy,vdotgradBz

 ! local variables
 real rcyl,drcyl,rintorus2,rintorus,ra2,term
 real costheta,sintheta,cosphi,sinphi,drintorus
 real Bextr,Bexttheta,Bextphi,dBthetadr
 real currJextr,currJexttheta,currJextphi
 real currJextx,currJexty,currJextz
 real frtorus,drhoi,valfven2,v2
 real currJintBextx,currJintBexty,currJintBextz
 real currJextBintx,currJextBinty,currJextBintz
!
!--get coordinate factors in torus coordinate system
!
 call get_torus_factors(xi,yi,zi,costheta,sintheta,cosphi,sinphi, &
                        rcyl,drcyl,rintorus,rintorus2,drintorus)
!
!--get 1/rho
!
 if (rhoi > epsilon(rhoi)) then
    drhoi = 1./rhoi
 else
    drhoi = 0.
    call warning('externB',' rho <= tiny in externalBfield !!!')
    return
 endif


!============ choose external B field profile ===============
!
!     iBext is a parameter which allows one to have different
!     choices for the external field profile
!
 select case(itype_externB)
 case(1)
!
!   these are the profiles as described in Wesson, Tokamaks
!
    ra2 = rintorus2*da2
    term = 1. - ra2
!
!   external currents
!
    currJextr = 0.
    currJexttheta = 0.
    currJextphi = currj0*term**nutorus
!
!   Bfield in torus "theta" direction
!
    Bextr = 0.
    Bexttheta = 0.5*currj0*atorus**2/(nutorus+1)*(1. - term**(nutorus+1))*drintorus
    Bextphi = Bphi*Rtorus*drcyl
!
!   derivative of Btheta with respect to torus 'r'
!
    dBthetadr = 0.5*currj0*atorus**2/(nutorus+1)* &
                (((2*nutorus + 1)*ra2 + 1.)*term**nutorus - 1.)*drintorus**2
!
!   external force is in torus "r" direction  (J X B)
!
    frtorus = -Bexttheta*currJextphi*drhoi
!
!   estimate of alfven speed used to give right order
!   of magnitude to boundary force
!
    valfven2 = (Bexttheta*Bexttheta + Bextphi*Bextphi)*drhoi
    v2 = valfven2

 case(2)
!
!   this is a Gaussian profile for J designed so that
!   J does not go to zero at the edge but rather tails off gently
!
    ra2 = rintorus2*da2
    term = exp(-4.*ra2)
    currJextr = 0.
    currJexttheta = 0.
    currJextphi = currj0*term

    Bextr = 0.
    Bexttheta = 0.125*currj0*atorus**2*(1.-term)*drintorus
    Bextphi = 0.
    dBthetadr = -0.125*currj0*atorus**2*(1.-term)*drintorus**2 + currj0*term

    frtorus = -Bexttheta*currJextphi*drhoi

    !v2 = gamma*2./3.*polyk*(rhozero)**(gamma-1.)
    valfven2 = (Bexttheta*Bexttheta + Bextphi*Bextphi)*drhoi
    v2 = valfven2 !+ v2
 case default
!
!   Default: no external B field; do nothing and return
!
    Bextx = 0.
    Bexty = 0.
    Bextz = 0.
    fextx = 0.
    fexty = 0.
    fextz = 0.
    vdotgradBx = 0.
    vdotgradBy = 0.
    vdotgradBz = 0.
    return
 end select

!============================================================!
 select case(trim(string))
 case('Bfield')
!
!--in this case return only the external B field
!  is returned in cartesian co-ordinates
!
    call vec_rthetaphi_to_xyz(Bextr,Bexttheta,Bextphi,Bextx,Bexty,Bextz, &
                              costheta,sintheta,cosphi,sinphi)

 case('fext')
!
!--add boundary force to force in torus 'r' direction
!
    !frtorus = frtorus + fbound(rintorus,atorus,hzero,v2)
!
!--get J_ext x B_ext as the external force
!
    call vec_rthetaphi_to_xyz(frtorus,0.,0.,fextx,fexty,fextz, &
                              costheta,sintheta,cosphi,sinphi)

 case('all')
!
!--in this case return the external force
!  including mixed Jint x Bext terms and also
!  the term -v.grad Bext needed in the B evolution equation
!
    call vec_rthetaphi_to_xyz(Bextr,Bexttheta,Bextphi,Bextx,Bexty,Bextz, &
                              costheta,sintheta,cosphi,sinphi)
!
!--external currents
!
    call vec_rthetaphi_to_xyz(currJextr,currJexttheta,currJextphi,&
                              currJextx,currJexty,currJextz,&
                              costheta,sintheta,cosphi,sinphi)
!
!--Add  J_int x B_ext
!
    currJintBextx = currJinty*Bextz - currJintz*Bexty
    currJintBexty = currJintz*Bextx - currJintx*Bextz
    currJintBextz = currJintx*Bexty - currJinty*Bextx
!
!--Add  J_ext x B_int
!
    currJextBintx = currJexty*Bintz - currJextz*Binty
    currJextBinty = currJextz*Bintx - currJextx*Bintz
    currJextBintz = currJextx*Binty - currJexty*Bintx
!
!--add boundary force to force in torus 'r' direction
!
    !frtorus = frtorus + fbound(rintorus,atorus,hzero,v2)
!
!--get J_ext x B_ext
!
    call vec_rthetaphi_to_xyz(frtorus,0.,0.,fextx,fexty,fextz, &
                              costheta,sintheta,cosphi,sinphi)
!
!--construct total external force
!
    fextx = fextx + (currJintBextx + currJextBintx)*drhoi
    fexty = fexty + (currJintBexty + currJextBinty)*drhoi
    fextz = fextz + (currJintBextz + currJextBintz)*drhoi
!
!--next calculate the advection terms
!
    call get_advectionterm(Bextr,Bexttheta,Bextphi,dBthetadr, &
             vxi,vyi,vzi,vdotgradBx,vdotgradBy,vdotgradBz, &
             drintorus,drcyl,costheta,sintheta,cosphi,sinphi)

 case default
    call error('externB','unknown string in call to externBfield')
 end select

 return
end subroutine externBfield

!--------------------------------------------------------
!+
!  This subroutine computes the v.grad B term
!  assuming only Btheta has a gradient
!+
!--------------------------------------------------------
subroutine get_advectionterm(br,btheta,bphi,dbthetadr, &
            vx,vy,vz,vdotgradbx,vdotgradby,vdotgradbz, &
            drintorus,drcyl,costheta,sintheta,cosphi,sinphi)
 real br,btheta,bphi,dbthetadr,vx,vy,vz
 real drintorus,drcyl,costheta,sintheta,cosphi,sinphi
 real vdotgradbx,vdotgradby,vdotgradbz

 real dbxdr,dbxdtheta,dbxdphi
 real dbydr,dbydtheta,dbydphi
 real dbzdr,dbzdtheta,dbzdphi
 real drdx,dthetadx,dphidx
 real drdy,dthetady,dphidy
 real drdz,dthetadz !,dphidz
 real dbxdx,dbxdy,dbxdz
 real dbydx,dbydy,dbydz
 real dbzdx,dbzdy,dbzdz

!
!--get gradients of cartesian components of B
!  with respect to toroidal coordinates
!
 dBxdr     = -dBthetadr*sintheta*cosphi
 dBxdtheta = -Br*sintheta*cosphi - Btheta*costheta*cosphi
 dBxdphi   = -Br*costheta*sinphi +Btheta*sintheta*sinphi -Bphi*cosphi

 dBydr     = -dBthetadr*sintheta*sinphi
 dBydtheta = -Br*sintheta*sinphi - Btheta*costheta*sinphi
 dBydphi   = Br*costheta*cosphi - Btheta*sintheta*cosphi -Bphi*sinphi

 dBzdr     = dBthetadr*costheta
 dBzdtheta = Br*costheta - Btheta*sintheta
 dBzdphi   = 0.
!
!--set transformation factors
!
 drdx = costheta*cosphi
 drdy = costheta*sinphi
 drdz = sintheta

 dthetadx = -sintheta*cosphi*drintorus
 dthetady = -sintheta*sinphi*drintorus
 dthetadz = costheta*drintorus

 dphidx = -sinphi*drcyl
 dphidy = cosphi*drcyl
 !dphidz = 0.
!
!--translate to get grad B in cartesians
!
 dBxdx = dBxdr*drdx + dBxdtheta*dthetadx + dBxdphi*dphidx
 dBxdy = dBxdr*drdy + dBxdtheta*dthetady + dBxdphi*dphidy
 dBxdz = dBxdr*drdz + dBxdtheta*dthetadz !+ dBxdphi*dphidz

 dBydx = dBydr*drdx + dBydtheta*dthetadx + dBydphi*dphidx
 dBydy = dBydr*drdy + dBydtheta*dthetady + dBydphi*dphidy
 dBydz = dBydr*drdz + dBydtheta*dthetadz !+ dBydphi*dphidz

 dBzdx = dBzdr*drdx + dBzdtheta*dthetadx + dBzdphi*dphidx
 dBzdy = dBzdr*drdy + dBzdtheta*dthetady + dBzdphi*dphidy
 dBzdz = dBzdr*drdz + dBzdtheta*dthetadz !+ dBzdphi*dphidz
!
!--get v.grad B in cartesians
!
 vdotgradBx = vx*dBxdx + vy*dBxdy + vz*dBxdz
 vdotgradBy = vx*dBydx + vy*dBydy + vz*dBydz
 vdotgradBz = vx*dBzdx + vy*dBzdy + vz*dBzdz

 return
end subroutine get_advectionterm

!------------------------------------------------------------
!+
!  this subroutine transforms a vector in torus coordinates
!  back to cartesian co-ordinates
!+
!------------------------------------------------------------
subroutine vec_rthetaphi_to_xyz(vr,vtheta,vphi,vx,vy,vz,costheta,sintheta,cosphi,sinphi)
 real, intent(in)  :: vr,vtheta,vphi
 real, intent(out) :: vx,vy,vz
 real, intent(in)  :: costheta,sintheta,cosphi,sinphi

 vx = vr*costheta*cosphi - vtheta*sintheta*cosphi - vphi*sinphi
 vy = vr*costheta*sinphi - vtheta*sintheta*sinphi + vphi*cosphi
 vz = vr*sintheta + vtheta*costheta

end subroutine vec_rthetaphi_to_xyz

!------------------------------------------------------------
!+
!  this subroutine transforms a vector in cartesian
!  co-ordinates into torus coordinates
!
!  note: interface is the same as for the previous function
!  so that both functions are called in exactly the same
!  way.
!+
!------------------------------------------------------------
subroutine vec_xyz_to_rthetaphi(vr,vtheta,vphi,vx,vy,vz,costheta,sintheta,cosphi,sinphi)
 real, intent(in)  :: vx,vy,vz
 real, intent(out) :: vr,vtheta,vphi
 real, intent(in)  :: costheta,sintheta,cosphi,sinphi

 vr     =  vx*costheta*cosphi + vy*costheta*sinphi + vz*sintheta
 vtheta = -vx*sintheta*cosphi - vy*sintheta*sinphi + vz*costheta
 vphi   = -vx*sinphi + vy*cosphi

end subroutine vec_xyz_to_rthetaphi

!------------------------------------------------------------
!+
!  this subroutine deals with the coordinate part of the
!  transformation (returns costheta,sintheta,cosphi,sinphi)
!+
!------------------------------------------------------------
subroutine get_torus_factors(xi,yi,zi,costheta,sintheta,cosphi,sinphi,&
           rcyl,drcyl,rintorus,rintorus2,drintorus)
 real, intent(in) :: xi,yi,zi
 real, intent(out) :: rcyl,drcyl,rintorus,rintorus2,drintorus
 real, intent(out) :: costheta,sintheta,cosphi,sinphi

 rcyl = sqrt(xi**2 + yi**2)
 if (rcyl > epsilon(rcyl)) then
    drcyl = 1./rcyl
 else
    drcyl = 0.
 endif
! rintorus is radius from centre of torus
 rintorus2 = (rcyl - rtorus)**2 + zi**2
 rintorus = sqrt(rintorus2)
 if (rintorus > epsilon(rintorus)) then
    drintorus = 1./rintorus
 else
    drintorus = 0.
 endif
 sintheta = zi*drintorus
 costheta = (rcyl-rtorus)*drintorus
 cosphi = xi*drcyl
 sinphi = yi*drcyl

 return
end subroutine get_torus_factors

!------------------------------------------------------------
!+
!  This function acts as an interface to externBfield
!  cutting out the dummy arguments needed on the first call
!+
!------------------------------------------------------------
real function Bexternal(xcoord,ycoord,zcoord,icomponent)
 use io, only:fatal
 integer, intent(in) :: icomponent
 real,    intent(in) :: xcoord,ycoord,zcoord
 real :: dumx,dumy,dumz,Bextx,Bexty,Bextz
 real :: dumh,dumrhoi,dumgx,dumgy,dumgz

 dumrhoi = 1.
 call externBfield(xcoord,ycoord,zcoord,dumh,0.,0.,0.,dumrhoi, &
                   0.,0.,0.,0.,0.,0., &
                   Bextx,Bexty,Bextz, &
                   dumx,dumy,dumz,dumgx,dumgy,dumgz,'Bfield')

 select case(icomponent)
 case(1)
    Bexternal= Bextx
 case(2)
    Bexternal= Bexty
 case(3)
    Bexternal= Bextz
 case default
    Bexternal = 0.
    call fatal('Bexternal','error in Bexternal call')
 end select

end function Bexternal

!------------------------------------------------------------
!+
!  This subroutine acts as an interface to externBfield
!  which returns only the external force component
!  (i.e. can be called from externf for doing
!   hydro relaxation runs in the external tokamak potential)
!+
!------------------------------------------------------------
subroutine get_externalB_force(xi,yi,zi,hi,rhoi,fextx,fexty,fextz)
 real, intent(in)  :: xi,yi,zi,hi,rhoi
 real, intent(out) :: fextx,fexty,fextz
 real :: dumx,dumy,dumz,dumgx,dumgy,dumgz

 call externBfield(xi,yi,zi,hi,0.,0.,0.,rhoi, &
                   0.,0.,0.,0.,0.,0.,dumx,dumy,dumz, &
                   fextx,fexty,fextz,dumgx,dumgy,dumgz,'fext')

end subroutine get_externalB_force


!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_externB(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to force from external magnetic field'
 call write_inopt(itype_externB,'itype_externB','type of external B field (0=none,1=Wesson torus,2=Gaussian torus)',iunit)
 call write_inopt(Rtorus,'Rtorus','radius of Torus containing tokamak B field',iunit)
 call write_inopt(a_on_Rtorus,'a_on_R','aspect ratio (a/R) of Torus containing tokamak B field',iunit)
 call write_inopt(nutorus,'nutorus','winding number of Torus (for itype_externB = 1)',iunit)
 call write_inopt(currJ0,'currJ0','normalisation of magnetic current at R=0',iunit)
 call write_inopt(Bphi,'Bphi','toroidal magnetic field',iunit)

end subroutine write_options_externB

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_externB(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use physcon, only:pi
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 logical, save :: read_Rtorus = .false.
 character(len=30), parameter :: label = 'read_options_externB'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('itype_externB')
    read(valstring,*,iostat=ierr) itype_externB
    if (itype_externB < 0 .or. itype_externB > 2) then
       call fatal(label,'invalid external field setting (itype_externB)')
    endif
    ngot = ngot + 1
 case('Rtorus')
    read(valstring,*,iostat=ierr) Rtorus
    if (Rtorus <= 0.) then
       call fatal(label,'invalid value for Rtorus (must be >0)')
    endif
    read_Rtorus = .true.
    ngot = ngot + 1
 case('a_on_R')
    read(valstring,*,iostat=ierr) a_on_Rtorus
    if (a_on_Rtorus <= 0. .or. a_on_Rtorus >= 1.) then
       call fatal(label,'invalid value for a_on_R (must be 0-1)')
    endif
    if (.not.read_Rtorus) call fatal(label,'must read Rtorus before a_on_R in input file')
    atorus = a_on_Rtorus*Rtorus
    da2 = 1./atorus**2
    ngot = ngot + 1
 case('nutorus')
    read(valstring,*,iostat=ierr) nutorus
    if (nutorus < 1) call fatal(label,'invalid setting for nutorus (must be >0)')
 case('currJ0')
    read(valstring,*,iostat=ierr) currJ0
 case('Bphi')
    read(valstring,*,iostat=ierr) Bphi
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 3)

end subroutine read_options_externB

!-----------------------------------------------------------------------
!+
!  checks that input file options are reasonable
!+
!-----------------------------------------------------------------------
subroutine check_externB_settings(ierr)
 use io, only:error
 integer, intent(out) :: ierr

 ierr = 0
 atorus = a_on_Rtorus*Rtorus
 da2 = 1./atorus**2
 if (atorus <= 0.) then
    call error('externB','atorus not set properly for external Bfield force')
    ierr = 1
 endif

end subroutine check_externB_settings

end module extern_Bfield
