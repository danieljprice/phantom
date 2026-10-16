!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Assemble two isolated stellar full dumps on an incoming hyperbolic orbit.
! Retain internal particle states and spin; add only centre-of-mass motion.
!
! :References: None
!
! :Owner: Francesco Radica
!
! :Runtime parameters:
!   - arg_peri   : *argument of pericentre in degrees*
!   - ecc        : *hyperbolic eccentricity*
!   - pericentre : *pericentre in code units*
!   - secondary  : *secondary stellar full dump*
!   - separation : *initial separation in code units*
!
! :Dependencies: centreofmass, damping, dim, dump_utils, io, moddump_utils,
!   options, part, prompting, readwrite_dumps, readwrite_dumps_common,
!   setbinary, units
!
 use io, only:fatal
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''
 character(len=1024) :: secondary = ''
 real :: separation = 4.,pericentre = 1.25,ecc = 1.6679,arg_peri = 180.

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 use centreofmass, only:get_centreofmass
 use dim, only:maxp,maxvxyzu,mpi,use_apr,gr,mhd,use_dust,ind_timesteps,inject_parts,do_radiation
 use damping, only:idamp
 use options, only:ieos,iexternalforce
 use part, only:igas,particle_mass,copy_particle,iorig,norig,npartoftypetot
 use readwrite_dumps, only:read_dump
 use setbinary, only:set_binary
 use moddump_utils, only:get_moddump_options
 use units, only:udist,umass,utime
 integer, intent(inout) :: npart,npartoftype(:)
 real, intent(inout) :: massoftype(:),xyzh(:,:),vxyzu(:,:)
 integer :: n1,n2,i,ierr,ntemporary,capacity
 integer(kind=8) :: idoffset,idsecondary
 real :: m1,m2,x1(3),v1(3),x2(3),v2(3),units1(3),time2,hfact2
 real :: a,f,cosf,pi,radius1,radius2
 real :: centres(6,2),velocities(3,2)

 if (mpi.or.use_apr.or.gr.or.mhd.or.use_dust.or.ind_timesteps.or.inject_parts.or.do_radiation) &
    call fatal('stellar_encounter','requires Newtonian hydro, shared timesteps, MPI=no and APR=no')
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)
 if (maxvxyzu /= 4 .or. .not.any(ieos==[2,10,12])) &
    call fatal('stellar_encounter','requires an ideal gas, MESA or gas plus radiation EOS')
 call check_star(npart,npartoftype,xyzh,vxyzu)
 if (len_trim(secondary)==0) call fatal('stellar_encounter','missing secondary full dump')
 if (.not.all(ieee_is_finite([separation,pericentre,ecc,arg_peri]))) &
    call fatal('stellar_encounter','orbital parameters must be finite')
 if (.not.(ecc > 1. .and. pericentre > 0. .and. separation > pericentre)) &
    call fatal('stellar_encounter','need ecc>1 and separation>pericentre>0')

 n1 = npart
 units1 = [udist,umass,utime]
 call check_secondary_header(secondary,units1,n2)
 if (int(n1,kind=8)+int(n2,kind=8) > int(maxp,kind=8)) &
    call fatal('stellar_encounter','insufficient capacity: pass --maxp >= combined particle count')
 capacity = maxp
 idoffset = maxval(iorig(1:n1))
 call get_centreofmass(x1,v1,n1,xyzh,vxyzu)
 m1 = sum(particle_mass(1:n1))
 radius1 = 0.
 ! Keep the first star beyond the entire second dump. The capacity preflight
 ! ensures read_dump cannot reallocate or invalidate the caller's arrays.
 do i=1,n1
    radius1 = max(radius1,norm2(xyzh(1:3,i)-x1))
 enddo
 do i=n1,1,-1
    call copy_particle(i,n2+i,new_part=.false.)
 enddo
 call read_dump(trim(secondary),time2,hfact2,93,6,0,1,ierr)
 if (ierr /= 0) call fatal('stellar_encounter','cannot read secondary full dump')
 if (maxp /= capacity .or. npart /= n2) call fatal('stellar_encounter','secondary header or capacity changed unexpectedly')
 call check_star(npart,npartoftype,xyzh,vxyzu)
 if (any(abs([udist,umass,utime]/units1-1.) > 1.e-12)) call fatal('stellar_encounter','incompatible units')
 call get_centreofmass(x2,v2,n2,xyzh,vxyzu)
 m2 = sum(particle_mass(1:n2))
 radius2 = 0.
 do i=1,n2
    radius2 = max(radius2,norm2(xyzh(1:3,i)-x2))
 enddo
 if (separation <= radius1+radius2) call fatal('stellar_encounter','stars overlap at the initial separation')
 idsecondary = maxval(iorig(1:n2))
 if (idsecondary > huge(idoffset)-idoffset) call fatal('stellar_encounter','particle ID offset would overflow')

 pi = acos(-1.)
 a = -pericentre/(ecc-1.)
 cosf = (pericentre*(1.+ecc)/separation-1.)/ecc
 f = -acos(max(-1.,min(1.,cosf)))*180./pi
 ntemporary = 0
 call set_binary(m1,m2,a,ecc,0.,0.,centres,velocities,ntemporary,ierr, &
                 posang_ascnode=-90.,arg_peri=arg_peri,incl=0.,f=f)
 if (ierr /= 0) call fatal('stellar_encounter','could not construct hyperbolic orbit')
 if (dot_product(centres(1:3,2)-centres(1:3,1),velocities(:,2)-velocities(:,1)) >= 0.) &
    call fatal('stellar_encounter','orbit is not incoming')

 do i=1,n2
    iorig(i) = iorig(i)+idoffset
    xyzh(1:3,i) = xyzh(1:3,i)-x2+centres(1:3,2)
    vxyzu(1:3,i) = vxyzu(1:3,i)-v2+velocities(:,2)
 enddo
 do i=n2+1,n2+n1
    xyzh(1:3,i) = xyzh(1:3,i)-x1+centres(1:3,1)
    vxyzu(1:3,i) = vxyzu(1:3,i)-v1+velocities(:,1)
 enddo
 npart = n1+n2
 npartoftype = 0
 npartoftype(igas) = npart
 npartoftypetot = int(npartoftype,kind=8)
 massoftype = 0.
 massoftype(igas) = (m1+m2)/npart
 norig = maxval(iorig(1:npart))
 idamp = 0
 iexternalforce = 0
 print*, 'Encounter stars (particles, mass, maximum particle radius):'
 print*, n1,m1,radius1
 print*, n2,m2,radius2
 print*, 'Maximum-radius overlap at point-mass pericentre:',radius1+radius2-pericentre
 print*, 'Individual particle masses and internal velocities retained; damping disabled.'
 print*, 'Primary IDs retained; secondary IDs offset by ',idoffset
end subroutine modify_dump

subroutine check_star(n,npartoftype,xyzh,vxyzu)
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 use part, only:variable_mass,nptmass,igas,particle_mass,iorig,rho,eos_vars,igasP,itemp,imu
 use options, only:ieos
 integer, intent(in) :: n,npartoftype(:)
 real, intent(in) :: xyzh(:,:),vxyzu(:,:)
 if (n <= 0 .or. nptmass /= 0 .or. npartoftype(igas) /= n) &
    call fatal('stellar_encounter','input must contain one isolated gas star without sinks')
 if (.not.variable_mass) call fatal('stellar_encounter','input must store individual particle masses')
 if (.not.all(ieee_is_finite(xyzh(:,1:n))) .or. .not.all(ieee_is_finite(vxyzu(:,1:n))) .or. &
     .not.all(ieee_is_finite(particle_mass(1:n))) .or. .not.all(ieee_is_finite(eos_vars(:,1:n)))) &
    call fatal('stellar_encounter','nonfinite particle fields')
 if (.not.all(ieee_is_finite(rho(1:n)))) call fatal('stellar_encounter','nonfinite density')
 if (any(particle_mass(1:n) <= 0.).or.any(xyzh(4,1:n) <= 0.).or.any(vxyzu(4,1:n) <= 0.).or. &
     any(rho(1:n) <= 0.).or.any(iorig(1:n) <= 0)) &
    call fatal('stellar_encounter','nonpositive particle fields or IDs')
 ! The gamma-law EOS does not store pressure or temperature in these dumps.
 if (ieos /= 2) then
    if (any(eos_vars(igasP,1:n) <= 0.).or.any(eos_vars(itemp,1:n) <= 0.)) &
       call fatal('stellar_encounter','nonpositive stored pressure or temperature')
 endif
 if (ieos==12 .and. any(eos_vars(imu,1:n) <= 0.)) &
    call fatal('stellar_encounter','ideal gas plus radiation requires positive molecular weights')
end subroutine check_star

subroutine check_secondary_header(path,units1,n)
 use dump_utils, only:dump_h,lenid,open_dumpfile_r,read_header,extract,free_header
 use readwrite_dumps_common, only:get_options_from_fileid
 use part, only:hfact
 use options, only:ieos
 character(len=*), intent(in) :: path
 real, intent(in) :: units1(3)
 integer, intent(out) :: n
 type(dump_h) :: hdr
 character(len=lenid) :: identifier
 logical :: tagged,phantomdump,smalldump,onefluid
 integer :: lu,ierr,ierrfree,ieos2,massmode,nblocks,ierrs(3)
 integer(kind=8) :: n8
 real :: units2(3),hfact2
 lu = 94
 call open_dumpfile_r(lu,path,identifier,ierr)
 if (ierr /= 0) call fatal('stellar_encounter','cannot open secondary full dump')
 call get_options_from_fileid(identifier,tagged,phantomdump,smalldump,onefluid,ierr)
 if (ierr /= 0 .or. .not.tagged .or. .not.phantomdump .or. onefluid) &
    call fatal('stellar_encounter','secondary must be a tagged Phantom hydro full dump')
 call read_header(lu,hdr,ierr,tagged=tagged)
 if (ierr /= 0) call fatal('stellar_encounter','cannot read secondary header')
 close(lu)
 call extract('nparttot',n8,hdr,ierr)
 if (ierr /= 0 .or. n8 <= 0 .or. n8 > huge(n)) call fatal('stellar_encounter','invalid secondary count')
 n = int(n8)
 call extract('ieos',ieos2,hdr,ierr)
 if (ierr /= 0 .or. ieos2 /= ieos) call fatal('stellar_encounter','secondary EOS differs from primary EOS')
 call extract('variable_mass',massmode,hdr,ierr)
 if (ierr /= 0 .or. massmode /= 1) call fatal('stellar_encounter','secondary does not store individual masses')
 call extract('nblocks',nblocks,hdr,ierr)
 if (ierr /= 0 .or. nblocks /= 1) call fatal('stellar_encounter','secondary must be a single-block dump')
 call extract('udist',units2(1),hdr,ierrs(1))
 call extract('umass',units2(2),hdr,ierrs(2))
 call extract('utime',units2(3),hdr,ierrs(3))
 if (any(ierrs /= 0).or.any(abs(units2/units1-1.) > 1.e-12)) call fatal('stellar_encounter','incompatible units')
 call extract('hfact',hfact2,hdr,ierr)
 if (ierr /= 0 .or. abs(hfact2/hfact-1.) > 1.e-12) call fatal('stellar_encounter','incompatible smoothing factor')
 call free_header(hdr,ierrfree)
end subroutine check_secondary_header

subroutine moddump_interactive()
 use prompting, only:prompt
 call prompt('Secondary full dump',secondary)
 call prompt('Initial separation in code units',separation)
 call prompt('Pericentre in code units',pericentre)
 call prompt('Eccentricity',ecc)
 call prompt('Argument of pericentre in degrees',arg_peri)
end subroutine moddump_interactive

subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit
 call write_inopt(secondary,'secondary','secondary stellar full dump',iunit)
 call write_inopt(separation,'separation','initial separation in code units',iunit)
 call write_inopt(pericentre,'pericentre','pericentre in code units',iunit)
 call write_inopt(ecc,'ecc','hyperbolic eccentricity',iunit)
 call write_inopt(arg_peri,'arg_peri','argument of pericentre in degrees',iunit)
end subroutine write_options_moddump

subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer, intent(inout) :: nerr
 call read_inopt(secondary,'secondary',db,errcount=nerr)
 call read_inopt(separation,'separation',db,errcount=nerr,min=0.)
 call read_inopt(pericentre,'pericentre',db,errcount=nerr,min=0.)
 call read_inopt(ecc,'ecc',db,errcount=nerr,min=1.)
 call read_inopt(arg_peri,'arg_peri',db,errcount=nerr)
end subroutine read_options_moddump

end module moddump
