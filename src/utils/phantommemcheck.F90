!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantommemcheck
!
!  DESCRIPTION: This program estimates the memory footprint of a phantom run
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantommemcheck [no arguments]
!
!  DEPENDENCIES: dim, part
!+
!--------------------------------------------------------------------------
program phantommemcheck
 implicit none
 integer, parameter :: megabyte = 1048576, gigabyte = 1073741824

 print "(/,a,/)",' Phantom memcheck: we calculate what you desire'

 call estimate_memusage()

 print "(/,a,/)",' Phantom memcheck: may your memories remain'

contains

subroutine estimate_memusage
 use dim
 use part, only:mhd,maxBevol
 integer(kind=8) :: memtot
 integer :: nreal

 print "(1x,a,/,/,a)",trim(modid),' Estimated memory usage for phantom (per MPI thread): '

 nreal = kind(1.0)
 !
 !--kind usually refers to number of bytes
 !  except on NAG, where kind=1 is 4 byte, kind=2 is 8 byte:
 !  handle this case below
 !
 if (nreal < 4) then
    nreal = 4*nreal
 endif

 memtot = 0_8
 memtot = memtot + memusage('xyzh',nreal,4,maxp)
 memtot = memtot + memusage('vxyzu',nreal,maxvxyzu,maxp)
 memtot = memtot + memusage('vxyzuhalf',4,maxvxyzu,maxp)
 memtot = memtot + memusage('fxyzu',nreal,maxvxyzu,maxp)
 memtot = memtot + memusage('div/curl v',4,ndivcurlv,maxp)
 memtot = memtot + memusage('alphaind',4,1,maxalpha)
 if (mhd) then
    print "(a)",' required for MHD: '
    memtot = memtot + memusage('Bevol',4,maxBevol,maxp)
    memtot = memtot + memusage('Bxyz',4,3,maxp)
    memtot = memtot + memusage('dBevol',4,maxBevol,maxp)
    memtot = memtot + memusage('div B',4,1,maxp)
    memtot = memtot + memusage('Bevolhalf',4,maxBevol,maxp)
 endif
 if (maxstrain==maxp) print "(a)",' physical viscosity: '
 memtot = memtot + memusage('straintensor',4,6,maxstrain)
#ifdef IND_TIMESTEPS
 print "(a)",' required with IND_TIMESTEP=yes: '
 memtot = memtot + memusage('ibin',1,1,maxp)
 memtot = memtot + memusage('iphase',1,1,maxp)
 memtot = memtot + memusage('gradh',1,1,maxp)
#else
 memtot = memtot + memusage('iphase',1,1,maxp)
#endif
 print "(a)",' link list storage: '
 memtot = memtot + memusage('ll',4,1,maxp)
 memtot = memtot + memusage('neighlist',4,1,maxneigh)
 memtot = memtot + memusage('ifirstincell',4,1,ncellsmax)
 memtot = memtot + memusage('hmaxcell',4,1,ncellsmax)
!$ memtot = memtot + memusage('ipart_omp_lock (approx)',8,1,maxp/10)

 if (memtot > gigabyte) then
    write(*,"(/,a50,' [',f7.1,'Gb ]')") 'TOTAL',memtot/real(gigabyte)
 else
    write(*,"(/,a50,' [',f7.1,'Mb ]')") 'TOTAL',memtot/real(megabyte)
 endif

end subroutine estimate_memusage

integer(kind=8) function memusage(var,nbytes,nrank,nsize)
 character(len=*), intent(in) :: var
 integer,          intent(in) :: nbytes,nrank,nsize

 memusage = nbytes*nrank*nsize
 !
 !--print out memory usage in Mb
 !
 if (memusage > 0) then
    write(*,"(a50,' [',f7.1,'Mb ]')") trim(var),memusage/real(megabyte)
 endif
end function memusage

end program phantommemcheck
