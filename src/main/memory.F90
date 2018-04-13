module memory
   use io,  only:fatal
 implicit none

 public :: allocate_memory

contains

subroutine allocate_memory
 use dim,   only:maxp

 use part,  only:xyzh
 use part,  only:xyzh_soa

 use dim,   only:maxvxyzu
 use part,  only:vxyzu

 use dim,   only:nalpha, maxalpha
 use part,  only:alphaind

 use dim,   only:ndivcurlv
 use part,  only:divcurlv

 use dim,   only:ndivcurlB
 use part,  only:divcurlB

 use dim,   only:maxBevol,maxmhd
 use part,  only:Bevol

 use part,  only:Bxyz

 integer :: allocstat

 allocate(xyzh(4, maxp), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','xyzh allocation error')

 allocate(xyzh_soa(maxp, 4), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','xyzh_soa allocation error')

 allocate(vxyzu(maxvxyzu,maxp), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','vxyzu allocation error')

 allocate(alphaind(nalpha,maxalpha), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','alphaind allocation error')

 allocate(divcurlv(ndivcurlv,maxp), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','divcurlv allocation error')

 allocate(divcurlB(ndivcurlB,maxp), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','divcurlB allocation error')

 allocate(Bevol(maxBevol,maxmhd), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','Bevol allocation error')

 allocate(Bxyz(3,maxmhd), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','Bxyz allocation error')

end subroutine allocate_memory

end module memory
