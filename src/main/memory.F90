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

 use dim,   only:maxp_growth
 use part,  only:dustprop

 use dim,   only:maxstrain
 use part,  only:straintensor

 use dim,   only:nabundances,maxp_h2
 use part,  only:abundance

 use dim,   only:maxtemp
 use part,  only:temperature

 use dim,   only:maxp_dustfrac
 use part,  only:dustfrac,dustevol,deltav

 use dim,   only:nsinkproperties,maxptmass
 use part,  only:xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink

 use dim,   only:maxgrav
 use part,  only:poten

 use dim,   only:maxmhdni
 use part,  only:n_R

 use dim,   only:maxne
 use part,  only:n_electronT

 use part,  only:eta_nimhd

 use dim,   only:maxlum
 use part,  only:luminosity

 use part,  only:fxyzu

 use part,  only:dBevol

 use part, only:divBsymm

 use part, only:fext

 use part, only:ddustfrac

 use part, only:ddustprop

 use part, only:vpred

 use part, only:dustpred

 use part, only:Bpred

 use part, only:dustproppred

#ifdef IND_TIMESTEPS
 use part, only:ibin

 use part, only:ibin_old

 use part, only:ibin_wake

 use part, only:dt_in

 use part, only:twas
#endif

 use part, only:iphase, iphase_soa

 use dim,  only:ngradh
 use part, only:gradh

 use part, only:tstop

 use part, only:ll

 !
 !--for analysis routines, do not allocate any more storage
 !  than is strictly necessary. This will eventually be deprecated by the
 !  memory manager.
 !
#ifdef ANALYSIS
 integer, parameter :: maxan = 0
 integer, parameter :: maxmhdan = 0
 integer, parameter :: maxdustan = 0
#else
 integer, parameter :: maxan = maxp
 integer, parameter :: maxmhdan = maxmhd
 integer, parameter :: maxdustan = maxp_dustfrac
#endif

 integer, parameter :: maxphase = maxan
 integer, parameter :: maxgradh = maxan

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

 allocate(dustprop(5,maxp_growth), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','dustprop allocation error')

 allocate(straintensor(6,maxstrain), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','straintensor allocation error')

 allocate(abundance(nabundances,maxp_h2), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','abundance allocation error')

 allocate(temperature(maxtemp), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','temperature allocation error')

 allocate(dustfrac(maxp_dustfrac), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','dustfrac allocation error')

 allocate(dustevol(maxp_dustfrac), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','dustevol allocation error')

 allocate(deltav(3,maxp_dustfrac), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','deltav allocation error')

 allocate(xyzmh_ptmass(nsinkproperties,maxptmass), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','xyzmh_ptmass allocation error')

 allocate(vxyz_ptmass(3,maxptmass), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','xyzmh_ptmass allocation error')

 allocate(fxyz_ptmass(4,maxptmass), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','fxyz_ptmass allocation error')

 allocate(fxyz_ptmass_sinksink(4,maxptmass), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','fxyz_ptmass_sinksink allocation error')

 allocate(poten(maxgrav), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','poten allocation error')

 allocate(n_R(4,maxmhdni), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','n_R allocation error')

 allocate(n_electronT(maxne), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','maxne allocation error')

 allocate(eta_nimhd(4,maxmhdni), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','eta_nimhd allocation error')

 allocate(luminosity(maxlum), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','luminosity allocation error')

 allocate(fxyzu(maxvxyzu,maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','fxyzu allocation error')

 allocate(dBevol(maxBevol,maxmhdan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','dBevol allocation error')

 allocate(divBsymm(maxmhdan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','divBsymm allocation error')

 allocate(fext(3,maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','fext allocation error')

 allocate(ddustfrac(maxdustan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','ddustfrac allocation error')

 allocate(ddustprop(5,maxp_growth), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','ddustprop allocation error')

 allocate(vpred(maxvxyzu,maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','vpred allocation error')

 allocate(dustpred(maxdustan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','dustpred allocation error')

 allocate(Bpred(maxBevol,maxmhdan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','Bpred allocation error')

 allocate(dustproppred(5,maxp_growth), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','dustproppred allocation error')

#ifdef IND_TIMESTEPS
 allocate(ibin(maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','ibin allocation error')

 allocate(ibin_old(maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','ibin_old allocation error')

 allocate(ibin_wake(maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','ibin_wake allocation error')

 allocate(dt_in(maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','dt_in allocation error')

 allocate(twas(maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','twas allocation error')
#endif

 allocate(iphase(maxphase), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','iphase allocation error')

 allocate(iphase_soa(maxphase), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','iphase_soa allocation error')

 allocate(gradh(ngradh,maxgradh), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','gradh allocation error')

 allocate(tstop(maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','tstop allocation error')

 allocate(ll(maxan), stat = allocstat)
 if (allocstat /= 0) call fatal('memory','ll allocation error')

end subroutine allocate_memory

end module memory
