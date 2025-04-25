!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_subgroup
!
! utils_subgroup
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 integer, parameter :: ck_size = 8
 real,dimension(8),parameter :: cks=(/0.3922568052387800,0.5100434119184585,-0.4710533854097566,&
                                      0.0687531682525181,0.0687531682525181,-0.4710533854097566,&
                                      0.5100434119184585,0.3922568052387800/)
 real,dimension(8),parameter :: cck_sorted=(/0.0976997828427615,0.3922568052387800,0.4312468317474820,&
                                              0.5000000000000000,0.5687531682525181,0.6077431947612200,&
                                              0.9023002171572385,1.0000000000000000/)
 real,dimension(8),parameter :: dks=(/0.7845136104775600,0.2355732133593570,-1.1776799841788701,&
                                      1.3151863206839063,-1.1776799841788701,0.2355732133593570,&
                                      0.7845136104775600,0.0000000000000000/)
 integer,dimension(8),parameter :: cck_sorted_id=(/6,1,3,4,5,7,2,8/)


contains

end module utils_subgroup
