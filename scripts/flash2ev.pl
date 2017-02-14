#!/usr/bin/perl
#
# swaps column order in FLASH evolution files to match Phantom .ev files
#
use strict;
if ($#ARGV != 0) {
   print "flash2ev: script to convert flash evolution file to Phantom .ev file\n";
   die "Usage: $0 flash.dat \n";
}

my ($flashfile) = @ARGV;
open(FLASHFILE,"< $flashfile");
open(EVFILE,"> $flashfile.ev");

while (<FLASHFILE>) {
my(
$time,
$mass,
$x_momentum,
$y_momentum,
$z_momentum,
$E_total,
$E_kinetic,
$E_internal,
$rms_Mach,
$rms_Mach_mw,
$rms_Mach_netto,
$rms_Mach_netto_mw,
$rms_Forcing,
$mean_temperature,
$mean_pressure,
$min_density,
$max_density,
$abs_vorticity,
$abs_vorticity_sqr,
$abs_vorticity_mw,
$abs_div_vel,
$abs_div_vel_sqr,
$abs_div_vel_mw,
$abs_rot_accel_norm,
$abs_rot_accel_norm_sq,
$abs_rot_accel_norm_mw,
$abs_div_accel,
$abs_div_accel_sqr,
$abs_div_accel_mw,
$E_magnetic,
$mean_Bx,
$rms_Bx,
$mean_By,
$rms_By,
$mean_Bz,
$rms_Bz,
$plasma_beta,
$mean_divB,
$rms_divB,
$mean_dens,
$rms_dens,
$sigma_dens,
$mean_ln_dens,
$rms_ln_dens,
$sigma_ln_dens ) = split ' ',$_;
print EVFILE "$time $E_kinetic $E_internal $E_magnetic 0.0 $E_total $x_momentum 0.0 $max_density 0.0 0.0 $rms_Mach \n";   
}
close(FLASHFILE);
close(EVFILE);
print "data written to $flashfile.ev\n";
