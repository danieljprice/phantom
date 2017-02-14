#!/bin/bash
if [ $# -ne 1 ]; then
   echo "Usage: $0 128 or 256";
else
   res=$1;
   for grid in 256; do
   for x in 'ekin' 'vx' 'vy' 'vz' 'rho0.33ekin' 'rho0.50ekin' 'rho'; do
       case $x in
       'ekin' )
          comp='2.0';;
       'vx' )
          comp='2.0';;
       'vy' )
          comp='2.0';;
       'vz' )
          comp='2.0';;
       'rho0.33ekin' )
          #comp='1.66666666666667';;
          comp='1.33333333333333';;
       'rho0.50ekin' )
          #comp='1.33333333333333';;
          comp='1.0';;
       'rho' )
          comp='0.33333333333333';;
       * )
          exit;;
       esac
       ~/phantom/utils/time_average_pspec $comp turb0[2-9]*_grid$grid$x.pow turb101_grid$grid$x.pow;
       outfile='averaged_pspec'$res'grid'$grid$x'_all.pow';
       echo $outfile;
       mv averaged_pspec.pow $outfile;
       ~/phantom/utils/time_average_pspec $comp turb0[2-9]1_grid$grid$x.pow turb101_grid$grid$x.pow;
       outfile='averaged_pspec'$res'grid'$grid$x'.pow';
       echo $outfile;
       mv averaged_pspec.pow $outfile;
   done
   done
fi
