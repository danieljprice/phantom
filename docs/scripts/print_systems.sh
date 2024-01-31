#!/bin/bash
#
# script to print out information about all the SETUP= options
# from the phantom Makefile
#
phantomdir='../../'
url="https://github.com/danieljprice/phantom/blob/master/"
echo ""
echo ".. tabularcolumns:: |p{3cm}|p{15cm}|"
echo ""
echo ".. table:: List of possible SYSTEM configurations"
echo "   :widths: auto"
echo ""
printf "   +"
printf -- '-%.0s' {1..18}
printf "+"
printf -- '-%.0s' {1..123}
printf "+\n"
printf "   | %-16s | %-121s | \n" "SYSTEM=" "description"
printf "   +"
printf -- '=%.0s' {1..18}
printf "+"
printf -- '=%.0s' {1..123}
printf "+\n"
print_system()
{
  system=$1;
  descript=`grep -A 1 "ifeq (\\$(SYSTEM), $system)" $phantomdir/build/Makefile_systems | grep '#' | cut -d'#' -f 2 | tail -1 | xargs`
  #lineno=`grep -n "ifeq (\\$(SETUP), $setup)" $phantomdir/build/Makefile_setups | cut -d':' -f 1`
  printf "   | %-16s | %-121s | \n" "$system" "$descript"
  printf "   +"
  printf -- '-%.0s' {1..18}
  printf "+"
  printf -- '-%.0s' {1..123}
  printf "+\n"
}
listofsystems=`grep 'ifeq ($(SYSTEM)' $phantomdir/build/Makefile_systems | grep -v skip | cut -d, -f 2 | cut -d')' -f 1 | sort`
for system in $listofsystems; do
    print_system $system
done
echo
