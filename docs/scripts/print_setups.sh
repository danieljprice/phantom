#!/bin/bash
#
# script to print out information about all the SETUP= options
# from the phantom Makefile
#
phantomdir='../../'
url="https://github.com/danieljprice/phantom/blob/master/"
echo "phantom compile-time configurations"
echo "==================================="
echo ""
printf "+"
printf -- '-%.0s' {1..105}
printf "+"
printf -- '-%.0s' {1..63}
printf "+"
printf -- '-%.0s' {1..63}
printf "+"
printf -- '-%.0s' {1..63}
printf "+\n"
printf "| %-103s | %-61s | %-61s | %-61s |  \n" "setup"  "description"  "compile-time options" "initial conditions file"
printf "+"
printf -- '=%.0s' {1..105}
printf "+"
printf -- '=%.0s' {1..63}
printf "+"
printf -- '=%.0s' {1..63}
printf "+"
printf -- '=%.0s' {1..63}
printf "+\n"
allsetups=`grep 'ifeq ($(SETUP)' $phantomdir/build/Makefile_setups | grep -v skip | cut -d, -f 2 | cut -d')' -f 1 | sort`
for setup in $allsetups; do
    descript=`grep -A 1 "ifeq (\\$(SETUP), $setup)" $phantomdir/build/Makefile_setups | grep '#' | cut -d'#' -f 2 | tail -1 | xargs`
    lineno=`grep -n "ifeq (\\$(SETUP), $setup)" $phantomdir/build/Makefile_setups | cut -d':' -f 1`
    options=`cd $phantomdir; make SETUP=$setup get_setup_opts`
    setupfile=`cd $phantomdir; make SETUP=$setup get_setup_file`
    lastfile='';
    for x in $setupfile; do
        lastfile=$x;
    done
    printf "| %-103s | %-61s | %-61s | %-61s |  \n" "\`$setup <$url/build/Makefile_setups#L$lineno>\`__" "$descript" "$options" "\`$lastfile <$url/src/setup/$lastfile>\`__"
    printf "+"
    printf -- '-%.0s' {1..105}
    printf "+"
    printf -- '-%.0s' {1..63}
    printf "+"
    printf -- '-%.0s' {1..63}
    printf "+"
    printf -- '-%.0s' {1..63}
    printf "+\n"
done
echo
