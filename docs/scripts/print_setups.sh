#!/bin/bash
#
# script to print out information about all the SETUP= options
# from the phantom Makefile
#
phantomdir='../../'
url="https://github.com/danieljprice/phantom/blob/master/"
echo ""
echo ".. tabularcolumns:: |p{3cm}|p{10cm}|p{5cm}|p{5cm}|"
echo ""
echo ".. table:: List of pre-cooked SETUP configurations"
echo "   :widths: auto"
echo ""
printf "   +"
printf -- '-%.0s' {1..18}
printf "+"
printf -- '-%.0s' {1..63}
printf "+"
printf -- '-%.0s' {1..52}
printf "+"
printf -- '-%.0s' {1..123}
printf "+\n"
printf "   | %-16s | %-61s | %-50s | %-121s |  \n" "SETUP=" "description" "compile-time options"  "initial conditions file"
printf "   +"
printf -- '=%.0s' {1..18}
printf "+"
printf -- '=%.0s' {1..63}
printf "+"
printf -- '=%.0s' {1..52}
printf "+"
printf -- '=%.0s' {1..123}
printf "+\n"
print_setup()
{
  setup=$1;
  descript=$(grep -A 3 "ifeq (\$(SETUP), $setup)" "$phantomdir/build/Makefile_setups" 2>/dev/null | grep '^#   ' | head -1 | sed 's/^#   *//')
  if [ -z "$descript" ]; then
    descript=$(echo "$setup" | sed 's/_/ /g')
  fi
  options=$(cd "$phantomdir" && make SETUP=$setup get_setup_opts 2>/dev/null | sed -e 's/, no\b//g' -e 's/\bno, *//g' -e 's/^no$//' | xargs)
  setupfile=$(cd "$phantomdir" && make SETUP=$setup get_setup_file 2>/dev/null)
  lastfile='';
  for x in $setupfile; do
    lastfile=$x;
  done
  printf "   | %-16s | %-61s | %-50s | %-121s |  \n" "$setup" "$descript" "$options" "\`$lastfile <$url/src/setup/$lastfile>\`__"
  printf "   +"
  printf -- '-%.0s' {1..18}
  printf "+"
  printf -- '-%.0s' {1..63}
  printf "+"
  printf -- '-%.0s' {1..52}
  printf "+"
  printf -- '-%.0s' {1..123}
  printf "+\n"
}
if [ "$1" == "best" ]; then
   listofsetups='binary cluster disc grtde jet turb wind'
else
   listofsetups=$(grep 'ifeq (\$(SETUP)' "$phantomdir/build/Makefile_setups" | grep -v skip | cut -d, -f 2 | cut -d')' -f 1 | sort)
fi
for setup in $listofsetups; do
    print_setup $setup
done
echo
