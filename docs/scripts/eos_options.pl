#!/usr/bin/perl
#
# @(#) Perl script to extract sink particle properties from src/main/part.F90
#
my $start = 0;
open(FILE,'../../src/main/eos.f90');
while (<FILE>) {
  my $line = $_;
  if ( m/select case\(eos_type\)/) {
     print ".. table:: Equations of state implemented in phantom\n";
     print "   :widths: auto\n\n";
     print "   +-----------+","-" x 82,"+\n";
     printf("   | ieos      | %-80s | \n","Description");
     print "   +===========+","=" x 82,"+\n";
     $start = 1;
  } elsif ($start == 1) {
     # last entry, close on matching ---- line
     if (m/^\s*case\(([\d,\,]+)\).*/) {
        if ($printed_case_num) {
           print "   +-----------+","-" x 82,"+\n";
        }
        $case_num = $1;
        $printed_case_num = 0;
        #exit();
     } elsif (m/^\!\s*$/ and !$printed_case_num) {
        next;
     } elsif (m/^\!--\s*(.*)/ or m/^\!\s*(.*)/) {
        if (!$printed_case_num) {
           printf("   | %-2d        | %-80s |\n",$case_num,substr("**$1**", 0, 80));
           $printed_case_num = 1;
        } else {
           printf("   |           | %-80s |\n",substr($1, 0, 80)); # additional comments
        }
     } elsif (m/end\s+select/) {
        print "   +-----------+","-" x 82,"+\n";
        exit();
     }
  }
}
