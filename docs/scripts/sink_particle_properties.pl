#!/usr/bin/perl
#
# @(#) Perl script to extract release notes from src/splash.f90
#
my $start = 0;
open(FILE,'../../src/main/part.F90');
while (<FILE>) {
  my $line = $_;
  if ( m/(^\!--sink particles$)/) {
     print "+------------------------------------------------------+\n";
     print "| Index     | Description                              | \n";
     print "+------------------------------------------------------+\n";
     $start = 1;
  } elsif ($start == 1) {
     # last entry, close on matching ---- line
     if (m/integer, parameter :: (.*) = (\d+)\s+\!\s+(.*)/) {
        printf "| %-9s | %-40s |\n",$1,$3;
        #exit();
     } elsif (m/^\!.*/) {
        # skip
     } else {
        print "+------------------------------------------------------+\n";
        exit();
     }
  }
}
