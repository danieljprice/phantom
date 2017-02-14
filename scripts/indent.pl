#!/usr/bin/perl
#
# @(#) Perl script to automatically indent Fortran 90 code
#
use strict;
use Text::Wrap;
my $indent = 0;
my $level = 0;
my $mod_indent = 1;
my $loop_indent = 3;
my $continuation;
my $incase = 0;
while (<STDIN>) {
  my $line = "$_";
  chomp ($line);
  #$line =~ m/^\s+(\S*.*)/;
  my $str = $line;
  $str =~ s/^\s*(.*$)/$1/;
  my $increment_next = 0;
  #
  # module, subroutine or program indent
  #
  if (not($str =~ m/^\s*\!/)) {
     if (($str =~ m/^\s*(module|subroutine|program|pure subroutine)\s\S+/) 
         and not($str =~ m/module procedure/)) {
        $increment_next = $mod_indent;
     } elsif ($str =~ m/^\s*end (module|subroutine|program|pure subroutine)\s*.*/) {
        $indent = $indent - $mod_indent;
     }
     #
     # function
     if ($str =~ m/^.*\s(function)\s.*\(.*\)/) {
        $increment_next = $mod_indent;
     } elsif ($str =~ m/^\s*end function/) {
        $indent = $indent - $mod_indent;
     }
     #
     # contains statement
     #
     if ($str =~ m/^\s*contains\s*$/) {
        $indent = 0;
     }
     #
     # if () then .. endif
     #
     if ($str =~ m/(^\s*if|\:\s*if)\s*.*then/) {
        #print "IF\n";
        $increment_next = $loop_indent;
     } elsif ($str =~ m/^\s*(else|elseif)/) {
        #print "ELSE\n";
        $indent = $indent - $loop_indent;
        $increment_next = $loop_indent;
     } elsif ($str =~ m/^\s*end\s*if.*/) {
        #print "ENDIF\n";
        $indent = $indent - $loop_indent;
     }
     #
     # do .. enddo
     #
     if ($str =~ m/(^\s*do|\:\s*do)(\s+|$)/) {
        #print "DO\n";
        $increment_next = $loop_indent;
     } elsif ($str =~ m/^\s*end\s*do/) {
        #print "ENDDO\n";
        $indent = $indent - $loop_indent;
     }
     #
     # type ... end type
     #
     if ($str =~ m/^\s*(type)\s*\w+$/) {
        #print "TYPE\n";
        $increment_next = $mod_indent;
     } elsif ($str =~ m/^\s*end type/) {
        #print "END TYPE\n";
        $indent = $indent - $mod_indent;
     }

     #
     # interface...end interface
     #
     if ($str =~ m/^\s*interface/) {
        $increment_next = $loop_indent;
        #print "INTERFACE $indent $increment_next\n";
     } elsif ($str =~ m/^\s*end\s+interface\s*/) {
        #print "END INTERFACE INDENT=$indent - $loop_indent\n";
        $indent = $indent - $loop_indent;
     }

     #
     # select case
     #
     if ($str =~ m/(^\s*case)/) {
        #print "CASE\n";
        if ($incase == 1) {
           $indent = $indent - $loop_indent;
        } else {
           $incase = 1;
        }
        $increment_next = $loop_indent;  
     } elsif ($str =~ m/(^\s*end select)/) {
        $indent = $indent - $loop_indent;
        $incase = 0;
        #print "END SELECT\n";
     }
  }

  #
  # comment
  #
  my $comment = 0;
  if ($line =~ m/^!/) { $comment = 1 };
  
  # replace my old-style comments with new-style comments
  $_ = $str;
  $str =~ s/^!--((\w|\s).*)/! $1/g;
  
  # preprocessing statements or numbered lines
  my $cpp = 0;
  if ($line =~ m/^(\#|\d+)/) { $cpp = 1 }; 

  #
  # now print the line
  #
  if ($continuation == 1 or $comment == 1 or $cpp == 1) {
     print "$line\n";
  } else {
     if ($indent < 0) { $indent = 0; }
     #print "INDENT $indent \n";
     my $newline = sprintf("%-*s%s\n",$indent,'',$str);
     print $newline;
  }
  $indent = $indent + $increment_next;

  if (($str =~ m/.*&\s*$/) and not ($str =~ m/^\s*\!/)) {
     $continuation = 1;
  } else {
     $continuation = 0;
  }
}
