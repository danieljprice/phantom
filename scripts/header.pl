#!/usr/bin/perl
#
# Perl script to generate headers for Phantom source files
#
use Getopt::Long;
use File::Basename;
my ( %opt );
my $replace=undef;
my $show=undef;
my $headerfile;
my $programfile;
#
# OPTIONS HANDLING
#
my $status = GetOptions("replace" => \$opt{"replace"},
                        "show" => \$opt{"show"},
                        "headerfile=s" => \$opt{"headerfile"},
                        "programfile=s" => \$opt{"programfile"});
if (defined $opt{headerfile}) {
  $headerfile=$opt{headerfile};
} else {
  $headerfile="HEADER-module";
  #
  # if headerfile not in current directory
  # then look for it in the same directory
  # as the script was run from
  #
  if (!(-e $headerfile)) {
     $headerfile=dirname($0)."/HEADER-module";
  }
}
if (defined $opt{programfile}) {
  $programfile=$opt{programfile};
} else {
  $programfile="HEADER-program";
  #
  # see above
  #
  if (!(-e $programfile)) {
     $programfile=dirname($0)."/HEADER-program";
  }

}
my $replace=$opt{replace} if defined $opt{replace};
my $show=$opt{show} if defined $opt{show};

#
# Usage
#
if ($#ARGV != 0) {
   die "Usage: $0 file.f90\nOptions:\n --headerfile=blah\n --replace\n";
}

#
# utility to get maximum in hash table
#
sub max_contrib (\%) {
    my $hash = shift;
    keys %$hash;       # reset the each iterator

    my ($large_key, $large_val) = each %$hash;

    while (my ($key, $val) = each %$hash) {
        if ($val > $large_val) {
            $large_val = $val;
            $large_key = $key;
        }
    }
    return ($large_key, $large_val);
}

#
# utility to get maximum length of key in hash table
#
sub max_key_len (\%) {
    my $hash = shift;
    keys %$hash;
    my $maxlen = 0;
    for $key ( keys %$hash ) {
       my $len = length($key);
       if ($len > $maxlen) {
          $maxlen = $len;
       }
    }
    return $maxlen;
}

#
# get the main author using git blame
#
sub getauthor {
   my $file = shift;
   my @blamelist = `git blame $file`;
   my %ncontribs = ();
   foreach (@blamelist) {
      my ($author) = m/.*\((.+\p{IsAlpha}+)\s+\d\d\d\d/;
      #printf("GOT AUTHOR $author\n");
      if ( length($author) > 0 ) {
         $ncontribs{$author}++;
      }
   }
   my ($authorm, $maxval) = max_contrib(%ncontribs);
   my $author = $authorm;
   # add names to author list if they have contributed
   # more than 20% of the code
   #for my $key ( keys %ncontribs ) {
   #   my $value = $ncontribs{$key};
   #   #print "$key => $value\n";
   #   if ( !($key =~ m/$authorm/) && $value > 0.2*$maxval ) {
   #      $author = "$author, $key";
   #   }
   #}
   return $author;
}

#
# get description from old-style header
#
sub get_descript {
   my $file = shift;
   my $inheader = 0;
   my $descript = '';
   my $n = 0;
   my $indescript = 0;
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      #
      # a header is anything between !+ and !+ in the file
      #
      if ( m/\!\+\s+$/ ) {
         $inheader = !$inheader;
         if ($inheader) {
            $newstyle = 0;
            $n++;
         }
      }
      #
      # header ends when a line that is not a comment
      # is encountered
      #
      if ( !m/\!.*/ ) { $inheader = 0 };

      if ($inheader && !m/\!\+/ && !m/\!-+/ && ($n == 1) ) {
         $descript = "$descript$_";
      }
   }
   # remove last empty line
   $descript =~ s/\n\!\s*$/\n/;
   return $descript;
}
#
# get description from free-form program file headers
#
sub get_program_descript_free {
   my $file = shift;
   my $indescript=0;
   my $descript='';
   my $text='';
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      while ( m/^\!.*/ ) {;
         $line = $_;
         if ( $line =~ m/!\s*(This program.*|program to.*)/ and
              not $line =~ m/This program is part of/) {
            $text=$1;
            $indescript = 1;
         }
         if ( m/^\!-+\s*$/ or m/^\!\+.*/ ) {
            $indescript = 0;
         }
         if ($indescript==1) {
            $descript = "$descript$_";
         }
         $_=<$fh>;
      }
   }
   # add DESCRIPTION keyword
   $descript =~ s/!\s*$text/!  $text/;
   # remove last empty line
   $descript =~ s/\n\!\s*$/\n/;
   #print "GOT DESCRIPTION=$descript";
   return $descript;
}
#
# get description from free-form program files
#
sub get_descript_free {
   my $file = shift;
   my $indescript=0;
   my $descript='';
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      while ( m/^\!.*/ ) {;
         $line = $_;
         if ( $line =~ m/!-+\s$/ or m/!\s*$/ ) {
            $indescript = 1;
            $_=<$fh>;
         }
         if ( m/^\!-+\s*$/ or m/^\!\+.*/ or m/^\!\s*\:.*/ or (not m/^!/) ) {
            $indescript = 0;
         }
         if ($indescript==1) {
            $descript = "$descript$_";
         }
         $_=<$fh>;
      }
   }
   # remove last empty line
   $descript =~ s/\n\!\s*$/\n/;
   #print "GOT DESCRIPTION=$descript";
   return $descript;
}

#
# get description from free-form module descriptions
#
sub get_module_descript_free {
   my $header = shift;
   my $descript = '';
   # match everything between !+ and !+
   if ($header =~ m/\!\+\s*\n(.*?\n)\!\+\s*/ ) {
      $descript=$1;
   # match everything starting with a blank comment line
   # and ending with another tag e.g. REFERENCES:
   } elsif ($header =~ m/^!\s*\n(.*?\n)^!(\s*?\n)!\s*[A-Z|a-z|\:]+\:/sm ) {
      $descript=$1;
      $descript =~ s/!\s+(\w)/\! $1/;  # single indentation in description text
      #print "GOT DESCRIPT=$descript\n";
   }
   return $descript;
}


#
# get header
#
sub get_header {
   my $file = shift;
   my $inheader = 1;
   my $header = '';
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      # header is defined as lines that are comments or blank
      # and include the top-level module or program definition statement
      if ( !m/!.*|^\s*$|^\s*module\s+\w.*$|\s*program\s+\w.*$/ ) { $inheader = 0; }
      if ($inheader==1) {
         $header = "$header$_";
      }
   }
   #print "GOT HEADER =$header\nDONE\n";
   return $header;

}

#
# get content of header field, e.g. DESCRIPTION:
#
sub get_field {
   my $fieldname = shift;
   my $header = shift;
   my $field='';
   $fieldname = ucfirst($fieldname);
   #print "finding $fieldname in $header";
   # following two match case where fieldname is NOT preceded by colon
   if ($header =~ m/(\!\s*$fieldname:.*?\n)\!\s*\n\!\s*[\:|A-Z|a-z]+\:/si ) {
      $field=$1;
      $field =~ s/$fieldname/\:$fieldname/i;  # add preceding colon
   } elsif ($header =~ m/(\!\s*$fieldname:.*?\n)\!\+/si ) {
      $field=$1;
      $field =~ s/$fieldname/\:$fieldname/i;  # add preceding colon
   # following two match case where fieldname is preceded by colon
   } elsif ($header =~ m/(\!\s*\:$fieldname:.*?\n)\!\s*\n\!\s*[\:|A-Z|a-z]+\:/si ) {
      $field=$1;
   } elsif ($header =~ m/(\!\s*\:$fieldname:.*?\n)\!\+/si ) {
      $field=$1;
   }
   # correct capitalisation
   $field =~ s/$fieldname/$fieldname/i;
   # correct indentation
   $field =~ s/!\s*\:/! \:/;

   if ($field =~ m/DESCRIPTION/i ) {
      $field =~ s/^! \:$fieldname\:\s*\n//i;  # do not add the word "DESCRIPTION"
      $field =~ s/\:$fieldname\://i;  # do not add the word "DESCRIPTION"
      $field =~ s/^!\s+(\w)/\! $1/;  # single indentation in description text
   }
   #print "GOT $fieldname=$field\nDONE\n";
   return $field;
}

#
# get current git $Id$ tag
#
sub get_id {
   my $header = shift;
   my $id='$Id$';
   if ($header =~ m/(\$Id.*?\$)/ ) {
      $id=$1;
   }
   return $id;
}

#
# get Usage from line matching "Usage:"
#
sub get_usage {
   my $file = shift;
   my $programname = $file;
   $programname =~ s/\..90//;
   $programname =~ s/.*\///;
   my $usage = '';
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      if ( m/sage:(.*)[\'|\"]/ ) {
         $usage=$1;
      }
   }
   $usage =~ s/\'\/\/trim\((.*?)\)\/\/\'/$programname/;
   $usage =~ s/\s+$//;
   return $usage;
}

#
# get module or program name from the Fortran 90 statement
#
sub get_module_name {
   my $file = shift;
   my $type = '';
   my $name = '';
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      if ( m/^\s*(module|program)\s(\w+)\s+$/ ) {
         ($type,$name) = ($1,$2);
      }
   }
   return ($type,$name);
}

#
# get module variables and their descriptions from
# the write_inopt calls
#
sub get_module_vars {
   my $file = shift;
   my %vars = ();
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      if ( m/\s*call write_inopt/ and not m/!/) {
         #print "MATCHING $_";
         my ($key)  = m/write_inopt\(.*\s*\,\s*\'(\w+)\'\s*\,\s*/;
         my ($value) = m/write_inopt\(.*\s*\,\s*\'\w+\'\s*\,\s*\'(.*)\'/;
         if  ( length($value) > 0 ) {
            $vars{$key} = $value;
         }
      }
   }
   #for my $key ( keys %vars ) {
   #  print "GET_VARS:  $key -- $vars{$key}\n";
   #}
   return %vars;
}

#
# get module dependencies from use statements
#
sub get_dependencies {
   my $file = shift;
   my %deps = ();
   open my $fh, '<', $file or die "Can't open $file\n";
   while ( <$fh> ) {
      if ( m/^\s*use\s\w+[\s,\,]/ ) {
         #print "MATCHING $_";
         my ($key)   = m/use\s(\w+)/;
         $deps{$key}++;
      }
   }
   return %deps;
}


#
# subroutine that writes the new header
#
sub write_module_header {
   use Text::Wrap;
   my $file = shift;
   my $author = shift;
   my $descript = shift;
   my $modulename = shift;
   my $refs = shift;
   my $id = shift;
   my $gen = shift;
   my $vartmp = shift;
   my %vars = %$vartmp;
   my $depstmp = shift;
   my %deps = %$depstmp;
   my $usage = shift;
   open my $fh, '<', $file or die "Can't open $file\nSpecify location as follows:\n\n$0 --headerfile=../../scripts/$file\n\n";
   while ( <$fh> ) {
      if ( m/OWNER:/i ) {
         print "! :Owner: $author\n";
      } elsif ( m/DESCRIPTION:/i ) {
         print "$descript";
      } elsif ( m/MODULE:/i ) {
         print "module $modulename\n";
      } elsif ( m/PROGRAM:/i ) {
         print "program $modulename\n";
      } elsif ( m/REFERENCES:/i ) {
         print "$refs";
      } elsif ( m/USAGE:/i ) {
         print "! :Usage:$usage\n";
#      } elsif ( m/\$Id/ ) {
#         print "!  $id\n";
      } elsif ( m/RUNTIME PARAMETERS:/i ) {
         if ( (keys %vars) > 0 ) {
            print "! :Runtime parameters:\n";
            my $len = max_key_len(%vars);
            for my $key ( sort keys %vars ) {
#              print "!  $key -- $vars{$key}\n";
              printf("!   - %-*s : \*%s\*\n",$len,$key,$vars{$key});
            }
         } else {
            print "! :Runtime parameters: None\n";
         }
      } elsif ( m/DEPENDENCIES:/i ) {
         if ( (keys %deps) > 0 ) {
            my $n = 0;
            my $string = ':Dependencies: ';
            for my $key ( sort keys(%deps) ) {
#            for my $key ( sort { $deps{$b} <=> $deps{$a} } keys(%deps) ) {
              if ($n == 0) {
                 $string = sprintf("%s%s",$string,$key);
              } else {
                 $string = sprintf("%s, %s",$string,$key);
              }
              $n++;
            }
            $string = "$string\n";
            my @tmp = split(/ /,$string);
            print wrap("! ",'!   ',@tmp);
         } else {
            print "! :Dependencies: None\n";
         }
      } elsif ( m/GENERATED:/i ) {
         if (length($gen) > 0) {
            print "!\n$gen";
         }
      } else {
         my $thisyear = 1900 + (localtime)[5];;
         $_ =~ s/THISYEAR/$thisyear/;
         print $_;
      }
   }
}

sub parsefile {
   my $file = shift;
   my $headerfile = shift;
   #print stderr "Parsing $file\n";
   my $author = getauthor($file);
   my ($filetype,$modulename) = get_module_name($file);
   my $header = get_header($file);
   my $descript = get_field('Description',$header);
   my $refs = get_field('References',$header);
   my $gen = get_field('Generated',$header);
   if (length($refs)==0) {
      $refs = "! :References: None\n";
   }
   my $id = get_id($header);
   my $inheader = 0;
   if ( $filetype =~ m/module/ ) {
      if (length($modulename) > 0) {
         my $usage = '';
         my %vars = get_module_vars($file);
         my %deps = get_dependencies($file);

         # find old-style description if current description is empty
         if (length($descript)==0) {
            $descript=get_module_descript_free($header);
            if (length($descript)==0) {
               $descript="!  $modulename\n";
            } else {
               $descript="$descript";
            }
         }
         write_module_header( $headerfile, $author, $descript, $modulename,
                              $refs, $id, $gen, \%vars, \%deps, $usage );
         $inheader = 1;
      }
   } elsif ( $filetype =~ m/program/ ) {
      if (length($modulename) > 0) {
         if (length($descript)==0) {
            $descript=get_module_descript_free($header);
         }
         if (length($descript)==0) {
            $descript = get_descript($file);
            if (length($descript)==0) {
               $descript = get_descript_free($file);
            }
            if (length($descript)==0) {
               $descript="! $modulename\n";
            } else {
               $descript="! $descript";
            }
         }
         my $usage = get_usage($file);
         if (length($usage)==0) {
            $usage = " $modulename [no arguments]"
         }
         my %vars = (); #get_module_vars($file);
         my %deps = get_dependencies($file);
         write_module_header( $programfile, $author, $descript, $modulename,
                              $refs, $id, $gen, \%vars, \%deps, $usage );
         $inheader = 1;
      }
   }
   if ($replace) {
      open my $fh, '<', $file or die "Can't open $file\n";
      while ( <$fh> ) {
         # end of header is when no longer a comment
         # but includes the "module blah" or "program blah" statement
         if ( !m/!.*|^\s*$|^\s*module\s+\w.*$|\s*program\s+\w.*$/ ) { $inheader = 0; }
         if ($inheader==0) {
            print $_;
         }
      }
   }
}

my $file = @ARGV[0];

if ( $file =~ m/.*\.[f,F]90/ ) {
   if ($show) {
      print get_header($file);
   } else {
     parsefile($file,$headerfile);
   }
} else {
   print stderr "skipping $file\n";
}
