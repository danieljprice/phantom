#!/usr/bin/perl
#
# @(#) Perl script to remove dimension() statements in Fortran 90 code
# @(#) Also aligns the intent() :: statements so they line up neatly
#
# Examples:
# 
# real, intent(in) :: x,y,z
# real, intent(out) :: u,v
# real, intent(inout) :: a,b,c
#
# becomes
#
# real, intent(in)    :: x,y,z
# real, intent(out)   :: u,v
# real, intent(inout) :: a,b,c
#
# For optional arguments, we align the optional keyword with the variable name, i.e.:
#
# real, intent(in) :: x,y,z
# real, intent(out) :: u,v
# real,    intent(inout), optional :: a,b,c
# integer, intent(in), optional :: i,j,k
#
# becomes
#
# real,    intent(in)  :: x,y,z
# real,    intent(out) :: u,v
# real,    intent(inout), optional :: a,b,c
# integer, intent(in),    optional :: i,j,k
#
use Text::Balanced qw<extract_bracketed>;
my $type_pattern = qr/(?:real|integer|logical|character|complex|double\s+precision|type|class)\b/i;
my @file = '';
while (<STDIN>) {
  my @block = ();
  my $gotblock = 0;
  #
  # match variable declaration block
  #
  while ( m/^(\s*)($type_pattern.*)\:\:\s*(.*)/ ) {
     my $pad = $1; # spaces
     my $argline = $2; # everything to left of  ::
     my $varline = $3; # everything to right of ::
     $argline =~ s/(\S)$/$1 /g; # ensure there is a space between variables and ::
     my $dims = '';
     my $newargline = '';
     my $newvarline = '';
     my $newline = '';
     my $comment = '';
     # only apply to intent/dimension lines and do not touch line
     # if it contains a line continuation
    if ( ($argline =~ m/intent|dimension/g) and not ($varline =~ m/\&|\=/g) ) {
        if ($argline =~ m/dimension\s*(\(.*\))/) {
           my $tmp = $1;
           $dims = extract_bracketed( $tmp, '()' );
        }

        # remove any dimension(...) attribute from the argument list,
        # regardless of exact spacing, so we don't leave a bare (:,:)
        # behind in $argline.
        $argline =~ s/\s*,\s*dimension\s*\([^)]*\)//ig;  # trailing ", dimension(:,:)"
        $argline =~ s/\bdimension\s*\([^)]*\)\s*,\s*//ig; # leading  "dimension(:,:),"
        
        # make sure there are spaces between variables
        $argline =~ s/(.*?),([a-z]+)/$1, $2/g;

        # get remaining arguments from splitting argline
        my @args = split(', ',$argline);
        my $count=0;
        foreach $arg (@args) {
           if ( !($arg =~ m/^\s*$/) ) { # skip blank arguments
              $arg =~ s/^\s+//;
              $arg =~ s/\s+$//;
              if ($count==0) {
                 $newargline = "$arg";
              } else {
                 $newargline = "$newargline, $arg";
              }
              $count++;
           }
        }
        
        if ( $varline =~ m/(.*?)(\s*\!.*)/ ) {
           $varline = $1;
           $comment = $2;
        }
        # split variable list on commas that are NOT inside parentheses,
        # so that vxyzu_add(:,:) stays as a single entry
        my @vars = ();
        my $current = '';
        my $depth = 0;
        foreach my $ch (split //, $varline) {
           if ($ch eq '(') { $depth++; }
           elsif ($ch eq ')' && $depth > 0) { $depth--; }
           if ($ch eq ',' && $depth == 0) {
              push @vars, $current;
              $current = '';
              next;
           }
           $current .= $ch;
        }
        push @vars, $current if $current =~ /\S/;

        $count=0;
        foreach $var (@vars) {
           $var =~ s/\s*$//g; # strip spaces at end of variable
           my $var_with_dims = $var;
           # Only append extracted dimensions if the variable does NOT
           # already have an explicit shape. This avoids cases like
           # vxyzu_add(:,:)(:,:) when one var already has (:,:).
           if ( $dims ne '' && $var !~ /\(/ ) {
              $var_with_dims = "$var$dims";
           }
           if ($count==0) {
              $newvarline = $var_with_dims;
           } else {
              if ($var =~ m/\w+/) {
                 # only add space if argument contains word character
                 $newvarline = "$newvarline,$var_with_dims";
              } else {
                 $newvarline = "$newvarline,$var_with_dims";
              }
           }
           $count++;
        }
        $newline = "$pad$newargline :: $newvarline$comment\n";
     } else {
        # make sure there are spaces between variables
        $argline =~ s/(.*?),([a-z]+)/$1, $2/g;

        # regurgitate line, with spacing adjusted
        $newline = "$pad$argline\:\: $varline\n";
     }
     push @block,$newline;
     $_ = <>;
     $gotblock = 1;
  }
  if ($gotblock==1) {
     # normalise declarations then pad to consistent lengths
     my @newblock = ();
     my @parsed = ();
    my $padlen  = 0;
    my $padleni = 0;
    my $padleni_nonopt = 0;

     foreach $line (@block) {
       if ( $line =~ m/(.*)(intent.*\))(.*)(\:\:.*)/ ) {
          my $vari = $1;
          my $in   = $2;
          my $ex   = $3;
          my $rest = $4;
          my $has_optional = 0;

          $vari =~ s/\s+$//;
          $in   =~ s/\s+$//;
          $ex   =~ s/\s+$//;

          $has_optional = 1 if $vari =~ s/\s*,\s*optional\s*,?//ig;
          $has_optional = 1 if $ex   =~ s/\s*,\s*optional\s*//ig;

          my @extras = ();
          if ( $ex =~ /\S/ ) {
             my @parts = split(/,/, $ex);
             foreach my $part (@parts) {
                $part =~ s/^\s+//;
                $part =~ s/\s+$//;
                next unless length $part;
                push @extras, $part;
             }
          }

          $vari =~ s/\s+$//;
          my $vari_clean = $vari;
          $vari_clean =~ s/,\s*$//;
          if (@extras) {
             my $extras_str = join(', ', @extras);
             if ( $vari_clean =~ /\S/ ) {
                $vari_clean = "$vari_clean, $extras_str";
             } else {
                $vari_clean = $extras_str;
             }
          }
          $vari = $vari_clean;
          $vari =~ s/\s+$//;
          if ( $vari =~ /\S/ and $vari !~ /,\s*$/ ) {
             $vari = "$vari,";
          }

          $rest =~ s/^\s*::/::/;

          my $base = $in;

          my %info = (
             vari     => $vari,
             base     => $base,
             rest     => $rest,
             optional => $has_optional,
          );
          push @parsed, \%info;

          my $len  = length($vari);
          my $leni = length($base);
          if ($len  > $padlen ) { $padlen  = $len; }
          if ($leni > $padleni) { $padleni = $leni; }
          if (!$has_optional && $leni > $padleni_nonopt) { $padleni_nonopt = $leni; }
       } else {
          push @parsed, { raw => $line };
       }
     }

    $padleni_nonopt = $padleni if $padleni_nonopt == 0;

    my $optlen = 0;
    foreach my $entry (@parsed) {
      next if exists $entry->{raw} || !$entry->{optional};
       my $base_full = sprintf("%-*s %-*s",$padlen,$entry->{vari},$padleni,$entry->{base});
       (my $base_trim = $base_full) =~ s/\s+$//;
       my $baselen = length($base_trim);
       if ($baselen > $optlen) { $optlen = $baselen; }
     }

     my @nonopt = grep { !exists $_->{raw} && !$_->{optional} } @parsed;
     my @opt    = grep { !exists $_->{raw} &&  $_->{optional} } @parsed;
     my @param  = grep {  exists $_->{raw} && $_->{raw} =~ /\bparameter\b/i } @parsed;
     my @raw    = grep {  exists $_->{raw} && $_->{raw} !~ /\bparameter\b/i } @parsed;

     foreach my $entry (@param, @nonopt, @opt, @raw) {
       my $newline = '';
       if ( exists $entry->{raw} ) {
          $newline = $entry->{raw};
       } else {
          my $vari = $entry->{vari};
          my $base = $entry->{base};
          my $rest = $entry->{rest};
          my $intent_width = $entry->{optional} ? $padleni : $padleni_nonopt;
          my $base_str = sprintf("%-*s %-*s",$padlen,$vari,$intent_width,$base);
          (my $base_trim = $base_str) =~ s/\s+$//;
          if ( $entry->{optional} ) {
             my $pad_spaces = $optlen - length($base_trim);
             $pad_spaces = 0 if $pad_spaces < 0;
             my $opt_str = ', ' . (' ' x $pad_spaces) . 'optional';
             $newline = sprintf("%s%s %s\n",$base_trim,$opt_str,$rest);
          } else {
             $newline = sprintf("%s %s\n",$base_str,$rest);
          }
       }
       push @newblock, $newline;
     }

     push @file, @newblock;
     $gotblock = 0;
  }
  push @file, $_;
}
print @file;
