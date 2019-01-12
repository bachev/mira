#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Copyright 2008, by Eurofins MWG/Medigenomix Ebersberg,  All rights reserved.
# Name:      454pairedEnd2caf.pl 
# Author(s): Dr. Jacqueline Weber-Lehmann 
#
#     Status.......: Preliminary
#     Status Date..: 01.05.2008
#     Purpose......: transform fasta 454 paired read info to CAF format
#
# Changes   Date(DD.MM.YYYY)   NAME   DESCRIPTION
#
#-------------------------------------------------------------------------------

my ($script) = ($0 =~ m|([^/]*)$|);

my $fasta_file;
my $qual_file;
my $vector_file;

my $minmatch = 8;
my $minscore = 10;

my $minreadlength = 20;

#-------------------------------------------------------------------------------

my $usage = "$script 
  -i <filename>      fasta file of 454 reads
  -q <filename>      quality file of 454 reads
  -v <filename>      fasta file of spacer sequence
  -minmatch <0..>    minimum match length for cross_match 
                     default: 8    
  -minscore <0..>    minimum match score for cross_match
                     default: 10
  -minreadlength <0..> minimum length of forward and reverse reads
                     default: 20
";

#-------------------------------------------------------------------------------
# read parameters
#------------------------------------------------------------------------------- 

while ( @ARGV ) {
  $arg = shift;
  if    ( $arg eq "-i" ) { $fasta_file  = shift; } # enter 454 read fasta file here
  elsif ( $arg eq "-q" ) { $qual_file = shift; }
  elsif ( $arg eq "-v" ) { $vector_file = shift; } # enter fasta file of vector sequence
  elsif ( $arg eq "-minmatch" ) { $minmatch = shift; }
  elsif ( $arg eq "-minscore" ) { $minscore = shift; }
  elsif ( $arg eq "-minreadlength" ) { $minreadlength = shift; }
  else { die "Unknown option $arg\n\n$usage"; }
}

#------------------------------------------------------------------------------- 

die "files missing" unless ( -f $fasta_file || -f $qual_file || -f $vector_file );

my $caf_file = $fasta_file;
$caf_file =~ s/\.fasta/\.caf/;

#------------------------------------------------------------------------------- 
# Running cross_match

printf STDERR "Running cross_match.";

system("cross_match.manyreads $fasta_file $vector_file -minmatch $minmatch -minscore $minscore -screen > screen.out");

printf STDERR "Running cross_match finished.";

my $screen_fasta_file  = $fasta_file . ".screen";

die "screenfile $screen_fasta_file missing" unless ( -f $screen_fasta_file );

#------------------------------------------------------------------------------- 

open(READ, $fasta_file);

while (<READ>) {
  if ( /^>(\S+)/ ) {
    $name = $1;
    $Read{$name}{seq} = "";
    next;
  }
  $Read{$name}{seq} .= $_;
}

close(READ);

#------------------------------------------------------------------------------- 

open(QUAL, $qual_file);

while (<QUAL>) {
  if ( /^>(\S+)/ ) {
    $name = $1;
    $Read{$name}{qual} = "";
    next;
  }
  $Read{$name}{qual} .= $_;
}

close(READ);

#------------------------------------------------------------------------------- 

my $name;
my $len;
my $ql;
my $qr;

open(CAF, ">$caf_file");
open(SCREEN, $screen_fasta_file) || die "cannot open $screen_fasta_file";

while (<SCREEN>) {
  chomp;
  
  if ( /^>/ ) {
    
    if ( /^>(\S+)\s+LEN\=(\d+)\s+QL\=(\d+)\s+QR\=(\d+)/ ) {
      
      if ( $sequence ne "" ) {

	if ( $sequence =~ /^([^X]*)([X]{1,})([^X]*)$/ ) {
	  $ql_read = $1;
	  $spacer  = $2;
	  $qr_read = $3;
	  
	  
          $length_fwd_read = length($ql_read) - $ql + 1;
          $length_rvs_read = $qr - ( $len-length($qr_read)+1 ) + 1;
	  
          if ( $length_fwd_read >= $minreadlength && $length_rvs_read >= $minreadlength ) {
	    
	    # reverse read
	    
	    $readname = $name . ".r";
	    
            printf CAF "Sequence : %s\n", $readname;
            printf CAF "Is_read\n";
            printf CAF "Padded\n";
            printf CAF "Strand Reverse\n";
            printf CAF "Template \"%s\"\n", $name;
            printf CAF "Insert_size 2000 3000\n";
	    
	    ($rvs_ql, $rvs_qr) = &reverseClipping($ql, length($ql_read), $len);
            printf CAF "Clipping QUAL %d %d\n", $rvs_ql, $rvs_qr;
	    
            printf CAF "Align_to_SCF 1 %d 1 %d\n", $len, $len;
            printf CAF "Tag MINF 1 1 \"ST=454GS\"";
            printf CAF "\n\n";
	    
            printf CAF "DNA : %s\n", $readname;
            printf CAF "%s\n\n", &reverseSequence($Read{$name}{seq});
	    
            printf CAF "BaseQuality : %s\n", $readname;
            printf CAF "%s\n\n", &reverseQuality($Read{$name}{qual});
	    
	    
	    # forward read

	    $readname = $name . ".f";
	    
	    printf CAF "Sequence : %s\n", $readname;
            printf CAF "Is_read\n";
            printf CAF "Padded\n";
            printf CAF "Strand Forward\n";
            printf CAF "Template \"%s\"\n", $name;
            printf CAF "Insert_size 2000 3000\n";
	    
            # ($rvs_ql, $rvs_qr) = &reverseClipping($len-length($qr_read)+1, $qr, $len);
	    
            printf CAF "Clipping QUAL %d %d\n", $len-length($qr_read)+1, $qr;
            printf CAF "Align_to_SCF 1 %d 1 %d\n", $len, $len;
            printf CAF "Tag MINF 1 1 \"ST=454GS\"";
            printf CAF "\n\n";
            printf CAF "DNA : %s\n", $readname;
            printf CAF "%s\n\n", $Read{$name}{seq};
            printf CAF "BaseQuality : %s\n", $readname;
            printf CAF "%s\n\n", $Read{$name}{qual};
	    
	  }
	}
      }
      
      $name = $1;
      $len  = $2;
      $ql   = $3;
      $qr   = $4;
      $sequence = "";
    }
    else {
      die "Wrong format in $screenfile\n";
    }
    next;
  }
  
  $sequence .= $_;

}

if ( $sequence ne "" ) {
  if ( $sequence =~ /^([^X]*)([X]{1,})([^X]*)$/ ) {
    $ql_read = $1;
    $spacer  = $2;
    $qr_read = $3;
    
    $length_fwd_read = length($ql_read) - $ql + 1;
    $length_rvs_read = $qr - ( $len-length($qr_read)+1 ) + 1;
    
    if ( $length_fwd_read >= $minreadlength && $length_rvs_read >= $minreadlength ) {
      
      # reverse read
      
      $readname = $name . ".r";
      
      printf CAF "Sequence : %s\n", $readname;
      printf CAF "Is_read\n";
      printf CAF "Padded\n";
      printf CAF "Strand Reverse\n";
      printf CAF "Template \"%s\"\n", $name;
      printf CAF "Insert_size 2000 3000\n";
      
      ($rvs_ql, $rvs_qr) = &reverseClipping($ql, length($ql_read), $len);
      printf CAF "Clipping QUAL %d %d\n", $rvs_ql, $rvs_qr;
      
      printf CAF "Align_to_SCF 1 %d 1 %d\n", $len, $len;
      printf CAF "Tag MINF 1 1 \"ST=454GS\"";
      printf CAF "\n\n";
      
      printf CAF "DNA : %s\n", $readname;
      printf CAF "%s\n\n", &reverseSequence($Read{$name}{seq});
      
      printf CAF "BaseQuality : %s\n", $readname;
      printf CAF "%s\n\n", &reverseQuality($Read{$name}{qual});
      
      
      # forward read
      
      $readname = $name . ".f";
      
      printf CAF "Sequence : %s\n", $readname;
      printf CAF "Is_read\n";
      printf CAF "Padded\n";
      printf CAF "Strand Forward\n";
      printf CAF "Template \"%s\"\n", $name;
      printf CAF "Insert_size 2000 3000\n";
      
      printf CAF "Clipping QUAL %d %d\n", $len-length($qr_read)+1, $qr;
      printf CAF "Align_to_SCF 1 %d 1 %d\n", $len, $len;
      printf CAF "Tag MINF 1 1 \"ST=454GS\"";
      printf CAF "\n\n";
      printf CAF "DNA : %s\n", $readname;
      printf CAF "%s\n\n", $Read{$name}{seq};

      printf CAF "BaseQuality : %s\n", $readname;
      printf CAF "%s\n\n", $Read{$name}{qual};
      
    }    
  }
}



close(SCREEN);

close(CAF);

printf "output written to $caf_file\n";


sub reverseQuality {
  my ($quality) = @_;
  $quality =~ s/\n//g;
  my $rvs_quality = "";
  my @q = ();

  @q = split / /, $quality;
  $n = $#q;  

  for ( $i=0; $i<$n ; $i++) { 
    $rvs_quality .= $q[$n-$i] . " ";
  }

  $rvs_quality .= $q[0];

  return $rvs_quality;
}

sub reverseClipping {
  my ($ql, $qr, $len) = @_;

  $ql_new = $len - $qr + 1;

  $qr_new = $len - $ql + 1;

  return ( $ql_new, $qr_new );

  

}

sub reverseSequence {
  my $sequence = shift;
  $sequence =~ s/\n//g;
  $sequence =~ tr/ACGTUYRMKSWHBVDNacgtuyrmkswhbvdn/TGCAARYKMSWDVBHNtgcaarykmswdvbhn/d;

  $sequence = reverse($sequence);
  return $sequence;
}
