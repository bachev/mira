#!/usr/bin/perl -w 
#
# bin_fasta_on_mid_primers.pl
#
# 	$Id$
#
#
# This perl script bins the reads from a file of fasta reads
# labeled with MID primers according to the MID primer sequence.
# MID primers are used in conjunction with a 454 sequencer on
# multiplexed samples. Samples labled with different MID primers
# can be run simultaneously on the same 454 sequencing run.
#
# This script requires as input
#  -r <fasta file of reads to extract MID primers from>
#  -q <quality file associated with fasta file of reads to extract MID primers from>
#  -m <fasta file of MID primers>
#  -b <base name of project>
#
#
# This script will output .fasta and .qual files for each and
# every MID primer present in both the MID fasta file
# and the fasta file of 454 reads.
# A MID primer must have an exact match at the begining of a read
# to be binned into the output files.
#
# If a read doesn't have an exact match to at least one primer in
# the MID fasta file, that read is put into junk file.

# Written by Greg Harhay, USDA-ARS-MARC gregory.harhay@ars.usda.gov
# United States Department of Agricultural, Agricultural Research Service
# Animal Health Research Unit
# Meat Animal Research Center, Clay Center, NE

# Unconditional Release - This code is now being released at no cost to the public for
# further development. ARS is releasing the code so that interested parties can use it for their
# own needs and purposes. ARS does not foresee providing monetary or technical support to refine,
# adapt, or use this code, and provides no warranty for its use for any purpose.
# ARS does not reserve any rights or interests in the work that may be performed by
# others to refine or adapt it. ARS does reserve the right to continue its own refinement
# of the current version of the code at a later date, should program needs require it.

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Perl;
use Getopt::Std;

my ($fasta_reads, $qual_reads, $fasta_midi, $base, %opts, @fasta_files_out,
   @qual_files_out);

getopts('r:q:m:b:', \%opts);
chomp(%opts);

if ($opts{h} || $opts{'?'}) {
   print "\n\nUsage:  binMidPrimers.pl \n\n";
   print "-r <fasta file of reads to extract MID primers from> \n";
   print
     "-q <quality file associated with fasta file of reads to extract MID primers from>\n";
   print "-m <fasta file of MID primers>\n";
   print "-b <base name of project files this script creates>\n\n";
   print "This perl script bins the reads from a file of fasta reads\n";
   print "labeled with MID primers according to the MID primer sequence.\n";
   print "MID primers are used in conjunction with a 454 sequencer on\n";
   print "multiplexed samples. Samples labled with different MID primers\n";
   die "can be run simultaneously on the same 454 sequencing run.\n\n";
} ## end if ($opts{h} || $opts{...

if ($opts{r}) {
   $fasta_reads = $opts{r};
   chomp $fasta_reads;
}

if ($opts{q}) {
   $qual_reads = $opts{q};
   chomp $qual_reads;
}

if ($opts{m}) {
   $fasta_midi = $opts{m};
   chomp $fasta_midi;
}

if ($opts{b}) {
   $base = $opts{b};
   chomp $base;
}

unless (defined $fasta_reads
   && defined $qual_reads
   && defined $fasta_midi
   && defined $base)
{
   print "\n\nUsage:  bin_fasta_on_mid_primers.pl \n\n";
   print "-r <fasta file of reads to extract MID primers from> \n";
   print
     "-q <quality file associated with fasta file of reads to extract MID primers from>\n";
   print "-m <fasta file of MID primers>\n";
   print "-b <base name of project files this script creates>\n\n";
   print "This perl script bins the reads from a file of fasta reads\n";
   print "labeled with MID primers according to the MID primer sequence.\n";
   print "MID primers are used in conjunction with a 454 sequencer on\n";
   print "multiplexed samples. Samples labled with different MID primers\n";
   die "can be run simultaneously on the same 454 sequencing run.\n\n";
} ## end unless (defined $fasta_reads...

my $seq_reads = Bio::SeqIO->new(-file => $fasta_reads, -format => "fasta");
my $seq_quals = Bio::SeqIO->new(-file => $qual_reads,  -format => "qual");
my $seq_midis = Bio::SeqIO->new(-file => $fasta_midi,  -format => "fasta");

undef my %midi_fasta;
undef my %midi_qual;

# For each  screen sequence, upper case sequence, open up a file to stuff data into - file name is
# concatination of project base name and midi-screen sequence.

while (my $inseq = $seq_midis->next_seq) {
   my $midi_seq = uc($inseq->seq);

   my $out_seq_name = $base . "_midi_" . $midi_seq . ".fasta";
   $midi_fasta{$midi_seq} =
     Bio::SeqIO->new(-file => ">$out_seq_name", -format => "fasta");

   my $out_qual_name = $base . "_midi_" . $midi_seq . ".fasta.qual";
   $midi_qual{$midi_seq} =
     Bio::SeqIO->new(-file => ">$out_qual_name", -format => "qual");

   push(@fasta_files_out, $out_seq_name);
   push(@qual_files_out,  $out_qual_name);

} ## end while (my $inseq = $seq_midis...

# Create a base_name.junk file that has reads that don't match the midi-screen sequence

my $out_junk_name = $base . "_junk.fasta";
my $out_junk_qual = $base . "_junk.fasta.qual";
my $midi_junk = Bio::SeqIO->new(-file => ">$out_junk_name", -format => "fasta");
my $midi_junk_qual =
  Bio::SeqIO->new(-file => ">$out_junk_qual", -format => "qual");

# Loop over file of fasta read file

while (my $read = $seq_reads->next_seq) {

   my $qual     = $seq_quals->next_seq;
   my $read_seq = uc($read->seq);

   my $flag = 0;

   foreach my $midi_seq (keys(%midi_fasta)) {
      $_ = $midi_seq;
      my $count_midi_bases = tr/A-Z//;
      if ($read_seq =~ /^$midi_seq/) {    # does read begin with midi seq

         my $length_seq   = length($read_seq);
         my $new_read_seq = $read->subseq($count_midi_bases + 1, $length_seq);
         my $new_qual_seq = $qual->subqual($count_midi_bases + 1, $length_seq);

         $read -> seq($new_read_seq);
         $qual -> qual($new_qual_seq);

         $midi_fasta{$midi_seq}->write_seq($read);
         $midi_qual{$midi_seq}->write_seq($qual);

         $flag = 1;

      } ## end if ($read_seq =~ /^$midi_seq/)
   } ## end foreach my $midi_seq (keys(...
   if ($flag == 0) {  # seq doesn't match any  midi in midi_screen,  put in junk
      $midi_junk->write_seq($read);
      $midi_junk_qual->write_seq($qual);
   }
} ## end while (my $read = $seq_reads...

## Remove MID bin files that are empty

foreach my $file (@fasta_files_out) {

   if (-z $file) {
      unlink($file);
   }

}

foreach my $file (@qual_files_out) {

   if (-z $file) {
      unlink($file);
   }

}

# 1) upper case bases in sequence
# 2) Loop over all files in screen sequence
# 3) Use index to look for hits of screen sequence against read sequence
# 4) IF there is a hit,remove mid, write new sequence to appropriate file
# 5) IF there is no hit to any screen sequence, stuff that read into a junk file
