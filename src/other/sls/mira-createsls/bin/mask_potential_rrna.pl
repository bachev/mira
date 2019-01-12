#!/usr/bin/env perl

#    Mask potential rRNAs in file
#    Output written to stdout
#
#    Warning: dumb coding, do nout use on very large files
#      TODO: correct that
#
#    Copyright (C) 2015- Bastien Chevreux
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::GFF;


my %glob_seq;
my %glob_desc;
my @glob_seqids;  # seq ids (names) in the order as read in the fasta file

my $infile=shift @ARGV or err("Please supply a contig fasta file on the command line.\n");
(-r $infile and !-d $infile and -s $infile) or err("'$infile' is not a readable non-empty FASTA file");

read_fastaseqs2annotate();

open my $BARRNAP, '-|', "barrnap --kingdom bac --reject 0.0001 --quiet $infile" or err("Please install the barrnap tool.\n");
my $gff = Bio::Tools::GFF->new(-fh => $BARRNAP, -gff_version => 3);
while (my $feat = $gff->next_feature) {
  if(exists $glob_seq{$feat->seq_id}){
    my $flen=$feat->end-$feat->start;
    substr($glob_seq{$feat->seq_id}, $feat->start, $flen) = 'n' x $flen;
  }
}

# things not found by barrnap
for my $sid (@glob_seqids) {
  $glob_seq{$sid} =~ s/GGCAGTTCCCTACTCTCGCATGGGGAGACCCCACACTACCATCGGCGCT/NNN/g;
  $glob_seq{$sid} =~ s/AGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCC/NNN/g;

  $glob_seq{$sid} =~ s/GTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGC/NNN/g;
  $glob_seq{$sid} =~ s/GCCTGGCAGTTCCCTACTCTCGCATGGGGAGACCCCACACTACCATCGGCGCTACGGCGTTTCACTTCTGAGTTCGGCATGGGGTCAGGTGGGACCACCGCGCTAC/NNN/g;

  $glob_seq{$sid} =~ s/GTTCGGCATGGGGTCAGGTGGGACCACCGCGCTAC/NNN/g;
  $glob_seq{$sid} =~ s/GTAGCGCGGTGGTCCCACCTGACCCCATGCCGAAC/NNN/g;

  $glob_seq{$sid} =~ s/TGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCCCCCATGCGAGAGTAGGGAACTGCCAG/NNN/g;
  $glob_seq{$sid} =~ s/CTGGCAGTTCCCTACTCTCGCATGGGGGACCCCACACTACCATCGGCGCTACGGCGTTTCACTTCTGAGTTCGGCA/NNN/g;

  $glob_seq{$sid} =~ s/CCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGC/NNN/g;
  $glob_seq{$sid} =~ s/GCCTGGCAGTTCCCTACTCTCGCATGGGGAGACCCCACACTACCATCGG/NNN/g;

  $glob_seq{$sid} =~ s/CTGATAAAACAGAATTTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTG/NNN/g;
  $glob_seq{$sid} =~ s/CAGGTGGGACCACCGCGCTACTGCCGCCAGGCAAATTCTGTTTTATCAG/NNN/g;

  $glob_seq{$sid} =~ s/AAGGGCGAATTCTGCAGATATCCATCACA/NNN/g;
  $glob_seq{$sid} =~ s/TGTGATGGATATCTGCAGAATTCGCCCTT/NNN/g;

  $glob_seq{$sid} =~ s/CCTGGCGGCAGTAGCGCGGTGGTCCCACCTG/NNN/g;
  $glob_seq{$sid} =~ s/CAGGTGGGACCACCGCGCTACTGCCGCCAGG/NNN/g;

  $glob_seq{$sid} =~ s/CCGGCCGCCATGGCGGCCGCGGGAATTCGATT/NNN/g;
  $glob_seq{$sid} =~ s/AATCGAATTCCCGCGGCCGCCATGGCGGCCGG/NNN/g;

###
# These do not look like regular rRNA, at least NCBI BLAST makes me suspicious
#AAGGGCGAATTCGCGGCCGCTAAATTCAATT
#CCGGCCGCCATGGCGGCCGCGGGAATTCGATT
#ATTCTGGATCCGATACGTAACGCGTCTGCAGCATGCGTGGTACCGAGCTTTCCCTAT
#ATAGAATACTCAAGCTATGCATCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCGAATTCACTAGTGATT
#AAGGGCGAATTCCAGCACACTGGCGGCCGTTACTAGTGGATCCGAGCTCGGTACCAAGCTTGATGCAT
#TAAGCCGAATTCCAGCACACTGGCGGCCGTTACTAGTGGATCCGAGCTCGGTACCAAGCTTGGG
#
}

#write
for my $sid (@glob_seqids) {
  print ">$sid ",$glob_desc{$sid},"\n",$glob_seq{$sid},"\n";
}


sub read_fastaseqs2annotate {
  my $fin = Bio::SeqIO->new(-file => $infile, -format=>'fasta');

  while (my $seq = $fin->next_seq) {
    $glob_seq{$seq->id} = $seq->seq;
    $glob_desc{$seq->id} = $seq->desc;
    push @glob_seqids, $seq->id;
  }
}

sub err {
  print STDERR @_;
  exit(2);
}
