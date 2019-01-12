#!/usr/bin/env perl

#    deduplicate FASTA file by sequence name
#    Only first sequence of every name gets written
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


my %glob_ids;

my $fin;

if(scalar(@ARGV)){
  my $infile=shift @ARGV or err("Please supply a contig fasta file on the command line.\n");
  (-r $infile and !-d _ and -s _) or err("'$infile' is not a readable non-empty FASTA file\n");
  $fin = Bio::SeqIO->new(-file => $infile, -format=>'fasta');
}else{
  $fin = Bio::SeqIO->new(-fh => \*STDIN, -format=>'fasta');
}

while (my $seq = $fin->next_seq) {
  next if(exists $glob_ids{$seq->id});
  $glob_ids{$seq->id} = 1;
  print ">",$seq->id," ",$seq->desc,"\n",$seq->seq,"\n",
}

sub err {
  print STDERR @_;
  exit(2);
}
