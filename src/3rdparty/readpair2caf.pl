#!/usr/bin/perl -w
use strict;
# configure this to your specific environment
use lib '/package/bioperl/default';
use Getopt::Long qw(:config auto_version);
use Bio::SeqIO;
use Bio::Seq;
use Pod::Usage;

=head1 NAME

readpair2caf.pl - convert FLX 454 long read pair concatenated fasta file to CAF output

=head1 SYNOPSIS

readpair2caf,pl <read pair fasta file> <read pair quality file> <adaptor blastoutput> <caf outputfile>

=head1 DESCRIPTION

Generates a CAF file of read pair information suitable for sequence
assembly programs such as MIRA

=item B<-v>

verbose output

=item B<-seqsize>

size cut off for forward and reverse reads

=item B<-c>

size cut off for blast alignment

=head1 EXAMPLE

1) Generate fasta, qual and xml files using of the readpair sff files  
sffinfo test.sff | sffinfo2mirafiles.tcl -project test
generates:
test_in.454.fasta
test_in.454.fasta.qual
test_traceinfo_in.454.xml

the generated fasta files should have lowercase and uppercase letters correspinding to low quality/high quality respectively.

2) with the readpair fasta file make a blast database of the fasta
formatdb -i test_in.454.fasta -p F -n test
generates:

test.nhr
test.nin
test.nsq
formatdb.log

3) Blast the sequences sequence vs the adaptor sequence using the -m8 option
blastall -p blastn -d adaptor -i test_in.454.fasta -o adaptor.blo -m8 -e 1e-10 -S 1
generates:

adaptor.blo

which give the coordinated about where the adaptor starts and ends

The adaptor sequence is:

>FLX Adaptor
GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC


4) Generate the CAF format file

readpair2caf.pl test_in.454.fasta test_in.454.fasta.qual adaptor.blo test.caf

=head1 AUTHOR

Written by Stephen Taylor, Computational Biology Research Group, Oxford 

=head1 REPORTING BUGS

Report bugs to genmail@molbiol.ox.ac.uk

=cut

my $help;
my $man;

# global hashes
my %seq_info;
my %qual_info;
my %hit;

my $INSERT_SIZE=2000;
my $INSERT_SIZE_LOWER=1400;
my $INSERT_SIZE_UPPER=2600;
# how long the crossmatch/blast alignment needs to be before rejection
my $adaptor_align_cutoff=44;
my $seqsize=10;
my $verbose;

my ($sff_fasta_file, $sff_qual_file, $adaptor_blast_file, $caf_output_file)=@ARGV;

GetOptions(    
        'c:i'=>\$adaptor_align_cutoff,	
	'v'=>\$verbose,
	'seqsize:i'=>\$seqsize,
        'h|help'=>\$help,
        'm|man'=>\$man

); 


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
if ($#ARGV<3) {
	pod2usage(1);
}


print STDERR "Loading quality files...\n";
&LoadQual($sff_qual_file);
print STDERR "Done\n";
print STDERR "Loading sequences...\n";
&LoadSeqs($sff_fasta_file);
print STDERR "Done\n";
print STDERR "Reading blast output...\n";
&BLORead($adaptor_blast_file);
print STDERR "Done\n";
print STDERR "Writing forward and reverse reads to CAF file...\n";
WriteFandR($caf_output_file);
print STDERR "Done\n";



exit;

# extract forward and reverse reads

sub BLORead {
	my ($adaptor_blast_file)=@_;
	my $count=0;
	open(BLO,$adaptor_blast_file) or die $!;
	while(<BLO>) {
		my @fields=split(/\t/);
		my $id=$fields[0];
		my $start_match=$fields[6];
		my $end_match=$fields[7];
		my $align_length=$fields[3];

		#print STDERR "$id $start_match $end_match\n";
		# if we have found the adaptor twice in the sequence need to merge it (if its a consecutive hit) 
		# for example
		#   44  0.00 0.00 0.00  E0VUAYN02JXVCL       38    81 (226)    FLX        1    44 (0)  
  		#   44  0.00 0.00 0.00  E0VUAYN02JXVCL       82   125 (182)  C FLX   (0)    44     1  
		if (exists($hit{$id})) {
			if ($verbose) {
				print "COORDS:  start_match=".$start_match." end_match=".$end_match." hit_id_start_match=".$hit{$id}{start_match}." hit_id_end_match=".$hit{$id}{end_match}."\n";
			}
			if ($hit{$id}{end_match} == ($start_match-1)) {
				$hit{$id}{end_match}=$end_match;
				if ($verbose) {					
					print "FIXED:$id hit_id_start_match=".$hit{$id}{start_match}." hit_id_end_match=".$hit{$id}{end_match}."\n";
				}
			} else {
			# not sure what to do with it so delete it from the entries	
				delete($hit{$id});
				if ($verbose) {
					print "DELETED:$id\n";
				}
			}
			next;
			
		}
		$hit{$id}{start_match}=$start_match;
		$hit{$id}{end_match}=$end_match;
		$hit{$id}{align_length}=$align_length;
		if ($verbose) {
			print "BLO:$id $start_match; $end_match; $align_length\n";
		}
	}

}

sub WriteFandR {
	my ($output_file)=@_;
	open(CAF,">$output_file") or die $!;
	my $count=0;	
	foreach my $id (keys %hit) {		 
		my $start_match=$hit{$id}{start_match};
		# the actual pos in the sequence 0 indexed we want to substr to i.e. 0 to start of match (not including the match to adaptor)
		my $l_end_0_index=$start_match-2;
		# the actual pos in the sequence 1 indexed that the alignment ends i.e 1 to start of match (not including the match to the adaptor) 	
		my $l_end_1_index=$start_match-1;
		my $end_match=$hit{$id}{end_match}; 
		# the actual pos in the sequence 0 indexed we want to substr from ie after the match
		my $r_start_0_index=$end_match;
		my $align_length=$hit{$id}{align_length};
		#Must have a good alignment for the adaptor sequence
		if ($align_length < $adaptor_align_cutoff) {
			if ($verbose) {
				print "ADAPTORCUTOFF:$id alignment=$align_length < adaptor_size=$adaptor_align_cutoff..rejecting\n";
			}
			next;
		} 
		my $l=length($seq_info{$id});
		my $l_read=substr($seq_info{$id},0,$l_end_1_index);	
		
		#print "r_start_0 = $r_start_0_index, l = $l\n".$seq_info{$id}."\n";	
		my $r_read=substr($seq_info{$id},$r_start_0_index,$l);
		
		
		# l_seq = seq on left of adaptor
		# r_seq = seq on right hand side of adaptor
		my $l_seq = Bio::Seq->new(-seq => $l_read, -alphabet => 'dna');
		my $r_seq = Bio::Seq->new(-seq => $r_read, -alphabet => 'dna');
		
		my $l_revcom_seq=$l_seq->revcom()->seq();
		my $r_for_seq=$r_seq->seq();
		
		my $l_seq_size=length($l_read);
		my $r_seq_size=length($r_read);

		#reject if length < $seqsize for forward or reverse sequence
		if ($l_seq_size < $seqsize or  $r_seq_size < $seqsize) {
			if ($verbose) {
				print "REJECTING:$id r size=$l f size=$r_seq_size\n"; 
			}
			next;
		}
		
		# print forward and reverse sizes
		if ($verbose) {
			print "READLENGTHS:$id	$l_seq_size	$r_seq_size\n";
		}
		
		# l reads
		my $clip_size_l=$l_seq_size-&GetClipSize($l_read);
		my @l_qual=@{$qual_info{$id}}[0..$l_end_0_index];
		my $l_rev_qual=reverse(@l_qual);
		my $l_qual=join(" ",@l_qual);
		
		if ( $clip_size_l > 0 ) {
			print CAF "DNA : $id.r\n$l_revcom_seq\n\n";
			print CAF "BaseQuality : $id.r\n";
			$l_qual=~s/^\s+//;
			print CAF "$l_qual\n\n";
			print CAF "Sequence: $id.r\n";
			print CAF "Is_read\nStrand Reverse\nPadded\nTemplate \"$id\"\n";
			print CAF "Insert size $INSERT_SIZE_LOWER $INSERT_SIZE_UPPER\n";
			print CAF "Align_to_SCF 1 $l_end_1_index 1 $l_end_1_index\n";
			print CAF "Clipping QUAL 1 $clip_size_l\n";
			print CAF "Tag MINF 1 1 \"ST=454GS\"\n\n";	
			$count++;
		}
		
		my $clip_size_r=$r_seq_size-&GetClipSize($r_for_seq);
		my @rev_qual=@{$qual_info{$id}}[$r_start_0_index..$l-1];
		my $rev_qual=join(" ",@rev_qual);  

		if ($clip_size_r > 0 ) {		
			# reverse reads
			print CAF "DNA : $id.f\n$r_for_seq\n\n";
			print CAF "BaseQuality : $id.f\n";
			# get rid of initial space
			$rev_qual=~s/^\s+//;
			print CAF "$rev_qual\n\n";
			print CAF "Sequence: $id.f\n";
			print CAF "Is_read\nStrand Forward\nPadded\nTemplate \"$id\"\n";
			print CAF "Insert size $INSERT_SIZE_LOWER $INSERT_SIZE_UPPER\n";
			print CAF "Align_to_SCF 1 $r_seq_size 1 $r_seq_size\n";
			print CAF "Clipping QUAL 1 ".$clip_size_r."\n";
			print CAF "Tag MINF 1 1 \"ST=454GS\"\n\n";	
			$count++;

		}
		
		if ($verbose) {
			print "SIZE:$id fseq=$l_seq_size fqual=$#l_qual rseq=$r_seq_size rqual=$#rev_qual\n";
		}
		
	}
	print STDERR "$count read pairs passed.\n";
}

sub LoadSeqs {
	my ($input_file)=@_;
	my $seq_input = Bio::SeqIO->new( '-format' => 'fasta' , -file => $input_file);
	while ((my $seqobj = $seq_input->next_seq())) 
	{    
        	my $id = $seqobj->display_id;
        	my $seq = $seqobj->seq;
		$seq_info{$id}=$seq;
	}       

}

sub GetClipSize {
	my ($seq)=@_;
	if ($seq=~ /([a-z]+)/) {
		#print STDERR "extracted = ".$1;
		#print STDERR "LENGTH ".(length($1))."\n";
		return (length($1));
	}
}

sub LoadQual {
	my ($input_file)=@_;
	my $seq_input = Bio::SeqIO->new( '-format' => 'qual' , -file => $input_file);
	while ((my $seqobj = $seq_input->next_seq())) 
	{    
        	my $id = $seqobj->display_id;
        	my @qual =@{$seqobj->qual()};
		$qual_info{$id}=\@qual;
	}       

}
