#!/usr/bin/perl -w

=head1 SYNOPSIS

caf2aceMiraConsed.pl - A script to use mira-produced caf and ace files with consed.

=head1 USAGE

caf2aceMiraConsed.pl --caf_file <string> --ace_file <string> [--sff_files <string>] [--solexa_fof <string>] [--solexa_paired][--phd_balls <string>] [--default_date <date>] [--consed_bin <string>] 

=head1 INPUT

=head2 --caf_file

CAF file containing the assembly

=head2 --ace_file

ACE file containing the same assembly as the CAF file. Obtained (for example) by running mira program convert_project -f caf -t ace project.caf project. Mandatory.

=head2 --sff_files

sff files containing (some of) the reads included in the caf file. Use commas to separate file names. File should be in a sff_dir of a consed assembly. consed -sff2PhdBall will be used on this files.

=head2 --solexa_files

For the moment, it is not recommended to use that option, since mira silently edit the solexa reads, for the sake of memory. A better solution is to let the script build fake phdballs, using if applicable --solexa_paired.

fof file ("file of file" in consed lingo) containing names of fastq files containing (some of) the reads included in the caf file. consed -solexa2PhdBall will be used on this files. File includes file names (without the path, files are assumed to be in ../solexa_dir). Two filenames on the same line separated by a space means paired ends.

=head2 --solexa_paired

Are solexa reads paired-ended? If stated, the script will assume that solexa reads are paired-end and will expect to have reads ending in /1 and /2, and will link them.

=head2 --phdballs

Existing phd_balls containing reads included in the caf file. Use commas to separate file names. Files should be in a phdball_dir of a consed assembly. These files will be modified to adapt them to the ace file.

=head2 --default_date

Sets the date for the reads. By defaults, sets the default date put by MIRA (Sat Jan 1 11:11:11 2000). You can either specify a particular date or use "current" to use the date at the launch of the program. Be careful that the reads need to have the same date as in the ace file to be read by consed.

=head2 --consed_bin

Location of the consed program. Can be the full path or only the name of the program used if available on the command line. consed is used by default.

=head1 OUTPUT

The main output is a modified ace file, with tags transferred at the end. Phballs are also created/modified.

To be read by consed, tags need to be recognized by consed. To achieve that, put tag definition in a consed-compatible file and link to it in your distribution's .consedrc file. For example, in the consed home folder, add that line to the .consedrc file (or create one if there is none): 

consed.fileOfTagTypes:/path/to/consed/home/mira_tags.conf

mira_tags.conf looks like:

MIRA OrangeRed both yes
SRMr OrangeRed both yes
SRMc OrangeRed both yes
WRMr Orange both yes
WRMc Orange both yes
ORMB orange1 both yes
CRMr palegoldenrod1 both yes
MNRr black both yes
UNSr Yellow both yes
UNSc Yellow both yes
IUPc royalblue3 both yes
MCVc firebrick3 both yes
SROr DarkTurquoise both yes
SROc DarkTurquoise both yes
SAOr SeaGreen both yes
SAOc SeaGreen both yes
SIOr PaleGreen both yes
SIOc PaleGreen both yes
ED_I pink both yes
ED_D pink both yes
ED_C pink both yes
R454 green both yes
PSHP orange1 both yes
DGPc orange1 both yes
STMS lightblue both yes
STMU lightblue both yes
HAF0 grey both yes
HAF1 white both yes
HAF2 palegreen1 both yes
HAF3 green both yes
HAF4 yellow2 both yes
HAF5 red both yes
HAF6 indianred4 both yes
HAF7 black both yes
Q454 yellow both yes
H454 red both yes

Since MIRA tags is a bit of a moving target, be sure to read the manual (e.g. http://mira-assembler.sourceforge.net/docs/mira.html#section_40) to be sure to include all possible tags. Colors used in consed are X11 colors, a list of which is available at http://en.wikipedia.org/wiki/X11_color_names

=head1 CAVEAT

This script produces huge phdball files, because it writes one line per base, plus headers. For example, for 350'000 454 reads and some long Sanger ones, you get a file which is 1.4 Gb...

=head1 AUTHOR

Lionel Guy (lionel.guy@ebc.uu.se)

=head1 DATE

Mon Feb  1 11:49:46 CET 2010

=cut

# libraries
$|++; # to turn autoflush on
use strict;
use Getopt::Long;
use File::Basename;
use File::Copy;
use Date::Manip;

#############
# options and arguments
#############
my $caf_file;
my $ace_file;
my $sff_files;
my $solexa_fof;
my $solexa_paired;
my $phdballs;
my $consed_bin = "consed";
my $default_date = "Sat Jan 1 11:11:11 2000";
my $debug;
GetOptions(
    'c|caf_file=s' =>   \$caf_file,
    'a|ace_file=s' =>   \$ace_file,
    'sff_files=s' =>    \$sff_files,
    'solexa_fof=s' =>   \$solexa_fof,
    'solexa_paired' =>  \$solexa_paired,
    'p|phdballs=s' =>   \$phdballs,
    'x|consed_bin=s' => \$consed_bin,
    'default_date=s' => \$default_date,
    'd|debug' =>        \$debug,
);
usage() unless ($caf_file && $ace_file);

# globals
my %reads;

# checks opt_d
$default_date = localtime()
    if ($default_date eq "actual" || $default_date eq "current");

# checks that we are in a consed edit_dir (or assembly) folder)
die "Not in a consed folder or no phdball_dir folder available\n" 
    unless (-d "../phdball_dir");

# arrays of files
my (@sffs, @balls);
@sffs =  split(/,/, $sff_files) if ($sff_files);
@balls = split(/,/, $phdballs) if ($phdballs);
my %balls;

# sets output ace file
my $ace_out = basename($ace_file);
my $version = 1;
while (-e $ace_out){
    if ($ace_out =~ /(.*)\.(\d+)$/){
	$ace_out = $1 . "." . ($2+1);
    }
    else {
	$ace_out .= ".1";
    }
}

###############
# Convert sff and fastq to phdballs
###############
if ($sff_files || $solexa_fof){
    test_program_presence($consed_bin, 
			  "$consed_bin: consed could not be found.");
}
# sff file
foreach my $sff (@sffs){
    die "$sff file unreadable\n" unless (-r $sff);
    print "Reading SFF files, making phdballs\n";
    my $sff_name = basename($sff);
    my $ball_file = "../phdball_dir/$sff_name.ball";
    my $cmd = $consed_bin . " -sff2PhdBall " . $sff . " -phdBall $ball_file";
    print "Executing:\n$cmd\n";
    # check that the temp files don't exist
    if ($debug && -e "$ball_file.tmp"){
	move("$ball_file.tmp", $ball_file);
    }
    else {
	die "Temp file $ball_file.tmp exists" if (-e "$ball_file.tmp");
    	!system($cmd) || die "Couldn't run:\n$cmd\n";
    }
    # check that the ball actually exist
    die "Ball $ball_file doesn't exist\n" unless (-e $ball_file);
    $balls{$sff_name}{'ball'} = $ball_file;
    $balls{$sff_name}{'chem'} = '454';
}
# solexa data
if ($solexa_fof){
    die "$solexa_fof file unreadable\n" unless (-r $solexa_fof);
    print "Reading solexa files, making phdballs\n";
    my $cmd = $consed_bin . " -solexa2PhdBall " . $solexa_fof
	. " -phdBallFOF solexa_phdballs.fof";
    print "Executing:\n$cmd\n";
    unless ($debug){
	!system($cmd) || die "Couldn't run:\n$cmd\n";
    }
    # parse phdball fof
    open(FOF, "solexa_phdballs.fof") || 
	die "Couldn't open solexa_phdballs.fof\n";
    while (<FOF>){
	chomp;
	# check that the temp files don't exist
	if (-e "$_.tmp" && !$debug){
	    move("$_.tmp", $_);
	}
	else {
	    die "temp file $_.tmp exists";
	}
	my $ball_name = basename($_);
	# check that the ball actually exist
	die "Ball $_ doesn't exist\n" unless (-e $_);
	$balls{$ball_name}{'ball'} = $_;
	$balls{$ball_name}{'chem'} = 'solexa';
    }
}
# already existing phdballs
foreach (@balls){
    die "$_ file unreadable\n" unless (-r $_);
    my $ball_name = basename($_);
    $balls{$ball_name}{'ball'} = $_;
    $balls{$ball_name}{'chem'} = 'sanger';
}
###############
# Read phdballs a first time to retrieve what reads are in it,
# and make fake phdfiles for the other ones.
###############
# changing the endofblock operator
$/="BEGIN_SEQUENCE ";
print "\nListing reads in phdballs\n";
my %reads_in_balls;
my $n_reads_in_balls = 0;
foreach my $ball (sort keys %balls){
    my $n_reads_in_ball = 0;
    open(BALL_IN, $balls{$ball}{'ball'});
    print "  Reading ball $ball:\n";
    # read every "chunk", ie each read
    while (my $read = <BALL_IN>){
	chomp $read;
	# get id
	$read =~ /^(\S+)/;
	my $id = $1;
	# next if starts with begin_sequence (void...)
	unless ($id && $id =~ /^[^BEGIN_SEQUENCE]/){
	    next;
	}
	# check that the read is unique
	die "Read $id already present in another phdball\n" 
	    if $reads_in_balls{$id};
	# else, store the read name
	$reads_in_balls{$id} = 1;
	$n_reads_in_balls++;
	$n_reads_in_ball++;
	# keep user updated
	print "\r  " . ($n_reads_in_ball/1000) . "k reads found" 
	    unless ($n_reads_in_ball % 1000);
    }
    print "\r  $n_reads_in_ball reads found\n";
}
print "Found $n_reads_in_balls reads in " . scalar(keys %balls) . " ball(s).\n";
# changing back the endofblock operator
$/ = "\n";

###############
# Read CAF
###############
# open caf file
open(CAF, $caf_file);
# open an (eventual) phd ball for reads not present in the ones given
my $fakeball = "fake_ball";
open(PHD, ">../phdball_dir/$fakeball");

# initialize
my ($seq_id, $qual_id, $id);
my (@nts, @qs);
my $n_reads = 0;
my $n_gapped_reads = 0;
my $n_reads_not_in_balls = 0;
my $chem;
my %contigs;
print "\nReading CAF file\n";
my @leadnts;

# loop
while (<CAF>){
    chomp;
    # read sequence
    if (/^DNA : (.+)$/){
	# keeping the user updated
	$n_reads++;
	print "\r  " . ($n_reads/1000) . "k reads read" 
	    unless ($n_reads % 1000);
	# DEBUG
	#last if ($n_reads > 100000)
	# parsing id
	$seq_id = $1;
	# if read not in balls, parse seq and qual to create phd
	#unless ($reads_in_balls{$seq_id} || $contigs{$seq_id}){
	unless ($contigs{$seq_id}){
	    #print "Not in balls: $seq_id\n";
	    my $line = <CAF>;
	    chomp $line;
	    @nts = split(//, $line);
	    # as long as there is something different than a blank line, go on
	    while ($line =~ /^\S+$/ && ($line = <CAF>)){
		chomp $line;
		push @nts, split(//, $line);
	    }
	    # detect lower-case charcters in the beginning of the read, 
	    # meaning crap that need to be reported in the phd file
	    # reset leadnts
	    @leadnts = ();
	    #
	    foreach (@nts){
		if (/[a-z]/){
		    push @leadnts, $_;
		}
		else {
		    last;
		}
	    }
	}
    }
    # read base quality
    elsif (/^BaseQuality : (.+)$/){
    	$qual_id = $1;
 	#unless ($reads_in_balls{$qual_id} || $contigs{$seq_id}){ 
 	unless ($contigs{$seq_id}){
	    my $line = <CAF>;
	    chomp $line;
	    @qs = split(/ /, $line);
	    # as long as there is something different than a blank line, go on
	    while ($line =~ /^[0-9]/ && ($line = <CAF>)){
		chomp $line;
		push @qs, split(/ /, $line);
	    }
	}
    }
    # parse (eventually) other infos, then print read
    elsif (/^Sequence : (.+)$/){
	$id = $1;
	#DEBUG
	#if ($1 eq "FRKZZ3R02H7MIT"){
	#   print "hello $id\n";
	#}
	my $line = <CAF>;
	chomp $line;
	# if is contig, store it (do avoid parsing them in the DNA section)
	if ($line =~ /^Is_contig/){
	    $contigs{$id} = 1;
	    next;
	}
	# next unless seq type is read
	elsif (!$line =~ /^Is_read/){
	    next;
	}
	# if is read, check ids and seq lengths
	die "IDs non matching: $id, $seq_id, $qual_id\n" 
	    unless ($id eq $seq_id && $id eq $qual_id);
	die "Non equal lengths\n" unless (scalar @nts == scalar @qs);
	# initialize seqtype
	$chem = "Sanger";
	# initialize gaps
	my @gaps;
	my %gaps;
	my @local_leadnts = @leadnts;
	# loop through tags
	while ($line =~ /^\S+/ && ($line = <CAF>)){
	    chomp $line;
	    if ($line =~ /^Align_to_SCF (.+)$/){
		my @arr = split(/ /, $1);
		die "Not 4 elements at line $line\n" unless $#arr == 3;
		push @gaps, [ @arr ];
	    }
	    if ($line =~ /^Tag MINF.+ST=(.+)"$/){
		$chem = lc($1);
	    }
	}
	$reads{$id}{'chem'} = $chem;
	# if read not in balls, don't store gaps and leadnts
	if ($reads_in_balls{$id}){
	    # prepare gap information
	    if ($#gaps > 0){
		for my $i ( 1 .. $#gaps ){
		    $gaps{$gaps[$i-1][3]+1} = $gaps[$i][2] - $gaps[$i-1][3] - 1;
		}
		$n_gapped_reads++;
		$reads{$id}{'gaps'} = \%gaps;
	    }
	    # push leadnts info if not solexa
	    $reads{$id}{'leadnts'} = \@local_leadnts
		unless (lc($chem) eq "solexa");
	}
	# print PHD if not alrady in phdballs
	else {
	    $n_reads_not_in_balls++;
	    # print header
	    $reads{$id}{'ball'} = $fakeball;
	    $reads{$id}{'date'} = $default_date;
	    print_fake_phd($id, $chem, \@nts, \@qs, \%gaps);
	}
	# reset leadnts
	@leadnts = ();
    }
}
print "\r  $n_reads reads read from $caf_file.\n";
print "  $n_gapped_reads gapped reads\n";
if ($n_reads_not_in_balls){
    $balls{$fakeball}{'ball'} = "../phdball_dir/$fakeball";
    $balls{$fakeball}{'chem'} = 'sanger';
    $balls{$fakeball}{'date'} = $default_date;
    $balls{$fakeball}{'fake'} = 1;
    print "  Found $n_reads_not_in_balls reads that were not in phdballs\n";
    print "  A phdball (../phd_ball/$fakeball) has been created for these\n";
}
else {
    print "  All reads were present in phdballs.\n";
    unlink("../phdball_dir/$fakeball");
}
###############
# Read & modify phdballs
###############
# changing the endofblock operator
$/="BEGIN_SEQUENCE ";

print "\nReading phdballs\n";
my $n_tot_reads = 0;
foreach my $ball (sort keys %balls){
    # skip fake balls: already OK
    next if ($ball eq $fakeball);
    my $n_reads = 0;
    my $file = $balls{$ball}{'ball'};
    my ($date, $chem);
    die "File $file not found\n" unless (-r $file);
    move($file, "$file.tmp");
    open(BALL_IN, "$file.tmp");
    print "  Reading $file\n";
    open(BALL_OUT, ">$file");
    # read every "chunk", ie each read
    while (my $read = <BALL_IN>){
	chomp $read;
	# get id
	$read =~ /^(\S+)/;
	my $id = $1;
	# next if starts with begin_sequence (void...)
	unless ($id && $id =~ /^[^BEGIN_SEQUENCE]/){
	    print BALL_OUT "$read";
	    next;
	}
	# next if read here not in list (assembly)
	next unless $reads{$id};
	#die "Read $id not in read list" unless $reads{$id};
	$n_reads++;
	# keep user updated
	print "\r  " . ($n_reads/1000) . "k reads read" 
	    unless ($n_reads % 1000);
	# parse date in the beginning only
	if ($n_reads < 20){
	    $read =~ /\nTIME: (.+)\n/;
	    $date = $1;
	    $balls{$ball}{'date'} = $date;
	    $read =~ /\nCHEM: (.+)\n/;
	    $chem = $1;
	    die "Chem $chem and date $date must be defined in $read\n"
		unless ($chem && $date);
	}
	# but attribute it to all reads
	$reads{$id}{'chem'} = $chem;
	$reads{$id}{'date'} = $date;
	$reads{$id}{'ball'} = $ball;
	die "Chem $chem and date $date must be defined in $read\n"
	    unless ($chem && $date);
	# DEBUG
	#last if $n_reads > 100000;
	#print "$id\n" if ($id);
	# if this read is listed in the gapped ones; else just print
	my %gaps;
	%gaps = %{ $reads{$id}{'gaps'} } if ($reads{$id}{'gaps'});
	my @leadnts;
	@leadnts = @{ $reads{$id}{'leadnts'} } if ($reads{$id}{'leadnts'});
	# edit read if gaps or solexa (add that n in front of the read)
	if (%gaps || @leadnts || (lc($chem) eq "solexa")){
	    #print "$id\n";
	    print BALL_OUT "BEGIN_SEQUENCE ";
	    my @lines = split(/\n/, $read);
	    my $pos = 0;
	    my $skip = 0;
	    foreach my $line (@lines){
		# start counting at BEGIN_DNA
		if ($line =~ /BEGIN_DNA/){
		    $pos = 1;
		    print BALL_OUT "$line\n";
		    # add n at pos 1 if solexa: trick from mira to save mem
		    print BALL_OUT "n 0\n" if lc($chem) eq "solexa";
		    # add leadnts if ever needed
		    if (@leadnts){
			$pos = 1;
			my $peak_pos = 15;
			foreach (@leadnts){
			    print BALL_OUT "$_ 0 $peak_pos\n";
			    $peak_pos += 19;
			    $pos++;
			}
		    }
		}
		elsif ($pos > 0 && %gaps){
		    # gap at that position? then increase skip
		    if (%gaps && $gaps{$pos}){
			#print "pos: $pos, gap: $gaps{$pos}\n";
			$skip = $gaps{$pos};
		    }
		    # skip reads until skip = 0; else print
		    if ($skip){
			$skip--;
		    }
		    else {
			print BALL_OUT "$line\n";
		    }
		    # in any case, increase pos
		    $pos++;
		}
		else {
		    print BALL_OUT "$line\n";
		}
	    }
	}
	else {
	    print BALL_OUT "$/$read";
	}
    }
    print "\r  Read $n_reads reads\n";
    $n_tot_reads += $n_reads;
    #unlink("$file.tmp");
}
print "Read $n_tot_reads reads in " . 
    scalar(grep {"^fakeball"} keys %balls) . " phdballs\n";
# changing back the endofblock operator
$/ = "\n";

###############
# ACE file
###############
# modify all CHROMAT_FILE tags in ace file, if ace file is given
if ($ace_file){
    my $ct_tags = 0;
    my $rt_tags = 0;
    die "ACE out file $ace_out already exists. Choose another name or remove it"
	if (-e $ace_out);
    open(ACEO, ">$ace_out");
    open(TAGS, ">$ace_out.tags");
    open(ACE, $ace_file);
    print "\nReading ACE file\n";
    my $c = 0;
    my ($id, $chem);
    while (<ACE>){
	# keeping the user updated
	$c++;
	print "\r  " . ($c/1000) ."k lines read" unless ($c % 10000);
	# in the RD line
	if (/^RD\s(\S+)\s/){
	    $id = $1;
	    die "Read $id not in phdballs\n" unless ($reads{$id});
	    $chem = $reads{$id}{'chem'};
	    die "No chem defined for $id\n" unless ($chem);
	}
	# in the DS line
	if (/^DS\s/){
	    # hack to prevent date not being set
	    my $date =  $reads{$id}{'date'};
	    unless ($date){
		$date = $default_date;
		print "\nNo date set for $id, in $reads{$id}{'ball'}\n";
	    }
	    # we want different headers depending on chemistry
	    # for 454, we want CHROMAT_FILE with the sff, PHD_FILE, TIME
	    # CHEM and DYE (maybe)
	    if ($chem eq "454"){
		$_ = "DS CHROMAT_FILE: sff:$reads{$id}{'ball'}:$id ";
		$_ .= "PHD_FILE: $id.phd.1 TIME: $date ";
		$_ .= "CHEM: $chem DYE: unknown\n";
	    }
	    # for solexa, only the VERSION, TIME and CHEM
	    elsif ($chem eq "solexa"){
		$_ = "DS VERSION: 1 TIME: $date CHEM: $chem\n";
	    }
	    # else just parse the time correctly
	    else {
		$_ =~ s/TIME:.+\s[A-Z]+:/TIME: $date/;
		$_ =~ s/TIME:.+$/TIME: $date/;
	    }
	}
	# move tags around
	if (/^CT\{/) { # Tag found
	    $ct_tags++;
	    print TAGS $_;
	    while (<ACE>) {
		print TAGS $_;
		if ($_ =~ /^\}/) { # End of tag, continue
		    print TAGS "\n";
		    last;
		}
	    }
	}
	elsif (/^RT\{/) { # Tag found
	    $rt_tags++;
	    print TAGS $_;
	    while (<ACE>) {
		print TAGS $_;
		if ($_ =~ /^\}/) { # End of tag, continue
		    print TAGS "\n";
		    last;
		}
	    }
	}
	else {
	    print ACEO $_;
	}
    }
    # print links to phdballs
    foreach my $ball (sort keys %balls){
	my $file = $balls{$ball}{'ball'};
	my $date = $balls{$ball}{'date'};
	$date = parse_phd_date($date);
	#print "$ball: $date";
	print ACEO "WA{\nphdBall consed $date \n$file\n}\n\n";
    }
    # done
    print "\r  Done reading ACE file. Read $c lines\n";
    # append tags at the end of the file 
    print "  Transferring tags...";
    system("cat $ace_out.tags >> $ace_out");
    # print final information
    print "\r  Transferred $ct_tags CT tags and $rt_tags RT tags ";
    print "at the end of file.\n";
    print "\nRun consed: \n";
    print "$consed_bin -ace $ace_out &\n";
}

exit 1;

sub usage{
    system("perldoc $0");
    exit;
}

sub test_program_presence{
    my ($prg, $message) = @_;
    $message = "" unless $message;
    system("which $prg &> /dev/null") == 0 
	or die "Failed to execute $prg. $message\n";
}

# parses Wed Dec 24 11:21:50 2008 to 081224:112150
# using Date::Manip
sub parse_phd_date{
    my ($str) = @_;
    my $date = Date::Manip::Date->new;
    $date->parse($str);
    my $res = $date->printf("%y%m%d:%H%M%S");
    return $res;
}

# prints the header of the phd file
sub print_fake_phd{
    my ($name, $chem, $nt_ref, $q_ref, $gaps_ref) = @_;
    # print header
    print PHD "BEGIN_SEQUENCE $name\n\n";
    print PHD "BEGIN_COMMENT\n\n";
    print PHD "CHROMAT_FILE: none\n";
    print PHD "CALL_METHOD: " . basename($0) . "\n";
    print PHD "QUALITY_LEVELS: 99\n";
    print PHD "CHEM: $chem\n";
    print PHD "TIME: $default_date\n";
    print PHD "END_COMMENT\n\n";
    print PHD "BEGIN_DNA\n";
    # loop through nts
    my $gap_counter = 0;
    my $peak_pos = 15;
    foreach my $i (0 .. scalar(@{$nt_ref})-1){
	# do not print gaps
	next if ($nt_ref->[$i] =~ /-/);
	# print to PHD
	print PHD $nt_ref->[$i] . " " . $q_ref->[$i];
	# print peak position only if chem is not solexa
	print PHD " " . $peak_pos unless ($chem eq "solexa");
	print PHD "\n";
	$peak_pos += 19;
	# go over gaps if any
	$peak_pos += 19 * $gaps_ref->{$i+1} if ($gaps_ref->{$i+1});
    }
    print PHD "END_DNA\n\n";
    print PHD "END_SEQUENCE\n\n";
    # if solexa and paired-end, print wr tags
    if ($solexa_paired && $chem){
	$id =~ /([^\/]+)\/([12])/;
	die "$id not a valid paired-end read id. Should end with /1 or /2.\n" 
	    unless ($1 && $2);
	my $template = $1;
	my $dir = 'fwd';
	$dir = 'rev' if ($2 == 2);
	print PHD "WR{\ntemplate " . basename($0) . " $default_date\n";
	print PHD "name: $template\n\}\n\n";
	print PHD "WR{\nprimer " . basename($0) . " $default_date\n";
	print PHD "type: univ $dir\n}\n";
    }
    # loop through nts
    # my $gap_counter = 0;
    # my $peak_pos = 15;
    # foreach my $i (0 .. $#nts){
    # 	# do not print gaps
    # 	next if ($nts[$i] =~ /-/);
    # 	# calculate peak position
    # 	my $q = $qs[$i];
    # 	print PHD "$nts[$i] $q $peak_pos\n";
    # 	$peak_pos += 19;
    # 	# go over gaps if any
    # 	$peak_pos += 19 * $gaps{$i+1} if ($gaps{$i+1});
    # }
    # print PHD "END_DNA\n\n";
    # print PHD "END_SEQUENCE\n\n";
}

# ###############
# # SFF file & chromat_dir
# ###############
# # if opt_f stated, check availability of sffinfo and mktrace
# # then run sffinfo to get the names of the reads into the sff file
# my %sff_reads;
# my %old_chromats;
# if ($opt_f){
#     die "Must use -s in conjunction with -f\n" unless $sff_files;
#     test_program_presence("sffinfo", 
# 			  "It is included in the GS assembler package");
#     test_program_presence("mktrace", "It is included in consed package");
#     # sff
#     my $call = "sffinfo -accno $sff_files |";
#     open(SFFINFO, $call) or die "Invalid call to sffinfo: $call\n";
#     while (<SFFINFO>){
# 	die "Invalid output of sffinfo: $_\n" unless /^\S+$/;
# 	$sff_reads{$_}++;
#     }
#     # chromat dir
#     opendir(CHROMATDIR, $opt_f);
#     while (<CHROMATDIR>){
# 	next if (!-e $_ || -d $_ || -z $_);
# 	$old_chromats{$_}++;
#     }
# }
# fix a list of reads in sff and in chromat_dir
