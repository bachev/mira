#!/usr/bin/perl
#########################################################################################
# lucy2xml.pl: A perl script to convert output from the Lucy quality & vector
# trimming program into an XML file according to the NCBI Traceinfo RFC found at 
# http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=rfc&m=main&s=rfc
#				 ____________________
#				|	           |
#				|     Sept 2008      |
#				|  Ross W. Whetten   |
#				|    version 0.3       |
#				|____________________|
#
#########################################################################################
#
# Usage: lucy2xml.pl -i (path/to/lucyoutput) -o (path/to/tracefile.xml) -p (y/n, will
#        paired-end reads be processed) [if yes, provide additional switches] -s (position
#        of 1st character of template_id in tracename string, starting with > as 0) -l (length
#        of template_id string) -t (position of a character specifying end identity for trace_end
#        field) -f (the character denoting FORWARD) -r (the character denoting REVERSE)
#########################################################################################
#
# Output from script:
# XML-tagged file with individual trace data required by MIRA for assembly.  Other information 
#   required by NCBI for submission must be added to create a Traceinfo XML file that meets the
#   requirements of the Trace Archive.
# TRACE_NAME, CLIP_VECTOR_LEFT, and CLIP_VECTOR_RIGHT fields are extracted from the lucy 
#   output file.
# Lucy does not provide separate CLIP_QUALITY_LEFT and CLIP_QUALITY_RIGHT values; it uses
#   both base quality values and presence of vector sequences to derive a single pair of 
#   clip values. The NCBI RFC specifies, if only a single set of clip values are provided, 
#   that those values occupy the CLIP_VECTOR_LEFT and CLIP_VECTOR_RIGHT fields of the 
#   traceinfo file.
#
#########################################################################################
use strict;
use Getopt::Std; 

# collect the options to variables; warn if not defined
my($lucyinfile, $xmloutfile, $paired, $startid, $lengthid, $traceloc, $forw, $rev, %opt);
getopts('i:o:p:c:l:t:f:r:h?',\%opt);
chomp(%opt);
if( $opt{h} || $opt{'?'} ) {
	print "\n\nUsage: lucy2xml.pl -i (path/to/lucyoutput) -o (path/to/tracefile.xml)\n";           
	print "-p y/n - (will paired-end reads be processed?)\n\t [if yes, provide additional switches]\n";       
	print "-c value (position of 1st character of template_id in\n\t tracename string, starting with > as 0)\n";   
	print "-l value (length of template_id string)\n";
	print "-t value (position of a character specifying clone end\n\t identity for trace_end field)\n";          
	print "-f char (the character in the tracename string that denotes FORWARD)\n";  
	die "-r char (the character in the tracename string that denotes REVERSE)\n\n";	
}
if($opt{i}) {
	$lucyinfile = $opt{i}; 
}
if($opt{o}) {
	$xmloutfile = $opt{o}; 
}
if($opt{p}) {
	$paired = $opt{p}; 
}
if($opt{c}) {
	$startid = $opt{c}; 
}
if($opt{l}) {
	$lengthid = $opt{l}; 
}
if($opt{t}) {
	$traceloc = $opt{t}; 
}
if($opt{f}) {
	$forw = $opt{f}; 
}
if($opt{r}) {
	$rev = $opt{r}; 
}
if($paired eq "n" || $paired eq "no") {
	$main::p = "n";
	unless (defined $lucyinfile && defined $xmloutfile) {
		print "usage: lucy2xml.pl -i path/to/lucyoutput -o path/to/traceinfo.xml \n" ;
		die "\tfor non-paired-end data. Enter lucy2xml.pl -h for full usage information.\n";
	}
}
elsif($paired eq "y" || $paired eq "yes") {
	$main::p = "y";
	unless (defined $lucyinfile && defined $xmloutfile && defined $startid && defined $lengthid && defined $traceloc && defined $forw && defined $rev) {
		print "\n\nUsage: lucy2xml.pl -i (path/to/lucyoutput) -o (path/to/tracefile.xml)\n";           
		print "-p y/n - (will paired-end reads be processed?)\n\t [if yes, provide additional switches]\n";       
		print "-c value (position of 1st character of template_id in\n\t tracename string, starting with > as 0)\n";   
		print "-l value (length of template_id string)\n";
		print "-t value (position of a character specifying clone end\n\t identity for trace_end field)\n";          
		print "-f char (the character in the tracename string that denotes FORWARD)\n";  
		die "-r char (the character in the tracename string that denotes REVERSE)\n\n";	
	}
}
else {
	print "\n\nUsage: lucy2xml.pl -i (path/to/lucyoutput) -o (path/to/tracefile.xml)\n";           
	print "-p y/n - (will paired-end reads be processed?)\n\t [if yes, provide additional switches]\n";       
	print "-c value (position of 1st character of template_id in\n\t tracename string, starting with > as 0)\n";   
	print "-l value (length of template_id string)\n";
	print "-t value (position of a character specifying clone end\n\t identity for trace_end field)\n";          
	print "-f char (the character in the tracename string that denotes FORWARD)\n";  
	die "-r char (the character in the tracename string that denotes REVERSE)\n\n";	
}

	
#########################################################################################
#
# MODIFY THESE VARIABLES TO INCLUDE THE APPROPRIATE VALUES IN THE TRACEINFO OUTPUT FIELDS
#
#   If these variables are to be defined, uncomment these lines and the print lines later in the  
#   script to print the appropriate values to the traceinfo.xml output file
#
# my $centername = "CENTER_NAME_ACRONYM"; #must be obtained through NCBI
# my $submission = "SUBMISSION_TYPE"; #new, update, updateinfo, withdraw
# my $centerproject = "PROJECT_NAME"; #for internal use within center
# my $strategy = "STRATEGY_TYPE"; #wgs, chip, cloneend, est - see RFC for complete list
# my $tracetype = "TRACE_TYPE_CODE"; #often the same as STRATEGY
# my $sourcetype = "SOURCE_TYPE"; #G (=genomic), N (=non-genomic), VIRAL RNA, SYNTHETIC
# my $genusspecies = "GENUS_SPECIES";
#
#########################################################################################
#
#  THESE VALUES SHOULD BE CHANGED; SET TO THESE DEFAULTS FOR TESTING ONLY
#
#  Change to appropriate values and uncomment the lines if necessary
# my $tracepath = 'local/path/to/traces/';
# my $traceformat = "scf"; # HERE SET TO SCF - change to correct value (scf, sff, ztr, or abi)
# my $insertsize = "10000"; # HERE SET TO 10 KB - change to correct value
# my $insertstddev = "500"; # HERE SET TO 500 BP - change to correct value
#########################################################################################
#
# Open the output file and print desired fields at the top of the file

open(my $out, ">", "$xmloutfile") or die "The named output file could not be opened\n";
print $out "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
print $out "<trace_volume>\n";

#########################################################################################
# UNCOMMENT THE FOLLOWING LINES IF THESE FIELDS ARE TO BE PRINTED TO THE XML FILE

#print $out "   <common_fields>\n";
#print $out "      <center_name>$centername</center_name>\n";
#print $out "      <submission_type>$submission</submission_type>\n";
#print $out "      <strategy>$strategy</strategy>\n";
#print $out "      <trace_type_code>$tracetype</trace_type_code>\n"; #often same as STRATEGY
#print $out "      <center_project>$centerproject</center_project>\n";
#print $out "      <source_type>$sourcetype</source_type>\n";
#print $out "      <species_code>$genusspecies</species_code>\n";
#print $out "      <insert_size>$insertsize</insert_size>\n";
#print $out "      <insert_stddev>$insertstddev</insert_stddev>\n";
#print $out "   </common_fields>\n";

#########################################################################################


# open the input file and read in one line at a time (to keep memory use low)
# discard any line that is not a header; split header lines on spaces into @data array
# use the leading > from the tracename to complete the <trace_name> tag


open(my $in,  "<",  "$lucyinfile")  or die "Can't open $lucyinfile";
if ($main::p eq "n") {
	##################################################################################
	#
	# FOR NON-PAIRED-END DATA
	while (<$in>){
		if (/^>/) {
			my @data = split(/\s+/);
			print $out "   <trace>\n";
			print $out "      <trace_name$data[0]</trace_name>\n";
			print $out "      <clip_vector_left>$data[4]</clip_vector_left>\n";
			print $out "      <clip_vector_right>$data[5]</clip_vector_right>\n";
			print $out "   </trace>\n";
			next;
		}
		else {next};
	}
	print $out "</trace_volume>\n";
}
elsif ($main::p eq "y") {

	##################################################################################
	#
	# FOR PAIRED-END DATA 
	while (<$in>) {	
		if (/^>/) {
			my @data = split(/\s+/);
	 		my $templateid = substr($data[0], $startid, $lengthid); #start at designated char of tracename, take designated length
	 		my $testchar = substr($data[0],$traceloc,1); # take the specified char of tracename, test against provided chars
			print $out "   <trace>\n";
			print $out "      <trace_name$data[0]</trace_name>\n";
			print $out "      <template_id>$templateid</template_id>\n";
			if($testchar eq $forw) {
				print $out "      <trace_end>FORWARD</trace_end>\n";
			}
			elsif($testchar eq $rev) {
				print $out "      <trace_end>REVERSE</trace_end>\n";
			}
			else {
				print "$testchar fails to match $forw or $rev for trace $data[0]\n";
			}				
			print $out "      <clip_vector_left>$data[4]</clip_vector_left>\n";
			print $out "      <clip_vector_right>$data[5]</clip_vector_right>\n";
			print $out "   </trace>\n";
			next;
		}
		else {next};
	}
	print $out "</trace_volume>\n";
}
close $in;
close $out;
exit;




