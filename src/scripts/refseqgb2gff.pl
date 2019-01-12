#!/usr/bin/env perl -w

eval 'exec /usr/bin/env perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell



use strict;
use warnings;
use diagnostics;

use lib "$ENV{HOME}/bioperl-live";
# chad put this here to enable situations when this script is tested
# against bioperl compiled into blib along with other programs using blib
BEGIN {
        unshift(@INC,'blib/lib');
};
use Pod::Usage;
use Bio::Root::RootI;
use Bio::SeqIO;
use File::Spec;
use Bio::SeqFeature::Tools::Unflattener;
use Bio::SeqFeature::Tools::TypeMapper;
use Bio::SeqFeature::Tools::IDHandler;
use Bio::Location::SplitLocationI;
use Bio::Location::Simple;
use Bio::Tools::GFF;
use Getopt::Long;


use vars qw/$split @filter $zip $outdir $help $ethresh
            $file @files $dir $summary $nolump 
            $source_type %proteinfa %exonpar $didheader $verbose $DEBUG $GFF_VERSION 
            $gene_id $rna_id $tnum $ncrna_id $rnum %method %id/;


use vars qw(%glob_seenbaseids
	    %glob_namepriorities
	    %glob_unneededtags
	    %glob_allowtagsonlythere
	    %glob_unsearchabletags
	    %glob_gb2gfftagmap
	    %glob_ptag2idpostfix
	    %glob_ptagotherpostfix
	    $glob_ptopcounter
	    @glob_fasta_sname
	    @glob_fasta_sdata
	    $glob_gffio 
	    $glob_outfh
	    $glob_tmpdir4tmpgff
	 );


$glob_namepriorities{'stdtag'} = [qw/gene standard_name gene_synonym label locus_tag protein_id old_locus_tag/];

# disallow certain tags in all primary tags but the cited ones
# e.g., 'Parent' may occur only in 'gene', 'CDS', etc.pp and will be deleted from all other ones
# e.g., 'gene_synonym' may occur in none (is merged into 'Alias' by the script
$glob_allowtagsonlythere{'Parent'} = [qw/gene CDS mRNA ncRNA rRNA tRNA tmRNA/];
$glob_allowtagsonlythere{'gene_synonym'} = [qw//];

# this removes unneeded (or even not allowed) tags from elements before they get written as GFF
# performs the same task of cleanup like $allowtagsonlythere, but allows a more fine grained control
$glob_unneededtags{'region'} = [qw/gene standard_name gene_synonym label locus_tag protein_id old_locus_tag/];
#$glob_unneededtags{'gene'} = [qw//];
#$glob_unneededtags{'mRNA'} = [qw//];
#$glob_unneededtags{'exon'} = [qw//];
#$glob_unneededtags{'CDS'} = [qw//];
#$glob_unneededtags{'mature_protein_region'} = [qw//];

# instead of removing 'Name' and 'Alias' in certain elements, rewrite them to 'name' and 'alias' so that
#  the info is still viewable by the user in GBrowse2, but not searched anymore in the
#  GBrowse2 search function
$glob_unsearchabletags{'gene'} = 1;
$glob_unsearchabletags{'mRNA'} = 1;
$glob_unsearchabletags{'exon'} = 1;
$glob_unsearchabletags{'mature_protein_region'} = 1;


# Tranform some GenBank tags to the canonical GFF3 versions
my %glob_gb2gfftagmap = (
  db_xref => 'Dbxref',
  note => 'Note'
);


# choose which postfix for which genbank(gff) feature
# "gene" and "contig" don't get postfixes!
# if not in this list, will chose a generix ".o<type>"
my %glob_ptag2idpostfix = (
  "gene" => "",
  "contig" => "",
  "CDS" => ".p",
  "protein" => ".p",
  "polypeptide" => ".p",
  "mature_protein_region" => ".mp",
  "region" => ".f",
  "mRNA" => ".t",
  "RNA" => ".r",
  "transcript" => ".r",
  "exon" => ".e",
  "signal_peptide" => ".sp"
);

my %glob_ptagotherpostfix;
my $glob_ptopcounter=1;

my $opt_informat="GenBank"; # swiss ; embl; genbank ; ** guess from SOURCEID **
my $opt_regionssearchable=1;


$GFF_VERSION = 3; # allow v2 ...
$verbose = 1; # right default? -nov to turn off

my $FORMAT="GenBank"; # swiss ; embl; genbank ; ** guess from SOURCEID **
my $SOURCEID= $FORMAT; # "UniProt" "GenBank"  "EMBL" should work  
   # other Bio::SeqIO formats may work.  TEST: EntrezGene < problematic tags; InterPro  KEGG 



$| = 1;
my $quiet= !$verbose;
my $ok= GetOptions( 'd|dir|input:s'   => \$dir,
            'z|zip'     => \$zip, 
            'h|help'    => \$help,
            's|summary' => \$summary,
            'o|outdir|output:s'=> \$outdir,
            'x|filter:s'=> \@filter,
            'y|split'   => \$split,
            "ethresh|e=s"=>\$ethresh,
            'f|format=s' => \$FORMAT,
            'typesource=s' => \$source_type,
            'GFF_VERSION=s' => \$GFF_VERSION,
            'quiet!'    => \$quiet, # swap quiet to verbose
            'DEBUG!'    => \$DEBUG,
            'n|nolump'  => \$nolump);

my $lump = 1 unless $nolump || $split;
$verbose= !$quiet;

# look for help request
pod2usage(2) if $help || !$ok;

# keep SOURCEID as-is and change FORMAT for SeqIO types; 
# note SeqIO uses file.suffix to guess type; not useful here
$SOURCEID= $FORMAT; 
$FORMAT  = "swiss" if $FORMAT =~/UniProt|trembl/;
$verbose =1 if($DEBUG);

# initialize handlers
my $unflattener = Bio::SeqFeature::Tools::Unflattener->new; # for ensembl genomes (-trust_grouptag=>1);
$unflattener->error_threshold($ethresh) if $ethresh;
$unflattener->verbose(1) if($DEBUG);
# $unflattener->group_tag('gene') if($FORMAT =~ /embl/i) ; #? ensembl only? 
# ensembl parsing is still problematic, forget this

my $tm  = Bio::SeqFeature::Tools::TypeMapper->new;
my $idh = Bio::SeqFeature::Tools::IDHandler->new;

# dgg
$source_type ||= "region"; # should really parse from FT.source contents below

my $FTSOmap = $tm->FT_SO_map();
# #convert $FTSOmap undefined to valid SO : moved to TypeMapper->map_types( -undefined => "region")

# stringify filter list if applicable
my $filter = join ' ', @filter  if @filter;

# determine input files
my $stdin=0; # dgg: let dir == stdin == '-' for pipe use
if ($dir && ($dir eq '-' || $dir eq 'stdin')) {
  $stdin=1;  $dir=''; @files=('stdin');
  
} elsif ( $dir ) {
    if ( -d $dir ) {
        opendir DIR, $dir or die "could not open $dir for reading: $!";
        @files = map { "$dir/$_";} grep { /\.gb.*/ } readdir DIR;  
        closedir DIR;
    }
    else {
        die "$dir is not a directory\n";
    }
}
else {
    @files = @ARGV;
    $dir = '';
}

# we should have some files by now
pod2usage(2) unless @files;


my $stdout=0; # dgg: let outdir == stdout == '-' for pipe use
if($outdir && ($outdir eq '-' || $outdir eq 'stdout')) {
  warn("std. output chosen: cannot split\n") if($split);
  warn("std. output chosen: cannot zip\n") if($zip);
  warn("std. output chosen: cannot nolump\n") if($nolump);
  $stdout=1; $lump=1; $split= 0; $zip= 0; # unless we pipe stdout thru gzip
  
} elsif ( $outdir && !-e $outdir ) {
    mkdir($outdir) or die "could not create directory $outdir: $!\n";        
}
elsif ( !$outdir ) {
    $outdir = $dir || '.';
}

$outdir .= '/' unless $outdir =~ m|/$|;

my $glob_tmpdir4tmpgff = File::Temp->newdir(DIR => '/tmp/', CLEANUP => 1);

my ($volume,$directories,$fn) = File::Spec->splitpath($files[0]);
my $outfile = $outdir . $fn . ".gff";
open $glob_outfh, ">$outfile";
writeHeader();


for my $file ( @files ) {
  process_file($file);
}

exit;

sub writeHeader {
  print $glob_outfh "##gff-version 3\n#conversion-by mygb2gff.pl\n"
}


sub process_file {
  my $filename = shift;

  # dgg ; allow 'stdin' / '-' input ?
  chomp $filename;
  die "$! $filename" unless($stdin || -e $filename);
  print "# Input: $filename\n" if($verbose);

  # open input file, unzip if req'd
  if ($stdin) {
    *FH= *STDIN;   
  } elsif ( $filename =~ /\.gz/ ) {
    open FH, "gunzip -c $filename |";
  }
  else {
    open FH, "<$filename";
  }

  my $in = Bio::SeqIO->new(-fh => \*FH, -format => $FORMAT, -debug=>$DEBUG);
  $glob_gffio = Bio::Tools::GFF->new( -noparse => 1, -gff_version => $GFF_VERSION );

  munchFile($in);
  appendTmpGFFs();

  writeFasta();

  close FH;

}

sub appendTmpGFFs {
  my @direntries;
  opendir(HU, $glob_tmpdir4tmpgff) || return;
  @direntries=sort(readdir(HU));
  closedir(HU);

  foreach my $tmpf (@direntries) {
    if($tmpf =~ /gff3/){
      open FH, "$glob_tmpdir4tmpgff/$tmpf" or next;
      while(my $line=<FH>){
	print $glob_outfh $line
      }
      close FH;
    }
  }

  system("rm -f $glob_tmpdir4tmpgff/*");
}


sub writeFasta {
  print "want to write FASTA ...\n";
  if(scalar(@glob_fasta_sname)){
    print $glob_outfh "##FASTA\n";
    for my $sname (@glob_fasta_sname){
      my $seq = shift(@glob_fasta_sdata);
      $seq  =~ s/(\S{1200})/$1\n/g;
      print $glob_outfh ">$sname\n";
      print $glob_outfh "$seq\n";
    }
    @glob_fasta_sdata= ();
    @glob_fasta_sname= ();
  }
}

sub munchFile {
  my $in = shift;

  while ( my $seq = $in->next_seq ) {
    my $seq_name = $seq->accession_number;
    my $end = $seq->length;

    #my $outfile = $outdir . $seq_name . ".gff";
    #open $glob_outfh, ">$outfile";

    print "sn: $seq_name\n";

    my $dna=$seq->seq;
    if($dna && $seq_name){
      print "Pushing ... $seq_name (",length($dna),")\n";
      push(@glob_fasta_sdata,$dna);
      push(@glob_fasta_sname,$seq_name);
    }

    # unflatten gene graphs, apply SO types, etc; this also does TypeMapper ..
    unflatten_seq($seq);


    my %elemcounter;
    for my $f ($seq->get_SeqFeatures) {
      if(exists($elemcounter{$f->primary_tag})){
	++$elemcounter{$f->primary_tag};
      }else{
	$elemcounter{$f->primary_tag}=1;
      }

      traverse($f,0,$elemcounter{$f->primary_tag},"","");
    }
  }
}


sub traverse {
  my ($sf,$level,$elemnum,$baseid,$parentid) = @_;

  $sf->source_tag("GenBank");
   #  $sf isa Bio::SeqfeatureI

  # get a base id for this whole thing, take care no doubles are made
  #if(length($baseid)==0){
    $baseid=generateBaseID($sf);
    if(exists $glob_seenbaseids{$baseid}){
      my $newid=$baseid;
      while(exists $glob_seenbaseids{$newid}){
	++$glob_seenbaseids{$baseid};
	$newid.="_$glob_seenbaseids{$baseid}";
      }
      $baseid=$newid;
    }else{
      $glob_seenbaseids{$baseid}=1;
    }
  #}


  # make a feature id
  my $fid="";
  if(length($parentid)){
    $fid=generateFeatureID($sf,$elemnum,$parentid);
  }else{
    $fid=generateFeatureID($sf,$elemnum,$baseid);
  }

  if(0){
    print $level; for(my $i=0; $i<$level;++$i) {print "  ";}
    print "PT: ", $sf->primary_tag, "\n";
    my @tags= $sf->get_all_tags();
    for my $t (@tags) {
      print $level; for(my $i=0; $i<$level;++$i) {print "  ";}
      my ($tv)=$sf->get_tag_values($t);
      print "TN: $t\t$tv\n";
    }
    
    print "BID: $baseid\n";
    print "FID: $fid\n";
    print "PID: $parentid\n";
  }

  maptags2gff($sf);

  addStandardTags($sf, $baseid, $fid, $parentid);
  reworkAliasTagValues($sf);

  leaveOnlyAllowedTags($sf);
  removeUnneededTags($sf);

  makeTagsUnsearchable($sf);
  renameTagValues($sf);

  reworkNote($sf);
  reworkFunction($sf);
  addNoteIfNotExisting($sf);


  if($sf->has_tag('translation')){
    moveTranslationToFasta($sf);
  }

  my $gffstr= $glob_gffio->gff_string($sf);
  print $glob_outfh "$gffstr\n" if $gffstr;

  my @subs = ();
  if ($sf->isa("Bio::FeatureHolderI")) {
    @subs = $sf->get_SeqFeatures;
  } elsif ($sf->isa("Bio::SeqFeatureI")) {
    @subs = $sf->sub_SeqFeature;
  }
  
  my %elemcounter;
  for my $subsf (@subs){
    if(exists($elemcounter{$subsf->primary_tag})){
      ++$elemcounter{$subsf->primary_tag};
    }else{
      $elemcounter{$subsf->primary_tag}=1;
    }

    traverse($subsf,$level+1,$elemcounter{$subsf->primary_tag},$baseid,$fid);
  }
  #traverse($_,$level+1,"") foreach @subs;

}


sub leaveOnlyAllowedTags {
  my $sf = shift;
  my $ptag=$sf->primary_tag;

  foreach my $tag (keys %glob_allowtagsonlythere) {
    if ($sf->has_tag($tag)) {
      my $keeptag=0;
      foreach my $allowed (@{$glob_allowtagsonlythere{$tag}}) {
	if($allowed eq $ptag){
	  $keeptag=1;
	  last;
	}
      }
      if(!$keeptag){
	$sf->remove_tag($tag);
      }
    }
  }

  if(exists($glob_unneededtags{$ptag})){
    foreach my $tag (@{$glob_unneededtags{$ptag}}) {
      if ($sf->has_tag($tag)) {
	$sf->remove_tag($tag);
      }
    }
  }
}

sub removeUnneededTags {
  my $sf = shift;
  my $ptag=$sf->primary_tag;

  if(exists($glob_unneededtags{$ptag})){
    foreach my $tag (@{$glob_unneededtags{$ptag}}) {
      if ($sf->has_tag($tag)) {
	$sf->remove_tag($tag);
      }
    }
  }
}


# In gene, mRNA, etc., rewrite the tags "Name" to "name" and "Alias" to "alias"
sub makeTagsUnsearchable {
  my $sf = shift;
  my $ptag=$sf->primary_tag;

  if(exists($glob_unsearchabletags{$ptag})){
    foreach my $tag (qw[Name Alias]) {
      if ($sf->has_tag($tag)) {
	my @val= $sf->get_tag_values($tag);
	$sf->remove_tag($tag);
	my $lctag=lc($tag);
	foreach my $v (@val) {
	  $sf->add_tag_value($lctag,$v);
	}
      }
    }
  }
}

# genbank may have "bla1; bla2; bla3" etc. in the names
# make it a regular multiple entry
sub reworkAliasTagValues {
  my $sf = shift;

  if ($sf->has_tag('Alias')) {
    my $val= join("; ",$sf->get_tag_values('Alias'));
    $sf->remove_tag('Alias');
    foreach my $a (split(/; /,$val)) {
      $sf->add_tag_value('Alias',$a);
    }
  }
}

# Some things are not written in notes and no note is existing
# in those case, make a simple note by taking the most probable entry
# so that GBrowse2 has something to display
sub addNoteIfNotExisting {
  my $sf = shift;

  if (!$sf->has_tag('Note')) {
    my $val="";
    if ($sf->has_tag('product')) {
      $val=join("; ",$sf->get_tag_values('product'));
    } elsif ($sf->has_tag('mobile_element')) {
      $val=join("; ",$sf->get_tag_values('mobile_element'));
    } elsif ($sf->has_tag('function')) {
      ($val)=$sf->get_tag_values('function');
    }
    if(length($val)){
      $sf->add_tag_value('Note',$val);
    }
  }
}

# for regions etc., rename "Name=xxx" to "Name=xxx_region"
sub renameTagValues {
  my $sf = shift;
  my $ptag=$sf->primary_tag;

  if($ptag eq "region"){
    if ($sf->has_tag('Name')) {
      my ($val)= $sf->get_tag_values('Name');
      $val.="_region";
      $sf->remove_tag('Name');
      $sf->add_tag_value('Name',$val);
    }
  }
}


# for regions, remove "; other site" at end of Note for misc_feature
# for regions, rewrite "Note=active sites" (and variations with upper/lower cases) to "Note=active site"
sub reworkNote {
  my $sf = shift;
  my $ptag=$sf->primary_tag;

  if($ptag eq "region"){
    if ($sf->has_tag('Note')) {
      my $val= join("; ",$sf->get_tag_values('Note'));
      #print "didu: $val\n";
      $val =~ s/; other site$//;
      #print "didi: $val\n";
      if (lc($val) eq "active site"){
	$val="active site";
      }
      if (lc($val) eq "active sites"){
	$val="active site";
      }
      $sf->remove_tag('Note');
      $sf->add_tag_value('Note',$val);
    }
  }
}

# 1) remove "; not classified" and end of "function" tags
sub reworkFunction {
  my $sf = shift;
  my $ptag=$sf->primary_tag;

  if ($sf->has_tag('function')) {
    my $val= join("; ",$sf->get_tag_values('function'));
    $val =~ s/; not classified$//i;
    $sf->remove_tag('function');
    $sf->add_tag_value('function',$val);
  }
}

sub addStandardTags {
  my ($sf, $baseid, $fid, $parentid) = @_;

  my $ptag=$sf->primary_tag;

  $sf->add_tag_value( ID => $fid);
  if(length($parentid)) {$sf->add_tag_value( Parent => $parentid);}

  my $hasname=$sf->has_tag('Name');
  my $hasalias=$sf->has_tag('Alias');

  if(! exists($glob_namepriorities{$ptag})){
    $ptag="stdtag";
  }
    
  for my $ntag (@{$glob_namepriorities{$ptag}}) {
    if($sf->has_tag($ntag)){
      my @names=$sf->get_tag_values($ntag);
      #print "Tag: ", $ntag, "\tvalues: ", @names,"\n";
      if(!$hasname){
	#print "Create name: ", $names[0],"\n";
	$hasname=1;
	$sf->add_tag_value( Name => $names[0]);
	shift @names;
      }
      if(scalar(@names)){
	if(!$hasalias){
	  #print "Create alias: ", @names,"\n";
	  $hasalias=1;
	}
	$sf->add_tag_value( Alias => @names);
      }
    }
  }

  if(!$hasname){
    $sf->add_tag_value( Name => $fid);
  }
}

sub generateBaseID {
  # modelled after get_persistent_id() of BioPerl
  # but using locus_tag if possible before the generic ID
  
  my $sf = shift;

  my $id='';
  if (!$sf->isa("Bio::SeqFeatureI")) {
    $sf->throw("not a Bio::SeqFeatureI");
  }
  my $seq_id = $sf->seq_id || $sf->throw("seq_id must be set: ".$sf->display_name);
  if ($sf->primary_tag eq "contig"){
    # let's pray that the contig names are unique
    $id=$seq_id;
  } elsif ($sf->has_tag('transcript_id')) {
    ($id) = $sf->get_tag_values('transcript_id');
  } elsif ($sf->has_tag('protein_id')) {
    ($id) = $sf->get_tag_values('protein_id');
  } elsif ($sf->has_tag('locus_tag')) {
    ($id) = $sf->get_tag_values('locus_tag');
  }

  # if no id was made yet or if primary tag is not one of contig|gene|mRNA|CDS|exon
  #  then generate a generic ID
  #print "dubidi: $id \t", $sf->primary_tag, "\n";
  if(length($id)==0 || $sf->primary_tag !~ /contig|gene|mRNA|CDS|exon/){
    #print "dubidu\n";
    my $source = $sf->source_tag || $sf->throw("source tag must be set: ".$sf->display_name);
    my $start = $sf->start || $sf->throw("start must be set or is zero: ".$sf->display_name);
    my $end = $sf->end || $sf->throw("end must be set");
    my $type = $sf->primary_tag || $sf->throw("primary_tag/type must be set: ".$sf->display_name);

    $id = "$source:$type:$seq_id:$start:$end";
  }
  return $id;
}

sub generateFeatureID {
  my ($sf,$elemnum,$parentid)=@_;

  my $fid=$parentid;
  my $ptag=$sf->primary_tag;

  if(exists($glob_ptag2idpostfix{$ptag})){
    # if there's a translation tag existing, two possibilities remain
    # 1) it's got a length, e.g. ".p", then append
    # 2) no lenth (like for 'gene' and 'contig'), then do nothing
    if(length($glob_ptag2idpostfix{$ptag})){
      $fid.=$glob_ptag2idpostfix{$ptag};
      $fid.=sprintf("%02d", $elemnum);
    }
  }else{
    if(!exists($glob_ptagotherpostfix{$ptag})){
      $glob_ptagotherpostfix{$ptag}=".o${glob_ptopcounter}x";
      $glob_ptopcounter+=1;
    }
    $fid.=$glob_ptagotherpostfix{$ptag};
    $fid.=sprintf("%02d", $elemnum);
  }

  return $fid;
}

sub moveTranslationToFasta {
  my ($sf)= @_;

  if($sf->has_tag('translation') ) {
    my ($aa) = $sf->get_tag_values("translation");
    $aa =~ s/\s//g;
    if($aa) {
      if($sf->has_tag('Name')) {
	push(@glob_fasta_sdata,$aa);
	push(@glob_fasta_sname,$sf->get_tag_values("Name"));
      }
      $sf->remove_tag("translation");
      $sf->add_tag_value("translation",$aa);
    }
  }
}

sub unflatten_seq {
  # from bp_genbank2gff3
  
  my $seq = shift;
  
  ## print "# working on $source_type:", $seq->accession, "\n"; 
  my $uh_oh = "Possible gene unflattening error with" .  $seq->accession_number .
    ": consult STDERR\n";
  
  eval {
    #        $unflattener->unflatten_seq( -seq => $seq,
    #                                     -use_magic => 1 );
    $unflattener->unflatten_seq( -seq => $seq,
				 -group_tag=>'locus_tag',
				 -use_magic => 1 );
    
  };
  
  # deal with unflattening errors
  if ( $@ ) {
    warn $seq->accession_number . " Unflattening error:\n";
    warn "Details: $@\n";
    print "# ".$uh_oh;
  }
  
  return 0 if !$seq || !$seq->all_SeqFeatures;
  
  # map feature types to the sequence ontology
  ## $tm->map_types_to_SO( -seq => $seq );
  $tm->map_types( -seq => $seq, -type_map => $FTSOmap, -undefined => "region" ); #dgg
  
  1;
}


sub maptags2gff {
  # verbatim from bp_genbank2gff3, except that rewriting note is done elsewhere
  my $f = shift;

  foreach my $tag (keys %glob_gb2gfftagmap) {
    if ($f->has_tag($tag)) {
      my $newtag= $glob_gb2gfftagmap{$tag};
      my @v= $f->get_tag_values($tag);
      $f->remove_tag($tag);
      $f->add_tag_value($newtag,@v);
    }
  }
}

