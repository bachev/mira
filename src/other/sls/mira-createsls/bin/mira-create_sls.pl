#!/usr/bin/env perl

#
#    Copyright (C) 2015- Bastien Chevreux
#
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
use FindBin;
use Cwd;

my $glob_version="1";
my $glob_author="Bastien Chevreux";

my(@glob_options,
   $opt_source,
   $opt_kmersize,
   $opt_filter,
   $opt_force,
   $opt_datadir,
   $opt_sednameaway,
   $opt_cleanup,

   $opt_curlfU,
   $opt_curlfx,
 );

my $glob_exe = $FindBin::RealScript;
my $glob_bindir = $FindBin::RealBin;
my $glob_dbdir = "$FindBin::RealBin/../db";
my $glob_cachedir = "$FindBin::RealBin/../db/cache";
my $glob_defaultdatadir="$FindBin::RealBin/../db/data";

my ($glob_univecraw,
    $glob_univecclean,
    $glob_rrnadataraw,
    $glob_rrnadatavect,
    $glob_rrnadatanovect,
    $glob_rrnadataclean,
    $glob_kmerstats,
    $glob_kmerstatsf,
    $glob_kmerstatsf2,
    $glob_dbg1,
    $glob_dbg2,
    $glob_tmpsls,
    $glob_finalsls
);

my @glob_rfamfiles= (
  'RF00001.fa.gz',
  'RF00002.fa.gz',
  'RF00177.fa.gz',
  'RF01959.fa.gz',
  'RF01960.fa.gz',
  'RF02540.fa.gz',
  'RF02541.fa.gz',
  'RF02542.fa.gz',
  'RF02543.fa.gz'
);

#="/scratch1/tmp/bachtmp/RFAM12/fasta_files";

my $glob_step=0;


##############################################################################################

##############################################################################################

setOptions();

checkTools();

if($opt_cleanup){
  `rm -rf $glob_defaultdatadir`;
  `rm -rf $glob_cachedir`;
  `rm -rf $glob_dbdir/*`;
  exit(0);
}

setFileNames();
ensureDataFiles();
#exit(0);
main();
exit(0);


sub checkTools {
  my %tools = (
    'mirabait' => {
      NEEDED => 1,
      DOWNLOAD => "http://sourceforge.net/projects/mira-assembler/"
    },
    'barrnap' => {
      NEEDED => 1,
      DOWNLOAD => "http://www.vicbioinformatics.com/software.barrnap.shtml"
    }
  );

  for my $tn (sort keys %tools) {
    my $exe=findExe($tn);
    if(!$exe) {
      err("Could not find $tn in your path. Please download and install from\n  ".$tools{$tn}->{DOWNLOAD}."\nAborting.\n");
    }
  }
}

sub findExe {
  my($bin) = shift;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);
    return $exe if -x $exe;
  }
  return;
}

sub ensureDataFiles {
  my $missing=0;

  mkdir $glob_dbdir unless -d $glob_dbdir;
  mkdir $glob_cachedir unless -d $glob_cachedir;
  mkdir $glob_defaultdatadir unless -d $glob_defaultdatadir;
  mkdir $opt_datadir unless -d $opt_datadir;

  if(-e $glob_cachedir) {
    if(!-d $glob_cachedir) {
      die "$glob_cachedir already exists but is not a directory.";
    }
  }else{
    unless(mkdir $glob_cachedir) {
      die "Unable to create $glob_cachedir\n";
    }
  }

  my @dldata=walkNeededFiles();
  if(scalar(@dldata)){
    downloadFiles(@dldata);
  }

  @dldata=walkNeededFiles();
  if(scalar(@dldata)){
    print STDERR "Some files are missing, cannot continue.\n";
    exit(10);
  }
}


sub walkNeededFiles {
  my @dldata;

  if(! -e $glob_univecraw){
    print "Need to download NCBI Univec database\n";
    push @dldata, "ftp", "ftp.ncbi.nlm.nih.gov", "ftp:ftp", "pub/UniVec/UniVec", $glob_univecraw;
  }

  if($opt_source eq "rfam"){
    foreach my $rff (@glob_rfamfiles){
      my $fpath="$opt_datadir/${rff}";
      my $fileok=0;

      if(-e "${fpath}.ok"){
	$fileok=1;
      }else{
	if(-e $fpath){
	  print "Checking integrity of $opt_source $rff\n";
	  if(system("gzip -t $fpath")==0){
	    `touch ${fpath}.ok`;
	    $fileok=1;
	  }else{
	    print "Integrity check of $fpath failed.\n";
	  }
	}
      }

      if(!$fileok) {
	print "Need to download $opt_source $rff\n";
	# protocol, username:passwd at site, dlpath, localpath
	#push @dldata, "ftp", "ftp:ftp", "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/$rff", $fpath;
	# protocol, site, username:passwd at site, dlpath, localpath
	push @dldata, "ftp", "ftp.ebi.ac.uk", "ftp:ftp", "pub/databases/Rfam/CURRENT/fasta_files/$rff", $fpath;
      }
    }
  }

  return @dldata;
}

sub downloadFiles {
  my @dldata=@_;

  while (my @dlset = splice(@dldata, 0, 5)) {
    my ($protocol,$site,$auth,$source,$dltarget)=@dlset;
    my $curlcmd="curl";
    if($protocol eq "ftp"){
      if(length($opt_curlfU)){
	my $tmp=$opt_curlfU;
	$tmp =~ s/%h/$site/;
	$curlcmd.=" -U $tmp";
      }
      $curlcmd.=" -u $auth";
      if(length($opt_curlfx)){
	$curlcmd.=" -x $opt_curlfx";
      }
      $curlcmd.=" ftp://$site/$source -o $dltarget";
    }else{
      print "Protocol $protocol not supported yet.\n";
      exit 2;
    }

    runcmd($curlcmd);
  }

}

sub setFileNames {
  $glob_univecraw = "$opt_datadir/UniVec.fasta";
  $glob_univecclean = "$glob_cachedir/UniVec_rrnamasked.fasta";
  $glob_rrnadataraw = "$glob_cachedir/${opt_source}_rrna_ref.fasta";
  $glob_rrnadatavect = "$glob_cachedir/${opt_source}_rrna_ref_vect.fasta";
  $glob_rrnadatanovect = "$glob_cachedir/${opt_source}_rrna_ref_novect.fasta";
  $glob_rrnadataclean = "$glob_cachedir/${opt_source}_rrna_ref_clean.fasta";
  $glob_kmerstats = "$glob_cachedir/${opt_source}_rrnakmerstats-$opt_kmersize.mhs.gz";
  $glob_kmerstatsf = "$glob_cachedir/${opt_source}_rrnakmerstats-$opt_kmersize-$opt_filter.mhs.gz";
  $glob_kmerstatsf2 = "$glob_cachedir/${opt_source}_rrnakmerstats-$opt_kmersize-$opt_filter.2.mhs.gz";
  $glob_dbg1 = "$glob_cachedir/${opt_source}_dbg1-$opt_kmersize-$opt_filter.fasta";
  $glob_dbg2 = "$glob_cachedir/${opt_source}_dbg2-$opt_kmersize-$opt_filter.fasta";
  $glob_tmpsls = "$glob_cachedir/${opt_source}_rrna-$opt_kmersize-$opt_filter.tmp.sls";

  my $cwdir = getcwd;
  $glob_finalsls = "$cwdir/${opt_source}_rrna-$opt_kmersize-$opt_filter.sls.gz";
}

sub main {
  my $sn=$opt_source;

  my $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    if($opt_source eq "rfam"){
      if($opt_sednameaway){
	runcmd("gunzip -c $opt_datadir/RF00001.fa.gz $opt_datadir/RF00002.fa.gz $opt_datadir/RF00177.fa.gz $opt_datadir/RF01959.fa.gz $opt_datadir/RF01960.fa.gz $opt_datadir/RF02540.fa.gz $opt_datadir/RF02541.fa.gz $opt_datadir/RF02542.fa.gz $opt_datadir/RF02543.fa.gz | $glob_bindir/dedup.pl | sed -e '/^>/c\\>x' >$glob_rrnadataraw");
      }else{
	runcmd("gunzip -c $opt_datadir/RF00001.fa.gz $opt_datadir/RF00002.fa.gz $opt_datadir/RF00177.fa.gz $opt_datadir/RF01959.fa.gz $opt_datadir/RF01960.fa.gz $opt_datadir/RF02540.fa.gz $opt_datadir/RF02541.fa.gz $opt_datadir/RF02542.fa.gz $opt_datadir/RF02543.fa.gz | $glob_bindir/dedup.pl >$glob_rrnadataraw");
      }
    }else{
      print "Source $opt_source not supported yet\n";
      exit;
    }
    `touch  $okfile`;
  }
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    if (!-r $glob_univecraw) {
      err("Missing file $glob_univecraw\n");
    }
    runcmd("$glob_bindir/mask_potential_rrna.pl $glob_univecraw >$glob_univecclean");
    `touch  $okfile`;
  }
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    runcmd("mirabait -d -l 0 -k 31 -n 1 -b $glob_univecclean -o $glob_rrnadatavect -O $glob_rrnadatanovect $glob_rrnadataraw");
    `touch $okfile`;
  }
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    runcmd("sed -e 's/[aA]\\{10,\\}/aaaaaaanaaaaaaa/g' -e 's/[cC]\\{10,\\}/cccccccnccccccc/g' -e 's/[gG]\\{10,\\}/gggggggcggggggg/g' -e 's/[tT]\\{10,\\}/tttttttnttttttt/g' $glob_rrnadatanovect >$glob_rrnadataclean");
    `touch $okfile`;
  }
  $sn.="-$opt_kmersize";
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    runcmd("mirabait -d -k $opt_kmersize -K $glob_kmerstats -b $glob_rrnadataclean");
    `touch $okfile`;
  }
  $sn.="-$opt_filter";
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    if($opt_filter == 1) {
      runcmd("ln -s $glob_kmerstats $glob_kmerstatsf");
    }else{
      runcmd("miramer -j filter -o $glob_kmerstatsf $opt_filter $glob_kmerstats");
    }
    `touch $okfile`;
  }
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    runcmd("miramer -j dbg -o $glob_dbg1 $glob_kmerstatsf");
    `touch $okfile`;
  }
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    runcmd("grep -v '>' ${glob_dbg1}.bymhc | sort >$glob_tmpsls");
    `touch $okfile`;
  }
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    runcmd("mirabait -k $opt_kmersize -K $glob_kmerstatsf2 -b fna::${glob_dbg1}.bymhc");
    `touch $okfile`;
  }
  $okfile=newStep($sn);
  if (!-e $okfile || $opt_force) {
    runcmd("miramer -j dbg -o $glob_dbg2 $glob_kmerstatsf2");
    `touch $okfile`;
  }

  runcmd("grep -v '>' ${glob_dbg2}.bymhc | sort | gzip -9 >$glob_finalsls");
}

###########################

sub newStep {
  my $retval="$glob_cachedir/step_".$glob_step;
  if(scalar(@_)){
    $retval.="-".join("-",@_);
  }
  $retval.=".ok";
  $glob_step+=1;
  print "### $retval\n";
  return $retval;
}

###########################

sub msg {
  print @_;
}

sub err {
  print @_;
  exit(2);
}

sub runcmd {
  msg("Running:", @_, "\n");
  system(@_)==0 or err("Could not run command:", @_);
}


###########################

sub setOptions {
  use Getopt::Long;

  @glob_options = (
    'Info:',
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"version", VAR=>\&version, DESC=>"Print version and exit"},
    'Options:',
    {OPT=>"source=s",   VAR=>\$opt_source,   DEFAULT=>'rfam',  DESC=>"rRNA database source. Currently only 'rfam'"},
    {OPT=>"datadir=s",  VAR=>\$opt_datadir,  DEFAULT=>'',      DESC=>"Directory to find rRNA database files in."},
    {OPT=>"kmer=i",     VAR=>\$opt_kmersize, DEFAULT=>21,      DESC=>"kmer size"},
    {OPT=>"filter=i",   VAR=>\$opt_filter,   DEFAULT=>12,      DESC=>"kmer minimum occurrence filter."},
    {OPT=>"force=i",    VAR=>\$opt_force,    DEFAULT=>0,       DESC=>"Force recalculation, do not use cached values."},
    'Proxy options:',
    {OPT=>"curlfU=s",   VAR=>\$opt_curlfU,   DEFAULT=>'',      DESC=>"for curl FTP: -U variable to pass proxy"},
    {OPT=>"curlfx=s",   VAR=>\$opt_curlfx,   DEFAULT=>'',      DESC=>"for curl FTP: -x variable to pass proxy"},
    'Cleanup options:',
    {OPT=>"cleanup!",   VAR=>\$opt_cleanup,  DEFAULT=>0,       DESC=>"cleanup temporary files and exit"},
  );

#    {OPT=>"proxyu=s",   VAR=>\$opt_proxyu,   DEFAULT=>'',      DESC=>"Proxy user (HTTP and FTP)"},
#    {OPT=>"proxyp=s",   VAR=>\$opt_proxyp,   DEFAULT=>'',      DESC=>"Proxy password (HTTP and FTP)"},
#    {OPT=>"hproxyu=s",  VAR=>\$opt_hproxyu,  DEFAULT=>'',      DESC=>"HTTP proxy user"},
#    {OPT=>"hproxyp=s",  VAR=>\$opt_hproxyp,  DEFAULT=>'',      DESC=>"HTTP proxy password"},
#    {OPT=>"fproxyu=s",  VAR=>\$opt_fproxyu,  DEFAULT=>'',      DESC=>"FTP proxy user"},
#    {OPT=>"fproxyp=s",  VAR=>\$opt_fproxyp,  DEFAULT=>'',      DESC=>"FTP proxy password"},

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} grep { ref } @glob_options) || usageerr();

  # Now setup default values.
  foreach (@glob_options) {
    if (ref $_
	  && defined($_->{DEFAULT})
	    && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }

  # set default data dir
  if(length($opt_datadir) == 0){
    $opt_datadir=$glob_defaultdatadir;
  }

}

sub usageerr {
  print STDERR "Use '$glob_exe  -h' for help.\n";
  exit(10);
}

sub usage {
  print
    "Name:\n  ", $glob_exe, " $glob_version by $glob_author\n",
    "Synopsis:\n  Create SLS files for installation with MIRA / MIRABAIT\n",
    "\nUsage:\n  $glob_exe [options]\n\n";
  foreach (@glob_options) {
    if (ref) {
      my $def = defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
      $def = ($def ? ' (default OFF)' : '(default ON)') if $_->{OPT} =~ m/!$/;
      my $opt = $_->{OPT};
      $opt =~ s/!$//;
      $opt =~ s/=s$/ [X]/;
      $opt =~ s/=i$/ [N]/;
      $opt =~ s/=f$/ [n.n]/;
      printf "  --%-15s %s%s\n", $opt, $_->{DESC}, $def;
    }
    else {
      print "$_\n";
    }
  }

  print
    "\nSubstitution strings for --curlfU and --curlfx\n",
    "$glob_exe uses curl to download data. Two curl-options are important\n",
    "for this: -X and -u (see curl manpage for more information on these).\n",
    "As there is no standard (at least for FTP proxies), both --curlfU and --curlfx\n",
    "understand substitution strings which are generated on the fly:\n",
    "  \%h host to be connected to\n",
    "  (others?)\n",
    "\nE.g., if you wanted to download, on the command line with curl, this file:\n",
    "  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF02540.fa.gz\n",
    "and if your proxy user account is '12345', the password is '67890',\n",
    "the FTP proxy is at 'http://myproxy:8080' and the proxy wants the FTP address\n",
    "as part of the login:password authentication like in the following curl\n",
    "command line:\n",
    "  curl -U 12345:\@ftp.ebi.ac.uk:67890\n",
    "       -x http://myproxy:8080\n",
    "       -u ftp:ftp\n",
    "       ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF02540.fa.gz\n",
    "(where -u ftp:ftp is the login:password of ftp.ebi.ac.uk), then you need to\n",
    "call $glob_exe like this:\n",
    "  $glob_exe --curlfU=12345:@%h:67890 --curlfx=http://myproxy:8080\n\n";

  exit(0);
}
