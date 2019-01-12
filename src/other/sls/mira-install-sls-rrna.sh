#!/bin/bash

BAITBIN=mirabait

# No user serviceable part below this line

function usage {
  echo "Installs files for filtering rRNA reads in MIRA."
  echo "Usage: $0 [targetdir] file.sls";
  echo
  echo "If no target directory is given, the files are installed in:"
  echo "<path-to-$BAITBIN>/../share/mira/mhs/"
}


# from http://stackoverflow.com/questions/7665/how-to-resolve-symbolic-links-in-a-shell-script
# but pretty much does what I implemented in C++ in MIRA :-)

function realpath {
    local r=$1; local t=$(readlink $r)
    while [ $t ]; do
        r=$(cd $(dirname $r) && cd $(dirname $t) && pwd -P)/$(basename $t)
        t=$(readlink $r)
    done
    echo $r
}

# Start here

# get the sed binary right: on OSX, we might need GNU sed, but let's try without first
SEDBIN=sed
if [ `uname` = 'Darwin' ] ; then
  SEDBIN=sed
fi

if [ $# -eq 0 ] ; then
  usage
  echo "Missing argument, need at least the name of the sls file." >& 2
  exit 2
fi

if [ x`which $BAITBIN` = x ] ; then
  BAITBIN=$(dirname $(realpath $0))/../bin/mirabait
  echo "Could not find the mirabait executable, trying $BAITBIN"
  if [ x`which $BAITBIN` = x ] ; then
    echo "Could not find any mira executable, this is needed for the installation."
    echo "Please make sure mirabait is in your path or change the BAITBIN variable in this script."
    exit 2;
  else
    echo "Cool, found it."
  fi
fi

if [ $# -eq 1 ] ; then
  if [ $1 = '-h' ] ; then
    usage
    exit 0
  fi
  SLSFILE=$1
  SHAREDIR=$(dirname $(realpath `which $BAITBIN`))
  SHAREDIR=${SHAREDIR}/../share/mira/mhs
fi

if [ $# -gt 2 ] ; then
  usage
  echo "Need exactly one or two arguments" >& 2
  exit 2
fi

if [ $# -eq 2 ] ; then
  if [ $1 = '-h' ] ; then
    usage
    exit 0
  fi
  SHAREDIR=$1
  SLSFILE=$2
fi

if [ ! -e $SLSFILE ]; then
  echo "Data file '$SLSFILE' not found in this directory?";
  exit 2;
fi

echo "using mirabait from $BAITBIN, installing data into $SHAREDIR"

if [ ! -d $SHAREDIR ]; then
  mkdir -p $SHAREDIR;
  if [ $? != 0 ]; then
    echo "Could not create $SHAREDIR. Missing permissions?";
    exit 2;
  fi
fi

# get away path (if any)
MSTMP=`echo $SLSFILE | $SEDBIN -e 's/.\+\/\([^/]\+\)$/\1/'`

SLSLIBNAME=`echo $MSTMP | $SEDBIN -e 's/\.sls.*//'`
KMERSIZE=`echo $SLSLIBNAME | $SEDBIN -e 's/.\+-\([0-9]\+\)-[0-9]\+$/\1/'`

echo "Will install $SLSLIBNAME with kmer size $KMERSIZE as MIRA default rRNA filter in"
echo "  $SHAREDIR"
echo
echo "Installing. This can take a minute or two, please be patient."

gunzip -c $SLSFILE | $SEDBIN -e 's/^/>x\'$'\n/' >$SLSLIBNAME.fasta
$BAITBIN -k $KMERSIZE -K $SHAREDIR/$SLSLIBNAME.mhs.gz -b $SLSLIBNAME.fasta >&mb.log
if [ $? != 0 ]; then
  echo "Some error occurred during execution of mirabait."
  echo "Please consult mb.log"
  exit 2
fi

rm $SLSLIBNAME.fasta mb.log

cd $SHAREDIR
ln -sf $SLSLIBNAME.mhs.gz filter_default_rrna.mhs.gz

echo "Done."
echo
echo "MIRA can now use the functionality to filter for rRNA during assemblies."
echo "MIRABAIT can now use the '-j rrna' option."
