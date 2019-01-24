#!/bin/bash

# No user serviceable part below this line

function usage {
  echo "Downloads SLS file for for MIRA."
  echo "Usage: $0 file.sls";
}


if [ $# -eq 0 ] ; then
  usage
  echo "Missing argument, need at least the name of the sls file." >& 2
  exit 2
fi

if [ $# -eq 1 ] ; then
  if [ $1 = '-h' ] ; then
    usage
    exit 0
  fi
  SLSFILE=$1
fi

if [ $# -gt 1 ] ; then
  usage
  echo "Need exactly one argument" >& 2
  exit 2
fi

loopi="0"
until gzip -t ${SLSFILE} ; do
  rm ${SLSFILE}
  echo "Need to download '$SLSFILE'";
  loopi=$[${loopi}+1]
  if [ ${loopi} -gt 1 ]; then
    echo "Maximum number of retries reached, couldn't download SLS $1"
    exit 2;
  fi
  wget https://sourceforge.net/projects/mira-assembler/files/slsstore/${SLSFILE}
done
