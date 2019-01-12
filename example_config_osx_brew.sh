#!/bin/sh

if [ ! -e ./configure ]; then
  echo "./configure does not exist. You might need to run ./bootstrap.sh"
  exit 1
fi

./configure --with-brew
rc=$?
if [ $rc != 0 ] ; then
  echo "Oooops, ./configure failed? Not good."
  exit $rc
fi
