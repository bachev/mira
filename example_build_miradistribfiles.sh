#!/bin/sh

# This script works out of the box for a GitHub clone on
#  (K)Ubuntu 18.04 and MacOS X once all dependencies are met!
#
# It will create distribution source and binary tar files

./bootstrap.sh

systemuname=`uname`
if test "${systemuname}" = "Darwin"; then
  ./configure --enable-mirastatic --with-brew
else
  ./configure --enable-mirastatic
fi

make clean
make -j 2
make docs
make -C distribution clean
make dist
make distrib
mv mira*.tar.gz distribution
mv mira*.tar.bz2 distribution

echo
echo
echo "Done. If all went right, you will find tar files ready in the"
echo "  distribution"
echo "directory. Enjoy."
echo
echo
