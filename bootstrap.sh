#!/bin/sh

rm -rf Makefile Makefile.in aclocal.m4 config.status config.h config.h.in config.log configure autom4te.cache config.guess config.sub missing libtool
rm -rf compile depcomp install-sh ltmain.sh
rm -rf config m4

# Gaaaaaaaaaaaaah. Something in the auto* system breaks if
#  the 'required' README file is missing

touch README

if [ -n "`which glibtoolize 2>/dev/null`" ]
then
  GLIBTOOLIZE="glibtoolize"
else
  GLIBTOOLIZE="libtoolize"
fi

echo "Copying m4"
cp -a ownm4 m4

echo "Running aclocal (can take a while)"
aclocal -I m4
rc=$?
if [ $rc != 0 ] ; then
  echo
  echo "Oooops, aclocal failed? Maybe a wrong libtool version?"
  exit $rc
fi

echo "Running ${GLIBTOOLIZE}"
$GLIBTOOLIZE --force --install

echo "Running autoconf"
autoconf
echo "Running autoheader"
autoheader
echo "Running automake"
automake --add-missing

chmod +x config/install-sh

echo "Done. You can now './configure'"
