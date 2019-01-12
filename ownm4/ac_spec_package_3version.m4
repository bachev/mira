dnl @synopsis AC_SPEC_PACKAGE_3VERSION(rpmspecfile)
dnl
dnl @obsoleted Use AX_SPEC_PACKAGE_VERSION.
dnl
dnl set PACKAGE and VERSION from the defines in the given specfile
dnl default to basename and currentdate if rpmspecfile is not found if
dnl the VERSION-number from the spec-file is shorter than 3 digits then
dnl additional numbers are taken from the current date using entries
dnl from `date` as %y and %W%w since libtool version numbers can be as
dnl max only be 3 digits. The year counts from 1900.
dnl
dnl spec example:
dnl
dnl     Name: testprog
dnl     Version: 2
dnl
dnl result: (on 1. April 2002, being monday of 13th week)
dnl
dnl     VERSION="2.102.131"
dnl     PACKAGE="testprog"
dnl
dnl See also AC_SET_RELEASEINFO_VERSIONINFO for the use of a 3VERSION.
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_SPEC_PACKAGE_3VERSION],[dnl
  pushdef([specfile], ac_spec_package_version_file)
  specfile=`basename $1`
  AC_MSG_CHECKING( $specfile package version)
  if test -z "$1"; then
    AC_MSG_ERROR( no rpm spec file given )
  else
    # find specfile
    for i in $1 $srcdir/$1 $srcdir/../$1 ; do
      if test -f "$i" ; then
        specfile="$i"
        break
      fi
    done
    if test ! -f $specfile ; then
      k="w/o spec"
    else
      if test -z "$PACKAGE" ; then
        i=`grep -i '^name:' $specfile | head -1 | sed -e 's/.*://'`
	PACKAGE=`echo $i | sed -e 's/ /-/'`
      fi
      if test -z "$VERSION" ; then
        i=`grep -i '^version:' $specfile | head -1 | sed -e 's/.*://'`
	VERSION=`echo $i | sed -e 's/ /-/'`
      fi
    fi
    if test -z "$PACKAGE" ; then
      PACKAGE=`basename $specfile .spec`
    fi
    if test -z "$VERSION" ; then
      VERSION=`date +0.1%y.%W%w`
    fi
    case "$VERSION" in
    *.*.*) ;;
    *.*) VERSION="$VERSION."`date +1%y%W%w` ;;
    *) VERSION="$VERSION."`date +1%y.%W%w` ;;
    esac
    VERSION=`echo $VERSION | sed -e 's/[[.]]0/./g'`
    AC_MSG_RESULT( $PACKAGE $VERSION $k )
  fi
])
