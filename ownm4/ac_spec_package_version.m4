dnl @synopsis AC_SPEC_PACKAGE_VERSION(rpmspecfile)
dnl
dnl @obsoleted Use AX_SPEC_PACKAGE_VERSION.
dnl
dnl set PACKAGE and VERSION from the defines in the given specfile
dnl default to basename and currentde if rpmspecfile is not found
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_SPEC_PACKAGE_VERSION],[dnl
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
      VERSION=`date +%Y.%m.%d`
    fi
    AC_MSG_RESULT( $PACKAGE $VERSION $k )
  fi
])
