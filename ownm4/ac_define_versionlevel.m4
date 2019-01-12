dnl @synopsis AC_DEFINE_VERSIONLEVEL(VARNAME [,VERSION [, DESCRIPTION]])
dnl
dnl if the VERSION is ommitted, shellvar $VERSION is used as defined by
dnl AM_INIT_AUTOMAKE's second argument.
dnl
dnl The versionlevel is the numeric representation of the given version
dnl string, thereby assuming the inputversion is a string with
dnl (maximal) three decimal numbers seperated by "."-dots. A "-patch"
dnl adds a percent.
dnl
dnl typical usage: AM_INIT_AUTOMAKE(mypkg,4.12.3)
dnl AC_DEFINE_VERSIONLEVEL(MYPKG_VERSION)
dnl
dnl the config.h created from autoheader's config.h.in will contain...
dnl /* mypkg versionlevel */ #define MYPKG_VERSION 4120003
dnl
dnl the MYKG_VERSION will be defined as both a shell-variable and
dnl AC_DEFINE
dnl
dnl examples:
dnl
dnl        3.0-beta1     3000001
dnl        3.1           3010000
dnl        3.11          3110000
dnl        3.11-dirpatch 3111000
dnl        3.11-patch6   3110006
dnl        2.2.18        2020018
dnl        2.0.112       2000112
dnl        2.4.2         2040002
dnl        2.4.2-pre     2040003
dnl        2.4.2-pre5    2040003
dnl        5.0-build125  5000125
dnl        5.0           5000000
dnl        0.30.17       30017
dnl
dnl @category Misc
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_DEFINE_VERSIONLEVEL],
[
ac_versionlevel_strdf=`echo ifelse($2, , $VERSION, $2) | sed -e 's:[[A-Z-]]*:.:' -e 's:[[^0-9.]]::g' -e 's:^[[.]]*::'`
AC_MSG_CHECKING(versionlevel $ac_versionlevel_strdf)
case $ac_versionlevel_strdf in
 *.*.*.|*.*.*.*) :
 ac_versionlevel_major=`echo $ac_versionlevel_strdf`
 ac_versionlevel_minor=`echo $ac_versionlevel_major | sed -e 's/[[^.]]*[[.]]//'`
 ac_versionlevel_patch=`echo $ac_versionlevel_minor | sed -e 's/[[^.]]*[[.]]//'`
 ac_versionlevel_major=`echo $ac_versionlevel_major | sed -e 's/[[.]].*//'`
 ac_versionlevel_minor=`echo $ac_versionlevel_minor | sed -e 's/[[.]].*//'`
 ac_versionlevel_patch=`echo $ac_versionlevel_patch | sed -e 's/[[.]].*//'`
 $1=`expr $ac_versionlevel_major '*' 1000000 \
        + $ac_versionlevel_minor '*'   10000 \
        + $ac_versionlevel_patch \
	+ 1` ;;
 *.*.*) :
 ac_versionlevel_major=`echo $ac_versionlevel_strdf`
 ac_versionlevel_minor=`echo $ac_versionlevel_major | sed -e 's/[[^.]]*[[.]]//'`
 ac_versionlevel_patch=`echo $ac_versionlevel_minor | sed -e 's/[[^.]]*[[.]]//'`
 ac_versionlevel_major=`echo $ac_versionlevel_major | sed -e 's/[[.]].*//'`
 ac_versionlevel_minor=`echo $ac_versionlevel_minor | sed -e 's/[[.]].*//'`
 ac_versionlevel_patch=`echo $ac_versionlevel_patch | sed -e 's/[[.]].*//'`
 $1=`expr $ac_versionlevel_major '*' 1000000 \
        + $ac_versionlevel_minor '*'   10000 \
        + $ac_versionlevel_patch`               ;;
 *.*.) :
 ac_versionlevel_major=`echo $ac_versionlevel_strdf`
 ac_versionlevel_minor=`echo $ac_versionlevel_major | sed -e 's/[[^.]]*[[.]]//'`
 ac_versionlevel_major=`echo $ac_versionlevel_major | sed -e 's/[[.]].*//'`
 ac_versionlevel_minor=`echo $ac_versionlevel_minor | sed -e 's/[[.]].*//'`
 ac_versionlevel_patch=0
 $1=`expr $ac_versionlevel_major '*' 1000000 \
        + $ac_versionlevel_minor '*'   10000 \
	+ 1000 \
        + $ac_versionlevel_patch`               ;;
 *.*) :
 ac_versionlevel_major=`echo $ac_versionlevel_strdf`
 ac_versionlevel_minor=`echo $ac_versionlevel_major | sed -e 's/[[^.]]*[[.]]//'`
 ac_versionlevel_major=`echo $ac_versionlevel_major | sed -e 's/[[.]].*//'`
 ac_versionlevel_minor=`echo $ac_versionlevel_minor | sed -e 's/[[.]].*//'`
 ac_versionlevel_patch=0
 $1=`expr $ac_versionlevel_major '*' 1000000 \
        + $ac_versionlevel_minor '*'   10000 \
        + $ac_versionlevel_patch`               ;;
 *.) :
 ac_versionlevel_major=0
 ac_versionlevel_minor=`echo $ac_versionlevel_strdf`
 ac_versionlevel_minor=`echo $ac_versionlevel_minor | sed -e 's/[[.]].*//'`
 ac_versionlevel_patch=0
 $1=`expr $ac_versionlevel_major '*' 1000000 \
        + $ac_versionlevel_minor '*'   10000 \
	+ 1000 \
        + $ac_versionlevel_patch`               ;;
 *) :
 ac_versionlevel_major=0
 ac_versionlevel_minor=`echo $ac_versionlevel_strdf`
 ac_versionlevel_minor=`echo $ac_versionlevel_minor | sed -e 's/[[.]].*//'`
 ac_versionlevel_patch=0
 $1=`expr $ac_versionlevel_major '*' 1000000 \
        + $ac_versionlevel_minor '*'   10000 \
        + $ac_versionlevel_patch`               ;;
esac
AC_MSG_RESULT($[$1])
AC_DEFINE_UNQUOTED( $1, $[$1], ifelse( $3, , $PACKAGE versionlevel, $3))
])
