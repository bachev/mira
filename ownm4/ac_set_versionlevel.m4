dnl @synopsis AC_SET_VERSIONLEVEL(VARNAME [,VERSION])
dnl
dnl if the VERSION is ommitted, shellvar $VERSION is used as defined by
dnl AM_INIT_AUTOMAKE's second argument.
dnl
dnl The versionlevel is the numeric representation of the given version
dnl string, thereby assuming the inputversion is a string with
dnl (maximal) three decimal numbers seperated by "."-dots. A "-patch"
dnl adds a percent.
dnl
dnl typical usage:
dnl
dnl  AM_INIT_AUTOMAKE(mypkg,4.12.3)
dnl  AC_SET_VERSIONLEVEL(MYPKG_VERSION)
dnl  AC_DEFINE_UNQUOTED(MYPKG_VERSION, $MYPKG_VERSION, [package version])
dnl
dnl (this macro shall superced AC_DEFINE_VERSIONLEVEL at some day)
dnl
dnl the version code has three digits per part which I feel is the most
dnl natural encoding - it makes it easier to be printf'd anyway
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

AC_DEFUN([AC_SET_VERSIONLEVEL],
[dnl
m4_pushdef(LVL, $1_LEVEL)
m4_pushdef(MJR, $1_MAJOR)
m4_pushdef(MNR, $1_MINOR)
m4_pushdef(MCR, $1_MICRO)
LVL=`echo ifelse($2, , $VERSION, $2) | sed -e 's:[[A-Z-]]*:.:g' -e 's:[[^0-9.]]::g' -e 's:[[.]]*:.:g' -e 's:^[[.]]*::'`
AC_MSG_CHECKING( $1 versionlevel $LVL)
case $LVL in
 *.*.*.|*.*.*.*|*.*.*) :
 MJR=`echo $LVL`
 MNR=`echo $MJR | sed -e 's/[[^.]]*[[.]]//'`
 MCR=`echo $MNR | sed -e 's/[[^.]]*[[.]]//'`
 MJR=`echo $MJR | sed -e 's/[[.]].*//'`
 MNR=`echo $MNR | sed -e 's/[[.]].*//'`
 MCR=`echo $MCR | sed -e 's/[[.]].*//'`
 ;;
 *.*.|*.*) :
 MJR=`echo $LVL`
 MNR=`echo $MJR | sed -e 's/[[^.]]*[[.]]//'`
 MJR=`echo $MJR | sed -e 's/[[.]].*//'`
 MNR=`echo $MNR | sed -e 's/[[.]].*//'`
 MCR=0
 ;;
 *.) :
 MJR=0
 MNR=`echo $LVL`
 MNR=`echo $MNR | sed -e 's/[[.]].*//'`
 MCR=0
 ;;
esac
# we trust sed greedy-match backtracking to extract the last three digits from each part, forming a nine-digit
$1=`echo 000$MJR.000$MNR.000$MCR | sed -e 's:\\(...\\)[[.]][[^.]]*\\((...\\))[[.]][[^.*]]\\((...\\)):\\1\\2\\3 -e 's:^0*::''
AC_MSG_RESULT($[$1] ($MJR,$MNR,$MCR)
dnl AC_DEFINE_UNQUOTED( $1, $[$1], ifelse( $3, , $PACKAGE versionlevel, $3))
m4_popdef(MCR)
m4_popdef(MNR)
m4_popdef(MJR)
m4_popdef(LVL)
])
