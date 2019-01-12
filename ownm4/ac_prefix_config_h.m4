dnl @synopsis AC_PREFIX_CONFIG_H [(PREFIX [,ORIG-HEADER [,OUTPUT-HEADER]])]
dnl
dnl @obsoleted Use AX_PREFIX_CONFIG_H.
dnl
dnl takes the usual config.h generated header file; looks for each of
dnl the generated "#define SOMEDEF" lines, and prefixes the defined
dnl name (ie. makes it "#define PREFIX_SOMEDEF". The result is written
dnl to the output config.header file. The PREFIX is converted to
dnl uppercase for the conversions. If PREFIX is absent, $PACKAGE will
dnl be assumed. If the ORIG-HEADER is absent, "config.h" will be
dnl assumed. If the OUTPUT-HEADER is absent, "PREFIX-config.h" will be
dnl assumed.
dnl
dnl In most cases, the configure.in will contain a line saying
dnl
dnl         AC_CONFIG_HEADER(config.h)
dnl
dnl somewhere *before* AC_OUTPUT and a simple line saying
dnl
dnl        AC_PREFIX_CONFIG_HEADER
dnl
dnl somewhere *after* AC_OUTPUT.
dnl
dnl example:
dnl
dnl   AC_INIT(config.h.in)        # config.h.in as created by "autoheader"
dnl   AM_INIT_AUTOMAKE(testpkg, 0.1.1)   # "#undef VERSION" and "PACKAGE"
dnl   AM_CONFIG_HEADER(config.h)         #                in config.h.in
dnl   AC_MEMORY_H                        # "#undef NEED_MEMORY_H"
dnl   AC_C_CONST_H                       # "#undef const"
dnl   AC_OUTPUT(Makefile)                # creates the "config.h" now
dnl   AC_PREFIX_CONFIG_H                 # creates "testpkg-config.h"
dnl         and the resulting "testpkg-config.h" contains lines like
dnl   #define TESTPKG_VERSION "0.1.1"
dnl   #define TESTPKG_NEED_MEMORY_H 1
dnl   #define TESTPKG_const const
dnl
dnl   and this "testpkg-config.h" can be installed along with other
dnl   header-files, which is most convenient when creating a shared
dnl   library (that has some headers) where some functionality is
dnl   dependent on the OS-features detected at compile-time. No
dnl   need to invent some "testpkg-confdefs.h.in" manually. :-)
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_PREFIX_CONFIG_H],
[changequote(<<, >>)dnl
ac_prefix_conf_PKG=`echo ifelse($1, , $PACKAGE, $1)`
ac_prefix_conf_PRE=`echo $ac_prefix_conf_PKG | tr 'abcdefghijklmnopqrstuvwxyz-' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_'`
ac_prefix_conf_PRE=`echo $ac_prefix_conf_PRE | sed -e '/^[0-9]/s/^/_/'`
ac_prefix_conf_INP=`echo ifelse($2, , config.h, $2)`
ac_prefix_conf_OUT=`echo ifelse($3, , $ac_prefix_conf_PKG-$ac_prefix_conf_INP, $3)`
ac_prefix_conf_DEF=`echo _$ac_prefix_conf_OUT | tr 'abcdefghijklmnopqrstuvwxyz./,-' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ____'`
changequote([, ])dnl
if test -z "$ac_prefix_conf_PKG" ; then
   AC_MSG_ERROR([no prefix for _PREFIX_CONFIG_H])
else
  AC_MSG_RESULT(creating $ac_prefix_conf_OUT - prefix $ac_prefix_conf_PRE for $ac_prefix_conf_INP defines)
  if test -f $ac_prefix_conf_INP ; then
    echo '#ifndef '$ac_prefix_conf_DEF >$ac_prefix_conf_OUT
    echo '#define '$ac_prefix_conf_DEF' 1' >>$ac_prefix_conf_OUT
    echo ' ' >>$ac_prefix_conf_OUT
    echo /'*' $ac_prefix_conf_OUT. Generated automatically at end of configure. '*'/ >>$ac_prefix_conf_OUT

    echo 's/#undef  */#undef '$ac_prefix_conf_PRE'_/' >conftest.sed
    echo 's/#define  *\([A-Za-z0-9_]*\)\(.*\)/#ifndef '$ac_prefix_conf_PRE"_\\1 \\" >>conftest.sed
    echo '#define '$ac_prefix_conf_PRE"_\\1 \\2 \\" >>conftest.sed
    echo '#endif/' >>conftest.sed
    sed -f conftest.sed $ac_prefix_conf_INP >>$ac_prefix_conf_OUT
    echo ' ' >>$ac_prefix_conf_OUT
    echo '/*' $ac_prefix_conf_DEF '*/' >>$ac_prefix_conf_OUT
    echo '#endif' >>$ac_prefix_conf_OUT
  else
    AC_MSG_ERROR([input file $ac_prefix_conf_IN does not exist, dnl
    skip generating $ac_prefix_conf_OUT])
  fi
  rm -f conftest.*
fi])
