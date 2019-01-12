dnl @synopsis AC_DEFINE_DIR_(VARNAME, DIR [, DESCRIPTION])
dnl
dnl @obsoleted Use AC_DEFINE_DIR instead.
dnl
dnl This macro _AC_DEFINEs VARNAME to the expansion of the DIR
dnl variable, taking care of fixing up ${prefix} and such.
dnl
dnl Note that the 3 argument form is only supported with autoconf 2.13
dnl and later (i.e. only where _AC_DEFINE supports 3 arguments).
dnl
dnl Examples:
dnl
dnl    AC_DEFINE_DIR_(DATADIR, datadir)
dnl    AC_DEFINE_DIR_(PROG_PATH, bindir, [Location of installed binaries])
dnl
dnl This macro is based on Alexandre Oliva's AC_DEFINE_DIR.
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_DEFINE_DIR_], [
  test "x$prefix" = xNONE && prefix="$ac_default_prefix"
  test "x$exec_prefix" = xNONE && exec_prefix='${prefix}'
  ac_define_dir=`eval "echo [$]$2"`
  ac_define_dir=`eval "echo [$]ac_define_dir"`
  ifelse($3, ,dnl
    AC_DEFINE_UNQUOTED($1, "$ac_define_dir"),dnl
    AC_DEFINE_UNQUOTED($1, "$ac_define_dir", $3))
])
