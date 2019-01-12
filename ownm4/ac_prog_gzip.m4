dnl @synopsis AC_PROG_GZIP
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if gzip is installed. If gzip is installed, it set
dnl $gzip to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_GZIP],[
AC_CHECK_PROGS(gzip,[gzip],no)
export gzip;
if test $gzip = "no" ;
then
	AC_MSG_ERROR([Unable to find the gzip application]);
fi
AC_SUBST(gzip)
])
