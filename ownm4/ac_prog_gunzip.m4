dnl @synopsis AC_PROG_GUNZIP
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if gunzip is installed. If gunzip is installed, it
dnl set $gunzip to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_GUNZIP],[
AC_CHECK_PROGS(gunzip,[gunzip],no)
export gunzip;
if test $gunzip = "no" ;
then
	AC_MSG_ERROR([Unable to find the gunzip application]);
fi
AC_SUBST(gunzip)
])
