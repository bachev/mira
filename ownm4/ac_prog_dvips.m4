dnl @synopsis AC_PROG_DVIPS
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if dvips is installed. If dvips is installed, it
dnl set $dvips to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_DVIPS],[
AC_CHECK_PROGS(dvips,dvips,no)
export dvips;
if test $dvips = "no" ;
then
	AC_MSG_ERROR([Unable to find a dvips application]);
fi;
AC_SUBST(dvips)
])
