dnl @synopsis AC_PROG_DVIPDF
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if dvipdf is installed. If dvipdf is installed, it
dnl set $dvipdf to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_DVIPDF],[
AC_CHECK_PROGS(dvipdf,dvipdf,no)
export dvipdf;
if test $dvipdf = "no" ;
then
	AC_MSG_ERROR([Unable to find a dvipdf application]);
fi;
AC_SUBST(dvipdf)
])
