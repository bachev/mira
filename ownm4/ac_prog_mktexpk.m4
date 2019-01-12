dnl @synopsis AC_PROG_MKTEXPK
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if mktexpk is installed. If mktexpk is installed,
dnl it set $mktexpk to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_MKTEXPK],[
AC_CHECK_PROGS(mktexpk,mktexpk,no)
export mktexpk;
if test $mktexpk = "no" ;
then
	AC_MSG_ERROR([Unable to find a mktexpk application]);
fi;
AC_SUBST(mktexpk)
])
