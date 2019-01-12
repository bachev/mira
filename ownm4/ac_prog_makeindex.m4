dnl @synopsis AC_PROG_MAKEINDEX
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if makeindex is installed. If makeindex is
dnl installed, it set $makeindex to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_MAKEINDEX],[
AC_CHECK_PROGS(makeindex,makeindex,no)
export makeindex;
if test $makeindex = "no" ;
then
	AC_MSG_ERROR([Unable to find a MakeIndex application]);
fi
AC_SUBST(makeindex)
])
