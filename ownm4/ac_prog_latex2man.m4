dnl @synopsis AC_PROG_LATEX2MAN
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if latex2man is installed. If latex2man is
dnl installed, it set $latex2man to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_LATEX2MAN],[
AC_CHECK_PROGS(latex2man,[latex2man],no)
export latex2man;
if test $latex2man = "no" ;
then
	AC_MSG_ERROR([Unable to find a LaTeX2man application]);
fi
AC_SUBST(latex2man)
])
