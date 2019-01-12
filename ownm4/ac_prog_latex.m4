dnl @synopsis AC_PROG_LATEX
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if latex is installed. If latex is installed, it
dnl set $latex to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_LATEX],[
AC_CHECK_PROGS(latex,[latex elatex lambda],no)
export latex;
if test $latex = "no" ;
then
	AC_MSG_ERROR([Unable to find a LaTeX application]);
fi
AC_SUBST(latex)
])
